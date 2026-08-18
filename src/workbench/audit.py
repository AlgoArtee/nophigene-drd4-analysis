"""Tamper-evident audit-chain utilities."""

from __future__ import annotations

import hashlib
import json
import uuid
from datetime import datetime, timezone
from typing import Any

from sqlalchemy import select
from sqlalchemy.orm import Session

from .models import AuditEvent


def _canonical_event(
    *,
    sequence: int,
    occurred_at: datetime,
    event_type: str,
    actor: str,
    entity_type: str,
    entity_id: str,
    payload: dict[str, Any],
    previous_hash: str,
) -> bytes:
    timestamp = occurred_at
    if timestamp.tzinfo is None:
        timestamp = timestamp.replace(tzinfo=timezone.utc)
    timestamp_text = timestamp.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")
    return json.dumps(
        {
            "sequence": sequence,
            "occurred_at": timestamp_text,
            "event_type": event_type,
            "actor": actor,
            "entity_type": entity_type,
            "entity_id": entity_id,
            "payload": payload,
            "previous_hash": previous_hash,
        },
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        default=str,
    ).encode("utf-8")


def append_audit_event(
    session: Session,
    event_type: str,
    *,
    entity_type: str = "",
    entity_id: str = "",
    payload: dict[str, Any] | None = None,
    actor: str = "local-user",
) -> AuditEvent:
    previous = session.scalar(select(AuditEvent).order_by(AuditEvent.sequence.desc()).limit(1))
    previous_hash = previous.event_hash if previous else ""
    event = AuditEvent(
        event_type=event_type,
        actor=actor,
        entity_type=entity_type,
        entity_id=entity_id,
        payload=dict(payload or {}),
        previous_hash=previous_hash,
        # A unique placeholder lets the database assign the append-only local
        # sequence before the canonical hash is calculated.
        event_hash=hashlib.sha256(uuid.uuid4().bytes).hexdigest(),
    )
    session.add(event)
    session.flush()
    event.event_hash = hashlib.sha256(
        _canonical_event(
            sequence=event.sequence,
            occurred_at=event.occurred_at,
            event_type=event_type,
            actor=actor,
            entity_type=entity_type,
            entity_id=entity_id,
            payload=dict(payload or {}),
            previous_hash=previous_hash,
        )
    ).hexdigest()
    session.flush()
    return event


def verify_audit_chain(session: Session) -> dict[str, Any]:
    previous_hash = ""
    count = 0
    for event in session.scalars(select(AuditEvent).order_by(AuditEvent.sequence)):
        expected = hashlib.sha256(
            _canonical_event(
                sequence=event.sequence,
                occurred_at=event.occurred_at,
                event_type=event.event_type,
                actor=event.actor,
                entity_type=event.entity_type,
                entity_id=event.entity_id,
                payload=event.payload or {},
                previous_hash=previous_hash,
            )
        ).hexdigest()
        if event.previous_hash != previous_hash or event.event_hash != expected:
            return {"valid": False, "event_count": count, "failed_sequence": event.sequence}
        previous_hash = event.event_hash
        count += 1
    return {"valid": True, "event_count": count, "head_hash": previous_hash}
