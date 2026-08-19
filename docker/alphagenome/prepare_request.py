"""Canonical, side-effect-free AlphaGenome request validation.

This is the callable worker form of the reviewed AlphaGenome skill helper.  It
does not contact the provider and it deliberately leaves reference validation
to the independently mounted GRCh38 FASTA check in ``worker.py``.
"""

from __future__ import annotations

import re
from typing import Any, Iterable


VARIANT = re.compile(
    r"^(chr(?:[1-9]|1\d|2[0-2]|X|Y)):(\d+):([ACGTN]+)>([ACGTN]+)$",
    re.IGNORECASE,
)
MODALITIES = {
    "RNA_SEQ",
    "DNASE",
    "ATAC",
    "CAGE",
    "CHIP_HISTONE",
    "CHIP_TF",
    "SPLICE_SITES",
    "SPLICE_JUNCTIONS",
}


def prepare_request(
    *,
    variant: str,
    assembly: str,
    modalities: Iterable[str],
    ontology_terms: Iterable[str],
    sequence_length: int,
    external_transfer_acknowledged: bool,
) -> dict[str, Any]:
    """Validate and normalize one request without starting inference."""
    if assembly != "GRCh38":
        raise ValueError("AlphaGenome requests require GRCh38.")
    match = VARIANT.fullmatch(str(variant).strip())
    if not match:
        raise ValueError("variant must look like chr16:27344882:G>A")
    chromosome, position_text, reference, alternate = match.groups()
    position = int(position_text)
    if position < 1 or reference.upper() == alternate.upper():
        raise ValueError("position must be positive and REF must differ from ALT")
    requested_modalities = list(dict.fromkeys(str(item).upper() for item in modalities))
    unsupported = sorted(set(requested_modalities) - MODALITIES)
    if not requested_modalities or unsupported:
        raise ValueError(f"unsupported or missing modalities: {', '.join(unsupported)}")
    length = int(sequence_length)
    if length < 65_536 or length > 1_048_576 or length & (length - 1):
        raise ValueError("sequence length must be a power of two from 65,536 through 1,048,576")
    if not external_transfer_acknowledged:
        raise ValueError("explicit external-transfer acknowledgement is required")
    center_zero = position - 1
    start = max(0, center_zero - length // 2)
    return {
        "assembly": assembly,
        "variant": {
            "chromosome": chromosome,
            "position_1_based": position,
            "reference": reference.upper(),
            "alternate": alternate.upper(),
        },
        "model_interval_0_based_half_open": {
            "start": start,
            "end": start + length,
            "length": length,
        },
        "modalities": requested_modalities,
        "ontology_terms": list(dict.fromkeys(str(item) for item in ontology_terms if str(item))),
        "external_transfer_acknowledged": True,
        "reference_allele_verified": False,
        "inference_started": False,
    }
