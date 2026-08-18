"""Typed multi-source gene interaction graph construction."""

from __future__ import annotations

from collections import defaultdict, deque
from typing import Any, Iterable

MAX_HOPS = 3
DEFAULT_NODE_CAP = 150
SOURCE_TIERS = {
    "intact": 5,
    "biogrid": 5,
    "ensembl_interactions": 5,
    "reactome": 4,
    "unibind": 4,
    "encode": 4,
    "screen": 4,
    "gtex": 3,
    "string": 2,
}


def _rank(edge: dict[str, Any], corroboration: int) -> tuple[float, float, str]:
    source = str(edge.get("source_key") or "").casefold()
    source_tier = SOURCE_TIERS.get(source, 1)
    native = edge.get("native_score")
    native_score = float(native) if native is not None else -1.0
    return float(source_tier + min(corroboration, 4)), native_score, source


def build_interaction_graph(
    query_gene: str,
    edges: Iterable[dict[str, Any]],
    *,
    max_hops: int = 1,
    node_cap: int = DEFAULT_NODE_CAP,
    per_hop_caps: tuple[int, int, int] = (30, 60, 59),
) -> dict[str, Any]:
    query = str(query_gene or "").upper()
    hops = max(1, min(MAX_HOPS, int(max_hops)))
    cap = max(2, min(DEFAULT_NODE_CAP, int(node_cap)))
    normalized: list[dict[str, Any]] = []
    pair_sources: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    adjacency: dict[str, list[int]] = defaultdict(list)
    for raw in edges:
        edge = dict(raw)
        source = str(edge.get("source_gene") or "").upper()
        target = str(edge.get("target_gene") or "").upper()
        if not source or not target or source == target:
            continue
        edge["source_gene"] = source
        edge["target_gene"] = target
        edge["edge_type"] = str(edge.get("edge_type") or "functional_association")
        edge["source_key"] = str(edge.get("source_key") or "unknown")
        index = len(normalized)
        normalized.append(edge)
        adjacency[source].append(index)
        if not edge.get("directed"):
            adjacency[target].append(index)
        pair_sources[(source, target, edge["edge_type"])].add(edge["source_key"])

    selected_nodes = {query}
    selected_edges: list[dict[str, Any]] = []
    visited_edges: set[int] = set()
    queue = deque([(query, 0)])
    hop_counts = defaultdict(int)
    while queue and len(selected_nodes) < cap:
        node, depth = queue.popleft()
        if depth >= hops:
            continue
        candidates: list[tuple[tuple[float, float, str], int, str]] = []
        for index in adjacency.get(node, []):
            if index in visited_edges:
                continue
            edge = normalized[index]
            neighbor = edge["target_gene"] if edge["source_gene"] == node else edge["source_gene"]
            if edge.get("directed") and edge["source_gene"] != node:
                continue
            corroboration = len(pair_sources[(edge["source_gene"], edge["target_gene"], edge["edge_type"])])
            candidates.append((_rank(edge, corroboration), index, neighbor))
        candidates.sort(reverse=True)
        hop_limit = min(per_hop_caps[depth], cap - len(selected_nodes))
        added_at_hop = 0
        for _score, index, neighbor in candidates:
            if added_at_hop >= hop_limit:
                break
            edge = dict(normalized[index])
            visited_edges.add(index)
            edge["hop"] = depth + 1
            edge["corroborating_source_count"] = len(
                pair_sources[(edge["source_gene"], edge["target_gene"], edge["edge_type"])]
            )
            edge["combined_rank_policy"] = "source tier + independent corroboration; not a probability"
            selected_edges.append(edge)
            if neighbor not in selected_nodes:
                selected_nodes.add(neighbor)
                queue.append((neighbor, depth + 1))
                hop_counts[depth + 1] += 1
                added_at_hop += 1
            if len(selected_nodes) >= cap:
                break
    return {
        "query_gene": query,
        "max_hops": hops,
        "node_cap": cap,
        "node_count": len(selected_nodes),
        "edge_count": len(selected_edges),
        "hop_counts": dict(hop_counts),
        "nodes": [{"id": gene, "type": "gene", "query": gene == query} for gene in sorted(selected_nodes)],
        "edges": selected_edges,
        "truncated": len(selected_nodes) >= cap,
        "hop_semantics": "gene_to_gene",
    }
