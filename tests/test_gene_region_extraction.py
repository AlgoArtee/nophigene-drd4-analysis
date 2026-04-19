from __future__ import annotations

from typing import Any

from src.gene_region_extraction import fetch_ucsc_region


class FakeResponse:
    def __init__(self, payload: dict[str, Any]) -> None:
        self.payload = payload

    def raise_for_status(self) -> None:
        return None

    def json(self) -> dict[str, Any]:
        return self.payload


def test_ucsc_lookup_uses_search_endpoint_for_any_chromosome(monkeypatch) -> None:
    """UCSC lookup should not be locked to the original chr11 DRD4 prototype."""

    def fake_get(url: str, *, params: dict[str, str], timeout: int) -> FakeResponse:
        assert url == "https://api.genome.ucsc.edu/search"
        assert params == {"search": "SIRT6", "genome": "hg19"}
        assert timeout > 0
        return FakeResponse(
            {
                "positionMatches": [
                    {
                        "trackName": "knownGene",
                        "matches": [
                            {
                                "position": "chr19:4174106-4182560",
                                "posName": "SIRT6 (ENST00000337491.6)",
                            },
                            {
                                "position": "chr11:637269-640706",
                                "posName": "DRD4 (ENST00000176183.6)",
                            },
                        ],
                    }
                ]
            }
        )

    monkeypatch.setattr("src.gene_region_extraction.requests.get", fake_get)

    assert fetch_ucsc_region("sirt6") == "19:4174106-4182560"


def test_ucsc_lookup_prefers_exact_symbol_matches(monkeypatch) -> None:
    """Search hits for interacting or similarly named genes should not leak in."""

    def fake_get(url: str, *, params: dict[str, str], timeout: int) -> FakeResponse:
        return FakeResponse(
            {
                "positionMatches": [
                    {
                        "trackName": "knownGene",
                        "matches": [
                            {
                                "position": "chrX:101906411-101914011",
                                "posName": "GPRASP1 (DRD4-associated sorting protein)",
                            },
                            {
                                "position": "chr11:637269-640706",
                                "posName": "DRD4 (ENST00000176183.6)",
                            },
                        ],
                    }
                ]
            }
        )

    monkeypatch.setattr("src.gene_region_extraction.requests.get", fake_get)

    assert fetch_ucsc_region("DRD4") == "11:637269-640706"
