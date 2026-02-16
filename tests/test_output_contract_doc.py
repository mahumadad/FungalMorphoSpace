from pathlib import Path


def test_output_contract_doc_present_and_mentions_canonical_tree():
    repo = Path(__file__).resolve().parents[1]
    doc = repo / "docs" / "OUTPUT_CONTRACT.md"
    assert doc.exists(), "Missing docs/OUTPUT_CONTRACT.md"

    text = doc.read_text(encoding="utf-8")
    # Must mention canonical directories
    for required in ["patterns/", "figures/", "tables/", "metrics/", "logs/", "analysis/", "sensitivity/"]:
        assert required in text
