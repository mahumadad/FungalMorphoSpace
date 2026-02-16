from pathlib import Path


def test_output_contract_exists_and_mentions_canonical_paths():
    repo = Path(__file__).resolve().parents[1]
    doc = repo / "docs" / "OUTPUT_CONTRACT.md"
    assert doc.exists()

    txt = doc.read_text(encoding="utf-8")
    # Canonical dirs
    assert "patterns/" in txt
    assert "figures/" in txt
    assert "tables/" in txt
    # Canonical files
    assert "COMPREHENSIVE_COMPARISON.png" in txt
    assert "validation_summary_machine.csv" in txt
