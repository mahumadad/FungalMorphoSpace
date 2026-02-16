from pathlib import Path

import sys

# Ensure src is importable when running tests from repo root
REPO_ROOT = Path(__file__).resolve().parents[1]
SRC = REPO_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from fungalmorphospace.contracts.output_contract import (
    CONTRACT_VERSION,
    MACHINE_COLUMNS_CANONICAL,
    LEDGER_CORE_COLUMNS,
)


def test_output_contract_doc_mentions_version_and_schema():
    doc = REPO_ROOT / "docs" / "OUTPUT_CONTRACT.md"
    assert doc.exists()
    txt = doc.read_text(encoding="utf-8")

    assert f"**Version:** {CONTRACT_VERSION}" in txt
    assert "output_contract.py" in txt  # single source of truth

    # Machine columns must all be documented (as backticked names in the table)
    for col in MACHINE_COLUMNS_CANONICAL:
        assert f"`{col}`" in txt

    # Ledger core columns must be mentioned somewhere
    for col in LEDGER_CORE_COLUMNS:
        assert col in txt
