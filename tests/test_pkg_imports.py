import sys
from pathlib import Path


def test_package_imports_and_species_db():
    repo = Path(__file__).resolve().parents[1]
    src = repo / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))

    from fungalmorphospace.runners import SPECIES_DATABASE, CALIBRATION_UM_PER_PX, IntegratedSimulationRunner

    assert isinstance(SPECIES_DATABASE, dict)
    assert len(SPECIES_DATABASE) >= 3
    assert CALIBRATION_UM_PER_PX > 0
    assert IntegratedSimulationRunner is not None
