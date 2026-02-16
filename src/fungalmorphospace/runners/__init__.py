"""Runners (pipelines) for FungalMorphoSpace.

These modules provide importable, testable orchestration logic used by the
`/scripts` entry points.

Design goals:
- Single source of truth for species presets and calibration.
- Scripts remain thin wrappers (CLI only).
- Outputs follow the documented Output Contract.
"""

from __future__ import annotations

from .species_database import (
    SPECIES_DATABASE,
    EXCLUDED_SPECIES,
    CALIBRATION_UM_PER_PX,
    load_species_database,
    load_calibration_um_per_px,
    get_excluded_species,
)

from .integrated_validation import IntegratedSimulationRunner

__all__ = [
    "SPECIES_DATABASE",
    "EXCLUDED_SPECIES",
    "CALIBRATION_UM_PER_PX",
    "load_species_database",
    "load_calibration_um_per_px",
    "get_excluded_species",
    "IntegratedSimulationRunner",
]
