"""Canonical species presets and microscopy calibration for FungalMorphoSpace.

This module is the **single source of truth** for fungal-species simulation
parameters and the pixel-to-micrometer calibration factor used throughout the
pipeline.  All data are loaded at import time from ``data/species_data.json``
(located at the repository root) and exposed as module-level constants.

Key design decisions:
    - Parameters live in a nested ``"parameters"`` object inside each species
      entry in the JSON file (see **CRITICAL FIX v0.6.2** below).
    - Species are keyed by short canonical identifiers (e.g. ``"fomes"``,
      ``"brumalis"``, ``"squamosus"``).  Aliases and full scientific names are
      resolved by :func:`resolve_species_key`.
    - Species marked ``"role": "EXCLUDED"`` in the JSON are omitted from the
      active database but are available via :func:`get_excluded_species` for
      documentation and audit purposes.

Module-level constants:
    ``CALIBRATION_UM_PER_PX``
        Micrometer-per-pixel scale factor (float).
    ``SPECIES_DATABASE``
        ``Dict[str, Dict[str, Any]]`` of active species presets.
    ``EXCLUDED_SPECIES``
        ``Dict[str, Dict[str, Any]]`` of excluded species with rationale.

See Also:
    ``docs/OUTPUT_CONTRACT.md`` for the downstream contract that consumes
    these parameters.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Any, Optional


def _find_repo_root(start: Path | None = None) -> Path:
    """Locate the repository root by walking up the directory tree.

    The root is identified as the first ancestor directory that contains
    both ``pyproject.toml`` and a ``data/`` subdirectory.  If no match is
    found, a fallback three levels above ``__file__`` is returned (which
    corresponds to the expected package layout).

    Args:
        start: Starting path for the upward search.  Defaults to the
            path of this module file.

    Returns:
        Absolute ``Path`` to the inferred repository root.
    """
    here: Path = (start or Path(__file__)).resolve()
    for p in [here] + list(here.parents):
        if (p / "pyproject.toml").exists() and (p / "data").exists():
            return p
    # Fallback: assume standard package nesting (src/fungalmorphospace/runners/).
    return Path(__file__).resolve().parents[3]


def _species_json_path() -> Path:
    """Return the absolute path to ``data/species_data.json``.

    Returns:
        ``Path`` object pointing to the canonical species data file.
    """
    root: Path = _find_repo_root()
    return root / "data" / "species_data.json"


def load_calibration_um_per_px() -> float:
    """Load the micrometer-per-pixel scale factor from ``species_data.json``.

    The calibration value is read from the ``"calibration"`` section of the
    JSON file.  If the key is missing, a default of ``400 / 46`` um/px is
    used (derived from the original microscope calibration slide).

    Returns:
        Scale factor in micrometers per pixel.
    """
    json_path: Path = _species_json_path()
    with open(json_path, "r", encoding="utf-8") as f:
        data: dict = json.load(f)
    cal: dict = data.get("calibration", {}) or {}
    return float(cal.get("scale_factor_um_per_px", 400 / 46))


def load_species_database() -> Dict[str, Dict[str, Any]]:
    """Load active species presets from ``species_data.json``.

    Returns a dictionary keyed by canonical short key (e.g. ``"fomes"``),
    where each value is a flat parameter dictionary suitable for direct
    consumption by :class:`IntegratedSimulationRunner`.

    **CRITICAL FIX (v0.6.2):** Simulation parameters (``D_v_D_u``, ``b``,
    ``rho``, etc.) are now read from the nested ``"parameters"`` object
    inside each species entry, **not** from the species top-level.  Earlier
    versions silently used default values because the top-level keys did
    not exist.

    Returns:
        Dictionary mapping canonical species keys to parameter dictionaries.
        Each parameter dictionary contains:

        - ``full_name`` (``str``): Scientific binomial name.
        - ``common_name`` (``str``): Vernacular name.
        - ``aliases`` (``list[str]``): Alternative lookup keys.
        - ``D_v_D_u`` (``float``): Inhibitor/activator diffusion ratio.
        - ``b`` (``float``): Gierer-Meinhardt saturation parameter.
        - ``rho`` (``float``): Source-density parameter.
        - ``a`` (``float``): Gierer-Meinhardt activation rate.
        - ``grid_size`` (``int``): Spatial grid side length (pixels).
        - ``T_target`` (``float``): Target simulation time.
        - ``dt_initial`` (``float``): Requested initial time step.
        - ``density_observed`` (``tuple[float, float]``): Observed pore
          density range (min, max) in pores/mm.
        - ``expected_wavelength_px`` (``float``): Expected pattern
          wavelength in pixels.
        - ``pore_spacing_um`` (``float``): Expected pore spacing in
          micrometers.
    """
    json_path: Path = _species_json_path()
    with open(json_path, "r", encoding="utf-8") as f:
        data: dict = json.load(f)

    species_db: Dict[str, Dict[str, Any]] = {}

    for species_key, species_data in data.get("species", {}).items():
        # Skip species explicitly marked as excluded from simulation.
        if species_data.get("role") == "EXCLUDED":
            continue

        # CRITICAL: Read from nested 'parameters' object (v0.6.2 fix).
        params: dict = species_data.get("parameters", {})

        species_db[species_key] = {
            "full_name": species_data.get("scientific_name", species_key),
            "common_name": species_data.get("common_name", ""),
            "aliases": species_data.get("aliases", []),
            # Parameters from nested object (THE FIX)
            "D_v_D_u": float(params.get("D_v_D_u", 150.0)),
            "b": float(params.get("b", 1.0)),
            "rho": float(params.get("rho", 0.2)),
            "a": float(params.get("a", 0.1)),
            "grid_size": int(params.get("grid_size", 512)),
            "T_target": float(params.get("T_target", 5.0)),
            "dt_initial": float(params.get("dt_initial", 0.0005)),
            # Observational data from species level
            "density_observed": tuple(species_data.get("pore_density_per_mm", [0, 0])),
            "expected_wavelength_px": float(species_data.get("expected_wavelength_px", 0.0)),
            "pore_spacing_um": float(species_data.get("pore_spacing_um", 0.0)),
        }

    return species_db


def resolve_species_key(key: str) -> Optional[str]:
    """Resolve an alias or scientific name to its canonical species key.

    Accepts any of the following forms (case-insensitive, stripped):

    - Canonical key: ``"fomes"``
    - Underscore form: ``"fomes_fomentarius"``
    - Scientific name: ``"Fomes fomentarius"``
    - Any string listed in the species ``"aliases"`` array.

    Args:
        key: Species identifier to resolve.

    Returns:
        The canonical short key (e.g. ``"fomes"``) if a match is found,
        or ``None`` if the key cannot be resolved.
    """
    key_lower: str = key.lower().strip()

    # Direct match against canonical keys.
    if key_lower in SPECIES_DATABASE:
        return key_lower

    # Search aliases and scientific names in every species entry.
    for canonical, data in SPECIES_DATABASE.items():
        aliases: list[str] = [a.lower() for a in data.get("aliases", [])]
        if key_lower in aliases:
            return canonical
        # Also check scientific name.
        if key_lower == data.get("full_name", "").lower():
            return canonical

    return None


def get_all_species_keys() -> list[str]:
    """Return a list of all valid species lookup keys.

    Includes both canonical short keys and all registered aliases.
    Useful for CLI tab-completion or input validation.

    Returns:
        Flat list of strings that :func:`resolve_species_key` would
        successfully resolve.
    """
    keys: list[str] = list(SPECIES_DATABASE.keys())
    for data in SPECIES_DATABASE.values():
        keys.extend(data.get("aliases", []))
    return keys


def get_excluded_species() -> Dict[str, Dict[str, Any]]:
    """Return species that are excluded from simulation, with rationale.

    Excluded species have ``"role": "EXCLUDED"`` in the JSON.  They are
    kept in the data file for documentation and reproducibility but are
    not loaded into ``SPECIES_DATABASE``.

    Returns:
        Dictionary mapping species keys to dicts containing:

        - ``full_name`` (``str``): Scientific binomial name.
        - ``exclusion_reason`` (``str``): Why the species was excluded.
        - ``density_observed`` (``tuple[float, float]``): Observed pore
          density range (min, max) in pores/mm.
    """
    json_path: Path = _species_json_path()
    with open(json_path, "r", encoding="utf-8") as f:
        data: dict = json.load(f)

    excluded: Dict[str, Dict[str, Any]] = {}
    for species_key, species_data in data.get("species", {}).items():
        if species_data.get("role") == "EXCLUDED":
            excluded[species_key] = {
                "full_name": species_data.get("scientific_name", species_key),
                "exclusion_reason": species_data.get("exclusion_rationale", ""),
                "density_observed": tuple(species_data.get("pore_density_per_mm", [0, 0])),
            }
    return excluded


# ---------------------------------------------------------------------------
# Module-level constants (loaded once at import time)
# ---------------------------------------------------------------------------

#: Micrometer-per-pixel calibration factor used for density predictions.
CALIBRATION_UM_PER_PX: float = load_calibration_um_per_px()

#: Active species presets keyed by canonical short identifier.
SPECIES_DATABASE: Dict[str, Dict[str, Any]] = load_species_database()

#: Excluded species with exclusion rationale (for audit/documentation).
EXCLUDED_SPECIES: Dict[str, Dict[str, Any]] = get_excluded_species()
