"""Output contract definitions for FungalMorphoSpace canonical artifacts.

This module is the **single source of truth** for:

- Canonical output file names and directory layout.
- Per-experiment snapshot naming conventions.
- Stable machine-readable CSV schema (column names, ordering, and types).
- Ledger core schema for the append-only experiment index.

Any change to the output contract must be accompanied by a version bump in
all three locations:

1. :data:`CONTRACT_VERSION` (this file).
2. Package version in ``pyproject.toml`` and ``fungalmorphospace.__version__``.
3. Documentation in ``docs/OUTPUT_CONTRACT.md``.

The contract guarantees downstream consumers (CI pipelines, notebooks,
dashboards) a stable interface even as the simulation engine evolves.

See Also:
    :func:`enforce_machine_schema` for runtime schema enforcement.
    :func:`validate_machine_schema` for non-mutating schema validation.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Any

import pandas as pd


# ---------------------------------------------------------------------------
# Contract version -- keep aligned with docs/OUTPUT_CONTRACT.md
# ---------------------------------------------------------------------------

#: Semantic version string for the output contract.  Bumped whenever
#: canonical file names, column schemas, or directory layout change.
CONTRACT_VERSION: str = "0.7.3"


# ---------------------------------------------------------------------------
# Canonical output paths (relative to the experiment root directory)
# ---------------------------------------------------------------------------

#: Mapping of artifact roles to their canonical relative paths.
#: These files are always present after a successful ``run_all_species`` call.
CANONICAL_FILES: dict[str, Path] = {
    "figure": Path("figures") / "COMPREHENSIVE_COMPARISON.png",
    "human_csv": Path("tables") / "validation_summary.csv",
    "machine_csv": Path("tables") / "validation_summary_machine.csv",
    "json": Path("tables") / "validation_summary.json",
}

#: Filename used for the append-only experiment ledger in each snapshot
#: subdirectory.
LEDGER_FILENAME: str = "_index.csv"

#: Per-experiment snapshot paths, parameterized by ``{exp_id}``.
#: Two modes are supported: ``"single"`` (one species) and ``"all"``
#: (full multi-species campaign).
SNAPSHOT_FILES: dict[str, dict[str, Path]] = {
    "single": {
        "human_csv": Path("tables") / "singles" / "validation_summary_{exp_id}.csv",
        "machine_csv": Path("tables") / "singles" / "validation_summary_machine_{exp_id}.csv",
        "json": Path("tables") / "singles" / "validation_summary_{exp_id}.json",
        "ledger": Path("tables") / "singles" / LEDGER_FILENAME,
    },
    "all": {
        "human_csv": Path("tables") / "all" / "validation_summary_{exp_id}.csv",
        "machine_csv": Path("tables") / "all" / "validation_summary_machine_{exp_id}.csv",
        "json": Path("tables") / "all" / "validation_summary_{exp_id}.json",
        "figure": Path("figures") / "all" / "COMPREHENSIVE_COMPARISON_{exp_id}.png",
        "ledger": Path("tables") / "all" / LEDGER_FILENAME,
    },
}

#: Legacy root-level copies maintained for backward compatibility with
#: older notebooks and scripts that expect artifacts at the project root.
LEGACY_COPIES: dict[str, Path] = {
    "figure": Path("COMPREHENSIVE_COMPARISON.png"),
    "human_csv": Path("validation_summary.csv"),
}


# ---------------------------------------------------------------------------
# Machine-readable CSV schema
# ---------------------------------------------------------------------------

#: Ordered column names for the canonical machine-readable CSV
#: (``validation_summary_machine.csv``).  The schema is strict: no extra
#: columns are permitted, and missing columns are filled with ``pd.NA``.
MACHINE_COLUMNS_CANONICAL: list[str] = [
    "species",
    "species_key",
    "D_v_D_u",
    "b",
    "rho",
    "grid_size",
    "wavelength_best_px",
    "wavelength_method",
    "wavelength_qc_pass",
    "wavelength_fft_px",
    "wavelength_fft_qc_pass",
    "wavelength_fft_peak_ratio",
    "wavelength_autocorr_px",
    "wavelength_autocorr_qc_pass",
    "spots",
    "holes",
    "euler_chi",
    "density_predicted_pores_per_mm",
    "density_observed_min_pores_per_mm",
    "density_observed_max_pores_per_mm",
    "dt_requested",
    "dt_final",
    "dt_adjusted",
    "dt_max_diffusion",
    "safety_factor",
    "T_target",
    "T_actual",
    "steps_computed",
    "had_warning",
    "pattern_path",
]

#: Snapshot machine CSVs prepend ``exp_id`` to the canonical column list,
#: keeping the rest identical and identically ordered.
MACHINE_COLUMNS_SNAPSHOT: list[str] = ["exp_id"] + MACHINE_COLUMNS_CANONICAL

#: Coarse type expectations for each machine CSV column.  Used by
#: :func:`validate_machine_schema` for best-effort type checking.
#: Values are one of ``"str"``, ``"float"``, ``"int"``, or ``"bool"``.
MACHINE_COLUMN_TYPES: dict[str, str] = {
    "species": "str",
    "species_key": "str",
    "D_v_D_u": "float",
    "b": "float",
    "rho": "float",
    "grid_size": "int",
    "wavelength_best_px": "float",
    "wavelength_method": "str",
    "wavelength_qc_pass": "bool",
    "wavelength_fft_px": "float",
    "wavelength_fft_qc_pass": "bool",
    "wavelength_fft_peak_ratio": "float",
    "wavelength_autocorr_px": "float",
    "wavelength_autocorr_qc_pass": "bool",
    "spots": "int",
    "holes": "int",
    "euler_chi": "int",
    "density_predicted_pores_per_mm": "float",
    "density_observed_min_pores_per_mm": "float",
    "density_observed_max_pores_per_mm": "float",
    "dt_requested": "float",
    "dt_final": "float",
    "dt_adjusted": "bool",
    "dt_max_diffusion": "float",
    "safety_factor": "float",
    "T_target": "float",
    "T_actual": "float",
    "steps_computed": "int",
    "had_warning": "bool",
    "pattern_path": "str",
    "exp_id": "str",
}

#: Core columns that must appear in every ledger (``_index.csv``) row.
#: Additional per-species wide columns may be appended dynamically via
#: union schema expansion.
LEDGER_CORE_COLUMNS: list[str] = ["exp_id", "mode", "generated_at", "version", "n_species", "species_json"]


def enforce_machine_schema(
    df: pd.DataFrame,
    *,
    snapshot: bool,
) -> pd.DataFrame:
    """Enforce the canonical machine CSV schema on a DataFrame.

    Ensures that the resulting DataFrame has exactly the expected columns
    in the correct order.  Missing columns are added as ``pd.NA``;
    unexpected extra columns are dropped.

    This function is the primary mechanism for maintaining schema stability
    across pipeline versions.

    Args:
        df: Input DataFrame, potentially with missing or extra columns.
        snapshot: If ``True``, use :data:`MACHINE_COLUMNS_SNAPSHOT`
            (includes ``exp_id``); otherwise use
            :data:`MACHINE_COLUMNS_CANONICAL`.

    Returns:
        A new DataFrame with exactly the contracted columns in the
        contracted order.
    """
    cols: list[str] = MACHINE_COLUMNS_SNAPSHOT if snapshot else MACHINE_COLUMNS_CANONICAL

    # Add missing columns as NA.
    for c in cols:
        if c not in df.columns:
            df[c] = pd.NA

    # Drop any columns not in the contract.
    extras: list[str] = [c for c in df.columns if c not in cols]
    if extras:
        df = df.drop(columns=extras)

    # Reorder to canonical column sequence.
    df = df[cols]
    return df


def validate_machine_schema(
    df: pd.DataFrame,
    *,
    snapshot: bool,
) -> tuple[bool, list[str]]:
    """Validate a DataFrame against the machine CSV schema without mutation.

    Checks for missing columns, unexpected extra columns, and performs
    best-effort type validation (numeric convertibility for float/int
    columns, boolean value range for bool columns).

    Args:
        df: DataFrame to validate.
        snapshot: If ``True``, validate against the snapshot schema;
            otherwise validate against the canonical schema.

    Returns:
        A 2-tuple ``(ok, messages)`` where *ok* is ``True`` if no issues
        were found, and *messages* is a list of human-readable diagnostic
        strings describing each detected violation.
    """
    msgs: list[str] = []
    cols: list[str] = MACHINE_COLUMNS_SNAPSHOT if snapshot else MACHINE_COLUMNS_CANONICAL

    missing: list[str] = [c for c in cols if c not in df.columns]
    extra: list[str] = [c for c in df.columns if c not in cols]

    if missing:
        msgs.append(f"Missing columns: {missing}")
    if extra:
        msgs.append(f"Unexpected columns: {extra}")

    # Best-effort type sanity checking.
    for c in cols:
        if c not in df.columns:
            continue
        expected: str | None = MACHINE_COLUMN_TYPES.get(c)
        if expected in {"float", "int"}:
            # Allow NA; check numeric convertibility for non-NA values.
            s = df[c].dropna()
            if not s.empty:
                try:
                    pd.to_numeric(s, errors="raise")
                except Exception:
                    msgs.append(f"Column {c} not numeric-convertible")
        elif expected == "bool":
            # Accept True/False/0/1/"True"/"False" (case-insensitive).
            s = df[c].dropna().astype(str).str.lower().unique().tolist()
            allowed: set[str] = {"true", "false", "0", "1"}
            if any(v not in allowed for v in s):
                msgs.append(f"Column {c} has non-bool values: {s[:8]}")

    return (len(msgs) == 0), msgs
