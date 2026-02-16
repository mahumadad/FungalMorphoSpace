#!/usr/bin/env python3
"""Parallel-safe validation CLI.

Purpose
-------
Run a *single* species and write outputs to species-specific files so that
multiple shells/processes can execute different species concurrently without
collisions.

Key outputs
-----------
- results/tables/<species_key>_validation.csv
- results/patterns/<species_key>/run_*.png (per-run analysis figures)
- results/patterns/<Full_Name_With_Underscores>.png (canonical copy of last run)

Note: Orchestration logic uses the package runner (v0.6+).
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path
import shutil
import pandas as pd

# Ensure local package is used
BASE_PATH = Path(__file__).resolve().parents[1]
SRC_PATH = BASE_PATH / "src"
if str(SRC_PATH) not in sys.path:
    sys.path.insert(0, str(SRC_PATH))

from fungalmorphospace.runners.species_database import (  # noqa: E402
    SPECIES_DATABASE,
    resolve_species_key,
    get_all_species_keys,
)
from fungalmorphospace.runners.integrated_validation import IntegratedSimulationRunner  # noqa: E402


class ParallelSafeValidator(IntegratedSimulationRunner):
    """Runner that writes species-specific summary files."""

    def run_single_species(self, species_key: str, n_runs: int = 1, verbose: bool = True):
        result = super().run_single_species(species_key=species_key, n_runs=n_runs, verbose=verbose)
        self._save_species_csv(species_key, result)
        self._refresh_canonical_pattern_copy(species_key)
        return result

    def _save_species_csv(self, species_key: str, result: dict):
        out_csv = self.paths.tables / f"{species_key}_validation.csv"

        # Backup if exists
        if out_csv.exists():
            ts = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
            shutil.copy(out_csv, out_csv.with_name(out_csv.stem + f"_{ts}" + out_csv.suffix))

        row = {
            "species_key": species_key,
            "species_name": SPECIES_DATABASE[species_key]["full_name"],
            "D_v_D_u": result.get("D_v_D_u"),
            "b": result.get("b"),
            "rho": result.get("rho"),
            "grid_size": result.get("grid_size"),
            "wavelength_best_px": result.get("wavelength_best_px", result.get("wavelength_px")),
            "wavelength_method": result.get("wavelength_method"),
            "wavelength_qc_pass": result.get("wavelength_qc_pass"),
            "spots": result.get("spots"),
            "holes": result.get("holes"),
            "euler_chi": result.get("euler_chi"),
            "density_predicted_pores_per_mm": result.get("density_predicted"),
            "dt_final": result.get("dt_final"),
            "T_target": result.get("T_target"),
            "T_actual": result.get("T_actual"),
            "steps_computed": result.get("steps_computed"),
            "pattern_path": result.get("pattern_path"),
            "timestamp": datetime.now().isoformat(),
        }

        df = pd.DataFrame([row])
        df.to_csv(out_csv, index=False)

    def _refresh_canonical_pattern_copy(self, species_key: str) -> Path | None:
        """Create/refresh canonical copy at patterns/<Full_Name>.png."""
        species_folder = self.paths.patterns / species_key
        if not species_folder.exists():
            return None

        pngs = sorted(species_folder.glob("*.png"), key=lambda p: p.stat().st_mtime)
        if not pngs:
            return None

        latest = pngs[-1]
        canonical_name = SPECIES_DATABASE[species_key]["full_name"].replace(" ", "_") + ".png"
        canonical_path = self.paths.patterns / canonical_name
        shutil.copy(latest, canonical_path)
        return canonical_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parallel-safe single-species validation",
        epilog=f"Valid species: {', '.join(SPECIES_DATABASE.keys())} (aliases also accepted)"
    )
    parser.add_argument("--species", required=True, help="Species key (or alias)")
    parser.add_argument("--n_runs", type=int, default=1)
    parser.add_argument("--output", type=str, default="results")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    # Resolve alias to canonical key
    canonical_key = resolve_species_key(args.species)
    if canonical_key is None:
        all_keys = get_all_species_keys()
        raise SystemExit(
            f"Unknown species: '{args.species}'\n"
            f"Valid options: {', '.join(SPECIES_DATABASE.keys())}\n"
            f"Aliases also accepted: {', '.join([k for k in all_keys if k not in SPECIES_DATABASE])}"
        )

    validator = ParallelSafeValidator(output_dir=args.output)
    validator.run_single_species(canonical_key, n_runs=args.n_runs, verbose=not args.quiet)

    # Print contract hints
    if not args.quiet:
        print("\nOutputs:")
        print(f"  - {Path(args.output) / 'tables' / (canonical_key + '_validation.csv')}")
        print(f"  - {Path(args.output) / 'patterns' / canonical_key}/run_*.png")
        print(f"  - {Path(args.output) / 'patterns' / (SPECIES_DATABASE[canonical_key]['full_name'].replace(' ', '_') + '.png')}")


if __name__ == "__main__":
    main()
