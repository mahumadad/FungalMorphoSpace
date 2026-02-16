#!/usr/bin/env python3
"""Integrated Validation CLI (wrapper).

This script is intentionally thin.
All orchestration logic lives in:
- `fungalmorphospace.runners.integrated_validation.IntegratedSimulationRunner`

Run from repo root:
  python3 scripts/run_integrated_validation.py --species all --n_runs 3
  python3 scripts/run_integrated_validation.py --species fomes
  python3 scripts/run_integrated_validation.py --species fomes_fomentarius  # alias works too

Outputs follow:
  docs/OUTPUT_CONTRACT.md
"""

from __future__ import annotations

import argparse
import sys
import datetime
import json
import hashlib
from pathlib import Path

# Ensure local package is used (avoid shadowing by any older site-packages installation)
BASE_PATH = Path(__file__).resolve().parents[1]
SRC_PATH = BASE_PATH / "src"
if str(SRC_PATH) not in sys.path:
    sys.path.insert(0, str(SRC_PATH))

from fungalmorphospace.runners.species_database import (  # noqa: E402
    SPECIES_DATABASE, 
    EXCLUDED_SPECIES,
    resolve_species_key,
    get_all_species_keys,
)
from fungalmorphospace.runners.integrated_validation import IntegratedSimulationRunner  # noqa: E402
from fungalmorphospace.utils.sensitivity_analysis import SensitivityAnalyzer  # noqa: E402



def _make_exp_id(mode: str, species: str, args_dict: dict) -> str:
    """Generate a compact, unique experiment id."""
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    payload = json.dumps({"mode": mode, "species": species, **args_dict}, sort_keys=True, default=str)
    h = hashlib.sha1(payload.encode("utf-8")).hexdigest()[:8]
    return f"{ts}_{mode}_{species}_{h}"

def _linspace(a: float, b: float, n: int):
    if n <= 1:
        return [float(a)]
    step = (b - a) / (n - 1)
    return [a + i * step for i in range(n)]


def main() -> None:
    # Build valid choices: canonical keys + "all"
    valid_species = list(SPECIES_DATABASE.keys()) + ["all"]
    
    parser = argparse.ArgumentParser(
        description="Integrated Polypore Validation (v0.7.3.post3)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python3 scripts/run_integrated_validation.py --species fomes\n"
            "  python3 scripts/run_integrated_validation.py --species all --n_runs 5\n"
            "  python3 scripts/run_integrated_validation.py --species fomes --b 1.5 --rho 0.3\n"
            "  python3 scripts/run_integrated_validation.py --sensitivity --param1 b --param2 rho --n_points 7\n"
            "\n"
            f"Valid species: {', '.join(valid_species)}\n"
            "Aliases also accepted: fomes_fomentarius, polyporus_brumalis, etc.\n"
        ),
    )

    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument(
        "--species",
        type=str,
        help="Species preset to use (short key or alias). Use 'all' for full validation.",
    )
    mode_group.add_argument("--sensitivity", action="store_true", help="Run sensitivity analysis")

    # Parameter overrides
    parser.add_argument("--D_v_D_u", type=float, help="Diffusion ratio")
    parser.add_argument("--b", type=float, help="Kinetic parameter b")
    parser.add_argument("--rho", type=float, help="Metabolic parameter rho")
    parser.add_argument("--grid", type=int, help="Grid size")
    parser.add_argument("--T_target", type=float, help="Target physical time")

    # Multiple runs
    parser.add_argument("--n_runs", type=int, default=1, help="Replicate runs for statistics")

    # Sensitivity params
    parser.add_argument("--param1", type=str, help="First parameter for sweep")
    parser.add_argument("--param2", type=str, help="Second parameter for sweep")
    parser.add_argument("--n_points", type=int, default=5, help="Points per parameter in sensitivity")    # Output
    parser.add_argument("--exp_id", type=str, default=None, help="Experiment ID (auto-generated if omitted)")
    parser.add_argument("--output", type=str, default="results", help="Output directory")
    parser.add_argument("--quiet", action="store_true", help="Minimal output")
    parser.add_argument(
        "--include_fields",
        action="store_true",
        help=(
            "Also generate per-species PNGs for the raw activator (u) and inhibitor (v) fields, "
            "annotated with species name and parameters. These are optional artifacts and do not "
            "replace the canonical COMPREHENSIVE_COMPARISON.png."
        ),
    )

    args = parser.parse_args()
    verbose = not args.quiet

    if args.sensitivity:
        if not args.param1 or not args.param2:
            raise SystemExit("--sensitivity requires --param1 and --param2")

        param1_range = _linspace(0.1, 2.0, args.n_points)
        param2_range = _linspace(0.05, 1.0, args.n_points)

        analyzer = SensitivityAnalyzer(
            base_params={
                "D_u": 1.0,
                "grid_size": args.grid or 256,
                "dx": 1.0,
                "dt": 0.0005,
                "steps": 5000,
                "model": "gierer_meinhardt",
            }
        )
        df = analyzer.sweep_2d(args.param1, param1_range, args.param2, param2_range, n_replicates=args.n_runs, verbose=verbose)

        out_dir = Path(args.output) / "sensitivity"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_csv = out_dir / f"sensitivity_{args.param1}_x_{args.param2}.csv"
        df.to_csv(out_csv, index=False)
        if verbose:
            print(f"\n✓ Sensitivity table saved: {out_csv}")
        return

    # Preset mode
    runner = IntegratedSimulationRunner(output_dir=args.output, include_field_images=args.include_fields)

    if args.species.lower() == "all":
        # Experiment id
        exp_id = args.exp_id or _make_exp_id(
            mode="all",
            species="all",
            args_dict={
                "n_runs": args.n_runs,
                "grid": args.grid,
                "T_target": args.T_target,
                "D_v_D_u": args.D_v_D_u,
                "b": args.b,
                "rho": args.rho,
            },
        )

        runner.run_all_species(n_runs=args.n_runs, grid_size=args.grid, T_target=args.T_target, verbose=verbose)
        runner.export_experiment_bundle(
            exp_id=exp_id,
            mode="all",
            include_figure=True,
            include_field_images=args.include_fields,
            include_json=True,
        )
        if verbose:
            print(f"✓ Experiment snapshot exported (exp_id={exp_id})")
        return

    # Single species - resolve alias to canonical key
    canonical_key = resolve_species_key(args.species)
    if canonical_key is None:
        all_keys = get_all_species_keys()
        raise SystemExit(
            f"Unknown species: '{args.species}'\n"
            f"Valid options: {', '.join(list(SPECIES_DATABASE.keys()) + ['all'])}\n"
            f"Aliases also accepted: {', '.join([k for k in all_keys if k not in SPECIES_DATABASE])}"
        )

    # Experiment id
    exp_id = args.exp_id or _make_exp_id(
        mode="single",
        species=canonical_key,
        args_dict={
            "n_runs": args.n_runs,
            "grid": args.grid,
            "T_target": args.T_target,
            "D_v_D_u": args.D_v_D_u,
            "b": args.b,
            "rho": args.rho,
        },
    )

    result = runner.run_single_species(
        species_key=canonical_key,
        D_v_D_u=args.D_v_D_u,
        b=args.b,
        rho=args.rho,
        grid_size=args.grid,
        T_target=args.T_target,
        n_runs=args.n_runs,
        verbose=verbose,
    )
    runner.export_experiment_bundle(
        exp_id=exp_id,
        mode="single",
        include_figure=False,
        include_field_images=args.include_fields,
        include_json=True,
    )
    if verbose:
        print(f"✓ Experiment snapshot exported (exp_id={exp_id})")

    if verbose:
        print("\nExcluded species (documented):")
        for k, info in EXCLUDED_SPECIES.items():
            print(f"  - {k}: {info.get('full_name')} :: {info.get('exclusion_reason')}")


if __name__ == "__main__":
    main()
