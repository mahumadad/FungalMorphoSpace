#!/usr/bin/env python3
"""
Smoke Test for FungalMorphoSpace Pipeline
==========================================
Runs minimal validation to ensure the pipeline is functional.

This test:
1. Runs --species all with small grid (256) for speed
2. Runs a single-species snapshot (brumalis) to exercise the snapshots/ledger path
3. Verifies that expected canonical output files are created
4. Verifies that per-experiment snapshots + append-only ledgers are created
5. Checks that validation_summary_machine.csv has correct columns
6. CRITICAL: Verifies that different species have different parameters
   (guards against the v0.6.0 bug where all species used defaults)

Usage:
    python3 scripts/smoke_test.py
    python3 scripts/smoke_test.py --grid 128  # faster but less accurate

Exit codes:
    0: All tests passed
    1: Tests failed
"""



from __future__ import annotations

# Prevent pytest from collecting this script as a test module.
__test__ = False

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path
import pandas as pd

# Ensure src/ is importable when running this script directly
REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_PATH = REPO_ROOT / "src"
if str(SRC_PATH) not in sys.path:
    sys.path.insert(0, str(SRC_PATH))

from fungalmorphospace.contracts.output_contract import MACHINE_COLUMNS_CANONICAL
from fungalmorphospace.utils.cleanup import cleanup_repo_temp_dirs

# Expected columns in machine-readable CSV (contract-defined)
EXPECTED_MACHINE_COLUMNS = MACHINE_COLUMNS_CANONICAL

# Expected species in output
EXPECTED_SPECIES = ["Fomes fomentarius", "Lentinus brumalis", "Polyporus squamosus"]


def run_command(cmd: list[str], cwd: Path) -> tuple[int, str, str]:
    """Run a command and return (returncode, stdout, stderr)."""
    result = subprocess.run(
        cmd,
        cwd=str(cwd),
        capture_output=True,
        text=True,
    )
    return result.returncode, result.stdout, result.stderr


def check_pipeline(output_dir: Path, grid_size: int, verbose: bool) -> bool:
    """Run the pipeline and check outputs."""

    repo_root = Path(__file__).resolve().parents[1]

    # Run integrated validation (ALL)
    cmd_all = [
        sys.executable,
        "scripts/run_integrated_validation.py",
        "--species", "all",
        "--grid", str(grid_size),
        "--n_runs", "1",
        "--T_target", "0.25",
        "--output", str(output_dir),
    ]

    if verbose:
        print(f"Running: {' '.join(cmd_all)}")

    returncode, stdout, stderr = run_command(cmd_all, repo_root)

    if verbose:
        print(stdout)
        if stderr:
            print(f"STDERR: {stderr}")

    if returncode != 0:
        print(f"❌ Pipeline (ALL) failed with return code {returncode}")
        print(stderr)
        return False

    # Run a single-species snapshot to exercise the per-run bundle + ledger path
    cmd_single = [
        sys.executable,
        "scripts/run_integrated_validation.py",
        "--species", "brumalis",
        "--grid", str(grid_size),
        "--n_runs", "1",
        "--T_target", "0.15",
        "--output", str(output_dir),
        "--quiet",
    ]

    if verbose:
        print(f"Running: {' '.join(cmd_single)}")

    returncode2, stdout2, stderr2 = run_command(cmd_single, repo_root)

    if verbose:
        print(stdout2)
        if stderr2:
            print(f"STDERR: {stderr2}")

    if returncode2 != 0:
        print(f"❌ Pipeline (single brumalis) failed with return code {returncode2}")
        print(stderr2)
        return False

    print("✓ Pipeline completed (ALL + single snapshot)")
    return True


def check_output_files(output_dir: Path) -> bool:
    """Verify expected output files exist."""
    
    expected_files = [
        output_dir / "figures" / "COMPREHENSIVE_COMPARISON.png",
        output_dir / "tables" / "validation_summary.csv",
        output_dir / "tables" / "validation_summary_machine.csv",
        output_dir / "tables" / "validation_summary.json",
    ]
    
    all_exist = True
    for f in expected_files:
        if f.exists():
            print(f"✓ Found: {f.relative_to(output_dir)}")
        else:
            print(f"❌ Missing: {f.relative_to(output_dir)}")
            all_exist = False
    
    return all_exist


def check_snapshot_outputs(output_dir: Path) -> bool:
    """Verify per-experiment snapshots and append-only ledgers exist.

    This checks both:
    - ALL snapshots: tables/all/* + figures/all/*
    - single snapshots: tables/singles/*
    """

    ok = True

    def _check_mode(mode: str, expect_figure: bool) -> bool:
        nonlocal ok
        base = output_dir / "tables" / ("all" if mode == "all" else "singles")
        ledger = base / "_index.csv"
        if not ledger.exists():
            print(f"❌ Missing ledger: {ledger.relative_to(output_dir)}")
            return False

        df = pd.read_csv(ledger)
        if df.empty or "exp_id" not in df.columns:
            print(f"❌ Ledger malformed/empty: {ledger.relative_to(output_dir)}")
            return False

        # pick last (most recent) exp_id
        exp_id = str(df.iloc[-1]["exp_id"])
        human = base / f"validation_summary_{exp_id}.csv"
        machine = base / f"validation_summary_machine_{exp_id}.csv"
        js = base / f"validation_summary_{exp_id}.json"

        for f in [human, machine, js]:
            if f.exists():
                print(f"✓ Found snapshot: {f.relative_to(output_dir)}")
            else:
                print(f"❌ Missing snapshot: {f.relative_to(output_dir)}")
                ok = False

        if expect_figure:
            fig = output_dir / "figures" / "all" / f"COMPREHENSIVE_COMPARISON_{exp_id}.png"
            if fig.exists():
                print(f"✓ Found snapshot figure: {fig.relative_to(output_dir)}")
            else:
                print(f"❌ Missing snapshot figure: {fig.relative_to(output_dir)}")
                ok = False

        # sanity: exp_id should appear in ledger (last row already), just confirm
        if exp_id not in df["exp_id"].astype(str).tolist():
            print(f"❌ exp_id not present in ledger (unexpected): {exp_id}")
            ok = False

        return True

    # ALL mode snapshots + figure
    if not _check_mode("all", expect_figure=True):
        ok = False

    # single mode snapshots (no figure copy)
    if not _check_mode("single", expect_figure=False):
        ok = False

    return ok


def check_csv_columns(output_dir: Path) -> bool:
    """Verify machine CSV has expected columns."""
    
    csv_path = output_dir / "tables" / "validation_summary_machine.csv"
    if not csv_path.exists():
        print(f"❌ Cannot check columns: {csv_path} not found")
        return False
    
    df = pd.read_csv(csv_path)
    
    missing_cols = []
    for col in EXPECTED_MACHINE_COLUMNS:
        if col not in df.columns:
            missing_cols.append(col)
    
    if missing_cols:
        print(f"❌ Missing columns: {missing_cols}")
        print(f"   Found columns: {list(df.columns)}")
        return False
    
    print(f"✓ All expected columns present ({len(EXPECTED_MACHINE_COLUMNS)} columns)")
    return True


def check_species_diversity(output_dir: Path) -> bool:
    """CRITICAL: Verify that different species have different parameters.
    
    This guards against the v0.6.0 bug where species_database.py
    failed to read nested parameters and used defaults for all species.
    """
    
    csv_path = output_dir / "tables" / "validation_summary_machine.csv"
    if not csv_path.exists():
        print(f"❌ Cannot check diversity: {csv_path} not found")
        return False
    
    df = pd.read_csv(csv_path)
    
    # Check that we have all expected species
    species_col = "species"
    if species_col not in df.columns:
        print(f"❌ No 'species' column found")
        return False
    
    found_species = df[species_col].tolist()
    for sp in EXPECTED_SPECIES:
        if sp not in found_species:
            print(f"❌ Missing species: {sp}")
            return False
    
    print(f"✓ All {len(EXPECTED_SPECIES)} species found")
    
    # CRITICAL CHECK: Verify parameters are different
    # If all D_v_D_u values are the same, something is wrong
    d_values = df["D_v_D_u"].unique()
    rho_values = df["rho"].unique()
    
    if len(d_values) == 1:
        print(f"❌ CRITICAL: All species have same D_v_D_u = {d_values[0]}")
        print("   This indicates the parameters bug is present!")
        return False
    
    if len(rho_values) == 1:
        print(f"❌ CRITICAL: All species have same rho = {rho_values[0]}")
        print("   This indicates the parameters bug is present!")
        return False
    
    print(f"✓ Parameter diversity verified:")
    print(f"   D_v_D_u values: {sorted(d_values)}")
    print(f"   rho values: {sorted(rho_values)}")
    
    # Verify wavelengths are different (strong indicator that simulation worked)
    wl_values = df["wavelength_best_px"].unique()
    if len(wl_values) == 1:
        print(f"⚠️ Warning: All species have same wavelength = {wl_values[0]}")
        print("   This is suspicious but not necessarily an error")
    else:
        print(f"✓ Wavelength diversity: {sorted(wl_values)}")
    
    return True


def main():
    parser = argparse.ArgumentParser(description="Smoke test for FungalMorphoSpace")
    parser.add_argument("--grid", type=int, default=256, help="Grid size for test (default: 256)")
    parser.add_argument("--output", type=str, default=None, help="Output directory (default: temp dir)")
    parser.add_argument("--no-clean-dev", action="store_true", help="Do not remove dev temp dirs (tmp_contract_test/tmp_check/tmp_check2) from repo root")
    parser.add_argument("--keep", action="store_true", help="Keep output directory after test")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    args = parser.parse_args()
    
    # Determine output directory
    if args.output:
        output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
        cleanup = False
    else:
        temp_dir = tempfile.mkdtemp(prefix="fms_smoke_")
        output_dir = Path(temp_dir)
        cleanup = not args.keep

    # Best-effort: remove known dev temp artifacts in repo root
    if not args.no_clean_dev:
        cleanup_repo_temp_dirs(REPO_ROOT)
    
    print("=" * 60)
    print("FungalMorphoSpace Smoke Test")
    print("=" * 60)
    print(f"Output dir: {output_dir}")
    print(f"Grid size: {args.grid}")
    print()
    
    all_passed = True
    
    # Test 1: Run pipeline
    print("TEST 1: Run pipeline")
    if not check_pipeline(output_dir, args.grid, args.verbose):
        all_passed = False
    print()
    
    # Test 2: Check output files (canonical)
    print("TEST 2: Check output files (canonical)")
    if not check_output_files(output_dir):
        all_passed = False
    print()
    
    # Test 2B: Check per-experiment snapshots + ledgers
    print("TEST 2B: Check per-experiment snapshots + ledgers")
    if not check_snapshot_outputs(output_dir):
        all_passed = False
    print()
    
    # Test 3: Check CSV columns
    print("TEST 3: Check CSV columns")
    if not check_csv_columns(output_dir):
        all_passed = False
    print()
    
    # Test 4: Check species diversity (CRITICAL)
    print("TEST 4: Check species diversity (CRITICAL)")
    if not check_species_diversity(output_dir):
        all_passed = False
    print()
    
    # Cleanup
    if cleanup:
        import shutil
        shutil.rmtree(output_dir)
        print(f"Cleaned up: {output_dir}")
    else:
        print(f"Output preserved: {output_dir}")

    # Final sweep of known dev temp dirs
    if not args.no_clean_dev:
        cleanup_repo_temp_dirs(REPO_ROOT)
    
    # Final result
    print("=" * 60)
    if all_passed:
        print("✅ ALL TESTS PASSED")
        sys.exit(0)
    else:
        print("❌ SOME TESTS FAILED")
        sys.exit(1)


if __name__ == "__main__":
    main()
