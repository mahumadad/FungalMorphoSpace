"""Integrated validation runner for the FungalMorphoSpace simulation pipeline.

This module provides the main orchestration logic for running Turing
reaction-diffusion simulations across multiple fungal species, collecting
topology metrics, and producing publication-ready comparison reports.

Originally extracted from ``scripts/run_integrated_validation.py`` to enable
programmatic (library) access.

Key guarantees (v0.6+):
    - No cross-imports between scripts.
    - Canonical outputs follow ``docs/OUTPUT_CONTRACT.md``.
    - Temporal comparability across species via ``T_target`` with
      ``steps = ceil(T_target / dt_final)``.
    - Wavelength selection is robust:
      FFT with QC --> autocorrelation with QC --> finite-positive fallback.

Typical usage::

    runner = IntegratedSimulationRunner("results")
    runner.run_all_species(n_runs=3, verbose=True)
    runner.create_comparison_report()

See Also:
    ``fungalmorphospace.contracts.output_contract`` for the canonical CSV
    schema enforced by :func:`enforce_machine_schema`.
"""

from __future__ import annotations

import json
import hashlib
import logging
import math
import shutil
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional, Tuple

import numpy as np
import pandas as pd

from fungalmorphospace.contracts.output_contract import enforce_machine_schema

from ..analysis.topology_analyzer import TopologyAnalyzer
from ..analysis.visualization import EnhancedVisualizer
from ..core.kinetics import create_kinetics
from ..core.turing_simulator import TuringSimulator
from fungalmorphospace import __version__
from .species_database import SPECIES_DATABASE, CALIBRATION_UM_PER_PX


def select_wavelength_px(
    metrics: Dict[str, Any],
) -> Tuple[float, str, bool, str]:
    """Select the best wavelength estimator using a QC cascade.

    Implements a three-tier quality-control cascade to choose the most
    reliable spatial wavelength estimate from a topology-metrics dictionary.

    Cascade order:
        1. **FFT radial peak** -- preferred when its QC flag is ``True``
           and the value is finite and positive.
        2. **Autocorrelation radial** -- used when FFT fails QC but
           autocorrelation passes.
        3. **Fallback** -- first finite positive among ``[FFT, autocorr]``
           is returned with ``qc_pass=False`` so downstream code can flag
           the estimate as unreliable.

    Args:
        metrics: Dictionary of topology/spectral metrics as returned by
            :meth:`TopologyAnalyzer.compute_all_metrics`.  Expected keys
            include ``wavelength_fft``, ``wavelength_fft_qc_pass``,
            ``wavelength_autocorr``, and ``wavelength_autocorr_qc_pass``.

    Returns:
        A 4-tuple ``(wavelength_px, method, qc_pass, qc_reason)`` where

        - *wavelength_px* (``float``): best wavelength in pixels (0.0 if
          no valid estimate exists).
        - *method* (``str``): one of ``"fft"``, ``"autocorr"``, or
          ``"none"``.
        - *qc_pass* (``bool``): whether the returned estimate passed its
          own quality-control check.
        - *qc_reason* (``str``): human-readable rationale for the
          selection decision.
    """
    # --- Tier 1: FFT radial peak (preferred) ---
    fft = metrics.get("wavelength_fft", np.nan)
    fft_ok = bool(metrics.get("wavelength_fft_qc_pass", False)) and np.isfinite(fft) and fft > 0

    # --- Tier 2: Autocorrelation radial ---
    ac = metrics.get("wavelength_autocorr", np.nan)
    ac_ok = bool(metrics.get("wavelength_autocorr_qc_pass", False)) and np.isfinite(ac) and ac > 0

    if fft_ok:
        return float(fft), "fft", True, "FFT QC pass"
    if ac_ok:
        return float(ac), "autocorr", True, "Autocorr QC pass (FFT failed)"

    # --- Tier 3: Fallback (QC failed) ---
    # Keep something finite if available for debugging, but mark QC fail.
    for v, m in [(fft, "fft"), (ac, "autocorr")]:
        if np.isfinite(v) and v > 0:
            return float(v), m, False, "QC fail; fallback used"

    return 0.0, "none", False, "No valid wavelength estimate"


@dataclass
class OutputPaths:
    """Canonical directory structure for simulation outputs.

    Each field maps to a subdirectory under the experiment root.  Directories
    are created eagerly by :func:`_ensure_dirs` at runner initialization.

    Attributes:
        root: Top-level output directory (e.g. ``results/``).
        patterns: Per-species pattern PNG images.
        figures: Publication figures (e.g. ``COMPREHENSIVE_COMPARISON.png``).
        tables: CSV / JSON summary tables (human and machine-readable).
        metrics: Raw numerical metrics dumps.
        logs: Error logs and diagnostic traces.
        analysis: Extended analysis artifacts.
        sensitivity: Sensitivity-sweep outputs.
    """

    root: Path
    patterns: Path
    figures: Path
    tables: Path
    metrics: Path
    logs: Path
    analysis: Path
    sensitivity: Path


def _ensure_dirs(output_dir: Path) -> OutputPaths:
    """Create the canonical output directory tree and return an ``OutputPaths``.

    All directories are created with ``parents=True, exist_ok=True`` so the
    function is idempotent and safe to call multiple times.

    Args:
        output_dir: Root directory under which subdirectories will be created.

    Returns:
        An :class:`OutputPaths` dataclass whose fields point to the freshly
        ensured subdirectories.
    """
    output_dir.mkdir(exist_ok=True, parents=True)

    paths = OutputPaths(
        root=output_dir,
        patterns=output_dir / "patterns",
        figures=output_dir / "figures",
        tables=output_dir / "tables",
        metrics=output_dir / "metrics",
        logs=output_dir / "logs",
        analysis=output_dir / "analysis",
        sensitivity=output_dir / "sensitivity",
    )

    # Materialize every subdirectory on disk.
    for p in [
        paths.patterns,
        paths.figures,
        paths.tables,
        paths.metrics,
        paths.logs,
        paths.analysis,
        paths.sensitivity,
    ]:
        p.mkdir(exist_ok=True, parents=True)

    return paths


class IntegratedSimulationRunner:
    """Orchestrator for multi-species Turing-pattern simulation campaigns.

    Manages parameter resolution from ``SPECIES_DATABASE``, replicate
    execution, topology analysis, statistical aggregation, and
    contract-compliant report generation.

    Attributes:
        output_dir: Root ``Path`` for all outputs.
        paths: :class:`OutputPaths` with canonical subdirectories.
        visualizer: :class:`EnhancedVisualizer` used for figure generation.
        results: Accumulated per-species result dictionaries.
        include_field_images: When ``True``, activator/inhibitor field PNGs
            are saved alongside the standard analysis figure.
    """

    def __init__(
        self,
        output_dir: str | Path = "results",
        *,
        include_field_images: bool = False,
    ) -> None:
        """Initialize the simulation runner and create output directories.

        Args:
            output_dir: Filesystem path (string or ``Path``) for the root
                output directory.  Created if it does not exist.
            include_field_images: If ``True``, per-run activator (u) and
                inhibitor (v) concentration fields are saved as additional
                PNG artifacts under the species pattern directory.
        """
        self.output_dir: Path = Path(output_dir)
        self.paths: OutputPaths = _ensure_dirs(self.output_dir)
        self.visualizer: EnhancedVisualizer = EnhancedVisualizer()
        self.results: list[Dict[str, Any]] = []
        # Optional artifact family: per-species activator/inhibitor field PNGs.
        self.include_field_images: bool = bool(include_field_images)

    def run_single_species(
        self,
        species_key: Optional[str] = None,
        D_v_D_u: Optional[float] = None,
        b: Optional[float] = None,
        rho: Optional[float] = None,
        grid_size: Optional[int] = None,
        T_target: Optional[float] = None,
        n_runs: int = 1,
        verbose: bool = True,
    ) -> Dict[str, Any]:
        """Run simulation for a single species with optional parameter overrides.

        Parameters are first resolved from ``SPECIES_DATABASE`` (when
        *species_key* is given) and then selectively overridden by any
        non-``None`` keyword arguments.  Multiple replicates are executed
        when ``n_runs > 1``; results are aggregated via
        :meth:`_aggregate_runs`.

        Args:
            species_key: Canonical key in ``SPECIES_DATABASE`` (e.g.
                ``"fomes"``).  If ``None``, a default "Custom" parameter
                set is used.
            D_v_D_u: Override for the inhibitor/activator diffusion ratio.
            b: Override for the Gierer-Meinhardt saturation parameter.
            rho: Override for the source-density parameter.
            grid_size: Override for the spatial grid side length (pixels).
            T_target: Override for the target simulation time.
            n_runs: Number of replicate simulations.  Seeds are
                deterministically derived as ``42 + run_idx``.
            verbose: If ``True``, emit progress bars (via *tqdm* if
                available) and summary statistics to stdout.

        Returns:
            A dictionary of aggregated (if ``n_runs > 1``) or single-run
            results including wavelength estimates, topology counts,
            predicted pore density, simulation metadata, and artifact paths.

        Raises:
            ValueError: If *species_key* is not found in
                ``SPECIES_DATABASE``.
        """
        # --- Resolve base parameters from database or defaults ---
        if species_key:
            if species_key not in SPECIES_DATABASE:
                raise ValueError(f"Unknown species: {species_key}")
            params = dict(SPECIES_DATABASE[species_key])
            params["species_key"] = species_key
            species_name = params["full_name"]
        else:
            # Sensible defaults when no species is specified.
            params = {
                "full_name": "Custom",
                "D_v_D_u": 150.0,
                "b": 1.0,
                "rho": 0.2,
                "grid_size": 512,
                "T_target": 5.0,
                "density_observed": (0, 0),
                "expected_wavelength_px": 0.0,
            }
            species_name = "Custom"

        # --- Apply per-parameter overrides ---
        if D_v_D_u is not None:
            params["D_v_D_u"] = float(D_v_D_u)
        if b is not None:
            params["b"] = float(b)
        if rho is not None:
            params["rho"] = float(rho)
        if grid_size is not None:
            params["grid_size"] = int(grid_size)
        if T_target is not None:
            params["T_target"] = float(T_target)

        if verbose:
            print(f"\n{'='*70}")
            print(f"SPECIES: {species_name}")
            print(f"{'='*70}")
            print("Parameters:")
            print(f"  D_v/D_u: {params['D_v_D_u']}")
            print(f"  b: {params['b']}, ρ: {params['rho']}")
            print(f"  Grid: {params['grid_size']}×{params['grid_size']}")
            print(f"  T_target: {params.get('T_target', 5.0)} (steps computed from dt_final)")
            print(f"  Runs: {n_runs}")

        # --- Execute replicates ---
        if n_runs > 1:
            run_results: list[Dict[str, Any]] = []

            iterable = range(n_runs)
            _tqdm = None
            if verbose:
                try:
                    from tqdm import tqdm as _tqdm  # type: ignore
                except Exception:
                    _tqdm = None

            # Wrap the iterable with tqdm if available for progress tracking.
            if _tqdm is not None:
                iterable = _tqdm(
                    iterable,
                    total=n_runs,
                    desc=f"Runs ({species_name})",
                    unit="run",
                    dynamic_ncols=True,
                    bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]",
                    leave=True,
                )

            for run_idx in iterable:
                # Deterministic seed per replicate for reproducibility.
                seed = 42 + int(run_idx)
                r = self._execute_simulation(
                    params,
                    verbose=verbose,
                    random_seed=seed,
                    run_idx=int(run_idx),
                    species_key=species_key or species_name,
                )
                r["species_name"] = species_name
                r["run_idx"] = int(run_idx)
                r["seed"] = seed
                run_results.append(r)

                # Update tqdm postfix with latest run's key metrics (best-effort).
                if _tqdm is not None and hasattr(iterable, "set_postfix"):
                    try:
                        iterable.set_postfix(
                            wl_px=f"{float(r.get('wavelength_px', float('nan'))):.1f}",
                            dens=f"{float(r.get('density_predicted', float('nan'))):.2f}",
                        )
                    except Exception:
                        pass

            if _tqdm is not None:
                try:
                    iterable.close()
                except Exception:
                    pass

            # Aggregate replicate statistics.
            aggregated = self._aggregate_runs(run_results)
            aggregated["species_key"] = species_key
            aggregated["species_name"] = species_name
            aggregated["n_runs"] = n_runs
            aggregated["params"] = params
            self.results.append(aggregated)
            return aggregated

        # --- Single run (no aggregation) ---

        result = self._execute_simulation(
            params,
            verbose=verbose,
            random_seed=42,
            run_idx=0,
            species_key=species_key or species_name,
        )
        result["species_key"] = species_key
        result["species_name"] = species_name
        result["params"] = params
        result["n_runs"] = n_runs
        self.results.append(result)
        return result

    def _execute_simulation(
        self,
        params: Dict[str, Any],
        verbose: bool = True,
        random_seed: int = 42,
        run_idx: int = 0,
        species_key: str | None = None,
    ) -> Dict[str, Any]:
        """Execute a single Turing-pattern simulation and analyse the result.

        This is the core inner loop: it sets up the kinetics model,
        constructs a :class:`TuringSimulator`, runs the time integration
        with CFL-adaptive time stepping, invokes topology analysis, saves
        pattern visualizations, and assembles a flat result dictionary.

        Args:
            params: Merged parameter dictionary (species defaults +
                user overrides).  Must contain at least ``D_v_D_u``,
                ``b``, ``rho``, ``grid_size``, and ``T_target``.
            verbose: Emit progress to stdout.
            random_seed: Seed for the simulator's PRNG.
            run_idx: Zero-based replicate index (used in file naming).
            species_key: Canonical species identifier for directory naming.

        Returns:
            Flat dictionary of simulation outputs including pattern arrays,
            wavelength estimates, topology counts, time-step metadata,
            and artifact file paths.
        """
        # --- Set up physics and kinetics ---
        D_u: float = 1.0
        D_v: float = D_u * float(params["D_v_D_u"])

        dt_requested: float = float(params.get("dt_initial", 0.0005))

        kinetics = create_kinetics(
            "gierer_meinhardt",
            rho=float(params["rho"]),
            a=0.1,
            b=float(params["b"]),
        )

        sim = TuringSimulator(
            kinetics_model=kinetics,
            D_u=D_u,
            D_v=D_v,
            grid_size=int(params["grid_size"]),
            dx=1.0,
            dt=dt_requested,
            random_seed=int(random_seed),
        )

        # The simulator may adjust dt for CFL stability.
        dt_final: float = float(sim.dt)
        T_target: float = float(params.get("T_target", 5.0))
        # Compute number of steps to reach T_target (ceiling ensures >= T_target).
        steps_computed: int = int(math.ceil(T_target / dt_final))
        T_actual: float = float(steps_computed * dt_final)

        if verbose and dt_final != dt_requested:
            print(f"    [dt adjusted: {dt_requested:.6f} → {dt_final:.6f} (CFL)]")
            print(f"    [T_target={T_target:.2f}, steps={steps_computed}, T_actual={T_actual:.4f}]")

        # Initialize activator/inhibitor fields with small random perturbation.
        sim.initialize(perturbation_amplitude=0.1)

        # --- Progress reporting ---
        use_tqdm: bool = bool(verbose)
        pbar = None
        _tqdm = None
        if use_tqdm:
            try:
                from tqdm import tqdm as _tqdm  # type: ignore
            except Exception:
                _tqdm = None

        if use_tqdm and _tqdm is not None:
            pbar = _tqdm(
                total=steps_computed,
                desc="Simulating",
                unit="step",
                dynamic_ncols=True,
                bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]",
                leave=False,
            )
        elif verbose:
            print("    Simulating...", end="", flush=True)

        # --- Blocked stepping ---
        # Updates progress at a coarse cadence to minimize per-step overhead.
        block_size: int = min(2000, max(1, steps_computed // 50))
        total_blocks: int = max(1, steps_computed // block_size)
        remaining: int = steps_computed

        for block in range(total_blocks):
            n: int = min(block_size, remaining)
            sim.run(num_steps=n, check_convergence=False)
            remaining -= n

            if pbar is not None:
                pbar.update(n)
            elif verbose:
                progress = int((block + 1) / total_blocks * 100)
                print(f"\r    Simulating... {progress}%", end="", flush=True)

        # Flush any residual steps that didn't fill a full block.
        if remaining > 0:
            sim.run(num_steps=remaining, check_convergence=False)
            if pbar is not None:
                pbar.update(remaining)

        if pbar is not None:
            try:
                pbar.close()
            except Exception:
                pass
        elif verbose:
            print(" ✓")

        # --- Topology analysis ---
        if verbose:
            print("    Analyzing topology...", end="", flush=True)

        analyzer = TopologyAnalyzer(sim.u, dx=1.0)
        metrics: Dict[str, Any] = analyzer.compute_all_metrics()

        if verbose:
            print(" ✓")

        # --- Save pattern visualization ---
        sp_key: str = (species_key or params.get("full_name", "Custom")).replace(" ", "_")
        pattern_dir: Path = self.paths.patterns / sp_key
        pattern_dir.mkdir(parents=True, exist_ok=True)

        timestamp: str = datetime.now().strftime("%Y%m%d_%H%M%S")
        pattern_name: str = f"run_{run_idx+1:02d}_seed{random_seed}_{timestamp}.png"
        pattern_path: Path = pattern_dir / pattern_name

        # Provenance annotation for PNGs (species + parameter set).
        params_txt: str = self.visualizer._format_params(
            {
                "D_v_D_u": params.get("D_v_D_u"),
                "b": params.get("b"),
                "rho": params.get("rho"),
                "grid_size": params.get("grid_size"),
                "T_target": params.get("T_target", T_target),
            }
        )
        species_label: str = str(params.get("full_name", species_key or "Custom"))

        # Main per-run multi-panel analysis figure.
        analyzer.visualize_analysis(
            save_path=str(pattern_path),
            show=False,
            suptitle=f"{species_label}\n{params_txt}" if params_txt else species_label,
        )

        # --- Optional raw field images ---
        activator_path: str | None = None
        inhibitor_path: str | None = None
        if self.include_field_images:
            try:
                from matplotlib import pyplot as _plt

                def _save(field: np.ndarray, suffix: str, label: str) -> str:
                    """Save a single concentration-field heatmap to disk."""
                    fig, ax = _plt.subplots(figsize=(6, 6))
                    im = ax.imshow(field, cmap="hot", interpolation="bilinear")
                    ax.set_title(f"{species_label}\n{label}\n{params_txt}", fontsize=10, fontweight="bold")
                    ax.axis("off")
                    cbar = _plt.colorbar(im, ax=ax, fraction=0.046)
                    cbar.set_label(label, fontsize=9)
                    out_path = pattern_dir / pattern_name.replace(".png", f"_{suffix}.png")
                    _plt.savefig(out_path, dpi=300, bbox_inches="tight")
                    _plt.close(fig)
                    return str(out_path)

                activator_path = _save(sim.u, "activator_u", "Activator (u)")
                inhibitor_path = _save(sim.v, "inhibitor_v", "Inhibitor (v)")
            except Exception as e:
                logging.debug(f"Failed to write field images for {sp_key}: {e}")

        # --- Select best wavelength via QC cascade ---
        wl_best, wl_method, wl_qc_pass, wl_qc_reason = select_wavelength_px(metrics)

        # --- Assemble flat result dictionary ---
        result: Dict[str, Any] = {
            # Keep arrays for in-memory ops only; downstream writers must exclude.
            "pattern": sim.u.copy(),
            "inhibitor_field": sim.v.copy(),

            # Wavelength estimates
            "wavelength_px": float(wl_best),  # backward-compatible primary key
            "wavelength_best_px": float(wl_best),
            "wavelength_method": wl_method,
            "wavelength_qc_pass": bool(wl_qc_pass),
            "wavelength_qc_reason": wl_qc_reason,
            "wavelength_fft_px": float(metrics.get("wavelength_fft", np.nan)) if np.isfinite(metrics.get("wavelength_fft", np.nan)) else np.nan,
            "wavelength_fft_qc_pass": bool(metrics.get("wavelength_fft_qc_pass", False)),
            "wavelength_fft_peak_ratio": float(metrics.get("wavelength_fft_peak_ratio", np.nan)) if np.isfinite(metrics.get("wavelength_fft_peak_ratio", np.nan)) else np.nan,
            "wavelength_autocorr_px": float(metrics.get("wavelength_autocorr", np.nan)) if np.isfinite(metrics.get("wavelength_autocorr", np.nan)) else np.nan,
            "wavelength_autocorr_qc_pass": bool(metrics.get("wavelength_autocorr_qc_pass", False)),

            # Topology / morphometrics
            "spots": int(metrics.get("n_components", 0)),
            "euler_chi": int(metrics.get("euler_characteristic", 0)),
            "holes": int(metrics.get("n_holes", 0)),
            "mean_spacing_px": float(metrics.get("mean_spacing", np.nan)) if np.isfinite(metrics.get("mean_spacing", np.nan)) else np.nan,
            "cv_spacing": float(metrics.get("cv_spacing", np.nan)) if np.isfinite(metrics.get("cv_spacing", np.nan)) else np.nan,

            # Energetics
            "final_energy": float(sim.compute_pattern_energy()),

            # Simulation parameters (echoed for traceability)
            "D_v_D_u": float(params["D_v_D_u"]),
            "b": float(params["b"]),
            "rho": float(params["rho"]),
            "grid_size": int(params["grid_size"]),
            "density_observed": tuple(params.get("density_observed", (0, 0))),

            # Time-step / stability metadata
            "dt_requested": dt_requested,
            "dt_final": dt_final,
            "dt_adjusted": bool(sim.stability_info.get("dt_adjusted", False)),
            "dt_max_diffusion": float(sim.stability_info.get("dt_max_diffusion", np.nan)),
            "safety_factor": float(sim.stability_info.get("safety_factor", np.nan)),
            "T_target": T_target,
            "T_actual": T_actual,
            "steps_computed": steps_computed,
            "had_warning": bool(sim.stability_info.get("had_stability_warning", False)),

            # Artifact paths
            "pattern_path": str(pattern_path),
            "activator_field_path": activator_path or "",
            "inhibitor_field_path": inhibitor_path or "",
        }

        # --- Predicted pore density from wavelength + global calibration ---
        result["density_predicted"] = 0.0
        if result["wavelength_px"] and result["wavelength_px"] > 0:
            # Convert wavelength (px) to physical spacing (um) via calibration.
            spacing_um: float = float(result["wavelength_px"]) * float(CALIBRATION_UM_PER_PX)
            if spacing_um > 0:
                # Density = 1 / spacing, converted from um to mm (factor 1000).
                result["density_predicted"] = float(1000.0 / spacing_um)

        if verbose:
            print("\n  Results:")
            print(f"    Wavelength(best): {result['wavelength_px']:.2f} px ({result['wavelength_method']}, QC={result['wavelength_qc_pass']})")
            print(f"    Spots: {result['spots']}")
            print(f"    Density predicted: {result['density_predicted']:.2f} pores/mm")
            print(f"    Pattern saved: {pattern_path}")

        return result

    def _aggregate_runs(
        self, run_results: list[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """Compute summary statistics across replicate simulation runs.

        Aggregates wavelength, spot count, Euler characteristic, and hole
        count using ``nanmean`` / ``nanstd``.  Representative values
        (first run) are kept for diagnostics that are not meaningfully
        averaged (e.g., time-step metadata and artifact paths).

        Args:
            run_results: List of per-run result dictionaries produced by
                :meth:`_execute_simulation`.

        Returns:
            A merged dictionary containing mean/std statistics, the
            first-run pattern array (for visualization), and representative
            diagnostic fields.
        """
        wavelengths: list[float] = [float(r.get("wavelength_px", np.nan)) for r in run_results]
        spots: list[int] = [int(r.get("spots", 0)) for r in run_results]

        aggregated: Dict[str, Any] = {
            # Wavelength statistics across replicates
            "wavelength_mean": float(np.nanmean(wavelengths)),
            "wavelength_std": float(np.nanstd(wavelengths)),
            "wavelength_px": float(np.nanmean(wavelengths)),
            # Spot-count statistics
            "spots_mean": float(np.nanmean(spots)),
            "spots_std": float(np.nanstd(spots)),
            "spots": int(np.nanmean(spots)) if np.isfinite(np.nanmean(spots)) else 0,
            "euler_chi": int(np.nanmean([r.get("euler_chi", 0) for r in run_results])),
            "holes": int(np.nanmean([r.get("holes", 0) for r in run_results])),
            # Representative pattern (first replicate) for visualization.
            "pattern": run_results[0]["pattern"],
            "inhibitor_field": run_results[0].get("inhibitor_field"),
            # Parameters (same across replicates)
            "D_v_D_u": float(run_results[0]["D_v_D_u"]),
            "b": float(run_results[0]["b"]),
            "rho": float(run_results[0]["rho"]),
            "grid_size": int(run_results[0]["grid_size"]),
            "density_observed": run_results[0].get("density_observed", (0, 0)),
            "density_obs_mean": float(np.mean(run_results[0].get("density_observed", (np.nan, np.nan))))
            if run_results[0].get("density_observed") is not None
            else np.nan,
            "density_predicted": float(np.nanmean([r.get("density_predicted", 0.0) for r in run_results])),

            # Wavelength diagnostics (representative from first run)
            "wavelength_method": run_results[0].get("wavelength_method", ""),
            "wavelength_qc_pass": bool(run_results[0].get("wavelength_qc_pass", False)),
            "wavelength_fft_px": float(run_results[0].get("wavelength_fft_px", np.nan)),
            "wavelength_fft_qc_pass": bool(run_results[0].get("wavelength_fft_qc_pass", False)),
            "wavelength_fft_peak_ratio": float(run_results[0].get("wavelength_fft_peak_ratio", np.nan)),
            "wavelength_autocorr_px": float(run_results[0].get("wavelength_autocorr_px", np.nan)),
            "wavelength_autocorr_qc_pass": bool(run_results[0].get("wavelength_autocorr_qc_pass", False)),
            # Time-step metadata (same across replicates)
            "dt_requested": float(run_results[0].get("dt_requested", np.nan)),
            "dt_final": float(run_results[0].get("dt_final", np.nan)),
            "dt_max_diffusion": float(run_results[0].get("dt_max_diffusion", np.nan)),
            "safety_factor": float(run_results[0].get("safety_factor", np.nan)),
            "dt_adjusted": bool(run_results[0].get("dt_adjusted", False)),
            "T_target": float(run_results[0].get("T_target", np.nan)),
            "T_actual": float(run_results[0].get("T_actual", np.nan)),
            "steps_computed": int(run_results[0].get("steps_computed", 0)),
            # True if *any* replicate experienced a stability warning.
            "had_warning": bool(any(r.get("had_warning", False) for r in run_results)),
            # Representative artifact paths (first run)
            "pattern_path": run_results[0].get("pattern_path", ""),
            "activator_field_path": run_results[0].get("activator_field_path", ""),
            "inhibitor_field_path": run_results[0].get("inhibitor_field_path", ""),
        }
        return aggregated

    def run_all_species(
        self,
        n_runs: int = 1,
        grid_size: Optional[int] = None,
        T_target: Optional[float] = None,
        verbose: bool = True,
    ) -> None:
        """Run simulations for every species in ``SPECIES_DATABASE``.

        Iterates over all canonical species keys, delegates to
        :meth:`run_single_species` for each, and finishes by calling
        :meth:`create_comparison_report`.

        Args:
            n_runs: Number of replicate runs per species.
            grid_size: Optional global grid-size override applied to all
                species.
            T_target: Optional global simulation-time override.
            verbose: Emit progress information to stdout.
        """
        if verbose:
            print(f"\n{'='*70}")
            print(f"RUNNING ALL SPECIES (n={n_runs} runs each)")
            if grid_size:
                print(f"Grid size override: {grid_size}×{grid_size}")
            if T_target is not None:
                print(f"T_target override: {T_target}")
            print(f"{'='*70}")

        for species_key in SPECIES_DATABASE.keys():
            self.run_single_species(species_key=species_key, n_runs=n_runs, grid_size=grid_size, T_target=T_target, verbose=verbose)

        self.create_comparison_report()

    def create_comparison_report(self) -> None:
        """Create and save the final (canonical) comparison report.

        Produces four canonical artifacts:

        1. **Human-readable CSV** -- ``tables/validation_summary.csv``
        2. **Machine-readable CSV** -- ``tables/validation_summary_machine.csv``
           (schema enforced by :func:`enforce_machine_schema`)
        3. **JSON summary** -- ``tables/validation_summary.json``
        4. **Comprehensive figure** -- ``figures/COMPREHENSIVE_COMPARISON.png``

        Contract hardening (fail-open design):
            - Canonical CSVs must *always* exist on disk after this call.
            - Canonical JSON must be JSON-serializable (no raw NumPy arrays).
            - Canonical figure is best-effort; a placeholder is written on
              failure.
            - If JSON or figure generation fails, previously written valid
              CSVs are **not** overwritten.
            - All errors are logged to ``logs/report_errors.log``.
        """

        import csv
        import traceback
        import datetime
        import json

        from matplotlib import pyplot as plt
        from fungalmorphospace.contracts.output_contract import CONTRACT_VERSION

        log_file: Path = self.paths.logs / "report_errors.log"

        def _log(stage: str, exc: BaseException) -> None:
            """Append a timestamped error entry to the report error log."""
            try:
                log_file.parent.mkdir(parents=True, exist_ok=True)
                with open(log_file, "a", encoding="utf-8") as f:
                    ts = datetime.datetime.now().isoformat(timespec="seconds")
                    f.write(f"[{ts}] stage={stage}\n")
                    f.write(traceback.format_exc())
                    f.write("\n\n")
            except Exception:
                # Absolute last resort: never crash because logging failed.
                pass

        def _write_placeholder_csv(path: Path, headers: list[str], message: str) -> None:
            """Write a minimal CSV with a single error-message row."""
            path.parent.mkdir(parents=True, exist_ok=True)
            with open(path, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f)
                w.writerow(headers)
                row = [message] + ["" for _ in range(max(0, len(headers) - 1))]
                w.writerow(row)

        def _write_json(path: Path, payload: dict) -> None:
            """Atomically write a JSON file."""
            path.parent.mkdir(parents=True, exist_ok=True)
            with open(path, "w", encoding="utf-8") as f:
                json.dump(payload, f, indent=2, ensure_ascii=False)

        def _write_placeholder_figure(path: Path, message: str) -> None:
            """Write a minimal placeholder PNG with an error message."""
            path.parent.mkdir(parents=True, exist_ok=True)
            try:
                fig = plt.figure(figsize=(10, 6))
                fig.text(0.5, 0.5, message, ha="center", va="center", wrap=True)
                plt.axis("off")
                fig.savefig(path, dpi=200, bbox_inches="tight")
                plt.close(fig)
            except Exception as e:
                _log("placeholder_figure", e)

        # --- Canonical output paths ---
        table_human: Path = self.paths.tables / "validation_summary.csv"
        table_machine: Path = self.paths.tables / "validation_summary_machine.csv"
        json_path: Path = self.paths.tables / "validation_summary.json"
        fig_path: Path = self.paths.figures / "COMPREHENSIVE_COMPARISON.png"

        # --- 0) Order species results by SPECIES_DATABASE key order ---
        # Runner stores results with keys like `species_name`, `D_v_D_u`, etc.
        results_by_key: dict[str, dict] = {}
        extras: list[dict] = []
        for r in list(self.results or []):
            sk = r.get("species_key")
            if sk is not None:
                results_by_key[str(sk)] = r
            else:
                extras.append(r)

        # Maintain canonical ordering from SPECIES_DATABASE, then extras.
        species_results: list[dict] = []
        for k in SPECIES_DATABASE.keys():
            if k in results_by_key:
                species_results.append(results_by_key[k])
        for r in extras:
            if r not in species_results:
                species_results.append(r)

        # --- 1) Backup previous canonical artifacts (best-effort) ---
        try:
            if fig_path.exists():
                self.paths.backups.mkdir(parents=True, exist_ok=True)
                shutil.copy2(fig_path, self.paths.backups / fig_path.name)
        except Exception as e:
            _log("backup_figure", e)

        try:
            if table_human.exists():
                self.paths.backups.mkdir(parents=True, exist_ok=True)
                shutil.copy2(table_human, self.paths.backups / table_human.name)
        except Exception as e:
            _log("backup_table_human", e)

        try:
            if table_machine.exists():
                self.paths.backups.mkdir(parents=True, exist_ok=True)
                shutil.copy2(table_machine, self.paths.backups / table_machine.name)
        except Exception as e:
            _log("backup_table_machine", e)

        try:
            if json_path.exists():
                self.paths.backups.mkdir(parents=True, exist_ok=True)
                shutil.copy2(json_path, self.paths.backups / json_path.name)
        except Exception as e:
            _log("backup_json", e)

        # --- 2) Canonical tables (contract-critical) ---
        machine_df: pd.DataFrame | None = None
        try:
            human_rows: list[dict] = []
            machine_rows: list[dict] = []

            for r in species_results:
                species_name = r.get("species_name") or r.get("full_name") or r.get("species") or ""
                species_key = r.get("species_key") or (r.get("params", {}) or {}).get("species_key") or ""

                # Extract observed density range, defaulting to NaN.
                d_obs = r.get("density_observed", (float("nan"), float("nan")))
                try:
                    dmin = float(d_obs[0])
                except Exception:
                    dmin = float("nan")
                try:
                    dmax = float(d_obs[1])
                except Exception:
                    dmax = float("nan")

                human_rows.append(
                    {
                        "Species": species_name,
                        "D_v/D_u": float(r.get("D_v_D_u", float("nan"))),
                        "b": float(r.get("b", float("nan"))),
                        "rho": float(r.get("rho", float("nan"))),
                        "Grid": int(r.get("grid_size", 0) or 0),
                        "Wavelength_best_px": round(float(r.get("wavelength_best_px", float("nan"))), 2)
                        if r.get("wavelength_best_px") is not None
                        else float("nan"),
                        "Wavelength_method": r.get("wavelength_method", ""),
                        "Wavelength_QC": bool(r.get("wavelength_qc_pass", False)),
                        "Spots": int(r.get("spots", 0) or 0),
                        "Density_predicted_pores_per_mm": round(float(r.get("density_predicted", float("nan"))), 3)
                        if r.get("density_predicted") is not None
                        else float("nan"),
                        "Biological_range_pores_per_mm": f"{dmin:g}-{dmax:g}" if (dmin == dmin and dmax == dmax) else "",
                    }
                )

                machine_rows.append(
                    {
                        "species": species_name,
                        "species_key": species_key,
                        "D_v_D_u": float(r.get("D_v_D_u", float("nan"))),
                        "b": float(r.get("b", float("nan"))),
                        "rho": float(r.get("rho", float("nan"))),
                        "grid_size": int(r.get("grid_size", 0) or 0),
                        "wavelength_best_px": float(r.get("wavelength_best_px", float("nan"))),
                        "wavelength_method": r.get("wavelength_method", ""),
                        "wavelength_qc_pass": bool(r.get("wavelength_qc_pass", False)),
                        "wavelength_fft_px": float(r.get("wavelength_fft_px", float("nan"))),
                        "wavelength_fft_qc_pass": bool(r.get("wavelength_fft_qc_pass", False)),
                        "wavelength_fft_peak_ratio": float(r.get("wavelength_fft_peak_ratio", float("nan"))),
                        "wavelength_autocorr_px": float(r.get("wavelength_autocorr_px", float("nan"))),
                        "wavelength_autocorr_qc_pass": bool(r.get("wavelength_autocorr_qc_pass", False)),
                        "spots": int(r.get("spots", 0) or 0),
                        "holes": int(r.get("holes", 0) or 0),
                        "euler_chi": int(r.get("euler_chi", 0) or 0),
                        "density_predicted_pores_per_mm": float(r.get("density_predicted", float("nan"))),
                        "density_observed_min_pores_per_mm": dmin,
                        "density_observed_max_pores_per_mm": dmax,
                        "dt_requested": float(r.get("dt_requested", float("nan"))),
                        "dt_final": float(r.get("dt_final", float("nan"))),
                        "dt_adjusted": bool(r.get("dt_adjusted", False)),
                        "dt_max_diffusion": float(r.get("dt_max_diffusion", float("nan"))),
                        "safety_factor": float(r.get("safety_factor", float("nan"))),
                        "T_target": float(r.get("T_target", float("nan"))),
                        "T_actual": float(r.get("T_actual", float("nan"))),
                        "steps_computed": int(r.get("steps_computed", 0) or 0),
                        "had_warning": bool(r.get("had_warning", False)),
                        "pattern_path": r.get("pattern_path", ""),
                    }
                )

            import pandas as pd

            human_df = pd.DataFrame(human_rows)
            # Enforce canonical column order and completeness.
            machine_df = enforce_machine_schema(pd.DataFrame(machine_rows), snapshot=False)

            table_human.parent.mkdir(parents=True, exist_ok=True)
            table_machine.parent.mkdir(parents=True, exist_ok=True)

            human_df.to_csv(table_human, index=False)
            machine_df.to_csv(table_machine, index=False)

        except Exception as e:
            _log("canonical_tables", e)
            _write_placeholder_csv(
                table_human,
                ["Species", "D_v/D_u", "b", "rho"],
                "ERROR: failed to write canonical table (see report_errors.log)",
            )
            try:
                import pandas as pd

                empty_machine = enforce_machine_schema(pd.DataFrame([]), snapshot=False)
                empty_machine.to_csv(table_machine, index=False)
            except Exception as e2:
                _log("canonical_tables_machine_fallback", e2)
                _write_placeholder_csv(table_machine, ["error"], "ERROR: failed to write machine CSV (see report_errors.log)")

        # --- 3) Canonical JSON (must be serializable; do not overwrite CSVs on failure) ---
        try:
            rows_for_json: list[dict] = []
            if machine_df is not None and not machine_df.empty:
                rows_for_json = machine_df.to_dict(orient="records")
            else:
                # Fallback: build minimal JSON from raw results.
                for r in species_results:
                    species_name = r.get("species_name") or r.get("full_name") or ""
                    rows_for_json.append(
                        {
                            "species": species_name,
                            "species_key": r.get("species_key") or "",
                            "D_v_D_u": float(r.get("D_v_D_u", float("nan"))),
                            "b": float(r.get("b", float("nan"))),
                            "rho": float(r.get("rho", float("nan"))),
                            "wavelength_best_px": float(r.get("wavelength_best_px", float("nan"))),
                            "density_predicted_pores_per_mm": float(r.get("density_predicted", float("nan"))),
                        }
                    )

            payload: dict = {
                "version": __version__,
                "generated_at": datetime.datetime.now().isoformat(timespec="seconds"),
                "calibration_um_per_px": float(CALIBRATION_UM_PER_PX),
                "species_results": rows_for_json,
                "contract": {
                    "contract_version": CONTRACT_VERSION,
                    "canonical_outputs": {
                        "figure": str(fig_path.as_posix()),
                        "validation_summary_human": str(table_human.as_posix()),
                        "validation_summary_machine": str(table_machine.as_posix()),
                        "validation_summary_json": str(json_path.as_posix()),
                    },
                },
            }
            _write_json(json_path, payload)

        except Exception as e:
            _log("canonical_json", e)
            _write_json(
                json_path,
                {
                    "version": __version__,
                    "generated_at": datetime.datetime.now().isoformat(timespec="seconds"),
                    "error": "failed to write canonical JSON",
                    "see": str(log_file.as_posix()),
                },
            )

        # --- 4) Main comparison figure (best-effort; placeholder on failure) ---
        try:
            fig = self.visualizer.create_comprehensive_figure(species_results)
            fig_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(fig_path, dpi=300, bbox_inches="tight")
            plt.close(fig)
        except Exception as e:
            _log("canonical_figure", e)
            _write_placeholder_figure(
                fig_path,
                "ERROR: failed to generate COMPREHENSIVE_COMPARISON.png\nSee results/logs/report_errors.log",
            )

        # --- 5) Optional per-species field visuals for the report ---
        if self.include_field_images:
            try:
                fields_dir: Path = self.paths.figures / "fields"
                fields_dir.mkdir(parents=True, exist_ok=True)
                self.visualizer.save_species_field_images(species_results, out_dir=fields_dir)
            except Exception as e:
                _log("field_images_report", e)

        # --- 6) Final contract sanity (ensure canonical files exist) ---
        try:
            if not table_human.exists():
                _write_placeholder_csv(table_human, ["Species"], "ERROR: canonical human table missing")
            if not table_machine.exists():
                _write_placeholder_csv(table_machine, ["error"], "ERROR: canonical machine table missing")
            if not json_path.exists():
                _write_json(json_path, {"error": "canonical JSON missing"})
            if not fig_path.exists():
                _write_placeholder_figure(fig_path, "ERROR: canonical figure missing")
        except Exception as e:
            _log("final_contract_sanity", e)



    def export_experiment_bundle(
        self,
        exp_id: str,
        mode: str,
        *,
        include_figure: bool = False,
        include_field_images: bool = False,
        include_json: bool = True,
    ) -> dict[str, str]:
        """Export per-experiment snapshot artifacts and append a ledger row.

        Writes human-readable and machine-readable CSVs, an optional JSON
        summary, and an optional copy of the comprehensive comparison figure
        into a mode-specific subdirectory (``tables/singles/`` or
        ``tables/all/``).  A row is appended to the corresponding
        ``_index.csv`` ledger file.

        This method follows a **fail-open** design for non-essential
        artifacts: if any stage fails, placeholder artifacts are written
        and the full traceback is appended to
        ``results/logs/report_errors.log``.

        Contract notes:
            - Canonical outputs (latest ALL run) remain unchanged by
              single-species runs.
            - Snapshots are written under ``tables/singles`` or
              ``tables/all`` with the *exp_id* prefix.
            - Machine CSV schema is enforced via
              :func:`enforce_machine_schema`.

        Args:
            exp_id: Unique experiment identifier used in file names and
                the ledger (e.g. ``"20240101_120000_fomes"``).
            mode: Either ``"single"`` or ``"all"``, determining the
                output subdirectory.
            include_figure: If ``True``, copy the canonical comparison
                figure into the snapshot directory.
            include_field_images: If ``True``, save per-species
                activator/inhibitor field images.
            include_json: If ``True`` (default), write a JSON summary
                alongside the CSVs.

        Returns:
            Dictionary mapping artifact names to their filesystem paths,
            e.g. ``{"human_csv": "...", "machine_csv": "...", ...}``.

        Raises:
            ValueError: If *mode* is not ``"single"`` or ``"all"``.
            RuntimeError: If no simulation results are available.
        """

        if mode not in {"single", "all"}:
            raise ValueError("mode must be 'single' or 'all'")

        if not self.results:
            raise RuntimeError("No results available. Run simulations before exporting.")

        import traceback
        import csv

        log_path: Path = self.paths.logs / "report_errors.log"

        def _log(stage: str, e: Exception) -> None:
            """Append a timestamped export error to the log file."""
            try:
                log_path.parent.mkdir(parents=True, exist_ok=True)
                with open(log_path, "a", encoding="utf-8") as f:
                    f.write(f"[{datetime.now().isoformat()}] export stage={stage} exp_id={exp_id} mode={mode}\n")
                    f.write(traceback.format_exc())
                    f.write("\n\n")
            except Exception:
                pass

        def _placeholder_csv(p: Path, headers: list[str], message: str) -> None:
            """Write a placeholder CSV with a single error-message row."""
            p.parent.mkdir(parents=True, exist_ok=True)
            with open(p, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f)
                w.writerow(headers)
                w.writerow([message] + [""] * (max(0, len(headers) - 1)))

        def _placeholder_json(p: Path, payload: dict) -> None:
            """Write a placeholder JSON file."""
            p.parent.mkdir(parents=True, exist_ok=True)
            with open(p, "w", encoding="utf-8") as f:
                json.dump(payload, f, indent=2, ensure_ascii=False)

        def _placeholder_figure(p: Path, message: str) -> None:
            """Write a placeholder PNG with an error message."""
            try:
                from matplotlib import pyplot as _plt

                p.parent.mkdir(parents=True, exist_ok=True)
                fig = _plt.figure(figsize=(10, 6))
                _plt.axis("off")
                _plt.text(0.5, 0.5, message, ha="center", va="center", wrap=True)
                fig.savefig(p, dpi=200, bbox_inches="tight")
                _plt.close(fig)
            except Exception:
                # As a last resort, leave it missing; caller will validate.
                pass

        # --- Build results index keyed by species name ---
        results_dict: Dict[str, Dict[str, Any]] = {}
        for r in self.results:
            results_dict[r.get("species_name", r.get("full_name", "Unknown"))] = r

        # --- Ensure subdirectories exist ---
        tables_root: Path = self.paths.tables
        singles_dir: Path = tables_root / "singles"
        all_dir: Path = tables_root / "all"
        singles_dir.mkdir(parents=True, exist_ok=True)
        all_dir.mkdir(parents=True, exist_ok=True)

        figures_all_dir: Path = self.paths.figures / "all"
        figures_all_dir.mkdir(parents=True, exist_ok=True)
        figures_singles_dir: Path = self.paths.figures / "singles"
        figures_singles_dir.mkdir(parents=True, exist_ok=True)

        dest_dir: Path = singles_dir if mode == "single" else all_dir

        human_path: Path = dest_dir / f"validation_summary_{exp_id}.csv"
        machine_path: Path = dest_dir / f"validation_summary_machine_{exp_id}.csv"
        json_path: Path = dest_dir / f"validation_summary_{exp_id}.json"

        # --- Build snapshot rows ---
        human_rows: list[dict[str, Any]] = []
        machine_rows: list[dict[str, Any]] = []

        for species_name, r in results_dict.items():
            human_rows.append(
                {
                    "exp_id": exp_id,
                    "Species": species_name,
                    "D_v/D_u": r.get("D_v_D_u"),
                    "b": r.get("b"),
                    "rho": r.get("rho"),
                    "Grid": r.get("grid_size"),
                    "λ(best) (px)": f"{float(r.get('wavelength_px', 0.0)):.2f}" if r.get("wavelength_px") else "",
                    "λ(method)": r.get("wavelength_method", ""),
                    "λ(QC)": bool(r.get("wavelength_qc_pass", False)),
                    "Spots": r.get("spots"),
                    "Euler_chi": r.get("euler_chi"),
                    "Holes": r.get("holes"),
                    "Density Pred. (pores/mm)": f"{float(r.get('density_predicted', 0.0)):.3f}",
                    "Density Obs. (pores/mm)": f"{r.get('density_observed', (0,0))[0]}-{r.get('density_observed', (0,0))[1]}",
                    "dt_final": f"{float(r.get('dt_final', float('nan'))):.6f}" if np.isfinite(r.get("dt_final", np.nan)) else "",
                    "T_target": r.get("T_target", ""),
                    "T_actual": f"{float(r.get('T_actual', float('nan'))):.4f}" if np.isfinite(r.get("T_actual", np.nan)) else "",
                    "dt_adjusted": bool(r.get("dt_adjusted", False)),
                    "had_warning": bool(r.get("had_warning", False)),
                    "pattern_path": r.get("pattern_path", ""),
                }
            )

            machine_rows.append(
                {
                    "exp_id": exp_id,
                    "species": species_name,
                    "species_key": r.get("species_key") or (r.get("params", {}) or {}).get("species_key"),
                    "D_v_D_u": float(r.get("D_v_D_u", np.nan)),
                    "b": float(r.get("b", np.nan)),
                    "rho": float(r.get("rho", np.nan)),
                    "grid_size": int(r.get("grid_size", 0)),
                    "wavelength_best_px": float(r.get("wavelength_px", np.nan)),
                    "wavelength_method": r.get("wavelength_method", ""),
                    "wavelength_qc_pass": bool(r.get("wavelength_qc_pass", False)),
                    "wavelength_fft_px": float(r.get("wavelength_fft_px", np.nan)),
                    "wavelength_fft_qc_pass": bool(r.get("wavelength_fft_qc_pass", False)),
                    "wavelength_fft_peak_ratio": float(r.get("wavelength_fft_peak_ratio", np.nan)),
                    "wavelength_autocorr_px": float(r.get("wavelength_autocorr_px", np.nan)),
                    "wavelength_autocorr_qc_pass": bool(r.get("wavelength_autocorr_qc_pass", False)),
                    "spots": int(r.get("spots", 0)),
                    "holes": int(r.get("holes", 0)),
                    "euler_chi": int(r.get("euler_chi", 0)),
                    "density_predicted_pores_per_mm": float(r.get("density_predicted", np.nan)),
                    "density_observed_min_pores_per_mm": float(r.get("density_observed", (np.nan, np.nan))[0]),
                    "density_observed_max_pores_per_mm": float(r.get("density_observed", (np.nan, np.nan))[1]),
                    "dt_requested": float(r.get("dt_requested", np.nan)),
                    "dt_final": float(r.get("dt_final", np.nan)),
                    "dt_adjusted": bool(r.get("dt_adjusted", False)),
                    "dt_max_diffusion": float(r.get("dt_max_diffusion", np.nan)),
                    "safety_factor": float(r.get("safety_factor", np.nan)),
                    "T_target": float(r.get("T_target", np.nan)) if r.get("T_target") is not None else np.nan,
                    "T_actual": float(r.get("T_actual", np.nan)),
                    "steps_computed": int(r.get("steps_computed", 0)),
                    "had_warning": bool(r.get("had_warning", False)),
                    "pattern_path": r.get("pattern_path", ""),
                }
            )

        # --- Write snapshot CSVs (fail-open) ---
        try:
            pd.DataFrame(human_rows).to_csv(human_path, index=False)
        except Exception as e:
            _log('human_csv', e)
            headers = list(human_rows[0].keys()) if human_rows else ['exp_id', 'Species']
            _placeholder_csv(human_path, headers, f'ERROR writing human CSV; see {log_path.name}')

        try:
            dfm: pd.DataFrame = enforce_machine_schema(pd.DataFrame(machine_rows), snapshot=True)
            dfm.to_csv(machine_path, index=False)
        except Exception as e:
            _log('machine_csv', e)
            dfm = enforce_machine_schema(pd.DataFrame([]), snapshot=True)
            dfm.to_csv(machine_path, index=False)

        # --- Optional JSON snapshot ---
        payload: dict[str, Any] = {}
        if include_json:
            try:
                payload = {
                    'exp_id': exp_id,
                    'mode': mode,
                    'version': __version__,
                    'generated_at': datetime.now().isoformat(),
                    'calibration_um_per_px': float(CALIBRATION_UM_PER_PX),
                    'species': {
                        name: {
                            'species_key': (r.get('species_key') or (r.get('params', {}) or {}).get('species_key')),
                            'D_v_D_u': float(r.get('D_v_D_u', np.nan)),
                            'b': float(r.get('b', np.nan)),
                            'rho': float(r.get('rho', np.nan)),
                            'grid_size': int(r.get('grid_size', 0)),
                            'wavelength_best_px': float(r.get('wavelength_px', np.nan)),
                            'wavelength_method': r.get('wavelength_method', ''),
                            'wavelength_qc_pass': bool(r.get('wavelength_qc_pass', False)),
                            'spots': int(r.get('spots', 0)),
                            'holes': int(r.get('holes', 0)),
                            'euler_chi': int(r.get('euler_chi', 0)),
                            'density_predicted_pores_per_mm': float(r.get('density_predicted', np.nan)),
                            'dt_final': float(r.get('dt_final', np.nan)),
                            'T_target': float(r.get('T_target', np.nan)) if r.get('T_target') is not None else np.nan,
                            'T_actual': float(r.get('T_actual', np.nan)),
                            'steps_computed': int(r.get('steps_computed', 0)),
                            'pattern_path': r.get('pattern_path', ''),
                        }
                        for name, r in results_dict.items()
                    },
                }
                _placeholder_json(json_path, payload)
            except Exception as e:
                _log('json', e)
                _placeholder_json(json_path, {
                    'exp_id': exp_id,
                    'mode': mode,
                    'version': __version__,
                    'generated_at': datetime.now().isoformat(),
                    'error': 'ERROR writing snapshot JSON; see report_errors.log',
                })

        # --- Optional figure copy ---
        fig_copy_path: Path | None = None
        if include_figure:
            try:
                canonical_fig: Path = self.paths.figures / 'COMPREHENSIVE_COMPARISON.png'
                fig_copy_path = figures_all_dir / f'COMPREHENSIVE_COMPARISON_{exp_id}.png'
                if canonical_fig.exists():
                    shutil.copy(canonical_fig, fig_copy_path)
                else:
                    _placeholder_figure(fig_copy_path, 'ERROR: canonical comparison figure missing')
            except Exception as e:
                _log('figure_copy', e)
                if fig_copy_path is None:
                    fig_copy_path = figures_all_dir / f'COMPREHENSIVE_COMPARISON_{exp_id}.png'
                _placeholder_figure(fig_copy_path, f'ERROR copying figure; see {log_path.name}')

        # --- Optional field images ---
        if include_field_images:
            try:
                fields_dest: Path = self.paths.figures / ('all' if mode == 'all' else 'singles')
                fields_dest.mkdir(parents=True, exist_ok=True)
                self.visualizer.save_species_field_images(results_dict, out_dir=fields_dest, exp_id=exp_id)
            except Exception as e:
                _log('field_images', e)

        # --- Append ledger row ---
        ledger_path: Path = dest_dir / '_index.csv'
        try:
            self._append_ledger_row(ledger_path, exp_id=exp_id, mode=mode, results_dict=results_dict)
        except Exception as e:
            _log('ledger', e)
            # Minimal ledger placeholder.
            if not ledger_path.exists():
                _placeholder_csv(ledger_path, ['exp_id', 'mode', 'generated_at', 'version', 'n_species', 'species_json'], f'ERROR writing ledger; see {log_path.name}')

        # --- Final validation: ensure snapshot artifacts exist ---
        if not human_path.exists():
            _placeholder_csv(human_path, ['exp_id', 'Species'], f'ERROR: missing human CSV; see {log_path.name}')
        if not machine_path.exists():
            dfm = enforce_machine_schema(pd.DataFrame([]), snapshot=True)
            dfm.to_csv(machine_path, index=False)
        if include_json and not json_path.exists():
            _placeholder_json(json_path, {
                'exp_id': exp_id,
                'mode': mode,
                'version': __version__,
                'generated_at': datetime.now().isoformat(),
                'error': 'ERROR: missing snapshot JSON; see report_errors.log',
            })
        if include_figure and fig_copy_path is not None and not fig_copy_path.exists():
            _placeholder_figure(fig_copy_path, f'ERROR: missing snapshot figure; see {log_path.name}')

        # --- Return paths of all written artifacts ---
        out: dict[str, str] = {
            'human_csv': str(human_path),
            'machine_csv': str(machine_path),
            'ledger': str(ledger_path),
        }
        if include_json:
            out['json'] = str(json_path)
        if include_figure and fig_copy_path is not None:
            out['figure'] = str(fig_copy_path)
        return out


    @staticmethod
    def _append_ledger_row(
        ledger_path: Path,
        *,
        exp_id: str,
        mode: str,
        results_dict: Dict[str, Dict[str, Any]],
    ) -> None:
        """Append a single experiment row to the append-only ledger CSV.

        The ledger uses a **union schema expansion** strategy: new columns
        introduced by a row are appended to the existing schema, and
        pre-existing rows receive ``pd.NA`` for those new columns.  Rows
        are de-duplicated by *exp_id* (keeping the latest entry).

        Core columns are always placed first in deterministic order, followed
        by per-species wide columns sorted alphabetically.

        Args:
            ledger_path: Filesystem path to the ``_index.csv`` ledger file.
                Created if it does not exist.
            exp_id: Unique experiment identifier.
            mode: ``"single"`` or ``"all"``.
            results_dict: Species-name-keyed results dictionary.
        """
        # --- Build the new row ---
        row: dict[str, Any] = {
            "exp_id": exp_id,
            "mode": mode,
            "generated_at": datetime.now().isoformat(),
            "version": __version__,
            "n_species": int(len(results_dict)),
        }

        # Compact JSON summary (portable; avoids schema explosion for long-term use).
        compact: dict[str, dict[str, Any]] = {
            name: {
                "species_key": (r.get("species_key") or (r.get("params", {}) or {}).get("species_key")),
                "D_v_D_u": float(r.get("D_v_D_u", np.nan)),
                "b": float(r.get("b", np.nan)),
                "rho": float(r.get("rho", np.nan)),
                "grid_size": int(r.get("grid_size", 0)),
                "wavelength_best_px": float(r.get("wavelength_px", np.nan)),
                "density_predicted_pores_per_mm": float(r.get("density_predicted", np.nan)),
                "dt_final": float(r.get("dt_final", np.nan)),
                "pattern_path": r.get("pattern_path", ""),
            }
            for name, r in results_dict.items()
        }
        row["species_json"] = json.dumps(compact, ensure_ascii=False)

        # Convenience wide columns (one set per species for quick filtering).
        for name, r in results_dict.items():
            sk: str = (r.get("species_key") or (r.get("params", {}) or {}).get("species_key") or name).lower().replace(" ", "_")
            row[f"{sk}_rho"] = float(r.get("rho", np.nan))
            row[f"{sk}_D_v_D_u"] = float(r.get("D_v_D_u", np.nan))
            row[f"{sk}_b"] = float(r.get("b", np.nan))
            row[f"{sk}_wavelength_best_px"] = float(r.get("wavelength_px", np.nan))
            row[f"{sk}_density_predicted_pores_per_mm"] = float(r.get("density_predicted", np.nan))
            row[f"{sk}_grid_size"] = int(r.get("grid_size", 0))
            row[f"{sk}_dt_final"] = float(r.get("dt_final", np.nan))

        # --- Load or create the ledger DataFrame ---
        if ledger_path.exists():
            try:
                df: pd.DataFrame = pd.read_csv(ledger_path)
            except Exception:
                df = pd.DataFrame()
        else:
            df = pd.DataFrame()

        # --- Build union schema (existing columns + any new columns in this row) ---
        existing_cols: list[str] = list(df.columns)
        for c in row.keys():
            if c not in existing_cols:
                existing_cols.append(c)

        # Align row to union schema (missing columns get pd.NA).
        row_aligned: dict[str, Any] = {c: row.get(c, pd.NA) for c in existing_cols}

        if df.shape[0] == 0:
            # Avoid pandas concat deprecation warning for empty frames.
            df = pd.DataFrame([row_aligned], columns=existing_cols)
        else:
            # Ensure the existing DataFrame has all union columns, then append.
            for c in existing_cols:
                if c not in df.columns:
                    df[c] = pd.NA
            df = pd.concat([df, pd.DataFrame([row_aligned], columns=existing_cols)], ignore_index=True)

        # De-duplicate by exp_id, keeping the most recent entry.
        if "exp_id" in df.columns:
            df = df.drop_duplicates(subset=["exp_id"], keep="last")

        # Stable column ordering: core columns first, then sorted extras.
        core: list[str] = ["exp_id", "mode", "generated_at", "version", "n_species", "species_json"]
        others: list[str] = [c for c in df.columns if c not in core]
        df = df[core + sorted(others)]

        df.to_csv(ledger_path, index=False)
