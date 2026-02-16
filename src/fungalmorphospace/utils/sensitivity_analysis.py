#!/usr/bin/env python3
"""Sensitivity analysis module for systematic parameter-space exploration.

Provides tools for performing two-dimensional parameter sweeps over
Turing reaction-diffusion simulations, enabling quantitative assessment
of how pattern morphometrics (wavelength, spot count, Euler characteristic)
respond to variations in model parameters such as the diffusion ratio
(``D_v/D_u``), saturation constant (``b``), and source density (``rho``).

Typical usage::

    analyzer = SensitivityAnalyzer(base_params={
        'D_u': 1.0,
        'grid_size': 256,
        'dx': 1.0,
        'dt': 0.0005,
        'steps': 5000,
        'model': 'gierer_meinhardt',
    })
    df = analyzer.sweep_2d('b', [0.8, 1.0, 1.5], 'rho', [0.2, 0.5, 1.0])
    analyzer.save_results('sensitivity_results.csv')

See Also:
    :class:`fungalmorphospace.core.turing_simulator.TuringSimulator`
    :class:`fungalmorphospace.analysis.topology_analyzer.TopologyAnalyzer`
"""

import numpy as np
import pandas as pd
from pathlib import Path
import json
from datetime import datetime
from typing import Any, Dict, Optional, Sequence

from ..core.turing_simulator import TuringSimulator
from ..core.kinetics import create_kinetics
from ..analysis.topology_analyzer import TopologyAnalyzer



class SensitivityAnalyzer:
    """Two-dimensional parameter sensitivity analysis for Turing patterns.

    Runs a grid of simulations over two user-specified parameters,
    collecting topology metrics for each combination.  Supports
    statistical replication with unique per-run random seeds.

    Attributes:
        base_params: Default simulation configuration used as the
            starting point for each sweep point.
        results: List of per-run result dictionaries accumulated across
            all calls to :meth:`sweep_2d`.
    """

    def __init__(self, base_params: Optional[Dict[str, Any]] = None) -> None:
        """Initialize the analyzer with a base parameter configuration.

        Args:
            base_params: Base configuration dictionary for simulations.
                If ``None``, sensible defaults for the Gierer-Meinhardt
                model are used (``D_u=1.0``, ``grid_size=256``,
                ``dx=1.0``, ``dt=0.0005``, ``steps=5000``).
        """
        self.base_params: Dict[str, Any] = base_params or {
            'D_u': 1.0,
            'grid_size': 256,
            'dx': 1.0,
            'dt': 0.0005,
            'steps': 5000,
            'model': 'gierer_meinhardt'
        }

        self.results: list[Dict[str, Any]] = []

    def sweep_2d(
        self,
        param1_name: str,
        param1_range: Sequence[float],
        param2_name: str,
        param2_range: Sequence[float],
        n_replicates: int = 1,
        verbose: bool = True,
    ) -> pd.DataFrame:
        """Execute a 2D parameter sweep and collect topology metrics.

        Iterates over every combination of *param1_range* x *param2_range*,
        running *n_replicates* independent simulations per combination.
        Each run uses a unique random seed (``42 + run_counter``) for true
        statistical replication.

        Args:
            param1_name: Name of the first parameter to vary.  Supported
                values: ``"b"``, ``"rho"``, ``"D_v_D_u"``.
            param1_range: Sequence of values for the first parameter.
            param2_name: Name of the second parameter to vary (same
                supported values as *param1_name*).
            param2_range: Sequence of values for the second parameter.
            n_replicates: Number of independent replicate runs per
                (param1, param2) combination.
            verbose: If ``True``, print progress to stdout.

        Returns:
            A :class:`pandas.DataFrame` containing all accumulated results
            (including results from previous ``sweep_2d`` calls on the
            same instance).
        """

        if verbose:
            print(f"\n{'='*70}")
            print(f"SENSITIVITY ANALYSIS: {param1_name} × {param2_name}")
            print(f"{'='*70}")
            print(f"Grid: {len(param1_range)} × {len(param2_range)} = {len(param1_range)*len(param2_range)} combinations")
            print(f"Replicates: {n_replicates}")
            print(f"Total runs: {len(param1_range)*len(param2_range)*n_replicates}")

        total_runs: int = len(param1_range) * len(param2_range) * n_replicates
        run_counter: int = 0

        for p1_val in param1_range:
            for p2_val in param2_range:
                for rep in range(n_replicates):
                    run_counter += 1

                    if verbose:
                        print(f"\r  Run {run_counter}/{total_runs}: "
                              f"{param1_name}={p1_val:.3f}, {param2_name}={p2_val:.3f}, "
                              f"rep={rep+1}/{n_replicates}", end='')

                    # Build merged params for this specific (p1, p2) combination.
                    run_params: Dict[str, Any] = self._build_params(param1_name, p1_val,
                                                    param2_name, p2_val)

                    # Unique seed per run for true statistical replication.
                    run_seed: int = 42 + run_counter

                    # Run simulation with unique seed.
                    result: Dict[str, Any] = self._run_single(run_params, random_seed=run_seed)

                    # Tag the result with sweep coordinates and metadata.
                    result[param1_name] = p1_val
                    result[param2_name] = p2_val
                    result['replicate'] = rep
                    result['seed'] = run_seed
                    self.results.append(result)

        if verbose:
            print("\n✓ Sweep complete")

        return pd.DataFrame(self.results)

    def _build_params(
        self,
        param1_name: str,
        param1_val: float,
        param2_name: str,
        param2_val: float,
    ) -> Dict[str, Any]:
        """Build a complete parameter dictionary for a single sweep point.

        Starts from :attr:`base_params` and overlays the two swept
        parameter values.  Kinetics parameters (``b``, ``rho``) are placed
        in a nested ``kinetics_params`` sub-dict, while diffusion-related
        parameters (``D_v_D_u``) update ``D_v`` directly.

        Args:
            param1_name: Name of the first swept parameter.
            param1_val: Value for the first swept parameter at this point.
            param2_name: Name of the second swept parameter.
            param2_val: Value for the second swept parameter at this point.

        Returns:
            Merged parameter dictionary ready for :meth:`_run_single`.
        """

        params: Dict[str, Any] = self.base_params.copy()

        # Map parameter names to kinetics kwargs vs. simulator kwargs.
        kinetics_params: Dict[str, float] = {}

        if param1_name == 'b':
            kinetics_params['b'] = param1_val
        elif param1_name == 'rho':
            kinetics_params['rho'] = param1_val
        elif param1_name == 'D_v_D_u':
            params['D_v'] = params['D_u'] * param1_val

        if param2_name == 'b':
            kinetics_params['b'] = param2_val
        elif param2_name == 'rho':
            kinetics_params['rho'] = param2_val
        elif param2_name == 'D_v_D_u':
            params['D_v'] = params['D_u'] * param2_val

        params['kinetics_params'] = kinetics_params

        return params

    def _run_single(
        self,
        params: Dict[str, Any],
        random_seed: int = 42,
    ) -> Dict[str, Any]:
        """Execute a single simulation and return topology metrics.

        Creates a :class:`TuringSimulator`, runs the time integration,
        and analyses the resulting activator field with
        :class:`TopologyAnalyzer`.

        Args:
            params: Simulation parameters including ``D_u``, ``D_v`` (or
                default ratio 150), ``grid_size``, ``dx``, ``dt``,
                ``steps``, ``model``, and optional ``kinetics_params``.
            random_seed: Seed for the simulator's PRNG.  Each run in a
                sweep receives a unique seed for true replication.

        Returns:
            Dictionary containing:

            - ``wavelength`` (``float``): Autocorrelation-based wavelength.
            - ``spots`` (``int``): Number of connected components.
            - ``euler_chi`` (``int``): Euler characteristic.
            - ``final_energy`` (``float``): Pattern energy at final time.
            - ``D_v_D_u_actual`` (``float``): Realized diffusion ratio.
        """

        # Create kinetics model with swept parameters.
        kinetics_params: Dict[str, float] = params.get('kinetics_params', {})
        kinetics = create_kinetics(params['model'], **kinetics_params)

        # Compute inhibitor diffusion coefficient (default ratio 150).
        D_v: float = params.get('D_v', params['D_u'] * 150)

        sim = TuringSimulator(
            kinetics_model=kinetics,
            D_u=params['D_u'],
            D_v=D_v,
            grid_size=params['grid_size'],
            dx=params['dx'],
            dt=params['dt'],
            random_seed=random_seed  # CRITICAL: unique seed per run
        )

        # Run the full time integration.
        sim.initialize(perturbation_amplitude=0.1)
        sim.run(num_steps=params['steps'], check_convergence=False)

        # Analyze the final activator field.
        analyzer = TopologyAnalyzer(sim.u, dx=params['dx'])
        metrics: Dict[str, Any] = analyzer.compute_all_metrics()

        # Compile topology metrics and energetics.
        result: Dict[str, Any] = {
            'wavelength': metrics.get('wavelength_autocorr', 0),
            'spots': metrics['n_components'],
            'euler_chi': metrics['euler_characteristic'],
            'final_energy': sim.compute_pattern_energy(),
            'D_v_D_u_actual': D_v / params['D_u']
        }

        return result

    def save_results(self, filename: str | Path) -> pd.DataFrame:
        """Persist accumulated sweep results to a CSV file.

        Args:
            filename: Output file path (string or ``Path``).

        Returns:
            The :class:`pandas.DataFrame` that was written.
        """
        df: pd.DataFrame = pd.DataFrame(self.results)
        df.to_csv(filename, index=False)
        print(f"✓ Results saved: {filename}")
        return df

    def get_summary_stats(self) -> Dict[str, float]:
        """Compute summary statistics over all accumulated results.

        Returns:
            Dictionary with keys ``wavelength_mean``, ``wavelength_std``,
            ``spots_mean``, ``spots_std``, and ``n_runs``.
        """
        df: pd.DataFrame = pd.DataFrame(self.results)

        summary: Dict[str, float] = {
            'wavelength_mean': float(df['wavelength'].mean()),
            'wavelength_std': float(df['wavelength'].std()),
            'spots_mean': float(df['spots'].mean()),
            'spots_std': float(df['spots'].std()),
            'n_runs': float(len(df)),
        }

        return summary


def quick_sensitivity_test() -> pd.DataFrame:
    """Run a small demonstration sensitivity sweep for smoke-testing.

    Creates a :class:`SensitivityAnalyzer` with a reduced grid
    (128 x 128, 2000 steps) and sweeps ``b`` x ``rho`` over 5 x 5
    values with 2 replicates each (50 total simulations).

    Returns:
        DataFrame of sweep results.
    """

    analyzer = SensitivityAnalyzer(base_params={
        'D_u': 1.0,
        'grid_size': 128,  # Small grid for speed
        'dx': 1.0,
        'dt': 0.0005,
        'steps': 2000,
        'model': 'gierer_meinhardt'
    })

    # Sweep: b (saturation) vs rho (source density).
    b_range: np.ndarray = np.linspace(0.8, 2.0, 5)
    rho_range: np.ndarray = np.linspace(0.2, 1.0, 5)

    df: pd.DataFrame = analyzer.sweep_2d('b', b_range, 'rho', rho_range, n_replicates=2)

    print("\n" + "="*70)
    print("RESULTS SUMMARY")
    print("="*70)
    print(df.describe())

    return df


if __name__ == "__main__":
    quick_sensitivity_test()
