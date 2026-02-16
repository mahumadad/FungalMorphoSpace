# API Reference

FungalMorphoSpace is organized into five modules.

## Core: Simulation Engine

### `fungalmorphospace.core.TuringSimulator`

The main simulation class. Solves the Gierer-Meinhardt reaction-diffusion system on a 2D grid with periodic boundary conditions.

```python
from fungalmorphospace.core import TuringSimulator, GiererMeinhardtKinetics

kinetics = GiererMeinhardtKinetics(rho=0.2, a=0.1, b=1.0)
sim = TuringSimulator(
    kinetics_model=kinetics,
    D_u=0.1, D_v=15.0,
    grid_size=512
)
sim.initialize_fields(noise=0.1)
sim.run(T_target=5.0)
```

### `fungalmorphospace.core.kinetics`

Available kinetics models:
- `GiererMeinhardtKinetics` — Saturating autocatalysis (used for hymenophore patterns)
- `SchnakenbergKinetics` — Cross-catalytic model
- `GrayScottKinetics` — Feed/kill model

Factory function: `create_kinetics(model_name, **params)`

## Analysis: Pattern Quantification

### `fungalmorphospace.analysis.TopologyAnalyzer`

Extracts morphometric and topological metrics from simulated patterns.

```python
from fungalmorphospace.analysis import TopologyAnalyzer

analyzer = TopologyAnalyzer(pattern_2d, dx=1.0)
metrics = analyzer.compute_all_metrics()
# Returns: wavelength_fft, wavelength_autocorr, n_components,
#          euler_characteristic, mean_spacing, cv_spacing, ...
```

### `fungalmorphospace.analysis.EnhancedVisualizer`

Publication-quality figure generation: species comparisons, scaling law plots, sensitivity heatmaps.

## Runners: Pipeline Orchestration

### `fungalmorphospace.runners.IntegratedSimulationRunner`

End-to-end validation pipeline for one or multiple species.

```python
from fungalmorphospace.runners import IntegratedSimulationRunner

runner = IntegratedSimulationRunner(output_dir="results/", n_runs=5, grid_size=512)
runner.run_all_species()
```

### `fungalmorphospace.runners.SPECIES_DATABASE`

Dict of calibrated species parameters loaded from `data/species_data.json`.

### `fungalmorphospace.runners.CALIBRATION_UM_PER_PX`

Spatial calibration constant (8.7 um/px, derived from *Fomes fomentarius* SEM data).

## Contracts: Output Schema

### `fungalmorphospace.contracts.output_contract`

Defines canonical output directory structure, CSV column schemas, and validation functions. See [OUTPUT_CONTRACT.md](OUTPUT_CONTRACT.md) for the full specification.

## Utils: Utilities

### `fungalmorphospace.utils.SensitivityAnalyzer`

2D parameter sweep tool for exploring the morphospace.

```python
from fungalmorphospace.utils.sensitivity_analysis import SensitivityAnalyzer

sa = SensitivityAnalyzer(base_params={...})
results = sa.sweep_2d("D_v_D_u", [50, 100, 200], "rho", [0.1, 0.2, 0.5])
```
