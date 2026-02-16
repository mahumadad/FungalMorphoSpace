# Examples

## 1. Single Species Simulation

Run a simulation for *Fomes fomentarius* (baseline species, small pores):

```bash
python scripts/run_integrated_validation.py --species fomes --n_runs 3
```

Output is saved to `results/<timestamp>/` with pattern images, wavelength metrics, and CSV summaries.

## 2. Full Multi-Species Validation

Validate all three calibrated species and generate comparison figures:

```bash
python scripts/run_integrated_validation.py --species all --n_runs 5 --grid 1024
```

This produces:
- `figures/COMPREHENSIVE_COMPARISON.png` — Side-by-side species comparison
- `tables/validation_summary_machine.csv` — Machine-readable metrics
- `tables/validation_summary.json` — JSON summary

## 3. Scaling Law Analysis

After running the validation, generate the allometric scaling plot:

```bash
python scripts/plot_scaling_law.py --results-dir results/<timestamp>/
```

## 4. Topological Continuity Test

Demonstrate that pores and gills are topological variants of the same mechanism:

```bash
python scripts/test_laminillas.py
```

This varies D_v/D_u while tracking the Euler characteristic to show smooth transitions.

## 5. Sensitivity Analysis

Explore how wavelength depends on diffusion ratio and reaction intensity:

```python
from fungalmorphospace.utils.sensitivity_analysis import SensitivityAnalyzer

analyzer = SensitivityAnalyzer(base_params={
    "D_u": 0.1, "D_v_D_u": 150, "rho": 0.2,
    "a": 0.1, "b": 1.0, "grid_size": 256
})
results = analyzer.sweep_2d(
    "D_v_D_u", [50, 100, 150, 200, 300],
    "rho", [0.1, 0.2, 0.4, 0.8]
)
```

## 6. Programmatic API Usage

```python
from fungalmorphospace.core import TuringSimulator, GiererMeinhardtKinetics
from fungalmorphospace.analysis import TopologyAnalyzer

# Configure Gierer-Meinhardt kinetics for Fomes fomentarius
kinetics = GiererMeinhardtKinetics(rho=0.2, a=0.1, b=1.0)

# Create simulator with calibrated diffusion ratio
sim = TuringSimulator(
    kinetics_model=kinetics,
    D_u=0.1,
    D_v=15.0,   # D_v/D_u = 150
    grid_size=512
)

# Initialize with small random perturbation around steady state
sim.initialize_fields(noise=0.1)

# Run to target physical time
sim.run(T_target=5.0)

# Analyze resulting pattern
analyzer = TopologyAnalyzer(sim.u, dx=1.0)
metrics = analyzer.compute_all_metrics()

print(f"Wavelength (FFT):    {metrics['wavelength_fft']:.1f} px")
print(f"Wavelength (autocorr): {metrics['wavelength_autocorr']:.1f} px")
print(f"Components:          {metrics['n_components']}")
print(f"Euler characteristic: {metrics['euler_characteristic']}")
```

## Species Keys

| Short key | Aliases | Scientific name |
|-----------|---------|-----------------|
| `fomes` | `fomes_fomentarius` | *Fomes fomentarius* |
| `brumalis` | `polyporus_brumalis`, `lentinus_brumalis` | *Lentinus brumalis* |
| `squamosus` | `polyporus_squamosus` | *Polyporus squamosus* |
