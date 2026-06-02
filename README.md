# FungalMorphoSpace

**Reaction-diffusion morphospace simulator for fungal hymenophores**

[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)
[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)
[![Version](https://img.shields.io/badge/version-0.7.3.post6-green.svg)]()

---

FungalMorphoSpace is a Python package for simulating and quantifying fungal hymenophore patterns using Gierer-Meinhardt reaction-diffusion systems. It explores the hypothesis that diverse hymenophore morphologies (pores, gills, labyrinths) emerge from parametric variation within a single self-organizing mechanism ([Kuhar et al., 2022](https://doi.org/10.1007/s12064-022-00363-z)).

## Key Features

- **Turing pattern simulation** with selectable kinetics (Gierer-Meinhardt, Schnakenberg, Gray-Scott)
- **Topological and morphometric analysis**: wavelength estimation (FFT + autocorrelation), Euler characteristic, connected components, spacing statistics
- **Calibrated species database** for 3 polypore species spanning an 8x scale range
- **Reproducible output contracts** with machine-readable CSV schemas and append-only ledgers
- **Sensitivity analysis** tools for parameter space exploration

## Validated Species

| Species | Key | lambda (px) | lambda (um) | D_v/D_u | rho | b |
|---------|-----|------------|------------|---------|-----|---|
| *Fomes fomentarius* | `fomes` | 46 | 400 | 150 | 0.2 | 1.0 |
| *Lentinus brumalis* | `brumalis` | 33 | 261 | 250 | 0.4 | 1.0 |
| *Polyporus squamosus* | `squamosus` | 235 | 2043 | 3750 | 5.0 | 3.0 |

Spatial calibration: 8.7 um/px (derived from *F. fomentarius* SEM data; [Klemm et al., 2024](https://doi.org/10.1371/journal.pone.0303122)).

## Installation

```bash
# Clone the repository
git clone https://github.com/mahumadad/FungalMorphoSpace.git
cd FungalMorphoSpace

# Install dependencies
pip install -r requirements.txt

# Or install as a package
pip install -e .

# Verify
python scripts/test_imports.py
```

### Requirements

- Python 3.9+
- NumPy, SciPy, Matplotlib, scikit-image, pandas, tqdm, PyYAML
- ~4 GB RAM for 1024x1024 grids

## Quick Start

```bash
# Validate all 3 species (1 run each, fast)
python scripts/run_integrated_validation.py --species all --n_runs 1

# Single species with higher resolution
python scripts/run_integrated_validation.py --species fomes --n_runs 5 --grid 1024

# Smoke test (verifies full pipeline)
python scripts/smoke_test.py --grid 256
```

### Python API

```python
from fungalmorphospace.core import TuringSimulator, GiererMeinhardtKinetics
from fungalmorphospace.analysis import TopologyAnalyzer

# Setup and run simulation
kinetics = GiererMeinhardtKinetics(rho=0.2, a=0.1, b=1.0)
sim = TuringSimulator(kinetics_model=kinetics, D_u=0.1, D_v=15.0, grid_size=512)
sim.initialize_fields(noise=0.1)
sim.run(T_target=5.0)

# Analyze pattern
analyzer = TopologyAnalyzer(sim.u, dx=1.0)
metrics = analyzer.compute_all_metrics()
print(f"Wavelength: {metrics['wavelength_fft']:.1f} px")
```

## Project Structure

```
FungalMorphoSpace/
├── src/fungalmorphospace/          # Main Python package
│   ├── core/                       # Turing simulator engine
│   │   ├── turing_simulator.py     # RD solver (explicit Euler, periodic BC)
│   │   └── kinetics.py             # Reaction kinetics models
│   ├── analysis/                   # Pattern analysis
│   │   ├── topology_analyzer.py    # Morphometric + topological metrics
│   │   └── visualization.py        # Publication-quality figures
│   ├── runners/                    # Pipeline orchestration
│   │   ├── integrated_validation.py # Main validation runner
│   │   └── species_database.py     # Species parameter loading
│   ├── contracts/                  # Output schema enforcement
│   └── utils/                      # Sensitivity analysis, cleanup
├── scripts/                        # CLI entry points
├── data/species_data.json          # Calibrated species parameters
├── docs/                           # Documentation
├── tests/                          # Unit and smoke tests
├── paper/                          # SoftwareX manuscript (LaTeX)
└── results/                        # Generated outputs (gitignored)
```

## Documentation

- [API Reference](docs/API.md)
- [Usage Examples](docs/EXAMPLES.md)
- [Output Contract](docs/OUTPUT_CONTRACT.md)
- [Changelog](CHANGELOG.md)
- [Contributing](CONTRIBUTING.md)

Additional thesis-specific documentation is available in `docs/archive/`.

## Scientific Background

The software implements the Gierer-Meinhardt reaction-diffusion system:

```
du/dt = D_u * laplacian(u) + rho * (u^2/v - u + a)
dv/dt = D_v * laplacian(v) + rho * (u^2 - b*v)
```

A key insight is that effective diffusion D_eff = D/rho jointly determines pattern wavelength, explaining cases where increasing the diffusion ratio *decreases* wavelength when reaction intensity increases simultaneously.

**Important note:** Parameters are currently calibrated empirically (fitted to observed morphologies). The software serves as an exploration and hypothesis-generation tool rather than an ab initio predictor. See `docs/archive/RUTAS_ANCLAJE_NO_TAUTOLOGICO.md` for a discussion of non-tautological anchoring strategies, and `docs/ESTRATEGIA_TESIS.md` for the publication roadmap.

**Grid resolution matters.** Large-wavelength species need a domain that fits several wavelengths. *Polyporus squamosus* defaults to a 1024 grid (where it forms ~28 genuine spots); forcing a smaller grid (e.g. `--grid 512`) collapses it to a single domain-scale blob (1 spot) that is *not* a genuine periodic pattern. The validator now flags this: results report `pattern_genuine` (≥4 resolved spots) and `under_resolved` (fewer than ~3 measured wavelengths fit the domain), and `validation_pass` combines the wavelength QC with the genuine-pattern gate. Keep *squamosus* at its default grid or larger. Note: even at adequate resolution the model's *squamosus* wavelength (~188 px) undershoots the biological target (235 px ≈ 2000 µm) by ~18% — a known limitation, not corrected by tuning the target.

## Citation

If you use this software in your research, please cite:

```bibtex
@software{ahumada2026fungalmorphospace,
  author       = {Ahumada Dur{\'a}n, Mario},
  title        = {FungalMorphoSpace: Reaction-diffusion morphospace
                  simulator for fungal hymenophores},
  year         = {2026},
  version      = {0.7.3.post6},
  url          = {https://github.com/mahumadad/FungalMorphoSpace}
}
```

## License

**Dual license:**

- **Academic/Research:** [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/) — Free for non-commercial use with attribution.
- **Commercial:** Contact the author for licensing terms.

## References

1. Turing, A.M. (1952). The chemical basis of morphogenesis. *Phil. Trans. R. Soc. B*, 237, 37-72.
2. Gierer, A. & Meinhardt, H. (1972). A theory of biological pattern formation. *Kybernetik*, 12, 30-39.
3. Kuhar, F. et al. (2022). Pattern formation features might explain homoplasy: fertile surfaces in higher fungi as an example. *Theory in Biosciences*, 141(1), 1-11.
4. Klemm, D. et al. (2024). Hierarchical structure of *Fomes fomentarius*. *PLOS ONE*.

---

(c) 2025-2026 Mario Ahumada Duran
