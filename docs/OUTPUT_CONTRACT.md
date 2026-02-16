# Output Contract (Canonical)

**Version:** 0.7.3

**Source of truth:** `src/fungalmorphospace/contracts/output_contract.py` (paths + machine CSV schema).

This document defines the **canonical, stable** file and directory layout produced by the FungalMorphoSpace CLI pipelines.

## Scope

Applies to:
- `python3 scripts/run_integrated_validation.py --species all ...`
- `python3 scripts/run_parallel_validation.py --species <key> ...`

Default output root: `results/` (can be overridden with `--output <dir>`)

## Species Keys

| Canonical Key | Aliases | Scientific Name |
|---------------|---------|-----------------|
| `fomes` | `fomes_fomentarius` | *Fomes fomentarius* |
| `brumalis` | `polyporus_brumalis`, `lentinus_brumalis` | *Lentinus brumalis* (syn. *Polyporus brumalis*) |
| `squamosus` | `polyporus_squamosus` | *Polyporus squamosus* |

Both canonical keys and aliases are accepted by all scripts.

## Canonical Directory Tree

Under `<OUTPUT_ROOT>` the pipelines create:

```
<OUTPUT_ROOT>/
  patterns/
  figures/
  tables/
  metrics/
  logs/
  analysis/
  sensitivity/
```

### 1) `patterns/`
Per-run analysis figures and canonical "last-run" copies.

**Per-run figures:**
```
<OUTPUT_ROOT>/patterns/<species_key>/run_<NN>_seed<SEED>_<YYYYMMDD_HHMMSS>.png
```

**Optional per-run field images (enabled with `--include_fields`):**
```
<OUTPUT_ROOT>/patterns/<species_key>/run_<NN>_seed<SEED>_<YYYYMMDD_HHMMSS>_activator_u.png
<OUTPUT_ROOT>/patterns/<species_key>/run_<NN>_seed<SEED>_<YYYYMMDD_HHMMSS>_inhibitor_v.png
```

**Canonical copy:**
```
<OUTPUT_ROOT>/patterns/<Full_Name_With_Underscores>.png
```

### 2) `figures/`
Cross-species or publication-ready figures.

**Canonical (only generated with `--species all`):**
```
<OUTPUT_ROOT>/figures/COMPREHENSIVE_COMPARISON.png
```

**Optional per-species field images (enabled with `--include_fields`):**
```
<OUTPUT_ROOT>/figures/fields/<species_key>_activator_u.png
<OUTPUT_ROOT>/figures/fields/<species_key>_inhibitor_v.png
```
These PNGs include the species name and the parameter set used (e.g., `D_v/D_u`, `b`, `ρ`, grid, and `T_target`).

### 3) `tables/`
Tabular outputs.

**Canonical (only generated with `--species all`):**
```
<OUTPUT_ROOT>/tables/validation_summary.csv         # Human-readable
<OUTPUT_ROOT>/tables/validation_summary_machine.csv # Machine-readable
<OUTPUT_ROOT>/tables/validation_summary.json        # JSON summary
```



### Per-experiment snapshots and ledgers (v0.7.0+)

The integrated validation CLI (`scripts/run_integrated_validation.py`) now emits **optional per-experiment snapshots**
(and an **append-only ledger**) without changing the canonical tables.

- Use `--exp_id <ID>` to set the ID explicitly (recommended for reproducibility).
- If omitted, an `exp_id` is auto-generated.

**Single-species snapshots (`--species <key>`):**
```
<OUTPUT_ROOT>/tables/singles/validation_summary_<exp_id>.csv
<OUTPUT_ROOT>/tables/singles/validation_summary_machine_<exp_id>.csv
<OUTPUT_ROOT>/tables/singles/validation_summary_<exp_id>.json
<OUTPUT_ROOT>/tables/singles/_index.csv                # append-only ledger (one row per experiment)
```

**ALL snapshots (`--species all`):**
```
<OUTPUT_ROOT>/tables/all/validation_summary_<exp_id>.csv
<OUTPUT_ROOT>/tables/all/validation_summary_machine_<exp_id>.csv
<OUTPUT_ROOT>/tables/all/validation_summary_<exp_id>.json
<OUTPUT_ROOT>/tables/all/_index.csv                    # append-only ledger (one row per experiment)

<OUTPUT_ROOT>/figures/all/COMPREHENSIVE_COMPARISON_<exp_id>.png

# Optional (when `--include_fields`)
<OUTPUT_ROOT>/figures/all/<exp_id>_<species_key>_activator_u.png
<OUTPUT_ROOT>/figures/all/<exp_id>_<species_key>_inhibitor_v.png
```

**Single-species snapshots (`--species <key>`, optional when `--include_fields`):**
```
<OUTPUT_ROOT>/figures/singles/<exp_id>_<species_key>_activator_u.png
<OUTPUT_ROOT>/figures/singles/<exp_id>_<species_key>_inhibitor_v.png
```

Snapshot CSVs include an `exp_id` column; the canonical CSVs do not.

### Ledgers (`_index.csv`) schema

Ledgers are **append-only** CSV indexes (with de-duplication by `exp_id`, keeping the most recent row).  
They are the long-term, audit-friendly record of what ran, when, and with which parameters.

**Core stable columns (guaranteed):**

| Column | Type | Description |
|--------|------|-------------|
| `exp_id` | str | Experiment identifier (primary key; used for de-duplication). |
| `mode` | str | `single` or `all`. |
| `generated_at` | str | ISO-8601 timestamp of export. |
| `version` | str | Package version used for the run. |
| `n_species` | int | Number of species included in this experiment bundle. |
| `species_json` | str | Compact JSON mapping (species name → `{species_key, D_v_D_u, b, rho, grid_size}`). |

**Optional “wide” convenience columns (may expand):**

For each species in the experiment, the exporter may add columns following the pattern:

- `<species_key>_rho`
- `<species_key>_D_v_D_u`
- `<species_key>_b`
- `<species_key>_wavelength_best_px`
- `<species_key>_density_predicted_pores_per_mm`
- `<species_key>_grid_size`
- `<species_key>_dt_final`

These columns are intended for quick spreadsheet-style comparisons. They are **not** the primary durable record; `species_json` + the per-experiment snapshots are.



**Parallel-safe single-species (run_parallel_validation.py):**
```
<OUTPUT_ROOT>/tables/<species_key>_validation.csv
```

### Machine-readable CSV columns

The `validation_summary_machine.csv` contains:

| Column | Type | Description |
|--------|------|-------------|
| `species` | str | Scientific name (human-readable). |
| `species_key` | str | Canonical key (e.g., `fomes`, `brumalis`, `squamosus`). |
| `D_v_D_u` | float | Diffusion ratio. |
| `b` | float | Kinetic parameter (Gierer–Meinhardt). |
| `rho` | float | Temporal rescaling factor (metabolic scaling). |
| `grid_size` | int | Simulation grid size (px). |
| `wavelength_best_px` | float | Best wavelength estimate in pixels (FFT if QC pass, else autocorr). |
| `wavelength_method` | str | Method used for `wavelength_best_px` (`fft` or `autocorr`). |
| `wavelength_qc_pass` | bool | QC flag for the chosen wavelength method. |
| `wavelength_fft_px` | float | FFT-based wavelength estimate (px). |
| `wavelength_fft_qc_pass` | bool | QC flag for FFT estimate. |
| `wavelength_fft_peak_ratio` | float | FFT spectral peak ratio (diagnostic). |
| `wavelength_autocorr_px` | float | Autocorrelation-based wavelength estimate (px). |
| `wavelength_autocorr_qc_pass` | bool | QC flag for autocorrelation estimate. |
| `spots` | int | Count of detected spots (pores). |
| `holes` | int | Count of holes (topological). |
| `euler_chi` | int | Euler characteristic (spots − holes). |
| `density_predicted_pores_per_mm` | float | Predicted pore density (pores/mm). |
| `density_observed_min_pores_per_mm` | float | Lower bound of observed biological range (pores/mm). |
| `density_observed_max_pores_per_mm` | float | Upper bound of observed biological range (pores/mm). |
| `dt_requested` | float | Requested initial time step. |
| `dt_final` | float | Final time step after CFL adjustment. |
| `dt_adjusted` | bool | Whether CFL forced a reduction of `dt_requested`. |
| `dt_max_diffusion` | float | Maximum stable dt from diffusion CFL. |
| `safety_factor` | float | Safety factor applied to diffusion CFL. |
| `T_target` | float | Requested simulation horizon. |
| `T_actual` | float | Actual simulated time (may differ slightly). |
| `steps_computed` | int | Computed steps (`ceil(T_target/dt_final)`). |
| `had_warning` | bool | Whether warnings were raised during the run. |
| `pattern_path` | str | Path to representative pattern file for the species/run. |

**Snapshots:** per-experiment machine CSVs add `exp_id` as the first column; all other columns are identical and ordered.


## Legacy Copies (backward compatibility)

When using `--species all`, legacy copies are created at:

```
<OUTPUT_ROOT>/COMPREHENSIVE_COMPARISON.png
<OUTPUT_ROOT>/validation_summary.csv
```

These are copies of the canonical files under `figures/` and `tables/`.

## Important Notes

1. **`--species all` required for comparison reports**: The comprehensive figure and validation summary tables are only generated when running `--species all`. Single-species runs generate only pattern files and species-specific CSVs.

2. **Time and stability metadata**: For every simulation, the system records:
   - `T_target` and `T_actual`
   - `dt_final` (after CFL stability adjustment)
   - `steps_computed = ceil(T_target / dt_final)`

3. **Wavelength sources**: The machine CSV contains both FFT and autocorrelation wavelength estimates. Use `wavelength_best_px` for the recommended value (FFT if QC pass, else autocorr).
