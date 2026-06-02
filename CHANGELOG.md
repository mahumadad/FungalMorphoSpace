# CHANGELOG — FungalMorphoSpace

## Unreleased (2026-06-02)

### Fixed — validation integrity (audit follow-up)
- **Genuine-pattern gate:** the validator now flags a result as a genuine periodic
  pattern only when it resolves at least `MIN_GENUINE_SPOTS` (=4) components. A single
  domain-scale blob (e.g. *Polyporus squamosus* under-resolved on a 512 grid) previously
  passed the FFT QC and fell inside a wide density band, producing a spurious
  "validation". New per-result fields: `pattern_genuine`, `validation_pass`.
- **Sub-resolution warning:** the validator now computes `wavelengths_in_domain` and warns
  (`under_resolved=True`) when fewer than `MIN_WAVELENGTHS_IN_DOMAIN` (=3) heuristic
  wavelengths fit the domain, i.e. when the grid is too small to resolve the pattern.
  *squamosus* on a 512 grid fits only ~1.3 wavelengths and is now flagged.
- Both flags surface in the verbose output and the human-readable CSVs. The canonical
  machine schema is unchanged (no contract bump).
- **Density gating:** `density_predicted` is now reported only for a genuine, QC-passing
  pattern (`validation_pass`), so a degenerate/under-resolved result cannot spuriously fall
  inside the observed density band.
- **Mid-run stability guard:** the solver re-checks for non-finite fields during the run
  (the timestep is fixed after init) and stops cleanly with a flag instead of emitting
  NaN-filled output. No-op for stable runs.
- *Attempted and reverted:* a "no-peak spectrum" QC criterion (reject FFT peaks at the
  frequency floor). It also rejected genuine large-wavelength patterns (*P. squamosus* at
  1024: a real 28-spot pattern fell back to a worse autocorrelation estimate), so it was
  reverted. The spot-count gate above covers the degenerate case robustly.

### Fixed — references
- Corrected the central **Kuhar et al. (2022)** citation, which carried a non-existent DOI
  (`10.1007/s12064-022-00380-8`) and incorrect title/authors/pages. The verified reference
  is `10.1007/s12064-022-00363-z`, *Theory in Biosciences* 141(1):1-11 (Kuhar, Terzzoli,
  Nouhra, Robledo, Mercker). Fixed in `paper/references.bib` and `README.md`.
- Added the standard modelling review **Davidson et al. (2011)**, *IMA Fungus* 2(1):33-37
  (`10.5598/imafungus.2011.02.01.06`).

## v0.7.3.post6 (2026-01-03)

### Docs
- Added an explicit tuning log and thesis-oriented calibration narrative for *Lentinus brumalis* ("The Tuning Journey"): `docs/CALIBRACION_BRUMALIS_TUNING_JOURNEY.md`.
- Updated thesis technical report to reference the tuning journey and to embed a LaTeX-ready summary table for the meso-poroide candidate calibration.
- Updated the hardening experiment guide with reproducible CLI commands for the key tuning runs (A–E) without modifying defaults.
- Updated the technical context handoff document to record the existence of a thesis candidate calibration (while clarifying contractual implications if defaults are changed).

## v0.7.3.post5 (2026-01-03)

### Fixed
- Fixed a `NameError` during snapshot export caused by an inconsistent pandas alias (`_pd` vs `pd`) in `export_experiment_bundle()`. Snapshot human CSV export is now reliable.
- Eliminated a pandas `FutureWarning` during ledger appends by avoiding `DataFrame` concatenation against empty/all-NA frames (forward-compatible with upcoming pandas behavior changes).
- Fixed an invalid `pyproject.toml` line break that could prevent TOML parsing in packaging/tooling.


## v0.7.3.post4 — 2026-01-03

### Dev hygiene: cleanup of leftover temp folders
- CLEAN: `scripts/smoke_test.py` now removes known dev temp folders in repo root (when they look like FMS output) by default.
  - Disable with `--no-clean-dev`.
- ADD: `scripts/clean_dev_artifacts.py` helper to manually remove `tmp_contract_test/`, `tmp_check/`, `tmp_check2/` (conservative deletion).
- ADD: `.gitignore` entries for the above temp folders and common Python caches.


## v0.7.3.post3 — 2026-01-03

### UX: Progress bars with elapsed time / ETA
- UX: integrated validator now uses `tqdm` to show progress with elapsed time and ETA during simulation stepping (per run) and across replicate runs (`--n_runs > 1`) when not in `--quiet` mode.
- QUIET: `--quiet` suppresses progress bars and keeps output minimal (useful for smoke tests and scripting).

## v0.7.3.post2 — 2026-01-03

### 🔴 CRITICAL FIX: Canonical report generation (smoke test)
- FIX: `IntegratedSimulationRunner.create_comparison_report()` now builds canonical CSV/JSON from the runner’s real result schema (`species_name`, `D_v_D_u`, `wavelength_best_px`, etc.).
- FIX: Canonical machine CSV now writes the contract column names (via `enforce_machine_schema`) without silently dropping data due to mismatched keys.
- FIX: Canonical JSON is now guaranteed JSON-serializable (no raw numpy arrays), eliminating the failure mode where an exception triggered placeholder outputs and made the smoke test report missing species.
- HARDEN: Fail-open behavior is now stage-isolated; JSON/figure failures no longer overwrite already-written canonical CSVs.
- NOTE: `scripts/smoke_test.py` deletes the temp output directory by default; use `--keep` to inspect generated artifacts after the run.

## v0.7.3.post1 — 2026-01-03

- CONTRACT: hardened validator report writer (fail-open). Canonical outputs are always created; if a stage fails, placeholders are emitted and full tracebacks are appended to `results/logs/report_errors.log`.
- FIX: snapshot exporter now correctly calls `enforce_machine_schema` (previously referenced `_enforce_machine_schema`), preventing missing per-experiment snapshot CSVs.
- VIS: per-run pattern analysis PNGs now include a suptitle with species name and parameter provenance.

## v0.7.3 — 2026-01-03

- VIS: added `--include_fields` option to the integrated validator. When enabled, the pipeline also writes two additional PNGs per species: raw activator field (u) and raw inhibitor field (v), both annotated with species name and the parameter set used.
- VIS: `COMPREHENSIVE_COMPARISON.png` now includes parameter annotations in each species panel (`D_v/D_u`, `b`, `ρ`, `grid`, `T_target`) for easier provenance.
- CONTRACT: output contract updated to document the optional `figures/fields/` artifacts and per-experiment field images under `figures/all/` and `figures/singles/`.

## v0.7.2 — 2026-01-03

- CONTRACT: hardened output contract via `fungalmorphospace.contracts.output_contract` (single source of truth for paths + machine CSV schema).
- CONTRACT: machine-readable CSV schema is now strict (stable columns + ordering); exporters enforce schema before writing.
- DOCS/TESTS: updated `docs/OUTPUT_CONTRACT.md` to enumerate machine columns and ledger core schema; smoke test reads schema from contract.

## v0.7.1 — 2026-01-03

- FIX: `IntegratedSimulationRunner.export_experiment_bundle()` and `_append_ledger_row()` were unintentionally defined at module scope; moved into the class to restore the snapshot/ledger workflow.
- FIX: `scripts/smoke_test.py` now verifies per-experiment snapshots + append-only ledgers (ALL + single brumalis).
- PARAMS: `squamosus` `dt_initial` lowered to `5e-05` for stability on 1024² grids.
- SCI: `scripts/test_laminillas.py` aligned to core Gierer–Meinhardt kinetics and includes a conservative anisotropy-aware CFL guard.
- LOGS: clarified that `lambda_theoretical` is a heuristic scale proxy (not a full dispersion prediction).
- DOCS: corrected Brumalis diffusion-effective values (`D_v/ρ=625`) and updated dt table.
- META: removed placeholder email from `pyproject.toml`.

## v0.7.0 — 2026-01-02

- Added per-experiment snapshots and append-only ledgers for integrated validation:
  - `results/tables/singles/validation_summary_<timestamp>_single_<species>_<exp_id>.csv` + machine + JSON
  - `results/tables/all/validation_summary_<timestamp>_all_all_<exp_id>.csv` + machine + JSON
  - `results/figures/all/COMPREHENSIVE_COMPARISON_<timestamp>_all_all_<exp_id>.png`
  - Ledgers: `results/tables/singles/_index.csv` and `results/tables/all/_index.csv`
- Canonical outputs remain unchanged and continue to be overwritten by the latest `--species all` run.
- Added repo-root wrapper `run_integrated_validation.py` (fixes common path error).
- Updated species database:
  - `brumalis`: ρ=0.4 (validated n=5), scientific name set to *Lentinus brumalis* with synonyms preserved as aliases.
  - `trametes` remains EXCLUDED (documented), with aliases retained.
- Hardened `scripts/smoke_test.py` to verify snapshot + ledger outputs.
- Documentation refresh (output contract, methods/technical thesis notes) aligned to the new parameters and workflow.

