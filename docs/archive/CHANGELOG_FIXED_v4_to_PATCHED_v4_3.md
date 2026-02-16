# CHANGELOG (auditoría) — FIXED_v4 → PATCHED_v4_3
Fecha: 2026-01-01 17:47:31

Este documento es un **diff narrado** desde `FungalMorphoSpace_FIXED_v4.zip` hacia el paquete completo `PATCHED_v4_3` (v0.4.3).

---

## Cambios clave (alto impacto)
1) **Refactor a paquete Python (`src/fungalmorphospace/`)**  
   - Reduce fragilidad por imports, facilita instalación editable y publicación futura.
2) **Blindaje de métricas**  
   - Spots/componentes por CC + centroides (sin sobreconteo por máximos locales).
   - λ(FFT) radial + QC; λ(autocorr) radial robusta.
3) **Comparabilidad temporal**  
   - Control por `T_target` (steps derivados de `dt_final`).
4) **Fix crítico v4.3: shadowing de imports**  
   - Scripts priorizan `./src` para evitar que una instalación vieja cause errores (ej. `compute_all_metrics`).
5) **Fix v4.3: `plot_scaling_law.py` sin hardcodes**  
   - Lee λ(px) desde `results/validation_summary.csv` y escala biológica desde `data/species_data.json`.

---

## Archivos modificados (comparado con FIXED_v4) — 9
- README.md
- data/species_data.json
- docs/GUIA_EXPERIMENTOS_BLINDAJE.md
- docs/INFORME_TECNICO_TESIS.md
- scripts/plot_scaling_law.py
- scripts/run_integrated_validation.py
- scripts/run_parallel_validation.py
- scripts/run_robustness_analysis.py

## Archivos agregados — 15
- CHANGELOG.md
- LICENSE
- data/_archive/species_data.json.old
- docs/INFORME_CONTEXTO_PROBLEMAS_RESUELTOS.md
- pyproject.toml
- src/fungalmorphospace/__init__.py
- src/fungalmorphospace/__main__.py
- src/fungalmorphospace/analysis/__init__.py
- src/fungalmorphospace/analysis/topology_analyzer.py
- src/fungalmorphospace/analysis/visualization.py
- src/fungalmorphospace/core/__init__.py
- src/fungalmorphospace/core/kinetics.py
- src/fungalmorphospace/core/turing_simulator.py
- src/fungalmorphospace/utils/__init__.py
- src/fungalmorphospace/utils/sensitivity_analysis.py

## Archivos removidos — 9
- data/species_data.json.old
- src/analysis/__init__.py
- src/analysis/topology_analyzer.py
- src/analysis/visualization.py
- src/core/__init__.py
- src/core/kinetics.py
- src/core/turing_simulator.py
- src/utils/__init__.py
- src/utils/sensitivity_analysis.py

---

## Nota sobre archivos de cache
En FIXED_v4 y/o el entorno de empaquetado pueden aparecer `__pycache__/` y `.pyc`.  
En `PATCHED_v4_3` **se excluyen del ZIP final** (ver limpieza pre-empaquetado).
