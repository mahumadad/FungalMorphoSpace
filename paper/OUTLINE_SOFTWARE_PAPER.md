# Software paper — Outline (FungalMorphoSpace)

**Target: SoftwareX** (Elsevier). Razón: acepta CC-BY-NC (preserva la licencia comercial dual);
es LaTeX; revisa software + significancia. JOSS descartado: exige licencia OSI-abierta, que
rompería la opción comercial.

Encuadre HONESTO: herramienta de exploración/generación de hipótesis. NO predicción ab initio,
NO mecanismo nuevo. El mecanismo nodo-inmóvil/equal-D ya está publicado (Nesterenko 2017,
Korvasová/Klika 2015, Marciniak-Czochra 2016) → se cita como fundamento.

## Estructura SoftwareX

1. **Motivation and significance** — qué hace, para quién; hipótesis de Kuhar; nicho del software.
2. **Software description** — arquitectura (core/analysis/runners/contracts/three_node);
   funcionalidades; 3 cinéticas; base de 3 especies; gates de validación; 16 tests.
   - Tabla 1: especies calibradas (con nota de calibración empírica).
   - Code metadata table (SoftwareX obligatoria).
3. **Illustrative examples** — (a) morfo-espacio 3 especies; (b) 3-nodos equal-D + ley √D
   (Fig. `figure_scale_law.png`); (c) diagnóstico de estabilidad lineal.
4. **Impact** — reproducibilidad, docencia, base para extensiones; sin claim predictivo.
5. **Limitations** — calibración empírica/tautológica; mecanismo no novedoso; caveat RD-ODE.
6. **Conclusions**.
7. Obligatorias: Data Availability (Zenodo DOI), AI disclosure, CRediT, CoI, Funding.

## Claims permitidos vs prohibidos

| ✅ | ❌ |
|---|---|
| herramienta de exploración/hipótesis | predice tamaño de poro ab initio |
| reproduce cualitativamente el morfo-espacio | explica la morfogénesis fúngica real |
| el ratio extremo es innecesario (mecanismo conocido) | descubrimos mecanismo nuevo |
| calibración empírica (limitación) | validado contra biología independiente |

## Figuras/tablas
- Fig 1: COMPREHENSIVE_COMPARISON (3 especies).
- Fig 2: figure_scale_law.png (ley √D, equal-D).
- Tabla 1: especies calibradas. Tabla code-metadata (SoftwareX).
