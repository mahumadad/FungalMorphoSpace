# Informe de Contexto Técnico Definitivo

## Objetivo del documento

Este informe consolida el contexto técnico del proyecto **FungalMorphoSpace** y la decisión metodológica clave (selección de especies, calibración de parámetros y criterios de validación). Se usa para mantener coherencia en el repositorio y como insumo directo para redacción académica (tesis/paper).

---

## 1. Resumen del proyecto

**FungalMorphoSpace** es una infraestructura computacional para explorar y validar la hipótesis de que ciertos patrones morfológicos de himenóforos (superficies fértiles con poros) pueden explicarse mediante **auto-organización físico-química**, en particular por mecanismos tipo **reacción–difusión** (patrones de Turing / LALIP).

El objetivo práctico del repositorio es:

- simular patrones en mallas 2D con cinéticas tipo Schnakenberg/Gierer–Meinhardt u otras,
- cuantificar **longitud de onda (λ)**, **topología** (Euler χ, holes) y métricas de textura,
- mapear resultados a magnitudes biológicas (densidad de poros, tamaño de poro) bajo una calibración explícita (µm/px),
- y construir un “morfoespacio” donde parámetros físicos actúan como variables causales.

---

## 2. El problema principal que se resolvió

La versión previa del proyecto usaba *Trametes versicolor* como “límite denso”. Esa elección generó una ambigüedad crítica:

1. **Ambigüedad en literatura:** *T. versicolor* tiene reportes variados de densidad (por ejemplo 3–5 poros/mm), pero no se comporta como extremo limpio sin solapamiento con especies intermedias.
2. **Ambigüedad morfométrica vs. modelo:** el rango real de *Trametes* se solapaba con especies intermedias, lo que hacía frágil la narrativa de escalamiento y la calibración de parámetros.
3. **Riesgo argumental:** si el “límite denso” no es extremo, cualquier conclusión sobre restricciones morfogenéticas pierde potencia, porque el morfoespacio queda mal anclado.

---

## 3. La solución metodológica adoptada

Se reemplazó *Trametes* por un esquema de **tres especies ancla** con separación clara de escala:

- **Baseline:** *Fomes fomentarius*.
- **Intermedio:** *Lentinus brumalis* (sin. *Polyporus brumalis*).
- **Extremo grande:** *Polyporus squamosus*.

Se mantiene *Trametes versicolor* en el dataset como **EXCLUDED**, con alias preservados, pero no participa en validación canónica.

---

## 4. Racional físico (por qué el ajuste de ρ es correcto)

Además de la razón de difusividades \(D_v/D_u\), existe un parámetro cinético \(\rho\) que escala la dinámica temporal. El argumento clave es el reescalamiento:

\[
\tau = \rho t
\]

Bajo este cambio, la difusión efectiva se comporta, en términos prácticos, como:

\[
D_{\mathrm{ef}} \propto \frac{D}{\rho}
\]

Consecuencia: aumentar \(\rho\) es geométricamente equivalente a disminuir la difusión efectiva, lo que permite ajustar la **densidad de spots/poros** sin introducir contradicciones internas. Esto habilita calibrar densidad (pores/mm) preservando topología.

---

## 5. Calibración vigente (v0.7.3)

### 5.1 Especies canónicas y parámetros

| Especie | Key | D_v/D_u | b | ρ | Grid | Objetivo biológico |
|---|---:|---:|---:|---:|---:|---|
| *Fomes fomentarius* | fomes | 150 | 1.0 | 0.2 | 512 | baseline |
| *Lentinus brumalis* (sin. *Polyporus brumalis*) | brumalis | 250 | 1.0 | 0.4 | 512 | 2–4 poros/mm |
| *Polyporus squamosus* | squamosus | 3750 | 3.0 | 5.0 | 1024 | 1–2 poros/mm |

### 5.2 Validación de brumalis (corrida de referencia)

Para fijar la ancla intermedia se calibró \(\rho=0.4\) y se corrió n=5 semillas:

- λ(best): media ≈ 32.6 px (SD ≈ 1.6; CV ≈ 4.8%)
- densidad predicha: media ≈ 3.54 poros/mm (SD ≈ 0.18; CV ≈ 5.1%)
- topología: régimen de spots estable con variación moderada entre semillas



### 5.3 Calibración alternativa (tesis): *L. brumalis* meso-poroide

Además del set vigente del software (ρ=0.4, λ≈33 px), se documentó una trayectoria de calibración orientada a ubicar a *Lentinus brumalis* como transición **meso-poroide** (objetivo λ≈60–70 px) entre *F. fomentarius* y *P. squamosus*.

- Run final de tesis (candidato): D_v/D_u=1000, ρ=0.15, b=1.0, con λ≈62.9 px (grid 1024) y alta población de spots.
- Bitácora completa y discusión física: `docs/CALIBRACION_BRUMALIS_TUNING_JOURNEY.md`.

**Nota contractual:** este candidato no reemplaza automáticamente los defaults del software; adoptarlo como preset canónico requiere re-baselinear validaciones y documentar el cambio de “gold standard”.
---

## 6. Outputs y trazabilidad

La salida del proyecto se organiza con un **contrato de outputs** estable, con dos capas:

1. **Canónico (último ALL):** sobre-escribe los CSV/figuras canónicas (estado actual).
2. **Por experimento (exp_id):** snapshot por corrida + ledger append-only.

En v0.7.0 se incorporan snapshots por corrida (human CSV, machine CSV, JSON) y un ledger `_index.csv` en:

- `results/tables/singles/` para corridas single,
- `results/tables/all/` para corridas ALL,
- con copia de figura por corrida en `results/figures/all/` cuando aplica.

---

## 7. Inserto bibliográfico (2025) y posicionamiento

Se agrega como contexto (estado del arte) el paper:

**“Deep Ensemble Learning and Explainable AI for Multi-Class Classification of Earthstar Fungal Species” (Biology, 2025).**

Motivo de inclusión: refuerza que la taxonomía computacional moderna (DL + XAI) se apoya en rasgos morfológicos y sufre por “small data”. FungalMorphoSpace se posiciona como un generador de datos sintéticos físicamente informados, útil para data augmentation no trivial, análisis explicable con variables causales (parámetros) y mitigación de desbalance de clases.

---

## 8. Prompt corto reutilizable (para LaTeX)

Puedes copiar lo siguiente como prompt de trabajo:

**Redacta en LaTeX (sin primera persona) una sección de métodos y resultados para una tesis de biofísica que valida patrones de Turing en himenóforos. Usa como anclas tres especies: Fomes fomentarius (D_v/D_u=150, b=1, ρ=0.2), Lentinus brumalis (D_v/D_u=250, b=1, ρ=0.4, validado n=5 con λ≈33 px y densidad≈3.54 poros/mm), y Polyporus squamosus (D_v/D_u=3750, b=3, ρ=5). Explica el reescalamiento temporal τ=ρt y la interpretación D_ef∝D/ρ. Describe outputs canónicos y snapshots por exp_id con ledger append-only.**
