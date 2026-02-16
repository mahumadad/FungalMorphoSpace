# Test plan ideal (blindaje científico + ingeniería)

Este documento define el set **ideal** de pruebas para elevar el proyecto a un estándar publicable (paper + repo público).

## 1) Reproducibilidad numérica

### 1.1 Replicación estadística por especie
- Ejecutar **n=20** replicas por especie (mismas condiciones salvo seed).
- Reportar:
  - media±sd de `λ(best)`
  - media±sd de `spots`, `holes`, `euler_chi`
  - tasa de QC-pass para FFT y autocorr

**Criterio de aprobación (propuesto):**
- CV(λ) < 10% para especies con patrón estable.
- Si FFT QC falla > 50% en una especie, justificar por régimen extremo y usar autocorr como primario.

### 1.2 Invariancia a resolución (grid y dx)
- Repetir para `grid_size` ∈ {256, 512, 1024}.
- Verificar que el patrón converge en métricas adimensionales y que la conversión a µm/px preserva el orden entre especies.

## 2) Validez de medición de λ

### 2.1 Banco de patrones sintéticos
Construir un set de patrones con λ conocido:
- seno/coseno 2D
- retículas hexagonales generadas
- patrones con ruido controlado

**Expectativa:** error relativo < 5%.

### 2.2 Robustez FFT vs autocorr
- Forzar casos extremos (λ muy grande o muy pequeño vs el tamaño de imagen).
- Confirmar:
  - FFT usa QC y se descarta cuando el pico no es dominante.
  - el selector `select_wavelength_px()` elige el método correcto.

## 3) Contrato de outputs

### 3.1 Estructura canónica
- El pipeline debe crear siempre:
  - `patterns/ figures/ tables/ metrics/ logs/ analysis/ sensitivity/`
- Verificar que los archivos canónicos existan tras el run:
  - `figures/COMPREHENSIVE_COMPARISON.png`
  - `tables/validation_summary_machine.csv`

### 3.2 Compatibilidad hacia atrás
- Verificar que:
  - `<OUTPUT_ROOT>/validation_summary.csv` y `<OUTPUT_ROOT>/COMPREHENSIVE_COMPARISON.png` existan como copias.

## 4) Blindaje de tiempo físico (M-04)

### 4.1 Invarianza temporal
- En cada run, verificar:
  - `steps_computed = ceil(T_target / dt_final)`
  - `T_actual >= T_target` (por el ceil)

### 4.2 Test de regresión
- Si un cambio en `_check_stability()` o `dt` altera `dt_final`, el pipeline debe seguir igualando tiempo físico vía steps.

## 5) Ingeniería (calidad de repo)

### 5.1 Imports y packaging
- Scripts no deben importar entre sí.
- La lógica de pipelines vive en `src/fungalmorphospace/runners/`.

### 5.2 Lint / formatting (opcional)
- `ruff`/`black` para consistencia.

### 5.3 Performance budget
- Smoke tests < 30 s.
- Runs de paper (n=20, grid=512) pueden demorar, pero deben ser reproducibles.

## 6) (Opcional) TDA y taxonomía
Cuando implementes TDA:
- tests de invariancia topológica (deformaciones suaves no cambian Betti/persistencia)
- tests de separación de clases (poros vs láminas) con métricas de distancia en espacio de persistencia.
