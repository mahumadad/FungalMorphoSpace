# GUÍA COMPLETA DEL PROYECTO FUNGALMORPHOSPACE
## Manual Operativo y de Referencia

**Autor:** Mario Ahumada Durán  
**Proyecto:** FungalMorphoSpace v0.5.1  
**Fecha:** Enero 2026

## Hoja de ruta

Para planificación estratégica (papers, derivados, productos y dependencias), ver `docs/HOJA_DE_RUTA.md` (también disponible como `docs/HOJA_DE_RUTA.pdf`).


---

## 📋 RESUMEN EJECUTIVO

Sistema de reacción-difusión (Gierer-Meinhardt) que valida computacionalmente la diversidad de himenóforos en hongos políporos mediante variación paramétrica.

| Campo | Valor |
|-------|-------|
| **Especies validadas** | 3 (Fomes, Brumalis, Squamosus) |
| **Rango de escala** | 8× (0.25-2.0 mm) |
| **Modelo** | Gierer-Meinhardt (1972) |
| **Versión** | 0.5.1 (PATCHED_v4_3) |
| **Licencia** | Dual: CC BY-NC 4.0 / Comercial |

---

## 🚀 INICIO RÁPIDO

### Instalación

```bash
# 1. Descomprimir
unzip FungalMorphoSpace_PATCHED_v4_3.zip
cd FungalMorphoSpace

# 2. Crear entorno virtual (recomendado)
python -m venv venv
source venv/bin/activate  # Linux/Mac
# venv\Scripts\activate   # Windows

# 3. Instalar dependencias
pip install -r requirements.txt

# 4. Instalar paquete
pip install -e .

# 5. Verificar instalación
python scripts/test_imports.py
```

### Primera Ejecución

```bash
# Prueba rápida (1 especie, 1 run)
python scripts/run_integrated_validation.py --species fomes --n_runs 1

# Validación completa (todas las especies)
python scripts/run_integrated_validation.py --species all --n_runs 5 --grid 1024
```

---

## 📁 ESTRUCTURA DEL PROYECTO

```
FungalMorphoSpace/
├── src/fungalmorphospace/          # Paquete Python principal
│   ├── __init__.py
│   ├── core/                       # Motor de simulación
│   │   ├── turing_simulator.py     # Solver Gierer-Meinhardt
│   │   └── kinetics.py             # Funciones cinéticas
│   ├── analysis/                   # Análisis de patrones
│   │   ├── topology_analyzer.py    # Métricas topológicas
│   │   └── visualization.py        # Visualización
│   └── utils/                      # Utilidades
│       └── sensitivity_analysis.py
├── scripts/                        # Scripts ejecutables
│   ├── run_integrated_validation.py    # Principal ⭐
│   ├── run_robustness_analysis.py      # Robustez n≥20
│   ├── run_parallel_validation.py      # Ejecución paralela
│   ├── plot_scaling_law.py             # Escalamiento
│   ├── test_laminillas.py              # Topología
│   └── test_imports.py                 # Verificación
├── data/
│   └── species_data.json           # Parámetros (FUENTE CANÓNICA)
├── docs/                           # Documentación
│   ├── METODOS_TESIS.md
│   ├── INFORME_TECNICO_TESIS.md
│   ├── GUIA_COMPLETA_PROYECTO.md   # Este archivo
│   └── GUIA_EXPERIMENTOS_BLINDAJE.md
├── config/
│   └── turing_params.yaml          # Config alternativa
├── results/                        # Outputs (generado)
├── LICENSE                         # Licencia dual
├── README.md
├── requirements.txt
└── pyproject.toml
```

---

## 🔬 ESPECIES Y PARÁMETROS

### Especies Validadas

| Especie | Key | D_v/D_u | ρ | Grid | Rol |
|---------|-----|---------|---|------|-----|
| **Fomes fomentarius** | `fomes` | 150 | 0.2 | 512 | Baseline |
| **Lentinus brumalis** | `brumalis` | 250 | 0.7 | 512 | Denso |
| **Polyporus squamosus** | `squamosus` | 3750 | 5.0 | 1024 | Disperso |

### Fuente de Parámetros

**Archivo canónico:** `data/species_data.json`

Todos los scripts leen parámetros de este archivo. Para modificar parámetros:
1. Editar `species_data.json`, o
2. Usar flags de override en línea de comandos

---

## 💻 COMANDOS DE EJECUCIÓN

### Script Principal: `run_integrated_validation.py`

```bash
# Sintaxis básica
python scripts/run_integrated_validation.py --species ESPECIE --n_runs N

# Opciones:
--species {fomes|brumalis|squamosus|all}  # Especie(s)
--n_runs N                                 # Número de repeticiones
--grid SIZE                                # Tamaño de grid (512, 1024, 2048)
--T_target TIME                            # Tiempo físico objetivo (default: 5.0)
--D_v_D_u RATIO                           # Override ratio difusión
--b VALUE                                  # Override parámetro b
--rho VALUE                                # Override parámetro ρ
```

### Ejemplos de Uso

```bash
# 1. Prueba rápida
python scripts/run_integrated_validation.py --species fomes --n_runs 1

# 2. Validación completa
python scripts/run_integrated_validation.py --species all --n_runs 5 --grid 1024

# 3. Alta resolución
python scripts/run_integrated_validation.py --species squamosus --n_runs 3 --grid 2048

# 4. Exploración paramétrica
python scripts/run_integrated_validation.py --species fomes --rho 0.5 --n_runs 3

# 5. Override múltiple
python scripts/run_integrated_validation.py --species brumalis --D_v_D_u 300 --rho 0.8
```

### Otros Scripts

```bash
# Análisis de robustez (n≥20 recomendado)
python scripts/run_robustness_analysis.py --species brumalis --n_runs 20

# Escalamiento alométrico
python scripts/plot_scaling_law.py

# Continuidad topológica
python scripts/test_laminillas.py

## 📊 OUTPUTS Y RESULTADOS

### Estructura de Outputs

```
results/
├── validation_summary.csv          # Tabla principal
├── validation_TIMESTAMP/           # Carpeta por ejecución
│   ├── fomes_run_0.png
│   ├── fomes_run_0_analysis.png
│   ├── brumalis_run_0.png
│   └── ...
├── robustness/                     # Si se ejecuta robustness
│   └── brumalis_n20/
└── scaling_validation_final.png    # Si se ejecuta plot_scaling
```

### Formato del CSV

`validation_summary.csv` contiene **1 fila por especie** con estadísticas agregadas:

| Columna | Descripción |
|---------|-------------|
| species | Nombre de especie |
| n_runs | Número de repeticiones |
| lambda_mean | Media de λ (px) |
| lambda_std | Desviación estándar |
| spots_mean | Media de componentes |
| density_mean | Densidad media (poros/mm) |
| euler_mean | Característica de Euler media |
| ... | Otras métricas |

---

## ⏱️ TIEMPOS ESTIMADOS

| Configuración | Tiempo aproximado |
|---------------|-------------------|
| `--species fomes --n_runs 1 --grid 512` | ~5 min |
| `--species all --n_runs 1 --grid 512` | ~20 min |
| `--species all --n_runs 5 --grid 512` | ~90 min |
| `--species all --n_runs 5 --grid 1024` | ~3 horas |
| `--species squamosus --n_runs 1 --grid 2048` | ~45 min |

**Nota:** Squamosus requiere más tiempo debido a mayor grid y dt más pequeño.

---

## 🔧 SOLUCIÓN DE PROBLEMAS

### Error: ImportError

```bash
# Verificar instalación
python scripts/test_imports.py

# Reinstalar paquete
pip install -e .
```

### Error: "species not found"

```bash
# Especies válidas:
--species fomes      # O
--species brumalis   # O
--species squamosus  # O
--species all        # Todas
```

### Memoria insuficiente

```bash
# Reducir grid
--grid 512  # En lugar de 1024/2048

# O ejecutar una especie a la vez
--species fomes
# Luego:
--species brumalis
```

### λ(FFT) QC FAIL

Esto es normal para patrones con periodicidad imperfecta. El sistema automáticamente usa λ(autocorr) como fallback.

### Muchos "spots" (>1000)

Verificar que estés usando v4.3+. Versiones anteriores tenían bug de sobreconteo.

---

## 📈 MÉTRICAS DISPONIBLES

### Métricas Primarias

| Métrica | Método | Unidad |
|---------|--------|--------|
| **λ (wavelength)** | FFT radial + autocorr | px o µm |
| **N_spots** | Componentes conexas | count |
| **Densidad** | N_spots / área | poros/mm |

### Métricas Secundarias

| Métrica | Método | Interpretación |
|---------|--------|----------------|
| **χ (Euler)** | Topología | >0: spots, <0: red |
| **CV_spacing** | Variación de distancias | <0.3: regular |
| **mean_spacing** | KD-tree centroides | px o µm |

### Control de Calidad

| Flag | Significado |
|------|-------------|
| `wavelength_fft_qc_pass` | True si FFT es confiable |
| `wavelength_autocorr_qc_pass` | True si autocorr es confiable |

---

## 🎓 FLUJO DE TRABAJO PARA TESIS

### Semana 1: Validación Base

```bash
# Ejecutar validación completa
python scripts/run_integrated_validation.py --species all --n_runs 5 --grid 1024

# Verificar resultados
cat results/validation_summary.csv
```

### Semana 2: Robustez

```bash
# Análisis de robustez para especie crítica
python scripts/run_robustness_analysis.py --species brumalis --n_runs 20

# Verificar CV < 5%
```

### Semana 3: Figuras

```bash
# Escalamiento alométrico
python scripts/plot_scaling_law.py

# Continuidad topológica
python scripts/test_laminillas.py
```

### Semana 4: Escritura

Usar documentos de `docs/`:
- `METODOS_TESIS.md` → Sección de Métodos
- `INFORME_TECNICO_TESIS.md` → Marco teórico y discusión

---

## 🔐 LICENCIAMIENTO

### Uso Académico (Gratuito)

```
Licencia: CC BY-NC 4.0
Requisito: Citar apropiadamente

Cita sugerida:
Ahumada Durán, M. (2026). FungalMorphoSpace: Validación computacional de 
morfologías de himenóforos mediante sistemas de reacción-difusión.
```

### Uso Comercial

Contactar al autor para licencia comercial.

---

## 📚 DOCUMENTACIÓN RELACIONADA

| Documento | Propósito |
|-----------|-----------|
| `README.md` | Inicio rápido y citación |
| `METODOS_TESIS.md` | Ecuaciones para copiar a tesis |
| `INFORME_TECNICO_TESIS.md` | Marco teórico completo |
| `GUIA_EXPERIMENTOS_BLINDAJE.md` | Protocolos de validación |
| `INFORME_CONTEXTO_PROBLEMAS_RESUELTOS.md` | Historial de correcciones |
| `LICENSE` | Términos de licencia |

---

## ❓ PREGUNTAS FRECUENTES

### ¿Por qué Brumalis tiene λ menor que Fomes si D_v/D_u es mayor?

Porque λ depende de D/ρ (difusión efectiva), no solo de D_v/D_u:
- Fomes: D_v/ρ = 150/0.2 = 750
- Brumalis: D_v/ρ = 250/0.4 = 625
- 625 < 750 → λ_Brumalis < λ_Fomes

### ¿Por qué el CSV tiene 3 filas y no 15?

El script agrega los runs por especie. Si ejecutas `--species all --n_runs 5`, obtienes:
- 5 runs de Fomes → 1 fila con media±std
- 5 runs de Brumalis → 1 fila con media±std
- 5 runs de Squamosus → 1 fila con media±std
- **Total: 3 filas**

### ¿Qué grid usar?

| Grid | Uso recomendado |
|------|-----------------|
| 512 | Pruebas rápidas, exploración |
| 1024 | Validación para tesis |
| 2048 | Publicación de alta calidad |

### ¿Cómo agrego una nueva especie?

1. Editar `data/species_data.json`
2. Agregar entrada con estructura similar
3. Asignar key corta en el mapeo del script

---

## 📞 SOPORTE

**Autor:** Mario Ahumada Durán  
**Email:** [tu-email]

Para reportar bugs o solicitar features, crear issue en el repositorio.

---

**© 2025-2026 Mario Ahumada Durán. Todos los derechos reservados.**

**Status:** ✅ Documentación completa y actualizada
