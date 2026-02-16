# Tests mínimos (locales)

Objetivo
--------
Tener un set de tests **rápidos** (segundos, no minutos) que permitan detectar:
- roturas de imports / empaquetado
- cambios en el contrato de output
- regresiones en métricas básicas (Topología / λ)

Este repositorio apunta a uso local. Por ahora no se configura CI.

## Dependencias

Instalar deps de desarrollo:

```bash
python3 -m pip install -r requirements.txt
python3 -m pip install -r requirements-dev.txt
```

## Ejecutar

Desde la raíz de la repo:

```bash
pytest -q
```

## Cobertura mínima definida

1. **Imports y presets**
   - importa `fungalmorphospace.runners` sin errores
   - `SPECIES_DATABASE` cargada desde `data/species_data.json`

2. **Topología / métricas (smoke)**
   - `TopologyAnalyzer.compute_all_metrics()` corre en un patrón sintético
   - retorna claves esperadas (wavelength_fft, wavelength_autocorr, QC flags)

3. **Contrato de output (docs)**
   - existe `docs/OUTPUT_CONTRACT.md`
   - contiene el árbol canónico y los paths esperados

4. **CLIs importables (smoke)**
   - `scripts/run_integrated_validation.py` y `scripts/run_parallel_validation.py` son ejecutables
   - `scripts/test_imports.py` no falla
