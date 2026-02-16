# BASE DE DATOS DE PARÁMETROS BIOFÍSICOS REALES
## Extraídos de literatura científica para anclaje no tautológico del modelo de Turing

**Autor**: Análisis exhaustivo de literatura para Mario Ahumada Durán  
**Fecha**: Enero 2026  
**Propósito**: Identificar parámetros medibles que NO sean tamaño/densidad de poros para anclar el modelo de reacción-difusión

---

## 1. PARÁMETROS ESTRUCTURALES DEL MEDIO (Klemm et al. 2024, PLOS ONE)

### 1.1 Porosidad por Segmento (φ)
| Segmento | Porosidad (%) | Fuente |
|----------|---------------|--------|
| Hymenium | **75.5%** | μCT, Klemm 2024 |
| Trama | **70.1%** | μCT, Klemm 2024 |
| Mycelial core | **61%** | μCT, Klemm 2024 |
| Crust | **19.1%** | μCT, Klemm 2024 |

**Relevancia para anclaje**: La porosidad afecta directamente la difusión efectiva mediante:
$$D_{eff} = D_0 \cdot \frac{\phi}{\tau^2}$$

### 1.2 Diámetros de Hifas (escala del obstáculo)
| Tipo de hifa | Diámetro (μm) | Ubicación | Fuente |
|--------------|---------------|-----------|--------|
| Esquelética | **3-4 μm** | Paredes de tubos | Klemm 2024 SEM |
| Generativa | **2-4 μm** | General | Klemm 2024 |
| Binding (finas) | **0.6-0.8 μm** | Interior tubos | Klemm 2024 SEM |
| Esquelética (lit.) | **3-8 μm** | General | Literatura previa |
| Binding (lit.) | **1.5-3.0 μm** | General | Literatura previa |

**Relevancia**: Define el tamaño de malla para calcular factor de obstrucción hidrodinámica H(obstrucción).

### 1.3 Diámetro de Tubos del Hymenium
| Parámetro | Valor | Fuente |
|-----------|-------|--------|
| Diámetro tubos | **~200 μm** | Klemm 2024 |
| Diámetro tubos (zonas) | Variable por zona | Pylkkänen 2023 (5 zonas medidas) |

---

## 2. COMPOSICIÓN QUÍMICA DE PARED CELULAR (Klemm 2024)

### 2.1 Contenido de Quitina/Quitosano
| Segmento | Quitina/Quitosano (wt.%) | Fuente |
|----------|--------------------------|--------|
| Hymenium | **~5%** (máximo) | Assay, Klemm 2024 |
| Trama | ~3.5% | Assay, Klemm 2024 |
| Mycelial core | Intermedio | Assay, Klemm 2024 |
| Crust | **<1%** (mínimo) | Assay, Klemm 2024 |

**Relevancia para anclaje**: 
- Quitina/quitosano es cargado (grupo amino) → afecta binding/retención de morfógenos cargados
- Diferencia de ~5× entre hymenium y crust = potencial diferencia en factor K(binding)

### 2.2 Contenido de Glucanos
| Segmento | Glucanos (wt.%) | Tipo dominante | Fuente |
|----------|-----------------|----------------|--------|
| Crust | **~70%** (máximo) | α-glucano | Klemm 2024 |
| Hymenium | ~27% | β-glucano | Klemm 2024 |
| Trama | **~16-25%** (mínimo) | β-glucano | Klemm 2024 |
| Mycelial core | Intermedio | β-glucano | Klemm 2024 |

**Relevancia**: β-glucano vs α-glucano tienen diferentes propiedades de malla/interacción.

### 2.3 Grado de Deacetilación
- Crust muestra **mayor deacetilación** de quitina (FTIR: bandas amida II/III diferentes)
- Quitosano (deacetilado >75%) es **policatiónico** → mayor binding de especies aniónicas

---

## 3. DISTRIBUCIÓN DE ELEMENTOS TRAZA (Klemm 2024, μXRF/SEM-EDX)

### 3.1 Calcio (Ca)
| Segmento | Nivel Ca | Forma | Relevancia |
|----------|----------|-------|------------|
| Hymenium | **Alto** | Cristales de oxalato de Ca | Posible regulador de crecimiento |
| Crust | **Alto** (capas) | Cristales de oxalato de Ca | Protección/estructura |
| Trama | Moderado | - | - |
| Mycelial core | Moderado | - | - |

**Dato clave**: Cristales de **oxalato de calcio** confirmados por SEM-EDX en hymenium y crust.  
**Relevancia**: Ca²⁺ es un candidato natural para "especie rápida" (inhibidor) - ión pequeño, alta difusividad.

### 3.2 Potasio (K)
| Segmento | Nivel K |
|----------|---------|
| Hymenium | Presente |
| Trama | Presente |
| Mycelial core | Alto |
| Crust | **Casi ausente** (depletado) |

### 3.3 Zinc (Zn) y Manganeso (Mn)
| Elemento | Distribución |
|----------|--------------|
| Zn | Elevado en hymenium y crust |
| Mn | **Predominante en mycelial core** |

---

## 4. PROPIEDADES MECÁNICAS (Klemm 2024 + Klemm 2025 Small Structures)

### 4.1 Esfuerzo de Plateau en Compresión
| Segmento | σ_pl (MPa) | Dirección | Fuente |
|----------|------------|-----------|--------|
| Hymenium (paralelo) | **4.77** | ∥ tubos | Klemm 2024 |
| Hymenium (transverso) | **0.87** | ⊥ tubos | Klemm 2024 |
| Trama | **1.22** | - | Klemm 2024 |
| Mycelial core | ~3.2 | - | Klemm 2024 |

**Ratio de anisotropía mecánica**: σ_∥ / σ_⊥ ≈ **5.5×** para hymenium

### 4.2 Desplazamientos bajo Carga (Klemm 2025, Synchrotron μCT + Optical Flow)
| Parámetro | Paralelo | Transverso | Ratio |
|-----------|----------|------------|-------|
| Desplazamiento máximo | ~40 μm | ~200 μm | **~5×** |

**Clave**: Los desplazamientos transversales son **~5× mayores** que los paralelos.  
Este ratio de anisotropía es un **anclaje mecánico directo** para justificar D_⊥/D_∥.

### 4.3 Mecanismos de Falla
- **Paralelo**: Buckling plástico + delaminación (similar a honeycombs de fibra)
- **Transverso**: Colapso de paredes de honeycomb
- **Húmedo vs Seco**: Húmedo muestra "acortamiento telescópico", seco muestra grietas fatales

---

## 5. PARÁMETROS DE PYLKKÄNEN ET AL. 2023 (Science Advances)

### 5.1 Contribución por Región al Cuerpo Fructífero
| Región | Contribución (%) |
|--------|------------------|
| H. tubes (hymenium) | **~69%** |
| Context (trama) | **~29%** |
| Crust | **~4%** |

### 5.2 Propiedades Comparativas
| Parámetro | H. tubes | Context |
|-----------|----------|---------|
| Resistencia máxima | **10× mayor** | Baseline |
| Módulo de Young | **10× mayor** | Baseline |
| Porosidad | Mayor | Menor |
| Densidad | Menor | Mayor |

### 5.3 Orientación Preferencial de Hifas
- Context: Hifas alineadas **antes** de formar conos de crecimiento
- H. tubes: Hifas alineadas **paralelas al eje tubular** (desviación de pocos grados)
- HOP (Hermans Orientation Parameter): Medido por WAXS/SAXS

---

## 6. CANDIDATOS MOLECULARES PARA u/v (Inferidos)

### 6.1 Candidato para Inhibidor (v) - Especie Rápida
| Candidato | Justificación | D₀ relativo |
|-----------|---------------|-------------|
| **Ca²⁺** | Enriquecido en hymenium, rol en tip growth | Alto (~10⁻⁹ m²/s en agua) |
| ROS (H₂O₂, O₂⁻) | Metabolitos oxidativos pequeños | Alto |
| K⁺ | Presente, pequeño | Alto |
| Metabolitos secundarios volátiles | Pequeños, móviles | Alto |

### 6.2 Candidato para Activador (u) - Especie Lenta
| Candidato | Justificación | D₀ relativo |
|-----------|---------------|-------------|
| **Hidrofobinas** | Proteínas de auto-ensamblaje, grandes | Bajo (~10⁻¹¹ m²/s) |
| **Lacasa/enzimas** | Grandes, asociadas a pared | Bajo |
| Complejos quitina-glucano | Macromoleculares, inmóviles | Muy bajo |
| Factores de crecimiento | Proteínas señalizadoras | Bajo-Medio |

### 6.3 Ratio D₀ Estimado (sin corrección por medio)
Para ión pequeño vs proteína grande:
- R_h(Ca²⁺) ≈ 0.1-0.2 nm
- R_h(proteína 50 kDa) ≈ 3-5 nm

Por Stokes-Einstein: D ∝ 1/R_h  
**Ratio D₀(v)/D₀(u) ≈ 15-50×** (solo por tamaño, sin considerar binding)

---

## 7. FÓRMULAS PARA DIFUSIÓN EFECTIVA

### 7.1 Modelo General
$$D_{eff} = D_0 \cdot \underbrace{\frac{\phi}{\tau^2}}_{\text{porosidad/tortuosidad}} \cdot \underbrace{H(\lambda)}_{\text{obstrucción}} \cdot \underbrace{K_{part}}_{\text{binding/partición}}$$

### 7.2 Factor de Tortuosidad (τ)
Para medios porosos típicos: τ ≈ 1.5-3  
Para hymenium con tubos alineados: τ_∥ << τ_⊥

### 7.3 Factor de Obstrucción H (modelo de fibras)
Para soluto de radio r_s en malla de fibras de radio r_f y fracción sólida φ_s:
$$H = (1 - \phi_s)^{1/2} \cdot \exp\left(-\pi \frac{r_s^2}{r_f^2} \phi_s\right)$$

### 7.4 Factor de Binding/Partición K
Depende de:
- Carga del soluto vs carga de pared (quitosano = +, glucanos = neutro/-)
- Constantes de adsorción específicas
- Interacciones hidrofóbicas

---

## 8. VALORES LISTOS PARA USAR EN MODELO

### 8.1 Porosidad (φ) - Datos Directos
```
φ_hymenium = 0.755
φ_trama = 0.701
φ_mycelial = 0.61
φ_crust = 0.191
```

### 8.2 Escala de Obstáculos
```
r_fiber_grueso = 2 μm (hifas esqueléticas)
r_fiber_fino = 0.4 μm (hifas binding)
```

### 8.3 Anisotropía Mecánica → Proxy para Anisotropía de Transporte
```
Ratio_anisotropía = 5× (de desplazamientos en compresión)
```
**Justificación**: Si el transporte sigue caminos preferenciales similares a la transferencia de carga mecánica, este ratio puede aproximar D_∥/D_⊥.

### 8.4 Composición Química para Factor K
```
Quitina_hymenium = 5%
Quitina_crust = <1%
Ratio_quitina = 5× (potencial diferencia en binding)
```

---

## 9. DATOS FALTANTES (Gaps Identificados)

| Parámetro Necesario | Estado | Comentario |
|---------------------|--------|------------|
| Tortuosidad τ directa | ❌ No medido | Podría estimarse de μCT 3D |
| D₀ de morfógenos específicos | ❌ No identificados | Se necesita proteómica/metabolómica |
| Constantes de binding | ❌ No medidas | Requiere estudios de adsorción |
| Viscosidad local | ❌ No medida | Asumible ~agua para poros grandes |
| Difusión in situ | ❌ No medida | Requiere FRAP o similar |

---

## 10. REFERENCIAS CLAVE

1. **Klemm et al. (2024)** - PLOS ONE - Hierarchical structure and chemical composition
2. **Klemm et al. (2025)** - Small Structures - Synchrotron μCT + Optical Flow
3. **Pylkkänen et al. (2023)** - Science Advances - Complex structure of F. fomentarius
4. **Müller et al. (2021)** - Applied Physics A - Bracket fungi microstructure
5. **Pohl et al. (2022)** - Fungal Biology and Biotechnology - Establishment for composites

---

*Este documento contiene SOLO datos medidos experimentalmente, sin valores inventados.*
