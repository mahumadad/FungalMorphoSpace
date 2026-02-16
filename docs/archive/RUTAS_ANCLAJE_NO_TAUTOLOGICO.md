# RUTAS DE ANCLAJE NO TAUTOLÓGICO
## Estrategias para conectar el modelo de Turing con datos biológicos reales

**Documento de trabajo para tesis - Mario Ahumada Durán**  
**Enero 2026**

---

## PRINCIPIO FUNDAMENTAL

> **Anclaje no tautológico** = El parámetro del modelo debe derivarse de una medición independiente del patrón que se quiere predecir.

```
❌ TAUTOLOGÍA: Observo poro de 2mm → Ajusto D hasta que sale 2mm → "Predigo" 2mm
✅ ANCLAJE:    Mido porosidad/química → Calculo D → Modelo predice 2mm → Comparo con real
```

---

## RUTA A: DIFUSIÓN EFECTIVA EN MEDIO POROSO FIBROSO
### (Recomendada como anclaje principal)

### Fundamento
Klemm et al. miden las propiedades del **medio** (porosidad, composición, estructura de hifas), no del tamaño de poro. Usamos estas propiedades para **restringir** los coeficientes de difusión efectiva.

### Cadena Causal
```
┌─────────────────────────────────────────────────────────────────┐
│  MEDICIONES (Klemm 2024)                                        │
│  ├── Porosidad φ = 75.5% (hymenium)                            │
│  ├── Diámetro hifas: 3-4 μm (esqueléticas), 0.6-0.8 μm (finas) │
│  ├── Composición: quitina 5%, glucanos 27%                      │
│  └── Ca²⁺ enriquecido en hymenium                              │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  MODELO DE TRANSPORTE                                           │
│  D_eff = D₀ × (φ/τ²) × H(obstrucción) × K(binding)             │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  SELECCIÓN DE AGENTES u/v                                       │
│  v = Ca²⁺ (pequeño, rápido, presente)                          │
│  u = Hidrofobina o similar (grande, lento, estructural)         │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  CÁLCULO DE RATIO                                               │
│  D_v/D_u = [D₀(v)/D₀(u)] × [H_v/H_u] × [K_v/K_u]               │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  PREDICCIÓN DEL MODELO                                          │
│  λ_teórico = f(D_v/D_u, otros parámetros)                      │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  VALIDACIÓN EXTERNA                                             │
│  Comparar λ_teórico con λ_observado (NO USADO EN CALIBRACIÓN)   │
└─────────────────────────────────────────────────────────────────┘
```

### Implementación Cuantitativa

**Paso 1: Ratio D₀ por tamaño molecular**
```
Para Ca²⁺ (R_h ≈ 0.1 nm) vs Hidrofobina (~15 kDa, R_h ≈ 2 nm):
D₀(Ca²⁺) ≈ 7×10⁻¹⁰ m²/s (en agua)
D₀(Hidrofobina) ≈ 7×10⁻¹¹ m²/s (estimado)
Ratio D₀ ≈ 10×
```

**Paso 2: Factor de obstrucción H**
```
Para malla de hifas (r_f = 2 μm, φ_sólido = 0.245):
H(Ca²⁺) ≈ 0.95 (casi sin obstrucción)
H(Hidrofobina) ≈ 0.7 (obstrucción moderada)
Ratio H ≈ 1.4×
```

**Paso 3: Factor de binding K**
```
Quitosano (cargado +) retiene aniones, repele cationes
Ca²⁺ es catión → K(Ca²⁺) ≈ 0.8-1.0 (poca retención)
Proteína puede tener afinidad por pared → K(prot) ≈ 0.1-0.5
Ratio K ≈ 2-8×
```

**Resultado combinado:**
```
D_v/D_u ≈ 10 × 1.4 × 4 = 56× (orden de magnitud: 10-100)
```

### Fortalezas
- ✅ Usa datos medidos de Klemm (φ, diámetros hifas, composición)
- ✅ No mira el tamaño de poro como input
- ✅ Físicamente justificable
- ✅ Candidatos moleculares razonables (Ca²⁺ documentado en hymenium)

### Debilidades
- ⚠️ No se han medido D₀ de morfógenos específicos in situ
- ⚠️ Requiere asumir identidad de u y v
- ⚠️ Factor K es el más incierto

---

## RUTA B: ANCLAJE MECÁNICO (POROELASTICIDAD)
### (Alternativa para modelo mecánico de Kuhar)

### Fundamento
Si interpretamos u/v como campos mecánicos (estrés/relajación), los coeficientes "D" representan difusividad de tensiones o tiempos de relajación viscoelástica.

### Cadena Causal
```
┌─────────────────────────────────────────────────────────────────┐
│  MEDICIONES MECÁNICAS (Klemm 2024 + 2025)                       │
│  ├── σ_∥ = 4.77 MPa, σ_⊥ = 0.87 MPa                            │
│  ├── Desplazamiento: ∥ ~40 μm, ⊥ ~200 μm (ratio 5×)            │
│  └── Mecanismos: buckling, delaminación                         │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  MODELO POROELÁSTICO/VISCOELÁSTICO                              │
│  "D" = difusividad de relajación = E/(η·τ_relax)               │
│  Anisotropía de D refleja anisotropía mecánica                  │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  INTERPRETACIÓN DE u/v                                          │
│  u = acumulación local de estrés (buckling iniciador)           │
│  v = relajación/redistribución de carga                         │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  RATIO DESDE ANISOTROPÍA                                        │
│  D_⊥/D_∥ ≈ 5× (del ratio de desplazamientos)                   │
│  D_v/D_u derivable de tiempos de relajación                     │
└─────────────────────────────────────────────────────────────────┘
```

### Implementación

**Datos disponibles de Klemm 2025:**
```
Ratio anisotropía de desplazamiento = 200/40 = 5×
Ratio anisotropía de esfuerzo = 4.77/0.87 = 5.5×
```

**Modelo de difusividad poroelástica:**
```
D_poro = k/(μ·S)
donde:
  k = permeabilidad (relacionada con φ y estructura)
  μ = viscosidad del fluido intersticial
  S = coeficiente de almacenamiento
```

**Para anisotropía:**
```
D_∥/D_⊥ ≈ k_∥/k_⊥ ≈ 5× (de anisotropía estructural)
```

### Fortalezas
- ✅ Datos mecánicos directos de compresión in situ
- ✅ Compatible con visión de Kuhar (modelo mecánico)
- ✅ Ratio 5× está medido experimentalmente

### Debilidades
- ⚠️ Requiere derivación explícita de parámetros de Turing desde mecánica
- ⚠️ Mapeo tensión→morfógeno es abstracto
- ⚠️ No distingue claramente activador de inhibidor

---

## RUTA C: ANCLAJE POR ESCALA JERÁRQUICA
### (Comparación multi-especie)

### Fundamento
Usar la variación entre especies (Fomes, Brumalis, Squamosus) como "experimento natural" donde las diferencias en estructura/composición predicen diferencias en tamaño de poro.

### Cadena Causal
```
┌─────────────────────────────────────────────────────────────────┐
│  DATOS POR ESPECIE (Literatura + Tu base de datos)              │
│  Fomes: poro ~400 μm, hymenium denso, quitina alta              │
│  Brumalis: poro ~250 μm, estructura intermedia                  │
│  Squamosus: poro ~2000 μm, estructura laxa                      │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  CORRELACIÓN ESTRUCTURA-PATRÓN                                  │
│  Si φ, composición o grosor de hifas varían entre especies,     │
│  ¿predicen el ratio de tamaños de poro?                         │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  TEST: Ley de escala                                            │
│  λ ∝ √(D_v/D_u) y D_eff ∝ f(φ, estructura)                     │
│  ¿La variación en φ entre especies explica λ_Fomes/λ_Squamosus? │
└─────────────────────────────────────────────────────────────────┘
```

### Problema
- ⚠️ Solo tienes datos detallados de Fomes (Klemm)
- ⚠️ Brumalis y Squamosus no tienen caracterización equivalente
- ⚠️ Sería correlacional, no predictivo

### Potencial
Si pudieras obtener/estimar φ para las tres especies, podrías hacer:
```
λ ∝ √(D_v/D_u) ∝ √(φ/τ²)

Si φ_Fomes = 0.755 y φ_Squamosus = 0.90 (hipotético):
λ_Squamosus/λ_Fomes ∝ √(0.90/0.755) ≈ 1.09 (no explica 5× diferencia)
```
→ La porosidad sola no explica la variación; se necesita otro factor.

---

## RUTA D: MODELO MULTIESCALA
### (Integración de escalas micro-meso-macro)

### Fundamento
Conectar explícitamente las escalas:
1. **Nano**: Composición de pared celular (quitina, glucanos)
2. **Micro**: Red de hifas (porosidad, orientación)
3. **Meso**: Tubos del hymenium (estructura honeycomb)
4. **Macro**: Patrón de poros

### Implementación Propuesta

```
┌─────────────────────────────────────────────────────────────────┐
│  ESCALA NANO (composición)                                      │
│  Input: quitina%, glucanos%, Ca                                 │
│  Output: K_binding, carga superficial                           │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  ESCALA MICRO (red de hifas)                                    │
│  Input: diámetro hifas, densidad, orientación                   │
│  Output: τ (tortuosidad), H (obstrucción)                       │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  ESCALA MESO (tubos honeycomb)                                  │
│  Input: φ_hymenium, anisotropía mecánica                        │
│  Output: D_eff anisótropo                                       │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  ESCALA MACRO (patrón)                                          │
│  Input: D_u_eff, D_v_eff, cinética                              │
│  Output: λ (longitud de onda del patrón)                        │
└─────────────────────────────────────────────────────────────────┘
```

### Ecuaciones de Acoplamiento

**Micro → Meso:**
```
D_eff(micro) = D₀ × (φ_micro/τ²_micro) × H_micro × K_nano
```

**Meso → Macro:**
```
D_eff(meso) = D_eff(micro) × f(φ_meso, anisotropía_meso)
```

**Patrón:**
```
λ = 2π × √(D_v_eff/ρ) × g(D_v_eff/D_u_eff)
```

### Fortalezas
- ✅ Integra todos los datos disponibles
- ✅ Físicamente más realista
- ✅ Permite identificar cuál escala domina

### Debilidades
- ⚠️ Complejidad alta
- ⚠️ Más parámetros = más incertidumbre
- ⚠️ Requiere implementación numérica cuidadosa

---

## RUTA E: VALIDACIÓN CRUZADA CON MORFOLOGÍA SECUNDARIA
### (Usar múltiples outputs como test)

### Fundamento
Si el modelo está bien anclado, debería predecir no solo el tamaño de poro sino también:
- Regularidad del patrón (hexagonalidad)
- Gradiente de tamaño entre zonas
- Respuesta a perturbaciones

### Implementación
```
PREDICCIONES MÚLTIPLES:
1. λ (tamaño de poro) → Comparar con Klemm
2. Regularidad (CV de tamaños) → Comparar con imágenes
3. Orientación → Comparar con HOP de WAXS
4. Respuesta a compresión → Comparar con optical flow
```

### Ventaja
Si múltiples predicciones coinciden, el anclaje es más robusto que si solo coincide una.

---

## RECOMENDACIÓN FINAL

### Estrategia Óptima: RUTA A + E

1. **Usar Ruta A** como anclaje principal:
   - Derivar D_v/D_u desde φ, estructura de hifas, y candidatos moleculares
   - NO usar tamaño de poro como input

2. **Usar Ruta E** para validación múltiple:
   - Comparar λ predicho vs observado
   - Verificar regularidad del patrón
   - Chequear consistencia con anisotropía

3. **Mencionar Ruta B** como interpretación alternativa:
   - "El modelo es matemáticamente isomorfo a poroelasticidad"
   - Los datos mecánicos de Klemm 2025 son consistentes

### Frase Clave para Tesis/Paper

> "Los coeficientes de difusión efectiva se derivaron de las propiedades medidas del medio poroso (porosidad φ = 75.5%, estructura de red hifal con diámetros de 0.6-4 μm, y composición química rica en quitina) siguiendo modelos establecidos de transporte en medios fibrosos, **sin utilizar la morfología del patrón como input**. La longitud de onda predicha se comparó entonces con las observaciones como **validación externa** del anclaje biofísico."

---

## TABLA RESUMEN DE RUTAS

| Ruta | Anclaje desde | Dato clave | Fortaleza | Riesgo |
|------|---------------|------------|-----------|--------|
| **A** | Medio poroso | φ, hifas, composición | Físicamente sólido | K incierto |
| **B** | Mecánica | σ, desplazamientos | Datos directos | Mapeo abstracto |
| **C** | Multi-especie | Comparación | Natural | Pocos datos |
| **D** | Multiescala | Integración | Completo | Complejo |
| **E** | Validación múltiple | Varios outputs | Robusto | Complementario |

**Recomendación**: A + E (con mención de B)
