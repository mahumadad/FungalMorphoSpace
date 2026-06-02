# Estrategia de tesis — FungalMorphoSpace y derivados

**Mario Ahumada Durán**
**Junio 2026**

> Documento de decisión estratégica. Fija la hipótesis central, la secuencia de
> papers, los criterios de falsabilidad y los hallazgos de auditoría de código que
> condicionan el plan. No es un roadmap de features; es la respuesta a la pregunta
> "¿qué publico, en qué orden, y qué tengo que arreglar antes de afirmarlo?".

---

## 0. TL;DR

1. El **riesgo dominante es narrativo**, no matemático: hay 4 papers latentes mezclados.
   Hay que comprometerse con UNA hipótesis central por artículo.
2. La **tautología es triple** (`D_v/D_u`, `rho`, `b` ajustados por especie) y **no hay
   ningún dato biofísico independiente** para ninguna especie. Por lo tanto la "predicción
   leave-one-out biofísica" **está bloqueada a nivel de datos**, no de código.
3. El **código de simulación y medición es de buena calidad**; el problema es de **diseño
   experimental** (de dónde salen los parámetros), no de implementación.
4. Secuencia recomendada: **(P1) tool paper honesto YA** → **(P2) reformulación de nodo
   no-difusible** (mata la tautología sin datos nuevos) → **(P3) MCAT mecano-químico** →
   parkear TDA y biomimética/TPMS como líneas separadas.

---

## 1. El problema narrativo

Del material acumulado emergen cuatro artículos potenciales, cada uno con hipótesis distinta:

| # | Paper latente | Hipótesis central | Estado |
|---|---------------|-------------------|--------|
| A | FungalMorphoSpace como herramienta | "Un único mecanismo RD genera el morfo-espacio de himenóforos" | Código listo, paper honesto escrito |
| B | Modelo MCAT mecano-químico | "Turing fija periodicidad; la mecánica fija la escala absoluta" | Solo diseño conceptual |
| C | TDA de redes miceliales | "Las especies tienen firmas topológicas (H1) distintas" | Datos + resultados preliminares (proyecto aparte) |
| D | Biomimética / TPMS | "El patrón de poros mapea a superficies útiles para diseño" | Solo idea |

**Regla de oro:** un revisor que no pueda responder *"¿cuál es la hipótesis falsable de este
paper?"* en una frase, rechaza. Mezclar A+B+C+D garantiza esa pregunta. **Un paper, una hipótesis.**

C y D son artículos independientes con su propio dataset y su propia pregunta. Se parkean.
No se mencionan como "contribución" de A/B salvo una línea de "future work".

---

## 2. El bloqueo de datos (lo que condiciona todo)

Verificado en `data/species_data.json`: las 4 especies (`fomes`, `brumalis`, `squamosus`,
`trametes`) tienen `D_v_D_u`, `rho` y `b` puestos a mano. **Ningún campo de porosidad,
tortuosidad, %quitina, %glucano ni radio hifal.** El `source: "Klemm 2024"` de *fomes* es
una justificación post-hoc, no un parámetro derivado.

**Consecuencias:**

- La tautología no es de un parámetro (`D_v/D_u`) sino de **tres** (`D_v/D_u`, `rho`, `b`),
  los tres fiteados por especie a su morfología observada. El tuning-log de *brumalis*
  (`docs/archive/CALIBRACION_BRUMALIS_TUNING_JOURNEY.md`) lo documenta explícitamente:
  búsqueda manual hasta la "Goldilocks zone".
- **`BiophysicalAnchor` (proyecto #9) NO se debe construir todavía.** Derivar `D_v/D_u`
  desde φ/τ/quitina es inútil si no hay φ/τ/quitina medidos para ≥2 especies contra las
  cuales validar. Construirlo ahora produce un anclaje sin testear → otra forma de tautología.
- La "predicción leave-one-out entre especies" requiere **adquisición de datos**
  (wet-lab o arqueología bibliográfica de microestructura), no programación. Es un
  sub-proyecto con su propio cronograma e incertidumbre.

---

## 3. Hallazgos de auditoría de código (junio 2026)

Auditoría de `turing_simulator.py`, `kinetics.py`, `topology_analyzer.py`,
`integrated_validation.py`. **Separar dos planos:**

### 3.1 Calidad de código: BUENA

- **Solver:** Euler explícito con chequeo CFL correcto (`dt_max = dx²/(4·D_max)`,
  factor de seguridad 0.8). Laplaciano 5-puntos con BC periódicas vía `np.roll`,
  sin bugs de indexado.
- **Medición de λ:** FFT radial + perfil de potencia suavizado + detección de pico con
  QC real (prominencia `peak_ratio ≥ 3.0`, rechazo de λ > 0.5·dominio). Cascada
  FFT→autocorr→fallback bien diseñada. **Esta capa es sólida y reproducible.**
- **GM kinetics:** protección de división `u²/(v+ε)` presente.

→ **No hay que reescribir el simulador ni el medidor para publicar.** El problema no es ahí.

### 3.2 Validez científica: CIRCULAR (severidad CRÍTICA)

| Hallazgo | Severidad | Evidencia |
|----------|-----------|-----------|
| Validación contra el mismo dato usado para tunear | CRÍTICA | `visualization.py` compara `density_predicted` contra `density_observed`, ambos del mismo JSON |
| Parámetros fiteados por matching de λ | CRÍTICA | `CALIBRACION_BRUMALIS_TUNING_JOURNEY.md` documenta la búsqueda manual |
| Sin held-out / leave-one-out | CRÍTICA | Cada especie se evalúa con sus propios parámetros fiteados |
| Tolerancia de "match" hardcodeada (±1.5 poros/mm) | ALTA | `visualization.py` ~L445, sin justificación publicada |
| Calibración 8.7 µm/px derivada de 1 especie, aplicada a todas | ALTA | `species_data.json` bloque `calibration` |
| `qc_pass` se calcula pero NO gatea la conclusión | MODERADA | `integrated_validation.py` L621-628: `density_predicted` se computa con `wavelength_px` ignorando `qc_pass` |
| *squamosus* (`D_v/D_u=3750`) probablemente sub-resuelto en grilla | MODERADA | λ teórica ≈ 385 px > 0.5·512; QC FFT falla y cae a fallback |
| dt no re-evaluado durante el run | MODERADA | `self.dt` fijo tras `__init__`; reacciones stiff no re-chequean CFL |
| Sin test numérico de regresión | MODERADA | 5 tests, todos de plumbing/schema; ninguno fija λ contra caso conocido |

**Lectura:** el código *funciona y mide bien*; lo que es circular es el **experimento**.
Eso es una buena noticia: arreglar la ciencia no exige reescribir el motor.

---

## 4. Hipótesis central por paper (falsables)

### P1 — Tool paper (SoftwareX). Hipótesis NO predictiva.

> "Un sistema RD de Gierer-Meinhardt con variación paramétrica reproduce
> cualitativamente el morfo-espacio de himenóforos poliporo (poros→laberinto), y
> proveemos una herramienta reproducible para explorarlo."

- **Falsabilidad:** que exista una morfología de himenóforo en el rango que el modelo
  NO pueda producir con ningún parámetro razonable.
- **Lo que NO debe afirmar:** que predice tamaño absoluto de poro. El paper actual ya es
  honesto en esto (L305 del .tex). Mantenerlo así.
- **Riesgo:** bajo. Es un software paper. La tautología está disclosed, no oculta.

### P2 — Reformulación de nodo no-difusible. Hipótesis teórica.

> "El morfo-espacio de himenóforos puede generarse con una red de Turing de 3 nodos con
> difusividades biológicamente plausibles (D_v≈D_u) más un nodo inmóvil (densidad hifal
> local), eliminando la necesidad del ratio D_v/D_u extremo (hasta 3750)."

- **Base:** Marcon et al., eLife 2016 (Turing con señales de igual difusión + nodo fijo);
  Raspopovic/Sheth Science 2012/2014 (Bmp-Sox9-Wnt, Sox9 no-difusible).
- **Falsabilidad:** que NINGUNA topología de 3 nodos con D plausibles reproduzca el rango
  de λ observado → la reformulación falla y el ratio extremo es irreducible.
- **Por qué es la jugada fuerte:** mata la crítica "ese D_v/D_u es un parámetro inventado"
  **sin requerir datos biofísicos nuevos**. Es puro modelado computacional sobre el motor
  existente.
- **Riesgo:** medio. Hay que demostrar que la reformulación cubre el mismo morfo-espacio.

### P3 — MCAT mecano-químico. Hipótesis de seguimiento.

> "Turing fija la periodicidad del patrón; el acoplamiento mecánico (turgor + tensión de
> pared, función de quitina/glucano) fija la escala absoluta de poro."

- **Base:** Mercker & Brinkmann 2016 ("biomechanical forces may replace the long-range
  inhibitor"); framework poroelástico bifásico; mecánica de punta hifal (Bartnicki-García,
  Cell wall dynamics 2023); datos mecánicos de micelio (Islam 2017).
- **Falsabilidad:** que la escala de poro NO correlacione con parámetros mecánicos medibles
  (σ_wall, P_turgor) entre especies.
- **Requisito:** parámetros mecánicos (parcialmente disponibles). Más trabajo que P2.
- **Riesgo:** alto, pero alineado con la dirección actual del campo (active matter,
  mechanochemical patterning).

---

## 5. Secuencia recomendada

```
AHORA           P1 — Tool paper SoftwareX (honesto, NO predictivo)
  │              ├─ arreglar huecos de bib (§6)
  │              ├─ verificar que ninguna afirmación reclame predicción de escala
  │              └─ (opcional) gatear density_predicted sobre qc_pass
  │
SIGUIENTE       P2 — Reformulación nodo no-difusible
  │              ├─ implementar red 3-nodos sobre el motor existente
  │              ├─ barrer D plausibles + nodo hifal inmóvil
  │              └─ test: ¿reproduce las 3 especies sin D_v/D_u extremo?
  │
FOLLOW-UP       P3 — MCAT mecano-químico
  │              └─ requiere σ_wall; usa P2 como base no-tautológica
  │
PARK            C (TDA) y D (biomimética/TPMS) — papers independientes,
                 dataset y pregunta propios. No mezclar con A/B/P.
```

**Diferencia clave con el plan intuitivo previo:** no invertir en `BiophysicalAnchor`
todavía (data-blocked). Invertir en P2 (reformulación computacional), que da el mismo
golpe no-tautológico sin el bloqueo de datos.

---

## 6. Huecos bibliográficos a cerrar antes de enviar P1

`paper/references.bib` tiene solo **15 entradas**. Faltan (mínimo):

- **Kuhar et al. 2022, `10.1007/s12064-022-00363-z`** — el segundo paper de Kuhar (simula
  el morfo-espacio agaricomycete completo). El .bib solo tiene `-00380-8`. Para un paper
  que dice "validar la hipótesis de Kuhar" esto es obligatorio.
- **Davidson, *Mathematical modelling of fungal growth and function*, IMA Fungus 2011** —
  review estándar de modelado fúngico. Su ausencia es bandera roja para un revisor de
  micología matemática.
- Para P2/P3 (cuando toque): **Marcon et al. eLife 2016**, **Raspopovic/Sheth Science
  2012/2014**, **Mercker & Brinkmann Biol. Direct 2016**, **Islam et al. Sci. Rep. 2017**.

---

## 7. Acciones concretas inmediatas (P1)

- [ ] Añadir las 2 referencias faltantes al `.bib` y citarlas en Intro/Discusión.
- [ ] Releer `paper/main.tex` y confirmar que NINGUNA frase reclama predicción de escala
      absoluta. El framing actual ("hypothesis-generation tool") es correcto — protegerlo.
- [ ] (Opcional, defensivo) Gatear `density_predicted` sobre `qc_pass` en
      `integrated_validation.py` para que *squamosus* (QC-fail por sub-resolución) no
      reporte una densidad "validada" silenciosamente. O documentar la limitación.
- [ ] (Opcional) Añadir 1 test numérico de regresión: fijar parámetros, afirmar λ medida
      contra un valor conocido ±tolerancia, para blindar reproducibilidad.

---

## 8. Criterios de "listo para enviar" (P1)

1. Hipótesis central enunciable en 1 frase, no predictiva. ✅ (ya lo es)
2. Tautología disclosed explícitamente como limitación. ✅ (L305, L418-420)
3. Bib sin huecos evidentes para un revisor de micología matemática. ❌ (cerrar §6)
4. Ninguna afirmación de predicción de escala. ⚠️ (verificar §7)
5. Reproducibilidad: smoke test pasa + (deseable) 1 test numérico. ⚠️

---

## 9. Verificación empírica (2026-06-02)

Se instaló el entorno (Python 3.14, numpy 2.4) y se corrió el pipeline real. Hallazgos:

### 9.1 El programa corre; 2 de 3 especies son genuinas

| Especie | grid | λ medida | esperada | spots | patrón |
|---------|------|----------|----------|-------|--------|
| *Fomes* | 512 | 43.5 px | 46 | 164 | genuino ✓ |
| *Lentinus* | 512 | 33.3 px | 33 | 508 | genuino ✓ |
| *Polyporus* | 512 | 188.6 px | 235 | **1** | **degenerado** ✗ |
| *Polyporus* | 1024 | 188.6 px | 235 | **28** | **genuino** ✓ |

### 9.2 *squamosus* — dos problemas reales descubiertos

1. **Sub-resolución (numérico).** Forzando grilla 512, *squamosus* produce un único blob de
   escala-dominio (1 spot), NO un patrón periódico. En su grilla **por defecto (1024)** forma
   28 spots genuinos. **Corrección de un hallazgo previo:** el default de *squamosus* ya es
   1024 (no 512); el caso degenerado solo aparece si se fuerza una grilla menor. El config
   por defecto es adecuado. Los gates (fix 1/4) protegen contra el uso de grilla chica.
2. **λ del modelo subestima el target biológico (datos).** La λ intrínseca del modelo para
   *squamosus* es **~188 px, estable entre grillas**, frente a los 235 px del target
   biológico (`expected_wavelength_px`; poro 2000 µm ÷ 8.7 µm/px ≈ 230 px). Es un
   **undershoot de ~18%**, independiente de la grilla. **NO se debe "reconciliar" bajando el
   target a 188** — eso ajustaría el objetivo a la salida del modelo (tautología). Se
   documenta como limitación conocida; el 235 se mantiene como target biológico legítimo.

### 9.3 Fixes aplicados (auditoría)

- **Fix 1 — gate de patrón genuino** (`MIN_GENUINE_SPOTS=4`): *squamosus*@512 pasa a
  `validation_pass=False`; ya no valida un blob único.
- **Fix 4 — aviso de sub-resolución** (`MIN_WAVELENGTHS_IN_DOMAIN=3`, usando la λ medida):
  *squamosus*@512 → `under_resolved=True` (2.7 λ en dominio); @1024 → genuino y resuelto.
- **Fix 3 — `validation_pass` gatea la densidad**: no se reporta densidad para un patrón
  no-genuino.
- **Fix 5 — re-chequeo de dt durante el run**: guarda defensiva que detecta blow-up
  (campos no-finitos) y para limpio; no-op para runs estables.
- **Fix 2 — QC rechaza espectros sin pico: INTENTADO y REVERTIDO.** Rechazaba el pico FFT
  cuando caía en el piso de frecuencia, pero eso también descarta patrones genuinos de λ
  grande (*squamosus*@1024, 28 spots, λ≈188px correcta → caía a un autocorr peor de 519px).
  Fix 1 (gate de spots) ya cubre el caso degenerado de forma robusta, así que fix 2 era
  redundante y causaba regresión. Lección: el gate por nº de spots es mejor señal que la
  forma del pico espectral.
- Solo se tocó `integrated_validation.py` y `topology_analyzer.py`; schema canónico intacto.

### 9.4 Integridad bibliográfica

La cita central **Kuhar 2022** tenía DOI inexistente (`-00380-8`) y metadatos incorrectos
(título/autores/páginas). Verificado vía Crossref y corregido a `-00363-z` (Theory in
Biosciences 141(1):1-11). Añadido Davidson 2011 (IMA Fungus). Ver CHANGELOG.

---

## 10. Prototipo 3-nodos (P2) — resultados y diagnóstico (2026-06-02)

Implementado `core/three_node.py` (solver + `ThreeNodeGM` + `MarconNetwork` +
estabilidad lineal), `scripts/explore_three_node.py` (cobertura + barrido + control con
grilla adaptativa), 14 tests TDD. Ejecutado a grilla adaptativa (256–512).

### 10.1 Resultados de simulación

- **PROBE (GM + nodo inmóvil):** GM puro (κ=0) patrona en TODO el rango de ratios (2–3750)
  cuando la grilla resuelve. → **El ratio extremo (3750) nunca fue necesario para patronar**
  (umbral de Turing ~2); era solo para fijar la escala de λ (λ ∝ √ratio). El nodo inmóvil
  tiene efecto **no-monótono** sobre el umbral (κ=0.5 lo empeora) → afirmación débil
  (nodo baja el umbral) **NO soportada**.
- **MARCON (D_u=D_v):** la simulación reportó 4/8 "genuinos" y CONTROL 0/8.

### 10.2 Diagnóstico por estabilidad lineal — corrección importante

Para `MarconNetwork` la reacción es lineal salvo un cúbico que se anula en el punto fijo, así
que el Jacobiano **es M**. La dispersión es `eigvals(M − k²·diag(D))`. Resultado:

- **NINGÚN punto del barrido MARCON es Turing-inestable** a D igual. `growth(k=0)` pasa de
  estable a inestable *en k=0* al subir el parámetro — el sistema **se salta el régimen de
  Turing**. → Los "4/8 patrones genuinos" eran **falsos positivos**: estructuras no-lineales
  de una inestabilidad en k=0 + saturación cúbica, contadas como spots. **El análisis lineal
  cazó un falso positivo de la simulación + conteo de spots.**
- **No es bug** (solver validado a 1e-9), **no es metodología** (persiste analíticamente):
  es la **topología ad-hoc** elegida a ojo.

### 10.3 ¿Se salva la afirmación fuerte?

Búsqueda aleatoria (20 000 matrices 3×3, D_u=D_v): **~10% SON Turing-inestables** a difusión
igual con nodo inmóvil → el mecanismo de Marcon/eLife 2016 **es real y alcanzable**; la
topología ad-hoc era mala elección. **Pero** las λ grandes (~188 px, escala *squamosus*) son
**raras** a D igual (distribución sesgada a λ chica) — la disparidad de difusión es lo que
facilita escala grande. Punto científico en sí.

### 10.4 Estado de P2 y siguiente paso

- **Soportado:** equal-D Turing-con-nodo-inmóvil existe (10% de topologías).
- **NO demostrado aún:** que *esta* familia genere el morfo-espacio fúngico con control de
  escala por un parámetro. El prototipo actual NO lo demuestra.
- **Siguiente paso riguroso:** (1) buscar M Turing-inestable a D igual vía
  `fastest_growing_wavelength`; (2) elegir el parámetro que mueve λ* monótono —verificado
  **analíticamente antes de simular**; (3) simular solo puntos con banda de Turing confirmada;
  (4) caracterizar si las λ grandes son accesibles a D igual o requieren disparidad moderada.
- **Caveat técnico:** `fastest_growing_wavelength` usa un `k_max` finito; λ* pequeñas pueden
  estar limitadas por el borde del scan (ampliar `k_max` para λ chicas).

### 10.6 Búsqueda de topología y ley de escala (resultado fuerte de P2)

Tras arreglar `fastest_growing_wavelength` (auto-ensanchado del scan; NaN honesto si no hay
selección de escala finita), se exploró analíticamente (sin simular) el espacio de topologías:

- **La biología "obvia" NO funciona.** Estructura de signos plausible (u activador, v
  inhibidor, w=densidad hifal *dirigida por u* e inhibiendo u): **0/8000** Turing a D igual.
- **Qué cableado SÍ funciona.** Enumerando los 3⁵=243 patrones de signo del nodo inmóvil:
  **46/243** admiten equal-D Turing. Los más robustos acoplan w al **inhibidor** (`w↔v`), con
  **w desacoplado de u**. Interpretación: el nodo inmóvil debe modular el inhibidor para crear
  la inhibición de largo alcance que normalmente aporta la disparidad de difusión.
- **Control de escala monótono: SÍ existe.** En una topología `w↔v` Turing-inestable, varios
  parámetros (M[0,1], M[1,0], M[1,2], M[2,1]) mueven λ* de forma monótona (verificado
  analíticamente). → afirmación de control por un parámetro **soportada** para la topología
  correcta.
- **Ley de escala — el hallazgo clave.** λ* sigue `√D` a difusión IGUAL:
  λ*≈3.8 px (D=0.5) → 41.8 px (D=50) → 125 px (D=500). **La escala biológica (33–235 px) se
  alcanza a D_u=D_v subiendo la D absoluta; la disparidad NO es necesaria para λ grande.**
- **Dos mecanismos distintos.** En la topología `w↔v`, subir la disparidad D_v/D_u **destruye**
  la banda de Turing (lo opuesto a Gierer-Meinhardt, que la requiere). GM fija escala por
  `√(ratio)`; la red de nodo inmóvil la fija por `√(D absoluta)`.

**Conclusión que disuelve la tautología:** la longitud de onda del himenóforo la fija la
**longitud de difusión `√(D/reacción)`, no el *ratio* de difusión.** El `D_v/D_u=3750` del
modelo original no es fundamental — es un proxy de longitud de difusión, alcanzable también sin
disparidad (con un nodo inmóvil que module el inhibidor). Esto convierte el "parámetro
sospechoso" en uno de (al menos) dos mecanismos que fijan la misma escala. **Es el resultado
publicable de P2** (pendiente: confirmar por simulación no-lineal a D grande / grilla grande;
el análisis lineal es el predictor riguroso y el solver está validado contra él).

### 10.7 Confirmación no-lineal (opción C) — bucle cerrado

La predicción lineal (§10.6) se confirmó **simulando** la topología `w↔v` a difusión igual
(`scripts/confirm_three_node.py`). Para que la simulación produjera patrones acotados se
corrigió `MarconNetwork`: la saturación cúbica ahora actúa en **los tres** nodos (antes solo
en u, dejando v/w crecer sin cota). El cúbico se anula en el punto fijo, así que el análisis
lineal (J=M) no cambia.

Resultado (grilla escalada ∝√D para comparar el patrón adimensional):

| D (igual) | spots | spacing (px) | λ* lineal (px) | spacing/√D |
|-----------|-------|--------------|----------------|------------|
| 5  | 16 | 24–32 | 12.5 | ~14 |
| 20 | 14–16 | 51–64 | 25.1 | ~14 |

- **`spacing/√D` ≈ constante** → el patrón no-lineal escala como **√D**, igual que la
  predicción lineal. Bucle lineal→no-lineal **cerrado**.
- El spacing no-lineal está un factor constante (~2.5×) por encima de λ* lineal (coarsening);
  irrelevante para la ley de escala (ambos ∝√D).

**Factor de coarsening (caracterizado):** spacing/λ* ≈ **2.4, constante** — independiente de D
(2.36→2.48 de D=5 a 20) **y de γ** (2.32 vs 2.36 a γ=1.0 vs 0.3). Que no dependa de la fuerza
de la no-linealidad indica que la selección no-lineal de longitud de onda la fija la geometría
de la banda de dispersión, no la magnitud del cúbico. → la ley √D transfiere lineal→no-lineal
con solo un offset multiplicativo constante.

**Regla unificadora de topología (condición necesaria):** enumerando las estructuras de signo
Turing-capaces a D igual, **el 100% acopla el nodo inmóvil al inhibidor v** (`w→v` o `v→w`≠0);
ninguna funciona sin ese acoplamiento. El acoplamiento al activador u NO es necesario (90%; hay
casos que patronan sin tocar u) y `w_self` es irrelevante (~51%). Mecánicamente: el nodo
inmóvil suple la inhibición de largo alcance que normalmente aporta el inhibidor rápido vía
disparidad, así que **debe actuar a través del inhibidor**. (Esto explica por qué la "biología
obvia" —w dirigido por u, w inhibe u— estaba en la clase 0%.)

**Conclusión P2 (analítica + no-lineal):** una red de 3 nodos con difusión igual y el nodo
inmóvil acoplado al inhibidor genera patrones de himenóforo genuinos cuya escala se fija por
`√(D absoluta)` — alcanzando el rango biológico **sin disparidad de difusión**. La
longitud de onda la fija la longitud de difusión, no el ratio; el `D_v/D_u=3750` del modelo
original es un proxy de esa longitud, no un requisito fundamental.

> **Pendiente (run para figura del paper):** la confirmación no-lineal se hizo a escala
> *brumalis/fomes* (grid ≤640). La validación a **escala *squamosus* (~230 px ≈ 2000 µm)**
> exige **grilla grande (~1408) y D≈270** — un hongo real necesita grid grande. Es el run
> pesado pendiente (`confirm_three_node.py --D 270 --grid 1408`, ~1–2 h por el CFL explícito).
> La ley √D predice ~16 spots con spacing ~230 px; ese punto cierra la figura spacing-vs-√D
> abarcando las 3 especies.

### 10.5 Lección de método

La simulación + conteo de spots puede producir **falsos positivos de patronamiento**. El
análisis de estabilidad lineal es el filtro barato y decisivo: **verificar la banda de Turing
analíticamente ANTES de declarar un patrón genuino.**

---

*Documento de estrategia. Informado por auditoría de código de junio 2026, verificación
empírica (ejecución del pipeline + análisis de estabilidad lineal) y revisión bibliográfica
vía academic-research-skills.*
