# Diseño: reformulación de Turing de 3 nodos (P2)

**Fecha:** 2026-06-02
**Autor:** Mario Ahumada Durán (con Claude Code)
**Contexto:** `docs/ESTRATEGIA_TESIS.md` §4 (P2)

---

## 1. Objetivo y afirmación

**Problema que resuelve:** la calibración actual de FungalMorphoSpace es tautológica
(3 parámetros — `D_v/D_u`, `rho`, `b` — tuneados por especie hasta calzar la λ observada,
y validados contra esa misma observación). Además exige `D_v/D_u` extremo (hasta 3750),
que un revisor lee como "parámetro inventado".

**Hipótesis central (falsable):**

> El morfo-espacio de himenóforos poliporo (rango de λ ≈ 46–235 px, topologías
> spots→laberinto) puede generarse con una red de Turing de 3 nodos —dos morfógenos
> difusibles con difusividades biológicamente plausibles más un nodo **no-difusible**
> (densidad hifal)— recorriendo el rango con **menos grados de libertad** que el modelo
> GM de 2 campos, y sin requerir disparidad de difusión extrema.

**Falsabilidad:** si NINGUNA topología de 3 nodos con `D_v≈D_u` reproduce el rango de λ,
o si reproducirlo exige tantos o más parámetros libres que el GM de 2 campos, la
reformulación falla y la disparidad extrema es irreducible (resultado negativo igualmente
publicable).

**Lo que P2 NO afirma:** predicción de escala absoluta de poro a partir de biofísica
(bloqueado por datos — ver estrategia §2). Eso es P3/MCAT.

---

## 2. Criterio de éxito — grados de libertad, con control

El éxito NO es "cobertura de rango" a secas (sería suficiencia débil: "otro modelo que
también ajusta"). El éxito se define por **reducción de grados de libertad**, con un
control que da contraste:

1. **Cobertura de morfo-espacio (modelo vs modelo).** El 3-nodos alcanza el mismo rango de
   (λ, topología) que el GM de 2 campos, con `D_v≈D_u`. Comparación modelo-vs-modelo →
   sortea el bloqueo de datos; la crítica de tautología no aplica (no afirmamos predecir
   biología, sino suficiencia mecanística).

2. **Barrido de 1 parámetro → mapa λ(p).** Mostrar que **un solo** parámetro biológicamente
   plausible recorre el rango de λ observado (≈46→235 px) de forma monótona, frente a los
   **3** parámetros que GM tunea por especie. Estilo Sheth/Hox (un parámetro tunea la
   wavelength).

3. **Control (topología barajada).** Una red de 3 nodos con la misma estructura pero signos
   de interacción aleatorizados/barajados **NO** cubre el morfo-espacio. Sin este contraste,
   un resultado positivo es trivial. El control es parte del experimento, no opcional.

**Métrica de "grados de libertad":** nº de parámetros que hay que variar para recorrer el
rango de λ. GM = 3 (por especie). Meta 3-nodos = 1.

---

## 3. Arquitectura — aislada del motor validado

El motor de 2 campos (`TuringSimulator` + `GiererMeinhardtKinetics`) **no se modifica**:
P1/tool paper depende de que siga estable. El prototipo vive en módulos nuevos.

### 3.1 Archivos nuevos

- **`src/fungalmorphospace/core/three_node.py`**
  - `ThreeNodeKinetics(ABC)` — interfaz paralela a `TuringKinetics`, con
    `f(u,v,w)`, `g(u,v,w)`, `h(u,v,w)` (reacciones de activador, inhibidor, nodo inmóvil).
  - `ThreeNodeGM(ThreeNodeKinetics)` — **sonda**: GM existente + nodo hifal `w`.
  - `MarconNetwork(ThreeNodeKinetics)` — **fuerte**: red parametrizada por matriz de
    interacción, diseñada para `D_v≈D_u`.
  - `ThreeNodeSimulator` — gestiona 3 campos `(u, v, w)`. `w` tiene `D_w=0` (sin término
    difusivo). Reimplementa el stencil Laplaciano de 5 puntos **inline** (es una sola
    expresión que solo depende de `dx`); NO se refactoriza `TuringSimulator` para no
    desestabilizar el path de envío.
  - CFL: `dt_max = dx²/(4·max(D_u,D_v))` (D_w=0 no restringe). **Re-chequeo de dt durante
    el run** (ver §5, modificación E) — no heredar el dt-fijo del motor viejo.

- **`scripts/explore_three_node.py`** — corre cobertura + barrido + control, vuelca
  figuras/tablas a `results/three_node/`.

- **`tests/test_three_node_solver.py`** — test de regresión numérico (ver §4).

### 3.2 Reuso de medición

λ por FFT y topología (Euler/spots) se miden con `TopologyAnalyzer` existente sobre el
campo `u` — sin cambios. La cascada QC FFT→autocorr ya está validada.

---

## 4. Validación del solver (correctitud, independiente de biología)

Antes de construir ciencia encima, el solver de 3 nodos debe demostrarse correcto:

- **Validación contra relación de dispersión lineal.** Para un set de parámetros, calcular
  analíticamente la longitud de onda del modo de crecimiento más rápido (linealización en
  torno al estado estacionario → autovalores de la matriz Jacobiana + difusión). Verificar
  que la λ medida por FFT coincide con la λ teórica dentro de tolerancia (p.ej. ±15%).
- **Test de regresión** (`test_three_node_solver.py`): fija parámetros y semilla, corre a
  estado estacionario, afirma λ_medida ≈ λ_teórica. Esto valida el motor nuevo Y empieza a
  cerrar el hueco de "cero tests numéricos" que halló la auditoría.
- **Conservación/sanidad:** sin NaN/Inf, campos acotados, `w≥0`.

---

## 5. Los dos modelos (detalle)

### 5.1 Sonda — `ThreeNodeGM` (time-boxed)

Sobre el GM existente, `w` = densidad hifal como **nodo dinámico lento real** (NO esclavo
algebraico — modificación D):

```
∂u/∂t = D_u ∇²u + ρ(a − b·u + u²/(v+ε))           (activador, difunde)
∂v/∂t = D_v ∇²v + ρ(u² − v)                         (inhibidor, difunde)
∂w/∂t = (1/τ_w)·(σ(u) − w)                          (densidad hifal, D_w=0, τ_w grande)
```

con **feedback de `w` hacia la dinámica** (sin esto `w` es decorativo): `w` modula la
saturación/inhibición local, p.ej. el término de inhibición efectivo pasa a `(v + κ·w)`.
`τ_w` grande = escala lenta de crecimiento de tejido.

**Experimento (time-boxed):** barrer `D_v/D_u` hacia abajo desde 3750 variando `κ` (fuerza
de acople de `w`) y ver hasta dónde baja el ratio manteniendo patrón. **Si en un barrido
pequeño el ratio no baja de forma significativa, pivotar de inmediato a Marcon** — la sonda
es diagnóstica, no un estudio completo. Resultado negativo ("GM es irreducible") es
informativo y se reporta.

### 5.2 Fuerte — `MarconNetwork`

Topología de 3 nodos del catálogo eLife 2016 (Marcon et al.) con **`D_u = D_v`** (difusión
igual) + `w` inmóvil. Parametrizada por matriz de interacción 3×3 (signos/fuerzas), con un
estado estacionario homogéneo y condiciones de Turing-con-nodo-inmóvil satisfechas.
Realización biológica análoga a Bmp-Sox9-Wnt (Raspopovic 2014), con `w`=hifa como el nodo
inmóvil (rol de Sox9).

**Barrido:** identificar el único parámetro plausible que tunea λ y mapear λ(p) sobre
46→235 px.

**Control:** misma estructura con signos de interacción barajados → verificar que NO
patterna o no cubre el rango.

---

## 6. Alcance (YAGNI)

**Incluido:** dos modelos 3-nodos, simulador aislado, validación de solver + test,
experimentos de cobertura/barrido/control, figuras/tablas en `results/three_node/`.

**Excluido explícitamente:**
- 4º nodo / matriz de pared → P3/MCAT.
- Anclaje biofísico (`BiophysicalAnchor`) → bloqueado por datos.
- Modificar el path de 2 campos (`TuringSimulator`, `GiererMeinhardtKinetics`).
- Integrar 3-nodos al pipeline de validación principal → se decide DESPUÉS de ver
  resultados; el prototipo es autónomo vía su script.

---

## 7. Relación con los bugs de auditoría

Modificación E: los fixes de correctitud del solver que afectan también al motor nuevo se
abordan **al construir** el 3-nodos (re-chequeo de dt, robustez de convergencia), no se
heredan. Los fixes del path viejo (`qc_pass` no gatea, *squamosus* sub-resuelto en grilla)
son independientes y se tratan aparte, después.

---

## 8. Entregables

1. `core/three_node.py` (simulador + 2 kinetics + ABC).
2. `scripts/explore_three_node.py` (cobertura + barrido + control).
3. `tests/test_three_node_solver.py` (regresión vs dispersión lineal).
4. Figuras/tablas en `results/three_node/`.
5. Un breve `results/three_node/FINDINGS.md` con el veredicto: ¿bajó el ratio? ¿un
   parámetro recorre el rango? ¿el control falla como se espera?
