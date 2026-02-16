# Calibración del Modelo para *Lentinus brumalis* ("The Tuning Journey")

**Proyecto:** FungalMorphoSpace  
**Propósito:** Anexo metodológico para tesis y documentación de experimentos  
**Fecha:** Enero 2026

---

## 1. Contexto del problema

El objetivo fue encontrar una configuración de parámetros que representara a *Lentinus brumalis* como una especie de transición ("meso-poro") entre dos anclas morfológicas del pipeline:

- **Micro-referencia:** *Fomes fomentarius* (patrón denso; longitud de onda simulada típica λ ≈ 46 px).
- **Macro-referencia:** *Polyporus squamosus* (poros gigantes; λ ≈ 230 px).

La hipótesis operativa fue que *L. brumalis* debería ubicarse en una zona intermedia del morfoespacio, con λ objetivo en el rango ≈ 60–70 px, manteniendo un **régimen de spots** con alta población (textura de esponja) y evitando el colapso a pocos poros gigantes.

---

## 2. Principio físico de diseño

La calibración se guía por un balance entre dos fuerzas efectivas:

- **Fuerza expansiva (Difusión):** controlada principalmente por el ratio \(D_v/D_u\). Al incrementarlo, el campo inhibidor se propaga más lejos y aumenta la escala espacial típica del patrón.
- **Presión de confinamiento (Reacción/Metabolismo):** controlada por \(\rho\) (reescalamiento temporal). Aumentar \(\rho\) acelera la reacción frente a la difusión efectiva, favoreciendo spots más pequeños y densos; disminuir \(\rho\) permite que los spots “se relajen” y expandan.

En esta interpretación, la zona “Goldilocks” de *L. brumalis* ocurre cuando la difusión es suficiente para superar la barrera de ~50 px (observada en varios intentos), pero la cinética mantiene una alta población de spots (evitando gigantismo).

---

## 3. Bitácora de experimentación (Chronology of Tuning)

### Fase A: Intento conservador (fallo por exceso de densidad)

**Run A:** \(D_v/D_u = 250\), \(\rho = 0.4\), \(b = 1.0\)  
**Resultado:** \(λ \approx 33\) px (menor que *Fomes*).  
**Diagnóstico físico:** \(\rho\) demasiado alto respecto a la difusión (la reacción “gana”), forzando spots muy pequeños y apretados.  
**Conclusión:** morfología biológicamente poco plausible para una especie intermedia (riesgo de solapamiento visual con *Fomes*).

### Fase B: Salto a alta difusión (fallo por gigantismo)

**Run B:** \(D_v/D_u = 2000\), \(\rho = 0.1\), \(b = 2.0\)  
**Resultado:** \(λ \approx 113\) px (pocos spots).  
**Diagnóstico físico:** el aumento extremo de difusión, combinado con mayor degradación (\(b\)), elimina poros intermedios y deja pocos poros gigantes aislados.  
**Conclusión:** deriva hacia un fenotipo tipo *P. squamosus*; se pierde la textura de “esponja”.

### Fase C: Barrera de estabilidad (estancamiento)

**Run C:** \(D_v/D_u = 1400\), \(\rho = 0.2\), \(b = 1.0\)  
**Resultado:** \(λ \approx 51\) px.  
**Diagnóstico físico:** al mantener \(\rho\) alto (igual que *Fomes*), el sistema tiende a un mínimo local alrededor de ~50 px: los poros se resisten a crecer más aun cuando aumenta la difusión.  
**Conclusión:** confirma una región de estabilidad del régimen de spots cerca de 50 px bajo \(\rho\) relativamente alto.

### Fase D: Contradicción biológica (retroceso)

**Run D:** \(D_v/D_u = 250\), \(\rho = 0.2\), \(b = 1.0\)  
**Resultado:** \(λ \approx 43.5\) px.  
**Diagnóstico:** al volver a difusión baja, incluso con metabolismo más lento, el patrón no alcanza poros medios.  
**Conclusión:** con difusión baja, el modelo no produce *meso-poros*; se revalida la necesidad de aumentar \(D_v/D_u\).

### Fase E: Convergencia final (éxito: “Goldilocks zone”)

**Run Final:** \(D_v/D_u = 1000\), \(\rho = 0.15\), \(b = 1.0\)  
**Resultado:** \(λ = 62.85\) px (\(~1000+\) spots en grid 1024).  

**Justificación del éxito:**

- **Ratio 1000:** provee la fuerza expansiva necesaria para romper la barrera de ~50 px sin colapsar a gigantismo.
- **\(\rho = 0.15\):** un metabolismo ligeramente más lento que *Fomes* (0.2) permite que los poros se expandan sin romper el régimen de spots.
- **\(b = 1.0\):** mantiene la supervivencia de una población alta de spots, conservando la textura densa.

Interpretación central: \(\rho\) opera como un **punto de inflexión metabólico**; al disminuir bajo cierto umbral, el sistema deja de favorecer poros pequeños y pasa a estructurar poros medios con estabilidad.

---

## 4. Tabla comparativa final (para insertar en LaTeX)

> Nota: los valores se presentan como el set “definitivo” para la narrativa de tesis. En el software, cambiar defaults implica re-baselinear validaciones y tests; ver §5.

```latex
\begin{table}[ht]
\centering
\caption{Par\'ametros calibrados y longitud de onda simulada (\emph{wavelength}) para las tres anclas del morfoespacio.}
\label{tab:calibracion_morfoespacio}
\begin{tabular}{lccccc}
\toprule
\textbf{Especie} & $D_v/D_u$ & $\rho$ & $b$ & $\lambda$ (px) & \textbf{Interpretaci\'on} \\
\midrule
\textit{F. fomentarius} & 150.0 & 0.20 & 1.0 & 46.0 & Referencia Micro \\
\textit{L. brumalis} & 1000.0 & 0.15 & 1.0 & 62.9 & Transici\'on Meso \\
\textit{P. squamosus} & 3750.0 & 5.00 & 3.0 & 230.0 & Referencia Macro \\
\bottomrule
\end{tabular}
\end{table}
```

---

## 5. Nota metodológica: ¿conviene cambiar los defaults del software?

Cambiar los defaults de *L. brumalis* a \(D_v/D_u=1000\), \(\rho=0.15\) es una decisión **de alto impacto** porque:

1. Re-define el “gold standard” y modifica la salida canónica (CSV/figuras) del pipeline.
2. Exige re-ejecutar robustez (múltiples semillas) para reportar media/SD en el nuevo punto.
3. Puede tensionar los rangos biológicos declarados en documentación previa (p.ej. densidad esperada), por lo que requiere una justificación explícita (incertidumbre de medición, variabilidad intra-específica, o criterio de separabilidad visual en tesis).

Recomendación práctica:

- Tratar este set como **preset de tesis** (o “candidate calibration”) hasta completar robustez y contrastación con micrografías.
- Si se adopta como default en el repo, hacerlo con:
  - bump de versión,
  - actualización de documentación contractual, y
  - actualización de smoke tests / expectativas.

