# TEXTO PARA TESIS — MÉTODOS Y LIMITACIONES
## Material Listo para Copiar y Adaptar

**Autor:** Mario Ahumada Durán  
**Proyecto:** FungalMorphoSpace  
**Fecha:** Enero 2026  
**Versión:** 5.0 (PATCHED_v4_3)

---

## 📋 1) MÉTODOS

### 1.1. Modelo Reacción-Difusión y Parametrización

La morfogénesis del himenóforo se modeló mediante un sistema de reacción-difusión tipo Gierer-Meinhardt en dos campos u(x,t) y v(x,t), donde u representa un activador local y v un inhibidor de largo alcance. El sistema general adopta la forma:

```
∂u/∂t = D_u∇²u + ρF(u,v)
∂v/∂t = D_v∇²v + ρG(u,v)
```

El parámetro ρ escala la cinética (intensidad efectiva de reacción), mientras D_u y D_v gobiernan la difusión. Se exploró explícitamente el cociente D_v/D_u como descriptor geométrico del régimen difusivo, manteniendo b=1.0 normalizado en las comparaciones finales (cuando aplica).

**Funciones cinéticas específicas (Gierer-Meinhardt):**
```
F(u,v) = u²/v - u + a
G(u,v) = u² - bv
```

---

### 1.2. Reescalamiento Temporal y Difusión Efectiva

Dado que ρ multiplica los términos cinéticos, el sistema puede analizarse mediante el cambio de variable temporal τ = ρt, obteniéndose:

```
∂u/∂τ = (D_u/ρ)∇²u + F(u,v)
∂v/∂τ = (D_v/ρ)∇²v + G(u,v)
```

En consecuencia, la escala espacial seleccionada (y por tanto la longitud de onda λ) depende del balance completo reacción-difusión a través de las difusiones efectivas:

```
D_u* = D_u/ρ
D_v* = D_v/ρ
```

**Este punto es central para interpretar casos donde λ disminuye al aumentar D_v/D_u si, simultáneamente, ρ aumenta.** Esta es la resolución de la "paradoja de Brumalis".

---

### 1.3. Control Temporal y Comparabilidad

#### Problema identificado

El ajuste automático de dt por estabilidad numérica (criterio CFL) producía tiempos físicos finales distintos entre especies si el número de pasos (steps) se mantenía fijo. Esto introducía un sesgo de "stopping-time" en las comparaciones.

#### Solución implementada (v4.3)

Se introdujo el parámetro `T_target` (tiempo físico objetivo) y se deriva el número de pasos desde:

```
steps = ceil(T_target / dt_final)
```

Donde `dt_final` es el paso temporal ajustado por estabilidad. Esto garantiza que todas las especies corran el mismo tiempo físico (T_target = 5.0 por defecto), eliminando el sesgo temporal.

---

### 1.4. Estimación de la Longitud de Onda λ

La longitud de onda λ se estimó mediante dos métodos complementarios:

#### Método 1: Espectro de Potencia Radial (FFT)

1. Calcular transformada de Fourier 2D del campo u
2. Obtener espectro de potencia |FFT|²
3. Promediar radialmente para obtener P(k) vs |k|
4. Excluir frecuencias ultra-bajas (menos de 2 ciclos por dominio)
5. Suavizar perfil radial
6. Detectar pico dominante k_peak
7. λ_FFT = 1 / k_peak

**Control de calidad:** Se verifica que el pico tenga prominencia suficiente (ratio ≥ 3.0 respecto al fondo) y que λ no sea comparable al tamaño del dominio. Si no pasa QC, λ_FFT se marca como no confiable.

#### Método 2: Autocorrelación Radial

1. Calcular autocorrelación 2D normalizada del campo u
2. Promediar radialmente para obtener C(r)
3. Suavizar perfil radial
4. Identificar primer mínimo local (anti-correlación)
5. λ_autocorr = posición del primer máximo local después del primer mínimo

Este método es más robusto para patrones con ruido o periodicidad imperfecta.

**Unidades:** Cuando dx=1.0, se reporta en píxeles (px). Cuando dx≠1.0, se aplica calibración a µm.

---

### 1.5. Métricas Topológicas y Morfométricas

#### Definición de "spots" (componentes)

**Problema anterior:** La detección de "picos" como máximos locales del campo continuo producía sobreconteo masivo (n≈4000-5000) en patrones ruidosos, haciendo las métricas de spacing no interpretables.

**Solución implementada (v4.3):** Los spots se definen como **componentes conexas** del patrón binarizado:

1. Normalizar campo u al rango [0, 1]
2. Binarizar con umbral (por defecto 0.5)
3. Etiquetar componentes conexas (8-conectividad)
4. N_spots = número de componentes con área > área_mínima

#### Espaciamiento (spacing)

El espaciamiento medio se calcula sobre los **centroides** de las componentes conexas:

1. Calcular centroide de cada componente
2. Construir KD-tree con centroides
3. Para cada centroide, encontrar distancia al vecino más cercano
4. mean_spacing = promedio de distancias NN
5. CV_spacing = std / mean (coeficiente de variación)

#### Característica de Euler (χ)

```
χ = N_componentes - N_agujeros
```

Calculada sobre el patrón binarizado. χ > 0 indica estructura tipo "spots"; χ ≈ 0 indica laberinto; χ < 0 indica red percolada con múltiples agujeros.

---

### 1.6. Calibración Espacial

La conversión de píxeles a unidades físicas se basa en la calibración con Fomes fomentarius:

```
λ_simulado(Fomes) = 46 px
λ_reportado(literatura) ≈ 400 µm

Factor de calibración = 400 / 46 = 8.695652 µm/px ≈ 8.7 µm/px
```

**Fuente de calibración:** Klemm et al. (2024), mediciones SEM de espaciamiento de poros.

---

### 1.7. Robustez: Semillas y Sensibilidad a Umbral

#### Análisis de robustez por semillas

Para descartar artefactos por inicialización, se ejecutan n≥20 simulaciones con semillas aleatorias independientes, reportando:
- Media ± desviación estándar de λ
- Coeficiente de variación (CV)
- Histograma de distribución

**Criterio de aceptación:** CV < 5% para métrica primaria.

#### Sensibilidad a umbral de binarización

Se evalúa la estabilidad de N_spots y spacing frente a variaciones del umbral:
- Umbral bajo: 0.4
- Umbral medio: 0.5 (default)
- Umbral alto: 0.6

Las conclusiones principales deben ser invariantes al umbral dentro de este rango.

---

## ⚠️ 2) LIMITACIONES Y AMENAZAS A LA VALIDEZ

### 2.1. Identificabilidad Paramétrica (Sistema Multivariable)

En régimen no lineal, múltiples combinaciones de parámetros (D_v/D_u, b, ρ, constantes cinéticas) pueden producir patrones visualmente similares. Por ello:
- La inferencia causal estricta de un único "controlador" requiere controles adicionales
- Se reporta análisis de sensibilidad paramétrica
- Los parámetros se interpretan como "efectivos", no como valores biológicos directos

### 2.2. Sensibilidad al Umbral en Métricas Discretas

El conteo de componentes conexas depende de la binarización. Aunque la metodología v4.3 es más robusta que la detección de máximos locales, el umbral sigue siendo un hiperparámetro. Por esta razón:
- λ se prioriza como métrica primaria (menos sensible)
- N_spots se reporta con el umbral explícito
- Se incluye análisis de sensibilidad al umbral

### 2.3. Discretización y Condiciones de Borde

- λ reportado en píxeles depende de resolución espacial
- Condiciones de borde periódicas imponen restricciones a longitudes de onda permitidas
- Comparaciones absolutas con medidas físicas requieren calibración explícita px→µm

### 2.4. Traducción Patrón→Poro Biológico

El campo simulado representa una abstracción morfogenética. La correspondencia exacta entre máximos/mínimos del campo y geometría real del poro requiere:
- Procedimiento de segmentación morfológicamente consistente
- Calibración con imágenes empíricas (SEM, µCT)
- Validación cruzada con datos experimentales (Klemm et al., 2024)

### 2.5. Ley de Escalamiento con N=3 Especies

Con solo 3 puntos de datos, la regresión log-log (R²=0.9997) tiene alto poder predictivo pero bajo poder estadístico para inferir exponentes. Se interpreta como:
- Consistencia cualitativa con la teoría
- No evidencia fuerte de exponente específico
- Validación requiere más especies en trabajos futuros

### 2.6. Modelo 2D

Las simulaciones son estrictamente 2D. La extensión a 3D (tubos reales) introduciría:
- Efectos de curvatura
- Anisotropía estructural
- Mayor costo computacional

---

## 📊 3) PIES DE FIGURA

### Figura X. Comparación de Escalas Espaciales entre Especies

**Figura X.** Longitud de onda estimada λ (en píxeles) para patrones simulados de himenóforo en tres especies representativas. El baseline *Fomes fomentarius* presenta λ≈44–46 px. Para *Lentinus brumalis* se fijó la configuración D_v/D_u=250, b=1.0, ρ=0.4, obteniéndose λ≈33 px (n=5). La relación entre escala y parámetros se interpreta con el reescalamiento temporal τ=ρt, donde la difusión efectiva escala como D/ρ, permitiendo ajustar la densidad de poros sin alterar la topología (spots). *Polyporus squamosus* (D_v/D_u=3750, ρ=5.0) exhibe λ≈190–235 px, cerrando el rango de escalas.

---

### Figura S1 (Apéndice). Sensibilidad del Conteo de Componentes al Umbral

**Figura S1.** Sensibilidad del número de componentes conexas al umbral de binarización (0.4, 0.5, 0.6) aplicado al campo normalizado. La metodología v4.3 (componentes conexas) muestra estabilidad relativa comparada con detección de máximos locales. La conclusión principal se fundamenta en λ como métrica primaria.

---

### Figura S2 (Apéndice). Robustez por Semillas Aleatorias

**Figura S2.** Distribución de λ (mejor estimador: FFT con QC, respaldo por autocorrelación) para n=5 simulaciones independientes de *L. brumalis* con la configuración oficial (D_v/D_u=250, b=1.0, ρ=0.4). Se observa estabilidad de escala y variación moderada entre semillas, consistente con un régimen robusto de spots.

---

### Figura S3 (Apéndice). Validación de Control Temporal

**Figura S3.** Comparación de tiempo físico final (T_final) entre especies. Panel A: metodología anterior (steps fijo) resulta en tiempos desiguales. Panel B: metodología v4.3 (T_target) garantiza T_final ≈ 5.0 para todas las especies, eliminando sesgo de stopping-time.

---

## 📝 4) NOTAS DE USO

### Para Métodos
- Secciones 1.1-1.7 van en el capítulo de Métodos
- Mantener ecuaciones en formato LaTeX para el documento final
- Sección 1.2 (reescalamiento) es CLAVE para defender la paradoja de Brumalis
- Sección 1.3 (control temporal) es NUEVA en v4.3

### Para Limitaciones
- Incluir en sección de Discusión o Limitaciones
- Demuestra honestidad científica y pensamiento crítico
- Previene preguntas incómodas del tribunal
- La limitación 2.5 (N=3) es importante reconocer

### Para Figuras
- Usar pies de figura como template
- Generar Figuras S1, S2, S3 para apéndice
- S3 demuestra la mejora metodológica v4.3

---

## ✅ VENTAJAS DE ESTA METODOLOGÍA

1. ✅ **Métricas robustas** - Componentes conexas eliminan sobreconteo
2. ✅ **Control temporal** - T_target garantiza comparabilidad
3. ✅ **Doble estimación λ** - FFT + autocorrelación con QC
4. ✅ **Calibración documentada** - 8.7 µm/px trazable a Klemm et al.
5. ✅ **Limitaciones explícitas** - Blindaje contra críticas
6. ✅ **Reproducible** - Código versionado y documentado

---

## 📚 REFERENCIAS METODOLÓGICAS

1. Turing, A. M. (1952). The chemical basis of morphogenesis. *Phil. Trans. R. Soc. B.*
2. Gierer, A. & Meinhardt, H. (1972). A theory of biological pattern formation. *Kybernetik.*
3. Klemm, D. et al. (2024). Hierarchical structure... Fomes fomentarius. *PLOS ONE.*
4. Kuhar, F. et al. (2022). Pattern formation features might explain homoplasy. *Theory in Biosciences.*

---

**© 2025-2026 Mario Ahumada Durán. Todos los derechos reservados.**

**Preparado por:** Mario Ahumada Durán  
**Proyecto:** FungalMorphoSpace v0.5.1  
**Status:** ✅ Listo para incluir en tesis
