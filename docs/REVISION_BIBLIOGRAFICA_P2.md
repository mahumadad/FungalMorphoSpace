# Revisión bibliográfica — Nodo inmóvil / Turing a difusión igual en himenóforos fúngicos

**Mario Ahumada Durán — Junio 2026**
**Método:** academic-research-skills `deep-research`, modo lit-review. Todas las referencias
verificadas en Crossref (DOI + autores + venue); ninguna en zona gris. Asistido por IA.

> Documento de honestidad para el encuadre de P2. Conclusión incómoda pero necesaria:
> el **mecanismo** ya está publicado (incluso aplicado); el aporte real es más estrecho.

---

## Veredicto en una frase

El mecanismo de 3 nodos (nodo inmóvil de binding → Turing a difusión igual) **ya está
publicado, y ya fue aplicado a sistemas biológicos reales con validación** (Nesterenko et al.
2017). El aporte de este trabajo **no es el mecanismo** sino: (1) la primera aplicación al
**himenóforo fúngico**, (2) la **crítica cuantitativa** de que el `D_v/D_u=3750` del modelo
previo es innecesario, y (3) una **herramienta reproducible**. Es un aporte de *aplicación +
crítica + software*, no de teoría nueva.

---

## HILO 1 — Estado del arte y novedad

### Lo que YA está publicado (el mecanismo)

| Referencia | Tipo | Qué hace |
|---|---|---|
| [Marciniak-Czochra, Karch & Suzuki 2016, *J. Math. Biol.*](https://doi.org/10.1007/s00285-016-1035-z) | Teoría pura | RD-ODE (nodo inmóvil); teorema: con nodo inmóvil **autocatalítico no existen patrones estables** (solo spikes). |
| [Härting, Marciniak-Czochra & Takagi 2017, *DCDS*](https://doi.org/10.3934/dcds.2017032) | Teoría | Con **histéresis** sí hay patrones estables (discontinuos). |
| [Korvasová, Gaffney, Maini, Ferreira & Klika 2015, *J. Theor. Biol.*](https://doi.org/10.1016/j.jtbi.2014.11.024) | Teoría | Condiciones de Turing con **sustrato inmóvil de binding**; reducción quasi-steady → ratio de difusión efectivo → Turing con D iguales. |
| [Nesterenko, Kuznetsov & Korotkova 2017, *PLOS ONE*](https://doi.org/10.1371/journal.pone.0171212) | **Teoría + aplicación + validación** | **Precedente directo.** Gierer-Meinhardt extendido a **3 ecuaciones con sitios de binding inmóviles**, Turing con **difusión igual**, **aplicado y validado** contra folículo piloso WNT/DKK, manchas de pez gato y placa neural de Xenopus, con parámetros realistas. |
| Raspopovic, Marcon, Russo & Sharpe 2014, *Science* | Aplicación | Red Bmp-Sox9-Wnt en dígitos; **Sox9 = nodo no-difusible**. |
| Marcon, Diego, Sharpe & Müller 2016, *eLife* | Teoría | Catálogo de redes de 3 nodos que dan Turing con señales de igual difusión. |
| [Madzvamuse, Chung & Venkataraman 2015](https://royalsocietypublishing.org/) (bulk-surface) | Teoría | Intercambio diferencial por interfaz reemplaza la disparidad de difusión. |
| λ ∝ √(D/k) (longitud de difusión) | Libro de texto | Estándar; la tasa de reacción ajusta la escala en órdenes de magnitud. |

### Veredicto de novedad, por afirmación

| Afirmación de P2 | ¿Novedosa? |
|---|---|
| Nodo inmóvil habilita Turing a D igual | ❌ Korvasová 2015, Nesterenko 2017, Madzvamuse 2015 |
| Modelo de 3 ecuaciones con binding inmóvil | ❌ Nesterenko 2017 (idéntico armazón) |
| Nodo inmóvil debe acoplar al inhibidor | ❌ implícito en el mecanismo de binding |
| λ ∝ √(D absoluta) / longitud de difusión | ❌ libro de texto |
| Robustez a perturbación | ❌ [Ahmad Shaberi et al. 2024](https://doi.org/10.1101/2024.10.15.618426) ya estudia robustez de Turing |
| **Aplicación al himenóforo fúngico** | ✅ **sí (nicho abierto)** |
| **Crítica: el ratio extremo del modelo previo es innecesario** | ✅ **sí (aporte al modelado fúngico)** |
| **Herramienta reproducible (simulador + análisis)** | ✅ **sí (software)** |

### Caveat crítico (estabilidad)
El teorema de Marciniak-Czochra 2016 (no hay patrones estables) aplica si el nodo inmóvil es
**autocatalítico**. La topología base de este trabajo tiene `w_self=0` (no autocatalítica) → cae
en la clase de **binding estable** (Korvasová/Nesterenko), no en la no-estable. **La confirmación
no-lineal no queda invalidada, pero requiere un chequeo de estabilidad a tiempo largo** (que los
spots no sean transitorios que decaigan a spikes).

---

## HILO 2 — Datos para parametrizar/discriminar (parcial)

| Variable | Datos | Fuente |
|---|---|---|
| Inhibidor = Ca²⁺ | Propagación **41–50 µm/s, alcance ~1500 µm** en red micelial | [Itani, Masuo, Yamamoto & Serizawa 2023, *PNAS Nexus*](https://doi.org/10.1093/pnasnexus/pgad012) |
| Activador = hidrofobina | Ensamblaje dependiente de concentración (~300 µg/ml), promovido por glucanos; sin D limpio | [de Vocht et al., PMC2785318](https://pmc.ncbi.nlm.nih.gov/articles/PMC2785318/) |
| Microestructura poliporo (φ, τ, quitina/glucano) | Por segmento (Fomes) | Klemm et al. 2024 |
| D absolutos / tasas de reacción | No encontrados listos para usar | — |

→ **Suficiente para HIPOTETIZAR las variables** (v=Ca²⁺ rápido/medido, u=hidrofobina lento,
w=densidad hifal acoplada al inhibidor); **insuficiente para validar/romper la tautología** sin
mediciones nuevas. El dato de Ca²⁺ (velocidad → D efectiva vía Stokes-Einstein) es el ancla
independiente más concreta para empezar.

---

## HILO 3 — A qué otros estudios ayuda

- **Materiales vivos / biomateriales:** generar morfologías de poro (TPMS, micelio-composites)
  antes de cultivar.
- **Familia teórica nodo-inmóvil/bulk-surface:** aporta un sistema biológico real adicional
  (Korvasová, Madzvamuse, [Pelz & Ward 2023, *Phil. Trans. R. Soc. A*](https://doi.org/10.1098/rsta.2022.0089)).
- **Otros patrones periódicos:** vellosidades intestinales, plumas — marco "longitud de difusión,
  múltiples mecanismos" transferible.
- **Computación con micelio:** las medidas de Ca²⁺ conectan con reservoir computing fúngico.

---

## Citabilidad — veredicto

| Como… | ¿Citable? |
|---|---|
| Software paper (JOSS / SoftwareX) | **Sí, fuerte** — independiente de la novedad del mecanismo |
| Nota aplicada/crítica fúngica (Theory in Biosciences / Fungal Biology / J. Theor. Biol.) | **Sí, condicionado** — citando Nesterenko/Korvasová/Marciniak-Czochra como fundamento; aporta caso fúngico + crítica del ratio extremo + herramienta |
| Paper teórico de math-bio (mecanismo nuevo) | **No** — re-deriva 2013–2017 |
| Paper empírico que explique la morfogénesis fúngica | **No** — tautología + datos bloqueados; además Nesterenko ya validó con datos en otros sistemas |

---

## Caminos constructivos (no es trabajo perdido)

1. **Tool paper YA** — `FungalMorphoSpace` + la suite de análisis (dispersión, búsqueda de
   topología, confirmación). Citable sin depender de la novedad del mecanismo.
2. **Conseguir las mediciones discriminantes** (D de Ca²⁺ vía Itani 2023; hidrofobina) →
   parametrizar y validar como hizo Nesterenko, pero para hongos. *Eso* sí sería un paper
   aplicado fuerte (no solo crítica).
3. **Encuadrar P2 honestamente** sobre el fundamento publicado, no contra él.

---

## Limitaciones de esta revisión + divulgación
- Búsqueda web en inglés; no exhaustiva en Scopus/WoS. "No encontré follow-up de Kuhar" es
  indicativo, no prueba de inexistencia.
- Referencias verificadas en Crossref; sin vibe-citing.
- Asistido por herramienta de IA (deep-research, lit-review).
