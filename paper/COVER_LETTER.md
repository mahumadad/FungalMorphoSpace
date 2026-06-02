# Cover letter — SoftwareX submission

Mario Ahumada Durán
Independent Researcher
mahumadad@gmail.com

To the Editors of *SoftwareX*

Dear Editors,

I am submitting the original software article **"FungalMorphoSpace: A
reaction-diffusion simulator for exploring fungal hymenophore morphogenesis"**
for consideration in *SoftwareX*.

**What the software does and why it is needed.** FungalMorphoSpace is an
open-source Python package that simulates Gierer–Meinhardt reaction–diffusion
patterns and quantifies them with a morphometric and topological analysis suite
(FFT and autocorrelation wavelength estimation with quality control, Euler
characteristic, connected-component spacing statistics). It ships calibrated
parameter presets for three polypore species and a contract-enforced output
schema for reproducibility. It also includes a three-node (immobile-node) module
with an analytical linear-stability utility. To my knowledge no existing tool
provides an integrated, reproducible pipeline tailored to fungal hymenophore
morphogenesis — combining simulation, objective morphometry, and a fixed output
contract — which is what makes systematic morphospace exploration practical.

**Scope of the contribution (stated plainly).** This is a *tools* contribution,
not a claim of a new pattern-forming mechanism or of an *ab initio* predictor of
fungal morphology. The equal-diffusion / immobile-node mechanism the software can
explore is established in the literature (Korvasová et al. 2015; Nesterenko et
al. 2017; Marcon et al. 2016; Raspopovic et al. 2014; Marciniak-Czochra et al.
2016), and the manuscript cites it as foundation rather than as novelty. The
parameter presets are calibrated empirically, which the manuscript discloses as
its central limitation: the tool establishes *consistency* (observed pore
spacings are reproducible within the model), not that fungi employ this
mechanism. Its value is as a reproducible instrument for hypothesis generation
and morphospace exploration, and it makes explicit — analytically and in
simulation — that the large diffusion ratio used in earlier modelling is one
parameterization of a diffusion length rather than a necessary ingredient.

**Software quality and reproducibility.** The package is modular and documented,
ships with an automated test suite, and enforces a machine-readable output schema
for cross-version reproducibility. All figures and results in the manuscript are
reproducible from the included scripts. The code is openly available under a
CC BY-NC 4.0 license at https://github.com/mahumadad/FungalMorphoSpace and a
versioned release is archived with a DOI (see the Code metadata table).

**Declarations.** This manuscript is original, has not been published previously,
and is not under consideration elsewhere. The author declares no competing
interests. An AI assistant was used for code review, literature search and
verification, and drafting; all outputs, including every cited reference, were
verified by the author, who takes full responsibility for the content.

Thank you for considering this submission.

Sincerely,
Mario Ahumada Durán
