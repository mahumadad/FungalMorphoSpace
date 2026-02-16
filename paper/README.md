# Paper: FungalMorphoSpace (SoftwareX)

## Compilation

```bash
cd paper/
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

Or with latexmk:

```bash
latexmk -pdf main.tex
```

## Figures

Place generated figures in `figures/`. To generate them:

```bash
# From project root:
python scripts/run_integrated_validation.py --species all --n_runs 5 --grid 1024
python scripts/plot_scaling_law.py --results-dir results/<timestamp>/
python scripts/test_laminillas.py
```

## TODOs before submission

- [ ] Replace author email placeholder
- [ ] Verify institution/affiliation
- [ ] Add acknowledgements
- [ ] Generate and include figures
- [ ] Final word count check (~2500-3000 for SoftwareX)
- [ ] Proofread
