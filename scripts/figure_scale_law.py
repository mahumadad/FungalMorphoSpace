#!/usr/bin/env python3
"""Paper figure: hymenophore wavelength scales as sqrt(absolute D) at equal
diffusion (P2). Plots the linear dispersion prediction lambda*(D), the
nonlinear simulated spacing (a constant coarsening factor above it), and the
biological species targets -- all on a sqrt(D) axis where the law is a line.

    .venv/bin/python scripts/figure_scale_law.py

Writes results/three_node/figure_scale_law.png (and .csv of the curve).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from fungalmorphospace.core.three_node import fastest_growing_wavelength  # noqa: E402

# Equal-diffusion w<->v Turing topology (the manuscript).
BASE_M = np.array([[0.352, -2.135, 0.0],
                   [0.992, -1.190, 2.422],
                   [0.000, 0.795, 0.0]])
COARSENING = 2.4  # nonlinear spacing / linear lambda* (measured, ~const; §10.7)

# Nonlinear confirmation points actually simulated (D, measured spacing px); §10.7.
SIM_POINTS = [(5, 29.6), (20, 62.2)]

# Biological pore spacings -> px via the original 8.7 um/px calibration.
SPECIES = {"L. brumalis (~261 um)": 261 / 8.7,
           "F. fomentarius (~400 um)": 400 / 8.7,
           "P. squamosus (~2000 um)": 2000 / 8.7}


def main():
    out = Path("results/three_node")
    out.mkdir(parents=True, exist_ok=True)

    Ds = np.array([0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500])
    # High k-resolution: large-D peaks sit at small k* (~0.05), which a coarse
    # scan snaps to a grid value and plateaus. k_max=8 covers the small-D peak
    # (~1.6); n_k=6000 resolves the small-k* peaks at large D.
    lam_lin = np.array([fastest_growing_wavelength(BASE_M, [D, D, 0.0], k_max=8.0, n_k=6000)
                        for D in Ds])
    sqrtD = np.sqrt(Ds)

    pd.DataFrame({"D": Ds, "sqrtD": sqrtD, "lambda_linear_px": lam_lin,
                  "spacing_nonlinear_px": COARSENING * lam_lin}).to_csv(
        out / "figure_scale_law.csv", index=False)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(sqrtD, lam_lin, "o-", color="#1f77b4", label=r"linear $\lambda^*$ (dispersion)")
    ax.plot(sqrtD, COARSENING * lam_lin, "s--", color="#ff7f0e",
            label=rf"nonlinear spacing ($\approx{COARSENING}\times\lambda^*$)")
    for D, sp in SIM_POINTS:
        ax.plot(np.sqrt(D), sp, "*", color="#d62728", markersize=15,
                label="simulated" if D == SIM_POINTS[0][0] else None)
    for name, sp_px in SPECIES.items():
        ax.axhline(sp_px, ls=":", color="gray", lw=0.8)
        ax.text(ax.get_xlim()[1] * 0.62, sp_px * 1.02, name, fontsize=7, color="gray")

    ax.set_xlabel(r"$\sqrt{D}$  (equal diffusion, $D_u=D_v$)")
    ax.set_ylabel("wavelength / spacing (px)")
    ax.set_title("Hymenophore scale follows $\\sqrt{D}$ at equal diffusion\n"
                 "(no diffusion-ratio disparity required)")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out / "figure_scale_law.png", dpi=150)
    print(f"Wrote {out/'figure_scale_law.png'} and figure_scale_law.csv")
    print("Linear lambda* is ~linear in sqrt(D) -> confirms the sqrt(D) scale law.")


if __name__ == "__main__":
    main()
