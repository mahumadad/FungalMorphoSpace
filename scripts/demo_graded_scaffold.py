#!/usr/bin/env python3
"""DEMO (exploratory): a graded-density Turing pattern as a 2D "scaffold".

Idea (user): seed the pattern centrally and grade it so the centre is denser,
like shark denticles (Turing + initiator + positional bias). The sqrt(D) scale
law (the manuscript) gives the recipe: spacing ~ sqrt(D), so a
radial diffusion gradient with SMALL D at the centre yields a FINER (denser)
pattern there and a COARSER one at the rim -- a functionally graded porosity.

This is a standalone illustrative demo (not library code, not a claim of
novelty). It uses a Gierer-Meinhardt activator-inhibitor with a spatially
varying inhibitor diffusivity D_v(x,y) and a central seed, no-flux (Neumann)
boundaries, conservative div(D grad u) discretization.

    .venv/bin/python scripts/demo_graded_scaffold.py
-> results/three_node/demo_graded_scaffold.png
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

N = 256                      # grid
dx = 1.0
steps = 30000
rho, a, b = 0.2, 0.1, 1.0    # GM kinetics
D_u = 0.5                    # activator diffusivity (constant; larger -> cleaner spots)
D_v_min, D_v_max = 8.0, 60.0 # inhibitor diffusivity: small centre -> dense centre
                             # lambda ~ 2*pi*sqrt(D_v/D_u): ~25 px centre, ~69 px rim
                             # both resolved on a 256 grid (no under-resolution)


def neumann_lap_const(f):
    """Laplacian with no-flux (reflective) BC, constant coefficient."""
    fp = np.pad(f, 1, mode="edge")
    return (fp[:-2, 1:-1] + fp[2:, 1:-1] + fp[1:-1, :-2] + fp[1:-1, 2:] - 4 * f) / dx ** 2


def div_D_grad(u, Dface_x, Dface_y):
    """Conservative div(D grad u) with no-flux BC. Dface_* are face diffusivities."""
    up = np.pad(u, 1, mode="edge")
    # fluxes across x-faces (between i and i+1): D_{i+1/2} (u_{i+1}-u_i)
    flux_x = Dface_x * (up[1:-1, 2:] - up[1:-1, 1:-1])   # right faces, shape (N,N)
    flux_xl = Dface_x * (up[1:-1, 1:-1] - up[1:-1, :-2])  # left faces
    flux_y = Dface_y * (up[2:, 1:-1] - up[1:-1, 1:-1])
    flux_yl = Dface_y * (up[1:-1, 1:-1] - up[:-2, 1:-1])
    return ((flux_x - flux_xl) + (flux_y - flux_yl)) / dx ** 2


def main():
    out = Path("results/three_node")
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)

    # radial diffusivity field: small at centre (dense) -> large at rim (coarse)
    yy, xx = np.mgrid[0:N, 0:N]
    cx = cy = (N - 1) / 2
    r = np.sqrt((xx - cx) ** 2 + (yy - cy) ** 2)
    r_norm = np.clip(r / r.max(), 0, 1)
    D_v = D_v_min + (D_v_max - D_v_min) * r_norm ** 2  # smooth radial grade

    # face diffusivities (simple neighbour average; clamp by edge padding)
    Dp = np.pad(D_v, 1, mode="edge")
    Dface_x = 0.5 * (Dp[1:-1, 1:-1] + Dp[1:-1, 2:])
    Dface_y = 0.5 * (Dp[1:-1, 1:-1] + Dp[2:, 1:-1])

    # CFL from the largest diffusivity
    dt = 0.8 * dx ** 2 / (4.0 * max(D_u, D_v_max))

    # initial state: homogeneous steady state + central SEED (positional bias)
    u0, v0 = (a + 1) / b, ((a + 1) / b) ** 2
    u = u0 + 0.01 * rng.standard_normal((N, N))
    v = v0 + 0.01 * rng.standard_normal((N, N))
    seed = np.exp(-(r ** 2) / (2 * (0.06 * N) ** 2))   # central Gaussian seed
    u += 0.5 * seed

    eps = 1e-6
    for _ in range(steps):
        lap_u = D_u * neumann_lap_const(u)
        diff_v = div_D_grad(v, Dface_x, Dface_y)
        f = rho * (a - b * u + u ** 2 / (v + eps))
        g = rho * (u ** 2 - v)
        u = np.maximum(u + dt * (lap_u + f), 0.0)
        v = np.maximum(v + dt * (diff_v + g), 1e-6)
        if not np.all(np.isfinite(u)):
            print("blow-up; reduce dt"); break

    # --- spot detection: density proxy = spots per unit area (finer = denser) ---
    from scipy import ndimage
    binu = u > (u.mean() + 0.4 * u.std())
    lbl, nspots = ndimage.label(binu)
    cen = np.array(ndimage.center_of_mass(binu, lbl, range(1, nspots + 1))) if nspots else np.empty((0, 2))
    rc = np.sqrt((cen[:, 1] - cx) ** 2 + (cen[:, 0] - cy) ** 2) / r.max() if len(cen) else np.array([])

    nb = 6
    edges = np.linspace(0, 1, nb + 1)
    spacing = []
    for k in range(nb):
        lo, hi = edges[k], edges[k + 1]
        n_in = int(np.sum((rc >= lo) & (rc < hi)))
        area = np.sum((r_norm >= lo) & (r_norm < hi)) * dx ** 2
        spacing.append(np.sqrt(area / n_in) if n_in > 0 else np.nan)  # local spot spacing
    mids = 0.5 * (edges[:-1] + edges[1:])

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    im0 = ax[0].imshow(D_v, cmap="viridis"); ax[0].set_title("Inhibitor diffusivity $D_v(x,y)$\n(small centre = dense)")
    plt.colorbar(im0, ax=ax[0], fraction=0.046)
    ax[1].imshow(u, cmap="bone_r"); ax[1].set_title(f"Graded scaffold (activator $u$)\nseed + radial $D$ grade, {nspots} spots")
    ax[2].plot(mids, spacing, "o-")
    ax[2].set_xlabel("normalized radius"); ax[2].set_ylabel("local spot spacing (px)")
    ax[2].set_title("Spacing grows outward\n(denser/finer at centre)"); ax[2].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out / "demo_graded_scaffold.png", dpi=130)
    sp = [s for s in spacing if np.isfinite(s)]
    print(f"Wrote {out/'demo_graded_scaffold.png'}  ({nspots} spots)")
    if len(sp) >= 2:
        print(f"centre spacing {sp[0]:.1f}px vs rim {sp[-1]:.1f}px "
              f"(graded if rim > centre, per sqrt(D) law)")


if __name__ == "__main__":
    main()
