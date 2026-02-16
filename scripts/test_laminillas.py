#!/usr/bin/env python3
"""TOPOLOGICAL TRANSITION EXPERIMENT: PORES -> LAMELLAE -> LABYRINTHS
======================================================================

This script demonstrates Kuhar's LALIP hypothesis on **topological continuity**:
pored, gilled, and labyrinthine hymenophores can emerge from the *same* reaction–diffusion
system under **anisotropic diffusion** (a proxy for directional mechanical stress during growth).

Kinetics: Gierer–Meinhardt (aligned with the package core)
----------------------------------------------------------
The core implementation in `fungalmorphospace.core.kinetics.GiererMeinhardtKinetics` uses:

    du/dt = D_u ∇²u + ρ( a - b u + u²/(v + ε) )
    dv/dt = D_v ∇²v + ρ( u² - v )

with ε = 1e-6 for numerical robustness.

Anisotropy is introduced in the Laplacian operator:

    ∇²_aniso = α_y * ∂²/∂y² + α_x * ∂²/∂x²

Where α_y > α_x creates vertical elongation (lamellae-like bands).

Author: Mario Ahumada Durán
Date: January 2026
Reference: Kuhar et al. (2022) Theory in Biosciences. DOI: 10.1007/s12064-022-00363-z
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Use the same kinetics as the main simulator to avoid divergence.
import sys

# Ensure local package import works when running from the repo
REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_PATH = REPO_ROOT / "src"
if str(SRC_PATH) not in sys.path:
    sys.path.insert(0, str(SRC_PATH))

from fungalmorphospace.core.kinetics import GiererMeinhardtKinetics


# =============================================================================
# BASE PARAMETERS (Fomes fomentarius - Validated)
# =============================================================================

PARAMS_BASE = {
    "D_v_D_u": 150.0,
    "b": 1.0,
    "rho": 0.2,
    "a": 0.1,
    "grid": 512,
    "T_target": 5.0,   # Physical time
    "dt": 0.0005,      # Requested dt (may be reduced by CFL condition below)
    "dx": 1.0,
}


def _laplacian_aniso(Z: np.ndarray, alpha_y: float, alpha_x: float) -> np.ndarray:
    """Anisotropic Laplacian with periodic boundary conditions.

    ∇²_aniso = α_y * ∂²/∂y² + α_x * ∂²/∂x²
    """
    dy = (np.roll(Z, 1, axis=0) + np.roll(Z, -1, axis=0) - 2 * Z) * alpha_y
    dx = (np.roll(Z, 1, axis=1) + np.roll(Z, -1, axis=1) - 2 * Z) * alpha_x
    return dy + dx


def _cfl_dt_safe(D_u: float, D_v: float, dx: float, alpha_y: float, alpha_x: float, *, safety: float = 0.8) -> float:
    """Conservative CFL dt bound for explicit diffusion updates (2D)."""
    # Conservative: scale by the *largest* anisotropy factor.
    D_max = max(D_u, D_v) * max(alpha_y, alpha_x)
    dt_max = (dx ** 2) / (4.0 * D_max)
    return safety * dt_max


def run_anisotropic_simulation(mode_name: str, *, anisotropy_y: float = 1.0, anisotropy_x: float = 1.0) -> np.ndarray:
    """Simulate GM kinetics with anisotropic diffusion.

    Notes
    -----
    - αy = αx = 1.0: isotropic -> spots (pores)
    - αy > αx: vertical elongation -> lamellae-like bands
    - intermediate: labyrinth patterns
    """
    print(f"\nGenerating morphotype: {mode_name}")
    print(f"  Anisotropy  Y={anisotropy_y}  X={anisotropy_x}")

    size = int(PARAMS_BASE["grid"])
    dx = float(PARAMS_BASE["dx"])

    # Initialize around the *core* steady state (consistent with kinetics in src/)
    rho = float(PARAMS_BASE["rho"])
    a = float(PARAMS_BASE["a"])
    b = float(PARAMS_BASE["b"])

    kinetics = GiererMeinhardtKinetics(rho=rho, a=a, b=b)
    u0, v0 = kinetics.get_steady_state()

    np.random.seed(42)
    U = u0 + np.random.normal(0.0, 0.01, (size, size))
    V = v0 + np.random.normal(0.0, 0.01, (size, size))
    U = np.maximum(U, 0.0)
    V = np.maximum(V, 1e-6)  # prevent division issues

    # Diffusion
    D_u = 1.0
    D_v = D_u * float(PARAMS_BASE["D_v_D_u"])

    # Time discretization with conservative CFL guard (anisotropy-aware)
    dt_req = float(PARAMS_BASE["dt"])
    dt_safe = _cfl_dt_safe(D_u, D_v, dx, anisotropy_y, anisotropy_x, safety=0.8)
    dt = min(dt_req, dt_safe)
    if dt < dt_req:
        print(f"  dt reduced by CFL: {dt_req:.6g} -> {dt:.6g} (safe={dt_safe:.6g})")

    T_target = float(PARAMS_BASE["T_target"])
    steps = int(np.ceil(T_target / dt))

    for i in range(steps):
        lap_U = _laplacian_aniso(U, anisotropy_y, anisotropy_x)
        lap_V = _laplacian_aniso(V, anisotropy_y, anisotropy_x)

        # Core-aligned kinetics (vectorized)
        f_uv = kinetics.f(U, V)
        g_uv = kinetics.g(U, V)

        U += (D_u * lap_U + f_uv) * dt
        V += (D_v * lap_V + g_uv) * dt

        U = np.maximum(U, 0.0)
        V = np.maximum(V, 1e-6)

        if i % max(1, steps // 10) == 0:
            print(f"  progress: {100 * i / steps:>3.0f}%", end="\r")

    print("  progress: 100%")

    return U


if __name__ == "__main__":
    # 1) Pore state (isotropic)
    res_poro = run_anisotropic_simulation(
        "1. Pores (Fomes)",
        anisotropy_y=1.0,
        anisotropy_x=1.0,
    )

    # 2) Lamella state (strong vertical anisotropy)
    res_lamina = run_anisotropic_simulation(
        "2. Lamellae (Lenzites)",
        anisotropy_y=4.0,
        anisotropy_x=0.5,
    )

    # 3) Labyrinth state (intermediate anisotropy)
    res_maze = run_anisotropic_simulation(
        "3. Labyrinth (Daedalea)",
        anisotropy_y=2.0,
        anisotropy_x=0.7,
    )

    # Visualization
    print("\nGenerating comparative panel...")
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    titles = [
        "State 1: PORES (Isotropic)\n(Fomes fomentarius)\nGM: Hexagonal spots",
        "State 2: LAMELLAE (Anisotropic)\n(Lenzites betulina)\nDirectional stress proxy",
        "State 3: LABYRINTH (Hybrid)\n(Daedalea quercina)\nSymmetry-breaking transition",
    ]
    patterns = [res_poro, res_lamina, res_maze]

    for ax, pat, title in zip(axes, patterns, titles):
        ax.imshow(pat, cmap="magma", origin="lower")
        ax.set_title(title, fontsize=11)
        ax.axis("off")

    plt.suptitle(
        "Kuhar LALIP Hypothesis: Topological Continuity in Fungal Hymenophores\n"
        "(Gierer–Meinhardt kinetics with anisotropic diffusion)",
        fontsize=13,
        fontweight="bold",
    )

    repo_root = Path(__file__).resolve().parents[1]
    output_dir = repo_root / "results" / "topology_experiment"
    output_dir.mkdir(parents=True, exist_ok=True)
    save_path = output_dir / "continuum_topology_kuhar.png"

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {save_path}")
    plt.show()
