#!/usr/bin/env python3
"""Explore the 3-node reformulation (P2): does a non-diffusible node let the
hymenophore morphospace form without an extreme diffusion ratio?

Runs three experiments with a tqdm progress bar:
  1. PROBE   — ThreeNodeGM: sweep D_v/D_u downward at several immobile-coupling
               strengths (kappa). How low can the ratio go and still pattern?
  2. MARCON  — equal-diffusion 3-node network: sweep a single plausible
               parameter and map the achievable wavelength range.
  3. CONTROL — scrambled-topology network must NOT cover the morphospace.

Run it in your own terminal (NOT inside Claude) to watch progress:

    cd FungalMorphoSpace
    .venv/bin/python scripts/explore_three_node.py            # full sweep
    .venv/bin/python scripts/explore_three_node.py --quick    # fast smoke

Outputs CSVs + FINDINGS.md under results/three_node/.

See docs/superpowers/specs/2026-06-02-three-node-reformulation-design.md.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover - tqdm is a declared dependency
    def tqdm(it, **kw):
        return it

from fungalmorphospace.core.three_node import (  # noqa: E402
    ThreeNodeSimulator,
    ThreeNodeGM,
    MarconNetwork,
    spacing_from_spots,
)
from fungalmorphospace.analysis.topology_analyzer import TopologyAnalyzer  # noqa: E402

MIN_GENUINE_SPOTS = 4  # mirrors the validator's genuine-pattern gate


def run_pattern(kinetics, D_u, D_v, grid, steps, init_mean, seed=0, dt_safety=0.8):
    """Run a 3-node sim to (approximate) steady state; return (lambda_px, spots, qc)."""
    dx = 1.0
    d_max = max(D_u, D_v, 1e-9)
    dt = dt_safety * dx ** 2 / (4.0 * d_max)
    sim = ThreeNodeSimulator(kinetics, D_u=D_u, D_v=D_v, grid_size=grid, dx=dx, dt=dt, seed=seed)
    rng = np.random.default_rng(seed)
    noise = lambda: init_mean + 0.05 * rng.standard_normal((grid, grid))
    sim.u, sim.v, sim.w = noise(), noise(), noise()

    for _ in range(steps):
        sim.step()
        if not np.all(np.isfinite(sim.u)):  # blow-up guard
            return float("nan"), 0, False

    ta = TopologyAnalyzer(sim.u, dx=dx)
    m = ta.compute_all_metrics()
    lam = float(m.get("wavelength_fft", float("nan")))
    spots = int(m.get("n_components", 0))
    qc = bool(m.get("wavelength_fft_qc_pass", False))
    return lam, spots, qc


def grid_for_ratio(ratio, target_wavelengths=6, min_grid=128, max_grid=768):
    """Pick a grid large enough to resolve ~target_wavelengths of the heuristic
    wavelength lambda ~ 2*pi*sqrt(D_v/D_u). High-ratio combos (large lambda)
    therefore get a bigger domain; rounded to a multiple of 64 and clamped to
    [min_grid, max_grid] so compute stays bounded."""
    lam = 2.0 * np.pi * np.sqrt(ratio)
    g = int(np.ceil(target_wavelengths * lam / 64.0) * 64)
    return int(max(min_grid, min(max_grid, g)))


def experiment_probe(base_grid, steps, max_grid=768, adaptive=True):
    """GM + immobile node: sweep D_v/D_u down at several kappa.

    With ``adaptive=True`` each combo uses a per-ratio grid (grid_for_ratio) so
    large-wavelength (high-ratio) runs are not under-resolved; low-ratio runs
    stay at ``base_grid``. With ``adaptive=False`` (e.g. --quick) all runs use
    ``base_grid``.
    """
    ratios = [3750, 1000, 250, 80, 30, 10, 5, 2]
    kappas = [0.0, 0.5, 2.0]
    rows = []
    combos = [(r, k) for k in kappas for r in ratios]
    for ratio, kappa in tqdm(combos, desc="1/3 PROBE (GM+w)", unit="run"):
        g = grid_for_ratio(ratio, min_grid=base_grid, max_grid=max_grid) if adaptive else base_grid
        kin = ThreeNodeGM(rho=0.2, a=0.1, b=1.0, kappa=kappa, tau_w=10.0)
        lam, spots, qc = run_pattern(kin, D_u=0.1, D_v=0.1 * ratio, grid=g, steps=steps, init_mean=1.0)
        wl_in_domain = float(g / lam) if (np.isfinite(lam) and lam > 0) else float("nan")
        rows.append({"ratio": ratio, "kappa": kappa, "grid": g,
                     "lambda_fft_px": lam,
                     "spacing_from_spots_px": round(spacing_from_spots(g, spots), 2),
                     "wavelengths_in_domain": round(wl_in_domain, 2) if np.isfinite(wl_in_domain) else wl_in_domain,
                     "spots": spots, "qc": qc, "genuine": spots >= MIN_GENUINE_SPOTS})
    return pd.DataFrame(rows)


def experiment_marcon(grid, steps):
    """Equal-diffusion 3-node network: sweep one parameter, map lambda."""
    # Base topology: u self-activates, inhibited by v and immobile w; v,w driven by u.
    base = np.array([[0.5, -1.0, -0.8], [1.0, -0.6, 0.0], [1.0, 0.0, -0.4]])
    ss = (1.0, 0.5, 0.8)
    D = 0.5  # EQUAL diffusion for the two diffusible nodes
    strengths = np.linspace(0.2, 2.0, 8)  # single plausible parameter: w-inhibition gain
    rows = []
    for s in tqdm(strengths, desc="2/3 MARCON (D_u=D_v)", unit="run"):
        M = base.copy()
        M[0, 2] = -0.8 * s  # scale the immobile-node inhibition of u
        kin = MarconNetwork(M=M, steady_state=ss, gamma=0.3)
        lam, spots, qc = run_pattern(kin, D_u=D, D_v=D, grid=grid, steps=steps, init_mean=1.0)
        rows.append({"w_inhibition": round(float(0.8 * s), 3),
                     "lambda_fft_px": lam,
                     "spacing_from_spots_px": round(spacing_from_spots(grid, spots), 2),
                     "spots": spots, "qc": qc, "genuine": spots >= MIN_GENUINE_SPOTS})
    return pd.DataFrame(rows)


def experiment_control(grid, steps):
    """Scrambled topology must NOT cover the morphospace."""
    base = np.array([[0.5, -1.0, -0.8], [1.0, -0.6, 0.0], [1.0, 0.0, -0.4]])
    ss = (1.0, 0.5, 0.8)
    D = 0.5
    rng = np.random.default_rng(123)
    rows = []
    for i in tqdm(range(8), desc="3/3 CONTROL (scrambled)", unit="run"):
        M = base.copy()
        # Scramble the off-diagonal interaction signs/magnitudes.
        perm = rng.permutation(M.flatten())
        M = perm.reshape(3, 3)
        kin = MarconNetwork(M=M, steady_state=ss, gamma=0.3)
        lam, spots, qc = run_pattern(kin, D_u=D, D_v=D, grid=grid, steps=steps, init_mean=1.0, seed=i)
        rows.append({"trial": i, "lambda_fft_px": lam,
                     "spacing_from_spots_px": round(spacing_from_spots(grid, spots), 2),
                     "spots": spots, "qc": qc, "genuine": spots >= MIN_GENUINE_SPOTS})
    return pd.DataFrame(rows)


def write_findings(out: Path, probe, marcon, control):
    # PROBE: lowest ratio that still patterns genuinely, per kappa.
    lines = ["# 3-node reformulation — FINDINGS", "",
             "Auto-generated by `scripts/explore_three_node.py`.", ""]
    lines.append("## 1. PROBE (GM + immobile node)")
    if "grid" in probe.columns:
        lines.append(f"- per-ratio adaptive grid used: {sorted(int(g) for g in probe['grid'].unique())} "
                     f"(>=~6 wavelengths in domain; see `wavelengths_in_domain` column)")
    for kappa, grp in probe.groupby("kappa"):
        gen = grp[grp["genuine"]]
        lowest = int(gen["ratio"].min()) if len(gen) else None
        lines.append(f"- kappa={kappa}: lowest D_v/D_u with a genuine pattern = "
                     + (f"{lowest}" if lowest is not None else "none patterned"))
    lines += ["", "## 2. MARCON (equal diffusion D_u=D_v)"]
    gen = marcon[marcon["genuine"]]
    if len(gen):
        lines.append(f"- genuine patterns at EQUAL diffusion: {len(gen)}/{len(marcon)} sweep points")
        # Use the spot-count spacing (robust); the FFT lambda is artifact-prone
        # at low frequency and is NOT a reliable wavelength-control signal.
        sp = gen["spacing_from_spots_px"].dropna()
        if len(sp):
            lines.append(f"- pattern-scale (spacing from spots) range: {sp.min():.1f}–{sp.max():.1f} px "
                         f"as ONE parameter (w-inhibition) varies; spots {int(gen['spots'].min())}–{int(gen['spots'].max())}")
        lines.append("- NOTE: the FFT `lambda_fft_px` column is unreliable here (snaps to a low-frequency "
                     "floor bin, e.g. ~188.6 px); use `spacing_from_spots_px` for scale control.")
    else:
        lines.append("- no genuine pattern formed at equal diffusion (strong claim NOT supported by this topology)")
    lines += ["", "## 3. CONTROL (scrambled topology)"]
    cg = int(control["genuine"].sum())
    lines.append(f"- genuine patterns from scrambled topologies: {cg}/{len(control)} "
                 f"(expected ~0; non-zero weakens the specificity claim)")
    lines += ["", "## Verdict",
              "- If PROBE shows the ratio dropping far below 3750 as kappa grows, the immobile",
              "  node reduces the required disparity (weak claim).",
              "- If MARCON covers a lambda range at D_u=D_v while CONTROL does not, the strong",
              "  (equal-diffusion) claim is supported and the tautology is dissolved.",
              "- Read the CSVs for the full sweep."]
    (out / "FINDINGS.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    ap = argparse.ArgumentParser(description="Explore the 3-node reformulation (P2)")
    ap.add_argument("--grid", type=int, default=256, help="base grid size (default 256)")
    ap.add_argument("--max-grid", type=int, default=768,
                    help="cap for the PROBE's per-ratio adaptive grid (default 768)")
    ap.add_argument("--steps", type=int, default=8000, help="time steps per run (default 8000)")
    ap.add_argument("--quick", action="store_true", help="fast smoke (small fixed grid, few steps)")
    ap.add_argument("--out", type=str, default="results/three_node", help="output dir")
    args = ap.parse_args()

    grid, steps = (64, 1500) if args.quick else (args.grid, args.steps)
    adaptive = not args.quick
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    print(f"3-node exploration | base_grid={grid} max_grid={args.max_grid} steps={steps} "
          f"adaptive_probe={adaptive} | out={out}")
    probe = experiment_probe(grid, steps, max_grid=args.max_grid, adaptive=adaptive)
    probe.to_csv(out / "coverage_probe.csv", index=False)
    marcon = experiment_marcon(grid, steps)
    marcon.to_csv(out / "marcon_sweep.csv", index=False)
    control = experiment_control(grid, steps)
    control.to_csv(out / "control.csv", index=False)

    write_findings(out, probe, marcon, control)
    print(f"\nDone. CSVs + FINDINGS.md written to {out}/")
    print((out / "FINDINGS.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
