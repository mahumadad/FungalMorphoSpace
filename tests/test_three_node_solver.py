"""Tests for the 3-node reaction-diffusion solver (P2 reformulation).

See docs/superpowers/specs/2026-06-02-three-node-reformulation-design.md
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


class _NullKinetics:
    """Trivial 3-node kinetics with zero reactions (for pure-diffusion tests)."""

    def f(self, u, v, w):
        return np.zeros_like(u)

    def g(self, u, v, w):
        return np.zeros_like(v)

    def h(self, u, v, w):
        return np.zeros_like(w)


def test_immobile_node_does_not_diffuse():
    """The defining property: w (hyphal density) has D_w=0, so under zero
    reaction it stays spatially fixed while a diffusible field (u) smooths."""
    from fungalmorphospace.core.three_node import ThreeNodeSimulator

    n = 64
    sim = ThreeNodeSimulator(
        kinetics=_NullKinetics(), D_u=0.2, D_v=0.2, grid_size=n, dx=1.0, dt=0.1, seed=1
    )

    # Seed a localized bump in both u and the immobile node w.
    bump = np.zeros((n, n))
    bump[n // 2, n // 2] = 1.0
    sim.u = bump.copy()
    sim.v = np.zeros((n, n))
    sim.w = bump.copy()

    w0 = sim.w.copy()
    for _ in range(50):
        sim.step()

    # u must have diffused (peak spread out -> max strictly decreased).
    assert sim.u.max() < 1.0 - 1e-6, "diffusible activator u failed to diffuse"
    # w must be unchanged (no diffusion, no reaction).
    assert np.allclose(sim.w, w0), "immobile node w must not change under pure diffusion"


class _LinearDecayKinetics:
    """Linear kinetics f = -a*u (g=h=0), for exact dispersion-relation checks."""

    def __init__(self, a: float):
        self.a = float(a)

    def f(self, u, v, w):
        return -self.a * u

    def g(self, u, v, w):
        return np.zeros_like(v)

    def h(self, u, v, w):
        return np.zeros_like(w)


def test_fourier_mode_matches_discrete_dispersion_relation():
    """A single Fourier mode must grow/decay at exactly the rate predicted by
    the discrete linear dispersion relation: for u_t = D_u lap(u) - a u, a mode
    with grid phase phi = 2*pi*k/n has discrete Laplacian eigenvalue
    mu = -2(1-cos phi)/dx^2, so the exact per-step multiplier is
    m = 1 + dt*(D_u*mu - a), and after N steps the modal amplitude is m^N."""
    from fungalmorphospace.core.three_node import ThreeNodeSimulator

    n, dx, dt = 64, 1.0, 0.05
    D_u, a, k, N = 0.3, 0.1, 4, 40

    sim = ThreeNodeSimulator(
        kinetics=_LinearDecayKinetics(a), D_u=D_u, D_v=D_u, grid_size=n, dx=dx, dt=dt, seed=2
    )

    # Pure cosine mode in x: u(x,y) = cos(2*pi*k*x/n).
    x = np.arange(n)
    mode = np.cos(2.0 * np.pi * k * x / n)[None, :] * np.ones((n, 1))
    sim.u = mode.copy()
    sim.v = np.zeros((n, n))
    sim.w = np.zeros((n, n))

    # Theoretical exact discrete multiplier per step.
    phi = 2.0 * np.pi * k / n
    mu = -2.0 * (1.0 - np.cos(phi)) / dx ** 2
    m = 1.0 + dt * (D_u * mu - a)
    expected_amp = m ** N

    norm0 = float((sim.u * mode).sum())
    for _ in range(N):
        sim.step()
    amp = float((sim.u * mode).sum()) / norm0

    assert np.isclose(amp, expected_amp, rtol=1e-9, atol=1e-12), (
        f"modal amplitude {amp:.6g} != dispersion prediction {expected_amp:.6g}"
    )


def test_three_node_gm_reduces_to_gierer_meinhardt_when_kappa_zero():
    """ThreeNodeGM with kappa=0 (immobile node decoupled from u's inhibition)
    must reproduce the validated 2-field Gierer-Meinhardt activator term."""
    from fungalmorphospace.core.kinetics import GiererMeinhardtKinetics
    from fungalmorphospace.core.three_node import ThreeNodeGM

    rng = np.random.default_rng(3)
    u = rng.uniform(0.1, 1.0, (16, 16))
    v = rng.uniform(0.1, 1.0, (16, 16))
    w = rng.uniform(0.1, 1.0, (16, 16))

    gm = GiererMeinhardtKinetics(rho=0.2, a=0.1, b=1.0)
    tn = ThreeNodeGM(rho=0.2, a=0.1, b=1.0, kappa=0.0, tau_w=10.0)

    assert np.allclose(tn.f(u, v, w), gm.f(u, v)), "kappa=0 must reduce to GM activator"
    assert np.allclose(tn.g(u, v, w), gm.g(u, v)), "kappa=0 must reduce to GM inhibitor"


def test_three_node_gm_immobile_node_relaxes_toward_target():
    """The immobile node w relaxes toward its activator-driven target sigma(u):
    h>0 where w is below target, h<0 where w is above it."""
    from fungalmorphospace.core.three_node import ThreeNodeGM

    tn = ThreeNodeGM(rho=0.2, a=0.1, b=1.0, kappa=0.5, tau_w=10.0)
    u = np.full((4, 4), 1.0)
    v = np.full((4, 4), 1.0)
    target = u ** 2  # sigma(u) = u^2

    below = np.full((4, 4), 0.1)  # w below target (1.0)
    above = np.full((4, 4), 2.0)  # w above target (1.0)

    assert np.all(tn.h(u, v, below) > 0), "w below target must increase"
    assert np.all(tn.h(u, v, above) < 0), "w above target must decrease"
    # At target, relaxation is zero.
    assert np.allclose(tn.h(u, v, target), 0.0)


def test_spacing_from_spots_tracks_pattern_scale():
    """spacing_from_spots(grid, spots) = grid/sqrt(spots) is a robust pattern-
    scale proxy (independent of FFT low-frequency binning): more spots in the
    same domain => smaller spacing, monotonically."""
    from fungalmorphospace.core.three_node import spacing_from_spots

    grid = 256
    # A 4x4 array of spots (16) has spacing ~ grid/4 = 64.
    assert np.isclose(spacing_from_spots(grid, 16), 64.0)
    # More spots -> smaller spacing; strictly decreasing.
    assert spacing_from_spots(grid, 4) > spacing_from_spots(grid, 16) > spacing_from_spots(grid, 64)
    # Degenerate / no pattern -> not a number (cannot define a spacing).
    assert np.isnan(spacing_from_spots(grid, 0))
    assert np.isnan(spacing_from_spots(grid, 1))


def test_marcon_network_bounded_under_strong_instability():
    """The nonlinear kinetics must bound ALL three fields, not just u. A strongly
    Turing-unstable M would let a purely-linear v/w grow without limit; the
    saturation on every node must keep the simulated fields finite and bounded."""
    from fungalmorphospace.core.three_node import MarconNetwork, ThreeNodeSimulator

    M = np.array([[1.845, -1.771, 0.0],
                  [0.323, -0.887, 3.249],
                  [0.000, 3.384, 0.0]])  # strongly unstable (max growth ~2.6)
    kin = MarconNetwork(M=M, steady_state=(1.0, 1.0, 1.0), gamma=0.5)
    sim = ThreeNodeSimulator(kin, D_u=5.0, D_v=5.0, grid_size=32, dx=1.0, dt=0.01, seed=1)
    rng = np.random.default_rng(1)
    sim.u = 1.0 + 0.1 * rng.standard_normal((32, 32))
    sim.v = 1.0 + 0.1 * rng.standard_normal((32, 32))
    sim.w = 1.0 + 0.1 * rng.standard_normal((32, 32))

    for _ in range(3000):
        sim.step()

    for field in (sim.u, sim.v, sim.w):
        assert np.all(np.isfinite(field)), "fields must stay finite"
        assert np.abs(field).max() < 100.0, "fields must stay bounded (saturated), not run away"


def test_dispersion_growth_rate_matches_diagonal_eigenvalues():
    """For a diagonal (decoupled) system the growth rate at wavenumber k is the
    max of the diagonal eigenvalues shifted by -k^2 * D_i. Exact, hand-checkable."""
    from fungalmorphospace.core.three_node import dispersion_growth_rate

    M = np.diag([-1.0, -2.0, -3.0])
    D = [0.1, 0.1, 0.0]
    # At k=2: eigenvalues -1-0.4, -2-0.4, -3 -> max = -1.4 (w has no diffusion).
    assert np.isclose(dispersion_growth_rate(M, D, k=2.0), -1.4)
    # At k=0: max diagonal = -1.0.
    assert np.isclose(dispersion_growth_rate(M, D, k=0.0), -1.0)


def test_fastest_growing_wavelength_nan_when_no_turing_instability():
    """A system stable at every k>0 (no diffusion-driven instability) has no
    pattern-forming mode, so the fastest-growing wavelength is undefined."""
    from fungalmorphospace.core.three_node import fastest_growing_wavelength

    M = np.diag([-1.0, -2.0, -3.0])  # always stable; growth only decreases with k
    D = [0.1, 0.1, 0.0]
    assert np.isnan(fastest_growing_wavelength(M, D))


def test_fastest_growing_wavelength_detects_turing_band():
    """A genuine Turing case (stable at k=0, a positive growth band at finite k)
    returns a finite, positive wavelength at the peak of the dispersion curve."""
    from fungalmorphospace.core.three_node import (
        fastest_growing_wavelength,
        dispersion_growth_rate,
    )

    # Activator-inhibitor with unequal diffusion: classic Turing.
    # J = [[1, -1.5], [2, -2.2]] (trace<0, det>0 -> stable homogeneous),
    # embedded as a 3x3 with a decoupled stable immobile node.
    M = np.array([[1.0, -1.5, 0.0],
                  [2.0, -2.2, 0.0],
                  [0.0, 0.0, -1.0]])
    D = [0.05, 1.0, 0.0]  # inhibitor diffuses faster -> Turing band exists

    assert np.isclose(dispersion_growth_rate(M, D, k=0.0), max(np.linalg.eigvals(M).real))
    lam = fastest_growing_wavelength(M, D)
    assert np.isfinite(lam) and lam > 0


def test_fastest_growing_wavelength_invariant_to_initial_kmax():
    """The reported wavelength must not depend on the initial k-scan bound: if
    the true peak lies beyond k_max the function must widen the scan to find it,
    not report the bogus edge value 2*pi/k_max."""
    from fungalmorphospace.core.three_node import fastest_growing_wavelength

    M = np.array([[1.0, -1.5, 0.0],
                  [2.0, -2.2, 0.0],
                  [0.0, 0.0, -1.0]])
    D = [0.05, 1.0, 0.0]

    lam_ref = fastest_growing_wavelength(M, D, k_max=20.0)  # wide enough: ground truth
    assert np.isfinite(lam_ref)
    # Auto-widening recovers the true peak from a reasonable initial bound below k*.
    assert np.isclose(fastest_growing_wavelength(M, D, k_max=1.0), lam_ref, rtol=0.05)
    # A pathologically narrow scan must never fabricate the bogus edge value 2*pi/k_max:
    # it returns NaN (honest "undetermined") or the correct value, never the artifact.
    lam_tiny = fastest_growing_wavelength(M, D, k_max=0.3)
    assert np.isnan(lam_tiny) or np.isclose(lam_tiny, lam_ref, rtol=0.05)
    assert not np.isclose(lam_tiny, 2.0 * np.pi / 0.3, rtol=0.05)


def test_marcon_network_reactions_vanish_at_steady_state():
    """A linear 3-node interaction network must have zero net reaction at its
    declared homogeneous steady state (so it is a genuine fixed point about
    which Turing instability can be analyzed)."""
    from fungalmorphospace.core.three_node import MarconNetwork

    M = np.array([
        [0.5, -1.0, -0.8],   # u: self-activation, inhibited by v and w
        [1.0, -0.6, 0.0],    # v: driven by u, self-decay
        [1.0, 0.0, -0.4],    # w: driven by u, self-decay (immobile node)
    ])
    ss = (1.0, 0.5, 0.8)
    net = MarconNetwork(M=M, steady_state=ss, gamma=0.3)

    u = np.full((8, 8), ss[0])
    v = np.full((8, 8), ss[1])
    w = np.full((8, 8), ss[2])

    assert np.allclose(net.f(u, v, w), 0.0), "f must vanish at steady state"
    assert np.allclose(net.g(u, v, w), 0.0), "g must vanish at steady state"
    assert np.allclose(net.h(u, v, w), 0.0), "h must vanish at steady state"
