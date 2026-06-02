"""Three-node reaction-diffusion solver (P2 reformulation prototype).

Two morphogens (u, v) diffuse; a third node (w, hyphal density) is immobile
(D_w = 0). This is the substrate for testing whether the hymenophore morphospace
can be generated without an extreme diffusion-ratio (Marcon et al., eLife 2016:
Turing patterns with equally diffusing signals plus a non-diffusible node).

See the project manuscript (paper/).

This module is isolated from the validated 2-field engine (turing_simulator.py)
so the tool-paper path stays stable.
"""
from __future__ import annotations

import numpy as np

_EPS = 1e-6  # matches GiererMeinhardtKinetics activator-denominator offset


def dispersion_growth_rate(M, D, k):
    """Linear growth rate of a Fourier mode with wavenumber ``k`` about the
    homogeneous steady state: the largest real part of the eigenvalues of
    ``M - k**2 * diag(D)``, where ``M`` is the reaction Jacobian and ``D`` the
    diffusion coefficients (the immobile node has D=0)."""
    M = np.asarray(M, dtype=float)
    Dk = np.diag(np.asarray(D, dtype=float) * (float(k) ** 2))
    eig = np.linalg.eigvals(M - Dk)
    return float(np.max(eig.real))


def fastest_growing_wavelength(M, D, k_max=3.0, n_k=600, _cap=64.0):
    """Wavelength ``2*pi/k*`` of the fastest-growing finite-k mode IF the system
    is Turing-unstable: stable at k->0 but with a positive growth band peaking at
    a finite, interior k. Returns NaN when there is no diffusion-driven
    instability, or when the peak never settles at a finite k.

    Because the immobile node has D=0, high-k modes are not damped to -inf, so the
    growth curve can keep rising toward a plateau (a fine-scale instability with no
    selected wavelength). The scan auto-widens until the peak is interior; if it is
    still at the boundary at ``_cap`` the wavelength is undefined (NaN) -- this is
    the correct verdict, not the bogus edge value 2*pi/k_max."""
    g0 = dispersion_growth_rate(M, D, 0.0)
    k_hi = float(k_max)
    while True:
        ks = np.linspace(1e-6, k_hi, int(n_k))
        growth = np.array([dispersion_growth_rate(M, D, k) for k in ks])
        i = int(np.argmax(growth))
        if growth[i] <= 0.0 or growth[i] <= g0:
            return float("nan")  # no instability, or maximal growth is at k->0
        if i >= len(ks) - 1:     # peak at the scan edge -> widen and retry
            if k_hi >= _cap:
                return float("nan")  # still rising at the cap: no finite-scale selection
            k_hi *= 2.0
            continue
        return float(2.0 * np.pi / ks[i])


def spacing_from_spots(grid, spots):
    """Robust pattern-scale proxy: the inter-feature spacing of ``spots``
    roughly equidistributed features in a ``grid`` x ``grid`` domain is
    ``grid / sqrt(spots)``. Independent of FFT low-frequency binning (which
    snaps large wavelengths to a single floor bin). Returns NaN when there is
    no genuine periodic arrangement to measure (fewer than 2 features)."""
    if spots is None or spots < 2:
        return float("nan")
    return float(grid) / np.sqrt(float(spots))


class ThreeNodeGM:
    """Probe model: Gierer-Meinhardt augmented with an immobile hyphal node w.

    The immobile node modifies the *local* inhibition seen by the activator
    (``v + kappa*w``) and relaxes slowly toward an activator-driven target
    ``sigma(u) = u**2`` on timescale ``tau_w``. With ``kappa=0`` the activator
    and inhibitor reduce exactly to the validated 2-field Gierer-Meinhardt
    kinetics, so this is a strict generalization used to ask: does the immobile
    node lower the diffusion ratio required for patterning?
    """

    def __init__(self, rho: float = 0.2, a: float = 0.1, b: float = 1.0,
                 kappa: float = 0.5, tau_w: float = 10.0) -> None:
        self.rho = float(rho)
        self.a = float(a)
        self.b = float(b)
        self.kappa = float(kappa)
        self.tau_w = float(tau_w)

    def f(self, u, v, w):
        return self.rho * (self.a - self.b * u + u ** 2 / (v + self.kappa * w + _EPS))

    def g(self, u, v, w):
        return self.rho * (u ** 2 - v)

    def h(self, u, v, w):
        # Slow relaxation of hyphal density toward activator-driven target u^2.
        return (u ** 2 - w) / self.tau_w


class MarconNetwork:
    """Strong model: a 3-node interaction network for Turing patterns with
    *equal* diffusion (D_u = D_v) plus a non-diffusible node w.

    Following Marcon et al. (eLife 2016) and Raspopovic et al. (Science 2014,
    Bmp-Sox9-Wnt), a Turing instability can arise without a diffusion-rate
    disparity when one node is immobile and the interaction topology is right.
    Reactions are linear about a homogeneous steady state ``x*`` (so it is a
    genuine fixed point), with a cubic saturation on the activator ``u`` to
    bound growth without shifting the fixed point. The interaction matrix ``M``
    encodes the network topology; the experiment sweeps it for instability.
    """

    def __init__(self, M, steady_state, gamma: float = 0.3) -> None:
        self.M = np.asarray(M, dtype=float)
        if self.M.shape != (3, 3):
            raise ValueError("M must be a 3x3 interaction matrix")
        self.uss, self.vss, self.wss = (float(s) for s in steady_state)
        self.gamma = float(gamma)

    def _sat(self, d):
        # Cubic self-saturation. Vanishes at the fixed point (d=0) with zero
        # derivative, so neither the steady state nor the Jacobian (J=M) change;
        # it only bounds the nonlinear growth. Clipped to avoid float overflow.
        return self.gamma * np.clip(d, -50.0, 50.0) ** 3

    def f(self, u, v, w):
        du, dv, dw = u - self.uss, v - self.vss, w - self.wss
        return self.M[0, 0] * du + self.M[0, 1] * dv + self.M[0, 2] * dw - self._sat(du)

    def g(self, u, v, w):
        du, dv, dw = u - self.uss, v - self.vss, w - self.wss
        return self.M[1, 0] * du + self.M[1, 1] * dv + self.M[1, 2] * dw - self._sat(dv)

    def h(self, u, v, w):
        du, dv, dw = u - self.uss, v - self.vss, w - self.wss
        return self.M[2, 0] * du + self.M[2, 1] * dv + self.M[2, 2] * dw - self._sat(dw)


class ThreeNodeSimulator:
    """Explicit-Euler solver for a 3-node RD system with one immobile node.

    Fields ``u`` and ``v`` diffuse with coefficients ``D_u``, ``D_v``; the node
    ``w`` does not diffuse (``D_w = 0``). Periodic (toroidal) boundary
    conditions via ``numpy.roll``.
    """

    def __init__(
        self,
        kinetics,
        D_u: float = 0.1,
        D_v: float = 0.1,
        grid_size: int = 256,
        dx: float = 1.0,
        dt: float = 0.1,
        seed: int = 42,
    ) -> None:
        self.kinetics = kinetics
        self.D_u = float(D_u)
        self.D_v = float(D_v)
        self.grid_size = int(grid_size)
        self.dx = float(dx)
        self.dt = float(dt)
        self.rng = np.random.default_rng(seed)

        shape = (self.grid_size, self.grid_size)
        self.u = np.zeros(shape)
        self.v = np.zeros(shape)
        self.w = np.zeros(shape)

    def _laplacian(self, field: np.ndarray) -> np.ndarray:
        """5-point discrete Laplacian with periodic BCs (toroidal)."""
        return (
            np.roll(field, -1, axis=0)
            + np.roll(field, 1, axis=0)
            + np.roll(field, -1, axis=1)
            + np.roll(field, 1, axis=1)
            - 4.0 * field
        ) / self.dx ** 2

    def step(self) -> None:
        """Advance one forward-Euler step. ``w`` carries no diffusion term."""
        lap_u = self.D_u * self._laplacian(self.u)
        lap_v = self.D_v * self._laplacian(self.v)

        f = self.kinetics.f(self.u, self.v, self.w)
        g = self.kinetics.g(self.u, self.v, self.w)
        h = self.kinetics.h(self.u, self.v, self.w)

        self.u = self.u + self.dt * (lap_u + f)
        self.v = self.v + self.dt * (lap_v + g)
        self.w = self.w + self.dt * h  # immobile: no Laplacian
