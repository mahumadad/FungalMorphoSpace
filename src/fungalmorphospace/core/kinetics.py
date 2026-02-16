#!/usr/bin/env python3
"""Reaction kinetics modules for Turing pattern formation in fungal hymenophores.

This module implements the reaction (local-kinetics) component of two-species
activator--inhibitor reaction-diffusion (RD) systems used to model
hymenophore morphogenesis in Basidiomycota.  Three classical kinetic
schemes are provided:

* **Schnakenberg** (1979) -- cubic autocatalysis; best suited for
  labyrinths and connected stripe networks at low diffusivity ratios.
* **Gierer--Meinhardt** (1972) -- saturating autocatalysis with
  u squared / v nonlinearity; produces discrete spot (pore) arrays and is the
  recommended default for high D_v / D_u ratios encountered in fungal
  species such as *Fomes fomentarius*.
* **Gray--Scott** (1984) -- feed-and-kill dynamics; generates an
  exceptionally rich zoo of patterns (spots, stripes, mitosis, solitons)
  depending on the (F, k) control parameters.

All models derive from the abstract base class :class:`TuringKinetics`,
which enforces a uniform interface consumed by
:class:`~fungalmorphospace.core.turing_simulator.TuringSimulator`.

References:
    Gierer, A. & Meinhardt, H. (1972). A theory of biological pattern
        formation. *Kybernetik*, 12, 30--39.
    Schnakenberg, J. (1979). Simple chemical reaction systems with limit
        cycle behaviour. *J. Theor. Biol.*, 81, 389--400.
    Gray, P. & Scott, S.K. (1984). Autocatalytic reactions in the isothermal,
        continuous stirred tank reactor. *Chem. Eng. Sci.*, 39, 1087--1097.

Author:
    Atlas Morphospace Project

FungalMorphoSpace v0.4.0
"""

import numpy as np
from typing import Tuple, Dict
from abc import ABC, abstractmethod


# =============================================================================
# Abstract base class
# =============================================================================

class TuringKinetics(ABC):
    """Abstract base class defining the interface for all RD kinetic models.

    Every concrete kinetic model must implement four methods that together
    fully specify the local (non-diffusive) dynamics of a two-species
    activator--inhibitor system:

    * :meth:`f` -- activator production / degradation rate.
    * :meth:`g` -- inhibitor production / degradation rate.
    * :meth:`get_steady_state` -- spatially homogeneous fixed point
      (u0, v0) around which pattern-forming instabilities are analysed.
    * :meth:`get_name` -- human-readable model identifier used in logging
      and plot titles.

    The simulator calls ``f`` and ``g`` element-wise on NumPy arrays,
    so implementations must be vectorised (i.e., use array arithmetic
    rather than Python loops).
    """

    @abstractmethod
    def f(self, u: np.ndarray, v: np.ndarray) -> np.ndarray:
        """Compute the activator reaction term f(u, v).

        Args:
            u: Activator concentration field. Shape ``(N, N)`` where *N* is
                the spatial grid size.
            v: Inhibitor concentration field, same shape as *u*.

        Returns:
            Pointwise activator reaction rate, same shape as *u*.
        """
        pass

    @abstractmethod
    def g(self, u: np.ndarray, v: np.ndarray) -> np.ndarray:
        """Compute the inhibitor reaction term g(u, v).

        Args:
            u: Activator concentration field. Shape ``(N, N)``.
            v: Inhibitor concentration field, same shape as *u*.

        Returns:
            Pointwise inhibitor reaction rate, same shape as *u*.
        """
        pass

    @abstractmethod
    def get_steady_state(self) -> Tuple[float, float]:
        """Return the spatially homogeneous steady state (u0, v0).

        The steady state satisfies ``f(u0, v0) = 0`` and ``g(u0, v0) = 0``
        simultaneously and is used by the simulator to set up initial
        conditions (small perturbation around the fixed point) and to
        compute the pattern energy (deviation from the homogeneous state).

        Returns:
            A 2-tuple ``(u0, v0)`` of steady-state concentrations.
        """
        pass

    @abstractmethod
    def get_name(self) -> str:
        """Return a human-readable name for this kinetic model.

        Returns:
            Model name string, e.g. ``"Schnakenberg"`` or
            ``"Gierer-Meinhardt"``.
        """
        pass


# =============================================================================
# Schnakenberg kinetics
# =============================================================================

class SchnakenbergKinetics(TuringKinetics):
    """Schnakenberg (1979) cubic-autocatalysis kinetic model.

    The reaction terms are::

        f(u, v) = a - u + u^2 * v        (activator)
        g(u, v) = b - u^2 * v             (inhibitor)

    The homogeneous steady state is ``u0 = a + b``, ``v0 = b / (a + b)^2``.

    Pattern-formation condition (Turing instability):
        A necessary condition derived from linear stability analysis is

            b > (a + 1)^2

        When this inequality is violated no spatially heterogeneous mode is
        unstable and the system converges back to the homogeneous state (or
        diverges numerically).

    Numerical stability warning:
        The cubic ``u^2 * v`` coupling makes the Schnakenberg scheme
        prone to explosive growth when the diffusivity ratio ``D_v / D_u``
        exceeds approximately 50--100.  For fungal species that require
        high ratios (e.g. *Fomes fomentarius*, D_v/D_u ~ 150), the
        Gierer--Meinhardt model is strongly recommended instead.

    Attributes:
        a: Source term parameter for the activator.
        b: Source term parameter for the inhibitor.
        u0: Activator steady-state concentration.
        v0: Inhibitor steady-state concentration.

    Example:
        >>> kin = SchnakenbergKinetics(a=0.1, b=1.3)
        >>> kin.get_steady_state()
        (1.4, 0.6632653061224489)
    """

    def __init__(self, a: float = 0.1, b: float = 0.9) -> None:
        """Initialise Schnakenberg kinetics and validate the pattern condition.

        Args:
            a: Activator source rate.  Typical range: 0.05--0.2.
                Defaults to 0.1.
            b: Inhibitor source rate.  Must satisfy ``b > (a + 1)^2``
                for Turing patterns.  Defaults to 0.9 (which intentionally
                violates the condition to trigger the warning).

        Raises:
            UserWarning: If the Turing pattern condition ``b > (a + 1)^2``
                is not satisfied.  The warning includes a detailed
                diagnostic message with recommended parameter values.
        """
        import warnings

        self.a: float = a
        self.b: float = b

        # --- Turing pattern condition check ---
        # From linear stability analysis the diffusion-driven instability
        # requires b > (a + 1)^2.  Below this threshold the homogeneous
        # state is stable to all spatial perturbations.
        b_min: float = (a + 1)**2
        if b <= b_min:
            warnings.warn(
                f"\n{'='*70}\n"
                f"⚠️  SCHNAKENBERG PATTERN CONDITION VIOLATED ⚠️\n"
                f"{'='*70}\n"
                f"Parameter b={b:.3f} is TOO LOW to form patterns!\n"
                f"\n"
                f"EXPLANATION:\n"
                f"  For Turing patterns: b > (a+1)²\n"
                f"  With a={a:.3f}:\n"
                f"    Required: b > {b_min:.3f}\n"
                f"    You have: b = {b:.3f}\n"
                f"\n"
                f"CONSEQUENCES:\n"
                f"  • Numerical instability (overflow/NaN)\n"
                f"  • No pattern will form\n"
                f"\n"
                f"SOLUTION:\n"
                f"  SchnakenbergKinetics(a=0.1, b=1.3)  # ✓ Correct b\n"
                f"\n"
                f"⚠️  ADDITIONAL WARNING:\n"
                f"  Even with correct b, Schnakenberg fails with high D_v/D_u!\n"
                f"  Species like 'Fomes fomentarius' use D_v/D_u=150 which is\n"
                f"  TOO HIGH for numerical stability in Schnakenberg.\n"
                f"\n"
                f"For PORES (spots), use Gierer-Meinhardt instead:\n"
                f"  create_kinetics('gierer_meinhardt')\n"
                f"{'='*70}\n",
                UserWarning,
                stacklevel=2
            )

        # Cache the steady state once at construction time
        self.u0: float
        self.v0: float
        self.u0, self.v0 = self.get_steady_state()

    def f(self, u: np.ndarray, v: np.ndarray) -> np.ndarray:
        """Activator reaction: f(u, v) = a - u + u^2 * v.

        The ``u^2 * v`` term represents cubic autocatalysis: the activator
        catalyses its own production in the presence of the inhibitor
        substrate *v*.

        Args:
            u: Activator field, shape ``(N, N)``.
            v: Inhibitor field, shape ``(N, N)``.

        Returns:
            Reaction rate array, same shape as inputs.
        """
        # a: constant production, -u: linear decay, u^2*v: autocatalysis
        return self.a - u + u**2 * v

    def g(self, u: np.ndarray, v: np.ndarray) -> np.ndarray:
        """Inhibitor reaction: g(u, v) = b - u^2 * v.

        The inhibitor is consumed by the same cubic term that produces
        the activator, coupling the two species.

        Args:
            u: Activator field, shape ``(N, N)``.
            v: Inhibitor field, shape ``(N, N)``.

        Returns:
            Reaction rate array, same shape as inputs.
        """
        # b: constant production, -u^2*v: consumption by autocatalysis
        return self.b - u**2 * v

    def get_steady_state(self) -> Tuple[float, float]:
        """Compute the homogeneous steady state analytically.

        Setting f = g = 0 and solving yields::

            u0 = a + b
            v0 = b / (a + b)^2

        Returns:
            Tuple of ``(u0, v0)``.
        """
        u0: float = self.a + self.b
        # v0 derived by substituting u0 into g(u0, v0) = 0 => v0 = b / u0^2
        v0: float = self.b / u0**2
        return u0, v0

    def get_name(self) -> str:
        """Return the model name.

        Returns:
            The string ``"Schnakenberg"``.
        """
        return "Schnakenberg"


# =============================================================================
# Gierer-Meinhardt kinetics
# =============================================================================

class GiererMeinhardtKinetics(TuringKinetics):
    """Gierer--Meinhardt (1972) saturating-autocatalysis kinetic model.

    The reaction terms are::

        f(u, v) = rho * (a  -  b * u  +  u^2 / (v + eps))   (activator)
        g(u, v) = rho * (u^2  -  v)                          (inhibitor)

    where ``eps = 1e-6`` is a small safety constant that prevents the
    catastrophic division-by-zero singularity when the inhibitor
    concentration *v* approaches zero.

    The key biological feature is **saturating autocatalysis**: the
    ``u^2 / v`` term means that activator production increases with *u*
    but is simultaneously suppressed by high inhibitor levels.  This
    produces sharply localised peaks (spots/pores) that resist merging,
    making the model ideal for simulating the regular pore arrays found
    on the hymenophore of poroid fungi such as *Fomes fomentarius*.

    The homogeneous steady state is::

        u0 = (a + 1) / b
        v0 = u0^2

    Attributes:
        rho: Overall reaction rate scale factor.
        a: Basal activator production rate.
        b: Activator linear degradation rate.
        u0: Activator steady-state concentration.
        v0: Inhibitor steady-state concentration.

    Example:
        >>> kin = GiererMeinhardtKinetics(rho=0.2, a=0.1, b=1.0)
        >>> kin.get_steady_state()
        (1.1, 1.2100000000000002)
    """

    def __init__(self, rho: float = 1.0, a: float = 0.1, b: float = 1.0) -> None:
        """Initialise Gierer--Meinhardt kinetics.

        Args:
            rho: Global reaction rate scaling factor.  Lower values
                (e.g. 0.2) slow the kinetics relative to diffusion,
                which can improve numerical stability for large grids.
                Defaults to 1.0.
            a: Basal (concentration-independent) activator production
                rate.  Defaults to 0.1.
            b: Linear activator degradation rate coefficient.
                Defaults to 1.0.
        """
        self.rho: float = rho
        self.a: float = a
        self.b: float = b

        # Cache steady state at construction
        self.u0: float
        self.v0: float
        self.u0, self.v0 = self.get_steady_state()

    def f(self, u: np.ndarray, v: np.ndarray) -> np.ndarray:
        """Activator reaction: f(u, v) = rho * (a - b*u + u^2 / (v + eps)).

        A small epsilon (``1e-6``) is added to *v* in the denominator to
        prevent catastrophic floating-point overflow.  This value is
        coordinated with the ``v_floor = 1e-6`` clamp applied by the
        simulator after each time step: together they guarantee that the
        ``u^2 / v`` ratio remains bounded even during transient excursions.

        Args:
            u: Activator field, shape ``(N, N)``.
            v: Inhibitor field, shape ``(N, N)``.

        Returns:
            Reaction rate array, same shape as inputs.
        """
        # Safety epsilon prevents u^2/(v+eps) from diverging when v -> 0.
        # Must be >= v_floor in the simulator to be effective.
        # Previous value 1e-10 was too small and caused overflow; 1e-6
        # provides robust numerical stability.
        epsilon: float = 1e-6  # Increased from 1e-10 for numerical robustness
        # a: basal production, -b*u: linear decay, u^2/(v+eps): autocatalysis
        return self.rho * (self.a - self.b * u + u**2 / (v + epsilon))

    def g(self, u: np.ndarray, v: np.ndarray) -> np.ndarray:
        """Inhibitor reaction: g(u, v) = rho * (u^2 - v).

        The inhibitor is produced proportionally to ``u^2`` (activated by
        the activator) and undergoes first-order decay (``-v``).

        Args:
            u: Activator field, shape ``(N, N)``.
            v: Inhibitor field, shape ``(N, N)``.

        Returns:
            Reaction rate array, same shape as inputs.
        """
        # u^2: production driven by activator, -v: linear self-decay
        return self.rho * (u**2 - v)

    def get_steady_state(self) -> Tuple[float, float]:
        """Compute the homogeneous steady state analytically.

        Setting f = g = 0 (ignoring the epsilon correction) yields::

            u0 = (a + 1) / b
            v0 = u0^2

        Returns:
            Tuple of ``(u0, v0)``.
        """
        # From g=0: v0 = u0^2;  substituting into f=0: a - b*u0 + 1 = 0
        u0: float = (self.a + 1) / self.b
        v0: float = u0**2
        return u0, v0

    def get_name(self) -> str:
        """Return the model name.

        Returns:
            The string ``"Gierer-Meinhardt"``.
        """
        return "Gierer-Meinhardt"


# =============================================================================
# Gray-Scott kinetics
# =============================================================================

class GrayScottKinetics(TuringKinetics):
    """Gray--Scott (1984) feed-and-kill kinetic model.

    The reaction terms are::

        f(u, v) = -u * v^2  +  F * (1 - u)    (activator / substrate)
        g(u, v) =  u * v^2  -  (F + k) * v    (inhibitor / autocatalyst)

    where *F* is the feed rate (replenishment of *u* from an external
    reservoir) and *k* is the kill rate (removal of *v*).

    The trivial steady state is ``(u0, v0) = (1, 0)``.  The system must
    be seeded with a finite *v* perturbation to trigger pattern formation.
    The (F, k) parameter space is famously rich, producing spots,
    stripes, moving pulses, mitosis, and complex labyrinthine structures.

    Typical parameter ranges for interesting patterns:
        * F in [0.01, 0.10]
        * k in [0.045, 0.070]

    Attributes:
        F: Feed rate for the substrate *u*.
        k: Kill (removal) rate for the autocatalyst *v*.
        u0: Substrate steady-state concentration (always 1.0).
        v0: Autocatalyst steady-state concentration (always 0.0).

    Example:
        >>> kin = GrayScottKinetics(F=0.042, k=0.063)
        >>> kin.get_name()
        'Gray-Scott'
    """

    def __init__(self, F: float = 0.020, k: float = 0.050) -> None:
        """Initialise Gray--Scott kinetics.

        Args:
            F: Feed rate -- controls how quickly the substrate *u* is
                replenished from the external reservoir.  Higher *F*
                favours spot splitting and mitosis.  Defaults to 0.020.
            k: Kill rate -- controls how quickly the autocatalyst *v* is
                removed from the system.  Higher *k* favours spot
                extinction and can lead to travelling pulses.
                Defaults to 0.050.
        """
        self.F: float = F
        self.k: float = k
        # Trivial (pattern-free) steady state: u=1 everywhere, v=0 everywhere
        self.u0: float = 1.0
        self.v0: float = 0.0

    def f(self, u: np.ndarray, v: np.ndarray) -> np.ndarray:
        """Substrate reaction: f(u, v) = -u * v^2 + F * (1 - u).

        The first term ``-u * v^2`` represents consumption of the
        substrate by the autocatalytic reaction, while ``F * (1 - u)``
        models replenishment from an external reservoir at rate *F*.

        Args:
            u: Substrate (activator) field, shape ``(N, N)``.
            v: Autocatalyst (inhibitor) field, shape ``(N, N)``.

        Returns:
            Reaction rate array, same shape as inputs.
        """
        # -u*v^2: autocatalytic consumption, F*(1-u): reservoir feed
        return -u * v**2 + self.F * (1.0 - u)

    def g(self, u: np.ndarray, v: np.ndarray) -> np.ndarray:
        """Autocatalyst reaction: g(u, v) = u * v^2 - (F + k) * v.

        The ``u * v^2`` term is the autocatalytic production (the reverse
        of the substrate consumption in :meth:`f`), while ``(F + k) * v``
        combines dilution by feed (*F*) and chemical kill (*k*).

        Args:
            u: Substrate field, shape ``(N, N)``.
            v: Autocatalyst field, shape ``(N, N)``.

        Returns:
            Reaction rate array, same shape as inputs.
        """
        # u*v^2: autocatalytic production, -(F+k)*v: combined removal
        return u * v**2 - (self.F + self.k) * v

    def get_steady_state(self) -> Tuple[float, float]:
        """Return the trivial homogeneous steady state.

        For Gray--Scott the trivial fixed point is ``(1, 0)``; the
        non-trivial steady states (when they exist) depend on *F* and *k*
        and require numerical treatment.  The simulator uses the trivial
        state as the baseline for initial perturbation and energy
        calculations.

        Returns:
            Tuple ``(1.0, 0.0)``.
        """
        return self.u0, self.v0

    def get_name(self) -> str:
        """Return the model name.

        Returns:
            The string ``"Gray-Scott"``.
        """
        return "Gray-Scott"


# =============================================================================
# Factory function
# =============================================================================

def create_kinetics(model_name: str, **kwargs: float) -> TuringKinetics:
    """Factory function to instantiate a kinetic model by name.

    This is the recommended entry point for creating kinetic model objects.
    It normalises the *model_name* string (lower-case, hyphens to
    underscores) and supports short aliases for convenience.

    Supported names and aliases:
        * ``"schnakenberg"``, ``"sch"`` -- :class:`SchnakenbergKinetics`
        * ``"gierer_meinhardt"``, ``"gierer-meinhardt"``, ``"gm"``
          -- :class:`GiererMeinhardtKinetics`
        * ``"gray_scott"``, ``"gray-scott"``, ``"gs"``
          -- :class:`GrayScottKinetics`

    Any extra keyword arguments are forwarded to the corresponding
    constructor (e.g. ``a``, ``b``, ``rho``, ``F``, ``k``).
    Unrecognised keyword arguments for a given model will raise a
    ``TypeError`` from the underlying ``__init__``.

    Args:
        model_name: Case-insensitive model identifier.  Hyphens are
            automatically converted to underscores before matching.
        **kwargs: Model-specific parameters forwarded to the constructor.

    Returns:
        An instance of the requested :class:`TuringKinetics` subclass.

    Raises:
        ValueError: If *model_name* does not match any known model or
            alias.

    Example:
        >>> kin = create_kinetics('gierer_meinhardt', rho=0.2, a=0.1, b=1.0)
        >>> kin.get_name()
        'Gierer-Meinhardt'
        >>> kin = create_kinetics('gm', rho=0.5)
        >>> kin.get_name()
        'Gierer-Meinhardt'
    """
    # Normalise: lowercase, replace hyphens with underscores for matching
    model_name = model_name.lower().replace('-', '_')

    # Match by substring or exact alias
    if 'schnakenberg' in model_name or model_name == 'sch':
        return SchnakenbergKinetics(**kwargs)
    elif 'gierer' in model_name or model_name == 'gm':
        return GiererMeinhardtKinetics(**kwargs)
    elif 'gray' in model_name or model_name == 'gs':
        return GrayScottKinetics(**kwargs)
    else:
        raise ValueError(f"Unknown model: {model_name}")
