#!/usr/bin/env python3
"""Polymorphic Turing pattern simulator for fungal hymenophore morphogenesis.

This module provides :class:`TuringSimulator`, a forward-Euler solver for
two-species activator--inhibitor reaction-diffusion (RD) systems on a 2-D
periodic square lattice.  It is the computational core of the
**FungalMorphoSpace** project, which models the morphogenesis of
hymenophore structures (pores, lamellae, reticulate networks) observed in
Basidiomycota fruit bodies.

The simulator is *polymorphic*: all reaction kinetics are delegated to a
pluggable :class:`~fungalmorphospace.core.kinetics.TuringKinetics` object,
so the same solver can run Schnakenberg, Gierer--Meinhardt, Gray--Scott,
or any future kinetic scheme without modification.

Numerical scheme:
    * **Spatial discretisation** -- standard 5-point Laplacian stencil
      with periodic boundary conditions implemented via ``numpy.roll``.
    * **Temporal integration** -- explicit forward Euler with an adaptive
      CFL stability check (safety factor 0.8) that automatically reduces
      the user-requested time step if it would violate the von Neumann
      stability criterion ``dt <= dx^2 / (4 * D_max)``.
    * **Blow-up prevention** -- a floor clamp ``v >= 1e-6`` is applied
      after every step to prevent catastrophic divergence of the
      ``u^2 / v`` term in Gierer--Meinhardt kinetics.

Analysis tools:
    * Pattern energy (L2 deviation from the homogeneous steady state)
      with automatic convergence detection.
    * Dominant wavelength extraction via 2-D FFT peak detection.
    * Four-panel diagnostic visualisation (activator, inhibitor, energy
      history, power spectrum).

CLI:
    The module can be run as a script (``python -m fungalmorphospace.core.turing_simulator``)
    with YAML configuration, species-preset JSON look-up, and per-parameter
    command-line overrides.

References:
    Gierer, A. & Meinhardt, H. (1972). *Kybernetik*, 12, 30--39.
    Turing, A.M. (1952). *Phil. Trans. R. Soc. Lond. B*, 237, 37--72.

FungalMorphoSpace v0.4.0
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft2, fftfreq
import argparse
import yaml
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# =============================================================================
# PROJECT IMPORTS (PACKAGED)
# =============================================================================
from .kinetics import (
    TuringKinetics,
    SchnakenbergKinetics,
    GiererMeinhardtKinetics,
    GrayScottKinetics,
    create_kinetics,
)

# =============================================================================
# LOGGING CONFIGURATION
# =============================================================================

def setup_logging(log_level: str = "INFO", log_file: Optional[str] = None) -> None:
    """Configure the root logger for the simulation session.

    Sets up console output and, optionally, a file handler.  Parent
    directories for the log file are created automatically if they do not
    exist.

    Args:
        log_level: Logging verbosity.  Must be a standard ``logging``
            level name (``"DEBUG"``, ``"INFO"``, ``"WARNING"``,
            ``"ERROR"``, ``"CRITICAL"``).  Defaults to ``"INFO"``.
        log_file: Optional path to a log file.  If ``None`` (default),
            only console output is produced.
    """
    handlers: List[logging.Handler] = [logging.StreamHandler()]
    if log_file:
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(
        level=getattr(logging, log_level),
        format='[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=handlers
    )


# =============================================================================
# TURING SYSTEM CLASS (POLYMORPHIC)
# =============================================================================

class TuringSimulator:
    """Polymorphic reaction-diffusion solver on a 2-D periodic lattice.

    The simulator delegates all local reaction kinetics to a
    :class:`~fungalmorphospace.core.kinetics.TuringKinetics` instance,
    so the same integration loop supports any two-species
    activator--inhibitor model.  Diffusion is computed via a 5-point
    Laplacian stencil with periodic (toroidal) boundary conditions.

    The integration uses explicit forward Euler.  At construction time a
    CFL stability check automatically reduces the requested ``dt`` if it
    exceeds the von Neumann limit for the faster-diffusing species,
    including a configurable safety margin (default 0.8) to account for
    stiff reaction terms.

    Attributes:
        kinetics: The pluggable kinetic model instance.
        D_u: Activator diffusion coefficient.
        D_v: Inhibitor diffusion coefficient.
        grid_size: Number of grid points along each spatial axis (N x N).
        dx: Spatial step size (grid spacing).
        dt: Current time step (may differ from the requested value after
            CFL adjustment).
        D_v_D_u_ratio: Ratio ``D_v / D_u`` (controls pattern selection).
        u_steady: Activator homogeneous steady-state concentration.
        v_steady: Inhibitor homogeneous steady-state concentration.
        lambda_theoretical: Heuristic wavelength estimate
            ``2 * pi * sqrt(D_v / D_u)`` (rule-of-thumb, not from the
            full dispersion relation).
        u: Activator concentration field, shape ``(N, N)`` or ``None``
            before :meth:`initialize` is called.
        v: Inhibitor concentration field, same shape, or ``None``.
        energy_history: List of pattern-energy samples recorded during
            :meth:`run`.
        time_history: List of corresponding physical times.
        stability_info: Dictionary of CFL stability diagnostics populated
            by :meth:`_check_stability`.
    """

    def __init__(
        self,
        kinetics_model: Optional[TuringKinetics] = None,
        D_u: float = 0.1,
        D_v: float = 15.0,
        grid_size: int = 256,
        dx: float = 1.0,
        dt: float = 0.01,
        random_seed: Optional[int] = 42,
    ) -> None:
        """Initialise the Turing simulator.

        Args:
            kinetics_model: An instance of a
                :class:`~fungalmorphospace.core.kinetics.TuringKinetics`
                subclass that provides the reaction terms ``f(u, v)`` and
                ``g(u, v)``.  If ``None``, a default
                :class:`GiererMeinhardtKinetics` with ``rho=0.2, a=0.1,
                b=1.0`` is used.  (Schnakenberg was removed as the
                default because its default ``b=0.9`` violates the
                pattern condition ``b > (a+1)^2 = 1.21`` for ``a=0.1``.)
            D_u: Activator diffusion coefficient.  Defaults to 0.1.
            D_v: Inhibitor diffusion coefficient.  Must be larger than
                ``D_u`` for Turing instability (short-range activation,
                long-range inhibition).  Defaults to 15.0.
            grid_size: Number of grid points per spatial dimension.
                The total domain is ``grid_size x grid_size``.
                Defaults to 256.
            dx: Spatial step (grid spacing) in physical units.
                Defaults to 1.0.
            dt: Requested integration time step.  This value is
                automatically reduced if it violates the CFL stability
                condition (see :meth:`_check_stability`).  Defaults to
                0.01.
            random_seed: Seed for NumPy's legacy PRNG to ensure
                reproducible initial perturbations.  Pass ``None`` for
                non-deterministic initialisation.  Defaults to 42.

        Side Effects:
            * Sets ``np.random.seed(random_seed)`` when *random_seed* is
              not ``None``.
            * Invokes :meth:`_check_stability`, which may reduce
              ``self.dt`` below the requested value and populate
              ``self.stability_info``.
        """
        # ----- 1. Kinetic model delegation -----
        if kinetics_model is None:
            # Default to Gierer-Meinhardt which is stable for high D_v/D_u ratios.
            # Note: Schnakenberg default was removed because b=0.9 violates
            # the pattern condition b > (a+1)^2 = 1.21 for a=0.1.
            logging.info("No kinetics model provided, defaulting to Gierer-Meinhardt")
            self.kinetics: TuringKinetics = GiererMeinhardtKinetics(rho=0.2, a=0.1, b=1.0)
        else:
            self.kinetics = kinetics_model

        # ----- 2. Physical / numerical parameters -----
        self.D_u: float = D_u
        self.D_v: float = D_v
        self.grid_size: int = grid_size
        self.dx: float = dx
        self.dt: float = dt

        # ----- 3. Reproducibility seed -----
        if random_seed is not None:
            np.random.seed(random_seed)

        # ----- 4. Derived quantities -----
        self.D_v_D_u_ratio: float = D_v / D_u

        # Obtain the homogeneous steady state from the kinetic model
        self.u_steady: float
        self.v_steady: float
        self.u_steady, self.v_steady = self.kinetics.get_steady_state()

        # Heuristic wavelength proxy (rule-of-thumb):
        #   lambda ~ 2*pi * sqrt(D_v / D_u)
        # NOTE: This is NOT the exact most-unstable wavelength from the
        # full dispersion relation; it is a convenient order-of-magnitude
        # estimate used for sanity checks and logging.
        self.lambda_theoretical: float = 2 * np.pi * np.sqrt(self.D_v_D_u_ratio)

        # ----- 5. Field arrays (populated by initialize()) -----
        self.u: Optional[np.ndarray] = None
        self.v: Optional[np.ndarray] = None
        self.energy_history: List[float] = []
        self.time_history: List[float] = []

        # ----- 6. CFL stability check (may reduce self.dt) -----
        self._check_stability()

        logging.info(f"TuringSimulator initialized with {self.kinetics.get_name()} kinetics")
        logging.info(f"  Steady state: u0={self.u_steady:.3f}, v0={self.v_steady:.3f}")
        logging.info(f"  Heuristic λ~ {self.lambda_theoretical:.2f} (grid units; rule-of-thumb, not a full dispersion prediction)")

    def _check_stability(self) -> None:
        """Enforce the CFL stability condition with a safety margin.

        For a 2-D explicit forward-Euler diffusion scheme on a square
        lattice, the von Neumann stability criterion requires::

            dt  <=  dx^2 / (4 * D_max)

        where ``D_max = max(D_u, D_v)`` is the largest diffusion
        coefficient and the factor 4 arises from the 5-point stencil in
        two dimensions (2 neighbours per axis).

        Because Gierer--Meinhardt and other RD kinetics can be *stiff*
        (the reaction Jacobian has large eigenvalues), operating at the
        bare CFL limit often leads to transient blow-ups.  A
        multiplicative **safety factor** of 0.8 (i.e. 20 % margin) is
        therefore applied::

            dt_max_safe  =  0.8 * dx^2 / (4 * D_max)

        If the user-requested ``dt`` exceeds ``dt_max_safe``, it is
        silently reduced and a warning is logged.  All diagnostic
        quantities are stored in ``self.stability_info`` for
        reproducibility and audit trails.

        Side Effects:
            * May reduce ``self.dt``.
            * Populates ``self.stability_info`` dict with keys:
              ``dt_requested``, ``dt_max_diffusion``, ``safety_factor``,
              ``dt_max_safe``, ``dt_adjusted``, ``had_stability_warning``,
              ``D_max``, ``dt_final``.
        """
        D_max: float = max(self.D_u, self.D_v)

        # Bare CFL limit for 2-D 5-point Laplacian on a square grid:
        # dt_max = dx^2 / (4 * D_max)
        dt_max_diffusion: float = self.dx**2 / (4 * D_max)

        # Apply a multiplicative safety factor to guard against stiff
        # reaction terms that can destabilise the explicit scheme even
        # when the pure-diffusion CFL is satisfied.
        # 0.8 = 20% margin (use 0.5 for very stiff systems)
        safety_factor: float = 0.8
        dt_max_safe: float = safety_factor * dt_max_diffusion

        # Store original request for traceability
        self.dt_requested: float = self.dt

        # Collect stability metadata in a dict for downstream reporting
        self.stability_info: Dict[str, object] = {
            'dt_requested': self.dt_requested,
            'dt_max_diffusion': dt_max_diffusion,
            'safety_factor': safety_factor,
            'dt_max_safe': dt_max_safe,
            'dt_adjusted': False,
            'had_stability_warning': False,
            'D_max': D_max
        }

        if self.dt > dt_max_safe:
            # Requested dt exceeds safe limit -- clamp it down
            self.stability_info['dt_adjusted'] = True
            self.stability_info['had_stability_warning'] = True

            logging.warning(
                f"[CFL STABILITY] dt reduced for numerical stability:\n"
                f"  dt_requested: {self.dt:.6f}\n"
                f"  dt_max_safe:  {dt_max_safe:.6f} (CFL × {safety_factor})\n"
                f"  dt_final:     {dt_max_safe:.6f}\n"
                f"  D_max={D_max:.1f}, dx={self.dx}"
            )
            self.dt = dt_max_safe

        self.stability_info['dt_final'] = self.dt

        # Log for traceability (always, not just on warning)
        if self.stability_info['dt_adjusted']:
            logging.info(f"[dt] {self.dt_requested:.6f} → {self.dt:.6f} (adjusted)")
        else:
            logging.info(f"[dt] {self.dt:.6f} (within CFL bounds)")

    def initialize(self, perturbation_amplitude: float = 0.01) -> None:
        """Set up initial concentration fields as small random perturbations.

        Each field is initialised to its homogeneous steady-state value
        plus Gaussian white noise scaled by *perturbation_amplitude*.
        Negative values are clamped to zero so that concentrations remain
        physically meaningful.  The energy and time histories are reset.

        This must be called before :meth:`run` or :meth:`step`.

        Args:
            perturbation_amplitude: Standard deviation of the Gaussian
                noise added to each steady-state field.  A small value
                (e.g. 0.01) seeds the Turing instability without
                overwhelming the linear regime.  Defaults to 0.01.
        """
        # Draw i.i.d. standard-normal noise for both fields
        u_noise: np.ndarray = np.random.randn(self.grid_size, self.grid_size)
        v_noise: np.ndarray = np.random.randn(self.grid_size, self.grid_size)

        # Perturb around the homogeneous steady state
        self.u = self.u_steady + perturbation_amplitude * u_noise
        self.v = self.v_steady + perturbation_amplitude * v_noise

        # Clamp to non-negative concentrations (physical constraint)
        self.u = np.maximum(self.u, 0)
        self.v = np.maximum(self.v, 0)

        # Reset monitoring histories for fresh run
        self.energy_history = []
        self.time_history = []

    def _laplacian(self, field: np.ndarray) -> np.ndarray:
        """Compute the discrete Laplacian using a 5-point stencil.

        Implements the standard second-order finite-difference
        approximation on a uniform square grid with **periodic boundary
        conditions** (toroidal topology).  Periodicity is achieved via
        ``numpy.roll``, which wraps array elements that shift past an
        edge back to the opposite side.

        The stencil in 2-D is::

            nabla^2 f_{i,j} = (f_{i+1,j} + f_{i-1,j}
                             + f_{i,j+1} + f_{i,j-1}
                             - 4 * f_{i,j}) / dx^2

        Args:
            field: 2-D concentration array, shape ``(N, N)``.

        Returns:
            Discrete Laplacian of *field*, same shape.
        """
        # Shift field one cell in each of the four cardinal directions.
        # np.roll with periodic wrapping implements toroidal BC.
        u_right: np.ndarray = np.roll(field, shift=-1, axis=1)  # f_{i, j+1}
        u_left: np.ndarray = np.roll(field, shift=1, axis=1)    # f_{i, j-1}
        u_up: np.ndarray = np.roll(field, shift=-1, axis=0)     # f_{i+1, j}
        u_down: np.ndarray = np.roll(field, shift=1, axis=0)    # f_{i-1, j}

        # 5-point stencil: sum of 4 neighbours minus 4 times centre, / dx^2
        return (u_right + u_left + u_up + u_down - 4 * field) / self.dx**2

    def _reaction_u(self, u: np.ndarray, v: np.ndarray) -> np.ndarray:
        """Delegate the activator reaction term to the kinetic model.

        Args:
            u: Activator field, shape ``(N, N)``.
            v: Inhibitor field, shape ``(N, N)``.

        Returns:
            Activator reaction rate, same shape.
        """
        return self.kinetics.f(u, v)  # Delegates to Gierer-Meinhardt, Schnakenberg, etc.

    def _reaction_v(self, u: np.ndarray, v: np.ndarray) -> np.ndarray:
        """Delegate the inhibitor reaction term to the kinetic model.

        Args:
            u: Activator field, shape ``(N, N)``.
            v: Inhibitor field, shape ``(N, N)``.

        Returns:
            Inhibitor reaction rate, same shape.
        """
        return self.kinetics.g(u, v)

    def step(self) -> None:
        """Advance the system by one forward-Euler time step.

        The update rule for each species is::

            u^{n+1} = u^n + dt * (D_u * nabla^2 u  +  f(u, v))
            v^{n+1} = v^n + dt * (D_v * nabla^2 v  +  g(u, v))

        After the update, two critical clamps are applied:

        1. ``u`` is clamped to ``>= 0`` (physical non-negativity).
        2. ``v`` is clamped to ``>= v_floor = 1e-6`` -- **not** zero.
           This prevents the ``u^2 / v`` term in Gierer--Meinhardt
           kinetics from producing catastrophic floating-point overflow
           when *v* transiently drops to near-zero values.  The floor
           is coordinated with the ``epsilon = 1e-6`` safety offset in
           :meth:`GiererMeinhardtKinetics.f`.
        """
        # --- Diffusion contributions (Laplacian * diffusivity) ---
        diff_u: np.ndarray = self.D_u * self._laplacian(self.u)
        diff_v: np.ndarray = self.D_v * self._laplacian(self.v)

        # --- Reaction contributions (delegated to kinetic model) ---
        react_u: np.ndarray = self._reaction_u(self.u, self.v)
        react_v: np.ndarray = self._reaction_v(self.u, self.v)

        # --- Forward Euler update ---
        self.u += self.dt * (diff_u + react_u)
        self.v += self.dt * (diff_v + react_v)

        # --- Post-step safety clamps ---
        # v_floor prevents division by near-zero in the u^2/(v+eps) term.
        # If v -> 0, the Gierer-Meinhardt activator term u^2/(v+eps)
        # explodes even with eps=1e-6, so we enforce v >= 1e-6 globally.
        v_floor: float = 1e-6  # Much safer than 0
        self.u = np.maximum(self.u, 0)
        self.v = np.maximum(self.v, v_floor)

    def compute_pattern_energy(self) -> float:
        """Compute the total pattern energy (L2 deviation from steady state).

        The pattern energy is defined as::

            E = sum((u - u0)^2) + sum((v - v0)^2)

        where ``u0`` and ``v0`` are the homogeneous steady-state
        concentrations.  This scalar quantifies how far the system has
        evolved from the uniform state and is used for convergence
        detection in :meth:`run`.

        Returns:
            Non-negative scalar energy value.
        """
        # Deviation of each field from its homogeneous steady state
        u_dev: np.ndarray = self.u - self.u_steady
        v_dev: np.ndarray = self.v - self.v_steady
        # Sum of squared deviations (L2 norm squared, not normalised)
        return float(np.sum(u_dev**2) + np.sum(v_dev**2))

    def run(
        self,
        num_steps: int = 5000,
        check_convergence: bool = True,
        convergence_threshold: float = 1e-6,
        check_interval: int = 100,
    ) -> Dict[str, object]:
        """Run the simulation for a specified number of time steps.

        Optionally monitors the pattern energy at regular intervals and
        stops early when the relative energy change falls below a
        threshold for three consecutive checks (to avoid false positives
        from transient plateaux).

        Args:
            num_steps: Maximum number of forward-Euler steps to execute.
                Defaults to 5000.
            check_convergence: Whether to monitor energy and stop early
                on convergence.  Defaults to ``True``.
            convergence_threshold: Relative energy change
                ``|E_new - E_old| / (E_old + 1e-12)`` below which a
                step counts as "converged".  Three consecutive converged
                checks trigger early termination.  Defaults to 1e-6.
            check_interval: Number of steps between energy evaluations.
                Defaults to 100.

        Returns:
            A dictionary with keys:

            * ``'converged'`` (*bool*) -- whether early stopping was
              triggered.
            * ``'final_step'`` (*int*) -- zero-based index of the last
              executed step.
            * ``'u_final'`` (*np.ndarray*) -- copy of the final activator
              field.
            * ``'v_final'`` (*np.ndarray*) -- copy of the final inhibitor
              field.
        """
        logging.info(f"Running {num_steps} steps...")

        prev_energy: float = self.compute_pattern_energy()
        converged: bool = False
        conv_cnt: int = 0  # Consecutive convergence-check passes

        for step in range(num_steps):
            self.step()

            if step % check_interval == 0:
                # Fix 5 (audit 2026-06-02): re-check stability during the run.
                # The timestep is fixed after init, so a stiff reaction regime
                # could drive a blow-up that the initial CFL check cannot
                # foresee. Detect non-finite fields and stop cleanly with a
                # flag instead of silently producing NaN-filled output. This is
                # a no-op for stable runs (fields remain finite).
                if not (np.all(np.isfinite(self.u)) and np.all(np.isfinite(self.v))):
                    logging.warning(
                        f"Non-finite field at step {step}; stopping (possible stiffness/CFL blow-up). "
                        f"Consider a smaller dt or grid."
                    )
                    self.stability_info["had_stability_warning"] = True
                    self.stability_info["blowup_step"] = int(step)
                    break

                curr: float = self.compute_pattern_energy()
                self.energy_history.append(curr)
                self.time_history.append(step * self.dt)

                if check_convergence and step > 500:
                    # Relative energy change (eps in denominator avoids 0/0)
                    rel: float = abs(curr - prev_energy) / (prev_energy + 1e-12)
                    if rel < convergence_threshold:
                        conv_cnt += 1
                        # Require 3 consecutive passes to confirm convergence
                        # (guards against transient energy plateaux)
                        if conv_cnt >= 3:
                            logging.info(f"Converged at step {step}")
                            converged = True
                            break
                    else:
                        conv_cnt = 0  # Reset counter on any non-converged check
                    prev_energy = curr

        return {
            'converged': converged,
            'final_step': step,
            'u_final': self.u.copy(),
            'v_final': self.v.copy()
        }

    def measure_wavelength(self) -> float:
        """Estimate the dominant pattern wavelength via 2-D FFT peak detection.

        Procedure:
            1. Subtract the steady-state mean from the activator field to
               isolate the pattern signal.
            2. Compute the 2-D FFT and its power spectrum.
            3. Zero the DC component (k=0) to ignore the mean.
            4. Locate the peak in the radial frequency domain.
            5. Convert the peak spatial frequency to a wavelength
               ``lambda = 1 / f_peak``.

        Returns:
            Dominant wavelength in spatial units determined by ``dx``.
            Returns ``numpy.nan`` if no dominant peak is found (i.e. the
            peak frequency is zero, indicating a flat field).
        """
        # 1. Compute 2-D FFT of the activator deviation from steady state
        fft_result: np.ndarray = fft2(self.u - self.u_steady)

        # 2. Power spectrum (squared magnitude)
        power: np.ndarray = np.abs(fft_result)**2

        # 3. Zero out DC (mean) component so it does not dominate the peak search
        power[0, 0] = 0

        # 4. Build a radial-frequency grid matching the FFT layout
        freqs: np.ndarray = fftfreq(self.grid_size, d=self.dx)
        fx: np.ndarray
        fy: np.ndarray
        fx, fy = np.meshgrid(freqs, freqs)
        fr: np.ndarray = np.sqrt(fx**2 + fy**2)  # Radial spatial frequency

        # 5. Find the 2-D index of the maximum power and read its frequency
        idx: Tuple[int, ...] = np.unravel_index(np.argmax(power), power.shape)
        peak_freq: float = float(fr[idx])

        if peak_freq > 0:
            wl: float = 1.0 / peak_freq
            # Choose a sensible unit label depending on grid spacing
            unit: str = "px" if self.dx == 1.0 else "μm"
            logging.info(f"Measured λ: {wl:.2f} {unit} (Heuristic λ~ {self.lambda_theoretical:.2f} grid units)")
            return wl
        return np.nan

    def visualize(
        self,
        save_path: Optional[str] = None,
        show: bool = True,
    ) -> None:
        """Generate a four-panel diagnostic figure of the simulation state.

        The panels are:

        1. **Top-left** -- Activator field *u* (``viridis`` colourmap).
        2. **Top-right** -- Inhibitor field *v* (``plasma`` colourmap).
        3. **Bottom-left** -- Pattern energy vs. physical time on a
           logarithmic y-axis (only if energy samples exist).
        4. **Bottom-right** -- Log10 power spectrum of the activator
           deviation, centred via ``fftshift`` (``inferno`` colourmap).

        Args:
            save_path: If provided, the figure is saved at 200 DPI to
                this path.  Parent directories are **not** auto-created
                (call ``Path.mkdir`` beforehand).  Defaults to ``None``.
            show: If ``True``, display the figure interactively via
                ``plt.show()``.  If ``False``, the figure is closed
                immediately after optional saving (useful in batch mode).
                Defaults to ``True``.
        """
        fig, axes = plt.subplots(2, 2, figsize=(10, 8))

        # --- Panel 1: Activator concentration ---
        axes[0, 0].imshow(self.u, cmap='viridis', origin='lower')
        axes[0, 0].set_title(f'Activator ({self.kinetics.get_name()})')

        # --- Panel 2: Inhibitor concentration ---
        axes[0, 1].imshow(self.v, cmap='plasma', origin='lower')
        axes[0, 1].set_title('Inhibitor')

        # --- Panel 3: Pattern energy history (log scale) ---
        if self.energy_history:
            axes[1, 0].plot(self.time_history, self.energy_history)
            axes[1, 0].set_yscale('log')
            axes[1, 0].set_title('Energy')

        # --- Panel 4: Centred log-power spectrum ---
        fft_result: np.ndarray = fft2(self.u - self.u_steady)
        # fftshift centres the zero-frequency component for display
        pwr: np.ndarray = np.abs(np.fft.fftshift(fft_result))**2
        # log10(pwr + 1) avoids log(0) while preserving dynamic range
        axes[1, 1].imshow(np.log10(pwr + 1), cmap='inferno')
        axes[1, 1].set_title('Power Spectrum')

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=200)
            logging.info(f"Saved: {save_path}")
        if show:
            plt.show()
        else:
            plt.close()


# =============================================================================
# CLI ENTRY POINT (supports YAML config + species presets + CLI overrides)
# =============================================================================

def main() -> int:
    """Command-line entry point for running Turing pattern simulations.

    The configuration is assembled from four layers, applied in order of
    increasing priority:

    1. **Hard-coded defaults** -- safe baseline parameters that always
       produce a valid (if generic) simulation.
    2. **YAML configuration file** -- loaded from ``--config`` (default
       ``config/turing_params.yaml``); merged on top of the defaults.
    3. **Species preset** -- loaded from ``data/species_data.json`` via
       ``--species`` (e.g. ``fomes``, ``brumalis``); merged on top.
    4. **Explicit CLI flags** -- ``--a``, ``--b``, ``--rho``, ``--F``,
       ``--k``, ``--D_v_D_u``, ``--grid_size``, ``--time_steps``,
       ``--T_target``; these override everything.

    The ``--model`` flag selects the kinetic scheme (``auto`` delegates to
    whatever the config/species specifies; defaults to
    ``gierer_meinhardt``).

    Returns:
        Exit code ``0`` on success.
    """
    parser = argparse.ArgumentParser(description='Simulador de Patrones de Turing para Himenóforos Fúngicos')
    parser.add_argument('--config', type=str, default='config/turing_params.yaml')
    parser.add_argument('--species', type=str, help='Species preset key (e.g., fomes, brumalis, squamosus)')
    parser.add_argument('--model', type=str, default='auto',
                        choices=['auto', 'schnakenberg', 'gierer_meinhardt', 'gray_scott'],
                        help='Kinetic model to use (auto = from config/species; fallback GM)')
    parser.add_argument('--a', type=float)
    parser.add_argument('--b', type=float)
    parser.add_argument('--rho', type=float)
    parser.add_argument('--F', type=float)
    parser.add_argument('--k', type=float)
    parser.add_argument('--D_v_D_u', type=float)
    parser.add_argument('--grid_size', type=int)
    parser.add_argument('--time_steps', type=int, help='Total integration steps (overrides T_target)')
    parser.add_argument('--T_target', type=float, help='Target physical time; steps computed from dt_final when total_steps not provided')
    parser.add_argument('--output', type=str, default='results/patterns/simulation.png')
    parser.add_argument('--test', action='store_true')
    parser.add_argument('--no-show', action='store_true')

    args: argparse.Namespace = parser.parse_args()
    setup_logging()

    # -----------------------------------------------------------------
    # Helper: shallow deep-merge of nested dicts
    # -----------------------------------------------------------------
    def _deep_merge(base: dict, incoming: dict) -> dict:
        """Shallow deep-merge: merges first-level nested dicts.

        For each key in *incoming*, if both the base and incoming values
        are dicts, their contents are merged (incoming wins on
        conflicts).  Otherwise the incoming value replaces the base.

        Args:
            base: The base configuration dictionary.
            incoming: Overriding configuration dictionary.

        Returns:
            A new merged dictionary (neither input is mutated).
        """
        out: dict = dict(base)
        for k, v in (incoming or {}).items():
            if isinstance(v, dict) and isinstance(out.get(k), dict):
                out[k] = {**out[k], **v}
            else:
                out[k] = v
        return out

    # -----------------------------------------------------------------
    # Helper: resolve a relative path from CWD or repository root
    # -----------------------------------------------------------------
    def _resolve_repo_path(rel: str) -> Path:
        """Resolve a relative path, trying CWD first then the repo root.

        Args:
            rel: Relative file path string.

        Returns:
            The resolved :class:`Path` (may or may not exist on disk).
        """
        p: Path = Path(rel)
        if p.exists():
            return p
        # Fall back to repository root (3 levels above this file)
        repo_root: Path = Path(__file__).resolve().parents[3]
        alt: Path = repo_root / rel
        return alt if alt.exists() else p

    # =================================================================
    # 1) Default config (safe baseline)
    # =================================================================
    config: dict = {
        'diffusion': {'D_u': 1.0, 'D_v_D_u_ratio': 150.0},
        'grid': {'size': 512, 'dx': 1.0},
        'time': {'dt': 0.0005, 'T_target': 5.0},
        'initial': {'random_seed': 42, 'perturbation_amplitude': 0.01},
        'kinetics': {'model': 'gierer_meinhardt', 'params': {'a': 0.1, 'b': 1.0, 'rho': 0.2}},
    }

    # =================================================================
    # 2) Load YAML config if present (merges on top of defaults)
    # =================================================================
    cfg_path: Path = _resolve_repo_path(args.config)
    if cfg_path.exists():
        with open(cfg_path, 'r') as f:
            loaded: dict = yaml.safe_load(f) or {}
        config = _deep_merge(config, loaded)
        logging.info(f"Loaded config: {cfg_path}")
    else:
        logging.warning(f"Config not found: {args.config} (using defaults)")

    # =================================================================
    # 3) Apply species preset (if provided via --species)
    # =================================================================
    if args.species:
        sp_path: Path = _resolve_repo_path('data/species_data.json')
        if not sp_path.exists():
            logging.warning(f"Species data not found: {sp_path} (ignoring --species)")
        else:
            with open(sp_path, 'r') as f:
                sdata: dict = json.load(f)

            # Canonical keys are short (fomes, brumalis, squamosus).
            # We also accept aliases listed under each species entry.
            raw_key: str = args.species.strip()
            key_norm: str = raw_key.lower().strip().replace("-", "_")

            species_key: Optional[str] = None
            species_map: dict = sdata.get('species', {})

            # Direct canonical match
            if key_norm in species_map:
                species_key = key_norm
            else:
                # Alias / scientific name match
                for canon, entry in species_map.items():
                    aliases: list = [str(a).lower().strip().replace("-", "_") for a in entry.get("aliases", [])]
                    sci: str = str(entry.get("scientific_name", "")).lower().strip()
                    if key_norm in aliases or key_norm == sci:
                        species_key = canon
                        break

            if species_key in species_map:
                species_data: dict = sdata['species'][species_key]
                preset: dict = species_data.get('parameters', species_data.get('turing_parameters', {})) or {}

                # --- Diffusion / grid overrides from species preset ---
                if 'D_v_D_u' in preset:
                    config['diffusion']['D_v_D_u_ratio'] = float(preset['D_v_D_u'])
                if 'grid_size' in preset:
                    config['grid']['size'] = int(preset['grid_size'])

                # --- Time overrides ---
                if 'dt_initial' in preset:
                    config['time']['dt'] = float(preset['dt_initial'])
                if 'T_target' in preset:
                    config['time']['T_target'] = float(preset['T_target'])

                # --- Kinetics parameter overrides (generic; unused params are ignored) ---
                config.setdefault('kinetics', {'model': 'gierer_meinhardt', 'params': {}})
                config['kinetics'].setdefault('params', {})
                for key in ('a', 'b', 'rho', 'F', 'k'):
                    if key in preset:
                        config['kinetics']['params'][key] = float(preset[key])

                logging.info(f"Loaded species preset: {args.species} -> {species_key}")
            else:
                logging.warning(f"Unknown species preset: {args.species}")

    # =================================================================
    # 4) CLI overrides (explicit flags take highest priority)
    # =================================================================
    if args.D_v_D_u is not None:
        config['diffusion']['D_v_D_u_ratio'] = float(args.D_v_D_u)
    if args.grid_size is not None:
        config['grid']['size'] = int(args.grid_size)
    if args.T_target is not None:
        config['time']['T_target'] = float(args.T_target)
    if args.time_steps is not None:
        config['time']['total_steps'] = int(args.time_steps)

    # --test shortcut: small grid + few steps for quick validation
    if args.test:
        config['grid']['size'] = 64
        config['time']['total_steps'] = 500

    # =================================================================
    # 5) Kinetics model creation (factory dispatch)
    # =================================================================
    cfg_kin: dict = config.get('kinetics', {}) or {}
    cfg_model: str = (cfg_kin.get('model') or 'gierer_meinhardt')
    cfg_params: dict = (cfg_kin.get('params') or {}).copy()

    # CLI param overrides (highest priority)
    if args.a is not None:
        cfg_params['a'] = args.a
    if args.b is not None:
        cfg_params['b'] = args.b
    if args.rho is not None:
        cfg_params['rho'] = args.rho
    if args.F is not None:
        cfg_params['F'] = args.F
    if args.k is not None:
        cfg_params['k'] = args.k

    model_name: str = cfg_model if args.model == 'auto' else args.model
    logging.info(f"Creating model: {model_name} with params {cfg_params}")
    kinetics: TuringKinetics = create_kinetics(model_name, **cfg_params)

    # =================================================================
    # 6) Instantiate the simulator
    # =================================================================
    D_u: float = float(config['diffusion']['D_u'])
    D_v: float = D_u * float(config['diffusion']['D_v_D_u_ratio'])

    sim: TuringSimulator = TuringSimulator(
        kinetics_model=kinetics,
        D_u=D_u,
        D_v=D_v,
        grid_size=int(config['grid']['size']),
        dx=float(config['grid']['dx']),
        dt=float(config['time']['dt']),
        random_seed=int(config['initial']['random_seed']),
    )

    # =================================================================
    # 7) Execute the simulation
    # =================================================================
    sim.initialize(perturbation_amplitude=float(config['initial']['perturbation_amplitude']))

    # Time control: prefer explicit total_steps; otherwise derive from
    # T_target using the (possibly CFL-adjusted) dt_final.
    time_cfg: dict = (config.get('time') or {})
    total_steps: Optional[int] = time_cfg.get('total_steps')
    if total_steps is None:
        T_target: Optional[float] = time_cfg.get('T_target')
        if T_target is None:
            total_steps = 5000
            logging.warning("No time.total_steps or time.T_target provided; defaulting to 5000 steps.")
        else:
            # Derive step count: ceil(T_target / dt_final)
            total_steps = int(np.ceil(float(T_target) / float(sim.dt)))
            logging.info(f"[time] Derived steps from T_target={T_target} with dt_final={sim.dt:.6g}: steps={total_steps}")
    else:
        logging.info(f"[time] Using explicit total_steps={total_steps} (dt_final={sim.dt:.6g})")

    logging.info(f"[time] Physical time simulated: T = {float(total_steps) * float(sim.dt):.6g}")

    sim.run(num_steps=int(total_steps))
    sim.measure_wavelength()

    # Save and optionally display results
    out_path: Path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    sim.visualize(save_path=str(out_path), show=not args.no_show)

    return 0

if __name__ == '__main__':
    main()
