#!/usr/bin/env python3
"""Topology Analyzer -- Quantitative Morphometry of Turing Patterns.

.. module:: fungalmorphospace.analysis.topology_analyzer
   :synopsis: Extract topological and morphometric features from 2D
              reaction-diffusion patterns for biological validation.

Overview
--------
This module provides the :class:`TopologyAnalyzer` class, which extracts
quantitative morphometric features from simulated (or experimental) Turing
patterns to enable objective comparison between simulated and biological
structures.  It is designed for the *FungalMorphoSpace* pipeline but is
general enough for any 2D spatial pattern analysis.

Topological Metrics
~~~~~~~~~~~~~~~~~~~
Topology studies properties preserved under continuous deformation. For
spatial patterns, key topological features include:

1. **Euler Characteristic** (chi = V - E + F)

   * chi > 0 : Spots / islands (disconnected regions).
   * chi = 0 : Labyrinth / maze (connected network).
   * chi < 0 : Holes / voids (connected with enclosed regions).

2. **Connectivity** -- number of connected components, useful for
   distinguishing spot arrays from labyrinthine networks.

3. **Regularity** -- coefficient of variation of nearest-neighbor
   spacing among connected-component centroids.  Low CV indicates a
   highly ordered lattice; high CV indicates disorder.

Morphometric Metrics
~~~~~~~~~~~~~~~~~~~~
These quantify the *geometry* (not just topology) of patterns:

1. **Wavelength** (lambda) -- dominant spatial periodicity measured via
   two independent methods (radial FFT power spectrum and Wiener-Khinchin
   autocorrelation).  Primary validation metric against empirical pore /
   gill / tooth spacing in polypore fungi.

2. **Isotropy index** -- angular uniformity of the 2D power spectrum.
   Values near 1.0 indicate no preferred direction; values near 0.0
   indicate strong stripe-like anisotropy.

Biological Relevance
~~~~~~~~~~~~~~~~~~~~
For fungal hymenophores:

* lambda  maps to pore spacing, gill spacing, or tooth spacing.
* Regularity (CV) maps to developmental precision.
* Connectivity maps to structural integrity of the fruiting body.
* Euler characteristic discriminates pattern type (spots vs. maze).

By quantifying these metrics the pipeline can:

1. Validate simulation accuracy against published morphological data.
2. Compare species objectively in a shared morphospace.
3. Identify morphological constraints imposed by reaction-diffusion
   dynamics.
4. Detect pattern-type transitions as parameters vary.

Author: Mario Ahumada Duran
Date: December 2025
License: MIT
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from scipy.fft import fft2, fftfreq
from scipy.spatial import distance, cKDTree
from skimage import measure, morphology, filters, feature
import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from datetime import datetime
import pandas as pd
from tqdm import tqdm

# Import logging helper from the packaged simulator (kept lightweight to avoid hard dependency cycles)
try:
    from ..core.turing_simulator import setup_logging
except Exception:
    # Minimal fallback (allows running this module in isolation in constrained environments)
    def setup_logging(log_level="INFO", log_file=None):
        import logging
        logging.basicConfig(
            level=getattr(logging, log_level),
            format='[%(asctime)s] %(levelname)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

# =============================================================================
# TOPOLOGY ANALYZER CLASS
# =============================================================================

class TopologyAnalyzer:
    """Comprehensive morphometric and topological analysis of Turing patterns.

    Provides a suite of methods to quantify spatial patterns and extract
    biologically relevant metrics for comparison with empirical data.
    Designed for the FungalMorphoSpace pipeline but applicable to any
    2D scalar field exhibiting periodic or quasi-periodic structure.

    Attributes:
        pattern: Original 2D array representing the pattern (activator
            field), stored unmodified.
        dx: Spatial resolution factor.  When ``dx=1.0`` (default), all
            spatial metrics are reported in pixel units.  Otherwise, *dx*
            is interpreted as micrometers per pixel and all lengths are
            converted accordingly.
        metrics: Dictionary of computed scalar metrics, progressively
            populated as individual analysis methods are called.  Fully
            JSON-serializable (no NumPy arrays).
        pattern_normalized: Pattern rescaled to ``[0, 1]`` (min-max
            normalization), used for binarization and display.
    """

    def __init__(self, pattern: np.ndarray, dx: float = 1.0) -> None:
        """Initialize the topology analyzer.

        Args:
            pattern: 2D NumPy array to analyze.  Typically the activator
                concentration field *u(x, y)* from a reaction-diffusion
                simulation.  Must be at least 3x3.
            dx: Spatial resolution.  If ``1.0`` (default), all spatial
                metrics are in pixels.  Otherwise, *dx* represents the
                micrometer-per-pixel conversion factor and all reported
                lengths carry physical units.
        """
        self.pattern = pattern
        self.dx = dx
        self.metrics = {}
        
        # Normalize pattern to [0, 1] for consistent analysis
        self.pattern_normalized = (pattern - pattern.min()) / (pattern.max() - pattern.min())
        
        # Determine unit label based on dx
        unit = "px" if dx == 1.0 else "μm"
        
        logging.info(f"TopologyAnalyzer initialized:")
        logging.info(f"  Pattern size: {pattern.shape}")
        logging.info(f"  Spatial resolution: dx={dx} ({'pixels' if dx == 1.0 else 'μm/pixel'})")
        logging.info(f"  Physical size: {pattern.shape[0]*dx:.1f} × {pattern.shape[1]*dx:.1f} {unit}²")
    
    # -------------------------------------------------------------------------
    # Internal helpers
    # -------------------------------------------------------------------------
    def _radial_average_pixels(self, image: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Compute the radial average of a 2D image around its center.

        Each pixel is assigned an integer radius (Euclidean distance from the
        image center, truncated to ``int``).  The mean value at each integer
        radius is then computed via ``np.bincount``, which is O(N) and avoids
        an explicit loop over radial bins.

        Args:
            image: 2D array to be radially averaged (e.g., a power spectrum
                or autocorrelation map).

        Returns:
            A tuple ``(r_px, profile)`` where:

            * ``r_px`` -- 1D array of integer radii ``[0, 1, ..., max_r]``
              in pixel units.
            * ``profile`` -- 1D array of the same length giving the mean
              image value at each integer radius.  Radii with zero pixel
              count are set to 0.0.
        """
        ny, nx = image.shape
        cy, cx = ny // 2, nx // 2

        # Build coordinate grids and compute Euclidean distance from center
        y, x = np.indices(image.shape)
        r = np.sqrt((x - cx) ** 2 + (y - cy) ** 2)
        r_int = r.astype(np.int32)  # Truncate to integer radial bin index

        max_r = int(r_int.max())

        # bincount accumulates the sum and count of pixels in each radial bin
        # in a single vectorized pass -- much faster than a Python loop.
        sums = np.bincount(r_int.ravel(), weights=image.ravel(), minlength=max_r + 1)
        counts = np.bincount(r_int.ravel(), minlength=max_r + 1)

        # Safe division: bins with no pixels get 0.0 instead of NaN
        profile = np.divide(sums, counts, out=np.zeros_like(sums, dtype=float), where=counts > 0)

        r_px = np.arange(max_r + 1, dtype=float)
        return r_px, profile

    @staticmethod
    def _smooth_1d(arr: np.ndarray, sigma: float = 2.0) -> np.ndarray:
        """Apply 1D Gaussian smoothing to a profile array.

        Used as a pre-processing step before peak detection on radial power
        spectra and autocorrelation profiles to suppress high-frequency noise
        that would otherwise produce spurious local extrema.

        Args:
            arr: 1D array to smooth (e.g., radial power profile).
            sigma: Standard deviation of the Gaussian kernel in array-index
                units.  Default of 2.0 provides moderate smoothing that
                preserves broad spectral peaks while removing pixel-scale
                noise.

        Returns:
            Smoothed copy of *arr* (same length), cast to ``float64``.
        """
        return ndimage.gaussian_filter1d(arr.astype(float), sigma=sigma)

    @staticmethod
    def _local_extrema_indices(arr: np.ndarray, kind: str = "max") -> np.ndarray:
        """Detect indices of local maxima or minima in a 1D array.

        A local maximum at index *i* satisfies
        ``arr[i] > arr[i-1]`` **and** ``arr[i] > arr[i+1]`` (strict
        inequality); analogously for minima.  Boundary elements (first and
        last) are never returned.

        This simple three-point comparison is intentionally used instead of
        ``scipy.signal.find_peaks`` to avoid introducing an extra dependency
        and to keep full control over the detection logic (no prominence /
        distance filtering is applied here -- callers handle that).

        Args:
            arr: 1D numeric array (typically a smoothed radial profile).
            kind: ``"max"`` for local maxima, ``"min"`` for local minima.

        Returns:
            Sorted 1D integer array of indices where local extrema occur.
            Empty array if ``arr`` has fewer than 3 elements.

        Raises:
            ValueError: If *kind* is not ``"max"`` or ``"min"``.
        """
        if arr.size < 3:
            return np.array([], dtype=int)
        if kind == "max":
            # Strict inequality: interior point must exceed both neighbors
            mask = (arr[1:-1] > arr[:-2]) & (arr[1:-1] > arr[2:])
        elif kind == "min":
            mask = (arr[1:-1] < arr[:-2]) & (arr[1:-1] < arr[2:])
        else:
            raise ValueError("kind must be 'max' or 'min'")
        # +1 offset because mask is indexed relative to arr[1:-1]
        return np.where(mask)[0] + 1

    def measure_wavelength_fft(self) -> float:
        """Measure characteristic wavelength via radially-averaged FFT power spectrum.

        This replaces the older "global max in 2D power" heuristic, which
        can fail when DC / very-low-frequency structure dominates (gradients,
        boundary effects), the pattern is weakly periodic, or noise produces
        many comparable local maxima.

        Algorithm:
            1. Subtract the spatial mean to remove the DC component.
            2. Compute the 2D FFT and its power spectrum ``|F(k)|^2``.
            3. Build a 1D radial power profile by binning pixels in
               frequency space into 200 uniform annuli of ``|k|`` and
               averaging the power within each annulus (bincount method).
            4. Smooth the radial profile (Gaussian, sigma=2 bins) and find
               the dominant peak at frequencies above a physically motivated
               floor (>= 2 full cycles across the smallest domain dimension).
            5. Convert to wavelength: ``lambda = 1 / f_peak``.

        Quality control (QC):
            The peak must satisfy **two** criteria to pass QC:

            * **Prominence ratio >= 3.0** -- the peak power must be at
              least 3x the median background power across all valid
              frequency bins.  This rejects patterns with no clear
              periodicity.
            * **Wavelength <= 0.5 * domain** -- rejects spurious
              domain-scale modes that are artifacts of finite-size
              effects rather than genuine Turing instabilities.

        Returns:
            Dominant wavelength in the spatial units defined by ``self.dx``
            (pixels when ``dx=1.0``, micrometers otherwise).  Returns
            ``np.nan`` if no valid peak is found or QC fails.

        Side Effects:
            Populates ``self.metrics`` with keys:

            * ``wavelength_fft`` -- the measured wavelength (or NaN).
            * ``wavelength_fft_qc_pass`` -- boolean QC flag.
            * ``wavelength_fft_peak_ratio`` -- prominence ratio.
            * ``wavelength_fft_peak_freq`` -- frequency of dominant peak.

            Also caches ``_fft_freq_centers``, ``_fft_radial_power``, and
            ``_fft_radial_power_smooth`` for downstream visualization.
        """
        # Step 1: Remove spatial mean (zeroes the DC / k=0 bin)
        pattern_centered = self.pattern - np.mean(self.pattern)

        # Step 2: 2D FFT and power spectrum |F(k)|^2
        fft = fft2(pattern_centered)
        power = np.abs(fft) ** 2

        ny, nx = self.pattern.shape
        freqs_x = fftfreq(nx, d=self.dx)  # Spatial frequencies along x
        freqs_y = fftfreq(ny, d=self.dx)  # Spatial frequencies along y
        fx, fy = np.meshgrid(freqs_x, freqs_y)
        freq_radial = np.sqrt(fx ** 2 + fy ** 2)  # Radial frequency |k|

        # Explicitly zero DC to prevent it from dominating the radial average
        power[0, 0] = 0.0

        # Step 3: Radial binning in frequency space (200 uniform annuli)
        n_bins = 200
        freq_bins = np.linspace(0.0, float(freq_radial.max()), n_bins + 1)
        bin_idx = np.digitize(freq_radial.ravel(), freq_bins) - 1
        bin_idx = np.clip(bin_idx, 0, n_bins - 1)

        # bincount gives O(N) accumulation of power per radial bin
        sums = np.bincount(bin_idx, weights=power.ravel(), minlength=n_bins).astype(float)
        counts = np.bincount(bin_idx, minlength=n_bins).astype(float)
        radial_power = np.divide(sums, counts, out=np.zeros_like(sums), where=counts > 0)

        freq_centers = 0.5 * (freq_bins[:-1] + freq_bins[1:])
        radial_power_s = self._smooth_1d(radial_power, sigma=2.0)

        # Step 4: Minimum frequency floor -- require at least ~2 full
        # periods across the smallest domain dimension to reject
        # domain-scale artifacts.
        domain = min(nx, ny) * self.dx
        f_min = 2.0 / max(domain, 1e-12)

        valid = (freq_centers >= f_min) & (counts > 0) & np.isfinite(radial_power_s)
        if np.count_nonzero(valid) < 5:
            wavelength = np.nan
            qc_pass = False
            peak_freq = np.nan
            peak_ratio = np.nan
            logging.warning("FFT wavelength: insufficient valid frequency bins after filtering")
        else:
            # Locate global maximum in the valid portion of the smoothed profile
            idx_valid = np.where(valid)[0]
            peak_i = idx_valid[int(np.argmax(radial_power_s[valid]))]
            peak_freq = float(freq_centers[peak_i])

            # Step 5: Convert dominant frequency to wavelength
            wavelength = (1.0 / peak_freq) if peak_freq > 0 else np.nan

            # QC criterion 1: prominence ratio (peak / median background)
            background = float(np.median(radial_power_s[valid]))
            peak_val = float(radial_power_s[peak_i])
            peak_ratio = peak_val / (background + 1e-12)

            # QC criterion 2: reject wavelengths that span > half the domain
            max_reasonable = 0.5 * domain
            qc_pass = bool(np.isfinite(wavelength) and (wavelength <= max_reasonable) and (peak_ratio >= 3.0))

            if not qc_pass:
                logging.warning(
                    f"FFT wavelength QC failed: λ={wavelength:.2f}, peak_ratio={peak_ratio:.2f}, domain={domain:.2f}"
                )

        # Store metrics for JSON serialization
        self.metrics['wavelength_fft'] = float(wavelength) if np.isfinite(wavelength) else np.nan
        self.metrics['wavelength_fft_qc_pass'] = bool(qc_pass)
        self.metrics['wavelength_fft_peak_ratio'] = float(peak_ratio) if np.isfinite(peak_ratio) else np.nan
        self.metrics['wavelength_fft_peak_freq'] = float(peak_freq) if np.isfinite(peak_freq) else np.nan

        # Cache profiles for debugging/visualization (not JSON-serialized)
        self._fft_freq_centers = freq_centers
        self._fft_radial_power = radial_power
        self._fft_radial_power_smooth = radial_power_s

        unit = "px" if self.dx == 1.0 else "μm"
        logging.info(f"Wavelength (FFT): {wavelength:.2f} {unit} (QC={'PASS' if qc_pass else 'FAIL'})")

        return float(wavelength)

    def measure_wavelength_autocorr(self) -> float:
        """Measure characteristic wavelength via radially-averaged 2D autocorrelation.

        Uses the **Wiener-Khinchin theorem**: the autocorrelation of a
        real-valued spatial signal equals the inverse FFT of its power
        spectrum, i.e. ``C(r) = IFFT(|FFT(u)|^2)``.  This avoids an
        expensive direct spatial convolution.

        Compared to the earlier "first local maximum" heuristic, this
        version is more robust because:

        * Radial averaging is computed on integer pixel radii via
          ``_radial_average_pixels`` (bincount, O(N)).
        * The characteristic wavelength is defined as the **first
          maximum after the first minimum** of the smoothed radial
          autocorrelation profile.  For an oscillatory correlation
          function ``C(r)``, this corresponds to one full period of the
          dominant spatial frequency and is the standard definition used
          in soft-matter physics (e.g., block-copolymer morphology
          analysis).

        Returns:
            Dominant wavelength in the spatial units defined by ``self.dx``
            (pixels when ``dx=1.0``, micrometers otherwise).  Returns
            ``np.nan`` if no oscillation is detected.

        Side Effects:
            Populates ``self.metrics`` with keys:

            * ``wavelength_autocorr`` -- the measured wavelength (or NaN).
            * ``wavelength_autocorr_qc_pass`` -- boolean indicating
              whether a valid oscillation peak was found.

            Also caches ``_autocorr_image``, ``_autocorr_r_px``,
            ``_autocorr_profile``, and ``_autocorr_profile_smooth`` for
            downstream visualization.
        """
        pattern_centered = self.pattern - np.mean(self.pattern)

        # Wiener-Khinchin theorem: C(r) = IFFT( |FFT(u)|^2 )
        # This computes the full 2D autocorrelation in O(N log N) time.
        fft = fft2(pattern_centered)
        autocorr = np.real(np.fft.ifft2(fft * np.conj(fft)))
        autocorr = np.fft.fftshift(autocorr)  # Center the zero-lag peak

        # Normalize to [-1, 1] range (zero-lag = 1.0)
        maxv = float(np.max(np.abs(autocorr)))
        if maxv > 0:
            autocorr = autocorr / maxv

        # Compute radial profile (in pixels), then convert to physical units
        r_px, profile = self._radial_average_pixels(autocorr)
        profile_s = self._smooth_1d(profile, sigma=2.0)

        # Skip the origin neighborhood (r < 3 px) to avoid the trivial
        # zero-lag peak dominating the extrema search.
        skip = 3
        mins = self._local_extrema_indices(profile_s, kind="min")
        mins = mins[mins >= skip]
        # The first minimum marks the anti-correlation trough
        start_idx = int(mins[0]) if mins.size > 0 else skip

        # First maximum *after* the first minimum = one full wavelength
        maxs = self._local_extrema_indices(profile_s, kind="max")
        maxs = maxs[maxs > start_idx]
        if maxs.size > 0:
            peak_idx = int(maxs[0])
            wavelength = float(r_px[peak_idx] * self.dx)  # Convert px -> physical
            qc_pass = True
        else:
            wavelength = np.nan
            qc_pass = False
            logging.warning("Autocorr wavelength: no suitable peak found (after first minimum).")

        self.metrics['wavelength_autocorr'] = float(wavelength) if np.isfinite(wavelength) else np.nan
        self.metrics['wavelength_autocorr_qc_pass'] = bool(qc_pass)

        # Cache profiles for debugging/visualization (not JSON-serialized)
        self._autocorr_image = autocorr
        self._autocorr_r_px = r_px
        self._autocorr_profile = profile
        self._autocorr_profile_smooth = profile_s

        unit = "px" if self.dx == 1.0 else "μm"
        logging.info(f"Wavelength (autocorrelation): {wavelength:.2f} {unit} (QC={'PASS' if qc_pass else 'FAIL'})")

        return float(wavelength)

    def compute_euler_characteristic(self, threshold: float = 0.5) -> int:
        """Compute the Euler characteristic of the binarized pattern.

        The Euler characteristic chi is a topological invariant defined as:

            chi = #connected_components - #holes

        For a 2D binary image this equals ``V - E + F`` in the CW-complex
        sense, but ``skimage.measure.euler_number`` computes it efficiently
        using the quad-tree bit-pattern method of Gray (1971), which
        correctly handles edge pixels without padding artifacts.

        Args:
            threshold: Binarization threshold applied to the ``[0, 1]``
                normalized pattern.  Default is 0.5 (Otsu-like midpoint).

        Returns:
            Integer Euler characteristic chi.

        Side Effects:
            Populates ``self.metrics`` with keys ``euler_characteristic``,
            ``n_components``, and ``n_holes``.  Also caches the binary
            mask (``_binary``), the labeled-component array (``_labeled``),
            and the threshold value (``_binary_threshold``) for downstream
            use by ``measure_regularity`` and ``visualize_analysis``.

        Notes:
            Different pattern types have characteristic chi values:

            * **Spots**: chi > 0 (many isolated components, few holes).
            * **Labyrinths**: chi ~ 0 (single connected network).
            * **Inverse spots / holes**: chi < 0 (connected background
              with many enclosed voids).

            **Connectivity = 2 (8-connectivity)** is used deliberately.
            In 2D, ``connectivity=1`` (4-connectivity) treats diagonal
            neighbors as disconnected, which would artificially fragment
            labyrinthine Turing patterns along diagonal ridges.  Using
            ``connectivity=2`` ensures that diagonally touching pixels
            belong to the same connected component, matching the
            biological interpretation where continuous tissue ridges are
            structurally connected regardless of orientation.
        """
        # Binarize pattern
        binary = self.pattern_normalized > threshold

        # Cache for downstream metrics/visualization
        self._binary = binary
        self._binary_threshold = float(threshold)
        
        # Use skimage's robust euler_number calculation
        # connectivity=2 means 8-connectivity in 2D (diagonal neighbors count)
        chi = measure.euler_number(binary, connectivity=2)
        
        # Also count components for reporting (still useful metric)
        labeled = measure.label(binary, connectivity=2)
        n_components = labeled.max()

        # Cache labeled components for spacing/regularity
        self._labeled = labeled
        
        # Derive holes from Euler formula: χ = components - holes
        # Note: This is now consistent with the euler_number calculation
        n_holes = max(0, n_components - chi)
        
        self.metrics['euler_characteristic'] = int(chi)
        self.metrics['n_components'] = int(n_components)
        self.metrics['n_holes'] = int(n_holes)
        
        logging.info(f"Euler characteristic: χ = {chi}")
        logging.info(f"  Components: {n_components}, Holes: {n_holes}")
        
        return chi
    
    def measure_regularity(self) -> Dict[str, float]:
        """Quantify the spatial regularity and uniformity of the pattern.

        This is a central metric for evaluating developmental precision in
        biological Turing patterns: a perfectly regular hexagonal lattice
        has CV(spacing) = 0, while a random Poisson point pattern has
        CV ~ 0.5.

        Design rationale (key change from earlier versions):
            The previous implementation detected **all local maxima** in the
            raw concentration field via ``maximum_filter`` equality, which
            wildly overcounts in noisy or weakly periodic patterns (e.g.,
            4000--5000 spurious peaks on a 256x256 grid).  This made the
            spacing statistics meaningless.

            The current implementation uses **connected-component centroids**
            from the binarized pattern (the same binary mask used for the
            Euler characteristic), which:

            * Aligns ``n_peaks`` with the biologically interpretable
              component count.
            * Avoids noise-driven overcounting.
            * Produces stable nearest-neighbor spacing estimates.

        Algorithm:
            1. Extract centroids of labeled connected components.
            2. Build a ``scipy.spatial.cKDTree`` from the centroid
               coordinates (O(n log n) construction).
            3. Query the tree for k=2 nearest neighbors (k=1 is self),
               yielding one nearest-neighbor distance per component.
            4. Compute summary statistics: mean, std, CV, and Shannon
               entropy of the NN distance distribution.

        Returns:
            Dictionary with keys:

            * ``mean_spacing`` -- mean nearest-neighbor distance.
            * ``std_spacing`` -- standard deviation of NN distances.
            * ``cv_spacing`` -- coefficient of variation (std / mean).
              Values < 0.2 indicate high regularity.
            * ``n_peaks`` -- number of connected components used.
            * ``entropy`` -- Shannon entropy of the 20-bin NN distance
              histogram (nats).  Higher values indicate broader, less
              regular distributions.
            * ``binarization_threshold`` -- the threshold used for the
              underlying binary mask.

        Side Effects:
            Updates ``self.metrics`` with the returned dictionary.  Caches
            ``_peak_coords`` (Nx2 array of centroid row/col positions) and
            ``_nn_distances`` (1D array of NN distances in physical units)
            for visualization.
        """
        # Ensure we have a cached binary/labels (computed in compute_euler_characteristic)
        if not hasattr(self, "_labeled") or self._labeled is None:
            # Fallback (should be rare): compute with default threshold
            self.compute_euler_characteristic(threshold=0.5)

        labeled = self._labeled
        n_components = int(labeled.max())

        if n_components < 2:
            logging.warning("Regularity: insufficient components for spacing analysis")
            regularity_metrics = {
                'mean_spacing': np.nan,
                'std_spacing': np.nan,
                'cv_spacing': np.nan,
                'n_peaks': n_components,
                'entropy': np.nan
            }
            self.metrics.update(regularity_metrics)
            return regularity_metrics

        # Extract centroid (row, col) of each connected component
        props = measure.regionprops(labeled)
        peak_coords = np.array([p.centroid for p in props], dtype=float)

        # Build a KD-tree for O(n log n) nearest-neighbor lookup instead
        # of the naive O(n^2) pairwise distance matrix (scipy cdist).
        tree = cKDTree(peak_coords)
        # k=2 because k=1 always returns self (distance 0); k=2 gives the
        # true nearest other component.
        dists, _ = tree.query(peak_coords, k=2)
        nn_distances = dists[:, 1] * self.dx  # Convert pixel distances to physical units

        # Coefficient of variation: CV = sigma / mu
        mean_d = float(np.mean(nn_distances))
        std_d = float(np.std(nn_distances))
        cv = float(std_d / (mean_d + 1e-12))  # Epsilon prevents division by zero

        # Shannon entropy of the NN distance distribution.
        # A 20-bin histogram with density=True gives a probability density;
        # entropy quantifies how spread out the distribution is.
        hist, _ = np.histogram(nn_distances, bins=20, density=True)
        hist = hist[hist > 0]  # Drop empty bins before log
        entropy = float(-np.sum(hist * np.log(hist + 1e-12)))

        regularity_metrics = {
            'mean_spacing': mean_d,
            'std_spacing': std_d,
            'cv_spacing': cv,
            'n_peaks': int(len(peak_coords)),
            'entropy': entropy,
            'binarization_threshold': float(getattr(self, "_binary_threshold", 0.5))
        }

        self.metrics.update(regularity_metrics)

        # Cache for visualization/debug
        self._peak_coords = peak_coords
        self._nn_distances = nn_distances

        unit = "px" if self.dx == 1.0 else "μm"
        logging.info("Regularity metrics (component centroids):")
        logging.info(f"  Peaks/components: {regularity_metrics['n_peaks']}")
        logging.info(f"  Mean spacing: {mean_d:.2f} {unit}")
        logging.info(f"  CV spacing: {cv:.3f}")
        logging.info(f"  Entropy: {entropy:.3f}")

        return regularity_metrics

    def measure_isotropy(self) -> float:
        """Measure directional bias (isotropy) of the spatial pattern.

        Isotropy quantifies whether the pattern has a preferred orientation.
        Perfectly isotropic patterns (e.g., hexagonal spot arrays) distribute
        their spectral power uniformly across all directions, while stripe
        patterns concentrate power along a single axis.

        Algorithm:
            1. Compute the 2D FFT power spectrum and shift zero-frequency
               to the center.
            2. For each pixel in the shifted spectrum, compute the polar
               angle ``theta = arctan2(dy, dx)`` relative to the center.
            3. Partition ``[-pi, +pi]`` into 36 angular bins (10 degrees
               each) and average the power within each wedge.
            4. Normalize the angular power distribution to ``[0, 1]``.
            5. Compute the isotropy index as
               ``I = 1 - CV(angular_power)`` where
               ``CV = std / mean``.  A uniform distribution (all bins
               equal) gives ``I = 1``; a single dominant direction gives
               ``I -> 0``.

        Returns:
            Isotropy index in ``[0, 1]``.  Values near 1 indicate no
            preferred direction; values near 0 indicate strong
            anisotropy (e.g., parallel stripes).

        Side Effects:
            Stores ``isotropy_index`` in ``self.metrics``.
        """
        # Compute centered 2D FFT power spectrum
        pattern_centered = self.pattern - np.mean(self.pattern)
        fft = fft2(pattern_centered)
        power = np.abs(np.fft.fftshift(fft))**2

        # Center pixel of the shifted spectrum
        cy, cx = np.array(power.shape) // 2

        # Polar angle of each frequency-domain pixel relative to center
        y, x = np.indices(power.shape)
        dy = y - cy
        dx = x - cx
        angle = np.arctan2(dy, dx)

        # 36 angular bins spanning [-pi, +pi] => 10-degree resolution
        n_bins = 36
        angle_bins = np.linspace(-np.pi, np.pi, n_bins + 1)
        angular_power = []

        for i in range(n_bins):
            mask = (angle >= angle_bins[i]) & (angle < angle_bins[i+1])
            if np.any(mask):
                angular_power.append(np.mean(power[mask]))
            else:
                angular_power.append(0)

        angular_power = np.array(angular_power)

        # Normalize to [0, 1] for a scale-invariant CV calculation
        if angular_power.max() > 0:
            angular_power /= angular_power.max()

        # Isotropy index: I = 1 - CV(angular_power)
        # Perfect isotropy => all bins equal => std = 0 => I = 1
        # Strong anisotropy => one dominant bin => high CV => I -> 0
        isotropy_index = 1.0 - (np.std(angular_power) / (np.mean(angular_power) + 1e-12))
        isotropy_index = np.clip(isotropy_index, 0, 1)

        self.metrics['isotropy_index'] = float(isotropy_index)

        logging.info(f"Isotropy index: {isotropy_index:.3f} (1 = perfectly isotropic)")

        return isotropy_index
    
    def compute_all_metrics(self) -> Dict[str, object]:
        """Orchestrate computation of all morphometric and topological metrics.

        Calls each analysis method in a fixed order that respects internal
        dependencies (``compute_euler_characteristic`` must run before
        ``measure_regularity`` because the latter re-uses the labeled
        connected-component array).

        Execution order:
            1. ``measure_wavelength_fft``
            2. ``measure_wavelength_autocorr``
            3. ``compute_euler_characteristic`` (produces ``_labeled``)
            4. ``measure_regularity`` (consumes ``_labeled``)
            5. ``measure_isotropy``

        After all sub-analyses complete, pattern metadata (shape, dx,
        ISO-8601 timestamp) is appended.

        Returns:
            Reference to ``self.metrics`` -- a dictionary containing all
            computed scalar metrics.  Array-valued diagnostics (radial
            profiles, images) are stored as private attributes but are
            **not** included in this dict to keep it JSON-serializable.
        """
        logging.info("="*70)
        logging.info("COMPUTING ALL METRICS")
        logging.info("="*70)
        
        # Wavelength
        self.measure_wavelength_fft()
        self.measure_wavelength_autocorr()
        
        # Topology
        self.compute_euler_characteristic()
        
        # Regularity
        self.measure_regularity()
        
        # Isotropy
        self.measure_isotropy()
        
        # Add metadata
        self.metrics['pattern_shape'] = self.pattern.shape
        self.metrics['dx'] = self.dx
        self.metrics['timestamp'] = datetime.now().isoformat()
        
        logging.info("="*70)
        logging.info("METRICS COMPUTATION COMPLETE")
        logging.info("="*70)
        
        return self.metrics
    
    def save_metrics(self, output_path: str) -> None:
        """Serialize computed metrics to a JSON file.

        Parent directories are created automatically if they do not exist.
        The output is pretty-printed with 2-space indentation for human
        readability.

        Args:
            output_path: Filesystem path for the output ``.json`` file.
                Can be absolute or relative.
        """
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w') as f:
            json.dump(self.metrics, f, indent=2)
        
        logging.info(f"Metrics saved to {output_path}")
    
    def visualize_analysis(self, save_path: Optional[str] = None, show: bool = True, suptitle: Optional[str] = None) -> None:
        """Create a comprehensive 3x3 diagnostic figure of all analyses.

        The figure is arranged as a 3-row, 3-column grid:

        +------------------+------------------+---------------------+
        | (0,0) Original   | (0,1) Binarized  | (0,2) Log power     |
        | pattern (viridis)| mask (gray) + chi| spectrum (inferno)  |
        +------------------+------------------+---------------------+
        | (1,0) 2D auto-   | (1,1) Component  | (1,2) Metrics       |
        | correlation      | centroids on     | summary text box    |
        | (RdBu_r)         | pattern (red pts)|                     |
        +------------------+------------------+---------------------+
        | (2,0) NN spacing | (2,1) Radial     | (2,2) [empty]       |
        | histogram        | power spectrum   |                     |
        +------------------+------------------+---------------------+

        Each panel is self-contained and labeled.  If a particular
        analysis has not been run yet, the corresponding panel shows a
        placeholder message instead of raising an error.

        Args:
            save_path: If provided, the figure is saved as a 300 dpi PNG
                at this path.  Parent directories are created if needed.
            show: If ``True`` (default), calls ``plt.show()``; otherwise
                closes the figure after (optional) saving.  Set to
                ``False`` in batch / headless mode.
            suptitle: Optional top-level title string placed above the
                grid.  If ``None``, a default title is used.
        """
        fig = plt.figure(figsize=(15, 10))
        gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

        if suptitle:
            fig.suptitle(suptitle, fontsize=12, fontweight='bold', y=0.98)
        
        # 1. Original pattern
        ax1 = fig.add_subplot(gs[0, 0])
        im1 = ax1.imshow(self.pattern, cmap='viridis', origin='lower')
        ax1.set_title('Original Pattern', fontweight='bold')
        plt.colorbar(im1, ax=ax1, fraction=0.046)
        
        # 2. Binary threshold
        ax2 = fig.add_subplot(gs[0, 1])
        binary = self.pattern_normalized > 0.5
        ax2.imshow(binary, cmap='gray', origin='lower')
        ax2.set_title(f'Binarized (χ = {self.metrics.get("euler_characteristic", "?")})', 
                     fontweight='bold')
        
        # 3. Power spectrum
        ax3 = fig.add_subplot(gs[0, 2])
        fft = fft2(self.pattern - np.mean(self.pattern))
        power = np.log10(np.abs(np.fft.fftshift(fft))**2 + 1)
        im3 = ax3.imshow(power, cmap='inferno', origin='lower')
        ax3.set_title('Power Spectrum (log)', fontweight='bold')
        plt.colorbar(im3, ax=ax3, fraction=0.046)
        
        # 4. Autocorrelation
        ax4 = fig.add_subplot(gs[1, 0])
        pattern_centered = self.pattern - np.mean(self.pattern)
        fft_ac = fft2(pattern_centered)
        autocorr = np.real(np.fft.ifft2(fft_ac * np.conj(fft_ac)))
        autocorr = np.fft.fftshift(autocorr)
        autocorr /= autocorr.max()
        im4 = ax4.imshow(autocorr, cmap='RdBu_r', origin='lower', vmin=-0.5, vmax=1)
        ax4.set_title('Autocorrelation', fontweight='bold')
        plt.colorbar(im4, ax=ax4, fraction=0.046)
        
        # 5. Feature points (component centroids)
        ax5 = fig.add_subplot(gs[1, 1])
        ax5.imshow(self.pattern, cmap='gray', origin='lower', alpha=0.5)

        # Use cached centroids from measure_regularity(), fall back to components if needed
        peak_coords = getattr(self, "_peak_coords", None)
        if peak_coords is None:
            # Ensure components exist
            if not hasattr(self, "_labeled") or self._labeled is None:
                self.compute_euler_characteristic(threshold=0.5)
            props = measure.regionprops(self._labeled)
            peak_coords = np.array([p.centroid for p in props], dtype=float)

        if peak_coords.size > 0:
            ax5.scatter(peak_coords[:, 1], peak_coords[:, 0], c='red', s=12, alpha=0.7)

        ax5.set_title(f'Peak Detection (n={len(peak_coords)})', fontweight='bold')
        
        # 6. Metrics summary
        ax6 = fig.add_subplot(gs[1, 2])
        ax6.axis('off')
        metrics_text = "MORPHOMETRIC METRICS\n" + "="*25 + "\n\n"
        
        # Determine unit label based on dx
        # If dx=1.0, we're working in pixels; otherwise assume dx is μm/pixel
        unit = "px" if self.dx == 1.0 else "μm"
        
        if 'wavelength_fft' in self.metrics:
            metrics_text += f"λ (FFT): {self.metrics['wavelength_fft']:.2f} {unit}\n"
        if 'wavelength_autocorr' in self.metrics:
            metrics_text += f"λ (autocorr): {self.metrics['wavelength_autocorr']:.2f} {unit}\n"
        metrics_text += "\n"
        
        if 'euler_characteristic' in self.metrics:
            metrics_text += f"Euler char (χ): {self.metrics['euler_characteristic']}\n"
            metrics_text += f"Components: {self.metrics['n_components']}\n"
            metrics_text += f"Holes: {self.metrics['n_holes']}\n"
        metrics_text += "\n"
        
        if 'mean_spacing' in self.metrics:
            metrics_text += f"Mean spacing: {self.metrics['mean_spacing']:.2f} {unit}\n"
            metrics_text += f"CV spacing: {self.metrics['cv_spacing']:.3f}\n"
        
        if 'isotropy_index' in self.metrics:
            metrics_text += f"Isotropy: {self.metrics['isotropy_index']:.3f}\n"
        
        ax6.text(0.1, 0.9, metrics_text, transform=ax6.transAxes, 
                fontsize=10, verticalalignment='top', family='monospace',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        # 7. Spacing histogram
        if 'mean_spacing' in self.metrics:
            ax7 = fig.add_subplot(gs[2, 0])

            nn_distances = getattr(self, "_nn_distances", None)
            if nn_distances is None:
                # Attempt to reconstruct by re-running regularity analysis
                try:
                    self.measure_regularity()
                    nn_distances = getattr(self, "_nn_distances", None)
                except Exception:
                    nn_distances = None

            if nn_distances is not None and len(nn_distances) > 0 and np.all(np.isfinite(nn_distances)):
                ax7.hist(nn_distances, bins=30, edgecolor='black', alpha=0.7)
                ax7.axvline(self.metrics.get('mean_spacing', np.nan), color='red',
                            linestyle='--', linewidth=2, label='Mean')
                spacing_unit = "px" if self.dx == 1.0 else "μm"
                ax7.set_xlabel(f'Spacing ({spacing_unit})')
                ax7.set_ylabel('Frequency')
                ax7.set_title('Spacing Distribution', fontweight='bold')
                ax7.legend()
            else:
                ax7.axis('off')
                ax7.text(0.1, 0.5, "No spacing data available", transform=ax7.transAxes)
# 8. Radial power spectrum
        ax8 = fig.add_subplot(gs[2, 1])

        # Prefer cached FFT radial profile (computed in measure_wavelength_fft)
        freq_centers = getattr(self, "_fft_freq_centers", None)
        radial_power_s = getattr(self, "_fft_radial_power_smooth", None)

        if freq_centers is None or radial_power_s is None:
            try:
                self.measure_wavelength_fft()
                freq_centers = getattr(self, "_fft_freq_centers", None)
                radial_power_s = getattr(self, "_fft_radial_power_smooth", None)
            except Exception:
                freq_centers, radial_power_s = None, None

        if freq_centers is not None and radial_power_s is not None:
            ax8.plot(freq_centers, radial_power_s, 'b-', linewidth=2)
            freq_unit = "1/px" if self.dx == 1.0 else "1/μm"
            ax8.set_xlabel(f'Spatial Frequency ({freq_unit})')
            ax8.set_ylabel('Power (radial mean)')
            ax8.set_title('Radial Power Spectrum', fontweight='bold')
            ax8.set_yscale('log')
            ax8.grid(True, alpha=0.3)
        else:
            ax8.axis('off')
            ax8.text(0.1, 0.5, "No FFT profile available", transform=ax8.transAxes)
# Overall title
        fig.suptitle('Topological and Morphometric Analysis', 
                    fontsize=16, fontweight='bold')
        
        if save_path:
            Path(save_path).parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logging.info(f"Analysis visualization saved to {save_path}")
        
        if show:
            plt.show()
        else:
            plt.close()


# =============================================================================
# BATCH ANALYSIS
# =============================================================================

def analyze_directory(input_dir: str, output_dir: str, dx: float = 1.0) -> None:
    """Batch-analyze all Turing-pattern arrays in a directory.

    Scans *input_dir* for ``.npy`` files, runs the full
    ``TopologyAnalyzer.compute_all_metrics`` pipeline on each, and
    produces three artifact types per pattern:

    * ``<stem>_metrics.json`` -- scalar metrics dictionary.
    * ``<stem>_analysis.png`` -- 3x3 diagnostic figure.

    After all patterns are processed, a summary CSV
    (``batch_metrics_summary.csv``) is written to *output_dir* with one
    row per pattern and columns for every metric, facilitating
    downstream statistical analysis in R, pandas, or Excel.

    Args:
        input_dir: Path to directory containing ``.npy`` pattern files.
            Each file must be a 2D NumPy array.
        output_dir: Path to directory where results are written.  Created
            automatically if it does not exist.
        dx: Spatial resolution (micrometers per pixel).  Passed to
            ``TopologyAnalyzer.__init__``.  Use ``1.0`` for pixel-unit
            analysis.
    """
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Find all .npy files
    pattern_files = list(input_path.glob("*.npy"))
    
    if not pattern_files:
        logging.error(f"No .npy files found in {input_dir}")
        return
    
    logging.info(f"Found {len(pattern_files)} patterns to analyze")
    
    all_metrics = []
    
    for pattern_file in tqdm(pattern_files, desc="Analyzing patterns"):
        # Load pattern
        pattern = np.load(pattern_file)
        
        # Analyze
        analyzer = TopologyAnalyzer(pattern, dx=dx)
        metrics = analyzer.compute_all_metrics()
        
        # Add file info
        metrics['filename'] = pattern_file.name
        all_metrics.append(metrics)
        
        # Save individual metrics
        metrics_file = output_path / f"{pattern_file.stem}_metrics.json"
        analyzer.save_metrics(str(metrics_file))
        
        # Save visualization
        viz_file = output_path / f"{pattern_file.stem}_analysis.png"
        analyzer.visualize_analysis(save_path=str(viz_file), show=False)
    
    # Create summary CSV
    df = pd.DataFrame(all_metrics)
    df.to_csv(output_path / "batch_metrics_summary.csv", index=False)
    
    logging.info(f"Batch analysis complete. Results in {output_path}")


# =============================================================================
# COMMAND LINE INTERFACE
# =============================================================================

def main() -> None:
    """Command-line entry point for single-file or batch topology analysis.

    Usage examples::

        # Single pattern
        python topology_analyzer.py --input pattern.npy --output results/ --dx 5.0

        # Batch mode (all .npy in a directory)
        python topology_analyzer.py --input data/ --output results/ --batch --dx 5.0
    """
    parser = argparse.ArgumentParser(
        description='Topology and Morphometry Analyzer for Turing Patterns'
    )
    
    parser.add_argument('--input', type=str, required=True,
                       help='Input pattern file (.npy) or directory')
    parser.add_argument('--output', type=str, default='results/analysis',
                       help='Output directory')
    parser.add_argument('--dx', type=float, default=1.0,
                       help='Spatial resolution (μm/pixel)')
    parser.add_argument('--batch', action='store_true',
                       help='Batch mode: analyze entire directory')
    parser.add_argument('--show', action='store_true',
                       help='Display visualizations')
    parser.add_argument('--log-level', type=str, default='INFO')
    
    args = parser.parse_args()
    
    # Setup logging
    log_file = Path(args.output) / "logs" / "analysis.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)
    setup_logging(args.log_level, str(log_file))
    
    if args.batch:
        # Batch mode: analyze directory
        analyze_directory(args.input, args.output, dx=args.dx)
    else:
        # Single file mode
        pattern = np.load(args.input)
        analyzer = TopologyAnalyzer(pattern, dx=args.dx)
        metrics = analyzer.compute_all_metrics()
        
        # Save results
        output_path = Path(args.output)
        output_path.mkdir(parents=True, exist_ok=True)
        
        analyzer.save_metrics(str(output_path / "metrics.json"))
        analyzer.visualize_analysis(
            save_path=str(output_path / "analysis.png"),
            show=args.show
        )
    
    logging.info("Analysis complete!")


if __name__ == '__main__':
    main()
