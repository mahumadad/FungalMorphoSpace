#!/usr/bin/env python3
"""Enhanced Visualization Module for Polypore Morphogenesis Validation.

.. module:: fungalmorphospace.analysis.visualization
   :synopsis: Publication-quality figures comparing simulated Turing
              patterns with empirical polypore morphological data.

Overview
--------
This module provides the :class:`EnhancedVisualizer` class, which generates
multi-panel, publication-ready figures for the FungalMorphoSpace validation
pipeline.  The primary outputs are:

* **Per-species field images** -- individual activator and inhibitor
  concentration maps saved as high-resolution PNGs.
* **Comprehensive comparison figure** -- a 4-row, 4-column layout that
  juxtaposes spatial patterns, zoomed-in detail, quantitative scaling-law
  analysis, and a tabular summary for up to four polypore species.
* **Sensitivity heatmaps** -- 2D parameter-sweep visualizations showing
  how a chosen metric (e.g., wavelength) varies across two model parameters.
* **Convergence plots** -- time-series of pattern energy for monitoring
  simulation steady-state convergence.

Design Notes
~~~~~~~~~~~~
All figures follow a consistent aesthetic: bold axis labels (12--14 pt),
300 dpi raster output, and a species-specific color palette inspired by
the natural pigmentation of the target polypore fungi.

Author: Mario Ahumada Duran
Date: December 2025
License: MIT
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path
from scipy.fft import fft2, fftshift
from scipy import ndimage
import pandas as pd

class EnhancedVisualizer:
    """Publication-quality visualizations for polypore morphogenesis validation.

    Generates multi-panel figures that compare simulated Turing-pattern
    outputs against empirical pore-density data for four reference polypore
    species.  The class encapsulates a species-specific color palette and
    several figure-building methods that share consistent formatting
    conventions (font sizes, DPI, colorbar placement).

    Attributes:
        species_colors: Mapping from binomial species name to a hex color
            string used for scatter plots and table highlights.  The
            palette is inspired by the natural pigmentation of each
            fungus:

            * *Fomes fomentarius* -- dark saddle brown (``#8B4513``).
            * *Polyporus brumalis* -- sienna (``#A0522D``).
            * *Trametes versicolor* -- peru / warm tan (``#CD853F``).
            * *Polyporus squamosus* -- burlywood / pale tan (``#DEB887``).
    """

    def __init__(self) -> None:
        """Initialize the visualizer with the default species color palette."""
        self.species_colors: dict[str, str] = {
            'Fomes fomentarius': '#8B4513',
            'Polyporus brumalis': '#A0522D',
            'Trametes versicolor': '#CD853F',
            'Polyporus squamosus': '#DEB887'
        }
    
    @staticmethod
    def _format_params(res: dict) -> str:
        """Build a compact, human-readable parameter string for figure annotations.

        Extracts up to five key simulation parameters from *res* and
        formats each with Python's ``g`` specifier (shortest meaningful
        representation) to keep the annotation short enough to embed in
        PNG titles without clipping.

        Parameters that are absent (``None``) are silently skipped; values
        that cannot be cast to ``float`` are included as-is to avoid
        crashing on unexpected types.

        Args:
            res: Simulation result dictionary.  Recognized keys:

                * ``D_v_D_u`` -- diffusivity ratio.
                * ``b`` -- kinetic parameter.
                * ``rho`` -- reaction-rate scaling.
                * ``grid_size`` -- spatial grid dimension.
                * ``T_target`` -- simulation target time.

        Returns:
            Comma-separated string, e.g.
            ``"D_v/D_u=10, b=0.6, rho=13, grid=256, T=500"``.
            Returns an empty string if no recognized keys are present.
        """
        dvdu = res.get("D_v_D_u", None)
        b = res.get("b", None)
        rho = res.get("rho", None)
        grid = res.get("grid_size", None)
        t_target = res.get("T_target", None)

        parts: list[str] = []
        if dvdu is not None:
            try:
                parts.append(f"D_v/D_u={float(dvdu):g}")
            except Exception:
                parts.append(f"D_v/D_u={dvdu}")
        if b is not None:
            try:
                parts.append(f"b={float(b):g}")
            except Exception:
                parts.append(f"b={b}")
        if rho is not None:
            try:
                parts.append(f"ρ={float(rho):g}")
            except Exception:
                parts.append(f"ρ={rho}")
        if grid is not None:
            parts.append(f"grid={grid}")
        if t_target is not None:
            try:
                parts.append(f"T={float(t_target):g}")
            except Exception:
                parts.append(f"T={t_target}")

        return ", ".join(parts)

    def save_species_field_images(
        self,
        results_dict: dict,
        *,
        out_dir: "str | Path",
        exp_id: str | None = None,
        dpi: int = 300,
    ) -> dict[str, dict[str, str]]:
        """Save per-species activator and inhibitor concentration field PNGs.

        Produces one or two standalone images per species: the activator
        field ``u(x,y)`` (colormap ``hot``) and, when available, the
        inhibitor field ``v(x,y)`` (colormap ``plasma``).  Each image
        includes the species name, field label, and a compact parameter
        annotation as its title.

        This is an *optional* artifact family that can be enabled via the
        validator CLI (``--save-fields``).  Failures on individual images
        are caught and logged so that one broken species does not abort
        the entire batch.

        Args:
            results_dict: Mapping from species name (str) to a result
                dictionary that must contain at least a ``"pattern"`` key
                (2D ``np.ndarray``).  An optional ``"inhibitor_field"``
                key triggers a second image.  Recognized annotation keys
                are those accepted by :meth:`_format_params`.
            out_dir: Output directory.  Created automatically (including
                parents) if it does not exist.
            exp_id: Optional experiment identifier prepended to filenames
                to disambiguate multiple runs sharing the same output
                directory.
            dpi: Resolution for saved PNGs (default 300).

        Returns:
            Dictionary mapping each successfully processed species name
            to a sub-dictionary with keys ``"activator"`` and
            ``"inhibitor"``, each holding the absolute filesystem path
            as a string (or an empty string if the corresponding field
            was absent or failed to render).
        """
        from pathlib import Path

        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        written: dict[str, dict[str, str]] = {}

        for species, res in results_dict.items():
            # Activator is the canonical 'pattern' array.
            u: np.ndarray | None = res.get("pattern", None)
            v: np.ndarray | None = res.get("inhibitor_field", None)
            if u is None:
                # No activator field means we cannot produce any field images.
                continue

            # Build a filesystem-safe base name from species key or name
            species_key = (res.get("species_key") or species.replace(" ", "_")).replace("/", "_")
            base = f"{species_key}" if exp_id is None else f"{exp_id}_{species_key}"

            params_txt: str = self._format_params(res)

            def _save(field: np.ndarray, kind: str, label: str, *, cmap: str) -> Path:
                """Render a single field image and write it to disk."""
                fig, ax = plt.subplots(figsize=(6, 6))
                im = ax.imshow(field, cmap=cmap, interpolation="bilinear")
                ax.set_title(f"{species}\n{label}\n{params_txt}", fontsize=10, fontweight="bold")
                ax.axis("off")
                cbar = plt.colorbar(im, ax=ax, fraction=0.046)
                cbar.set_label(label, fontsize=9)
                path = out_dir / f"{base}_{kind}.png"
                plt.savefig(path, dpi=dpi, bbox_inches="tight")
                plt.close(fig)
                return path

            a_path: Path | None = None
            i_path: Path | None = None

            try:
                a_path = _save(u, "activator_u", "Activator (u)", cmap="hot")
            except Exception:
                # Do not allow a single image failure to break the full validator run.
                a_path = None

            if v is not None:
                try:
                    i_path = _save(v, "inhibitor_v", "Inhibitor (v)", cmap="plasma")
                except Exception:
                    i_path = None

            if a_path is not None or i_path is not None:
                written[species] = {
                    "activator": str(a_path) if a_path is not None else "",
                    "inhibitor": str(i_path) if i_path is not None else "",
                }

        return written

    def create_comprehensive_figure(
        self,
        results_dict: dict[str, dict],
        save_path: str | None = None,
        *,
        annotate_params: bool = True,
    ) -> plt.Figure:
        """Create the main 4-row comparison figure for multi-species validation.

        This is the central visualization artifact of the validation
        pipeline.  It places up to four species side by side and
        combines qualitative pattern inspection with quantitative
        scaling-law analysis in a single, publication-ready figure.

        Layout (4 rows x 4 columns, ``GridSpec(4, 4)``):

        * **Row 0 -- Spatial patterns**: Full activator field per species
          (``hot`` colormap) with wavelength and spot-count annotations.
        * **Row 1 -- Zoomed center**: 128x128 px crop from the field
          center (``nearest`` interpolation) to inspect fine structure.
        * **Row 2 -- Quantitative panels** (two wide sub-panels):

          * *Left (cols 0-1)* -- Wavelength vs. pore density scatter
            with a power-law fit ``lambda = a * rho^b`` (log y-scale).
          * *Right (cols 2-3)* -- Predicted vs. observed density: per-
            species horizontal bars show the empirical range while
            diamond markers show the model prediction.  Green / yellow /
            pink color coding indicates agreement quality.

        * **Row 3 -- Results table**: A ``matplotlib`` table summarizing
          wavelength, spot count, predicted and observed densities,
          diffusivity ratio, and kinetic parameters for quick reference.

        Args:
            results_dict: Mapping from binomial species name to a result
                dictionary.  Required keys per species:

                * ``pattern`` (``np.ndarray``) -- 2D activator field.
                * ``wavelength_px`` (``float``) -- dominant wavelength
                  in pixels.
                * ``density_predicted`` (``float``) -- model pore density.
                * ``density_observed`` (``tuple[float, float]``) --
                  empirical (min, max) pore density range.
                * ``D_v_D_u`` (``float``) -- diffusivity ratio.
                * ``b`` (``float``) -- kinetic parameter.
                * ``rho`` (``float``) -- reaction-rate scaling.
                * ``spots`` (``int``) -- detected spot count.

            save_path: If provided, the figure is saved as a 300 dpi PNG.
            annotate_params: If ``True`` (default), each pattern tile
                includes a compact parameter string below the species
                name.

        Returns:
            The ``matplotlib.figure.Figure`` object (not closed), allowing
            the caller to perform further customization if desired.
        """

        fig = plt.figure(figsize=(20, 14))
        gs = GridSpec(4, 4, figure=fig, hspace=0.35, wspace=0.3)
        
        species_list: list[str] = list(results_dict.keys())

        # ---- ROW 0: Full spatial patterns (one column per species) ----
        for i, species in enumerate(species_list):
            ax = fig.add_subplot(gs[0, i])
            pattern = results_dict[species]['pattern']
            im = ax.imshow(pattern, cmap='hot', interpolation='bilinear')
            
            short_name = ' '.join(species.split()[:2])
            λ = results_dict[species]['wavelength_px']
            spots = results_dict[species].get('spots', 'N/A')

            params_txt = self._format_params(results_dict[species]) if annotate_params else ""
            title_lines = [short_name]
            if params_txt:
                title_lines.append(params_txt)
            title_lines.append(f"λ={λ:.1f} px | Spots={spots}")
            
            ax.set_title("\n".join(title_lines), fontsize=10, fontweight='bold')
            ax.axis('off')
            cbar = plt.colorbar(im, ax=ax, fraction=0.046)
            cbar.set_label('Activator u', fontsize=9)
        
        # ---- ROW 1: Zoomed center crops (128x128 px) ----
        for i, species in enumerate(species_list):
            ax = fig.add_subplot(gs[1, i])
            pattern = results_dict[species]['pattern']
            center = pattern.shape[0] // 2
            zoom_size = 64
            zoom = pattern[center-zoom_size:center+zoom_size, 
                          center-zoom_size:center+zoom_size]
            ax.imshow(zoom, cmap='hot', interpolation='nearest')
            ax.set_title('Zoom (128×128 px)', fontsize=9)
            ax.axis('off')
        
        # ---- ROW 2: Quantitative analyses (two wide sub-panels) ----

        # Panel A (cols 0-1): Wavelength vs. pore density scatter + power-law fit
        ax1 = fig.add_subplot(gs[2, :2])
        
        densities = [results_dict[sp].get('density_obs_mean', 
                    np.mean(results_dict[sp]['density_observed'])) 
                    for sp in species_list]
        wavelengths = [results_dict[sp]['wavelength_px'] for sp in species_list]
        colors = [self.species_colors.get(sp, '#888888') for sp in species_list]
        
        ax1.scatter(densities, wavelengths, s=300, c=colors, 
                   edgecolors='black', linewidths=2.5, alpha=0.8, zorder=3)
        
        # Annotations
        for i, sp in enumerate(species_list):
            genus = sp.split()[0][:4]  # First 4 letters
            ax1.annotate(genus, (densities[i], wavelengths[i]), 
                        xytext=(8, 0), textcoords='offset points', 
                        fontsize=10, fontweight='bold')
        
        # Fit a power-law scaling relation: lambda = a * rho^b
        # Biological expectation: b ~ -0.5 (inverse-square-root) for
        # hexagonal packing of pores on a surface.
        from scipy.optimize import curve_fit

        def power_law(x: np.ndarray, a: float, b: float) -> np.ndarray:
            """Power-law model: y = a * x^b."""
            return a * x**b
        
        try:
            popt, _ = curve_fit(power_law, np.array(densities), 
                               np.array(wavelengths), p0=[100, -1])
            x_fit = np.linspace(min(densities)*0.8, max(densities)*1.2, 100)
            y_fit = power_law(x_fit, *popt)
            ax1.plot(x_fit, y_fit, 'r--', linewidth=2.5, alpha=0.7, 
                    label=f'Fit: λ = {popt[0]:.1f} × ρ^{popt[1]:.2f}', zorder=2)
            ax1.legend(fontsize=11, loc='best')
        except (RuntimeError, ValueError) as e:
            import logging
            logging.debug(f"curve_fit failed in scaling law fit: {e}")
        
        ax1.set_xlabel('Pore Density (pores/mm)', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Wavelength λ (px)', fontsize=12, fontweight='bold')
        ax1.set_title('Scaling Law Validation', fontsize=13, fontweight='bold')
        ax1.grid(True, alpha=0.3, zorder=1)
        ax1.set_yscale('log')
        
        # Panel B (cols 2-3): Predicted vs. observed density comparison
        ax2 = fig.add_subplot(gs[2, 2:])
        
        species_names = [sp.split()[1] for sp in species_list]
        pred_densities = [results_dict[sp].get('density_predicted', 0) 
                         for sp in species_list]
        obs_min = [results_dict[sp]['density_observed'][0] for sp in species_list]
        obs_max = [results_dict[sp]['density_observed'][1] for sp in species_list]
        
        x_pos = np.arange(len(species_list))
        
        # Observed range bars
        ax2.barh(x_pos, np.array(obs_max) - np.array(obs_min), left=obs_min, 
                height=0.4, color='lightblue', alpha=0.6, 
                label='Observed range', edgecolor='blue', linewidth=2)
        
        # Predicted points
        ax2.scatter(pred_densities, x_pos, s=200, c=colors, marker='D', 
                   edgecolors='black', linewidths=2, label='Model prediction', zorder=3)
        
        ax2.set_yticks(x_pos)
        ax2.set_yticklabels(species_names, fontsize=10)
        ax2.set_xlabel('Pore Density (pores/mm)', fontsize=12, fontweight='bold')
        ax2.set_title('Prediction vs Observation', fontsize=13, fontweight='bold')
        ax2.legend(fontsize=10)
        ax2.grid(True, alpha=0.3, axis='x')
        ax2.invert_yaxis()
        
        # ---- ROW 3: Tabular summary ----
        ax3 = fig.add_subplot(gs[3, :])
        ax3.axis('tight')
        ax3.axis('off')
        
        table_data = []
        for sp in species_list:
            res = results_dict[sp]
            row = [
                sp.split()[1],  # Genus
                f"{res['wavelength_px']:.1f}",
                f"{res.get('spots', 'N/A')}",
                f"{res.get('density_predicted', 0):.2f}",
                f"{res['density_observed'][0]}-{res['density_observed'][1]}",
                f"{res['D_v_D_u']:.0f}",
                f"b={res['b']:.2f}, ρ={res['rho']:.2f}"
            ]
            table_data.append(row)
        
        table = ax3.table(
            cellText=table_data,
            colLabels=['Species', 'λ (px)', 'Spots', 'Density\nPred.', 
                      'Density\nObs.', 'D_v/D_u', 'Parameters'],
            cellLoc='center',
            loc='center',
            bbox=[0, 0, 1, 1]
        )
        
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 2.5)
        
        # Color-code the "Density Pred." column:
        #   Green  -- prediction falls within the observed range.
        #   Yellow -- prediction is within 1.5 pores/mm of the range midpoint.
        #   Pink   -- prediction is far outside the observed range.
        for i in range(len(species_list)):
            res = results_dict[species_list[i]]
            pred: float = res.get('density_predicted', 0)
            obs_range: tuple[float, float] = res['density_observed']

            if obs_range[0] <= pred <= obs_range[1]:
                color = '#90EE90'  # Light green -- within observed range
            elif abs(pred - np.mean(obs_range)) < 1.5:
                color = '#FFFFE0'  # Light yellow -- near observed range
            else:
                color = '#FFB6C1'  # Light pink -- far from observed range
            table[(i+1, 3)].set_facecolor(color)
        
        plt.suptitle('INTEGRATED VALIDATION: Polypore Morphogenesis via Reaction-Diffusion', 
                    fontsize=16, fontweight='bold', y=0.995)
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"✓ Figure saved: {save_path}")
        
        return fig
    
    def create_sensitivity_heatmap(
        self,
        sweep_results: dict[tuple[float, float], float],
        param1_name: str,
        param2_name: str,
        metric_name: str = 'wavelength',
        save_path: str | None = None,
    ) -> plt.Figure:
        """Create a 2D heatmap visualizing how a metric varies over a parameter sweep.

        Given a dictionary of ``{(param1, param2): metric_value}`` entries
        produced by a two-parameter sweep, this method reconstructs the 2D
        grid by enumerating the sorted unique values of each parameter,
        fills in the grid (using ``np.nan`` for missing combinations), and
        renders it as an ``imshow`` heatmap with a ``viridis`` colormap.

        Args:
            sweep_results: Mapping from ``(param1_value, param2_value)``
                tuples to the scalar metric value measured at that
                parameter combination.  Missing combinations are rendered
                as transparent (NaN).
            param1_name: Human-readable label for the x-axis (first
                element of each key tuple).
            param2_name: Human-readable label for the y-axis (second
                element of each key tuple).
            metric_name: Label for the colorbar and figure title
                (default ``"wavelength"``).
            save_path: If provided, the figure is saved as a 300 dpi PNG.

        Returns:
            The ``matplotlib.figure.Figure`` object.
        """

        # Extract and sort the unique parameter values along each axis
        param1_vals: list[float] = sorted(set(k[0] for k in sweep_results))
        param2_vals: list[float] = sorted(set(k[1] for k in sweep_results))

        # Reconstruct a 2D grid from the sparse dictionary
        grid = np.zeros((len(param2_vals), len(param1_vals)))

        for i, p2 in enumerate(param2_vals):
            for j, p1 in enumerate(param1_vals):
                grid[i, j] = sweep_results.get((p1, p2), np.nan)

        # Render heatmap
        fig, ax = plt.subplots(figsize=(12, 10))

        im = ax.imshow(grid, aspect='auto', origin='lower', cmap='viridis',
                      extent=[param1_vals[0], param1_vals[-1],
                             param2_vals[0], param2_vals[-1]])

        ax.set_xlabel(param1_name, fontsize=14, fontweight='bold')
        ax.set_ylabel(param2_name, fontsize=14, fontweight='bold')
        ax.set_title(f'Sensitivity Analysis: {metric_name}',
                    fontsize=16, fontweight='bold')

        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label(metric_name, fontsize=12)

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Heatmap saved: {save_path}")

        return fig
    
    def plot_convergence_history(
        self,
        energy_history: list[float] | np.ndarray,
        time_history: list[float] | np.ndarray,
        save_path: str | None = None,
    ) -> plt.Figure:
        """Plot the pattern energy convergence over simulation time.

        Renders a single-axis line plot of ``energy(t)`` on a
        logarithmic y-scale.  A well-converged simulation should show
        the energy curve plateauing (flattening) at late times,
        indicating that the reaction-diffusion system has reached a
        stationary or slowly evolving state.

        The y-axis uses a log scale because pattern energy typically
        decays exponentially during the linear instability phase before
        saturating at nonlinear steady state.

        Args:
            energy_history: 1D sequence of scalar energy values sampled
                at each output time step.  Typically computed as
                ``sum((u - u_mean)^2)`` or a free-energy functional.
            time_history: 1D sequence of simulation times corresponding
                to *energy_history*.  Must have the same length.  Units
                are arbitrary (``a.u.``).
            save_path: If provided, the figure is saved at 200 dpi.

        Returns:
            The ``matplotlib.figure.Figure`` object.
        """

        fig, ax = plt.subplots(figsize=(10, 6))

        ax.plot(time_history, energy_history, 'b-', linewidth=2)
        ax.set_xlabel('Time (a.u.)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Pattern Energy', fontsize=12, fontweight='bold')
        ax.set_title('Convergence History', fontsize=14, fontweight='bold')
        ax.set_yscale('log')
        ax.grid(True, alpha=0.3)

        if save_path:
            plt.savefig(save_path, dpi=200, bbox_inches='tight')

        return fig
