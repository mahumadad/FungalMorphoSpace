#!/usr/bin/env python3
"""
FUNGAL MORPHOSPACE - Scaling Law Validation Plot
================================================
Generates log-log plot to validate allometric consistency between:

- X: biological scale (mm), from `data/species_data.json` (pore_spacing_um)
- Y: simulated wavelength λ (px), from `results/tables/validation_summary_machine.csv`

v0.6.2: Reads from machine-readable CSV with correct column names.
        Supports --wavelength-source {best,fft,autocorr}.

Author: Mario Ahumada Durán
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress


# Short canonical keys -> display labels
SPECIES_ORDER = [
    ("brumalis", "P. brumalis"),
    ("fomes", "F. fomentarius"),
    ("squamosus", "P. squamosus"),
]


def load_species_json(json_path: Path) -> dict:
    if not json_path.exists():
        raise FileNotFoundError(f"Species JSON not found: {json_path}")
    with open(json_path, "r", encoding="utf-8") as f:
        return json.load(f)


def load_simulated_wavelengths(summary_csv: Path, wavelength_source: str) -> dict[str, float]:
    """Return dict mapping canonical species key -> wavelength_px (float).
    
    Args:
        summary_csv: Path to validation_summary_machine.csv (preferred) or validation_summary.csv
        wavelength_source: One of 'best', 'fft', 'autocorr'
    """
    if not summary_csv.exists():
        return {}

    df = pd.read_csv(summary_csv)
    
    # Map column names based on source preference
    # Machine-readable CSV uses these columns:
    column_map = {
        "best": "wavelength_best_px",
        "fft": "wavelength_fft_px",
        "autocorr": "wavelength_autocorr_px",
    }
    
    # Human-readable CSV fallback columns
    human_column_map = {
        "best": "λ(best) (px)",
        "fft": "λ(FFT) (px)",
        "autocorr": "λ(autocorr) (px)",
    }
    
    # Determine which column to use
    target_col = column_map.get(wavelength_source, "wavelength_best_px")
    human_col = human_column_map.get(wavelength_source, "λ(best) (px)")
    
    # Try machine column first, then human column
    if target_col in df.columns:
        col_to_use = target_col
    elif human_col in df.columns:
        col_to_use = human_col
    else:
        # Try to find any wavelength column
        for c in ["wavelength_best_px", "wavelength_px", "λ(best) (px)", "Wavelength (px)"]:
            if c in df.columns:
                col_to_use = c
                print(f"⚠️ Using fallback column: {c}")
                break
        else:
            print(f"⚠️ No wavelength column found in {summary_csv}")
            return {}
    
    # Map scientific names to canonical keys
    name_to_key = {
        "Fomes fomentarius": "fomes",
        "Polyporus brumalis": "brumalis",
        "Polyporus squamosus": "squamosus",
        # Also support legacy long keys
        "fomes_fomentarius": "fomes",
        "polyporus_brumalis": "brumalis",
        "polyporus_squamosus": "squamosus",
    }

    out: dict[str, float] = {}
    
    # Try to find species column
    species_col = None
    for c in ["species", "Species", "species_key"]:
        if c in df.columns:
            species_col = c
            break
    
    if species_col is None:
        print(f"⚠️ No species column found in {summary_csv}")
        return {}
    
    for _, row in df.iterrows():
        species_name = str(row.get(species_col, "")).strip()
        key = name_to_key.get(species_name, species_name.lower())
        
        # Normalize to short key
        if key in name_to_key.values():
            pass
        elif key in name_to_key:
            key = name_to_key[key]
        else:
            continue

        try:
            raw = row.get(col_to_use)
            wl = float(str(raw).strip())
            out[key] = wl
        except (ValueError, TypeError):
            continue

    return out


def main():
    parser = argparse.ArgumentParser(description="Generate scaling law validation plot")
    parser.add_argument(
        "--output",
        type=str,
        default="results/figures/scaling_validation.png",
        help="Output path for figure",
    )
    parser.add_argument(
        "--summary-csv",
        type=str,
        default=None,
        help="CSV with pipeline results. Default: tries tables/validation_summary_machine.csv, then validation_summary.csv",
    )
    parser.add_argument(
        "--species-json",
        type=str,
        default="data/species_data.json",
        help="Canonical species data file",
    )
    parser.add_argument(
        "--wavelength-source",
        type=str,
        choices=["best", "fft", "autocorr"],
        default="best",
        help="Which wavelength estimator to use (default: best)",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Don't display the figure (only save)",
    )

    args = parser.parse_args()
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    json_path = Path(args.species_json)
    
    # Find summary CSV
    if args.summary_csv:
        summary_csv = Path(args.summary_csv)
    else:
        # Try machine-readable first, then human-readable
        candidates = [
            Path("results/tables/validation_summary_machine.csv"),
            Path("results/tables/validation_summary.csv"),
            Path("results/validation_summary.csv"),
        ]
        summary_csv = None
        for c in candidates:
            if c.exists():
                summary_csv = c
                print(f"Using summary CSV: {summary_csv}")
                break
        if summary_csv is None:
            summary_csv = candidates[0]  # Will fail gracefully

    species_data = load_species_json(json_path)["species"]
    sim_from_csv = load_simulated_wavelengths(summary_csv, args.wavelength_source)

    bio_sizes_mm = []
    sim_sizes_px = []
    labels = []
    source_used = []
    missing = []

    for key, label in SPECIES_ORDER:
        sp = species_data.get(key, {})
        pore_spacing_um = sp.get("pore_spacing_um")
        
        if pore_spacing_um is None:
            print(f"⚠️ Missing pore_spacing_um for {key}, skipping")
            continue

        bio_mm = float(pore_spacing_um) / 1000.0  # um -> mm

        # Prefer measured wavelengths from CSV (pipeline output)
        if key in sim_from_csv:
            sim_px = float(sim_from_csv[key])
            source_used.append(f"csv({args.wavelength_source})")
        else:
            # Fallback: expected wavelength stored in canonical JSON
            exp_wl = sp.get("expected_wavelength_px")
            if exp_wl is None:
                print(f"⚠️ No wavelength data for {key}, skipping")
                continue
            sim_px = float(exp_wl)
            source_used.append("expected_json")
            missing.append(key)

        bio_sizes_mm.append(bio_mm)
        sim_sizes_px.append(sim_px)
        labels.append(label)

    if len(bio_sizes_mm) < 2:
        raise RuntimeError("Not enough data points for regression")

    bio_sizes_mm = np.array(bio_sizes_mm, dtype=float)
    sim_sizes_px = np.array(sim_sizes_px, dtype=float)

    # Log-log regression
    log_bio = np.log10(bio_sizes_mm)
    log_sim = np.log10(sim_sizes_px)
    slope, intercept, r_value, p_value, std_err = linregress(log_bio, log_sim)

    # Plot
    plt.figure(figsize=(8, 6))
    colors = plt.cm.Set1(np.linspace(0, 1, len(labels)))
    
    for i in range(len(labels)):
        plt.scatter(bio_sizes_mm[i], sim_sizes_px[i], s=120, c=[colors[i]], 
                   edgecolors="k", zorder=5, label=labels[i])

    x_fit = np.linspace(0.1, 3.0, 200)
    y_fit = 10 ** (intercept + slope * np.log10(x_fit))
    plt.plot(x_fit, y_fit, "k--", alpha=0.6, label=f"Fit ($R^2={r_value**2:.3f}$, slope={slope:.2f})")

    plt.xscale("log")
    plt.yscale("log")
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.xlabel("Biological scale (pore spacing, mm)")
    plt.ylabel(f"Simulated wavelength $\\lambda$ (px) [{args.wavelength_source}]")
    plt.title("Allometric Scaling Validation\n(Model vs Biology)")
    plt.legend()
    plt.tight_layout()

    plt.savefig(output_path, dpi=300)
    print(f"✅ Figure saved: {output_path}")
    print(f"   R²: {r_value**2:.4f}")
    print(f"   Slope: {slope:.2f}")
    print(f"   Wavelength source: {args.wavelength_source}")
    print(f"   Data sources: {', '.join(set(source_used))}")
    
    if missing:
        print("⚠️  Warning: some λ(px) came from expected_wavelength_px (JSON), not from CSV:")
        for k in missing:
            print(f"   - {k}")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
