#!/usr/bin/env python3
"""
FUNGAL MORPHOSPACE - Robustness Analysis
==========================================
Statistical validation through multiple independent runs

Author: Mario Ahumada Durán
Project: FungalMorphoSpace
Date: January 2026
"""

import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from datetime import datetime

# Package-first imports
base_path = Path(__file__).resolve().parents[1]

# Ensure local package is used (avoid shadowing by any older site-packages installation)
src_path = base_path / "src"
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))


try:
    from fungalmorphospace.core.turing_simulator import TuringSimulator
    from fungalmorphospace.core.kinetics import create_kinetics
    from fungalmorphospace.analysis.topology_analyzer import TopologyAnalyzer
except ImportError:
    sys.path.insert(0, str(base_path / "src"))
    from fungalmorphospace.core.turing_simulator import TuringSimulator
    from fungalmorphospace.core.kinetics import create_kinetics
    from fungalmorphospace.analysis.topology_analyzer import TopologyAnalyzer

# =============================================================================
# SPECIES PARAMETERS - LOADED FROM CANONICAL SOURCE (data/species_data.json)
# =============================================================================

def load_species_params():
    """Load species parameters from canonical JSON source."""
    json_path = base_path / "data" / "species_data.json"
    
    if not json_path.exists():
        raise FileNotFoundError(f"Species data not found: {json_path}")
    
    with open(json_path, 'r') as f:
        data = json.load(f)
    
    # Convert JSON format to robustness format
    params = {}
    key_mapping = {
        'fomes_fomentarius': 'fomes',
        'polyporus_brumalis': 'brumalis',
        'polyporus_squamosus': 'squamosus'
    }
    
    for species_key, species_data in data['species'].items():
        if species_data.get('role') == 'EXCLUDED':
            continue
            
        short_key = key_mapping.get(species_key, species_key)
        sp_params = species_data.get('parameters', {})
        
        params[short_key] = {
            'name': species_data.get('scientific_name', species_key),
            'D_v_D_u': sp_params.get('D_v_D_u', 150.0),
            'b': sp_params.get('b', 1.0),
            'rho': sp_params.get('rho', 0.2),
            'grid': sp_params.get('grid_size', 1024),  # Use higher resolution for robustness
            'T_target': sp_params.get('T_target', 5.0),
            'dt_initial': sp_params.get('dt_initial', 0.0005)
        }
    
    return params

SPECIES_PARAMS = load_species_params()

# =============================================================================
# CALIBRATION (px → μm) - LOADED FROM CANONICAL SOURCE (data/species_data.json)
# =============================================================================
def load_calibration_um_per_px():
    json_path = base_path / "data" / "species_data.json"
    with open(json_path, 'r') as f:
        data = json.load(f)
    cal = data.get('calibration', {}) or {}
    return float(cal.get('scale_factor_um_per_px', 400/46))

CALIBRATION_UM_PER_PX = load_calibration_um_per_px()


# =============================================================================
# ROBUSTNESS ANALYSIS FUNCTION
# =============================================================================

def run_robustness_analysis(species_key, n_runs=20, save_patterns=False):
    """
    Run multiple simulations with different random seeds to assess stability
    
    Parameters:
    -----------
    species_key : str
        Species identifier ('fomes', 'brumalis', 'squamosus')
    n_runs : int
        Number of independent runs (default: 20)
    save_patterns : bool
        Whether to save all pattern images (default: False)
    
    Returns:
    --------
    DataFrame with results from all runs
    """
    
    if species_key not in SPECIES_PARAMS:
        raise ValueError(f"Unknown species: {species_key}")
    
    params = SPECIES_PARAMS[species_key]
    species_name = params['name']
    
    print("="*70)
    print(f"ROBUSTNESS ANALYSIS: {species_name}")
    print("="*70)
    print(f"Number of runs: {n_runs}")
    print(f"Parameters: D_v/D_u={params['D_v_D_u']}, b={params['b']}, rho={params['rho']}")
    print()
    
    results = []
    
    for run in range(n_runs):
        print(f"Run {run+1}/{n_runs}...", end=' ', flush=True)
        
        # Unique seed per run for true statistical replication
        seed = 42 + run
        
        # Physics
        D_u = 1.0
        D_v = D_u * params['D_v_D_u']
        dt_requested = params.get('dt_initial', 0.0005)  # From JSON or default
        
        # Kinetics
        kinetics = create_kinetics('gierer_meinhardt',
                                   rho=params['rho'],
                                   a=0.1,
                                   b=params['b'])
        
        # Simulate with unique seed passed to constructor
        sim = TuringSimulator(
            kinetics_model=kinetics,
            D_u=D_u,
            D_v=D_v,
            grid_size=params['grid'],
            dx=1.0,
            dt=dt_requested,
            random_seed=seed  # CRITICAL: unique seed per run
        )
        
        # Calculate steps from T_target using dt_final (after CFL adjustment)
        import math
        dt_final = sim.dt
        T_target = params.get('T_target', 5.0)
        steps_computed = math.ceil(T_target / dt_final)
        
        sim.initialize(perturbation_amplitude=0.1)
        sim.run(num_steps=steps_computed, check_convergence=False)
        
        # Analyze
        analyzer = TopologyAnalyzer(sim.u, dx=1.0)
        metrics = analyzer.compute_all_metrics()
        
        # Save pattern if requested
        if save_patterns:
            output_dir = base_path / 'results' / 'robustness' / species_key
            output_dir.mkdir(parents=True, exist_ok=True)
            pattern_path = output_dir / f'run_{run+1:02d}.png'
            analyzer.visualize_analysis(save_path=str(pattern_path), show=False)
        
        # Store results
        wavelength = metrics.get('wavelength_autocorr', 0)
        spots = metrics['n_components']
        euler = metrics['euler_characteristic']
        
        # Calculate density
        scale_um_per_px = CALIBRATION_UM_PER_PX  # Calibration from Fomes
        spacing_um = wavelength * scale_um_per_px
        density = 1000 / spacing_um if spacing_um > 0 else 0
        
        results.append({
            'run': run + 1,
            'seed': seed,
            'wavelength_px': wavelength,
            'spots': spots,
            'euler_chi': euler,
            'density_pores_mm': density
        })
        
        print(f"λ={wavelength:.1f} px, spots={spots}")
    
    # Convert to DataFrame
    df = pd.DataFrame(results)
    
    # Calculate statistics
    print("\n" + "="*70)
    print("STATISTICAL SUMMARY")
    print("="*70)
    
    stats = {
        'wavelength_mean': df['wavelength_px'].mean(),
        'wavelength_std': df['wavelength_px'].std(),
        'wavelength_cv': df['wavelength_px'].std() / df['wavelength_px'].mean() * 100,
        'spots_mean': df['spots'].mean(),
        'spots_std': df['spots'].std(),
        'density_mean': df['density_pores_mm'].mean(),
        'density_std': df['density_pores_mm'].std()
    }
    
    print(f"\nWavelength (λ):")
    print(f"  Mean ± SD: {stats['wavelength_mean']:.2f} ± {stats['wavelength_std']:.2f} px")
    print(f"  CV: {stats['wavelength_cv']:.2f}%")
    
    print(f"\nSpots:")
    print(f"  Mean ± SD: {stats['spots_mean']:.1f} ± {stats['spots_std']:.1f}")
    
    print(f"\nDensity:")
    print(f"  Mean ± SD: {stats['density_mean']:.2f} ± {stats['density_std']:.2f} pores/mm")
    
    # Save results
    output_dir = base_path / 'results' / 'robustness'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    csv_path = output_dir / f'{species_key}_robustness_n{n_runs}.csv'
    df.to_csv(csv_path, index=False)
    print(f"\n✓ Data saved: {csv_path}")
    
    # Create visualization
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    # Wavelength distribution
    axes[0].hist(df['wavelength_px'], bins=15, color='steelblue', edgecolor='black', alpha=0.7)
    axes[0].axvline(stats['wavelength_mean'], color='red', linestyle='--', linewidth=2, 
                    label=f"Mean: {stats['wavelength_mean']:.1f} px")
    axes[0].set_xlabel('Wavelength (px)', fontsize=11)
    axes[0].set_ylabel('Frequency', fontsize=11)
    axes[0].set_title(f'λ Distribution (CV={stats["wavelength_cv"]:.1f}%)', fontsize=12)
    axes[0].legend()
    axes[0].grid(alpha=0.3)
    
    # Spots distribution
    axes[1].hist(df['spots'], bins=15, color='forestgreen', edgecolor='black', alpha=0.7)
    axes[1].axvline(stats['spots_mean'], color='red', linestyle='--', linewidth=2,
                    label=f"Mean: {stats['spots_mean']:.0f}")
    axes[1].set_xlabel('Number of Spots', fontsize=11)
    axes[1].set_ylabel('Frequency', fontsize=11)
    axes[1].set_title('Spots Distribution', fontsize=12)
    axes[1].legend()
    axes[1].grid(alpha=0.3)
    
    # Wavelength vs Spots scatter
    axes[2].scatter(df['wavelength_px'], df['spots'], alpha=0.6, s=50, color='purple')
    axes[2].set_xlabel('Wavelength (px)', fontsize=11)
    axes[2].set_ylabel('Spots', fontsize=11)
    axes[2].set_title('λ vs Spots Correlation', fontsize=12)
    axes[2].grid(alpha=0.3)
    
    plt.suptitle(f'{species_name} - Robustness Analysis (n={n_runs})', 
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    fig_path = output_dir / f'{species_key}_robustness_n{n_runs}.png'
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    print(f"✓ Figure saved: {fig_path}")
    
    plt.close()
    
    return df, stats

# =============================================================================
# MAIN
# =============================================================================

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Robustness Analysis for FungalMorphoSpace')
    parser.add_argument('--species', type=str, required=True,
                       choices=['fomes', 'brumalis', 'squamosus'],
                       help='Species to analyze (fomes, brumalis, squamosus)')
    parser.add_argument('--n_runs', type=int, default=20,
                       help='Number of independent runs (default: 20)')
    parser.add_argument('--save_patterns', action='store_true',
                       help='Save all pattern images')
    
    args = parser.parse_args()
    
    # Run analysis
    df, stats = run_robustness_analysis(
        species_key=args.species,
        n_runs=args.n_runs,
        save_patterns=args.save_patterns
    )
    
    print("\n" + "="*70)
    print("✅ ROBUSTNESS ANALYSIS COMPLETED")
    print("="*70)
