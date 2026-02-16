#!/usr/bin/env python3
"""
Quick test for parallel validation scripts
"""

import sys
from pathlib import Path

print("="*70)
print("TESTING PARALLEL VALIDATION IMPORTS")
print("="*70)

# Ensure local package is importable
base_path = Path(__file__).resolve().parent.parent
src_path = base_path / "src"
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

# Also allow importing CLI modules by name
scripts_path = base_path / "scripts"
if str(scripts_path) not in sys.path:
    sys.path.insert(0, str(scripts_path))

try:
    print("\n1. Testing import of IntegratedSimulationRunner (package)...")
    from fungalmorphospace.runners import IntegratedSimulationRunner, SPECIES_DATABASE
    print("   ✓ IntegratedSimulationRunner imported successfully")
    print(f"   ✓ SPECIES_DATABASE has {len(SPECIES_DATABASE)} species:")
    for key in SPECIES_DATABASE.keys():
        print(f"      - {key}: {SPECIES_DATABASE[key]['full_name']}")
    
    print("\n2. Testing ParallelSafeValidator (CLI module)...")
    from run_parallel_validation import ParallelSafeValidator
    print("   ✓ ParallelSafeValidator imported successfully")
    
    print("\n3. Creating validator instance...")
    validator = ParallelSafeValidator(output_dir='test_results')
    print("   ✓ Validator created successfully")
    
    print("\n4. Checking output directories...")
    print(f"   - Output dir: {validator.output_dir}")
    print(f"   - Patterns dir exists: {(validator.output_dir / 'patterns').exists()}")
    
    print("\n" + "="*70)
    print("✅ ALL IMPORTS SUCCESSFUL - Scripts are ready to use!")
    print("="*70)
    print("\nYou can now run:")
    print("  python scripts/run_parallel_validation.py --species fomes --n_runs 1")
    print("\n" + "="*70)
    
except Exception as e:
    print(f"\n❌ ERROR: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
