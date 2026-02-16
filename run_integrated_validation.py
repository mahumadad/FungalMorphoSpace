#!/usr/bin/env python3
"""Convenience wrapper for the integrated validation CLI.

This exists because many users naturally run:

    python3 run_integrated_validation.py --species all --n_runs 1

The canonical entrypoint remains:

    python3 scripts/run_integrated_validation.py ...
"""

from __future__ import annotations

import runpy
from pathlib import Path

if __name__ == "__main__":
    script = Path(__file__).resolve().parent / "scripts" / "run_integrated_validation.py"
    runpy.run_path(str(script), run_name="__main__")
