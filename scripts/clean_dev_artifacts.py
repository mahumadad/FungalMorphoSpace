#!/usr/bin/env python3
"""Clean up dev/test artifacts created by local runs.

Removes known temporary folders that can be created during local checks:
- tmp_contract_test
- tmp_check
- tmp_check2

This script is intentionally conservative: it only deletes these directories if they
look like FungalMorphoSpace output folders.
"""

from pathlib import Path
import sys

# Ensure src/ is importable
REPO_ROOT = Path(__file__).resolve().parents[1]
SRC = REPO_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from fungalmorphospace.utils.cleanup import cleanup_repo_temp_dirs, FMS_TEMP_DIR_NAMES


def main() -> None:
    removed_n, removed = cleanup_repo_temp_dirs(REPO_ROOT, FMS_TEMP_DIR_NAMES)
    if removed_n == 0:
        print("No dev temp dirs removed.")
        return
    print(f"Removed {removed_n} dev temp dir(s): {', '.join(removed)}")


if __name__ == "__main__":
    main()
