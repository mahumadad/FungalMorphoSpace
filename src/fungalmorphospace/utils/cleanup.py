"""Safe cleanup utilities for FungalMorphoSpace temporary artifacts.

Provides functions to remove known development/test temporary directories
from the repository root without risking deletion of user data or
non-FMS directories.

Safety model:
    Deletion is guarded by :func:`_looks_like_fms_output_dir`, a heuristic
    that checks for characteristic FMS subdirectory names (``tables``,
    ``figures``, ``patterns``, etc.) or canonical output files
    (``COMPREHENSIVE_COMPARISON.png``, ``validation_summary.csv``, etc.).
    Only directories passing this check are removed.

Typical usage::

    from fungalmorphospace.utils.cleanup import cleanup_repo_temp_dirs
    from pathlib import Path

    removed_count, removed_names = cleanup_repo_temp_dirs(Path("/path/to/repo"))
"""

from __future__ import annotations

from pathlib import Path
import shutil
from typing import Iterable, Tuple


#: Names of temporary directories created by CI, contract tests, and
#: development scripts.  These are located directly under the repository
#: root and are safe to remove when they match the FMS output heuristic.
FMS_TEMP_DIR_NAMES: tuple[str, ...] = ("tmp_contract_test", "tmp_check", "tmp_check2")


def _looks_like_fms_output_dir(p: Path) -> bool:
    """Heuristic safety check: only return ``True`` for FMS-like directories.

    A directory is considered an FMS artifact if it either:

    1. Contains at least one of the characteristic subdirectories used by
       the simulation runner (``tables``, ``figures``, ``patterns``,
       ``logs``, ``results``, ``metrics``, ``analysis``), **or**
    2. Contains at least one of the known canonical output files
       (``COMPREHENSIVE_COMPARISON.png``, ``validation_summary.csv``,
       ``validation_summary_machine.csv``, ``validation_summary.json``).

    This prevents accidental deletion of unrelated directories that happen
    to share a name with an FMS temp directory.

    Args:
        p: Path to evaluate.

    Returns:
        ``True`` if the path is an existing directory that matches the
        FMS output heuristic; ``False`` otherwise.
    """
    if not p.is_dir():
        return False

    # Check for characteristic FMS subdirectory structure.
    markers: tuple[str, ...] = ("tables", "figures", "patterns", "logs", "results", "metrics", "analysis")
    for m in markers:
        if (p / m).exists():
            return True

    # Also treat as safe if it contains known canonical output files.
    known_files: tuple[str, ...] = (
        "COMPREHENSIVE_COMPARISON.png",
        "validation_summary.csv",
        "validation_summary_machine.csv",
        "validation_summary.json",
    )
    for f in known_files:
        if (p / f).exists():
            return True

    return False


def cleanup_repo_temp_dirs(
    repo_root: Path,
    names: Iterable[str] = FMS_TEMP_DIR_NAMES,
) -> Tuple[int, list[str]]:
    """Delete known development temporary directories from the repository root.

    Iterates over the given directory *names*, and for each one that
    exists under *repo_root* and passes the
    :func:`_looks_like_fms_output_dir` safety heuristic, recursively
    removes the entire directory tree.

    This function is **best-effort**: individual deletion failures are
    silently suppressed so that a single permission error does not prevent
    cleanup of other directories.

    Args:
        repo_root: Path to the repository root directory.
        names: Iterable of directory names to consider for removal.
            Defaults to :data:`FMS_TEMP_DIR_NAMES`.

    Returns:
        A 2-tuple ``(count_removed, removed_names)`` where
        *count_removed* is the number of directories successfully deleted
        and *removed_names* is a list of their base names.
    """
    removed: list[str] = []
    for name in names:
        p: Path = repo_root / name
        if p.exists() and _looks_like_fms_output_dir(p):
            try:
                shutil.rmtree(p)
                removed.append(name)
            except Exception:
                # Best-effort: do not raise in cleanup.
                pass
    return len(removed), removed
