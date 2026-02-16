"""Package entry point.

Allows running:
  python -m fungalmorphospace ...

Defaults to the simulator CLI.
"""

from .core.turing_simulator import main

if __name__ == '__main__':
    raise SystemExit(main())
