# Contributing to FungalMorphoSpace

Thank you for your interest in contributing to FungalMorphoSpace! This document provides guidelines for contributing to the project.

## Reporting Bugs

1. Check the [issue tracker](https://github.com/mahumada408/FungalMorphoSpace/issues) to see if the bug has already been reported.
2. If not, open a new issue with:
   - A clear, descriptive title
   - Steps to reproduce the bug
   - Expected behavior vs. actual behavior
   - Python version, OS, and dependency versions (`pip freeze`)
   - Error tracebacks if applicable

## Suggesting Features

Open an issue with the `enhancement` label describing:
- The problem your feature would solve
- How you envision the solution
- Whether you're willing to implement it

## Development Setup

```bash
# Clone the repository
git clone https://github.com/mahumada408/FungalMorphoSpace.git
cd FungalMorphoSpace

# Create a virtual environment
python -m venv .venv
source .venv/bin/activate  # Linux/macOS
# .venv\Scripts\activate   # Windows

# Install in development mode
pip install -e .
pip install -r requirements-dev.txt

# Verify installation
python scripts/test_imports.py
```

## Running Tests

```bash
# Unit tests
pytest tests/

# Integration smoke test
python scripts/smoke_test.py --grid 256
```

## Code Style

- **Docstrings:** Google style (see existing code for examples)
- **Type hints:** Use type annotations on all public function signatures
- **Language:** Docstrings and code comments in English; scientific documentation may be bilingual (English/Spanish)

## Pull Request Process

1. Fork the repository and create your branch from `main`.
2. Add or update tests for any new functionality.
3. Ensure all tests pass (`pytest tests/`).
4. Update documentation if you changed public APIs.
5. Write a clear PR description explaining the changes and motivation.

## Adding New Species

To add a new validated species:

1. Add calibrated parameters to `data/species_data.json` following the existing schema.
2. Document the calibration process and biological references.
3. Run the validation pipeline and verify metrics.
4. Update the species table in `README.md`.

## License

By contributing, you agree that your contributions will be licensed under the project's dual license (CC BY-NC 4.0 for academic use / Commercial license for commercial use).
