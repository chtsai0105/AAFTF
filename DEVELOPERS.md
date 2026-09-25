# Developer Guide

## Running Tests

Tests live in `tests/` and use pytest. Run from the repo root.

```bash
# Unit tests only — no external bioinformatics tools required
pixi run -e dev pytest tests/ -m unit -v

# All tests — integration tests require the full tool environment
pixi run -e dev pytest tests/ -v
```

## Compile Check

Quick syntax check of all modules without running anything:

```bash
python -m py_compile aaftf/*.py
```

## Lint

```bash
pixi run -e dev pre-commit run --all-files   # ruff, ruff-format, codespell, pydocstyle, ...
```

## Code Style and Architecture

See [AGENTS.md](AGENTS.md) for full development guidelines, subcommand conventions, and code patterns.
