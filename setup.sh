#!/usr/bin/env sh
# Recreate the isolated, pinned runtime; tests are data independent.
set -eu
cd "$(dirname "$0")"
uv sync --extra test --frozen
if [ "${1:-}" = "--tests" ]; then
  .venv/bin/python -m pytest tests -q
fi
