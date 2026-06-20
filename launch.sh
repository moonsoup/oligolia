#!/bin/bash
# Oligolia — launch standalone GUI app
# Usage: ./launch.sh
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

VENV="backend/.venv"
PYTHON="$VENV/bin/python"

if [ ! -f "$PYTHON" ]; then
    echo "Creating virtual environment…"
    python3 -m venv "$VENV"
    "$VENV/bin/pip" install --quiet -r backend/requirements.txt PyQt6
fi

exec "$PYTHON" oligolia.py "$@"
