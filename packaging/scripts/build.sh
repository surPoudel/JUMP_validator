#!/usr/bin/env bash
set -euo pipefail

# Build a standalone executable for the *current* OS using PyInstaller.
#
# Usage:
#   ./packaging/scripts/build.sh

python -m pip install -U pip
python -m pip install -r requirements.txt
python -m pip install pyinstaller

pyinstaller -y packaging/pyinstaller/jump-validator.spec
pyinstaller -y packaging/pyinstaller/jump-validator-gui.spec

echo "Build complete. See dist/jump-validator/ and dist/jump-validator-gui/"
