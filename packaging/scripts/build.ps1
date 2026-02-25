# Build a standalone executable for the *current* OS using PyInstaller.
#
# Usage (PowerShell):
#   powershell -ExecutionPolicy Bypass -File .\packaging\scripts\build.ps1

python -m pip install -U pip
python -m pip install -r requirements.txt
python -m pip install pyinstaller

pyinstaller -y packaging/pyinstaller/jump-validator.spec
pyinstaller -y packaging/pyinstaller/jump-validator-gui.spec

Write-Host "Build complete. See dist\jump-validator\ and dist\jump-validator-gui\" 
