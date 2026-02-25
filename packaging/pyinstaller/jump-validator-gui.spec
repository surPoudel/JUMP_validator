# -*- mode: python ; coding: utf-8 -*-

# PyInstaller spec to build a standalone `jump-validator-gui` executable.
#
# Build (current OS only):
#   pyinstaller -y packaging/pyinstaller/jump-validator-gui.spec

from PyInstaller.utils.hooks import collect_data_files, collect_submodules

hiddenimports = []

# These libraries often need explicit hidden-import collection.
for pkg in [
    "pyteomics",
    "seaborn",
    "matplotlib",
    "PIL",
    "dataframe_image",
    "docx",
    "openpyxl",
    "xlsxwriter",
]:
    hiddenimports += collect_submodules(pkg)

# Matplotlib needs data files (fonts, etc.) for consistent rendering.
datas = collect_data_files("matplotlib")

block_cipher = None


a = Analysis(
    ["jump_validator/gui.py"],
    pathex=["."],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name="jump-validator-gui",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,  # Windowed app (no console)
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name="jump-validator-gui",
)
