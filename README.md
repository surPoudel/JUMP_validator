# JUMP Validator (spectral library / MS2 validation)

This repo packages the legacy **JUMP spectral validator** into a small, cross‑platform app:

- **CLI**: batch/automation friendly
- **GUI**: pick a `.params` file, run validation, and **browse the per‑peptide plots** inside the app

✅ The **core algorithm + output files** are unchanged from the legacy scripts (workflow/packaging improvements only).

---

## What it produces

Given a `jump_validator.params`, the validator generates (same as legacy):

- `ManualValidation/Report_of_manual_validation.xlsx`
- `ManualValidation/all_ms2_input.pkl`
- Per‑spectrum plots under `ManualValidation/` (e.g. `__intensityPlot.png/.pdf`, `__ToleranceMap.png/.pdf`, etc.)

The GUI “Results Browser” simply previews the PNG plots already written to disk.

---

## Requirements

- Python **3.9+**
- Packages in `requirements.txt`
- **Tkinter** (for the GUI)
  - macOS: included with python.org / conda Python
  - Windows: included with python.org / conda Python
  - Linux: you may need `python3-tk` (see Linux section)

---

## Install

You can install either with **pip/venv** (recommended) or **conda/mamba**.

### macOS / Linux (pip + venv)

```bash
cd JUMP_validator_app

python3 -m venv .venv
source .venv/bin/activate

python -m pip install -U pip
python -m pip install -r requirements.txt
```

(Optional) install as an editable package to get the `jump-validator` command:

```bash
python -m pip install -e .
```

### Windows (pip + venv, PowerShell)

```powershell
cd JUMP_validator_app

py -m venv .venv
.\.venv\Scripts\Activate.ps1

python -m pip install -U pip
python -m pip install -r requirements.txt
```

(Optional) editable install:

```powershell
python -m pip install -e .
```

### Linux note: Tkinter

If the GUI won’t launch on Linux and you see a Tk error, install Tkinter:

- Debian/Ubuntu:
  ```bash
  sudo apt-get update
  sudo apt-get install -y python3-tk
  ```

(You can always use the CLI on headless servers.)

### Conda / mamba (all OS)

```bash
cd JUMP_validator_app
conda env create -f environment.yml
conda activate jump_validator
```

---

## Run

### CLI

If you did the editable install:

```bash
jump-validator /path/to/jump_validator.params
```

Without installing:

```bash
python -m jump_validator /path/to/jump_validator.params
```

### GUI (with plot viewer)

If you did the editable install:

```bash
jump-validator-gui
```

Without installing:

```bash
python -m jump_validator.gui
```

In the GUI:
1. Choose the `.params` file
2. Click **Run**
3. When finished, click **View results** to browse the plots per spectrum

---

## Params file notes

The validator reads configuration from the `.params` file (INI format).

Common keys under `[caseEvaluation]`:

- `ms2_fileType = mzXML` **or** `ms2`
- If `ms2_fileType = mzXML` → provide `mzXML_path = ...`
- If `ms2_fileType = ms2` → provide `ms2_path = ...`

`peptideFile` can be Excel or tab‑delimited text (the tool will try both).

### If you see: “No matching spectra found…”

This means spectra keys in `peptideFile` (run/scan/charge) do not match spectra keys in the ID file.

Common causes:
- scan numbers were **renumbered during mzXML concatenation**
- scan/charge were parsed as floats (e.g. `12345.0`)
- run name differs between `peptideFile` and the ID file

This build includes an **automatic fallback** (only when there are zero matches) to infer
per‑run scan offsets and align scans where possible.

---

## Building standalone executables (optional)

This repo includes PyInstaller specs and a GitHub Actions workflow to build binaries for:

- macOS
- Windows
- Linux

### Local build (current OS only)

```bash
python -m pip install -r requirements.txt
python -m pip install pyinstaller

pyinstaller -y packaging/pyinstaller/jump-validator.spec
pyinstaller -y packaging/pyinstaller/jump-validator-gui.spec
```

Outputs appear in `dist/`.

### CI build (recommended)

See: `.github/workflows/build-binaries.yml`

When you push to GitHub, Actions can build and upload artifacts for each OS.

---

## Pushing this project to your GitHub

I can’t directly push to your GitHub account from here, but you can do it in a minute:

### Option A: using git (works everywhere)

1) Create an empty repo on GitHub (e.g. `JUMP_validator_app`)

2) From the folder that contains `JUMP_validator_app/`:

```bash
cd JUMP_validator_app
git init
git add .
git commit -m "Add JUMP Validator (CLI+GUI) cross-platform package"

# replace with your repo URL
git remote add origin https://github.com/<YOUR_USER>/<YOUR_REPO>.git
git branch -M main
git push -u origin main
```

### Option B: using GitHub CLI (`gh`)

```bash
cd JUMP_validator_app
git init
git add .
git commit -m "Add JUMP Validator (CLI+GUI) cross-platform package"

gh repo create <YOUR_REPO> --private --source=. --remote=origin --push
```

> Tip: if you keep large example files, consider Git LFS.

---

## Revised scripts

- `jump_validator/validator.py`
- `jump_validator/gui.py`
- `jump_validator/JUMP_l_modules.py`
