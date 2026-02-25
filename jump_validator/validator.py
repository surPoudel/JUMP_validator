"""JUMP Validator

This module contains the *unchanged* core validation pipeline that was previously
implemented in the notebook-derived `JUMP-validator.py` script.

Goals of this refactor:
- Keep the same algorithm and output files.
- Make the code easier to run cross-platform (macOS/Linux/Windows).
- Centralize the pipeline in a callable function (`run`).

The original script was converted from a Jupyter notebook and relied on
`from JUMP_l_modules import *` to pull in helper functions and even some
imports (e.g. `attrgetter`). For backwards compatibility and to avoid any
behavioral changes, this module still imports the helper module at import time.
"""

from __future__ import annotations

import configparser
import os
import pickle
import sys
from pathlib import Path

import matplotlib

# IMPORTANT: Use a non-interactive backend so the tool works on headless systems
# (e.g., servers/CI) and behaves consistently across platforms.
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Third-party scientific parsing
import pyteomics as pyteomics  # noqa: F401  (kept for compatibility with legacy code)
from pyteomics import auxiliary, mzxml, pepxml  # noqa: F401

# Legacy helper functions used throughout the pipeline.
# NOTE: This is intentionally a star import to preserve the legacy namespace.
from .JUMP_l_modules import *  # noqa: F403,F401


def _read_id_file(jump_f_id: str) -> pd.DataFrame:
    """Read the peptide identification (ID) file.

    Supported formats (same as legacy script):
    - JUMP ID text file (semicolon-delimited, with a variable header size)
    - Tab-delimited text file
    - Excel file

    The try/except structure and messages are kept to avoid changing
    user-visible behavior.
    """
    try:
        # Legacy behavior: `return_skiprows` scans the file to find the header.
        df = pd.read_csv(jump_f_id, delimiter=";", skiprows=return_skiprows(jump_f_id, ";", "Peptide"))  # noqa: F405
    except Exception:
        print("The input file is not JUMP idtxt. Trying for ; delimited file")
        try:
            df = pd.read_csv(jump_f_id, delimiter="\t")
        except Exception:
            print(
                "Please check the input file format. Either JUMP id file(separated by ;) "
                "or tab delimited file is accepted. Trying for excel file"
            )
            try:
                df = pd.read_excel(jump_f_id)
            except Exception:
                print("Please make sure your ID file is either excel or Tab delimited file")
                raise

    # Legacy quirk: some JUMP outputs have mis-indexed columns.
    if "IDwDecoy" in jump_f_id:
        df.reset_index(inplace=True)
        df.columns = list(df.columns[1:]) + [str(np.nan)]

    return df


def _select_ms_files(
    input_df: pd.DataFrame,
    ms2_file_type: str,
    ms2_path: str,
    mzxml_path: str,
) -> list[str]:
    """Select the MS/MS files required to extract spectra.

    This preserves the legacy behavior:
    - For MZXML mode: prefer <mzxml_path>/<sample>/*.mzXML, then fallback to <mzxml_path>/<sample>.mzXML
    - For MS2 mode: prefer <ms2_path>/<sample>/*.ms2*

    Returns a list of file paths.
    """
    mzxml_list: list[str] = []
    samples = list(set(input_df.Exp))

    for sample in samples:
        if ms2_file_type.upper() == "MZXML":
            try:
                mzxml_file = glob.glob(str(Path(mzxml_path) / sample / "*.mzXML"))[0]  # noqa: F405
            except Exception:
                print(
                    "The mzxml file is not present inside the fraction subfolder of given path. "
                    "Trying for directly checking mzXML in the given path"
                )
                try:
                    mzxml_file = glob.glob(str(Path(mzxml_path) / f"{sample}.mzXML"))[0]  # noqa: F405
                except Exception:
                    print("No mzXML file present in the given path")
                    raise
        else:
            mzxml_file = glob.glob(str(Path(ms2_path) / sample / "*.ms2*"))[0]  # noqa: F405

        mzxml_list.append(mzxml_file)

    return mzxml_list


def _intlike_to_str(v) -> str:
    """Normalize integer-like values to clean strings.

    This mirrors the normalization used when reading the peptide/spectrum list,
    and helps with spectrum-key matching when scan/charge were parsed as floats.
    """

    if pd.isna(v):
        return ""
    try:
        if isinstance(v, (int, np.integer)):
            return str(int(v))
        if isinstance(v, (float, np.floating)):
            if float(v).is_integer():
                return str(int(v))
            return str(v)
    except Exception:
        pass
    s = str(v).strip()
    # Common Excel artifact: 12345.0
    if s.endswith('.0') and s.replace('.', '', 1).isdigit():
        try:
            f = float(s)
            if f.is_integer():
                return str(int(f))
        except Exception:
            pass
    return s


def _spectrum_from_columns(df: pd.DataFrame, run_col: str, scan_col: str, charge_col: str) -> pd.Series:
    run = df[run_col].astype(str).str.strip()
    scan = df[scan_col].apply(_intlike_to_str)
    charge = df[charge_col].apply(_intlike_to_str)
    return run + "." + scan.astype(str) + "." + charge.astype(str)


def _spectrum_from_outfile_series(out_series: pd.Series, which_scan: str = 'first') -> pd.Series:
    """Robustly derive run/scan/charge from an Outfile-like string.

    JUMP/SEQUEST-style outfiles often look like:
      run.scan.scan.charge.dta
    but variants exist. This parser is only used as a fallback when
    legacy parsing yields zero matches.
    """

    import re

    rx4 = re.compile(r"^(?P<run>.+?)\.(?P<scan1>\d+)\.(?P<scan2>\d+)\.(?P<charge>\d+)(?:\..*)?$")
    rx3 = re.compile(r"^(?P<run>.+?)\.(?P<scan>\d+)\.(?P<charge>\d+)(?:\..*)?$")

    def parse_one(v):
        if pd.isna(v):
            return None
        base = os.path.basename(str(v).replace('\\', '/')).strip()
        m = rx4.match(base)
        if m:
            scan = m.group('scan1') if which_scan == 'first' else m.group('scan2')
            return f"{m.group('run')}.{scan}.{m.group('charge')}"
        m = rx3.match(base)
        if m:
            return f"{m.group('run')}.{m.group('scan')}.{m.group('charge')}"
        return None

    return out_series.apply(parse_one)


def _maybe_fix_spectrum_keys(df_pho2: pd.DataFrame, input_spectra: set[str]) -> pd.DataFrame:
    """If legacy spectrum construction yields zero matches, try alternatives.

    This is a workflow robustness improvement only. When legacy parsing already
    yields >=1 match, we keep it unchanged to avoid altering behavior.
    """

    legacy_matches = int(df_pho2['spectrum'].isin(input_spectra).sum())
    if legacy_matches > 0:
        return df_pho2

    candidates: list[tuple[str, pd.Series, int]] = []

    cols = set(df_pho2.columns)

    # 1) Column-based spectrum keys (common variants across JUMP exports).
    run_cols = ["Run#", "Run", "Exp", "exp"]
    scan_cols = ["Scan#", "Scan", "scan"]
    charge_cols = ["z", "Charge", "charge", "Z"]

    for r in run_cols:
        for s in scan_cols:
            for c in charge_cols:
                if {r, s, c}.issubset(cols):
                    spec = _spectrum_from_columns(df_pho2, r, s, c)
                    candidates.append((f"{r}/{s}/{c}", spec, int(spec.isin(input_spectra).sum())))


    # 3) Robust Outfile parsing (both scan positions).
    if 'Outfile' in cols:
        spec_first = _spectrum_from_outfile_series(df_pho2['Outfile'], which_scan='first')
        spec_second = _spectrum_from_outfile_series(df_pho2['Outfile'], which_scan='second')

        # Keep original spectrum for rows we couldn't parse.
        spec_first = spec_first.fillna(df_pho2['spectrum'])
        spec_second = spec_second.fillna(df_pho2['spectrum'])

        candidates.append(('Outfile (scan1)', spec_first, int(spec_first.isin(input_spectra).sum())))
        candidates.append(('Outfile (scan2)', spec_second, int(spec_second.isin(input_spectra).sum())))

    if not candidates:
        return df_pho2

    best_name, best_series, best_count = max(candidates, key=lambda t: t[2])
    if best_count > 0:
        print(f"No matches using legacy spectrum keys; switching to: {best_name} (matched {best_count} spectra).")
        df_pho2 = df_pho2.copy()
        df_pho2['spectrum'] = best_series.astype(str)

    return df_pho2


def _maybe_align_input_spectra(input_df: pd.DataFrame, id_df: pd.DataFrame) -> pd.DataFrame:
    """Align peptideFile spectra keys to the ID file when there are zero matches.

    In some workflows, multiple mzXML files are concatenated into a single file.
    Depending on the concatenation method, scan numbers in one source may be
    shifted by a constant *offset* per run/fraction (e.g., second file scans are
    incremented by the number of scans in the first file). When the peptide/spectrum
    list was generated using the *pre*-concatenation scan numbers, but the ID file
    (and extracted spectra) use the *post*-concatenation scan numbers, the strict
    spectrum key intersection becomes empty.

    This helper only runs when there are **zero** matches after all ID-side spectrum
    parsing strategies have been attempted. It tries to infer an integer scan offset
    per run by comparing scan lists, applies it to the peptideFile scans, and then
    (only when unambiguous) reconciles charge based on the ID file.

    This is a workflow robustness improvement and does **not** alter behavior when
    inputs already match.
    """

    try:
        input_specs = set(input_df["spectrum"].astype(str))
        id_specs = set(id_df["spectrum"].astype(str))
    except Exception:
        return input_df

    # If anything matches already, do nothing.
    if input_specs & id_specs:
        return input_df

    # Require the expected columns from prepareInput().
    if not {"Exp", "scan", "charge"}.issubset(set(input_df.columns)):
        return input_df

    # Parse ID spectra into run/scan/charge using rsplit so run names can contain dots.
    id_parts = id_df["spectrum"].astype(str).str.rsplit(".", n=2, expand=True)
    if id_parts.shape[1] != 3:
        return input_df

    id_runs = id_parts[0].astype(str).str.strip()
    id_scan = pd.to_numeric(id_parts[1].apply(_intlike_to_str), errors="coerce")
    id_charge = pd.to_numeric(id_parts[2].apply(_intlike_to_str), errors="coerce")

    from collections import Counter, defaultdict

    run_to_scans: dict[str, np.ndarray] = {}
    runscan_to_charges: dict[tuple[str, int], set[int]] = defaultdict(set)

    # Build scan arrays per run for fast nearest-neighbor lookups.
    tmp: dict[str, list[int]] = defaultdict(list)
    for run, scan, ch in zip(id_runs, id_scan, id_charge):
        if pd.isna(scan):
            continue
        s = int(scan)
        tmp[run].append(s)
        if not pd.isna(ch):
            runscan_to_charges[(run, s)].add(int(ch))

    for run, scans in tmp.items():
        run_to_scans[run] = np.array(sorted(set(scans)), dtype=int)

    # Prepare a working copy of the input file with integer scans.
    aligned = input_df.copy()
    aligned["_run"] = aligned["Exp"].astype(str).str.strip()
    aligned["_scan_i"] = pd.to_numeric(aligned["scan"].apply(_intlike_to_str), errors="coerce")
    aligned["_charge_i"] = pd.to_numeric(aligned["charge"].apply(_intlike_to_str), errors="coerce")

    offsets: dict[str, int] = {}

    # Infer an offset per run by looking for the most common shift between
    # peptideFile scans and the nearest ID scan within that run.
    for run in sorted(set(aligned["_run"])):

        id_scans = run_to_scans.get(run)
        if id_scans is None or id_scans.size == 0:
            continue

        pep_scans = aligned.loc[aligned["_run"] == run, "_scan_i"].dropna().astype(int).tolist()
        if not pep_scans:
            continue

        diffs: list[int] = []
        for s in pep_scans:
            idx = int(np.searchsorted(id_scans, s))
            cand: list[int] = []
            if idx > 0:
                cand.append(int(id_scans[idx - 1]))
            if idx < id_scans.size:
                cand.append(int(id_scans[idx]))
            if not cand:
                continue
            nearest = min(cand, key=lambda x: abs(x - s))
            diffs.append(int(nearest - s))

        if not diffs:
            continue

        cnt = Counter(diffs)
        # Evaluate the top few candidate offsets by how many scans become exact matches.
        id_set = set(map(int, id_scans.tolist()))
        best_off = None
        best_matches = -1

        for off, _ in cnt.most_common(5):
            m = sum((s + int(off)) in id_set for s in pep_scans)
            if m > best_matches:
                best_matches = m
                best_off = int(off)

        if best_off is not None and best_matches > 0 and best_off != 0:
            offsets[run] = best_off

    if offsets:
        print("No exact spectrum matches found; attempting scan-offset alignment (common after mzXML concatenation).")

        # Apply per-run offsets.
        for run, off in offsets.items():
            mask = aligned["_run"] == run

            aligned.loc[mask, "_scan_i"] = aligned.loc[mask, "_scan_i"] + off

            # Report how many scans are now present in the ID set for that run.
            id_set = set(map(int, run_to_scans[run].tolist()))
            after = aligned.loc[mask, "_scan_i"].dropna().astype(int).tolist()
            matches = sum(int(s) in id_set for s in after)
            denom = max(1, len(after))
            print(f"  {run}: inferred scan offset {off:+d} (matches {matches}/{denom} scans after shift)")

        # Update scan strings from the shifted integer scans.
        aligned["scan"] = aligned["_scan_i"].apply(lambda v: "" if pd.isna(v) else str(int(v)))

        # Rebuild spectrum with the *current* charge column.
        aligned["charge"] = aligned["charge"].apply(_intlike_to_str)
        aligned["spectrum"] = aligned["_run"] + "." + aligned["scan"].astype(str) + "." + aligned["charge"].astype(str)
    # Try reconciling charge when unambiguous (unique charge per run/scan in the ID file).
    aligned_specs = set(aligned["spectrum"].astype(str))
    current_matches = len(aligned_specs & id_specs)

    aligned_charge = aligned.copy()
    fixes = 0
    for i, row in aligned_charge.iterrows():
        run = str(row["_run"]).strip()
        scan_i = row["_scan_i"]
        if pd.isna(scan_i):
            continue

        charges = runscan_to_charges.get((run, int(scan_i)))
        if not charges:
            continue

        if len(charges) == 1:
            ch = str(next(iter(charges)))
            if str(row["charge"]).strip() != ch:
                aligned_charge.at[i, "charge"] = ch
                fixes += 1

    if fixes:
        aligned_charge["charge"] = aligned_charge["charge"].apply(_intlike_to_str)
        aligned_charge["spectrum"] = (
            aligned_charge["_run"] + "." + aligned_charge["scan"].astype(str) + "." + aligned_charge["charge"].astype(str)
        )
        new_matches = len(set(aligned_charge["spectrum"].astype(str)) & id_specs)
        if new_matches > current_matches:
            print(f"  Adjusted charge for {fixes} spectra based on ID file (unique charge per run/scan).")
            aligned = aligned_charge

    # Clean up helper columns.
    for c in ["_run", "_scan_i", "_charge_i"]:
        if c in aligned.columns:
            aligned.drop(columns=[c], inplace=True)

    # Keep the aligned version only if we gained at least one match.
    if set(aligned["spectrum"].astype(str)) & id_specs:
        return aligned

    return input_df



def run(params_file: str) -> None:
    """Run the JUMP spectral validation pipeline.

    Parameters
    ----------
    params_file:
        Path to a `.params` INI-style configuration file.

    Output
    ------
    The function generates the same output folder structure and report
    as the legacy script.
    """

    config = configparser.ConfigParser()
    config.read(params_file)

    case = config["caseEvaluation"]

    peptideFile = case["peptideFile"]
    jump_f_id = case["jump_f_id"]
    jumpl = case["jumpl"]

    # NOTE:
    # The legacy params file includes BOTH `ms2_path` and `mzXML_path`, even though
    # only one is needed depending on `ms2_fileType`.
    #
    # Users sometimes omit the unused path. In the first GUI refactor we accessed
    # both unconditionally, which caused a KeyError even when the missing path
    # wasn't required. To preserve the algorithm and outputs while improving
    # usability, we now only require the relevant key.
    ms2_fileType = case["ms2_fileType"]
    if ms2_fileType.upper() == "MZXML":
        # Required for mzXML mode; MS2 path is optional.
        mzXML_path = case.get("mzXML_path", fallback="").strip()
        if not mzXML_path:
            # Keep the same exception type users would have seen when a required
            # key is missing.
            raise KeyError("mzXML_path")
        ms2_path = case.get("ms2_path", fallback="").strip()
    else:
        # Required for MS2 mode; mzXML path is optional.
        ms2_path = case.get("ms2_path", fallback="").strip()
        if not ms2_path:
            raise KeyError("ms2_path")
        mzXML_path = case.get("mzXML_path", fallback="").strip()

    ion_types_str = case["ion_types"]
    ionLoss_str = case["ionLoss"]
    tol = case["tol"]
    out_fol = case["out_fol"]

    # Example pepxml file to parse the modification information
    pepxml_file = case["pepxml"]

    ion_types = ion_types_str.split(",")
    ionLoss = ionLoss_str.split(",")

    # 1) Read the ID file (JUMP output / tab / excel)
    df_pho22 = _read_id_file(jump_f_id)
    df_pho2 = df_pho22.copy()

    # 2) Read the list of spectra to validate (excel/tab)
    inputCheck = peptideFile
    inputDF = prepareInput(inputCheck)  # noqa: F405

    # 3) Create the JUMP spectrum key and normalize peptide columns based on jumpl mode
    df_pho2["spectrum"] = df_pho2.apply(createOutfile, df=df_pho2, axis=1)  # noqa: F405

    # If the spectrum keys derived from the ID file do not match the
    # peptide/spectrum list (common when Outfile formatting varies),
    # try a few alternative strategies. This is only applied when
    # legacy parsing yields *zero* matches to avoid changing behavior.
    input_spectra = set(inputDF.spectrum.astype(str))
    df_pho2 = _maybe_fix_spectrum_keys(df_pho2, input_spectra)

    # If there are still zero matches, try aligning the peptideFile scan/charge.
    # This commonly occurs when mzXML files were concatenated and scan numbers
    # in the peptide list refer to the pre-concatenation scan indices.
    if int(df_pho2["spectrum"].isin(input_spectra).sum()) == 0:
        inputDF = _maybe_align_input_spectra(inputDF, df_pho2)
        input_spectra = set(inputDF.spectrum.astype(str))
        matched_after = int(df_pho2["spectrum"].isin(input_spectra).sum())
        if matched_after > 0:
            print(f"Matched {matched_after} spectra after aligning peptideFile scan/charge.")


    if jumpl.upper() == "YES":
        df_pho2.rename(columns={"Peptide": "JUMP-f_Peptide", "JUMPl_site": "Peptides"}, inplace=True)
        scoreDict = dict(zip(df_pho2.spectrum, df_pho2.JUMPl_score))
    else:
        df_pho2.rename(columns={"Peptide": "Peptides"}, inplace=True)
        scoreDict = {}

    df_pho2[["exp", "scan", "charge"]] = df_pho2["spectrum"].str.split(".", expand=True)
    df_pho2["JUMPl_score"] = df_pho2.spectrum.map(scoreDict)

    df_pho2.drop_duplicates(subset=["spectrum", "Peptides"], inplace=True, keep="first")

    df_pho3 = df_pho2.loc[df_pho2.spectrum.isin(list(inputDF.spectrum))]
    df_pho = df_pho3.copy()

    # If no spectra match, later steps fail with confusing pandas errors.
    # This guard does not change the algorithm; it simply provides a clearer
    # message when the ID file and peptide list refer to different runs/scans
    # (or when scan/charge formatting differs, e.g., 12345 vs 12345.0).
    if df_pho.empty:
        try:
            input_examples = list(inputDF.spectrum.head(5))
        except Exception:
            input_examples = []
        try:
            id_examples = list(df_pho2.spectrum.head(5))
        except Exception:
            id_examples = []

        print("No matching spectra were found between the ID file and the peptide/spectrum list.")
        if input_examples:
            print("Example spectra from peptideFile:", input_examples)
        if id_examples:
            print("Example spectra from ID file:", id_examples)

        raise ValueError(
            "No matching spectra found. Check that peptideFile Exp/scan/charge match the ID file run/scan/charge. "
            "A common issue is scan/charge being read as floats (e.g., 12345.0)."
        )
    df_pho["expScan"] = df_pho.exp + "." + df_pho.scan.astype("str")

    # 4) Collect the MS2 / mzXML files needed
    mzXML_list = _select_ms_files(inputDF, ms2_fileType, ms2_path, mzXML_path)

    # 5) Extract MS2 spectra into a dataframe
    if ms2_fileType.upper() == "MZXML":
        testID_mzXML = mzmlDF(mzXML_list, list(set(df_pho["expScan"])))  # noqa: F405
        mz_spectrumDict = dict(zip(testID_mzXML.expScan, testID_mzXML["m/z"]))
        int_spectrumDict = dict(zip(testID_mzXML.expScan, testID_mzXML["intensity"]))

        df_pho["m/z"] = df_pho["expScan"].map(mz_spectrumDict)
        df_pho["intensity"] = df_pho["expScan"].map(int_spectrumDict)
        ms2DF = df_pho[["spectrum", "exp", "scan", "charge", "expScan", "m/z", "intensity", "JUMPl_score"]]
    else:
        ms2DF = msToDF(mzXML_list, list(inputDF.spectrum))  # noqa: F405
        ms2DF["JUMPl_score"] = ms2DF.spectrum.map(scoreDict)

    # 6) Parse modification info from pepXML
    jump_modAA_dict, jump_mod_dict, sta_AA = getDynStatModsInfoPepXml(pepxml_file)  # noqa: F405
    print(jump_modAA_dict, jump_mod_dict, sta_AA)

    df_pho[["plain_peptide", "modifications"]] = df_pho.apply(
        computeModifications,  # noqa: F405
        sta_AA=sta_AA,
        jump_mod_dict=jump_mod_dict,
        axis=1,
    )

    reqdCols = ["spectrum", "plain_peptide", "modifications", "XCorr", "JUMPl_score"]
    df_pho2 = df_pho[reqdCols]
    inputFileDf = df_pho2.copy()

    # NOTE: This string is used later to support cases where a single spectrum is
    # evaluated against multiple hypothetical modification patterns.
    inputFileDf["spectrum_modifications"] = inputFileDf.spectrum + "_" + inputFileDf.modifications
    inputFileDf["combined"] = inputFileDf.apply(lambda row: "\t".join(row.values[0:3].astype(str)), axis=1)

    combinedDict = dict(zip(inputFileDf.spectrum_modifications, inputFileDf.combined))
    plainPepDict = dict(zip(inputFileDf.spectrum_modifications, inputFileDf.plain_peptide))
    modSiteDict = dict(zip(inputFileDf.spectrum_modifications, inputFileDf.modifications))

    # 7) Merge spectra with peptide/mod information
    ms2DF_2 = ms2DF.merge(inputFileDf, how="inner", on="spectrum")
    ms2DF_final = ms2DF_2.copy()

    ms2DF2 = ms2DF_final.dropna(subset=["spectrum_modifications"]).reset_index()
    newDF = ms2DF2.rename(columns={"m/z": "exp_mz_list", "intensity": "intensity_list"})

    # 8) Prepare output folders/files
    mkdir(out_fol)  # noqa: F405
    newDir = str(Path(out_fol) / "ManualValidation")
    mkdir(newDir)  # noqa: F405

    workbookName = str(Path(newDir) / "Report_of_manual_validation.xlsx")
    writer = pd.ExcelWriter(workbookName, engine="xlsxwriter")

    updater = 1
    msviewer_input_dict: dict[str, list[list[list[float]]]] = {}

    # 9) Process each peptide/spectrum and write to the report
    for index_sp, specs in enumerate(list(inputFileDf["combined"])):
        xcorr = "missing"
        prob = "missing"
        lscoreSite = "missing"
        absDelRT = "missing"

        spectrumSplit = specs.split("\t")

        # NOTE: Legacy file naming uses a double-underscore join of the combined fields.
        filename = str(Path(newDir) / ("__".join(spectrumSplit)))

        if "XCorr" in inputFileDf.columns:
            xcorr = inputFileDf.loc[inputFileDf["combined"] == specs].XCorr.values[0]

        if "JUMPl_score" in inputFileDf.columns:
            lscoreSite = inputFileDf.loc[inputFileDf["combined"] == specs].JUMPl_score.values[0]

        spectrum_DF = newDF.loc[newDF.combined == specs]

        peptide = spectrum_DF.loc[spectrum_DF.combined == specs].plain_peptide.values[0]
        modsOri = spectrum_DF.loc[spectrum_DF.combined == specs].modifications.values[0]
        maxCharge = int(spectrum_DF.loc[spectrum_DF.combined == specs].charge.values[0])

        mods = modsForReport(modsOri, peptide)  # noqa: F405
        massPosDict1 = spectrumToDict(spectrumSplit[-1])  # noqa: F405

        # Generate theoretical fragment series
        if maxCharge == 1:
            df_pep = ionSeriesIonLossSpeRes(peptide, maxcharge=2, massPosDict=massPosDict1, useMod="Yes")  # noqa: F405
        else:
            df_pep = ionSeriesIonLossSpeRes(
                peptide,
                maxcharge=maxCharge,
                massPosDict=massPosDict1,
                useMod="Yes",
            )  # noqa: F405

        # Select columns requested by ion types / neutral losses
        reqdCols2 = ["Seq"]
        for x in df_pep.columns:
            for y in ion_types:
                if y + "+" in x:
                    reqdCols2.append(x)
                for z in ionLoss:
                    if y + "-" + z in x:
                        reqdCols2.append(x)

        df_pep2 = df_pep.copy()[reqdCols2]
        df_pep2["b-series"] = [*range(1, len(df_pep2) + 1, 1)]
        df_pep2["y-series"] = [*range(len(df_pep2), 0, -1)]

        exp_mz_list = spectrum_DF.exp_mz_list.values[0]
        exp_int_list = spectrum_DF.intensity_list.values[0]

        match_list, match_int_list, true_ions_ppm_list = ionMatches(  # noqa: F405
            exp_mz_list,
            exp_int_list,
            df_pep,
            ion_types=ion_types,
            ionLoss=ionLoss,
            tol=int(tol),
        )

        msviewer_input_dict[specs] = [[exp_mz_list, exp_int_list]]

        df_pep2["Peptide_Mod_Seq"] = df_pep2.apply(displaySeqTable, massPosDict=massPosDict1, axis=1)  # noqa: F405

        # Column order for display table
        displayTableCols: list[str] = []
        for cols in df_pep2.columns:
            if "b+" in cols and "series" not in cols:
                displayTableCols.append(cols)

        displayTableCols = displayTableCols[::-1]
        displayTableCols.append("b-series")
        displayTableCols.append("Peptide_Mod_Seq")
        displayTableCols.append("y-series")

        for cols in df_pep2.columns:
            if "y+" in cols and "series" not in cols:
                displayTableCols.append(cols)

        df = df_pep2.copy()[displayTableCols]

        # Create human-readable ion labels for matched peaks
        seriesName: list[str] = []
        for vals in match_list:
            val = float(vals)
            row_column_pair = df_pep2[df_pep2.isin([val])].stack().index[0]
            row = row_column_pair[0]
            column = row_column_pair[1]

            if "-" in column:
                columnSplit = column.split("-")
                if ("x" in column) or ("y" in column) or ("z" in column):
                    ionNu = columnSplit[0] + str(df_pep2["y-series"][row]) + "-" + columnSplit[1]
                if ("a" in column) or ("b" in column) or ("c" in column):
                    ionNu = columnSplit[0] + str(df_pep2["b-series"][row]) + "-" + columnSplit[1]
            else:
                columnSplit = column.split("+")
                if ("x" in column) or ("y" in column) or ("z" in column):
                    ionNu = columnSplit[0] + str(df_pep2["y-series"][row]) + "+" + columnSplit[1]
                if ("a" in column) or ("b" in column) or ("c" in column):
                    ionNu = columnSplit[0] + str(df_pep2["b-series"][row]) + "+" + columnSplit[1]

            seriesName.append(ionNu)

        seriesName = list(map(lambda x: x.replace("+1", ""), seriesName))

        # -------------------------
        # Plot A: Spectrum (matched ions highlighted)
        # -------------------------
        fig, ax = plt.subplots(figsize=(7, 2.5))
        plt.style.use("ggplot")

        plt.rcParams["axes.edgecolor"] = "#010101"
        plt.rcParams["axes.facecolor"] = "#FFFFFF"

        plt.rcParams.update({"font.size": 8, "figure.max_open_warning": 0})

        plt.bar(
            list(spectrum_DF.exp_mz_list)[0],
            list(spectrum_DF.intensity_list)[0],
            width=0.1,
            linewidth=0.5,
            edgecolor="black",
        )

        plt.bar(
            match_list,
            match_int_list,
            width=0.1,
            linewidth=0.5,
            edgecolor="red",
            alpha=0.7,
        )

        ax.set_ylabel("Absolute Intensity", color="black")
        ax.set_xlabel("m/z", color="black")
        ax.set_xlim(0, max(match_list) + 100)
        ax.tick_params(axis="x", colors="black")
        ax.tick_params(axis="y", colors="black")
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

        for i, txt in enumerate(seriesName):
            for label in ["b", "y"]:
                if label in txt and "-" not in txt:
                    plt.annotate(txt, (match_list[i], match_int_list[i]), ha="center")

        figurename = filename + "__intensityPlot.pdf"
        figurename1 = filename + "__intensityPlot.png"
        plt.savefig(figurename, bbox_inches="tight", dpi=600)
        plt.savefig(figurename1, bbox_inches="tight", dpi=600)

        # -------------------------
        # Plot B: Publication-style spectrum (no annotations, normalized)
        # -------------------------
        fig, ax = plt.subplots(figsize=(7, 2.5))
        plt.style.use("ggplot")

        plt.rcParams["axes.edgecolor"] = "#010101"
        plt.rcParams["axes.facecolor"] = "#FFFFFF"

        plt.rcParams.update({"font.size": 8, "figure.max_open_warning": 0})

        new_int_list = list(spectrum_DF.intensity_list)[0]
        max_value = max(new_int_list)
        normalized_int = normalize_intensity(new_int_list, max_value)  # noqa: F405

        plt.bar(
            list(spectrum_DF.exp_mz_list)[0],
            normalized_int,
            width=0.1,
            linewidth=0.5,
            edgecolor="black",
        )

        normalized_int_matched = normalize_intensity(match_int_list, max_value)  # noqa: F405

        plt.bar(
            match_list,
            normalized_int_matched,
            width=0.1,
            linewidth=0.5,
            edgecolor="black",
            alpha=0.7,
        )

        ax.set_ylabel("Relative Abundance", color="black")
        ax.set_xlabel("m/z", color="black")
        ax.set_xlim((min(list(spectrum_DF.exp_mz_list)[0]) - 50), max(match_list) + 100)
        ax.tick_params(axis="x", colors="black")
        ax.tick_params(axis="y", colors="black")
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

        figurename_noannot2 = filename + "__intensityPlot_no_annotation_pub.pdf"
        figurename_noannot = filename + "__intensityPlot_no_annotation_pub.png"
        plt.savefig(figurename_noannot2, bbox_inches="tight", dpi=600)
        plt.savefig(figurename_noannot, bbox_inches="tight", dpi=600)

        # -------------------------
        # Plot C: Tolerance map (ppm vs m/z)
        # -------------------------
        ionDFmatch = pd.DataFrame(columns=["matched_ions", "matched_intensity", "matched_ppm"])
        ionDFmatch["matched_ions"] = match_list
        ionDFmatch["matched_intensity"] = match_int_list
        ionDFmatch["matched_ppm"] = true_ions_ppm_list

        bins = [0, 10000, 100000, 1000000, 10000000, 100000000, 1000000000, 10000000000]
        ionDFmatch["intRange"] = pd.cut(ionDFmatch.matched_intensity, bins=bins, include_lowest=True)
        ionDFmatch["right"] = ionDFmatch["intRange"].map(attrgetter("right"))  # noqa: F405
        ionDFmatch["intensity"] = ionDFmatch.apply(rightInterval, axis=1)  # noqa: F405
        ionDFmatch["ionType"] = seriesName

        fig, ax = plt.subplots(figsize=(6, 2.5))
        plt.style.use("ggplot")

        plt.rcParams["axes.edgecolor"] = "#010101"
        plt.rcParams["axes.facecolor"] = "#FFFFFF"

        plt.rcParams.update({"font.size": 8, "figure.max_open_warning": 0})

        ionDFmatch2 = ionDFmatch.sort_values(by="right", ascending=True)
        ionDFmatch2.to_excel(str(Path(newDir) / "{}_match_info.xlsx".format(spectrumSplit[0])))

        sns.scatterplot(data=ionDFmatch2, x="matched_ions", y="matched_ppm", s=25, hue="intensity")

        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)
        plt.hlines(0, 50, 2200, color="green", linestyles="dotted", linewidth=0.5)

        ax.set_ylabel("mass error (ppm)", color="black")
        ax.set_xlabel("all matched ions", color="black")
        ax.set_xlim(0, max(match_list) + 100)
        ax.tick_params(axis="x", colors="black")
        ax.tick_params(axis="y", colors="black")
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

        for i, txt in enumerate(seriesName):
            for label in ["b", "y"]:
                if label in txt and "-" not in txt:
                    plt.annotate(
                        txt,
                        (list(ionDFmatch.matched_ions)[i], list(ionDFmatch.matched_ppm)[i] + 0.3),
                        ha="center",
                    )

        figurename_tolmap_pdf = filename + "__ToleranceMap.pdf"
        figurename_tolmap_png = filename + "__ToleranceMap.png"

        plt.savefig(figurename_tolmap_pdf, bbox_inches="tight", dpi=600)
        plt.savefig(figurename_tolmap_png, bbox_inches="tight", dpi=600)

        matched = str(len(ionDFmatch)) + "/" + str(len(list(spectrum_DF.exp_mz_list)[0]))

        # -------------------------
        # Prepare the formatted fragment table and write to Excel
        # -------------------------
        finalDF, seqSymbol = reformat_dataframe2(df, massPosDict1, jump_mod_dict, sta_AA, match_list)  # noqa: F405
        finalDF2 = finalDF.style.applymap(lambda x: "color: red" if x in match_list else "color: black")

        massSeriesLength = df.shape[0]

        try:
            lscoreSite = lscoreSite.astype("str")
            xcorr = xcorr.astype("str")
            prob = prob.astype("str")
            absDelRT = absDelRT.asype("str")
        except Exception:
            print("String has no astype! Warning")

        writer, updater = excelWriter3(  # noqa: F405
            writer,
            finalDF2,
            "Sheet1",
            figurename1,
            spectrumSplit,
            xcorr,
            massSeriesLength,
            peptide,
            updater,
            index_sp,
        )

    # Save the workbook
    writer.save()

    # Save msviewer input as a pickle for downstream tooling
    filenamepkl = str(Path(newDir) / "all_ms2_input.pkl")
    print("...Finalizing the ms2 peaks for msviewer input generation and saving as a pickle file")

    with open(filenamepkl, "wb") as handle:
        pickle.dump(msviewer_input_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
