"""JUMP Validator GUI

This GUI serves two purposes:

1) Run the validator by selecting a `.params` file.
2) Browse the generated outputs inside the GUI.

The backend algorithm/output is **not** changed: the GUI simply calls
:func:`jump_validator.validator.run`.

Result browsing
---------------
The validator generates multiple plots per entry (PNG/PDF). By default the
pipeline writes these to:

    <out_fol>/ManualValidation/

The GUI scans that folder and shows the plots for the selected entry.

Relative paths
--------------
The legacy scripts were usually executed from the folder containing the params
file, so relative paths inside the params were resolved against that directory.

To preserve that behavior, this GUI offers a checkbox to temporarily `chdir` to
`params.parent` while running. This is ON by default.
"""

from __future__ import annotations

import configparser
import contextlib
import os
import queue
import subprocess
import sys
import threading
import traceback
from pathlib import Path
from typing import Optional

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from .validator import run

# Pillow is optional. If available, we use it for smoother image scaling.
try:  # pragma: no cover
    from PIL import Image, ImageTk  # type: ignore

    _HAVE_PIL = True
except Exception:  # pragma: no cover
    Image = None  # type: ignore
    ImageTk = None  # type: ignore
    _HAVE_PIL = False


class _QueueWriter:
    """File-like object that pushes text into a queue."""

    def __init__(self, q: "queue.Queue[str]") -> None:
        self._q = q

    def write(self, s: str) -> int:
        if s:
            self._q.put(s)
        return len(s)

    def flush(self) -> None:  # pragma: no cover
        return


def _open_folder(path: Path) -> None:
    """Open a folder in the system file explorer."""
    try:
        if sys.platform.startswith("win"):
            os.startfile(str(path))  # type: ignore[attr-defined]
        elif sys.platform == "darwin":
            subprocess.run(["open", str(path)], check=False)
        else:
            subprocess.run(["xdg-open", str(path)], check=False)
    except Exception:
        # Non-fatal.
        pass


def _open_file(path: Path) -> None:
    """Open a file using the system default application."""
    try:
        if sys.platform.startswith("win"):
            os.startfile(str(path))  # type: ignore[attr-defined]
        elif sys.platform == "darwin":
            subprocess.run(["open", str(path)], check=False)
        else:
            subprocess.run(["xdg-open", str(path)], check=False)
    except Exception:
        pass


def _resolve_path_maybe_relative(path_str: str, base_dir: Path) -> Path:
    """Resolve a possibly-relative path against a base directory."""
    p = Path(path_str).expanduser()
    if p.is_absolute():
        return p
    return (base_dir / p).resolve()


def _find_manual_validation_dir(out_fol: Path) -> Optional[Path]:
    """Return the ManualValidation folder if present (best effort)."""
    if not out_fol.exists():
        return None

    mv = out_fol / "ManualValidation"
    if mv.exists() and mv.is_dir():
        return mv

    # Some users may set out_fol directly to ManualValidation.
    if out_fol.name.lower() == "manualvalidation":
        return out_fol

    # Fallback: search one level deep.
    try:
        for child in out_fol.iterdir():
            if child.is_dir() and child.name.lower() == "manualvalidation":
                return child
    except Exception:
        return None

    return None


class ResultsViewer(tk.Toplevel):
    """Browse validator plot outputs per peptide/spectrum."""

    # These suffixes are created by `jump_validator.validator.run()`.
    PLOT_SUFFIXES = {
        "Intensity": "__intensityPlot.png",
        "Publication": "__intensityPlot_no_annotation_pub.png",
        "Tolerance": "__ToleranceMap.png",
    }

    def __init__(self, master: tk.Tk, out_fol: Path, manual_dir: Path) -> None:
        super().__init__(master)
        self.title("JUMP Validator - Results")
        self.geometry("1120x720")

        self.out_fol = out_fol
        self.manual_dir = manual_dir

        # Base filename (without suffix)
        self._entries: list[str] = []
        self._filtered_entries: list[str] = []

        # Keep references to images to avoid garbage-collection.
        self._photo_refs: dict[str, object] = {}

        self._current_base: Optional[str] = None
        self._search_var = tk.StringVar(value="")

        self._build_ui()
        self._load_entries()

    def _build_ui(self) -> None:
        outer = ttk.Frame(self, padding=10)
        outer.pack(fill="both", expand=True)

        # Top bar
        top = ttk.Frame(outer)
        top.pack(fill="x")

        ttk.Label(top, text="ManualValidation:").pack(side="left")
        ttk.Label(top, text=str(self.manual_dir), foreground="#444").pack(side="left", padx=(8, 0))

        ttk.Button(top, text="Open folder", command=lambda: _open_folder(self.manual_dir)).pack(side="right")

        report_path = self.manual_dir / "Report_of_manual_validation.xlsx"
        ttk.Button(
            top,
            text="Open report.xlsx",
            command=lambda p=report_path: _open_file(p)
            if p.exists()
            else messagebox.showinfo("Report", "Report workbook not found."),
        ).pack(side="right", padx=(0, 8))

        # Search
        search_row = ttk.Frame(outer)
        search_row.pack(fill="x", pady=(10, 0))
        ttk.Label(search_row, text="Filter:").pack(side="left")
        ent = ttk.Entry(search_row, textvariable=self._search_var)
        ent.pack(side="left", fill="x", expand=True, padx=(8, 0))
        ent.bind("<KeyRelease>", lambda _e: self._apply_filter())

        # Main split
        main = ttk.Frame(outer)
        main.pack(fill="both", expand=True, pady=(10, 0))

        # Left: entry list
        left = ttk.Frame(main)
        left.pack(side="left", fill="y")

        ttk.Label(left, text="Entries").pack(anchor="w")

        self.listbox = tk.Listbox(left, width=52, height=25)
        self.listbox.pack(side="left", fill="y", expand=False)
        sb = ttk.Scrollbar(left, orient="vertical", command=self.listbox.yview)
        sb.pack(side="right", fill="y")
        self.listbox.configure(yscrollcommand=sb.set)
        self.listbox.bind("<<ListboxSelect>>", self._on_select)

        # Right: plot previews
        right = ttk.Frame(main)
        right.pack(side="left", fill="both", expand=True, padx=(12, 0))

        self.plot_notebook = ttk.Notebook(right)
        self.plot_notebook.pack(fill="both", expand=True)

        self._plot_canvases: dict[str, tk.Canvas] = {}
        self._plot_inners: dict[str, ttk.Frame] = {}

        for tab_name in ["Intensity", "Publication", "Tolerance"]:
            frame = ttk.Frame(self.plot_notebook, padding=6)
            self.plot_notebook.add(frame, text=tab_name)

            canvas = tk.Canvas(frame, highlightthickness=0)
            canvas.pack(side="left", fill="both", expand=True)

            vs = ttk.Scrollbar(frame, orient="vertical", command=canvas.yview)
            hs = ttk.Scrollbar(frame, orient="horizontal", command=canvas.xview)
            vs.pack(side="right", fill="y")
            hs.pack(side="bottom", fill="x")
            canvas.configure(yscrollcommand=vs.set, xscrollcommand=hs.set)

            inner = ttk.Frame(canvas)
            canvas.create_window((0, 0), window=inner, anchor="nw")

            ttk.Label(inner, text="Select an entry to preview plots.").pack(anchor="center", pady=10)

            def _on_inner_configure(event, c=canvas):
                c.configure(scrollregion=c.bbox("all"))

            inner.bind("<Configure>", _on_inner_configure)

            self._plot_canvases[tab_name] = canvas
            self._plot_inners[tab_name] = inner

        # Bottom actions
        bottom = ttk.Frame(outer)
        bottom.pack(fill="x", pady=(10, 0))

        ttk.Button(bottom, text="Open current plot externally", command=self._open_current_plot).pack(side="left")
        ttk.Button(bottom, text="Close", command=self.destroy).pack(side="right")

    def _load_entries(self) -> None:
        """Scan the ManualValidation folder and populate the list."""
        self._entries.clear()

        suffix = self.PLOT_SUFFIXES["Intensity"]
        for p in sorted(self.manual_dir.glob(f"*{suffix}")):
            base = p.name[: -len(suffix)]
            self._entries.append(base)

        self._apply_filter(initial=True)

    def _format_entry_label(self, base: str) -> str:
        """Make a human-readable label from the legacy filename."""
        parts = base.split("__")
        spectrum = parts[0] if parts else ""
        peptide = parts[1] if len(parts) > 1 else ""
        mods = "__".join(parts[2:]) if len(parts) > 2 else ""

        if mods:
            return f"{spectrum} | {peptide} | {mods}"
        if peptide:
            return f"{spectrum} | {peptide}"
        return spectrum or base

    def _apply_filter(self, initial: bool = False) -> None:
        term = self._search_var.get().strip().lower()
        if not term:
            self._filtered_entries = list(self._entries)
        else:
            self._filtered_entries = [b for b in self._entries if term in self._format_entry_label(b).lower()]

        self.listbox.delete(0, "end")
        for base in self._filtered_entries:
            self.listbox.insert("end", self._format_entry_label(base))

        if initial and self._filtered_entries:
            self.listbox.selection_set(0)
            self.listbox.event_generate("<<ListboxSelect>>")

    def _on_select(self, _event=None) -> None:
        if not self.listbox.curselection():
            return
        idx = int(self.listbox.curselection()[0])
        if idx < 0 or idx >= len(self._filtered_entries):
            return

        base = self._filtered_entries[idx]
        self._current_base = base
        self._render_all_plots(base)

    def _open_current_plot(self) -> None:
        base = self._current_base
        if not base:
            return

        tab = self.plot_notebook.tab(self.plot_notebook.select(), "text")
        suffix = self.PLOT_SUFFIXES.get(tab)
        if not suffix:
            return

        path = self.manual_dir / f"{base}{suffix}"
        if path.exists():
            _open_file(path)
        else:
            messagebox.showinfo("Plot", "Plot file not found on disk.")

    def _render_all_plots(self, base: str) -> None:
        for tab_name, suffix in self.PLOT_SUFFIXES.items():
            self._render_plot(tab_name, self.manual_dir / f"{base}{suffix}")

    def _render_plot(self, tab_name: str, path: Path) -> None:
        inner = self._plot_inners[tab_name]
        canvas = self._plot_canvases[tab_name]

        # Clear previous content.
        for child in inner.winfo_children():
            child.destroy()

        if not path.exists():
            ttk.Label(inner, text=f"Missing file: {path.name}").pack(anchor="center", pady=10)
            return

        # Get a target size based on the canvas.
        max_w = max(480, int(canvas.winfo_width()) - 30)
        max_h = max(320, int(canvas.winfo_height()) - 30)

        # Pillow path (best quality)
        if _HAVE_PIL:
            try:
                img = Image.open(path)  # type: ignore[call-arg]
                scale = min(max_w / img.width, max_h / img.height, 1.0)
                if scale < 1.0:
                    img = img.resize((int(img.width * scale), int(img.height * scale)), Image.LANCZOS)  # type: ignore[attr-defined]
                photo = ImageTk.PhotoImage(img)  # type: ignore[call-arg]
                lbl = ttk.Label(inner, image=photo)
                lbl.image = photo  # type: ignore[attr-defined]
                lbl.pack(anchor="nw")
                self._photo_refs[tab_name] = photo
                return
            except Exception:
                # Fall back to Tk loader.
                pass

        # Tk fallback: integer subsampling.
        try:
            photo0 = tk.PhotoImage(file=str(path))
            w, h = photo0.width(), photo0.height()

            factor = 1
            if w > max_w or h > max_h:
                factor = max(1, int(min(w / max_w, h / max_h)))

            photo = photo0.subsample(factor, factor)
            lbl = ttk.Label(inner, image=photo)
            lbl.image = photo  # type: ignore[attr-defined]
            lbl._orig = photo0  # type: ignore[attr-defined]
            lbl.pack(anchor="nw")
            self._photo_refs[tab_name] = photo
        except Exception as e:
            ttk.Label(inner, text=f"Failed to load image: {e}").pack(anchor="center", pady=10)


class ValidatorGUI(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("JUMP Validator")
        self.geometry("820x520")

        self._q: "queue.Queue[str]" = queue.Queue()
        self._worker: Optional[threading.Thread] = None

        self.params_path = tk.StringVar(value="")
        self.use_params_dir_as_cwd = tk.BooleanVar(value=True)

        # CWD when the GUI was launched.
        self._launch_cwd = Path.cwd()

        # Cached output folder from the most recent successful run.
        self._last_out_fol: Optional[Path] = None

        self._build_ui()
        self._poll_log()

    def _build_ui(self) -> None:
        outer = ttk.Frame(self, padding=10)
        outer.pack(fill="both", expand=True)

        # Top controls
        top = ttk.Frame(outer)
        top.pack(fill="x")

        ttk.Label(top, text="Params file:").pack(side="left")
        ttk.Entry(top, textvariable=self.params_path).pack(side="left", fill="x", expand=True, padx=(8, 8))

        ttk.Button(top, text="Browse…", command=self._browse).pack(side="left")
        self.run_btn = ttk.Button(top, text="Run", command=self._start_run, state="disabled")
        self.run_btn.pack(side="left", padx=(8, 0))

        # Options row
        opt = ttk.Frame(outer)
        opt.pack(fill="x", pady=(10, 0))

        ttk.Checkbutton(
            opt,
            text="Use params file folder as working directory (recommended)",
            variable=self.use_params_dir_as_cwd,
        ).pack(side="left")

        self.open_out_btn = ttk.Button(opt, text="Open output folder", command=self._open_output, state="disabled")
        self.open_out_btn.pack(side="right")

        self.view_results_btn = ttk.Button(opt, text="View results", command=self._view_results, state="disabled")
        self.view_results_btn.pack(side="right", padx=(0, 8))

        # Log panel
        log_frame = ttk.LabelFrame(outer, text="Log", padding=6)
        log_frame.pack(fill="both", expand=True, pady=(10, 0))

        self.log_text = tk.Text(log_frame, wrap="word", height=20)
        self.log_text.pack(side="left", fill="both", expand=True)

        scroll = ttk.Scrollbar(log_frame, command=self.log_text.yview)
        scroll.pack(side="right", fill="y")
        self.log_text.configure(yscrollcommand=scroll.set)

        # Bottom bar
        bottom = ttk.Frame(outer)
        bottom.pack(fill="x", pady=(10, 0))

        ttk.Button(bottom, text="Copy log", command=self._copy_log).pack(side="left")
        ttk.Button(bottom, text="Clear", command=self._clear_log).pack(side="left", padx=(8, 0))
        ttk.Button(bottom, text="Quit", command=self.destroy).pack(side="right")

    def _browse(self) -> None:
        path = filedialog.askopenfilename(
            title="Select jump_validator.params",
            filetypes=[("Params files", "*.params"), ("All files", "*")],
        )
        if path:
            self.params_path.set(path)
            self.run_btn.configure(state="normal")
            self.open_out_btn.configure(state="disabled")
            self.view_results_btn.configure(state="disabled")
            self._last_out_fol = None

    def _append_log(self, s: str) -> None:
        self.log_text.insert("end", s)
        self.log_text.see("end")

    def _poll_log(self) -> None:
        try:
            while True:
                s = self._q.get_nowait()
                self._append_log(s)
        except queue.Empty:
            pass
        self.after(100, self._poll_log)

    def _clear_log(self) -> None:
        self.log_text.delete("1.0", "end")

    def _copy_log(self) -> None:
        text = self.log_text.get("1.0", "end").strip()
        self.clipboard_clear()
        self.clipboard_append(text)
        self.update()  # keep clipboard after exit on some platforms

    def _set_running(self, running: bool) -> None:
        self.run_btn.configure(state="disabled" if running else ("normal" if self.params_path.get() else "disabled"))

    def _get_out_folder_resolved(self, params: Path) -> Optional[Path]:
        """Read out_fol and resolve relative paths consistently with the run."""
        try:
            config = configparser.ConfigParser()
            config.read(str(params))
            out_str = config.get("caseEvaluation", "out_fol", fallback="").strip()
            if not out_str:
                return None

            base = params.parent if self.use_params_dir_as_cwd.get() else self._launch_cwd
            return _resolve_path_maybe_relative(out_str, base)
        except Exception:
            return None

    def _open_output(self) -> None:
        params = Path(self.params_path.get()).expanduser()
        out_fol = self._get_out_folder_resolved(params)
        if out_fol and out_fol.exists():
            _open_folder(out_fol)
        else:
            messagebox.showinfo("Output folder", "Output folder not found (check 'out_fol' in params).")

    def _view_results(self) -> None:
        params = Path(self.params_path.get()).expanduser()

        out_fol = self._last_out_fol or self._get_out_folder_resolved(params)
        if not out_fol or not out_fol.exists():
            messagebox.showinfo(
                "Results",
                "Output folder not found. Run the validator first, or check 'out_fol' in params.",
            )
            return

        manual_dir = _find_manual_validation_dir(out_fol)
        if not manual_dir:
            messagebox.showinfo("Results", "ManualValidation folder not found inside output folder.")
            return

        ResultsViewer(self, out_fol=out_fol, manual_dir=manual_dir)

    def _start_run(self) -> None:
        if self._worker and self._worker.is_alive():
            return

        params = Path(self.params_path.get()).expanduser()
        if not params.exists():
            messagebox.showerror("Error", "Params file does not exist.")
            return

        self._set_running(True)
        self.open_out_btn.configure(state="disabled")
        self.view_results_btn.configure(state="disabled")

        # Run in a background thread so the UI stays responsive.
        self._worker = threading.Thread(target=self._run_worker, args=(params,), daemon=True)
        self._worker.start()

    def _run_worker(self, params: Path) -> None:
        writer = _QueueWriter(self._q)

        self._q.put("=== JUMP Validator GUI ===\n")
        self._q.put(f"Params: {params}\n\n")

        old_cwd = Path.cwd()
        try:
            # Optional: temporarily set working dir to params folder (legacy behavior).
            if self.use_params_dir_as_cwd.get():
                os.chdir(str(params.parent))

            with contextlib.redirect_stdout(writer), contextlib.redirect_stderr(writer):
                run(str(params))

            self._q.put("\n✅ Finished successfully.\n")

            out_fol = self._get_out_folder_resolved(params)
            if out_fol:
                self._last_out_fol = out_fol
                self._q.put(f"Output folder (out_fol): {out_fol}\n")

            # Enable post-run buttons.
            self.after(0, lambda: self.open_out_btn.configure(state="normal"))
            self.after(0, lambda: self.view_results_btn.configure(state="normal"))

            self.after(
                0,
                lambda: messagebox.showinfo(
                    "Done",
                    "Validation finished successfully.\n\nClick 'View results' to preview plots, or 'Open output folder' to browse files.",
                ),
            )
        except KeyError as e:
            # Common user error: missing a required key in the params file.
            # We keep the traceback in the log panel for debugging, but show a
            # friendlier message in the dialog.
            tb = traceback.format_exc()
            missing = str(e).strip("\"'")

            hint = ""
            if missing == "ms2_path":
                hint = (
                    "\n\nYour params file is missing `ms2_path` under [caseEvaluation].\n"
                    "- If `ms2_fileType = ms2`, `ms2_path` is required (folder containing fraction subfolders).\n"
                    "- If you are using mzXML instead, set `ms2_fileType = mzXML` and provide `mzXML_path`."
                )
            elif missing in {"mzXML_path", "mzxml_path"}:
                hint = (
                    "\n\nYour params file is missing `mzXML_path` under [caseEvaluation].\n"
                    "- If `ms2_fileType = mzXML`, `mzXML_path` is required (folder containing mzXML files).\n"
                    "- If you are using ms2 instead, set `ms2_fileType = ms2` and provide `ms2_path`."
                )
            elif missing == "ms2_fileType":
                hint = (
                    "\n\nYour params file is missing `ms2_fileType` under [caseEvaluation].\n"
                    "Set it to either `ms2` or `mzXML`."
                )

            self._q.put("\n❌ ERROR\n")
            self._q.put(tb + "\n")
            self.after(
                0,
                lambda: messagebox.showerror(
                    "Error",
                    f"Validation failed (missing key: {missing}).{hint}",
                ),
            )
        except Exception as e:
            tb = traceback.format_exc()
            self._q.put("\n❌ ERROR\n")
            self._q.put(tb + "\n")
            msg = str(e)
            self.after(0, lambda msg=msg: messagebox.showerror("Error", f"Validation failed:\n{msg}"))
        finally:
            try:
                os.chdir(str(old_cwd))
            except Exception:
                pass
            self.after(0, lambda: self._set_running(False))


def main(argv: list[str] | None = None) -> int:  # argv is ignored (GUI)
    """Launch the GUI."""
    app = ValidatorGUI()
    app.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
