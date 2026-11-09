"""Shared look for figures and the one way tables are written.

Figures use one font, one palette and one set of sizes, so the panels of
every figure match. Tables are written twice: as CSV, and as an Excel file in
which every cell is text, so that Excel cannot reformat a value such as
"0.0 (0.0–0.5)" or drop a trailing zero.

Nothing here estimates anything: functions receive finished results and only
format or draw them, rounding where a display needs it. Figures 1 and 2 are
drawn by the accepted shared renderer (settings.RENDERER); the other figures
are drawn by their section scripts.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")                     # draw to files, never open a window
import matplotlib.pyplot as plt          # noqa: E402
import pandas as pd                      # noqa: E402
from openpyxl import Workbook            # noqa: E402
from openpyxl.styles import Alignment, Font  # noqa: E402

import settings as S                     # noqa: E402

INK, MUTED, GRID = "#222222", "#666666", "#E5E5E5"
MEN, WOMEN = "#0072B2", "#E69F00"         # Okabe-Ito blue and orange (colour-blind safe)
GREEN, PURPLE = "#009E73", "#5B4B9A"
MM = 1 / 25.4                             # millimetres to inches

plt.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 9,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.edgecolor": MUTED,
    "axes.linewidth": .65, "text.color": INK, "axes.labelcolor": INK, "xtick.color": MUTED,
    "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
    "svg.fonttype": "none", "pdf.fonttype": 42, "svg.hashsalt": "nhanes",
})


def save_figure(fig, folder: Path, stem: str) -> None:
    """Write an editable SVG, a 600 dpi PNG for the journal and a small preview PNG."""
    svg = folder / f"{stem}.svg"
    fig.savefig(svg, facecolor="white", metadata={"Date": None})
    # Matplotlib leaves trailing spaces in SVG paths; removing them keeps reruns diff-clean.
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text(encoding="utf-8").splitlines()) + "\n",
                   encoding="utf-8")
    fig.savefig(folder / f"{stem}.png", dpi=600, facecolor="white", metadata={"Software": None})
    fig.savefig(folder / f"{stem}_preview.png", dpi=150, facecolor="white", metadata={"Software": None})
    plt.close(fig)


def write_table(table: pd.DataFrame, folder: Path, stem: str) -> None:
    """Write a display table as CSV and as an all-text Excel sheet."""
    table.to_csv(folder / f"{stem}.csv", index=False, encoding="utf-8-sig")
    wb = Workbook()
    ws = wb.active
    ws.append([str(c) for c in table.columns])
    for row in table.itertuples(index=False):
        ws.append(["" if pd.isna(v) else str(v) for v in row])
    for cell in ws[1]:
        cell.font = Font(bold=True)
    for row in ws.iter_rows():
        for cell in row:
            cell.number_format = "@"              # text: Excel keeps the string as written
            cell.alignment = Alignment(vertical="top", wrap_text=True)
    ws.column_dimensions["A"].width = 44
    for letter in "BCDEFG":
        ws.column_dimensions[letter].width = 18
    wb.save(folder / f"{stem}.xlsx")


def publication_renderer():
    """Load the accepted drawing routines shared with the saved-display rebuild.

    The caller supplies its newly computed values. The renderer does not fit
    models or load stored estimates, keeping the drawing layout separate from
    the analysis and avoiding duplicate figure code drifting out of sync.
    """
    import importlib.util
    import sys

    name = "nhanes_current_publication_figures"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, S.RENDERER)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module            # registered, so runs.code_fingerprint includes the renderer
    spec.loader.exec_module(module)
    return module
