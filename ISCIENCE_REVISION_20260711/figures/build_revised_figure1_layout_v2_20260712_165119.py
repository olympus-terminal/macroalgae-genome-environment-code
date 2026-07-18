#!/usr/bin/env python3
"""Compose Figure 1 layout V4 with the revised map and original panels B–I.

The map is inserted as one panel A. The manually assembled specimen and
habitat panels are copied directly from the original PDF without alteration.
"""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime
from pathlib import Path

import fitz


SCRIPT = Path(__file__).resolve()
ROOT = SCRIPT.parents[2]
OUTDIR = SCRIPT.parent
GENERATION_ID = "20260712_180210"
MAP_PDF = (
    ROOT
    / "ISCIENCE_REVISION_20260711/spatial/"
    / "Figure_global_sampling_map_layout_v4_20260712_180210.pdf"
)
MAP_MANIFEST = (
    ROOT
    / "ISCIENCE_REVISION_20260711/spatial/"
    / "sampling_map_layout_v4_manifest_20260712_180210.json"
)
ORIGINAL_FIGURE = ROOT / "FIGURES/FIGURE_1-25SEP2025-maac2-drn.pdf"
OUT_PDF = OUTDIR / f"Figure1_revised_map_and_specimens_layout_v4_{GENERATION_ID}.pdf"
OUT_SVG = OUTDIR / f"Figure1_revised_map_and_specimens_layout_v4_{GENERATION_ID}.svg"
OUT_PNG = OUTDIR / f"Figure1_revised_map_and_specimens_layout_v4_{GENERATION_ID}.png"
MANIFEST = OUTDIR / f"Figure1_revised_layout_v4_integrity_{GENERATION_ID}.json"

DPI = 300
DOCX_ASPECT_WIDTH = 1028
DOCX_ASPECT_HEIGHT = 1348
ORIGINAL_CROP_TOP_PT = 338.0
VISIBLE_GAP_BEFORE_PT = 51.75
GAP_REDUCTION_FRACTION = 0.75
LOWER_TRANSLATE_UP_PT = VISIBLE_GAP_BEFORE_PT * GAP_REDUCTION_FRACTION


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    inputs = (MAP_PDF, MAP_MANIFEST, ORIGINAL_FIGURE)
    for path in inputs:
        if not path.is_file():
            raise FileNotFoundError(path)
    if any(path.exists() for path in (OUT_PDF, OUT_SVG, OUT_PNG, MANIFEST)):
        raise FileExistsError("Refusing to overwrite an existing Figure 1 layout-v4 output")

    with fitz.open(MAP_PDF) as map_document, fitz.open(ORIGINAL_FIGURE) as original_document:
        if map_document.page_count != 1 or original_document.page_count != 1:
            raise RuntimeError("Expected one-page map and original Figure 1 PDFs")
        map_rect = map_document[0].rect
        original_rect = original_document[0].rect
        output_width = original_rect.width
        baseline_output_height = output_width * DOCX_ASPECT_HEIGHT / DOCX_ASPECT_WIDTH
        lower_height = original_rect.y1 - ORIGINAL_CROP_TOP_PT
        scaled_map_height = map_rect.height * output_width / map_rect.width
        canvas_gap_height = baseline_output_height - scaled_map_height - lower_height
        if not 0 <= canvas_gap_height <= 12.0:
            raise RuntimeError(
                f"Unexpected canvas gap {canvas_gap_height:.3f} pt; expected 0–12 pt"
            )
        lower_y0 = scaled_map_height + canvas_gap_height - LOWER_TRANSLATE_UP_PT
        output_height = lower_y0 + lower_height
        visible_gap_target = VISIBLE_GAP_BEFORE_PT * (1.0 - GAP_REDUCTION_FRACTION)

        output_document = fitz.open()
        page = output_document.new_page(width=output_width, height=output_height)
        map_target = fitz.Rect(0, 0, output_width, scaled_map_height)
        page.show_pdf_page(map_target, map_document, 0, keep_proportion=True)

        original_clip = fitz.Rect(
            original_rect.x0,
            ORIGINAL_CROP_TOP_PT,
            original_rect.x1,
            original_rect.y1,
        )
        lower_target = fitz.Rect(
            0,
            lower_y0,
            output_width,
            output_height,
        )
        if abs(lower_target.width - original_clip.width) > 1e-6 or abs(
            lower_target.height - original_clip.height
        ) > 1e-6:
            raise RuntimeError(
                "Panels B–I target dimensions differ from the original clip; refusing distortion"
            )
        page.show_pdf_page(
            lower_target,
            original_document,
            0,
            clip=original_clip,
            keep_proportion=False,
        )
        output_document.save(OUT_PDF, garbage=4, deflate=True)
        OUT_SVG.write_text(page.get_svg_image(text_as_path=False), encoding="utf-8")
        pixmap = page.get_pixmap(matrix=fitz.Matrix(DPI / 72, DPI / 72), alpha=False)
        pixmap.set_dpi(DPI, DPI)
        pixmap.save(OUT_PNG)
        output_document.close()

    payload = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "status": "MAP_WORDING_LAYOUT_AND_COLOR_UPDATE_PLUS_UNALTERED_ORIGINAL_SPECIMEN_PANELS",
        "command": sys.argv,
        "generator": {
            "path": str(SCRIPT.relative_to(ROOT)),
            "sha256": sha256(SCRIPT),
        },
        "inputs": {
            str(path.relative_to(ROOT)): sha256(path) for path in inputs
        },
        "operations": {
            "map_changes": "wording, placement, typography, and display colors only; no coordinate or count changes",
            "dpi_for_docx_raster": DPI,
            "vector_submission_outputs": [str(OUT_PDF.name), str(OUT_SVG.name)],
            "baseline_target_aspect_pixels": [DOCX_ASPECT_WIDTH, DOCX_ASPECT_HEIGHT],
            "original_crop_top_points": ORIGINAL_CROP_TOP_PT,
            "canvas_gap_points": canvas_gap_height,
            "visible_gap_before_points": VISIBLE_GAP_BEFORE_PT,
            "visible_gap_target_points": visible_gap_target,
            "visible_gap_reduction_fraction": GAP_REDUCTION_FRACTION,
            "panels_B_to_I_translation_up_points": LOWER_TRANSLATE_UP_PT,
            "panels_B_to_I_target_points": [
                lower_target.x0,
                lower_target.y0,
                lower_target.x1,
                lower_target.y1,
            ],
            "output_page_points": [output_width, output_height],
            "panels_B_to_I_dimensions_unchanged": True,
            "preserved_photo_panel_labels": "B-I",
            "inset_panel_letter_redaction": False,
        },
        "outputs": {
            str(path.relative_to(ROOT)): {
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in (OUT_PDF, OUT_SVG, OUT_PNG)
        },
    }
    MANIFEST.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(OUT_PDF)
    print(OUT_PNG)
    print(MANIFEST)


if __name__ == "__main__":
    main()
