#!/usr/bin/env python3
"""Create submission-facing Figure S3 from completed analysis outputs.

The script performs no scientific analysis and creates no synthetic values. It
loads the completed 20260711_131706 result tables and invokes the existing plot
function after substituting two display titles in memory. All other SVG text,
geometry, embedded image data, and numerical labels are required to match the
previous vector rendering exactly. A clean one-page caption PDF is then built
from numerical values read by the established Figure S3 package validator.

Outputs are versioned and are never overwritten.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import platform
import re
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from lxml import etree


STAMP = "20260713_100151"
ROOT = Path(__file__).resolve().parents[2]
SUPP = ROOT / "ISCIENCE_REVISION_20260711" / "supplement"
STATS = ROOT / "ISCIENCE_REVISION_20260711" / "analysis_stats"

SOURCE_SCRIPT = STATS / "run_robustness_20260711_085930.py"
SOURCE_RESULTS = STATS / "archived_AEF_priority_robustness_20260711_131706.csv"
SOURCE_NULL = STATS / "archived_AEF_within_phylum_null_20260711_131706.csv"
SOURCE_BUILDER = SUPP / "build_figure_s3_supporting_package_20260711_175934.py"
REFERENCE_PDF = SUPP / "Figure_S3_revision_stage_sensitivity_20260712_074152.pdf"
REFERENCE_SVG = SUPP / "Figure_S3_revision_stage_sensitivity_20260712_074152.svg"

FIGURE_STEM = f"Figure_S3_posthoc_sensitivity_checks_V2_{STAMP}"
PAGE_STEM = f"Figure_S3_posthoc_sensitivity_checks_supporting_page_V2_{STAMP}"
FIG_PDF = SUPP / f"{FIGURE_STEM}.pdf"
FIG_SVG = SUPP / f"{FIGURE_STEM}.svg"
CAPTION_TXT = SUPP / f"Figure_S3_posthoc_sensitivity_checks_caption_V2_{STAMP}.txt"
PAGE_TEX = SUPP / f"{PAGE_STEM}.tex"
PAGE_PDF = SUPP / f"{PAGE_STEM}.pdf"
MANIFEST = SUPP / f"Figure_S3_posthoc_sensitivity_checks_V2_integrity_{STAMP}.json"
SCRIPT = Path(__file__).resolve()

OLD_PANEL_A = "A  Priority archived latent-feature associations across robustness checks"
REFERENCE_PANEL_A = "A  Revision-stage latent-feature associations across robustness checks"
NEW_PANEL_A = "A  Selected Pfam–AEF associations across sensitivity checks"
OLD_PANEL_C = "C  Retention by check"
NEW_PANEL_C = "C  Consistency across checks"

EXPECTED_INPUT_SHA256 = {
    SOURCE_SCRIPT: "ac5c3dc9676dfc69bddca373679445b8ecbde095e88489f17e34b34e0ba226ca",
    SOURCE_RESULTS: "99b54015007ef2c642772ed2bbe57f1bfa0572f36bc66fed1e1ce3f15441b39e",
    SOURCE_NULL: "7b49c7e2535cfb4a69f24a31f45b785d0467441b77e24ea92470db27ffa5d972",
    SOURCE_BUILDER: "15582c1d15876b4549e52e00bf808040939bfffe77e899f70f2cfc0e61b27779",
    REFERENCE_PDF: "0aeadc82e970909225af08d36521f4cb74e11838fd93f74ac7297fd5a30edbe6",
    REFERENCE_SVG: "307668f7dcc74c509d16533225ad31ce30a4f259efb2afe0320e532c53311501",
}

SVG_NS = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"
SVG_NAMESPACES = {"svg": SVG_NS, "xlink": XLINK_NS}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def numerical_tokens(text: str) -> Counter[str]:
    return Counter(re.findall(r"\d+(?:[.,]\d+)*(?:e[+-]?\d+)?", text, flags=re.IGNORECASE))


def load_module(path: Path, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def svg_root(path: Path) -> etree._Element:
    parser = etree.XMLParser(resolve_entities=False, remove_blank_text=False)
    return etree.parse(str(path), parser=parser).getroot()


def svg_texts(root: etree._Element) -> list[str]:
    return ["".join(element.itertext()) for element in root.xpath(".//svg:text", namespaces=SVG_NAMESPACES)]


def svg_geometry_signatures(root: etree._Element) -> list[tuple[str, tuple[tuple[str, str], ...]]]:
    """Capture ordered drawable geometry while ignoring randomized SVG IDs."""
    drawable = {"path", "rect", "use", "image", "line", "circle", "polygon", "polyline"}
    signatures: list[tuple[str, tuple[tuple[str, str], ...]]] = []
    for element in root.iter():
        local_name = etree.QName(element).localname
        if local_name not in drawable:
            continue
        attributes: list[tuple[str, str]] = []
        for key, value in element.attrib.items():
            attribute_name = etree.QName(key).localname if key.startswith("{") else key
            if attribute_name in {"id", "href", "clip-path"}:
                continue
            attributes.append((attribute_name, value))
        signatures.append((local_name, tuple(sorted(attributes))))
    return signatures


def svg_embedded_images(root: etree._Element) -> list[str]:
    return [
        element.get(f"{{{XLINK_NS}}}href", "")
        for element in root.xpath(".//svg:image", namespaces=SVG_NAMESPACES)
    ]


def validate_source_files() -> None:
    for path, expected_hash in EXPECTED_INPUT_SHA256.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        observed_hash = sha256(path)
        if observed_hash != expected_hash:
            raise RuntimeError(
                f"Input hash mismatch for {path}: expected {expected_hash}, observed {observed_hash}"
            )
    for output in (FIG_PDF, FIG_SVG, CAPTION_TXT, PAGE_TEX, PAGE_PDF, MANIFEST):
        if output.exists():
            raise FileExistsError(f"Refusing to overwrite versioned output: {output}")


def render_plot() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    source_text = SOURCE_SCRIPT.read_text(encoding="utf-8")
    substitutions = {
        OLD_PANEL_A: NEW_PANEL_A,
        OLD_PANEL_C: NEW_PANEL_C,
    }
    modified_source = source_text
    for old_text, new_text in substitutions.items():
        if modified_source.count(old_text) != 1:
            raise RuntimeError(f"Expected one plotting-label occurrence: {old_text}")
        modified_source = modified_source.replace(old_text, new_text)

    namespace = {"__name__": "figure_s3_plot_source", "__file__": str(SOURCE_SCRIPT)}
    exec(compile(modified_source, str(SOURCE_SCRIPT), "exec"), namespace)
    plotter = namespace.get("create_robustness_figure")
    if plotter is None:
        raise RuntimeError("Existing Figure S3 plot function was not defined")

    results = pd.read_csv(SOURCE_RESULTS)
    null_results = pd.read_csv(SOURCE_NULL)
    if len(results) != 2856 or len(null_results) != 204:
        raise RuntimeError(
            f"Unexpected completed-output row counts: {len(results)} association rows, "
            f"{len(null_results)} structured-null rows"
        )
    plotter(results, null_results, FIG_PDF, FIG_SVG)
    return results, null_results, {
        "from": substitutions,
        "count": len(substitutions),
    }


def clean_caption(source_caption: str) -> tuple[str, list[dict[str, str]]]:
    replacements = [
        (
            "Post hoc selected-set sensitivity analysis of archived Pfam–AEF latent-feature associations",
            "Post hoc selected-set sensitivity analysis of Pfam–AEF latent-feature associations",
        ),
        ("the four previously highlighted discovery pairs", "four priority discovery pairs"),
        ("four highlighted pairs plus PF00092", "four priority discovery pairs plus PF00092"),
        ("These 12 are a display subset of 68 pairs tested", "These 12 are a display subset of 68 pairs evaluated"),
        (
            "Only PF01638/A52 and PF01638/A53 met direction consistency and all seven required selected-set gates",
            "PF01638–A52 and PF01638–A53 were directionally consistent and met the selected-set criterion across all seven required checks",
        ),
        (
            "PF10988/A36 (6/7) and PF13411/A18 (5/7) did not pass all gates.",
            "PF10988–A36 met six of seven checks (6/7), and PF13411–A18 met five of seven checks (5/7).",
        ),
    ]
    caption = source_caption
    applied: list[dict[str, str]] = []
    for old_text, new_text in replacements:
        if caption.count(old_text) != 1:
            raise RuntimeError(f"Expected one caption-language occurrence: {old_text}")
        caption = caption.replace(old_text, new_text)
        applied.append({"from": old_text, "to": new_text})

    if numerical_tokens(source_caption) != numerical_tokens(caption):
        raise RuntimeError(
            "Numerical caption tokens changed:\n"
            f"SOURCE {numerical_tokens(source_caption)}\nCLEAN {numerical_tokens(caption)}"
        )

    banned = re.compile(
        r"\b(?:revision|revision-stage|previously|gate|gates|retention|archived|retained|provenance)\b",
        flags=re.IGNORECASE,
    )
    banned_matches = sorted(set(match.group(0) for match in banned.finditer(caption)))
    if banned_matches:
        raise RuntimeError(f"Submission-facing caption language remains: {banned_matches}")
    return caption, applied


def run_checked(command: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, cwd=cwd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed ({result.returncode}): {' '.join(command)}\n"
            f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )
    return result


def pdf_pages(path: Path) -> int:
    info = run_checked(["pdfinfo", str(path)]).stdout
    match = re.search(r"^Pages:\s+(\d+)$", info, flags=re.MULTILINE)
    if not match:
        raise RuntimeError(f"Could not read PDF page count: {path}")
    return int(match.group(1))


def pdf_page_size(path: Path) -> str:
    info = run_checked(["pdfinfo", str(path)]).stdout
    match = re.search(r"^Page size:\s+(.+)$", info, flags=re.MULTILINE)
    if not match:
        raise RuntimeError(f"Could not read PDF page size: {path}")
    return match.group(1).strip()


def audit_svg_equivalence() -> dict[str, object]:
    reference_root = svg_root(REFERENCE_SVG)
    output_root = svg_root(FIG_SVG)

    reference_text = svg_texts(reference_root)
    output_text = svg_texts(output_root)
    if len(reference_text) != len(output_text):
        raise RuntimeError("SVG text element count changed")
    expected_text = [
        NEW_PANEL_A if text == REFERENCE_PANEL_A else NEW_PANEL_C if text == OLD_PANEL_C else text
        for text in reference_text
    ]
    if output_text != expected_text:
        differences = [
            {"index": index, "expected": expected, "observed": observed}
            for index, (expected, observed) in enumerate(zip(expected_text, output_text))
            if expected != observed
        ]
        raise RuntimeError(f"Unexpected SVG text changes: {differences}")

    reference_geometry = svg_geometry_signatures(reference_root)
    output_geometry = svg_geometry_signatures(output_root)
    if reference_geometry != output_geometry:
        raise RuntimeError("SVG drawable geometry changed")

    reference_images = svg_embedded_images(reference_root)
    output_images = svg_embedded_images(output_root)
    if reference_images != output_images:
        raise RuntimeError("Embedded image payloads changed")

    canvas_attributes = ("width", "height", "viewBox")
    reference_canvas = {key: reference_root.get(key) for key in canvas_attributes}
    output_canvas = {key: output_root.get(key) for key in canvas_attributes}
    if reference_canvas != output_canvas:
        raise RuntimeError(f"SVG canvas changed: {reference_canvas} != {output_canvas}")

    font_style_failures = []
    for index, element in enumerate(output_root.xpath(".//svg:text", namespaces=SVG_NAMESPACES)):
        style = element.get("style", "")
        if "6px 'Arial', 'Helvetica', sans-serif" not in style:
            font_style_failures.append({"index": index, "style": style, "text": "".join(element.itertext())})
    if font_style_failures:
        raise RuntimeError(f"Non-compliant SVG text styles: {font_style_failures}")

    svg_string = FIG_SVG.read_text(encoding="utf-8")
    stroke_widths = re.findall(r"stroke-width:\s*([0-9.]+)", svg_string)
    if not stroke_widths or set(stroke_widths) != {"0.25"}:
        raise RuntimeError(f"Unexpected SVG stroke widths: {sorted(set(stroke_widths))}")

    return {
        "reference_text_elements": len(reference_text),
        "output_text_elements": len(output_text),
        "intended_text_changes": [
            {"from": REFERENCE_PANEL_A, "to": NEW_PANEL_A},
            {"from": OLD_PANEL_C, "to": NEW_PANEL_C},
        ],
        "all_other_text_elements_identical": True,
        "drawable_geometry_elements": len(output_geometry),
        "drawable_geometry_identical": True,
        "embedded_image_count": len(output_images),
        "embedded_image_payloads_identical": True,
        "canvas_identical": True,
        "canvas": output_canvas,
        "all_text_6px_arial": True,
        "stroke_width_values": sorted(set(stroke_widths)),
    }


def audit_pdf_fonts(path: Path) -> dict[str, object]:
    output = run_checked(["pdffonts", str(path)]).stdout
    font_rows = [line for line in output.splitlines()[2:] if line.strip()]
    if not font_rows:
        raise RuntimeError(f"No fonts reported for {path}")
    if any("Arial" not in row or " yes " not in f" {row} " for row in font_rows):
        raise RuntimeError(f"Figure PDF font audit failed:\n{output}")
    return {"font_rows": font_rows, "arial_only": True, "fonts_embedded": True}


def main() -> None:
    validate_source_files()

    package = load_module(SOURCE_BUILDER, "figure_s3_package_source")
    validated = package.validate_inputs(ROOT)
    source_caption = package.build_caption(validated)
    caption, caption_replacements = clean_caption(source_caption)

    results, null_results, plot_substitutions = render_plot()
    CAPTION_TXT.write_text(caption + "\n", encoding="utf-8")
    PAGE_TEX.write_text(package.build_tex(caption, FIG_PDF.name), encoding="utf-8")
    run_checked(
        [
            "tectonic",
            "--only-cached",
            "--chatter",
            "minimal",
            "--outdir",
            str(SUPP),
            str(PAGE_TEX),
        ],
        cwd=SUPP,
    )

    for path in (FIG_PDF, PAGE_PDF):
        run_checked(["qpdf", "--check", str(path)])
        if pdf_pages(path) != 1:
            raise RuntimeError(f"Expected a one-page PDF: {path}")

    svg_audit = audit_svg_equivalence()
    font_audit = audit_pdf_fonts(FIG_PDF)

    supporting_text = run_checked(["pdftotext", "-raw", str(PAGE_PDF), "-"]).stdout
    supporting_normalized = " ".join(supporting_text.split())
    for required_text in (
        "Selected Pfam–AEF associations across sensitivity checks",
        "Consistency across checks",
        "four priority discovery pairs",
        "all seven required checks",
    ):
        if required_text not in supporting_normalized:
            raise RuntimeError(f"Supporting-page PDF lacks required text: {required_text}")
    banned_supporting = re.compile(
        r"\b(?:revision|revision-stage|previously|gate|gates|retention|archived|retained|provenance)\b",
        flags=re.IGNORECASE,
    )
    supporting_banned_matches = sorted(
        set(match.group(0) for match in banned_supporting.finditer(supporting_text))
    )
    if supporting_banned_matches:
        raise RuntimeError(f"Supporting-page process language remains: {supporting_banned_matches}")

    outputs = (FIG_PDF, FIG_SVG, CAPTION_TXT, PAGE_TEX, PAGE_PDF)
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": str(SCRIPT),
        "script_sha256": sha256(SCRIPT),
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "matplotlib": matplotlib.__version__,
        "inputs": {str(path): sha256(path) for path in EXPECTED_INPUT_SHA256},
        "completed_analysis_rows": {
            "association_rows_loaded": len(results),
            "structured_null_rows_loaded": len(null_results),
            "selected_pairs_validated": len(validated["candidates"]),
            "displayed_pairs_validated": len(validated["displayed"]),
        },
        "outputs": {str(path): sha256(path) for path in outputs},
        "operation": {
            "scientific_analysis_rerun": False,
            "scientific_values_recomputed": False,
            "synthetic_or_placeholder_data_used": False,
            "plotting_from_completed_result_rows": True,
            "plot_title_substitutions": plot_substitutions,
            "caption_replacements": caption_replacements,
            "caption_numerical_tokens_preserved": True,
        },
        "validations": {
            "source_hashes_match_pinned_inputs": True,
            "candidate_set_equals_68_pairs": True,
            "display_subset_equals_12_pairs": True,
            "all_check_consistent_pairs_equal_PF01638_A52_A53": True,
            "svg": svg_audit,
            "figure_pdf_pages": pdf_pages(FIG_PDF),
            "figure_pdf_page_size": pdf_page_size(FIG_PDF),
            "supporting_page_pdf_pages": pdf_pages(PAGE_PDF),
            "supporting_page_pdf_page_size": pdf_page_size(PAGE_PDF),
            "figure_pdf_fonts": font_audit,
            "qpdf_checks_passed": True,
            "supporting_page_required_wording_present": True,
            "supporting_page_banned_terms": supporting_banned_matches,
        },
        "interpretation_boundary": (
            "Post hoc selected-set sensitivity analysis of priority Pfam–AEF associations; "
            "AEF axes are unitless latent descriptors and do not identify physical variables, "
            "adaptation, or mechanism."
        ),
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print(json.dumps({
        "figure_pdf": str(FIG_PDF.resolve()),
        "figure_svg": str(FIG_SVG.resolve()),
        "caption": str(CAPTION_TXT.resolve()),
        "supporting_page_tex": str(PAGE_TEX.resolve()),
        "supporting_page_pdf": str(PAGE_PDF.resolve()),
        "manifest": str(MANIFEST.resolve()),
        "figure_pdf_sha256": sha256(FIG_PDF),
        "figure_svg_sha256": sha256(FIG_SVG),
        "supporting_page_pdf_sha256": sha256(PAGE_PDF),
        "geometry_elements_preserved": svg_audit["drawable_geometry_elements"],
        "text_elements_preserved_except_two_titles": svg_audit["output_text_elements"],
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
