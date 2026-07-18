#!/usr/bin/env python3
"""Render the global distribution map for the 126-genome analysis set.

Only placement and typography are changed. Coordinates, site aggregation,
marker areas, phylum colors, and the original non-jittered data remain those
of the audited map renderer.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import platform
import re
import sys
from datetime import datetime
from pathlib import Path

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shapereader
import fitz
import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import pyproj
import scipy
import shapely


SCRIPT = Path(__file__).resolve()
OUTDIR = SCRIPT.parent
SOURCE = OUTDIR / "spatial_revision_analysis_20260711_080029.py"
DEFAULT_RUN_TAG = "20260712_180210"
FINAL_WIDTH_PT = 595.276001
FIGURE_WIDTH_IN = FINAL_WIDTH_PT / 72.0
FIGURE_HEIGHT_IN = 3.18 * FIGURE_WIDTH_IN / 7.0
DISPLAY_PHYLUM_COLORS = {
    "Rhodophyta": "#004488",
    "Ochrophyta": "#DDAA33",
    "Chlorophyta": "#228833",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-tag", default=DEFAULT_RUN_TAG)
    return parser.parse_args()


def relative(path: Path, root: Path) -> str:
    return str(path.resolve().relative_to(root.resolve()))


def plot_site_pies(
    ax,
    site_counts: pd.DataFrame,
    phylum_colors: dict[str, str],
    wedge_marker,
    base_size: float,
    label_min_n: int,
    label_offset_points: tuple[float, float] = (6.0, 4.0),
) -> None:
    """Plot the audited site pies with labels clear of marker boundaries."""
    order = list(phylum_colors)
    for site in site_counts.itertuples(index=False):
        total = int(site.total)
        start = 0.0
        marker_size = base_size * (0.75 + 0.45 * math.sqrt(total))
        for phylum in order:
            count = int(getattr(site, phylum))
            if count == 0:
                continue
            end = start + 2.0 * np.pi * count / total
            ax.scatter(
                [site.longitude],
                [site.latitude],
                s=marker_size,
                marker=wedge_marker(start, end),
                facecolor=phylum_colors[phylum],
                edgecolor="none",
                transform=ccrs.PlateCarree(),
                zorder=5,
            )
            start = end
        ax.scatter(
            [site.longitude],
            [site.latitude],
            s=marker_size,
            facecolor="none",
            edgecolor="#202020",
            linewidth=0.25,
            transform=ccrs.PlateCarree(),
            zorder=6,
        )
        if total >= label_min_n:
            ax.annotate(
                str(total),
                xy=(site.longitude, site.latitude),
                xytext=label_offset_points,
                xycoords=ccrs.PlateCarree()._as_mpl_transform(ax),
                textcoords="offset points",
                ha="left",
                va="bottom",
                fontsize=6,
                color="#111111",
                bbox={
                    "boxstyle": "round,pad=0.08",
                    "fc": "white",
                    "ec": "none",
                    "alpha": 0.85,
                },
                zorder=8,
            )


def main() -> None:
    args = parse_args()
    spec = importlib.util.spec_from_file_location("audited_spatial_source", SOURCE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load audited map source: {SOURCE}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    output_pdf = OUTDIR / f"Figure_global_sampling_map_layout_v4_{args.run_tag}.pdf"
    output_svg = OUTDIR / f"Figure_global_sampling_map_layout_v4_{args.run_tag}.svg"
    manifest = OUTDIR / f"sampling_map_layout_v4_manifest_{args.run_tag}.json"
    if any(path.exists() for path in (output_pdf, output_svg, manifest)):
        raise FileExistsError("Refusing to overwrite an existing map-layout output")

    map_data = module.load_map_cohort()
    counts = (
        map_data.groupby(["DD latitude", "DD longitude", "Phylum"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=list(DISPLAY_PHYLUM_COLORS), fill_value=0)
        .reset_index()
        .rename(columns={"DD latitude": "latitude", "DD longitude": "longitude"})
    )
    counts["total"] = counts[list(DISPLAY_PHYLUM_COLORS)].sum(axis=1)
    uae = counts[
        counts["latitude"].between(module.UAE_BBOX["lat_min"], module.UAE_BBOX["lat_max"])
        & counts["longitude"].between(module.UAE_BBOX["lon_min"], module.UAE_BBOX["lon_max"])
    ].copy()

    mpl.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 6,
            "axes.labelsize": 6,
            "axes.titlesize": 6,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "legend.fontsize": 6,
            "axes.linewidth": 0.25,
            "xtick.major.width": 0.25,
            "ytick.major.width": 0.25,
            "xtick.major.size": 2,
            "ytick.major.size": 2,
            "axes.labelpad": 1,
            "xtick.major.pad": 1,
            "ytick.major.pad": 1,
        }
    )

    fig = plt.figure(figsize=(FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN))
    projection = ccrs.Robinson(central_longitude=0)
    world_rect = [0.025, 0.08, 0.66, 0.87]
    inset_rect = [0.74, 0.69, 0.24, 0.24]
    legend_rect = [0.74, 0.39, 0.24, 0.22]
    note_rect = [0.74, 0.19, 0.24, 0.13]

    ax = fig.add_axes(world_rect, projection=projection)
    ax.set_global()
    ax.add_feature(
        cfeature.LAND.with_scale("110m"),
        facecolor="#E7E3DA",
        edgecolor="none",
        zorder=0,
    )
    ax.coastlines(resolution="110m", color="#444444", linewidth=0.25, zorder=2)
    ax.add_feature(
        cfeature.BORDERS.with_scale("110m"),
        edgecolor="#8C8C8C",
        linewidth=0.15,
        zorder=1,
    )
    ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=False,
        xlocs=np.arange(-180, 181, 30),
        ylocs=np.arange(-60, 91, 30),
        linewidth=0.2,
        color="#777777",
        alpha=0.5,
        linestyle=":",
        zorder=1,
    )
    plot_site_pies(
        ax,
        counts,
        DISPLAY_PHYLUM_COLORS,
        module.wedge_marker,
        base_size=8.0,
        label_min_n=10_000,
    )
    ax.set_title(
        "A  Global distribution of 126 macroalgal genomes",
        loc="left",
        fontweight="bold",
        pad=2,
    )

    inset = fig.add_axes(inset_rect, projection=ccrs.PlateCarree())
    inset.set_extent([52.85, 54.55, 24.15, 24.62], crs=ccrs.PlateCarree())
    inset.add_feature(
        cfeature.LAND.with_scale("10m"),
        facecolor="#E7E3DA",
        edgecolor="none",
        zorder=0,
    )
    inset.coastlines(resolution="10m", color="#444444", linewidth=0.25, zorder=2)
    grid = inset.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        xlocs=np.arange(53.0, 54.6, 0.5),
        ylocs=np.arange(24.2, 24.7, 0.2),
        linewidth=0.2,
        color="#777777",
        alpha=0.5,
        linestyle=":",
    )
    grid.top_labels = False
    grid.right_labels = False
    grid.xlabel_style = {"size": 6}
    grid.ylabel_style = {"size": 6}
    plot_site_pies(
        inset,
        uae,
        DISPLAY_PHYLUM_COLORS,
        module.wedge_marker,
        base_size=55.0,
        label_min_n=2,
        label_offset_points=(7.0, 5.0),
    )
    inset.set_title(
        "UAE: 9 genomes at 3 coordinate pairs",
        loc="left",
        fontweight="bold",
        pad=2,
    )

    scale_lat = 24.18
    scale_lon_start = 53.0
    scale_km = 50.0
    scale_deg_lon = scale_km / (111.32 * np.cos(np.radians(scale_lat)))
    inset.plot(
        [scale_lon_start, scale_lon_start + scale_deg_lon],
        [scale_lat, scale_lat],
        transform=ccrs.PlateCarree(),
        color="#111111",
        linewidth=0.8,
        solid_capstyle="butt",
        zorder=10,
    )
    for x_value in (scale_lon_start, scale_lon_start + scale_deg_lon):
        inset.plot(
            [x_value, x_value],
            [scale_lat - 0.008, scale_lat + 0.008],
            transform=ccrs.PlateCarree(),
            color="#111111",
            linewidth=0.5,
            zorder=10,
        )
    inset.text(
        scale_lon_start + scale_deg_lon / 2,
        scale_lat + 0.016,
        "50 km",
        ha="center",
        va="bottom",
        transform=ccrs.PlateCarree(),
        fontsize=6,
    )
    # The arrow is placed in the central-water gap, away from all three UAE sites.
    inset.text(
        53.72,
        24.55,
        "N",
        transform=ccrs.PlateCarree(),
        ha="center",
        va="center",
        fontsize=6,
        fontweight="bold",
        zorder=10,
    )
    inset.annotate(
        "",
        xy=(53.72, 24.50),
        xytext=(53.72, 24.42),
        xycoords=ccrs.PlateCarree()._as_mpl_transform(inset),
        textcoords=ccrs.PlateCarree()._as_mpl_transform(inset),
        arrowprops={"arrowstyle": "-|>", "lw": 0.4, "color": "#111111"},
        zorder=10,
    )

    legend_ax = fig.add_axes(legend_rect)
    legend_ax.axis("off")
    phylum_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor=color,
            markeredgecolor="#222222",
            markeredgewidth=0.25,
            markersize=4.5,
            label=f"{phylum} (n = {int((map_data['Phylum'] == phylum).sum())})",
        )
        for phylum, color in DISPLAY_PHYLUM_COLORS.items()
    ]
    phylum_legend = legend_ax.legend(
        handles=phylum_handles,
        loc="upper left",
        frameon=False,
        title="Phylum",
        title_fontsize=6,
        handletextpad=0.4,
        borderaxespad=0,
        labelspacing=0.35,
        alignment="left",
    )
    legend_ax.add_artist(phylum_legend)
    size_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor="#BDBDBD",
            markeredgecolor="#222222",
            markeredgewidth=0.25,
            markersize=size,
            label=f"{n} genome{'s' if n > 1 else ''}",
        )
        for n, size in [(1, 3.0), (4, 4.5), (7, 5.5)]
    ]
    legend_ax.legend(
        handles=size_handles,
        loc="upper right",
        frameon=False,
        title="Genomes per coordinate",
        title_fontsize=6,
        handletextpad=0.4,
        borderaxespad=0,
        labelspacing=0.35,
        alignment="left",
    )

    note_ax = fig.add_axes(note_rect)
    note_ax.axis("off")
    note_ax.text(
        0.0,
        1.0,
        "Pie segments: phylum composition\n"
        "UAE numbers: shared-coordinate counts\n"
        "Markers use recorded coordinates; no jitter.",
        transform=note_ax.transAxes,
        ha="left",
        va="top",
        fontsize=6,
        linespacing=1.30,
    )

    raw_pdf = output_pdf.with_name(f".{output_pdf.stem}.matplotlib-render.pdf")
    if raw_pdf.exists():
        raise FileExistsError(f"Refusing to overwrite temporary render: {raw_pdf}")
    fig.savefig(
        raw_pdf,
        format="pdf",
        transparent=True,
        edgecolor="none",
    )
    fig.savefig(
        output_svg,
        format="svg",
        transparent=True,
        edgecolor="none",
    )
    plt.close(fig)

    # Matplotlib 3.9 quantizes the physical canvas to 0.01 inch. Embed that
    # unscaled page on the exact A4-width specimen canvas so the 6-pt text is
    # not enlarged by the downstream composer.
    with fitz.open(raw_pdf) as raw_document:
        raw_rect = raw_document[0].rect
        normalized_document = fitz.open()
        normalized_page = normalized_document.new_page(
            width=FINAL_WIDTH_PT,
            height=raw_rect.height,
        )
        normalized_page.show_pdf_page(
            fitz.Rect(0, 0, raw_rect.width, raw_rect.height),
            raw_document,
            0,
            keep_proportion=False,
        )
        normalized_document.save(output_pdf, garbage=4, deflate=True)
        normalized_document.close()
    raw_pdf.unlink()

    svg_text = output_svg.read_text(encoding="utf-8")
    svg_text, width_hits = re.subn(
        r'width="594\.72pt" height="270pt" viewBox="0 0 594\.72 270"',
        f'width="{FINAL_WIDTH_PT:.6f}pt" height="270pt" viewBox="0 0 {FINAL_WIDTH_PT:.6f} 270"',
        svg_text,
        count=1,
    )
    if width_hits != 1:
        raise RuntimeError("Could not normalize the SVG canvas to exact A4 width")
    output_svg.write_text(svg_text, encoding="utf-8")

    figure_stats = {
        "samples": int(len(map_data)),
        "unique_coordinate_pairs": int(len(counts)),
        "shared_coordinate_pairs": int((counts["total"] > 1).sum()),
        "samples_at_shared_coordinates": int(
            counts.loc[counts["total"] > 1, "total"].sum()
        ),
        "max_genomes_at_one_coordinate": int(counts["total"].max()),
        "uae_samples": int(
            map_data["DD latitude"]
            .between(module.UAE_BBOX["lat_min"], module.UAE_BBOX["lat_max"])
            .mul(
                map_data["DD longitude"].between(
                    module.UAE_BBOX["lon_min"], module.UAE_BBOX["lon_max"]
                )
            )
            .sum()
        ),
        "uae_coordinate_pairs": int(len(uae)),
    }
    if figure_stats != {
        "samples": 126,
        "unique_coordinate_pairs": 90,
        "shared_coordinate_pairs": 19,
        "samples_at_shared_coordinates": 55,
        "max_genomes_at_one_coordinate": 7,
        "uae_samples": 9,
        "uae_coordinate_pairs": 3,
    }:
        raise RuntimeError(f"Map statistics changed unexpectedly: {figure_stats}")

    natural_earth_specs = [
        ("land", "physical", "110m"),
        ("coastline", "physical", "110m"),
        ("admin_0_boundary_lines_land", "cultural", "110m"),
        ("land", "physical", "10m"),
        ("coastline", "physical", "10m"),
    ]
    natural_earth_assets = []
    for name, category, resolution in natural_earth_specs:
        shapefile = Path(
            shapereader.natural_earth(
                resolution=resolution,
                category=category,
                name=name,
            )
        )
        natural_earth_assets.append(
            {
                "dataset": name,
                "category": category,
                "resolution": resolution,
                "components": [
                    {"filename": component.name, "sha256": sha256(component)}
                    for component in sorted(shapefile.parent.glob(f"{shapefile.stem}.*"))
                    if component.is_file()
                ],
            }
        )

    payload = {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "status": "MAP_LAYOUT_ONLY_NO_DATA_OR_COORDINATE_CHANGE",
        "command": sys.argv,
        "generator": {"path": relative(SCRIPT, module.ROOT), "sha256": sha256(SCRIPT)},
        "imported_audited_source": {
            "path": relative(SOURCE, module.ROOT),
            "sha256": sha256(SOURCE),
        },
        "inputs": [
            {"path": relative(path, module.ROOT), "sha256": sha256(path)}
            for path in (module.AEF_EMBEDDINGS, module.CLEAN_METADATA)
        ],
        "outputs": [
            {"path": relative(path, module.ROOT), "sha256": sha256(path)}
            for path in (output_pdf, output_svg)
        ],
        "layout_parameters": {
            "figure_inches": [FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN],
            "final_width_points": FINAL_WIDTH_PT,
            "world_axes": world_rect,
            "uae_inset_axes": inset_rect,
            "legend_axes": legend_rect,
            "note_axes": note_rect,
            "north_arrow_longitude": 53.72,
            "north_arrow_latitude_range": [24.42, 24.55],
            "uae_label_offset_points": [7.0, 5.0],
            "fixed_canvas_no_tight_bbox": True,
            "matplotlib_raw_canvas_points": [594.72, 270.0],
            "normalized_canvas_points": [FINAL_WIDTH_PT, 270.0],
            "normalized_without_content_scaling": True,
            "display_phylum_colors": DISPLAY_PHYLUM_COLORS,
            "color_change_scope": "presentation only; genome-set and phylum assignments unchanged",
            "coordinates_jittered": False,
        },
        "figure_statistics": figure_stats,
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "cartopy": cartopy.__version__,
            "matplotlib": matplotlib.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "pyproj": pyproj.__version__,
            "scipy": scipy.__version__,
            "shapely": shapely.__version__,
        },
        "natural_earth_assets": natural_earth_assets,
    }
    manifest.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(output_pdf)
    print(output_svg)
    print(manifest)


if __name__ == "__main__":
    main()
