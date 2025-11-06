from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple

import geopandas as gpd
import numpy as np
import pandas as pd
import rioxarray as rxr
from rasterstats import zonal_stats
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

import utils


@dataclass(slots=True)
class AspectData:
    bin_edges: np.ndarray
    background_counts: np.ndarray
    event_angles: np.ndarray
    historical_angles: np.ndarray


@dataclass(slots=True)
class FrequencyResults:
    event_d8: pd.DataFrame
    historical_d8: pd.DataFrame
    coverage_table: pd.DataFrame


def _circular_mean_deg(values: Iterable[float]) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    arr = arr[(arr >= 0) & (arr <= 360)]
    if arr.size == 0:
        return float("nan")
    ang = np.deg2rad(arr)
    s = np.sin(ang).mean()
    c = np.cos(ang).mean()
    if np.isclose(s, 0.0) and np.isclose(c, 0.0):
        return float("nan")
    return float((np.degrees(np.arctan2(s, c)) % 360.0))


def _load_aspect_angles(cfg, inventories: utils.Inventories) -> AspectData:
    aspect_da = rxr.open_rasterio(cfg.aspect_file, masked=True).squeeze(drop=True)
    aspect_vals = aspect_da.values
    aspect_vals = aspect_vals[np.isfinite(aspect_vals)]
    aspect_vals = aspect_vals[(aspect_vals >= 0) & (aspect_vals <= 360)]

    bin_edges = np.linspace(0, 360, 37)
    background_counts, _ = np.histogram(aspect_vals, bins=bin_edges)

    custom_stats = {"aspect_circ_mean": _circular_mean_deg}

    event_stats = []
    for geom in tqdm(inventories.event.geometry, desc="Event aspect means"):
        out = zonal_stats(
            geom,
            cfg.aspect_file,
            stats=[],
            add_stats=custom_stats,
            nodata=None,
            geojson_out=False,
        )
        event_stats.append(out[0])

    hist_stats = []
    for geom in tqdm(inventories.historical.geometry, desc="Historical aspect means"):
        out = zonal_stats(
            geom,
            cfg.aspect_file,
            stats=[],
            add_stats=custom_stats,
            nodata=None,
            geojson_out=False,
        )
        hist_stats.append(out[0])

    event_angles = np.asarray(
        [row.get("aspect_circ_mean", np.nan) for row in event_stats], dtype=float
    )
    historical_angles = np.asarray(
        [row.get("aspect_circ_mean", np.nan) for row in hist_stats], dtype=float
    )

    event_angles = event_angles[np.isfinite(event_angles)] % 360.0
    historical_angles = historical_angles[np.isfinite(historical_angles)] % 360.0

    return AspectData(
        bin_edges=bin_edges,
        background_counts=background_counts.astype(float),
        event_angles=event_angles,
        historical_angles=historical_angles,
    )


def _angles_in_range(theta: np.ndarray, start: float, end: float) -> np.ndarray:
    start %= 360.0
    end %= 360.0
    if start <= end:
        return (theta >= start) & (theta < end)
    return (theta >= start) | (theta < end)


def _proportion_from_hist(bin_edges: np.ndarray, counts: np.ndarray, start: float, end: float) -> float:
    total = counts.sum()
    if total <= 0:
        return 0.0
    start %= 360.0
    end %= 360.0
    prop = 0.0
    for (low, high), count in zip(zip(bin_edges[:-1], bin_edges[1:]), counts):
        if count <= 0:
            continue
        width = high - low
        segments: List[Tuple[float, float]]
        if start <= end:
            segments = [(start, end)]
        else:
            segments = [(start, 360.0), (0.0, end)]
        overlap = 0.0
        for s, e in segments:
            lo = max(low, s)
            hi = min(high, e)
            overlap = max(overlap, hi - lo)
        if overlap > 0:
            prop += count * (overlap / width)
    return float(prop / total)


def _build_d8_table(
    angles: np.ndarray,
    label: str,
    background_counts: np.ndarray,
    bin_edges: np.ndarray,
) -> pd.DataFrame:
    d8_bins = [
        ("N", 337.5, 22.5),
        ("NE", 22.5, 67.5),
        ("E", 67.5, 112.5),
        ("SE", 112.5, 157.5),
        ("S", 157.5, 202.5),
        ("SW", 202.5, 247.5),
        ("W", 247.5, 292.5),
        ("NW", 292.5, 337.5),
    ]

    landscape_total = background_counts.sum()
    group_total = angles.size

    rows = []
    for name, start, end in d8_bins:
        p_bg = _proportion_from_hist(bin_edges, background_counts, start, end)
        mask = _angles_in_range(angles, start, end)
        p_group = mask.sum() / group_total if group_total else 0.0
        fr = (p_group / p_bg) if p_bg > 0 else np.nan
        rows.append(
            {
                "Bin": name,
                "Landscape %": p_bg * 100.0,
                f"{label} %": p_group * 100.0,
                "FR": fr,
            }
        )
    return pd.DataFrame(rows)


def _coverage_table(angles_event: np.ndarray, angles_hist: np.ndarray, background_counts: np.ndarray, bin_edges: np.ndarray) -> pd.DataFrame:
    ranges = [
        ("N–E", 337.5, 112.5),
        ("E–S", 112.5, 202.5),
        ("S–W", 202.5, 292.5),
        ("W–N", 292.5, 337.5),
    ]

    rows = []
    for label, start, end in ranges:
        # Landscape percentage from background
        landscape_pct = _proportion_from_hist(bin_edges, background_counts, start, end) * 100.0
        event_pct = _angles_in_range(angles_event, start, end).mean() * 100.0
        hist_pct = _angles_in_range(angles_hist, start, end).mean() * 100.0
        fr = (event_pct / landscape_pct) if landscape_pct > 0 else np.nan
        rows.append(
            {
                "Range": label,
                "Landscape %": landscape_pct,
                "Event %": event_pct,
                "Historical %": hist_pct,
                "FR": fr,
            }
        )
    return pd.DataFrame(rows)


def compute_aspect_frequency_results(cfg) -> FrequencyResults:
    inventories = utils.load_inventories(cfg)
    aspect_data = _load_aspect_angles(cfg, inventories)

    event_table = _build_d8_table(
        aspect_data.event_angles,
        label="Event",
        background_counts=aspect_data.background_counts,
        bin_edges=aspect_data.bin_edges,
    )
    hist_table = _build_d8_table(
        aspect_data.historical_angles,
        label="Historical",
        background_counts=aspect_data.background_counts,
        bin_edges=aspect_data.bin_edges,
    )
    coverage = _coverage_table(
        aspect_data.event_angles, 
        aspect_data.historical_angles,
        aspect_data.background_counts,
        aspect_data.bin_edges
    )

    return FrequencyResults(
        event_d8=event_table,
        historical_d8=hist_table,
        coverage_table=coverage,
    )


def plot_rose_diagram(aspect_data: AspectData, cfg) -> plt.Figure:
    """Create polar rose diagram matching notebook style"""
    landscape_color = "#6B8E23"  # olive green
    landslide_color = "#762A83"   # reddish
    
    # Compute proportions for each bin
    bin_centers = (aspect_data.bin_edges[:-1] + aspect_data.bin_edges[1:]) / 2
    bin_widths = np.diff(aspect_data.bin_edges)
    
    # Background (landscape) proportions
    bg_p = aspect_data.background_counts / aspect_data.background_counts.sum()
    
    # Landslide proportions
    ls_counts, _ = np.histogram(aspect_data.event_angles, bins=aspect_data.bin_edges)
    ls_p = ls_counts / ls_counts.sum()
    
    # Convert to radians
    centers_rad = np.deg2rad(bin_centers)
    bin_widths_rad = np.deg2rad(bin_widths)
    
    # Create polar figure
    fig, ax = plt.subplots(figsize=(60/25.4, 60/25.4), subplot_kw=dict(polar=True))
    
    # Plot landscape bars (full bin width)
    ax.bar(
        centers_rad,
        bg_p,
        width=bin_widths_rad,
        bottom=0.0,
        color=landscape_color,
        alpha=0.35,
        edgecolor="none",
        label="Landscape"
    )
    
    # Plot landslide bars (narrower, on top)
    ax.bar(
        centers_rad,
        ls_p,
        width=bin_widths_rad * 0.7,
        bottom=0.0,
        color=landslide_color,
        alpha=0.75,
        edgecolor="white",
        linewidth=0.5,
        label="Landslide"
    )
    
    # Orientation: North up, clockwise
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    
    # Ticks and labels
    xtick_angles = [0, 45, 90, 135, 180, 225, 270, 315]
    xtick_labels = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]
    
    ax.set_xticks(np.deg2rad(xtick_angles))
    
    for label, angle in zip(xtick_labels, xtick_angles):
        ax.text(
            np.deg2rad(angle), 
            ax.get_rmax() * 1.05, 
            label, 
            ha="center", 
            va="center",
        )
    # Remove default xticklabels to avoid overlap
    ax.set_xticklabels([])
    
    ax.set_yticklabels([])  # hide radial tick labels
    
    # Add radial tick labels as percentages in the SW quadrant
    ticks = np.arange(0.01, 0.051, 0.02)
    ax.set_yticks(ticks)
    
    for label, tick in zip([f"{int(t*100)}%" for t in ticks], ticks):
        ax.text(
            np.deg2rad(247),  # angle in radians (SW)
            tick,
            label,
            ha="center",
            va="center",
            color="k",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.3, boxstyle="round,pad=0.1")
        )
    
    ax.set_yticklabels([])
    ax.grid(alpha=0.8)
    
    plt.tight_layout()
    
    return fig


def plot_coverage_bars(results: FrequencyResults, cfg) -> plt.Figure:
    """Create N-E coverage bar plot matching notebook style"""
    import matplotlib as mpl
    import seaborn as sns
    
    mpl.rcParams['font.size'] = 10
    mpl.rcParams['font.family'] = 'Arial'
    font_size = mpl.rcParams['font.size']
    font_family = mpl.rcParams['font.family']
    
    landscape_color = "#6B8E23"
    landslide_color = "#762A83"
    
    fig, ax = plt.subplots(figsize=(40/25.4, 50/25.4), sharex=False, constrained_layout=True)
    
    # Prepare data for grouped bar plot - only N-E range
    cust_long = results.coverage_table.copy()
    # Filter to only N-E range
    cust_long = cust_long[cust_long["Range"] == "N–E"]
    # Rename columns to match notebook
    cust_long = cust_long.rename(columns={"Event %": "Landslide %"})
    cust_long = cust_long.melt(
        id_vars=["Range", "FR"], 
        value_vars=["Landscape %", "Landslide %"],
        var_name="Group", 
        value_name="Percent"
    )
    
    palette = {"Landscape %": landscape_color, "Landslide %": landslide_color}
    
    barplot = sns.barplot(
        data=cust_long, x="Range", y="Percent", hue="Group", palette=palette,
        edgecolor="black", linewidth=0.5, ax=ax
    )
    
    # Set alpha for each bar according to group
    for patch, (_, row) in zip(ax.patches, cust_long.iterrows()):
        if row["Group"] == "Landscape %":
            patch.set_alpha(0.35)
        elif row["Group"] == "Landslide %":
            patch.set_alpha(0.75)
    
    ax.set_ylabel("coverage [%]", fontsize=font_size, fontfamily=font_family)
    ax.set_xlabel("", fontsize=font_size, fontfamily=font_family)
    ax.tick_params(axis='both', which='major', labelsize=font_size)
    ax.get_legend().remove()
    ymax2 = float(cust_long["Percent"].max())
    ax.set_ylim(0, ymax2 * 1.30)
    
    xticklabels = [label.get_text() for label in ax.get_xticklabels()]
    ax.set_xticklabels(xticklabels)
    
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(font_size)
        label.set_fontfamily(font_family)
    ax.set_ylabel(ax.get_ylabel(), fontsize=font_size, fontfamily=font_family)
    ax.set_xlabel(ax.get_xlabel(), fontsize=font_size, fontfamily=font_family)
    
    ax.grid(False)
    sns.despine(ax=ax)
    
    return fig


def plot_frequency_ratio_comparison(results: FrequencyResults, cfg) -> plt.Figure:
    """Create side-by-side D8 frequency ratio comparison matching notebook style"""
    import matplotlib as mpl
    import seaborn as sns
    from matplotlib.colors import to_rgba, to_hex
    
    mpl.rcParams['font.size'] = 6
    mpl.rcParams['font.family'] = 'Arial'
    font_size = mpl.rcParams['font.size']
    font_family = mpl.rcParams['font.family']
    
    def blend_with_white(color, alpha):
        """Blend a color with white at the given alpha, return as hex."""
        rgba = to_rgba(color)
        r, g, b, _ = rgba
        r_blend = (1 - alpha) * 1.0 + alpha * r
        g_blend = (1 - alpha) * 1.0 + alpha * g
        b_blend = (1 - alpha) * 1.0 + alpha * b
        return to_hex((r_blend, g_blend, b_blend))
    
    landslide_color = "#762A83"
    landslide_color_blend = blend_with_white(landslide_color, 0.75)
    
    fig, ax = plt.subplots(figsize=(165/25.4, 50/25.4), sharex=False, constrained_layout=True)
    
    # Prepare dataframes
    hist_d8_df = results.historical_d8.copy()
    event_d8_df = results.event_d8.copy()
    
    # Bar width and positions
    bar_width = 0.4
    bins = hist_d8_df["Bin"]
    x = np.arange(len(bins))
    
    # Plot historical (gray) bars
    bars_hist = ax.bar(
        x - bar_width/2, hist_d8_df["FR"], width=bar_width, 
        color=str(0.7), edgecolor="black", linewidth=0.5, 
        label="Historical", alpha=1, zorder=2
    )
    
    # Plot event (blended purple) bars
    bars_event = ax.bar(
        x + bar_width/2, event_d8_df["FR"], width=bar_width, 
        color=landslide_color_blend, edgecolor="black", linewidth=0.5, 
        label="Event (2022)", alpha=1, zorder=2
    )
    
    ax.set_ylabel("frequency ratio", fontsize=font_size, fontfamily=font_family)
    ax.set_xlabel("", fontsize=font_size, fontfamily=font_family)
    ax.tick_params(axis='both', which='major', labelsize=font_size)
    ymax = float(max(hist_d8_df["FR"].max(), event_d8_df["FR"].max()))
    ax.set_ylim(0, ymax * 1.25)
    ax.set_xticks(x)
    ax.set_xticklabels(bins, fontsize=font_size, fontfamily=font_family)
    sns.despine(ax=ax)
    
    # Set grid opacity to 0.3 and move grid to back
    ax.grid(True, alpha=0.3, zorder=0)
    
    # Add dashed red line at y=1
    ax.axhline(1.0, color="red", lw=1.2, ls="--", alpha=0.7, zorder=3)
    
    # Add FR value as text at the bottom of each bar for both datasets
    for i, (patch, fr) in enumerate(zip(bars_hist, hist_d8_df["FR"])):
        x_text = patch.get_x() + patch.get_width() / 2
        y = patch.get_y()
        text_y = y + 0.05 * ymax
        ax.text(
            x_text, text_y, f"{fr:.2f}",
            ha="center", va="bottom",
            fontsize=font_size, fontfamily=font_family,
            color="black", zorder=4
        )
    for i, (patch, fr) in enumerate(zip(bars_event, event_d8_df["FR"])):
        x_text = patch.get_x() + patch.get_width() / 2
        y = patch.get_y()
        text_y = y + 0.05 * ymax
        ax.text(
            x_text, text_y, f"{fr:.2f}",
            ha="center", va="bottom",
            fontsize=font_size, fontfamily=font_family,
            color="black", zorder=4
        )
    
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(font_size)
        label.set_fontfamily(font_family)
    
    return fig


def save_tables(results: FrequencyResults, cfg) -> None:
    results.event_d8.to_csv(cfg.figure_path / "aspect_frequency_event.csv", index=False)
    results.historical_d8.to_csv(cfg.figure_path / "aspect_frequency_historical.csv", index=False)
    results.coverage_table.to_csv(cfg.figure_path / "aspect_frequency_coverage.csv", index=False)


def main() -> None:
    import sys
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    if str(root) not in sys.path:
        sys.path.append(str(root))

    import config as cfg

    print("Computing aspect frequency results...")
    inventories = utils.load_inventories(cfg)
    aspect_data = _load_aspect_angles(cfg, inventories)
    results = compute_aspect_frequency_results(cfg)
    
    # Save CSV tables
    save_tables(results, cfg)
    
    # Figure 1: Rose diagram
    print("Creating rose diagram...")
    fig_rose = plot_rose_diagram(aspect_data, cfg)
    output_rose = cfg.figure_path / "aspect_rose_bars.pdf"
    fig_rose.savefig(output_rose, dpi=300)
    plt.close(fig_rose)
    print(f"  Rose diagram saved to: {output_rose}")
    
    # Figure 2: Coverage bars (N-E vs coverage)
    print("Creating coverage bar plot...")
    fig_coverage = plot_coverage_bars(results, cfg)
    output_coverage = cfg.figure_path / "aspect_bars_select.pdf"
    fig_coverage.savefig(output_coverage, dpi=300)
    plt.close(fig_coverage)
    print(f"  Coverage bars saved to: {output_coverage}")
    
    # Figure 3: Frequency ratio comparison
    print("Creating frequency ratio comparison...")
    fig_fr = plot_frequency_ratio_comparison(results, cfg)
    output_fr = cfg.figure_path / "aspect_bars_d8_historical_event.pdf"
    fig_fr.savefig(output_fr, dpi=300)
    plt.close(fig_fr)
    print(f"  Frequency ratio comparison saved to: {output_fr}")
    
    print("\nAspect analysis complete!")


if __name__ == "__main__":
    main()
