from dataclasses import dataclass

import geopandas as gpd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

import utils


@dataclass(slots=True)
class ReactivationResults:
    cumulative: pd.DataFrame


def _nearest_distances(hist: gpd.GeoDataFrame, event: gpd.GeoDataFrame) -> np.ndarray:
    """Compute nearest distance from each historical landslide to any event landslide"""
    try:
        union = event.geometry.union_all()
    except AttributeError:
        # Fallback for older geopandas versions
        union = event.geometry.unary_union
    return np.array([geom.distance(union) for geom in hist.geometry], dtype=float)


def _cumulative_table(distances: np.ndarray, areas: np.ndarray) -> pd.DataFrame:
    """Build cumulative reactivation table with mean/median historical areas"""
    thresholds = np.array([0, 1, 2, 3, 5, 7, 10, 15, 20, 30, 50, 75, 100, 150, 200, 300, 400, 500, 700, 1000])
    total = distances.size
    rows = []
    for d in thresholds:
        mask = distances <= d
        count = mask.sum()
        pct = (count / total) * 100.0 if total else 0.0
        areas_masked = areas[mask]
        mean_area = float(areas_masked.mean()) if count else float("nan")
        median_area = float(np.median(areas_masked)) if count else float("nan")
        rows.append(
            {
                "distance_m": float(d),
                "reactivated_count": int(count),
                "total_historical": int(total),
                "cumulative_percentage": pct,
                "mean_hist_area_m2": mean_area,
                "median_hist_area_m2": median_area,
            }
        )
    return pd.DataFrame(rows)


def compute_reactivation_results(cfg) -> ReactivationResults:
    inventories = utils.load_inventories(cfg)
    distances = _nearest_distances(inventories.historical, inventories.event)
    cumulative = _cumulative_table(distances, inventories.historical["area_m2"].to_numpy())
    return ReactivationResults(cumulative=cumulative)


def plot_reactivation_curve(results: ReactivationResults, cfg) -> plt.Figure:
    """Create reactivation distance curve matching notebook Cell 7 style"""
    mpl.rcParams['font.size'] = 8
    mpl.rcParams['font.family'] = 'Arial'
    
    fig, ax1 = plt.subplots(figsize=(150/25.4, 60/25.4))
    
    # Add shaded zones for different reactivation types
    ax1.axvspan(0, 1, alpha=0.2, color='#a6dba0')
    ax1.axvspan(1, 10, alpha=0.2, color='#80cdc1')
    ax1.axvspan(10, 1000, alpha=0.2, color='#bf812d')
    
    # Plot the main data line, using 0.00001 instead of 0 for log scale
    distances_for_plot = results.cumulative['distance_m'].copy()
    distances_for_plot.iloc[0] = 0.00001  # Replace 0 with a small positive value for log scale
    ax1.semilogx(
        distances_for_plot, 
        results.cumulative['cumulative_percentage'], 
        '-', linewidth=1.5, color='#762A83', zorder=5, label='cumulative % reactivated'
    )
    
    # Add boundary lines
    ax1.axvline(1, color='#a6dba0', linestyle='-', linewidth=1, alpha=0.8, zorder=4)
    ax1.axvline(10, color='#80cdc1', linestyle='-', linewidth=1, alpha=0.8, zorder=4)
    ax1.axvline(1000, color='#bf812d', linestyle='-', linewidth=1, alpha=0.8, zorder=4)
    
    # Add zone labels
    label_y = 8
    ax1.text(
        0.49, label_y, 'DIRECT', fontsize=8, fontweight='bold', ha='center', va='center',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#a6dba0', alpha=0.3)
    )
    ax1.text(
        3, label_y, 'LOCAL', fontsize=8, fontweight='bold', ha='center', va='center',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#80cdc1', alpha=0.3)
    )
    ax1.text(
        40, label_y, 'REMOTE', fontsize=8, fontweight='bold', ha='center', va='center',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#bf812d', alpha=0.3)
    )
    
    # Add annotations for % reactivated at x=1 and x=10
    pct_1 = results.cumulative[results.cumulative['distance_m'] == 1]['cumulative_percentage'].iloc[0]
    pct_10 = results.cumulative[results.cumulative['distance_m'] == 10]['cumulative_percentage'].iloc[0]
    
    ax1.annotate(
        f"{pct_1:.1f}%", 
        xy=(1, pct_1), 
        xytext=(1.2, pct_1 + 7),
        textcoords='data',
        ha='left', va='bottom',
        fontsize=9, fontweight='bold',
        color='red',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='red', alpha=0.8),
        arrowprops=dict(arrowstyle='->', color='red', lw=1)
    )
    ax1.annotate(
        f"{pct_10:.1f}%", 
        xy=(10, pct_10), 
        xytext=(13, pct_10 + 7),
        textcoords='data',
        ha='left', va='bottom',
        fontsize=9, fontweight='bold',
        color='orange',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='orange', alpha=0.8),
        arrowprops=dict(arrowstyle='->', color='orange', lw=1)
    )
    
    ax1.set_xlabel('distance from historical landslide (m)', fontsize=8)
    ax1.set_ylabel('reactivated historical landslides (%)', fontsize=8)
    ax1.grid(True, which="both", alpha=0.3)
    ax1.set_xlim(0.25, 1000)
    ax1.set_ylim(0, 105)
    
    # Right axis: historical landslide size stats (m²)
    ax2 = ax1.twinx()
    mean_series = results.cumulative['mean_hist_area_m2']
    median_series = results.cumulative['median_hist_area_m2']
    ax2.set_ylabel('historical landslide size (m²)', fontsize=8)
    ax2.semilogx(distances_for_plot, mean_series, '-', color='#1b9e77', linewidth=1.2, label='mean area')
    ax2.semilogx(distances_for_plot, median_series, '--', color='#1b9e77', linewidth=1.2, label='median area')
    
    # y-limit based on both series
    try:
        ymax = np.nanmax([np.nanmax(mean_series.values), np.nanmax(median_series.values)])
        ax2.set_ylim(0, ymax * 1.1 if np.isfinite(ymax) else 1)
    except Exception:
        pass
    
    plt.tight_layout()
    
    return fig


def main() -> None:
    import sys
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    if str(root) not in sys.path:
        sys.path.append(str(root))

    import config as cfg

    print("Computing reactivation distance metrics...")
    results = compute_reactivation_results(cfg)
    
    # Save cumulative table
    results.cumulative.to_csv(cfg.figure_path / "reactivation_cumulative.csv", index=False)
    
    # Create and save figure
    fig = plot_reactivation_curve(results, cfg)
    output = cfg.figure_path / "cumulative_distance_curve_with_size.pdf"
    fig.savefig(output, format="pdf", bbox_inches='tight', dpi=300)
    plt.close(fig)
    
    # Print key insights
    print(f"\nKey Insights:")
    print(f"• Direct reactivation (0m): {results.cumulative.iloc[0]['cumulative_percentage']:.1f}%")
    print(f"• Within 10m: {results.cumulative[results.cumulative['distance_m'] == 10]['cumulative_percentage'].iloc[0]:.1f}%")
    print(f"• Within 100m: {results.cumulative[results.cumulative['distance_m'] == 100]['cumulative_percentage'].iloc[0]:.1f}%")
    print(f"• Within 1000m: {results.cumulative[results.cumulative['distance_m'] == 1000]['cumulative_percentage'].iloc[0]:.1f}%")
    
    print(f"\nReactivation distance figure saved to: {output}")


if __name__ == "__main__":
    main()
