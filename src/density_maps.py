from dataclasses import dataclass
from typing import Dict

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
import shapely
from scipy.stats import gaussian_kde

import utils


@dataclass(slots=True)
class DensitySurfaces:
    grid_x: np.ndarray
    grid_y: np.ndarray
    count_event: np.ndarray
    count_hist: np.ndarray
    area_event: np.ndarray
    area_hist: np.ndarray


def _build_grid(mapping_area: gpd.GeoDataFrame, grid_size: int = 500) -> Dict[str, np.ndarray]:
    xmin, ymin, xmax, ymax = mapping_area.total_bounds
    grid_x, grid_y = np.mgrid[xmin:xmax:complex(0, grid_size), ymin:ymax:complex(0, grid_size)]
    return {
        "X": grid_x,
        "Y": grid_y,
        "coords": np.vstack([grid_x.ravel(), grid_y.ravel()]),
    }


def _apply_area_mask(data: np.ndarray, mask: np.ndarray) -> np.ndarray:
    out = data.copy()
    out[~mask] = np.nan
    return out


def compute_density_surfaces(cfg, grid_size: int = 500) -> DensitySurfaces:
    inv = utils.load_inventories(cfg)
    mapping_area = utils.load_mapping_area(cfg)

    grid = _build_grid(mapping_area, grid_size=grid_size)
    coords = grid["coords"]

    def _kde(points: np.ndarray, weights=None):
        kde = gaussian_kde(points, bw_method=cfg.kde_method, weights=weights)
        vals = kde(coords).reshape(grid["X"].shape)
        return vals

    def _extract_xy(gdf: gpd.GeoDataFrame) -> np.ndarray:
        centroids = gdf.geometry.centroid
        return np.vstack([centroids.x.values, centroids.y.values])

    xy_event = _extract_xy(inv.event)
    xy_hist = _extract_xy(inv.historical)

    count_event = _kde(xy_event)
    count_hist = _kde(xy_hist)

    weights_event = inv.event["area_m2"].values / 1e6
    weights_hist = inv.historical["area_m2"].values / 1e6
    area_event = _kde(xy_event, weights=weights_event)
    area_hist = _kde(xy_hist, weights=weights_hist)

    # Scale density surfaces
    # grid spacing (assumes regular grid)
    dx = (grid["X"].max() - grid["X"].min()) / (grid["X"].shape[0] - 1)
    dy = (grid["Y"].max() - grid["Y"].min()) / (grid["Y"].shape[1] - 1)
    pixel_area_km2 = (dx * dy) / 1e6

    count_event_scaled = count_event * len(inv.event) / (np.nansum(count_event) * pixel_area_km2)
    count_hist_scaled = count_hist * len(inv.historical) / (np.nansum(count_hist) * pixel_area_km2)

    total_area_event = weights_event.sum()
    total_area_hist = weights_hist.sum()
    area_event_scaled = area_event * total_area_event / (np.nansum(area_event) * pixel_area_km2)
    area_hist_scaled = area_hist * total_area_hist / (np.nansum(area_hist) * pixel_area_km2)

    # Build mask
    try:
        mask_geom = mapping_area.geometry.union_all()
    except AttributeError:
        # Fallback for older geopandas versions
        mask_geom = mapping_area.geometry.unary_union
    pts = shapely.points(coords[0], coords[1])
    mask_flat = shapely.contains(mask_geom, pts)
    flat_mask = mask_flat.reshape(grid["X"].shape)

    count_event_scaled = _apply_area_mask(count_event_scaled, flat_mask)
    count_hist_scaled = _apply_area_mask(count_hist_scaled, flat_mask)
    area_event_scaled = _apply_area_mask(area_event_scaled, flat_mask)
    area_hist_scaled = _apply_area_mask(area_hist_scaled, flat_mask)

    return DensitySurfaces(
        grid_x=grid["X"],
        grid_y=grid["Y"],
        count_event=count_event_scaled,
        count_hist=count_hist_scaled,
        area_event=area_event_scaled,
        area_hist=area_hist_scaled,
    )


def _add_scalebar(ax, length_m, x_pos, y_pos, map_width_m):
    """Add a simple scalebar to the axis"""
    from matplotlib.patches import Rectangle
    from matplotlib.lines import Line2D
    
    # Calculate scalebar width in axis coordinates
    scalebar_width = length_m / map_width_m
    
    # Add black background rectangle
    rect = Rectangle(
        (x_pos, y_pos), scalebar_width, 0.015,
        transform=ax.transAxes,
        facecolor='black',
        edgecolor='white',
        linewidth=1,
        zorder=10
    )
    ax.add_patch(rect)
    
    # Add text label
    ax.text(
        x_pos + scalebar_width / 2, y_pos + 0.025,
        f'{length_m/1000:.0f} km',
        transform=ax.transAxes,
        ha='center', va='bottom',
        fontsize=8, color='white',
        fontweight='bold',
        zorder=11
    )


def plot_density_figures(surfaces: DensitySurfaces, cfg) -> tuple[plt.Figure, plt.Figure]:
    """Create two 1x2 figures: one for event data, one for historical data"""
    mapping_area = utils.load_mapping_area(cfg)
    inv = utils.load_inventories(cfg)
    
    xmin, ymin, xmax, ymax = mapping_area.total_bounds
    map_width_m = xmax - xmin
    
    # Normalize density surfaces to 0-1
    def normalize(data):
        return (data - np.nanmin(data)) / (np.nanmax(data) - np.nanmin(data))
    
    count_event_norm = normalize(surfaces.count_event)
    count_hist_norm = normalize(surfaces.count_hist)
    area_event_norm = normalize(surfaces.area_event)
    area_hist_norm = normalize(surfaces.area_hist)
    
    # Figure 1: Event (2022) density maps
    fig_event, axs_event = plt.subplots(1, 2, figsize=(170/25.4, 125/25.4))
    
    # Event count density
    ax = axs_event[0]
    im = ax.imshow(
        np.rot90(count_event_norm),
        cmap='inferno',
        extent=[xmin, xmax, ymin, ymax],
        alpha=0.8,
        zorder=2,
        vmin=0,
        vmax=1
    )
    inv.event.boundary.plot(ax=ax, edgecolor='white', linewidth=0.2, alpha=0.4, zorder=3, rasterized=True)
    fig_event.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
    ax.set_axis_off()
    _add_scalebar(ax, 5000, 0.03, 0.05, map_width_m)
    
    # Event area fraction
    ax = axs_event[1]
    im = ax.imshow(
        np.rot90(area_event_norm),
        cmap='inferno',
        extent=[xmin, xmax, ymin, ymax],
        alpha=0.8,
        zorder=2,
        vmin=0,
        vmax=1
    )
    inv.event.boundary.plot(ax=ax, edgecolor='white', linewidth=0.2, alpha=0.4, zorder=3, rasterized=True)
    fig_event.colorbar(im, ax=ax, fraction=0.03, pad=0.04, label='percentile')
    ax.set_axis_off()
    _add_scalebar(ax, 5000, 0.03, 0.05, map_width_m)
    
    # Figure 2: Historical density maps
    fig_hist, axs_hist = plt.subplots(1, 2, figsize=(170/25.4, 125/25.4))
    
    # Historical count density
    ax = axs_hist[0]
    im = ax.imshow(
        np.rot90(count_hist_norm),
        cmap='inferno',
        extent=[xmin, xmax, ymin, ymax],
        alpha=0.8,
        zorder=2,
        vmin=0,
        vmax=1
    )
    inv.event.boundary.plot(ax=ax, edgecolor='white', linewidth=0.2, alpha=0.4, zorder=3, rasterized=True)
    fig_hist.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
    ax.set_axis_off()
    _add_scalebar(ax, 5000, 0.03, 0.05, map_width_m)
    
    # Historical area fraction
    ax = axs_hist[1]
    im = ax.imshow(
        np.rot90(area_hist_norm),
        cmap='inferno',
        extent=[xmin, xmax, ymin, ymax],
        alpha=0.8,
        zorder=2,
        vmin=0,
        vmax=1
    )
    inv.event.boundary.plot(ax=ax, edgecolor='white', linewidth=0.2, alpha=0.4, zorder=3, rasterized=True)
    fig_hist.colorbar(im, ax=ax, fraction=0.03, pad=0.04, label='percentile')
    ax.set_axis_off()
    _add_scalebar(ax, 5000, 0.03, 0.05, map_width_m)
    
    return fig_event, fig_hist


def main() -> None:
    import sys
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    if str(root) not in sys.path:
        sys.path.append(str(root))

    import config as cfg

    print("Computing density surfaces...")
    surfaces = compute_density_surfaces(cfg)
    
    print("Creating density figures...")
    fig_event, fig_hist = plot_density_figures(surfaces, cfg)
    
    # Save event figure
    output_event = cfg.figure_path / "density-maps-event.pdf"
    fig_event.savefig(output_event, format="pdf", bbox_inches=None, pad_inches=0, dpi=300)
    plt.close(fig_event)
    print(f"Event density figure saved to: {output_event}")
    
    # Save historical figure
    output_hist = cfg.figure_path / "density-maps-historical.pdf"
    fig_hist.savefig(output_hist, format="pdf", bbox_inches=None, pad_inches=0, dpi=300)
    plt.close(fig_hist)
    print(f"Historical density figure saved to: {output_hist}")


if __name__ == "__main__":
    main()


