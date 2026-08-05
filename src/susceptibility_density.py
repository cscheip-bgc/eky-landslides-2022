"""Susceptibility vs. observed landslide density across the mapping area.

Tests whether the spatial density patterns shown in the density-map figures can
be explained by spatial variation in terrain predisposition. An independent,
pre-existing susceptibility model and the mapped event inventory are aggregated
onto a common grid of fixed-size cells; if density patterns were driven by
predisposing factors, cells with high landslide density would also carry
systematically higher susceptibility.

Produces tabular output only - no figure. A susceptibility map at the scale of
the 525 km2 mapping area would ask the reader to perceive the *absence* of a
correlation by eye, which maps communicate poorly.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import rasterio
from scipy.stats import pearsonr, spearmanr
import sys

# Project config
sys.path.append(str(Path(__file__).resolve().parent.parent))
import config as cfg

import utils


# ---------------------------------------------------------------------------
# Grid definition
# ---------------------------------------------------------------------------

def _grid_spec(src, cell_size):
    """
    Define a cell grid aligned to round projected coordinates.

    Anchoring to multiples of ``cell_size`` (rather than to the raster corner)
    means the grid is reproducible from the CRS alone and can be described in
    the manuscript as a plain 1-km UTM grid.

    Returns
    -------
    dict
        Grid origin, cell counts, and pixel size
    """
    b = src.bounds
    x_anchor = np.floor(b.left / cell_size) * cell_size
    y_anchor = np.ceil(b.top / cell_size) * cell_size

    n_x = int(np.ceil((b.right - x_anchor) / cell_size))
    n_y = int(np.ceil((y_anchor - b.bottom) / cell_size))

    return {
        'x_anchor': x_anchor,
        'y_anchor': y_anchor,
        'n_x': n_x,
        'n_y': n_y,
        'n_cells': n_x * n_y,
        'pixel_size': abs(src.transform.a),
        'cell_size': cell_size,
    }


def _cell_ids(xs, ys, grid):
    """Map projected coordinates to flat cell indices (row-major)."""
    cx = np.floor((xs - grid['x_anchor']) / grid['cell_size']).astype(np.int64)
    cy = np.floor((grid['y_anchor'] - ys) / grid['cell_size']).astype(np.int64)
    valid = (cx >= 0) & (cx < grid['n_x']) & (cy >= 0) & (cy < grid['n_y'])
    return cy * grid['n_x'] + cx, valid


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def aggregate_susceptibility(cell_size=None, thresholds=None,
                             susceptibility_file=None):
    """
    Accumulate susceptibility statistics per grid cell at native resolution.

    The raster is read block by block rather than in full: at 9,173 x 16,182
    float32 a whole-array read is ~594 MB before any float64 promotion.

    Parameters
    ----------
    cell_size : float, optional
        Grid cell size in metres (defaults to ``config.SUSC_CELL_SIZE_M``)
    thresholds : sequence of float, optional
        Susceptibility values above which to count pixels (defaults to
        ``config.SUSC_HIGH_THRESHOLDS``)
    susceptibility_file : str or Path, optional
        Path to the susceptibility raster (defaults to
        ``config.susceptibility_file``)

    Returns
    -------
    tuple
        (DataFrame indexed by flat cell id, grid spec dict)
    """
    if cell_size is None:
        cell_size = cfg.SUSC_CELL_SIZE_M
    if thresholds is None:
        thresholds = cfg.SUSC_HIGH_THRESHOLDS
    if susceptibility_file is None:
        susceptibility_file = cfg.susceptibility_file

    with rasterio.open(susceptibility_file) as src:
        grid = _grid_spec(src, cell_size)
        n = grid['n_cells']

        n_px = np.zeros(n, dtype=np.int64)
        sum_v = np.zeros(n, dtype=np.float64)
        hi = {t: np.zeros(n, dtype=np.int64) for t in thresholds}

        px = grid['pixel_size']
        x0, y1 = src.transform.c, src.transform.f

        for _, window in src.block_windows(1):
            block = src.read(1, window=window)

            col0, row0 = window.col_off, window.row_off
            h, w = block.shape
            xs = x0 + (col0 + np.arange(w) + 0.5) * px
            ys = y1 - (row0 + np.arange(h) + 0.5) * px

            ok = np.isfinite(block) & (block > -1e30)
            if src.nodata is not None:
                ok &= block != src.nodata
            if not ok.any():
                continue

            xx = np.broadcast_to(xs, (h, w))[ok]
            yy = np.broadcast_to(ys[:, None], (h, w))[ok]
            vals = block[ok].astype(np.float64)

            cid, in_grid = _cell_ids(xx, yy, grid)
            cid, vals = cid[in_grid], vals[in_grid]

            n_px += np.bincount(cid, minlength=n)
            sum_v += np.bincount(cid, weights=vals, minlength=n)
            for t in thresholds:
                hi[t] += np.bincount(cid, weights=(vals > t), minlength=n).astype(np.int64)

    out = pd.DataFrame({'n_px': n_px, 'sum_susc': sum_v})
    for t in thresholds:
        out[f'n_hi{int(round(t * 100))}'] = hi[t]
    return out, grid


def aggregate_landslides(event, grid):
    """
    Accumulate landslide count and area per grid cell.

    Uses polygon centroids, matching how ``density_maps.compute_density_surfaces``
    builds the KDE surfaces, so this compares against the same quantity the
    density figures display.

    Parameters
    ----------
    event : gpd.GeoDataFrame
        Event inventory in the target CRS, with ``area_m2`` populated
    grid : dict
        Grid spec from :func:`aggregate_susceptibility`

    Returns
    -------
    pd.DataFrame
        Per-cell landslide counts and summed areas
    """
    centroids = event.geometry.centroid
    cid, in_grid = _cell_ids(centroids.x.values, centroids.y.values, grid)

    n = grid['n_cells']
    counts = np.bincount(cid[in_grid], minlength=n)
    areas = np.bincount(cid[in_grid], weights=event['area_m2'].values[in_grid],
                        minlength=n)

    return pd.DataFrame({'n_landslides': counts, 'landslide_area_m2': areas})


def build_cell_table(susc, landslides, grid, min_coverage=None):
    """
    Join the aggregations and derive per-cell rates, dropping edge cells.

    Parameters
    ----------
    susc, landslides : pd.DataFrame
        Outputs of the two aggregation functions
    grid : dict
        Grid spec
    min_coverage : float, optional
        Minimum fraction of a cell that must carry valid raster data
        (defaults to ``config.SUSC_MIN_CELL_COVERAGE``)

    Returns
    -------
    pd.DataFrame
        One row per retained cell
    """
    if min_coverage is None:
        min_coverage = cfg.SUSC_MIN_CELL_COVERAGE

    px_per_cell = (grid['cell_size'] / grid['pixel_size']) ** 2

    d = susc.join(landslides)
    d['coverage'] = d.n_px / px_per_cell
    d = d[d.coverage >= min_coverage].copy()

    cell_area_km2 = d.n_px * grid['pixel_size'] ** 2 / 1e6
    d['mean_susc'] = d.sum_susc / d.n_px
    for col in [c for c in d.columns if c.startswith('n_hi')]:
        d[f'frac_hi{col[4:]}'] = d[col] / d.n_px
    d['landslides_km2'] = d.n_landslides / cell_area_km2
    d['landslide_area_frac'] = d.landslide_area_m2 / (cell_area_km2 * 1e6)

    d.index.name = 'cell_id'
    return d


# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------

def _susc_metrics(cells):
    return ['mean_susc'] + sorted(c for c in cells.columns if c.startswith('frac_hi'))


DENSITY_METRICS = ['landslides_km2', 'landslide_area_frac']


def summarize_dispersion(cells):
    """
    Report spread of each metric across cells.

    The coefficient of variation and the 90:10 ratio are the comparison of
    interest: susceptibility should be far less variable across the mapping
    area than observed landslide density.
    """
    rows = []
    for m in _susc_metrics(cells) + DENSITY_METRICS:
        x = cells[m]
        p10, p90 = np.percentile(x, [10, 90])
        rows.append({
            'Metric': m,
            'Mean': round(float(x.mean()), 4),
            'SD': round(float(x.std()), 4),
            'CV': round(float(x.std() / x.mean()), 3),
            'Min': round(float(x.min()), 4),
            'P10': round(float(p10), 4),
            'P90': round(float(p90), 4),
            'Max': round(float(x.max()), 4),
            'P90_over_P10': round(float(p90 / p10), 1) if p10 > 0 else np.inf,
        })
    return pd.DataFrame(rows)


def correlate(cells):
    """Correlate every susceptibility metric against every density metric."""
    rows = []
    for m in _susc_metrics(cells):
        for t in DENSITY_METRICS:
            r, p = pearsonr(cells[m], cells[t])
            rho, p_s = spearmanr(cells[m], cells[t])
            rows.append({
                'Susceptibility_metric': m,
                'Density_metric': t,
                'Pearson_r': round(float(r), 3),
                'R2': round(float(r ** 2), 4),
                'Pearson_p': float(f'{p:.2e}'),
                'Spearman_rho': round(float(rho), 3),
                'Spearman_p': float(f'{p_s:.2e}'),
            })
    return pd.DataFrame(rows)


def summarize_by_density_quartile(cells, n_groups=4):
    """
    Mean susceptibility within landslide-density quartiles of cells.

    This is the response-to-reviewers table: if predisposition drove the density
    pattern, susceptibility would rise monotonically and materially from the
    lowest to the highest density quartile.
    """
    labels = ['Q1 (low)', 'Q2', 'Q3', 'Q4 (high)'][:n_groups]
    d = cells.copy()
    d['Density_quartile'] = pd.qcut(d.landslides_km2, n_groups, labels=labels,
                                    duplicates='drop')

    agg = {'N_cells': ('mean_susc', 'size'),
           'Landslides_km2': ('landslides_km2', 'mean'),
           'Mean_susceptibility': ('mean_susc', 'mean'),
           'SD_susceptibility': ('mean_susc', 'std')}
    for c in sorted(c for c in cells.columns if c.startswith('frac_hi')):
        agg[c] = (c, 'mean')

    out = d.groupby('Density_quartile', observed=True).agg(**agg).reset_index()
    return out.round(4)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def run_for_cell_size(event, cell_size, label):
    """Run the full comparison at one cell size and print the tables."""
    susc, grid = aggregate_susceptibility(cell_size=cell_size)
    ls = aggregate_landslides(event, grid)
    cells = build_cell_table(susc, ls, grid)

    valid_km2 = susc.n_px.sum() * grid['pixel_size'] ** 2 / 1e6
    retained = int(cells.n_landslides.sum())
    print(f"\n{'=' * 78}\n{label}: {len(cells)} retained cells "
          f"({cell_size / 1000:.0f} km), {retained} of {len(event)} landslides")
    print(f"valid raster extent: {valid_km2:.1f} km2\n{'=' * 78}")

    dispersion = summarize_dispersion(cells)
    corr = correlate(cells)
    quartiles = summarize_by_density_quartile(cells)

    with pd.option_context('display.width', 200, 'display.max_columns', None):
        print('\nDispersion across cells')
        print(dispersion.to_string(index=False))
        print('\nSusceptibility vs. density correlations')
        print(corr.to_string(index=False))
        print('\nSusceptibility by landslide-density quartile')
        print(quartiles.to_string(index=False))

    return cells, dispersion, corr, quartiles


def main():
    """Run the susceptibility-density comparison and write the tables."""
    output_dir = utils.ensure_output_dir(cfg.figure_path)

    inv = utils.load_inventories(cfg)
    event = inv.event
    print(f"Loaded {len(event)} event landslides")

    mapping_area = utils.load_mapping_area(cfg)
    print(f"Mapping area: {mapping_area.area.sum() / 1e6:.1f} km2")

    cells, dispersion, corr, quartiles = run_for_cell_size(
        event, cfg.SUSC_CELL_SIZE_M, 'PRIMARY')

    print('\n\nSensitivity check at a coarser cell size')
    run_for_cell_size(event, cfg.SUSC_SENSITIVITY_CELL_SIZE_M, 'SENSITIVITY')

    summary = pd.concat([
        dispersion.assign(Table='dispersion'),
        corr.assign(Table='correlation'),
    ], ignore_index=True)

    cells_path = output_dir / 'susceptibility_density_cells.csv'
    summary_path = output_dir / 'susceptibility_density_summary.csv'
    quartile_path = output_dir / 'susceptibility_by_density_quartile.csv'
    cells.to_csv(cells_path)
    summary.to_csv(summary_path, index=False)
    quartiles.to_csv(quartile_path, index=False)

    print(f"\nSaved per-cell table to {cells_path}")
    print(f"Saved summary statistics to {summary_path}")
    print(f"Saved density-quartile table to {quartile_path}")


if __name__ == '__main__':
    main()
