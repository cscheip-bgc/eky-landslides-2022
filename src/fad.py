from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

import utils


@dataclass(slots=True)
class FADCurves:
    bin_centers_event: np.ndarray
    bin_centers_hist: np.ndarray
    event_density: np.ndarray
    historical_density: np.ndarray
    event_pdf: np.ndarray
    historical_pdf: np.ndarray
    event_params: tuple
    historical_params: tuple
    powerlaw_slope: float
    powerlaw_intercept: float
    powerlaw_beta: float
    n_event: int
    n_hist: int
    max_area_event: float
    max_area_hist: float


def _compute_fad(values: np.ndarray, n_bins: int) -> tuple[np.ndarray, np.ndarray]:
    """Compute frequency-area distribution matching kgs.compute_fad"""
    min_val = values.min()
    max_val = values.max()
    bins = np.logspace(np.log10(min_val), np.log10(max_val), n_bins)
    counts, bin_edges = np.histogram(values, bins=bins, density=True)
    bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2
    return counts, bin_centers


def _fit_invgamma_from_raw(values: np.ndarray, n_bins: int) -> tuple[tuple, np.ndarray, np.ndarray]:
    """Fit inverse gamma to raw data matching kgs.fit_invgamma_from_raw"""
    params = stats.invgamma.fit(values)
    bins = np.logspace(np.log10(values.min()), np.log10(values.max()), n_bins)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    pdf = stats.invgamma.pdf(bin_centers, *params)
    return params, pdf, bin_centers


def compute_fad_curves(cfg) -> FADCurves:
    inventories = utils.load_inventories(cfg)

    event_area = inventories.event["area_m2"].values
    hist_area = inventories.historical["area_m2"].values

    # Compute FADs
    counts_event, bin_centers_event = _compute_fad(event_area, cfg.n_bins)
    counts_hist, bin_centers_hist = _compute_fad(hist_area, cfg.n_bins)

    # Fit inverse gamma distributions
    params_event, pdf_event, _ = _fit_invgamma_from_raw(event_area, cfg.n_bins)
    params_hist, pdf_hist, _ = _fit_invgamma_from_raw(hist_area, cfg.n_bins)

    # Compute power law fit on event data
    mask = (
        (bin_centers_event >= cfg.AREA_MIN)
        & (bin_centers_event <= cfg.AREA_MAX)
        & (counts_event > 0)
    )
    
    if mask.sum() >= 3:
        log_areas = np.log10(bin_centers_event[mask])
        log_freq = np.log10(counts_event[mask])
        slope, intercept, r_value, _, _ = stats.linregress(log_areas, log_freq)
        beta = -slope
    else:
        slope = intercept = beta = np.nan

    # Find max frequency areas
    max_idx_event = np.argmax(pdf_event)
    max_area_event = bin_centers_event[max_idx_event]
    
    max_idx_hist = np.argmax(pdf_hist)
    max_area_hist = bin_centers_hist[max_idx_hist]

    return FADCurves(
        bin_centers_event=bin_centers_event,
        bin_centers_hist=bin_centers_hist,
        event_density=counts_event,
        historical_density=counts_hist,
        event_pdf=pdf_event,
        historical_pdf=pdf_hist,
        event_params=params_event,
        historical_params=params_hist,
        powerlaw_slope=slope,
        powerlaw_intercept=intercept,
        powerlaw_beta=beta,
        n_event=len(event_area),
        n_hist=len(hist_area),
        max_area_event=max_area_event,
        max_area_hist=max_area_hist,
    )


def plot_fad(curves: FADCurves, cfg) -> plt.Figure:
    """Create FAD plot matching notebook Cell 6 style"""
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 8
    
    fig, ax = plt.subplots(figsize=(80/25.4, 60/25.4))
    
    # Scatter plots (semi-transparent)
    ax.scatter(curves.bin_centers_event, curves.event_density, 
               color="#762A83", marker="o", zorder=2.5, s=15, alpha=0.3)
    ax.scatter(curves.bin_centers_hist, curves.historical_density, 
               color="#1B7837", marker="o", zorder=2.5, s=15, alpha=0.3)
    
    # Inverse gamma fit lines
    ax.plot(curves.bin_centers_event, curves.event_pdf, 
            color="#762A83", linestyle="-", zorder=2.5, 
            label=f"2022 (n={curves.n_event})")
    ax.plot(curves.bin_centers_hist, curves.historical_pdf, 
            color="#1B7837", linestyle="-", zorder=2.5, 
            label=f"Historical (n={curves.n_hist})")
    
    # Power law line
    if not np.isnan(curves.powerlaw_slope):
        areas_powerlaw = np.logspace(np.log10(cfg.AREA_MIN), np.log10(cfg.AREA_MAX), 100)
        freq_powerlaw = 10**(curves.powerlaw_intercept) * areas_powerlaw**(curves.powerlaw_slope)
        ax.plot(areas_powerlaw, freq_powerlaw, 
                color="black", linestyle="--", linewidth=1, zorder=5, 
                label=f"Power Law: β = {curves.powerlaw_beta:.2f}")
    
    # Vertical lines at max frequency areas
    ax.axvline(x=curves.max_area_event, color="#762A83", linestyle="--", zorder=2.5, alpha=0.5)
    ax.axvline(x=curves.max_area_hist, color="#1B7837", linestyle="--", zorder=2.5, alpha=0.5)
    
    # Formatting
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("area (m²)")
    ax.set_ylabel("frequency (m$^{-2}$)")
    ax.grid(True, which="major", ls="-", zorder=0)
    ax.set_ylim(top=1e-1)
    ax.set_xlim([1e0, 1e5])
    
    plt.tight_layout()
    
    return fig


def main() -> None:
    import sys
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    if str(root) not in sys.path:
        sys.path.append(str(root))

    import config as cfg

    print("Computing frequency-area curves...")
    curves = compute_fad_curves(cfg)
    
    # Print summary statistics
    print(f"\n2022: shape={curves.event_params[0]:.2f}, loc={curves.event_params[1]:.2f}, scale={curves.event_params[2]:.2f}")
    print(f"Historical: shape={curves.historical_params[0]:.2f}, loc={curves.historical_params[1]:.2f}, scale={curves.historical_params[2]:.2f}")
    print(f"\n2022: {curves.n_event}")
    print(f"Historical: {curves.n_hist}")
    print(f"\n2022 max frequency area: {curves.max_area_event:.2f} m²")
    print(f"Historical max frequency area: {curves.max_area_hist:.2f} m²")
    print(f"\nPower Law: β = {curves.powerlaw_beta:.2f}")
    
    # Save figure
    fig = plot_fad(curves, cfg)
    output_pdf = cfg.figure_path / "fad_2022_vs_historical.pdf"
    fig.savefig(output_pdf, bbox_inches='tight')
    plt.close(fig)
    print(f"\nFAD figure saved to: {output_pdf}")


if __name__ == "__main__":
    main()
