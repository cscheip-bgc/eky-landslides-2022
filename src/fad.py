"""Frequency-area distribution (FAD) analysis and visualization.

Produces a two-panel figure:
  (a) Event and historical FADs with inverse-gamma fits and power-law line
  (b) Mapping deficit detail: extrapolated power law vs mapped data in the
      rollover-to-cutoff range, with shaded deficit region
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

import utils


@dataclass(slots=True)
class FADCurves:
    bin_centers_event: np.ndarray
    bin_centers_hist: np.ndarray
    bins_event: np.ndarray          # raw bin edges for count recovery
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
    event_areas: np.ndarray         # raw area values for deficit calc


def _compute_fad(values: np.ndarray, n_bins: int):
    """Binned frequency-area distribution (density-normalized)."""
    min_val = values.min()
    max_val = values.max()
    bins = np.logspace(np.log10(min_val), np.log10(max_val), n_bins)
    counts, bin_edges = np.histogram(values, bins=bins, density=True)
    bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2
    return counts, bin_centers, bins


def _fit_invgamma_from_raw(values: np.ndarray, n_bins: int):
    """Fit inverse gamma to raw data."""
    params = stats.invgamma.fit(values)
    bins = np.logspace(np.log10(values.min()), np.log10(values.max()), n_bins)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    pdf = stats.invgamma.pdf(bin_centers, *params)
    return params, pdf, bin_centers


def compute_fad_curves(cfg) -> FADCurves:
    inventories = utils.load_inventories(cfg)

    event_area = inventories.event["area_m2"].values
    hist_area = inventories.historical["area_m2"].values

    counts_event, bin_centers_event, bins_event = _compute_fad(event_area, cfg.n_bins)
    counts_hist, bin_centers_hist, _ = _compute_fad(hist_area, cfg.n_bins)

    params_event, pdf_event, _ = _fit_invgamma_from_raw(event_area, cfg.n_bins)
    params_hist, pdf_hist, _ = _fit_invgamma_from_raw(hist_area, cfg.n_bins)

    # Power law fit on event data
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

    # Max frequency areas (rollover)
    max_idx_event = np.argmax(pdf_event)
    max_area_event = bin_centers_event[max_idx_event]

    max_idx_hist = np.argmax(pdf_hist)
    max_area_hist = bin_centers_hist[max_idx_hist]

    return FADCurves(
        bin_centers_event=bin_centers_event,
        bin_centers_hist=bin_centers_hist,
        bins_event=bins_event,
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
        event_areas=event_area,
    )


# ---------------------------------------------------------------------------
# Panel (a): Full FAD comparison
# ---------------------------------------------------------------------------

def _plot_panel_a(ax, curves: FADCurves, cfg):
    """Event + historical FADs with inverse-gamma and power-law fits."""
    # Scatter plots
    ax.scatter(curves.bin_centers_event, curves.event_density,
               color="#762A83", marker="o", zorder=2.5, s=15, alpha=0.3)
    ax.scatter(curves.bin_centers_hist, curves.historical_density,
               color="#1B7837", marker="o", zorder=2.5, s=15, alpha=0.3)

    # Inverse gamma fits
    ax.plot(curves.bin_centers_event, curves.event_pdf,
            color="#762A83", linestyle="-", zorder=2.5,
            label=f"2022 (n={curves.n_event})")
    ax.plot(curves.bin_centers_hist, curves.historical_pdf,
            color="#1B7837", linestyle="-", zorder=2.5,
            label=f"Historical (n={curves.n_hist})")

    # Power law line
    if not np.isnan(curves.powerlaw_slope):
        areas_pl = np.logspace(np.log10(cfg.AREA_MIN), np.log10(cfg.AREA_MAX), 100)
        freq_pl = 10**curves.powerlaw_intercept * areas_pl**curves.powerlaw_slope
        ax.plot(areas_pl, freq_pl,
                color="black", linestyle="--", linewidth=1, zorder=5,
                label=f"Power Law: \u03b2 = {curves.powerlaw_beta:.2f}")

    # Vertical rollover lines
    ax.axvline(x=curves.max_area_event, color="#762A83", linestyle="--",
               zorder=2.5, alpha=0.5)
    ax.axvline(x=curves.max_area_hist, color="#1B7837", linestyle="--",
               zorder=2.5, alpha=0.5)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("area (m\u00b2)")
    ax.set_ylabel("frequency (m$^{-2}$)")
    ax.grid(True, which="major", ls="-", alpha=0.3, zorder=0)
    ax.set_ylim(top=1e-1)
    ax.set_xlim([1e0, 1e5])


# ---------------------------------------------------------------------------
# Panel (b): Mapping deficit detail
# ---------------------------------------------------------------------------

def _compute_deficit(curves: FADCurves, cfg):
    """Estimate theoretically missed landslides between rollover and AREA_MIN."""
    bin_widths = np.diff(curves.bins_event)
    n_total = curves.n_event
    rollover = curves.max_area_event

    missed_bins = []
    for i, (bc, dens) in enumerate(zip(curves.bin_centers_event, curves.event_density)):
        if bc < rollover or bc > cfg.AREA_MIN:
            continue
        # Mapped count from raw data
        lo, hi = curves.bins_event[i], curves.bins_event[i + 1]
        mapped = int(np.sum((curves.event_areas >= lo) & (curves.event_areas < hi)))
        # Theoretical count from power law
        theo_dens = 10**curves.powerlaw_intercept * bc**curves.powerlaw_slope
        theo_count = theo_dens * n_total * bin_widths[i]
        missed = max(0.0, theo_count - mapped)
        missed_bins.append({
            "center": bc, "mapped_density": dens,
            "theoretical_density": theo_dens, "missed_count": missed,
        })

    total_missed = sum(b["missed_count"] for b in missed_bins)
    return missed_bins, total_missed


def _plot_panel_b(ax, curves: FADCurves, cfg):
    """Mapping deficit: extrapolated power law vs mapped data."""
    missed_bins, total_missed = _compute_deficit(curves, cfg)
    rollover = curves.max_area_event

    # Mapped data scatter (event only)
    ax.scatter(curves.bin_centers_event, curves.event_density,
               color="#762A83", marker="o", zorder=3, s=25, alpha=0.5,
               label="Mapped Data")

    # Solid power-law line in fit range
    areas_pl = np.logspace(np.log10(cfg.AREA_MIN), np.log10(cfg.AREA_MAX), 100)
    freq_pl = 10**curves.powerlaw_intercept * areas_pl**curves.powerlaw_slope
    ax.plot(areas_pl, freq_pl, color="black", linestyle="-", linewidth=2, zorder=5)

    # Dashed extrapolation below AREA_MIN
    extrap = np.logspace(np.log10(rollover), np.log10(cfg.AREA_MIN), 100)
    extrap_freq = 10**curves.powerlaw_intercept * extrap**curves.powerlaw_slope
    ax.plot(extrap, extrap_freq, color="black", linestyle="--", linewidth=2, zorder=5)

    # Shaded deficit
    if missed_bins:
        centers = [b["center"] for b in missed_bins]
        mapped = [b["mapped_density"] for b in missed_bins]
        theo = [b["theoretical_density"] for b in missed_bins]
        ax.fill_between(centers, mapped, theo,
                        color="red", alpha=0.1, zorder=1,
                        label=f"Mapping Deficit (n={int(round(total_missed))})")
        # Count labels
        for b in missed_bins:
            if b["missed_count"] > 10:
                ax.text(b["center"] * 1.4, b["theoretical_density"],
                        f"+{int(round(b['missed_count']))}",
                        color="#762A83", fontsize=8, ha="center",
                        fontweight="bold", alpha=0.7)

    ax.axvline(x=rollover, color="#762A83", linestyle="--", zorder=2, alpha=0.5)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(True, which="major", ls="-", alpha=0.3, zorder=0)
    ax.set_ylim([1e-3, 1e-1])
    ax.set_xlim([1e1, 1e3])

    return total_missed


# ---------------------------------------------------------------------------
# Combined figure
# ---------------------------------------------------------------------------

def plot_fad(curves: FADCurves, cfg) -> tuple[plt.Figure, float]:
    """Two-panel FAD figure: (a) full comparison, (b) mapping deficit."""
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 8

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(165 / 25.4, 60 / 25.4))

    _plot_panel_a(ax_a, curves, cfg)
    total_missed = _plot_panel_b(ax_b, curves, cfg)

    fig.tight_layout()
    return fig, total_missed


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    import sys as _sys

    root = Path(__file__).resolve().parents[1]
    if str(root) not in _sys.path:
        _sys.path.append(str(root))

    import config as cfg

    print("Computing frequency-area curves...")
    curves = compute_fad_curves(cfg)

    print(f"\n2022: shape={curves.event_params[0]:.2f}, "
          f"loc={curves.event_params[1]:.2f}, scale={curves.event_params[2]:.2f}")
    print(f"Historical: shape={curves.historical_params[0]:.2f}, "
          f"loc={curves.historical_params[1]:.2f}, scale={curves.historical_params[2]:.2f}")
    print(f"\n2022 count: {curves.n_event}")
    print(f"Historical count: {curves.n_hist}")
    print(f"\n2022 max frequency area (rollover): {curves.max_area_event:.2f} m\u00b2")
    print(f"Historical max frequency area: {curves.max_area_hist:.2f} m\u00b2")
    print(f"\nPower Law: \u03b2 = {curves.powerlaw_beta:.2f}")

    fig, total_missed = plot_fad(curves, cfg)
    print(f"Estimated missed landslides (rollover to {cfg.AREA_MIN} m\u00b2): "
          f"~{int(round(total_missed))}")

    output_pdf = cfg.figure_path / "fad_2022_vs_historical.pdf"
    fig.savefig(output_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"\nFAD figure saved to: {output_pdf}")


if __name__ == "__main__":
    main()
