"""Area-volume scaling, cumulative erosion analysis, and DoD vs M3C2 comparison.

Reads precomputed erosion volumes from CSV (landslide_erosion_volumes.csv) and
the event inventory shapefile for planform areas. Produces three publication
figures:
  1. Area-volume power-law scaling (with Larsen et al. 2010 reference)
  2. Cumulative erosion semilog fit
  3. DoD vs M3C2 erosion volume scatter comparison
"""

from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import ticker as mticker
from scipy.stats import linregress, pearsonr, spearmanr
from scipy import stats
import warnings
import sys

# Project config
sys.path.append(str(Path(__file__).resolve().parent.parent))
import config as cfg

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_volumes():
    """Load areas, erosion volumes, and slope from the precomputed CSV."""
    return pd.read_csv(cfg.erosion_volumes_csv)


# ---------------------------------------------------------------------------
# Area-volume scaling
# ---------------------------------------------------------------------------

def _filter_dod_volumes(inv):
    """Return filtered DataFrame with valid DoD erosion volumes for scaling."""
    df = inv[
        (inv["dod_erosion_volume_m3"].notnull())
        & (inv["dod_erosion_volume_m3"] < -1)
        & (inv["dod_erosion_volume_m3"] > -1e4)
        & (np.isfinite(inv["dod_erosion_volume_m3"]))
        & (inv["area_m2"] > 0)
    ].copy()
    df["erosion_abs"] = -1 * df["dod_erosion_volume_m3"]
    return df


def fit_area_volume(df):
    """Fit V = alpha * A^gamma in log-log space. Returns dict of results."""
    area = df["area_m2"].values
    volume = df["erosion_abs"].values

    log_a = np.log10(area)
    log_v = np.log10(volume)
    slope, intercept, r_value, _, _ = linregress(log_a, log_v)

    return {
        "alpha": 10**intercept,
        "gamma": slope,
        "r_squared": r_value**2,
        "n": len(area),
        "area": area,
        "volume": volume,
    }


def plot_area_volume(fit, output_dir):
    """Area-volume scaling figure with Larsen et al. (2010) reference."""
    _set_rc()
    fig, ax = plt.subplots(figsize=(80 / 25.4, 60 / 25.4))

    area, volume = fit["area"], fit["volume"]
    alpha, gamma, r2 = fit["alpha"], fit["gamma"], fit["r_squared"]

    # Scatter
    sns.scatterplot(x=area, y=volume, s=10, color="#762A83",
                    edgecolor="black", alpha=0.4, ax=ax, zorder=1)

    # Best-fit line
    x_vals = np.logspace(np.log10(area.min()), np.log10(area.max()), 100)
    ax.plot(x_vals, alpha * x_vals**gamma, color="black", linewidth=2,
            label=f"This study (r\u00b2 = {r2:.3f})", zorder=3)

    # Larsen et al. (2010) reference
    larsen_alpha = 10**-0.44
    ax.fill_between(x_vals,
                    larsen_alpha * x_vals**1.1,
                    larsen_alpha * x_vals**1.4,
                    color="gray", alpha=0.4, zorder=2)
    ax.plot(x_vals, larsen_alpha * x_vals**1.145,
            linestyle="--", color=(80/255, 80/255, 80/255),
            linewidth=2, zorder=2.5)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, which="major", ls="-", zorder=0)
    ax.set_axisbelow(True)
    ax.set_xlabel("area (m\u00b2)")
    ax.set_ylabel("volume (m\u00b3)")
    fig.tight_layout()

    path = output_dir / "area_volume_scaling_dod.pdf"
    fig.savefig(path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved: {path}")
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Cumulative erosion semilog fit
# ---------------------------------------------------------------------------

def compute_cumulative_erosion(df):
    """Sort by area, compute cumulative erosion percentage."""
    df_sorted = df.sort_values("area_m2").reset_index(drop=True)
    df_sorted["cumulative_erosion"] = df_sorted["erosion_abs"].cumsum()
    total = df_sorted["cumulative_erosion"].iloc[-1]
    df_sorted["cumulative_erosion_pct"] = df_sorted["cumulative_erosion"] / total * 100
    return df_sorted, total


def fit_semilog_erosion(areas, cum_pct, area_min=None, area_max=None):
    """Semilog regression: cum_pct = m * log10(area) + b."""
    area_min = area_min or cfg.EROSION_AREA_MIN
    area_max = area_max or cfg.EROSION_AREA_MAX

    mask = (areas >= area_min) & (areas <= area_max) & (cum_pct >= 0)
    a_fit, e_fit = areas[mask], cum_pct[mask]
    if len(a_fit) < 3:
        return None

    m, b, r, _, se = stats.linregress(np.log10(a_fit), e_fit)
    return {"m": m, "b": b, "r_squared": r**2, "std_err": se,
            "areas_fit": a_fit, "erosion_pct_fit": e_fit}


def bootstrap_semilog_erosion(areas, cum_pct, area_min=None, area_max=None,
                              n_bootstrap=1000):
    """Bootstrap CIs for semilog cumulative erosion fit."""
    area_min = area_min or cfg.EROSION_AREA_MIN
    area_max = area_max or cfg.EROSION_AREA_MAX

    mask = (areas >= area_min) & (areas <= area_max) & (cum_pct >= 0)
    a_fit, e_fit = areas[mask], cum_pct[mask]
    n = len(a_fit)
    if n < 3:
        return None

    slopes, intercepts, r2s = [], [], []
    for _ in range(n_bootstrap):
        idx = np.random.choice(n, size=n, replace=True)
        try:
            m, b, r, _, _ = stats.linregress(np.log10(a_fit[idx]), e_fit[idx])
            slopes.append(m)
            intercepts.append(b)
            r2s.append(r**2)
        except Exception:
            continue

    return {
        "slopes": np.array(slopes),
        "intercepts": np.array(intercepts),
        "r_squared": np.array(r2s),
    }


def plot_cumulative_erosion(df_sorted, total_erosion, semilog_fit, boot,
                            output_dir):
    """Cumulative erosion figure with semilog fit and secondary y-axis."""
    _set_rc()
    fig, ax = plt.subplots(figsize=(82 / 25.4, 63 / 25.4))

    areas = df_sorted["area_m2"].values
    cum_pct = df_sorted["cumulative_erosion_pct"].values

    # Shaded fit range
    ax.axvspan(cfg.EROSION_AREA_MIN, cfg.EROSION_AREA_MAX,
               color="#762A83", alpha=0.15, zorder=1,
               label=f"Semilog fit range ({cfg.EROSION_AREA_MIN}\u2013{cfg.EROSION_AREA_MAX} m\u00b2)")

    # Cumulative curve
    ax.plot(areas, cum_pct, color="#762A83", linewidth=1.5,
            label="Cumulative Erosion", zorder=3)

    # Semilog fit line
    if semilog_fit is not None:
        x_line = np.logspace(np.log10(cfg.EROSION_AREA_MIN),
                             np.log10(cfg.EROSION_AREA_MAX), 100)
        m_mean = np.mean(boot["slopes"]) if boot else semilog_fit["m"]
        y_line = semilog_fit["m"] * np.log10(x_line) + semilog_fit["b"]
        ax.plot(x_line, y_line, color="black", linewidth=1.5, linestyle="--",
                label=f"Semilog Fit: m = {m_mean:.1f} %/log(m\u00b2)", zorder=4)

    ax.set_xscale("log")
    ax.set_xlabel("landslide area (m\u00b2)")
    ax.set_ylabel("cumulative erosion (%)")
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 100)

    # Secondary y-axis: raw volume
    secax = ax.secondary_yaxis(
        "right",
        functions=(lambda p: total_erosion * p / 100,
                   lambda v: (v / total_erosion) * 100),
    )
    secax.set_ylabel("m\u00b3")
    secax.yaxis.set_major_locator(mticker.MultipleLocator(1e5))
    sf = mticker.ScalarFormatter(useMathText=True)
    sf.set_scientific(True)
    sf.set_powerlimits((5, 5))
    secax.yaxis.set_major_formatter(sf)

    fig.tight_layout()
    path = output_dir / "cumulative_erosion_semilog_fit.pdf"
    fig.savefig(path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved: {path}")
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# DoD vs M3C2 scatter comparison
# ---------------------------------------------------------------------------

def plot_dod_vs_m3c2(inv, output_dir):
    """Scatter plot comparing DoD and M3C2 slope-normal erosion volumes."""
    # Filter valid pairs
    df = inv[
        inv["dod_erosion_volume_m3"].notnull()
        & inv["m3c2_sn_erosion_volume_m3"].notnull()
        & (inv["dod_erosion_volume_m3"] != 0)
        & (inv["m3c2_sn_erosion_volume_m3"] != 0)
    ].copy()
    df = df[
        np.isfinite(df["dod_erosion_volume_m3"])
        & np.isfinite(df["m3c2_sn_erosion_volume_m3"])
    ]

    df["dod_abs"] = np.abs(df["dod_erosion_volume_m3"])
    df["m3c2_abs"] = np.abs(df["m3c2_sn_erosion_volume_m3"])
    df["ratio"] = df["dod_abs"] / df["m3c2_abs"]

    # Remove extreme outliers for correlation
    cdf = df[(df["ratio"] < 10) & (df["ratio"] > 0.1)].copy()
    print(f"Valid comparison pairs: {len(cdf)} landslides")

    # Correlations
    r_pearson, p_pearson = pearsonr(np.log10(cdf["dod_abs"]),
                                    np.log10(cdf["m3c2_abs"]))
    r_spearman, p_spearman = spearmanr(cdf["dod_abs"], cdf["m3c2_abs"])

    print(f"Pearson  r = {r_pearson:.3f}, p < 0.001")
    print(f"Spearman \u03c1 = {r_spearman:.3f}, p < 0.001")

    # Plot
    _set_rc()
    fig, ax = plt.subplots(figsize=(80 / 25.4, 60 / 25.4))

    sc = ax.scatter(cdf["m3c2_abs"], cdf["dod_abs"],
                     s=10, c=cdf["mean_slope_deg"], cmap="viridis",
                     alpha=0.4, edgecolors="black", linewidth=0.3)

    # 1:1 identity line
    lo = min(cdf["m3c2_abs"].min(), cdf["dod_abs"].min())
    hi = max(cdf["m3c2_abs"].max(), cdf["dod_abs"].max())
    ax.plot([lo, hi], [lo, hi], linestyle="--", color=(0, 0, 0, 0.7),
            linewidth=1)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("M3C2 slope-normal erosion volume (m\u00b3)")
    ax.set_ylabel("DoD erosion volume (m\u00b3)")
    ax.set_xlim([1e0, 1e4])
    ax.set_ylim([1e0, 1e4])
    ax.grid(True, alpha=0.3)

    # Colorbar
    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label("mean slope (\u00b0)")

    # Annotation box
    fs = mpl.rcParams["font.size"]
    table_text = (
        f"Pearson r: {r_pearson:.3f}\n"
        f"p-value: <0.001\n\n"
        f"Spearman \u03c1: {r_spearman:.3f}\n"
        f"p-value: <0.001"
    )
    ax.text(0.03, 0.97, table_text, transform=ax.transAxes,
            fontsize=fs * 0.8, verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.2,rounding_size=0.05",
                      facecolor="white", edgecolor="black", alpha=0.9))

    fig.tight_layout()
    path = output_dir / "dod_vs_m3c2_comparison.pdf"
    fig.savefig(path, format="pdf", dpi=300, bbox_inches="tight")
    print(f"Saved: {path}")
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _set_rc():
    """Set matplotlib rcParams for publication figures."""
    mpl.rcParams.update(mpl.rcParamsDefault)
    mpl.rcParams["font.family"] = "Arial"
    mpl.rcParams["font.size"] = 8


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("Loading inventory and precomputed volumes...")
    inv = load_volumes()
    df_dod = _filter_dod_volumes(inv)

    # 1. Area-volume scaling
    print("\n--- Area-Volume Scaling ---")
    av_fit = fit_area_volume(df_dod)
    print(f"V = {av_fit['alpha']:.3f} * A^{av_fit['gamma']:.3f}  "
          f"(r\u00b2 = {av_fit['r_squared']:.3f}, n = {av_fit['n']})")
    print(f"log\u2081\u2080(\u03b1) = {np.log10(av_fit['alpha']):.4f}")
    plot_area_volume(av_fit, cfg.figure_path)

    # 2. Cumulative erosion
    print("\n--- Cumulative Erosion ---")
    df_sorted, total = compute_cumulative_erosion(df_dod)
    areas_e = df_sorted["area_m2"].values
    cum_pct = df_sorted["cumulative_erosion_pct"].values
    print(f"Total erosion volume: {total:,.0f} m\u00b3")

    semilog = fit_semilog_erosion(areas_e, cum_pct)
    boot = None
    if semilog:
        print(f"Semilog fit: Ec = {semilog['m']:.1f} * log10(A) + ({semilog['b']:.1f})  "
              f"(r\u00b2 = {semilog['r_squared']:.4f})")
        boot = bootstrap_semilog_erosion(areas_e, cum_pct)
        if boot and len(boot["slopes"]) > 100:
            ci = np.percentile(boot["slopes"], [2.5, 97.5])
            print(f"Bootstrap m: {np.mean(boot['slopes']):.1f} "
                  f"[{ci[0]:.1f}, {ci[1]:.1f}] (95% CI)")

    plot_cumulative_erosion(df_sorted, total, semilog, boot, cfg.figure_path)

    # 3. DoD vs M3C2 comparison
    print("\n--- DoD vs M3C2 Comparison ---")
    plot_dod_vs_m3c2(inv, cfg.figure_path)

    print("\nDone.")


if __name__ == "__main__":
    main()
