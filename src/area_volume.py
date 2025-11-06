from pathlib import Path
import geopandas as gpd
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import rasterio as rio
import rasterio.features
import rioxarray as rxr
import xarray as xr
from tqdm import tqdm
from scipy.stats import linregress
import warnings

# --
import sys
# Add the parent path of this script to sys.path so imports from the project root work
sys.path.append(str(Path(__file__).resolve().parent.parent))
import config as cfg

# ----------------------------
# I/O Setup
# ----------------------------
# Project folder
project_path = cfg.project_path

# Input files
inv_2022_path = cfg.inv_2022_path
inv_historical_path = cfg.inv_historical_path
aspect_file = cfg.aspect_file
dod_file = cfg.dod_aligned_file

# Output directory
figure_path = cfg.figure_path

# ----------------------------
# Parameters
# ----------------------------
target_crs = cfg.target_crs

# ----------------------------
# Read the inventories
# ----------------------------

# Read the inventory
inv_2022 = gpd.read_file(inv_2022_path).to_crs(target_crs)
inv_historical = gpd.read_file(inv_historical_path).to_crs(target_crs)

# # add a guid to the dataset that maintains after overlay
# inv_2022["guid"] = inv_2022.index
# inv_historical = inv_historical.reset_index(drop=True)
# inv_historical["guid"] = inv_historical.index

# for both, compute the area in square meters and save to field `area_m2`
inv_2022["area_m2"] = inv_2022.geometry.area
inv_historical["area_m2"] = inv_historical.geometry.area

# =============================================================================
# DOD VOLUME ANALYSIS
# =============================================================================
print("Setting up DOD volume computation...")
dod_erosion_volumes = []
dod_deposition_volumes = []
dod_net_volumes = []

with rio.open(dod_file) as dod:
    print(f"DOD file info:")
    print(f"  Shape: {dod.width} x {dod.height}")
    print(f"  Resolution: {dod.res[0]:.2f} x {dod.res[1]:.2f} m")
    print(f"  CRS: {dod.crs}")
    print(f"  Bounds: {dod.bounds}")
    
    pixel_area = dod.res[0] * dod.res[1]
    
    for landslide in tqdm(inv_2022.itertuples(), total=len(inv_2022), desc="Computing DOD volumes"):
        try:
            # Get bounding window for landslide
            window = rio.features.geometry_window(dod, [landslide.geometry], pad_x=100, pad_y=100)
            dod_data = dod.read(1, window=window)
            window_transform = dod.window_transform(window)
            
            # Create mask for landslide polygon
            mask_array = rio.features.geometry_mask(
                [landslide.geometry],
                transform=window_transform,
                invert=True,
                out_shape=dod_data.shape
            )
            
            # Extract DOD values within landslide
            dod_masked = dod_data[mask_array]
            
            # Remove NaN and infinite values
            dod_valid = dod_masked[np.isfinite(dod_masked)]
            
            if len(dod_valid) > 0:
                # Calculate volumes (positive = deposition, negative = erosion)
                # DOD represents elevation change, so negative = loss/erosion
                erosion_vol = np.sum(dod_valid[dod_valid < 0]) * pixel_area
                deposition_vol = np.sum(dod_valid[dod_valid > 0]) * pixel_area
                net_vol = np.sum(dod_valid) * pixel_area
            else:
                erosion_vol = np.nan
                deposition_vol = np.nan
                net_vol = np.nan
                
        except Exception as e:
            print(f"Error processing landslide {landslide.Index}: {e}")
            erosion_vol = np.nan
            deposition_vol = np.nan
            net_vol = np.nan
        
        dod_erosion_volumes.append(erosion_vol)
        dod_deposition_volumes.append(deposition_vol)
        dod_net_volumes.append(net_vol)

# Add DOD volumes to GeoDataFrame
inv_2022["dod_erosion_volume_m3"] = dod_erosion_volumes
inv_2022["dod_deposition_volume_m3"] = dod_deposition_volumes
inv_2022["dod_net_volume_m3"] = dod_net_volumes

print(f"✅ DOD volume computation complete!")
print(f"DOD volume statistics:")
print(f"  Erosion (negative): {len([v for v in dod_erosion_volumes if v < 0])} landslides")
print(f"  Deposition (positive): {len([v for v in dod_deposition_volumes if v > 0])} landslides")
print(f"  Valid volumes: {len([v for v in dod_net_volumes if not np.isnan(v)])} landslides")


# =============================================================================
# CREATE FIGURE
# =============================================================================
# Set matplotlib font to Arial and font size to 8 for all text elements
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)  # Reset to defaults
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.size'] = 8

# Use DoD volumes for area-volume scaling
# Filter: only landslides with valid DoD volumes (not null, not nan, not inf, and >0)
df_dod = inv_2022[
    (inv_2022['dod_erosion_volume_m3'].notnull()) &
    (inv_2022['dod_erosion_volume_m3'] < -1) &
    (inv_2022['dod_erosion_volume_m3'] >-1e4) &
    (np.isfinite(inv_2022['dod_erosion_volume_m3'])) &
    (inv_2022['area_m2'] > 0)
].copy()

area = df_dod['area_m2']
volume = -1*df_dod['dod_erosion_volume_m3']

# Fit in log-log space
log_area = np.log10(area)
log_volume = np.log10(volume)
slope, intercept, r_value, _, _ = linregress(log_area, log_volume)
r_squared = r_value ** 2
alpha = 10 ** intercept
gamma = slope

# Generate best-fit line in linear space
x_vals = np.logspace(np.log10(area.min()), np.log10(area.max()), 100)
y_fit = alpha * x_vals ** gamma

# Larsen parameters from literature
larsen_alpha = 10 ** -0.44  # From log₁₀(α) = -0.44
larsen_gamma = 1.145

# Larsen shallow range bounds (1.1 to 1.4 gamma range)
y_gamma_1_1 = larsen_alpha * x_vals ** 1.1
y_gamma_1_4 = larsen_alpha * x_vals ** 1.4

# Larsen specific reference line
y_larsen = larsen_alpha * x_vals ** larsen_gamma

# Plot
fig, ax = plt.subplots(figsize=(80/25.4, 60/25.4))

# Scatter with Purple (#762A83)
sns.scatterplot(
    x=area, y=volume, 
    s=10, color='#762A83', edgecolor='black', alpha=0.4, ax=ax, zorder=1
)

# Best-fit line 
line, = ax.plot(
    x_vals, y_fit, color='black', linewidth=2, 
    label=f"This study (r² = {r_squared:.3f})",
    zorder=3
)
line.set_gid('best_fit')  # optional, for clarity
from matplotlib.lines import Line2D

# Shaded Larsen shallow landslide range (updated to 1.1-1.4 gamma range)
ax.fill_between(
    x_vals, y_gamma_1_1, y_gamma_1_4, 
    color='gray', alpha=0.4, 
    zorder=2
)

# Larsen reference line (dark gray dashed)
ax.plot(x_vals, y_larsen, linestyle='--', color=(80/255, 80/255, 80/255), linewidth=2, zorder=2.5)

# Log scale axes
ax.set_xscale('log')
ax.set_yscale('log')

# add grid
ax.grid(True, which="major", ls="-", zorder=0)
ax.set_axisbelow(True)

# Labels and title
ax.set_xlabel("area (m²)")
ax.set_ylabel("volume (m³)")
fig.tight_layout()
plt.show()

# Summary
print(f"Power-law model: V = α · A^γ")
print(f"Estimated parameters: α = {alpha:.4e} (log₁₀(α) = {np.log10(alpha):.4f}), γ = {gamma:.4f}")
print(f"R² = {r_squared:.4f}")
print(f"n = {len(area)}")

# Export to disk
pdf_path = figure_path / "area_volume_scaling_dod.pdf"
fig.savefig(pdf_path, format="pdf")
print(f"Figure saved to: {pdf_path}")

# =============================================================================
# CREATE TABLE
# =============================================================================

# Import great_tables
import great_tables as gt
from great_tables import GT, html, md


# Create teh Larsen vars
larsen_log_alpha = -0.44  # ± 0.02 
larsen_gamma = 1.145  # ± 0.008
larsen_r_squared = 0.9
larsen_n = 1741

# Recreate the comparison data
table_data_gt = {
    '': ['Larsen et al. (2010)', 'This study'],
    'log₁₀(α)': [larsen_log_alpha, np.log10(alpha)],
    'γ': [larsen_gamma, gamma],
    'r²': [larsen_r_squared, r_squared],
    'n': [larsen_n, len(area)]
}

# Create DataFrame
df_gt = pd.DataFrame(table_data_gt)

# Create the great_tables object (compact version)
gt_table = (
    GT(df_gt)
    .fmt_number(
        columns=['log₁₀(α)'],
        decimals=3
    )
    .fmt_number(
        columns=['γ'],
        decimals=3
    )
    .fmt_number(
        columns=['r²'],
        decimals=3
    )
    .fmt_integer(
        columns=['n']
    )
    .tab_style(
        style=[
            gt.style.text(weight="bold"),
            gt.style.text(align="left")
        ],
        locations=gt.loc.body(columns="Study")
    )
    .tab_style(
        style=gt.style.text(align="center"),
        locations=gt.loc.body(columns=['log₁₀(α)', 'γ', 'r²', 'n'])
    )
    .tab_style(
        style=gt.style.fill(color="#f9f9f9"),
        locations=gt.loc.body(rows=[0])  # Larsen et al. row
    )
    .tab_style(
        style=gt.style.fill(color="#ffffff"),
        locations=gt.loc.body(rows=[1])  # This study row
    )
    .opt_table_font(
        font="Arial"
    )
    .tab_style(
        style=gt.style.text(size="7px"),
        locations=gt.loc.body()
    )
    .tab_options(
        # table_width="auto",  # Auto-size to content
        table_width="52mm",
        column_labels_border_top_width="2px",
        column_labels_border_bottom_width="2px",
        table_border_top_width="0px",
        table_border_bottom_width="2px",
        data_row_padding="4px",  # Compact row padding
        column_labels_padding="4px",  # Compact column padding
        column_labels_font_size="8px",
        column_labels_font_weight="bold",
        column_labels_background_color="#f0f0f0"
    )
)

# # Display the table
# gt_table.show()

pdf_gt_file = figure_path / "power_law_comparison_great_tables_dod.pdf"

gt_table.save(
    file=pdf_gt_file,
    scale=8.0,  # High resolution for publication quality
    web_driver='chrome'  # Use Chrome for rendering
)


# =============================================================================
# SEMILOG AREA VS CUMULATIVE EROSION FITTING AND PLOT
# =============================================================================

from scipy import stats

def fit_semilog_to_cumulative_erosion(areas, cumulative_erosion_pct, area_min=None, area_max=None):
    """
    Fit semilog relationship to area-cumulative erosion percentage data:
    cumulative_erosion_pct = m × log₁₀(area) + b
    
    This is the correct semilog relationship where:
    - x-axis: log₁₀(area) 
    - y-axis: cumulative_erosion_pct (0-100%)
    
    Args:
        areas: Landslide areas (m²)
        cumulative_erosion_pct: Cumulative erosion percentage (0-100%)
        area_min, area_max: Fitting range bounds
    
    Returns:
        slope_m: Logarithmic scaling coefficient (m) [%/log(m²)]
        intercept_b: Intercept (b) [%]
        r_squared: Coefficient of determination
        fit_details: Tuple with fitting details
    """
    # Use global parameters if not specified
    if area_min is None:
        area_min = cfg.EROSION_AREA_MIN
    if area_max is None:
        area_max = cfg.EROSION_AREA_MAX

    # Filter to the fitting range
    mask = (areas >= area_min) & (areas <= area_max) & (cumulative_erosion_pct >= 0)
    areas_fit = areas[mask]
    erosion_pct_fit = cumulative_erosion_pct[mask]

    if len(areas_fit) < 3:
        return None, None, None, None

    # Semilog regression: cumulative_erosion_pct = m × log₁₀(area) + b
    log_areas = np.log10(areas_fit)
    slope_m, intercept_b, r_value, p_value, std_err = stats.linregress(log_areas, erosion_pct_fit)

    r_squared = r_value**2

    return slope_m, intercept_b, r_squared, (slope_m, intercept_b, std_err, areas_fit, erosion_pct_fit)

def bootstrap_semilog_cumulative_erosion(areas, cumulative_erosion_pct, area_min=None, area_max=None, n_bootstrap=1000):
    """
    Bootstrap confidence intervals for area-cumulative erosion percentage semilog fit:
    cumulative_erosion_pct = m × log₁₀(area) + b
    
    Args:
        areas: Landslide areas (m²)
        cumulative_erosion_pct: Cumulative erosion percentage (0-100%)
        area_min, area_max: Fitting range bounds
        n_bootstrap: Number of bootstrap samples
    
    Returns:
        slopes_m: Array of bootstrap slope coefficients [%/log(m²)]
        intercepts_b: Array of bootstrap intercepts [%]
        r_squared_vals: Array of bootstrap R² values
    """
    # Use global parameters if not specified
    if area_min is None:
        area_min = cfg.EROSION_AREA_MIN
    if area_max is None:
        area_max = cfg.EROSION_AREA_MAX

    # Filter to fitting range
    mask = (areas >= area_min) & (areas <= area_max) & (cumulative_erosion_pct >= 0)
    areas_fit = areas[mask]
    erosion_pct_fit = cumulative_erosion_pct[mask]
    n_points = len(areas_fit)

    if n_points < 3:
        return None, None, None

    slopes = []
    intercepts = []
    r_squared_vals = []

    for _ in range(n_bootstrap):
        # Bootstrap sample
        bootstrap_indices = np.random.choice(n_points, size=n_points, replace=True)
        boot_areas = areas_fit[bootstrap_indices]
        boot_erosion_pct = erosion_pct_fit[bootstrap_indices]

        try:
            # Semilog regression: cumulative_erosion_pct = m × log₁₀(area) + b
            log_areas = np.log10(boot_areas)
            slope_m, intercept_b, r_value, _, _ = stats.linregress(log_areas, boot_erosion_pct)
            slopes.append(slope_m)
            intercepts.append(intercept_b)
            r_squared_vals.append(r_value**2)
        except:
            continue

    return np.array(slopes), np.array(intercepts), np.array(r_squared_vals)

print("✅ Semilog fitting functions for area-cumulative erosion data loaded successfully!")


# Use the DoD (Digital Elevation Model of Difference) erosion volumes, filter outliers
df_cdf = inv_2022[
    (inv_2022['dod_erosion_volume_m3'].notnull()) &
    (np.isfinite(inv_2022['dod_erosion_volume_m3'])) &
    (inv_2022['dod_erosion_volume_m3'] < -1) &    # Valid erosion (negative values)
    (inv_2022['dod_erosion_volume_m3'] > -1e4) &  # Remove extreme outliers
    (inv_2022['area_m2'] > 0)
].copy()

# Convert to positive (absolute) erosion volumes as is convention
df_cdf['erosion_volume_abs'] = -1 * df_cdf['dod_erosion_volume_m3']

# Sort by landslide area (size)
df_cdf_sorted = df_cdf.sort_values('area_m2').reset_index(drop=True)

# Calculate cumulative erosion by DoD volume
df_cdf_sorted['cumulative_erosion'] = df_cdf_sorted['erosion_volume_abs'].cumsum()
total_erosion = df_cdf_sorted['cumulative_erosion'].iloc[-1]
df_cdf_sorted['cumulative_erosion_pct'] = (df_cdf_sorted['cumulative_erosion'] / total_erosion) * 100

# Data summary for reference
areas_erosion = df_cdf_sorted['area_m2'].values
erosion_volumes = df_cdf_sorted['erosion_volume_abs'].values
cumulative_erosion = df_cdf_sorted['cumulative_erosion'].values
cumulative_erosion_pct = df_cdf_sorted['cumulative_erosion_pct'].values

print(f"✅ Cumulative erosion data prepared (using DoD volumes):")
print(f"  Total landslides: {len(df_cdf_sorted):,}")
print(f"  Total erosion (DoD): {total_erosion:,.0f} m³")
print(f"  Area range: {areas_erosion.min():.0f} - {areas_erosion.max():,.0f} m²")
print(f"  Cumulative erosion %: {cumulative_erosion_pct.min():.1f}% - {cumulative_erosion_pct.max():.1f}%")

def bootstrap_erosion_semilog(areas, cumulative_erosion, area_min=None, area_max=None, n_bootstrap=1000):
    """
    Bootstrap confidence intervals for area-cumulative erosion semilog (log-linear) fit:
    log10(cumulative_erosion) ~ slope * area + intercept
    """
    # Use global parameters if not specified
    if area_min is None:
        area_min = cfg.EROSION_AREA_MIN
    if area_max is None:
        area_max = cfg.EROSION_AREA_MAX

    # Filter to fitting range
    mask = (areas >= area_min) & (areas <= area_max) & (cumulative_erosion > 0)
    areas_fit = areas[mask]
    erosion_fit = cumulative_erosion[mask]
    n_points = len(areas_fit)

    if n_points < 3:
        return None, None, None

    slopes = []
    intercepts = []
    r_squared_vals = []

    for _ in range(n_bootstrap):
        # Bootstrap sample
        bootstrap_indices = np.random.choice(n_points, size=n_points, replace=True)
        boot_areas = areas_fit[bootstrap_indices]
        boot_erosion = erosion_fit[bootstrap_indices]

        try:
            log_erosion = np.log10(boot_erosion)
            slope, intercept, r_value, _, _ = stats.linregress(boot_areas, log_erosion)
            slopes.append(slope)
            intercepts.append(intercept)
            r_squared_vals.append(r_value**2)
        except:
            continue

    return np.array(slopes), np.array(intercepts), np.array(r_squared_vals)

print("✅ Power law fitting functions for area-erosion data loaded successfully (using DoD volumes)!")

# =============================================================================
# PERFORM SEMILOG ANALYSIS ON AREA-CUMULATIVE EROSION DATA (USING DOD VOLUMES)
# =============================================================================

print("🔍 FITTING SEMILOG RELATIONSHIP TO AREA-CUMULATIVE EROSION DATA (USING DOD VOLUMES)")
print("="*60)
print("   (Confirmed: cumulative erosion is computed from differencing DoD volumes)")

# Fit semilog to cumulative erosion PERCENTAGE data (0-100%)
# This uses DoD-derived cumulative erosion volumes (expressed as percentages)
slope_m, intercept_b, r_squared_semilog, fit_details_semilog = fit_semilog_to_cumulative_erosion(
    areas_erosion, cumulative_erosion_pct  # <-- This is from DoD volumes as percent!
)

if slope_m is not None:
    slope, intercept, std_err, areas_fit, erosion_pct_fit = fit_details_semilog

    print(f"Semilog Fitting Range: {cfg.EROSION_AREA_MIN} - {cfg.EROSION_AREA_MAX} m²")
    print(f"Data points in range: {len(areas_fit)}")
    print(f"")
    print(f"Semilog Fit Results (using DoD volumes):")
    print(f"  Cumulative_Erosion_% = m × log₁₀(Area) + b")
    print(f"  m (logarithmic scaling coefficient) = {slope_m:.1f} %/log(m²)")
    print(f"  b (intercept) = {intercept_b:.1f} %")
    print(f"  R² = {r_squared_semilog:.4f}")
    print(f"  Standard error = {std_err:.2f}")

    # Bootstrap confidence intervals
    print(f"\nBootstrapping confidence intervals (using DoD volumes)...")
    boot_slopes, boot_intercepts, boot_r2_semilog = bootstrap_semilog_cumulative_erosion(
        areas_erosion, cumulative_erosion_pct, n_bootstrap=1000
    )

    if boot_slopes is not None and len(boot_slopes) > 100:
        slope_mean = np.mean(boot_slopes)
        slope_std = np.std(boot_slopes)
        slope_ci = np.percentile(boot_slopes, [2.5, 97.5])

        intercept_mean = np.mean(boot_intercepts)
        intercept_std = np.std(boot_intercepts)
        intercept_ci = np.percentile(boot_intercepts, [2.5, 97.5])
        


# =============================================================================
# PLOT THE SEMILOG FIT
# =============================================================================

# =============================================================================
# PUBLICATION FIGURE: CUMULATIVE EROSION WITH SEMILOG FIT
# =============================================================================

# Set matplotlib style for publication quality
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.size'] = 8

# Create figure
fig, ax = plt.subplots(figsize=(82/25.4, 63/25.4))

# Add shading for the fitting range
ax.axvspan(cfg.EROSION_AREA_MIN, cfg.EROSION_AREA_MAX, color='#762A83', alpha=0.15, zorder=1, 
           label=f"Semilog fit range ({cfg.EROSION_AREA_MIN}-{cfg.EROSION_AREA_MAX} m²)")

# Plot the cumulative erosion curve
ax.plot(areas_erosion, cumulative_erosion_pct, color='#762A83', linewidth=1.5, 
        label='Cumulative Erosion', zorder=3)

# Plot the semilog fit line
if slope_m is not None:
    areas_fit_line = np.logspace(np.log10(cfg.EROSION_AREA_MIN), np.log10(cfg.EROSION_AREA_MAX), 100)
    cumulative_fit_line = slope_m * np.log10(areas_fit_line) + intercept_b
    
    ax.plot(areas_fit_line, cumulative_fit_line, color='black', linewidth=1.5, linestyle='--', zorder=4,
            label=f'Semilog Fit: m = {slope_mean:.1f} %/log(m²)')

# Formatting
ax.set_xscale('log')
ax.set_xlabel('landslide area (m²)')
ax.set_ylabel('cumulative erosion (%)')
ax.grid(True, alpha=0.3)
ax.set_ylim(0, 100)
# Right y-axis: raw cumulative erosion volume (m³)
from matplotlib import ticker as mticker
secax = ax.secondary_yaxis(
    'right',
    functions=(lambda p: total_erosion * p / 100,
               lambda v: (v / total_erosion) * 100)
)
secax.set_ylabel('m³')
# Major ticks every 1e5 m³ → labels 0–5 with ×10^5
secax.yaxis.set_major_locator(mticker.MultipleLocator(1e5))
sf = mticker.ScalarFormatter(useMathText=True)
sf.set_scientific(True)
sf.set_powerlimits((5, 5))
secax.yaxis.set_major_formatter(sf)
# ax.legend(fontsize=7)

plt.tight_layout()
plt.show()

# Export the figure
pdf_path = figure_path / "cumulative_erosion_semilog_fit.pdf"
fig.savefig(pdf_path, format="pdf", dpi=300, bbox_inches='tight')
print(f"\n💾 Figure saved to: {pdf_path}")