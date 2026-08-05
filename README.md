# Lidar change detection applied to a convective storm-triggered landslide event in Kentucky, USA highlights limitations of conventional inventories

**Supporting code for the manuscript submitted to _Landslides_ by Scheip et al.**

This repository contains the code to reproduce figures and tables for the above-mentioned manuscript to support open research standards.

## Overview

This codebase generates figures and tables used in the manuscript from raw raster or vector data inputs. The analysis includes:

- Landslide density and area fraction maps (KDE surfaces)
- Aspect frequency analysis and rose diagrams
- Frequency-area distributions with power law fitting
- Area-volume scaling relationships from lidar change detection
- Reactivation distance analysis
- Cumulative erosion analysis
- Precipitation analysis and rolling sum plots

## Requirements

### Software Dependencies

This project uses [uv](https://docs.astral.sh/uv/) for environment and dependency
management. Install uv (see the [install docs](https://docs.astral.sh/uv/getting-started/installation/)
for other options):

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Then create the environment from the repository root:

```bash
uv sync
```

This reads `pyproject.toml` and `uv.lock`, downloads Python 3.11 if needed, and
installs the exact pinned dependency versions into a local `.venv/`. No separate
activation step is required — prefix commands with `uv run` as shown below.

### Data Requirements

Download the following additional datasets and place them in the `data/` subdirectory:

| Filename | Description | Source |
|----------|-------------|--------|
| `event-inventory.shp` | 2022 event landslide inventory (vector) | [OSF Project](https://osf.io/tuvqw/overview?view_only=79501ea561984e198188d3cf4dd7f1c7) |
| `historical-inventory.shp` | Historical landslide inventory (vector) | [OSF Project](https://osf.io/tuvqw/overview?view_only=79501ea561984e198188d3cf4dd7f1c7) |
| `mapping-area.shp` | Study/mapping area boundary polygon | [Zenodo Record](https://zenodo.org/uploads/16813914) |
| `aspect.tif` | DEM-derived aspect raster | [Kentucky from Above](https://kyfromabove.ky.gov) |
| `m3c2_lcd.tif` | Lidar change detection (LCD) raster (M3C2, Surface-Normal) | [Zenodo Record](https://zenodo.org/uploads/16813914) |
| `dod_aligned.tif` | Aligned DEM-of-difference raster | [Zenodo Record](https://zenodo.org/uploads/16813914) |
| `landslide_erosion_volumes.csv` | Precomputed M3C2 and DoD erosion volumes with mean slope per landslide | Included in repository |
| `eky-gauge-data.csv` | Precipitation data for four gauges near mapping area | Included in repository |

## Configuration

Edit `config.py` to adjust paths if your data is stored elsewhere. By default, the configuration assumes:
- Data files are in `./data/`
- Outputs will be written to `./outputs/`
- Figures will be saved to `./outputs/figures/`

## Usage

Run the analysis scripts in the following order from the repository root:

```bash
# 1. Density maps (event and historical)
uv run python src/density_maps.py

# 2. Aspect frequency analysis
uv run python src/aspect_frequency.py

# 3. Frequency-area distributions
uv run python src/fad.py

# 4. Area-volume scaling and cumulative erosion
uv run python src/area_volume.py

# 5. Reactivation distance analysis
uv run python src/reactivation_distance.py

# 6. Precipitation analysis and plotting
uv run python src/precip_plotting.py
```

`uv run` syncs the environment before executing, so `uv sync` is optional if you
go straight to running a script. All scripts resolve their input and output paths
through `config.py`, so they can be invoked from any working directory.

### Expected Outputs

After running all scripts, the `outputs/figures/` directory will contain:

**Density Maps:**
- `density-maps-event.pdf` - Event landslide density and area fraction
- `density-maps-historical.pdf` - Historical landslide density and area fraction

**Aspect Analysis:**
- `aspect_rose_bars.pdf` - Polar rose diagram of aspect distributions
- `aspect_bars_select.pdf` - N-E coverage frequency ratios
- `aspect_bars_d8_historical_event.pdf` - D8 frequency ratio comparison
- `aspect_frequency_*.csv` - Frequency tables

**Frequency-Area Distribution:**
- `fad_2022_vs_historical.pdf` - FAD comparison with inverse gamma fits and mapping deficit detail

**Area-Volume Analysis:**
- `area_volume_scaling_dod.pdf` - Area vs. volume scatter with power law fit
- `cumulative_erosion_semilog_fit.pdf` - Cumulative erosion curve
- `dod_vs_m3c2_comparison.pdf` - DoD vs M3C2 erosion volume scatter comparison

**Reactivation Analysis:**
- `cumulative_distance_curve_with_size.pdf` - Reactivation distance curve
- `reactivation_cumulative.csv` - Cumulative reactivation statistics

**Precipitation Analysis:**
- `Month_AverageFullSpan.png` - 30-day rolling sum precipitation (full time span)
- `FourDay_GaugeInset.png` - 4-day rolling sum precipitation (June-August 2022)
- `precip_daily_percentiles.csv` - 1-day total distribution per gauge (90th/95th/99th percentiles over all days and over wet days, plus record daily total)
- `precip_daily_top_events.csv` - 20 wettest days in the record ranked by highest single-gauge total
- `precip_daily_exceedances.csv` - counts of days meeting 75/100/150/200 mm daily thresholds, per gauge, pooled, and network-wide

## Analysis Parameters

Key parameters are defined in `config.py`:

- `target_crs`: Coordinate reference system (EPSG:32617 - UTM Zone 17N)
- `kde_method`: Kernel density estimation bandwidth method ("scott")
- `n_bins`: Number of bins for frequency-area distributions (50)
- `AREA_MIN`, `AREA_MAX`: Power law fitting range (150-10000 m²)
- `EROSION_AREA_MIN`, `EROSION_AREA_MAX`: Erosion analysis range (150-4000 m²)

## Code Structure

```
eky-landslides-2022/
├── config.py              # Configuration and paths
├── pyproject.toml         # Project metadata and dependencies (uv)
├── uv.lock                # Pinned dependency versions (uv)
├── .python-version        # Python version used by uv (3.11)
├── README.md             # This file
├── data/                 # Place raw data files here (not tracked)
├── outputs/              # Generated outputs (not tracked)
│   └── figures/          # Generated figures
└── src/                  # Source code modules
    ├── utils.py          # Utility functions for loading data
    ├── density_maps.py   # KDE density surface generation
    ├── aspect_frequency.py # Aspect analysis and frequency ratios
    ├── fad.py            # Frequency-area distributions
    ├── area_volume.py    # Volume calculations and scaling
    ├── reactivation_distance.py # Reactivation analysis
    └── precip_plotting.py # Precipitation analysis and plotting
```

## Notes

- All scripts can be run independently after data is downloaded
- Progress bars will display for long-running computations (e.g., aspect analysis)
- Figures are saved as PDFs
- CSV tables are saved alongside figures for reference

## Citation

If you use this code, please cite:

```
Scheip, C., Crawford, M., Koch, H., Bibbins, E., (2025). RReconsidering the magnitude of convective storm-triggered landslide events in the Appalachian Plateau, USA. Journal of Geophysical Research: Earth Surface. [Submitted]
```

## License

This project is licensed under the Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0) - see the [LICENSE](LICENSE) file for details.

[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)

**You are free to:**
- **Share** — copy and redistribute the material in any medium or format
- **Adapt** — remix, transform, and build upon the material

**Under the following terms:**
- **Attribution** — You must give appropriate credit, provide a link to the license, and indicate if changes were made
- **NonCommercial** — You may not use the material for commercial purposes

For commercial use, please contact cscheip@bgcengineering.com

## Contact

For questions or issues, please contact cscheip@bgcengineering.com
