# Reconsidering the magnitude of convective storms in triggering landslide events in the Appalachian Plateau, USA

**Supporting code for the manuscript submitted to _Geophysical Research Letters_ by Scheip et al.**

This repository contains the code to reproduce figures and tables for the above-mentioned manuscript to support open research standards.

## Overview

This codebase generates figures and tables used in the manuscript from raw raster or vector data inputs. The analysis includes:

- Landslide density and area fraction maps (KDE surfaces)
- Aspect frequency analysis and rose diagrams
- Frequency-area distributions with power law fitting
- Area-volume scaling relationships from lidar change detection
- Reactivation distance analysis
- Cumulative erosion analysis

## Requirements

### Software Dependencies

Create the conda environment from the provided `environment.yml`:

```bash
conda env create -f environment.yml
conda activate landslide-figures
```

### Data Requirements

Download the following datasets and place them in the `data/` subdirectory:

| Filename | Description | Source |
|----------|-------------|--------|
| `event_inventory_2022.shp` | 2022 event landslide inventory (vector) | [OSF Project](https://osf.io/tuvqw/overview?view_only=79501ea561984e198188d3cf4dd7f1c7) |
| `historical_inventory.shp` | Historical landslide inventory (vector) | [OSF Project](https://osf.io/tuvqw/overview?view_only=79501ea561984e198188d3cf4dd7f1c7) |
| `mapping_area.shp` | Study/mapping area boundary polygon | [Zenodo Record](https://zenodo.org/uploads/16813914) |
| `aspect.tif` | DEM-derived aspect raster | [Kentucky from Above](https://kyfromabove.ky.gov) |
| `M3C2_ICP_2017_vs_2023.tif` | Change-detection raster | [Zenodo Record](https://zenodo.org/uploads/16813914) |
| `DoD_ICP_2023_minus_2017.tif` | Aligned DEM-of-difference raster | [Zenodo Record](https://zenodo.org/uploads/16813914) |

## Configuration

Edit `config.py` to adjust paths if your data is stored elsewhere. By default, the configuration assumes:
- Data files are in `./data/`
- Outputs will be written to `./outputs/`
- Figures will be saved to `./outputs/figures/`

## Usage

Run the analysis scripts in the following order from the `release-artifacts` directory:

```bash
# 1. Density maps (event and historical)
python src/density_maps.py

# 2. Aspect frequency analysis
python src/aspect_frequency.py

# 3. Frequency-area distributions
python src/fad.py

# 4. Area-volume scaling and cumulative erosion
python src/area_volume.py

# 5. Reactivation distance analysis
python src/reactivation_distance.py
```

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
- `fad_2022_vs_historical.pdf` - FAD comparison with inverse gamma fits

**Area-Volume Analysis:**
- `area_volume_scaling.pdf` - Area vs. volume scatter with power law fit
- `area_volume_table.csv` - Volume statistics table
- `cumulative_erosion_semilog_fit.pdf` - Cumulative erosion curve

**Reactivation Analysis:**
- `cumulative_distance_curve_with_size.pdf` - Reactivation distance curve
- `reactivation_cumulative.csv` - Cumulative reactivation statistics

## Analysis Parameters

Key parameters are defined in `config.py`:

- `target_crs`: Coordinate reference system (EPSG:32617 - UTM Zone 17N)
- `kde_method`: Kernel density estimation bandwidth method ("scott")
- `n_bins`: Number of bins for frequency-area distributions (50)
- `AREA_MIN`, `AREA_MAX`: Power law fitting range (150-10000 m²)
- `EROSION_AREA_MIN`, `EROSION_AREA_MAX`: Erosion analysis range (150-4000 m²)

## Code Structure

```
release-artifacts/
├── config.py              # Configuration and paths
├── environment.yml        # Conda environment specification
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
    └── reactivation_distance.py # Reactivation analysis
```

## Notes

- All scripts can be run independently after data is downloaded
- Progress bars will display for long-running computations (e.g., aspect analysis)
- Figures are saved as PDFs
- CSV tables are saved alongside figures for reference

## Citation

If you use this code, please cite:

```
Scheip, C., Crawford, M., Koch, H., Bibbins, E., (2025). Reconsidering the magnitude of convective storms in triggering landslide events in the Appalachian Plateau, USA. Geophysical Research Letters. [Submitted]
```

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.

You are free to use, modify, and distribute this code for any purpose, including commercial use, with minimal restrictions.

## Contact

For questions or issues, please contact cscheip@bgcengineering.com
