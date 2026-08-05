"""Project configuration for the release-artifacts figure kit.

Edit the paths below to point to your local data files. By default, this config
assumes data files are placed in the 'data' subdirectory and outputs will be
written to the 'outputs' subdirectory within the release-artifacts folder.

Adjust these paths as needed for your local setup."""

from pathlib import Path

# Root locations -----------------------------------------------------------------
# Get the release-artifacts directory (parent of this config file)
release_root = Path(__file__).resolve().parent

# Data directory (place your downloaded data files here)
data_path = release_root / "data"

# Output directory (figures and tables will be saved here)
output_path = release_root / "outputs"
output_path.mkdir(parents=True, exist_ok=True)

figure_path = output_path / "figures"
figure_path.mkdir(parents=True, exist_ok=True)

# Inventory + study boundary inputs ---------------------------------------------
inv_2022_path = data_path / "event-inventory.shp"
inv_historical_path = data_path / "historical-inventory.shp"

mapping_area_file = data_path / "mapping-area.shp"

# Terrain rasters ----------------------------------------------------------------
aspect_file = data_path / "aspect.tif"

# LCD/DoD rasters ----------------------------------------------------------------
lcd_file = data_path / "m3c2_lcd.tif"
dod_aligned_file = data_path / "dod_aligned.tif"

# Independent landslide susceptibility model ------------------------------------
susceptibility_file = data_path / "event-area-susceptibility.tif"

# Precomputed erosion volumes (from raster extraction) -------------------------
erosion_volumes_csv = data_path / "landslide_erosion_volumes.csv"

# Precipitation gauge records ----------------------------------------------------
precip_gauge_csv = data_path / "eky-gauge-data.csv"


# Analysis parameters ------------------------------------------------------------
target_crs = "EPSG:32617"
kde_method = "scott"

n_bins = 50

AREA_MIN = 151
AREA_MAX = 10000

EROSION_AREA_MIN = 151
EROSION_AREA_MAX = 4000

# Susceptibility vs. landslide density comparison --------------------------------
SUSC_CELL_SIZE_M = 1000.0        # primary aggregation cell size
SUSC_SENSITIVITY_CELL_SIZE_M = 2000.0
SUSC_MIN_CELL_COVERAGE = 0.90    # drop partially covered edge cells
SUSC_HIGH_THRESHOLDS = (0.5, 0.6)

