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
inv_2022_path = data_path / "event_inventory_2022.shp"
inv_historical_path = data_path / "historical_inventory.shp"

mapping_area_file = data_path / "mapping_area.shp"

# Terrain rasters ----------------------------------------------------------------
aspect_file = data_path / "aspect.tif"

# LCD/DoD rasters ----------------------------------------------------------------
lcd_file = data_path / "M3C2_ICP_2017_vs_2023.tif"
dod_aligned_file = data_path / "DoD_ICP_2023_minus_2017.tif"


# Analysis parameters ------------------------------------------------------------
target_crs = "EPSG:32617"
kde_method = "scott"

n_bins = 50

AREA_MIN = 150
AREA_MAX = 10000

EROSION_AREA_MIN = 150
EROSION_AREA_MAX = 4000

