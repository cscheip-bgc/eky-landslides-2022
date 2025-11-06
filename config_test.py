"""Project configuration for the release-artifacts figure kit.

Edit the absolute paths below so they point to the local copies of each data
source before running any scripts.  This file intentionally performs no path
validation to keep things simple."""

from pathlib import Path

# Root locations -----------------------------------------------------------------
project_path = Path("/Users/cscheip/Library/CloudStorage/OneDrive-BGCEngineeringInc/P24236_Document_Set/project_folder/writing/paper-1-intro/data-for-repos/code/test")

output_path = project_path / "outputs"
output_path.mkdir(parents=True, exist_ok=True)

figure_path = output_path / "figures"
figure_path.mkdir(parents=True, exist_ok=True)

# Inventory + study boundary inputs ---------------------------------------------
inv_2022_path = project_path / "data" / "Event_mapped_landslides.shp"
inv_historical_path = project_path / "data" / "AllHistoricPolygons_InStudy.shp"

mapping_area_file = project_path / "data" / "usgs_study_area_rev1.shp"

# Terrain rasters ----------------------------------------------------------------
aspect_file = project_path / "data" / "2023-01-01_aspect.tif"

# LCD rasters ----------------------------------------------------------------
lcd_file = project_path / "data" / "128_ALS_2017-01-01_vs_2023-01-01_mapping_area.tif"
dod_aligned_file = project_path / "data" / "2023-minus-2017-icp.tif"


# Analysis parameters ------------------------------------------------------------
target_crs = "EPSG:32617"
kde_method = "scott"

n_bins = 50

AREA_MIN = 150
AREA_MAX = 10000

EROSION_AREA_MIN = 150
EROSION_AREA_MAX = 4000

