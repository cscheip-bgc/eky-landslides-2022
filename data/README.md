# Data Directory

This directory contains some raw input data files required for the analysis. Additional large data files need to be downloaded separately.

## Files Included in Repository

The following file is included in the repository:

- **`eky-gauge-data.csv`** - Daily precipitation data from eastern Kentucky gauges (Buckhorn, Jackson, Oneida, Whitesburg)
  - Source: [National Weather Service [NWS] NOWData—NOAA online weather data: Jackson, KY (WFO JKL)](https://www.weather.gov/wrh/Climate?wfo=jkl) (Accessed July 7, 2025)
  - This file is small enough to be version controlled

## Required Files to Download

Place the following files in this directory before running the analysis scripts:

### Vector Data (Shapefiles)

1. **`event_inventory_2022.shp`** (+ associated .shx, .dbf, .prj files)
   - 2022 landslide event inventory
   - Source: [OSF Project](https://osf.io/tuvqw/overview?view_only=79501ea561984e198188d3cf4dd7f1c7)
   
2. **`historical_inventory.shp`** (+ associated .shx, .dbf, .prj files)
   - Historical landslide inventory
   - Source: [OSF Project](https://osf.io/tuvqw/overview?view_only=79501ea561984e198188d3cf4dd7f1c7)
   
3. **`mapping_area.shp`** (+ associated .shx, .dbf, .prj files)
   - Study area boundary polygon
   - Source: [Zenodo Record](https://zenodo.org/uploads/16813914)

### Raster Data (GeoTIFF)

4. **`aspect.tif`**
   - DEM-derived aspect raster
   - Source: [Kentucky from Above](https://kyfromabove.ky.gov)
   - Navigate to the study area and download the elevation raster to produce a D8 aspect raster.
   
5. **`lcd_baseline.tif`**
   - Landslide change-detection (LCD) baseline raster
   - Source: [Zenodo Record](https://zenodo.org/uploads/16813914)
   
6. **`dod_aligned.tif`**
   - Aligned DEM-of-difference raster
   - Source: [Zenodo Record](https://zenodo.org/uploads/16813914)

## Data Sources

### OSF (Open Science Framework)
Preliminary landslide inventories are published on OSF. These will eventually be superseded by official Kentucky Geological Survey publications with DOIs.

**OSF Project:** https://osf.io/tuvqw/overview?view_only=79501ea561984e198188d3cf4dd7f1c7

Download the event and historical inventory shapefiles from this project.

### Zenodo
Lidar change detection data, processed DEM products, and study area boundary are archived on Zenodo.

**Zenodo Record:** https://zenodo.org/uploads/16813914

Download the mapping area shapefile, LCD baseline raster, and aligned DEM-of-difference raster from this record.

### Kentucky from Above
Base lidar data and derived products (aspect, slope, etc.) are available from the Kentucky Division of Geographic Information.

**Portal:** https://kyfromabove.ky.gov

Navigate to the study area in eastern Kentucky and download the aspect raster. The study area is centered approximately at coordinates 37.4°N, 83.5°W.

### National Weather Service (NWS)
Daily precipitation data from weather stations in eastern Kentucky is available from the NWS Climate Data portal.

**Portal:** https://www.weather.gov/wrh/Climate?wfo=jkl

The precipitation data includes daily measurements from gauges at Buckhorn, Jackson, Oneida, and Whitesburg.

## File Naming

If your data files have different names, you can either:
1. Rename them to match the expected filenames above, or
2. Edit the paths in `../config.py` to point to your actual filenames

## Coordinate Reference System

All data should be in or will be reprojected to **EPSG:32617** (WGS 84 / UTM Zone 17N).

## Data Not Included

Most data files are not included in the repository. Users must download the spatial data (shapefiles and rasters) from the sources listed above. These large data files are not included due to:
- File size constraints (rasters can be several GB)
- Ongoing data publication process
- Licensing considerations

The precipitation CSV (`eky-gauge-data.csv`) is small (90KB) and is included in the repository for convenience.

## Verification

After downloading the required files, verify you have everything:

```bash
ls -1 data/
```

Expected output (note that `eky-gauge-data.csv` and `README.md` are already included in the repository):
```
aspect.tif
dod_aligned.tif
eky-gauge-data.csv          # included in repo
event_inventory_2022.dbf
event_inventory_2022.prj
event_inventory_2022.shp
event_inventory_2022.shx
historical_inventory.dbf
historical_inventory.prj
historical_inventory.shp
historical_inventory.shx
lcd_baseline.tif
mapping_area.dbf
mapping_area.prj
mapping_area.shp
mapping_area.shx
README.md                   # included in repo
```

## Questions?

For data access issues or questions, please refer to the main README.md or contact [contact information to be added].

