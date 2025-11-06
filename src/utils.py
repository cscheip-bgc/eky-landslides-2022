from dataclasses import dataclass
from pathlib import Path
from typing import Tuple

import geopandas as gpd
import numpy as np


@dataclass(slots=True)
class Inventories:
    event: gpd.GeoDataFrame
    historical: gpd.GeoDataFrame


def _prep_inventory(df: gpd.GeoDataFrame, *, drop_min_area: float | None = None) -> gpd.GeoDataFrame:
    df = df.copy()
    df["guid"] = df.index.to_numpy()
    df["area_m2"] = df.geometry.area
    df = df[~df.area_m2.isna() & (df.area_m2 > 0)]
    if drop_min_area is not None:
        df = df[df.area_m2 > drop_min_area]
    return df.reset_index(drop=True)


def load_inventories(cfg) -> Inventories:
    target_crs = cfg.target_crs

    event = gpd.read_file(cfg.inv_2022_path).to_crs(target_crs)
    historical = gpd.read_file(cfg.inv_historical_path).to_crs(target_crs)

    if "modified_d" in event.columns:
        event = event.rename(columns={"modified_d": "modified_descriptor"})

    event = _prep_inventory(event)
    historical = _prep_inventory(historical, drop_min_area=125)

    event["react"] = 0
    try:
        hits = historical.sindex.query_bulk(event.geometry, predicate="intersects")
        if hits.size:
            event.loc[np.unique(hits[0]), "react"] = 1
    except Exception:
        pass

    return Inventories(event=event, historical=historical)


def load_mapping_area(cfg) -> gpd.GeoDataFrame:
    target_crs = cfg.target_crs
    return gpd.read_file(cfg.mapping_area_file).to_crs(target_crs)


def ensure_output_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


