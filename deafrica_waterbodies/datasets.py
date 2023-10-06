import multiprocessing
from functools import partial

import datacube
import geopandas as gpd
import tqdm


def check_ds_intersects_polygons(
    polygons_gdf: gpd.GeoDataFrame | None, ds: datacube.model.Dataset
) -> str:
    """
    Check if the extent of a dataset intersects with a set of polygons.

    Parameters
    ----------
    polygons_gdf : gpd.GeoDataFrame | None
    ds : datacube.model.Dataset

    Returns
    -------
    str | None
        Dataset id if the dataset's extent intersects with the polygons.
    """
    if polygons_gdf is not None:
        # Get the extent of the dataset.
        ds_extent = ds.extent
        # Reproject the extent of the dataset to match the polygons.
        ds_extent = ds_extent.to_crs(polygons_gdf.crs)
        # Get the shapely geometry of the reprojected extent of the dataset.
        ds_extent_geom = ds_extent.geom
        # Check if the dataset's extent intersects with any of the polygons.
        check_intersection = polygons_gdf.geometry.intersects(ds_extent_geom).any()
        if check_intersection:
            return str(ds.id)
        else:
            return ""
    else:
        return str(ds.id)


def filter_datasets(
    dss: list[datacube.model.Dataset], polygons_gdf: gpd.GeoDataFrame | None, num_workers: int = 8
) -> list[str]:
    """
    Filter out datasets that do not intersect with a set of polygons, using a
    multi-process approach to run the `check_ds_intersects_polygons` function.

    Parameters
    ----------
    dss : list[datacube.model.Dataset]
        A list of Datasets to filter.
    polygons_gdf : gpd.GeoDataFrame | None
        A set of polygons in a GeoDataFrame
    num_workers : int, optional
        Number of worker processes to use during multi-processing, by default 8

    Returns
    -------
    list[str]
        A list of the filtered datasets ids.
    """
    with multiprocessing.Pool(processes=num_workers) as pool:
        filtered_datasets_ids_ = list(
            tqdm.tqdm(pool.imap(partial(check_ds_intersects_polygons, polygons_gdf), dss))
        )

    # Remove empty strings.
    filtered_datasets_ids = [item for item in filtered_datasets_ids_ if item]

    return filtered_datasets_ids


def get_datasets_ids(
    aoi_gdf: gpd.GeoDataFrame | None, dc: datacube.Datacube | None = None, num_workers: int = 8
):
    """
    Get the dataset ids of the WOfS All Time Summary datasets whose extents intersect
    with any of the area of interest polygons.

    Parameters
    ----------
    aoi_gdf : gpd.GeoDataFrame | None
        Area of interest
    dc : datacube.Datacube | None, optional
        Datacube connection, by default None
    num_workers : int, optional
        Number of worker processes to use when filtering datasets, by default 8
    Returns
    -------
    str
        Dataset ids of the WOfS All Time Summary datasets whose extents intersect
        with any of the area of interest polygons.
    """
    # Connect to the datacube.
    if dc is None:
        dc = datacube.Datacube(app="WaterbodiesPolygons")

    # Find all datasets available for the WOfS All Time summary product.
    dss = dc.index.datasets.search(product=["wofs_ls_summary_alltime"])
    dss = list(dss)

    # Filter the datasets to the area of interest.
    filtered_datasets_ids = filter_datasets(dss, aoi_gdf, num_workers)

    return filtered_datasets_ids
