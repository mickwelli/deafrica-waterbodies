import multiprocessing
from functools import partial

import datacube
import geopandas as gpd
import tqdm
from odc.dscache.tools.tiling import parse_gridspec_with_name


def check_tile_intersects_polygons(
    polygons_gdf: gpd.GeoDataFrame | None,
    tile: tuple[tuple[int, int], datacube.api.grid_workflow.Tile],
) -> tuple[int, int]:
    """
    Check if the extent of the geobox of a tile intersects with a set of polygons.

    Parameters
    ----------
    polygons_gdf : gpd.GeoDataFrame | None
    tile : tuple[tuple[int,int], datacube.api.grid_workflow.Tile]

    Returns
    -------
    tuple
        Tile id if the extent of the geobox of a tile intersects with the polygons.
    """
    tile_id = tile[0]
    tile_extent = tile[1].geobox.extent

    if polygons_gdf is not None:
        # Reproject the extent of the geobox of a tile to match the polygons.
        tile_extent = tile_extent.to_crs(polygons_gdf.crs)
        # Get the shapely geometry of the reprojected extent of the tile's geobox.
        tile_extent_geom = tile_extent.geom
        # Check if the extent intersects with any of the polygons.
        check_intersection = polygons_gdf.geometry.intersects(tile_extent_geom).any()
        if check_intersection:
            return tile_id
        else:
            return ()
    else:
        return tile_id


def filter_tiles(
    tiles: dict[tuple[int, int], datacube.api.grid_workflow.Tile],
    polygons_gdf: gpd.GeoDataFrame | None,
    num_workers: int = 8,
) -> list[tuple[int, int]]:
    """
    Filter out tiles that do not intersect with a set of polygons, using a
    multi-process approach to run the `check_tile_intersects_polygons` function.

    Parameters
    ----------
    tiles: dict[tuple[int,int], datacube.api.grid_workflow.Tile]
        A dictionary of tiles.
    polygons_gdf : gpd.GeoDataFrame | None
        A set of polygons in a GeoDataFrame
    num_workers : int, optional
        Number of worker processes to use during multi-processing, by default 8

    Returns
    -------
    list[tuple[int, int]]
        A list of the filtered tiles ids.
    """
    with multiprocessing.Pool(processes=num_workers) as pool:
        filtered_tile_ids_ = list(
            tqdm.tqdm(
                pool.imap(partial(check_tile_intersects_polygons, polygons_gdf), tiles.items())
            )
        )

    # Remove empty tuples.
    filtered_tile_ids = [item for item in filtered_tile_ids_ if item]

    return filtered_tile_ids
