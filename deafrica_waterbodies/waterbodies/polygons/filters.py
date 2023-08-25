"""
Filter waterbody polygons based on different criteria.
"""
import math
import logging
import warnings
import numpy as np
import pandas as pd
import geopandas as gpd

_log = logging.getLogger(__name__)


def filter_geodataframe_by_intersection(gpd_data,
                                        gpd_filter,
                                        filtertype="intersects",
                                        invert_mask=True,
                                        return_inverse=False):
    """
    Filter out polygons that intersect with another polygon shapefile.

    Parameters
    ----------
    gpd_data: geopandas GeoDataFrame
        Polygon data to be filtered.
    gpd_filter: geopandas GeoDataFrame
        Polygon dataset to be used as a filter.

    Optional
    --------
    filtertype: default = 'intersects'
        Options = ['intersects', 'contains', 'within']
    invert_mask: boolean
        Default = 'True'. This determines whether you want areas that
        DO ( = 'False') or DON'T ( = 'True') intersect with the filter dataset.
    return_inverse: boolean
        Default = 'False'. If True, then return both parts of the intersection:
        - those that intersect AND
        - those that don't as two GeoDataFrames.

    Returns
    -------
    gpd_data_filtered: geopandas GeoDataFrame
        If invert_mask==True, `gpd_data_filtered` is a filtered polygon set,
        with polygons that DO intersect with gpd_filter removed.
        If invert_mask==False, `gpd_data_filtered` is a filtered polygon set,
        with polygons that DON'T intersect with gpd_filter removed.
    intersect_index: list of indices of gpd_data that intersect with gpd_filter.

    Optional
    --------
    if 'return_inverse = True'
    gpd_data_inverse: geopandas GeoDataFrame
        If invert_mask==True, `gpd_data_inverse` is a filtered polygon set,
        with polygons that DON'T intersect with gpd_filter removed (inverse of gpd_data_filtered).
        If invert_mask==False, `gpd_data_inverse` is a filtered polygon set,
        with polygons that DO intersect with gpd_filter removed (inverse of gpd_data_filtered).

    """

    # Check that the coordinate reference systems of both GeoDataFrames are the same.
    assert gpd_data.crs == gpd_filter.crs

    # Find the index of all the polygons in gpd_data that intersect with gpd_filter.
    intersections = gpd_filter.sjoin(gpd_data,
                                     how="inner",
                                     predicate=filtertype)
    intersect_index = np.sort(intersections["index_right"].unique())

    if invert_mask:
        # Grab only the polygons that are NOT in the intersect_index.
        gpd_data_filtered = gpd_data.loc[~gpd_data.index.isin(intersect_index)]
    else:
        # Grab only the polygons that ARE in the intersect_index.
        gpd_data_filtered = gpd_data.loc[gpd_data.index.isin(intersect_index)]

    if return_inverse:
        # We need to use the indices from intersect_index to find the inverse dataset, so we
        # will just swap the '~'.

        if invert_mask:
            # Grab only the polygons that ARE in the intersect_index.
            gpd_data_inverse = gpd_data.loc[gpd_data.index.isin(
                intersect_index)]
        else:
            # Grab only the polygons that are NOT in the intersect_index.
            gpd_data_inverse = gpd_data.loc[~gpd_data.index.isin(intersect_index)]

        return gpd_data_filtered, intersect_index, gpd_data_inverse
    else:
        return gpd_data_filtered, intersect_index


def pp_test_gdf(gdf):
    """
    Function to calculate the Polsby–Popper test values on a
    geopandas GeoDataFrame.
    """
    gdf["area"] = gdf["geometry"].area
    gdf["perimeter"] = gdf["geometry"].length
    gdf["pp_test"] = (4 * math.pi * gdf["area"]) / (gdf["perimeter"] ** 2)

    return gdf


def split_large_polygons(input_gdf, pp_thresh: int = 0.005, method="nothing"):

    warnings.simplefilter(action='ignore', category=FutureWarning)

    # Calculate the Polsby–Popper values.
    input_gdf = pp_test_gdf(input_gdf)

    # Confirm the option to use.
    valid_options = ['erode-dilate-v1', 'erode-dilate-v2', 'nothing']
    if method not in valid_options:
        _log.info(f"{method} method not implemented to handle large polygons. Defaulting to not splitting large polygons.")
        method = "nothing"

    # Split large polygons.
    if method == 'nothing':
        info_msg = 'You have chosen not to split large polygons. If you meant to use this option, please ' \
            f'select one of the following methods: {valid_options[:2]}.'
        _log.info(info_msg)
        large_polygons_handled = input_gdf
    else:
        _log.info(f"Splitting large polygons using the `{method}` method, using the threshold {pp_thresh}.")
        if method == 'erode-dilate-v1':
            splittable_polygons = input_gdf[input_gdf.pp_test <= pp_thresh]
            not_splittable_polygons = input_gdf[input_gdf.pp_test > pp_thresh]

            splittable_polygons_buffered = splittable_polygons.buffer(-50)
            split_polygons = splittable_polygons_buffered.explode(index_parts=True).reset_index(drop=True).buffer(50)
            split_polygons_gdf = gpd.GeoDataFrame(geometry=split_polygons, crs=input_gdf.crs)
            split_polygons_gdf = pp_test_gdf(split_polygons_gdf)

            large_polygons_handled = pd.concat([not_splittable_polygons, split_polygons_gdf], ignore_index=True)
        elif method == 'erode-dilate-v2':
            splittable_polygons = input_gdf[input_gdf.pp_test <= pp_thresh]
            not_splittable_polygons = input_gdf[input_gdf.pp_test > pp_thresh]

            if len(splittable_polygons) >= 1:
                splittable_polygons_buffered = splittable_polygons.buffer(-100)
                splittable_polygons_buffered = splittable_polygons_buffered.buffer(125)

                splittable_polygons_buffered_union = gpd.GeoDataFrame(geometry=[splittable_polygons_buffered.unary_union], crs=input_gdf.crs)
                subtracted = gpd.overlay(splittable_polygons, splittable_polygons_buffered_union, how='difference').explode(index_parts=True).reset_index(drop=True)
                resubtracted = gpd.overlay(splittable_polygons, subtracted, how='difference').explode(index_parts=True).reset_index(drop=True)

                # Assign each chopped-off bit of the polygon to its nearest big
                # neighbour.
                unassigned = np.ones(len(subtracted), dtype=bool)
                recombined = []

                for row in resubtracted.itertuples():
                    mask = subtracted.exterior.intersects(row.geometry.exterior) & unassigned
                    neighbours = subtracted[mask]
                    unassigned[mask] = False
                    poly = row.geometry.union(neighbours.unary_union)
                    recombined.append(poly)

                recombined_gdf = gpd.GeoDataFrame(geometry=recombined, crs=input_gdf.crs)
                # Get only the actual geometry objects that are neither missing nor empty
                recombined_gdf_masked = recombined_gdf[~(recombined_gdf.geometry.is_empty | recombined_gdf.geometry.isna())]
                # All remaining polygons are not part of a big polygon.
                results = pd.concat([recombined_gdf_masked, subtracted[unassigned], not_splittable_polygons], ignore_index=True)

                results = pp_test_gdf(results)

                large_polygons_handled = results.explode(index_parts=True).reset_index(drop=True)

            else:
                info_msg = f"There are no polygons with a Polsby–Popper score above the {pp_thresh}. " \
                    "No polygons were split."
                _log.info(info_msg)
                large_polygons_handled = input_gdf

    large_polygons_handled.drop(columns=['area', 'perimeter', 'pp_test'], inplace=True)

    return large_polygons_handled


def filter_waterbodies(
        primary_threshold_polygons,
        secondary_threshold_polygons,
        min_polygon_size: float = 4500,
        max_polygon_size: float = math.inf,
        filter_out_ocean_polygons: bool = False,
        land_sea_mask_fp: str = None,
        filter_out_major_rivers_polygons: bool = False,
        major_rivers_mask_fp: str = None,
        filter_out_urban_polygons: bool = False,
        urban_mask_fp: str = None,
        handle_large_polygons: str = "nothing",
        pp_test_threshold: float = 0.005,
        ):

    # Assert both polygons have the same crs
    assert primary_threshold_polygons.crs == secondary_threshold_polygons.crs

    crs = primary_threshold_polygons.crs

    _log.info(f"Filtering the waterbody polygons using the minimum area {min_polygon_size} and the maximum area {max_polygon_size}.")
    primary_threshold_polygons["area"] = primary_threshold_polygons.area
    secondary_threshold_polygons["area"] = secondary_threshold_polygons.area

    area_filtered_primary = primary_threshold_polygons.loc[((primary_threshold_polygons['area'] > min_polygon_size) & (primary_threshold_polygons['area'] <= max_polygon_size))]
    area_filtered_secondary = secondary_threshold_polygons.loc[(secondary_threshold_polygons['area'] <= max_polygon_size)]

    # Filter out polygons that intersect with the ocean.
    if filter_out_ocean_polygons:
        _log.info("Filtering out ocean waterbody polygons")
        if land_sea_mask_fp is None:
            error_msg = 'Please provide a file path to the dataset you wish to use ' \
                'to filter out ocean polygons. The dataset needs to be a vector dataset, ' \
                'and able to be read in by the fiona python library. '
            raise ValueError(error_msg)
        else:
            try:
                land_sea_mask = gpd.read_file(land_sea_mask_fp).to_crs(crs)
            except Exception as error:
                _log.error(error)
                raise
            ocean_filtered_primary, _ = filter_geodataframe_by_intersection(area_filtered_primary, land_sea_mask, invert_mask=True)
            ocean_filtered_secondary, _ = filter_geodataframe_by_intersection(area_filtered_secondary, land_sea_mask, invert_mask=True)
    else:
        ocean_filtered_primary = area_filtered_primary
        ocean_filtered_secondary = area_filtered_secondary

    # Filter out CBD areas.
    if filter_out_urban_polygons:
        _log.info("Filtering out urban/CBD areas from waterbody polygons")
        if urban_mask_fp is None:
            error_msg = 'Please provide a file path to the dataset you wish to use ' \
                'to filter out urban/CBD area polygons. The dataset needs to be a vector dataset, ' \
                'and able to be read in by the fiona python library. '
            raise ValueError(error_msg)
        else:
            try:
                urban_mask = gpd.read_file(urban_mask_fp).to_crs(crs)
            except Exception as error:
                _log.error(error)
                raise
            cbd_filtered_primary = filter_geodataframe_by_intersection(ocean_filtered_primary, urban_mask)
            cbd_filtered_secondary = ocean_filtered_secondary
    else:
        info_msg = 'You have chosen not to filter out waterbodies within CBDs. If you meant to use this option, please ' \
            'set `filter_out_urban_areas = True`, and set the path to your urban filter shapefile in `urban_mask_fp`.'
        _log.info(info_msg)
        cbd_filtered_primary = ocean_filtered_primary
        cbd_filtered_secondary = ocean_filtered_secondary

    # Find the polygons identified using the secondary threshold that intersect with those identified
    # using the primary threshold.
    _log.info("Finding polygons identified using the secondary threshold that intersect with those identified using the primary threshold.")

    _, intersect_indices = filter_geodataframe_by_intersection(cbd_filtered_secondary, cbd_filtered_primary)

    do_intersect_with_primary = cbd_filtered_secondary.loc[cbd_filtered_secondary.index.isin(intersect_indices)]

    # Concat the two polygon datasets together.
    combined_polygons = gpd.GeoDataFrame(pd.concat([do_intersect_with_primary, cbd_filtered_primary], ignore_index=True))
    # Merge overlapping polygons.
    merged_combined_polygons_geoms = combined_polygons.unary_union

    # `Explode` the multipolygons back out into individual polygons.
    merged_combined_polygons = gpd.GeoDataFrame(crs=crs, geometry=[merged_combined_polygons_geoms])
    merged_combined_polygons = merged_combined_polygons.explode(index_parts=True).reset_index(drop=True)

    # Filter out polygons that intersect with major rivers.
    if filter_out_major_rivers_polygons:
        _log.info("Filtering out major rivers from the waterbodies.")
        if major_rivers_mask_fp is None:
            error_msg = 'Please provide a file path to the dataset you wish to use ' \
                'to filter out major rivers. The dataset needs to be a vector dataset, ' \
                'and able to be read in by the fiona python library. '
            raise ValueError(error_msg)
        else:
            try:
                major_rivers = gpd.read_file(major_rivers_mask_fp).to_crs(crs)
            except Exception as error:
                _log.error(error)
                raise
            major_rivers_filtered, _ = filter_geodataframe_by_intersection(merged_combined_polygons, major_rivers)
    else:
        info_msg = 'You have chosen not to filter out major rivers from the waterbodies. If you meant to use this option, please ' \
            'set `filter_out_rivers = True`, and set the path to your major rivers shapefile in `major_rivers_fp`.'
        _log.info(info_msg)
        major_rivers_filtered = merged_combined_polygons

    # Handling large polygons.
    large_polygons_handled = split_large_polygons(major_rivers_filtered, pp_thresh=pp_test_threshold, method=handle_large_polygons)

    # Reapply the size filtering, just to check that all of the split and filtered waterbodies are
    # still in the size range we want.
    large_polygons_handled["area"] = large_polygons_handled.area
    filtered_polygons = large_polygons_handled.loc[((large_polygons_handled['area'] > min_polygon_size) & (large_polygons_handled['area'] <= max_polygon_size))]

    filtered_polygons.drop(columns=['area'], inplace=True)

    return filtered_polygons
