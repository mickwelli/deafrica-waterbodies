import os
import math
import logging
import numpy as np
import geohash as gh
import geopandas as gpd

_log = logging.getLogger(__name__)


def check_wetness_thresholds(minimum_wet_thresholds):
    # Test whether the wetness threshold has been correctly set.

    if minimum_wet_thresholds[0] > minimum_wet_thresholds[-1]:
        error_msg = 'We will be running a hybrid wetness threshold. ' \
            'Please ensure that the primary threshold has a higher value than the ' \
            'secondary threshold. \n'
        raise ValueError(error_msg)
    else:
        print_msg = 'We will be running a hybrid wetness threshold. \n' \
            f'**You have set {minimum_wet_thresholds[-1]} as the ' \
            'primary threshold, which will define the location of the waterbody ' \
            f'polygons \n with {minimum_wet_thresholds[0]} set as the supplementary ' \
            'threshold, which will define the extent/shape of the waterbody polygons.**'
        return (print_msg)


def pp_test_gdf(gdf):
    """
    Function to calculate the Polsby–Popper test values on a
    geopandas GeoDataFrame.
    """
    gdf["area"] = gdf["geometry"].area
    gdf["perimeter"] = gdf["geometry"].length
    gdf["pp_test"] = (4 * math.pi * gdf["area"]) / (gdf["perimeter"] ** 2)

    return gdf


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


def assign_unique_ids(polygons):

    crs = polygons.crs

    # Generate a unique id for each polygon.
    polygons_with_unique_ids = polygons.to_crs(epsg=4326)
    polygons_with_unique_ids['UID'] = polygons_with_unique_ids.apply(lambda x: gh.encode(x.geometry.centroid.y, x.geometry.centroid.x, precision=9), axis=1)

    # Check that our unique ID is in fact unique
    assert polygons_with_unique_ids['UID'].is_unique

    # Make an arbitrary numerical ID for each polygon. We will first sort the dataframe by geohash
    # so that polygons close to each other are numbered similarly.
    polygons_with_unique_ids_sorted = polygons_with_unique_ids.sort_values(by=['UID']).reset_index()
    polygons_with_unique_ids_sorted['WB_ID'] = polygons_with_unique_ids_sorted.index

    # The step above creates an 'index' column, which we don't actually want, so drop it.
    polygons_with_unique_ids_sorted.drop(columns=['index'], inplace=True)

    # Reproject to the same crs as the input polygons.
    polygons_with_unique_ids_sorted = polygons_with_unique_ids_sorted.to_crs(crs)

    return polygons_with_unique_ids_sorted


def get_product_regions(
    product: str,
        ):
    """
    Returns a GeoDataFrame of all the tiles/regions from the product.

    Parameters
    ----------
    product: str
        DE Africa product to get regions for.

    Returns
    -------
    geopandas.geodataframe.GeoDataFrame

    """
    base_url = "https://explorer.digitalearth.africa/api/regions/"
    regions_url = f"{base_url}{product}"

    regions = gpd.read_file(regions_url).drop('count', axis=1)

    return regions


def get_product_tiles(product="wofs_ls_summary_alltime", aoi_gdf=None):

    # Load the product regions.
    regions = get_product_regions(product=product)

    if aoi_gdf is None:
        tiles = regions
    else:
        # Reproject the regions to match the area of interest.
        crs = aoi_gdf.crs
        regions = regions.to_crs(crs)

        tiles = filter_geodataframe_by_intersection(aoi_gdf,
                                                    regions,
                                                    filtertype="intersects",
                                                    invert_mask=False,
                                                    return_inverse=False)
    return tiles


def setup_output_fp(
        ouptut_folder: str = os.getcwd(),
        output_base_filename: str = "waterbodies",
        output_file_type: str = "GeoJSON",):

    valid_output_file_type = ["GeoJSON", "Shapefile"]
    if output_file_type not in valid_output_file_type:
        raise ValueError(f"{output_file_type} is not implemented. Please select from {valid_output_file_type}.")

    if output_file_type == "GeoJSON":
        output_file_extension = ".geojson"
    else:
        output_file_extension = ".shp"

    output_dir_suffix = "waterbodies_outputs"
    output_dir_fp = os.path.join(ouptut_folder, output_dir_suffix)

    if not os.path.exists(output_dir_fp):
        os.makedirs(output_dir_fp)

    final_output_fn = f"{output_base_filename}{output_file_extension}"
    final_output_fp = os.path.join(output_dir_fp, final_output_fn)

    return final_output_fp
