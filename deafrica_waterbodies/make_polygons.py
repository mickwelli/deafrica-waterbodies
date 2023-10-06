"""
Make waterbody polygons from the Water Observations from Space all-time
summary.

Geoscience Australia - 2021
    Claire Krause
    Matthew Alger
"""

import logging

import datacube
import datacube.model
import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from deafrica_tools.spatial import xr_vectorize

from deafrica_waterbodies.filters import filter_by_intersection

_log = logging.getLogger(__name__)


def check_wetness_thresholds(minimum_wet_thresholds: list) -> str:
    """
    Function to validate the wetness thresholds.

    Parameters
    ----------
    minimum_wet_thresholds : list
        A list containing the primary and secondary thresholds, with the secondary
        threshold listed first.

    Returns
    -------
    str
        Validation message

    """
    # Test whether the wetness threshold has been correctly set.

    if minimum_wet_thresholds[0] > minimum_wet_thresholds[-1]:
        _log.error("Primary threshold value is less than the secondary threshold.")
        error_msg = (
            "We will be running a hybrid wetness threshold. "
            "Please ensure that the primary threshold has a higher value than the "
            "secondary threshold. \n"
        )
        raise ValueError(error_msg)
    else:
        print_msg = (
            "We will be running a hybrid wetness threshold. \n"
            f"**You have set {minimum_wet_thresholds[-1]} as the "
            "primary threshold, which will define the location of the waterbody "
            f"polygons \n with {minimum_wet_thresholds[0]} set as the supplementary "
            "threshold, which will define the extent/shape of the waterbody polygons.**"
        )
        return print_msg


def merge_polygons_at_dataset_boundaries(waterbody_polygons: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Function to merge waterbody polygons located at WOfS All Time Summary scene/dataset boundaries.

    Parameters
    ----------
    waterbody_polygons : gpd.GeoDataFrame
        The waterbody polygons.

    Returns
    -------
    gpd.GeoDataFrame
        Waterbody polygons with polygons located at WOfS All Time Summary scene/dataset boundaries merged.
    """
    # Get the dataset extents/regions for the WOfS All Time Summary product.
    ds_extents = gpd.read_file(
        "https://explorer.digitalearth.africa/api/regions/wofs_ls_summary_alltime"
    ).to_crs(waterbody_polygons.crs)

    # Add a 1 pixel (30 m) buffer to the dataset extents.
    buffered_30m_ds_extents_geom = ds_extents.boundary.buffer(
        30, cap_style="flat", join_style="mitre"
    )
    buffered_30m_ds_extents = gpd.GeoDataFrame(
        geometry=buffered_30m_ds_extents_geom, crs=waterbody_polygons.crs
    )

    # Get the polygons at the dataset boundaries.
    boundary_polygons, not_boundary_polygons = filter_by_intersection(
        gpd_data=waterbody_polygons,
        gpd_filter=buffered_30m_ds_extents,
        invert_mask=False,
        return_inverse=True,
    )

    # Now combine overlapping polygons in boundary_polygons.
    merged_boundary_polygons_geoms = shapely.ops.unary_union(boundary_polygons["geometry"])

    # `Explode` the multipolygon back out into individual polygons.
    merged_boundary_polygons = gpd.GeoDataFrame(
        crs=waterbody_polygons.crs, geometry=[merged_boundary_polygons_geoms]
    )
    merged_boundary_polygons = merged_boundary_polygons.explode(index_parts=True).reset_index(
        drop=True
    )

    # Then combine our merged_boundary_polygons with the not_boundary_polygons.
    all_polygons = gpd.GeoDataFrame(
        pd.concat([not_boundary_polygons, merged_boundary_polygons], ignore_index=True, sort=True)
    ).set_geometry("geometry")

    return all_polygons


def get_polygons_from_dataset(
    dataset_id: str,
    dask_chunks: dict[str, int] = {"x": 3200, "y": 3200, "time": 1},
    resolution: tuple[int, int] = (-30, 30),
    output_crs: str = "EPSG:6933",
    min_valid_observations: int = 128,
    primary_threshold: float = 0.1,
    secondary_threshold: float = 0.05,
    dc: datacube.Datacube | None = None,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Generate water body polygons by thresholding a WOfS All Time Summary scene.

    Parameters
    ----------
    dataset_id : str
        The dataset id of a WOfs All Time summary scene for which to
        generate waterbody polygons for.
    dask_chunks : dict, optional
        dask_chunks to use to load WOfS data, by default {"x": 3200, "y": 3200, "time": 1}
    resolution : tuple[int, int], optional
        Resolution to use to load WOfS data, by default (-30, 30)
    output_crs : str, optional
        CRS to load data and for the output waterbody polygons, by default "EPSG:6933"
    min_valid_observations : int, optional
        Threshold to use to mask out pixels based on the number of valid WOfS observations for each pixel, by default 128
    primary_threshold : float, optional
        Threshold to use to determine the location of the waterbody polygons, by default 0.1
    secondary_threshold : float, optional
        Threshold to use to determine the extent / shape of the waterbodies polygons, by default 0.05
    dc : datacube.Datacube | None, optional
        Datacube connection, by default None

    Returns
    -------
    tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]
        A tuple containing GeoDataFrames of waterbody polygons generated from thresholding WOfS All Time Summary data
        using the primary and secondary thresholds.

    """
    # Set up the primary and secondary thresholds.
    minimum_wet_thresholds = [secondary_threshold, primary_threshold]

    # Create a datacube query object.
    query = dict(
        dask_chunks=dask_chunks,
        resolution=resolution,
        output_crs=output_crs,
    )

    # Connect to the datacube.
    if dc is None:
        dc = datacube.Datacube(app="GenerateWaterbodyPolygons")

    # Get the dataset.
    dataset = dc.index.datasets.get(dataset_id)

    # Generate the waterbody polygons using the primary and secondary thresholds,
    # from the dataset/scene.
    try:
        _log.info(f"Generating water body polygons for dataset {dataset_id}")

        # Load the WOfS All-Time Summary dataset.
        wofs_alltime_summary = dc.load(datasets=[dataset], **query).squeeze()

        # Set the no-data values to nan.
        # Masking here is done using the frequency measurement because for multiple
        # areas NaN values are present in the frequency measurement but the
        # no data value -999 is not present in the count_clear and
        # count_wet measurements.
        # Note: it seems some pixels with NaN values in the frequency measurement
        # have a value of zero in the count_clear and/or the count_wet measurements.
        wofs_alltime_summary = wofs_alltime_summary.where(~np.isnan(wofs_alltime_summary.frequency))

        # Mask pixels not observed at least min_valid_observations times.
        wofs_alltime_summary_valid_clear_count = (
            wofs_alltime_summary.count_clear >= min_valid_observations
        )

        # Generate the polygons.
        generated_polygons = {}
        for threshold in minimum_wet_thresholds:
            # Mask any pixels whose frequency of water detection is less than the threshold.
            wofs_alltime_summary_valid_wetness = wofs_alltime_summary.frequency > threshold

            # Now find pixels that meet both the minimum valid observations
            # and minimum wet threshold criteria.
            wofs_alltime_summary_valid = wofs_alltime_summary_valid_wetness.where(
                wofs_alltime_summary_valid_wetness & wofs_alltime_summary_valid_clear_count
            )

            # Convert the raster to polygons.
            # We use a mask of '1' to only generate polygons around values of '1' (not NaNs).
            polygons_mask = wofs_alltime_summary_valid == 1

            polygons = xr_vectorize(
                wofs_alltime_summary_valid,
                mask=polygons_mask,
                crs=wofs_alltime_summary.geobox.crs,
            )

            # Combine any overlapping polygons.
            merged_polygon_geoms = shapely.ops.unary_union(polygons["geometry"])

            # Turn the combined multipolygon back into a GeoDataFrame.
            try:
                merged_polygons = gpd.GeoDataFrame(geometry=list(merged_polygon_geoms.geoms))
            except AttributeError:
                merged_polygons = gpd.GeoDataFrame(geometry=[merged_polygon_geoms])

            # We need to add the crs back onto the GeoDataFrame.
            merged_polygons.crs = wofs_alltime_summary.geobox.crs

            generated_polygons[threshold] = merged_polygons

    except Exception as error:
        _log.exception(
            f"\nDataset {str(dataset_id)} did not run. \n"
            "This is probably because there are no waterbodies present in this scene."
        )
        _log.exception(error)

    primary_threshold_polygons = generated_polygons[primary_threshold]
    secondary_threshold_polygons = generated_polygons[secondary_threshold]

    return primary_threshold_polygons, secondary_threshold_polygons
