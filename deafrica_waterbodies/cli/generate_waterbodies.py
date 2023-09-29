import logging
import math

import click
import geopandas as gpd
import shapely

from deafrica_waterbodies.cli.group_options import MutuallyExclusiveOption
from deafrica_waterbodies.cli.io import write_waterbodies_to_file
from deafrica_waterbodies.cli.logs import logging_setup
from deafrica_waterbodies.waterbodies.polygons.attributes import (
    add_area_and_perimeter_attributes,
    add_timeseries_attribute,
)
from deafrica_waterbodies.waterbodies.polygons.make_polygons import get_waterbodies


@click.command("generate-waterbodies", no_args_is_help=True)
@click.option(
    "--run-type",
    type=click.Choice(["continental", "custom"]),
    default=["custom"],
    help="If run-type is continental, generate waterbodies for all the "
    " WOfS All Time Summary regions. If run-type is custom, pass a vector "
    "file defining the area of interest to --vector-file-fp OR the coordinates "
    "of a bounding box for the area of interest to --bbox.",
)
@click.option(
    "--vector-file-fp",
    cls=MutuallyExclusiveOption,
    mutually_exclusive=["bbox"],
    help="The path to a vector defining the area of interest.",
)
@click.option(
    "--bbox",
    cls=MutuallyExclusiveOption,
    mutually_exclusive=["vector_file_fp"],
    help="Coordinates of the area of interest's bounding box in the format 'minx,miny,maxx,maxy'.",
)
@click.option(
    "--bbox-crs",
    default="EPSG:4326",
    help="CRS of the bounding box coordinates.",
    show_default=True,
)
@click.option(
    "--primary-threshold",
    default=0.1,
    type=click.FLOAT,
    help="Threshold to define the location of the waterbody polygons.",
    show_default=True,
)
@click.option(
    "--secondary-threshold",
    default=0.05,
    type=click.FLOAT,
    help="Threshold to define the extent/shape of the waterbody polygons.",
    show_default=True,
)
@click.option(
    "--min-polygon-size",
    default=4500,
    show_default=True,
    type=click.FLOAT,
    help="Minimum area in m2 of the waterbody polygons to be included.",
)
@click.option(
    "--max-polygon-size",
    type=click.FLOAT,
    default=math.inf,
    show_default=True,
    help="Maximum area in m2 of the waterbody polygons to be included.",
)
@click.option(
    "--remove-ocean-polygons",
    is_flag=True,
    default=False,
    show_default=True,
    help="Filter out ocean polygons from the waterbody polygons.",
)
@click.option(
    "--land-sea-mask-fp",
    type=click.Path(),
    default=None,
    help="File path to vector dataset to use to filter out ocean polygons.",
)
@click.option(
    "--remove-major-rivers",
    is_flag=True,
    default=False,
    help="Filter out major rivers from the waterbody polygons.",
)
@click.option(
    "--major-rivers-mask-fp",
    type=click.Path(),
    default=None,
    help="File path to vector dataset to use to filter out major rivers.",
)
@click.option(
    "--remove-cbd", is_flag=True, default=False, help="Filter out CBDs from the waterbody polygons."
)
@click.option(
    "--urban-mask-fp",
    type=click.Path(),
    default=None,
    help="File path to vector dataset to use to filter out urban/CBD areas.",
)
@click.option(
    "--handle-large-polygons",
    default="nothing",
    type=click.Choice(["erode-dilate-v1", "erode-dilate-v2", "nothing"]),
    show_default=True,
    help="Method to use to split large polygons above the Polsby–Popper test value.",
)
@click.option(
    "--pp-test-threshold",
    type=click.FLOAT,
    default=0.005,
    show_default=True,
    help="Polsby–Popper test value to use to split large polygons.",
)
@click.option("-v", "--verbose", count=True)
@click.option(
    "--product-version",
    type=str,
    default="0.0.1",
    show_default=True,
    help="Product version for the DE Africa Waterbodies product.",
)
@click.option("--s3", "storage_location", flag_value="s3", help="Save the output to an s3 bucket.")
@click.option(
    "--local",
    "storage_location",
    flag_value="local",
    default=True,
    help="Save the output to a local folder.",
)
@click.option(
    "--output-bucket-name",
    type=str,
    show_default=True,
    cls=MutuallyExclusiveOption,
    mutually_exclusive=["output_local_folder"],
    help="The s3 bucket to write the output to.",
)
@click.option(
    "--output-local-folder",
    type=click.Path(),
    cls=MutuallyExclusiveOption,
    mutually_exclusive=["output_bucket_name"],
    help="Local directory to write the waterbody polygons to.",
)
@click.option(
    "--output-file-name",
    default="waterbodies",
    show_default=True,
    help="File name for the output waterbodies.",
)
@click.option(
    "--output-file-type",
    default="GeoJSON",
    show_default=True,
    help="File type for the outputs waterbodies.",
    type=click.Choice(["GeoJSON", "Shapefile"], case_sensitive=False),
)
def generate_waterbodies(
    run_type,
    vector_file_fp,
    bbox,
    bbox_crs,
    primary_threshold,
    secondary_threshold,
    min_polygon_size,
    max_polygon_size,
    remove_ocean_polygons,
    land_sea_mask_fp,
    remove_major_rivers,
    major_rivers_mask_fp,
    remove_cbd,
    urban_mask_fp,
    handle_large_polygons,
    pp_test_threshold,
    verbose,
    product_version,
    storage_location,
    output_bucket_name,
    output_local_folder,
    output_file_name,
    output_file_type,
):
    """
    Generate waterbodies for all the WOfS All Time Summary regions.
    """
    logging_setup(verbose)
    _log = logging.getLogger(__name__)

    dask_chunks = {"x": 3000, "y": 3000, "time": 1}
    resolution = (-30, 30)
    output_crs = "EPSG:6933"
    min_valid_observations = 128

    if run_type == "custom":
        if not vector_file_fp and not bbox:
            _log.error("No vector file or bounding box for the area of interest specified.")
            error_msg = (
                "If not running for all of Africa please specify an "
                "area of interest by EITHER passing a vector file path to --vector-file-fp "
                "OR passing bounding box coordinates to --bbox"
            )
            raise ValueError(error_msg)
        elif vector_file_fp and bbox:
            _log.error("Both vector file and bounding box for the area of interest specified.")
            error_msg = (
                "Please specify EITHER a vector file path using --vector-file-fp OR  bounding box "
                "coordinates using --bbox ."
            )
            raise ValueError(error_msg)
        elif vector_file_fp and not bbox:
            # Read the vector file.
            try:
                aoi_gdf = gpd.read_file(vector_file_fp)
            except Exception as error:
                _log.exception(f"Could not read the vector file {vector_file_fp}")
                raise error
            else:
                aoi_gdf = aoi_gdf.to_crs(output_crs)
                continental_run = False
        elif bbox and not vector_file_fp:
            try:
                # Convert bounding box to GeoDataFrame.
                bbox_ = [float(i.strip()) for i in bbox.split(",")]
                aoi_gdf = gpd.GeoDataFrame(geometry=[shapely.geometry.box(*bbox_)], crs=bbox_crs)
            except Exception as error:
                _log.exception(
                    f"Error in defining the bounding box {bbox} using the crs {bbox_crs}."
                )
                raise error
            else:
                aoi_gdf = aoi_gdf.to_crs(output_crs)
                continental_run = False
    elif run_type == "continental":
        if vector_file_fp or bbox:
            _log.error(
                "Vector file or bounding box defining an area of interest specified, yet run type is continental."
            )
            error_msg = (
                "If running for all of Africa do not pass a vector file to "
                "--vector-file-fp or bounding box coordinates to --bbox"
            )
            raise ValueError(error_msg)
        else:
            aoi_gdf = None
            continental_run = True

    if remove_ocean_polygons:
        filter_out_ocean_polygons = True
    else:
        filter_out_ocean_polygons = False

    if remove_major_rivers:
        filter_out_major_rivers_polygons = True
    else:
        filter_out_major_rivers_polygons = False

    if remove_cbd:
        filter_out_urban_polygons = True
    else:
        filter_out_urban_polygons = False

    waterbodies_gdf = get_waterbodies(
        aoi_gdf=aoi_gdf,
        continental_run=continental_run,
        dask_chunks=dask_chunks,
        resolution=resolution,
        output_crs=output_crs,
        min_valid_observations=min_valid_observations,
        primary_threshold=primary_threshold,
        secondary_threshold=secondary_threshold,
        min_polygon_size=min_polygon_size,
        max_polygon_size=max_polygon_size,
        filter_out_ocean_polygons=filter_out_ocean_polygons,
        land_sea_mask_fp=land_sea_mask_fp,
        filter_out_major_rivers_polygons=filter_out_major_rivers_polygons,
        major_rivers_mask_fp=major_rivers_mask_fp,
        filter_out_urban_polygons=filter_out_urban_polygons,
        urban_mask_fp=urban_mask_fp,
        handle_large_polygons=handle_large_polygons,
        pp_test_threshold=pp_test_threshold,
    )

    waterbodies_gdf = add_area_and_perimeter_attributes(waterbodies_gdf)
    waterbodies_gdf = add_timeseries_attribute(waterbodies_gdf, product_version, output_bucket_name)

    # Reproject to EPSG:4326
    waterbodies_gdf_4326 = waterbodies_gdf.to_crs("EPSG:4326")

    write_waterbodies_to_file(
        waterbodies_gdf_4326,
        product_version,
        storage_location,
        output_bucket_name,
        output_local_folder,
        output_file_name,
        output_file_type,
    )
