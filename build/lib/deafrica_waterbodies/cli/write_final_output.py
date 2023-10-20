import logging
import os

import click
import geopandas as gpd

from deafrica_waterbodies.attributes import (
    add_area_and_perimeter_attributes,
    add_timeseries_attribute,
    assign_unique_ids,
)
from deafrica_waterbodies.cli.logs import logging_setup
from deafrica_waterbodies.io import write_waterbodies_to_file


@click.command("write-final-output", no_args_is_help=True)
@click.option("-v", "--verbose", count=True)
@click.option(
    "--output-directory",
    type=click.Path(),
    help="Directory containing the waterbody polygons.",
)
@click.option(
    "--product-version",
    type=str,
    default="0.0.1",
    show_default=True,
    help="Product version for the DE Africa Waterbodies product.",
)
@click.option(
    "--timeseries-bucket",
    type=str,
    show_default=True,
    help="The s3 bucket to containing the timeseries for the polygons.",
)
def write_final_output(
    verbose,
    output_directory,
    product_version,
    timeseries_bucket,
):
    # Set up logger.
    logging_setup(verbose=verbose)
    _log = logging.getLogger(__name__)

    # Support pathlib paths.
    output_directory = str(output_directory)

    _log.info("Loading filtered waterbody polygons...")
    filtered_polygons_fp = os.path.join(output_directory, "filtered_polygons.parquet")
    filtered_polygons = gpd.read_parquet(filtered_polygons_fp)
    _log.info(f"Waterbody polygons count {len(filtered_polygons)}.")

    waterbodies_gdf = assign_unique_ids(polygons=filtered_polygons)
    waterbodies_gdf = add_area_and_perimeter_attributes(polygons=waterbodies_gdf)
    waterbodies_gdf = add_timeseries_attribute(
        polygons=waterbodies_gdf,
        product_version=product_version,
        timeseries_bucket=timeseries_bucket,
    )

    # Reproject to EPSG:4326
    waterbodies_gdf_4326 = waterbodies_gdf.to_crs("EPSG:4326")

    # Write to disk.
    write_waterbodies_to_file(
        waterbodies_gdf=waterbodies_gdf_4326,
        product_version=product_version,
        output_directory=output_directory,
    )
