import logging
import math
import os

import click
import geopandas as gpd

from deafrica_waterbodies.cli.logs import logging_setup
from deafrica_waterbodies.filters import filter_waterbodies


@click.command("filter-waterbody-polygons", no_args_is_help=True)
@click.option("-v", "--verbose", count=True)
@click.option(
    "--output-directory",
    type=click.Path(),
    help="Directory containing the waterbody polygons.",
)
@click.option(
    "--min-polygon-size",
    default=4500,
    show_default=True,
    help="Minimum area in m2 of the waterbody polygons to be included.",
)
@click.option(
    "--max-polygon-size",
    default=math.inf,
    show_default=True,
    help="Maximum area in m2 of the waterbody polygons to be included.",
)
@click.option(
    "--land-sea-mask-fp",
    default="",
    help="File path to vector dataset to use to filter out ocean polygons.",
)
@click.option(
    "--urban-mask-fp",
    type=click.Path(),
    default="",
    help="File path to vector dataset to use to filter out urban/CBD areas.",
)
@click.option(
    "--major-rivers-mask-fp",
    type=click.Path(),
    default="",
    help="File path to vector dataset to use to filter out major rivers.",
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
    default=0.005,
    show_default=True,
    help="Polsby–Popper test value to use to split large polygons.",
)
def filter_waterbody_polygons(
    verbose,
    output_directory,
    min_polygon_size,
    max_polygon_size,
    land_sea_mask_fp,
    urban_mask_fp,
    major_rivers_mask_fp,
    handle_large_polygons,
    pp_test_threshold,
):
    """
    Filter the primary and secondary threshold waterbody polygons.
    """
    # Set up logger.
    logging_setup(verbose=verbose)
    _log = logging.getLogger(__name__)

    # Support pathlib paths.
    output_directory = str(output_directory)

    _log.info("Loading primary and secondary threshold polygons...")

    primary_threshold_polygons_fp = os.path.join(
        output_directory, "primary_threshold_polygons_merged_at_ds_boundaries.parquet"
    )
    secondary_threshold_polygons_fp = os.path.join(
        output_directory, "secondary_threshold_polygons_merged_at_ds_boundaries.parquet"
    )
    primary_threshold_polygons = gpd.read_parquet(primary_threshold_polygons_fp)
    secondary_threshold_polygons = gpd.read_parquet(secondary_threshold_polygons_fp)

    _log.info(f"Primary threshold polygons count {len(primary_threshold_polygons)}.")
    _log.info(f"Secondary threshold polygons count {len(secondary_threshold_polygons)}.")

    filtered_polygons = filter_waterbodies(
        primary_threshold_polygons=primary_threshold_polygons,
        secondary_threshold_polygons=secondary_threshold_polygons,
        min_polygon_size=min_polygon_size,
        max_polygon_size=max_polygon_size,
        land_sea_mask_fp=land_sea_mask_fp,
        urban_mask_fp=urban_mask_fp,
        major_rivers_mask_fp=major_rivers_mask_fp,
        handle_large_polygons=handle_large_polygons,
        pp_test_threshold=pp_test_threshold,
    )

    filtered_polygons_fp = os.path.join(output_directory, "filtered_polygons.parquet")
    filtered_polygons.to_parquet(filtered_polygons_fp)
    _log.info(f"Filtered waterbody polygons written to {filtered_polygons_fp}")
