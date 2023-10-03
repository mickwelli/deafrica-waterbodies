import logging
import os

import click
import geopandas as gpd
import pandas as pd

import deafrica_waterbodies.io
import deafrica_waterbodies.make_polygons
from deafrica_waterbodies.cli.logs import logging_setup


@click.command("merge-polygons-at-ds-boundaries", no_args_is_help=True)
@click.option("-v", "--verbose", count=True)
@click.option(
    "--output-directory",
    type=click.Path(),
    help="Directory containing the waterbody polygons.",
)
def merge_polygons_at__ds_boundaries(verbose, output_directory):
    # Set up logger.
    logging_setup(verbose=verbose)
    _log = logging.getLogger(__name__)

    # Support pathlib paths.
    output_directory = str(output_directory)

    polygons_from_thresholds_dir = os.path.join(output_directory, "polygons_from_thresholds")

    # Find all parquet files for the primary threshold.
    primary_threshold_polygons_paths = deafrica_waterbodies.io.find_parquet_files(
        path=polygons_from_thresholds_dir, pattern=".*primary.*"
    )
    _log.info(f"Found {len(primary_threshold_polygons_paths)} primary threshold polygons.")

    # Load all the primary threshold polygons into a single GeoDataFrame.
    primary_threshold_polygons_list = []
    for path in primary_threshold_polygons_paths:
        gdf = gpd.read_parquet(path)
        primary_threshold_polygons_list.append(gdf)

    primary_threshold_polygons = pd.concat(primary_threshold_polygons_list, ignore_index=True)

    # Find all parquet files for the secondary threshold.
    secondary_threshold_polygons_paths = deafrica_waterbodies.io.find_parquet_files(
        path=polygons_from_thresholds_dir, pattern=".*secondary.*"
    )
    _log.info(f"Found {len(secondary_threshold_polygons_paths)} secondary threshold polygons.")

    # Load all the secondary threshold polygons into a single GeoDataFrame.
    secondary_threshold_polygons_list = []
    for path in secondary_threshold_polygons_paths:
        gdf = gpd.read_parquet(path)
        secondary_threshold_polygons_list.append(gdf)

    secondary_threshold_polygons = pd.concat(secondary_threshold_polygons_list, ignore_index=True)

    # Merge polygons at WOfS All Time Summary scene/dataset boundaries.
    _log.info("Merging waterbody polygons located at dataset/scene boundaries.")
    primary_threshold_polygons_merged = (
        deafrica_waterbodies.make_polygons.merge_polygons_at_dataset_boundaries(
            primary_threshold_polygons
        )
    )
    secondary_threshold_polygons_merged = (
        deafrica_waterbodies.make_polygons.merge_polygons_at_dataset_boundaries(
            secondary_threshold_polygons
        )
    )

    # Write the merged polygons to disk.
    primary_threshold_polygons_output_fp = os.path.join(
        output_directory, "primary_threshold_polygons_merged_at_ds_boundaries.parquet"
    )
    secondary_threshold_polygons_output_fp = os.path.join(
        output_directory, "secondary_threshold_polygons_merged_at_ds_boundaries.parquet"
    )

    primary_threshold_polygons_merged.to_parquet(primary_threshold_polygons_output_fp)
    secondary_threshold_polygons_merged.to_parquet(secondary_threshold_polygons_output_fp)

    _log.info(
        f"Polygons written to {primary_threshold_polygons_output_fp} and {secondary_threshold_polygons_output_fp}"
    )
