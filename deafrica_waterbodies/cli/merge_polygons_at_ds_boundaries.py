import logging
import os

import click
import geopandas as gpd
import pandas as pd

from deafrica_waterbodies.cli.logs import logging_setup
from deafrica_waterbodies.io import find_parquet_files
from deafrica_waterbodies.make_polygons import merge_polygons_at_dataset_boundaries


@click.command("merge-polygons-at-ds-boundaries", no_args_is_help=True)
@click.option("-v", "--verbose", count=True)
@click.option(
    "--output-directory",
    type=click.Path(),
    help="Directory containing the waterbody polygons.",
)
def merge_polygons_at_ds_boundaries(verbose, output_directory):
    """
    Merge polygons at dataset boundaries.
    """
    # Set up logger.
    logging_setup(verbose=verbose)
    _log = logging.getLogger(__name__)

    # Support pathlib paths.
    output_directory = str(output_directory)

    # Directory containing the water body polygons generated from
    # thresholding WOfS All time summary datasets.
    polygons_from_thresholds_dir = os.path.join(output_directory, "polygons_from_thresholds")

    # Find all parquet files for the primary threshold.
    primary_threshold_polygons_paths = find_parquet_files(
        path=polygons_from_thresholds_dir, pattern=".*primary.*"
    )
    _log.info(f"Found {len(primary_threshold_polygons_paths)} primary threshold polygons.")

    # Load all the primary threshold polygons into a single GeoDataFrame.
    _log.info("Loading the primary threshold polygons parquet files..")
    primary_threshold_polygons_list = []
    for path in primary_threshold_polygons_paths:
        gdf = gpd.read_parquet(path)
        primary_threshold_polygons_list.append(gdf)

    primary_threshold_polygons = pd.concat(primary_threshold_polygons_list, ignore_index=True)
    _log.info(f"Found {len(primary_threshold_polygons)} primary threshold polygons.")

    _log.info("Merging primary threshold waterbody polygons located at dataset/scene boundaries...")
    primary_threshold_polygons_merged = merge_polygons_at_dataset_boundaries(
        primary_threshold_polygons
    )
    _log.info(f"Primary threshold polygons count {len(primary_threshold_polygons_merged)}.")

    _log.info("Writing primary threshold polygons merged at dataset boundaries to disk..")
    primary_threshold_polygons_output_fp = os.path.join(
        output_directory, "primary_threshold_polygons_merged_at_ds_boundaries.parquet"
    )

    primary_threshold_polygons_merged.to_parquet(primary_threshold_polygons_output_fp)
    _log.info(f"Polygons written to {primary_threshold_polygons_output_fp}")

    # Find all parquet files for the secondary threshold.
    secondary_threshold_polygons_paths = find_parquet_files(
        path=polygons_from_thresholds_dir, pattern=".*secondary.*"
    )
    _log.info(
        f"Found {len(secondary_threshold_polygons_paths)} parquet files for the secondary threshold polygons."
    )

    # Load all the secondary threshold polygons into a single GeoDataFrame.
    _log.info("Loading the secondary threshold polygons parquet files...")
    secondary_threshold_polygons_list = []
    for path in secondary_threshold_polygons_paths:
        gdf = gpd.read_parquet(path)
        secondary_threshold_polygons_list.append(gdf)

    secondary_threshold_polygons = pd.concat(secondary_threshold_polygons_list, ignore_index=True)
    _log.info(f"Found {len(secondary_threshold_polygons)} secondary threshold polygons.")

    _log.info(
        "Merging secondary threshold waterbody polygons located at dataset/scene boundaries..."
    )
    secondary_threshold_polygons_merged = merge_polygons_at_dataset_boundaries(
        secondary_threshold_polygons
    )
    _log.info(f"Secondary threshold polygons count {len(secondary_threshold_polygons_merged)}.")

    _log.info("Writing secondary threshold polygons merged at dataset boundaries to disk..")
    secondary_threshold_polygons_output_fp = os.path.join(
        output_directory, "secondary_threshold_polygons_merged_at_ds_boundaries.parquet"
    )

    secondary_threshold_polygons_merged.to_parquet(secondary_threshold_polygons_output_fp)

    _log.info(f"Polygons written to {secondary_threshold_polygons_output_fp}")
