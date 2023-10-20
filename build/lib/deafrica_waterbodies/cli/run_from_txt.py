import logging
import os

import click
import datacube
import fsspec

from deafrica_waterbodies.cli.logs import logging_setup
from deafrica_waterbodies.io import check_dir_exists, check_file_exists, check_if_s3_uri
from deafrica_waterbodies.make_polygons import check_wetness_thresholds, get_polygons_from_dataset


@click.command("run-from-txt", no_args_is_help=True)
@click.option("-v", "--verbose", count=True)
@click.option(
    "--primary-threshold",
    default=0.1,
    type=float,
    help="Threshold to define the location of the waterbody polygons.",
    show_default=True,
)
@click.option(
    "--secondary-threshold",
    default=0.05,
    type=float,
    help="Threshold to define the extent/shape of the waterbody polygons.",
    show_default=True,
)
@click.option(
    "--minimum-valid-observations",
    default=128,
    type=int,
    help="Minimum number of observations for a pixel to be considered valid.",
    show_default=True,
)
@click.option(
    "--output-directory",
    type=click.Path(),
    help="Directory to write the waterbody polygons to.",
)
@click.option(
    "--dataset-ids-text-file",
    type=click.Path(),
    required=True,
    help="Path to dataset ids text file.",
)
@click.option(
    "--overwrite/--no-overwrite",
    default=False,
    help="Rerun scenes that have already been processed.",
)
def run_from_txt(
    verbose,
    primary_threshold,
    secondary_threshold,
    minimum_valid_observations,
    output_directory,
    dataset_ids_text_file,
    overwrite,
):
    """
    Generate waterbody polygons from WOfS All Time Summary scenes whose ids are
    in a text file.
    """
    # Set up logger.
    logging_setup(verbose=verbose)
    _log = logging.getLogger(__name__)

    # Support pathlib paths.
    output_directory = str(output_directory)
    dataset_ids_text_file = str(dataset_ids_text_file)

    dask_chunks = {"x": 3200, "y": 3200, "time": 1}
    resolution = (-30, 30)
    output_crs = "EPSG:6933"

    # Read the dataset ids from the text file.
    if not check_file_exists(dataset_ids_text_file):
        _log.error(f"Could not find text file {dataset_ids_text_file}!")
        raise FileNotFoundError(f"Could not find text file {dataset_ids_text_file}!")
    else:
        if check_if_s3_uri(dataset_ids_text_file):
            fs = fsspec.filesystem("s3")
        else:
            fs = fsspec.filesystem("file")
        with fs.open(dataset_ids_text_file, "r") as file:
            lines = file.readlines()
            dataset_ids = [line.strip() for line in lines]

    # Directory to write generated waterbody polygons to.
    polygons_from_thresholds_dir = os.path.join(output_directory, "polygons_from_thresholds")

    # Set the filesystem to use.
    if check_if_s3_uri(polygons_from_thresholds_dir):
        fs = fsspec.filesystem("s3")
    else:
        fs = fsspec.filesystem("file")

    # Check if the directory exists. If it does not, create it.
    if not check_dir_exists(polygons_from_thresholds_dir):
        fs.mkdirs(polygons_from_thresholds_dir, exist_ok=True)
        _log.info(f"Created directory {polygons_from_thresholds_dir}")

    # Check if the wetness thresholds have been set correctly.
    minimum_wet_thresholds = [secondary_threshold, primary_threshold]
    _log.info(check_wetness_thresholds(minimum_wet_thresholds))

    # Connect to the datacube.
    dc = datacube.Datacube(app="GenerateWaterbodyPolygons")

    # For each dataset id, threshold the scene to generate the primary and secondary threshold
    # waterbody polygons.
    for dataset_id in dataset_ids:
        primary_threshold_polygons_fp = os.path.join(
            polygons_from_thresholds_dir, f"{dataset_id}_primary_threshold_polygons.parquet"
        )
        secondary_threshold_polygons_fp = os.path.join(
            polygons_from_thresholds_dir, f"{dataset_id}_secondary_threshold_polygons.parquet"
        )

        if not overwrite:
            _log.info(
                f"Checking existence of {primary_threshold_polygons_fp} and {secondary_threshold_polygons_fp}"
            )
            exists = check_file_exists(primary_threshold_polygons_fp) and check_file_exists(
                secondary_threshold_polygons_fp
            )

        if overwrite or not exists:
            (
                primary_threshold_polygons,
                secondary_threshold_polygons,
            ) = get_polygons_from_dataset(
                dataset_id=dataset_id,
                dask_chunks=dask_chunks,
                resolution=resolution,
                output_crs=output_crs,
                min_valid_observations=minimum_valid_observations,
                primary_threshold=primary_threshold,
                secondary_threshold=secondary_threshold,
                dc=dc,
            )
            # Write the polygons to parquet files.
            primary_threshold_polygons.to_parquet(primary_threshold_polygons_fp)
            secondary_threshold_polygons.to_parquet(secondary_threshold_polygons_fp)
