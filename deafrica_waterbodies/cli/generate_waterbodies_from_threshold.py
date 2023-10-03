import logging
import os

import click
import datacube
import fsspec

import deafrica_waterbodies.io
import deafrica_waterbodies.make_polygons
from deafrica_waterbodies.cli.logs import logging_setup


@click.command("generate-waterbodies-from-thresholds", no_args_is_help=True)
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
def generate_waterbodies_from_thresholds(
    verbose,
    primary_threshold,
    secondary_threshold,
    minimum_valid_observations,
    output_directory,
):
    # Set up logger.
    logging_setup(verbose=verbose)
    _log = logging.getLogger(__name__)

    # Support pathlib paths.
    output_directory = str(output_directory)

    dask_chunks = {"x": 3200, "y": 3200, "time": 1}
    resolution = (-30, 30)
    output_crs = "EPSG:6933"

    # Instanstiate the filesystem to use.
    if deafrica_waterbodies.io.check_if_s3_uri(output_directory):
        fs = fsspec.filesystem("s3")
    else:
        fs = fsspec.filesystem("file")

    # Check if the dataset ids text file exists in the output directory.
    dataset_ids_txt = os.path.join(output_directory, "dataset_ids.txt")

    if not deafrica_waterbodies.io.check_file_exists(dataset_ids_txt):
        _log.error(f"Could not find text file {dataset_ids_txt}!")
        raise FileNotFoundError(f"Could not find text file {dataset_ids_txt}!")

    # Read the dataset ids from the text file.
    with fs.open(dataset_ids_txt, "r") as file:
        lines = file.readlines()
        dataset_ids = [line.strip() for line in lines]

    # Directory to write generated waterbody polygons to.
    polygons_from_thresholds_dir = os.path.join(output_directory, "polygons_from_thresholds")

    # Check if the directory exists. If it does not, create it.
    if not deafrica_waterbodies.io.check_dir_exists(polygons_from_thresholds_dir):
        fs.mkdirs(polygons_from_thresholds_dir, exist_ok=True)
        _log.info(f"Created directory {polygons_from_thresholds_dir}")

    # Check if the wetness thresholds have been set correctly.
    minimum_wet_thresholds = [secondary_threshold, primary_threshold]
    _log.info(deafrica_waterbodies.make_polygons.check_wetness_thresholds(minimum_wet_thresholds))

    # Connect to the datacube.
    dc = datacube.Datacube(app="GenerateWaterbodyPolygons")

    # For each dataset id, threshold the scene to generate the primary and secondary threshold
    # waterbody polygons.
    for dataset_id in dataset_ids:
        (
            primary_threshold_polygons,
            secondary_threshold_polygons,
        ) = deafrica_waterbodies.make_polygons.get_polygons_using_thresholds(
            dataset_id=dataset_id,
            dask_chunks=dask_chunks,
            resolution=resolution,
            output_crs=output_crs,
            min_valid_observations=minimum_valid_observations,
            primary_threshold=primary_threshold,
            secondary_threshold=secondary_threshold,
            dc=dc,
        )
        # Write the polygons to parquet.
        primary_threshold_polygons_fp = os.path.join(
            polygons_from_thresholds_dir, f"{dataset_id}_primary_threshold_polygons.parquet"
        )
        secondary_threshold_polygons_fp = os.path.join(
            polygons_from_thresholds_dir, f"{dataset_id}_secondary_threshold_polygons.parquet"
        )

        primary_threshold_polygons.to_parquet(primary_threshold_polygons_fp)
        secondary_threshold_polygons.to_parquet(secondary_threshold_polygons_fp)
