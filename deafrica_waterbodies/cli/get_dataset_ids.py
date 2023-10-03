import logging
import os

import click
import fsspec
import geopandas as gpd

import deafrica_waterbodies.io
import deafrica_waterbodies.make_polygons
from deafrica_waterbodies.cli.logs import logging_setup


@click.command("get-dataset-ids", no_args_is_help=True)
@click.option("-v", "--verbose", count=True)
@click.option(
    "--aoi-vector-file",
    type=str,
    default=None,
    help="Path to the vector file defining the area of interest.",
)
@click.option(
    "--num-workers", type=int, help="Number of worker processes to use when filtering datasets."
)
@click.option(
    "--output-directory",
    type=click.Path(),
    help="Directory to write the dataset ids text file to.",
)
def get_dataset_ids(
    verbose,
    aoi_vector_file,
    num_workers,
    output_directory,
):
    """
    Get the dataset ids of the WOfS All Time summary datasets/scenes to generate
    the waterbody polygons for.
    """
    # Set up logger.
    logging_setup(verbose=verbose)
    _log = logging.getLogger(__name__)

    # Support pathlib Paths.
    aoi_vector_file = str(aoi_vector_file)
    output_directory = str(output_directory)

    # Load the area of interest as a GeoDataFrame.
    if aoi_vector_file is not None:
        try:
            aoi_gdf = gpd.read_file(aoi_vector_file)
        except Exception as error:
            _log.exception(f"Could not read the file {aoi_vector_file}")
            raise error
    else:
        aoi_gdf = None

    # Get the WOfS All Time Summary scene ids for the scenes whose extent
    # intersects with the area of interest.
    dataset_ids = deafrica_waterbodies.make_polygons.get_datasets_ids(
        aoi_gdf=aoi_gdf, num_workers=num_workers
    )

    # Instanstiate the filesystem to use.
    if deafrica_waterbodies.io.check_if_s3_uri(output_directory):
        fs = fsspec.filesystem("s3")
    else:
        fs = fsspec.filesystem("file")

    # Check if the output directory exists. If it does not, create it.
    if not deafrica_waterbodies.io.check_dir_exists(output_directory):
        fs.mkdirs(output_directory, exist_ok=True)
        _log.info(f"Created directory {output_directory}")

    # Write the dataset ids to a text file.
    output_file_path = os.path.join(output_directory, "dataset_ids.txt")

    with fs.open(output_file_path, "w") as file:
        for dataset_id in dataset_ids:
            file.write(f"{dataset_id}\n")

    _log.info(f"Dataset IDs written to: {output_file_path}.")
