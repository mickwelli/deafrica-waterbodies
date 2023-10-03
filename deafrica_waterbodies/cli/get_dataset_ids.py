import logging

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
    "--dataset-ids-text-file",
    type=click.Path(),
    help="File URI or S3 URI of the text file to write the dataset ids to.",
)
def get_dataset_ids(
    verbose,
    aoi_vector_file,
    num_workers,
    dataset_ids_text_file,
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
    dataset_ids_text_file = str(dataset_ids_text_file)

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
    if deafrica_waterbodies.io.check_if_s3_uri(dataset_ids_text_file):
        fs = fsspec.filesystem("s3")
    else:
        fs = fsspec.filesystem("file")

    # Write the dataset ids to the text file.
    with fs.open(dataset_ids_text_file, "w") as file:
        for dataset_id in dataset_ids:
            file.write(f"{dataset_id}\n")

    _log.info(f"Dataset IDs written to: {dataset_ids_text_file}.")
