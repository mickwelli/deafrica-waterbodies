import logging
import os
import urllib

import click
import fsspec
import geopandas as gpd

from deafrica_waterbodies.cli.logs import logging_setup
from deafrica_waterbodies.io import check_if_s3_uri
from deafrica_waterbodies.make_polygons import get_datasets_ids


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
    "--output-file-path",
    type=click.Path(),
    help="File URI or S3 URI of the text file to write the dataset ids to.",
)
def get_dataset_ids(
    verbose,
    aoi_vector_file,
    num_workers,
    output_file_path,
):
    # Set up logger.
    logging_setup(verbose=verbose)
    _log = logging.getLogger(__name__)

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
    dataset_ids = get_datasets_ids(aoi_gdf=aoi_gdf, num_workers=num_workers)

    # Write the dataset ids to a text file.
    if check_if_s3_uri(output_file_path):
        fs = fsspec.filesystem("s3")
    else:
        fs = fsspec.filesystem("file")
        parsed_output_fp = urllib.parse.urlparse(output_file_path).path
        absolute_output_fp = os.path.abspath(parsed_output_fp)
        path_head, path_tail = os.path.split(absolute_output_fp)
        if path_head:
            if not fs.exists(path_head):
                fs.mkdirs(path_head, exist_ok=True)
                _log.info(f"Local folder {path_head} created.")

    with fs.open(output_file_path, "w") as file:
        for dataset_id in dataset_ids:
            file.write(f"{dataset_id}\n")

    _log.info(f"Dataset IDs written to: {output_file_path}.")
