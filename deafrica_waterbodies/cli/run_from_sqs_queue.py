import logging
import os

import boto3
import click
import datacube
import fsspec
from rasterio.errors import RasterioIOError

from deafrica_waterbodies.cli.logs import logging_setup
from deafrica_waterbodies.io import check_dir_exists, check_file_exists, check_if_s3_uri
from deafrica_waterbodies.make_polygons import check_wetness_thresholds, get_polygons_from_dataset
from deafrica_waterbodies.queues import (
    delete_batch_with_retry,
    get_queue_url,
    move_to_dead_letter_queue,
    receive_a_message,
)


@click.command("run-from-sqs-queue", no_args_is_help=True)
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
    "--dataset-ids-queue",
    type=str,
    help="Name of the SQS queue to read dataset IDs from.",
)
@click.option(
    "--visibility-timeout",
    default=18 * 60,
    help="The duration in seconds that a received SQS msg is invisible.",
)
@click.option(
    "--max-retries",
    default=10,
    help="Maximum number of times to retry sending/receiving messages to/from a SQS queue.",
)
@click.option(
    "--overwrite/--no-overwrite",
    default=False,
    help="Rerun scenes that have already been processed.",
)
def run_from_sqs_queue(
    verbose,
    primary_threshold,
    secondary_threshold,
    minimum_valid_observations,
    output_directory,
    dataset_ids_queue,
    visibility_timeout,
    max_retries,
    overwrite,
):
    """
    Generate waterbody polygons from WOfS All Time Summary scenes whose ids are
    in a SQS queue.
    """
    # Set up logger.
    logging_setup(verbose=verbose)
    _log = logging.getLogger(__name__)

    # Support pathlib paths.
    output_directory = str(output_directory)

    dask_chunks = {"x": 3200, "y": 3200, "time": 1}
    resolution = (-30, 30)
    output_crs = "EPSG:6933"

    # Instanstiate the filesystem to use.
    if check_if_s3_uri(output_directory):
        fs = fsspec.filesystem("s3")
    else:
        fs = fsspec.filesystem("file")

    # Directory to write generated waterbody polygons to.
    polygons_from_thresholds_dir = os.path.join(output_directory, "polygons_from_thresholds")

    # Check if the directory exists. If it does not, create it.
    if not check_dir_exists(polygons_from_thresholds_dir):
        fs.mkdirs(polygons_from_thresholds_dir, exist_ok=True)
        _log.info(f"Created directory {polygons_from_thresholds_dir}")

    # Check if the wetness thresholds have been set correctly.
    minimum_wet_thresholds = [secondary_threshold, primary_threshold]
    _log.info(check_wetness_thresholds(minimum_wet_thresholds))

    # Connect to the datacube.
    dc = datacube.Datacube(app="GenerateWaterbodyPolygons")

    # Create the service client.
    sqs_client = boto3.client("sqs")

    dataset_ids_queue_url = get_queue_url(queue_name=dataset_ids_queue, sqs_client=sqs_client)
    # Get the dead-letter queue.
    dead_letter_queue_name = f"{dataset_ids_queue}-deadletter"
    dead_letter_queue_url = get_queue_url(queue_name=dead_letter_queue_name, sqs_client=sqs_client)

    retries = 0
    while retries <= max_retries:
        # Retrieve a single message from the dataset_ids_queue.
        message = receive_a_message(
            queue_url=dataset_ids_queue_url,
            max_retries=max_retries,
            visibility_timeout=visibility_timeout,
            sqs_client=sqs_client,
        )
        if message is None:
            retries += 1
        else:
            retries = 0  # reset the count

            # Process the ID.
            dataset_id = message["Body"]
            _log.info(f"Read dataset id {dataset_id} from queue {dataset_ids_queue_url}")

            entry_to_delete = [
                {"Id": message["MessageId"], "ReceiptHandle": message["ReceiptHandle"]}
            ]

            # Produce the primary and secondary threshold polygons.
            success_flag = True

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
                try:
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

                except KeyError as keyerr:
                    _log.exception(f"Found {dataset_id} has KeyError: {str(keyerr)}")
                    _log.error(f"Moving {dataset_id} to deadletter queue {dead_letter_queue_url}")
                    move_to_dead_letter_queue(
                        dead_letter_queue_url=dead_letter_queue_url,
                        message_body=dataset_id,
                        sqs_client=sqs_client,
                    )
                    success_flag = False
                except TypeError as typeerr:
                    _log.exception(f"Found {dataset_id} has TypeError: {str(typeerr)}")
                    _log.error(f"Moving {dataset_id} to deadletter queue {dead_letter_queue_url}")
                    move_to_dead_letter_queue(
                        dead_letter_queue_url=dead_letter_queue_url,
                        message_body=dataset_id,
                        sqs_client=sqs_client,
                    )
                    success_flag = False
                except RasterioIOError as ioerror:
                    _log.exception(f"Found {dataset_id} has RasterioIOError: {str(ioerror)}")
                    _log.error(f"Moving {dataset_id} to deadletter queue {dead_letter_queue_url}")
                    move_to_dead_letter_queue(
                        dead_letter_queue_url=dead_letter_queue_url,
                        message_body=dataset_id,
                        sqs_client=sqs_client,
                    )
                    success_flag = False
            else:
                _log.info(
                    f"{primary_threshold_polygons_fp} and {secondary_threshold_polygons_fp} already exist, skipping"
                )

            # Delete datased id from queue.
            if success_flag:
                _log.info(f"Successful, deleting {dataset_id}")
                (
                    successfully_deleted,
                    failed_to_delete,
                ) = delete_batch_with_retry(
                    queue_url=dataset_ids_queue_url,
                    entries=entry_to_delete,
                    max_retries=max_retries,
                    sqs_client=sqs_client,
                )
                if failed_to_delete:
                    _log.error(
                        f"Failed to delete dataset id {dataset_id} from queue {dataset_ids_queue_url}"
                    )
                    raise RuntimeError(f"Failed to delete dataset id: {dataset_id}")
                else:
                    _log.info(f"Deleted dataset id {dataset_id} from queue")

            else:
                _log.info(
                    f"Not successful, moved {dataset_id} to dead letter queue {dead_letter_queue_url}"
                )
