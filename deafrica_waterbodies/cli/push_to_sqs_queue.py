import logging
import os

import boto3
import click

import deafrica_waterbodies.queues
from deafrica_waterbodies.cli.logs import logging_setup


@click.command("push-to-sqs-queue", no_args_is_help=True)
@click.option(
    "--output-directory",
    type=click.Path(),
    required=True,
    help="Directory containing the dataset ids text file.",
)
@click.option(
    "--dataset-ids-queue", required=True, help="Name of the queue to push the dataset ids to."
)
@click.option(
    "--max-retries",
    default=10,
    help="Maximum number of times to retry sending/receiving messages to/from a SQS queue.",
)
@click.option("-v", "--verbose", count=True)
def push_to_sqs_queue(output_directory, dataset_ids_queue, max_retries, verbose):
    """
    Push dataset ids from the lines of a text file to a SQS queue.
    """
    logging_setup(verbose=verbose)
    _log = logging.getLogger(__name__)  # noqa F841

    # Create an sqs client.
    sqs_client = boto3.client("sqs")

    # Dataset ids text file.
    text_file_path = os.path.join(output_directory, "dataset_ids.txt")

    failed_to_push = deafrica_waterbodies.queues.push_dataset_ids_to_queue_from_txt(
        text_file_path=text_file_path,
        queue_name=dataset_ids_queue,
        max_retries=max_retries,
        sqs_client=sqs_client,
    )

    if failed_to_push:
        # Push the failed dataset ids to the deadletter queue.
        dead_letter_queue_name = f"{dataset_ids_queue}-deadletter"
        dead_letter_queue_url = deafrica_waterbodies.queues.get_queue_url(
            queue_name=dead_letter_queue_name, sqs_client=sqs_client
        )

        for idx in failed_to_push:
            deafrica_waterbodies.queues.move_to_dead_letter_queue(
                dead_letter_queue_url=dead_letter_queue_url,
                message_body=idx,
                max_retries=max_retries,
                sqs_client=sqs_client,
            )
