import logging

import boto3
import click

from deafrica_waterbodies.cli.logs import logging_setup
from deafrica_waterbodies.queues import (
    get_queue_url,
    move_to_dead_letter_queue,
    push_dataset_ids_to_queue_from_txt,
)


@click.command("push-to-sqs-queue", no_args_is_help=True)
@click.option("-v", "--verbose", count=True)
@click.option(
    "--dataset-ids-text-file",
    type=click.Path(),
    required=True,
    help="Path to dataset ids text file.",
)
@click.option(
    "--dataset-ids-queue", required=True, help="Name of the queue to push the dataset ids to."
)
@click.option(
    "--max-retries",
    default=10,
    help="Maximum number of times to retry sending/receiving messages to/from a SQS queue.",
)
def push_to_sqs_queue(verbose, dataset_ids_text_file, dataset_ids_queue, max_retries):
    """
    Push dataset ids from the lines of a text file to a SQS queue.
    """
    logging_setup(verbose=verbose)
    _log = logging.getLogger(__name__)  # noqa F841

    # Create an sqs client.
    sqs_client = boto3.client("sqs")

    # Support pathlib paths.
    dataset_ids_text_file = str(dataset_ids_text_file)

    failed_to_push = push_dataset_ids_to_queue_from_txt(
        text_file_path=dataset_ids_text_file,
        queue_name=dataset_ids_queue,
        max_retries=max_retries,
        sqs_client=sqs_client,
    )

    if failed_to_push:
        # Push the failed dataset ids to the deadletter queue.
        dead_letter_queue_name = f"{dataset_ids_queue}-deadletter"
        dead_letter_queue_url = get_queue_url(
            queue_name=dead_letter_queue_name, sqs_client=sqs_client
        )

        for idx in failed_to_push:
            move_to_dead_letter_queue(
                dead_letter_queue_url=dead_letter_queue_url,
                message_body=idx,
                max_retries=max_retries,
                sqs_client=sqs_client,
            )
