import click

import deafrica_waterbodies.__version__
from deafrica_waterbodies.cli.get_dataset_ids import get_dataset_ids
from deafrica_waterbodies.cli.merge_polygons_at_ds_boundaries import merge_polygons_at_ds_boundaries
from deafrica_waterbodies.cli.push_to_sqs_queue import push_to_sqs_queue
from deafrica_waterbodies.cli.run_from_sqs_queue import run_from_sqs_queue
from deafrica_waterbodies.cli.run_from_txt import run_from_txt


@click.version_option(package_name="deafrica_waterbodies", version=deafrica_waterbodies.__version__)
@click.group(help="Run deafrica-waterbodies.")
def main():
    pass


main.add_command(get_dataset_ids)
main.add_command(push_to_sqs_queue)
main.add_command(run_from_sqs_queue)
main.add_command(run_from_txt)
main.add_command(merge_polygons_at_ds_boundaries)
