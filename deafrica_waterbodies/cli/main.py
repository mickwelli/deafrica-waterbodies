import click

import deafrica_waterbodies.__version__

from .continental_run import waterbodies_continental_run
from .from_bbox import waterbodies_from_bbox
from .from_vector_file import waterbodies_from_vector_file


@click.version_option(package_name="deafrica_waterbodies", version=deafrica_waterbodies.__version__)
@click.group(help="Run deafrica-waterbodies.")
def main():
    pass


main.add_command(waterbodies_continental_run)
main.add_command(waterbodies_from_bbox)
main.add_command(waterbodies_from_vector_file)
