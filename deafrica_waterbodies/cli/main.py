import click

import deafrica_waterbodies.__version__
from deafrica_waterbodies.cli.generate_timeseries import generate_timeseries
from deafrica_waterbodies.cli.generate_waterbodies import generate_waterbodies


@click.version_option(package_name="deafrica_waterbodies", version=deafrica_waterbodies.__version__)
@click.group(help="Run deafrica-waterbodies.")
def main():
    pass


main.add_command(generate_waterbodies)
main.add_command(generate_timeseries)
