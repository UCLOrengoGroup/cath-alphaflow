import logging
import click

from .settings import get_default_settings
from .commands import create_dataset_uniprot_ids
from .commands import create_dataset_cath_files
from .commands import optimise_domain_boundaries
from .commands import convert_dssp_to_sse_summary
from .commands import convert_cif_to_dssp

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s"
)

LOG = logging.getLogger(__name__)


@click.group()
@click.version_option()
@click.option("--verbose", "-v", "verbosity", default=0, count=True)
@click.pass_context
def cli(ctx, verbosity):
    "Workflow tools to assign CATH structural domains to AlphaFold predictions"

    root_logger = logging.getLogger()
    log_level = root_logger.getEffectiveLevel() - (10 * verbosity)
    root_logger.setLevel(log_level)
    LOG.info(
        f"Starting logging... (level={logging.getLevelName(root_logger.getEffectiveLevel())})"
    )


@click.command()
def dump_config():
    """
    Dump the current settings
    """
    settings = get_default_settings()
    click.echo("Settings:")
    for key, val in settings.to_dict().items():
        click.echo(f"{key:25s} {val}")


cli.add_command(dump_config)
cli.add_command(create_dataset_uniprot_ids.create_dataset_uniprot_ids)
cli.add_command(create_dataset_cath_files.create_dataset_cath_files)
cli.add_command(optimise_domain_boundaries.optimise_domain_boundaries)
cli.add_command(convert_dssp_to_sse_summary.convert_dssp_to_sse_summary)
cli.add_command(convert_cif_to_dssp.convert_cif_to_dssp)
