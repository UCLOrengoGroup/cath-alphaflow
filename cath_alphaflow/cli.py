import logging
import click

from .settings import get_default_settings
from .commands import create_dataset_uniprot_ids
from .commands import create_dataset_cath_files
from .commands import optimise_domain_boundaries
from .commands import filter_domains_by_sse

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
    for key in dir(settings):
        if key.startswith("__"):
            continue
        val = getattr(settings, key)
        if "PASSWORD" in key:
            val = "******"
        click.echo(f"{key} {val}")


cli.add_command(dump_config)
cli.add_command(create_dataset_uniprot_ids.create_dataset_uniprot_ids)
cli.add_command(create_dataset_cath_files.create_dataset_cath_files)
cli.add_command(optimise_domain_boundaries.optimise_domain_boundaries)
cli.add_command(filter_domains_by_sse.filter_domains_by_sse)
