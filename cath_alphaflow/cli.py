import logging
import click

from .commands import create_dataset_uniprot_ids
from .commands import create_dataset_cath_files
from .commands import optimise_domain_boundaries

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


cli.add_command(create_dataset_uniprot_ids.create_dataset_uniprot_ids)
cli.add_command(create_dataset_cath_files.create_dataset_cath_files)
cli.add_command(optimise_domain_boundaries.optimise_domain_boundaries)
