import logging
import click


from .settings import get_default_settings
from .commands import create_dataset_uniprot_ids
from .commands import create_dataset_cath_files
from .commands import optimise_domain_boundaries
from .commands import convert_dssp_to_sse_summary
from .commands import convert_cif_to_dssp
from .commands import extract_plddt_and_lur
from .commands import convert_cif_to_foldseek_db
from .commands import run_foldseek
from .commands import convert_foldseek_output_to_summary
from .commands import chop_domain
from .commands import create_md5
from .commands import convert_cif_to_fasta
from .commands import load_mongo
from .commands import measure_globularity
from .commands import measure_globularity_af_cif

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
cli.add_command(create_dataset_cath_files.create_cath_dataset_from_db)
cli.add_command(create_dataset_cath_files.create_cath_dataset_from_files)
cli.add_command(optimise_domain_boundaries.optimise_domain_boundaries)
cli.add_command(convert_dssp_to_sse_summary.convert_dssp_to_sse_summary)
cli.add_command(convert_cif_to_dssp.convert_cif_to_dssp)
cli.add_command(extract_plddt_and_lur.convert_cif_to_plddt_summary)
cli.add_command(convert_cif_to_foldseek_db.convert_cif_to_foldseek_db)
cli.add_command(run_foldseek.run_foldseek)
cli.add_command(convert_foldseek_output_to_summary.convert_foldseek_output_to_summary)
cli.add_command(chop_domain.chop_cif_command)
cli.add_command(create_md5.create_md5)
cli.add_command(convert_cif_to_fasta.convert_cif_to_fasta)
cli.add_command(load_mongo.load_af_from_archive)
cli.add_command(measure_globularity.measure_globularity)
cli.add_command(measure_globularity_af_cif.measure_globularity_af_cif)
