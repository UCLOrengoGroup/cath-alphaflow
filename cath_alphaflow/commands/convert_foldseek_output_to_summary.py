import logging
from nis import match
from pathlib import Path
import click
from cath_alphaflow.io_utils import yield_first_col
from cath_alphaflow.settings import get_default_settings
from cath_alphaflow.constants import DEFAULT_FS_BITS_CUTOFF, DEFAULT_FS_OVERLAP

config = get_default_settings()

LOG = logging.getLogger()


@click.command()
@click.option(
    "--id_file",
    type=click.File("rt"),
    required=True,
    help="Input: CSV file containing list of ids to convert from CIF to DSSP",
)
@click.option(
    "--fs_input_file",
    type=click.Path(exists=True, file_okay=True, resolve_path=True),
    default="fs_query_results.m8",
    help=f"Foldseek tabular output as input",
)
@click.option(
    "--fs_results",
    type=click.Path(resolve_path=True),
    default="fs_hits.tsv",
    help=f"Foldseek hits file",
)
def convert_foldseek_output_to_summary(id_file, fs_input_file, fs_results):
    "Convert Foldseek tabular output to summary of best hits"
    seen_ids = set()
    with open(fs_input_file, "rt") as fs_fh:
        for line in fs_fh:
            (
                query,
                target,
                qstart,
                qend,
                qlen,
                tstart,
                tend,
                tlen,
                qcov,
                tcov,
                bits,
                evalue,
            ) = line.split()
            if query.endswith(".cif"):
                query = query[:-4]
            if query in seen_ids:
                continue
            if (
                float(tcov) >= DEFAULT_FS_OVERLAP
                and int(bits) >= DEFAULT_FS_BITS_CUTOFF
            ):
                seen_ids.add(query)

    for file_stub in yield_first_col(id_file):
        click.echo(file_stub)
