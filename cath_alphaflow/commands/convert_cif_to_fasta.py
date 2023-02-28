import logging
from pathlib import Path
import click
from cath_alphaflow.io_utils import yield_first_col
from cath_alphaflow.constants import DEFAULT_CIF_SUFFIX
from cath_alphaflow.settings import get_default_settings
from cath_alphaflow.seq_utils import cif_to_fasta, combine_fasta_files, write_fasta_file

config = get_default_settings()

LOG = logging.getLogger()


@click.command()
@click.option(
    "--cif_in_dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, resolve_path=True),
    required=True,
    help="Input: directory of CIF files",
)
@click.option(
    "--id_file",
    type=click.File("rt"),
    required=True,
    help="Input: CSV file containing list of ids to convert from CIF to FASTA",
)
@click.option(
    "--cif_suffix",
    type=str,
    default=DEFAULT_CIF_SUFFIX,
    help=f"Input: optional suffix to add to id when looking for cif file (default: {DEFAULT_CIF_SUFFIX})",
)
@click.option(
    "--fasta_out_dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, resolve_path=True),
    required=True,
    help="Output: FASTA Output Folder",
)
@click.option(
    "--combined_fasta_file",
    type=click.File("wt"),
    required=True,
    default="merged.fasta",
    help="Output multiFASTA file containing all FASTA sequences",
)
def convert_cif_to_fasta(
    cif_in_dir, id_file, cif_suffix, fasta_out_dir, combined_fasta_file
):
    "Convert CIF to FASTA"

    for file_stub in yield_first_col(id_file):
        cif_path = Path(cif_in_dir) / (file_stub + cif_suffix)
        header, sequence = cif_to_fasta(cif_path, fasta_out_dir)
        write_fasta_file(header=header, sequence=sequence, fasta_out_dir=fasta_out_dir)
    combine_fasta_files(fasta_in_dir=fasta_out_dir, fasta_out_file=combined_fasta_file)

    click.echo("DONE")
