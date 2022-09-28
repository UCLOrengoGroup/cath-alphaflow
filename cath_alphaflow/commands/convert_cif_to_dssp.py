import logging
from pathlib import Path
import click
import subprocess
from cath_alphaflow.io_utils import yield_first_col

from cath_alphaflow.constants import DEFAULT_CIF_SUFFIX, DEFAULT_DSSP_SUFFIX
from cath_alphaflow.settings import get_default_settings

config = get_default_settings()

DSSP_BINARY_PATH = config.DSSP_BINARY_PATH
DSSP_PDB_DICT = config.DSSP_PDB_DICT


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
    help="Input: CSV file containing list of ids to convert from CIF to DSSP",
)
@click.option(
    "--cif_suffix",
    type=str,
    default=DEFAULT_CIF_SUFFIX,
    help=f"Input: optional suffix to add to id when looking for cif file (default: {DEFAULT_CIF_SUFFIX})",
)
@click.option(
    "--dssp_suffix",
    type=str,
    default=DEFAULT_DSSP_SUFFIX,
    help=f"Input: optional suffix to add to id when writing dssp files (default: {DEFAULT_DSSP_SUFFIX})",
)
@click.option(
    "--dssp_out_dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, resolve_path=True),
    required=True,
    help="Output: DSSP Output Folder",
)
def convert_cif_to_dssp(cif_dir, id_file, cif_suffix, dssp_suffix, dssp_out_dir):
    "Converts CIF to DSSP files"

    for file_stub in yield_first_col(id_file):
        cif_path = Path(cif_dir) / file_stub + cif_suffix
        dssp_path = Path(dssp_out_dir) / file_stub + dssp_suffix
        click.echo(f"Running DSSP: {cif_path} {dssp_path}")
        run_dssp(cif_path, dssp_path)

    click.echo("DONE")


def run_dssp(cif_path: Path, dssp_path: Path):

    if not cif_path.exists():
        msg = f"failed to locate CIF input file {cif_path}"
        LOG.error(msg)
        raise FileNotFoundError(msg)

    if not dssp_path.parent.exists():
        msg = f"failed to find output DSSP directory {dssp_path.parent}"
        LOG.error(msg)
        raise FileNotFoundError(msg)

    subprocess.call(
        [
            DSSP_BINARY_PATH,
            "--mmcif-dictionary",
            DSSP_PDB_DICT,
            "--output-format",
            "dssp",
            f"{cif_path}",
            f"{dssp_path}",
        ],
        stderr=subprocess.DEVNULL,
        check=True,
    )
