import logging
from pathlib import Path
import os
import click
import subprocess
from cath_alphaflow.io_utils import yield_first_col

from cath_alphaflow.constants import (
    DEFAULT_FS_QUERYDB_SUFFIX,
    DEFAULT_CIF_SUFFIX,
)
from cath_alphaflow.settings import get_default_settings

config = get_default_settings()

FS_BINARY_PATH = config.FS_BINARY_PATH

LOG = logging.getLogger()


@click.command()
@click.option(
    "--cif_dir",
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
    "--fs_querydb_suffix",
    type=str,
    default=DEFAULT_FS_QUERYDB_SUFFIX,
    help=f"Output: optional suffix to add to Foldseek query database during creation (default: {DEFAULT_FS_QUERYDB_SUFFIX})",
)
@click.option(
    "--fs_querydb_dir",
    type=click.Path(file_okay=False, dir_okay=True, resolve_path=True),
    required=True,
    help=f"Output: directory to use for Foldseek query database during creation",
)
def convert_cif_to_foldseek_db(
    cif_dir, fs_querydb_dir, id_file, cif_suffix, fs_querydb_suffix
):
    "Create Foldseek query database from mmCIF folder"
    fs_querydb_path = Path(fs_querydb_dir)
    if not fs_querydb_path.exists():
        os.makedirs(fs_querydb_path)
    for file_stub in yield_first_col(id_file):
        cif_path = Path(cif_dir) / f"{file_stub}{cif_suffix}"
        click.echo(cif_path)
        if not cif_path.exists():
            msg = f"failed to locate CIF input file {cif_path}"
            LOG.error(msg)
            raise FileNotFoundError(msg)
        # Create symlinks to querydb_dir
        subprocess.call(
            ["ln", "-s", cif_path, f"{fs_querydb_path}/{cif_path.name}"],
        )
    subprocess.call(
        [
            FS_BINARY_PATH,
            "createdb",
            f"{fs_querydb_path}/",
            f"{fs_querydb_dir}{fs_querydb_suffix}",
        ],
        stderr=subprocess.DEVNULL,
    )
    click.echo("DONE")
    return Path.exists(Path(f"{fs_querydb_dir}{fs_querydb_suffix}"))
