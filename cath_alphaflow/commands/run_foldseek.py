from email.policy import default
import logging
from pathlib import Path
import os
import click
import subprocess

from cath_alphaflow.settings import get_default_settings

config = get_default_settings()

FS_BINARY_PATH = config.FS_BINARY_PATH
FS_DB_PATH = config.FS_DB_PATH
FS_TMP_PATH = config.FS_TMP_PATH

LOG = logging.getLogger()


@click.command()
@click.option(
    "--fs_querydb",
    type=click.Path(exists=True, file_okay=True, resolve_path=True),
    default="fs_query_structures.db",
    help=f"Input: Foldseek query database)",
)
@click.option(
    "--fs_targetdb",
    type=click.Path(exists=True, file_okay=True, resolve_path=True),
    default=FS_DB_PATH,
    help=f"Target Database for Foldseek. default:{FS_DB_PATH}",
)
@click.option(
    "--fs_rawdata",
    type=click.Path(resolve_path=True),
    default="fs_query_structures.raw",
    help=f"Raw output of Foldseek (before convertalis). default: fs_query_structures.raw",
)
@click.option(
    "--fs_results",
    type=click.Path(resolve_path=True),
    default="fs_query_results.m8",
    help=f"Foldseek tabular output",
)
@click.option(
    "--tmp_dir",
    type=click.Path(file_okay=False, dir_okay=True, resolve_path=True),
    required=True,
    default=FS_TMP_PATH,
    help=f"Output: Foldseek temp folder (default:{FS_TMP_PATH})",
)
def run_foldseek(fs_querydb, fs_targetdb, fs_rawdata, fs_results, tmp_dir):
    "Run Foldseek Query DB against Target DB"
    subprocess.call(
        [
            FS_BINARY_PATH,
            "search",
            fs_querydb,
            fs_targetdb,
            fs_rawdata,
            tmp_dir,
            "-s",
            "9",
        ],
        stderr=subprocess.DEVNULL,
    )
    subprocess.call(
        [
            FS_BINARY_PATH,
            "convertalis",
            fs_querydb,
            fs_targetdb,
            fs_rawdata,
            fs_results,
            "--format-output",
            "query,target,qstart,qend,qlen,tstart,tend,tlen,qcov,tcov,bits,evalue",
        ],
        stderr=subprocess.DEVNULL,
    )

    click.echo("DONE")
