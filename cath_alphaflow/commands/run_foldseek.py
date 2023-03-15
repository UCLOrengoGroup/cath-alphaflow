from email.policy import default
import logging
from pathlib import Path
import os
import glob
import click
import subprocess

from cath_alphaflow.settings import get_default_settings

config = get_default_settings()

FS_BINARY_PATH = config.FS_BINARY_PATH
FS_TMP_PATH = config.FS_TMP_PATH
FS_OVERLAP = config.FS_OVERLAP
DEFAULT_FS_COV_MODE = "0" # overlap over query and target
DEFAULT_FS_ALIGNER = "2" # 3di+AA (fast, accurate)
DEFAULT_FS_FORMAT_OUTPUT = "query,target,qstart,qend,qlen,tstart,tend,tlen,qcov,tcov,bits,evalue"

LOG = logging.getLogger()


@click.command()
@click.option(
    "--fs_querydb",
    type=click.Path(exists=True, file_okay=True, resolve_path=True),
    required=True,
    help=f"Input: Foldseek query database)",
)
@click.option(
    "--fs_targetdb",
    type=click.Path(exists=True, file_okay=True, resolve_path=True),
    required=True,
    help=f"Target Database for Foldseek",
)
@click.option(
    "--fs_rawdata",
    type=click.Path(resolve_path=True),
    default="./fs_query_structures.raw",
    help=f"Raw output of Foldseek (before convertalis). default: fs_query_structures.raw",
)
@click.option(
    "--fs_results",
    type=click.Path(resolve_path=True),
    default="./fs_query_results.m8",
    help=f"Foldseek tabular output",
)
@click.option(
    "--tmp_dir",
    type=click.Path(file_okay=False, dir_okay=True, resolve_path=True),
    required=True,
    default=FS_TMP_PATH,
    help=f"Output: Foldseek temp folder (default:{FS_TMP_PATH})",
)
@click.option(
    "--cov_mode",
    type=click.Choice(["0","1","2","3","4","5"]),
    default=DEFAULT_FS_COV_MODE,
    help=f"Foldseek coverage mode: 0: % of query and target 1: % target 2: % query 3: target seqlen is %of query 4: query seqlen is % of target 5: shortest seq is % of longest (default: 0)",
)
@click.option(
    "--coverage",
    type=float,
    default=FS_OVERLAP,
    help=f'Foldseek overlap (default:{FS_OVERLAP})'
)
@click.option(
    "--fs_bin_path",
    type=click.Path(file_okay=True, resolve_path=True),
    default=FS_BINARY_PATH,
    help=f"Option: directory containing the Foldseek executable. (default: {FS_BINARY_PATH})"
)
@click.option(
    "--alignment-type",
    type=click.Choice(["0","1","2"]),
    default=DEFAULT_FS_ALIGNER,
    help=f"Option: Foldseek alignment engine: 0: 3di alignment 1: TMalign 2: 3di+AA. (default: {DEFAULT_FS_ALIGNER})",
)
def run_foldseek(fs_querydb, fs_targetdb, fs_rawdata, fs_results, tmp_dir, cov_mode, coverage, alignment_type, fs_bin_path):
    "Run Foldseek Query DB against Target DB"
    assert str(fs_rawdata) != ''
    subprocess.run(
        [
            fs_bin_path,
            "search",
            fs_querydb,
            fs_targetdb,
            fs_rawdata,
            tmp_dir,
            "-s",
            "9",
            "--cov-mode",
            str(cov_mode),
            "-c",
            str(coverage),
            "--alignment-type",
            str(alignment_type)
        ],
        stderr=subprocess.DEVNULL,
        check=True,
    )
    subprocess.run(
        [
            fs_bin_path,
            "convertalis",
            fs_querydb,
            fs_targetdb,
            fs_rawdata,
            fs_results,
            "--format-output",
            DEFAULT_FS_FORMAT_OUTPUT,
        ],
        stderr=subprocess.DEVNULL,
        check=True,
    )
    files_to_remove = glob.glob(f"{fs_rawdata}*")
    if files_to_remove:
        for file in files_to_remove:
            os.unlink(file)
    click.echo("DONE")
