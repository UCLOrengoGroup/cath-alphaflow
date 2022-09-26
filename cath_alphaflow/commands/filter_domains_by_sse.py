from ast import Raise
import logging
from pathlib import Path
import click
import subprocess

from cath_alphaflow.io_utils import get_af_domain_id_reader
from cath_alphaflow.io_utils import get_csv_dictwriter
from cath_alphaflow.models import AFDomainID
from cath_alphaflow.settings import get_default_settings

config = get_default_settings()

DSSP_BINARY_PATH = config.DSSP_BINARY_PATH
DSSP_PDB_DICT = config.DSSP_PDB_DICT

LOG = logging.getLogger()

# Input:
# - `FILE:AF_CHAIN_MMCIF`
# - `FILE:AF_DOMAIN_LIST`
# - `FILE:AF_DSSP_FOLDER`

# Output:
# - `FILE:AF_DSSP_FILE`
# - `FILE:AF_DOMAIN_SSE` -- `(AF_domain_ID_tailchop, SSE_num, PERC_SS)`


@click.command()
@click.option(
    "--af_chain_mmcif_dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, resolve_path=True),
    required=True,
    help="Input: directory of mmCIF files",
)
@click.option(
    "--af_domain_list_post_tailchop",
    type=click.File("rt"),
    required=True,
    help="Input: CSV file for AF2 domain list after chopping",
)
@click.option(
    "--af_dssp_folder",
    "af_dssp_folder",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, resolve_path=True),
    required=True,
    help="Output: DSSP Output Folder",
)
def filter_domains_by_sse(
    af_chain_mmcif_dir, af_domain_list_post_tailchop, af_dssp_folder
):
    "Runs DSSP on AF2 MMCIF files to obtain percentage of secondary structure and Secondary Structure Elements (SSEs)"

    af_domain_list_reader = get_af_domain_id_reader(af_domain_list_post_tailchop)
    click.echo(
        f"Running DSSP on AF2 MMCIFs"
        f"(mmcif_dir={af_chain_mmcif_dir}, in_file={af_domain_list_post_tailchop.name}, "
        # f"out_file={af_domain_list_post_tailchop_writer.name} ) ..."
    )
    for af_domain_id in af_domain_list_reader:
        run_process_dssp(
            af_domain_id=af_domain_id,
            af_chain_mmcif_dir=af_chain_mmcif_dir,
            af_dssp_folder=af_dssp_folder,
        )
    click.echo("DONE")


def run_process_dssp(af_domain_id, af_chain_mmcif_dir, af_dssp_folder):
    input_dssp = f"{af_chain_mmcif_dir}/{af_domain_id.to_chain_str()}.cif"
    output_dssp = f"{af_dssp_folder}/{af_domain_id.to_chain_str()}.dssp"

    click.echo(f"Running DSSP: " f"{input_dssp}" f"{output_dssp}")
    subprocess.run(
        [
            DSSP_BINARY_PATH,
            "--mmcif-dictionary",
            DSSP_PDB_DICT,
            f"{af_chain_mmcif_dir}/{af_domain_id.to_chain_str()}.cif",
            f"{af_dssp_folder}/{af_domain_id.to_chain_str()}.dssp",
        ],
        stderr=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        check=True,
    )
    with open(output_dssp, "r") as output_dssp_fh:
        dssp_string = []
        print(output_dssp_fh)
        read_file = 0
        for line in output_dssp_fh:
            if line.startswith("  #"):
                read_file = 1
                continue
            if read_file == 1:
                dssp_code = line[16]
                dssp_string.append(dssp_code)
        LOG.info(f"{af_domain_id.to_str()} {dssp_string}")

    # subprocess.run(["rm", f"{af_dssp_folder}/{af_domain_id.to_chain_str()}.dssp"])
