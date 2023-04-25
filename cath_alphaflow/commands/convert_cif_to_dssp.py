import logging
from pathlib import Path
import click
import subprocess
from cath_alphaflow.io_utils import yield_first_col
from cath_alphaflow.models.domains import AFDomainID, AFChainID
from cath_alphaflow.constants import (
    DEFAULT_CIF_SUFFIX,
    DEFAULT_DSSP_SUFFIX,
    ID_TYPE_AF_DOMAIN,
    ID_TYPE_UNIPROT_DOMAIN,
    ID_TYPE_AF_CHAIN,
)
from cath_alphaflow.settings import get_default_settings
from cath_alphaflow.errors import ArgumentError


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
    "--id_type",
    type=click.Choice([ID_TYPE_AF_DOMAIN, ID_TYPE_UNIPROT_DOMAIN, ID_TYPE_AF_CHAIN]),
    default=ID_TYPE_AF_DOMAIN,
    help=f"Option: specify the type of ID to specify the chopping [{ID_TYPE_AF_DOMAIN}]",
)
@click.option(
    "--af_version",
    type=int,
    help=f"Option: specify the AF version when parsing uniprot ids",
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
def convert_cif_to_dssp(
    cif_in_dir, id_file, id_type, af_version, cif_suffix, dssp_suffix, dssp_out_dir
):
    "Converts CIF to DSSP files"

    for af_id_str in yield_first_col(id_file):
        af_chain_id = None
        chopping = None
        if id_type == ID_TYPE_UNIPROT_DOMAIN:
            af_domain_id = AFDomainID.from_uniprot_str(af_id_str, version=af_version)
            af_chain_id = af_domain_id.af_chain_id
            chopping = af_domain_id.chopping
        elif id_type == ID_TYPE_AF_DOMAIN:
            af_domain_id = AFDomainID.from_str(af_id_str)
            af_chain_id = af_domain_id.af_chain_id
            chopping = af_domain_id.chopping
        elif id_type == ID_TYPE_AF_CHAIN:
            af_chain_id = AFChainID.from_str(af_id_str)
        else:
            msg = f"failed to understand id_type '${id_type}'"
            raise ArgumentError(msg)

        file_stub = af_chain_id

        cif_path = Path(cif_in_dir) / f"{file_stub}{cif_suffix}"
        dssp_path = Path(dssp_out_dir) / f"{file_stub}{dssp_suffix}"

        LOG.debug(f"Running DSSP: {cif_path} {dssp_path}")
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

    args = [
        DSSP_BINARY_PATH,
    ]

    if not DSSP_PDB_DICT is None:
        args.extend(
            [
                "--mmcif-dictionary",
                DSSP_PDB_DICT,
            ]
        )

    args.extend(
        [
            "--output-format",
            "dssp",
            f"{cif_path}",
            f"{dssp_path}",
        ]
    )

    LOG.debug(f"Running: `{' '.join(args)}`")

    subprocess.run(
        args,
        stderr=subprocess.DEVNULL,
        check=True,
    )
