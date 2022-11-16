from pathlib import Path
import logging
import click

from cath_alphaflow.io_utils import (
    yield_first_col,
)
from cath_alphaflow.models import AFDomainID
from cath_alphaflow.constants import DEFAULT_CIF_SUFFIX
from cath_alphaflow.chopping import chop_cif

LOG = logging.getLogger()

ID_TYPE_AF_DOMAIN = "af"
ID_TYPE_UNIPROT_DOMAIN = "uniprot"
DEFAULT_AF_FRAGMENT_NUMBER = 1


@click.command("chop-cif")
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
    help="Input: CSV file containing list of ids to process (AF domain id in first col)",
)
@click.option(
    "--cif_out_dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, resolve_path=True),
    required=True,
    help="Output: directory of CIF files",
)
@click.option(
    "--id_type",
    type=click.Choice([ID_TYPE_AF_DOMAIN, ID_TYPE_UNIPROT_DOMAIN]),
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
    help=f"Option: suffix to use for mmCIF files (default: {DEFAULT_CIF_SUFFIX})",
)
def chop_cif_command(
    cif_in_dir,
    id_file,
    id_type,
    cif_out_dir,
    cif_suffix,
    af_version,
):
    "Apply chopping to CIF files"

    if id_type == ID_TYPE_UNIPROT_DOMAIN and af_version is None:
        raise click.UsageError(
            f"option --af_version must be specified when using id_type={id_type}"
        )

    for id_str in yield_first_col(id_file):
        if id_type == ID_TYPE_AF_DOMAIN:
            af_domain_id = AFDomainID.from_str(id_str)
        elif id_type == ID_TYPE_UNIPROT_DOMAIN:
            af_domain_id = AFDomainID.from_uniprot_str(
                id_str, fragment_number=DEFAULT_AF_FRAGMENT_NUMBER, version=af_version
            )
        else:
            raise UsageError(f"failed to recognise id_type={id_type}")

        af_chain_stub = af_domain_id.af_chain_id
        chain_cif_path = Path(cif_in_dir) / f"{af_chain_stub}{cif_suffix}"
        if not chain_cif_path.exists():
            msg = f"failed to locate CIF input file {chain_cif_path}"
            LOG.error(msg)
            raise FileNotFoundError(msg)

        domain_cif_path = Path(cif_out_dir) / af_domain_id.to_str() + cif_suffix
        if domain_cif_path.exists():
            msg = f"output domain CIF file already exists: {domain_cif_path}"
            LOG.error(msg)
            raise FileExistsError(msg)

        chop_cif(
            chain_cif_path=chain_cif_path,
            domain_cif_path=domain_cif_path,
            chopping=af_domain_id.chopping,
        )

    click.echo("DONE")
