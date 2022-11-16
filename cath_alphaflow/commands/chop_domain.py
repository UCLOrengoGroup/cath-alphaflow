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


@click.command("chop_cif")
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
    "--cif_suffix",
    type=str,
    default=DEFAULT_CIF_SUFFIX,
    help=f"Option: suffix to use for mmCIF files (default: {DEFAULT_CIF_SUFFIX})",
)
def chop_cif_command(
    cif_in_dir,
    id_file,
    cif_out_dir,
    cif_suffix,
):
    "Creates structure files corresponding to"

    for af_domain_id_str in yield_first_col(id_file):
        af_domain_id = AFDomainID.from_str(af_domain_id_str)
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
