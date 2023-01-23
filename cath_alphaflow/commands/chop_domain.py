from pathlib import Path
import logging
import click

from cath_alphaflow.io_utils import (
    yield_first_col,
)
from cath_alphaflow.models.domains import AFDomainID
from cath_alphaflow.constants import (
    DEFAULT_CIF_SUFFIX,
    ID_TYPE_AF_DOMAIN,
    ID_TYPE_UNIPROT_DOMAIN,
    VALID_CIF_SUFFIXES,
)
from cath_alphaflow.chopping import chop_cif
from cath_alphaflow.errors import ChoppingError

LOG = logging.getLogger()

FILE_POLICY_SKIP = "skip"
FILE_POLICY_ERROR = "error"
FILE_POLICY_OVERWRITE = "overwrite"
DEFAULT_INPUT_FILE_POLICY = FILE_POLICY_SKIP
DEFAULT_OUTPUT_FILE_POLICY = FILE_POLICY_SKIP

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
    "--input_file_policy",
    type=click.Choice(
        [
            FILE_POLICY_ERROR,
            FILE_POLICY_SKIP,
        ]
    ),
    default=DEFAULT_INPUT_FILE_POLICY,
    help=f"Option: specify behaviour on missing input file [{DEFAULT_INPUT_FILE_POLICY}]",
)
@click.option(
    "--output_file_policy",
    type=click.Choice(
        [
            FILE_POLICY_OVERWRITE,
            FILE_POLICY_ERROR,
            FILE_POLICY_SKIP,
        ]
    ),
    default=DEFAULT_OUTPUT_FILE_POLICY,
    help=f"Option: specify behaviour on encountering an existing output file [{DEFAULT_OUTPUT_FILE_POLICY}]",
)
@click.option(
    "--af_version",
    type=int,
    help=f"Option: specify the AF version when parsing uniprot ids",
)
def chop_cif_command(
    cif_in_dir,
    id_file,
    id_type,
    cif_out_dir,
    input_file_policy,
    output_file_policy,
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
            af_domain_id = AFDomainID.from_uniprot_str(id_str, version=af_version)
        else:
            raise click.UsageError(f"failed to recognise id_type={id_type}")

        af_chain_stub = af_domain_id.af_chain_id

        chain_cif_path = None
        for cif_suffix in VALID_CIF_SUFFIXES:
            chain_cif_path = Path(cif_in_dir) / f"{af_chain_stub}{cif_suffix}"
            if chain_cif_path.exists():
                break

        if not chain_cif_path.exists():
            msg = f"failed to locate CIF input file {cif_in_dir}/{af_chain_stub}{VALID_CIF_SUFFIXES}"
            if input_file_policy == FILE_POLICY_SKIP:
                LOG.warning(msg + " (skipping)")
                continue
            elif input_file_policy == FILE_POLICY_ERROR:
                LOG.error(msg)
                raise FileNotFoundError(msg)
            else:
                msg = f"unexpected input file policy {input_file_policy}"
                click.UsageError(msg)

        domain_cif_path = Path(cif_out_dir) / (af_domain_id.to_file_stub() + cif_suffix)
        if domain_cif_path.exists():
            msg = f"output domain CIF file already exists: {domain_cif_path}"
            if output_file_policy == FILE_POLICY_OVERWRITE:
                LOG.warning(msg + " (overwriting)")
            elif output_file_policy == FILE_POLICY_ERROR:
                LOG.error(msg)
                raise FileExistsError(msg)
            elif output_file_policy == FILE_POLICY_SKIP:
                LOG.info(msg + " (skipping)")
                continue
            else:
                msg = f"unexpected output file policy {output_file_policy}"
                raise click.UsageError(msg)

        try:
            chop_cif(
                domain_id=af_domain_id.to_str(),
                chain_cif_path=chain_cif_path,
                domain_cif_path=domain_cif_path,
                chopping=af_domain_id.chopping,
            )
        except ChoppingError as err:
            msg = f"failed to chop cif file {chain_cif_path} with chopping {af_domain_id.chopping}: {err}"
            LOG.error(msg)
            raise

    click.echo("DONE")
