import logging
from pathlib import Path
import click

from cath_alphaflow.errors import ParseError

FILE_POLICY_SKIP = "skip"
FILE_POLICY_ERROR = "error"
DEFAULT_DSSP_CHECK_POLICY = FILE_POLICY_ERROR

LOG = logging.getLogger(__name__)

from cath_alphaflow.io_utils import (
    yield_first_col,
    get_sse_summary_writer,
)
from cath_alphaflow.models.domains import SecStrSummary, AFDomainID
from cath_alphaflow.constants import (
    DEFAULT_DSSP_SUFFIX,
    ID_TYPE_SIMPLE,
    ID_TYPE_AF_DOMAIN,
    ID_TYPE_UNIPROT_DOMAIN,
)


@click.command()
@click.option(
    "--dssp_dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, resolve_path=True),
    required=True,
    help="Input: directory of DSSP files",
)
@click.option(
    "--id_file",
    type=click.File("rt"),
    required=True,
    help="Input: CSV file for containing list of ids to process",
)
@click.option(
    "--id_type",
    type=click.Choice([ID_TYPE_SIMPLE, ID_TYPE_AF_DOMAIN, ID_TYPE_UNIPROT_DOMAIN]),
    default=ID_TYPE_AF_DOMAIN,
    help=f"Option: specify the type of ID to specify the chopping [{ID_TYPE_AF_DOMAIN}]",
)
@click.option(
    "--sse_out_file",
    type=click.File("wt"),
    required=True,
    help="Output: SSE output file",
)
@click.option(
    "--dssp_suffix",
    type=str,
    default=DEFAULT_DSSP_SUFFIX,
    help="Option: suffix to use for dssp files (default: .dssp)",
)
@click.option(
    "--af_version",
    type=int,
    help=f"Option: specify the AF version when parsing uniprot ids",
)
@click.option(
    "--dssp_check_policy",
    type=click.Choice(
        [
            FILE_POLICY_ERROR,
            FILE_POLICY_SKIP,
        ]
    ),
    default=DEFAULT_DSSP_CHECK_POLICY,
    help=f"Option: specify behaviour on dssp check [{DEFAULT_DSSP_CHECK_POLICY}]",
)
def convert_dssp_to_sse_summary(
    dssp_dir, id_file, id_type, sse_out_file, dssp_suffix, af_version, dssp_check_policy
):
    "Creates summary of secondary structure elements (SSEs) from DSSP files"

    if id_type == ID_TYPE_UNIPROT_DOMAIN and af_version is None:
        raise click.UsageError(
            f"option --af_version must be specified when using id_type={id_type}"
        )

    sse_out_writer = get_sse_summary_writer(sse_out_file)

    for id_str in yield_first_col(id_file):

        click.echo(f"Processing '{id_str}' (id:{id_type}) ...")

        dssp_file_stub = None
        acc_id = None
        chopping = None
        if id_type == ID_TYPE_SIMPLE:
            dssp_file_stub = id_str
            acc_id = id_str
        elif id_type == ID_TYPE_AF_DOMAIN:
            _af_domain_id = AFDomainID.from_str(id_str)
            acc_id = _af_domain_id.to_file_stub()
            dssp_file_stub = _af_domain_id.af_chain_id
            chopping = _af_domain_id.chopping
        elif id_type == ID_TYPE_UNIPROT_DOMAIN:
            _af_domain_id = AFDomainID.from_uniprot_str(id_str, version=af_version)
            acc_id = _af_domain_id.to_file_stub()
            dssp_file_stub = _af_domain_id.af_chain_id
            chopping = _af_domain_id.chopping
        else:
            raise click.UsageError(f"failed to recognise id_type={id_type}")

        dssp_path = Path(dssp_dir) / f"{dssp_file_stub}{dssp_suffix}"

        try:
            ss_sum = get_sse_summary_from_dssp(
                dssp_path,
                acc_id=acc_id,
                chopping=chopping,
            )
        except ParseError:
            msg = f"failed to get SSE summary for entry {acc_id}"
            if dssp_check_policy == FILE_POLICY_SKIP:
                LOG.warning(f"{msg} (skipping)")
                continue
            else:
                LOG.error(f"{msg} (use --dssp_check_policy=skip to ignore)")
                raise

        sse_out_writer.writerow(ss_sum.to_dict())

    click.echo("DONE")


def get_sse_summary_from_dssp(
    dssp_path: Path, *, chopping=None, acc_id=None
) -> SecStrSummary:

    dssp_string = ""
    read_headers = False
    if acc_id is None:
        acc_id = dssp_path.stem

    with dssp_path.open("rt") as dssp_fh:
        for line in dssp_fh:
            if line.startswith("  #"):
                read_headers = True
                continue
            if read_headers:
                dssp_code = line[16]
                dssp_string += dssp_code

    segment_dssp = ""
    if chopping:
        for segment in chopping.segments:
            _seg_dssp = dssp_string[(segment.start - 1) : segment.end]
            segment_dssp += _seg_dssp
    else:
        segment_dssp = dssp_string

    try:
        secstr_summary = SecStrSummary.new_from_dssp_str(segment_dssp, acc_id)
    except ParseError:
        msg = f"failed to create SecStrSummary from DSSP {dssp_path} (dssp_length: {len(dssp_string)}, chopping: {chopping}, chopping_dssp: '{segment_dssp}')"
        LOG.error(msg)
        raise

    return secstr_summary
