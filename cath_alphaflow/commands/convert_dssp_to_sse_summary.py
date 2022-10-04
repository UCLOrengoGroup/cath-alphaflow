from pathlib import Path
import click

from cath_alphaflow.errors import ParseError

from cath_alphaflow.io_utils import (
    yield_first_col,
    get_sse_summary_writer,
)
from cath_alphaflow.models import SecStrSummary
from cath_alphaflow.constants import (
    DEFAULT_DSSP_SUFFIX,
    DEFAULT_HELIX_MIN_LENGTH,
    DEFAULT_STRAND_MIN_LENGTH,
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
def convert_dssp_to_sse_summary(dssp_dir, id_file, sse_out_file, dssp_suffix):
    "Creates summary of secondary structure elements (SSEs) from DSSP files"

    sse_out_writer = get_sse_summary_writer(sse_out_file)

    for file_stub in yield_first_col(id_file):
        click.echo(file_stub)
        dssp_path = Path(dssp_dir) / f"{file_stub}{dssp_suffix}"
        ss_sum = get_sse_summary_from_dssp(dssp_path)
        sse_out_writer.writerow(ss_sum.to_dict())

    click.echo("DONE")


def get_sse_summary_from_dssp(
    dssp_path: Path, *, chopping=None, acc_id=None
) -> SecStrSummary:

    dssp_string = []
    read_headers = False
    domain_length = 0
    ss_total = 0
    if acc_id is None:
        acc_id = dssp_path.stem

    with dssp_path.open("rt") as dssp_fh:
        for line in dssp_fh:
            if line.startswith("  #"):
                read_headers = True
                continue
            if read_headers:
                dssp_code = line[16]
                dssp_string.append(dssp_code)

    segment_dssp = ""
    if chopping:
        for segment in chopping.segments:
            segment_dssp += dssp_string[(segment.start - 1) : segment.end]
    else:
        segment_dssp = dssp_string

    return SecStrSummary.new_from_dssp_str(segment_dssp, acc_id)
