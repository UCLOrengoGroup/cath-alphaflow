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
        sse_out_writer.writerow(ss_sum.__dict__)

    click.echo("DONE")


def get_sse_summary_from_dssp(dssp_path, *, chopping=None) -> SecStrSummary:

    dssp_string = []
    read_headers = False
    domain_length = 0
    ss_total = 0

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

    # Calculate percentage of residues in secondary structures
    ss_total += segment_dssp.count("H") + segment_dssp.count("E")
    domain_length += len(segment_dssp)
    if domain_length == 0:
        msg = f"failed to find any SS data in DSSP file '{dssp_path}'"
        raise ParseError(msg)
    else:
        perc_not_in_ss = round(((domain_length - ss_total) / domain_length) * 100, 2)

    # Calculate number of SSEs
    sse_H_num = sse_H_res = sse_E_num = sse_E_res = 0
    sse_H = sse_E = False

    # Calculate number of alpha helices
    for residue in dssp_string:
        if residue == "H":
            sse_H_res += 1
            if sse_H_res > DEFAULT_HELIX_MIN_LENGTH and sse_H == False:
                sse_H = True
                sse_H_num += 1

        if residue == "E":
            sse_E_res += 1
            if sse_E_res > DEFAULT_STRAND_MIN_LENGTH and sse_E == False:
                sse_E = True
                sse_E_num += 1

        if residue != "H" and residue != "E":
            sse_H = sse_E = False
            sse_H_res = sse_E_res = 0

    sse_num = sse_H_num + sse_E_num

    ss_sum = SecStrSummary(
        af_domain_id=dssp_path,
        ss_res_total=ss_total,
        res_count=domain_length,
        perc_not_in_ss=perc_not_in_ss,
        sse_H_num=sse_H_num,
        sse_E_num=sse_E_num,
        sse_num=sse_num,
    )

    return ss_sum
