from pathlib import Path
import gzip
from Bio.PDB import MMCIF2Dict
import logging
import click
from cath_alphaflow.io_utils import (
    yield_first_col,
    get_plddt_summary_writer,
)
from cath_alphaflow.models import LURSummary, pLDDTSummary
from cath_alphaflow.constants import MIN_LENGTH_LUR

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
    help="Input: CSV file containing list of ids to process from CIF to pLDDT",
)
@click.option(
    "--plddt_stats_file",
    type=click.File("wt"),
    required=True,
    help="Output: pLDDT and LUR output file",
)
@click.option(
    "--cif_suffix",
    type=str,
    default=".cif",
    help="Option: suffix to use for mmCIF files (default: .cif)",
)
def convert_cif_to_plddt_summary(
    cif_in_dir,
    id_file,
    plddt_stats_file,
    cif_suffix,
):
    "Creates summary of secondary structure elements (SSEs) from DSSP files"

    plddt_out_writer = get_plddt_summary_writer(plddt_stats_file)

    for file_stub in yield_first_col(id_file):
        cif_path = Path(cif_in_dir) / f"{file_stub}{cif_suffix}"
        if not cif_path.exists():
            msg = f"failed to locate CIF input file {cif_path}"
            LOG.error(msg)
            raise FileNotFoundError(msg)

        avg_plddt = get_average_plddt_from_plddt_string(cif_path)
        perc_LUR_summary = get_LUR_residues_percentage(cif_path)
        plddt_stats = pLDDTSummary(
            af_domain_id=file_stub,
            avg_plddt=avg_plddt,
            perc_LUR=perc_LUR_summary.LUR_perc,
        )
        plddt_out_writer.writerow(plddt_stats.__dict__)

    click.echo("DONE")


def get_average_plddt_from_plddt_string(
    cif_path: Path, *, chopping=None, acc_id=None
) -> float:
    if acc_id is None:
        acc_id = cif_path.stem
    open_func = open
    if cif_path.name.endswith(".gz"):
        open_func = gzip.open
    with open_func(str(cif_path), mode="rt") as cif_fh:
        mmcif_dict = MMCIF2Dict.MMCIF2Dict(cif_fh)
    chain_plddt = mmcif_dict["_ma_qa_metric_global.metric_value"][0]
    plddt_strings = mmcif_dict["_ma_qa_metric_local.metric_value"]
    chopping_plddt = []
    if chopping:
        for segment in chopping.segments:
            segment_plddt = [
                float(plddt)
                for plddt in plddt_strings[int(segment.start) - 1 : int(segment.end)]
            ]
            chopping_plddt += segment_plddt
        domain_length = len(chopping_plddt)
        average_plddt = round((sum(chopping_plddt) / domain_length), 2)

    else:
        average_plddt = chain_plddt
    return average_plddt


def get_LUR_residues_percentage(cif_path: Path, *, chopping=None, acc_id=None):
    if acc_id is None:
        acc_id = cif_path.stem
    open_func = open
    if cif_path.name.endswith(".gz"):
        open_func = gzip.open
    with open_func(str(cif_path), mode="rt") as cif_fh:
        mmcif_dict = MMCIF2Dict.MMCIF2Dict(cif_fh)
    plddt_strings = mmcif_dict["_ma_qa_metric_local.metric_value"]
    chopping_plddt = []
    if chopping:
        for segment in chopping.segments:
            segment_plddt = [
                float(plddt)
                for plddt in plddt_strings[int(segment.start) - 1 : int(segment.end)]
            ]
            chopping_plddt += segment_plddt
    else:
        segment_plddt = plddt_strings
    # Calculate LUR
    LUR_perc = 0
    LUR_total = 0
    LUR_res = 0
    LUR_stretch = False
    min_res_lur = MIN_LENGTH_LUR
    for residue in segment_plddt:
        plddt_res = float(residue)
        if plddt_res < 70:
            LUR_res += 1
            if LUR_stretch:
                LUR_total += 1

            if LUR_res == min_res_lur and not LUR_stretch:
                LUR_stretch = True
                LUR_total += min_res_lur

        else:
            LUR_stretch = False
            LUR_res = 0
    LUR_perc = round(LUR_total / len(segment_plddt) * 100, 2)

    return LURSummary(
        LUR_perc=LUR_perc, LUR_total=LUR_total, residues_total=len(segment_plddt)
    )
