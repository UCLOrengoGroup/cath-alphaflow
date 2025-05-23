from pathlib import Path
import gzip
from Bio.PDB import MMCIF2Dict
import logging
import click
from cath_alphaflow.io_utils import (
    yield_first_col,
    get_plddt_summary_writer,
)
from cath_alphaflow.seq_utils import cif_to_md5
from cath_alphaflow.models.domains import AFDomainID, LURSummary, pLDDTSummary
from cath_alphaflow.constants import (
    MIN_LENGTH_LUR,
    ID_TYPE_AF_DOMAIN,
    ID_TYPE_UNIPROT_DOMAIN,
)
from cath_alphaflow.errors import ArgumentError

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
    "--id_type",
    type=click.Choice([ID_TYPE_AF_DOMAIN, ID_TYPE_UNIPROT_DOMAIN]),
    default=ID_TYPE_AF_DOMAIN,
    help=f"Option: specify the type of ID to specify the chopping [{ID_TYPE_AF_DOMAIN}]",
)
@click.option(
    "--plddt_stats_file",
    type=click.File("wt"),
    required=True,
    help="Output: pLDDT and LUR output file",
)
@click.option(
    "--af_version",
    type=int,
    help="Option: specify the AF version when parsing uniprot ids",
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
    id_type,
    plddt_stats_file,
    af_version,
    cif_suffix,
):
    "Creates summary of secondary structure elements (SSEs) from DSSP files"

    if id_type == ID_TYPE_UNIPROT_DOMAIN and af_version is None:
        raise click.UsageError(
            f"option --af_version must be specified when using id_type={id_type}"
        )

    plddt_out_writer = get_plddt_summary_writer(plddt_stats_file)

    for af_domain_id_str in yield_first_col(id_file):
        if id_type == ID_TYPE_UNIPROT_DOMAIN:
            af_domain_id = AFDomainID.from_uniprot_str(
                af_domain_id_str, version=af_version
            )
        elif id_type == ID_TYPE_AF_DOMAIN:
            af_domain_id = AFDomainID.from_str(af_domain_id_str)
        else:
            msg = f"failed to understand id_type '${id_type}'"
            raise ArgumentError(msg)

        file_stub = af_domain_id.af_chain_id
        chopping = af_domain_id.chopping

        cif_path = Path(cif_in_dir) / f"{file_stub}{cif_suffix}"
        if not cif_path.exists():
            msg = f"failed to locate CIF input file {cif_path}"
            LOG.error(msg)
            raise FileNotFoundError(msg)

        # get md5 hash of cif file

        md5 = cif_to_md5(cif_path, chopping=chopping)

        avg_plddt = get_average_plddt_from_plddt_string(cif_path, chopping=chopping)
        perc_LUR_summary = get_LUR_residues_percentage(cif_path, chopping=chopping)
        plddt_stats = pLDDTSummary(
            af_domain_id=str(af_domain_id),
            md5=md5,
            avg_plddt=avg_plddt,
            perc_LUR=perc_LUR_summary.LUR_perc,
            residues_total=perc_LUR_summary.residues_total,
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
        chopping_plddt = plddt_strings
    # Calculate LUR
    LUR_perc = 0
    LUR_total = 0
    LUR_res = 0
    LUR_stretch = False
    min_res_lur = MIN_LENGTH_LUR
    for residue in chopping_plddt:
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
    LUR_perc = round(LUR_total / len(chopping_plddt) * 100, 2)

    return LURSummary(
        LUR_perc=LUR_perc, LUR_total=LUR_total, residues_total=len(chopping_plddt)
    )
