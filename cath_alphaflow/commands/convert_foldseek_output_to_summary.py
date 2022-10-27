import logging
import click
from cath_alphaflow.io_utils import yield_first_col, get_foldseek_summary_writer
from cath_alphaflow.settings import get_default_settings
from cath_alphaflow.constants import DEFAULT_FS_BITS_CUTOFF, DEFAULT_FS_OVERLAP
from cath_alphaflow.models import FoldseekSummary

config = get_default_settings()

LOG = logging.getLogger()


@click.command()
@click.option(
    "--id_file",
    type=click.File("rt"),
    required=True,
    help="Input: CSV file containing list of ids to convert from CIF to DSSP",
)
@click.option(
    "--fs_input_file",
    type=click.Path(exists=True, file_okay=True, resolve_path=True),
    default="fs_query_results.m8",
    help=f"Foldseek tabular output as input",
)
@click.option(
    "--fs_results",
    type=click.File("wt"),
    default="fs_hits.tsv",
    help=f"Foldseek hits file",
)
def convert_foldseek_output_to_summary(id_file, fs_input_file, fs_results):
    fs_results_writer = get_foldseek_summary_writer(fs_results)
    best_hits_dict = extract_best_hits_foldseek(fs_input_file)
    for file_stub in yield_first_col(id_file):
        extract_hits_from_hits_dict_foldseek(file_stub, best_hits_dict)


def extract_best_hits_foldseek(fs_input_file):
    "Convert Foldseek tabular output to summary of best hits"
    results_dict = {}
    with open(fs_input_file, "rt") as fs_fh:
        for line in fs_fh:
            line = line.rstrip()
            (
                query,
                target,
                qstart,
                qend,
                qlen,
                tstart,
                tend,
                tlen,
                qcov,
                tcov,
                bits,
                evalue,
            ) = line.split()
            if query.endswith(".cif"):
                query = query[:-4]
            if query in results_dict:
                continue
            results_dict[query] = [
                target,
                qstart,
                qend,
                qlen,
                tstart,
                tend,
                tlen,
                qcov,
                tcov,
                bits,
                evalue,
            ]
    return results_dict


def extract_hits_from_hits_dict_foldseek(query_id, best_hits_dict):
    no_hits_placeholder = ""
    if query_id in best_hits_dict:
        (
            target,
            qstart,
            qend,
            qlen,
            tstart,
            tend,
            tlen,
            qcov,
            tcov,
            bits,
            evalue,
        ) = best_hits_dict[query_id]
    return FoldseekSummary(
        target,
        qstart,
        qend,
        qlen,
        tstart,
        tend,
        tlen,
        qcov,
        tcov,
        bits,
        evalue=no_hits_placeholder,
    )
