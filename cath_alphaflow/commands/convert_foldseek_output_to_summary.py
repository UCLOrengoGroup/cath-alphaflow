import logging
from unittest import result
import click
from cath_alphaflow.io_utils import (
    yield_first_col,
    get_foldseek_reader,
    get_foldseek_summary_writer,
)
from cath_alphaflow.settings import get_default_settings
from cath_alphaflow.constants import DEFAULT_FS_BITS_CUTOFF, DEFAULT_FS_OVERLAP
from cath_alphaflow.models.domains import FoldseekSummary

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
    type=click.File("rt"),
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
    unique_af_ids = set()
    unique_af_ids.add("NOHIT")
    best_hits = set()
    foldseek_results_writer = get_foldseek_summary_writer(fs_results)
    for file_stub in yield_first_col(id_file):
        unique_af_ids.add(file_stub)
    foldseek_reader = get_foldseek_reader(fs_input_file)
    for foldseek_result_as_dict in foldseek_reader:
        result = FoldseekSummary(**foldseek_result_as_dict)
        if result.query.endswith(".cif"):
            result.query = result.query[:-4]
        if (
            result.query not in best_hits
            and float(result.tcov) >= DEFAULT_FS_OVERLAP
            and int(result.bits) >= DEFAULT_FS_BITS_CUTOFF
            and result.query in unique_af_ids
        ):
            best_hits.add(result.query)
            foldseek_results_writer.writerow(result.__dict__)
