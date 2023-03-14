import logging
from unittest import result
import click
from cath_alphaflow.io_utils import (
    yield_first_col,
    get_foldseek_reader,
    get_foldseek_summary_writer,
)
from cath_alphaflow.settings import get_default_settings, DEFAULT_AF_VERSION
from cath_alphaflow.constants import DEFAULT_FS_BITS_CUTOFF, DEFAULT_FS_OVERLAP,ID_TYPE_AF_DOMAIN,ID_TYPE_UNIPROT_DOMAIN
from cath_alphaflow.models.domains import FoldseekSummary,AFDomainID
from cath_alphaflow.errors import ArgumentError

config = get_default_settings()

LOG = logging.getLogger()


@click.command()
@click.option(
    "--id_file",
    type=click.File("rt"),
    required=True,
    help="Input: CSV file containing list of ids to process.",
)
@click.option(
    "--id_type",
    type=click.Choice([ID_TYPE_AF_DOMAIN, ID_TYPE_UNIPROT_DOMAIN]),
    default=ID_TYPE_AF_DOMAIN,
    help=f"Option: specify the type of ID in id_file [{ID_TYPE_AF_DOMAIN}]",
)
@click.option(
    "--fs_input_file",
    type=click.File("rt"),
    required=True,
    help=f"Foldseek tabular output as input",
)
@click.option(
    "--fs_results",
    type=click.File("wt"),
    required=True,
    help=f"Foldseek results summary file",
)
@click.option(
    "--af_version",
    type=int,
    default=DEFAULT_AF_VERSION,
    help=f"Option: specify the AF version when parsing uniprot ids. (default:{DEFAULT_AF_VERSION})",
)

def convert_foldseek_output_to_summary(id_file, fs_input_file, fs_results, id_type, af_version):
    unique_af_ids = set()
    unique_af_ids.add("NOHIT")
    best_hit_by_query = {}
    foldseek_results_writer = get_foldseek_summary_writer(fs_results)
    # Build set of unique AF IDs
    for af_domain_id_str in yield_first_col(id_file):
        if id_type == ID_TYPE_UNIPROT_DOMAIN:
            af_domain_id = AFDomainID.from_uniprot_str(
                af_domain_id_str,version=af_version
            )
            af_domain_id_str = af_domain_id.to_file_stub()
        elif id_type == ID_TYPE_AF_DOMAIN:
            af_domain_id = AFDomainID.from_str(af_domain_id_str)
            af_domain_id_str = af_domain_id.to_file_stub()
        else:
            msg = f"failed to understand id_type '${id_type}'"
            raise ArgumentError(msg)
        unique_af_ids.add(af_domain_id_str)
   
    # Build Foldseek Reader
    foldseek_reader = get_foldseek_reader(fs_input_file)

    # Extract best hit per query for filtered ids
    for foldseek_result_as_dict in foldseek_reader:
        result = FoldseekSummary(**foldseek_result_as_dict)
        
        """
        result: 
        {'query': 'AF-A0A059CHW2-F1-model_v4-22-322.cif', 'target': '1xhlA00', 'qstart': '7', 'qend': '271', 
         'qlen': '301', 'tstart': '3', 'tend': '252', 'tlen': '274', 'qcov': '0.880', 'tcov': '0.912', 
         'bits': '509', 'evalue': '1.305E-11'}
        """
        
        af_query_id = AFDomainID.from_foldseek_query(result.query)
        
        
        af_query_id_str = af_query_id.from_str(str(af_query_id)).to_file_stub()
    
        if af_query_id_str not in unique_af_ids:
            continue
        
        if float(result.tcov) <= DEFAULT_FS_OVERLAP or int(result.bits) <= DEFAULT_FS_BITS_CUTOFF:
            continue

        if af_query_id_str not in best_hit_by_query:
            best_hit_by_query[af_query_id_str] = result

        if int(result.bits) > int(best_hit_by_query[af_query_id_str].bits):
            best_hit_by_query[af_query_id_str] = result
        


        
    # Write results to summary file
    for af_query_id, result in best_hit_by_query.items():   
        foldseek_results_writer.writerow(result.__dict__)
