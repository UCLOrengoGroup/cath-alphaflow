import logging
import click

from cath_alphaflow.io_utils import get_uniprot_id_dictreader
from cath_alphaflow.io_utils import get_csv_dictwriter
from cath_alphaflow.io_utils import chunked_iterable
from cath_alphaflow.db_utils import OraDB


LOG = logging.getLogger()

DEFAULT_CHUNK_SIZE = 1000


@click.command()
@click.option(
    "--csv_uniprot_ids",
    type=click.File("rt"),
    required=True,
    help="Input: CSV file containing UniProt IDs",
)
@click.option(
    "--csv_uniprot_md5",
    type=click.File("wt"),
    help="Output: CSV file of UniProt to MD5 mapping",
)
@click.option(
    "--gene3d_crh_output",
    type=click.File("wt"),
    help="Output: CRH output file for Gene3D domains",
)
@click.option(
    "--af_domainlist_ids",
    type=click.File("wt"),
    help="Output: CSV file of AF2 domain ids",
)
@click.option(
    "--af_chainlist_ids",
    type=click.File("wt"),
    help="Output: CSV file of AF2 chain ids",
)
@click.option(
    "--af_cath_annotations",
    type=click.File("wt"),
    help="Output: CSV file of CATH annotations",
)
@click.option(
    "--dbname",
    "gene3d_dbname",
    type=str,
    help="Param: database to use when querying sequences",
)
@click.option(
    "--chunk",
    "chunk_size",
    type=int,
    default=DEFAULT_CHUNK_SIZE,
    help="Param: size of chunk when processing data",
)
def create_dataset_cath_files(
    csv_uniprot_ids,
    csv_uniprot_md5,
    gene3d_crh_output,
    af_domainlist_ids,
    af_chainlist_ids,
    af_cath_annotations,
    gene3d_dbname,
    chunk_size,
):
    "Creates CATH data files for a given dataset"

    # setup reader
    uniprot_reader = get_uniprot_id_dictreader(csv_uniprot_ids)

    # setup writers
    crh_output_headers = [
        "sequence_md5",
        "cath_domain_id",
        "score",
        "boundaries",
        "resolved",
    ]
    crh_output_writer = get_csv_dictwriter(
        gene3d_crh_output, fieldnames=crh_output_headers
    )
    # NOTE: do not write headers to crh output

    csv_uniprot_md5_headers = ["uniprot_acc", "sequence_md5"]
    csv_uniprot_md5_writer = get_csv_dictwriter(
        csv_uniprot_md5, fieldnames=csv_uniprot_md5_headers
    )
    csv_uniprot_md5_writer.writeheader()

    af_domainlist_headers = ["af_domain_id"]
    af_domainlist_writer = get_csv_dictwriter(
        af_domainlist_ids, fieldnames=af_domainlist_headers
    )
    af_domainlist_writer.writeheader()

    af_chainlist_headers = ["af_chain_id"]
    af_chainlist_writer = get_csv_dictwriter(
        af_chainlist_ids, fieldnames=af_chainlist_headers
    )
    af_chainlist_writer.writeheader()

    af_cath_annotations_headers = [
        "cath_domain_id",
        "uniprot_acc",
        "md5",
        "bitscore",
        "chopping",
    ]
    af_cath_annotations_writer = get_csv_dictwriter(
        af_cath_annotations, fieldnames=af_cath_annotations_headers
    )
    af_cath_annotations_writer.writeheader()

    # how we are going to process a chunk of uniprot ids
    def process_uniprot_ids(uniprot_ids):

        db = OraDB()

        for entry in db.next_cath_dataset_entry(
            gene3d_dbname=gene3d_dbname,
            uniprot_ids=uniprot_ids,
        ):
            # sort out all variables that we are going to use in our data files
            uniprot_acc = entry.uniprot_acc
            sequence_md5 = entry.sequence_md5
            gene3d_domain_id = entry.gene3d_domain_id
            bitscore = entry.bitscore
            chopping = entry.chopping
            af_domain_id = "???"
            af_chain_id = "???"

            # write data
            csv_uniprot_md5_writer.writerow(
                {"uniprot_acc": uniprot_acc, "sequence_md5": sequence_md5}
            )
            af_domainlist_writer.writerow({"af_domain_id": af_domain_id})
            af_chainlist_writer.writerow({"af_chain_id": af_chain_id})
            af_cath_annotations_writer.writerow(
                {
                    "cath_domain_id": gene3d_domain_id,
                    "uniprot_acc": uniprot_acc,
                    "md5": sequence_md5,
                    "bitscore": bitscore,
                    "chopping": chopping,
                }
            )

            # use the 'correct' chopping in both cols
            crh_output_writer.writerow(
                {
                    "cath_domain_id": gene3d_domain_id,
                    "sequence_md5": sequence_md5,
                    "score": bitscore,
                    "boundaries": chopping,
                    "resolved": chopping,
                }
            )

    # process chunks of uniprot ids
    for chunked_uniprot_rows in chunked_iterable(uniprot_reader, chunk_size=chunk_size):
        uniprot_ids = [row.get("uniprot_acc") for row in chunked_uniprot_rows]
        LOG.info("Processing %s UniProtIDs %s...", len(uniprot_ids), uniprot_ids[:3])
        process_uniprot_ids(uniprot_ids)
