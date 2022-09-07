import logging
import click

from cath_alphaflow.io_utils import get_csv_dictreader
from cath_alphaflow.io_utils import get_csv_dictwriter
from cath_alphaflow.io_utils import chunked_iterable
from cath_alphaflow.db_utils import get_cathora_connection
from cath_alphaflow.db_utils import next_cath_dataset_entry


LOG = logging.getLogger()

DEFAULT_CHUNK_SIZE = 1000
SEQUENCE_TABLE_NAMES = ["sequences", "sequences_extra"]
UNIPROT_ACC_TABLE = ["uniprot_prim_acc"]


@click.command()
@click.option(
    "--uniprot_ids_csv",
    type=click.File("wt"),
    required=True,
    help="Output: CSV file containing UniProt IDs",
)
@click.option(
    "--dbname",
    "gene3d_dbname",
    type=str,
    help="Param: database to use when querying sequences",
)
@click.option(
    "--max_evalue",
    "max_evalue",
    type=str,
    help="Param: restrict UniProt IDs by independent evalue of Gene3D hits",
)
@click.option(
    "--max_records",
    "max_records",
    type=str,
    help="Param: maximum records to return",
)
def create_dataset_uniprot_ids(uniprot_ids_csv, gene3d_dbname, max_evalue, max_records):
    "Creates UniProt IDs for given dataset"

    conn = get_cathora_connection()
    headers = ["uniprot_acc"]

    csv_writer = get_csv_dictwriter(uniprot_ids_csv, fieldnames=headers)
    csv_writer.writeheader()

    click.echo(
        f"Querying Dataset UniProtIDs "
        f"(max_evalue={max_evalue}, max_records={max_records}, file={uniprot_ids_csv.name}) ..."
    )
    for entry in next_cath_dataset_entry(
        conn,
        max_independent_evalue=max_evalue,
        max_records=max_records,
        gene3d_dbname=gene3d_dbname,
    ):
        csv_writer.writerow({"uniprot_acc": entry.uniprot_acc})

    click.echo("DONE")


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
    "--af2_domainlist_ids",
    type=click.File("wt"),
    help="Output: CSV file of AF2 domain ids",
)
@click.option(
    "--af2_chainlist_ids",
    type=click.File("wt"),
    help="Output: CSV file of AF2 chain ids",
)
@click.option(
    "--af2_cath_annotations",
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
    af2_domainlist_ids,
    af2_chainlist_ids,
    af2_cath_annotations,
    gene3d_dbname,
    chunk_size,
):
    "Creates CATH data files for a given dataset"

    # setup reader
    uniprot_reader = get_csv_dictreader(csv_uniprot_ids)
    # skip headers
    next(uniprot_reader)

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

    af2_domainlist_headers = ["af2_domain_id"]
    af2_domainlist_writer = get_csv_dictwriter(
        af2_domainlist_ids, fieldnames=af2_domainlist_headers
    )
    af2_domainlist_writer.writeheader()

    af2_chainlist_headers = ["af2_chain_id"]
    af2_chainlist_writer = get_csv_dictwriter(
        af2_chainlist_ids, fieldnames=af2_chainlist_headers
    )
    af2_chainlist_writer.writeheader()

    af2_cath_annotations_headers = [
        "cath_domain_id",
        "uniprot_acc",
        "md5",
        "bitscore",
        "chopping",
    ]
    af2_cath_annotations_writer = get_csv_dictwriter(
        af2_cath_annotations, fieldnames=af2_cath_annotations_headers
    )
    af2_cath_annotations_writer.writeheader()

    # how we are going to process a chunk of uniprot ids
    def process_uniprot_ids(uniprot_ids):

        conn = get_cathora_connection()

        for entry in next_cath_dataset_entry(
            conn,
            gene3d_dbname=gene3d_dbname,
            uniprot_ids=uniprot_ids,
        ):
            # sort out all variables that we are going to use in our data files
            uniprot_acc = entry.uniprot_acc
            sequence_md5 = entry.sequence_md5
            gene3d_domain_id = entry.gene3d_domain_id
            bitscore = entry.bitscore
            chopping = entry.chopping
            af2_domain_id = "???"
            af2_chain_id = "???"

            # write data
            csv_uniprot_md5_writer.writerow(
                {"uniprot_acc": uniprot_acc, "sequence_md5": sequence_md5}
            )
            af2_domainlist_writer.writerow({"af2_domain_id": af2_domain_id})
            af2_chainlist_writer.writerow({"af2_chain_id": af2_chain_id})
            af2_cath_annotations_writer.writerow(
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
