import logging
import click

from cath_alphaflow.io_utils import get_csv_dictwriter
from cath_alphaflow.db_utils import OraDB


LOG = logging.getLogger()


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

    db = OraDB()
    headers = ["uniprot_acc"]

    csv_writer = get_csv_dictwriter(uniprot_ids_csv, fieldnames=headers)
    csv_writer.writeheader()

    click.echo(
        f"Querying Dataset UniProtIDs "
        f"(max_evalue={max_evalue}, max_records={max_records}, file={uniprot_ids_csv.name}) ..."
    )
    for entry in db.next_cath_dataset_entry(
        max_independent_evalue=max_evalue,
        max_records=max_records,
        gene3d_dbname=gene3d_dbname,
    ):
        csv_writer.writerow({"uniprot_acc": entry.uniprot_acc})

    click.echo("DONE")
