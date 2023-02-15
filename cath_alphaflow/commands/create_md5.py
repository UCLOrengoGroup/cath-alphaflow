import logging
from pathlib import Path
import click
import csv
import hashlib

from Bio import SeqIO

from cath_alphaflow.io_utils import yield_first_col, get_uniprot_md5_summary_writer
from cath_alphaflow.seq_utils import str_to_md5
from cath_alphaflow.models import AFDomainID
from cath_alphaflow.constants import (
    ID_TYPE_AF_DOMAIN,
    ID_TYPE_UNIPROT_DOMAIN,
    ID_TYPE_SIMPLE,
)

DEFAULT_CHUNK_SIZE = 1000000

LOG = logging.getLogger()


@click.command()
@click.option(
    "--id_file",
    type=click.File("rt"),
    required=True,
    help="Input: CSV file containing list of ids to convert from CIF to DSSP",
)
@click.option(
    "--fasta",
    "fasta_file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    required=True,
    help=f"Input: the fasta database containing all AF sequences",
)
@click.option(
    "--id_type",
    type=click.Choice([ID_TYPE_AF_DOMAIN, ID_TYPE_UNIPROT_DOMAIN, ID_TYPE_SIMPLE]),
    default=ID_TYPE_AF_DOMAIN,
    help=f"Option: specify the type of ID to specify the chopping [{ID_TYPE_AF_DOMAIN}]",
)
@click.option(
    "--uniprot_md5_csv",
    "uniprot_md5_csv_file",
    type=click.File("wt"),
    required=True,
    help="Output: UniProt to MD5 CSV file",
)
@click.option(
    "--chunk_size",
    type=int,
    default=DEFAULT_CHUNK_SIZE,
    help=f"Options: process uniprot accessions in chunks (if you run out of memory) [{DEFAULT_CHUNK_SIZE}]",
)
def create_md5(id_file, fasta_file, id_type, uniprot_md5_csv_file, chunk_size):
    "Calculate MD5 for FASTA sequences"
    with uniprot_md5_csv_file as out_fp:
        md5_out_writer = get_uniprot_md5_summary_writer(out_fp)

        # chunk uniprot ids (in case we have many millions)
        for uniprot_ids in yield_first_col_chunked(
            id_file, id_type, chunk_size=chunk_size
        ):
            # work through fasta file, calculate output for relevant records
            with open(fasta_file, "rt") as fasta_fp:
                for record in SeqIO.parse(fasta_fp, "fasta"):
                    click.echo(f"record: {record.id} seq='{record.seq[:10]}...'")

                    _db, uniprot_acc, _name = record.id.split("|")

                    if uniprot_acc not in uniprot_ids:
                        click.echo(f"skipping: {record.id}")
                        continue

                    row_data = {
                        "uniprot_acc": uniprot_acc,
                        "sequence_md5": str_to_md5(str(record.seq)),
                    }
                    md5_out_writer.writerow(row_data)

    click.echo("DONE")


def yield_first_col_chunked(id_file, id_type, chunk_size):
    uniprot_ids = set()
    for id_str in yield_first_col(id_file):
        uniprot_id = None
        if id_type == ID_TYPE_AF_DOMAIN:
            uniprot_id = AFDomainID.from_str(id_str).uniprot_acc
        elif id_type == ID_TYPE_SIMPLE:
            uniprot_id = id_str
        else:
            raise click.UsageError(f"failed to recognise id_type={id_type}")
        uniprot_ids.add(uniprot_id)
        if len(uniprot_ids) % chunk_size == 0:
            yield uniprot_ids
            uniprot_ids = set()

    yield uniprot_ids
