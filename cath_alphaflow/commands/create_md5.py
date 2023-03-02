import logging
from pathlib import Path
import click
import csv
import hashlib
import re

from Bio import SeqIO

from cath_alphaflow.errors import ParseError
from cath_alphaflow.io_utils import yield_first_col, get_csv_dictwriter
from cath_alphaflow.models.domains import AFDomainID
from cath_alphaflow.constants import (
    ID_TYPE_AF_DOMAIN,
    ID_TYPE_UNIPROT_DOMAIN,
    ID_TYPE_SIMPLE,
)

DEFAULT_CHUNK_SIZE = 1000000

LOG = logging.getLogger()

# AFDB:AF-A0A3B9RYK9-F1
RE_AF_MODEL_ID = re.compile("AFDB:AF-(?P<uniprot_acc>[A-Z0-9]+)-F(?P<frag_num>[0-9]+)$")


@click.command()
@click.option(
    "--id_file",
    type=click.File("rt"),
    help="Input: CSV file containing list of ids to convert MD5 (optional)",
)
@click.option(
    "--fasta",
    "fasta_file",
    type=click.File("rt"),
    required=True,
    help=f"Input: the fasta database containing all AF sequences",
)
@click.option(
    "--id_type",
    type=click.Choice([ID_TYPE_AF_DOMAIN, ID_TYPE_UNIPROT_DOMAIN, ID_TYPE_SIMPLE]),
    default=ID_TYPE_SIMPLE,
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

    csv_fieldnames = ["uniprot_acc", "sequence_md5", "af_model_id"]
    csv_writer = get_csv_dictwriter(uniprot_md5_csv_file, fieldnames=csv_fieldnames)
    csv_writer.writeheader()

    id_file_has_header = False

    LOG.info(f"FASTA_FILE:            {fasta_file.name}")
    LOG.info(f"UNIPROT_MD5_CSV_FILE:  {uniprot_md5_csv_file.name}")
    LOG.info(f"ID_FILE:               {id_file.name if id_file else None}")
    LOG.info(f"ID_TYPE:               {id_type}")
    LOG.info(f"ID_FILE_HAS_HEADER:    {id_file_has_header}")
    LOG.info(f"CHUNK_SIZE:            {chunk_size}")

    rows_written = 0
    if id_file:
        # chunk uniprot ids (in case we have many millions)
        for chunk_num, uniprot_ids in enumerate(
            yield_first_col_chunked(
                id_file, id_type, chunk_size=chunk_size, header=id_file_has_header
            ),
            1,
        ):
            LOG.info(
                f"Processing UniProts [chunk={chunk_num}]: uniprot_ids={len(uniprot_ids)} (e.g. {sorted(uniprot_ids)[:2]})"
            )
            rows_written += process_fasta_file(
                fasta_file, csv_writer, uniprot_ids=uniprot_ids
            )
    else:
        rows_written = process_fasta_file(fasta_file, csv_writer)

    LOG.info(f"Wrote {rows_written} rows")

    LOG.info("DONE")


def process_fasta_file(fasta_file, csv_writer, uniprot_ids=None):

    rows_written = 0

    fasta_file.seek(0)

    # work through fasta file, calculate output for relevant records
    for record_num, record in enumerate(SeqIO.parse(fasta_file, "fasta"), 1):

        # >AFDB:AF-A0A3B9RYK9-F1
        match = RE_AF_MODEL_ID.match(record.id)
        if not match:
            msg = f"failed to parse AF model id from fasta id '{record.id}' (seq: {record_num})"
            raise ParseError(msg)

        uniprot_acc = match.group("uniprot_acc")
        frag_num = match.group("frag_num")

        if uniprot_ids and uniprot_acc not in uniprot_ids:
            LOG.debug(
                f"id={record.id} seq='{record.seq[:10]}...' uniprot_acc={uniprot_acc} (skipping)"
            )
            continue

        rows_written += 1

        row_data = {
            "uniprot_acc": uniprot_acc,
            "sequence_md5": str_to_md5(record.seq),
            "af_model_id": record.id,
        }
        csv_writer.writerow(row_data)

    return rows_written


def str_to_md5(in_str):
    md5 = hashlib.md5(in_str.encode("utf-8")).hexdigest()
    return md5


def yield_first_col_chunked(id_file, id_type: str, chunk_size: int, header: bool):

    uniprot_ids = set()
    for id_str in yield_first_col(id_file, header=header):
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
