import logging
import click

from Bio import SeqIO

from cath_alphaflow.io_utils import get_af_uniprot_md5_summary_writer
from cath_alphaflow.seq_utils import str_to_md5
from cath_alphaflow.errors import ParseError

DEFAULT_CHUNK_SIZE = 1000000

LOG = logging.getLogger()

# AFDB:AF-A0A3B9RYK9-F1
RE_AF_MODEL_ID = re.compile("AFDB:AF-(?P<uniprot_acc>[A-Z0-9]+)-F(?P<frag_num>[0-9]+)$")


@click.command()
@click.option(
    "--fasta",
    "fasta_file",
    type=click.File("rt"),
    required=True,
    help=f"Input: the fasta database containing all AF sequences",
)
@click.option(
    "--uniprot_md5_csv",
    "uniprot_md5_csv_file",
    type=click.File("wt"),
    required=True,
    help="Output: UniProt to MD5 CSV file",
)
def create_md5(fasta_file, uniprot_md5_csv_file):
    "Calculate MD5 for FASTA sequences"
    with uniprot_md5_csv_file as out_fh:
        md5_out_writer = get_af_uniprot_md5_summary_writer(out_fh)
        # work through fasta file, calculate output for relevant records
        for record in SeqIO.parse(fasta_file, "fasta"):
            click.echo(f"record: {record.id} seq='{record.seq[:10]}...'")
            # >AFDB:AF-A0A2L2JPH6-F1
            try:
                af_chain_id = record.id.split(":")[1]
                af_uniprot_id = af_chain_id.split("-")[1]
            except:
                raise ParseError(f"Failed to parse {record.id} to AlphaFold chain ID")
            row_data = {
                "af_chain_id": af_chain_id,
                "uniprot_id": af_uniprot_id,
                "sequence_md5": str_to_md5(str(record.seq)),
            }
            md5_out_writer.writerow(row_data)

    click.echo("DONE")
