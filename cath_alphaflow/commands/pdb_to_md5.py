import logging
import csv
from pathlib import Path
import click

from cath_alphaflow.io_utils import get_csv_dictwriter
from cath_alphaflow.seq_utils import pdb_to_seq, str_to_md5
from cath_alphaflow.errors import ParseError

LOG = logging.getLogger()


# click command that takes a PDB file (arg: -i) or directory of PDB files (arg: -d) as input
# and outputs a TSV file with ID to MD5 mappings (arg: -o)
@click.command()
@click.option(
    "-i",
    "--input-file",
    type=click.Path(exists=True, dir_okay=False),
    help="Input PDB file",
)
@click.option(
    "-d",
    "--input-dir",
    type=click.Path(exists=True, file_okay=False),
    help="Input directory of PDB files",
)
@click.option(
    "-o",
    "--output-file",
    type=click.Path(dir_okay=False),
    required=True,
    help="Output TSV file for ID to MD5 mappings",
)
def pdb_to_md5(input_file, input_dir, output_file):
    """Convert PDB files to MD5 checksums of their sequences and write to a TSV file."""
    if (not input_file and not input_dir) or (input_file and input_dir):
        LOG.error("Either --input-file or --input-dir must be provided.")
        return

    with open(output_file, "w", newline="") as output_fh:
        fieldnames = ["id", "md5"]
        writer = get_csv_dictwriter(output_fh, fieldnames=fieldnames)
        writer.writeheader()

        def process_pdb_file(pdb_path):
            LOG.info(f"Processing PDB file: {pdb_path}")
            try:
                pdb_id = Path(pdb_path).stem
                sequence = pdb_to_seq(pdb_path)
                md5_checksum = str_to_md5(sequence)
                writer.writerow({"id": pdb_id, "md5": md5_checksum})
            except Exception as e:
                LOG.error(f"Error processing {pdb_path}: {e}")
                raise ParseError(f"Failed to parse PDB file {pdb_path}") from e

        if input_file:
            process_pdb_file(input_file)
        elif input_dir:
            import os

            for filename in os.listdir(input_dir):
                if filename.endswith(".pdb"):
                    pdb_path = os.path.join(input_dir, filename)
                    process_pdb_file(pdb_path)

    LOG.info(f"MD5 checksums written to {output_file}")
