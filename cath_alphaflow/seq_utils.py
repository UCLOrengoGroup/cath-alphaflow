import hashlib
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
import logging

LOG = logging.getLogger(__name__)


def str_to_md5(in_str):
    md5 = hashlib.md5(in_str.encode("utf-8")).hexdigest()
    return md5


# TODO add chain_id exception
def cif_to_fasta(cif_path=Path, chain_id=0):
    if not cif_path.exists():
        msg = f"failed to locate CIF input file {cif_path}"
        LOG.error(msg)
        # raise FileNotFoundError(msg)
    header = cif_path.stem
    structure = MMCIF2Dict(cif_path)
    if "_entity_poly.pdbx_seq_one_letter_code" in structure:
        sequence = structure["_entity_poly.pdbx_seq_one_letter_code"][0].replace(
            "\n", ""
        )
    else:
        parser = MMCIFParser()
        structure = parser.get_structure(header, cif_path)
        model = structure[0]
        chain = model[chain_id]
        sequence = ""
        for residue in chain:
            resname = residue.get_resname()
            if (
                resname
                not in {
                    "HOH",
                    "H2O",
                }
                and seq1(resname) != "X"
            ):
                sequence += seq1(resname)
    return header, sequence


# TODO Test that it's still working
def write_fasta_file(header, sequence, fasta_out_file):
    """Write a FASTA file with the given header and sequence to the specified directory,
    optionally appending to a multi-FASTA file if it already exists.
    """
    record = SeqRecord(Seq(sequence), id=header)
    SeqIO.write(record, fasta_out_file, "fasta")


def combine_fasta_files(fasta_in_dir, fasta_out_file):
    """
    Combines all FASTA files in fasta_in_dir into a single multiFASTA file at fasta_out_file
    """
    # Combine the sequences from each FASTA file into a single multiFASTA file
    with fasta_out_file as fasta_out_fh:
        for fasta_file in Path(fasta_in_dir).iterdir():
            if fasta_file.suffix in [".fasta", ".fa"]:
                for seq_record in SeqIO.parse(fasta_file, "fasta"):
                    SeqIO.write(seq_record, fasta_out_fh, "fasta")
