import hashlib
import logging
import re

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB import MMCIFParser
from Bio.PDB import Structure
from Bio.PDB.Residue import Residue
from Bio.SeqUtils import seq1
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path

LOG = logging.getLogger(__name__)

RE_RESIDUE_LABEL = re.compile(r"^(?P<residue_number>-?[0-9]+)(?P<insert_code>[A-Z])?$")


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
            if seq1(resname) != "X":
                sequence += seq1(resname)
    return header, sequence


def write_fasta_file(header, sequence, fasta_out_file):
    """Write a FASTA file with the given header and sequence to the specified directory,
    optionally appending to a multi-FASTA file if it already exists.
    """
    record = SeqRecord(Seq(sequence), id=header, description="")
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


def get_local_plddt_for_res(
    structure: Structure,
    residue: Residue,
    *,
    model_num: int = 0,
    chain_id: str = "A",
    atom_type: str = "CA",
):
    """
    Returns the local pLDDT score for a given residue
    """
    if isinstance(residue, Residue):
        res = residue.id
    elif isinstance(residue, int):
        res = (" ", residue, " ")
    elif isinstance(residue, str):
        match = RE_RESIDUE_LABEL.match(residue)
        if not match:
            raise ValueError(f"failed to parse '{residue}' as residue label")
        ins_code = match.group("insert_code")
        if not ins_code:
            ins_code = " "
        res = (" ", int(match.group("residue_number")), ins_code)

    return structure[model_num][chain_id][res][atom_type].get_bfactor()
