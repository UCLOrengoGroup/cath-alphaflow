import gzip
import hashlib
import logging
from pathlib import Path
import re

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB import MMCIFParser
from Bio.PDB import PDBParser
from Bio.PDB import Structure
from Bio.PDB.Residue import Residue
from Bio.SeqUtils import seq1
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from cath_alphaflow.models.domains import ChoppingPdbResLabel
from cath_alphaflow.models.domains import SegmentStr


LOG = logging.getLogger(__name__)

RE_RESIDUE_LABEL = re.compile(r"^(?P<residue_number>-?[0-9]+)(?P<insert_code>[A-Z])?$")


def str_to_md5(in_str):
    md5 = hashlib.md5(in_str.encode("utf-8")).hexdigest()
    return md5


def biostructure_to_md5(structure: Structure) -> str:
    """
    Convert the sequence from a PDB structure to an MD5 hash string.
    """
    sequence = ""
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname()
                if seq1(resname) != "X":
                    sequence += seq1(resname)
    return str_to_md5(sequence)


def pdb_to_seq(pdb_file, model=0, chain=0):
    """
    Convert a PDB file to a sequence string.
    """
    parser = PDBParser()
    structure = parser.get_structure("PDB", pdb_file)
    model = structure[model]
    chain = model[chain]
    sequence = ""
    for residue in chain:
        resname = residue.get_resname()
        if seq1(resname) != "X":
            sequence += seq1(resname)
    return sequence


def cif_to_md5(cif_path=Path, chain_id=0, chopping=None):
    """
    Convert a CIF file to a sequence string and return the MD5
    """

    if not cif_path.exists():
        msg = f"failed to locate CIF input file {cif_path}"
        LOG.error(msg)
        raise FileNotFoundError(msg)

    _hdr, seq = cif_to_fasta(cif_path, chain_id=chain_id)
    if chopping:
        # apply chopping to sequence
        seq = "".join(
            [
                seq[int(segment.start) - 1 : int(segment.end)]
                for segment in chopping.segments
            ]
        )

    md5 = str_to_md5(seq)
    return md5


# TODO add chain_id exception
def cif_to_fasta(cif_path=Path, chain_id=0):
    if not cif_path.exists():
        msg = f"failed to locate CIF input file {cif_path}"
        LOG.error(msg)
        raise FileNotFoundError(msg)

    if cif_path.name.endswith(".gz"):
        open_func = gzip.open
    else:
        open_func = open

    with open_func(str(cif_path), mode="rt") as cif_fh:
        header = cif_path.stem
        structure = MMCIF2Dict(cif_fh)

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


def guess_chopping_from_biostructure(
    structure: Structure, target_chain_id=None, assume_all_atom_breaks_are_segments=True
) -> ChoppingPdbResLabel:
    """
    Reverse engineer a `Chopping` object from a Bio Structure model

    'Real' PDB files can have weird residue numberings, so it is not possible to
    reliably know whether a PDB file has been chopped by an external process
    (i.e. a CATH domain) or whether the numbering in the original PDB is just weird.
    The best we can do is guess, i.e. check for cases where the sequential atom number
    is not contiguous.

    If `assume_all_atom_breaks_are_segments` is turned off then we just create a chopping
    with a single segment that spans the entire PDB file.
    """

    segments = []
    start_res = None
    end_res = None
    last_atom_num = None
    last_res_label = None

    # Iterate over all models and chains in the structure
    for model in structure:
        for chain in model:
            # Skip chains that are not the target chain
            if target_chain_id and chain.id != target_chain_id:
                continue

            # Iterate over all residues in the chain
            for residue in chain:

                # Concat the residue number and insertion code to create a unique residue label
                res_id = residue.get_id()
                res_label = f"{res_id[1]}{res_id[2]}".strip()

                # Iterate over all atoms in the residue
                for atom in residue.get_unpacked_list():

                    atom_num = atom.get_serial_number()
                    # LOG.info(f"atom_num: {atom_num}, res_label: {res_label}")

                    # If we haven't specified a chain, then we take the first one we see
                    if not target_chain_id:
                        target_chain_id = chain.id

                    if not last_atom_num:
                        last_atom_num = atom_num
                    if not last_res_label:
                        last_res_label = res_label

                    if not start_res:
                        start_res = res_label
                    end_res = res_label

                    if assume_all_atom_breaks_are_segments:
                        if (
                            res_label != last_res_label
                            and atom_num != last_atom_num + 1
                        ):
                            LOG.info(
                                f"found discontinuity in residue numbering: {last_res_label} -> {res_label} ({last_atom_num} -> {atom_num})"
                            )
                            segments.append(
                                SegmentStr(start=start_res, end=last_res_label)
                            )
                            start_res = None
                            end_res = None

                    last_res_label = res_label
                    last_atom_num = atom_num

    if start_res and end_res:
        segments.append(SegmentStr(start=start_res, end=end_res))

    return ChoppingPdbResLabel(segments=segments)


def guess_chopping_from_pdb_file(
    pdbfile, target_chain_id=None, assume_all_atom_breaks_are_segments=True
) -> ChoppingPdbResLabel:

    pdbpath = Path(str(pdbfile))
    if not pdbpath.exists():
        msg = f"failed to locate input file {pdbpath}"
        LOG.error(msg)
        raise FileNotFoundError(msg)
    parser = PDBParser()
    structure = parser.get_structure(pdbpath.stem, pdbpath)
    return guess_chopping_from_biostructure(
        structure, target_chain_id, assume_all_atom_breaks_are_segments
    )
