import logging
import gzip
from pathlib import Path
import pytest
import tempfile

from Bio import SeqIO
from Bio.PDB import MMCIFParser
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1

from cath_alphaflow.models.domains import (
    ChoppingSeqres,
    ChoppingPdbResLabel,
    AFDomainID,
)
from cath_alphaflow.chopping import ChoppingProcessor
from cath_alphaflow.chopping import chop_cif
from cath_alphaflow.chopping import chop_structure

FIXTURE_PATH = Path(__file__).parent / "fixtures"
PDB_DIR = FIXTURE_PATH / "pdb"
EXAMPLE_CIF_FILE = FIXTURE_PATH / "cif" / "AF-P00520-F1-model_v3.cif.gz"

MULTI_FRAG_EXAMPLE_CIF_PATH = FIXTURE_PATH / "cif" / "AF-Q15772-F11-model_v4.cif.gz"
MULTI_FRAG_EXAMPLE_FASTA_PATH = FIXTURE_PATH / "fasta" / "Q15772.fasta"

LOG = logging.getLogger(__name__)


@pytest.fixture
def example_cif_gzipped_path():
    return EXAMPLE_CIF_FILE


@pytest.fixture
def example_cif_tmpfile(example_cif_gzipped_path):
    return get_unzipped_cif_file(example_cif_gzipped_path)


def get_tmp_pdb_file(pdb_path):
    tmp_pdb_file = tempfile.NamedTemporaryFile(mode="wt", suffix=".pdb")
    with open(str(pdb_path), mode="rt") as fh:
        for line in fh:
            tmp_pdb_file.write(line)
    tmp_pdb_file.flush()
    return tmp_pdb_file


def get_unzipped_cif_file(zipped_cif_path):
    tmp_cif_file = tempfile.NamedTemporaryFile(mode="wt", suffix=".cif")
    with gzip.open(str(zipped_cif_path), mode="rt") as fh:
        for line in fh:
            tmp_cif_file.write(line)
    return tmp_cif_file


@pytest.fixture
def example_cif_chain(example_cif_tmpfile):
    parser = MMCIFParser()
    structure = parser.get_structure("1abc", example_cif_tmpfile.name)
    return list(structure.get_chains())[0]


@pytest.fixture
def example_chopping_str():
    return "12-23,34-45"


@pytest.fixture
def example_expected_resids(example_chopping_str):
    chopping = ChoppingPdbResLabel.from_str(example_chopping_str)
    expected_resids = [
        (" ", resnum, " ")
        for seg in chopping.segments
        for resnum in range(int(seg.start), int(seg.end) + 1)
    ]
    return expected_resids


def test_chopping_processor(
    example_chopping_str, example_cif_chain, example_expected_resids
):
    domain_resnames = []
    chopper = ChoppingProcessor(
        ChoppingPdbResLabel.from_str(example_chopping_str),
        on_segment_residue=lambda res: domain_resnames.append(res.id),
    )
    chopper.process_chain(example_cif_chain)
    assert domain_resnames == example_expected_resids


def test_chopping_processor_from_str(
    example_chopping_str, example_cif_chain, example_expected_resids
):
    domain_resnames = []
    chopper = ChoppingProcessor(
        example_chopping_str,
        on_segment_residue=lambda res: domain_resnames.append(res.id),
    )
    chopper.process_chain(example_cif_chain)
    assert domain_resnames == example_expected_resids


def test_chop_cif(example_cif_tmpfile, example_chopping_str, example_expected_resids):
    new_cif_tmpfile = tempfile.NamedTemporaryFile(mode="wt", suffix=".cif")
    example_chopping = ChoppingPdbResLabel.from_str(example_chopping_str)

    chop_cif(
        domain_id="1abc",
        chain_cif_path=Path(example_cif_tmpfile.name),
        domain_cif_path=Path(new_cif_tmpfile.name),
        chopping=example_chopping,
    )

    assert Path(new_cif_tmpfile.name).exists()

    parser = MMCIFParser()
    structure = parser.get_structure("1abc", new_cif_tmpfile.name)
    # this structure should only contain residues from the given chopping and nothing else
    resids = [res.get_id() for res in list(structure.get_residues())]

    assert resids == example_expected_resids


def test_chop_pdb():
    example_pdb_file = PDB_DIR / "5yh0J.pdb"
    example_pdb_tmpfile = get_tmp_pdb_file(example_pdb_file)
    example_chopping_str = "135-366"
    # expected_atom_labels = "135-159,198-366"
    example_expected_resids = [
        (" ", i, " ") for i in list(range(135, 159 + 1)) + list(range(198, 366 + 1))
    ]
    new_pdb_tmpfile = tempfile.NamedTemporaryFile(mode="wt", suffix=".pdb")
    example_chopping = ChoppingPdbResLabel.from_str(example_chopping_str)

    chop_structure(
        domain_id="1abc",
        chain_path=Path(example_pdb_tmpfile.name),
        domain_path=Path(new_pdb_tmpfile.name),
        chopping=example_chopping,
    )

    assert Path(new_pdb_tmpfile.name).exists()

    parser = PDBParser()
    structure = parser.get_structure("1abc", new_pdb_tmpfile.name)
    # this structure should only contain residues from the given chopping and nothing else
    resids = [res.get_id() for res in list(structure.get_residues())]

    assert resids == example_expected_resids


def test_chopping_equals_filter_bio_residues():
    example_pdb_file = PDB_DIR / "5yh0J.pdb"
    example_pdb_tmpfile = get_tmp_pdb_file(example_pdb_file)
    example_chopping_str = "135-366"
    # expected_atom_labels = "135-159,198-366"
    example_expected_resids = [
        (" ", i, " ") for i in list(range(135, 159 + 1)) + list(range(198, 366 + 1))
    ]
    new_pdb_tmpfile = tempfile.NamedTemporaryFile(mode="wt", suffix=".pdb")
    example_chopping = ChoppingPdbResLabel.from_str(example_chopping_str)

    # get the residues by chopping the file
    chop_structure(
        domain_id="1abc",
        chain_path=Path(example_pdb_tmpfile.name),
        domain_path=Path(new_pdb_tmpfile.name),
        chopping=example_chopping,
    )

    assert Path(new_pdb_tmpfile.name).exists()

    # get the residues by filtering the Bio.Structure residues
    parser = PDBParser()
    structure = parser.get_structure("1abc", new_pdb_tmpfile.name)
    filtered_residues = example_chopping.filter_bio_residues(structure.get_residues())

    filtered_resids = [res.get_id() for res in list(filtered_residues)]

    # this structure should only contain residues from the given chopping and nothing else
    resids = [res.get_id() for res in list(structure.get_residues())]

    assert resids == example_expected_resids
    assert resids == filtered_resids


def test_chop_multi_fragment():
    new_cif_tmpfile = tempfile.NamedTemporaryFile(mode="wt", suffix=".cif")
    expected_fragment_number = 11
    dom_start = 2944
    dom_end = 3260
    uniprot_with_chopping = f"Q15772/{dom_start}-{dom_end}"
    af_domain_id = AFDomainID.from_uniprot_str(uniprot_with_chopping, version=3)

    assert af_domain_id.fragment_number == expected_fragment_number
    example_chopping = af_domain_id.chopping
    assert example_chopping.map_to_uniprot_residue is not None

    chop_cif(
        domain_id="1abc",
        chain_cif_path=Path(MULTI_FRAG_EXAMPLE_CIF_PATH),
        domain_cif_path=Path(new_cif_tmpfile.name),
        chopping=example_chopping,
    )

    # check that the sequences match up
    seq_entry = next(SeqIO.parse(MULTI_FRAG_EXAMPLE_FASTA_PATH, "fasta"))
    assert len(seq_entry.seq) == 3267

    sequence_from_chopped_fasta = seq_entry.seq[dom_start - 1 : dom_end]

    parser = MMCIFParser()
    # check that the sifts residue numbering looks sensible
    with open(new_cif_tmpfile.name, mode="rt") as fp:
        structure = parser.get_structure(uniprot_with_chopping, fp)

    sequence_from_chopped_cif = "".join(
        [protein_letters_3to1[res.get_resname()] for res in structure.get_residues()]
    )

    assert sequence_from_chopped_fasta == sequence_from_chopped_cif


def test_chop_gzip_cif(example_chopping_str, example_expected_resids):
    """
    check we can transparently use .gz files to read and write
    """

    new_cif_tmpfile = tempfile.NamedTemporaryFile(suffix=".cif.gz")
    example_chopping = ChoppingPdbResLabel.from_str(example_chopping_str)

    chop_cif(
        domain_id="1abc",
        chain_cif_path=Path(EXAMPLE_CIF_FILE),
        domain_cif_path=Path(new_cif_tmpfile.name),
        chopping=example_chopping,
    )

    assert Path(new_cif_tmpfile.name).exists()

    parser = MMCIFParser()
    with gzip.open(new_cif_tmpfile.name, mode="rt") as fp:
        structure = parser.get_structure("1abc", fp)

    # this structure should only contain residues from the given chopping and nothing else
    resids = [res.get_id() for res in list(structure.get_residues())]

    assert resids == example_expected_resids
