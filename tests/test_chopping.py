import logging
import gzip
from pathlib import Path
import pytest
import tempfile

from Bio.PDB import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO

from cath_alphaflow.models import Chopping
from cath_alphaflow.chopping import ChoppingProcessor, chop_cif

FIXTURE_PATH = Path(__file__).parent / "fixtures"
EXAMPLE_CIF_FILE = FIXTURE_PATH / "cif" / "AF-P00520-F1-model_v3.cif.gz"

LOG = logging.getLogger(__name__)


@pytest.fixture
def example_cif_gzipped_path():
    return EXAMPLE_CIF_FILE


@pytest.fixture
def example_cif_tmpfile(example_cif_gzipped_path):
    tmp_cif_file = tempfile.NamedTemporaryFile(mode="wt", suffix=".cif")
    with gzip.open(str(example_cif_gzipped_path), mode="rt") as fh:
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
    chopping = Chopping.from_str(example_chopping_str)
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
        Chopping.from_str(example_chopping_str),
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
    example_chopping = Chopping.from_str(example_chopping_str)

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
