import os
from pathlib import Path
import csv
import tempfile

import pytest
from click.testing import CliRunner
from Bio.PDB import PDBParser

from cath_alphaflow.cli import cli
from cath_alphaflow.commands.measure_globularity import (
    measure_globularity,
    guess_chopping_from_pdb_file,
    calculate_normed_radius_of_gyration,
    calculate_packing_density,
)
from cath_alphaflow.models.domains import GeneralDomainID, Chopping, Segment
from cath_alphaflow.chopping import chop_structure

PDB_ID = "2xdqA"
AF_ID = "AF-Q96HM7-F1-model_v4"
PDB_DIR = Path(__file__).parent / "fixtures" / "pdb"
EXAMPLE_PDB_FILE = PDB_DIR / f"{PDB_ID}.pdb"
EXAMPLE_AF_FILE = PDB_DIR / f"{AF_ID}.pdb"

SUBCOMMAND = "measure-globularity"


def test_cli_usage():
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(cli, [SUBCOMMAND, "--help"])
        assert result.exit_code == 0
        assert "Usage:" in result.output


def write_ids_to_file(fh, headers, ids):
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(headers)
    for _id in ids:
        writer.writerow([_id])
    fh.flush()


def create_fake_pdb_dir(pdb_path, ids, pdb_src=EXAMPLE_PDB_FILE):
    pdb_path = Path(str(pdb_path))
    if not pdb_path.exists():
        pdb_path.mkdir()
    for _id in ids:
        pdb_path_dest = pdb_path / f"{_id}.pdb"
        # print(f"Creating symbolic link for test file {pdb_src} -> {pdb_path_dest}")
        os.symlink(pdb_src, f"{pdb_path_dest}")
    return pdb_path


def test_guess_chopping_from_pdb_file():
    pdb_path = EXAMPLE_PDB_FILE
    chopping = guess_chopping_from_pdb_file(pdb_path)
    assert chopping.to_str() == "6-459"

    tmp_pdb_file = tempfile.NamedTemporaryFile(mode="wt", suffix=".pdb")
    with pdb_path.open("rt") as fh:
        for line in fh:
            if line.startswith("ATOM"):
                res_num = int(line[22:26].strip())
                if res_num > 50 and res_num < 100:
                    continue
            tmp_pdb_file.write(line)
    tmp_pdb_file.flush()

    chopping = guess_chopping_from_pdb_file(tmp_pdb_file.name)
    assert chopping.to_str() == "6-50_100-459"


def test_globularity_for_domain_matches_chain_with_chopping():
    """
    Check that the globularity is the same for a chopped domain vs a chain with chopping
    """
    pdb_path = PDB_DIR / "5yh0J.pdb"
    model_id = pdb_path.stem

    chain_structure = PDBParser().get_structure(model_id, pdb_path)

    chopping = Chopping.from_str("135-366")
    domain_id = GeneralDomainID(raw_id=model_id, chopping=chopping)

    chain_with_chopping_packing_density = calculate_packing_density(
        domain_id, chain_structure, 5, include_all_atoms=False
    )
    assert chain_with_chopping_packing_density == 11.374

    tmp_domain_pdb = tempfile.NamedTemporaryFile("wt", suffix=".pdb")

    chop_structure(
        domain_id=model_id,
        chain_path=pdb_path,
        domain_path=tmp_domain_pdb.name,
        chopping=chopping,
    )

    domain_structure = PDBParser().get_structure(model_id, tmp_domain_pdb.name)
    domain_packing_density = calculate_packing_density(
        domain_id, domain_structure, 5, include_all_atoms=True
    )
    assert domain_packing_density == chain_with_chopping_packing_density


def test_calculate_normed_radius_of_gyration_af():
    pdb_path = EXAMPLE_AF_FILE
    model_id = pdb_path.stem
    expected_gyration_radius = 0.373

    model_structure = PDBParser().get_structure(model_id, pdb_path)

    chopping = Chopping.from_str("1-432")
    domain_id = GeneralDomainID(raw_id=model_id, chopping=chopping)
    gyration_radius = calculate_normed_radius_of_gyration(domain_id, model_structure, 5)
    assert gyration_radius == expected_gyration_radius

    del chopping
    del domain_id
    del gyration_radius

    # if we change the chopping then the gyration radius should also change
    chopping = Chopping.from_str("100-200")
    domain_id = GeneralDomainID(raw_id=model_id, chopping=chopping)
    gyration_radius = calculate_normed_radius_of_gyration(domain_id, model_structure, 5)
    assert gyration_radius != expected_gyration_radius


def test_calculate_normed_radius_of_gyration_cath():
    pdb_path = EXAMPLE_PDB_FILE
    model_id = pdb_path.stem
    expected_gyration_radius = 0.296

    model_structure = PDBParser().get_structure(model_id, pdb_path)

    chopping = Chopping.from_str("6-459")
    domain_id = GeneralDomainID(raw_id=model_id, chopping=chopping)
    gyration_radius = calculate_normed_radius_of_gyration(domain_id, model_structure, 5)
    assert gyration_radius == expected_gyration_radius

    del chopping
    del domain_id
    del gyration_radius

    # if we don't specify the chopping then the gyration radius should be the same
    domain_id = GeneralDomainID(raw_id=model_id)
    gyration_radius = calculate_normed_radius_of_gyration(domain_id, model_structure, 5)
    assert gyration_radius == expected_gyration_radius

    del domain_id
    del gyration_radius

    # if we change the chopping then the gyration radius should also change
    chopping = Chopping.from_str("100-200")
    domain_id = GeneralDomainID(raw_id=model_id, chopping=chopping)
    gyration_radius = calculate_normed_radius_of_gyration(domain_id, model_structure, 5)
    assert gyration_radius != expected_gyration_radius

    del chopping
    del domain_id
    del gyration_radius


def test_calculate_packing_density_af():
    pdb_path = EXAMPLE_AF_FILE
    model_id = pdb_path.stem

    model_structure = PDBParser().get_structure(model_id, pdb_path)

    chopping = Chopping.from_str("1-432")
    domain_id = GeneralDomainID(raw_id=model_id, chopping=chopping)
    packing_density = calculate_packing_density(domain_id, model_structure, 5)
    assert packing_density == 9.422


def test_calculate_packing_density_cath():
    pdb_path = EXAMPLE_PDB_FILE
    model_id = pdb_path.stem

    model_structure = PDBParser().get_structure(model_id, pdb_path)

    chopping = Chopping.from_str("6-459")
    domain_id = GeneralDomainID(raw_id=model_id, chopping=chopping)
    packing_density = calculate_packing_density(domain_id, model_structure, 5)
    assert packing_density == 12.363


def test_globularity_pdb():
    pdb_path = PDB_DIR / "5yh0J.pdb"
    model_id = pdb_path.stem

    model_structure = PDBParser().get_structure(model_id, pdb_path)

    chopping = Chopping.from_str("135-366")
    domain_id = GeneralDomainID(raw_id=model_id, chopping=chopping)

    packing_density_all = calculate_packing_density(domain_id, model_structure, 5)
    assert packing_density_all == 12.571

    packing_density_just_domain = calculate_packing_density(
        domain_id, model_structure, 5, include_all_atoms=False
    )
    assert packing_density_just_domain == 11.374
