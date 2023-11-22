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
from cath_alphaflow.models.domains import GeneralDomainID


PDB_IDS = ["2xdqA"]
FIXTURE_PATH = Path(__file__).parent / "fixtures"
EXAMPLE_PDB_FILE = FIXTURE_PATH / "pdb" / "2xdqA.pdb"

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


def test_calculate_normed_radius_of_gyration():
    pdb_path = EXAMPLE_PDB_FILE
    model_id = pdb_path.stem
    domain_id = GeneralDomainID(raw_id=model_id)
    model_structure = PDBParser().get_structure(model_id, pdb_path)
    gyration_radius = calculate_normed_radius_of_gyration(domain_id, model_structure, 5)
    assert gyration_radius == 0
