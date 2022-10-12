import os
from pathlib import Path
from click.testing import CliRunner
from cath_alphaflow.cli import cli
from cath_alphaflow.commands.extract_plddt_and_lur import (
    get_average_plddt_from_plddt_string,
    get_LUR_residues_percentage,
)
from cath_alphaflow.models import Chopping, LURSummary, Segment


UNIPROT_IDS = ["P00520"]
FIXTURE_PATH = Path(__file__).parent / "fixtures"
EXAMPLE_CIF_FILE = FIXTURE_PATH / "cif" / "AF-P00520-F1-model_v3.cif.gz"

SUBCOMMAND = "convert-cif-to-plddt-summary"


def test_cli_usage():
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(cli, [SUBCOMMAND, "--help"])
        assert result.exit_code == 0
        assert "Usage:" in result.output


def create_fake_cif_path(dirname, cif_id, cif_src=EXAMPLE_CIF_FILE):
    dir_path = Path(dirname)
    dir_path.mkdir()

    cif_path_dest = dir_path / f"{cif_id}.cif.gz"
    os.symlink(cif_src, f"{cif_path_dest}")
    return cif_path_dest


def test_extract_plddt_summary(tmp_path):
    acc_id = "test1"
    cif_path = create_fake_cif_path(tmp_path.name, acc_id)
    chopping = Chopping(segments=[Segment("10", "20")])

    average_plddt = get_average_plddt_from_plddt_string(
        cif_path, chopping=chopping, acc_id=acc_id
    )

    assert average_plddt == 33.71

    chopping = Chopping(segments=[Segment("10", "20"), Segment("20", "35")])

    average_plddt = get_average_plddt_from_plddt_string(
        cif_path, chopping=chopping, acc_id=acc_id
    )

    assert average_plddt == 32.88


def get_total_residues_from_chopping(chopping):
    return sum([int(seg.end) - int(seg.start) + 1 for seg in chopping.segments])


def test_extract_LUR_summary(tmp_path):
    acc_id = "test1"
    cif_path = create_fake_cif_path(tmp_path.name, acc_id)
    chopping = Chopping(segments=[Segment("10", "20")])

    lur_summary = get_LUR_residues_percentage(
        cif_path, chopping=chopping, acc_id=acc_id
    )

    assert lur_summary == LURSummary(
        LUR_perc=100.0,
        LUR_total=11,
        residues_total=get_total_residues_from_chopping(chopping),
    )
    # clean up after test
    del chopping
    del lur_summary

    chopping = Chopping(segments=[Segment("1", "200"), Segment("200", "1120")])

    lur_summary = get_LUR_residues_percentage(
        cif_path, chopping=chopping, acc_id=acc_id
    )

    assert lur_summary == LURSummary(
        LUR_perc=47.55,
        LUR_total=533,
        residues_total=get_total_residues_from_chopping(chopping),
    )
    del chopping
    del lur_summary
