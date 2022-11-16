import os
from pathlib import Path
import csv

import pytest
from click.testing import CliRunner
from cath_alphaflow.cli import cli
from cath_alphaflow.commands.optimise_domain_boundaries import (
    calculate_domain_id_post_tailchop,
)
from cath_alphaflow.models import AFDomainID


UNIPROT_IDS = ["P00520"]
FIXTURE_PATH = Path(__file__).parent / "fixtures"
EXAMPLE_AF_ID = "AF-P00520-F1-model_v3"
EXAMPLE_CIF_FILE = FIXTURE_PATH / "cif" / f"{EXAMPLE_AF_ID}.cif.gz"

SUBCOMMAND = "optimise-domain-boundaries"


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


def create_fake_cif_dir(cif_path, ids, cif_src=EXAMPLE_CIF_FILE):
    cif_path = Path(str(cif_path))
    if not cif_path.exists():
        cif_path.mkdir()
    for _id in ids:
        cif_path_dest = cif_path / f"{_id}.cif.gz"
        # print(f"Creating symbolic link for test file {cif_src} -> {cif_path_dest}")
        os.symlink(cif_src, f"{cif_path_dest}")
    return cif_path


def test_optimise_domain_boundaries(tmp_path):

    headers = ["header"]
    chain_ids = ["AF-P00520-F1-model_v3", "AF-P00521-F1-model_v3"]
    af_domain_ids = [
        "AF-P00520-F1-model_v3/1-100",
        "AF-P00521-F1-model_v3/800-1123",
    ]
    # TODO: manually check that these results are correct...
    new_af_domain_ids = [
        "AF-P00520-F1-model_v3/61-100",
        "AF-P00521-F1-model_v3/1018-1123",
    ]
    expected_mapping_rows = [
        [row[0], row[1]] for row in zip(af_domain_ids, new_af_domain_ids)
    ]

    runner = CliRunner()
    with runner.isolated_filesystem(temp_dir=tmp_path) as td:

        _dir = Path(f"{td}")

        tmp_cif_path = create_fake_cif_dir(_dir / "cif", chain_ids)
        tmp_id_path = _dir / "ids.csv"
        with tmp_id_path.open("wt") as fh:
            write_ids_to_file(fh, headers, af_domain_ids)

        tmp_list_post_chop_path = _dir / "domain_list_post_tailchop.csv"
        tmp_mapping_post_chop_path = _dir / "domain_mapping_post_tailchop.csv"

        args = (
            SUBCOMMAND,
            "--af_domain_list",
            f"{tmp_id_path}",
            "--af_chain_mmcif_dir",
            f"{tmp_cif_path}",
            "--af_domain_list_post_tailchop",
            f"{tmp_list_post_chop_path}",
            "--af_domain_mapping_post_tailchop",
            f"{tmp_mapping_post_chop_path}",
        )
        # print(f"Running: cath-af-cli {' '.join(args)}")
        result = runner.invoke(cli, args)
        assert result.exit_code == 0
        assert "DONE" in result.output

        assert_csv_matches(
            tmp_list_post_chop_path,
            ["af_domain_id"],
            # convert list of ids to list of list of ids (ie csv rows)
            [[_id] for _id in new_af_domain_ids],
        )
        assert_csv_matches(
            tmp_mapping_post_chop_path,
            ["af_domain_id_orig", "af_domain_id_post_tailchop"],
            expected_mapping_rows,
        )


def assert_csv_matches(csv_path, expected_headers, expected_rows, max_lines=None):

    assert csv_path.exists()
    with csv_path.open("rt") as fh:
        csvreader = csv.reader(fh, delimiter="\t")
        got_headers = next(csvreader)
        got_rows = list(csvreader)

    if max_lines:
        got_headers = got_headers[:max_lines]
        got_rows = got_rows[:max_lines]

    assert got_headers == expected_headers
    assert got_rows == expected_rows


@pytest.mark.parametrize(
    "old_chopping, cutoff, new_chopping",
    [
        # check single domain at the end of the chain using different cutoffs
        ("800-1123", 70, "1018-1123"),
        ("800-1123", 90, "1021-1122"),
        ("800-1123", 95, "1025-1058"),
        # (simple) check domain with multiple segments (default cutoff)
        ("1000-1050_1070-1123", 70, "1018-1050_1070-1123"),
        # TODO: add more tests
    ],
)
def test_calculate_domain_id_post_tailchop(
    tmp_path, old_chopping, cutoff, new_chopping
):
    """
    Test the main function separately from the CLI
    """

    af_chain_id = EXAMPLE_AF_ID
    af_domain_id = AFDomainID.from_str(f"{af_chain_id}/{old_chopping}")
    cif_dir = create_fake_cif_dir(tmp_path, ids=[af_chain_id])

    new_af_domain_id = calculate_domain_id_post_tailchop(
        af_domain_id=af_domain_id, af_chain_mmcif_dir=cif_dir, cutoff_plddt_score=cutoff
    )

    assert new_af_domain_id.chopping.to_str() == new_chopping
