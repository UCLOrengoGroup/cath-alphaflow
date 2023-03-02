import pytest
from pathlib import Path
import csv
import logging
from click.testing import CliRunner
from cath_alphaflow.cli import cli

LOG = logging.getLogger(__name__)

UNIPROT_IDS = ["P00520"]

COMMAND_DB = "create-cath-dataset-from-db"
COMMAND_FILES = "create-cath-dataset-from-files"

INPUT_OPTIONS = ("--csv_uniprot_ids",)
OUTPUT_OPTIONS = (
    "--csv_uniprot_md5",
    "--gene3d_crh_output",
    "--af_domainlist_ids",
    "--af_chainlist_ids",
    "--af_cath_annotations",
)

SHARED_OPTIONS = (
    *INPUT_OPTIONS,
    *OUTPUT_OPTIONS,
)


@pytest.mark.parametrize(
    "subcommand,missing_options",
    [
        (COMMAND_DB, ["--dbname"]),
        (COMMAND_FILES, ["--src_af_uniprot_md5", "--src_crh"]),
    ],
)
def test_create_dataset_usage(subcommand, missing_options):
    runner = CliRunner()
    with runner.isolated_filesystem():

        # provide fake arguments for the shared params and check we still get an error
        # about the missing arguments for the specific command
        shared_args = []
        for opt in INPUT_OPTIONS:
            infilename = opt.replace("--", "") + ".txt"
            shared_args.extend([opt, infilename])
            Path(infilename).touch()
        for opt in OUTPUT_OPTIONS:
            shared_args.extend([opt, "fake_argument"])
        result = runner.invoke(cli, [subcommand, *shared_args])
        assert result.exit_code != 0
        assert "Error: Missing option" in result.output
        # for missing_option in missing_options:
        #     assert missing_option in "".join(result.output)


def test_create_dataset_db_command(create_mock_query):

    expected_uniprot_id = "P00520"
    mock_rows = [
        {
            "uniprot_acc": expected_uniprot_id,
            "sequence_md5": "12acb12cab1c2ab2cab1c2ab12cab1",
            "gene3d_domain_id": "1abcA00__1.10.8.10/23-345",
            "bitscore": 123.4,
            "chopping": "23-345",
            "indp_evalue": 0,
        }
    ]
    expected_outfiles = {}
    shared_args = []
    for opt in SHARED_OPTIONS:
        name = opt.replace("--", "")
        filename = name + ".txt"
        expected_outfiles[name] = filename
        shared_args.extend([opt, filename])

    runner = CliRunner()
    create_mock_query(["uniprot_id"], mock_rows)
    with runner.isolated_filesystem():

        with open("csv_uniprot_ids.txt", "wt") as fp:
            for uniprot_id in UNIPROT_IDS:
                fp.write(uniprot_id + "\n")

        args = (
            "create-cath-dataset-from-db",
            *shared_args,
            "--dbname",
            "test_db_name",
        )

        result = runner.invoke(cli, args)
        if result.exit_code != 0:
            LOG.error(f"cmd: cath-af-cli {' '.join(args)}")
            LOG.error(f"exit_code: {result.exit_code}")
            LOG.error(f"stdout: {result.stdout}")

        assert result.exit_code == 0
        assert "DONE" in result.output

        for name, filename in expected_outfiles.items():
            expected_outfile_path = Path(filename)
            assert expected_outfile_path.exists()
