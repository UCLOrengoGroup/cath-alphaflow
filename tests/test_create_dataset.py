from pathlib import Path
import csv
from click.testing import CliRunner
from cath_alphaflow.cli import cli

UNIPROT_IDS = ["P00520"]


def test_cli_usage():
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(cli, [])
        assert result.exit_code == 0
        assert "Usage: cli [OPTIONS] COMMAND [ARGS]..." in result.output


def test_create_dataset_uniprot_ids_fails_without_args():
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(cli, ["create-dataset-uniprot-ids"])
        assert result.exit_code != 0
        assert "Error: Missing option" in result.output


def test_create_dataset_uniprot_ids(create_mock_query):

    expected_uniprot_id = "P00520"
    mock_rows = [
        {
            "uniprot_acc": expected_uniprot_id,
            "sequence_md5": "12acb12cab1c2ab2cab1c2ab12cab1",
            "gene3d_domain_id": "1abcA00__1.10.8.10/23-345",
            "bitscore": 123.4,
            "chopping": "23-345",
        }
    ]
    expected_outfile = "uniprot_ids.csv"

    runner = CliRunner()
    create_mock_query(["uniprot_id"], mock_rows)
    with runner.isolated_filesystem():
        result = runner.invoke(
            cli,
            [
                "create-dataset-uniprot-ids",
                "--uniprot_ids_csv",
                expected_outfile,
                "--max_evalue",
                1e-50,
                "--max_records",
                100,
            ],
        )
        assert result.exit_code == 0
        assert "DONE" in result.output

        expected_outfile_path = Path(expected_outfile)
        assert expected_outfile_path.exists()

        with expected_outfile_path.open("rt") as fh:
            csvreader = csv.reader(fh)
            headers = next(csvreader)
            rows = list(csvreader)

        assert headers == ["uniprot_acc"]
        assert rows == [[expected_uniprot_id]]
