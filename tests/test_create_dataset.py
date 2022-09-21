from click.testing import CliRunner
from cath_alphaflow.cli import cli

UNIPROT_IDS = ["P00520"]


def write_uniprot_ids(fh, uniprot_ids):
    fh.write("uniprot_id\n")
    for uniprot_id in uniprot_ids:
        fh.write(uniprot_id + "\n")


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
    runner = CliRunner()
    create_mock_query(["uniprot_id"], UNIPROT_IDS)
    with runner.isolated_filesystem():
        result = runner.invoke(
            cli,
            [
                "create-dataset-uniprot-ids",
                "--uniprot_ids_csv",
                "uniprot_ids.csv",
                "--max_evalue",
                1e-50,
                "--max_records",
                100,
            ],
        )
        assert result.exit_code == 0
        assert "foo" in result.output
