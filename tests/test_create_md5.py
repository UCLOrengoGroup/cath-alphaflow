import os
from tempfile import NamedTemporaryFile
from pathlib import Path
from click.testing import CliRunner
from cath_alphaflow.cli import cli
from cath_alphaflow.commands.create_md5 import (
    create_md5,
)


FIXTURE_PATH = Path(__file__).parent / "fixtures"
EXAMPLE_FASTA_FILE = FIXTURE_PATH / "fasta" / "fasta_test10.fasta"

SUBCOMMAND = "create-md5"


def test_cli_usage():
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(cli, [SUBCOMMAND, "--help"])
        assert result.exit_code == 0
        assert "Usage:" in result.output


def test_create_md5():
    expected_tsv_content = (
        """
af_chain_id	uniprot_id	sequence_md5
AF-A0A2L2JPH6-F1	A0A2L2JPH6	8d7bd48f425bed920868023f24865164
AF-A0A7Y2HFE6-F1	A0A7Y2HFE6	5247969083247b855d4cf4aff8e9e88a
AF-A0A1E4JP01-F1	A0A1E4JP01	9d1f4a47e6171de10d8cdde4cae5d146
AF-A0A543JGU7-F1	A0A543JGU7	cd2e171e2bf4cd6305c8f45cea3bbf6d
AF-A0A428DD76-F1	A0A428DD76	597ef8c3d835b2554d1a79dbaef616e1
    """.strip()
        + "\n"
    )

    uniprot_md5_csv_file = "output.csv"
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(
            cli,
            [
                SUBCOMMAND,
                "--fasta",
                str(EXAMPLE_FASTA_FILE),
                "--uniprot_md5_csv",
                uniprot_md5_csv_file,
            ],
        )
        assert result.exit_code == 0
        assert "DONE" in result.output

        output_path = Path(uniprot_md5_csv_file)
        actual_tsv_content = ""
        for line in output_path.open("rt"):
            actual_tsv_content += line
    assert actual_tsv_content == expected_tsv_content
