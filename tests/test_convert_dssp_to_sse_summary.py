from pathlib import Path
import csv
from click.testing import CliRunner
from cath_alphaflow.cli import cli


UNIPROT_IDS = ["P00520"]

SUBCOMMAND = "convert-dssp-to-sse-summary"


def test_cli_usage():
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(cli, [SUBCOMMAND, "--help"])
        assert result.exit_code == 0
        assert "Usage:" in result.output


def write_ids_to_file(fh, headers, ids):
    fh.write("\t".join([h for h in headers]) + "\n")
    for _id in ids:
        fh.write(f"{_id}\n")
    fh.flush()


def create_fake_dssp_dir(dirname, ids):
    dir_path = Path(dirname)
    dir_path.mkdir()
    for _id in ids:
        dssp_path = dir_path / f"{_id}.dssp"
        dssp_path.open("wt")
    return dir_path


def test_convert_dssp(tmp_path):

    headers = ["header"]
    ids = ["id1", "id2"]

    runner = CliRunner()
    with runner.isolated_filesystem(temp_dir=tmp_path):

        cwd_path = Path.cwd()

        tmp_dssp_path = create_fake_dssp_dir("dssp", ids)
        tmp_id_path = cwd_path / "ids.csv"
        with tmp_id_path.open("wt") as fh:
            write_ids_to_file(fh, headers, ids)
        tmp_sse_path = cwd_path / "sse_out.csv"

        args = (
            SUBCOMMAND,
            "--dssp_dir",
            tmp_dssp_path.name,
            "--id_file",
            tmp_id_path.name,
            "--sse_out_file",
            tmp_sse_path.name,
        )
        print(f"cd {Path.cwd()}")
        print(f"cath-af-cli {' '.join(args)}")
        result = runner.invoke(cli, args)
        assert result.exit_code == 0
        assert "DONE" in result.output

        expected_outfile_path = Path(tmp_sse_path.name)
        assert expected_outfile_path.exists()

        with expected_outfile_path.open("rt") as fh:
            csvreader = csv.reader(fh)
            got_headers = next(csvreader)
            got_rows = list(csvreader)

        assert got_headers == headers
        assert len(got_rows) == len(rows)
        assert got_rows == [
            "ss_res_total",
            "res_count",
            "perc_not_in_ss",
            "sse_H_num",
            "sse_E_num",
        ]
