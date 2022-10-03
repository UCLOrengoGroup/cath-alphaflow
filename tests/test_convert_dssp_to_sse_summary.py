import os
from pathlib import Path
import csv
from click.testing import CliRunner
from cath_alphaflow.cli import cli
from cath_alphaflow.models import SecStrSummary


UNIPROT_IDS = ["P00520"]
FIXTURE_PATH = Path(__file__).parent / "fixtures"
EXAMPLE_DSSP_FILE = FIXTURE_PATH / "dssp" / "AF-P00520-F1-model_v3.dssp"

SUBCOMMAND = "convert-dssp-to-sse-summary"


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


def create_fake_dssp_dir(dirname, ids, dssp_src=EXAMPLE_DSSP_FILE):
    dir_path = Path(dirname)
    dir_path.mkdir()
    for _id in ids:
        dssp_path_dest = dir_path / f"{_id}.dssp"
        os.symlink(dssp_src, f"{dssp_path_dest}")
    return dir_path


def test_convert_dssp(tmp_path):

    headers = ["header"]
    ids = ["id1", "id2"]

    # TODO: fix the last two cols
    expected_rows = [
        ["id1", "313", "1123", "72.13", "17", "23", "40"],
        ["id2", "313", "1123", "72.13", "17", "23", "40"],
    ]

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
            f"{tmp_dssp_path}",
            "--id_file",
            f"{tmp_id_path}",
            "--sse_out_file",
            f"{tmp_sse_path}",
        )
        result = runner.invoke(cli, args)
        assert result.exit_code == 0
        assert "DONE" in result.output

        expected_outfile_path = Path(tmp_sse_path.name)
        assert expected_outfile_path.exists()

        with expected_outfile_path.open("rt") as fh:
            csvreader = csv.reader(fh, delimiter="\t")
            got_headers = next(csvreader)
            got_rows = list(csvreader)

        assert got_headers == [
            "af_domain_id",
            "ss_res_total",
            "res_count",
            "perc_not_in_ss",
            "sse_H_num",
            "sse_E_num",
            "sse_num",
        ]
        assert len(got_rows) == len(ids)
        assert got_rows == expected_rows


def test_sse_summary():
    test_dssp_str = "HHHHHHHHEEEEEEEHEHEEEEEIIIIIIIIIIIEEH"
    sss = SecStrSummary.new_from_dssp_str(
        test_dssp_str,
        acc_id="test1",
    )
    assert sss.to_dict() == {
        "af_domain_id": "test1",
        "ss_res_total": len(test_dssp_str) - 11,
        "res_count": len(test_dssp_str),
        "perc_not_in_ss": round((100 * 11) / len(test_dssp_str), 2),
        "sse_H_num": 1,
        "sse_E_num": 2,
        "sse_num": 3,
    }
