import os
from pathlib import Path
import csv
import logging
from click.testing import CliRunner
from cath_alphaflow.cli import cli
from cath_alphaflow.settings import get_default_settings

config = get_default_settings()

FS_BINARY_PATH = Path(config.FS_BINARY_PATH)
FIXTURE_PATH = Path(__file__).parent / "fixtures"
EXAMPLE_CIF_FILE = FIXTURE_PATH / "cif" / "AF-P00520-F1-model_v3.cif.gz"
SUBCOMMAND = "convert-cif-to-foldseek-db"


if not FS_BINARY_PATH.exists():
    msg = f"cannot run tests as foldseek is not installed: {FS_BINARY_PATH}"


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


def create_fake_cif_dir(dirname, ids, cif_src=EXAMPLE_CIF_FILE):
    dir_path = Path(dirname)
    dir_path.mkdir()
    for _id in ids:
        path_dest = dir_path / f"{_id}.cif.gz"
        if path_dest.is_symlink():
            path_dest.unlink()
        os.symlink(cif_src, f"{path_dest}")
    return dir_path


def test_convert_cif_to_foldseek_db(tmp_path, create_cli_runner):
    headers = ["header"]
    ids = ["AF-P00520-F1-model_v3.cif.gz"]

    runner = create_cli_runner(extra_settings={"FS_BINARY_PATH": "foldseek-fake-path"})
    with runner.isolated_filesystem(temp_dir=tmp_path):
        cwd_path = Path.cwd()

        tmp_fs_path = create_fake_cif_dir("cif", ids)
        tmp_id_path = cwd_path / "ids.csv"
        with tmp_id_path.open("wt") as fh:
            write_ids_to_file(fh, headers, ids)
        tmp_foldseek_db_path = cwd_path / "tmp_foldseek_db"
        tmp_foldseek_db_path.mkdir()

        args = (
            SUBCOMMAND,
            "--cif_dir",
            f"{tmp_fs_path}",
            "--fs_querydb_dir",
            f"{tmp_foldseek_db_path}",
        )
        # print("CMD: " + " ".join(args))
        result = runner.invoke(cli, args)
        if result.exception:
            print(f"RESULT.exception: {result.exception}")
        assert result.exception is None
        assert result.exit_code == 0
        assert "DONE" in result.output
