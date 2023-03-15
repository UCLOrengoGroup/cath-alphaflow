import pytest
from pathlib import Path
import logging
from click.testing import CliRunner
from cath_alphaflow.cli import cli
from cath_alphaflow.io_utils import DecoratedCrhReader, Gene3DCrhReader

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

COUNT_UNIPROT_IDS = 5


def make_uniprot_id(n: int):
    return f"A000A0000{n}"


def make_cath_domain_id(n: int):
    letter = chr(96 + n)
    return f"{n}{letter * 3}{letter.upper()}00"


def make_af_chain_id(n: int):
    uniprot_id = make_uniprot_id(n)
    return f"AF-{uniprot_id}-F1-model_v4"


def make_md5(n: int):
    return f"{n}" * 32


def make_sfam_id(n: int):
    return f"1.1.1.{n}"


@pytest.fixture
def mock_uniprot_ids():
    return [make_uniprot_id(n) for n in range(1, COUNT_UNIPROT_IDS + 1)]


def mock_af_uniprot_md5_text():
    text = "\t".join(["af_chain_id", "uniprot_id", "sequence_md5"]) + "\n"
    for n in range(1, COUNT_UNIPROT_IDS + 1):
        text += "\t".join([make_af_chain_id(n), make_uniprot_id(n), make_md5(n)]) + "\n"

    return text


def mock_decorated_crh_text():
    # "1aaaA00","1.1.1.1","11111111111111111111111111111111","1aaaA00-i1","111.1","1-100","1-101","1-100","1.1e-111","1.2e-111",""
    text = ""
    for n in range(1, COUNT_UNIPROT_IDS + 1):
        cols = [
            make_cath_domain_id(n),
            make_sfam_id(n),
            make_md5(n),
            f"{n}{n}{n}.{n}",
            f"{n}-{n}00",
            f"{n}-{n}0{n}",
            f"{n}.{n}e-{n}{n}{n}",
            f"{n}.{n+1}e-{n}{n}{n}",
            "",
        ]
        text += ",".join([f'"{val}"' for val in cols]) + "\n"
    return text


def mock_gene3d_crh_text():
    # 3ce18771b4195d6aad287c3965a3c4f8        5ksdA03__3.40.50.1000/327-341_490-626   1054.6  327-341,490-626 327-341,490-626
    text = ""
    for n in range(1, COUNT_UNIPROT_IDS + 1):
        domain_id = make_cath_domain_id(n)
        sfam_id = make_sfam_id(n)
        cols = [
            make_md5(n),
            f"{domain_id}__{sfam_id}/{n}-{n}00_{n}0{n}-{n+1}0{n}",
            f"{n}{n}{n}.{n}",
            f"{n}-{n}00,{n}0{n}-{n+1}0{n}",
            f"{n}-{n}00,{n}0{n}-{n+1}0{n+1}",
            "",
        ]
        text += "\t".join(cols) + "\n"
    return text


@pytest.mark.parametrize(
    "subcommand,file_args,error",
    [
        (COMMAND_DB, (), "Error: Missing option '--dbname'"),
        (
            COMMAND_FILES,
            ("--src_gene3d_crh",),
            "Error: Missing option '--src_af_uniprot_md5'",
        ),
        (
            COMMAND_FILES,
            ("--src_decorated_crh",),
            "Error: Missing option '--src_af_uniprot_md5'",
        ),
        (
            COMMAND_FILES,
            ("--src_af_uniprot_md5",),
            "must specify at least one",
        ),
        (
            COMMAND_FILES,
            ("--src_af_uniprot_md5", "--src_gene3d_crh", "--src_decorated_crh"),
            "cannot specify more than one",
        ),
    ],
)
def test_create_dataset_unexpected_args_error(subcommand, file_args, error):
    """
    Checks we get sensible errors if we are missing required args
    """
    runner = CliRunner()
    with runner.isolated_filesystem():

        shared_args = []
        for opt in INPUT_OPTIONS + file_args:
            infilename = opt.replace("--", "") + ".txt"
            shared_args.extend([opt, infilename])
            Path(infilename).touch()
        for opt in OUTPUT_OPTIONS:
            shared_args.extend([opt, "fake_argument"])
        result = runner.invoke(cli, [subcommand, *shared_args])
        assert result.exit_code != 0
        assert error in "".join(result.output)


AF_UNIPROT_MD5_FILENAME = "af_uniprot_md5.csv"
GENE3D_CRH_FILENAME = "gene3d.crh"
DECORATED_CRH_FILENAME = "decorated.crh"


@pytest.mark.parametrize(
    "subcommand,extra_args,src_file,src_text_func",
    [
        (
            COMMAND_FILES,
            (
                "--src_af_uniprot_md5",
                AF_UNIPROT_MD5_FILENAME,
                "--src_gene3d_crh",
                GENE3D_CRH_FILENAME,
            ),
            GENE3D_CRH_FILENAME,
            mock_gene3d_crh_text,
        ),
        (
            COMMAND_FILES,
            (
                "--src_af_uniprot_md5",
                AF_UNIPROT_MD5_FILENAME,
                "--src_decorated_crh",
                DECORATED_CRH_FILENAME,
            ),
            DECORATED_CRH_FILENAME,
            mock_decorated_crh_text,
        ),
    ],
)
def test_create_dataset_crh_output(subcommand, extra_args, src_file, src_text_func):
    """
    Checks we get the same output for gene3d/decorated crh input
    """

    src_text = src_text_func()

    runner = CliRunner()
    with runner.isolated_filesystem():

        shared_args = []
        for opt in INPUT_OPTIONS:
            infilename = opt.replace("--", "") + ".txt"
            shared_args.extend([opt, infilename])
            Path(infilename).touch()
        expected_outfiles = []
        for opt in OUTPUT_OPTIONS:
            outfilename = opt.replace("--", "") + ".txt"
            expected_outfiles.extend([outfilename])
            shared_args.extend([opt, outfilename])

        with Path(AF_UNIPROT_MD5_FILENAME).open("wt") as fp:
            fp.write(mock_af_uniprot_md5_text())

        with Path(src_file).open("wt") as fp:
            fp.write(src_text)

        result = runner.invoke(cli, [subcommand, *shared_args, *extra_args])
        assert result.exit_code == 0

        annotations_path = Path("af_cath_annotations.txt")
        assert annotations_path.exists()
        with annotations_path.open("rt") as fp:
            annotation_lines = []
            for _line in fp:
                annotation_lines.extend([_line])

        if not len(annotation_lines) == COUNT_UNIPROT_IDS:
            print(f"OUTPUT: {''.join(result.output)}")
            print(f"LINES: {''.join(annotation_lines)}")

        assert len(annotation_lines) == COUNT_UNIPROT_IDS


@pytest.mark.parametrize(
    "subcommand,extra_args",
    [
        (COMMAND_DB, ["--dbname", "gene3d_21"]),
        (
            COMMAND_FILES,
            ["--src_af_uniprot_md5", "foo.csv", "--src_gene3d_crh", "bar.csv"],
        ),
        (
            COMMAND_FILES,
            ["--src_af_uniprot_md5", "foo.csv", "--src_decorated_crh", "bar.csv"],
        ),
    ],
)
def test_create_dataset_expected_args(subcommand, extra_args, create_mock_query):
    """
    Checks job runs to completion if we have expected args
    """
    runner = CliRunner()

    create_mock_query(["uniprot_id"], [])
    with runner.isolated_filesystem():

        for mock_file in ["foo.csv", "bar.csv"]:
            Path(mock_file).touch()

        shared_args = []
        for opt in INPUT_OPTIONS:
            infilename = opt.replace("--", "") + ".txt"
            shared_args.extend([opt, infilename])
            Path(infilename).touch()
        for opt in OUTPUT_OPTIONS:
            shared_args.extend([opt, "fake_argument"])
        result = runner.invoke(cli, [subcommand, *shared_args, *extra_args])
        assert result.exit_code == 0
        assert "DONE" in result.output


def test_create_dataset_db_command(create_mock_query):
    """
    Checks that the db mock works
    """

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
