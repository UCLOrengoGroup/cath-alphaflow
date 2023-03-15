import logging
import shutil
import filecmp
import difflib

LOG = logging.getLogger(__name__)


def write_uniprot_ids(fh, uniprot_ids):
    fh.write("uniprot_id\n")
    for uniprot_id in uniprot_ids:
        fh.write(uniprot_id + "\n")


def assert_files_match(*, test, expected, bootstrap_results: bool = False):
    if bootstrap_results:
        LOG.warning(f"BOOTSTRAP test results: {test} -> {expected}")
        shutil.copyfile(str(test), str(expected))

    files_match = filecmp.cmp(str(test), str(expected))
    if not files_match:
        with test.open("rt") as fp:
            test_lines = [line.rstrip() for line in fp]
        with expected.open("rt") as fp:
            expected_lines = [line.rstrip() for line in fp]

        LOG.error(
            f"Found differences between test result file ({test}) and expected file ({expected})"
        )
        for diffline in difflib.unified_diff(
            test_lines, expected_lines, fromfile=str(test), tofile=str(expected)
        ):
            LOG.error(diffline)

    assert files_match
