import gzip
import io
from pathlib import Path
import pytest

from cath_alphaflow.commands.load_mongo import (
    get_beacons_uniprot_summary_from_af_cif,
    beacons,
)

CIF_DIR = Path(__file__).parent / "fixtures" / "cif"


@pytest.mark.parametrize(
    "example_cif_fname,uniprot_id",
    [
        ["AF-P00520-F1-model_v3.cif.gz", "P00520"],
        ["AF-Q15772-F3-model_v4.cif.gz", "Q15772"],
        ["AF-Q15772-F11-model_v4.cif.gz", "Q15772"],
    ],
)
def test_get_beacons_uniprot_summary_from_af_cif(example_cif_fname, uniprot_id):

    cif_path = CIF_DIR / example_cif_fname
    with gzip.open(str(cif_path), "rt") as fh:
        uniprot_summary = get_beacons_uniprot_summary_from_af_cif(fh)
        assert isinstance(uniprot_summary, beacons.UniprotSummary)
        assert uniprot_summary.uniprot_entry.ac == uniprot_id
