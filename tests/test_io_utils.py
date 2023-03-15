import logging
import tempfile

from cath_alphaflow.models.domains import AFChainID, AFDomainID
from cath_alphaflow.io_utils import get_af_chain_id_reader, get_af_domain_id_reader
from cath_alphaflow.io_utils import Gene3DCrhReader, DecoratedCrhReader
from cath_alphaflow.models.domains import Gene3DCrh, DecoratedCrh, CrhBase

LOG = logging.getLogger(__name__)


def test_af_chain_id_reader():
    """
    Check that af_chain_id_reader returns AFChainID objects
    """

    expected_af_chain_ids = ["AF-P00520-F1-model_v3", "AF-P00521-F1-model_v3"]
    expected_af_chains = [
        AFChainID.from_str(chain_id) for chain_id in expected_af_chain_ids
    ]

    with tempfile.NamedTemporaryFile(mode="wt") as tmpfh:
        tmpfh.write("af_chain_id\n")
        for af_id in expected_af_chain_ids:
            tmpfh.write(af_id + "\n")
        tmpfh.flush()

        with open(tmpfh.name, "rt") as new_fh:
            reader = get_af_chain_id_reader(new_fh)
            af_chain_ids = list(reader)
            assert af_chain_ids == expected_af_chains


def test_af_domain_id_reader():
    """
    Check that af_domain_id_reader returns AFDomainID objects
    """

    expected_af_domain_ids = [
        "AF-P00520-F1-model_v3/12-345",
        "AF-P00521-F1-model_v3/234-456",
    ]
    expected_af_domains = [
        AFDomainID.from_str(domain_id) for domain_id in expected_af_domain_ids
    ]

    with tempfile.NamedTemporaryFile(mode="wt") as tmpfh:
        tmpfh.write("af_domain_id\n")
        for af_id in expected_af_domain_ids:
            tmpfh.write(af_id + "\n")
        tmpfh.flush()

        with open(tmpfh.name, "rt") as new_fh:
            reader = get_af_domain_id_reader(new_fh)
            af_domain_ids = list(reader)
            assert af_domain_ids == expected_af_domains


def test_gene3d_crh_reader():

    test_crh = (
        """
000122ad8c8fccfd2991bbd4a138d3c6	1c52A00__1.10.760.10/18-148	146.2	18-148	18-148
00015352b8446c4ccdb74461696a601f	6dxpC00__3.40.50.360/1-202	202.1	1-202	1-202
0002fb2a82c8a28bb6119ed72997a450	2eslA00__2.40.100.10/12-200	268.8	12-200	12-200
""".strip()
        + "\n"
    )

    with tempfile.NamedTemporaryFile(mode="wt") as tmpfh:
        tmpfh.write(test_crh)
        tmpfh.flush()

        with open(tmpfh.name, "rt") as new_fh:
            reader = Gene3DCrhReader(new_fh)
            crh_entries = list(reader)
            assert len(crh_entries) == 3
            assert isinstance(crh_entries[0], Gene3DCrh)
            crh_entries[0].sequence_md5 = "000122ad8c8fccfd2991bbd4a138d3c6"

            crh = crh_entries[0].to_crh()
            assert isinstance(crh, CrhBase)
            crh.domain_id = "1c52A00"
            crh.superfamily_id = "1.10.760.10"


def test_decorated_crh_reader():

    test_crh = (
        """
"1c52A00","1.10.760.10","000122ad8c8fccfd2991bbd4a138d3c6","1c52A00-i2","146.2","18-148","18-148","1-129,18-146","1.1e-44","1.5e-38",""
"6dxpC00","3.40.50.360","00015352b8446c4ccdb74461696a601f","6dxpC00-i2","202.1","1-202","1-202","2-61,2-61;62-147,63-148;148-197,152-201","1e-61","1.3e-55",""
"2eslA00","2.40.100.10","0002fb2a82c8a28bb6119ed72997a450","2eslA00-i2","268.8","12-200","12-200","7-186,17-196","1.2e-82","3.9e-76",""
""".strip()
        + "\n"
    )

    with tempfile.NamedTemporaryFile(mode="wt") as tmpfh:
        tmpfh.write(test_crh)
        tmpfh.flush()

        with open(tmpfh.name, "rt") as new_fh:
            reader = DecoratedCrhReader(new_fh)
            crh_entries = list(reader)
            assert len(crh_entries) == 3
            assert isinstance(crh_entries[0], DecoratedCrh)
            crh_entries[0].sequence_md5 = "000122ad8c8fccfd2991bbd4a138d3c6"

            crh = crh_entries[0].to_crh()
            assert isinstance(crh, CrhBase)
            crh.domain_id = "1c52A00"
            crh.superfamily_id = "1.10.760.10"
