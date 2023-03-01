import logging
import tempfile

from cath_alphaflow.models import AFChainID, AFDomainID
from cath_alphaflow.io_utils import get_af_chain_id_reader, get_af_domain_id_reader

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
