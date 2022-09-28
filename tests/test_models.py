import pytest

from cath_alphaflow.errors import ParseError
from cath_alphaflow.models import Segment, Chopping, AFChainID, AFDomainID


def test_af_ids():

    chain_id = "AF-P00520-F1-model_v3"

    chain = AFChainID.from_str(chain_id)
    assert f"{chain}" == chain_id
    assert chain.uniprot_acc == "P00520"
    with pytest.raises(AttributeError):
        chain.chopping  # only a domain has a chopping
    del chain

    dom_id = f"{chain_id}/12-234"

    dom = AFDomainID.from_str(dom_id)
    assert f"{dom}" == dom_id
    assert dom.af_domain_id == dom_id
    assert dom.af_chain_id == chain_id
    assert dom.chopping.segments == [Segment(start=12, end=234)]
    del dom


def test_chopping_parser():

    assert Chopping.from_str("123-456") == Chopping(
        segments=[Segment(start=123, end=456)]
    )

    assert Chopping.from_str("123-456_789-1230") == Chopping(
        segments=[Segment(start=123, end=456), Segment(start=789, end=1230)],
    )


def test_af_chain_parser():

    assert AFChainID.from_str("AF-P00520-F1-model_v3") == AFChainID(
        uniprot_acc="P00520", fragment_number=1, version=3
    )

    bad_ids = [
        "P00520-F1-v3",
        "AF_P00520_F1_v3",
        "AF-P00520_-F1-v",
        "AF-P00520-A-v3",
        "AF-P00520-F1-v",
    ]

    for bad_id in bad_ids:
        with pytest.raises(ParseError):
            AFChainID.from_str(bad_id)


def test_af_domain_parser():

    assert AFDomainID.from_str("AF-P00520-F1-model_v3/12-234") == AFDomainID(
        uniprot_acc="P00520",
        fragment_number=1,
        version=3,
        chopping=Chopping(segments=[Segment(start=12, end=234)]),
    )

    bad_ids = [
        "AF-P00520-F1-v3/12_234",
        "AF-P00520-F1-v3/12-234_234",
        "AF-P00520-F1-v3/-12-234",
        "AF-P00520-F1-v3/12A-234",
        "AF-P00520-F1-v3/12A-234",
    ]

    for bad_id in bad_ids:
        with pytest.raises(ParseError):
            AFDomainID.from_str(bad_id)
