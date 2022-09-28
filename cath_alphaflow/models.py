import logging
import re
from typing import List
from dataclasses import dataclass

from .errors import ParseError

RE_AF_CHAIN_ID = re.compile(
    r"AF-(?P<uniprot_acc>[0-9A-Z]+)-F(?P<frag_num>[0-9])-model_v(?P<version>[0-9]+)"
)

RE_AF_DOMAIN_ID = re.compile(
    r"AF-(?P<uniprot_acc>[0-9A-Z]+)-F(?P<frag_num>[0-9])-model_v(?P<version>[0-9]+)/(?P<chopping>[0-9\-_]+)"
)

LOG = logging.getLogger(__name__)


@dataclass
class PredictedCathDomain:
    """
    Holds data on a PredictedCathDomain (from Gene3D)
    """

    uniprot_acc: str
    sequence_md5: str
    gene3d_domain_id: str
    bitscore: float
    chopping: str


@dataclass
class Segment:
    start: str
    end: str


@dataclass
class Chopping:
    segments: List[Segment]

    RE_SEGMENT_SPLITTER = re.compile(r"[_,]")
    RE_SEGMENT_PARSER = re.compile(r"(?P<start>[0-9]+)-(?P<end>[0-9]+)")

    @classmethod
    def from_str(cls, chopping_str: str):
        segs = []
        for seg_str in re.split(cls.RE_SEGMENT_SPLITTER, chopping_str):
            match = cls.RE_SEGMENT_PARSER.match(seg_str)
            if not match:
                msg = f"failed to match segment '{seg_str}'"
                raise ParseError(msg)

            seg = Segment(start=int(match.group("start")), end=int(match.group("end")))
            segs.append(seg)
        return Chopping(segments=segs)

    def to_str(self):
        return "_".join([f"{seg.start}-{seg.end}" for seg in self.segments])


@dataclass
class AFChainID:
    uniprot_acc: str
    fragment_number: int
    version: int

    @classmethod
    def from_str(cls, raw_chainid: str):

        match = RE_AF_CHAIN_ID.match(raw_chainid)
        if not match:
            msg = f"failed to match AF chain id '{raw_chainid}'"
            raise ParseError(msg)

        chainid = AFChainID(
            uniprot_acc=match.group("uniprot_acc"),
            fragment_number=int(match.group("frag_num")),
            version=int(match.group("version")),
        )

        return chainid

    @property
    def af_chain_id(self):
        return f"AF-{self.uniprot_acc}-F{self.fragment_number}-model_v{self.version}"

    def to_str(self):
        return self.af_chain_id

    def __str__(self):
        return self.to_str()


@dataclass
class AFDomainID(AFChainID):

    chopping: Chopping

    @classmethod
    def from_str(cls, raw_domid: str):

        match = RE_AF_DOMAIN_ID.match(raw_domid)
        try:
            domid = AFDomainID(
                uniprot_acc=match.group("uniprot_acc"),
                fragment_number=int(match.group("frag_num")),
                version=int(match.group("version")),
                chopping=Chopping.from_str(match.group("chopping")),
            )
        except (KeyError, AttributeError):
            msg = f"failed to parse AFDomainId from {raw_domid}"
            LOG.error(msg)
            raise ParseError(msg)

        return domid

    @property
    def af_domain_id(self):
        return self.af_chain_id + "/" + self.chopping.to_str()

    def to_str(self):
        return self.af_domain_id
