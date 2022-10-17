import logging
import re
from typing import List
from dataclasses import dataclass, asdict

from .errors import ParseError
from .constants import DEFAULT_HELIX_MIN_LENGTH, DEFAULT_STRAND_MIN_LENGTH

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

    def deep_copy(self):
        return Segment(start=self.start, end=self.end)


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

    def deep_copy(self):
        new_segments = [s.deep_copy() for s in self.segments]
        return Chopping(segments=new_segments)


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

    def deep_copy(self):
        return AFChainID(asdict(self))


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

    def deep_copy(self):
        flds = asdict(self)
        flds["chopping"] = self.chopping.deep_copy()
        return AFDomainID(**flds)


@dataclass
class LURSummary:
    LUR_perc: float
    LUR_total: int
    residues_total: int


@dataclass
class SecStrSummary:
    af_domain_id: str
    ss_res_total: int
    res_count: int
    perc_not_in_ss: float
    sse_H_num: int
    sse_E_num: int

    @property
    def sse_num(self):
        return self.sse_E_num + self.sse_H_num

    def to_dict(self):
        d = self.__dict__
        d["sse_num"] = self.sse_num
        return d

    @classmethod
    def new_from_dssp_str(
        cls,
        dssp_str,
        acc_id,
        *,
        min_helix_length=DEFAULT_HELIX_MIN_LENGTH,
        min_strand_length=DEFAULT_STRAND_MIN_LENGTH,
    ):
        # Calculate percentage of residues in secondary structures
        ss_total = dssp_str.count("H") + dssp_str.count("E")
        domain_length = len(dssp_str)
        if domain_length == 0:
            msg = f"failed to find any SS data in DSSP string '{dssp_str}'"
            raise ParseError(msg)
        else:
            perc_not_in_ss = round(
                ((domain_length - ss_total) / domain_length) * 100, 2
            )

        # Calculate number of SSEs
        sse_H_num = sse_H_res = sse_E_num = sse_E_res = 0
        sse_H = sse_E = False

        # Calculate number of alpha helices and beta strands
        for residue in dssp_str:
            if residue == "H":
                sse_H_res += 1
                if sse_H_res >= min_helix_length and not sse_H:
                    sse_H = True
                    sse_H_num += 1

            if residue == "E":
                sse_E_res += 1
                if sse_E_res >= min_strand_length and not sse_E:
                    sse_E = True
                    sse_E_num += 1

            if residue != "H" and residue != "E":
                sse_H = sse_E = False
                sse_H_res = sse_E_res = 0

        ss_sum = SecStrSummary(
            af_domain_id=acc_id,
            ss_res_total=ss_total,
            res_count=domain_length,
            perc_not_in_ss=perc_not_in_ss,
            sse_H_num=sse_H_num,
            sse_E_num=sse_E_num,
        )

        return ss_sum


@dataclass
class pLDDTSummary:
    af_domain_id: str
    avg_plddt: float
    perc_LUR: float
    LUR_residues: int
    total_residues: int
