import logging
import re
from typing import List, Callable
from dataclasses import dataclass, asdict
from pydantic import BaseModel as PyBaseModel, ConfigDict

from ..errors import ParseError, NoMatchingFragmentError
from ..constants import (
    DEFAULT_HELIX_MIN_LENGTH,
    DEFAULT_STRAND_MIN_LENGTH,
    AF_FRAGMENT_MAX_RESIDUES,
    AF_FRAGMENT_OVERLAP_WINDOW,
)

RE_UNIPROT_ID = re.compile(r"(?P<uniprot_acc>[0-9A-Z]{6}|[0-9A-Z]{10})$")

RE_AF_CHAIN_ID = re.compile(
    r"^AF-(?P<uniprot_acc>[0-9A-Z]+)-F(?P<frag_num>[0-9])-model_v(?P<version>[0-9]+)$"
)

RE_AF_DOMAIN_ID = re.compile(
    r"^(?P<raw_id>AF-(?P<uniprot_acc>[0-9A-Z]+)-F(?P<frag_num>[0-9])-model_v(?P<version>[0-9]+))[/\-](?P<chopping>[0-9\-_]+)$"
)

RE_UNIPROT_DOMAIN_ID = re.compile(
    r"^(?P<raw_id>(?P<uniprot_acc>[0-9A-Z]+))[/\-](?P<chopping>[0-9\-_]+)$"
)

LOG = logging.getLogger(__name__)


class BaseModel(PyBaseModel):
    model_config = ConfigDict(protected_namespaces=())


class CrhBase(BaseModel):
    """
    Defines the common interface for "original" and "decorated" CRH rows
    """

    domain_id: str
    superfamily_id: str
    sequence_md5: str
    model_id: str
    bitscore: float
    chopping_raw: str
    chopping_final: str


class CrhProvider(BaseModel):
    def to_crh(self):
        return CrhBase(
            domain_id=self.domain_id,
            superfamily_id=self.superfamily_id,
            sequence_md5=self.sequence_md5,
            model_id=self.model_id,
            bitscore=self.bitscore,
            chopping_raw=self.chopping_raw,
            chopping_final=self.chopping_final,
        )


class Gene3DCrh(CrhProvider):
    """
    Holds data corresponding to an entry from a Gene3D CRH file
    """

    # 3ce18771b4195d6aad287c3965a3c4f8        5ksdA01__1.20.1110.10/95-132_218-326_627-816    1054.6  95-132,218-326,627-816  95-132,218-326,627-816

    sequence_md5: str
    domain_sfam_id: str
    bitscore: float
    chopping_raw: str
    chopping_final: str

    @property
    def domain_id(self):
        return self.domain_sfam_id.split("__")[0]

    @property
    def model_id(self):
        return self.domain_id

    @property
    def superfamily_id(self):
        return self.domain_sfam_id.split("__")[1].split("/")[0]


class DecoratedCrh(CrhProvider):
    """
    Holds data corresponding to an entry in a 'Decorated' CATH Resolve Hits file
    """

    domain_id: str
    superfamily_id: str
    sequence_md5: str
    model_id: str
    bitscore: float
    chopping_raw: str
    chopping_final: str
    alignment_regs: str
    cond_evalue: float
    indp_evalue: float
    reg_ostats: str


@dataclass
class StatusLog:
    """
    Holds data corresponding to an entry in a Status Log file
    """

    entry_id: str
    status: str
    error: str
    description: str


class PredictedCathDomain(BaseModel):
    """
    Holds data on a PredictedCathDomain (from Gene3D)
    """

    uniprot_acc: str
    sequence_md5: str
    gene3d_domain_id: str
    bitscore: float
    chopping: str
    indp_evalue: float


@dataclass
class Segment:
    start: int
    end: int

    def deep_copy(self):
        return Segment(start=self.start, end=self.end)


@dataclass
class Chopping:
    segments: List[Segment]
    map_to_uniprot_residue: Callable = None

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

    def residue_count(self):
        res_count = sum([seg.end - seg.start + 1 for seg in self.segments])
        return res_count

    @property
    def first_residue(self):
        return self.segments[0].start

    @property
    def last_residue(self):
        return self.segments[-1].end

    def __str__(self):
        _map_str = None
        if self.map_to_uniprot_residue:
            _map_str = f"map_to_uniprot_residue(1 => {self.map_to_uniprot_residue(1)})"

        return f"<Chopping str='{self.to_str()}' {_map_str}>"


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
class GeneralDomainID:
    raw_id: str
    chopping: Chopping
    acc: str = None
    version: str = None

    @classmethod
    def from_uniprot_str(cls, raw_id: str, *, version: str = None):
        match = RE_UNIPROT_DOMAIN_ID.match(raw_id)
        if not match:
            msg = f"failed to parse AFDomainID from Uniprot Domain ID {raw_id}"
            LOG.error(msg)
            raise ParseError(msg)

        chopping = Chopping.from_str(match.group("chopping"))

        domid = GeneralDomainID(
            raw_id=match.group("raw_id"),
            acc=match.group("uniprot_acc"),
            version=version,
            chopping=chopping,
        )

        return domid

    def get_domain_residues_from_numeric_chopping(self) -> List[int]:
        chopping = self.chopping
        segments = chopping.segments
        domain_residues = []
        for segment in segments:
            for res in range(segment.start, segment.end):
                domain_residues.append(res)
        return domain_residues


@dataclass
class AFDomainID(AFChainID):
    chopping: Chopping

    @classmethod
    def from_uniprot_str(
        cls, raw_id: str, *, version: int, fragment_number: int = None
    ):
        match = RE_UNIPROT_DOMAIN_ID.match(raw_id)
        if not match:
            msg = f"failed to parse AFDomainID from Uniprot Domain ID {raw_id}"
            LOG.error(msg)
            raise ParseError(msg)

        chopping = Chopping.from_str(match.group("chopping"))

        # https://alphafold.ebi.ac.uk/faq
        # For human proteins longer than 2,700 amino acids, check the whole proteome download.
        # This contains longer proteins predicted as overlapping fragments. For example, Titin
        # has predicted fragment structures named as
        # Q8WZ42-F1 (residues 1–1400), Q8WZ42-F2 (residues 201–1600), etc.
        if fragment_number is None:
            # F1=1-1400, F2=201-1600, etc
            for _frag_num in range(1, 1000):
                # 1:0, 2:200, 3:400, ...
                offset_to_pdb = (_frag_num - 1) * AF_FRAGMENT_OVERLAP_WINDOW
                frag_window_start = offset_to_pdb + 1
                frag_window_end = frag_window_start + AF_FRAGMENT_MAX_RESIDUES - 1
                if (
                    chopping.last_residue >= frag_window_start
                    and chopping.last_residue <= frag_window_end
                ):
                    fragment_number = _frag_num
                    chopping.map_to_uniprot_residue = lambda x: x + offset_to_pdb
                    break
            else:
                msg = (
                    f"failed to find any potential AF fragments within chopping window "
                    f"{chopping.first_residue}-{chopping.last_residue} (chopping={chopping})"
                )
                raise NoMatchingFragmentError(msg)

        if fragment_number > 1:
            LOG.warning(
                f"using AF fragment {fragment_number} for chopping {chopping} (id: {raw_id})"
            )

        domid = AFDomainID(
            uniprot_acc=match.group("uniprot_acc"),
            fragment_number=fragment_number,
            version=version,
            chopping=chopping,
        )

        return domid

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

    @classmethod
    def from_foldseek_query(cls, raw_query_id: str):
        if raw_query_id.endswith(".cif"):
            raw_query_id = raw_query_id.replace(".cif", "")
        return cls.from_str(raw_query_id)

    @property
    def af_domain_id(self):
        return self.af_chain_id + "/" + self.chopping.to_str()

    def to_str(self):
        return self.af_domain_id

    def to_file_stub(self):
        return self.af_domain_id.replace("/", "-")

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
class FoldseekSummary:
    query: str
    target: str
    qstart: int
    qend: int
    qlen: int
    tstart: int
    tend: int
    tlen: int
    qcov: float
    tcov: float
    bits: int
    evalue: float


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
        dssp_str: str,
        acc_id: str,
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
    residues_total: int
