import gzip
import logging
from pathlib import Path
import re
from typing import Callable, Union
from enum import Enum

from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB import PDBIO

from .models.domains import ChoppingSeqres, ChoppingPdbResLabel
from .errors import ChoppingError, MultipleModelsError, MultipleChainsError, ParseError

LOG = logging.getLogger(__name__)

RE_RESIDUE_MATCH = re.compile(r"(-?[0-9]+)([A-Z]?)")


class StructureFileType(Enum):
    PDB = "PDB"
    CIF = "CIF"


def default_map_to_pdb_resid(res):
    resnum, insert_code = RE_RESIDUE_MATCH.match(str(res)).groups()
    return (" ", int(resnum), " " if not insert_code else insert_code)


def chop_cif(
    *,
    domain_id: str,
    chain_cif_path: Path,
    domain_cif_path: Path,
    chopping: Union[ChoppingPdbResLabel, ChoppingSeqres],
    map_to_pdb_resid: Callable = default_map_to_pdb_resid,
):
    """
    Chops a CIF file into a domain based on the given chopping

    The optional `map_to_pdb_resid` parameter should provide a callable that accepts
    the start/end residue (`str`) from the chopping/segment and returns the equivalent
    residue name from the PDB/CIF structure. By default this is a 1:1 mapping, i.e. it
    is assumed that the start/end values in the chopping are equivalent to the residue
    labels in CIF.

    """

    chop_structure(
        domain_id=domain_id,
        chain_path=chain_cif_path,
        domain_path=domain_cif_path,
        chopping=chopping,
        map_to_pdb_resid=map_to_pdb_resid,
        file_type=StructureFileType.CIF,
    )


def chop_structure(
    *,
    domain_id: str,
    chain_path: Path,
    domain_path: Path,
    chopping: Union[ChoppingPdbResLabel, ChoppingSeqres],
    map_to_pdb_resid: Callable = default_map_to_pdb_resid,
    file_type: StructureFileType = StructureFileType.PDB,
):
    """
    Chops a PDB/CIF file into a domain based on the given chopping

    The optional `map_to_pdb_resid` parameter should provide a callable that accepts
    the start/end residue (`str`) from the chopping/segment and returns the equivalent
    residue name from the PDB/CIF structure. By default this is a 1:1 mapping, i.e. it
    is assumed that the start/end values in the chopping are equivalent to the residue
    labels in CIF.

    """

    if not file_type:
        if ".cif" in str(chain_path):
            file_type = StructureFileType.CIF
        elif ".pdb" in str(chain_path):
            file_type = StructureFileType.PDB
        else:
            raise ParseError(f"Unable to determine file type from {chain_path}")

    if file_type == StructureFileType.CIF:
        parser = MMCIFParser()
        io = MMCIFIO()
    elif file_type == StructureFileType.PDB:
        parser = PDBParser()
        io = PDBIO()
    else:
        raise ParseError(f"Unknown file type {file_type}")

    if str(chain_path).endswith(".gz"):
        with gzip.open(str(chain_path), mode="rt") as fp:
            structure = parser.get_structure(domain_id, fp)
    else:
        structure = parser.get_structure(domain_id, str(chain_path))

    models = structure.get_list()
    if len(models) != 1:
        msg = f"expected exactly 1 model, found {len(models)} in {chain_path}"
        raise MultipleModelsError(msg)
    model = models[0]

    chains = model.get_list()
    if len(chains) != 1:
        msg = f"expected exactly 1 chain, found {len(chains)} in {chain_path}"
        raise MultipleChainsError(msg)
    chain = chains[0]

    all_valid_resids = set()

    chopper = ChoppingProcessor(
        chopping,
        map_to_pdb_resid=map_to_pdb_resid,
        on_segment_residue=lambda res: all_valid_resids.add(res.get_id()),
    )
    chopper.process_chain(chain)

    if not all_valid_resids:
        msg = f"failed to find any valid residues when applying chopping {chopping} to chain {chain}"
        raise ChoppingError(msg)

    io.set_structure(structure)

    class ResSelector:
        def accept_model(self, _model):
            # LOG.info(f"accept_model: {_model} == {model} ({_model == model})")
            return 1 if _model == model else 0

        def accept_chain(self, _chain):
            # LOG.info(f"accept_chain: {_chain} == {chain} ({_chain == chain})")
            return 1 if _chain == chain else 0

        def accept_residue(self, _residue):
            return 1 if _residue.get_id() in all_valid_resids else 0

        def accept_atom(self, _atom):
            return 1

    res_selector = ResSelector()

    if str(domain_path).endswith(".gz"):
        with gzip.open(str(domain_path), mode="wt") as fp:
            io.save(fp, select=res_selector)
    else:
        io.save(str(domain_path), select=res_selector)


class ChoppingProcessor:
    """
    Applies a function to residues of `Bio.PDB.Chain` that fall within segments of a given chopping

    This is only necessary since we can't assume residue names within a PDB/CIF structure
    are numbered 1-n (residue labels are strings that often look like numbers but can be
    non-sequential and include optional insert codes).

    Also note that the chopping object can have a `map_to_uniprot_residue` function that maps
    the numbers in the chopping to the actual uniprot residue (e.g. fragment 2 of AF model is
    numbered 1-1400 in the PDB/CIF file, but this maps to UniProt resodies 201-1600)

    Typical usage example:

        # compile a list of all the residues within a domain "chopping"
        # where a chopping comprises of one or more segments, each with a start/end residue.

        domain_boundaries = "12-34A,56-78" # note the "insert code"
        domain_residues = []
        chopper = ChoppingProcessor(
            domain_boundaries,
            on_segment_residue=lambda res: domain_residues.append(res)
        )
        chopper.process_chain(structure.models[0].chains[0])

        # domain_residues = [12, 13, ..., 33, 34A, 56, 57, ..., 77, 78]

    """

    def __init__(
        self,
        chopping: Union[str, ChoppingPdbResLabel, ChoppingSeqres],
        on_segment_residue: Callable,
        *,
        map_to_pdb_resid: Callable = default_map_to_pdb_resid,
        chopping_class=ChoppingPdbResLabel,
    ):
        """
        Args:
            chopping: str|Chopping
            on_segment_residue: Callable
            map_to_pdb_resid: Callable = lambda x: (" ", x, " ")
        """

        if isinstance(chopping, str):
            chopping = chopping_class.from_str(chopping)
        self.chopping = chopping
        self.map_to_pdb_resid = map_to_pdb_resid
        self.on_segment_residue = on_segment_residue
        self.init()

    def init(self):
        self.next_segment_idx = None
        self.current_segment = None
        self.next_segment = None
        self.setup_next_segment()

    def setup_next_segment(self):
        if self.next_segment_idx is None:
            self.next_segment_idx = 0
        else:
            self.next_segment_idx += 1

        # if there are more segments to find ...
        if len(self.chopping.segments) >= self.next_segment_idx + 1:
            self.next_segment = self.chopping.segments[self.next_segment_idx]
        # otherwise stop looking ...
        else:
            self.next_segment_idx = None
            self.next_segment = None

    def enter_segment(self):
        self.current_segment = self.chopping.segments[self.next_segment_idx]

    def exit_segment(self):
        self.current_segment = None
        self.setup_next_segment()

    @property
    def current_start_resid(self):
        return self.map_to_pdb_resid(self.current_segment.start)

    @property
    def current_end_resid(self):
        return self.map_to_pdb_resid(self.current_segment.end)

    @property
    def next_start_resid(self):
        return self.map_to_pdb_resid(self.next_segment.start)

    @property
    def next_end_resid(self):
        return self.map_to_pdb_resid(self.next_segment.end)

    def process_chain(self, chain):
        self.init()
        for residue in chain.get_unpacked_list():
            residue_id = residue.id
            # if the chopping has set a mapping between PDB residue and chopping
            # (e.g due to AF fragments being offset), then apply that here
            if self.chopping.map_to_uniprot_residue:
                res_num = residue_id[1]
                uniprot_num = self.chopping.map_to_uniprot_residue(res_num)
                residue_id = (residue_id[0], uniprot_num, residue_id[2])

            # not currently in a segment and we don't have any more segments to look for
            if self.current_segment is None and self.next_segment is None:
                break

            # not currently in a segment, but we are looking for one ...
            if self.current_segment is None and self.next_segment is not None:
                # check for the start of the segment
                if residue_id == self.next_start_resid:
                    self.enter_segment()

            # in a segment
            if self.current_segment is not None:
                # do something with this residue
                self.on_segment_residue(residue)
                # check for the end of the segment
                if residue_id == self.current_end_resid:
                    self.exit_segment()
