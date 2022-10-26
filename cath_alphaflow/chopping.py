import logging
from pathlib import Path
import re
from typing import Callable

from Bio.PDB import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO

from .models import Chopping
from .errors import ChoppingError, MultipleModelsError, MultipleChainsError

LOG = logging.getLogger(__name__)

RE_RESIDUE_MATCH = re.compile(r"(-?[0-9]+)([A-Z]?)")


def default_map_to_pdb_resid(res):
    resnum, insert_code = RE_RESIDUE_MATCH.match(str(res)).groups()
    return (" ", int(resnum), " " if not insert_code else insert_code)


def chop_cif(
    *,
    domain_id: str,
    chain_cif_path: Path,
    domain_cif_path: Path,
    chopping: Chopping,
    map_to_pdb_resid: Callable = default_map_to_pdb_resid,
):
    """
    Chops a CIF file into a domain based on the given chopping

    The optional `map_to_pdb_resid` parameter should provide a callable that accepts
    the start/end residue (`str`) from the chopping/segment and returns the equivalent
    residue name from the PDB/CIF structure. By default this is a 1:1 mapping, i.e. it
    is assumed that the start/end in the chopping are based on PDB residue names.

    """

    parser = MMCIFParser()

    structure = parser.get_structure(domain_id, str(chain_cif_path))

    io = MMCIFIO()

    models = structure.get_list()
    if len(models) != 1:
        msg = f"expected exactly 1 model, found {len(models)} in {chain_cif_path}"
        raise MultipleModelsError(msg)
    model = models[0]

    chains = model.get_list()
    if len(chains) != 1:
        msg = f"expected exactly 1 chain, found {len(chains)} in {chain_cif_path}"
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

    class CifSelector:
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

    cif_selector = CifSelector()

    io.save(str(domain_cif_path), select=cif_selector)
    # io.save(str(domain_cif_path))


class ChoppingProcessor:
    """
    Applies a function to residues of `Bio.PDB.Chain` that fall within segments of a given chopping

    This is only necessary since we can't assume residue names within a PDB/CIF structure
    are numbered 1-n (residue labels are strings that often look like numbers but can be
    non-sequential and include optional insert codes).

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
        chopping: str | Chopping,
        on_segment_residue: Callable,
        *,
        map_to_pdb_resid: Callable = default_map_to_pdb_resid,
    ):
        """
        Args:
            chopping: str|Chopping
            on_segment_residue: Callable
            map_to_pdb_resid: Callable = lambda x: (" ", x, " ")
        """

        if isinstance(chopping, str):
            chopping = Chopping.from_str(chopping)
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

            # not currently in a segment and we don't have any more segments to look for
            if self.current_segment is None and self.next_segment is None:
                break

            # not currently in a segment, but we are looking for one ...
            if self.current_segment is None and self.next_segment is not None:
                # check for the start of the segment
                if residue.id == self.next_start_resid:
                    self.enter_segment()

            # in a segment
            if self.current_segment is not None:
                # do something with this residue
                self.on_segment_residue(residue)
                # check for the end of the segment
                if residue.id == self.current_end_resid:
                    self.exit_segment()
