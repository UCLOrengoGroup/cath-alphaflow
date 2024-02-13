from Bio.PDB import Structure, MMCIFParser, Selection, NeighborSearch
import logging
from pathlib import Path
import gzip
import numpy as np
import click

from cath_alphaflow.io_utils import get_csv_dictreader
from cath_alphaflow.io_utils import get_csv_dictwriter
from cath_alphaflow.io_utils import get_pdb_structure
from cath_alphaflow.models.domains import GeneralDomainID
from cath_alphaflow.models.domains import ChoppingPdbResLabel
from cath_alphaflow.models.domains import SegmentStr

from cath_alphaflow.constants import DEFAULT_GLOB_DISTANCE, DEFAULT_GLOB_VOLUME

DEFAULT_PDB_SUFFIX = ".pdb"

LOG = logging.getLogger()

# Path: cath_alphaflow/commands/measure_globularity.py

# --consensus_domain_list /cluster/project9/afdb_domain_ext/results/proteome_consensus_domains.domain_summary.tsv
# ==


@click.command()
@click.option(
    "--consensus_domain_list",
    type=click.File("rt"),
    required=False,
    help="Input: CSV file containing consensus domains",
)
@click.option(
    "--chainsaw_domain_list",
    type=click.File("rt"),
    required=False,
    help="Input: CSV file containing chainsaw domain choppings",
)
@click.option(
    "--pdb_dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, resolve_path=True),
    required=True,
    help="Input: directory of PDB files",
)
@click.option(
    "--gzipped",
    "chains_are_gzipped",
    type=bool,
    default=False,
    help="Whether chain files are stored in .gz files",
)
@click.option(
    "--pdb_suffix",
    type=str,
    default=DEFAULT_PDB_SUFFIX,
    help="Suffix to use when looking for the PDB model file",
)
@click.option(
    "--domain_globularity",
    type=click.File("wt"),
    required=True,
    help="Output: CSV file for domain results including both globularity parameters",
)
@click.option(
    "--distance_cutoff",
    type=int,
    required=False,
    default=DEFAULT_GLOB_DISTANCE,
    help=f"Distance cutoff for the packing density in Angstrom. (default: {DEFAULT_GLOB_DISTANCE})",
)
@click.option(
    "--volume_resolution",
    type=int,
    required=False,
    default=DEFAULT_GLOB_VOLUME,
    help=f"The voxel resolution for approximating the protein volume. (default: {DEFAULT_GLOB_VOLUME})",
)
def measure_globularity(
    consensus_domain_list,
    chainsaw_domain_list,
    pdb_dir,
    pdb_suffix,
    domain_globularity,
    distance_cutoff,
    volume_resolution,
    chains_are_gzipped,
):
    "Checks the globularity of the AF domain"

    if consensus_domain_list and chainsaw_domain_list:
        msg = "Cannot specify both --consensus_domain_list and --chainsaw_domain_list"
        raise click.BadParameter(msg)

    if consensus_domain_list:
        general_domain_provider = yield_domain_from_consensus_domain_list(
            consensus_domain_list
        )
    elif chainsaw_domain_list:
        general_domain_provider = yield_domain_from_chainsaw_domain_list_csv(
            chainsaw_domain_list
        )
    else:
        general_domain_provider = yield_domain_from_pdbdir(pdb_dir)

    globularity_writer = get_csv_dictwriter(
        domain_globularity,
        fieldnames=[
            "model_id",
            "chopping",
            "packing_density",
            "normed_radius_gyration",
        ],
    )
    globularity_writer.writeheader()

    click.echo(
        f"Checking domain globularity "
        f"(model_dir={pdb_dir}, pdb_suffix={pdb_suffix}, out_file={domain_globularity.name} ) ..."
    )

    for domain_id in general_domain_provider:
        LOG.info(f"Working on: {domain_id} ...")

        model_structure = get_pdb_structure(
            domain_id, pdb_dir, chains_are_gzipped, pdb_suffix=pdb_suffix
        )

        domain_packing_density = calculate_packing_density(
            domain_id, model_structure, distance_cutoff
        )

        domain_normed_radius_gyration = calculate_normed_radius_of_gyration(
            domain_id, model_structure, volume_resolution
        )
        LOG.info(
            f"Processed entry: {domain_id}\tpacking_density:{domain_packing_density}\tgyration:{domain_normed_radius_gyration}"
        )
        globularity_writer.writerow(
            {
                "model_id": domain_id.raw_id,
                "chopping": domain_id.chopping.to_str() if domain_id.chopping else "-",
                "packing_density": domain_packing_density,
                "normed_radius_gyration": domain_normed_radius_gyration,
            }
        )

    click.echo("DONE")


def check_domain_chopping_matches_model_residues(
    domain_id: GeneralDomainID, model_structure: Structure
):
    """
    Make sure the domain chopping matches the residues in the structure
    """

    chopping_residue_labels = set(domain_id.get_domain_residues_from_numeric_chopping())
    structure_residue_labels = set(
        [res.id[1] for res in model_structure.get_residues()]
    )

    if chopping_residue_labels != structure_residue_labels:
        raise Exception(
            f"Residues in chopping ({len(chopping_residue_labels)}) and structure ({len(structure_residue_labels)}) do not match "
            f"(mismatch={structure_residue_labels-chopping_residue_labels})"
        )


def calculate_packing_density(
    domain_id: GeneralDomainID,
    model_structure: Structure,
    distance_cutoff: int,
    include_all_atoms: bool = True,
) -> float:
    chain_residues = list(model_structure.get_chains())[0].get_residues()
    if domain_id.chopping:
        target_residues = domain_id.chopping.filter_bio_residues(chain_residues)
    else:
        target_residues = chain_residues

    # List of hydrophobic residues
    hydrophobic_list = [
        "GLY",
        "ALA",
        "LEU",
        "ILE",
        "VAL",
        "PRO",
        "PHE",
        "MET",
        "TRP",
    ]

    neighbor_list_protein = []
    if include_all_atoms:
        # nearby residues for all hydrophobic residues in the structure
        atom_list = [atom for atom in model_structure.get_atoms()]
    else:
        # nearby residues for all hydrophobic residues in the domain
        atom_list = [atom for res in target_residues for atom in res.get_atoms()]
    for res in target_residues:
        if res.get_resname() in hydrophobic_list and not res.is_disordered():
            center_atoms = Selection.unfold_entities(res, "A")
            ns = NeighborSearch(atom_list)

            nearby_residues = {
                _res
                for center_atom in center_atoms
                for _res in ns.search(center_atom.coord, distance_cutoff, "R")
            }
            neighbor_list_protein.append(len(nearby_residues))

    return round(np.mean(neighbor_list_protein), 3)


# Function to get the approximated Volume of the domain
def get_volume(residues, volume_resolution):
    reduced_coords = list()
    for residue in residues:
        for atm in residue.get_atoms():
            x_reduced = int(atm.coord.tolist()[0] / volume_resolution)
            y_reduced = int(atm.coord.tolist()[1] / volume_resolution)
            z_reduced = int(atm.coord.tolist()[2] / volume_resolution)
            reduced_coords.append([x_reduced, y_reduced, z_reduced])

    unique_coords = [list(x) for x in set(tuple(x) for x in reduced_coords)]
    volume = len(unique_coords) * volume_resolution**3

    return volume


# Function to get the normed radius of gyration
def calculate_normed_radius_of_gyration(
    domain_id: GeneralDomainID,
    model_structure: Structure,
    volume_resolution: int,
) -> float:
    chain_residues = list(model_structure.get_chains())[0].get_residues()
    if domain_id.chopping:
        target_residues = domain_id.chopping.filter_bio_residues(chain_residues)
    else:
        target_residues = list(chain_residues)

    coords = []
    masses = []
    # Fill lists for the coords and masses for all atoms in the domain
    for residue in target_residues:
        for atom in residue:
            coords.append(atom.coord.tolist())
            if atom.element == "C":
                masses.append(12.0107)
            elif (
                atom.element == "H"
            ):  # We do not consider hydrogens for the radius of gyration
                masses.append(0.0)
            elif atom.element == "O":
                masses.append(15.9994)
            elif atom.element == "N":
                masses.append(14.0067)
            elif atom.element == "S":
                masses.append(32.065)
            else:
                raise Exception(f"Protein contains unknown atom type {atom.element}")

    # Calculate the radius of gyration
    mass_coords = [(m * i, m * j, m * k) for (i, j, k), m in zip(coords, masses)]
    total_mass = sum(masses)
    r_tmp = sum(
        mi * i + mj * j + mk * k for (i, j, k), (mi, mj, mk) in zip(coords, mass_coords)
    )
    m_tmp = sum((sum(i) / total_mass) ** 2 for i in zip(*mass_coords))
    radius_gyration = np.sqrt(r_tmp / total_mass - m_tmp)

    # Calculate the radius of gyration for a sphere with the same volume as the protein
    volume = get_volume(target_residues, volume_resolution)
    radius = np.cbrt(volume) * 3 / 4 * np.pi
    radius_gyration_sphere = np.sqrt(3 / 5 * radius**2)

    # Return normed radius of gyration
    return round(radius_gyration / radius_gyration_sphere, 3)


def guess_chopping_from_pdb_file(
    pdbfile, target_chain_id=None, assume_all_atom_breaks_are_segments=True
) -> ChoppingPdbResLabel:
    """
    Reverse engineer a `Chopping` object from a PDB file

    'Real' PDB files can have weird residue numberings, so it is not possible to
    reliably know whether a PDB file has been chopped by an external process
    (i.e. a CATH domain) or whether the numbering in the original PDB is just weird.
    The best we can do is guess, i.e. check for cases where the sequential atom number
    is not contiguous.

    If `assume_all_atom_breaks_are_segments` is turned off then we just create a chopping
    with a single segment that spans the entire PDB file.
    """

    pdbpath = Path(str(pdbfile))
    segments = []
    start_res = None
    end_res = None
    last_atom_num = None
    last_res_label = None
    with pdbpath.open("rt") as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue

            atom_num = int(line[6:11])
            chain_id = line[21]
            res_label = line[22:27].strip()  # includes the (optional) insert code

            # if we haven't specified a chain, then we take the first one we see
            if not target_chain_id:
                target_chain_id = chain_id

            if not last_atom_num:
                last_atom_num = atom_num
            if not last_res_label:
                last_res_label = res_label

            if chain_id != target_chain_id:
                continue

            if not start_res:
                start_res = res_label
            end_res = res_label

            if assume_all_atom_breaks_are_segments:
                if res_label != last_res_label and atom_num != last_atom_num + 1:
                    segments.append(SegmentStr(start=start_res, end=last_res_label))
                    start_res = None
                    end_res = None

            last_res_label = res_label
            last_atom_num = atom_num

    if start_res and end_res:
        segments.append(SegmentStr(start=start_res, end=end_res))

    return ChoppingPdbResLabel(segments=segments)


def yield_domain_from_pdbdir(pdbdir, pdb_suffix=".pdb") -> GeneralDomainID:
    pdbdir = Path(str(pdbdir))
    for pdbpath in pdbdir.iterdir():
        if not str(pdbpath).endswith(pdb_suffix):
            continue
        model_id = pdbpath.stem
        yield GeneralDomainID(raw_id=model_id, chopping=None)


def yield_domain_from_consensus_domain_list(consensus_domain_list) -> GeneralDomainID:
    consensus_domain_list_reader = get_csv_dictreader(
        consensus_domain_list,
        fieldnames=[
            "domain_id",
            "domain_md5",
            "consensus_level",
            "chopping",
            "domain_nres",
            "no_segments",
            "plddt",
            "total_ss",
            "no_helix",
            "no_strand",
            "no_helix+strand",
            "no_turn",
            "proteome_id",
        ],
    )

    for domain_row in consensus_domain_list_reader:
        domain = GeneralDomainID(
            raw_id=domain_row["domain_id"],
            chopping=ChoppingPdbResLabel.from_str(domain_row["chopping"]),
        )
        yield domain


def yield_domain_from_chainsaw_domain_list_csv(
    consensus_domain_list,
) -> GeneralDomainID:
    """

    ```
    $ head combined_chainsaw_domains.tsv
    chain_id	sequence_md5	nres	ndom	chopping	uncertainty
    5yntA	5e88befc4122daa13edac258bd490776	300	1	3-292	0.00836
    5yh0J	dcb023981e453a818025bc8bc1789d51	378	2	135-366	0.171
    5yh0J	dcb023981e453a818025bc8bc1789d51	378	2	377-547	0.171
    ```
    """

    fieldnames = [
        "chain_id",
        "sequence_md5",
        "nres",
        "ndom",
        "chopping",
        "uncertainty",
    ]
    chainsaw_domain_list_reader = get_csv_dictreader(
        consensus_domain_list,
        fieldnames=fieldnames,
    )

    headers = next(chainsaw_domain_list_reader)

    assert list(headers.values()) == list(
        fieldnames
    ), f"Unexpected headers: {list(headers.values())} (expected {fieldnames})"

    domains_by_chain_id = {}

    for domain_row in chainsaw_domain_list_reader:
        chain_id = domain_row["chain_id"]
        if chain_id not in domains_by_chain_id:
            domains_by_chain_id[chain_id] = []
        domain_count = len(domains_by_chain_id[chain_id])
        domain_id = f"{chain_id}{domain_count:02}"
        domains_by_chain_id[chain_id].append(domain_id)
        domain = GeneralDomainID(
            # raw_id=domain_id,
            raw_id=chain_id,
            chopping=ChoppingPdbResLabel.from_str(domain_row["chopping"]),
        )
        yield domain
