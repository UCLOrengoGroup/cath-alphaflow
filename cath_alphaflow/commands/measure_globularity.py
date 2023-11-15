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
from cath_alphaflow.models.domains import Chopping

from cath_alphaflow.constants import DEFAULT_GLOB_DISTANCE, DEFAULT_GLOB_VOLUME

LOG = logging.getLogger()

# Path: cath_alphaflow/commands/measure_globularity.py

# --consensus_domain_list /cluster/project9/afdb_domain_ext/results/proteome_consensus_domains.domain_summary.tsv
# ==


@click.command()
@click.option(
    "--consensus_domain_list",
    type=click.File("rt"),
    required=True,
    help="Input: CSV file containing consensus domains",
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
    pdb_dir,
    domain_globularity,
    distance_cutoff,
    volume_resolution,
    chains_are_gzipped,
):
    "Checks the globularity of the AF domain"

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
        f"Checking globularity for AF domain"
        f"(model_dir={pdb_dir}, in_file={consensus_domain_list.name}, "
        f"out_file={domain_globularity.name} ) ..."
    )

    for domain_row in consensus_domain_list_reader:
        domain_id = GeneralDomainID(
            raw_id=domain_row["domain_id"],
            chopping=Chopping.from_str(domain_row["chopping"]),
        )

        LOG.debug(f"Working on: {domain_id} ...")

        model_structure = get_pdb_structure(domain_id, pdb_dir, chains_are_gzipped)

        check_domain_chopping_matches_model_residues(domain_id, model_structure)

        domain_packing_density = calculate_packing_density(
            domain_id, model_structure, distance_cutoff
        )

        domain_normed_radius_gyration = calculate_normed_radius_of_gyration(
            domain_id, model_structure, volume_resolution
        )
        LOG.debug(
            f"Processed entry: {domain_id}\tpacking_density:{domain_packing_density}\tgiration:{domain_normed_radius_gyration}"
        )
        globularity_writer.writerow(
            {
                "model_id": domain_id.raw_id,
                "chopping": domain_id.chopping,
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
            f"Residues in chopping ({len(chopping_residue_labels)}) and structure ({len(structure_residue_labels)}) do not match"
        )


def calculate_packing_density(
    domain_id: GeneralDomainID,
    chain_structure: Structure,
    distance_cutoff: int,
) -> float:
    domain_residue_numbers = domain_id.get_domain_residues_from_numeric_chopping()

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
    # Get the number of nearby residues for all hydrophobic residues in the domain
    atom_list = [atom for atom in chain_structure.get_atoms()]
    for residue_number in domain_residue_numbers:
        residue = chain_structure[0]["A"][residue_number]
        if residue.get_resname() in hydrophobic_list and not residue.is_disordered():
            center_atoms = Selection.unfold_entities(residue, "A")
            ns = NeighborSearch(atom_list)

            nearby_residues = {
                residue
                for center_atom in center_atoms
                for residue in ns.search(center_atom.coord, distance_cutoff, "R")
            }
            neighbor_list_protein.append(len(nearby_residues))

    return round(np.mean(neighbor_list_protein), 3)


# Function to get the approximated Volume of the domain
def get_volume(structure, domain_residue_numbers, volume_resolution):
    reduced_coords = list()
    for residue_number in domain_residue_numbers:
        residue = structure[0]["A"][residue_number]
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
    domain_residue_numbers = domain_id.get_domain_residues_from_numeric_chopping()

    coords = []
    masses = []
    # Fill lists for the coords and masses for all atoms in the domain
    for residue_number in domain_residue_numbers:
        residue = model_structure[0]["A"][residue_number]

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
    volume = get_volume(model_structure, domain_residue_numbers, volume_resolution)
    radius = np.cbrt(volume) * 3 / 4 * np.pi
    radius_gyration_sphere = np.sqrt(3 / 5 * radius**2)

    # Return normed radius of gyration
    return round(radius_gyration / radius_gyration_sphere, 3)
