from Bio.PDB import MMCIFParser, Selection, NeighborSearch
import logging
from pathlib import Path
import gzip
import numpy as np

import click


from cath_alphaflow.io_utils import get_af_domain_id_reader
from cath_alphaflow.io_utils import get_csv_dictwriter
from cath_alphaflow.models import AFDomainID, Segment, Chopping
from cath_alphaflow.errors import NoMatchingResiduesError

LOG = logging.getLogger()

# Input:
# - `FILE:AF_CHAIN_MMCIF`
# - `FILE:AF_DOMAIN_LIST`

# Output:
# - `FILE:AF_DOMAIN_LIST_POST_TAILCHOP` -- `(AF_domain_ID_tailchop)`
# - `FILE:AF_DOMAIN_MAPPING_POST_TAILCHOP` -- `(AF_domain_ID_orig, AF_domain_ID_tailchop)`


@click.command()
@click.option(
    "--af_domain_list",
    type=click.File("rt"),
    required=True,
    help="Input: CSV file containing AF2 domains",
)
@click.option(
    "--af_chain_mmcif_dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, resolve_path=True),
    required=True,
    help="Input: directory of mmCIF files",
)
@click.option(
    "--af_domain_globularity",
    type=click.File("wt"),
    required=True,
    help="Output: CSV file for AF2 domain list including both globularity parameters",
)
# @click.option(
#     "--af_domain_mapping_post_tailchop",
#     type=click.File("wt"),
#     required=True,
#     help="Output: CSV file for mapping of AF2 domain before/after chopping",
# )
@click.option(
    "--distance_cutoff",
    type=int,
    required=False,
    default=5,
    help="Distance cutoff for the packing density in Angstrom. If none is given, a cutoff of 5 Angstrom will be applied",
)
@click.option(
    "--volume_resolution",
    type=int,
    required=False,
    default=5,
    help="The voxel resolution for approximating the protein volume. If none is given, a resolution of 5 Angstrom will be applied",
)
def measure_globularity(
    af_domain_list,
    af_chain_mmcif_dir,
    af_domain_globularity,
    # af_domain_mapping_post_tailchop,
    distance_cutoff,
    volume_resolution,
):
    "Checks the globularity of the AF domain"

    af_domain_list_reader = get_af_domain_id_reader(af_domain_list)

    # af_domain_globularity_writer = get_csv_dictwriter(
        # af_domain_list_post_tailchop, fieldnames=["af_domain_id"]
    # )
    # af_domain_list_post_tailchop_writer.writeheader()

    af_globularity_writer = get_csv_dictwriter(
        af_domain_globularity,
        fieldnames=["af_domain_id", "packing_density", "normed_radius_gyration"],
    )
    af_globularity_writer.writeheader()

    click.echo(
        f"Checking globularity for AF domain"
        f"(mmcif_dir={af_chain_mmcif_dir}, in_file={af_domain_list.name}, "
        f"out_file={af_domain_globularity.name} ) ..."
    )

    for af_domain_id in af_domain_list_reader:
        click.echo(f"Working on: {af_domain_id} ...")

        af_domain_packing_density = calculate_packing_density(
            af_domain_id, af_chain_mmcif_dir, distance_cutoff
        )
        print(af_domain_packing_density)

        af_domain_normed_radius_gyration = calculate_normed_radius_of_gyration(
            af_domain_id, af_chain_mmcif_dir,volume_resolution
        )

        print(af_domain_normed_radius_gyration)
        # af_domain_list_post_tailchop_writer.writerow(
            # {"af_domain_id": af_domain_packing_density}
        # )

        af_globularity_writer.writerow(
            {
                "af_domain_id": af_domain_id,
                "packing_density": af_domain_packing_density,
                "normed_radius_gyration": af_domain_normed_radius_gyration,
            }
        )

    click.echo("DONE")


def get_local_plddt_for_res(
    structure,
    residue_num: int,
    *,
    model_num: int = 0,
    chain_id: str = "A",
    atom_type: str = "CA",
):
    """
    Returns the local pLDDT score for a given residue
    """
    return structure[model_num][chain_id][residue_num][atom_type].get_bfactor()



def calculate_packing_density(
    af_domain_id: AFDomainID,
    af_chain_mmcif_dir: Path,
    distance_cutoff: int,
    *,
    cif_filename=None,
) -> AFDomainID:

    chopping = af_domain_id.chopping
    segments = chopping.segments
    domain_residues = []
    for segment in segments:
        for res in range(segment.start,segment.end):
            domain_residues.append(res)


    # create default filename
    if cif_filename is None:
        cif_filename = af_domain_id.af_chain_id + ".cif.gz"

    cif_path = Path(af_chain_mmcif_dir, cif_filename)

    with gzip.open(f"{cif_path}", mode="rt") as cif_fh:

        structure = MMCIFParser(QUIET=1).get_structure("struc", cif_fh)

        # List of hydrophobic residues
        hydrophobic_list=['GLY','ALA','LEU','ILE','VAL','PRO','PHE','MET','TRP']
        
        neighbor_list_protein = []
        # Get the number of nearby residues for all hydrophobic residues in the domain
        for residue_number in domain_residues:
            residue = structure[0]["A"][residue_number]
            if residue.get_resname() in hydrophobic_list and not residue.is_disordered(): 
                center_atoms = Selection.unfold_entities(residue, 'A')
                atom_list = [atom for atom in structure.get_atoms()]
                ns = NeighborSearch(atom_list)

                nearby_residues = {residue for center_atom in center_atoms
                    for residue in ns.search(center_atom.coord, distance_cutoff, 'R')}
                neighbor_list_protein.append(len(nearby_residues))

        return round(np.mean(neighbor_list_protein),3)


# Function to get the approximated Volume of the domain
def get_volume(structure,domain_residues,volume_resolution):

    reduced_coords = list()
    for residue_number in domain_residues:
        residue = structure[0]["A"][residue_number]
        for atm in residue.get_atoms():
            x_reduced = int(atm.coord.tolist()[0]/volume_resolution)
            y_reduced = int(atm.coord.tolist()[1]/volume_resolution)
            z_reduced = int(atm.coord.tolist()[2]/volume_resolution)
            reduced_coords.append([x_reduced,y_reduced,z_reduced])

    unique_coords = [list(x) for x in set(tuple(x) for x in reduced_coords)]
    volume = len(unique_coords) * volume_resolution**3
   
    return volume



# Function to get the normed radius of gyration
def calculate_normed_radius_of_gyration(
    af_domain_id: AFDomainID,
    af_chain_mmcif_dir: Path,
    volume_resolution: int,
    *,
    cif_filename=None,
) -> AFDomainID:

    chopping = af_domain_id.chopping
    segments = chopping.segments
    domain_residues = []
    for segment in segments:
        for res in range(segment.start,segment.end):
            domain_residues.append(res)


    # create default filename
    if cif_filename is None:
        cif_filename = af_domain_id.af_chain_id + ".cif.gz"

    cif_path = Path(af_chain_mmcif_dir, cif_filename)

    with gzip.open(f"{cif_path}", mode="rt") as cif_fh:

        structure = MMCIFParser(QUIET=1).get_structure("struc", cif_fh)

    coords = []
    masses = []
    # Fill lists for the coords and masses for all atoms in the domain
    for residue_number in domain_residues:
        residue = structure[0]["A"][residue_number]

        for atom in residue:
            coords.append(atom.coord.tolist())
            if atom.element == 'C':
                masses.append(12.0107)
            elif atom.element == 'H':    # We do not consider hydrogens for the radius of gyration
                masses.append(0.0)
            elif atom.element == 'O':
                masses.append(15.9994)
            elif atom.element == 'N':
                masses.append(14.0067)
            elif atom.element == 'S':
                masses.append(32.065)    
            else:
                raise Exception(f'Protein contains unknown atom type {atom.element}')

    # Calculate the radius of gyration
    mass_coords = [(m*i, m*j, m*k) for (i, j, k), m in zip(coords, masses)]
    total_mass = sum(masses)
    r_tmp = sum(mi*i + mj*j + mk*k for (i, j, k), (mi, mj, mk) in zip (coords, mass_coords))
    m_tmp = sum((sum(i) / total_mass)**2 for i in zip (*mass_coords))
    radius_gyration = np.sqrt(r_tmp / total_mass-m_tmp)


    # Calculate the radius of gyration for a sphere with the same volume as the protein
    volume = get_volume(structure,domain_residues,volume_resolution)
    radius=np.cbrt(volume)*3/4*np.pi
    radius_gyration_sphere = np.sqrt(3/5*radius**2)

    # Return normed radius of gyration
    return(round(radius_gyration/radius_gyration_sphere,3))


