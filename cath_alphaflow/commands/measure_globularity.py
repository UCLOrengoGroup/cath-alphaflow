from Bio.PDB import MMCIFParser, Selection, NeighborSearch
import logging
from pathlib import Path
import gzip
import numpy as np
import click


from cath_alphaflow.io_utils import get_af_domain_id_reader
from cath_alphaflow.io_utils import get_csv_dictwriter
from cath_alphaflow.models.domains import AFDomainID
from cath_alphaflow.constants import DEFAULT_GLOB_DISTANCE, DEFAULT_GLOB_VOLUME

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
    "--gzipped_af_chains",
    type=bool,
    default=False,
    help="Set True if AF-chain files are stored in .gz files"
)
@click.option(
    "--af_domain_globularity",
    type=click.File("wt"),
    required=True,
    help="Output: CSV file for AF2 domain list including both globularity parameters",
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
    af_domain_list,
    af_chain_mmcif_dir,
    af_domain_globularity,
    distance_cutoff,
    volume_resolution,
    gzipped_af_chains,
):
    "Checks the globularity of the AF domain"

    af_domain_list_reader = get_af_domain_id_reader(af_domain_list)

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
        LOG.debug(f"Working on: {af_domain_id} ...")

        af_domain_packing_density = calculate_packing_density(
            af_domain_id, af_chain_mmcif_dir, distance_cutoff,gzipped_af_chains
        )

        af_domain_normed_radius_gyration = calculate_normed_radius_of_gyration(
            af_domain_id, af_chain_mmcif_dir,volume_resolution,gzipped_af_chains
        )
        LOG.debug(f"Processed entry: {af_domain_id}\tpacking_density:{af_domain_packing_density}\tgiration:{af_domain_normed_radius_gyration}")
        af_globularity_writer.writerow(
            {
                "af_domain_id": af_domain_id,
                "packing_density": af_domain_packing_density,
                "normed_radius_gyration": af_domain_normed_radius_gyration,
            }
        )

    click.echo("DONE")


def calculate_packing_density(
    af_domain_id: AFDomainID,
    af_chain_mmcif_dir: Path,
    distance_cutoff: int,
    gzipped_af_chains,
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
        if gzipped_af_chains == False:
            open_func = open
            cif_filename = af_domain_id.af_chain_id + ".cif"
        else:
            open_func = gzip.open
            cif_filename = af_domain_id.af_chain_id + ".cif.gz"

    cif_path = Path(af_chain_mmcif_dir, cif_filename)

    with open_func(f"{cif_path}", mode="rt") as cif_fh:

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
    gzipped_af_chains,
    *,
    cif_filename=None,
) -> AFDomainID:

    chopping = af_domain_id.chopping
    segments = chopping.segments
    domain_residues = []
    for segment in segments:
        for res in range(segment.start,segment.end):
            domain_residues.append(res)


    if cif_filename is None:
        if gzipped_af_chains == False:
            open_func = open
            cif_filename = af_domain_id.af_chain_id + ".cif"
        else:
            open_func = gzip.open
            cif_filename = af_domain_id.af_chain_id + ".cif.gz"

    cif_path = Path(af_chain_mmcif_dir, cif_filename)

    with open_func(f"{cif_path}", mode="rt") as cif_fh:

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


