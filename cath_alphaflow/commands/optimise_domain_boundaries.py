from Bio.PDB import MMCIFParser
import logging
from pathlib import Path

import click


from cath_alphaflow.io_utils import get_af_domain_id_reader
from cath_alphaflow.io_utils import get_csv_dictwriter
from cath_alphaflow.models import AFDomainID


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
    "--af_domain_list_post_tailchop",
    type=click.File("wt"),
    required=True,
    help="Output: CSV file for AF2 domain list after chopping",
)
@click.option(
    "--af_domain_mapping_post_tailchop",
    type=click.File("wt"),
    required=True,
    help="Output: CSV file for mapping of AF2 domain before/after chopping",
)
@click.option(
    "--cutoff_plddt_score",
    type=int,
    required=False,
    default=70,
    help="pLDDT cut-off score to apply to the domains. If not given, a cut-off score of 70 will be applied",
)
def optimise_domain_boundaries(
    af_domain_list,
    af_chain_mmcif_dir,
    af_domain_list_post_tailchop,
    af_domain_mapping_post_tailchop,
    cutoff_plddt_score,
):
    "Adjusts the domain boundaries of AF2 by removing unpacked tails"



    af_domain_list_reader = get_af_domain_id_reader(af_domain_list)

    af_domain_list_post_tailchop_writer = get_csv_dictwriter(
        af_domain_list_post_tailchop, fieldnames=["af_domain_id"]
    )
    af_domain_list_post_tailchop_writer.writeheader()

    af_mapping_list_post_tailchop_writer = get_csv_dictwriter(
        af_domain_mapping_post_tailchop,
        fieldnames=["af_domain_id_orig", "af_domain_id_post_tailchop"],
    )
    af_mapping_list_post_tailchop_writer.writeheader()

    click.echo(
        f"Chopping tails from AF domains"
        f"(mmcif_dir={af_chain_mmcif_dir}, in_file={af_domain_list.name}, "
        f"out_file={af_domain_list_post_tailchop.name} ) ..."
    )
    for af_domain_id in af_domain_list_reader:
        print(af_domain_id)

        af_domain_id_post_tailchop = calculate_domain_id_post_tailchop(
            af_domain_id, af_chain_mmcif_dir, cutoff_plddt_score
        )

        af_domain_list_post_tailchop_writer.writerow(
            {"af_domain_id": af_domain_id_post_tailchop}
        )

        af_mapping_list_post_tailchop_writer.writerow(
            {
                "af_domain_id_orig": af_domain_id,
                "af_domain_id_post_tailchop": af_domain_id_post_tailchop,
            }
        )

    click.echo("DONE")

def cut_boundary(structure,boundary_to_cut,cutoff_plddt_score,cut_start,cut_end):
    # Cuts one boundary based on the cutoff_plddt_score and returns the reduced boundary
    new_boundary=''
    boundary_lower_value = int(boundary_to_cut.split('-')[0])
    boundary_higher_value = int(boundary_to_cut.split('-')[1])
    new_boundary_lower_value = boundary_lower_value
    new_boundary_higher_value = boundary_higher_value

    if cut_start:
        for res in range(boundary_lower_value,boundary_higher_value):
            if structure[0]["A"][res]["CA"].get_bfactor() > cutoff_plddt_score:
                break
            else:
                new_boundary_lower_value += 1
                if new_boundary_lower_value == new_boundary_higher_value:
                    return('NaN')


    if cut_end:
        for res in range(boundary_higher_value,boundary_lower_value, -1):
            if structure[0]["A"][res]["CA"].get_bfactor() > cutoff_plddt_score:
                break
            else:
                new_boundary_higher_value -= 1
                if new_boundary_higher_value == new_boundary_lower_value:
                    return('NaN')
                

    if (new_boundary_higher_value - new_boundary_lower_value) + 0:
        new_boundary=str(new_boundary_lower_value)+"-"+str(new_boundary_higher_value)

    return new_boundary


def cut_first_boundary(structure,boundaries,cutoff_plddt_score):
    # Cut from the start of the first boundary until a residue has a pLDDT score > cutoff_plddt_score
    new_boundary = cut_boundary(structure,boundaries[0],cutoff_plddt_score,cut_start=True,cut_end=False)
    if new_boundary == 'NaN':
        boundaries.remove(boundaries[0])
        if len(boundaries) == 0:
            return []
        else:
            cut_first_boundary(structure,boundaries,cutoff_plddt_score)
    else:
        boundaries[0] = new_boundary
        return boundaries


def cut_last_boundary(structure,boundaries,cutoff_plddt_score):
    # Cut from the end of the last boundary until a residue has a pLDDT score > cutoff_plddt_score
    new_boundary = cut_boundary(structure,boundaries[-1],cutoff_plddt_score,cut_start=False,cut_end=True)
    if new_boundary == 'NaN':
        boundaries.remove(boundaries[-1])
        if len(boundaries) == 1:
            return boundaries
        else:
            cut_last_boundary(structure,boundaries,cutoff_plddt_score)
    else:
        boundaries[0] = new_boundary
        return boundaries



def calculate_domain_id_post_tailchop(
    af_domain_id: AFDomainID, af_chain_mmcif_dir: Path, cutoff_plddt_score: int
):
    af_domain_id_post_tailchop = None

    up_id = str(af_domain_id).split('/')[0]
    old_boundaries = str(af_domain_id).split('/')[1].split('_')
    structure = MMCIFParser(QUIET=1).get_structure('struc',af_chain_mmcif_dir+'/'+up_id+'.cif')
    if len(old_boundaries) == 1:
        # For single region domains just cut the one region
        new_boundary = cut_boundary(structure,old_boundaries[0],cutoff_plddt_score,cut_start=True,cut_end=True)
        af_domain_id_post_tailchop=up_id+'/'+new_boundary
        return af_domain_id_post_tailchop
    else:
        # For contigs, only cut the outer most parts and leave everything in the middle intact
        # We do this by first cutting from the start and then from the end. 
        cut_first_boundary(structure,old_boundaries,cutoff_plddt_score)
        if len(old_boundaries) == 0:
            return up_id+'/NaN'
        cut_last_boundary(structure,old_boundaries,cutoff_plddt_score)
        if len(old_boundaries) == 0:
            return up_id+'/NaN'

        af_domain_id_post_tailchop=up_id+'/'+"_".join(str(boundary) for boundary in old_boundaries)

        return af_domain_id_post_tailchop

