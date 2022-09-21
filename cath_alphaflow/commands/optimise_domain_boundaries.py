import logging
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
def optimise_domain_boundaries(
    af_domain_list,
    af_chain_mmcif_dir,
    af_domain_list_post_tailchop,
    af_domain_mapping_post_tailchop,
):
    "Adjusts the domain boundaries of AF2 by removing unpacked tails"

    af_domain_list_reader = get_af_domain_id_reader(af_domain_list)

    af_domain_list_post_tailchop_writer = get_csv_dictwriter(
        af_domain_list_post_tailchop, fieldnames=["af_domain_id"]
    )
    af_domain_list_post_tailchop_writer.writeheader()

    af_domain_list_post_tailchop_writer = get_csv_dictwriter(
        af_domain_mapping_post_tailchop,
        fieldnames=["af_domain_id_orig", "af_domain_id_post_tailchop"],
    )
    af_domain_list_post_tailchop_writer.writeheader()

    click.echo(
        f"Chopping tails from AF domains"
        f"(mmcif_dir={af_chain_mmcif_dir}, in_file={af_domain_list.name}, "
        f"out_file={af_domain_list_post_tailchop_writer.name} ) ..."
    )
    for af_domain_id in af_domain_list_reader:

        af_domain_id_post_tailchop = calculate_domain_id_post_tailchop(
            af_domain_id, af_chain_mmcif_dir
        )

        af_domain_list_post_tailchop_writer.writerow(
            {"af_domain_id": af_domain_id_post_tailchop}
        )

        af_domain_mapping_post_tailchop.writerow(
            {
                "af_domain_id": af_domain_id_post_tailchop,
                "af_domain_id_post_tailchop": af_domain_id_post_tailchop,
            }
        )

    click.echo("DONE")


def calculate_domain_id_post_tailchop(af_domain_id: AFDomainID, af_chain_mmcif_dir):

    af_domain_id_post_tailchop = None

    # parse AF2 chopping from original domain id

    # process residue coordinates from mmCIF file

    # calculate updated chopping

    # return updated AF2 domain id

    return af_domain_id_post_tailchop
