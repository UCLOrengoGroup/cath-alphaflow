import logging
import click

from cath_alphaflow.io_utils import get_csv_dictreader
from cath_alphaflow.io_utils import get_csv_dictwriter


LOG = logging.getLogger()

# Input:
# - `FILE:AF2_CHAIN_MMCIF`
# - `FILE:AF2_DOMAIN_LIST`

# Output:
# - `FILE:AF2_DOMAIN_LIST_POST_TAILCHOP` -- `(AF_domain_ID_tailchop)`
# - `FILE:AF2_DOMAIN_MAPPING_POST_TAILCHOP` -- `(AF_domain_ID_orig, AF_domain_ID_tailchop)`


@click.command()
@click.option(
    "--af2_domain_list",
    type=click.File("rt"),
    required=True,
    help="Input: CSV file containing AF2 domains",
)
@click.option(
    "--af2_chain_mmcif_dir",
    type=click.Path,
    required=True,
    help="Input: directory of mmCIF files",
)
@click.option(
    "--af2_domain_list_post_tailchop",
    type=click.File("wt"),
    required=True,
    help="Output: CSV file for AF2 domain list after chopping",
)
@click.option(
    "--af2_domain_mapping_post_tailchop",
    type=click.File("wt"),
    required=True,
    help="Output: CSV file for mapping of AF2 domain before/after chopping",
)
def create_dataset_uniprot_ids(
    af2_domain_list,
    af2_chain_mmcif_dir,
    af2_domain_list_post_tailchop,
    af2_domain_mapping_post_tailchop,
):
    "Adjusts the domain boundaries of AF2 by removing unpacked tails"

    af2_domain_list_reader = get_csv_dictreader(
        af2_domain_list,
    )
    next(af2_domain_list_reader)

    af2_domain_list_post_tailchop_writer = get_csv_dictwriter(
        af2_domain_list_post_tailchop, fieldnames=["af2_domain_id"]
    )
    af2_domain_list_post_tailchop_writer.writeheader()

    af2_domain_list_post_tailchop_writer = get_csv_dictwriter(
        af2_domain_mapping_post_tailchop,
        fieldnames=["af2_domain_id_orig", "af2_domain_id_post_tailchop"],
    )
    af2_domain_list_post_tailchop_writer.writeheader()

    click.echo(
        f"Chopping tails from AF2 domains"
        f"(mmcif_dir={af2_chain_mmcif_dir}, in_file={af2_domain_list.name}, "
        f"out_file={af2_domain_list_post_tailchop_writer.name} ) ..."
    )
    for rowdict in af2_domain_list_reader:

        af2_domain_id = rowdict["af2_domain_id"]
        af2_domain_id_post_tailchop = calculate_domain_id_post_tailchop(
            af2_domain_id, af2_chain_mmcif_dir
        )

        af2_domain_list_post_tailchop_writer.writerow(
            {"af2_domain_id": af2_domain_id_post_tailchop}
        )

        af2_domain_mapping_post_tailchop.writerow(
            {
                "af2_domain_id": af2_domain_id_post_tailchop,
                "af2_domain_id_post_tailchop": af2_domain_id_post_tailchop,
            }
        )

    click.echo("DONE")


def calculate_domain_id_post_tailchop(af2_domain_id, af2_chain_mmcif_dir):

    af2_domain_id_post_tailchop = None

    # parse AF2 chopping from original domain id

    # process residue coordinates from mmCIF file

    # calculate updated chopping

    # return updated AF2 domain id

    return af2_domain_id_post_tailchop
