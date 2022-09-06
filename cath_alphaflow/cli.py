import click

from .create_dataset import commands as create_dataset


@click.group()
@click.version_option()
def cli():
    "Workflow tools to assign CATH structural domains to AlphaFold predictions"
    pass


cli.add_command(create_dataset.create_dataset_uniprot_ids)

cli.add_command(create_dataset.create_dataset_cath_files)
