from pathlib import Path
import logging
import click

import pandas as pd
import numpy as np

LOG = logging.getLogger()

PK_COLNAME = "uniprot_domain_id"


@click.command("combine-results")
@click.option(
    "--id_file",
    type=click.File("rt"),
    required=True,
    help="Input: CSV file containing list of ids to process (AF domain id in first col)",
)
@click.option(
    "--uniprot_md5_file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
    required=True,
    help="Input: UniProt to MD5 mapping",
)
@click.option(
    "--cath_annotation_file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
    required=True,
    help="Input: CATH dataset file",
)
@click.option(
    "--results_file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
    required=True,
    help="Output: combined results file",
)
@click.option(
    "--af_version",
    type=int,
    help=f"Option: specify the AF version when parsing uniprot ids",
)
def merge_results_command(
    id_file,
    uniprot_md5_file,
    dataset_file,
    results_file,
    af_version,
):
    "Merge results into a single spreadsheet"

    merger = ResultsMerger(
        id_file=id_file, uniprot_md5_file=uniprot_md5_file, dataset_file=dataset_file
    )
    df = merger.get_merged_df()
    df.to_csv(results_file)

    click.echo("DONE")


def csv_to_df(csv_file, *, cols=None, index_col=PK_COLNAME, **kwargs):
    """
    Converts a CSV file to a pandas DataFrame with the expected column names
    """
    df_in = pd.read_csv(csv_file, index_col=index_col, **kwargs)
    if cols:
        df_out = df_in.loc[cols]
    return df_out


class ResultsMerger:
    """
    Merges columns from the equivalent rows in various files
    """

    def __init__(self, *, id_file, uniprot_md5_file, dataset_file):
        self.id_file = id_file
        self.uniprot_md5_file = uniprot_md5_file
        self.dataset_file = dataset_file

    def get_merged_df(self):
        id_df = self.get_id_df()
        uniprot_md5_df = self.get_uniprot_md5_df()
        dataset_df = self.get_dataset_df()
        df = id_df.join([dataset_df])
        return df

    def get_uniprot_md5_df(self):
        df = csv_to_df(
            self.uniprot_md5_file,
            sep=None,
            cols=["uniprot_acc", "md5"],
            index_col="uniprot_acc",
        )
        return df

    def get_id_df(self):
        df = csv_to_df(
            self.id_file,
            sep=None,
            names=["uniprot_acc"],
            cols=["uniprot_acc"],
            index_col="uniprot_acc",
        )
        return df

    def get_dataset_df(self):
        df = csv_to_df(self.dataset_file, sep="\t")
        df_clean = df.loc[
            ["cath_domain_id", "uniprot_acc", "md5", "bitscore", "chopping"]
        ]
        df_clean[PK_COLNAME] = (
            df["uniprot_acc"] + "/" + df["chopping"].str.replace(",", "_")
        )
        df_clean.set_index(PK_COLNAME, inplace=True)
        return df_clean
