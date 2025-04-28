from pathlib import Path
import logging
import click

import pandas as pd
import numpy as np

LOG = logging.getLogger()

PK_COLNAME = "uniprot_domain_id"


@click.command("merge-domain-results")
@click.option(
    "--id_file",
    "-i",
    type=click.File("rt"),
    required=True,
    help="Input: CSV file containing list of ids to process (AF domain id in first col)",
)
@click.option(
    "--results_file",
    "-o",
    type=click.Path(exists=False, file_okay=True, dir_okay=False, resolve_path=True),
    required=True,
    help="Output: combined results file",
)
@click.option(
    "--chainsaw_file",
    "-c",
    type=click.Path(exists=False, file_okay=True, dir_okay=False, resolve_path=True),
    required=True,
    help="Output: chainsaw results file",
)
@click.option(
    "--merizo_file",
    "-m",
    type=click.Path(exists=False, file_okay=True, dir_okay=False, resolve_path=True),
    required=True,
    help="Output: merizo results file",
)
@click.option(
    "--unidoc_file",
    "-u",
    type=click.Path(exists=False, file_okay=True, dir_okay=False, resolve_path=True),
    required=True,
    help="Output: unidox results file",
)
def merge_domain_results_command(
    id_file,
    results_file,
    chainsaw_file,
    merizo_file,
    unidoc_file,
):
    "Merge results into a single CSV"

    merger = DomainResultsMerger(
        id_file=id_file,
        chainsaw_file=chainsaw_file,
        merizo_file=merizo_file,
        unidoc_file=unidoc_file,
    )
    df = merger.get_merged_df()
    df.to_csv(results_file)

    click.echo("DONE")


def csv_to_df(csv_file, *, cols=None, index_col=None, **kwargs):
    """
    Converts a CSV file to a pandas DataFrame with the expected column names
    """
    LOG.info(f"Reading CSV file {csv_file} ... ({kwargs})")
    try:
        df_in = pd.read_csv(csv_file, **kwargs)
    except Exception as err:
        msg = f"failed to read CSV file '{csv_file}' index_col={index_col}, kwargs={kwargs} (err={err})"
        LOG.error(msg)
        raise

    if index_col:
        df_in = df_in.set_index(index_col)
        df_in[index_col] = df_in.index

    if cols:
        LOG.info(f"DF: {df_in}")
        LOG.info(f"Restricting columns to {cols} ... (cols: {df_in.columns})")
        try:
            df_in = df_in[cols]
        except KeyError as err:
            msg = f"failed to get subset of columns={cols} from {[str(col) for col in df_in.columns]} in CSV file '{csv_file}' (err={err})"
            LOG.error(msg)
            raise KeyError(msg)

    return df_in


class DomainResultsMerger:
    """
    Merges columns from the equivalent rows in various files
    """

    def __init__(
        self,
        *,
        id_file,
        chainsaw_file,
        merizo_file,
        unidoc_file,
    ):
        self.id_file = id_file
        self.chainsaw_file = chainsaw_file
        self.merizo_file = merizo_file
        self.unidoc_file = unidoc_file

    def get_merged_df(self):
        id_df = self.get_id_df()
        dfs = {
            "chainsaw": self.get_chainsaw_df(),
            "merizo": self.get_merizo_df(),
            "unidoc": self.get_unidoc_df(),
        }
        df = id_df
        for name, df_in in dfs.items():
            try:
                df = df.join(df_in, rsuffix=f"_{name}")
            except:
                LOG.error(f"failed to join {name} results to df")
                raise
        return df

    def get_id_df(self):
        df = csv_to_df(
            self.id_file,
            sep="\s+",
            names=["model_id"],
            cols=["model_id"],
            index_col="model_id",
        )
        df[PK_COLNAME] = df["model_id"]
        df.set_index(PK_COLNAME, inplace=True)
        df[PK_COLNAME] = df.index
        return df

    def get_chainsaw_df(self):
        df = csv_to_df(
            self.chainsaw_file,
            names=[
                "chain_id",
                "sequence_md5",
                "nres",
                "ndom",
                "chopping",
                "uncertainty",
            ],
            sep="\t",
        )
        df.rename({"sequence_md5": "md5"}, axis=1, inplace=True)
        df[PK_COLNAME] = df["chain_id"] + "/" + df["chopping"]
        df.set_index(PK_COLNAME, inplace=True)
        df[PK_COLNAME] = df.index
        return df

    def get_merizo_df(self):
        df = csv_to_df(
            self.merizo_file,
            names=[
                "chain_id",
                "nres",
                "nres_dom",
                "nres_ndr",
                "ndom",
                "pIoU",
                "runtime",
                "chopping",
            ],
            sep="\t",
        )
        df[PK_COLNAME] = df["chain_id"] + "/" + df["chopping"]
        df.set_index(PK_COLNAME, inplace=True)
        df[PK_COLNAME] = df.index
        return df

    def get_unidoc_df(self):
        df = csv_to_df(
            self.unidoc_file,
            names=["chain_id", "md5", "ndom", "chopping", "score"],
            header=0,
            sep="\t",
        )
        # standarise the chopping format from unidoc
        df["chopping"] = df["chopping"].str.replace(",", "_")
        df[PK_COLNAME] = df["chain_id"] + "/" + df["chopping"]
        df.set_index(PK_COLNAME, inplace=True)
        df[PK_COLNAME] = df.index
        return df
