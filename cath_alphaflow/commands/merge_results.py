from pathlib import Path
import logging
import click

import pandas as pd
import numpy as np

LOG = logging.getLogger()

PK_COLNAME = "af_domain_id"


@click.command("merge-results")
@click.option(
    "--id_file",
    type=click.File("rt"),
    required=True,
    help="Input: CSV file containing list of ids to process (AF domain id in first col)",
)
@click.option(
    "--uniprot_md5_file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
    required=False,
    help="Input: UniProt to MD5 mapping",
)
@click.option(
    "--cath_annotation_file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=True),
    required=False,
    help="Input: CATH dataset file",
)
@click.option(
    "--plddt_file",
    type=click.Path(exists=False, file_okay=True, dir_okay=False, resolve_path=True),
    required=True,
    help="Input: plddt file",
)
@click.option(
    "--glob_file",
    type=click.Path(exists=False, file_okay=True, dir_okay=False, resolve_path=True),
    required=True,
    help="Input: globularity summary file",
)
@click.option(
    "--results_file",
    type=click.Path(exists=False, file_okay=True, dir_okay=False, resolve_path=True),
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
    cath_annotation_file,
    plddt_file,
    glob_file,
    results_file,
    af_version,
):
    "Merge results into a single CSV"

    merger = ResultsMerger(
        id_file=id_file,
        uniprot_md5_file=uniprot_md5_file,
        cath_annotation_file=cath_annotation_file,
        plddt_file=plddt_file,
        glob_file=glob_file,
    )
    df = merger.get_merged_df()
    df.to_csv(results_file)

    click.echo("DONE")


def csv_to_df(csv_file, *, cols=None, index_col=None, **kwargs):
    """
    Converts a CSV file to a pandas DataFrame with the expected column names
    """
    try:
        LOG.info(f"read_csv({csv_file}, {kwargs})")
        df_in = pd.read_csv(csv_file, **kwargs)
        LOG.info(f"AFTER read_csv: df_in: {df_in}")
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


class ResultsMerger:
    """
    Merges columns from the equivalent rows in various files
    """

    def __init__(
        self,
        *,
        id_file,
        plddt_file,
        glob_file,
        uniprot_md5_file=None,
        cath_annotation_file=None,
    ):
        self.id_file = id_file
        self.plddt_file = plddt_file
        self.glob_file = glob_file
        self.uniprot_md5_file = uniprot_md5_file
        self.cath_annotation_file = cath_annotation_file

    def get_merged_df(self):
        df = self.get_id_df()

        if self.uniprot_md5_file:
            uniprot_md5_df = self.get_uniprot_md5_df()
            df = df.join([uniprot_md5_df])

        if self.cath_annotation_file:
            cath_annotation_df = self.get_cath_annotation_df()
            df = df.join([cath_annotation_df])

        plddt_df = self.get_plddt_df().drop(columns=[PK_COLNAME])
        glob_df = self.get_glob_df().drop(columns=[PK_COLNAME])

        df = df.join([plddt_df])
        df = df.join([glob_df])

        return df

    def get_uniprot_md5_df(self):
        df = csv_to_df(
            self.uniprot_md5_file,
            sep="\s+",
            index_col="uniprot_acc",
        )
        return df

    def get_id_df(self):
        df = csv_to_df(
            self.id_file,
            sep="\s+",
            names=["af_domain_id"],
            cols=["af_domain_id"],
            index_col="af_domain_id",
        )
        return df

    def get_cath_annotation_df(self):
        df = csv_to_df(
            self.cath_annotation_file,
            cols=["cath_domain_id", "uniprot_acc", "md5", "bitscore", "chopping"],
            sep="\t",
        )
        df[PK_COLNAME] = df["uniprot_acc"] + "/" + df["chopping"].str.replace(",", "_")
        df.set_index(PK_COLNAME, inplace=True)
        df[PK_COLNAME] = df.index
        return df

    def get_plddt_df(self):
        df = csv_to_df(
            self.plddt_file,
            cols=["af_domain_id", "avg_plddt", "perc_LUR", "residues_total"],
            sep="\t",
        )
        df[PK_COLNAME] = df["af_domain_id"]
        df.set_index(PK_COLNAME, inplace=True)
        df[PK_COLNAME] = df.index
        return df

    def get_glob_df(self):
        df = csv_to_df(
            self.glob_file,
            cols=["af_domain_id", "packing_density", "normed_radius_gyration"],
            sep="\t",
        )
        df[PK_COLNAME] = df["af_domain_id"]
        df.set_index(PK_COLNAME, inplace=True)
        df[PK_COLNAME] = df.index
        return df
