import csv
from dataclasses import dataclass
import itertools
import logging

import click
import cx_Oracle


logging.basicConfig(
    level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s"
)
LOG = logging.getLogger()

DEFAULT_CHUNK_SIZE = 1000
SEQUENCE_TABLE_NAMES = ["sequences", "sequences_extra"]
UNIPROT_ACC_TABLE = ["uniprot_prim_acc"]


@dataclass
class CathDataset:
    uniprot_acc: str
    sequence_md5: str
    gene3d_domain_id: str
    bitscore: str
    resolved: str


def get_cathora_connection():
    """Connect to the Oracle Database"""
    dsn = cx_Oracle.makedsn("localhost", 1521, sid="cathora1")
    conn = cx_Oracle.connect(user="orengoreader", password="orengoreader", dsn=dsn)
    return conn


def get_csv_dictwriter(csvfile, fieldnames, **kwargs):
    """Common CSV writer"""
    return csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter="\t", **kwargs)


def get_csv_dictreader(csvfile, **kwargs):
    """Common CSV reader"""
    return csv.DictReader(csvfile, delimiter="\t", **kwargs)


def chunked_iterable(iterable, *, size):
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, size))
        if not chunk:
            break
        yield chunk


@click.command()
@click.option(
    "--uniprot_ids_csv",
    type=click.File("wt"),
    required=True,
    help="Output: CSV file containing UniProt IDs",
)
@click.option(
    "--dbname",
    "gene3d_dbname",
    type=str,
    help="Param: database to use when querying sequences",
)
@click.option(
    "--max_evalue",
    "max_evalue",
    type=str,
    help="Param: restrict UniProt IDs by independent evalue of Gene3D hits",
)
@click.option(
    "--max_records",
    "max_records",
    type=str,
    help="Param: maximum records to return",
)
def create_dataset_uniprot_ids(uniprot_ids_csv, gene3d_dbname, max_evalue, max_records):
    "Creates UniProt IDs for given dataset"

    conn = get_cathora_connection()
    headers = ["uniprot_acc"]

    csv_writer = get_csv_dictwriter(uniprot_ids_csv, fieldnames=headers)
    csv_writer.writeheader()

    click.echo("Writing ")
    for entry in next_cath_dataset_entry(
        conn,
        max_independent_evalue=max_evalue,
        max_records=max_records,
        gene3d_dbname=gene3d_dbname,
    ):
        csv_writer.writerow([entry.uniprot_acc])


@click.command()
@click.option(
    "--csv_uniprot_ids",
    type=click.File("rt"),
    required=True,
    help="Input: CSV file containing UniProt IDs",
)
@click.option(
    "--csv_uniprot_md5",
    type=click.File("wt"),
    help="Output: CSV file of UniProt to MD5 mapping",
)
@click.option(
    "--gene3d_crh_output",
    type=click.File("wt"),
    help="Output: CRH output file for Gene3D domains",
)
@click.option(
    "--af2_domainlist_ids",
    type=click.File("wt"),
    help="Output: CSV file of AF2 domain ids",
)
@click.option(
    "--af2_chainlist_ids",
    type=click.File("wt"),
    help="Output: CSV file of AF2 chain ids",
)
@click.option(
    "--af2_cath_annotations",
    type=click.File("wt"),
    help="Output: CSV file of CATH annotations",
)
@click.option(
    "--dbname",
    "gene3d_dbname",
    type=str,
    help="Param: database to use when querying sequences",
)
@click.option(
    "--max_evalue",
    "max_evalue",
    type=str,
    help="Param: restrict UniProt IDs by independent evalue of Gene3D hits",
)
@click.option(
    "--chunk",
    "chunk_size",
    type=int,
    default=DEFAULT_CHUNK_SIZE,
    help="Param: size of chunk when processing data",
)
def create_dataset_cath_files(
    csv_uniprot_ids,
    csv_uniprot_md5,
    crh_output,
    af2_domainlist_ids,
    af2_chainlist_ids,
    af2_cath_annotations,
    gene3d_dbname,
    chunk_size,
):
    "Creates CATH data files for a given dataset"

    conn = get_cathora_connection()

    # setup reader
    uniprot_reader = get_csv_dictreader(csv_uniprot_ids)
    # skip headers
    next(uniprot_reader)

    # setup writers
    crh_output_headers = [
        "sequence_md5",
        "gene3d_domain_id",
        "score",
        "boundaries",
        "resolved",
    ]
    crh_output_writer = get_csv_dictwriter(crh_output, fieldnames=crh_output_headers)
    # NOTE: do not write headers to crh output

    csv_uniprot_md5_headers = ["uniprot_acc", "sequence_md5"]
    csv_uniprot_md5_writer = get_csv_dictwriter(
        csv_uniprot_md5, fieldnames=csv_uniprot_md5_headers
    )
    csv_uniprot_md5_writer.writeheader()

    af2_domainlist_headers = ["af2_domain_id"]
    af2_domainlist_writer = get_csv_dictwriter(
        af2_domainlist_ids, fieldnames=af2_domainlist_headers
    )
    af2_domainlist_writer.writeheader()

    af2_chainlist_headers = ["af2_chain_id"]
    af2_chainlist_writer = get_csv_dictwriter(
        af2_chainlist_ids, fieldnames=af2_chainlist_headers
    )
    af2_chainlist_writer.writeheader()

    af2_cath_annotations_headers = ["AF2_CHAIN_ID"]
    af2_cath_annotations_writer = get_csv_dictwriter(
        af2_cath_annotations, fieldnames=af2_cath_annotations_headers
    )
    af2_cath_annotations_writer.writeheader()

    # how we are going to process a chunk of uniprot ids
    def process_uniprot_ids(uniprot_ids):

        for entry in next_cath_dataset_entry(
            conn,
            max_independent_evalue=1e-50,
            gene3d_dbname=gene3d_dbname,
            uniprot_ids=uniprot_ids,
        ):
            # sort out all variables that we are going to use in our data files
            uniprot_acc = entry.uniprot_acc
            sequence_md5 = entry.sequence_md5
            gene3d_domain_id = entry.gene3d_domain_id
            bitscore = entry.bitscore
            chopping = entry.chopping
            af2_domain_id = "???"
            af2_chain_id = "???"

            # write data
            csv_uniprot_md5_writer.writerow([uniprot_acc, sequence_md5])
            af2_domainlist_writer.writerow([af2_domain_id])
            af2_chainlist_writer.writerow([af2_chain_id])
            af2_cath_annotations_writer.writerow([])

            # use the 'correct' chopping in both cols
            crh_output_writer.writerow(
                [sequence_md5, gene3d_domain_id, bitscore, chopping, chopping]
            )

    # process chunks of uniprot ids
    for chunked_uniprot_rows in chunked_iterable(uniprot_reader, chunk_size=chunk_size):
        uniprot_ids = [row.get("UNIPROT_ACC") for row in chunked_uniprot_rows]
        process_uniprot_ids(uniprot_ids)


def make_dict_factory(cursor):
    """Turn a row into a dict"""
    col_names = [d[0] for d in cursor.description]

    def create_row(*args):
        return dict(zip(col_names, args))

    return create_row


def next_cath_dataset_entry(
    conn,
    *,
    max_independent_evalue,
    gene3d_dbname,
    max_records=None,
    uniprot_ids=None,
) -> CathDataset:
    """
    Returns a generator that provides CathDataset entries
    """

    if not max_records and not uniprot_ids:
        raise RuntimeError("need to specify one of [max_records, uniprot_ids]")

    sql_args = {}
    sql_where = ""
    if uniprot_ids:
        uniprot_id_placeholders = {
            "u" + num: uniprot_id for num, uniprot_id in enumerate(uniprot_ids, 1)
        }
        sql_args.update(uniprot_id_placeholders)
        sql_where += (
            " upa.ACCESSION IN ("
            + ", ".join([":u" + num for num in range(1, len(uniprot_ids))])
            + ") "
        )

    if max_independent_evalue:
        sql_args.update({"max_independent_evalue": max_independent_evalue})
        sql_where += " INDEPENDENT_EVALUE <= :max_independent_evalue "

    if max_records:
        sql_args.update({"max_records": max_records})
        sql_where += " ROWNUM <= :max_records "

    with conn.cursor() as curs:
        try:
            sql = f"""
SELECT
    upa.ACCESSION                           AS uniprot_acc,
    upa.SEQUENCE_MD5                        AS sequence_md5,
    DOMAIN_ID || '__' || SUPERFAMILY || '/'
        || REPLACE(RESOLVED, ',', '_')      AS gene3d_domain_id,
    SCORE                                   AS bitscore,
    RESOLVED                                AS chopping
FROM
    {gene3d_dbname}.CATH_DOMAIN_PREDICTIONS cdp
    LEFT JOIN {gene3d_dbname}.UNIPROT_PRIM_ACC upa
        ON (cdp.SEQUENCE_MD5 = upa.SEQUENCE_MD5)
WHERE
    {sql_where}
ORDER BY
    INDEPENDENT_EVALUE ASC
"""
            curs.execute(sql, **sql_args)
            curs.rowfactory = make_dict_factory(curs)
        except cx_Oracle.Error as err:
            LOG.warning(f"SQL command failed: {err}")
            raise
        for rowdict in curs:
            entry = CathDataset(**rowdict)
            yield entry


def dump_md5_crh_output(
    conn, gene3d_dbname, output_fh, independent_evalue, max_records
):
    md5_list = []
    with conn.cursor() as cur:
        try:
            sql = f"""
SELECT
    SEQUENCE_MD5,
    DOMAIN_ID || '__' || SUPERFAMILY || '/' || REPLACE(RESOLVED,',','_') AS GENE3D_DOMAIN_ID,
    SCORE,
    BOUNDARIES,
    RESOLVED
FROM {gene3d_dbname}.CATH_DOMAIN_PREDICTIONS cdp
WHERE
    INDEPENDENT_EVALUE <= :independent_evalue
    AND
    ROWNUM <= :max_records
ORDER BY
    INDEPENDENT_EVALUE ASC
;"""
            cur.execute(
                sql, independent_evalue=independent_evalue, max_records=max_records
            )
        except cx_Oracle.Error as err:
            LOG.warning(f"SQL command failed: {err}")
            raise
        for sequence_md5, gene3d_domain_id, score, boundaries, resolved in cur:
            output_fh.write(
                "\t".join([sequence_md5, gene3d_domain_id, score, boundaries, resolved])
                + "\n"
            )
            md5_list.append(sequence_md5)
