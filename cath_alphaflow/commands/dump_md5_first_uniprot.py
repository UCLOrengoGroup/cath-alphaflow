import click
import logging
import cx_Oracle


logging.basicConfig(
    level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s"
)
LOG = logging.getLogger()

DEFAULT_DBNAME = "gene3d_21"
DEFAULT_CHUNK_SIZE = 1000
SEQUENCE_TABLE_NAMES = ["sequences", "sequences_extra"]
UNIPROT_ACC_TABLE = ["uniprot_prim_acc"]


@click.command()
@click.option(
    "--dbname",
    "gene3d_dbname",
    type=str,
    default=DEFAULT_DBNAME,
    help="database to use when querying sequences",
)
@click.option(
    "--chunk",
    "chunk_size",
    type=int,
    default=DEFAULT_CHUNK_SIZE,
    help="chunk size for db queries",
)
@click.option(
    "--md5_list",
    "md5_list_file",
    type=click.File("rt"),
    required=True,
    help="input file containing md5s to retrieve from DB",
)
@click.option(
    "--md5_uniprot_file",
    "md5_uniprot_file",
    type=click.File("wt"),
    required=True,
    help="output file containing md5 to uniprot mapping",
)
def run(gene3d_dbname, chunk_size, md5_list_file, md5_uniprot_file):
    md5_list = []
    try:
        dsn = cx_Oracle.makedsn("localhost", 1521, sid="cathora1")
        conn = cx_Oracle.connect(user="orengoreader", password="orengoreader", dsn=dsn)

    except cx_Oracle.DatabaseError as er:
        print("Error Oracle:", er)
        raise
    # Read md5_list_fh into a list of md5s

    for line in md5_list_file:
        line = line.rstrip()
        md5_list.append(line)

    # Go through md5 list in chunks of 1000 and dump to file
    chunk_from = 0
    chunk_num = 1
    while True:
        chunk_to = chunk_size + chunk_from
        chunked_md5_list = md5_list[chunk_from:chunk_to]

        if len(chunked_md5_list) == 0:
            break
        LOG.info(f"Processing chunk {chunk_num}: {chunk_from}->{chunk_to}")
        dump_md5_first_uniprot_mapping(
            conn, chunked_md5_list, gene3d_dbname, md5_uniprot_file
        )
        chunk_from = chunk_to
        chunk_num += 1


def dump_md5_first_uniprot_mapping(conn, md5_list, gene3d_dbname, output_fh):
    output_fh.write(f"UNIPROT_ID\tDOMAIN_MD5\n")
    with conn.cursor() as cur:
        md5_placeholder = " OR ".join(
            [f"SEQUENCE_MD5=:{md5_num}" for md5_num in range(1, len(md5_list) + 1)]
        )
        try:
            sql = f"""
SELECT MIN(ACCESSION),
    SEQUENCE_MD5  
FROM {gene3d_dbname}.UNIPROT_PRIM_ACC upa 
WHERE {md5_placeholder}
GROUP BY SEQUENCE_MD5
"""
            cur.execute(
                sql,
                md5_list,
            )
        except cx_Oracle.Error as err:
            LOG.warning(f"SQL command failed: {err}")
            raise
        for accession, seq_md5 in cur:
            output_fh.write("\t".join([accession, seq_md5]) + "\n")


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


if __name__ == "__main__":
    run()
