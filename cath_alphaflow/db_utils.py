import logging
import cx_Oracle

from cath_alphaflow.models import PredictedCathDomain

LOG = logging.getLogger(__name__)


def get_cathora_connection(
    host="localhost",
    port=1521,
    sid="cathora1",
    user="orengoreader",
    password="orengoreader",
):
    """
    Connect to the Oracle Database
    """
    dsn = cx_Oracle.makedsn(host, port, sid=sid)
    conn = cx_Oracle.connect(user=user, password=password, dsn=dsn)
    return conn


def yieldall(connection, sql_stmt, *args, return_type=None):
    """
    Run SQL query and yield row one-by-one
    """

    with connection as conn:
        with conn.cursor() as curs:

            LOG.debug("SQL: %s (args:%s)", sql_stmt, *args)
            curs.execute(sql_stmt, *args)

            if return_type:
                if issubclass(return_type, dict):
                    curs.rowfactory = make_dict_factory(curs)
                else:
                    msg = f"Do not know how to process return_type={return_type}"
                    raise NotImplementedError(msg)

            for record in curs:
                yield record


def make_dict_factory(cursor):
    """Turn a row into a dict"""
    col_names = [d[0].lower() for d in cursor.description]

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
) -> PredictedCathDomain:
    """
    Returns a generator that provides `PredictedCathDomain` entries
    """

    if not max_records and not uniprot_ids:
        raise RuntimeError("need to specify one of [max_records, uniprot_ids]")

    sql_args = {}
    sql_where_args = []
    if uniprot_ids:
        uniprot_id_placeholders = {
            "u" + num: uniprot_id for num, uniprot_id in enumerate(uniprot_ids, 1)
        }
        sql_args.update(uniprot_id_placeholders)
        sql_where_args += [
            (
                "upa.ACCESSION IN ("
                + ", ".join([":u" + num for num in range(1, len(uniprot_ids))])
                + ")"
            )
        ]

    if max_independent_evalue:
        sql_args.update({"max_independent_evalue": max_independent_evalue})
        sql_where_args += ["INDEPENDENT_EVALUE <= :max_independent_evalue"]

    if max_records:
        sql_args.update({"max_records": max_records})
        sql_where_args += ["ROWNUM <= :max_records"]

    sql_where = " AND ".join(sql_where_args)

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
    INNER JOIN {gene3d_dbname}.UNIPROT_PRIM_ACC upa
        ON (cdp.SEQUENCE_MD5 = upa.SEQUENCE_MD5)
WHERE
    {sql_where}
    AND ROWNUM <= 100
ORDER BY
    INDEPENDENT_EVALUE ASC, upa.ACCESSION ASC
"""

    for rowdict in yieldall(conn, sql, sql_args, return_type=dict):
        entry = PredictedCathDomain(**rowdict)
        yield entry
