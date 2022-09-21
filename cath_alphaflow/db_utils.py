import logging
import cx_Oracle

from cath_alphaflow.settings import get_default_settings
from cath_alphaflow.models import PredictedCathDomain

LOG = logging.getLogger(__name__)

config = get_default_settings()

DEFAULT_HOST = config.ORACLE_DB_HOST
DEFAULT_PORT = config.ORACLE_DB_PORT
DEFAULT_SID = config.ORACLE_DB_SID
DEFAULT_USER = config.ORACLE_DB_USERNAME
DEFAULT_PASSWORD = config.ORACLE_DB_PASSWORD


class OraDB(cx_Oracle.Connection):
    def __init__(
        self,
        host=DEFAULT_HOST,
        port=DEFAULT_PORT,
        sid=DEFAULT_SID,
        user=DEFAULT_USER,
        password=DEFAULT_PASSWORD,
    ):
        self._dsn = cx_Oracle.makedsn(host, port, sid=sid)
        self._conn = cx_Oracle.connect(user=user, password=password, dsn=self._dsn)

    @property
    def conn(self):
        return self._conn

    def yieldall(self, sql_stmt, *args, return_type=None):
        """
        Run SQL query and yield rows one-by-one
        """

        with self.conn as conn:
            with conn.cursor() as curs:
                LOG.debug("curs: %s", curs)

                LOG.debug("SQL: %s (args:%s)", sql_stmt, *args)
                curs.execute(sql_stmt, *args)

                if return_type:
                    if issubclass(return_type, dict):
                        curs.rowfactory = self.make_dict_factory(curs)
                    else:
                        msg = f"Do not know how to process return_type={return_type}"
                        raise NotImplementedError(msg)

                for record in curs:
                    yield record

    def make_dict_factory(self, cursor):
        """Turn a row into a dict"""
        col_names = [d[0].lower() for d in cursor.description]

        def create_row(*args):
            return dict(zip(col_names, args))

        return create_row

    def next_cath_dataset_entry(
        self,
        *,
        gene3d_dbname,
        max_independent_evalue=None,
        max_records=None,
        uniprot_ids=None,
    ) -> PredictedCathDomain:
        """
        Returns a generator that provides `PredictedCathDomain` entries
        """

        if not max_records and not uniprot_ids:
            raise RuntimeError("need to specify one of [max_records, uniprot_ids]")

        # organise filters in the where clause
        # (SQL placeholders and cooresponding substitutions key/values)
        sql_args = {}
        sql_where_args = []
        if uniprot_ids:

            uniprot_id_placeholders = {
                f"u{num}": uniprot_id for num, uniprot_id in enumerate(uniprot_ids, 1)
            }
            sql_args.update(uniprot_id_placeholders)
            sql_where_args += [
                (
                    "upa.ACCESSION IN ("
                    + ", ".join([f":u{num}" for num in range(1, len(uniprot_ids) + 1)])
                    + ")"
                )
            ]

        if max_independent_evalue:
            sql_args.update({"max_independent_evalue": max_independent_evalue})
            sql_where_args += ["INDEPENDENT_EVALUE <= :max_independent_evalue"]

        if max_records:
            sql_args.update({"max_records": max_records})
            sql_where_args += ["ROWNUM <= :max_records"]

        # join all where clauses together with 'AND' operator
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
    ORDER BY
        INDEPENDENT_EVALUE ASC, upa.ACCESSION ASC
    """

        # execute query and yield results (as dict) row by row
        for rowdict in self.yieldall(sql, sql_args, return_type=dict):
            # convert dict to data structure
            entry = PredictedCathDomain(**rowdict)
            yield entry
