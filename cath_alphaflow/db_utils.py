import logging
import cx_Oracle

from cath_alphaflow.settings import get_default_settings
from cath_alphaflow.predicted_domain_provider import OraclePredictedCathDomainProvider

LOG = logging.getLogger(__name__)

config = get_default_settings()

DEFAULT_HOST = config.ORACLE_DB_HOST
DEFAULT_PORT = config.ORACLE_DB_PORT
DEFAULT_SID = config.ORACLE_DB_SID
DEFAULT_USER = config.ORACLE_DB_USERNAME
DEFAULT_PASSWORD = config.ORACLE_DB_PASSWORD


class OraDB(OraclePredictedCathDomainProvider, cx_Oracle.Connection):
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
