import logging
import oracledb as Oracle

from cath_alphaflow import settings
from cath_alphaflow.predicted_domain_provider import OraclePredictedCathDomainProvider

LOG = logging.getLogger(__name__)


class OraDB(OraclePredictedCathDomainProvider):
    def __init__(
        self,
        host=None,
        port=None,
        sid=None,
        user=None,
        password=None,
    ):
        config = settings.get_default_settings()

        if not host:
            host = config.ORACLE_DB_HOST
        if not port:
            port = config.ORACLE_DB_PORT
        if not sid:
            sid = config.ORACLE_DB_SID
        if not user:
            user = config.ORACLE_DB_USERNAME
        if not password:
            password = config.ORACLE_DB_PASSWORD

        LOG.info("Connecting to Oracle DB: %s@%s:%s/%s", user, host, port, sid)
        self._dsn = Oracle.makedsn(host, port, sid=sid)
        self._conn = Oracle.connect(user=user, password=password, dsn=self._dsn)

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
