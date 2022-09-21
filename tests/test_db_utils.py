from unittest import mock
import logging

import cx_Oracle

from cath_alphaflow.db_utils import OraDB

LOG = logging.getLogger(__name__)


@mock.patch.object(cx_Oracle, "connect")
def test_mock_connection(mock_connect):

    OraDB()

    mock_connect.assert_called_with(
        user="orengoreader",
        password="orengoreader",
        dsn="(DESCRIPTION=(ADDRESS=(PROTOCOL=TCP)(HOST=odb.cs.ucl.ac.uk)(PORT=1521))(CONNECT_DATA=(SID=cathora1)))",
    )


def test_yieldall(create_mock_query):

    mock_rows = [
        {"cath_code": "1", "description": "Mainly Alpha"},
        {"cath_code": "2", "description": "Mainly Beta"},
        {"cath_code": "3", "description": "Alpha Beta"},
        {"cath_code": "4", "description": "Few Secondary Structures"},
        {"cath_code": "6", "description": "Special"},
    ]

    create_mock_query(description=["cath_code", "description"], rows=mock_rows)

    db = OraDB()
    sql = "select cath_code, description from gene3d_21.cath_names where rownum <= :max_records"
    sql_args = {"max_records": 5}

    entries = list(db.yieldall(sql, sql_args, return_type=dict))

    assert entries == mock_rows
