from unittest import mock
import logging

import cx_Oracle

from cath_alphaflow.db_utils import OraDB
from cath_alphaflow.settings import get_default_settings

LOG = logging.getLogger(__name__)


def test_mock_query(create_mock_query):
    expected_rows = [["P00520"]]
    create_mock_query(["uniprot_id"], expected_rows)
    db = OraDB()
    rows = list(db.yieldall("select * from foo"))
    assert rows == expected_rows


@mock.patch.object(cx_Oracle, "connect")
def test_mock_connection(mock_connect, mock_settings):

    settings = get_default_settings()

    OraDB()

    mock_connect.assert_called_with(
        user=settings.ORACLE_DB_USERNAME,
        password=settings.ORACLE_DB_PASSWORD,
        dsn=(
            f"(DESCRIPTION="
            f"(ADDRESS=(PROTOCOL=TCP)(HOST={settings.ORACLE_DB_HOST})(PORT={settings.ORACLE_DB_PORT}))"
            f"(CONNECT_DATA=(SID={settings.ORACLE_DB_SID}))"
            f")"
        ),
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
