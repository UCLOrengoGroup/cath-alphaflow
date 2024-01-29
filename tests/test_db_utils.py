from unittest import mock
import logging

import oracledb

from cath_alphaflow.db_utils import OraDB
from cath_alphaflow import settings

LOG = logging.getLogger(__name__)


def test_mock_query(create_mock_query, mock_settings):
    expected_rows = [["P00520"]]
    create_mock_query(["uniprot_id"], expected_rows)
    db = OraDB()
    rows = list(db.yieldall("select * from foo"))
    assert rows == expected_rows


@mock.patch.object(oracledb, "connect")
def test_mock_connection(mock_connect, mock_settings):
    OraDB()

    config = settings.TestSettings()

    mock_connect.assert_called_with(
        user=config.ORACLE_DB_USERNAME,
        password=config.ORACLE_DB_PASSWORD,
        dsn=(
            f"(DESCRIPTION="
            f"(ADDRESS=(PROTOCOL=TCP)(HOST={config.ORACLE_DB_HOST})(PORT={config.ORACLE_DB_PORT}))"
            f"(CONNECT_DATA=(SID={config.ORACLE_DB_SID}))"
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
