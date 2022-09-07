from cath_alphaflow.db_utils import get_cathora_connection
from cath_alphaflow.db_utils import yieldall


def test_yieldall():

    conn = get_cathora_connection()
    sql = "select * from gene3d_21.cath_names where rownum <= :max_records"
    sql_args = {"max_records": 5}

    entries = yieldall(conn, sql, sql_args, return_type=dict)

    assert len(entries) == 5

    assert entries == [
        {"cath_code": "1", "description": "Mainly Alpha"},
        {"cath_code": "2", "description": "Mainly Beta"},
        {"cath_code": "3", "description": "Alpha Beta"},
        {"cath_code": "4", "description": "Few Secondary Structures"},
        {"cath_code": "6", "description": "Special"},
    ]
