
"""
Check your connection to the oracle database here by running this file.
Instructions are in README.oracle.md
"""

use_new_lib = True

ORACLE_DB_HOST="localhost"
ORACLE_DB_PORT=1521
ORACLE_DB_SID="ask someone at UCL"
ORACLE_DB_USERNAME="ask someone at UCL"
ORACLE_DB_PASSWORD="ask someone at UCL"

dsn=(
            f"(DESCRIPTION="
            f"(ADDRESS=(PROTOCOL=TCP)(HOST={ORACLE_DB_HOST})(PORT={ORACLE_DB_PORT}))"
            f"(CONNECT_DATA=(SID={ORACLE_DB_SID}))"
            f")"
        )

if use_new_lib:
  import oracledb as db
  db.init_oracle_client() # this is for a thick client and needs local installation - for thin client comment out  
else:
  import cx_Oracle as db
  
conn = db.connect(user=ORACLE_DB_USERNAME, password=ORACLE_DB_PASSWORD, dsn=dsn)
print("Connected with oracle to",ORACLE_DB_USERNAME, ORACLE_DB_PASSWORD, ORACLE_DB_HOST, ORACLE_DB_PORT, ORACLE_DB_SID)
conn.close()

