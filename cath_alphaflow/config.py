from prettyconf import config


class Config:
    ORACLE_DB_HOST = config("ORACLE_DB_HOST")
    ORACLE_DB_PORT = config("ORACLE_DB_PORT", 1521)
    ORACLE_DB_SID = config("ORACLE_DB_SID")
    ORACLE_DB_USERNAME = config("ORACLE_DB_USERNAME")
    ORACLE_DB_PASSWORD = config("ORACLE_DB_PASSWORD")


def get_default_config():
    return Config()
