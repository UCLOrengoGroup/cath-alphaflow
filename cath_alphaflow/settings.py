from prettyconf import config


class Settings:
    ORACLE_DB_HOST = None
    ORACLE_DB_PORT = None
    ORACLE_DB_SID = None
    ORACLE_DB_USERNAME = None
    ORACLE_DB_PASSWORD = None


class ProductionSettings(Settings):
    ORACLE_DB_HOST = config("ORACLE_DB_HOST")
    ORACLE_DB_PORT = config("ORACLE_DB_PORT", default=1521)
    ORACLE_DB_SID = config("ORACLE_DB_SID")
    ORACLE_DB_USERNAME = config("ORACLE_DB_USERNAME")
    ORACLE_DB_PASSWORD = config("ORACLE_DB_PASSWORD")


class TestSettings(Settings):
    ORACLE_DB_HOST = "TEST_ORACLE_DB_HOST"
    ORACLE_DB_PORT = "TEST_ORACLE_DB_PORT"
    ORACLE_DB_SID = "TEST_ORACLE_DB_SID"
    ORACLE_DB_USERNAME = "TEST_ORACLE_DB_USERNAME"
    ORACLE_DB_PASSWORD = "TEST_ORACLE_DB_PASSWORD"


def get_default_settings():
    return ProductionSettings()
