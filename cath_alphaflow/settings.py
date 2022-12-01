from prettyconf import config

DEFAULT_AF_VERSION = 3
DEFAULT_AF_FRAGMENT = 1


class Settings:
    ORACLE_DB_HOST = config("ORACLE_DB_HOST", default=None)
    ORACLE_DB_PORT = config("ORACLE_DB_PORT", default=1521)
    ORACLE_DB_SID = config("ORACLE_DB_SID", default=None)
    ORACLE_DB_USERNAME = config("ORACLE_DB_USERNAME", default=None)
    ORACLE_DB_PASSWORD = config("ORACLE_DB_PASSWORD", default=None)
    DSSP_BINARY_PATH = config("DSSP_BINARY_PATH", default="mkdssp")
    DSSP_PDB_DICT = config("DSSP_PDB_DICT", default=None)

    def to_dict(self):
        dict = {}
        for key in dir(self):
            if key.startswith("__"):
                continue
            val = getattr(self, key)
            if callable(val):
                continue
            if "PASSWORD" in key:
                val = "******"
            dict[key] = val
        return dict


class ProductionSettings(Settings):
    pass


class TestSettings(Settings):
    ORACLE_DB_HOST = "TEST_ORACLE_DB_HOST"
    ORACLE_DB_PORT = "TEST_ORACLE_DB_PORT"
    ORACLE_DB_SID = "TEST_ORACLE_DB_SID"
    ORACLE_DB_USERNAME = "TEST_ORACLE_DB_USERNAME"
    ORACLE_DB_PASSWORD = "TEST_ORACLE_DB_PASSWORD"


def get_default_settings():
    return ProductionSettings()
