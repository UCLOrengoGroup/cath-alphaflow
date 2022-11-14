from pathlib import Path
from prettyconf import config

DEFAULT_AF_VERSION = 3
DEFAULT_AF_FRAGMENT = 1

PROJECT_ROOT_DIR = Path(__file__).parent.parent
DEFAULT_FS_BINARY_PATH = str(PROJECT_ROOT_DIR / "foldseek" / "bin" / "foldseek")


def resolve_path(raw_path_str):
    return str(Path(raw_path_str).resolve())


class Settings:
    ORACLE_DB_HOST = config("ORACLE_DB_HOST", default=None)
    ORACLE_DB_PORT = config("ORACLE_DB_PORT", default=1521)
    ORACLE_DB_SID = config("ORACLE_DB_SID", default=None)
    ORACLE_DB_USERNAME = config("ORACLE_DB_USERNAME", default=None)
    ORACLE_DB_PASSWORD = config("ORACLE_DB_PASSWORD", default=None)
    DSSP_BINARY_PATH = config("DSSP_BINARY_PATH", default=None)
    DSSP_PDB_DICT = config("DSSP_PDB_DICT", default=None)
    FS_BINARY_PATH = config(
        "FS_BINARY_PATH",
        default=DEFAULT_FS_BINARY_PATH,
        cast=resolve_path,
    )
    FS_DB_PATH = config("FS_DB_PATH", default="foldseek_db", cast=resolve_path)
    FS_TMP_PATH = config("FS_TMP_PATH", default="foldseek_tmp", cast=resolve_path)

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
