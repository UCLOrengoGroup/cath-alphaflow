import cx_Oracle
import sys
import os

DEFAULT_INSTANTCLIENT_VERSION = "instantclient_19_8"

try:
    if sys.platform.startswith("darwin"):
        lib_dir_guess = os.path.join(
            os.environ.get("HOME"), "Downloads", DEFAULT_INSTANTCLIENT_VERSION
        )
        lib_dir = os.environ.get("LD_LIBRARY_PATH", lib_dir_guess)
        cx_Oracle.init_oracle_client(lib_dir=lib_dir)
    elif sys.platform.startswith("win32"):
        lib_dir = f"C:\\oracle\\{DEFAULT_INSTANTCLIENT_VERSION}"
        cx_Oracle.init_oracle_client(lib_dir=lib_dir)
except Exception as err:
    print("Failed to find instantclient: bugging out as cx_Oracle will fail ...")
    print(err)
    sys.exit(1)
