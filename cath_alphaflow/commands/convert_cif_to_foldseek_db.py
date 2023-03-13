import logging
from pathlib import Path
import os
import click
import subprocess
from cath_alphaflow.io_utils import yield_first_col
from cath_alphaflow.models.domains import AFDomainID
from cath_alphaflow.constants import (
    DEFAULT_FS_QUERYDB_SUFFIX,
    DEFAULT_FS_QUERYDB_NAME,
    DEFAULT_CIF_SUFFIX,
    ID_TYPE_AF_DOMAIN,
    ID_TYPE_UNIPROT_DOMAIN,
)
from cath_alphaflow.settings import get_default_settings
from cath_alphaflow.errors import ArgumentError

config = get_default_settings()

FS_BINARY_PATH = config.FS_BINARY_PATH

LOG = logging.getLogger()


@click.command()
@click.option(
    "--cif_dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, resolve_path=True),
    required=True,
    help="Input: directory of CIF files",
)
@click.option(
    "--cif_suffix",
    type=str,
    default=DEFAULT_CIF_SUFFIX,
    help=f"Input: optional suffix to add to id when looking for cif file (default: {DEFAULT_CIF_SUFFIX})",
)
@click.option(
    "--fs_querydb_name",
    type=str,
    default=DEFAULT_FS_QUERYDB_NAME,
    help=f"Output: Foldseek Query Database Name (default: {DEFAULT_FS_QUERYDB_NAME})",
)
@click.option(
    "--fs_querydb_suffix",
    type=str,
    default=DEFAULT_FS_QUERYDB_SUFFIX,
    help=f"Option: optional suffix to add to Foldseek query database during creation (default: {DEFAULT_FS_QUERYDB_SUFFIX})",
)
@click.option(
    "--id_file",
    type=click.File("rt"),
    default=None,
    required=False,
    help='Optional id list file if generating a subset of a larger folder. (default: False)'
)
@click.option(
    "--id_type",
    type=click.Choice([ID_TYPE_AF_DOMAIN, ID_TYPE_UNIPROT_DOMAIN]),
    default=ID_TYPE_AF_DOMAIN,
    help=f"Option: specify the type of ID to specify the chopping [{ID_TYPE_AF_DOMAIN}]",
)
@click.option(
    "--fs_querydb_dir",
    type=click.Path(file_okay=False, dir_okay=True, resolve_path=True),
    required=True,
    help=f"Output: directory to use for Foldseek query database during creation",
)
@click.option(
    "--fs_bin_path",
    type=click.Path(file_okay=True, resolve_path=True),
    default=FS_BINARY_PATH,
    help=f"Option: directory containing the Foldseek executable. (default: {FS_BINARY_PATH})"
)
@click.option(
    "--af_version",
    type=int,
    default=4,
    help=f"Option: specify the AF version when parsing uniprot ids",
)
def convert_cif_to_foldseek_db(
    cif_dir, fs_querydb_dir, fs_querydb_name, id_file, id_type, cif_suffix, fs_querydb_suffix, fs_bin_path, af_version
):
    "Create Foldseek query database from mmCIF folder"

    if FS_BINARY_PATH is None:
        msg = "expected foldseek binary path (FS_BINARY_PATH) to be set"
        raise RuntimeError(msg)

    fs_querydb_path = Path(fs_querydb_dir).resolve()
    if not fs_querydb_path.exists():
        os.makedirs(fs_querydb_path)
    
    if id_file is not None:
        af_tmp_dir = 'af_tmp_dir'
        if Path(af_tmp_dir).is_dir==False:
            os.mkdir(af_tmp_dir)
        for af_domain_id_str in yield_first_col(id_file):
            if id_type == ID_TYPE_UNIPROT_DOMAIN:
                af_domain_id = AFDomainID.from_uniprot_str(
                    af_domain_id_str,version=af_version
                )
            elif id_type == ID_TYPE_AF_DOMAIN:
                af_domain_id = AFDomainID.from_str(af_domain_id_str)
            else:
                msg = f"failed to understand id_type '${id_type}'"
                raise ArgumentError(msg)

            file_stub = af_domain_id.to_file_stub()
            cif_path = Path(cif_dir) / f"{file_stub}{cif_suffix}"
            if not cif_path.exists():
                msg = f"failed to locate CIF input file {cif_path}"
                LOG.error(msg)
                raise FileNotFoundError(msg)
            # Create symlinks to querydb_dir
            dest_cif_path = Path(af_tmp_dir) / cif_path.name
            if not dest_cif_path.exists():
                os.symlink(str(cif_path), str(dest_cif_path))
        cif_input_dir = af_tmp_dir
        fs_querydb = Path(f"{fs_querydb_dir}/{fs_querydb_name}{fs_querydb_suffix}")
    else:
        cif_input_dir = cif_dir
        fs_querydb = Path(f"{fs_querydb_dir}/{fs_querydb_name}{fs_querydb_suffix}")
    LOG.info(f'{cif_input_dir} {fs_querydb_dir}')
    subprocess.run(
        [
            f"{fs_bin_path}",
            "createdb",
            f"{cif_input_dir}",
            f"{fs_querydb}",
        ],
        stderr=subprocess.DEVNULL,
        check=True,
    )
    if id_file is not None:
        for root, dirs, files in os.walk(af_tmp_dir, topdown=False):
            for name in files:
                file_path = os.path.join(root, name)
                if os.path.islink(file_path):
                    os.remove(file_path)
            for name in dirs:
                dir_path = os.path.join(root, name)
                if os.path.islink(dir_path):
                    os.remove(dir_path)
            if not fs_querydb.exists():
                msg = f"failed to create expected foldseek database file: {fs_querydb_db_path}"
                raise FileNotFoundError(msg)
        os.rmdir(af_tmp_dir)

    click.echo("DONE")
    return
