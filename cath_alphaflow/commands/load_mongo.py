import gzip
import json
import logging
import io
import os
import tempfile
import tarfile
import re
from typing import Collection, Dict, List
from urllib.parse import quote

import click
import pymongo
from pymongo import UpdateOne

from cath_alphaflow.settings import get_default_settings
from cath_alphaflow.models.mongo import AFFile, AFFileType

LOG = logging.getLogger(__name__)

MONGO_INDEXES = ("dataset", "fileType", "afVersion", "uniprotAccession")

AF_CIF_FILE_RE = re.compile(
    "AF-(?P<up_accession>[0-9A-Z]+)-F(?P<frag_num>[0-9]+)-model_v(?P<af_version>[0-9]+).cif.gz"
)


@click.command("load-mongo-archive")
@click.option(
    "-i",
    "--archive-path",
    help="Path to archive tar file.",
    required=True,
)
@click.option(
    "-h",
    "--mongo-db-url",
    help="Mongo DB URL",
    required=False,
)
@click.option(
    "-d",
    "--dataset",
    help="Dataset name",
    required=True,
)
@click.option(
    "-b",
    "--batch-size",
    help="Number of documents to load in a batch, default 1000",
    required=False,
    default=1000,
    type=int,
)
def load_af_from_archive(
    mongo_db_url: str, archive_path: str, dataset: str, batch_size: int
):  # pragma: no cover

    settings = get_default_settings()

    if not mongo_db_url:
        mongo_db_url = (
            f"mongodb://{quote(settings.MONGO_USERNAME, safe='')}:"
            f"{quote(settings.MONGO_PASSWORD, safe='')}"
            f"@{settings.MONGO_HOST}"
        )

    run(
        archive_path=archive_path,
        mongo_db_url=mongo_db_url,
        dataset=dataset,
        batch_size=batch_size,
    )


class MongoLoad:

    data: List[Dict]
    collection: Collection

    def __init__(self) -> None:
        self.data = []
        self.key_fields = list((x, "text") for x in MONGO_INDEXES)

    def init_collection(self, mongo_db_url):
        self.collection = pymongo.MongoClient(mongo_db_url).models.afCollection

    def load(self):
        self.collection.bulk_write(self.data)

    def create_index(self):
        LOG.info("Creating index")
        self.collection.create_index(self.key_fields)


def yield_archive_file(archive_path, *, dataset) -> AFFile:
    tar = tarfile.open(archive_path)
    for tarinfo in tar:
        if not tarinfo.isfile():
            LOG.debug(f"archive entry '{tarinfo.name}' is not a valid file (skipping)")
            continue

        cif_file_match = AF_CIF_FILE_RE.match(tarinfo.name)
        if not cif_file_match:
            LOG.debug(
                f"archive entry '{tarinfo.name}' does not match expected file name (skipping)"
            )
            continue

        tar_file = tar.extractfile(tarinfo)
        gzip_file = gzip.GzipFile(fileobj=tar_file)
        file_contents = gzip_file.read()

        af_file = AFFile(
            id=tarinfo.name,
            dataset=dataset,
            fileType=AFFileType.MODEL_CIF,
            afVersion=cif_file_match.group("af_version"),
            fragNum=cif_file_match.group("frag_num"),
            uniprotAccession=cif_file_match.group("up_accession"),
            contents=file_contents,
        )

        yield af_file


def run(*, archive_path: str, dataset: str, mongo_db_url: str, batch_size: int):
    """Load AlphaFold tar archive model files into MONGO

    Args:
        archive_path (str): Path to the tar archive file
        mongo_db_url (str): Mongo DB URL
        dataset (str): Name to associate with these files
        batch_size (int): Number of documents to batch in a single commit
    """

    lm = MongoLoad()

    LOG.info(f"Initiating Mongo collection {mongo_db_url}")
    lm.init_collection(mongo_db_url)

    LOG.info(f"Loading all model files from {archive_path}")
    for af_file in yield_archive_file(archive_path=archive_path, dataset=dataset):

        total = incr = 0

        LOG.info(f"AF Model file: {af_file}")

        af_file_dict: dict = af_file.dict()

        lm.data.append(
            UpdateOne(
                {"_id": af_file_dict.get("_id")}, {"$set": af_file_dict}, upsert=True
            )
        )
        incr += 1
        if incr == batch_size:
            total += incr
            incr = 0
            lm.load()
            lm.data.clear()
            LOG.info(f"Loading done: {incr} documents")

        if lm.data:
            lm.load()
            LOG.info(f"Loading done: {incr} documents")

    lm.create_index()

    return 0
