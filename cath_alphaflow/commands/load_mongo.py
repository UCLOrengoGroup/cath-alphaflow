from datetime import datetime
import gzip
import logging
import io
import tarfile
import re
from typing import Collection, Dict, List, Optional
from urllib.parse import quote
from pathlib import Path

import click
import pymongo
import pydantic
from pymongo import UpdateOne
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from cath_alphaflow.settings import get_default_settings
from cath_alphaflow.models.mongo import AFFile, AFFileType
from cath_alphaflow.models import beacons
from cath_alphaflow.errors import ParseError
from pydantic import ConfigDict

LOG = logging.getLogger(__name__)

MONGO_INDEXES = ("dataset", "fileType", "afVersion", "uniprotAccession")

AF_CIF_FILE_RE = re.compile(
    "AF-(?P<up_accession>[0-9A-Z]+)-F(?P<frag_num>[0-9]+)-model_v(?P<af_version>[0-9]+)(?P<cif_suffix>\.cif\.gz)"
)

DEFAULT_BATCH_SIZE = 50
PROVIDER_ALPHAFOLD = "AlphaFold DB"
MODEL_URL_STEM = "https://alphafold.ebi.ac.uk/files"
MODEL_PAGE_URL_STEM = "https://alphafold.ebi.ac.uk/entry"


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
    help=f"Number of documents to load in a batch [default: {DEFAULT_BATCH_SIZE}]",
    default=DEFAULT_BATCH_SIZE,
    type=int,
)
@click.option(
    "--force-overwrite",
    is_flag=True,
    help="Whether or not to overwrite an existing file MongoDB [default: False]",
    default=False,
    type=bool,
)
def load_af_from_archive(
    mongo_db_url: str,
    archive_path: str,
    dataset: str,
    batch_size: int,
    force_overwrite: bool,
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
        force_overwrite=force_overwrite,
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


class AFArchiveFile(pydantic.BaseModel):
    """
    Represents a single file in the AlphaFold tar archive
    """

    filename: str
    tarfileobj: tarfile.ExFileObject
    model_config = ConfigDict(arbitrary_types_allowed=True)


def yield_next_file_from_archive(archive_path) -> AFArchiveFile:
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

        tarfileobj = tar.extractfile(tarinfo)

        yield AFArchiveFile(filename=tarinfo.name, tarfileobj=tarfileobj)


def make_af_file(
    *,
    tarfileobj,
    filename,
    dataset,
) -> AFFile:

    cif_file_match = AF_CIF_FILE_RE.match(filename)
    if not cif_file_match:
        msg = f"archive entry '{filename}' does not match expected file name"
        raise ParseError(msg)

    gzip_file = gzip.GzipFile(fileobj=tarfileobj)
    file_contents = gzip_file.read()
    gzip_file.seek(0)

    text_fh = io.TextIOWrapper(gzip_file)
    try:
        uniprot_summary = get_beacons_uniprot_summary_from_af_cif(
            text_fh, filename=filename
        )
    except Exception as err:
        msg = f"failed to parse UniprotSummary from {filename}: {err}"
        LOG.error(msg)
        raise

    up_accession = cif_file_match.group("up_accession")

    af_file = AFFile(
        fileName=filename,
        fileType=AFFileType.MODEL_CIF,
        dataset=dataset,
        afVersion=cif_file_match.group("af_version"),
        fragNum=cif_file_match.group("frag_num"),
        uniprotAccession=up_accession,
        contents=file_contents,
        uniprot_summary=uniprot_summary,
    )

    return af_file


def get_beacons_uniprot_summary_from_af_cif(
    cif_fh: io.TextIOWrapper,
    filename: str,
) -> beacons.UniprotSummary:
    """
    Returns the 3D-Beacons compliant UniprotSummary object from AlphaFold CIF file

    Note: this function relies on metadata and that should be present in all v4 AF files,
    however they may not be present in files that have been processed.
    """

    cif_file_match = AF_CIF_FILE_RE.match(filename)
    if not cif_file_match:
        raise ParseError(f"failed to parse cif filename '{filename}'")

    cif_suffix = cif_file_match.group("cif_suffix")
    af_id = filename
    if af_id.endswith(cif_suffix):
        af_id = af_id[-(len(cif_suffix)) :]

    unp_accession = cif_file_match.group("up_accession")

    cif_dict = MMCIF2Dict(cif_fh)

    # db_accession                 P00520
    # db_code                      ABL1_MOUSE
    # db_name                      UNP
    # gene_name                    Abl1
    # ncbi_taxonomy_id             10090
    # organism_scientific          "Mus musculus"
    # seq_db_align_begin           1
    # seq_db_align_end             1123
    # seq_db_isoform               ?
    # seq_db_sequence_checksum     BD48ADE8557AE95C
    # seq_db_sequence_version_date 2005-02-15
    # target_entity_id             1

    # _ma_qa_metric_global.metric_value 64.96
    global_plddt_score = float(cif_dict["_ma_qa_metric_global.metric_value"][0])
    unp_start = int(cif_dict["_ma_target_ref_db_details.seq_db_align_begin"][0])
    unp_end = int(cif_dict["_ma_target_ref_db_details.seq_db_align_end"][0])
    unp_id = cif_dict["_ma_target_ref_db_details.db_code"][0]
    unp_checksum = cif_dict["_ma_target_ref_db_details.seq_db_sequence_checksum"][0]

    # not part of 3D Beacons, but might be useful...
    taxon_id = cif_dict["_ma_target_ref_db_details.ncbi_taxonomy_id"][0]
    organism = cif_dict["_ma_target_ref_db_details.organism_scientific"][0]

    pdbx_description = cif_dict["_entity.pdbx_description"][0]
    entry_date = cif_dict["_pdbx_database_status.recvd_initial_deposition_date"][0]

    beacons_summary = beacons.UniprotSummary(
        uniprot_entry=beacons.UniprotEntry(
            ac=unp_accession,
            id=unp_id,
            uniprot_checksum=unp_checksum,
            sequence_length=(unp_end - unp_start + 1),
            segment_start=unp_start,
            segment_end=unp_end,
        ),
        structures=[
            beacons.Overview(
                summary=beacons.SummaryItems(
                    model_identifier=af_id,
                    model_category=beacons.ModelCategory.AB_INITIO,
                    model_url=f"{MODEL_URL_STEM}/{af_id}.cif",
                    model_format=beacons.ModelFormat.MMCIF,
                    model_type=beacons.ModelType.ATOMIC,
                    model_page_url=f"{MODEL_PAGE_URL_STEM}/{af_id}",
                    provider=PROVIDER_ALPHAFOLD,
                    # number_of_conformers=None,
                    # ensemble_sample_url=None,
                    # ensemble_sample_format=None,
                    created=entry_date,
                    sequence_identity=100,
                    uniprot_start=unp_start,
                    uniprot_end=unp_end,
                    coverage=100,
                    experimental_method=beacons.ExperimentalMethod.THEORETICAL_MODEL,
                    # resolution=None,
                    confidence_type=beacons.ConfidenceType.pLDDT,
                    # confidence_version=None,
                    confidence_avg_local_score=global_plddt_score,
                    entities=[
                        beacons.Entity(
                            entity_type=beacons.EntityType.POLYMER,
                            entity_poly_type=beacons.EntityPolyType.POLYPEPTIDE_L_,
                            identifier=unp_accession,
                            identifier_category=beacons.IdentifierCategory.UNIPROT,
                            description=pdbx_description,
                            chain_ids=["A"],
                        )
                    ],
                )
            )
        ],
    )
    return beacons_summary


def run(
    *,
    archive_path: str,
    dataset: str,
    mongo_db_url: str,
    batch_size: int,
    force_overwrite: bool = False,
):
    """Load AlphaFold tar archive model files into MONGO

    Args:
        archive_path (str): Path to the tar archive file
        mongo_db_url (str): Mongo DB URL
        dataset (str): Name to associate with these files
        batch_size (int): Number of documents to batch in a single commit
        force_overwrite (bool): Whether to overwrite existing documents

    """

    lm = MongoLoad()

    LOG.info(f"Initiating Mongo collection {mongo_db_url}")
    lm.init_collection(mongo_db_url)

    LOG.info(f"Loading all model files from {archive_path}")
    for af_archive_file in yield_next_file_from_archive(archive_path=archive_path):

        af_filename = af_archive_file.filename
        af_tarfileobj = af_archive_file.tarfileobj

        total = incr = 0

        unique_criteria = {"dataset": dataset, "fileName": af_filename}

        if not force_overwrite:
            if lm.collection.find_one(unique_criteria):
                LOG.info(f"Skipping {af_filename} as it already exists")
                continue

        af_file = make_af_file(
            tarfileobj=af_tarfileobj, filename=af_filename, dataset=dataset
        )

        af_file_dict = af_file.dict()

        LOG.info(f"Adding AF Model: {af_file}")

        lm.data.append(
            UpdateOne(
                unique_criteria,
                {"$set": af_file_dict},
                upsert=True,
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

    # lm.create_index()

    return 0
