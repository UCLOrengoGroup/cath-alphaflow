import csv
import gzip
import logging
from pathlib import Path
import itertools
from typing import List, Type
import dataclasses

from Bio.PDB import PDBParser
from Bio.PDB import Structure

from .models.domains import AFChainID
from .models.domains import AFDomainID
from .models.domains import GeneralDomainID
from .models.domains import DecoratedCrh
from .models.domains import Gene3DCrh
from .models.domains import StatusLog
from .models.domains import RE_UNIPROT_ID
from .models.domains import FoldseekSummary
from .errors import CsvHeaderError

LOG = logging.getLogger(__name__)


class CsvReaderBase(csv.DictReader):
    """
    Generic CSV reader that can return rows as objects

    Args:
        fieldnames (`List[str]`): names of the CSV fields
        object_class (`Type[object]`): return each row as an instance of this class

    If `object_class` is set then the reader will attempt to yield an instance of
    the object class using the `dict` provided by the CSV row. If `object_class` is
    not set then the re will just return a `dict`.

    CSV fieldnames can be set explicitly via `fieldnames`. If this is not set
    then the fieldnames will be inferred from the attributes of `object_class`.

    """

    object_class: Type[object] = None
    fieldnames: List[str] = None
    has_header: bool = True

    DEFAULT_DELIMITER: str = ","

    def __init__(
        self,
        f,
        *args,
        **kwargs,
    ):
        if "delimiter" not in kwargs:
            kwargs["delimiter"] = self.get_default_delimiter()

        super().__init__(f, *args, **kwargs)
        self._seen_header = False

        if self.fieldnames is None:
            if "fieldnames" in kwargs:
                fieldnames = kwargs["fieldnames"]
            else:
                fieldnames = [
                    f.name
                    for f in self.get_result_object_fields()
                    if not f.name.startswith("_")
                ]
            self.fieldnames = fieldnames

    def get_default_delimiter(self):
        return self.DEFAULT_DELIMITER

    def get_result_object_fields(self):
        # pydantic model
        if hasattr(self.object_class, "__fields__"):
            fields = list(self.object_class.__fields__.values())
        # dataclasses
        else:
            fields = dataclasses.fields(self.object_class)
        return fields

    def __next__(self):
        """
        Checks the headers and handles converting the CSV row to object
        """

        dictrow = super().__next__()
        if not self._seen_header and self.has_header:
            self._seen_header = True
            if list(dictrow.keys()) != self.fieldnames:
                msg = (
                    f"expected first line of {self.__class__.__name__} "
                    f"to contain fieldnames {self.fieldnames}, "
                    f"but found {list(dictrow.keys())}) (reader: {self.reader})"
                )
                raise CsvHeaderError(msg)
            dictrow = super().__next__()

        if not self.object_class:
            return dictrow

        obj = self.dict_to_obj(dictrow)

        return obj

    def dict_to_obj(self, row: dict):
        return self.object_class(**row)


class StatusLogReader(CsvReaderBase):
    object_class = StatusLog


def get_status_log_dictwriter(csvfile, **kwargs):
    fieldnames = ["entry_id", "status", "error", "description"]
    writer = get_csv_dictwriter(csvfile, fieldnames=fieldnames, **kwargs)
    writer.writeheader()
    return writer


class AFDomainIDReader(CsvReaderBase):
    object_class = AFDomainID
    fieldnames = ["af_domain_id"]

    def dict_to_obj(self, row: dict):
        return self.object_class.from_str(row["af_domain_id"])


class GeneralDomainIDReader(CsvReaderBase):
    object_class = GeneralDomainID
    fieldnames = ["domain_id"]

    def dict_to_obj(self, row: dict):
        return self.object_class.from_str(row["domain_id"])


class AFChainIDReader(CsvReaderBase):
    object_class = AFChainID
    fieldnames = ["af_chain_id"]

    def dict_to_obj(self, row: dict):
        return self.object_class.from_str(row["af_chain_id"])


class DecoratedCrhReader(CsvReaderBase):
    object_class = DecoratedCrh
    has_header = False


class Gene3DCrhReader(CsvReaderBase):
    object_class = Gene3DCrh
    has_header = False

    def get_default_delimiter(self):
        return "\t"


class UniprotIDReader(CsvReaderBase):
    fieldnames = ["uniprot_id"]

    def dict_to_obj(self, row: dict):
        uniprot_id = row["uniprot_id"]
        assert RE_UNIPROT_ID.match(uniprot_id)
        return {"uniprot_id": uniprot_id}


def get_csv_dictwriter(csvfile, fieldnames, delimiter="\t", **kwargs):
    """Common CSV writer"""
    return csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=delimiter, **kwargs)


def get_csv_dictreader(csvfile, delimiter="\t", **kwargs):
    """Common CSV reader"""
    return csv.DictReader(csvfile, delimiter=delimiter, **kwargs)


def get_uniprot_id_dictreader(csvfile, **kwargs):
    reader = UniprotIDReader(csvfile)
    return reader


def get_uniprot_id_dictwriter(csvfile, **kwargs):
    writer = get_csv_dictwriter(csvfile, fieldnames=["uniprot_id"], **kwargs)
    writer.writeheader()
    return writer


def get_sse_summary_reader(csvfile):
    reader = get_csv_dictreader(csvfile)
    next(reader)
    return reader


def get_foldseek_reader(csvfile):
    foldseek_fieldnames = [
        "query",
        "target",
        "qstart",
        "qend",
        "qlen",
        "tstart",
        "tend",
        "tlen",
        "qcov",
        "tcov",
        "bits",
        "evalue",
    ]
    foldseek_reader = get_csv_dictreader(csvfile, fieldnames=foldseek_fieldnames)
    return foldseek_reader


def get_foldseek_summary_writer(csvfile):
    writer = get_csv_dictwriter(
        csvfile,
        fieldnames=[
            "query",
            "target",
            "qstart",
            "qend",
            "qlen",
            "tstart",
            "tend",
            "tlen",
            "qcov",
            "tcov",
            "bits",
            "evalue",
        ],
    )
    writer.writeheader()
    return writer


def get_sse_summary_writer(csvfile):
    writer = get_csv_dictwriter(
        csvfile,
        fieldnames=[
            "af_domain_id",
            "ss_res_total",
            "res_count",
            "perc_not_in_ss",
            "sse_H_num",
            "sse_E_num",
            "sse_num",
        ],
    )
    writer.writeheader()
    return writer


def get_af_uniprot_md5_summary_writer(csvfile):
    writer = get_csv_dictwriter(
        csvfile,
        fieldnames=[
            "af_chain_id",
            "uniprot_id",
            "sequence_md5",
        ],
    )
    writer.writeheader()
    return writer


def get_plddt_summary_writer(csvfile):
    writer = get_csv_dictwriter(
        csvfile,
        fieldnames=[
            "af_domain_id",
            "avg_plddt",
            "perc_LUR",
            "residues_total",
        ],
    )
    writer.writeheader()
    return writer


def yield_first_col(infile, *, header=True):
    if header:
        next(infile)
    for line in infile:
        first_col = line.split()[0].strip()  # take ID from first column
        yield first_col


def get_af_domain_id_reader(csvfile):
    reader = AFDomainIDReader(csvfile)
    return reader


def get_af_chain_id_reader(csvfile):
    reader = AFChainIDReader(csvfile)
    return reader


def get_general_domain_id_reader(csvfile):
    reader = GeneralDomainIDReader(csvfile)
    return reader


def chunked_iterable(iterable, *, chunk_size):
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, chunk_size))
        if not chunk:
            break
        yield chunk


def get_pdb_structure(
    model_id, chain_pdb_dir, chains_are_gzipped=False, model_filename=None
) -> Structure:
    """Return a Bio.PDB.Structure for a given domain (model_id / chopping)"""

    # create default filename
    if model_filename is None:
        if not chains_are_gzipped:
            open_func = open
            model_filename = model_id.raw_id + ".pdb"
        else:
            open_func = gzip.open
            model_filename = model_id.raw_id + ".pdb.gz"

    pdb_path = Path(chain_pdb_dir, model_filename)

    parser = PDBParser(QUIET=1)

    with open_func(str(pdb_path), mode="rt") as pdb_fh:
        structure = parser.get_structure(model_id.raw_id, pdb_fh)

    return structure
