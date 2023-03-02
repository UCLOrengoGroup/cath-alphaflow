import csv
import logging
import itertools
from typing import List, Type
from dataclasses import fields

from .models.domains import AFChainID
from .models.domains import AFDomainID
from .models.domains import DecoratedCrh
from .models.domains import RE_UNIPROT_ID
from .errors import CsvHeaderError

LOG = logging.getLogger(__name__)


class CsvReaderBase(csv.DictReader):
    """
    Generic CSV reader that maps rows to objects
    """

    object_class: Type[object] = None
    fieldnames: List[str] = None

    def __init__(
        self,
        f,
        *args,
        **kwargs,
    ):
        super().__init__(f, *args, **kwargs)
        self._seen_header = False

        if not self.fieldnames:
            if "fieldnames" in kwargs:
                fieldnames = kwargs["fieldnames"]
            else:
                fieldnames = [
                    f.name
                    for f in fields(self.object_class)
                    if not f.name.startswith("_")
                ]
            self.fieldnames = fieldnames

    def __next__(self):
        """
        Checks the headers and handles converting the CSV row to object
        """

        dictrow = super().__next__()
        if not self._seen_header:
            self._seen_header = True
            if list(dictrow.keys()) != self.fieldnames:
                msg = (
                    f"expected first line of {self.__class__.__name__} "
                    f"to contain fieldnames {self.fieldnames}, "
                    f"but found {dictrow.keys()}) (reader: {self.reader})"
                )
                raise CsvHeaderError(msg)
            dictrow = super().__next__()

        if not self.object_class:
            return dictrow

        return self.dict_to_obj(dictrow)

    def dict_to_obj(self, row: dict):
        return self.return_class(**row)


class AFDomainIDReader(CsvReaderBase):
    object_class = AFDomainID
    fieldnames = ["af_domain_id"]

    def dict_to_obj(self, row: dict):
        return self.object_class.from_str(row["af_domain_id"])


class AFChainIDReader(CsvReaderBase):
    object_class = AFChainID
    fieldnames = ["af_chain_id"]

    def dict_to_obj(self, row: dict):
        return self.object_class.from_str(row["af_chain_id"])


class DecoratedCrhReader(CsvReaderBase):
    object_class = DecoratedCrh


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


def chunked_iterable(iterable, *, chunk_size):
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, chunk_size))
        if not chunk:
            break
        yield chunk
