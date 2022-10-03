import csv
import logging
import itertools

from .models import AFChainID
from .models import AFDomainID

LOG = logging.getLogger(__name__)


def get_csv_dictwriter(csvfile, fieldnames, **kwargs):
    """Common CSV writer"""
    return csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter="\t", **kwargs)


def get_csv_dictreader(csvfile, **kwargs):
    """Common CSV reader"""
    return csv.DictReader(csvfile, delimiter="\t", **kwargs)


def get_uniprot_id_dictreader(csvfile, **kwargs):
    reader = get_csv_dictreader(csvfile, **kwargs)
    next(reader)  # header
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


class AFDomainIDReader(csv.DictReader):
    def __init__(self, *args):
        self._seen_header = False
        super().__init__(*args, fieldnames=["af_domain_id"])

    def __next__(self):
        dictrow = super().__next__()
        if not self._seen_header:
            self._seen_header = True
            return dictrow
        else:
            return AFDomainID.from_str(dictrow["af_domain_id"])


class AFChainIDReader(csv.DictReader):
    def __init__(self, *args):
        LOG.info("AFChainIDReader.args: %s", args)
        self._seen_header = False
        super().__init__(*args, fieldnames=["af_chain_id"])
        LOG.info("AFChainIDReader.self: %s", self)

    def __next__(self):
        dictrow = super().__next__()
        if not self._seen_header:
            self._seen_header = True
            return dictrow
        else:
            return AFChainID.from_str(dictrow["af_chain_id"])


def yield_first_col(infile, *, header=True):
    if header:
        next(infile)
    for line in infile:
        first_col = line.split()[0].strip()  # take ID from first column
        yield first_col


def get_af_domain_id_reader(csvfile):
    reader = AFDomainIDReader(csvfile)
    next(reader)
    return reader


def get_af_chain_id_reader(csvfile):
    reader = AFChainIDReader(csvfile)
    next(reader)
    return reader


def chunked_iterable(iterable, *, chunk_size):
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, chunk_size))
        if not chunk:
            break
        yield chunk
