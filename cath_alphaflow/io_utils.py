import csv
import itertools


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


def chunked_iterable(iterable, *, chunk_size):
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, chunk_size))
        if not chunk:
            break
        yield chunk
