import csv
import itertools


def get_csv_dictwriter(csvfile, fieldnames, **kwargs):
    """Common CSV writer"""
    return csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter="\t", **kwargs)


def get_csv_dictreader(csvfile, **kwargs):
    """Common CSV reader"""
    return csv.DictReader(csvfile, delimiter="\t", **kwargs)


def chunked_iterable(iterable, *, chunk_size):
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, chunk_size))
        if not chunk:
            break
        yield chunk
