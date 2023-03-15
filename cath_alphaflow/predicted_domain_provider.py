"""
Classes to generate `PredictedCathDomain` objects from databases or files
"""

import logging

import pydantic
from io import TextIOWrapper

from cath_alphaflow.models.domains import PredictedCathDomain
from cath_alphaflow.io_utils import DecoratedCrhReader, Gene3DCrhReader
from cath_alphaflow.errors import NoMatchingMd5Error

LOG = logging.getLogger()


class PredictedCathDomainProviderBase:
    """
    Interface that provides `PredictedCathDomain` for a given data source
    """

    def __init__(self, datasource, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if datasource is None:
            raise ValueError("datasource must be set")
        self.datasource = datasource

    def next_cath_dataset_entry(
        self,
        *,
        max_independent_evalue=None,
        max_records=None,
        uniprot_ids=None,
        **kwargs,
    ) -> PredictedCathDomain:
        raise NotImplementedError


class AfMd5Uniprot(pydantic.BaseModel):
    af_id: str
    md5: str
    uniprot_id: str


class OraclePredictedCathDomainProvider(PredictedCathDomainProviderBase):
    """
    Provides datasets from Oracle database
    """

    def next_cath_dataset_entry(
        self,
        *,
        dbname,
        max_independent_evalue=None,
        max_records=None,
        uniprot_ids=None,
    ) -> PredictedCathDomain:
        """
        Returns a generator that provides `PredictedCathDomain` entries
        """

        if not max_records and not uniprot_ids:
            raise RuntimeError("need to specify one of [max_records, uniprot_ids]")

        # organise filters in the where clause
        # (SQL placeholders and cooresponding substitutions key/values)
        sql_args = {}
        sql_where_args = []
        if uniprot_ids:

            uniprot_id_placeholders = {
                f"u{num}": uniprot_id for num, uniprot_id in enumerate(uniprot_ids, 1)
            }
            sql_args.update(uniprot_id_placeholders)
            sql_where_args += [
                (
                    "upa.ACCESSION IN ("
                    + ", ".join([f":u{num}" for num in range(1, len(uniprot_ids) + 1)])
                    + ")"
                )
            ]

        if max_independent_evalue:
            sql_args.update({"max_independent_evalue": max_independent_evalue})
            sql_where_args += ["INDEPENDENT_EVALUE <= :max_independent_evalue"]

        if max_records:
            sql_args.update({"max_records": max_records})
            sql_where_args += ["ROWNUM <= :max_records"]

        # join all where clauses together with 'AND' operator
        sql_where = " AND ".join(sql_where_args)

        sql = f"""
    SELECT DISTINCT
        upa.ACCESSION                           AS uniprot_acc,
        upa.SEQUENCE_MD5                        AS sequence_md5,
        DOMAIN_ID || '__' || SUPERFAMILY || '/'
            || REPLACE(RESOLVED, ',', '_')      AS gene3d_domain_id,
        SCORE                                   AS bitscore,
        RESOLVED                                AS chopping,
        INDEPENDENT_EVALUE                      AS indp_evalue
    FROM
        {dbname}.CATH_DOMAIN_PREDICTIONS cdp
        INNER JOIN {dbname}.UNIPROT_PRIM_ACC upa
            ON (cdp.SEQUENCE_MD5 = upa.SEQUENCE_MD5)
    WHERE
        {sql_where}
    ORDER BY
        INDEPENDENT_EVALUE ASC, upa.ACCESSION ASC
    """

        # execute query and yield results (as dict) row by row
        for rowdict in self.yieldall(sql, sql_args, return_type=dict):
            # convert dict to data structure
            entry = PredictedCathDomain(**rowdict)
            yield entry


class CrhPredictedCathDomainProviderBase(PredictedCathDomainProviderBase):
    """
    Provides datasets from CRH files ("original" or "decorated")
    """

    def __init__(self, *args, af_uniprot_md5_file, **kwargs):
        super().__init__(*args, **kwargs)
        self.af_uniprot_md5_file = af_uniprot_md5_file
        self.md5_to_af_md5_uniprot_mapping = None

    def get_datasource_reader(self):
        raise NotImplementedError

    def build_md5_to_af_md5_uniprot_mapping(self):
        md5_to_af_md5_uniprot_mapping = dict()
        LOG.info(f"Building UniProt to MD5 mapping ... {self.af_uniprot_md5_file}")

        def process_line(line):
            af_id, uniprot_id, md5 = line.strip().split()
            entry = AfMd5Uniprot(af_id=af_id, md5=md5, uniprot_id=uniprot_id)
            LOG.info(f"entry: {entry}")
            if md5 not in md5_to_af_md5_uniprot_mapping:
                md5_to_af_md5_uniprot_mapping[md5] = []
            md5_to_af_md5_uniprot_mapping[md5].extend([entry])

        if isinstance(self.af_uniprot_md5_file, TextIOWrapper):
            for line in self.af_uniprot_md5_file:
                process_line(line)
        else:
            with open(self.af_uniprot_md5_file) as fh:
                for line in self.af_uniprot_md5_file:
                    process_line(line)

        self.md5_to_af_md5_uniprot_mapping = md5_to_af_md5_uniprot_mapping

    def next_cath_dataset_entry(
        self,
        *,
        max_independent_evalue=None,
        max_records=None,
        uniprot_ids=None,
        **kwargs,
    ) -> PredictedCathDomain:
        if uniprot_ids is not None:
            uniprot_ids = set(uniprot_ids)

        if not self.md5_to_af_md5_uniprot_mapping:
            self.build_md5_to_af_md5_uniprot_mapping()

        records_counter = 0
        reader = self.get_datasource_reader()

        for crh in reader:
            if (
                max_independent_evalue is not None
                and crh.indp_evalue > max_independent_evalue
            ):
                continue

            if max_records and records_counter > max_records:
                break

            seq_md5 = crh.sequence_md5

            if seq_md5 not in self.md5_to_af_md5_uniprot_mapping:
                msg = f"failed to find sequence MD5 {seq_md5} in lookup"
                raise NoMatchingMd5Error(msg)

            for af_md5_uniprot in self.md5_to_af_md5_uniprot_mapping[seq_md5]:
                pred_dom = PredictedCathDomain(
                    uniprot_acc=af_md5_uniprot.uniprot_id,
                    sequence_md5=crh.sequence_md5,
                    gene3d_domain_id=crh.domain_id,
                    bitscore=crh.bitscore,
                    chopping=crh.chopping_final,
                    indp_evalue=crh.indp_evalue,
                )
                yield pred_dom


class Gene3DCrhPredictedCathDomainProvider(CrhPredictedCathDomainProviderBase):
    """
    Provides a CATH dataset from CRH / MD5 files
    """

    def get_datasource_reader(self):
        return Gene3DCrhReader(self.datasource)


class DecoratedCrhPredictedCathDomainProvider(CrhPredictedCathDomainProviderBase):
    """
    Provides a CATH dataset from "decorated" CRH / MD5 files
    """

    def get_datasource_reader(self):
        return DecoratedCrhReader(self.datasource)
