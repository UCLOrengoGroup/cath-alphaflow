import logging

from cath_alphaflow.models.domains import PredictedCathDomain
from cath_alphaflow.models.domains import DecoratedCrh
from cath_alphaflow.io_utils import DecoratedCrhReader
from cath_alphaflow.errors import NoMatchingMd5Error

LOG = logging.getLogger()


class DatasetProvider:
    """
    Interface that provides `PredictedCathDomain` for a given data source
    """

    def __init__(self, datasource, *args, **kwargs):
        super().__init__(*args, **kwargs)
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


class AnnotatedCrhDatasetProvider(DatasetProvider):
    def __init__(self, *args, af_uniprot_md5_file, **kwargs):
        super().__init__(*args, **kwargs)
        self.af_uniprot_md5_file = af_uniprot_md5_file
        self.md5_to_uniprot_mapping = None

    def build_af_md5_uniprot_mapping(self):
        md5_to_uniprot = dict()
        LOG.info("Building UniProt to MD5 mapping ...")
        with self.af_uniprot_md5_file as fh:
            for line in fh:
                af_id, uniprot_id, md5 = line.strip().split()
                if md5 not in md5_to_uniprot:
                    md5_to_uniprot[md5] = []
                md5_to_uniprot[md5].extend([uniprot_id])

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

        if not self.md5_to_uniprot_mapping:
            self.build_af_md5_uniprot_mapping()

        md5_to_uniprot = self.md5_to_uniprot_mapping

        records_counter = 0
        reader = DecoratedCrhReader(self.datasource)
        for crh in reader:
            if (
                max_independent_evalue is not None
                and crh.indp_evalue > max_independent_evalue
            ):
                continue

            if records_counter > max_records:
                break

            seq_md5 = crh.sequence_md5

            if seq_md5 not in md5_to_uniprot:
                msg = f"failed to find sequence MD5 {seq_md5} in lookup"
                raise NoMatchingMd5Error(msg)

            for uniprot_id in md5_to_uniprot[seq_md5]:
                pred_dom = PredictedCathDomain(
                    uniprot_acc=uniprot_id,
                    sequence_md5=crh.sequence_md5,
                    gene3d_domain_id=crh.domain_id,
                    bitscore=crh.bitscore,
                    chopping=crh.chopping,
                    indp_evalue=crh.indp_evalue,
                )
                yield pred_dom
