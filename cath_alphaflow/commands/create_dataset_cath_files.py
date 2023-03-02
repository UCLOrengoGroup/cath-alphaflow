from email.policy import default
import logging
import click
import typing
import pydantic

from cath_alphaflow.io_utils import get_uniprot_id_dictreader
from cath_alphaflow.io_utils import get_csv_dictwriter
from cath_alphaflow.io_utils import chunked_iterable
from cath_alphaflow.db_utils import OraDB
from cath_alphaflow.models.domains import Chopping, AFChainID, AFDomainID
from cath_alphaflow.settings import DEFAULT_AF_VERSION, DEFAULT_AF_FRAGMENT
from cath_alphaflow.dataset_provider import AnnotatedCrhDatasetProvider

LOG = logging.getLogger()

DEFAULT_CHUNK_SIZE = 1000

SHARED_OPTIONS = [
    click.core.Option(
        ("--csv_uniprot_ids",),
        type=click.File("rt"),
        required=True,
        help="Input: CSV file of UniProt IDs to process",
    ),
    click.core.Option(
        ("--csv_uniprot_md5",),
        type=click.File("wt"),
        required=True,
        help="Output: CSV file of UniProt to MD5 mapping",
    ),
    click.core.Option(
        ("--gene3d_crh_output",),
        type=click.File("wt"),
        required=True,
        help="Output: CRH output file for Gene3D domains",
    ),
    click.core.Option(
        ("--af_domainlist_ids",),
        type=click.File("wt"),
        required=True,
        help="Output: CSV file of AF2 domain ids",
    ),
    click.core.Option(
        ("--af_chainlist_ids",),
        type=click.File("wt"),
        required=True,
        help="Output: CSV file of AF2 chain ids",
    ),
    click.core.Option(
        ("--af_cath_annotations",),
        type=click.File("wt"),
        required=True,
        help="Output: CSV file of CATH annotations",
    ),
    click.core.Option(
        ("--chunk", "chunk_size"),
        type=int,
        default=DEFAULT_CHUNK_SIZE,
        help="Param: size of chunk when processing data",
    ),
    click.core.Option(
        ("--af_version",),
        type=str,
        default=DEFAULT_AF_VERSION,
        help="AlphaFoldDB version",
    ),
    click.core.Option(
        ("--fragment_number",),
        type=str,
        default=DEFAULT_AF_FRAGMENT,
        help="default fragment number: 1",
    ),
]


class BaseCommand(click.core.Command):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for opt in SHARED_OPTIONS:
            self.params.insert(0, opt)


@click.command(cls=BaseCommand)
@click.option(
    "--dbname",
    type=str,
    required=True,
    help="Param: database to use when querying sequences",
)
def create_cath_dataset_from_db(*args, **kwargs):
    """
    Creates CATH data files for a given dataset (based on DB)
    """

    generator = CathDatasetGeneratorFromDB(*args, **kwargs)
    generator.run()


@click.command(cls=BaseCommand)
@click.option(
    "--input_af_uniprot_md5",
    type=click.File("rt"),
    required=True,
    help="Input: CSV file containing AF chain ID, UniProt IDs, MD5 lookup",
)
@click.option(
    "--input_crh",
    type=click.File("rt"),
    required=True,
    help="Input: CRH file containing matches",
)
def create_cath_dataset_from_files(*args, **kwargs):
    """
    Creates CATH data files for a given dataset (based on flat files)
    """
    generator = CathDatasetGeneratorFromFiles(*args, **kwargs)
    generator.run()


class CathDatasetGeneratorBase(pydantic.BaseModel):

    csv_uniprot_ids: typing.Any  # click.File
    gene3d_crh_output: typing.Any  # click.File
    csv_uniprot_md5: typing.Any  # click.File
    af_domainlist_ids: typing.Any  # click.File
    af_chainlist_ids: typing.Any  # click.File
    af_cath_annotations: typing.Any  # click.File
    chunk_size: int
    af_version: int
    fragment_number: int
    crh_output_writer: typing.Any = None
    csv_uniprot_md5_writer: typing.Any = None
    af_domainlist_writer: typing.Any = None
    af_chainlist_writer: typing.Any = None
    af_cath_annotations_writer: typing.Any = None

    class Config:
        arbitrary_types_allowed: True

    def run(self):

        csv_uniprot_ids = self.csv_uniprot_ids
        chunk_size = self.chunk_size

        click.echo("Setting up output file writers ...")
        self.setup_writers()

        # setup reader
        click.echo(f"Setting up UniProt reader '{csv_uniprot_ids.name}'")
        uniprot_reader = get_uniprot_id_dictreader(csv_uniprot_ids)

        # process chunks of uniprot ids
        for chunked_uniprot_rows in chunked_iterable(
            uniprot_reader, chunk_size=chunk_size
        ):
            uniprot_ids = [row.get("uniprot_id") for row in chunked_uniprot_rows]
            click.echo(
                f"Processing {len(uniprot_ids)} UniProtIDs (e.g. {uniprot_ids} ...)",
            )
            self.process_uniprot_ids(uniprot_ids)

        click.echo("DONE")

    def setup_writers(self):
        """
        Create all the file writers that we are going to use
        """

        crh_output_headers = [
            "sequence_md5",
            "cath_domain_id",
            "score",
            "boundaries",
            "resolved",
        ]
        self.crh_output_writer = get_csv_dictwriter(
            self.gene3d_crh_output, fieldnames=crh_output_headers
        )
        # NOTE: do not write headers to crh output

        csv_uniprot_md5_headers = ["uniprot_acc", "sequence_md5"]
        self.csv_uniprot_md5_writer = get_csv_dictwriter(
            self.csv_uniprot_md5, fieldnames=csv_uniprot_md5_headers
        )
        self.csv_uniprot_md5_writer.writeheader()

        af_domainlist_headers = ["af_domain_id"]
        self.af_domainlist_writer = get_csv_dictwriter(
            self.af_domainlist_ids, fieldnames=af_domainlist_headers
        )
        self.af_domainlist_writer.writeheader()

        af_chainlist_headers = ["af_chain_id"]
        self.af_chainlist_writer = get_csv_dictwriter(
            self.af_chainlist_ids, fieldnames=af_chainlist_headers
        )
        self.af_chainlist_writer.writeheader()

        af_cath_annotations_headers = [
            "cath_domain_id",
            "uniprot_acc",
            "md5",
            "bitscore",
            "chopping",
        ]
        self.af_cath_annotations_writer = get_csv_dictwriter(
            self.af_cath_annotations, fieldnames=af_cath_annotations_headers
        )
        self.af_cath_annotations_writer.writeheader()

    def process_uniprot_ids(self, uniprot_ids):
        """
        General method to process uniprot ids
        """

        for entry in self.next_cath_dataset_entry(
            uniprot_ids=uniprot_ids,
        ):

            # sort out all variables that we are going to use in our data files
            uniprot_acc = entry.uniprot_acc
            sequence_md5 = entry.sequence_md5
            gene3d_domain_id = entry.gene3d_domain_id
            bitscore = entry.bitscore
            chopping = entry.chopping

            af_domain_id = AFDomainID(
                uniprot_acc=entry.uniprot_acc,
                fragment_number=self.fragment_number,
                version=self.af_version,
                chopping=Chopping.from_str(chopping_str=chopping),
            ).to_str()

            af_chain_id = AFChainID(
                uniprot_acc=uniprot_acc,
                fragment_number=self.fragment_number,
                version=self.af_version,
            ).to_str()

            # write data
            self.csv_uniprot_md5_writer.writerow(
                {"uniprot_acc": uniprot_acc, "sequence_md5": sequence_md5}
            )
            self.af_domainlist_writer.writerow({"af_domain_id": af_domain_id})
            self.af_chainlist_writer.writerow({"af_chain_id": af_chain_id})
            self.af_cath_annotations_writer.writerow(
                {
                    "cath_domain_id": gene3d_domain_id,
                    "uniprot_acc": uniprot_acc,
                    "md5": sequence_md5,
                    "bitscore": bitscore,
                    "chopping": chopping,
                }
            )

            # use the 'correct' chopping in both cols
            self.crh_output_writer.writerow(
                {
                    "cath_domain_id": gene3d_domain_id,
                    "sequence_md5": sequence_md5,
                    "score": bitscore,
                    "boundaries": chopping,
                    "resolved": chopping,
                }
            )


class CathDatasetGeneratorFromDB(CathDatasetGeneratorBase):
    dbname: str

    def next_cath_dataset_entry(self, uniprot_ids):
        db = OraDB()
        for entry in db.next_cath_dataset_entry(
            dbname=self.dbname,
            uniprot_ids=uniprot_ids,
        ):
            yield entry


class CathDatasetGeneratorFromFiles(CathDatasetGeneratorBase):

    src_crh_file: typing.Any  # click.File
    src_af_uniprot_md5: typing.Any  # click.File

    def next_cath_dataset_entry(self, uniprot_ids):
        provider = AnnotatedCrhDatasetProvider(
            datasource=self.src_crh_file,
            af_uniprot_md5_file=self.src_af_uniprot_md5,
        )

        for entry in provider.next_cath_dataset_entry(
            uniprot_ids=uniprot_ids,
        ):
            yield entry
