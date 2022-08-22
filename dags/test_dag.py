
# [START Test_DAG for downloading and flattening a FASTA file from UniProt]
# [START import_module]
from datetime import datetime, timedelta
from textwrap import dedent

# The DAG object; we'll need this to instantiate a DAG
from airflow import DAG

# Operators; we need this to operate!
from airflow.operators.bash import BashOperator

# [END import_module]


# [START instantiate_dag]
with DAG(
    'test_download_flatten_fasta',
    # [START default_args]
    # These args will get passed on to each operator
    # You can override them on a per-task basis during operator initialization
    default_args={
        'depends_on_past': False,
        'email': ['airflow@example.com'],
        'email_on_failure': False,
        'email_on_retry': False,
        'retries': 1,
        'retry_delay': timedelta(minutes=5),
        # 'queue': 'bash_queue',
        # 'pool': 'backfill',
        # 'priority_weight': 10,
        # 'end_date': datetime(2016, 1, 1),
        # 'wait_for_downstream': False,
        # 'sla': timedelta(hours=2),
        # 'execution_timeout': timedelta(seconds=300),
        # 'on_failure_callback': some_function,
        # 'on_success_callback': some_other_function,
        # 'on_retry_callback': another_function,
        # 'sla_miss_callback': yet_another_function,
        # 'trigger_rule': 'all_success'
    },
    # [END default_args]
    description='Fetches the Human UniProt Proteome, Unzips and Flattens the Fasta Files',
    schedule_interval=timedelta(days=1),
    start_date=datetime(2021, 1, 1),
    catchup=False,
    tags=['example'],
) as dag:
    # [END instantiate_dag]

    # t1, t2 and t3 are examples of tasks created by instantiating operators
    # [START basic_task]
    t1 = BashOperator(
        task_id='wget_proteome_seq',
        bash_command='wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz -P /opt/airflow/sequences/',
    )

    t2 = BashOperator(
        task_id='wget_proteome_struct',
        bash_command='wget https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v3.tar -P /opt/airflow/structures/',
    )
    
    t3 = BashOperator(
        task_id='unzip_seq',
        depends_on_past=False,
        bash_command='gunzip -f /opt/airflow/sequences/UP000005640_9606.fasta.gz',
        retries=3,
    )

    t4 = BashOperator(
        task_id='unzip_struct',
        bash_command='tar xzf opt/airflow/structures/UP000005640_9606_HUMAN_v3.tar',
    )

    t5 = BashOperator(
        task_id='flatten_fasta',
        depends_on_past=False,
        bash_command='sleep 5',
    )
    
    t6 = BashOperator(
        task_id='list_structures',
        depends_on_past=False,
        bash_command='sleep 5',
        retries=3,
    )
    # [END basic_task]

    # [START documentation]
    t1.doc_md = dedent(
        """\
    #### Task 1
    Downloads the proteome of Homo sapiens using the UniProt FTP
    """
    )

    t2.doc_md = dedent(
        """\
        ### Task 2
        Download the AF2 human proteome using AFDB
        """
    )

    t3.doc_md = dedent(
        """\
        ### Task 3
        Unzip File with human proteome sequence
        """
    )

    t4.doc_md = dedent(
        """\
        ### Task 4
        Unzip File with Human Proteome structures
        """
    )

    t5.doc_md = dedent(
        """\
        ### Task 5
        Flatten Fasta File
        """
    )

    t6.doc_md = dedent(
        """\
        ### Task 6
        Create list of structures
        """
    )

    dag.doc_md = __doc__  # providing that you have a docstring at the beginning of the DAG
    dag.doc_md = """
    This is a documentation placed anywhere
    """  # otherwise, type it like this
    # [END documentation]

    # [START jinja_template]
    templated_command = dedent(
        """
    {% for i in range(5) %}
        echo "{{ ds }}"
        echo "{{ macros.ds_add(ds, 7)}}"
    {% endfor %}
    """
    )

    # [END jinja_template]

    t1 >> t3 >> t5
    t2 >> t4 >> t6
# [END tutorial]