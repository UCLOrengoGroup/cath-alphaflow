#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Params defined before include are also used in the processes
params.dataset_name = 'kinfams'

// Processes to include from the shared module
include { uniprot_domain_to_uniprot } from './cath-shared-core'
include { create_af_manifest_file } from './cath-shared-core'
include { retrieve_af_chain_cif_files } from './cath-shared-core'
include { chop_cif as CHOP_CIF } from './cath-shared-core'
include { cif_paths_to_uniprot_accessions } from './cath-shared-core'
include { create_missing_uniprot_domain_ids } from './cath-shared-core'
include { ids_from_cif_files as domain_ids_from_cif_files} from './cath-shared-core'
include { ids_from_cif_files as chain_ids_from_cif_files } from './cath-shared-core'


// Params defined after include are used locally
params.dataset_dir = "$workflow.launchDir/tests/fixtures/${params.dataset_name}"
params.uniprot_domain_ids_csv_fn = "${params.dataset_name}.uniprot_domain_ids.csv"
params.cath_version = 'v4_3_0'

// Existing Params
params.gs_bucket = 'gs://public-datasets-deepmind-alphafold'
params.chunk_size = 1000

params.cath_af_cli = 'cath-af-cli'
params.cath_data_root = "/cath/data/${params.cath_version}"
params.all_af_chain_fasta_url = "${params.gs_bucket}/sequences.fasta"
params.all_af_chain_fasta_file_name = 'all_af2_chain_sequences.fasta'
params.cath_s95_pdb_dir_name = 'cath_s95_pdb'
params.cath_s95_foldseek_library_name = 'foldseek_s95_lib'
params.cath_odb_name = 'GENE3D_21'
params.dataset_max_records = '100000'
params.dataset_max_evalue = '1E-50'

params.publish_dir = "$workflow.launchDir/kinfams-results"

params.af_version = 4
params.af_download_stem = "gs://public-datasets-deepmind-alphafold-v${params.af_version}"

params.af_to_cluster_fn = "human_kinfams_uniprot_list_full"

params.uniprot_ids_csv_fn = "uniprot_ids.csv"
params.uniprot_md5_csv_fn = "uniprot_md5.csv"
params.gene3d_crh_output_fn = "gene3d_crh_output.csv"
params.af_domainlist_ids_csv_fn = "af_domainlist_ids.csv"
params.af_chainlist_ids_csv_fn = "af_chainlist_ids.csv"
params.af_cath_orig_annotations_csv_fn = "af_cath_orig_annotations.csv"
params.af_manifest_file_fn = "af_manifest_file"
params.uniprot_domain_ids_fn = "uniprot_domain_id.txt"
params.af_manifest_fn = "af_manifest.txt"
params.af_cif_raw_dir = "cif_raw"
params.af_cif_chopped_dir = "cif_chopped"
params.af_dssp_raw_dir = "dssp_raw"
params.plddt_stats_fn = "plddt_summary.csv"
params.sse_summary_fn = "sse_summary.csv"

def uniprot_ids_file = file("${params.publish_dir}/${params.uniprot_ids_csv_fn}")

process kinfam_clusters_to_uniprot_domain_ids {
    publishDir params.publish_dir, mode: 'copy'

    input:
    path cluster_path

    output:
    path params.uniprot_domain_ids_fn

    """
    cat ${cluster_path} | awk '{print \$1}' > ${params.uniprot_domain_ids_fn}
    """
}

process create_cath_dataset_from_db {
    publishDir params.publish_dir, mode: 'copy'

    input:
    path 'uniprot_ids_csv'

    output:
    path params.uniprot_md5_csv_fn
    path params.gene3d_crh_output_fn
    path params.af_domainlist_ids_csv_fn
    path params.af_chainlist_ids_csv_fn
    path params.af_cath_orig_annotations_csv_fn

    """
    cath-af-cli create-cath-dataset-from-db \
        --csv_uniprot_ids uniprot_ids_csv \
        --csv_uniprot_md5 ${params.uniprot_md5_csv_fn} \
        --gene3d_crh_output ${params.gene3d_crh_output_fn} \
        --af_domainlist_ids ${params.af_domainlist_ids_csv_fn} \
        --af_chainlist_ids ${params.af_chainlist_ids_csv_fn} \
        --af_cath_annotations ${params.af_cath_orig_annotations_csv_fn} \
        --dbname ${params.cath_odb_name} \
    """
}


process create_cath_dataset_from_files {
    publishDir params.publish_dir, mode: 'copy'

    input:
    path 'uniprot_ids_csv'
    path 'all_crh_csv'
    path 'all_af_uniprot_md5_csv'

    output:
    path params.uniprot_md5_csv_fn
    path params.gene3d_crh_output_fn
    path params.af_domainlist_ids_csv_fn
    path params.af_chainlist_ids_csv_fn
    path params.af_cath_orig_annotations_csv_fn

    """
    grep -F -f uniprot_ids_csv all_crh_csv > filtered_crh_csv
    grep -F -f uniprot_ids_csv all_af_uniprot_md5_csv > filtered_af_uniprot_md5_csv
    cath-af-cli create-cath-dataset-from-files \
        --input_af_uniprot_md5 filtered_crh_csv \
        --input_crh filtered_af_crh \
        --csv_uniprot_ids uniprot_ids_csv \
        --csv_uniprot_md5 ${params.uniprot_md5_csv_fn} \
        --gene3d_crh_output ${params.gene3d_crh_output_fn} \
        --af_domainlist_ids ${params.af_domainlist_ids_csv_fn} \
        --af_chainlist_ids ${params.af_chainlist_ids_csv_fn} \
        --af_cath_annotations ${params.af_cath_orig_annotations_csv_fn}
    """
}

process create_dssp {
    publishDir "${params.publish_dir}/${params.af_dssp_raw_dir}", mode: 'copy', overwrite: false

    input:
    path 'ids.csv'
    path cif_files

    output:
    path "*.dssp"

    """
    find . -name "*.cif" | sed 's#^./##' | sed 's/.cif\$//' > manually_created_ids.csv
    cath-af-cli convert-cif-to-dssp \
        --cif_in_dir . \
        --id_file manually_created_ids.csv \
        --cif_suffix .cif \
        --dssp_out_dir .
    """
}


process create_sse_summary {
    input:
    val domain_ids
    path dssp_files

    output:
    path 'dssp_summary.csv'

    """
    echo ${domain_ids} > af_ids_csv
    cath-af-cli convert-dssp-to-sse-summary \
        --dssp_dir . \
        --id_file af_ids_csv \
        --id_type af \
        --sse_out_file dssp_summary.csv \
        --dssp_check_policy=error
    """
}

process create_plddt_summary {

    // HACK: 
    // This currently expects to find CIF files in $publish_dir which goes
    // against best practices. If anyone can think of a neat way to coordinate 
    // the domain/chain id streams without this hack, then please change.  

    input:
    path 'af_ids_csv'
    path _chain_cif_files_list

    output:
    path params.plddt_stats_fn

    """
    # HACK
    CIF_DIR='${params.publish_dir}/${params.af_cif_raw_dir}'
    cath-af-cli convert-cif-to-plddt-summary \
        --cif_in_dir \$CIF_DIR \
        --id_file af_ids_csv \
        --plddt_stats_file ${params.plddt_stats_fn}
    """
}

process ids_from_cif_dir {
    input:
    path 'cif_dir'

    output:
    path 'ids.txt'

    """
    ls cif_dir/ | grep ".cif" | sed 's/.cif\$//' > ids.txt
    """
}

process af_domain_ids_from_cif_dir {
    input:
    path 'cif_dir'

    output:
    path 'af_models.txt'

    """
    ls cif_dir/ | grep ".cif" | sed 's/.cif\$//' > af_models.txt
    """
}

process cif_to_pdb {
    input:
    path 'tmp.cif'

    output:
    path 'tmp.pdb'

    """
    pdb_fromcif tmp.cif > tmp.pdb
    """
}

process uniprot_csv_from_af_domains {
    input:
    path af_domain_id_file

    output:
    path 'uniprot_ids.csv'

    """
    echo uniprot_acc > uniprot_ids.csv
    cat ${af_domain_id_file} | tr '-' ' ' | awk '{print \$2}' >> uniprot_ids.csv
    """
}

process filter_uniprot_id_file {
    input:
    path 'full_af_domain_id_file'
    path 'required_uniprot_ids'

    output:
    path 'filtered_af_domain_id_file'

    """
    grep -F -f required_uniprot_ids full_af_domain_id_file > filtered_af_domain_id_file
    """
}

workflow AF_DOWNLOAD_CIF {
    take:
        uniprot_domain_ids_ch
    main:
        def _cif_files = uniprot_domain_ids_ch.splitText(by: 100, file: true)
            | uniprot_domain_to_uniprot
            | create_af_manifest_file
            | retrieve_af_chain_cif_files
    emit:
        cif_files = _cif_files
}

workflow AF_CHOP_CIF {

    take:
        chain_cif_files
        uniprot_domain_ids_ch
    main:
        CHOP_CIF( chain_cif_files, uniprot_domain_ids_ch )
    emit:
        cif_files = CHOP_CIF.out
}

workflow AF_ANNOTATE_CHAINS_CIF {
    take:
        af_chain_cif_files_ch
        af_domain_
}

workflow AF_ANNOTATE_DOMAINS_CIF {

    take:
        af_chain_cif_files_ch
        af_domain_cif_files_ch

    main:
        def _plddt_file = file("${params.publish_dir}/${params.plddt_stats_fn}")
        def _sse_summary_file = file("${params.publish_dir}/${params.sse_summary_fn}")

        def af_chain_ids_ch = chain_ids_from_cif_files(af_chain_cif_files_ch)
        def af_domain_ids_ch = domain_ids_from_cif_files(af_domain_cif_files_ch)

        def uniprot_csv_ch = uniprot_csv_from_af_domains(af_domain_ids_ch)
        def uniprot_dataset = create_dataset_cath_files(uniprot_csv_ch)

        def chain_dssp_files_ch = create_dssp(af_chain_ids_ch, af_chain_cif_files_ch)

        // HACK:
        // simplest way to proceed here is to finish creating all the DSSP files
        // then assume that DSSP chain files can be found in publish_dir (against best practices)
        // so collect all these file names and make that a dependency of the next step
        def all_chain_dssp_files = chain_dssp_files_ch.collectFile(name: "chain_dssp_files.txt")

        // The pLDDT process requires meta data that is present in the chain AF CIF file.
        // However, this metadata is not present in the domain CIF file after chopping. 
        // So, while this process needs to take the AF domain IDs as input, it also needs 
        // to be given the CHAIN CIF files and a chopping (rather than the DOMAIN CIF files)
        def plddt_summary_ch = create_plddt_summary(af_domain_ids_ch, all_chain_dssp_files)
        
        // we now need to combine the channels so that we get dssp files with the same
        // uniprot ids as the dssp files

        def sse_summary_ch = create_sse_summary(af_domain_ids_ch, all_chain_dssp_files)
        sse_summary_ch.collectFile(name: _sse_summary_file)

        // TODO:
        // combine these outputs into a single spreadsheet...
    emit:
        plddt_summary_file = _plddt_file
        sse_summary_file = _sse_summary_file 
}

workflow {
    def kinfams_ch = Channel.fromPath(params.af_to_cluster_fn, checkIfExists: true)

    def uniprot_domain_ids_ch = kinfam_clusters_to_uniprot_domain_ids(kinfams_ch)

    def uniprot_csv_ch = uniprot_csv_from_af_domains(uniprot_domain_ids_ch)
    def uniprot_dataset = create_dataset_cath_files(uniprot_csv_ch)

    // def af_full_fasta_ch = retrieve_af_fasta_database()
    // def af_uniprot_md5_ch = create_af_uniprot_md5_from_fasta(af_full_fasta_ch, uniprot_domain_ids_ch)
    AF_DOWNLOAD_CIF(uniprot_domain_ids_ch)
    def chain_cif_files = AF_DOWNLOAD_CIF.out.cif_files

    AF_CHOP_CIF(chain_cif_files, uniprot_domain_ids_ch)
    def domain_cif_files = AF_CHOP_CIF.out.cif_files

    AF_ANNOTATE_DOMAINS_CIF(chain_cif_files, domain_cif_files)

    // Collect results together
    // AF_COLLECT_RESULTS(
    //     AF_ANNOTATE_DOMAINS_CIF.out.plddt_file,
    //     AF_ANNOTATE_DOMAINS_CIF.out.sse_summary_file,
    // )
}