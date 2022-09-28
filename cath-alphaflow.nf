#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.cath_af_cli = 'cath-af-cli'
params.cath_version = 'v4_3_0'
params.cath_data_root = "/cath/data/${params.cath_version}"
params.all_af_chain_fasta_url = 'gs://public-datasets-deepmind-alphafold/sequences.fasta'
params.all_af_chain_fasta_file_name = 'all_af2_chain_sequences.fasta'
params.cath_s95_pdb_dir_name = 'cath_s95_pdb'
params.cath_s95_foldseek_library_name = 'foldseek_s95_lib'
params.cath_odb_name = 'GENE3D_21'
params.dataset_name = 'dataset100k'
params.dataset_max_records = '100000'
params.dataset_max_evalue = '1E-50'

params.uniprot_ids_csv = "${params.dataset_name}.uniprot_ids.csv"
params.uniprot_md5_csv = "${params.dataset_name}.uniprot_md5.csv"
params.gene3d_crh_output = "${params.dataset_name}.gene3d_crh_output.csv"
params.af_domainlist_ids_csv = "${params.dataset_name}.af_domainlist_ids.csv"
params.af_chainlist_ids_csv = "${params.dataset_name}.af_chainlist_ids.csv"
params.af_cath_orig_annotations_csv = "${params.dataset_name}.af_cath_orig_annotations.csv"


process create_all_af_chain_fasta {
    output:
        file all_af2_chain_fasta

    """
    gsutil cp ${params.all_af2_chain_fasta_url} ${params.all_af2_chain_fasta_file_name}
    """
}

process create_domain_s95_pdb_dir {
    input:
        file cath_domainlist_s95
        path cath_pdb_domain_dir
    output:
        path "${params.cath_s95_pdb_dir_name}/*"
    
    """
    rsync --files-from=$cath_domainlist_s95 $cath_pdb_domain_dir/ ${params.cath_s95_pdb_dir_name}/
    """
}

process create_foldseek_s95_library {
    input:
        path ${params.cath_s95_pdb_dir_name}
    
    output:
        file foldseek_s95_library

    """
    foldseek createdb ${params.cath_s95_pdb_dir_name} ${params.cath_s95_foldseek_library_name}
    """
}

process create_dataset_uniprot_ids {
    publishDir "${params.dataset_name}", mode: 'symlink'

    output:
        path "${params.uniprot_ids_csv}", emit: uniprot_ids_csv

    // export LD_LIBRARY_PATH=/Users/ian/Downloads/instantclient_19_8
    """
    cath-af-cli create-dataset-uniprot-ids \
        --uniprot_ids_csv ${params.uniprot_ids_csv} \
        --dbname ${params.cath_odb_name} \
        --max_evalue ${params.dataset_max_evalue} \
        --max_records ${params.dataset_max_records}
    """
}

process create_dataset_cath_files {
    publishDir "${params.dataset_name}", mode: 'symlink'

    input:
        path uniprot_ids_csv

    output:
        path params.uniprot_ids_csv, emit: uniprot_md5_csv
        path params.uniprot_ids_csv, emit: gene3d_crh_output
        path params.af_domainlist_ids_csv, emit: af_domainlist_ids_csv
        path params.af_chainlist_ids_csv, emit: af_chainlist_ids_csv
        path params.af_cath_orig_annotations_csv, emit: af_cath_orig_annotations_csv

    """
    cath-af-cli create-dataset-cath-files \
        --csv_uniprot_ids ${params.uniprot_ids_csv} \
        --csv_uniprot_md5 ${params.uniprot_md5_csv} \
        --gene3d_crh_output ${params.gene3d_crh_output} \
        --af_domainlist_ids ${params.af_domainlist_ids_csv} \
        --af_chainlist_ids ${params.af_chainlist_ids_csv} \
        --af_cath_annotations ${params.af_cath_orig_annotations_csv} \
        --dbname ${params.cath_odb_name} \

    """
}

// process create_dataset_af_files {}
// process create_annotation_chain_disorder {}
// process filter_domainlist_by_chain_disorder {}
// process optimise_domain_boundaries {}
// process create_annotation_domain_disorder {}
// process filter_domainlist_by_disorder {}
// process filter_domainlist_by_af2_quality {}
// process filter_domainlist_by_af2_packing {}
// process filter_domainlist_by_sse {}
// process chop_af_domains {}
// process create_annotation_foldseek_s95 {}
// process create_results_table {}


workflow {

    create_dataset_uniprot_ids | create_dataset_cath_files
    
}
