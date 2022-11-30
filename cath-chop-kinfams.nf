#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.gs_bucket = 'gs://public-datasets-deepmind-alphafold'
params.chunk_size = 1000

params.cath_af_cli = 'cath-af-cli'
params.cath_version = 'v4_3_0'
params.cath_data_root = "/cath/data/${params.cath_version}"
params.all_af_chain_fasta_url = "${params.gs_bucket}/sequences.fasta"
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
params.af_manifest_file = "${params.dataset_name}.af_manifest_file"

params.publish_dir = "$workflow.projectDir/kinfams-results"
params.af_to_cluster_file = "${params.publish_dir}/human_kinfams_uniprot_list_full"
params.uniprot_ids_filename = "uniprot_id.txt"
params.uniprot_domain_ids_filename = "uniprot_domain_id.txt"
params.af_manifest_filename = "af_manifest.txt"
params.af_cif_raw_dir = "${params.publish_dir}/cif_raw"
params.af_cif_chopped_dir = "${params.publish_dir}/cif_chopped"
params.af_version = 4
params.af_download_stem = "gs://public-datasets-deepmind-alphafold-v${params.af_version}"


// A0A0S2Z4D1/43-337 kinases_4.3-FF-000306.faa
// Q15831/43-337 kinases_4.3-FF-000306.faa

process create_af_manifest_file {
    publishDir params.publish_dir, mode: 'copy', overwrite: true

    input:
    path uniprot_id_file
    
    output:
    path params.af_manifest_filename

    """
    cat ${uniprot_id_file} | awk '{print "${params.af_download_stem}/AF-"\$1"-F*.cif"}' > ${params.af_manifest_filename}
    """
}

process retrieve_af_chain_cif_files {
    publishDir params.af_cif_raw_dir, mode: 'copy', overwrite: true

    input:
    path af_model_urls_file

    output:
    path "AF-*.cif"

    // If Google returns 401 errors then make sure you have logged in:
    //   gcloud auth application-default login
    //
    // see: https://www.nextflow.io/docs/latest/google.html

    """
    cat ${af_model_urls_file} | (gsutil -m cp -I . || echo "Ignoring non-zero exit code: \$?")
    """
}

process kinfam_clusters_to_uniprot_domain_ids {
    publishDir params.publish_dir, mode: 'copy', overwrite: true

    input:
    path cluster_path

    output:
    path params.uniprot_domain_ids_filename

    """
    cat ${cluster_path} | awk '{print \$1}' > ${params.uniprot_domain_ids_filename}
    """

}

process uniprot_domain_to_uniprot {
    publishDir params.publish_dir, mode: 'copy', overwrite: true

    input:
    path uniprot_domain_id

    output:
    path params.uniprot_ids_filename

    """
    cat ${uniprot_domain_id} | tr '/' ' ' | awk '{print \$1}' > ${params.uniprot_ids_filename}
    """
   
}

process chop_cif {
    publishDir params.publish_dir, mode: 'copy', overwrite: true

    input:
    path uniprot_ids_file
    path raw_cif_files

    output:
    path 'cif_chopped/*.cif'

    """
    mkdir cif_chopped
    cath-af-cli chop-cif \
        --cif_in_dir ${params.af_cif_raw_dir} \
        --id_file ${uniprot_ids_file} \
        --cif_out_dir cif_chopped/ \
        --id_type uniprot \
        --input_file_policy skip \
        --output_file_policy skip \
        --af_version ${params.af_version}
    """
}


workflow AF_CHOP_CIF {

    def kinfams_ch = Channel.fromPath(params.af_to_cluster_file, checkIfExists: true)

    def uniprot_domain_ids_ch = kinfam_clusters_to_uniprot_domain_ids(kinfams_ch)

    def cif_files = uniprot_domain_ids_ch.splitText(by: 100, file: true)
        | uniprot_domain_to_uniprot
        | create_af_manifest_file
        | retrieve_af_chain_cif_files

    cif_files.collect()

    chop_cif( uniprot_domain_ids_ch.splitText(by: 100, file: true), cif_files )
}

workflow {
    AF_CHOP_CIF ()
}