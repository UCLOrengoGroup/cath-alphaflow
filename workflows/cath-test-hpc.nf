#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.chunk_size = 10

params.cath_af_cli = 'cath-af-cli'
params.cath_version = 'v4_3_0'
params.dataset_name = 'dataset100'

params.dataset_dir = "$workflow.launchDir/tests/fixtures/dataset100"
params.publish_dir = "$workflow.launchDir/cath-test-hpc-${params.dataset_name}"
params.af_version = 4
params.af_download_stem = "gs://public-datasets-deepmind-alphafold-v${params.af_version}"

params.uniprot_ids_csv_fn = "${params.dataset_name}.uniprot_ids.csv"
params.uniprot_domain_ids_csv_fn = "${params.dataset_name}.uniprot_domain_ids.csv"

params.af_manifest_fn = "af_manifest.txt"
params.af_cif_raw_dir = "cif_raw"
params.af_cif_chopped_dir = "cif_chopped"


process create_af_manifest_file {
    publishDir params.publish_dir, mode: 'copy'

    input:
    path uniprot_id_file
    
    output:
    path params.af_manifest_fn

    """
    cat ${uniprot_id_file} | awk '{print "${params.af_download_stem}/AF-"\$1"-F*.cif"}' > ${params.af_manifest_fn}
    """
}

process retrieve_af_chain_cif_files {
    publishDir params.af_cif_raw_dir, mode: 'copy'

    input:
    path af_model_urls_file

    output:
    path "AF-*.cif", optional: true

    // If Google returns 401 errors then make sure you have logged in:
    //   gcloud auth application-default login
    //
    // see: https://www.nextflow.io/docs/latest/google.html

    """
    cat ${af_model_urls_file} | (gsutil -m cp -I . || echo "Ignoring non-zero exit code: \$?")
    """
}

process chop_cif {
    publishDir params.publish_dir, mode: 'copy'

    input:
    path uniprot_domain_ids_file
    path raw_cif_files

    output:
    path 'cif_chopped/*.cif'

    """
    mkdir cif_chopped
    cath-af-cli chop-cif \
        --cif_in_dir . \
        --id_file ${uniprot_domain_ids_file} \
        --cif_out_dir cif_chopped/ \
        --id_type uniprot \
        --input_file_policy skip \
        --output_file_policy skip \
        --af_version ${params.af_version}
    """
}

process uniprot_domain_to_uniprot {
    publishDir params.publish_dir, mode: 'copy'

    input:
    path uniprot_domain_id

    output:
    path params.uniprot_ids_csv_fn

    """
    cat ${uniprot_domain_id} | tr '/' ' ' | awk '{print \$1}' > ${params.uniprot_ids_csv_fn}
    """
}

workflow {
    def uniprot_domain_ids_ch = Channel.fromPath("${params.dataset_dir}/${params.uniprot_domain_ids_csv_fn}", checkIfExists: true)

    def cif_files = uniprot_domain_ids_ch.splitText(by: 10, file: true)
        | uniprot_domain_to_uniprot
        | create_af_manifest_file
        | retrieve_af_chain_cif_files

    cif_files.collect()

    // not all of the uniprot accessions will have AF models
    // of the AF models we do have, we would expect to be able to chop all of them without error
    def existing_uniprot_ids = cif_files.flatten().map { n -> n.name.split('-')[1] }.collect().unique()

    def existing_uniprot_domain_ids_ch = uniprot_domain_ids_ch.splitText(by: 100, file: true).filter {
        n -> existing_uniprot_ids.contains(n.name.split('/')[0])
    }

    chop_cif( existing_uniprot_domain_ids_ch.splitText(by: 10, file: true), cif_files )
}
