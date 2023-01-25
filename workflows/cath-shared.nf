#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.chunk_size = 10

params.cath_af_cli = 'cath-af-cli'
params.cath_version = 'v4_3_0'

params.publish_dir = "$workflow.launchDir/cath-test-hpc-${params.dataset_name}"
params.af_version = 4
params.af_download_stem = "gs://public-datasets-deepmind-alphafold-v${params.af_version}"
params.uniprot_ids_csv_fn = "${params.dataset_name}.uniprot_ids.csv"
params.af_manifest_fn = "af_manifest.txt"
params.af_cif_raw_dir = "cif_raw"
params.af_cif_chopped_dir = "cif_chopped"

process create_af_manifest_file {
    input:
    path uniprot_id_file
    
    output:
    path params.af_manifest_fn

    """
    cat ${uniprot_id_file} | awk '{print "${params.af_download_stem}/AF-"\$1"-F*.cif"}' > ${params.af_manifest_fn}
    """
}

process retrieve_af_chain_cif_files {
    publishDir "${params.publish_dir}/${params.af_cif_raw_dir}", mode: 'copy'

    input:
    path af_model_urls_file

    output:
    path "AF-*.cif", optional: true

    // If Google returns 401 errors then make sure you have logged in:
    // 
    // gcloud auth application-default login
    //
    // see: https://www.nextflow.io/docs/latest/google.html

    """
    cat ${af_model_urls_file} | (gsutil -o GSUtil:parallel_process_count=1 -m cp -I . || echo "Ignoring non-zero exit code: \$?")
    """
}

process chop_cif {
    publishDir params.publish_dir, mode: 'copy'

    input:
    path raw_cif_files
    path full_uniprot_domain_ids_file
 
    output:
    path "${params.af_cif_chopped_dir}/*.cif"

    """
    mkdir ${params.af_cif_chopped_dir}

    # get a list of uniprot accessions from CIF files
    find . -name "*.cif" | tr '-' ' ' | awk '{print \$2}' | sort | uniq > uniprot_ids.txt

    # extract the uniprot domain ids for which we have CIF files
    grep -F -f uniprot_ids.txt ${full_uniprot_domain_ids_file} > uniprot_domain_ids.txt

    # do the chopping
    cath-af-cli chop-cif \
        --cif_in_dir . \
        --id_file uniprot_domain_ids.txt \
        --cif_out_dir ${params.af_cif_chopped_dir}/ \
        --id_type uniprot \
        --input_file_policy skip \
        --output_file_policy skip \
        --af_version ${params.af_version}
    """
}

process uniprot_domain_to_uniprot {
    input:
    path uniprot_domain_id

    output:
    path params.uniprot_ids_csv_fn

    """
    cat ${uniprot_domain_id} | tr '/' ' ' | awk '{print \$1}' | sort | uniq  > ${params.uniprot_ids_csv_fn}
    """
}

process cif_paths_to_uniprot_accessions {
    input: 
    path 'cif_paths.txt'

    output:
    path 'uniprot_ids.txt'

    """
    cat cif_paths.txt | perl -ne '/AF-(\\w+)-/ && print "\$1\\n"' | sort | uniq > uniprot_ids.txt
    """
}

process create_missing_uniprot_domain_ids {
    input:
    path 'found_uniprot_ids.txt'
    path 'all_uniprot_domain_ids.txt'

    output:
    path 'missing_uniprot_domain_ids.txt'

    """
    grep -F -v -f found_uniprot_ids.txt all_uniprot_domain_ids.txt > missing_uniprot_domain_ids.txt
    """
}