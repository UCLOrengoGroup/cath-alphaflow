#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Params defined before include are also used in the processes
params.dataset_name = 'dataset100'

// Processes to include from the shared module
include { uniprot_domain_to_uniprot } from './cath-shared-core'
include { create_af_manifest_file } from './cath-shared-core'
include { retrieve_af_chain_cif_files } from './cath-shared-core'
include { chop_cif } from './cath-shared-core'
include { cif_paths_to_uniprot_accessions } from './cath-shared-core'
include { create_missing_uniprot_domain_ids } from './cath-shared-core'

// Params defined after include are used locally
params.dataset_dir = "$workflow.launchDir/tests/fixtures/${params.dataset_name}"
params.uniprot_domain_ids_csv_fn = "${params.dataset_name}.uniprot_domain_ids.csv"

workflow {
    def uniprot_domain_ids_path = "${params.dataset_dir}/${params.uniprot_domain_ids_csv_fn}"

    def all_uniprot_domain_ids_ch = Channel.fromPath(uniprot_domain_ids_path, checkIfExists: true)

    // reduce the list of uniprot domains to a list of unique uniprot accessions 
    def all_uniprot_ids_ch = all_uniprot_domain_ids_ch
        | uniprot_domain_to_uniprot

    // split this into smaller chunks, retrieve files (download from AF), collect file paths
    def all_cif_files = all_uniprot_ids_ch.splitText(by: 10, file: true)
        | create_af_manifest_file
        | retrieve_af_chain_cif_files

    // actually chop the CIF files according to the CATH domain definitions 
    chop_cif( all_cif_files, uniprot_domain_ids_path )

    // find out the uniprot ids that have been successfully downloaded
    def downloaded_uniprot_ids_ch = all_cif_files
        | cif_paths_to_uniprot_accessions

    // create a list of all the uniprot domain ids that have not been downloaded 
    def undownloaded_uniprot_domain_ids_ch = create_missing_uniprot_domain_ids(
        downloaded_uniprot_ids_ch,
        all_uniprot_domain_ids_ch
    )
}
