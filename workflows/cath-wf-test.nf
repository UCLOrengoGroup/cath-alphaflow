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

// Submodule workflows to include
include { AF_CIF_FILES } from './cath-module-core'
include { CHOP_CIF } from './cath-module-core'
include { OUTPUT_UNIPROT } from './cath-module-core'

// Params defined after include are used locally
params.dataset_dir = "$workflow.launchDir/tests/fixtures/${params.dataset_name}"
params.uniprot_domain_ids_csv_fn = "${params.dataset_name}.uniprot_domain_ids.csv"

workflow {
    def uniprot_domain_ids_path = "${params.dataset_dir}/${params.uniprot_domain_ids_csv_fn}"

    def all_uniprot_domain_ids_ch = Channel.fromPath(uniprot_domain_ids_path, checkIfExists: true)
                    
    def all_cif_files = AF_CIF_FILES(all_uniprot_domain_ids_ch)
    
    CHOP_CIF(all_cif_files,uniprot_domain_ids_path)
    
    OUTPUT_UNIPROT(all_cif_files,all_uniprot_domain_ids_ch)
}
