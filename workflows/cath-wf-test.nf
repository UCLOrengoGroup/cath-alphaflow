#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Params defined before include are also used in the processes
params.dataset_name = 'dataset100'

// Processes to include from the shared module
include { create_cath_dataset_from_db } from './cath-shared-database'

// Submodule workflows to include
include { GATHER_CIF } from './cath-module-core'
include { AF_ANNOTATE_DOMAINS_CIF } from './cath-module-kinfams'

// Params defined after include are used locally
params.dataset_dir = "$workflow.launchDir/tests/fixtures/${params.dataset_name}"
params.uniprot_domain_ids_csv_fn = "${params.dataset_name}.uniprot_domain_ids.csv"

workflow {
    
    def uniprot_domain_ids_path = "${params.dataset_dir}/${params.uniprot_domain_ids_csv_fn}"

    GATHER_CIF(uniprot_domain_ids_path)
        
    def uniprot_csv_ch = AF_ANNOTATE_DOMAINS_CIF()            
}
