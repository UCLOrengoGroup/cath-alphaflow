#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Params defined before include are also used in the processes
params.dataset_name = 'dataset10'//'dataset200000'
params.dataset_dir = "$workflow.launchDir/tests/fixtures/${params.dataset_name}/${params.dataset_name}"
params.publish_dir = "$workflow.launchDir/results-${params.dataset_name}"


params.uniprot_md5_csv_fn = "uniprot_md5.csv"
params.uniprot_ids_csv_fn = "uniprot_ids.csv"
params.gene3d_crh_output_fn = "gene3d_crh_output.csv"
params.af_domainlist_ids_csv_fn = "af_domainlist_ids.csv"
params.af_chainlist_ids_csv_fn = "af_chainlist_ids.csv"
params.af_cath_orig_annotations_csv_fn = "af_cath_orig_annotations.csv"

params.src_decorated_crh = "${params.dataset_dir}.decorated_crh.csv"
params.src_af_uniprot_md5 = "${params.dataset_dir}.af_uniprot_md5.csv"
params.csv_uniprot_md5 = "${params.dataset_dir}.uniprot_ids.csv"
params.gene3d_crh_output = "${params.dataset_dir}.gene_3d_crh.csv"
params.af_domainlist_ids = "${params.dataset_dir}.csv"
params.af_chainlist_ids = "${params.dataset_dir}.csv"



// Processes to include from the shared module
include { create_cath_dataset_from_files } from './cath-shared-files'

// Submodule workflows to include
include { GATHER_CIF } from './cath-module-core'
include { AF_ANNOTATE_DOMAINS_CIF } from './cath-module-kinfams'

// Params defined after include are used locally
params.uniprot_domain_ids_csv_fn = "uniprot_domain_ids.csv"
params.all_crh_csv = "${params.dataset_dir}.af_uniprot_md5.csv"
params.all_af_uniprot_md5 = "${params.dataset_dir}.uniprot_ids.csv"




workflow {
    
    //The first 2 are commented out due to missing domain file in dataset10 - dataset10.uniprot_domain_ids.csv 
    
    def uniprot_domain_ids_path = "${params.dataset_dir}.${params.uniprot_domain_ids_csv_fn}"

    GATHER_CIF(uniprot_domain_ids_path)

    def uniprot_csv_ch = AF_ANNOTATE_DOMAINS_CIF()

    def uniprot_dataset = create_cath_dataset_from_files(uniprot_csv_ch, params.all_crh_csv, params.all_af_uniprot_md5)            
}
