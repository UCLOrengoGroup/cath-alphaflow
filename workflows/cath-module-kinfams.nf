#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Params at this level
params.af_cif_chopped_dir = "cif_chopped"
//params.publish_dir = "$workflow.launchDir/kinfams-results"
params.publish_dir = "$workflow.launchDir/results-${params.dataset_name}"
params.af_cif_raw_dir = "cif_raw"
params.af_dssp_raw_dir = "dssp_raw"
params.plddt_stats_fn = "plddt_summary.csv"

// Params in modules
params.uniprot_md5_csv_fn = "uniprot_md5.csv"
params.gene3d_crh_output_fn = "gene3d_crh_output.csv"
params.af_domainlist_ids_csv_fn = "af_domainlist_ids.csv"
params.af_chainlist_ids_csv_fn = "af_chainlist_ids.csv"
params.af_cath_orig_annotations_csv_fn = "af_cath_orig_annotations.csv"
params.cath_odb_name = 'GENE3D_21'

// Processes to include from the shared module
include { af_domain_ids_from_cif_dir } from './cath-shared-kinfams'
include { ids_from_cif_dir } from './cath-shared-kinfams'
include { uniprot_csv_from_af_domains } from './cath-shared-kinfams'
include { create_dataset_cath_files } from './cath-shared-kinfams'
include { create_dssp } from './cath-shared-kinfams'
include { create_plddt_summary } from './cath-shared-kinfams'
include { create_sse_summary } from './cath-shared-kinfams'


// ********** WORKFLOWS ********** //



workflow AF_ANNOTATE_DOMAINS_CIF {

    def cif_chopped_dir = file("${params.publish_dir}/${params.af_cif_chopped_dir}")
    def cif_raw_dir = file("${params.publish_dir}/${params.af_cif_raw_dir}")
    def dssp_raw_dir = file("${params.publish_dir}/${params.af_dssp_raw_dir}")
    def plddt_file = file("${params.publish_dir}/${params.plddt_stats_fn}")

    def af_domain_ids_ch = af_domain_ids_from_cif_dir(cif_chopped_dir)

    def af_domain_ids_chunked_ch = af_domain_ids_ch.splitText(by: 100, file: true)

    def plddt_summary_ch = create_plddt_summary(af_domain_ids_chunked_ch, cif_raw_dir)
    
    plddt_summary_ch.collectFile(name: plddt_file)
    
    def cif_ids_ch = ids_from_cif_dir(cif_raw_dir).splitText(by: 100, file: true)
        
    def uniprot_csv_ch = uniprot_csv_from_af_domains(af_domain_ids_ch)

    // This gets files from the oracle database
    def uniprot_dataset = create_dataset_cath_files(uniprot_csv_ch)
                    
    //def dssp_files_ch = create_dssp(cif_ids_ch, cif_raw_dir, dssp_raw_dir)                
    //dssp_files_ch.collect()
    
    // what does this do? why do I get errors?
    //create_sse_summary(af_domain_ids_ch, dssp_raw_dir)

}

