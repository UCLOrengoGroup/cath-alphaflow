#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.chunk_size = 10

params.cath_af_cli = 'cath-af-cli'
params.cath_version = 'v4_3_0'

params.publish_dir = "$workflow.launchDir/results-${params.dataset_name}"
params.af_version = 4
params.af_download_stem = "gs://public-datasets-deepmind-alphafold-v${params.af_version}"
params.af_cif_raw_dir = "cif_raw"
params.af_cif_chopped_dir = "cif_chopped"

// Processes to include from the shared module
include { uniprot_domain_to_uniprot } from './cath-shared-core'
include { create_af_manifest_file } from './cath-shared-core'
include { retrieve_af_chain_cif_files } from './cath-shared-core'
include { chop_cif } from './cath-shared-core'
include { cif_paths_to_uniprot_accessions } from './cath-shared-core'
include { create_missing_uniprot_domain_ids } from './cath-shared-core'

///////// SUBMODULE WORKFLOWS /////////////////////////////////
workflow AF_CIF_FILES {

    take:
        uniprot_domain_ids_ch        
            
    main:
        def uniprot_ids_ch = uniprot_domain_ids_ch
        | uniprot_domain_to_uniprot

        // split this into smaller chunks, retrieve files (download from AF), collect file paths
        uniprot_ids_ch.splitText(by: 10, file: true)
        | create_af_manifest_file
        | retrieve_af_chain_cif_files
    
    emit:
        retrieve_af_chain_cif_files.out

}

workflow CHOP_CIF {

    take:
        cif_files
        uniprot_domain_ids_path
            
    main:
        // actually chop the CIF files according to the CATH domain definitions 
        chop_cif( cif_files, uniprot_domain_ids_path )      
                          
}

workflow OUTPUT_UNIPROT {
    
    take:
        cif_files
        uniprot_domain_ids_ch
    
    main:        
        // find out the uniprot ids that have been successfully downloaded
        def downloaded_uniprot_ids_ch = cif_files
            | cif_paths_to_uniprot_accessions

        // create a list of all the uniprot domain ids that have not been downloaded 
        def undownloaded_uniprot_domain_ids_ch = create_missing_uniprot_domain_ids(
            downloaded_uniprot_ids_ch,
            uniprot_domain_ids_ch
        )

        
}

