#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Params in modules
//params.af_cif_chopped_dir = "cif_chopped"
params.publish_dir = "$workflow.launchDir/results-${params.dataset_name}"
//params.af_cif_raw_dir = "cif_raw"
//params.af_dssp_raw_dir = "dssp_raw"
//params.plddt_stats_fn = "plddt_summary.csv"
params.uniprot_md5_csv_fn = "uniprot_md5.csv"
params.gene3d_crh_output_fn = "gene3d_crh_output.csv"
params.af_domainlist_ids_csv_fn = "af_domainlist_ids.csv"
params.af_chainlist_ids_csv_fn = "af_chainlist_ids.csv"
params.af_cath_orig_annotations_csv_fn = "af_cath_orig_annotations.csv"
params.cath_odb_name = 'GENE3D_21'





// ********** PROCESSES ********** //
process create_cath_dataset_from_db {
    publishDir params.publish_dir, mode: 'copy'

    input:
    path 'uniprot_ids_csv'

    output:
    path params.uniprot_md5_csv_fn
    path params.gene3d_crh_output_fn
    path params.af_domainlist_ids_csv_fn
    path params.af_chainlist_ids_csv_fn
    path params.af_cath_orig_annotations_csv_fn

    """
    cath-af-cli create-cath-dataset-from-db \
        --csv_uniprot_ids uniprot_ids_csv \
        --csv_uniprot_md5 ${params.uniprot_md5_csv_fn} \
        --gene3d_crh_output ${params.gene3d_crh_output_fn} \
        --af_domainlist_ids ${params.af_domainlist_ids_csv_fn} \
        --af_chainlist_ids ${params.af_chainlist_ids_csv_fn} \
        --af_cath_annotations ${params.af_cath_orig_annotations_csv_fn} \
        --dbname ${params.cath_odb_name} \
    """
}

