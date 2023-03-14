#!/usr/bin/env nextflow
nextflow.enable.dsl=2





// ********** PROCESSES ********** //
process create_cath_dataset_from_files {
    publishDir params.publish_dir, mode: 'copy'

    input:
    path 'uniprot_ids_csv'
    path 'all_crh_csv'
    path 'all_af_uniprot_md5_csv'

    output:
    path params.uniprot_md5_csv_fn
    path params.gene3d_crh_output_fn
    path params.af_domainlist_ids_csv_fn
    path params.af_chainlist_ids_csv_fn
    path params.af_cath_orig_annotations_csv_fn

    """
    grep -F -f uniprot_ids_csv all_crh_csv > filtered_crh_csv
    grep -F -f uniprot_ids_csv all_af_uniprot_md5_csv > filtered_af_uniprot_md5_csv
    cath-af-cli create-cath-dataset-from-files \
        --src_af_uniprot_md5 filtered_crh_csv \
        --src_crh filtered_af_crh \
        --csv_uniprot_ids uniprot_ids_csv \
        --csv_uniprot_md5 ${params.uniprot_md5_csv_fn} \
        --gene3d_crh_output ${params.gene3d_crh_output_fn} \
        --af_domainlist_ids ${params.af_domainlist_ids_csv_fn} \
        --af_chainlist_ids ${params.af_chainlist_ids_csv_fn} \
        --af_cath_annotations ${params.af_cath_orig_annotations_csv_fn}
    """
}

