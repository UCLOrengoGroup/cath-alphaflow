#!/usr/bin/env nextflow
nextflow.enable.dsl=2





// ********** PROCESSES ********** //
process af_domain_ids_from_cif_dir {
    input:
    path 'cif_dir'

    output:
    path 'af_models.txt'

    """
    ls cif_dir/ | grep ".cif" | sed 's/.cif\$//' > af_models.txt
    """
}

process ids_from_cif_dir {
    input:
    path 'cif_dir'

    output:
    path 'ids.txt'

    """
    ls cif_dir/ | grep ".cif" | sed 's/.cif\$//' > ids.txt
    """
}

process uniprot_csv_from_af_domains {
    input:
    path af_domain_id_file

    output:
    path 'uniprot_ids.csv'

    """
    echo uniprot_acc > uniprot_ids.csv
    cat ${af_domain_id_file} | tr '-' ' ' | awk '{print \$2}' >> uniprot_ids.csv
    """
}

process create_dataset_cath_files {
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
    cath-af-cli create-dataset-cath-files \
        --csv_uniprot_ids uniprot_ids_csv \
        --csv_uniprot_md5 ${params.uniprot_md5_csv_fn} \
        --gene3d_crh_output ${params.gene3d_crh_output_fn} \
        --af_domainlist_ids ${params.af_domainlist_ids_csv_fn} \
        --af_chainlist_ids ${params.af_chainlist_ids_csv_fn} \
        --af_cath_annotations ${params.af_cath_orig_annotations_csv_fn} \
        --dbname ${params.cath_odb_name} \
    """
}

process create_dssp {
    publishDir params.publish_dir, mode: 'copy'

    input:
    path 'uniprot_ids_csv'
    path 'cif_dir'
    path dssp_dir

    output:
    path "${dssp_dir}/*.dssp"

    """
    mkdir dssp
    cath-af-cli convert-cif-to-dssp \
        --cif_in_dir cif_dir \
        --id_file uniprot_ids_csv \
        --cif_suffix .cif \
        --dssp_out_dir ${dssp_dir}
    """
}

process create_plddt_summary {
    input:
    path 'uniprot_ids_csv'
    path 'cif_chopped'

    output:
    path params.plddt_stats_fn

    """
    cath-af-cli convert-cif-to-plddt-summary \
        --cif_in_dir cif_chopped \
        --id_file uniprot_ids_csv \
        --plddt_stats_file ${params.plddt_stats_fn}
    """
}

process create_sse_summary {
    publishDir params.publish_dir, mode: 'copy'

    input:
    path 'uniprot_ids_csv'
    path 'dssp_dir'

    output:
    path 'dssp_summary.csv'

    """
    cath-af-cli convert-dssp-to-sse-summary \
        --dssp_dir dssp_dir \
        --id_file uniprot_ids_csv \
        --sse_out_file dssp_summary.csv
    """
}