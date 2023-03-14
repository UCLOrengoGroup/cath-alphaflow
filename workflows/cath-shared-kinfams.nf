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