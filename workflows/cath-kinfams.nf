#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.gs_bucket = 'gs://public-datasets-deepmind-alphafold'
params.chunk_size = 1000

params.cath_af_cli = 'cath-af-cli'
params.cath_version = 'v4_3_0'
params.cath_data_root = "/cath/data/${params.cath_version}"
params.all_af_chain_fasta_url = "${params.gs_bucket}/sequences.fasta"
params.all_af_chain_fasta_file_name = 'all_af2_chain_sequences.fasta'
params.cath_s95_pdb_dir_name = 'cath_s95_pdb'
params.cath_s95_foldseek_library_name = 'foldseek_s95_lib'
params.cath_odb_name = 'GENE3D_21'
params.dataset_name = 'dataset100k'
params.dataset_max_records = '100000'
params.dataset_max_evalue = '1E-50'

params.publish_dir = "$workflow.launchDir/kinfams-results"

params.af_version = 4
params.af_download_stem = "gs://public-datasets-deepmind-alphafold-v${params.af_version}"

params.af_to_cluster_fn = "human_kinfams_uniprot_list_full"

params.uniprot_ids_csv_fn = "uniprot_ids.csv"
params.uniprot_md5_csv_fn = "uniprot_md5.csv"
params.gene3d_crh_output_fn = "gene3d_crh_output.csv"
params.af_domainlist_ids_csv_fn = "af_domainlist_ids.csv"
params.af_chainlist_ids_csv_fn = "af_chainlist_ids.csv"
params.af_cath_orig_annotations_csv_fn = "af_cath_orig_annotations.csv"
params.af_manifest_file_fn = "af_manifest_file"
params.uniprot_domain_ids_fn = "uniprot_domain_id.txt"
params.af_manifest_fn = "af_manifest.txt"
params.af_cif_raw_dir = "cif_raw"
params.af_cif_chopped_dir = "cif_chopped"
params.af_dssp_raw_dir = "dssp_raw"
params.plddt_stats_fn = "plddt_summary.csv"


// A0A0S2Z4D1/43-337 kinases_4.3-FF-000306.faa
// Q15831/43-337 kinases_4.3-FF-000306.faa

process create_af_manifest_file {
    publishDir params.publish_dir, mode: 'copy', overwrite: true

    input:
    path uniprot_id_file
    
    output:
    path params.af_manifest_fn

    """
    cat ${uniprot_id_file} | awk '{print "${params.af_download_stem}/AF-"\$1"-F*.cif"}' > ${params.af_manifest_fn}
    """
}

process retrieve_af_chain_cif_files {
    publishDir params.af_cif_raw_dir, mode: 'copy', overwrite: true

    input:
    path af_model_urls_file

    output:
    path "AF-*.cif"

    // If Google returns 401 errors then make sure you have logged in:
    //   gcloud auth application-default login
    //
    // see: https://www.nextflow.io/docs/latest/google.html

    """
    cat ${af_model_urls_file} | (gsutil -m cp -I . || echo "Ignoring non-zero exit code: \$?")
    """
}

process kinfam_clusters_to_uniprot_domain_ids {
    publishDir params.publish_dir, mode: 'copy', overwrite: true

    input:
    path cluster_path

    output:
    path params.uniprot_domain_ids_fn

    """
    cat ${cluster_path} | awk '{print \$1}' > ${params.uniprot_domain_ids_fn}
    """
}

process uniprot_domain_to_uniprot {
    publishDir params.publish_dir, mode: 'copy', overwrite: true

    input:
    path uniprot_domain_id

    output:
    path params.uniprot_ids_fn

    """
    cat ${uniprot_domain_id} | tr '/' ' ' | awk '{print \$1}' > ${params.uniprot_ids_fn}
    """
}

process chop_cif {
    publishDir params.publish_dir, mode: 'copy', overwrite: true

    input:
    path uniprot_ids_file
    path raw_cif_files

    output:
    path 'cif_chopped/*.cif'

    """
    mkdir cif_chopped
    cath-af-cli chop-cif \
        --cif_in_dir ${params.af_cif_raw_dir} \
        --id_file ${uniprot_ids_file} \
        --cif_out_dir cif_chopped/ \
        --id_type uniprot \
        --input_file_policy skip \
        --output_file_policy skip \
        --af_version ${params.af_version}
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

process ids_from_cif_dir {
    input:
    path 'cif_dir'

    output:
    path 'ids.txt'

    """
    ls cif_dir/ | grep ".cif" | sed 's/.cif\$//' > ids.txt
    """
}

process af_domain_ids_from_cif_dir {
    input:
    path 'cif_dir'

    output:
    path 'af_models.txt'

    """
    ls cif_dir/ | grep ".cif" | sed 's/.cif\$//' > af_models.txt
    """
}

process cif_to_pdb {
    input:
    path 'tmp.cif'

    output:
    path 'tmp.pdb'

    """
    pdb_fromcif tmp.cif > tmp.pdb
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

workflow AF_CHOP_CIF {

    def kinfams_ch = Channel.fromPath(params.af_to_cluster_file, checkIfExists: true)

    def uniprot_domain_ids_ch = kinfam_clusters_to_uniprot_domain_ids(kinfams_ch)

    def cif_files = uniprot_domain_ids_ch.splitText(by: 100, file: true)
        | uniprot_domain_to_uniprot
        | create_af_manifest_file
        | retrieve_af_chain_cif_files

    cif_files.collect()

    chop_cif( uniprot_domain_ids_ch.splitText(by: 100, file: true), cif_files )
}

workflow AF_ANNOTATE_DOMAINS_CIF {

    def cif_chopped_dir = file("${params.publish_dir}/${params.af_cif_chopped_dir}")
    def cif_raw_dir = file("${params.publish_dir}/${params.af_cif_raw_dir}")
    def dssp_raw_dir = file("${params.publish_dir}/${params.af_dssp_raw_dir}")
    def plddt_file = file("${params.publish_dir}/${params.plddt_stats_fn}")

    def af_domain_ids_ch = af_domain_ids_from_cif_dir(cif_chopped_dir)
    def cif_ids_ch = ids_from_cif_dir(cif_raw_dir).splitText(by: 100, file: true)

    def uniprot_csv_ch = uniprot_csv_from_af_domains(af_domain_ids_ch)

    def uniprot_dataset = create_dataset_cath_files(uniprot_csv_ch)

    def dssp_files_ch = create_dssp(cif_ids_ch, cif_raw_dir, dssp_raw_dir)

    def af_domain_ids_chunked_ch = af_domain_ids_ch.splitText(by: 100, file: true)

    def plddt_summary_ch = create_plddt_summary(af_domain_ids_chunked_ch, cif_raw_dir)
    
    plddt_summary_ch.collectFile(name: plddt_file)

    dssp_files_ch.collect()

    create_sse_summary(af_domain_ids_ch, dssp_raw_dir)

}

workflow {
    // AF_CHOP_CIF ()
    AF_ANNOTATE_DOMAINS_CIF ()
}