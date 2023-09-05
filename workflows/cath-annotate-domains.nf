#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// the number of uniprot ids processed in each chunk of work 
params.chunk_size = 3

params.cath_af_cli = 'cath-af-cli'
params.cath_version = 'v4_3_0'

params.af_version = 4

// params.af_archive_root = "/cluster/ref3/alphafolddb/proteomes"
params.proteome_zip_file = "${workflow.launchDir}/tests/fixtures/zipfiles/proteome-tax_id-2025159-0_v4.zip"
params.uniprot_csv_file = "${workflow.launchDir}/tests/fixtures/dataset100/dataset100.uniprot_ids.csv"
params.tools_root = "${workflow.launchDir}/tools"
params.af_archive_root = "${workflow.launchDir}/datasets/chainsaw"
params.chainsaw_dir = "${params.tools_root}/chainsaw"
params.chainsaw_repo = "https://github.com/JudeWells/chainsaw"
params.merizo_dir = "${params.tools_root}/Merizo"
params.merizo_repo = "https://github.com/psipred/Merizo"
params.unidoc_dir = "${params.tools_root}/UniDoc"
params.unidoc_script = "${params.unidoc_dir}/Run_UniDoc_from_scratch_structure.py"
params.alphafold_url_stem = "https://alphafold.ebi.ac.uk/files"

process cif_files_from_web {
    input:
    path 'af_ids.txt'

    output:
    path 'AF-*.cif'

    """
    awk '{print "${params.alphafold_url_stem}/"\$1".cif"}' af_ids.txt > af_model_urls.txt
    wget -i af_model_urls.txt || true
    """
}

process cif_files_from_gs {
    input:
    path 'af_ids.txt'

    output:
    path "AF-*.cif", optional: true

    // If Google returns 401 errors then make sure you have logged in:
    // 
    // gcloud auth application-default login
    // 
    // see: https://www.nextflow.io/docs/latest/google.html

    """
    awk '{print \$1".cif"}' uniprot_ids.txt > af_ids.txt
    cat af_model_urls.txt | (gsutil -o GSUtil:parallel_process_count=1 -m cp -I . || echo "Ignoring non-zero exit code: \$?")
    """
}

process pdb_files_from_zip {
    input:
    path 'proteome.zip'

    output:
    path '*.pdb'

    """
    unzip proteome.zip
    """
}

process cif_to_pdb {
    input:
    path '*'

    output:
    path '*.pdb'

    """
    for cif_file in *.cif; do
        pdb_file=\${cif_file%.cif}.pdb
        pdb_fromcif \$cif_file > \$pdb_file
    done
    """
}

process install_chainsaw {
    // do this in a docker container? 

    """
    mkdir ${params.tools_root} && cd ${params.tools_root}
    git clone ${params.chainsaw_repo}
    cd chainsaw
    python3 -m venv venv
    . venv/bin/activate
    pip install --upgrade pip wheel
    pip install -r requirements.txt
    cd stride && tar -zxvf stride.tgz && make
    """
}

process install_merizo {

    """
    mkdir ${params.tools_root} && cd ${params.tools_root}
    git clone ${params.chainsaw_repo}
    cd chainsaw
    python3 -m venv venv
    . venv/bin/activate
    pip install --upgrade pip wheel
    pip install -r requirements.txt
    cd stride && tar -zxvf stride.tgz && make
    """
}

process run_chainsaw {
    input:
    path '*'

    output:
    path 'chainsaw_results.csv'

    """
    ${params.chainsaw_dir}/venv/bin/python3 ${params.chainsaw_dir}/get_predictions.py \
        --structure_directory . -o chainsaw_results.csv
    """
}

process run_merizo {
    input:
    path '*'

    output:
    path 'merizo_results.csv'

    """
    ${params.merizo_dir}/venv/bin/python3 ${params.merizo_dir}/predict.py \
        -i *.pdb > merizo_results.csv
    """
}

process run_measure_globularity {
    input:
    path 'af_domain_list.csv'
    path '*.cif'

    output:
    path 'af_domain_globularity.csv'

    """
    cath-af-cli measure-globularity \
        --af_domain_list af_domain_list.csv \
        --af_chain_mmcif_dir . \
        --af_domain_globularity af_domain_globularity.csv \
    """
}

process run_unidoc {
    input:
    path '*'

    output:
    path 'unidoc_results.csv'

    """
    # UniDoc expects to be run from its own directory
    ln -s ${params.unidoc_dir}/bin bin

    for pdb_file in *.pdb; do
        file_stem=\${pdb_file%.pdb}

        python3 ${params.unidoc_script} -i \$pdb_file -c A -o \${file_stem}.unidoc

        # extract the chopping from the last line of the unidoc output (possibly blank)
        # and change chopping string to the expected format
        chopping=\$(tail -n 1 \${file_stem}.unidoc | tr '~' '-' | tr ',' '_' | tr '/' ',' | tr -d '\\n')

        echo "\$file_stem\t\$chopping" >> unidoc_results.csv
    done
    """
}

process merge_domain_results {
    input:
    file 'id_file.tsv'
    file 'chainsaw_results.tsv'
    file 'merizo_results.tsv'
    file 'unidoc_results.tsv'

    output:
    file 'all_results.tsv'

    """
    cath-af-cli merge-domain-results \
        -i id_file.tsv \
        -c chainsaw_results.tsv \
        -m merizo_results.tsv \
        -u unidoc_results.tsv \
        -o all_results.tsv
    """
}

workflow {
    def zip_files_ch = Channel.fromPath( params.proteome_zip_file )

    def pdb_files_ch = pdb_files_from_zip( zip_files_ch )

    def af_ids_ch = pdb_files_ch.map { it.collect { it.getBaseName() } }

    def id_file_ch = af_ids_ch.collectFile(name: 'af_ids.txt', newLine: true)

    def chainsaw_results_ch = run_chainsaw( pdb_files_ch )
    def merizo_results_ch = run_merizo( pdb_files_ch )
    def unidoc_results_ch = run_unidoc( pdb_files_ch )

    def combine_results_ch = merge_domain_results( id_file_ch, chainsaw_results_ch, merizo_results_ch, unidoc_results_ch )

    combine_results_ch.collectFile(name: 'all_results.tsv', storeDir: workflow.launchDir)
        .subscribe {
            println "Entries are saved to file: $it"
        }
}

// workflow {

//     // Create a channel from the uniprot csv file
//     def uniprot_ids_ch = Channel.fromPath( params.uniprot_csv_file )
//         // process the file as a CSV with a header line
//         .splitCsv(header: true)
//         // only process a few ids when debugging
//         .take( 5 )

//     uniprot_ids_ch | view

//     // Generate files containing chunks of AlphaFold ids
//     // NOTE: this will only retrieve the first fragment in the AF prediction (F1)
//     def af_ids = uniprot_ids_ch
//         // make sure we don't have duplicate uniprot ids
//         .unique()
//         // map uniprot id (CSV row) to AlphaFold id
//         .map { up_row -> "AF-${up_row.uniprot_id}-F1-model_v${params.af_version}" }
//         // collect all ids into a single file
//         .collectFile(name: 'all_af_ids.txt', newLine: true)
//         // split into chunks and save to files
//         .splitText(file: 'chunked_af_ids.txt', by: params.chunk_size)

//     // download cif files from afdb server
//     def cif_ch = cif_files_from_web( af_ids )

//     // convert cif to pdb files
//     def pdb_ch = cif_to_pdb( cif_ch )

//     // run chainsaw on the pdb files
//     def chainsaw_results_ch = run_chainsaw( pdb_ch )

//     chainsaw_results_ch.view { println "chainsaw_results: '${it}'" }

//     def chainsaw_results = chainsaw_results_ch
//         .collectFile(name: 'results.chainsaw.csv', 
//             // skip: 1,
//             storeDir: workflow.launchDir)
//         .subscribe {
//             println "Entries are saved to file: $it"
//         }
// }