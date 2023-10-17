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
params.uniprot_domain_ids_file = "${workflow.launchDir}/tests/fixtures/dataset100/dataset100.uniprot_domain_ids.csv"
// params.afid_chopping_file = "${workflow.launchDir}/tests/fixtures/dataset100/dataset100.afid_chopping.csv"
params.afid_chopping_file = "${workflow.launchDir}/zipfiles.00000001.results.head.csv"
params.tools_root = "${workflow.launchDir}/tools"
params.af_archive_root = "${workflow.launchDir}/datasets/chainsaw"
params.chainsaw_dir = "${params.tools_root}/chainsaw"
params.chainsaw_repo = "https://github.com/JudeWells/chainsaw"
params.merizo_dir = "${params.tools_root}/Merizo"
params.merizo_repo = "https://github.com/psipred/Merizo"
params.unidoc_dir = "${params.tools_root}/UniDoc"
params.unidoc_script = "${params.unidoc_dir}/Run_UniDoc_from_scratch_structure.py"
params.alphafold_url_stem = "https://alphafold.ebi.ac.uk/files"


include { chop_pdb as CHOP_CIF } from './cath-shared-core'


process uniprot_domain_id_to_af_id {
    input:
    path 'uniprot_domain_ids.csv'

    output:
    path 'af_ids.csv'

    """
    cat uniprot_domain_ids.csv | tr '/' ' ' | \
        perl -ne '(\$id, \$ch) = split(/\s+/, \$_); print("AF-\$id-F1-model_v4 \$ch")' > af_ids.csv
    """
}

process cif_files_from_web {
    input:
    path 'af_ids.txt'

    output:
    path 'AF-*.cif'

    """
    cat af_ids.txt | tr '/' ' ' | awk '{print "${params.alphafold_url_stem}/"\$1".cif"}' > af_model_urls.txt
    wget -i af_model_urls.txt || true
    """
}

process pdb_files_from_web {
    input:
    path 'af_ids.txt'

    output:
    path 'AF-*.pdb'

    """
    awk '{print "${params.alphafold_url_stem}/"\$1".pdb"}' af_ids.txt > af_model_urls.txt
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

process run_plddt_summary {
    input:
    path 'af_domain_list.csv'
    path '*'

    output:
    path 'plddt_summary.csv'

    """
    find -L . -name '*.cif' | sed 's/^\\.\\///' | sed 's/\\.cif//g' > cif_ids.txt
    echo "af_domain_id" > af_domain_list_files_present.csv
    grep -F -f cif_ids.txt af_domain_list.csv >> af_domain_list_files_present.csv
    cath-af-cli convert-cif-to-plddt-summary \
        --id_file af_domain_list_files_present.csv \
        --id_type af \
        --cif_in_dir . \
        --plddt_stats_file plddt_summary.csv
    
    """
}


process run_globularity_summary {
    input:
    path 'af_domain_list.csv'
    path '*'

    output:
    path 'globularity_summary.csv'

    """
    find -L . -name '*.cif' | sed 's/^\\.\\///' | sed 's/\\.cif//g' > cif_ids.txt
    echo "af_domain_id" > af_domain_list_files_present.csv
    grep -F -f cif_ids.txt af_domain_list.csv >> af_domain_list_files_present.csv
    cath-af-cli measure-globularity \
        --af_domain_list af_domain_list_files_present.csv \
        --af_chain_mmcif_dir . \
        --af_domain_globularity globularity_summary.csv
    """
}

process merge_results {
    input:
    file 'id_file.tsv'
    // file 'uniprot_md5_file.tsv'
    file 'globularity_results.tsv'
    file 'plddt_results.tsv'

    output:
    file 'combined_results.tsv'

    """
    cath-af-cli merge-results \
        --id_file id_file.tsv \
        --glob_file globularity_results.tsv \
        --plddt_file plddt_results.tsv \
        --results_file combined_results.tsv
    """
}


process af_domain_ids_from_chopping_file {
    input:
    file 'chopping_file.csv'

    output:
    file 'af_domain_ids.csv'

    // chain_id        sequence_md5    nres    ndom    chopping        uncertainty
    // AF-A0A1N7SYI2-F1-model_v4       f26aa0fe84a8270c2a2c3dfa811e40a3        258     2       37-117,124-219  0.0178

    """
#!/usr/bin/env python3

IS_ZERO_INDEXED = True

with open("af_domain_ids.csv", "wt") as fp_out:
    fp_out.write("af_domain_id\\n")
    with open("chopping_file.csv", "rt") as fp_in:
        next(fp_in) # skip header
        for line in fp_in:
            cols = line.split()
            chain_id = cols[0]
            chopping = cols[4]
            doms = chopping.split(',')
            for dom in doms:
                if dom == 'NULL':
                    continue
                segs = dom.split('_')
                if IS_ZERO_INDEXED:
                    segs_new = []
                    for seg in segs:
                        print(f"dom: {dom}, segs: {segs}, seg: {seg}")
                        start, end = seg.split('-')
                        segs_new.append(f"{int(start)+1}-{int(end)+1}")
                    segs = segs_new
                fp_out.write(f"{chain_id}/{'_'.join(segs)}\\n")
    """
}

workflow {

    // Create a channel from AF id file
    def afid_chopping_file = Channel.fromPath( params.afid_chopping_file )

    def af_domain_ids = af_domain_ids_from_chopping_file( afid_chopping_file )

    // // get AF ids from the chopping file
    def af_ids = afid_chopping_file
                    .splitCsv( header: true, sep: '\t')
                    .map { it.chain_id }
                    .collectFile( newLine: true )

    def cif_ch = cif_files_from_web( af_ids )

    // calculate globularity
    def globularity_results_file = run_globularity_summary( af_domain_ids, cif_ch ).collectFile()

    // calculate plddt
    def plddt_summary_file = run_plddt_summary( af_domain_ids, cif_ch ).collectFile()

    // merge results
    def combine_results_ch = merge_results( 
        af_domain_ids,
        globularity_results_file,
        plddt_summary_file,
    )
    
    combine_results_ch.collectFile(name: 'combined_results.tsv', storeDir: './')
}