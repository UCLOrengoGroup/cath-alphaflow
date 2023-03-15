Testing chopping of 100,000 domains from AlphaFold

1) Install repo
2) Create .env file (use cathora1)
3) export instantclient using export LD_LIBRARY_PATH=~/instantclient_21_7
4) load python module - module load python/3.8.5
5) Retrieve list of 100k UniProt IDs from database 

start: 14.31 - end: 14.38 - time: 7minutes

Input: cath-af-cli create-dataset-uniprot-ids --uniprot_ids_csv af_uniprot_id_100k_from_g3d.list --dbname gene3d_21 --max_evalue 1e-50 --max_records 100000
Output: af_uniprot_id_100k_from_g3d.list

6) Run create-cath-dataset-from-db
start 14:40 - end: 14:44 time: 4 minutes
Input: af_uniprot_id_100k_from_g3d.list

cath-af-cli create-cath-dataset-from-db \
    --csv_uniprot_ids af_uniprot_id_100k_from_g3d.list \
    --dbname gene3d_21 \
    --af_cath_annotations af_100k_cath_annotations.tsv \
    --af_chainlist_ids af_100k_chainlist_ids.tsv \
    --af_domainlist_ids af_100k_domainlist_ids.tsv \
    --gene3d_crh_output af_100k.crh \
    --csv_uniprot_md5 af_100k_md5.tsv

Output: 
    af_100k_cath_annotations.tsv 
    af_100k_chainlist_ids.tsv 
    af_100k_chainlist_ids.tsv 
    af_100k.crh 
    af_100k_md5.tsv

7) Create GCloud Manifest File - less than 1 min
Input: af_100k_chainlist_ids.tsv
Output: af_100k_manifest.txt

grep -v 'chain_id' af_100k_chainlist_ids.tsv | sort | uniq | tr -d '\r' | awk '{print "gs://public-datasets-deepmind-alphafold-v4/"$1".cif"}' > af_100k_manifest.txt


8) Download structures from GCloud

start: 15.14 end: 17.35 - time: 2hrs and 20mins

Input: af_100k_manifest.txt
Output: af_raw_cif folder containing unchopped AlphaFold chains

mkdir af_cif_raw
cat af_100k_manifest.txt | (/share/apps/google-cloud-sdk/bin/gsutil -o GSUtil:parallel_process_count=1 -m cp -I af_cif_raw/ || echo "Ignoring non-zero exit code: \$?")

CommandException: 6532 files/objects could not be transferred. Total downloaded CIF: 62,212

9) Create chainlist after gsutil. (i.e. needed now as cif_to_fasta expects a list of files that exist, so we need to feed it just the chains that are in AFDB and are available on the filesystem).

cd af_cif_raw
find *.cif > ../af_100k_chainlist_after_gsutil.txt
sed -i 's/\.cif//g' af_100k_chainlist_after_gsutil.txt
echo af_chain_id | cat - af_100k_chainlist_after_gsutil.txt > temp && mv temp af_100k_chainlist_after_gsutil.txt



10) Convert CIF to FASTA (62,212 CIF files -> 62,212 FASTA files + merged.fasta)
Input: af_cif_raw folder, af_100k_chainlist_after_gsutil.txt, af_fasta_raw folder 
Output: af_fasta_raw folder containing a FASTA file for each chain available in the filesystem

Thoughts: Too slow on 1CPU (i.e. 30mins for 10k sequences, splitting and parallelizing)


Start: 12:11  
End:   13:20
Time:  1hr 9 mins

mkdir af_fasta_raw
split -l 10000 af_100k_chainlist_after_gsutil.txt splitlist
echo `date` start; find splitlista* | xargs -P4 -I XXX sh -c 'cath-af-cli convert-cif-to-fasta --cif_in_dir af_cif_raw/ --id_file XXX --fasta_out_dir af_fasta_raw/' ; echo `date` end

Careful: the module expects each chunk to have a header, so the first structure doesn't get processed if we use chunks instead of the first file.

11) Create MD5 from FASTA (seconds)

Input: merged.fasta
Output: af-100k_cif_raw_md5.txt

cath-af-cli create-md5 --fasta merged.fasta --uniprot_md5_csv af_100k_cif_raw_md5.txt

12) Filter CRH (get only matches between CIF-md5 and chopping MD5, remove all CRH entries that rely on a missing AFDB entry)

Issue: The files generated in step 6 have some carriage returns (like ^M) that needs removing using tr -d '\r'
Temporary fix: strip these characters before doing the grep command.

Command: 

echo `date` start; grep -v 'af_chain_id' af_100k_cif_raw_md5.txt | awk '{print $3}' | tr -d '\r' | grep -F -f - af_100k.crh > af_100k_crh_after_filtering.crh; echo `date` end
sort af_100k_crh_after_filtering.crh | uniq > af_100k_crh_after_filtering_uniq.crh
awk '{print $1}' af_100k_crh_after_filtering_uniq.crh | grep -F -f - af_100k_cif_raw_md5.txt | awk '{print $1}' | grep -F -f - af_100k_domainlist_ids.txt > af_100k_domainlist_after_md5_filter.txt
echo af_domain_id | cat - af_100k_domainlist_after_md5_filter.txt > temp && mv temp af_100k_domainlist_after_md5_filter.txt

Time: 1 second
Output: af_100k_crh_after_filtering.crh

13) CHOP CIF

Running on xargs with 12 processes with a uniq-ed crh file and domain list (120,129 domains)

Start: 15:15
End:   17.01
Time:  1hr 45 minutes

Command: echo `date` start; find split_crh_uniq_000* | xargs -P12 -I XXX sh -c 'cath-af-cli chop-cif --cif_in_dir af_cif_raw/ --id_file XXX --cif_out_dir chopped_cif_before_optimization/'; echo `date` end

Comments: remove af_domain_id from header, program can't parse it or add exception, as it skips the first entry but it crashes as it can't find af_domain_id

Chopped 120,129 domains

14) Optimise boundaries
On 16 threads, 7500 entries each
Start: 15.07
End: 16.08
Time: 1hr 1 min

Comments: include percentage of delta length (i.e. how much the boundaries are changing?)
Remove af_domain_id from files, it crashes the process.
Optimise boundaries works on chains, not on domains. The initial chopping is redundant, keeping it in this instance as it's interesting for debugging.

Commands:
split -l 7500 af_100k_domainlist_after_md5_filter_uniq.txt split_af_100k -a 5 -d

echo `date` start; find split_af_100k000* | xargs -P16 -I XXX sh -c 'cath-af-cli optimise-domain-boundaries --af_domain_list XXX --af_chain_mmcif_dir af_cif_raw/ --af_domain_list_post_tailchop af_100k_domainlist_post_tailchop_XXX.txt --af_domain_mapping_post_tailchop af_100k_domain_mapping_post_tailchop_XXX.tsv --status_log optimise_boundaries_status_log_XXX'; echo `date` end

Afterwards, concatenate logs, af_domain_list_post_tailchop and af_domain_mapping_post_tailchop and remove duplicated headers.

cat af_100k_domain_mapping_post_tailchop_split_af_100k000* | grep -v 'af_domain_id' > af_100k_domain_mapping_post_tailchop.tsv
cat af_100k_domainlist_post_tailchop_split_af_100k000* grep -v 'af_domain_id' > af_100k_domainlist_post_tailchop.txt
cat optimise_boundaries_status_log_split_af_100k000* > optimise_boundaries_status_log
rm af_100k_domain_mapping_post_tailchop_split_af_100k000*
rm optimise_boundaries_status_log_split_af_100k000*
rm af_100k_domainlist_post_tailchop_split_af_100k000*
rm split_af_100k000*


Stats: 84,604 have unchanged boundaries
       35,378 have adjusted boundaries
       129 have failed due to not having residues above pLDDT 70



15) CHOP CIF (after optimisation)

Start: 18:19 
End: 20.06
Time: 1hr 46 mins

split -l 7500 af_100k_domainlist_post_tailchop.txt af_100k_domainlist -a 5 -d
echo `date` start; find af_100k_domainlist00* | xargs -P12 -I XXX sh -c 'cath-af-cli chop-cif --cif_in_dir af_cif_raw/ --id_file XXX --cif_out_dir chopped_cif_after_optimisation/'; echo `date` end

TODO: Add StatusLog to ALL MODULES
Chopped 120,112 domains with optimised boundaries

16) Calculate pLDDT and LUR report

Start: 13:42
End: 14.41
Time:  59 mins

split -l 7500 af_100k_domainlist_post_tailchop.txt split_af_100k_plddt_ -a 5 -d
echo `date` start; find split_af_100k_plddt_000* | xargs -P16 -I XXX sh -c 'cath-af-cli convert-cif-to-plddt-summary --cif_in_dir af_cif_raw/ --id_file XXX --plddt_stats_file af_100k_plddt_summary_XXX.tsv'; echo `date` end
cat af_100k_plddt_summary_split_af_100k_plddt_*.tsv > af_100k_plddt_summary.tsv

17) Create DSSP Files

Start: 18.08
End: 22.11
Time: 4 hrs

module load boost/1.71.0
mkdir af_dssp_dir
echo `date` start; find split_af_100k_dssp_000* | xargs -P16 -I XXX sh -c 'cath-af-cli convert-cif-to-dssp --cif_in_dir af_cif_raw/ --id_file XXX --dssp_out_dir af_dssp_dir/'; echo `date` end

Implement a skip if the file is already present, maybe switch to IUPRED? Or some other secondary structure predictors based on sequence. 

18) Convert DSSP to SSE summary

Start: 12.03
End: 12.07
Time: 4 mins

Command: echo `date` start; find split_af_100k_dssp_000* | xargs -P16 -I XXX sh -c 'cath-af-cli convert-dssp-to-sse-summary --dssp_dir af_dssp_dir/ --id_file XXX --sse_out_file af_100k_sse_XXX.tsv'; echo `date` end


19) Convert optimised domains to Foldseek DB

Start: 17.33
End: 18.40
Time: 2 hrs 7 mins

cath-af-cli convert-cif-to-foldseek-db --cif_dir chopped_cif_after_optimisation/ --fs_querydb_dir af_foldseek_db/ --fs_bin_path /SAN/biosciences/alphafold/foldseek/bin/foldseek

Added an option to generate a database only from a few structures. Needs to be specified as a CLI option using id_list

20) Run Foldseek against query_db

Start: 15.27
End: 17.59
Time: 2hrs 30

Command: cath-af-cli run-foldseek --fs_querydb af_foldseek_db/af_query_foldseek.db --fs_targetdb /SAN/cath/cath_v4_3_0/databases/foldseek/cath_s95/cath_s95_db --fs_bin_path /SAN/biosciences/alphafold/foldseek/bin/foldseek

21) Extract results to Foldseek summary

Start: 13.02
End: 13.07
Time: 5 mins

Command: echo `date` start; cath-af-cli convert-foldseek-output-to-summary --id_file af_100k_domainlist_post_tailchop.txt --fs_input_file fs_query_results.m8 --fs_results fs_query_results.txt; echo `date` end;

Starting from 110,250 unique domain ids in foldseek output, resulting into 110,218 hits passing the thresholds.

22) Globularity prediction

Start: 12.11
End:
Time: 

Split into 16 processes. Processing 120,113 domains

Command:
split -l 7500 af_100k_domainlist_post_tailchop.txt split_af_100k_globularity_ -a 5 -d




















