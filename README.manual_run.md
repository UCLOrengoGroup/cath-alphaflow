# Running CATH-AlphaFlow-CLI - Step-by-step Tutorial
## UCL CS HPC (Test on 100,000 structures)

Initial Steps:

- load python module - `module load python/3.8.5`
- create virtual environment `python3 -m venv venv && . venv/bin/activate`
- Upgrade pip, wheel and pytest `pip install --upgrade pip wheel pytest`
- Install repo - `git clone git@github.com:UCLOrengoGroup/cath-alphaflow.git`
- Create .env file (use prod settings) - see `config.env`
- export instantclient for Oracle  `export LD_LIBRARY_PATH=~/instantclient_21_7`


## Retrieve list of 100k UniProt IDs from database 

Time: 7 minutes

Command: 

```
cath-af-cli create-dataset-uniprot-ids --uniprot_ids_csv af_uniprot_id_100k_from_g3d.list \
    --dbname gene3d_21 \
    --max_evalue 1e-50 \
    --max_records 100000 \
 ```

Output: `af_uniprot_id_100k_from_g3d.list`

## Run create-cath-dataset-from-db

Time: 4 minutes

Input: af_uniprot_id_100k_from_g3d.list

Command: 

```
cath-af-cli create-cath-dataset-from-db \
    --csv_uniprot_ids af_uniprot_id_100k_from_g3d.list \
    --dbname gene3d_21 \
    --af_cath_annotations af_100k_cath_annotations.tsv \
    --af_chainlist_ids af_100k_chainlist_ids.tsv \
    --af_domainlist_ids af_100k_domainlist_ids.tsv \
    --gene3d_crh_output af_100k.crh \
    --csv_uniprot_md5 af_100k_md5.tsv \
```

Output: 

    - af_100k_cath_annotations.tsv 
    - af_100k_chainlist_ids.tsv 
    - af_100k_chainlist_ids.tsv 
    - af_100k.crh 
    - af_100k_md5.tsv

## Create GCloud Manifest File

Time: less than 1 min

Input: `af_100k_chainlist_ids.tsv`

Output: `af_100k_manifest.txt`

```
grep -v 'chain_id' af_100k_chainlist_ids.tsv | sort | uniq | tr -d '\r' | awk '{print "gs://public-datasets-deepmind-alphafold-v4/"$1".cif"}' > af_100k_manifest.txt
```

## Download structures from GCloud

time: 2hrs and 20mins

Input: af_100k_manifest.txt

Output: `af_raw_cif` folder containing unchopped AlphaFold chains

```
mkdir af_cif_raw
cat af_100k_manifest.txt | (/share/apps/google-cloud-sdk/bin/gsutil -o GSUtil:parallel_process_count=1 -m cp -I af_cif_raw/ || echo "Ignoring non-zero exit code: \$?")
```

The step outputs the number of UniProt IDs not available in AlphaFoldDB.

```
CommandException: 6532 files/objects could not be transferred. Total downloaded CIF: 62,212
```

## Create new AFChainList after GCloud download. 

Rationale: 

The `cif_to_fasta` module in `cath-af-cli` expects a list of existing files, so we need to feed it just the chains that are in AFDB and are available on the filesystem.

```
cd af_cif_raw
find *.cif > ../af_100k_chainlist_after_gsutil.txt
sed -i 's/\.cif//g' af_100k_chainlist_after_gsutil.txt
echo af_chain_id | cat - af_100k_chainlist_after_gsutil.txt > temp && mv temp af_100k_chainlist_after_gsutil.txt
```


## Convert CIF to FASTA 

Process: 62,212 CIF files -> 62,212 FASTA files + merged.fasta
Thoughts: Too slow on 1CPU (i.e. 30mins for 10k sequences, splitting and parallelizing)
Time:  1hr 9 mins

Input: 

```
af_cif_raw/ 
af_100k_chainlist_after_gsutil.txt
```

Output: 
`af_fasta_raw`  - folder containing a FASTA file for each chain available in the filesystem
`merged.fasta`  - file containing a multiFASTA of all processed CIF files.

```
mkdir af_fasta_raw
split -l 10000 af_100k_chainlist_after_gsutil.txt splitlist
find splitlista* | xargs -P4 -I XXX sh -c 'cath-af-cli convert-cif-to-fasta --cif_in_dir af_cif_raw/ \
    --id_file XXX \
    --fasta_out_dir af_fasta_raw/'
```

Careful: the module expects each chunk to have a header, so the first structure doesn't get processed if we use chunks instead of the first file.

## Create MD5 from FASTA

Time: seconds

Input: `merged.fasta`

Output: `af_100k_cif_raw_md5.txt`

```
cath-af-cli create-md5 --fasta merged.fasta \
    --uniprot_md5_csv af_100k_cif_raw_md5.txt
```

## Filter CRH 

Rationale:

- Extract and further process only matches between CIF-md5 and chopping MD5
- Remove all CRH entries that rely on a missing AFDB entry.

Time: 1 second

Issue: The files generated in step 6 have some carriage returns (like ^M) that needs removing using tr -d '\r'

Temporary fix: Strip these characters before doing the grep command.

Command: 

```
grep -v 'af_chain_id' af_100k_cif_raw_md5.txt | awk '{print $3}' | tr -d '\r' | grep -F -f - af_100k.crh > af_100k_crh_after_filtering.crh
sort af_100k_crh_after_filtering.crh | uniq > af_100k_crh_after_filtering_uniq.crh
awk '{print $1}' af_100k_crh_after_filtering_uniq.crh | grep -F -f - af_100k_cif_raw_md5.txt | awk '{print $1}' | grep -F -f - af_100k_domainlist_ids.txt > af_100k_domainlist_after_md5_filter.txt
echo af_domain_id | cat - af_100k_domainlist_after_md5_filter.txt > temp && mv temp af_100k_domainlist_after_md5_filter.txt
```

Output: `af_100k_crh_after_filtering.crh`

## Chop CIF before optimisation

### WARNING!! 

### Unnecessary step as `optimise-domain-boundaries` acts on chains, not on domains, so we can optimise first and chop once.

Running on xargs with 12 processes with a uniq-ed crh file and domain list (120,129 domains)
Time:  1hr 45 minutes

Command: 

```
split -l 7500 af_100k_domainlist_after_md5_filter.txt split_crh_uniq_ -a 5 -d 
find split_crh_uniq_000* | xargs -P12 -I XXX sh -c 'cath-af-cli chop-cif --cif_in_dir af_cif_raw/ \
    --id_file XXX \
    --cif_out_dir chopped_cif_before_optimization/'
```

Issues Found: 
Remove af_domain_id from header, program can't parse it or add exception, as it skips the first entry but it crashes as it can't find af_domain_id

Stats: Chopped 120,129 domains

## Optimise boundaries
On 16 threads, 7500 entries each

Time: 1hr 1 min

Comments/TODO: 
 - Include percentage of delta length (i.e. how much the boundaries are changing?)
 - Remove af_domain_id from header, it crashes the process.
 - Optimise boundaries works on chains, not on domains. 

The initial chopping is redundant, keeping it in this instance as it's interesting for debugging.

Commands:

```
split -l 7500 af_100k_domainlist_after_md5_filter_uniq.txt split_af_100k -a 5 -d
find split_af_100k000* | xargs -P16 -I XXX sh -c 'cath-af-cli optimise-domain-boundaries --af_domain_list XXX \
    --af_chain_mmcif_dir af_cif_raw/ \
    --af_domain_list_post_tailchop af_100k_domainlist_post_tailchop_XXX.txt \
    --af_domain_mapping_post_tailchop af_100k_domain_mapping_post_tailchop_XXX.tsv \
    --status_log optimise_boundaries_status_log_XXX'
```

Afterwards:

- Concatenate logs
- Concatenate af_domain_list_post_tailchop 
- Concatenate af_domain_mapping_post_tailchop 
- Remove duplicated headers

```
cat af_100k_domain_mapping_post_tailchop_split_af_100k000* | grep -v 'af_domain_id' > af_100k_domain_mapping_post_tailchop.tsv
cat af_100k_domainlist_post_tailchop_split_af_100k000* grep -v 'af_domain_id' > af_100k_domainlist_post_tailchop.txt
cat optimise_boundaries_status_log_split_af_100k000* > optimise_boundaries_status_log
rm af_100k_domain_mapping_post_tailchop_split_af_100k000*
rm optimise_boundaries_status_log_split_af_100k000*
rm af_100k_domainlist_post_tailchop_split_af_100k000*
rm split_af_100k000*
```

Stats:

- 84,604 have unchanged boundaries
- 35,378 have adjusted boundaries
- 129 have failed due to not having residues above pLDDT 70



## Chop CIF files after boundaries optimisation

Time: 1hr 46 mins on 12 threads

```
split -l 7500 af_100k_domainlist_post_tailchop.txt af_100k_domainlist -a 5 -d
find af_100k_domainlist00* | xargs -P12 -I XXX sh -c 'cath-af-cli chop-cif --cif_in_dir af_cif_raw/ \
    --id_file XXX \
    --cif_out_dir chopped_cif_after_optimisation/'
```

TODO: Add StatusLog to ALL MODULES

Chopped 120,112 domains with optimised boundaries

## Calculate pLDDT and LUR report

Time:  59 mins on 16 threads

```
split -l 7500 af_100k_domainlist_post_tailchop.txt split_af_100k_plddt_ -a 5 -d
find split_af_100k_plddt_000* | xargs -P16 -I XXX sh -c 'cath-af-cli convert-cif-to-plddt-summary --cif_in_dir af_cif_raw/ \
    --id_file XXX \
    --plddt_stats_file af_100k_plddt_summary_XXX.tsv'
cat af_100k_plddt_summary_split_af_100k_plddt_*.tsv > af_100k_plddt_summary.tsv
```

## Create DSSP Files

Time: 4 hrs

Requirements: boost library. On CS HPC can be loaded as a module.

```
module load boost/1.71.0
mkdir af_dssp_dir
find split_af_100k_dssp_000* | xargs -P16 -I XXX sh -c 'cath-af-cli convert-cif-to-dssp --cif_in_dir af_cif_raw/ \
    --id_file XXX \
    --dssp_out_dir af_dssp_dir/'
```

Comments: Implement a skip if the file is already present, maybe switch to IUPRED? Or some other secondary structure predictors based on sequence. 

## Convert DSSP to SSE summary

Time: 4 mins

Command: 
```
find split_af_100k_dssp_000* | xargs -P16 -I XXX sh -c 'cath-af-cli convert-dssp-to-sse-summary --dssp_dir af_dssp_dir/ \
    --id_file XXX \
    --sse_out_file af_100k_sse_XXX.tsv'
```

## Convert optimised domains to Foldseek DB

Time: 2 hrs 7 mins on a single thread, RAM:1G

Command:
```
cath-af-cli convert-cif-to-foldseek-db --cif_dir chopped_cif_after_optimisation/ \
    --fs_querydb_dir af_foldseek_db/ \
    --fs_bin_path /SAN/biosciences/alphafold/foldseek/bin/foldseek
```

Comment:
- Added an option to generate a database only from a few structures. Needs to be specified as a CLI option using id_list
- Investigate why m_core on CS is not seen by Foldseek. Increasing RAM available speeds up the process significantly.



## Run Foldseek of Query_Domain_Structures against CATH_S95 Representatives

Time: 2hrs 30

Command: 
```
cath-af-cli run-foldseek --fs_querydb af_foldseek_db/af_query_foldseek.db \
    --fs_targetdb /SAN/cath/cath_v4_3_0/databases/foldseek/cath_s95/cath_s95_db \
    --fs_bin_path /SAN/biosciences/alphafold/foldseek/bin/foldseek 
```

## Convert Foldseek Results to Summary

Time: 5 mins

Command: 
```
cath-af-cli convert-foldseek-output-to-summary --id_file af_100k_domainlist_post_tailchop.txt \
    --fs_input_file fs_query_results.m8 \
    --fs_results fs_query_results.txt
```

Starting from 110,250 unique domain ids in foldseek output, resulting into 110,218 hits passing the thresholds.

## Globularity prediction

Start: 12.11
End: 16.17
Time: 2hrs 6 mins

Split into 16 processes. Processing 120,113 domains

Command:
```
split -l 7500 af_100k_domainlist_post_tailchop.txt split_af_100k_globularity_ -a 5 -d
find split_af_100k_globularity_000* | xargs -P16 -I XXX sh -c 'cath-af-cli measure-globularity --af_domain_list XXX \
    --af_chain_mmcif_dir af_cif_raw/ \
    --af_domain_globularity af_100k_globularity_XXX.txt'
```

## Combine results
























