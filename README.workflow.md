This is the initial overview for the CATH AlphaFold NextFlow workflow.

#### 1a. Global dataset (AF2 sequences)

Name: CREATE_ALL_AF2_CHAIN_FASTA

Input:

- `PARAM:AF2_DOWNLOAD_URL` -- AlphaFold Download (gs://public-datasets-deepmind-alphafold/sequences.fasta)

Output:

- `FILE:ALL_AF2_CHAIN_FASTA`

#### 1b. Global dataset (Foldseek S95 library)

Name: CREATE_DOMAIN_S95_PDB_DIR

Input:

- `FILE:DOMAIN_LIST_S95`
- `DIRECTORY:CATH_PDB_DOMAINS` -- v4.3
- `PARAM:S95_DIR` -- folder containing S95 PDB files

Output:

- `DIRECTORY:FOLDSEEK_S95_PDB_DB`

Process:

- `rsync --files-from <FILE:DOMAIN_LIST_S95> <DIRECTORY:CATH_PDB_DOMAINS> <DIRECTORY:FOLDSEEK_S95_PDB_DB>`

#### 1c. Global dataset (Foldseek S95 library)

Name: CREATE_FOLDSEEK_S95_LIBRARY

Input:

- `PARAM:S95_DIR` -- folder containing S95 PDB files

Output:

- `FILE:FOLDSEEK_S95_LIBRARY`

Process:

- `foldseek createdb <s95_dir> <s95_db>`

Notes:

- currently in `/SAN/cath/cath_v4_3_0/databases/foldseek/cath_s95/`

#### 2a. Local dataset (UniProt IDs)

Name: CREATE_DATASET_UNIPROT_IDS

Input:

- `SOURCE:OracleDB` -- `GENE3D_21.CATH_DOMAIN_PREDICTIONS{_EXTRA}`
- `FILTER` -- CONDITIONAL_EVALUE <= 1e-50, top 10k hits (ordered by increasing evalue)

Output:

- `FILE:CSV_UNIPROT_IDS` -- `UniProt_ID`

#### 2b. Local dataset (CATH domains predictions)

Name: CREATE_DATASET_CATH_FILES

Input:

- `SOURCE:OracleDB` -- `GENE3D_21.CATH_DOMAIN_PREDICTIONS{_EXTRA}`
- `FILE:CSV_UNIPROT_IDS` -- `(UniProt_ID)`

Output:

- `FILE:CSV_UNIPROT_MD5` -- `(UniProt_ID, md5)`
- `FILE:CRH_OUTPUT` -- `(md5, CATH_domain_ID, bitscore, boundaries, resolved_boundaries)`
- `FILE:AF2_DOMAIN_LIST` -- `(AF_domain_ID_orig)`
- `FILE:AF2_CHAIN_LIST` -- `(AF_chain_ID)`
- `FILE:AF2_CATH_ORIG_ANNOTATIONS` -- `(AF_domain_ID_orig, CATH_domain_ID, UniProt_ID, md5, bitscore, resolved_boundaries, sfam_id, class_id)`

#### 2c. Local dataset (AF2 Structure)

Name: CREATE_DATASET_AF2_FILES

Input:

- `FILE:CSV_UNIPROT_MD5`
- `SOURCE:HPC_ENVIRONMENT` -- Google Cloud Storage (GCS), Computer Science (CS)
- `FILE:AF2_CHAIN_LIST`

Output:

- `FILE:ALL_CHAIN_FASTA` -- one FASTA file containing all AF2 chains
- `DIRECTORY:AF2_CHAIN_MMCIF` -- one mmCIF file per AF2 chain

Process:

- `GCS` -- `cat [manifest file] | gsutil -m cp -I . `
- `CS` -- (local): symlinking (check if it exists)

#### 3a. Run SETH (Chain sequence)

Name: CREATE_ANNOTATION_CHAIN_DISORDER

Input:

- `FILE:ALL_CHAIN_FASTA`

Output:

- `FILE:SETH_CHAIN_OUTPUT_FILE` -- Per-residue SETH scores (disordered 0 -> 1 ordered)

Process:

- `SETH_1.py -i <your input fasta file name> -o <the desired name of your output file>`

#### 3b. Filter (Disordered Chain)

Name: FILTER_DOMAINLIST_BY_CHAIN_DISORDER

Input:

- `FILE:SETH_CHAIN_OUTPUT_FILE`
- `FILTER` -- `(SETH_GLOBAL_CHAIN_DISORDER, SETH_LOCAL_DOMAIN_DISORDER)`
- `FILE:AF2_DOMAIN_LIST`

Output:

- `FILE:AF2_DOMAIN_LIST`
- `FILE:SETH_GLOBAL_CHAIN_DISORDER` -- `(AF_domain_ID_orig, seth_global_disorder_score)`

Process:

- `python filter_disorder.py`

#### 4. Optimize domain boundaries (chop tails)

Name: OPTIMISE_DOMAIN_BOUNDARIES

Input:

- `FILE:AF2_CATH_ORIG_ANNOTATIONS`
- `FILE:AF2_CHAIN_MMCIF`
- `FILE:AF2_DOMAIN_LIST`

Output:

- `FILE:AF2_DOMAIN_LIST_POST_TAILCHOP` -- `(AF_domain_ID_tailchop)`
- `FILE:AF2_DOMAIN_MAPPING_POST_TAILCHOP` -- `(AF_domain_ID_orig, AF_domain_ID_tailchop)`

Process:

- `python optimize_boundaries.py`

#### 5a. Run SETH (Domain sequence)

Name: CREATE_ANNOTATION_DOMAIN_DISORDER

Input:

- `FILE:ALL_CHAIN_FASTA`
- `FILE:AF2_DOMAIN_LIST_POST_TAILCHOP`

Output:

- `FILE:SETH_DOMAIN_OUTPUT_FILE` -- Per-residue SETH scores (disordered 0 -> 1 ordered)

Process:

- `SETH_1.py -i <your input fasta file name> -o <the desired name of your output file>`

#### 5b. Filter (Disordered Domain)

Name: FILTER_DOMAINLIST_BY_DISORDER

Input:

- `FILE:SETH_DOMAIN_OUTPUT_FILE`
- `FILTER` -- `(SETH_GLOBAL_CHAIN_DISORDER, SETH_LOCAL_DOMAIN_DISORDER)`
- `FILE:AF2_DOMAIN_LIST_POST_TAILCHOP`

Output:

- `FILE:AF2_DOMAIN_LIST_POST_DOMAIN_DISORDER`
- `FILE:AF2_SETH_ANNOTATIONS` -- `(AF_domain_ID_tailchop, seth_global_disorder_score)`

Process:

- `python filter_disorder.py`

#### 6. Filter: AF2 Quality (pLDDT, LUR)

Name: FILTER_DOMAINLIST_BY_AF2_QUALITY

Input:

- `FILE:AF2_CHAIN_MMCIF`
- `FILE:AF2_DOMAIN_LIST_POST_DOMAIN_DISORDER`

Output:

- `FILE:AF2_DOMAIN_LIST_POST_AF_QUALITY`
- `FILE:AF2_QUALITY_ANNOTATIONS` -- `(AF_domain_ID_tailchop, PLDDT_average, LUR_score)`

Process:

- `python filter_af2_quality.py`

#### 7. Filter: AF2 Packing (packing, Surf/Vol)

Name: FILTER_DOMAINLIST_BY_AF2_PACKING

Input:

- `FILE:AF2_CHAIN_MMCIF`
- `FILE:AF2_DOMAIN_LIST_POST_AF_QUALITY`

Output:

- `FILE:AF2_DOMAIN_LIST_POST_AF_PACKING`
- `FILE:AF2_PACKING_ANNOTATIONS` -- `(AF_domain_ID_tailchop, packing_score, surf_vol_score)`

Process:

- `python filter_af2_packing.py`

#### 8. Filter: SSE (DSSP+secmake)

Name: FILTER_DOMAINLIST_BY_SSE

Input:

- `FILE:AF2_CHAIN_MMCIF`
- `FILE:AF2_DOMAIN_LIST_POST_AF_PACKING`

Output:

- `FILE:AF2_DOMAIN_LIST_POST_SSE`
- `FILE:AF2_SSE_ANNOTATIONS` -- `(AF_domain_ID_tailchop, sse_number)`

Process:

- `python3 filter_sse.py`

#### 9. Run Chopping

Name: CHOP_AF2_DOMAINS

Input:

- `FILE:AF2_DOMAIN_LIST_POST_SSE` -- `af_<Uniprot_ID>/<start>-<stop>`
- `DIRECTORY:AF2_CHAIN_MMCIF`

Output:

- `DIRECTORY:AF2_DOMAIN_MMCIF`
- `FILE:AF2_DOMAIN_MAPPING_POST_CHOPPING` -- `(AF_domain_ID_tailchop, AF_domain_ID_chopped)`

Process:

```
submit.sh
/home/ucbcisi/work/2022_06_29.alphafold_rechop_corrected_domains/submit.sh
```

#### 10. Foldseek S95

Name: CREATE_ANNOTATION_FOLDSEEK_S95

Input:

- `FILE:FOLDSEEK_S95_DB`
- `FILE:AF2_DOMAIN_LIST_POST_SSE`
- `DIRECTORY:AF2_DOMAIN_MMCIF`

Output:

- `FILE:AF2_FOLDSEEK_ANNOTATIONS` -- `(AF_domain_ID_tailchop, foldseek_bitscore, foldseek_overlap)`

Process:

```
foldseek createdb AF2_DOMAIN_MMCIF AF2_DOMAIN_MMCIF_DB
foldseek search AF2_DOMAIN_MMCIF_DB FOLDSEEK_S95_DB
foldseek convertalis AF2_DOMAIN_MMCIF_DB FOLDSEEK_S95_DB
```

#### 11. Create Table

Name: CREATE_RESULTS_TABLE

Input:

- `FILE:AF2_DOMAIN_LIST_POST_SSE`
- `FILE:AF2_SETH_ANNOTATIONS`
- `FILE:AF2_QUALITY_ANNOTATIONS`
- `FILE:AF2_PACKING_ANNOTATIONS`
- `FILE:AF2_CATH_ORIG_ANNOTATIONS`
- `FILE:AF2_SSE_ANNOTATIONS`
- `FILE:AF2_FOLDSEEK_ANNOTATIONS`

Output:

- `FILE:CATH_AF2_TABLE`

Process:

- `python collate_data_to_table.py`
