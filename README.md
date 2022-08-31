# CATH AlphaFlow

Workflow to annotate AlphaFold2 model structures with CATH domains.

#### 1. Global dataset

Input: 
- SOURCE: AlphaFold Download URL (gs://public-datasets-deepmind-alphafold/sequences.fasta)

Output:
- FILE:ALL_AF2_CHAIN_FASTA

#### 2. Local dataset (CATH domains predictions)

Input:
- SOURCE:OracleDB -- `GENE3D_21.CATH_DOMAIN_PREDICTIONS{_EXTRA}`
- FILTER:CONDITIONAL_EVALUE <= 1e-50, top 10k hits (ordered by increasing evalue) 

Output: 
- FILE:CSV_UNIPROT_MD5 -- `UniProt_ID, md5`
- FILE:CRH_OUTPUT -- `md5, CATH_domain_ID, bitscore,boundaries, resolved_boundaries`
- FILE:AF2_DOMAIN_LIST
- FILE:AF2_CHAIN_LIST
- FILE:AF2_SFAM_PREDICTION
- FILE:AF2_CLASS_PREDICTION

#### 3. Local dataset (AF2 Structure)

Input: 
- CSV_UNIPROT_MD5 (UniProt_ID)
- Source: GCS,CS
- AF2_CHAIN_LIST

Output:
- ALL_CHAIN_FASTA 
- AF2_CHAIN_MMCIF

Process:
- GCS: cat [manifest file] | gsutil -m cp -I .   
- CS (local): symlinking (check if it exists) 
     
#### 4. Run SETH (sequences)

Input:
- ALL_CHAIN_FASTA

Output:
- SETH_OUTPUT_FILE: Per-residue SETH scores (disordered 0 -> 1 ordered)

Process:
- SETH_1.py -i <your input fasta file name> -o <the desired name of your output file>

#### 5. Filter (Disorder Chain/Domain)

Input:
- SETH_OUTPUT_FILE 
- Filter criteria (SETH_GLOBAL_CHAIN_DISORDER,SETH_LOCAL_DOMAIN_DISORDER)
- AF2_DOMAIN_LIST

Output:
- AF2_DOMAIN_LIST
- SETH_GLOBAL_CHAIN_DISORDER
- SETH_LOCAL_DOMAIN_DISORDER

Process: 
- `python filter_disorder.py`

#### 6. Optimize domain boundaries

Input:
- CRH_OUTPUT
- AF2_CHAIN_MMCIF
- AF2_DOMAIN_ID

Output:
- AF2_DOMAIN_LIST_POST_TAILS
- AF2_DOMAIN_MAPPING_POST_TAILS
            
Process:
- `python optimize_boundaries.py`

#### 7. Filter: AF2 Quality (pLDDT,LUR)

Input:
- AF2_CHAIN_MMCIF
- AF2_DOMAIN_LIST_POST_TAILS

Output:
- AF2_DOMAIN_LIST_POST_TAILS
- AVG_PLDDT
- PERC_LUR

Process:
- `python filter_af2_quality.py`

#### 8. Filter: AF2 Packing (packing, Surf/Vol)

Input:
- AF2_CHAIN_MMCIF
- AF2_DOMAIN_LIST_POST_TAILS

Output:
- AF2_DOMAIN_LIST_POST_TAILS
- PACKING_SCORE
- SURF_VOL_SCORE

Process:
- `python filter_af2_packing.py`

#### 9. Filter: SSE (DSSP+secmake)

Input:
- AF2_CHAIN_MMCIF
- AF2_DOMAIN_LIST_POST_TAILS

Output:
- AF2_DOMAIN_LIST_POST_TAILS
- SSE_NUM

Process:
- `python3 filter_sse.py`

#### 10. Run Chopping 

Input:
- CRH_OUTPUT
- CSV_UNIPROT_MD5
- AF2_DOMAIN_LIST_POST_TAILS
- AF2_CHAIN_MMCIF

Output:
- AF2_DOMAIN_MMCIF
- AF2_DOMAIN_LIST_POST_TAILS

Process:
```
submit.sh 
/home/ucbcisi/work/2022_06_29.alphafold_rechop_corrected_domains/submit.sh
```

#### 11. Foldseek S95

Input:
- FOLDSEEK_S95_DB
- AF2_DOMAIN_LIST_POST_TAILS
- AF2_DOMAIN_MMCIF

Output:
- FOLDSEEK_BITSCORE
- FOLDSEEK_OVERLAP

Process:
```
foldseek createdb AF2_DOMAIN_MMCIF AF2_DOMAIN_MMCIF_DB 
foldseek search AF2_DOMAIN_MMCIF_DB FOLDSEEK_S95_DB
foldseek convertalis AF2_DOMAIN_MMCIF_DB FOLDSEEK_S95_DB
```

#### 12. Create Table

Input:
- AF2_DOMAIN_LIST_POST_TAILS
- SETH_LOCAL_DOMAIN_DISORDER
- AVG_PLDDT
- PERC_LUR
- PACKING_SCORE
- SURF_VOL_SCORE
- AF2_SFAM_PREDICTION
- AF2_CLASS_PREDICTION
- SSE_NUM
- FOLDSEEK_BITSCORE
- FOLDSEEK_OVERLAP
            
Output:
- CATH_AF2_TABLE

Process:
- `python collate_data_to_table.py`

