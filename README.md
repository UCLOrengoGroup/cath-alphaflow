# CATH AlphaFlow

NextFlow workflow to annotate AlphaFold model structures with CATH domains.

## Overview

This project aims to create a pipeline that provides CATH structural domain annotations for protein models predicted by AlphaFold.
The pipeline aims to take advantage of a modern workflow framework (Nextflow) to allow this process to be reusable, portable and scalable (e.g. parallel on HPC). The aim is to make these domain annotations externally available to the scientific community through a 3D-Beacons client.

## Useful Terms / Links

**Protein Domain** - proteins are often made of one or more protein domains, where each domain is defined as a globular, compact, semi-independently folding structural unit.

**CATH** - [CATH](https://www.cathdb.info) is a resource that "chops" domains from
3D coordinates from experimental data (wwPDB) and assigns these structural domains into homologous
superfamilies (i.e. domains that are related by evolution).

**AlphaFold** - [AlphaFold](https://alphafold.ebi.ac.uk/) is a protein structure prediction algorithm written by DeepMind that is capable of predicting protein 3D coordinates directly from amino acid sequence to a high degree of accuracy.

**3D Beacons** - [3D Beacons](https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/) is an open collaboration between providers of macromolecular structure models. This collaboration aims to provide model coordinates and meta-information from all the contributing data resources in a standardised data format, on a unified platform.

**NextFlow** - [NextFlow](https://www.nextflow.io/) enables scalable and reproducible scientific workflows using software containers. It allows the adaptation of pipelines written in the most common scripting languages.

## Project

The NextFlow workflows are kept here:

```
workflows/
```

Most of the non-trivial steps within the workflow are carried out via Python CLI commands, which can be found in:

```
cath_alphaflow/
```

Installing the Python dependencies in a virtual environment:

```
python3 -m venv venv
. venv/bin/activate
pip install --upgrade pip wheel pytest
pip install -e .
```

Accessing the Python CLI tools

```
. venv/bin/activate
cath-af-cli
```

Running the local tests

```
pytest
```

## Nextflow

Install NextFlow (more details [here](https://www.nextflow.io/index.html#GetStarted))

```bash
curl -s https://get.nextflow.io | bash
```

Run the basic workflow via:

```bash
./nextflow run -resume workflows/cath-test-workflow.nf
```

Notes:

- the `-resume` flag instructs NextFlow to use cached data where possible
- running this workflow will currently fail without the steps below

One of the steps in the basic workflow tries to download files from Google Storage. This requires you to be logged into Google via the `gcloud` utility (more details [here](https://cloud.google.com/sdk/docs/install)).

```bash
gcloud auth application-default login
```

Now running the workflow should run successfully...

```bash
$ ./nextflow run -resume workflows/cath-test-workflow.nf
N E X T F L O W ~ version 22.10.3
Launching `workflows/cath-test-workflow.nf` [maniac_bhabha] DSL2 - revision: 21b594e11a
executor > local (48)
[21/0dcefb] process > uniprot_domain_to_uniprot (51) [100%] 51 of 51, cached: 51 ✔
[4c/3b022d] process > create_af_manifest_file (31) [100%] 51 of 51, cached: 51 ✔
[da/f74788] process > retrieve_af_chain_cif_files (47) [100%] 51 of 51, cached: 49 ✔
[40/074465] process > chop_cif (42) [100%] 46 of 46 ✔
Completed at: 17-Jan-2023 15:25:49
Duration : 1m 8s
CPU hours : 0.3 (32.1% cached)
Succeeded : 48
Cached : 151

```
