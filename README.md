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

**dssp** - [dssp](https://github.com/PDB-REDO/dssp) applies secondary structure to protein sequence. 
Needs first installed: [libmcfp](https://github.com/mhekkel/libmcfp)

**oracle** - Some datasources come optionally from an oracle database, for which there is access only within UCL. Instructions for UCL access to the database are [here](README.oracle.md)

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

It's also possible to install Nextflow via a python wrapper package, this automatically adds it to the path

```bash
pip install --upgrade nextflow
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

## Configuring for your execution environment

The platform folder contains subfolders for each execution environment the pipeline has been tested with so far. You can create a new subfolder and copy and customise the files found in one of the existing subfolders. Each subfolder contains:

- an `include` file which sets environment variables
- an `install` script which performs one-off package and dependency installation
- a `source` script which should to sourced immediately prior to launching the pipeline
- a `nextflow.config` file containing overrides to the base nextflow.config

To run nextflow using a custom configuration use the following from the toplevel folder of the repository, ie the folder containing the base nextflow.config, such that it will be automatically picked up but selectively overridden by the explicitly specified platform specific config:

```
nextflow run -resume -c ./platforms/<your_platform>/nextflow.config ./workflows/cath-test-workflow.nf
```

You would typically run this from an HPC login node.

---

## Tracing and visualising timelines of the pipeline
https://www.nextflow.io/docs/latest/tracing.html

The production of a pipeline image needs the installation of graphviz
```
sudo apt install graphviz
```
Call nextflow with -with-dag, e.g.
```
./nextflow run workflows/cath-test-workflow.nf -with-dag cath-test-workflow.png
```
An html report that shows the process timings, gantt chart style
```
./nextflow run workflows/cath-test-workflow.nf -with-timeline cath-test-workflow-gantt.html
```
An html html report with plotly summaries of processes: 
```
./nextflow run workflows/cath-test-workflow.nf -with-report cath-test-workflow-plots.html
```