# CATH AlphaFlow

NextFlow workflow to annotate AlphaFold model structures with CATH domains.

## Overview

This project aims to create a pipeline that provides CATH structural domain annotations for protein models predicted by AlphaFold.
The pipeline aims to take advantage of a modern workflow framework (Nextflow) to allow this process to be reusable, portable and scalable (e.g. parallel on HPC). The aim is to make these domain annotations externally available to the scientific community through a 3D-Beacons client.

## Useful Terms / Links

**Protein Domain** - proteins are often made of one or more protein domains, where each domain is defined as a globular, compact, semi-independently folding structural unit.

**CATH** - [https://www.cathdb.info](CATH) is a resource that "chops" domains from
3D coordinates from experimental data (wwPDB) and assigns these structural domains into homologous
superfamilies (i.e. domains that are related by evolution).

**AlphaFold** - [https://alphafold.ebi.ac.uk/](AlphaFold) is a protein structure prediction algorithm written by DeepMind that is capable of predicting protein 3D coordinates directly from amino acid sequence to a high degree of accuracy.

**3D Beacons** - [https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/](3D-Beacons) is an open collaboration between providers of macromolecular structure models. This collaboration aims to provide model coordinates and meta-information from all the contributing data resources in a standardised data format, on a unified platform.

**NextFlow** - [https://www.nextflow.io/](NextFlow) enables scalable and reproducible scientific workflows using software containers. It allows the adaptation of pipelines written in the most common scripting languages.

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
