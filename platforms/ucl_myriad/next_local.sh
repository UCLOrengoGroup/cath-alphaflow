#!/bin/bash -l

SCRIPT_DIR=~/repos/cath-alphaflow/platforms/ucl_myriad
source ${SCRIPT_DIR}/include
. ${CATH_VENV_PATH}/bin/activate

gcloud auth login ${CATH_GCLOUD_USER}

source ~/.bashrc
module load python3
module load java/openjdk-11/11.0.1

 ./nextflow run -resume workflows/cath-test-workflow.nf



