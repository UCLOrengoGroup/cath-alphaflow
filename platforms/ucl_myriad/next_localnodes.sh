#!/bin/bash -l

SCRIPT_DIR=~/repos/cath-alphaflow/platforms/ucl_myriad
source ${SCRIPT_DIR}/include
module load python3
module load java/openjdk-11/11.0.1
./nextflow run workflows/cath-test-workflow.nf -c platforms/ucl_myriad/nextflow.config