#!/usr/bin/env bash

# Myriad version / Rachel version

#initial install of everything using pip, wget and git clone
#also currently runs cath-af pytest and nextflow hello example

#SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
SCRIPT_DIR=~/repos/cath-alphaflow/platforms/ucl_myriad
source ${SCRIPT_DIR}/include

#load java and python modules
module load python3
module load java/openjdk-11/11.0.1

#test nextflow
cd ${CATH_REPO_PATH}
./nextflow run hello
