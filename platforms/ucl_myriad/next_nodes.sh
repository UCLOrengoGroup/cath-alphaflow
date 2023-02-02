#!/bin/bash -l

# Batch script to run a serial job under SGE.

# Request ten minutes of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=0:10:0

# Request 1 gigabyte of RAM (must be an integer followed by M, G, or T)
#$ -l mem=1G

# Request 15 gigabyte of TMPDIR space (default is 10 GB - remove if cluster is diskless)
#$ -l tmpfs=15G

# Set the name of the job.
#$ -N cath_ctrl

# Set the working directory to somewhere in your scratch space.  
#  This is a necessary step as compute nodes cannot write to $HOME.
# Replace "<your_UCL_id>" with your UCL user ID.
#$ -wd /home/ucbtlcr/Scratch/workspace

cd ~
#path to create the venv
export CATH_VENV_PATH=~/cath_venv
#path to install gcloud
export CATH_GCLOUD=~/gcloud
export CATH_GCLOUD_USER='rachelalcraft@gmail.com'
#path to clone the repo to
export REPOS_PATH=~/repos
export REPO_NAME=cath-alphaflow
export CATH_REPO_PATH=${REPOS_PATH}/${REPO_NAME}
#path to store bulk data outside the repo path if required
#eg if repo is on home and bulk storage is a project folder somewhere else?
export CATH_DATA_PATH=~/cath-workspace-test
#put python, java and google cloud sdk in the path
module load python3
module load java/openjdk-11/11.0.1
export PATH=${CATH_GCLOUD}/google-cloud-sdk/bin:${PATH}
. ${CATH_VENV_PATH}/bin/activate
pip install -e .
cath-af-cli
cd ${CATH_REPO_PATH}

 ./nextflow run workflows/cath-test-workflow.nf -c platforms/ucl_myriad/nextflow.config
