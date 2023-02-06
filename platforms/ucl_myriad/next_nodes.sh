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
SCRIPT_DIR=~/repos/cath-alphaflow/platforms/ucl_myriad
source ${SCRIPT_DIR}/include
. ${CATH_VENV_PATH}/bin/activate

gcloud auth login ${CATH_GCLOUD_USER}

source ~/.bashrc
module load python3
module load java/openjdk-11/11.0.1

 ./nextflow run workflows/cath-test-workflow.nf -c platforms/ucl_myriad/nextflow.config
