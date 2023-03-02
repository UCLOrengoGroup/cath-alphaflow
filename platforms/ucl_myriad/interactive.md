
# Interactive kick off for cath-aplhaflow

## PATH INCLUDES and END VARIABLES
```
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
```

## INSTALLATIONS AND CLONES
Initial install of everything using pip, wget and git clone
```
#venv setup, install dependencies
python3 -m venv ${CATH_VENV_PATH}
. ${CATH_VENV_PATH}/bin/activate
pip install --upgrade pip wheel pytest nextflow
#get gcloud utils
mkdir -p ${CATH_GCLOUD}
cd ${CATH_GCLOUD}
wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-415.0.0-linux-x86_64.tar.gz
tar -xf google-cloud-cli-415.0.0-linux-x86_64.tar.gz
#get cath-af repo
mkdir -p ${REPOS_PATH}
cd ${REPOS_PATH}
git clone https://github.com/UCLOrengoGroup/cath-alphaflow.git
cd ${REPO_NAME}
pip install -e .
cath-af-cli
```

## Install NEXTFLOW
```
cd ${CATH_REPO_PATH}
curl -s https://get.nextflow.io | bash
```
# run locally on myriad
```
cd ${CATH_REPO_PATH}
 ./nextflow run -resume workflows/cath-test-workflow.nf
```
 # or on the nodes
 ```
 cd ${CATH_REPO_PATH}
 ./nextflow run workflows/cath-test-workflow.nf -c platforms/ucl_myriad/nextflow.config
```
 # or on the nodes from a node
 ```
 cd ${CATH_REPO_PATH}
 platforms/ucl_myriad/next_nodes.sh

```


