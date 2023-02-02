# Install cath-alphaflow on MYRIAD

### Version that is already in the repo
```
cd ~
mkdir repos
cd repos
git clone https://github.com/UCLOrengoGroup/cath-alphaflow.git
```


# Run cath-alphaflow on MYRIAD

### Locally
This can help you test that the provisioning and set up is working
```
cd ~
./nextflow run -resume workflows/cath-test-workflow.nf
```
### From a login node
This can help you test that the log in nodes are picking up the provisioning
```
cd ~
 ./nextflow run workflows/cath-test-workflow.nf -c platforms/ucl_myriad/nextflow.config
```
### From a cluster node
Preferred ordinary use as a login node is not tied up interactively
```
cd ~
platforms/ucl_myriad/next_nodes.sh
```
