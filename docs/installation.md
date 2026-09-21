# Install SWARM

## Download the code and models

Simply clone from github (install lfs to download large h5 files):

```bash
git lfs install
git clone https://github.com/comprna/SWARM/ && cd SWARM 
```

If git lfs cannot be installed, download the models from the Zenodo link:

```bash
git clone https://github.com/comprna/SWARM/ && cd SWARM 
rm -rf SWARM_models
wget 'https://zenodo.org/records/22123294/files/SWARM_models.tar.gz' -O SWARM_models.tar.gz
tar -xzf SWARM_models.tar.gz && rm -f SWARM_models.tar.gz
```

## Compile SWARM preprocessing

```bash
cd SWARM_scripts/preprocess/
#build and compile htslib, slow5tools, SWARM_preprocess
bash build.sh
```

## Dependencies

SWARM supports GPU inference with tensorflow, tested with versions 2.8.0 and 2.15.0. 

### Using pre-installed tensorflow
If your HPC has a tensorflow module, simply load tensorflow and use the loaded python path for creating venv:
```
module load tensorflow/2.15.0
python3 -m venv swarm_env
source swarm_env/bin/activate
python3 -m pip install pysam==0.22.1 numpy==1.26.2 pandas==2.2.0 scikit-learn==1.4.0

# make sure to activate the venv before running SWARM read-level and site-level prediction
module load tensorflow/2.15.0
source /PATH/TO/swarm_env/bin/activate
```

### Using containerised environment
If tensorflow with GPU configuration is not pre-installed or there are issues with dependencies, we provide a containerised environment with tensorflow and pysam:
https://zenodo.org/records/22123294
```
# You can use singularity to run read-level and site-level prediction:
singularity exec --nv tensorflow_24.01-tf2-py3-pysam.sif python3 script.py ...

# Environment tested on NCI gadi HPC using singularity version 3.11.3 and NVIDIA Volta GPUs
```

