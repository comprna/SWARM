# Install SWARM

## Using containerised environment
If tensorflow with GPU configuration is not pre-installed or there are issues with dependencies, we provide a containerised environment:
https://zenodo.org/records/22123294

You can run scripts from the SWARM repo using singularity and /opt/SWARM/path/to/script
```
# For example to run SWARM_read_level.py located at SWARM/SWARM_scripts/SWARM_read_level.py
singularity exec --nv SWARM.sif python3 /opt/SWARM/SWARM_scripts/SWARM_read_level.py --OPTIONS

#The image was built using singularity v3.11.0 and GO v1.18.2

#The image was tested on x86-64 Linux systems running CentOS and Ubuntu, and with NVIDIA Volta, Hopper, and Blackwell GPUs
```

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

