# Install SWARM

## Download the code and models

Simply clone from github (install lfs to download large h5 files):

```bash
git lfs install
git clone https://github.com/comprna/SWARM/ && cd SWARM 
```

If git lfs cannot be installed, download the models from the dropbox link:

```bash
git clone https://github.com/comprna/SWARM/ && cd SWARM 
rm -rf SWARM_models
wget 'https://zenodo.org/records/22123294/files/SWARM_models.tar.gz' -O SWARM_models.tgz
tar -xzf SWARM_models.tgz && rm -f SWARM_models.tgz
```

## Compile SWARM preprocessing

```bash
cd SWARM_scripts/preprocess/
#build and compile htslib, slow5tools, SWARM_preprocess
bash build.sh
```

## Dependencies

SWARM supports GPU inference with tensorflow, tested with versions 2.8.0 and 2.15.0

GPU-configured tensorflow should be available on most HPC systems. Otherwise, you can install tensorflow configured for GPU as per https://www.tensorflow.org/install/

python requirements:

```bash
python==3.11.7
tensorflow==2.15.0
numpy==1.26.2
pandas==2.2.0
scikit-learn==1.4.0
pysam==0.22.1
scipy==1.14.1
statsmodels==0.14.4
```

Example for setting up the SWARM environment with conda:

```bash
conda create -n SWARM python==3.11.7 numpy==1.26.2 pandas==2.2.0 scikit-learn==1.4.0 pysam==0.22.1 scipy==1.14.1 statsmodels==0.14.4
conda activate SWARM
```

##File tree

```
└── SWARM
    ├── README.md
    ├── SWARM_models
    │   ├── kmer_model
    │   │   ├── model_5-mer.RNA002.csv
    │   │   └── model_5-mer.RNA004.csv
    │   ├── Model1
    │   │   ├── RNA002
    │   │   │   ├── m5C
    │   │   │   │   └── Model_100_epoch_relu.h5
    │   │   │   ├── m6A
    │   │   │   │   └── Model_100_epoch_relu.h5
    │   │   │   └── pU
    │   │   │       └── Model_100_epoch_relu.h5
    │   │   └── RNA004
    │   │       ├── m5C
    │   │       │   └── Model_100_epoch_relu.h5
    │   │       ├── m6A
    │   │       │   └── Model_100_epoch_relu.h5
    │   │       └── pU
    │   │           └── Model_100_epoch_relu.h5
    │   └── Model2
    │       ├── RNA002
    │       │   ├── m5C
    │       │   │   └── Model_100_epoch_relu.h5
    │       │   ├── m6A
    │       │   │   └── Model_100_epoch_relu.h5
    │       │   └── pU
    │       │       └── Model_100_epoch_relu.h5
    │       └── RNA004
    │           ├── m5C
    │           │   └── Model_100_epoch_relu.h5
    │           ├── m6A
    │           │   └── Model_100_epoch_relu.h5
    │           └── pU
    │               └── Model_100_epoch_relu.h5
    └── SWARM_scripts
        ├── predict
        │   ├── DL_models.py
        │   ├── network_21122023.py
        │   ├── network_2132024.py
        │   ├── network_27082022.py
        │   ├── predict_model1_from_pickle.py
        │   ├── predict_model1_parallel_modbam.py
        │   └── predict_model1_parallel.py
        ├── preprocess
        │   ├── argagg.hpp
        │   ├── build.sh
        │   ├── check_RNA_kit.cpp
        │   ├── Makefile
        │   ├── split_bams.py
        │   ├── SWARM_preprocess.cpp
        │   ├── SWARM_preprocess.py
        │   ├── SWARM_preprocess_target_9mers.cpp
        │   └── SWARM_preprocess_targets.cpp
        ├── process_modbam.py
        ├── SWARM_diff.py
        ├── SWARM_read_level.py
        ├── SWARM_site_level.py
        └── train_models
            ├── assemble_data.py
            ├── network_27082022.py
            ├── split_training_by_9mers.py
            ├── train_model1.py
            └── trim_tsv_events.py
```
