#!/bin/bash
  
### download and compile htslib
echo -e "downloading htslib:\n"
wget https://github.com/samtools/htslib/releases/download/1.19/htslib-1.19.tar.bz2
tar -xf htslib-1.19.tar.bz2
rm -f htslib-1.19.tar.bz2

echo -e "\ncompiling htslib:\n"
cd htslib-1.19/ && make

### download and compile slow5lib
echo -e "\ndownloading slow5lib\n"
cd ../
git clone https://github.com/hasindu2008/slow5lib

echo -e "\ncompiling slow5lib\n"
cd slow5lib/ && make

### compile SWARM_preprocess
cd ../
echo -e "\ncompiling SWARM_preprocess\n"
make && echo -e "\nSWARM_preprocess is ready."
