#!/bin/bash

# the number of chromosomes
chrom=5

mkdir data
rm -r data/*

# specify simulation parameters in slim_code/make_vcf.py
cp slim_code/make_vcf.py ./data

# this is a SLiM code working under version 4
cp slim_code/simulator.slim ./data/simulator.slim

# use simulator_for_v5.slim for the current version of the SLiM
#cp slim_code/simulator_for_v5.slim ./data/simulator.slim

cd data

# specify sampling methods in the following python file
# this program randomly determines sampling locations
python3 ../codes/prepare_ini_pop.py

# specify an environmental map
envfile="env_file_map1.txt"

# run simulations
for((i=1; i<=${chrom}; i++)); do
  # simulate each chromosome
  # generate VCF-like files for both neutral and selected variants

  mkdir folder${i}
  cd folder${i}
  cp ../simulator.slim ./
  cp ../make_vcf.py ./
  cp ../../env_maps/${envfile} ./env_file.txt
  python3 make_vcf.py

  rm *.out
  rm *.py
  cd ../
done

# bind VCFs across chromosomes
# specify the number of chromosomes (per_run), MAF filter, and the number of 
# sampling strategies (sampling_method) in the python code
python3 ../codes/bind_vcf_files.py

# convert file format (BayPass and LFMM style)
python3 ../codes/make_input_file.py

cp folder1/env_paras.txt ./
rm *.py
rm -r folder*

samplingStrategies=5
for ((j=0; j<samplingStrategies; j++)); do
  gzip -9 bind_common_alleles${j}.vcf
done

## R code "codes/lea.R" was used for LFMM2 analysis
## "codes/bin.txt" was used to specify MAF bins in SimGEA

exit 0