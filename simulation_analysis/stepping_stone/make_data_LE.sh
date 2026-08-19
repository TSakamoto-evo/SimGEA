#!/bin/bash

# the number of chromosomes
chrom=5

mkdir data
rm -r data/*

# specify simulation parameters in linkage_equilibrium_simulator/make_vcf.py

cd linkage_equilibrium_simulator
g++ *.cpp -Wall -Wextra -std=c++17 -O3 -o unlinked_simu.out
mv unlinked_simu.out ../data
cp make_vcf.py ../data

cd ../data

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
  cp ../unlinked_simu.out ./
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
## SimGEA was run in mode 1 using the default parameters with 100000 sites per window

exit 0