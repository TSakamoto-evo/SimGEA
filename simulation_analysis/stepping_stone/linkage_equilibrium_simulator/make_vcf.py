import subprocess
import numpy as np
import random, math, statistics
import warnings

##################### set simulation parameters #####################
# diploid population size per deme
deme_size = 100
# migration rate between adjacent subpopulations
mig_rate = 0.01
# mutation rate at the QTL
qtl_mut_rate = 1e-7
# gene length
gene_len = 5000

# mutation rate at neutral sites
neutral_mut_rate = 1e-8
##################### end set simulation parameters #####################

env = []
with open("env_file.txt", mode="r") as f:
  for line in f:
    line = line.rstrip("\n")
    env.append(float(line))

env_mean = statistics.mean(env)
env_sd = statistics.stdev(env)

with open("env_paras.txt", mode="w") as f:
  for i in env:
    print((i - env_mean) / env_sd, file=f)

subprocess.run("./unlinked_simu.out env_paras.txt " + str(deme_size) + " " + str(mig_rate) + " " +
  str(qtl_mut_rate) + " " + str(neutral_mut_rate) + " " + str(gene_len), shell=True)

for reso in range(5):
  sampled_pop_list = []
  sample_ind_num_list = []

  with open("../sampled_pop_list" + str(reso) + ".txt", mode="r") as f:
    for line in f:
      line = line.rstrip("\n")
      tmplist = line.split()

      sampled_pop_list.append(int(tmplist[0]))
      sample_ind_num_list.append(int(tmplist[1]))

  initial_populations = open("initial_populations" + str(reso) + ".txt", mode="w")
  for i in range(len(sample_ind_num_list)):
    for j in range(2 * sample_ind_num_list[i]):
      print(sampled_pop_list[i], file=initial_populations)
  initial_populations.close()

  selected_vcf = open("output_sele" + str(reso) + ".vcf", mode="w")
  print("#pseudoVCF code", file=selected_vcf)
  print("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sep="\t", end="", file=selected_vcf)

  for i in range(sum(sample_ind_num_list)):
    print("\t", "tsk_" + str(i), sep="", end="", file=selected_vcf)
  print("", file=selected_vcf)

  with open("mean_freq_effect.txt", mode="r") as f:
    for line in f:
      line = line.rstrip("\n")
      tmplist = line.split("\t")

      retlist = [str(tmplist[0]), str(tmplist[1]), "", "", "", ".", "PASS", ".", "GT"]
      ac = 0

      for i in range(len(sample_ind_num_list)):
        for j in range(sample_ind_num_list[i]):
          pop = sampled_pop_list[i]
          freq = float(tmplist[3 + pop - 1])

          allele0 = int(random.random() < freq)
          allele1 = int(random.random() < freq)

          retlist.append(str(allele0) + "|" + str(allele1))
          ac += allele0 + allele1

      if ac != 0 and ac != 2 * sum(sample_ind_num_list):
        print("\t".join(retlist), file=selected_vcf)
  selected_vcf.close()


  neutral_vcf = open("output_neu" + str(reso) + ".vcf", mode="w")
  print("#pseudoVCF code", file=neutral_vcf)
  print("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sep="\t", end="", file=neutral_vcf)

  for i in range(sum(sample_ind_num_list)):
    print("\t", "tsk_" + str(i), sep="", end="", file=neutral_vcf)
  print("", file=neutral_vcf)

  with open("mean_freq_neutral.txt", mode="r") as f:
    for line in f:
      line = line.rstrip("\n")
      tmplist = line.split("\t")

      retlist = [str(tmplist[0]), str(tmplist[1]), "", "", "", ".", "PASS", ".", "GT"]
      ac = 0

      for i in range(len(sample_ind_num_list)):
        for j in range(sample_ind_num_list[i]):
          pop = sampled_pop_list[i]
          freq = float(tmplist[2 + pop - 1])

          allele0 = int(random.random() < freq)
          allele1 = int(random.random() < freq)

          retlist.append(str(allele0) + "|" + str(allele1))
          ac += allele0 + allele1

      if ac != 0 and ac != 2 * sum(sample_ind_num_list):
        print("\t".join(retlist), file=neutral_vcf)
  neutral_vcf.close()

print("vcf has been produced.")
exit()
