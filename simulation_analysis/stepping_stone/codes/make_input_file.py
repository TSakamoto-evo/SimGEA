import re, statistics, os
import numpy as np

folder = "."
sampling_methods = 5

for reso in range(sampling_methods):
  pop_dict = {}

  with open(os.path.join(folder, "sampled_pop_list" + str(reso) + ".txt"), mode="r") as f:
    for line in f:
      line = line.rstrip("\n")
      tmplist = line.split()

      if int(tmplist[0]) not in pop_dict.keys():
        pop_dict[int(tmplist[0])] = len(pop_dict)


  output_env = open("env_dist" + str(reso) + ".txt", mode="w")

  list_env = []
  with open(os.path.join(folder, "folder1/env_paras.txt"), mode="r") as f:
    for line in f:
      line = line.rstrip("\n")
      list_env.append(float(line))


  outlist = [0.0 for i in range(len(pop_dict))]
  for key, value in pop_dict.items():
    outlist[value] = list_env[key - 1]
  
  outlist2 = [str((outlist[i] - statistics.mean(outlist)) / statistics.stdev(outlist)) for i in range(len(outlist))]

  print(" ".join(outlist2), file=output_env)
  output_env.close()

  output_snp = open("snpsfile" + str(reso) + ".txt", mode="w")

  pop_list = []

  with open(os.path.join(folder, "folder1/initial_populations" + str(reso) + ".txt"), mode="r") as f:
    linenum = 0
    for line in f:
      line = line.rstrip("\n")

      if linenum % 2 == 0:
        pop_list.append(int(line))
      linenum += 1
  
  with open(os.path.join(folder, "bind_common_alleles" + str(reso) + ".vcf"), mode="r") as f:
    for line in f:
      line = line.rstrip("\n")

      if not re.match(r'#', line):
        tmplist = line.split("\t")

        list_allele = [0 for i in range(2 * len(pop_dict))]

        for i in range(9, len(tmplist)):
          genotype = tmplist[i].split("|")

          if genotype[0] == "0":
            list_allele[2 * pop_dict[pop_list[i - 9]]] += 1
          elif genotype[0] == "1":
            list_allele[2 * pop_dict[pop_list[i - 9]] + 1] += 1
            
          
          if genotype[1] == "0":
            list_allele[2 * pop_dict[pop_list[i - 9]]] += 1
          elif genotype[1] == "1":
            list_allele[2 * pop_dict[pop_list[i - 9]] + 1] += 1
        
        list_allele_str = [str(i) for i in list_allele]
        print(" ".join(list_allele_str), file=output_snp)
  output_snp.close()


  output_lfmm = open("genofile" + str(reso) + ".geno", mode="w")
  with open(os.path.join(folder, "bind_common_alleles" + str(reso) + ".vcf"), mode="r") as f:
    for line in f:
      line = line.rstrip("\n")

      if not re.match(r'#', line):
        tmplist = line.split("\t")
        ret = ""

        for i in range(9, len(tmplist)):
          genotype = tmplist[i].split("|")
          allele_count = int(genotype[0]) + int(genotype[1])
          ret = ret + str(allele_count)
        
        print(ret, file=output_lfmm)
  output_lfmm.close()

  output_env_lfmm = open("envfile" + str(reso) + ".env", mode="w")
  for i in pop_list:
    print(list_env[i - 1], file=output_env_lfmm)
  output_env_lfmm.close()

exit()  