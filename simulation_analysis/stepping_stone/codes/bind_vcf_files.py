import re, statistics, os
import numpy as np

folder = "."

cutoff = 0.05
per_run = 5
sampling_methods = 5

for reso in range(sampling_methods):
  output_answer = open("answers" + str(reso) + ".txt", mode="w")
  output_vcf = open("bind_common_alleles" + str(reso) + ".vcf", mode="w")
  output_effect = open("effect_freq" + str(reso) + ".txt", mode="w")

  for folder_num in range(1, per_run + 1):
    tmp_sele_pos = []
    tmp_sele_row = []

    tmp_dict = {}

    with open(os.path.join(folder, "folder" + str(folder_num), "mean_freq_effect.txt"), mode="r") as f:
      for line in f:
        line = line.rstrip("\n")
        tmplist = line.split()

        if float(tmplist[2]) > 0.0:
          tmp_dict[int(tmplist[1])] = float(tmplist[2])
        
        tmplist[0] = str(folder_num)
        print("\t".join(tmplist), file=output_effect)

    with open(os.path.join(folder, "folder" + str(folder_num), "output_sele" + str(reso) + ".vcf"), mode="r") as f:
      for line in f:
        line = line.rstrip("\n")

        if not re.match(r'#', line):
          tmplist = line.split("\t")

          biallelic = 1
          allele0 = 0
          allele1 = 0

          for i in range(9, len(tmplist)):
            genotype = tmplist[i].split("|")

            if genotype[0] == "0":
              allele0 += 1
            elif genotype[0] == "1":
              allele1 += 1
            else:
              biallelic = 0
            
            if genotype[1] == "0":
              allele0 += 1
            elif genotype[1] == "1":
              allele1 += 1
            else:
              biallelic = 0
          
          if biallelic == 1 and allele0 / (allele0 + allele1) >= cutoff and allele1 / (allele0 + allele1) >= cutoff:
            tmp_sele_pos.append(int(tmplist[1]))

            tmplist[0] = str(folder_num)
            tmp_sele_row.append("\t".join(tmplist))

    next_sele_pos = 0
    with open(os.path.join(folder, "folder" + str(folder_num), "output_neu" + str(reso) + ".vcf"), mode="r") as f:
      for line in f:
        line = line.rstrip("\n")

        if re.match(r'#', line):
          if folder_num == 1:
            print(line, file=output_vcf)
        else:
          tmplist = line.split("\t")

          biallelic = 1
          allele0 = 0
          allele1 = 0

          for i in range(9, len(tmplist)):
            genotype = tmplist[i].split("|")

            if genotype[0] == "0":
              allele0 += 1
            elif genotype[0] == "1":
              allele1 += 1
            else:
              biallelic = 0
            
            if genotype[1] == "0":
              allele0 += 1
            elif genotype[1] == "1":
              allele1 += 1
            else:
              biallelic = 0
          
          if biallelic == 1 and allele0 / (allele0 + allele1) >= cutoff and allele1 / (allele0 + allele1) >= cutoff:
            while next_sele_pos < len(tmp_sele_pos) and tmp_sele_pos[next_sele_pos] < int(tmplist[1]):
              print(tmp_sele_row[next_sele_pos], file=output_vcf)
              print(folder_num, tmp_sele_pos[next_sele_pos], tmp_dict[tmp_sele_pos[next_sele_pos]], file=output_answer)
              next_sele_pos += 1
            
            tmplist[0] = str(folder_num)
            print("\t".join(tmplist), file=output_vcf)
            print(folder_num, tmplist[1], 0.0, file=output_answer)
    
    while next_sele_pos < len(tmp_sele_pos):
      print(tmp_sele_row[next_sele_pos], file=output_vcf)
      print(folder_num, tmp_sele_pos[next_sele_pos], tmp_dict[tmp_sele_pos[next_sele_pos]], file=output_answer)
      next_sele_pos += 1

  output_answer.close()
  output_vcf.close()
  output_effect.close()
        
exit()  