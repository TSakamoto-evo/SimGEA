import subprocess, msprime, pyslim, tskit
import numpy as np
import random, math, statistics
import warnings

warnings.simplefilter('ignore', msprime.TimeUnitsMismatchWarning)

##################### set simulation parameters #####################
# diploid population size per deme
deme_size = 100
# migration rate between adjacent subpopulations
mig_rate = 0.01
# recombination rate within gene
r1 = 1e-7
# recombination rate between adjacent genes
r2 = 0.005
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

# run slim simulation
subprocess.run("slim" + " -m" + " -d" + " optimaFile=\\\"" + "env_paras.txt" + "\\\"" +
               " -d" + " N=" + str(deme_size) +
               " -d" + " m=" + str(mig_rate) +
               " -d" + " r1=" + str(r1) +
               " -d" + " r2=" + str(r2) +
               " -d" + " us=" + str(qtl_mut_rate) +
               " -d" + " gene_length=" + str(gene_len) +
               " simulator.slim", shell=True
)


for reso in range(5):
  ts = tskit.load("slim_local_adaptation.trees")

  pop_num = len(env)
  L = int(np.sqrt(pop_num))

  keep_nodes = []
  sampled_pop_list = []
  sample_ind_num_list = []

  with open("../sampled_pop_list" + str(reso) + ".txt", mode="r") as f:
    for line in f:
      line = line.rstrip("\n")
      tmplist = line.split()

      sampled_pop_list.append(int(tmplist[0]))
      sample_ind_num_list.append(int(tmplist[1]))

  for i in range(len(sampled_pop_list)):
    alive_inds = pyslim.individuals_alive_at(ts, 0, population=sampled_pop_list[i])
    keep_indivs = np.random.choice(alive_inds, sample_ind_num_list[i], replace=False)

    for j in keep_indivs:
      keep_nodes.extend(ts.individual(j).nodes)

  ts = ts.simplify(keep_nodes, keep_input_roots=True, filter_populations = False)

  print(f"Number of trees with only one root: {sum([t.num_roots == 1 for t in ts.trees()])}\n"
        f"Number with more than one root: {sum([t.num_roots > 1 for t in ts.trees()])}")

  demography = msprime.Demography.from_tree_sequence(ts)
  for pop in demography.populations:
    pop.initial_size = deme_size

  for i in range(L):
    for j in range(L):
      index = i * L + j + 1
      if i > 0:
        index1 = (i - 1) * L + j + 1
        demography.add_migration_rate_change(time=ts.metadata['SLiM']['tick'],
                rate=mig_rate, source=index, dest=index1)
      if i < (L - 1):
        index1 = (i + 1) * L + j + 1
        demography.add_migration_rate_change(time=ts.metadata['SLiM']['tick'],
                rate=mig_rate, source=index, dest=index1)
      if j > 0:
        index1 = i * L + (j - 1) + 1
        demography.add_migration_rate_change(time=ts.metadata['SLiM']['tick'],
                rate=mig_rate, source=index, dest=index1)
      if j < (L - 1):
        index1 = i * L + (j + 1) + 1
        demography.add_migration_rate_change(time=ts.metadata['SLiM']['tick'],
                rate=mig_rate, source=index, dest=index1)

  reco_pos = [0]
  reco_rate = []

  for i in range(200):
    reco_pos.append((i + 1) * gene_len)
    reco_rate.append(r1)

    if i != 199:
      reco_pos.append((i + 1) * gene_len + 1)
      reco_rate.append(r2)

  recomb_map = msprime.RateMap(position=reco_pos, rate=reco_rate)

  recap0 = pyslim.recapitate(ts, recombination_rate=recomb_map,
    demography=demography)

  print(f"Number of trees with only one root: {sum([t.num_roots == 1 for t in recap0.trees()])}\n"
        f"Number with more than one root: {sum([t.num_roots > 1 for t in recap0.trees()])}")

  recap0 = recap0.simplify(filter_populations = False)

  recap = msprime.sim_mutations(recap0, rate=neutral_mut_rate,
      model=msprime.SLiMMutationModel(type=0), keep = False)
  
  indivlist = []
  initial_populations = open("initial_populations" + str(reso) + ".txt", mode="w")
  for i in pyslim.individuals_alive_at(recap, 0):
    ind = recap.individual(i)
    if recap.node(ind.nodes[0]).is_sample():
      indivlist.append(i)
      pop = recap.node(ind.nodes[0]).population
      print(pop, file=initial_populations)
      print(pop, file=initial_populations)

      # if one node is a sample, the other should be also:
      assert recap.node(ind.nodes[1]).is_sample()

  with open("output_sele" + str(reso) + ".vcf", "w") as vcffile:
    recap0.write_vcf(vcffile, individuals=indivlist, allow_position_zero=True)
  with open("output_neu" + str(reso) + ".vcf", "w") as vcffile:
    recap.write_vcf(vcffile, individuals=indivlist, allow_position_zero=True)

print("vcf has been produced.")
exit()

