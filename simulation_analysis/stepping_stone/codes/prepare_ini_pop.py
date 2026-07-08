import numpy as np

pop_num = 14 * 14
sample_pop_num = [160, 80, 40, 20, 10]
ind_per_pop = [1, 2, 4, 8, 16]

for reso in range(len(sample_pop_num)):
  pop_list = list(range(1, pop_num + 1))
  sampled_pop_list = np.random.choice(pop_list, sample_pop_num[reso], replace=False)

  with open("sampled_pop_list" + str(reso) + ".txt", mode="w") as f:
    for i in sampled_pop_list:
      print(i, ind_per_pop[reso], sep=" ", file=f)

exit()