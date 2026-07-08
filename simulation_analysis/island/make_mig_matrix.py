import numpy as np
import subprocess

folder = "100pop_100_L100000_a8"
size = 100
numset = 20

exp_total = 100

for rep in range(numset):
  mig_mat = np.zeros((size, size))

  for i in range(size):
    for j in range(i + 1, size):
      mig_mat[i, j] = np.random.exponential(exp_total / size)
      mig_mat[j, i] = mig_mat[i, j]

  for i in range(size):
    sum_val = 0.0
    for j in range(size):
      if i != j:
        sum_val += mig_mat[i, j]

    mig_mat[i, i] = -sum_val

  with open("mig_matrix.txt", mode="w") as f:
    for i in range(size):
      out_list = [str(j) for j in mig_mat[i]]
      print(" ".join(out_list), file=f)
  
  subprocess.run("mkdir " + folder + "/data" + str(rep + 1), shell=True)
  subprocess.run("cp mig_matrix.txt " + folder + "/data" + str(rep + 1), shell=True)

