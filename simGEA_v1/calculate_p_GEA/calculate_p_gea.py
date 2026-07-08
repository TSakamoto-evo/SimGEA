import os
import argparse
import time

# Optional: uncomment these lines to limit each process to one thread
# os.environ["OPENBLAS_NUM_THREADS"] = "1"
# os.environ["OMP_NUM_THREADS"] = "1"

from scipy import stats
import numpy as np
from scipy.stats import norm
import copy
import sys

import utils

np.set_printoptions(precision=8, suppress=False)

def pct(x: str) -> float:
  v = float(x)
  if not (0.5 < v < 1.0):
      raise argparse.ArgumentTypeError("Specify q within (0.5, 1)")
  return v

parser = argparse.ArgumentParser(
  description="Calculate p-value using null distribution (Version 1.0.0).",
  formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument("-b", "--bin-file", required=True, metavar="<bin file>", help="bin file")
parser.add_argument("-e", "--env-file", required=True, metavar="<env file>", help="environment file")
parser.add_argument("-n", "--null-file", required=True, metavar="<null file>", help="neutral loci used for null distribution")
parser.add_argument("-d", "--data-file", required=True, metavar="<data file>", help="data for p-value calculation")

parser.add_argument("-l", "--log-file", metavar="<log file>", default="log.txt", help="log file")
parser.add_argument("-o", "--out-file", default="p_val.txt", metavar="<output file>", help="output file")
parser.add_argument("-q", "--tail-quantile", default=0.995, type=pct, metavar="<tail quantile>", help="definition of the tail")

args = parser.parse_args()

########## PARAMETERS ##########
q_tail = args.tail_quantile

bin_file   = args.bin_file
env_file   = args.env_file
null_file  = args.null_file
data_file  = args.data_file
output_file = args.out_file
log = open(args.log_file, mode="a")
################################

print("P-VALUE CALCULATION\n", file=log, flush=True)

start = time.time()

print("Start reading nullfile", file=log, flush=True)

bin_val = utils.read_bin_file(bin_file)
env_list = utils.read_env_file(env_file)
z_list, weights, ranks_list, tie_term_list, max_z_list = utils.calculate_null_corr(null_file, bin_val, env_list)

print("Start calculating smoothed distribution", file=log, flush=True)

for i in range(len(env_list)):
  for j in range(len(bin_val)):
    num_outlier = len(z_list[i][j]) * (1.0 - q_tail)

    if num_outlier < 250:
      print("Warning: Tail (q > " + str(q_tail) + ") has only " + str(int(num_outlier)) + "samples", file=log, flush=True)
    
    if len(z_list[i][j]) == 0:
      print("Error: No samples for Env " + str(i) + " , Bin " + str(j), file=log, flush=True)
      print("Error: No samples for Env " + str(i) + " , Bin " + str(j))
      sys.exit(1)


band_width = [[0.0 for j in range(len(bin_val))] for i in range(len(env_list))]
u_thre = [[0.0 for j in range(len(bin_val))] for i in range(len(env_list))]
cf_thre = [[0.0 for j in range(len(bin_val))] for i in range(len(env_list))]
c_para = [[0.0 for j in range(len(bin_val))] for i in range(len(env_list))]
s_para = [[0.0 for j in range(len(bin_val))] for i in range(len(env_list))]

print("fitted parameters for extreme z_values", file=log)
print("env bin z_thre cf_thre c scale max_z", file=log)
for i in range(len(env_list)):
  for j in range(len(bin_val)):
    arr = np.asarray(z_list[i][j], dtype=float)

    q1 = np.percentile(arr, 25)
    q3 = np.percentile(arr, 75)
    iqr = q3 - q1

    sig_est = np.min([np.std(arr, ddof=1), iqr/1.34])
    h = 0.9 * sig_est * (arr.size)**(-0.2)
    band_width[i][j] = h

    u_thre[i][j] = np.percentile(arr, 100 * q_tail)
    
    tmpvals = (u_thre[i][j] - arr) / band_width[i][j]
    cf_thre[i][j] = min(max(np.mean(norm.sf(tmpvals)), 0.0), 1.0)

    y = arr[arr > u_thre[i][j]] - u_thre[i][j]
    c, loc, scale = stats.genpareto.fit(y, floc=0.0)

    if c < 0 and max_z_list[i][j] - u_thre[i][j] > -scale / c:
      cs = np.linspace(c, 0.0, 100)  

      best = None
      best_ll = -np.inf

      for c_test in cs:
        scale_test = stats.genpareto.fit(y, floc=0, fc=c_test)[2]
        
        if c_test >= 0.0 or -scale_test / c_test >= max_z_list[i][j] - u_thre[i][j]:
          ll = np.sum(stats.genpareto.logpdf(y, c_test, loc=0, scale=scale_test))
          if ll > best_ll:
            best_ll = ll
            best = (c_test, scale_test)
      
      c = best[0]
      scale = best[1]

    c_para[i][j] = c
    s_para[i][j] = scale

    print(i, j, u_thre[i][j], cf_thre[i][j], c, scale, max_z_list[i][j], file=log)
    
print("Start calculating p-value", file=log, flush=True)

num_out_of_range = 0
out_p_val = open(output_file, mode="w")
total_sample = np.sum(weights)

SQRT2 = np.sqrt(2.0)
inv_h_s2 = [[1.0 / (band_width[i][j] * SQRT2) for j in range(len(bin_val))] for i in range(len(env_list))]
z_scaled = copy.deepcopy(z_list)

for i in range(len(env_list)):
  for j in range(len(bin_val)):
    z_scaled[i][j] = np.array(z_list[i][j], dtype=float) * inv_h_s2[i][j]

with open(data_file, mode="r") as f:
  for line in f:
    line = line.rstrip("\n")
    tmplist = line.split()

    count1 = np.array([int(tmplist[i]) for i in range(0, len(tmplist), 2)], dtype=float)
    count0 = np.array([int(tmplist[i]) for i in range(1, len(tmplist), 2)], dtype=float)

    if not np.array_equal(count1 + count0, weights):
      raise ValueError("The total number of samples are different")

    freq = np.sum(count1) / total_sample

    if freq > 0.5:
      freq = 1.0 - freq

    if freq > bin_val[0] and (np.sum(count0) > 0 and np.sum(count1) > 0):
      bin_idx = 0
      while bin_idx < len(bin_val) - 1 and freq > bin_val[bin_idx + 1]:
        bin_idx += 1

      tmp_pval = []
      for i in range(len(env_list)):
        U, N1, N0 = utils.u_from_ranks(count1, count0, ranks_list[i])
        varU = utils.varU_tie_corrected(N1, N0, tie_term_list[i])
        z_val = utils.z_from_U(U, N1, N0, varU)

        if np.abs(z_val) <= u_thre[i][bin_idx]:
          z_val_scaled = np.abs(z_val) * inv_h_s2[i][bin_idx]
          p_val = utils.calc_kernel_pval(z_val_scaled, z_scaled[i][bin_idx])
          tmp_pval.append(f"{p_val:.6e}")
        else:
          pval = 2.0 * cf_thre[i][bin_idx] * stats.genpareto.sf(np.abs(z_val) - u_thre[i][bin_idx],
                                                c_para[i][bin_idx], loc=0.0, scale=s_para[i][bin_idx])
          tmp_pval.append(f"{pval:.6e}")
        
      print(' '.join(tmp_pval), file=out_p_val)

    else:
      num_out_of_range += 1
      tmp_pval = []
      for i in range(len(env_list)):
        tmp_pval.append(str(np.nan))
      print(' '.join(tmp_pval), file=out_p_val)

out_p_val.close()

print("Finish run\n", flush=True, file=log)

if num_out_of_range >= 1:
  print(str(num_out_of_range) + " sites are out of bin range\n", file=log, flush=True)

end = time.time()

secs = int(end - start)
days = int(secs / (24 * 60 * 60))
secs -= days * 24 * 60 * 60
hours = int(secs / (60 * 60))
secs -= hours * 60 * 60
minutes = int(secs / 60)
secs -= minutes * 60

print("Total execution time: " + str(days) + "d " + str(hours) + "h " +str(minutes) + "m " + str(secs) + "s", file=log, flush=True)

print("END P-VALUE CALCULATION\n\n", file=log, flush=True)

log.close()