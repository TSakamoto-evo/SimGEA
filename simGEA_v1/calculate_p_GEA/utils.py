import numpy as np
from numba import njit
import math

def read_bin_file(bin_file):
  ret = []

  with open(bin_file, mode="r") as f:
    for line in f:
      line = line.rstrip("\n")
      ret.append(float(line))

  ret = sorted(ret)
  return ret

def read_env_file(env_file):
  ret = []

  with open(env_file, mode="r") as f:
    for line in f:
      line = line.rstrip("\n")
      tmplist = line.split()
      ret.append([float(i) for i in tmplist])
  
  return ret

def calculate_null_corr(null_file, bin_val, env_list, use_continuity=False):
  ret = [[[] for i in range(len(bin_val))] for j in range(len(env_list))]
  max_z_list = []

  ranks_list = []
  tie_term_list = []

  with open(null_file, mode="r") as f:
    for line in f:
      line = line.rstrip("\n")
      tmplist = line.split()

      count1 = np.array([int(tmplist[i]) for i in range(1, len(tmplist), 2)], dtype=float)
      count0 = np.array([int(tmplist[i]) for i in range(2, len(tmplist), 2)], dtype=float)

      if count0.size != count1.size:
        raise ValueError("Allele count list has size of odd number")

      weights = np.asarray(count1 + count0, dtype=float)
      break

  for i in range(len(env_list)):
    env = np.array(env_list[i])

    if weights.size != env.size:
      raise ValueError("Number of populations is different between nullfile and envfile")
    
    if np.all(env == env[0]):
      raise ValueError("No variation in environment " + str(i))

    order, starts, ends = prep_tie_groups(env)
    ranks, tie_term = midrank_from_weights(order, starts, ends, weights)
    ranks = np.asarray(ranks, dtype=float)

    ranks_list.append(ranks)
    tie_term_list.append(tie_term)
    max_z_list.append(max_u_val(weights, ranks, tie_term, bin_val, use_continuity))

  with open(null_file, mode="r") as f:
    for line in f:
      line = line.rstrip("\n")
      tmplist = line.split()

      bin_idx = int(tmplist[0])

      if bin_idx < 0 or bin_idx >= len(bin_val):
        print("Out of range of bin index in nullfile")
        continue

      count1 = np.array([int(tmplist[i]) for i in range(1, len(tmplist), 2)], dtype=float)
      count0 = np.array([int(tmplist[i]) for i in range(2, len(tmplist), 2)], dtype=float)

      if not np.array_equal(count1 + count0, weights):
        raise ValueError("The total number of samples are different")
      elif np.sum(count0) == 0 or np.sum(count1) == 0:
        print("There is monomorphic site in nullfile")
        continue

      for i in range(len(env_list)):
        U, N1, N0 = u_from_ranks(count1, count0, ranks_list[i])
        varU = varU_tie_corrected(N1, N0, tie_term_list[i])
        z_val = z_from_U(U, N1, N0, varU, use_continuity)

        ret[i][bin_idx].append(z_val)
        ret[i][bin_idx].append(-z_val)

  return ret, weights, ranks_list, tie_term_list, max_z_list

def prep_tie_groups(env):
  n = env.size

  order = np.argsort(env, kind="stable")
  xs = env[order]

  starts = []
  ends = []
  start = 0

  for i in range(1, n):
    if xs[i] != xs[i - 1]:
      starts.append(start)
      ends.append(i)
      start = i 

  starts.append(start)
  ends.append(n)
  
  return order, np.array(starts, dtype=int), np.array(ends, dtype=int)

def midrank_from_weights(order, starts, ends, weights):
  ws = weights[order]
  cw = np.cumsum(ws)

  r_sorted = np.empty_like(ws)

  for s, e in zip(starts, ends):
    left = (cw[s] - ws[s]) + 1.0
    right = cw[e - 1]

    r_sorted[s:e] = 0.5 * (left + right)

  ranks = np.empty_like(r_sorted)
  ranks[order] = r_sorted

  t = np.add.reduceat(ws, starts)
  N = weights.sum()
  tie_term = (np.sum(t**3 - t)) / (N * (N - 1))

  return ranks, tie_term


def max_u_val(weights, ranks, tie_term, bin_val, use_continuity=False):
  wsum = int(np.sum(weights))
  rank_order = np.argsort(ranks, kind="stable")
  ret = [0.0 for i in range(len(bin_val))]

  for i in range(1, wsum):
    p = 1.0 * i / wsum
    if p > 0.5:
      p = 1.0 - p
    
    if p > bin_val[0]:
      bin_idx = 0
      while bin_idx < len(bin_val) - 1 and p > bin_val[bin_idx + 1]:
        bin_idx += 1

      c1 = np.zeros_like(weights)
      c0 = np.copy(weights)
      remaining = i

      for j in range(len(weights)):
        add = min(weights[rank_order[j]], remaining)
        remaining -= add
        c1[rank_order[j]] = add
        c0[rank_order[j]] -= add

        if remaining == 0:
          U, N1, N0 = u_from_ranks(c1, c0, ranks)
          varU = varU_tie_corrected(N1, N0, tie_term)
          z = z_from_U(U, N1, N0, varU, use_continuity)

          if np.abs(z) > ret[bin_idx]:
            ret[bin_idx] = np.abs(z)
          
          break
  return ret

@njit(cache=True, fastmath=False)
def u_from_ranks(counts1, counts0, ranks):
  c1 = counts1
  c0 = counts0
  r  = ranks

  N1 = c1.sum()
  N0 = c0.sum()

  R1 = np.dot(c1, r)
  U = R1 - N1 * (N1 + 1.0) / 2.0

  return U, N1, N0

@njit(cache=True, fastmath=False)
def varU_tie_corrected(N1, N0, tie_term):
  varU = (N0 * N1 / 12.0) * (N1 + N0 + 1.0 - tie_term)

  return(float(varU))

@njit(cache=True, fastmath=False)
def z_from_U(U, N1, N0, varU, use_continuity=False):
  if varU <= 0:
    return np.nan
  
  meanU = N1 * N0 / 2.0
  num = U - meanU

  if use_continuity:
    cc = 0.5
    num = np.sign(num) * np.maximum((abs(num) - cc), 0.0)

  z = num / np.sqrt(varU)

  return z

@njit(cache=True, fastmath=False)
def calc_kernel_pval(z_focal, z_list):
  s = 0.0
  n = z_list.size

  for k in range(n):
    s += math.erfc(z_focal - z_list[k])
  return(s / n)