using Printf
using ROCAnalysis
using HypothesisTests
using MultipleTesting

include("kendalls_tau_b.jl")

folder = "slim_simu/m=0.04/map3/"
min_f = 1e-4
sep = 10001

gene_length = 5000
gene_per_chrom = 200
num_chrom = 5

function pr_curve(scores::Vector{Float64}, labels::Vector{Bool})
  n = length(scores)
  @assert length(labels) == n

  P = count(labels)
  if P == 0
    error("No positive samples in labels")
  end

  ord = sortperm(scores, rev=true)
  tp = 0
  fp = 0

  recall = Float64[]
  precision = Float64[]

  push!(recall, 0.0)
  push!(precision, 1.0)

  for k in 1:n
    i = ord[k]
    if labels[i]
      tp += 1
    else
      fp += 1
    end
    
    r = tp / P
    p = tp / (tp + fp)
      
    push!(recall, r)
    push!(precision, p)
  end

  return recall, precision
end

function auprc(recall::Vector{Float64}, precision::Vector{Float64})
  @assert length(recall) == length(precision)
  a = 0.0
  for i in 1:length(recall)-1
    dx = recall[i+1] - recall[i]
    if dx > 0
      a += dx * precision[i+1]
    end
  end
  return a
end

function step_eval_pr!(recall, precision, grid, tar_list)
  j = 1
  for (i, g) in pairs(grid)
    while j < length(recall) && recall[j+1] <= g
      j += 1
    end
    tar_list[i] += precision[j]
  end
end

function gene_level_val(val_list, ret_idx, len_idx)
  ret_list = fill(minimum(val_list) - 1.0, len_idx)

  for i in 1:length(val_list)
    if val_list[i] > ret_list[ret_idx[i]]
      ret_list[ret_idx[i]] = val_list[i]
    end
  end

  return ret_list
end

for reso in 0:4
  all_case = 0

  x_latent = LinRange(0.0, 1.0, sep)
  kendall = fill(0.0, sep)
  mwU = fill(0.0, sep)
  baypass = fill(0.0, sep)
  matrix = fill(0.0, sep)
  lfmm = fill(0.0, sep)

  kendall_auc = Float64[]
  mwU_auc = Float64[]
  baypass_auc = Float64[]
  matrix_auc = Float64[]
  lfmm_auc = Float64[]

  for i in 1:50
    chrom = fill(-1, gene_per_chrom * num_chrom)
    pos = fill(-1, gene_per_chrom * num_chrom)
    ret_idx = Int[]
    true_false = fill(false, gene_per_chrom * num_chrom)

    path = folder * "data$(i)/answers$(reso).txt"

    if !isfile(path)
      println("$(i) $(reso)")
    else
      open(path, "r") do f
        while !eof(f)
          line = strip(readline(f))
          parts = split(line)

          now_chrom = parse(Int, parts[1])
          now_pos = parse(Int, parts[2])

          if div(now_pos, gene_length) == gene_per_chrom
            error("invalid position")
          end

          idx = div(now_pos, gene_length) + (now_chrom - 1) * gene_per_chrom + 1
          push!(ret_idx, idx)

          chrom[idx] = now_chrom
          pos[idx] = now_pos

          if abs(parse(Float64, parts[3])) > 1e-8
            true_false[idx] = true
          end
        end
      end

      all_case += 1

      env = Float64[]
      open(folder * "data$(i)/envfile$(reso).env", "r") do f
        while !eof(f)
          line = strip(readline(f))
          push!(env, parse(Float64, line))
        end
      end

      # kendall's tau
      p_val = Float64[]
      open(folder * "data$(i)/genofile$(reso).geno", "r") do f
        while !eof(f)
          line = strip(readline(f))
          parts = collect(line)
          ns = [parse(Int, parts[j]) for j in 1:length(parts)]

          if length(ns) != length(env)
            error("ns length does not match with env length")
          end

          test = kendall_tau_b(ns, env)
          push!(p_val, test[3])
        end
      end

      scores = gene_level_val(-p_val, ret_idx, gene_per_chrom * num_chrom)
      labels = true_false

      recall, precision = pr_curve(scores, labels)
      step_eval_pr!(recall, precision, x_latent, kendall)
      push!(kendall_auc, auprc(recall, precision))

      # Mann-Whitney U test
      p_val = Float64[]
      open(folder * "data$(i)/genofile$(reso).geno", "r") do f
        while !eof(f)
          line = strip(readline(f))
          parts = collect(line)

          group0 = Float64[]
          group1 = Float64[]

          for j in 1:length(parts)
            ac = parse(Int, parts[j])

            if ac == 0
              push!(group0, env[j])
              push!(group0, env[j])
            elseif ac == 1
              push!(group0, env[j])
              push!(group1, env[j])
            elseif ac == 2
              push!(group1, env[j])
              push!(group1, env[j])
            else
              error("Invalid allele count")
            end
          end

          test = MannWhitneyUTest(group0, group1)
          push!(p_val, pvalue(test))
        end
      end

      scores = gene_level_val(-p_val, ret_idx, gene_per_chrom * num_chrom)
      labels = true_false

      recall, precision = pr_curve(scores, labels)
      step_eval_pr!(recall, precision, x_latent, mwU)
      push!(mwU_auc, auprc(recall, precision))


      # BAYPASS
      bf = Float64[]

      open(folder * "data$(i)/baypass/sample$(reso)/summary_betai_reg.out", "r") do f
        linenum = 0

        while !eof(f)
          line = strip(readline(f))

          if linenum > 0
            parts = split(line)
            push!(bf, 10^(parse(Float64, parts[7]) / 10.0))
          end

          linenum += 1
        end
      end

      scores = gene_level_val(bf, ret_idx, gene_per_chrom * num_chrom)
      labels = true_false

      recall, precision = pr_curve(scores, labels)
      step_eval_pr!(recall, precision, x_latent, baypass)
      push!(baypass_auc, auprc(recall, precision))


      # matrix
      p_val = Float64[]

      open(folder * "data$(i)/matrix/sample$(reso)/p_val.txt", "r") do f
        while !eof(f)
          line = strip(readline(f))
          push!(p_val, parse(Float64, line))
        end
      end

      scores = gene_level_val(-p_val, ret_idx, gene_per_chrom * num_chrom)
      labels = true_false

      recall, precision = pr_curve(scores, labels)
      step_eval_pr!(recall, precision, x_latent, matrix)
      push!(matrix_auc, auprc(recall, precision))

      # LFMM
      p_val = Float64[]

      open(folder * "data$(i)/LFMM/sample$(reso)/p-values.txt", "r") do f
        while !eof(f)
          line = strip(readline(f))
          push!(p_val, parse(Float64, line))
        end
      end

      scores = gene_level_val(-p_val, ret_idx, gene_per_chrom * num_chrom)
      labels = true_false

      recall, precision = pr_curve(scores, labels)
      step_eval_pr!(recall, precision, x_latent, lfmm)
      push!(lfmm_auc, auprc(recall, precision))
    end
  end

  kendall ./= all_case
  mwU ./= all_case
  baypass ./= all_case
  matrix ./= all_case
  lfmm ./= all_case

  open(folder * "/pr_data$(reso)_gene.txt", "w") do f
    println(f, "recall\tkendall\tmwU\tbaypass\tmatrix\tlfmm")

    for i in 1:sep
      @printf(f, "%.5f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", x_latent[i], kendall[i],
        mwU[i], baypass[i], matrix[i], lfmm[i])
    end
  end

  open(folder * "/aucpr_data$(reso)_gene.txt", "w") do f
    println(f, "kendall\tmwU\tbaypass\tmatrix\tlfmm")

    for i in 1:length(kendall_auc)
      @printf(f, "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", kendall_auc[i],
        mwU_auc[i], baypass_auc[i], matrix_auc[i], lfmm_auc[i])
    end
  end
end
