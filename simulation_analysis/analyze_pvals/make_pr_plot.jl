using Printf
using ROCAnalysis
using HypothesisTests
using MultipleTesting

include("kendalls_tau_b.jl")

folder = "LE_simu/m=0.04/map3/"
min_f = 1e-4
sep = 10001

q_val_crit = 0.05
baypass_crit = 100.0

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

output_fdr = open(folder * "output_power.txt", "w")
println(output_fdr, "reso\tmethod\tmindex\tpower\tfdr")

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
    chrom = String[]
    pos = Int[]
    true_false = Bool[]

    path = folder * "data$(i)/answers$(reso).txt"

    if !isfile(path)
      println("$(i) $(reso)")
    else
      open(path, "r") do f
        while !eof(f)
          line = strip(readline(f))
          parts = split(line)

          push!(chrom, parts[1])
          push!(pos, parse(Int, parts[2]))
          push!(true_false, abs(parse(Float64, parts[3])) > 1e-8)
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

      scores = -p_val
      labels = true_false

      recall, precision = pr_curve(scores, labels)
      step_eval_pr!(recall, precision, x_latent, kendall)
      push!(kendall_auc, auprc(recall, precision))

      qvals = adjust(p_val, BenjaminiHochberg())
      power = sum(true_false .& (qvals .< q_val_crit)) / sum(true_false)

      n_sig = sum(qvals .< q_val_crit)
      fdr = (n_sig == 0) ? 0.0 : sum((.!true_false) .& (qvals .< q_val_crit)) / n_sig

      @printf(output_fdr, "%d\tkendall\t0\t%.6e\t%.6e\n", reso, power, fdr)

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

      scores = -p_val
      labels = true_false

      recall, precision = pr_curve(scores, labels)
      step_eval_pr!(recall, precision, x_latent, mwU)
      push!(mwU_auc, auprc(recall, precision))

      qvals = adjust(p_val, BenjaminiHochberg())
      power = sum(true_false .& (qvals .< q_val_crit)) / sum(true_false)

      n_sig = sum(qvals .< q_val_crit)
      fdr = (n_sig == 0) ? 0.0 : sum((.!true_false) .& (qvals .< q_val_crit)) / n_sig

      @printf(output_fdr, "%d\tmwU\t1\t%.6e\t%.6e\n", reso, power, fdr)

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

      scores = bf
      labels = true_false

      recall, precision = pr_curve(scores, labels)
      step_eval_pr!(recall, precision, x_latent, baypass)
      push!(baypass_auc, auprc(recall, precision))

      power = sum(true_false .& (bf .> baypass_crit)) / sum(true_false)

      n_sig = sum(bf .> baypass_crit)
      fdr = (n_sig == 0) ? 0.0 : sum((.!true_false) .& (bf .> baypass_crit)) / n_sig

      @printf(output_fdr, "%d\tbaypass\t4\t%.6e\t%.6e\n", reso, power, fdr)

      # matrix
      p_val = Float64[]

      open(folder * "data$(i)/matrix/sample$(reso)/p_val.txt", "r") do f
        while !eof(f)
          line = strip(readline(f))
          push!(p_val, parse(Float64, line))
        end
      end

      scores = -p_val
      labels = true_false

      recall, precision = pr_curve(scores, labels)
      step_eval_pr!(recall, precision, x_latent, matrix)
      push!(matrix_auc, auprc(recall, precision))

      qvals = adjust(p_val, BenjaminiHochberg())
      power = sum(true_false .& (qvals .< q_val_crit)) / sum(true_false)

      n_sig = sum(qvals .< q_val_crit)
      fdr = (n_sig == 0) ? 0.0 : sum((.!true_false) .& (qvals .< q_val_crit)) / n_sig

      @printf(output_fdr, "%d\tmatrix\t2\t%.6e\t%.6e\n", reso, power, fdr)

      # LFMM
      p_val = Float64[]

      open(folder * "data$(i)/LFMM/sample$(reso)/p-values.txt", "r") do f
        while !eof(f)
          line = strip(readline(f))
          push!(p_val, parse(Float64, line))
        end
      end

      scores = -p_val
      labels = true_false

      recall, precision = pr_curve(scores, labels)
      step_eval_pr!(recall, precision, x_latent, lfmm)
      push!(lfmm_auc, auprc(recall, precision))

      qvals = adjust(p_val, BenjaminiHochberg())
      power = sum(true_false .& (qvals .< q_val_crit)) / sum(true_false)

      n_sig = sum(qvals .< q_val_crit)
      fdr = (n_sig == 0) ? 0.0 : sum((.!true_false) .& (qvals .< q_val_crit)) / n_sig

      @printf(output_fdr, "%d\tlfmm\t3\t%.6e\t%.6e\n", reso, power, fdr)
    end
  end

  kendall ./= all_case
  mwU ./= all_case
  baypass ./= all_case
  matrix ./= all_case
  lfmm ./= all_case

  open(folder * "/pr_data$(reso).txt", "w") do f
    println(f, "recall\tkendall\tmwU\tbaypass\tmatrix\tlfmm")

    for i in 1:sep
      @printf(f, "%.5f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", x_latent[i], kendall[i],
        mwU[i], baypass[i], matrix[i], lfmm[i])
    end
  end

  open(folder * "/aucpr_data$(reso).txt", "w") do f
    println(f, "kendall\tmwU\tbaypass\tmatrix\tlfmm")

    for i in 1:length(kendall_auc)
      @printf(f, "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", kendall_auc[i],
        mwU_auc[i], baypass_auc[i], matrix_auc[i], lfmm_auc[i])
    end
  end
end

close(output_fdr)
