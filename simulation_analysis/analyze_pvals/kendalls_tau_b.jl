using SpecialFunctions

function logaddexp(a::Float64, b::Float64)::Float64
  if isinf(a) && a < 0
    return b
  elseif isinf(b) && b < 0
    return a
  elseif isinf(a) || isinf(b)
    return Inf
  end

  m = max(a, b)
  return m + log1p(exp(-abs(a - b)))
end

function logsubexp(a::Float64, b::Float64)::Float64
  if b == -Inf
    return a
  end

  if a <= b
    if a < b
      error("a >= b is not satisfied")
    end

    return -Inf
  end

  return a + log1p(-exp(b - a))
end

function kendall_tau_b(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}) 
  if length(y) != length(x)
    error("length is different")
  end
  idx = .!(isnan.(x) .| isnan.(y))

  x = x[idx]
  y = y[idx]

  n = length(x)
  if n < 2
    return (0.0, 0.0, 1.0, 0.0)
  end

  c1 = Dict{Float64, Int}()
  c2 = Dict{Float64, Int}()

  @inbounds for i in 1:n
    xi = x[i]
    yi = y[i]

    c1[xi] = get(c1, xi, 0) + 1
    c2[yi] = get(c2, yi, 0) + 1
  end

  if length(c1) == 1 || length(c2) == 1
    return (0.0, 0.0, 1.0, 0.0)
  end

  n0 = div(n * (n - 1), 2)

  n1 = 0
  for (i, t) in c1
    if t >= 2
      n1 += div(t * (t - 1), 2)
    end
  end

  n2 = 0
  for (i, u) in c2
    if u >= 2
      n2 += div(u * (u - 1), 2)
    end
  end

  nc = Int64(0)
  nd = Int64(0)

  @inbounds for i in 1:(n-1)
    xi = x[i]
    yi = y[i]

    for j in (i+1):n
      xj = x[j]
      yj = y[j]

      if (xi > xj && yi > yj) || (xi < xj && yi < yj)
        nc += 1
      elseif (xi > xj && yi < yj) || (xi < xj && yi > yj)
        nd += 1
      end
    end
  end

  # tau-b
  denom = sqrt(Float64(n0 - n1) * Float64(n0 - n2))

  tau = 0.0
  if n0 != n1 && n0 != n2
    tau = Float64(nc - nd) / denom
  end

  d = abs(nc - nd)
  if d == 0
      return (tau, 0.0, 1.0, 0.0)
  end

  # z-value
  log_v0 = log(Float64(n)) + log(Float64(n - 1)) + log(2.0 * n + 5.0)

  log_vt = -Inf
  log_vu = -Inf
  log_v11 = -Inf
  log_v12 = -Inf
  log_v21 = -Inf
  log_v22 = -Inf

  for (i, t) in c1
    if t >= 2
      tf = Float64(t)
      log_vt = logaddexp(log_vt, log(tf) + log(tf - 1.0) + log(2.0 * tf + 5.0))
      log_v11 = logaddexp(log_v11, log(tf) + log(tf - 1.0))

      if t >= 3
        log_v12 = logaddexp(log_v12, log(tf) + log(tf - 1.0) + log(tf - 2.0))
      end
    end
  end

  for (i, u) in c2
    if u >= 2
      uf = Float64(u)
      log_vu = logaddexp(log_vu, log(uf) + log(uf - 1.0) + log(2.0 * uf + 5.0))
      log_v21 = logaddexp(log_v21, log(uf) + log(uf - 1.0))

      if u >= 3
        log_v22 = logaddexp(log_v22, log(uf) + log(uf - 1.0) + log(uf - 2.0))
      end
    end
  end

  log_v1 = log_v11 + log_v21 - log(Float64(n)) - log(Float64(n - 1)) - log(2.0)
  log_v2 = log_v12 + log_v22 - log(Float64(n)) - log(Float64(n - 1)) - log(Float64(n - 2)) - log(9.0)

  log_va = logaddexp(log_v0 - log(18.0), logaddexp(log_v1, log_v2))
  log_vb = logaddexp(log_vt, log_vu) - log(18.0)
  log_v = logsubexp(log_va, log_vb)

  log_z = log(Float64(d)) - 0.5 * log_v
  z = exp(log_z)

  zn = z / sqrt(2.0)
  logp = log(erfcx(zn)) - zn^2
  p = exp(logp)

  return (tau, z, p, logp)
end
  