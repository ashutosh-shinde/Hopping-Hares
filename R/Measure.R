# ======================================================================
# Measure.R   numbers the outcome page reports
# ======================================================================
#
#   From the parameters alone
#     F_crit()         the one-cell fold f_c (section 2.1)
#     Cell_states()    the settled and tipping counts of one cell
#     Exact_T_one()    exact mean years until one cell empties (1.3)
#     Exact_T_two()    the same for two linked cells (checks, 1.6)
#     Exact_P_empty()  exact chance one cell is empty by time T (E2's
#                      window theory, figure "exact one-cell curves")
#   From runs
#     Fit_logistic()   tipping curve and rho_50 (E2)
#     Boot_rho50()     bootstrap interval of rho_50 (E2)
#     Km()             Kaplan-Meier survival curve (E2)
#     Wilson()         95% interval of a share (E4)
#     Alive_at()       share of runs alive at a given year
# ======================================================================

# ----------------------------------------------------------------------
# F_crit: the single-cell fold, the lowest food at which a cell can
# hold hares in the average (deterministic) picture
#
#   f_c = mu / (b0 * g*),   g* = x / (x + h) * (1 - (x + 1) / K),
#   x = sqrt(h^2 + h (K - 1)) - h      (so n* = x + 1)
#
#   g* is the top of the hill phi(n) s(n): births per hare peak at
#   n* = 5.5 hares. A cell of grass q keeps hares (alone) only while
#   rain > f_c / q. Example: q = 0.8 needs rain above 0.62.
#
#   p   parameters             uses mu, b0, h, K
#
#   returns  f_c, one number (0.498 at the defaults)
# ----------------------------------------------------------------------
F_crit = function(p) {
  x = sqrt(p$h^2 + p$h * (p$K - 1)) - p$h
  g = x / (x + p$h) * (1 - (x + 1) / p$K)  # max of phi(n) * s(n)
  f_c = p$mu / (p$b0 * g)
  return(f_c)
}

# ----------------------------------------------------------------------
# Cell_states: where one sealed cell can rest, in the average picture
#
#   Settled counts solve births = deaths per hare:
#     b0 f phi(x) s(x) = mu
#   The hill crosses the line twice: x_minus (tipping count,
#   unstable) and x_plus (settled count, stable). Below the fold
#   (f < f_c) neither exists. The empty cell (x = 0) is always
#   stable. Example: at food 1.0, x_plus = 14.2 and x_minus = 1.8
#   (outcome page, section 2.1).
#
#   f   food, rho * q
#   p   parameters
#
#   returns  c(x_minus, x_plus); c(0, 0) below the fold
# ----------------------------------------------------------------------
Cell_states = function(f, p) {
  g = function(x) {
    return(p$b0 * f * (x - 1) / (x - 1 + p$h) * (1 - x / p$K) - p$mu)
  }
  x_top = 1 + sqrt(p$h^2 + p$h * (p$K - 1)) - p$h   # top of the hill
  if (f <= 0 || g(x_top) < 0) return(c(0, 0))
  lo = uniroot(g, c(1 + 1e-9, x_top), tol = 1e-10)$root
  hi = uniroot(g, c(x_top, p$K), tol = 1e-10)$root
  return(c(lo, hi))
}

# ----------------------------------------------------------------------
# Exact_T_one: exact mean years until one isolated cell empties
#
#   The master equation (section 1.3) gives, one event ahead:
#   (B(n) + D(n)) T(n) - B(n) T(n + 1) - D(n) T(n - 1) = 1,
#   T(0) = 0, for n = 1..K (births stop at K)
#
#   f   food in the cell, rho * q, 0 or more
#   p   parameters
#
#   Example: a full cell (13 hares) lives 544 years at food 1.0 and
#   74 years at food 0.8.
#
#   returns  T, length K: element n is T(n), in years
# ----------------------------------------------------------------------
Exact_T_one = function(f, p) {
  K = p$K
  A = matrix(0, K, K)                    # K x K: the linear system
  for (n in 1:K) {
    B = Birth_rate(n, f, p)
    D = Death_rate(n, p)
    A[n, n] = B + D
    if (n < K) A[n, n + 1] = -B
    if (n > 1) A[n, n - 1] = -D          # T(0) = 0 drops out
  }
  T_n = solve(A, rep(1, K))
  return(T_n)
}

# ----------------------------------------------------------------------
# Exact_T_two: exact mean years until both cells of a 2-cell map are
# empty, for every starting pair (n1, n2)
#
#   R T(n1, n2) - sum over the 6 events [rate x T(pair after)] = 1,
#   T(0, 0) = 0. Pairs with a count above n_max are left out (a
#   reflecting edge). At the defaults the answer does not change in
#   the first 4 digits between n_max = 35 and 55, up to rain 1.1.
#
#   map     2-cell map from Build_map()
#   rho     rain, one number
#   p       parameters
#   n_max   largest count kept per cell
#
#   returns  T, (n_max + 1) x (n_max + 1): T[n1 + 1, n2 + 1] in years
# ----------------------------------------------------------------------
Exact_T_two = function(map, rho, p, n_max = 35) {
  if (map$C != 2) stop("Exact_T_two: the map needs exactly 2 cells")
  f = rho * map$q                        # length 2: food
  a = Hop_table(map, p$alpha)            # 2 x 4: crossing chances
  hop12 = (p$m / 4) * sum(a[1, map$nbr[1, ] == 2])
  hop21 = (p$m / 4) * sum(a[2, map$nbr[2, ] == 1])
  S = (n_max + 1)^2                      # number of pairs
  A = matrix(0, S, S)                    # S x S: the linear system
  b = rep(1, S)
  for (n1 in 0:n_max) {
    for (n2 in 0:n_max) {
      s = n1 * (n_max + 1) + n2 + 1      # row of pair (n1, n2)
      if (n1 == 0 && n2 == 0) {
        A[s, s] = 1
        b[s] = 0
        next
      }
      # 6 x 3: each event's rate, then the pair it leads to
      ev = rbind(c(Birth_rate(n1, f[1], p), n1 + 1, n2),
                 c(Death_rate(n1, p),       n1 - 1, n2),
                 c(Birth_rate(n2, f[2], p), n1,     n2 + 1),
                 c(Death_rate(n2, p),       n1,     n2 - 1),
                 c(hop12 * n1,              n1 - 1, n2 + 1),
                 c(hop21 * n2,              n1 + 1, n2 - 1))
      keep = ev[, 1] > 0 & ev[, 2] <= n_max & ev[, 3] <= n_max
      ev = ev[keep, , drop = FALSE]
      A[s, s] = sum(ev[, 1])
      for (e in seq_len(nrow(ev))) {
        s2 = ev[e, 2] * (n_max + 1) + ev[e, 3] + 1
        A[s, s2] = A[s, s2] - ev[e, 1]
      }
    }
  }
  T_pair = matrix(solve(A, b), n_max + 1, n_max + 1, byrow = TRUE)
  return(T_pair)
}

# ----------------------------------------------------------------------
# Exact_P_empty: exact chance that one sealed cell is empty by time T
#
#   The master equation of one cell, dP/dt = M P, solved as
#   P(T) = exp(M T) P(0). M is the generator: column n holds the
#   rates out of n hares (birth to n + 1, death to n - 1). The cell
#   starts at its settled count of a normal year, round(x_plus(1.0)).
#
#   How exp(M T) is computed (no eigenvectors: those turn singular
#   for K = 100 and 200, because the empty state traps):
#     1. a small step dt = 2^-10 years; E = exp(M dt) from its Taylor
#        series (|M dt| stays near 1, so 30 terms are plenty)
#     2. square E again and again: E^2 = exp(M 2 dt), E^4, E^8, ...
#     3. write T / dt in binary and apply the matching powers to the
#        start vector, one after another
#   Every power holds chances (nothing negative, columns sum to 1), so
#   squaring stays accurate even for T = 1e8 years (47 squarings).
#   Example: f = 0.9, T = 200, Omega = 1 takes 47 small matrix
#   squarings of a 21 x 21 matrix.
#
#   System size Omega scales the cell: K = 20 Omega and h = 2 Omega.
#   The average picture (and the fold) barely moves while the hares
#   multiply, so the curves show what more hares do to the switch
#   (outcome page, E2 figure of exact one-cell curves: about 14, 28,
#   71 and 142 hares for Omega = 1, 2, 5, 10).
#
#   f       foods (vector)
#   T       windows in years (vector; whole multiples of 2^-10)
#   Omega   system size, a whole number
#   p       parameters (uses b0, mu)
#
#   returns  length(T) x length(f) matrix of P(empty by T)
# ----------------------------------------------------------------------
Exact_P_empty = function(f, T, Omega = 1, p = Set_params()) {
  q = p
  q$K = 20 * Omega
  q$h = 2 * Omega
  start = round(Cell_states(1.0, q)[2])
  dt = 2^-10                             # step 1: the small step
  steps = round(T / dt)                  # T in small steps
  if (any(abs(steps * dt - T) > 1e-9 * pmax(1, T))) {
    stop("Exact_P_empty: T must be whole multiples of 2^-10 years")
  }
  n_pow = max(1, ceiling(log2(max(steps) + 1)))
  out = matrix(NA_real_, length(T), length(f))
  for (j in seq_along(f)) {
    n = 0:q$K
    B = sapply(n, function(x) Birth_rate(x, f[j], q))
    D = q$mu * n
    M = matrix(0, q$K + 1, q$K + 1)      # generator: columns = from
    for (i in seq_along(n)) {
      if (i <= q$K) M[i + 1, i] = M[i + 1, i] + B[i]
      if (i > 1) M[i - 1, i] = M[i - 1, i] + D[i]
      M[i, i] = M[i, i] - (B[i] + D[i])
    }
    # step 1: E = exp(M dt) by Taylor series
    A = M * dt
    E = diag(q$K + 1)
    term = E
    for (k in 1:30) {
      term = (term %*% A) / k
      E = E + term
    }
    # step 2: the powers E^(2^k), k = 0, 1, ..., n_pow - 1
    pow = vector("list", n_pow)
    pow[[1]] = E
    for (k in seq_len(n_pow - 1) + 1) pow[[k]] = pow[[k - 1]] %*%
      pow[[k - 1]]
    # step 3: apply the powers named by the binary digits of steps
    for (k in seq_along(T)) {
      P = numeric(q$K + 1)
      P[start + 1] = 1                   # start: the settled count
      left = steps[k]
      for (b in seq_len(n_pow)) {
        if (left %% 2 == 1) P = as.vector(pow[[b]] %*% P)
        left = left %/% 2
      }
      out[k, j] = min(max(P[1], 0), 1)   # P[1]: zero hares
    }
  }
  return(out)
}

# ----------------------------------------------------------------------
# Fit_logistic: the tipping curve of a set of runs
#
#   P(collapse | rain x) = 1 / (1 + exp((x - rho50) / w))
#
#   rho50 is the rain at which half the runs collapse; w > 0 sets the
#   width of the switch (from 90% to 10% collapsed takes 2 w ln 9 of
#   rain). Maximum likelihood on one 0/1 outcome per run, from three
#   starting widths (0.005, 0.02, 0.08); the best fit wins.
#
#   x   rain of each run
#   y   1 = collapsed, 0 = survived (one per run)
#
#   returns  c(rho50 = ..., width = ...); width NA when every run
#            collapsed (rho50 = top rain) or none did (lowest rain)
# ----------------------------------------------------------------------
Fit_logistic = function(x, y) {
  y = as.numeric(y)
  if (all(y == 1)) return(c(rho50 = max(x), width = NA))
  if (all(y == 0)) return(c(rho50 = min(x), width = NA))
  xs = sort(unique(x))
  k = sapply(xs, function(v) sum(y[x == v]))       # collapsed
  n = sapply(xs, function(v) sum(x == v))          # runs
  nll = function(th) {
    z = (xs - th[1]) / exp(th[2])
    l1p = ifelse(z > 30, z, log1p(exp(pmin(z, 30))))   # log(1 + e^z)
    return(sum(k * l1p + (n - k) * (l1p - z)))
  }
  i = which.min(abs(k / n - 0.5))        # share nearest one half
  best = NULL
  for (lw0 in log(c(0.005, 0.02, 0.08))) {
    fit = optim(c(xs[i], lw0), nll, method = "Nelder-Mead",
                control = list(reltol = 1e-12, maxit = 4000))
    if (is.null(best) || fit$value < best$value) best = fit
  }
  return(c(rho50 = best$par[1], width = exp(best$par[2])))
}

# ----------------------------------------------------------------------
# Boot_rho50: bootstrap interval of rho50
#
#   Every rain reuses the same 20 start boards, so we resample the
#   start boards (with replacement), refit, and take the middle 95% of
#   the refitted rho50 (outcome page, E2 section 5).
#
#   x, y    as in Fit_logistic()
#   board   start board of each run (1-20)
#   n_boot  number of resamples
#   seed    random seed (NULL = no reset)
#
#   returns  c(lo = 2.5% point, hi = 97.5% point)
# ----------------------------------------------------------------------
Boot_rho50 = function(x, y, board, n_boot = 200, seed = 1) {
  if (!is.null(seed)) set.seed(seed)
  ids = sort(unique(board))
  v = numeric(n_boot)
  for (b in seq_len(n_boot)) {
    pick = sample(ids, length(ids), replace = TRUE)
    rows = unlist(lapply(pick, function(r) which(board == r)))
    v[b] = Fit_logistic(x[rows], y[rows])[["rho50"]]
  }
  out = quantile(v, c(0.025, 0.975), names = FALSE)
  return(c(lo = out[1], hi = out[2]))
}

# ----------------------------------------------------------------------
# Km: Kaplan-Meier survival curve
#
#   S(t) = product over collapse times t_i <= t of (1 - 1 / at risk)
#
#   times    collapse years of each run (NA = alive at the end)
#   censor   last year watched; runs alive then are censored
#
#   returns  data frame: t (0 and each collapse time), S just after
# ----------------------------------------------------------------------
Km = function(times, censor) {
  t = ifelse(is.na(times), Inf, times)
  ev = sort(t[is.finite(t) & t <= censor])
  at_risk = length(t)
  S = 1
  out_t = 0
  out_s = 1
  for (v in ev) {                        # no ties in continuous time
    S = S * (at_risk - 1) / at_risk
    at_risk = at_risk - 1
    out_t = c(out_t, v)
    out_s = c(out_s, S)
  }
  return(data.frame(t = out_t, S = out_s))
}

# ----------------------------------------------------------------------
# Wilson: 95% Wilson interval of a share k / n
#
#   k   successes
#   n   tries
#   z   1.96 for 95%
#
#   returns  c(lo, hi)
# ----------------------------------------------------------------------
Wilson = function(k, n, z = 1.96) {
  p = k / n
  den = 1 + z^2 / n
  mid = (p + z^2 / (2 * n)) / den
  half = z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  return(c(lo = mid - half, hi = mid + half))
}

# ----------------------------------------------------------------------
# Alive_at: share of runs that still had hares at a given year
#
#   t_ext   collapse years of the runs (NA = alive at the end)
#   year    time in years
#
#   returns  one number from 0 to 1
# ----------------------------------------------------------------------
Alive_at = function(t_ext, year) {
  alive = is.na(t_ext) | t_ext > year
  share = mean(alive)
  return(share)
}
