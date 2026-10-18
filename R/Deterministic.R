# ======================================================================
# Deterministic.R   the deterministic board: averages, no luck
# ======================================================================
#
#   Outcome page, section 1.4 and steps 5-9 of the build chart (1.1);
#   results in section 2.
#
#   Symbols
#     x_i     average hares in cell i (a real number, not whole)
#     q_i     grass of cell i; rho the rain; f_i = rho q_i the food
#     B       births per year, b0 f x phi(x) s(x) for 1 < x < K, else 0
#     mu      deaths per hare per year
#     a_ij    crossing chance from i into neighbour j (Map.R)
#     m / 4   hop tries per hare per year toward each side
#
#   The average equations (step 5): one per cell,
#     dx_i/dt = B(x_i, rho q_i) - mu x_i            births, deaths
#               - (m/4) x_i sum_j a_ij              hops out
#               + (m/4) sum_j a_ji x_j              hops in
#   written below as dx/dt = F(x). They drop the curvature term
#   (1/2) B'' Var(n) of the exact average, so they hold 16-26% more
#   hares than the simulator at rain 1.0 (E1).
#
#   What we compute
#     step 6a settle: from 15 hares per cell, run the equations
#             forward (RK4) until nothing changes, then polish with
#             Newton, x <- x - J^-1 F(x), J = dF/dx the Jacobian
#     step 6b stability: the eigenvalues of J. The largest, lambda_1,
#             is the slowest recovery rate; a settled board is stable
#             while lambda_1 < 0, and 1 / |lambda_1| is its recovery
#             time (section 2.4, Plot_recovery)
#     step 7  continuation: lower the rain in steps of 0.0025 and
#             follow the settled board, each Newton solve starting
#             from the board one step wetter
#     step 8  folds: where lambda_1 reaches 0 the settled board
#             vanishes (a saddle-node). Bisection pins the rain to
#             1e-5; the eigenvector of lambda_1 marks the cells that
#             stop recovering first (section 2.3, Plot_staircase)
#     step 9  below the fold: settle again; stop if the board is
#             empty, else repeat steps 7-9 on the lower branch
#
#   A shortcut for the eigenvalues. Hops obey detailed balance with
#   weights pi_i = q_i^alpha: a_ij pi_i = a_ji pi_j. So
#   S = diag(pi)^(-1/2) J diag(pi)^(1/2) is symmetric, has the same
#   eigenvalues as J (all real), and R's symmetric eigen() is fast and
#   exact enough. An eigenvector w of S gives v = pi^(1/2) w of J.
#
#   Det_board()  Det_birth()  Det_dbirth()  Det_change()
#   Det_jacobian()  Det_settle()  Det_newton()  Det_slowest()
#   Det_stable_near()  Det_follow()  Det_path()  Det_ramp()
#   Det_start()  Fast_hop_limit()
# ======================================================================

# ----------------------------------------------------------------------
# Det_board: everything the average equations need for one map
#
#   map   map from Build_map()
#   p     parameters
#
#   returns  a list: map, p, C, w_out (hops out per hare per year,
#            length C), w_in (C x 4: hops in per hare of each
#            neighbour), nb (C x 4 neighbour index, C + 1 = none),
#            s (sqrt(pi) for the symmetric form), sym (TRUE if every
#            q > 0, so the symmetric form exists)
# ----------------------------------------------------------------------
Det_board = function(map, p) {
  a = Hop_table(map, p$alpha)            # C x 4: crossing chances
  C = map$C
  back = c(2, 1, 4, 3)                   # opposite of up, down, l, r
  w_in = matrix(0, C, 4)                 # C x 4
  for (d in 1:4) {
    has = map$nbr[, d] > 0
    j = map$nbr[has, d]
    w_in[has, d] = (p$m / 4) * a[cbind(j, back[d])]  # j -> i
  }
  nb = map$nbr
  nb[nb == 0] = C + 1L                   # points at an extra zero
  pi_w = map$q^p$alpha
  db = list(map = map, p = p, C = C,
            w_out = (p$m / 4) * rowSums(a), w_in = w_in, nb = nb,
            s = sqrt(pi_w), sym = all(pi_w > 0))
  return(db)
}

# ----------------------------------------------------------------------
# Det_birth: births per year at average counts x (vector)
#
#   B(x) = b0 f x (x - 1) / (x - 1 + h) (1 - x / K)  for 1 < x < K,
#   0 otherwise (the same rule as Birth_rate() for whole numbers)
#
#   x, f   counts and foods, same length
#   p      parameters
#
#   returns  births per year, one per cell
# ----------------------------------------------------------------------
Det_birth = function(x, f, p) {
  ok = x > 1 & x < p$K
  xx = ifelse(ok, x, 2)
  B = p$b0 * f * xx * (xx - 1) / (xx - 1 + p$h) * (1 - xx / p$K)
  return(ifelse(ok, B, 0))
}

# ----------------------------------------------------------------------
# Det_dbirth: dB/dx, the slope of the birth curve (for the Jacobian)
#
#   with u = x (x - 1), v = x - 1 + h, w = 1 - x / K:
#   dB/dx = b0 f (u' w v - u v / K - u w) / v^2,  u' = 2 x - 1
#
#   returns  one slope per cell (0 outside 1 < x < K)
# ----------------------------------------------------------------------
Det_dbirth = function(x, f, p) {
  ok = x > 1 & x < p$K
  xx = ifelse(ok, x, 2)
  u = xx * (xx - 1)
  v = xx - 1 + p$h
  w = 1 - xx / p$K
  dg = ((2 * xx - 1) * w * v - u * v / p$K - u * w) / v^2
  return(ifelse(ok, p$b0 * f * dg, 0))
}

# ----------------------------------------------------------------------
# Det_change: F(x), the change per year of every average count
#
#   F_i = B(x_i, rho q_i) - mu x_i - w_out_i x_i + sum_d w_in[i, d] x_j
#
#   x     average counts, length C
#   rho   rain
#   db    from Det_board()
#
#   returns  dx/dt, length C
# ----------------------------------------------------------------------
Det_change = function(x, rho, db) {
  p = db$p
  hop_in = rowSums(db$w_in * matrix(c(x, 0)[db$nb], db$C, 4))
  dx = Det_birth(x, rho * db$map$q, p) - p$mu * x - db$w_out * x +
    hop_in
  return(dx)
}

# ----------------------------------------------------------------------
# Det_jacobian: J = dF/dx, C x C
#
#   J = diag(B'(x) - mu - w_out) + (hops in: J[i, j] = w_in from j)
#
#   returns  C x C matrix
# ----------------------------------------------------------------------
Det_jacobian = function(x, rho, db) {
  p = db$p
  C = db$C
  J = matrix(0, C, C)
  for (d in 1:4) {
    has = which(db$nb[, d] <= C)
    J[cbind(has, db$nb[has, d])] = db$w_in[has, d]
  }
  diag(J) = Det_dbirth(x, rho * db$map$q, p) - p$mu - db$w_out
  return(J)
}

# ----------------------------------------------------------------------
# Det_settle: run the average equations forward until nothing changes
# (step 6a; also step 9, below each fold)
#
#   Fourth-order Runge-Kutta with a step short enough for the fastest
#   hops (0.05 years at m = 10, shorter for fast hoppers). Stops early
#   once every |dx_i/dt| < tol.
#
#   db    from Det_board()
#   x0    starting counts, length C
#   rho   rain, one number
#   T     longest time to run, years
#   tol   stop when max |F| falls below this
#
#   returns  counts at the end (negative rounding set to 0)
# ----------------------------------------------------------------------
Det_settle = function(db, x0, rho, T = 3000, tol = 1e-8) {
  fast = 2 * max(db$w_out) + db$p$mu + db$p$b0 * rho
  dt = min(0.05, 1.5 / fast)             # RK4 is stable below 2.8
  x = x0
  t = 0
  while (t < T) {
    k1 = Det_change(x, rho, db)
    if (max(abs(k1)) < tol) break
    k2 = Det_change(x + dt / 2 * k1, rho, db)
    k3 = Det_change(x + dt / 2 * k2, rho, db)
    k4 = Det_change(x + dt * k3, rho, db)
    x = x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    t = t + dt
  }
  return(pmax(x, 0))
}

# ----------------------------------------------------------------------
# Det_newton: a settled board near x, by Newton's method (step 6a)
#
#   x <- x + lam dx,  J dx = -F(x); lam halves while a count would go
#   below 0 (counts stay non-negative)
#
#   db, x, rho   as above
#   it           most iterations
#   tol          stop when max |F| < tol
#
#   returns  list(x = counts, ok = TRUE if max |F| < 1e-8 at the end)
# ----------------------------------------------------------------------
Det_newton = function(db, x, rho, it = 60, tol = 1e-10) {
  for (k in seq_len(it)) {
    Fv = Det_change(x, rho, db)
    if (max(abs(Fv)) < tol) return(list(x = x, ok = TRUE))
    dx = tryCatch(solve(Det_jacobian(x, rho, db), -Fv),
                  error = function(e) NULL)
    if (is.null(dx)) return(list(x = x, ok = FALSE))
    lam = 1
    while (any(x + lam * dx < -1e-9) && lam > 1e-4) lam = lam / 2
    x = x + lam * dx
  }
  ok = max(abs(Det_change(x, rho, db))) < 1e-8
  return(list(x = x, ok = ok))
}

# ----------------------------------------------------------------------
# Det_slowest: the slowest recovery rate lambda_1 (the largest
# eigenvalue of J) and, if asked, its pattern (step 6b)
#
#   db, x, rho   as above
#   vector       TRUE = also return the eigenvector
#
#   returns  list(lambda = lambda_1, pattern = |v| / max |v| (length
#            C, with vector = TRUE)). lambda < 0: a nudge fades
#            within about 1 / |lambda| years.
# ----------------------------------------------------------------------
Det_slowest = function(db, x, rho, vector = FALSE) {
  J = Det_jacobian(x, rho, db)
  if (db$sym) {
    S = J * outer(1 / db$s, db$s)        # diag(1/s) J diag(s)
    S = (S + t(S)) / 2                   # symmetric up to rounding
    e = eigen(S, symmetric = TRUE, only.values = !vector)
    lam = e$values[1]
    if (vector) v = abs(db$s * e$vectors[, 1])
  } else {
    e = eigen(J, only.values = !vector)
    k = which.max(Re(e$values))
    lam = Re(e$values[k])
    if (vector) v = abs(Re(e$vectors[, k]))
  }
  out = list(lambda = lam)
  if (vector) out$pattern = v / max(v)
  return(out)
}

# ----------------------------------------------------------------------
# Det_stable_near: the stable settled board at rain rho close to
# x_prev, or NULL if it is gone (one continuation step, step 7)
#
#   Newton from x_prev must converge, no cell may move by more than
#   2 hares (else it jumped to another state), and lambda_1 < 0.
#
#   returns  list(x, lambda) or NULL
# ----------------------------------------------------------------------
Det_stable_near = function(db, x_prev, rho) {
  nw = Det_newton(db, x_prev, rho)
  if (!nw$ok || max(abs(nw$x - x_prev)) > 2) return(NULL)
  lam = Det_slowest(db, nw$x, rho)$lambda
  if (lam >= 0) return(NULL)
  return(list(x = nw$x, lambda = lam))
}

# ----------------------------------------------------------------------
# Det_follow: follow the settled board down in rain; record every fold
# (steps 7-9; outcome page, section 2.2)
#
#   db       from Det_board()
#   rho_hi   start: settle from 15 hares per cell at this rain
#   rho_lo   stop below this rain (or once the board is empty)
#   step     continuation step in rain
#   tol      bisection tolerance for the fold rain
#   x_start  optional starting counts instead of 15 per cell
#
#   returns  a list:
#     rows    data frame: rho, N (hares), occupied (cells with
#             x >= 1), lambda (slowest rate), one row per step
#     folds   data frame: fold, rho, N_before, N_after, lambda
#     maps    list per fold: x_before (board just above the fold),
#             pattern (eigenvector), x_after (board it falls to)
# ----------------------------------------------------------------------
Det_follow = function(db, rho_hi = 1.2, rho_lo = 0.3, step = 0.0025,
                      tol = 1e-5, x_start = NULL) {
  x0 = if (is.null(x_start)) rep(15, db$C) else x_start
  x = Det_newton(db, Det_settle(db, x0, rho_hi), rho_hi)$x
  rho = rho_hi
  lam = Det_slowest(db, x, rho)$lambda
  rows = list()
  folds = list()
  maps = list()
  while (rho > rho_lo && sum(x) > 1e-6) {
    rows[[length(rows) + 1]] = c(rho, sum(x), sum(x >= 1), lam)
    nxt = Det_stable_near(db, x, rho - step)
    if (!is.null(nxt)) {                 # the board survives the step
      x = nxt$x
      lam = nxt$lambda
      rho = rho - step
      next
    }
    # the settled board vanishes in (rho - step, rho): bisect
    lo = rho - step
    hi = rho
    xh = x
    while (hi - lo > tol) {
      mid = (lo + hi) / 2
      xm = Det_stable_near(db, xh, mid)
      if (is.null(xm)) {
        lo = mid
      } else {
        hi = mid
        xh = xm$x
      }
    }
    sl = Det_slowest(db, xh, hi, vector = TRUE)
    after = Det_settle(db, xh, lo - 1e-4, T = 5000)   # where it falls
    k = length(folds) + 1
    folds[[k]] = c(k, hi, sum(xh), sum(after), sl$lambda)
    maps[[k]] = list(x_before = xh, pattern = sl$pattern,
                     x_after = after)
    rows[[length(rows) + 1]] = c(hi, sum(xh), sum(xh >= 1), sl$lambda)
    x = after
    rho = lo - 1e-4
    lam = Det_slowest(db, x, rho)$lambda
  }
  rows = as.data.frame(do.call(rbind, rows))
  names(rows) = c("rho", "N", "occupied", "lambda")
  folds = if (length(folds) == 0) {
    data.frame(fold = integer(0), rho = numeric(0),
               N_before = numeric(0), N_after = numeric(0),
               lambda = numeric(0))
  } else {
    as.data.frame(do.call(rbind, folds))
  }
  names(folds) = c("fold", "rho", "N_before", "N_after", "lambda")
  return(list(rows = rows, folds = folds, maps = maps))
}

# ----------------------------------------------------------------------
# Det_path: the average equations along a rain that changes in time
#
#   Runge-Kutta steps of dt; rain_fun(t) gives the rain at time t.
#   Optionally stops when N first falls below stop_below.
#
#   db          from Det_board()
#   x0          starting counts
#   rain_fun    function of t returning the rain
#   t_out       times to report N (increasing, from 0)
#   dt          step in years (cut for fast hoppers)
#   stop_below  stop when N < this (NA = never)
#
#   returns  list(N = N at t_out (NA after a stop), t_stop = the
#            crossing time, interpolated (NA if none), x = counts at
#            the end)
# ----------------------------------------------------------------------
Det_path = function(db, x0, rain_fun, t_out, dt = 0.05,
                    stop_below = NA) {
  fast = 2 * max(db$w_out) + db$p$mu + db$p$b0 * 1.5
  dt = min(dt, 1.5 / fast)
  x = x0
  t = 0
  N_out = rep(NA_real_, length(t_out))
  k = 1
  t_stop = NA_real_
  N_now = sum(x)
  while (k <= length(t_out) && t_out[k] <= 1e-12) {
    N_out[k] = N_now
    k = k + 1
  }
  t_end = max(t_out)
  while (t < t_end - 1e-12) {
    h = min(dt, t_end - t)
    k1 = Det_change(x, rain_fun(t), db)
    k2 = Det_change(x + h / 2 * k1, rain_fun(t + h / 2), db)
    k3 = Det_change(x + h / 2 * k2, rain_fun(t + h / 2), db)
    k4 = Det_change(x + h * k3, rain_fun(t + h), db)
    x_new = x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    N_new = sum(x_new)
    # report N at output times inside (t, t + h], by straight lines
    while (k <= length(t_out) && t_out[k] <= t + h + 1e-12) {
      w = (t_out[k] - t) / h
      N_out[k] = (1 - w) * N_now + w * N_new
      k = k + 1
    }
    if (!is.na(stop_below) && N_new < stop_below) {
      t_stop = t + h * (N_now - stop_below) / (N_now - N_new)
      x = x_new
      break
    }
    x = x_new
    N_now = N_new
    t = t + h
  }
  return(list(N = N_out, t_stop = t_stop, x = pmax(x, 0)))
}

# ----------------------------------------------------------------------
# Det_ramp: the deterministic board on an E3 ramp (lag alone, no luck)
#
#   Start settled at rain `from`, let the rain fall along a straight
#   line to `to` over `years`, then hold it; report the rain when N
#   first drops below frac of its start (E3's dashed black line).
#   Example (Shuffled): 0.608 on the 50-year ramp, 0.754 on the
#   4,050-year ramp.
#
#   db      from Det_board()
#   x0      settled board at rain `from` (Det_start())
#   years   ramp length
#   from, to, frac, hold   ramp ends, the 10% rule, years held after
#
#   returns  list(rho10 = rain at the crossing (NA if none), t10)
# ----------------------------------------------------------------------
Det_ramp = function(db, x0, years, from = 1.0, to = 0.4, frac = 0.1,
                    hold = 50) {
  r = (from - to) / years
  rain_fun = function(t) max(to, from - r * t)
  path = Det_path(db, x0, rain_fun, t_out = years + hold,
                  stop_below = frac * sum(x0))
  rho10 = if (is.na(path$t_stop)) NA_real_ else rain_fun(path$t_stop)
  return(list(rho10 = rho10, t10 = path$t_stop))
}

# ----------------------------------------------------------------------
# Det_start: the settled board at one rain, from 15 hares per cell
#
#   returns  counts, length C
# ----------------------------------------------------------------------
Det_start = function(db, rho = 1.0) {
  x = Det_settle(db, rep(15, db$C), rho)
  return(Det_newton(db, x, rho)$x)
}

# ----------------------------------------------------------------------
# Fast_hop_limit: the fold of a board whose hares hop very fast
#
#   Very fast hops spread hares as pi_i = q_i^alpha / sum(q^alpha)
#   (detailed balance), so the board acts as one population N with
#   x_i = N pi_i. Its fold is the rain where
#     max over N of [ sum_i B(N pi_i, rho q_i) ] / N = mu
#   (outcome page, section 2.5: 0.742 on the Shuffled board).
#
#   q   grass of every cell
#   p   parameters
#
#   returns  the fold rain
# ----------------------------------------------------------------------
Fast_hop_limit = function(q, p) {
  w = q^p$alpha / sum(q^p$alpha)
  Ns = seq(1, p$K * length(q), length.out = 2001)
  best = function(rho) {
    v = sapply(Ns, function(N) sum(Det_birth(N * w, rho * q, p)) / N)
    return(max(v) - p$mu)
  }
  return(uniroot(best, c(0.4, 1.2), tol = 1e-7)$root)
}
