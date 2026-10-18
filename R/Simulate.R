# ======================================================================
# Simulate.R   the Gillespie simulator
# ======================================================================
#
#   Outcome page, section 1.6 (steps 10-14 of the build chart, 1.1).
#   Each cell has six events: a birth, a death, and a hop toward
#   each of its 4 sides (rates in Rates.R). R = sum of all rates.
#     10. total   R = sum of every rate on the board
#     11. wait    tau = -log(u1) / R, u1 uniform on (0, 1)
#     12. choose  cell i with chance tot_i / R, then event e of cell i
#                 with chance rate_e / tot_i
#     13. apply   birth: n_i + 1; death: n_i - 1;
#                 hop i -> j: n_i - 1, n_j + 1;
#                 then update only the rates of the cells that changed
#     14. many runs: Run_many(), one seed per run (Seed_of())
#   Rain changes in steps of dt_rain years (1 = yearly; 1/12 =
#   monthly, as in E3). When the next event would fall past a step,
#   the clock stops there, the birth rates take the new rain, and a
#   fresh wait is drawn (waits have no memory, so the run stays
#   exact). The run ends at t_max or when no hare is left.
#
#   Two engines run this algorithm:
#     "R"   Gillespie_R() below: plain R, easy to read, slow
#     "C"   src/hh_engine.c through R/Engine.R: the same steps, about
#           100 times faster; needs a C compiler once
#   Both draw from R's random numbers, so set.seed() makes either
#   repeatable. options(hh.engine = "R") forces plain R.
#
#   The simulator returns data only; it never prints or plots.
#
#   Run_sim()  Run_many()  Gillespie_R()  Seed_of()
# ======================================================================

# ----------------------------------------------------------------------
# Run_sim: one run of the model
#
#   p          parameters from Set_params()
#   map        map from Build_map(); Run_sim adds a from p$alpha
#   rho        rain, one value per step of dt_rain years, 0 or more
#              (see Rho.R); the last value holds after the last step
#   n0         starting hares: one number (the same in every cell) or
#              one per cell, length C (see Place_hares())
#   t_max      end time in years, finite; default: the rain's length
#   dt_rain    length of one rain step in years (1 = yearly)
#   dt_rec     recording step for the total N in years (NULL = none)
#   record     TRUE = also keep the board at every recording time
#   snap_t     times to save the whole board (for example when the
#              rain crosses 0.995, 0.990, ... in E3); NULL = none
#   occ_window c(from, to): years each cell holds hares inside this
#              window (E2: years 180 to 200); NULL = off
#   frac       record the first time N falls below frac * N(0)
#              (E3: 0.1); NA = off
#   groups     named list of cell numbers, e.g. list(poor = ...,
#              good = ...); totals kept per group at each record
#   seed       a number makes the run repeatable; NULL = no reset
#   engine     "R", "C" or NULL (= getOption("hh.engine", "auto"))
#
#   returns  a list:
#     t, N, rho_t     recording times; total hares; rain at each time
#     n               board at each record (n_rec x C) if record
#     G               group totals (n_rec x n_groups) if groups
#     snap            boards at snap_t (n_snap x C) if snap_t
#     occ_time        years each cell held hares in the window
#     t_frac          first time N < frac * N(0) (NA if never)
#     extinct         TRUE if every cell emptied
#     t_extinct       when that happened (NA if it never did)
#     events          counts of births, deaths and hops
#     n0, n_final     first and last board (length C)
#     rates_final     last rate table (C x 6; engine "R" only)
#     f_final         last food per cell (length C; engine "R" only)
#     p, map, rho, dt_rain, t_max, t_end, engine   inputs and end time
# ----------------------------------------------------------------------
Run_sim = function(p, map, rho, n0, t_max = NULL, dt_rain = 1,
                   dt_rec = 1 / 12, record = TRUE, snap_t = NULL,
                   occ_window = NULL, frac = NA, groups = NULL,
                   seed = NULL, engine = NULL) {
  # check the inputs (no random numbers used here)
  if (!is.numeric(rho) || length(rho) == 0 || any(!is.finite(rho)) ||
      any(rho < 0)) {
    stop("Run_sim: rain must be numbers, 0 or more")
  }
  if (length(dt_rain) != 1 || !is.finite(dt_rain) || dt_rain <= 0) {
    stop("Run_sim: dt_rain must be one number above 0")
  }
  if (is.null(t_max)) {
    if (length(rho) == 1) {
      stop("Run_sim: rain is one number; give t_max in years")
    }
    t_max = length(rho) * dt_rain
  }
  if (length(t_max) != 1 || !is.finite(t_max) || t_max <= 0) {
    stop("Run_sim: t_max must be one finite number above 0")
  }
  if (!is.null(dt_rec) && (length(dt_rec) != 1 ||
                           !is.finite(dt_rec) || dt_rec <= 0)) {
    stop("Run_sim: dt_rec must be one number above 0, or NULL")
  }
  if (record && is.null(dt_rec)) {
    stop("Run_sim: record = TRUE needs a recording step dt_rec")
  }
  if (record && !is.null(snap_t)) {
    stop("Run_sim: use record = TRUE or snap_t, not both")
  }
  if (!(length(n0) %in% c(1, map$C))) {
    stop("Run_sim: n0 needs length 1 or ", map$C,
         " (one per cell); see Place_hares()")
  }
  if (any(!is.finite(n0)) || any(n0 < 0) || any(n0 != round(n0))) {
    stop("Run_sim: n0 must be whole numbers, 0 or more")
  }
  if (!is.null(occ_window) && (length(occ_window) != 2 ||
                               occ_window[2] <= occ_window[1])) {
    stop("Run_sim: occ_window must be c(from, to) with from < to")
  }
  engine = Engine_choose(engine)
  map$a = Hop_table(map, p$alpha)        # C x 4: crossing chances

  # recording times, snapshots and groups
  rec_t = numeric(0)                     # recording times of N
  if (!is.null(dt_rec)) rec_t = seq(0, t_max, by = dt_rec)
  if (record) snap_t = rec_t
  snap_t = if (is.null(snap_t)) numeric(0) else sort(snap_t)
  gid = rep(0L, map$C)                   # length C: group of each cell
  for (g in seq_along(groups)) {
    if (any(gid[groups[[g]]] > 0)) {
      stop("Run_sim: a cell sits in two groups")
    }
    gid[groups[[g]]] = g
  }
  n0 = rep_len(as.integer(round(n0)), map$C)

  if (!is.null(seed)) set.seed(seed)
  core = if (engine == "C") {
    Gillespie_C(p, map, rho, dt_rain, n0, t_max, rec_t, snap_t,
                occ_window, frac, gid, length(groups))
  } else {
    Gillespie_R(p, map, rho, dt_rain, n0, t_max, rec_t, snap_t,
                occ_window, frac, gid, length(groups))
  }

  out = list(type = "run", p = p, map = map, rho = rho,
             dt_rain = dt_rain, t_max = t_max, t_end = core$t_end,
             extinct = sum(core$end) == 0, t_extinct = core$t_ext,
             events = core$events, n0 = n0, n_final = core$end,
             engine = engine)
  if (length(rec_t) > 0) {
    out$t = rec_t
    out$N = core$N
    out$rho_t = Rho_at(rec_t, rho, dt_rain)
  }
  if (record) out$n = core$snap
  if (!record && length(snap_t) > 0) {
    out$snap_t = snap_t
    out$snap = core$snap
  }
  if (length(groups) > 0) {
    out$G = core$G
    colnames(out$G) = names(groups)
  }
  if (!is.null(occ_window)) out$occ_time = core$occ
  if (!is.na(frac)) out$t_frac = core$t_frac
  if (engine == "R") {
    out$rates_final = core$rates_final
    out$f_final = core$f_final
  }
  return(out)
}

# ----------------------------------------------------------------------
# Gillespie_R: the algorithm in plain R (called by Run_sim)
#
#   p, map (with a), rho, dt_rain, n0 (length C), t_max   as in Run_sim
#   rec_t, snap_t   record and snapshot times (may be empty)
#   occ_window      c(from, to) or NULL
#   frac            share of N(0) for t_frac, or NA
#   gid, n_g        group of each cell (0 = none) and number of groups
#
#   returns  a list: N, G, snap, occ, t_frac, end, t_ext, t_end,
#            events, rates_final, f_final
# ----------------------------------------------------------------------
Gillespie_R = function(p, map, rho, dt_rain, n0, t_max, rec_t, snap_t,
                       occ_window, frac, gid, n_g) {
  C = map$C
  n = n0                                 # length C: hares per cell
  N_tot = sum(n)
  thresh = if (is.na(frac)) -Inf else frac * N_tot
  t_frac = NA_real_
  G_tot = integer(n_g)                   # length n_g: group totals
  for (g in seq_len(n_g)) G_tot[g] = sum(n[gid == g])

  # rain: k_rho = step in force; t_change = when it ends (never, once
  # the last value is in force)
  k_rho = 1
  t_change = if (length(rho) > 1) dt_rain else Inf
  f = rho[k_rho] * map$q                 # length C: food per cell

  rates = matrix(0, C, 6)                # C x 6: birth, death, 4 hops
  for (i in seq_len(C)) rates[i, ] = Cell_rates(i, n[i], f[i], map, p)
  tot = rowSums(rates)                   # length C: total per cell

  n_rec = length(rec_t)
  out_N = integer(n_rec)                 # n_rec: total hares
  out_G = matrix(0L, n_rec, n_g)         # n_rec x n_g: group totals
  k = 1                                  # next record to fill
  n_snap = length(snap_t)
  out_snap = matrix(0L, n_snap, C)       # n_snap x C: saved boards
  s = 1                                  # next snapshot to fill
  do_occ = !is.null(occ_window)
  occ = numeric(C)                       # years held inside the window
  last = numeric(C)                      # last change of each cell

  t = 0
  t_ext = NA_real_
  ev = c(birth = 0, death = 0, hop = 0)

  repeat {
    if (N_tot == 0) {                    # no hare left: the trap
      t_ext = t
      break
    }

    # 10-11. total rate, and the wait: exponential with mean 1 / R
    cs = cumsum(tot)                     # length C
    R = cs[C]
    tau = if (R > 0) -log(runif(1)) / R else Inf
    t_next = t + tau
    t_stop = min(t_change, t_max)

    # records and snapshots: the board holds until t_next (or t_stop)
    t_fill = min(t_next, t_stop)
    while (k <= n_rec && rec_t[k] < t_fill) {
      out_N[k] = N_tot
      if (n_g > 0) out_G[k, ] = G_tot
      k = k + 1
    }
    while (s <= n_snap && snap_t[s] < t_fill) {
      out_snap[s, ] = n
      s = s + 1
    }

    # 11. no event before the rain changes or the run ends: move the
    # clock there, take the new rain, draw a fresh wait
    if (t_next >= t_stop) {
      t = t_stop
      if (t_stop >= t_max) break
      k_rho = k_rho + 1
      t_change = if (k_rho < length(rho)) k_rho * dt_rain else Inf
      f = rho[k_rho] * map$q
      for (i in seq_len(C)) rates[i, 1] = Birth_rate(n[i], f[i], p)
      tot = rowSums(rates)
      next
    }
    t = t_next

    # 12. choose the cell, then one of its six events
    i = which(cs > runif(1) * R)[1]
    cr = cumsum(rates[i, ])              # length 6
    e = which(cr > runif(1) * cr[6])[1]
    j = if (e >= 3) map$nbr[i, e - 2] else 0L
    changed = if (j > 0) c(i, j) else i

    # close the occupancy spells of the cells about to change
    if (do_occ) {
      for (c_i in changed) {
        if (n[c_i] > 0) {
          occ[c_i] = occ[c_i] + max(0, min(t, occ_window[2]) -
                                      max(last[c_i], occ_window[1]))
        }
        last[c_i] = t
      }
    }

    # 13. apply the event
    if (e == 1) {                        # birth
      n[i] = n[i] + 1L
      N_tot = N_tot + 1
      ev[1] = ev[1] + 1
      if (n_g > 0 && gid[i] > 0) G_tot[gid[i]] = G_tot[gid[i]] + 1L
    } else if (e == 2) {                 # death
      n[i] = n[i] - 1L
      N_tot = N_tot - 1
      ev[2] = ev[2] + 1
      if (n_g > 0 && gid[i] > 0) G_tot[gid[i]] = G_tot[gid[i]] - 1L
    } else {                             # hop into neighbour j
      n[i] = n[i] - 1L
      n[j] = n[j] + 1L
      ev[3] = ev[3] + 1
      if (n_g > 0) {
        if (gid[i] > 0) G_tot[gid[i]] = G_tot[gid[i]] - 1L
        if (gid[j] > 0) G_tot[gid[j]] = G_tot[gid[j]] + 1L
      }
    }
    if (is.na(t_frac) && N_tot < thresh) t_frac = t

    # update only the cells that changed
    for (c_i in changed) {
      rates[c_i, ] = Cell_rates(c_i, n[c_i], f[c_i], map, p)
      tot[c_i] = sum(rates[c_i, ])
    }
  }

  # the last board holds until t_max (all zeros after extinction)
  if (k <= n_rec) {
    out_N[k:n_rec] = N_tot
    if (n_g > 0) out_G[k:n_rec, ] = rep(G_tot, each = n_rec - k + 1)
  }
  if (s <= n_snap) out_snap[s:n_snap, ] = rep(n, each = n_snap - s + 1)
  if (do_occ) {
    held = which(n > 0)
    occ[held] = occ[held] + pmax(0, min(t, occ_window[2]) -
                                   pmax(last[held], occ_window[1]))
  }
  out = list(N = out_N, G = out_G, snap = out_snap, occ = occ,
             t_frac = t_frac, end = n, t_ext = t_ext, t_end = t,
             events = ev, rates_final = rates, f_final = f)
  return(out)
}

# ----------------------------------------------------------------------
# Run_many: n_runs independent runs with the same rain (step 14)
#
#   p, map, rho, t_max, dt_rain, dt_rec, record, snap_t, occ_window,
#   frac, groups, engine   as in Run_sim()
#   n0       as in Run_sim(), or a C x n_runs matrix: column r is the
#            starting board of run r
#   n_runs   number of runs
#   seeds    one seed per run (each run then repeats on its own,
#            whatever else runs with it); see Seed_of()
#   seed     one seed for the whole set (used when seeds is NULL)
#
#   returns  a list:
#     t_ext    length n_runs: when each run emptied (NA = survived)
#     t, rho_t   recording times and rain (length n_rec)
#     N        n_rec x n_runs: total hares
#     end      C x n_runs: the board at the end of each run
#     G        n_rec x n_groups x n_runs, if groups
#     snap     n_snap x C x n_runs, if snap_t
#     occ_time C x n_runs, if occ_window
#     t_frac   length n_runs, if frac
#     n_mean, occ   with record = TRUE: board averaged over runs
#                   (n_rec x C) and occupied cells (n_rec x n_runs)
#     p, map, rho, t_max, n_runs, groups, engine   inputs
# ----------------------------------------------------------------------
Run_many = function(p, map, rho, n0, n_runs, t_max = NULL, dt_rain = 1,
                    dt_rec = 1 / 12, record = FALSE, snap_t = NULL,
                    occ_window = NULL, frac = NA, groups = NULL,
                    seeds = NULL, seed = NULL, engine = NULL) {
  if (!is.null(seeds) && length(seeds) != n_runs) {
    stop("Run_many: give one seed per run")
  }
  if (is.null(seeds) && !is.null(seed)) set.seed(seed)
  engine = Engine_choose(engine)
  per_run = is.matrix(n0)
  if (per_run && (nrow(n0) != map$C || ncol(n0) != n_runs)) {
    stop("Run_many: n0 as a matrix must be C x n_runs")
  }
  out = list(type = "many", p = p, rho = rho, n_runs = n_runs,
             groups = groups, engine = engine,
             t_ext = rep(NA_real_, n_runs))
  for (r in seq_len(n_runs)) {
    n0_r = if (per_run) n0[, r] else n0
    res = Run_sim(p, map, rho, n0_r, t_max = t_max, dt_rain = dt_rain,
                  dt_rec = dt_rec, record = record, snap_t = snap_t,
                  occ_window = occ_window, frac = frac, groups = groups,
                  seed = if (is.null(seeds)) NULL else seeds[r],
                  engine = engine)
    out$t_ext[r] = res$t_extinct
    if (r == 1) {                        # set up the tables
      out$map = res$map
      out$t_max = res$t_max
      out$end = matrix(0L, map$C, n_runs)          # C x n_runs
      if (!is.null(res$t)) {
        out$t = res$t
        out$rho_t = res$rho_t
        out$N = matrix(0L, length(res$t), n_runs)  # n_rec x n_runs
      }
      if (!is.null(res$G)) {
        out$G = array(0L, c(nrow(res$G), length(groups), n_runs),
                      dimnames = list(NULL, names(groups), NULL))
      }
      if (!is.null(res$snap)) {
        out$snap_t = res$snap_t
        out$snap = array(0L, c(nrow(res$snap), map$C, n_runs))
      }
      if (!is.null(res$occ_time)) {
        out$occ_time = matrix(0, map$C, n_runs)    # C x n_runs
      }
      if (!is.null(res$t_frac)) out$t_frac = rep(NA_real_, n_runs)
      if (record) {
        n_sum = matrix(0, length(res$t), map$C)    # n_rec x C
        out$occ = matrix(0L, length(res$t), n_runs)
      }
    }
    out$end[, r] = res$n_final
    if (!is.null(res$t)) out$N[, r] = res$N
    if (!is.null(res$G)) out$G[, , r] = res$G
    if (!is.null(res$snap)) out$snap[, , r] = res$snap
    if (!is.null(res$occ_time)) out$occ_time[, r] = res$occ_time
    if (!is.null(res$t_frac)) out$t_frac[r] = res$t_frac
    if (record) {
      n_sum = n_sum + res$n
      out$occ[, r] = rowSums(res$n > 0)
    }
  }
  if (record) out$n_mean = n_sum / n_runs
  return(out)
}

# ----------------------------------------------------------------------
# Seed_of: a whole-number seed from any key, the same on every computer
#
#   Seed_of("E2", 3, 0.785, 12) gives the seed of run 12 of board 3 at
#   rain 0.785. A string hash: h = (31 h + code) mod (2^31 - 1).
#
#   ...   parts of the key (numbers or text)
#
#   returns  one whole number from 0 to 2^31 - 2
# ----------------------------------------------------------------------
Seed_of = function(...) {
  key = paste(..., sep = "|")
  h = 0
  for (x in utf8ToInt(key)) h = (h * 31 + x) %% 2147483647
  return(as.integer(h))
}
