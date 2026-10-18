# ======================================================================
# Run_E2.R   E2. Fixed rain: where does each board tip?
# ======================================================================
#
#   Outcome page, section 4. Every run starts from one of E1's start
#   boards B_1 ... B_20 and switches to the test rain at year 0;
#   collapse means no hare left. Per board:
#     1. coarse pass  rains 0.4, 0.5, ..., 1.2 (10 runs each); widen
#                     the range until the ends are reached
#     2. fine pass    every 0.025 over the rains where 1-9 of 10 runs
#                     collapsed (or the 10 -> 0 jump), 0.1 to each
#                     side (20 runs each)
#     3. year 300     runs alive at year 200 continue 100 years at
#                     every rain up to rho_50(200) + 0.15, and higher
#                     while the top rain still loses over 10% of runs
#     4. refine       every 0.005 inside the switch (20 runs each;
#                     survivors continue to year 300)
#   Then: rho_50(T) for windows of 50-300 years with bootstrap
#   intervals (Fit_logistic(), Boot_rho50()); the occupancy of each
#   cell in years 180-200; the pullback runs at the wet and dry edges
#   rho_50 +/- w ln 9 (10% and 90% of runs die); 2,000-year runs on
#   Real and Shuffled; and the exact one-cell curves
#   (Exact_P_empty()).
#
#   Needs E1's start boards (results/E1/start_boards.csv.gz, or the
#   quick ones). Run from the Hopping_hare_final folder:
#     source("run/Run_E2.R")
#   quick = TRUE  fewer runs, no refine pass, short long runs: about
#                 6 minutes with the C engine; writes results_quick/
#   quick = FALSE the full design (2,890 runs, 478,000 run-years, plus
#                 400 runs of 2,000 years): about 1.5 hours with the C
#                 engine; writes results/ (replacing the tables)
#   Plain R takes about 100 times longer: compile the C engine first.
#
#   Writes  E2/runs.csv.gz       one row per run: N at years 0-300
#           E2/occupancy.csv.gz  share of years 180-200 occupied
#           E2/rho50.csv         rho_50, width and bootstrap bounds
#           E2/pullback.csv.gz   poor and good cells, year by year
#           E2/long_runs.csv     the 2,000-year runs
#           E2/window_theory.csv exact one-cell curves
# ======================================================================

source("Load.R")
quick = TRUE

S = if (quick) {
  list(n_coarse = 2, n_fine = 3, refine = FALSE, n_boot = 50,
       n_pull = 3, long_rains = 2, n_long = 3, long_years = 300,
       Omega = c(1, 2), f_step = 0.02)
} else {
  list(n_coarse = 10, n_fine = 20, refine = TRUE, n_boot = 200,
       n_pull = 20, long_rains = 10, n_long = 20, long_years = 2000,
       Omega = c(1, 2, 5, 10), f_step = 0.01)
}
YEARS = 200                              # main window
EXTRA = 100                              # continuation to year 300
OCC = c(180, 200)                        # occupancy window
old = options(hh.results = if (quick) "results_quick" else "results")
p = Set_params()
Bd = Read_boards()
cat("engine:", Engine_choose(), "\n")
t0 = Sys.time()

# ---- start boards B_r (E1), as one C x n matrix per board ------------
start_file = file.path(Res_dir(), "E1", "start_boards.csv.gz")
sb = if (file.exists(start_file)) {
  read.csv(start_file)
} else {
  cat("no", start_file, "- using the saved start boards\n")
  read.csv(file.path("results", "E1", "start_boards.csv.gz"))
}
Start_matrix = function(b, map) {
  s = sb[sb$board == b, ]
  rs = sort(unique(s$run))
  return(sapply(rs, function(r) {
    g = matrix(0, 20, 20)
    t = s[s$run == r, ]
    g[cbind(t$row, t$col)] = t$n
    return(Grid_to_cells(g, map))
  }))
}

Key = function(rain) round(rain * 1000)  # rain as a whole number

out_runs = list()
out_occ = list()
out_rho50 = list()
out_pull = list()
for (b in 1:7) {
  map = Build_map(Bd[[b]])
  Bst = Start_matrix(b, map)             # C x n_start
  n_start = ncol(Bst)
  runs = new.env()                       # key "rain run" -> one run

  # simulate: 200-year runs at one rain from start boards rs
  simulate = function(rain, rs) {
    rs = rs[rs <= n_start]
    rs = rs[!(paste(Key(rain), rs) %in% ls(runs))]
    if (length(rs) == 0) return(invisible(NULL))
    m = Run_many(p, map, rain, Bst[, rs, drop = FALSE],
                 n_runs = length(rs), t_max = YEARS, dt_rec = 1,
                 occ_window = OCC,
                 seeds = sapply(rs, function(r) {
                   return(Seed_of("E2", b, Key(rain), r))
                 }))
    for (j in seq_along(rs)) {
      assign(paste(Key(rain), rs[j]), list(
        rain = rain, run = rs[j], t_ext = m$t_ext[j], N = m$N[, j],
        end = m$end[, j], occ = m$occ_time[, j] / diff(OCC),
        t_ext300 = NA_real_, N300 = rep(NA_integer_, EXTRA)),
        envir = runs)
    }
    return(invisible(NULL))
  }
  at_rain = function(rain) {
    keys = ls(runs)
    keys = keys[sapply(keys, function(k) runs[[k]]$rain == rain)]
    return(keys)
  }
  n_dead = function(rain) {
    return(sum(sapply(at_rain(rain), function(k) {
      return(!is.na(runs[[k]]$t_ext))
    })))
  }
  # continue: survivors at one rain run 100 more years (to year 300);
  # returns the share of runs at that rain dead by year 300
  continue = function(rain) {
    keys = at_rain(rain)
    alive = keys[sapply(keys, function(k) {
      return(is.na(runs[[k]]$t_ext) && all(is.na(runs[[k]]$N300)))
    })]
    if (length(alive) > 0) {
      n0 = sapply(alive, function(k) runs[[k]]$end)
      m = Run_many(p, map, rain, n0, n_runs = length(alive),
                   t_max = EXTRA, dt_rec = 1,
                   seeds = sapply(alive, function(k) {
                     return(Seed_of("E2 cont", b, Key(rain),
                                    runs[[k]]$run))
                   }))
      for (j in seq_along(alive)) {
        r = runs[[alive[j]]]
        r$t_ext300 = YEARS + m$t_ext[j]
        r$N300 = m$N[-1, j]              # years 201-300
        assign(alive[j], r, envir = runs)
      }
    }
    dead = sapply(keys, function(k) {
      return(!is.na(runs[[k]]$t_ext) || !is.na(runs[[k]]$t_ext300))
    })
    return(mean(dead))
  }
  all_rains = function() {
    return(sort(unique(sapply(ls(runs), function(k) runs[[k]]$rain))))
  }

  # 1. coarse pass, widened until the ends are reached
  coarse = round(seq(0.4, 1.2, by = 0.1), 3)
  for (v in coarse) simulate(v, seq_len(S$n_coarse))
  lo = min(coarse)
  while (n_dead(lo) < 0.9 * S$n_coarse && lo > 0.25) {
    lo = round(lo - 0.1, 3)
    coarse = c(lo, coarse)
    simulate(lo, seq_len(S$n_coarse))
  }
  hi = max(coarse)
  while (n_dead(hi) > 0.1 * S$n_coarse && hi < 1.35) {
    hi = round(hi + 0.1, 3)
    coarse = c(coarse, hi)
    simulate(hi, seq_len(S$n_coarse))
  }
  cc = sapply(coarse, n_dead)

  # 2. fine pass over the switch
  mid = coarse[cc >= 1 & cc <= S$n_coarse - 1]
  if (length(mid) == 0) {                # the 10 -> 0 jump
    k = which(cc[-length(cc)] == S$n_coarse & cc[-1] == 0)[1]
    mid = coarse[c(k, k + 1)]
  }
  f_lo = max(min(mid) - 0.1, min(coarse))
  f_hi = min(max(mid) + 0.1, max(coarse))
  for (v in round(seq(f_lo, f_hi + 1e-9, by = 0.025), 3)) {
    simulate(v, seq_len(S$n_fine))
  }

  # 3. continuations to year 300
  table_now = function() {
    keys = ls(runs)
    return(data.frame(
      rain = sapply(keys, function(k) runs[[k]]$rain),
      dead = sapply(keys, function(k) !is.na(runs[[k]]$t_ext))))
  }
  tn = table_now()
  r50 = Fit_logistic(tn$rain, tn$dead)[["rho50"]]
  top = min(all_rains())
  share = 1
  for (v in all_rains()) {
    if (v <= r50 + 0.15 + 1e-9) {
      share = continue(v)
      top = v
    }
  }
  higher = all_rains()[all_rains() > top]
  while (share > 0.1 && length(higher) > 0) {
    top = higher[1]
    higher = higher[-1]
    share = continue(top)
  }

  # 4. refine: every 0.005 inside the switch
  if (S$refine) {
    rains = all_rains()
    dead_share = sapply(rains, function(v) n_dead(v) /
                          length(at_rain(v)))
    all_dead = rains[dead_share == 1 & rains < 1.0]
    lo = if (length(all_dead) > 0) max(all_dead) else min(rains)
    none_dead = rains[dead_share == 0 & rains > lo]
    hi = if (length(none_dead) > 0) min(none_dead) else max(rains)
    for (v in round(seq(lo + 0.005, hi - 1e-9, by = 0.005), 3)) {
      if (any(abs(rains - v) < 1e-9)) next
      simulate(v, seq_len(S$n_fine))
      continue(v)
    }
  }

  # ---- one row per run -----------------------------------------------
  keys = ls(runs)
  keys = keys[order(sapply(keys, function(k) runs[[k]]$rain),
                    sapply(keys, function(k) runs[[k]]$run))]
  pos = (map$row - 1) * 20 + map$col     # ties broken by position
  poor = order(map$q, pos)[1:200]        # the 200 poorest cells
  quad = which((map$row <= 10 & map$col > 10) |
                 (map$row > 10 & map$col <= 10))   # poor quadrants
  rows = lapply(keys, function(k) {
    r = runs[[k]]
    N = c(r$N, r$N300)
    names(N) = paste0("N_", 0:(YEARS + EXTRA))
    return(data.frame(board = b, rain = r$rain, run = r$run,
                      t_ext = r$t_ext, t_ext300 = r$t_ext300,
                      cont = as.integer(!all(is.na(r$N300))),
                      top300 = top, poor200 = sum(r$end[poor]),
                      quad200 = sum(r$end[quad]), t(N)))
  })
  d = do.call(rbind, rows)
  out_runs[[b]] = d

  # ---- occupancy in years 180-200 of the runs alive at year 200 ------
  for (v in sort(unique(d$rain))) {
    alive = keys[sapply(keys, function(k) {
      return(runs[[k]]$rain == v && is.na(runs[[k]]$t_ext))
    })]
    if (length(alive) == 0) next
    occ = rowMeans(sapply(alive, function(k) runs[[k]]$occ))
    g = Cells_to_grid(occ, map)
    out_occ[[length(out_occ) + 1]] = data.frame(
      board = b, rain = v, n_alive = length(alive),
      row = as.vector(row(g)), col = as.vector(col(g)),
      occ = round(as.vector(g), 4))
  }

  # ---- rho_50(T) with bootstrap intervals ----------------------------
  for (T in c(50, 100, 150, 200, 300)) {
    dd = E2_dead(d, T)
    u = dd$usable
    fit = Fit_logistic(d$rain[u], dd$dead[u])
    ci = Boot_rho50(d$rain[u], dd$dead[u], d$run[u], n_boot = S$n_boot,
                    seed = 1)
    out_rho50[[length(out_rho50) + 1]] = data.frame(
      board = b, window = T, rho50 = fit[["rho50"]],
      width = fit[["width"]], lo = ci[["lo"]], hi = ci[["hi"]])
  }

  # ---- pullback: poor and good cells at the wet and dry edges --------
  f200 = Fit_logistic(d$rain, E2_dead(d, 200)$dead)
  edge = c(wet = f200[["rho50"]] + f200[["width"]] * log(9),
           dry = f200[["rho50"]] - f200[["width"]] * log(9))
  groups = list(poor = poor, good = setdiff(seq_len(map$C), poor))
  for (e in names(edge)) {
    rs = seq_len(min(S$n_pull, n_start))
    m = Run_many(p, map, edge[[e]], Bst[, rs, drop = FALSE],
                 n_runs = length(rs), t_max = YEARS, dt_rec = 1,
                 groups = groups,
                 seeds = sapply(rs, function(r) {
                   return(Seed_of("E2 pull", b, e, r))
                 }))
    for (g in c("poor", "good")) {
      Y = t(m$G[, g, ])                  # runs x 201
      colnames(Y) = paste0("Y_", 0:YEARS)
      out_pull[[length(out_pull) + 1]] = data.frame(
        board = b, edge = e, rain = edge[[e]], run = rs,
        t_ext = m$t_ext, group = g, Y)
    }
  }
  cat(sprintf(paste0("%-15s %d runs, rho50(200) %.3f, followed to ",
                     "300 up to rain %.3f  [%.0f s]\n"),
              Board_names[b], nrow(d), f200[["rho50"]], top,
              as.numeric(Sys.time() - t0, units = "secs")))
}
Res_write(do.call(rbind, out_runs), "E2/runs", gzip = TRUE)
Res_write(do.call(rbind, out_occ), "E2/occupancy", gzip = TRUE)
Res_write(do.call(rbind, out_rho50), "E2/rho50")
Res_write(do.call(rbind, out_pull), "E2/pullback", gzip = TRUE)

# ---- long runs: Real and Shuffled for 2,000 years --------------------
#   rains every 0.005 from rho_50(200), rounded to that grid, upward
long = list()
r50_tab = do.call(rbind, out_rho50)
for (b in 1:2) {
  map = Build_map(Bd[[b]])
  Bst = Start_matrix(b, map)
  r50 = r50_tab$rho50[r50_tab$board == b & r50_tab$window == 200]
  rains = round(round(r50 * 200) / 200 + 0.005 * (seq_len(S$long_rains)
                                                  - 1), 3)
  for (v in rains) {
    rs = seq_len(min(S$n_long, ncol(Bst)))
    m = Run_many(p, map, v, Bst[, rs, drop = FALSE],
                 n_runs = length(rs),
                 t_max = S$long_years, dt_rec = NULL,
                 seeds = sapply(rs, function(r) {
                   return(Seed_of("E2 long", b, Key(v), r))
                 }))
    long[[length(long) + 1]] = data.frame(board = b, rain = v, run = rs,
                                          t_ext = m$t_ext)
  }
}
Res_write(do.call(rbind, long), "E2/long_runs")

# ---- the exact one-cell curves (no simulation) -----------------------
#   K = 20 Omega, h = 2 Omega; start at the settled count of food 1.0
Ts = c(10, 50, 200, 1e3, 1e4, 1e5, 1e6, 1e8)
fs = round(seq(0.40, 1.60, by = S$f_step), 2)
wt = list()
for (Om in S$Omega) {
  P = Exact_P_empty(fs, Ts, Omega = Om, p = p)
  fold = F_crit(Set_params(K = 20 * Om, h = 2 * Om))
  for (i in seq_along(Ts)) {
    wt[[length(wt) + 1]] = data.frame(K = 20 * Om, fold = fold,
                                      window = Ts[i], f = fs,
                                      P_empty = P[i, ])
  }
}
Res_write(do.call(rbind, wt), "E2/window_theory")

cat(sprintf("done in %.1f minutes; tables in %s/E2/\n",
            as.numeric(Sys.time() - t0, units = "mins"), Res_dir()))
options(old)
