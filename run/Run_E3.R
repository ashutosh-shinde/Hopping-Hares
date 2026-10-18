# ======================================================================
# Run_E3.R   E3. Falling rain: does the speed of the drought matter?
# ======================================================================
#
#   Outcome page, section 5. Per board and ramp length D (50, 150,
#   450, 1,350 and 4,050 years): rain falls along a straight line from
#   1.0 to 0.4 over D years in monthly steps (Rho_ramp(), each month
#   at the line's middle value), then holds 0.4 for up to 50 years.
#   20 runs from E1's start boards B_1 ... B_20.
#   Saved per run: N every month, the board each time the line passes
#   1.0, 0.995, 0.990, ..., 0.4, the first time N < 10% of its start
#   (t_frac) and the last hare (t_ext).
#   The deterministic ramps (E3's dashed black line) come from
#   run/Run_deterministic.R.
#
#   Needs E1's start boards. Run from the Hopping_hare_final folder:
#     source("run/Run_E3.R")
#   quick = TRUE  3 runs, ramps of 50 and 450 years: under a minute
#                 with the C engine; writes results_quick/
#   quick = FALSE the full design (700 runs): about 30 minutes with the
#                 C engine; writes results/ (replacing the tables)
#
#   Writes  E3/runs.csv            t_frac and t_ext of every run
#           E3/mean_N.csv.gz       mean N against the current rain
#           E3/empty_rain.csv.gz   rain at which each cell empties
# ======================================================================

source("Load.R")
quick = TRUE

S = list(runs = if (quick) 3 else 20,
         durs = if (quick) c(50, 450) else c(50, 150, 450, 1350, 4050),
         hold = 50, from = 1.0, to = 0.4)
old = options(hh.results = if (quick) "results_quick" else "results")
p = Set_params()
Bd = Read_boards()
cat("engine:", Engine_choose(), "\n")
t0 = Sys.time()

start_file = file.path(Res_dir(), "E1", "start_boards.csv.gz")
sb = if (file.exists(start_file)) {
  read.csv(start_file)
} else {
  cat("no", start_file, "- using the saved start boards\n")
  read.csv(file.path("results", "E1", "start_boards.csv.gz"))
}

grid_rain = seq(S$from, S$to, by = -0.005)       # saved boards
bin_rain = round(seq(S$from, S$to, by = -0.0025), 4)  # mean N bins
out_runs = list()
out_mean = list()
out_empty = list()
for (b in 1:7) {
  map = Build_map(Bd[[b]])
  s = sb[sb$board == b, ]
  rs = seq_len(min(S$runs, max(s$run)))
  Bst = sapply(rs, function(r) {
    g = matrix(0, 20, 20)
    t = s[s$run == r, ]
    g[cbind(t$row, t$col)] = t$n
    return(Grid_to_cells(g, map))
  })
  for (D in S$durs) {
    r = (S$from - S$to) / D              # rain lost per year
    rho = c(Rho_ramp(S$from, S$to, D, dt_rain = 1 / 12),
            Rho_const(S$to, S$hold, dt_rain = 1 / 12))
    m = Run_many(p, map, rho, Bst, n_runs = length(rs),
                 t_max = D + S$hold, dt_rain = 1 / 12, dt_rec = 1 / 12,
                 snap_t = (S$from - grid_rain) / r, frac = 0.1,
                 seeds = sapply(rs, function(k) {
                   return(Seed_of("E3", b, D, k))
                 }))
    out_runs[[length(out_runs) + 1]] = data.frame(
      board = b, dur = D, run = rs, t_frac = m$t_frac, t_ext = m$t_ext)

    # mean N (collapsed runs = 0) inside rain bins of 0.0025
    rain_t = S$from - r * m$t            # the line at each month
    N_mean = rowMeans(m$N)
    for (v in bin_rain) {
      sel = abs(rain_t - v) <= 0.00125 + 1e-12
      if (any(sel)) {
        out_mean[[length(out_mean) + 1]] = data.frame(
          board = b, dur = D, rain = v, N_mean = mean(N_mean[sel]))
      }
    }

    # rain at which each cell is empty in more than half the runs
    occ = apply(m$snap > 0, c(1, 2), mean)     # n_snap x C
    empty = apply(occ, 2, function(o) {
      k = which(o < 0.5)
      return(if (length(k) > 0) grid_rain[k[1]] else S$to)
    })
    g = Cells_to_grid(empty, map)
    out_empty[[length(out_empty) + 1]] = data.frame(
      board = b, dur = D, row = as.vector(row(g)),
      col = as.vector(col(g)), rain_empty = as.vector(g))
  }
  cat(sprintf("%-15s done  [%.0f s]\n", Board_names[b],
              as.numeric(Sys.time() - t0, units = "secs")))
}
Res_write(do.call(rbind, out_runs), "E3/runs")
Res_write(do.call(rbind, out_mean), "E3/mean_N", gzip = TRUE)
Res_write(do.call(rbind, out_empty), "E3/empty_rain", gzip = TRUE)
cat(sprintf("done in %.1f minutes; tables in %s/E3/\n",
            as.numeric(Sys.time() - t0, units = "mins"), Res_dir()))
options(old)
