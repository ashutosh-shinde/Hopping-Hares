# ======================================================================
# Run_E4.R   E4. Coming back: how big must a released colony be?
# ======================================================================
#
#   Outcome page, section 6. The Real board, empty apart from the
#   released colony: a k x k square of cells with 14 hares each (about
#   the settled count of a good cell at rain 1.0), placed where the
#   mean grass over the square is highest.
#     Arm A, colony size: k = 1, 2, 3, 4, 5, 6, 8 (14 to 896 hares)
#     Arm B, one big or several small, 224 hares each: one 4x4
#            colony, four 2x2 colonies (the best spot inside each
#            quadrant), sixteen single cells (the best cell inside
#            each 5 x 5 part)
#   Rains 1.0, 0.9, 0.85 and 0.785 (Real's rho_50 at 200 years);
#   50 runs of 50 years each. Success: hares alive at year 50.
#
#   Run from the Hopping_hare_final folder:
#     source("run/Run_E4.R")
#   quick = TRUE  5 runs per design and rain: under a minute with the
#                 C engine (about 30 minutes in plain R); writes
#                 results_quick/
#   quick = FALSE the full design (2,000 runs): about 3 minutes with
#                 the C engine; writes results/ (replacing the tables)
#
#   Writes  E4/runs.csv.gz   N at years 0-50 and occupied cells at 50
#           E4/sites.csv     where each colony sits
# ======================================================================

source("Load.R")
quick = TRUE

S = list(runs = if (quick) 5 else 50, years = 50, per_cell = 14,
         rains = c(1.0, 0.9, 0.85, 0.785))
old = options(hh.results = if (quick) "results_quick" else "results")
p = Set_params()
q2 = Read_boards()[[1]]                  # the Real board, 20 x 20
map = Build_map(q2)
cat("engine:", Engine_choose(), "\n")
t0 = Sys.time()

# ----------------------------------------------------------------------
# Best_block: the k x k square with the highest mean grass inside the
# rows and columns given (the first best one, scanning row by row)
#
#   returns  c(row, col, mean q): the square's top-left cell
# ----------------------------------------------------------------------
Best_block = function(k, rows = 1:20, cols = 1:20) {
  best = c(NA, NA, -1)
  for (r in rows[rows + k - 1 <= max(rows)]) {
    for (c in cols[cols + k - 1 <= max(cols)]) {
      v = mean(q2[r:(r + k - 1), c:(c + k - 1)])
      if (v > best[3]) best = c(r, c, v)
    }
  }
  return(best)
}

# ----------------------------------------------------------------------
# Colony: the starting board of one design and its sites
#
#   k       colony side
#   parts   1, 4 or 16: one colony in each of parts equal squares of
#           the board (the best spot inside each)
#
#   returns  list(n0 = hares per cell (map order), sites = data frame)
# ----------------------------------------------------------------------
Colony = function(k, parts = 1) {
  g = round(sqrt(parts))
  h = 20 %/% g
  n0 = matrix(0, 20, 20)
  sites = data.frame()
  for (i in seq_len(g)) {
    for (j in seq_len(g)) {
      bb = Best_block(k, (i - 1) * h + 1:h, (j - 1) * h + 1:h)
      n0[bb[1]:(bb[1] + k - 1), bb[2]:(bb[2] + k - 1)] = S$per_cell
      sites = rbind(sites, data.frame(row = bb[1], col = bb[2], k = k,
                                      mean_q = round(bb[3], 4)))
    }
  }
  return(list(n0 = Grid_to_cells(n0, map), sites = sites))
}

designs = rbind(
  data.frame(arm = "A", design = sprintf("%dx%d", c(1:6, 8), c(1:6, 8)),
             k = c(1:6, 8), parts = 1),
  data.frame(arm = "B", design = c("1 x 4x4", "4 x 2x2", "16 x 1x1"),
             k = c(4, 2, 1), parts = c(1, 4, 16)))

out_runs = list()
out_sites = list()
for (i in seq_len(nrow(designs))) {
  ds = designs[i, ]
  colony = Colony(ds$k, ds$parts)
  out_sites[[i]] = data.frame(arm = ds$arm, design = ds$design,
                              site = seq_len(nrow(colony$sites)),
                              colony$sites)
  for (rain in S$rains) {
    m = Run_many(p, map, rain, colony$n0, n_runs = S$runs,
                 t_max = S$years, dt_rec = 1,
                 seeds = sapply(seq_len(S$runs), function(r) {
                   return(Seed_of("E4", ds$arm, ds$design, rain, r))
                 }))
    N = t(m$N)                           # runs x 51
    colnames(N) = paste0("N_", 0:S$years)
    out_runs[[length(out_runs) + 1]] = data.frame(
      arm = ds$arm, design = ds$design, rain = rain,
      run = seq_len(S$runs), t_ext = m$t_ext,
      occupied50 = colSums(m$end > 0), N)
  }
  cat(sprintf("%s %-9s done  [%.0f s]\n", ds$arm, ds$design,
              as.numeric(Sys.time() - t0, units = "secs")))
}
Res_write(do.call(rbind, out_runs), "E4/runs", gzip = TRUE)
Res_write(do.call(rbind, out_sites), "E4/sites")
cat(sprintf("done in %.1f minutes; tables in %s/E4/\n",
            as.numeric(Sys.time() - t0, units = "mins"), Res_dir()))
options(old)
