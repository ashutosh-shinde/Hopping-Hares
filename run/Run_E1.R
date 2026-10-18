# ======================================================================
# Run_E1.R   E1. Start: how every run begins
# ======================================================================
#
#   Outcome page, section 3. Per board, rain (1.0 and 0.8) and start
#   (low, mid, high, random; plus 1 hare per cell at 0.8): 20 runs of
#   50 years from the same starting board, and the deterministic path
#   from the same start. The high start at rain 1.0 also gives the 20
#   start boards B_1 ... B_20 that E2, E3 and the pullback runs reuse:
#   run r's board at year 20.
#
#   Starts per cell, from the average equation at the run's rain
#   (Cell_states(): x_minus = tipping count, x_plus = settled count):
#     low     round(x_minus) + 2
#     mid     round((x_minus + x_plus) / 2) + 2
#     high    round(x_plus) + 2
#     random  0-20 hares per cell, one board per run
#     one  1 hare per cell
#   (2 hares where the cell's food lies below the fold.)
#
#   Run from the Hopping_hare_final folder:
#     source("run/Run_E1.R")
#   quick = TRUE  3 runs per start, all boards: about 1 minute with the
#                 C engine (an hour in plain R); writes results_quick/
#   quick = FALSE the full design (1,260 runs): about 6 minutes with
#                 the C engine; writes results/ (replacing the tables)
#
#   Writes  E1/runs.csv.gz          N at years 0-50 of every run
#           E1/det_paths.csv        the deterministic paths
#           E1/start_boards.csv.gz  B_1 ... B_n (board at year 20)
# ======================================================================

source("Load.R")
quick = TRUE

S = list(runs = if (quick) 3 else 20, years = 50, burn = 20,
         rains = c(1.0, 0.8))
old = options(hh.results = if (quick) "results_quick" else "results")
p = Set_params()
B = Read_boards()
cat("engine:", Engine_choose(), "\n")
t0 = Sys.time()

# ----------------------------------------------------------------------
# E1_start_board: one starting board (length C, map cell order)
#
#   kind   "low", "mid", "high", "random" or "one"
#   map    the board's map; rain   the run's rain; p   parameters
#
#   returns  hares per cell
# ----------------------------------------------------------------------
E1_start_board = function(kind, map, rain, p) {
  if (kind == "one") return(rep(1, map$C))
  if (kind == "random") return(sample(0:20, map$C, replace = TRUE))
  st = sapply(map$q, function(q) Cell_states(rain * q, p))  # 2 x C
  x = switch(kind, low = st[1, ], mid = (st[1, ] + st[2, ]) / 2,
             high = st[2, ])
  return(round(x) + 2)
}

runs = list()
dets = list()
starts = list()
for (b in 1:7) {
  map = Build_map(B[[b]])
  db = Det_board(map, p)
  for (rain in S$rains) {
    kinds = c("low", "mid", "high", "random",
              if (rain < 1) "one" else NULL)
    for (kind in kinds) {
      seeds = sapply(seq_len(S$runs), function(r) {
        return(Seed_of("E1", b, rain, kind, r))
      })
      if (kind == "random") {           # one random board per run
        n0 = sapply(seq_len(S$runs), function(r) {
          set.seed(Seed_of("E1 random board", b, rain, r))
          return(E1_start_board("random", map, rain, p))
        })
      } else {
        n0 = matrix(E1_start_board(kind, map, rain, p), map$C, S$runs)
      }
      many = Run_many(p, map, rain, n0, n_runs = S$runs,
                      t_max = S$years, dt_rec = 1, snap_t = S$burn,
                      seeds = seeds)
      N = t(many$N)                      # runs x 51
      colnames(N) = paste0("N_", 0:S$years)
      runs[[length(runs) + 1]] = data.frame(
        board = b, rain = rain, start = kind, run = seq_len(S$runs),
        t_ext = many$t_ext, N)
      # deterministic path from the same start (the mean random board)
      path = Det_path(db, rowMeans(n0), function(t) rain, 0:S$years)
      dp = matrix(path$N, 1)
      colnames(dp) = paste0("N_", 0:S$years)
      dets[[length(dets) + 1]] = data.frame(board = b, rain = rain,
                                            start = kind, dp)
      # B_r: the high start at rain 1.0, run r, at year 20
      if (rain == 1.0 && kind == "high") {
        for (r in seq_len(S$runs)) {
          g = Cells_to_grid(many$snap[1, , r], map)
          starts[[length(starts) + 1]] = data.frame(
            board = b, run = r, row = as.vector(row(g)),
            col = as.vector(col(g)), n = as.vector(g))
        }
      }
    }
  }
  cat(sprintf("%-15s done  [%.0f s]\n", Board_names[b],
              as.numeric(Sys.time() - t0, units = "secs")))
}
Res_write(do.call(rbind, runs), "E1/runs", gzip = TRUE)
Res_write(do.call(rbind, dets), "E1/det_paths")
Res_write(do.call(rbind, starts), "E1/start_boards", gzip = TRUE)
cat(sprintf("done in %.1f minutes; tables in %s/E1/\n",
            as.numeric(Sys.time() - t0, units = "mins"), Res_dir()))
options(old)
