# ======================================================================
# Run_deterministic.R   the deterministic board: folds of the seven
#                       boards, the hop sweep and the E3 ramps
# ======================================================================
#
#   Outcome page, section 2 (steps 5-9 of the build chart, 1.1); the
#   ramps feed E3's dashed black line (section 5). No random numbers:
#   every run of this file gives the same numbers.
#
#   Run from the Hopping_hare_final folder:
#     source("run/Run_deterministic.R")
#   quick = TRUE  three hop rates and two ramps: about 8 minutes;
#                 writes to results_quick/
#   quick = FALSE the full design of the outcome page: about 30
#                 minutes; writes to results/ (replacing the saved
#                 tables)
#   Keep the rain step at 0.0025: a coarser step can jump between
#   branches near folds that lie close together and report extra
#   folds. Checked against the saved tables: every fold of the seven
#   boards and every ramp agree to 1e-5 in rain. At hop rate m = 1
#   (dozens of tiny one-cell folds) the count can differ by one,
#   because two folds closer than the rain step merge.
#
#   Writes  deterministic/branch.csv     the settled branch per board
#           deterministic/folds.csv      every fold
#           deterministic/fold_maps.csv.gz  board and pattern at folds
#           deterministic/limits.csv     f_c and two hop limits
#           deterministic/hop_sweep.csv  folds against hop rate m
#           deterministic/ramps.csv      E3 ramps without luck
# ======================================================================

source("Load.R")
quick = TRUE

S = if (quick) {
  list(step = 0.0025, hop_m = c(10, 30, 100), durs = c(50, 450))
} else {
  list(step = 0.0025, hop_m = c(1, 3, 10, 30, 100, 400),
       durs = c(50, 150, 450, 1350, 4050))
}
old = options(hh.results = if (quick) "results_quick" else "results")
p = Set_params()
B = Read_boards()
t0 = Sys.time()

# Grid_long: per-cell values of a map as rows of (row, col, value)
Grid_long = function(x, map) {
  g = Cells_to_grid(x, map)
  return(data.frame(row = as.vector(row(g)), col = as.vector(col(g)),
                    value = as.vector(g)))
}

# ---- 1. the settled branch and its folds, board by board ------------
#   (steps 6-9: Newton, eigenvalues, continuation from rain 1.2 down
#   to 0.3 in steps of 0.0025, bisection at each fold)
branch = list()
folds = list()
maps = list()
for (b in 1:7) {
  map = Build_map(B[[b]])
  db = Det_board(map, p)
  fol = Det_follow(db, rho_hi = 1.2, rho_lo = 0.3, step = S$step)
  branch[[b]] = data.frame(board = b, fol$rows)
  folds[[b]] = data.frame(board = b, fol$folds[, c("fold", "rho",
                                                   "N_before",
                                                   "N_after")])
  for (k in seq_along(fol$maps)) {
    nb = Grid_long(fol$maps[[k]]$x_before, map)
    pt = Grid_long(fol$maps[[k]]$pattern, map)
    maps[[length(maps) + 1]] = data.frame(board = b, fold = k,
                                          row = nb$row, col = nb$col,
                                          n_before = nb$value,
                                          pattern = pt$value)
  }
  cat(sprintf("%-15s %d folds: %s  [%.0f s]\n", Board_names[b],
              nrow(fol$folds), paste(sprintf("%.4f", fol$folds$rho),
                                     collapse = " "),
              as.numeric(Sys.time() - t0, units = "secs")))
}
Res_write(do.call(rbind, branch), "deterministic/branch")
Res_write(do.call(rbind, folds), "deterministic/folds")
Res_write(do.call(rbind, maps), "deterministic/fold_maps", gzip = TRUE)

# ---- 2. limits: the one-cell fold and the two hop limits -------------
#   (section 2.5: no hopping, the best cell dies at f_c / max q; very
#   fast hopping pools the board into one population)
q_shuffled = Grid_to_cells(B[[2]], Build_map(B[[2]]))
lim = data.frame(name = c("f_c", "no_hop_last_cell", "fast_hop_limit"),
                 value = c(F_crit(p), F_crit(p) / max(q_shuffled),
                           Fast_hop_limit(q_shuffled, p)))
Res_write(lim, "deterministic/limits")

# ---- 3. the hop sweep: the Shuffled board at other hop rates --------
#   (section 2.5: rain 1.0 down to 0.45 in steps of 0.005)
hop = data.frame()
for (m in S$hop_m) {
  db = Det_board(Build_map(B[[2]]), Set_params(m = m))
  fol = Det_follow(db, rho_hi = 1.0, rho_lo = 0.45, step = 0.005)
  hop = rbind(hop, data.frame(m = m, n_folds = nrow(fol$folds),
                              first_fold = fol$folds$rho[1],
                              final_fold = fol$folds$rho[nrow(
                                fol$folds)]))
  cat(sprintf("hop sweep m = %g: %d folds, first %.4f, final %.4f\n",
              m, nrow(fol$folds), fol$folds$rho[1],
              fol$folds$rho[nrow(fol$folds)]))
}
Res_write(hop, "deterministic/hop_sweep")

# ---- 4. E3's ramps without luck --------------------------------------
#   (start settled at rain 1.0; rain falls in a straight line to 0.4
#   over D years; the rain when N first drops below 10% of its start)
ramps = data.frame()
for (b in 1:7) {
  db = Det_board(Build_map(B[[b]]), p)
  x0 = Det_start(db, 1.0)
  for (D in S$durs) {
    rho10 = Det_ramp(db, x0, D)$rho10
    ramps = rbind(ramps, data.frame(board = b, dur = D, rho10 = rho10))
  }
}
Res_write(ramps, "deterministic/ramps")

cat(sprintf("done in %.1f minutes; tables in %s/deterministic/\n",
            as.numeric(Sys.time() - t0, units = "mins"), Res_dir()))
options(old)
