# ======================================================================
# Plot_deterministic.R   figures of the deterministic board
# ======================================================================
# Plotting code written with AI assistance (Claude).
#
#   Outcome page, section 2 ("What the deterministic board tells
#   us"). Data: results/deterministic/ (run/Run_deterministic.R),
#   data/boards.csv, and E2's runs for the simulator's means.
#
#   Plot_boards()          the seven boards (sections 2.2, 3, 4)
#   Plot_one_cell_births() births per hare (2.1, A)
#   Plot_one_cell_states() settled states of one cell (2.1, B)
#   Plot_branches()        hares on the settled board as rain falls
#   Plot_staircase()       the board just above each fold (2.3)
#   Plot_recovery()        recovery time along the branch (2.4)
#   Plot_hop_sweep()       fold rain against hop rate (2.5)
#   nicknames at the end of the file
# ======================================================================

# ----------------------------------------------------------------------
# Det_rows, Det_folds: the saved branch and folds of board b
#
#   returns  data frames (rho, N, occupied, lambda) and (fold, rho,
#            N_before, N_after)
# ----------------------------------------------------------------------
Det_rows = function(b) {
  d = Res("deterministic/branch")
  return(d[d$board == b, ])
}

Det_folds = function(b) {
  d = Res("deterministic/folds")
  return(d[d$board == b, ])
}

# ----------------------------------------------------------------------
# Det_branch_lines: the settled branch of board b as lines, broken at
# each fold, with the drop at each fold dotted
#
#   b      board number
#   col    line colour
#   lwd    line width
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Det_branch_lines = function(b, col, lwd = 1.8) {
  r = Det_rows(b)
  f = Det_folds(b)
  # Det_follow() saves one row at each fold's rain (the board just
  # before the drop); break the line right after those rows. This
  # works for any rain step, unlike a cut on the size of the drop.
  cut = which(sapply(r$rho, function(v) {
    return(any(abs(v - f$rho) < 1e-6))
  }))
  seg = split(seq_len(nrow(r)), findInterval(seq_len(nrow(r)),
                                             cut + 1))
  for (s in seg) lines(r$rho[s], r$N[s], col = col, lwd = lwd)
  segments(f$rho, f$N_before, f$rho, f$N_after, col = col, lwd = 0.9,
           lty = 3)
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Sim_mean_N: the simulator's mean hares at a few fixed rains
#
#   Mean N over years 100-200 of E2's runs still alive at year 200
#   (the hollow markers of section 2.2).
#
#   b       board number
#   rains   rains to report (those E2 ran on this board)
#
#   returns  data frame: rain, N
# ----------------------------------------------------------------------
Sim_mean_N = function(b, rains = c(1.0, 0.95, 0.9, 0.875, 0.85)) {
  d = Res("E2/runs")
  d = d[d$board == b & is.na(d$t_ext), ]
  N = Wide_years(d, "N_")[, as.character(100:200), drop = FALSE]
  out = data.frame(rain = numeric(0), N = numeric(0))
  for (v in rains) {
    sel = abs(d$rain - v) < 1e-9
    if (any(sel)) out = rbind(out, data.frame(rain = v,
                                              N = mean(N[sel, ])))
  }
  return(out)
}

# ----------------------------------------------------------------------
# Plot_boards: the seven boards side by side (darker = better grass)
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_boards = function() {
  B = Read_boards()
  op = par(mfrow = c(1, 7), mar = c(0.3, 0.4, 2.0, 0.4), oma = c(0, 0,
                                                                  0, 0))
  on.exit(par(op))
  for (b in 1:7) {
    Map_panel(B[[b]], c(0.3, 0.96), Pal$q, main = Board_label[b],
              cex_main = 1.15)
  }
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Plot_one_cell_births: births per hare against hares in one cell, for
# four foods, with deaths per hare as a flat line (section 2.1, A)
#
#   The hill: b0 f phi(n) s(n). Where it crosses mu = 1 a cell can
#   rest; at food f_c = 0.498 the hill only touches the line (the
#   one-cell fold).
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_one_cell_births = function() {
  p = Set_params()
  f_c = F_crit(p)
  n = seq(1, 20, length.out = 600)
  g = (n - 1) / (n - 1 + p$h) * (1 - n / p$K)
  foods = c(1.0, 0.662, f_c, 0.40)
  cols = c(Board_col[1], Board_col[3], Board_col[2], Col$muted)
  op = par(mar = c(3.4, 4.2, 2.2, 1))
  on.exit(par(op))
  Panel(c(1, 20), c(0, 2.6),
        main = "A. Births per hare against crowding",
        xlab = "hares in the cell, n",
        ylab = "births per hare per year")
  for (k in seq_along(foods)) {
    lines(n, p$b0 * foods[k] * g, col = cols[k], lwd = 2.2)
  }
  abline(h = p$mu, col = Col$ink, lwd = 1.3, lty = 2)
  text(19.8, 1.08, expression("deaths per hare (" * mu * " = 1)"),
       adj = c(1, 0), cex = 0.85, col = Col$ink)
  legend("topright", sprintf("food %.3f", foods), col = cols, lwd = 2.2,
         bty = "n", cex = 0.9, inset = c(0.02, 0.06), seg.len = 2.2)
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Plot_one_cell_states: one cell's settled states against food
# (section 2.1, B)
#
#   Blue: the settled count x_plus (stable); dashed orange: the
#   tipping count x_minus (unstable); black: the empty cell (always
#   stable). The two counts meet and vanish at the fold f_c.
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_one_cell_states = function() {
  p = Set_params()
  f_c = F_crit(p)
  fs = c(f_c + 1e-9, seq(f_c + 0.0005, 1.4, length.out = 400))
  st = sapply(fs, function(f) Cell_states(f, p))     # 2 x length(fs)
  op = par(mar = c(3.4, 4.2, 2.2, 1))
  on.exit(par(op))
  Panel(c(0.3, 1.4), c(-0.8, 20),
        main = "B. One cell's settled states against food",
        xlab = "food in the cell, f = rain \u00d7 grass",
        ylab = "hares, n")
  lines(c(0.3, 1.4), c(0, 0), col = Col$ink, lwd = 3)
  lines(fs, st[2, ], col = Board_col[1], lwd = 2.2)
  lines(fs, st[1, ], col = Board_col[2], lwd = 2.2, lty = 2)
  abline(v = f_c, col = Col$ink2, lwd = 1, lty = 3)
  text(f_c + 0.01, 18.6, sprintf("fold\nfood %.3f", f_c), adj = c(0, 1),
       cex = 0.85, col = Col$ink2)
  legend("right", c("settled count (stable)",
                    "tipping count (unstable)",
                    "empty (always stable)"),
         col = c(Board_col[1], Board_col[2], Col$ink),
         lwd = c(2.2, 2.2, 3),
         lty = c(1, 2, 1), bty = "n", cex = 0.9, seg.len = 2.2)
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Plot_branches: hares on the settled board as rain falls, wet on the
# left (section 2.2)
#
#   Lines: the settled branch of each board; filled markers: its
#   folds; dotted lines: the drop at each fold; hollow markers: the
#   simulator's mean N at the same rain (always lower; section 2.6).
#
#   group   "fine" (Shuffled, Checker 1x1, Blocks 2x2 and 4x4; panel
#           A) or "clumped" (Real and the 10x10 boards; panel B)
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_branches = function(group = "fine") {
  fine = group == "fine"
  boards = if (fine) c(2, 3, 4, 5) else c(1, 6, 7)
  main = if (fine) "A. Fine-grained boards" else "B. Clumped boards"
  op = par(mar = c(3.4, 4.6, 2.2, 1))
  on.exit(par(op))
  Panel(c(1.2, 0.68), c(0, 5000), main = main, xlab = expression(
    "rain " * rho), ylab = "hares on the board, N")
  for (b in boards) {
    Det_branch_lines(b, Board_col[b])
    f = Det_folds(b)
    Board_points(f$rho, f$N_before, b, cex = 1.15)
    s = Sim_mean_N(b)
    Board_points(s$rain, s$N, b, cex = 1.15, hollow = TRUE)
  }
  legend("topright", Board_label[boards], col = Board_col[boards],
         pch = Board_pch[boards], pt.bg = Board_col[boards], lwd = 2,
         pt.lwd = ifelse(Board_pch[boards] >= 21, 0.8, 2),
         bty = "n", cex = 0.9, seg.len = 2)
  if (fine) {
    text(1.2, 200, paste0("filled markers: folds\n",
                           "hollow markers: simulator mean"),
         adj = c(0, 1), cex = 0.7, col = Col$ink2)
  }
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Plot_staircase: the board just above each fold on the Real and
# 10x10 random boards (section 2.3)
#
#   First column: the grass. Then one map per fold: hares on the
#   settled board just above it (white = none, dark orange = 16);
#   the purple line circles the cells where that fold's collapse
#   starts (the slowest eigenvector above half its peak).
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_staircase = function(boards = c(1, 6)) {
  B = Read_boards()
  maps = Res("deterministic/fold_maps")
  nf = max(sapply(boards, function(b) nrow(Det_folds(b))))
  op = par(mfrow = c(length(boards), nf + 1),
           mar = c(0.3, 0.4, 2.0, 0.4), oma = c(0, 2.2, 0, 0))
  on.exit(par(op))
  for (r in seq_along(boards)) {
    b = boards[r]
    f = Det_folds(b)
    Map_panel(B[[b]], c(0.3, 0.96), Pal$q, main = "grass", frame = TRUE)
    mtext(Board_label[b], side = 2, line = 0.4, col = Col$ink,
          cex = 0.8)
    for (k in seq_len(nf)) {
      if (k > nrow(f)) {
        plot.new()
        next
      }
      m = maps[maps$board == b & maps$fold == k, ]
      x = matrix(NA_real_, 20, 20)
      x[cbind(m$row, m$col)] = m$n_before
      pat = matrix(NA_real_, 20, 20)
      pat[cbind(m$row, m$col)] = m$pattern
      Map_panel(x, c(0, 16), Pal$n, frame = TRUE,
                main = sprintf("fold %d \u00b7 %.3f", k, f$rho[k]))
      contour(1:20, 1:20, t(pat[20:1, ]), levels = 0.5, add = TRUE,
              drawlabels = FALSE, col = Board_col[7], lwd = 1.4)
    }
  }
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Plot_recovery: recovery time 1 / |lambda_1| along each settled
# branch (section 2.4). Each spike is a fold: recovery slows down
# before every fold (critical slowing down).
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_recovery = function() {
  op = par(mar = c(3.4, 4.6, 2.2, 1))
  on.exit(par(op))
  y_at = 10^(0:3)
  Panel(c(1.2, 0.68), c(0.6, 1100), log = "y",
        main = "Recovery slows down before every fold",
        xlab = expression("rain " * rho),
        ylab = "recovery time, years (1 / slowest rate)",
        y_at = y_at, y_lab = Log_labels(y_at))
  for (b in 1:7) {
    r = Det_rows(b)
    lines(r$rho, 1 / abs(r$lambda), col = Board_col[b], lwd = 1.5)
  }
  legend("topleft", Board_label, col = Board_col, lwd = 1.8, ncol = 2,
         bty = "n", cex = 0.85, seg.len = 2)
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Plot_hop_sweep: the Shuffled board's folds against the hop rate m
# (section 2.5)
#
#   Slow hoppers keep their best patches longest; fast hoppers pool
#   into one population; m = 10 sits near the worst point. Dotted
#   lines: no hopping (the best cell dies at f_c / max q) and the
#   fast-hopping limit (Fast_hop_limit()).
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_hop_sweep = function() {
  h = Res("deterministic/hop_sweep")
  lim = Res("deterministic/limits")
  fast = lim$value[lim$name == "fast_hop_limit"]
  slow = lim$value[lim$name == "no_hop_last_cell"]
  op = par(mar = c(3.4, 4.2, 2.2, 1))
  on.exit(par(op))
  x_at = 10^(0:2)
  Panel(range(h$m), c(0.5, 0.97), log = "x",
        main = "Shuffled board: how hopping sets the fold",
        xlab = "hop tries per hare per year, m (default 10)",
        ylab = "rain at the fold", x_at = x_at,
        x_lab = Log_labels(x_at))
  abline(v = 10, col = Col$muted, lwd = 1)
  abline(h = c(fast, slow), col = Col$ink2, lwd = 0.9, lty = 3)
  text(12, fast - 0.012, sprintf("fast-hopping limit %.3f", fast),
       adj = c(0, 1), cex = 0.8, col = Col$ink2)
  text(max(h$m), slow + 0.008,
       sprintf("no hopping: last cell dies at %.3f", slow),
       adj = c(1, 0), cex = 0.8, col = Col$ink2)
  oc = Board_col[2]
  lines(h$m, h$final_fold, col = oc, lwd = 2.2)
  points(h$m, h$final_fold, pch = 22, bg = oc, col = oc, cex = 1.3)
  lines(h$m, h$first_fold, col = oc, lwd = 2.2, lty = 2)
  points(h$m, h$first_fold, pch = 22, bg = "white", col = oc, cex = 1.3,
         lwd = 1.5)
  legend("topright", c("final fold (board empty)",
                       "first fold (first region lost)"),
         col = oc, lwd = 2.2, lty = c(1, 2), pch = 22,
         pt.bg = c(oc, "white"), bty = "n", cex = 0.9, seg.len = 3)
  return(invisible(NULL))
}

# ---- nicknames: one line each, no input (README.Rmd calls these) -----
boards_plot = function() Plot_boards()
deterministic_one_cell_plot_1 = function() Plot_one_cell_births()
deterministic_one_cell_plot_2 = function() Plot_one_cell_states()
deterministic_fold_plot_1 = function() Plot_branches("fine")
deterministic_fold_plot_2 = function() Plot_branches("clumped")
deterministic_staircase_plot = function() Plot_staircase()
deterministic_recovery_plot = function() Plot_recovery()
deterministic_hop_plot = function() Plot_hop_sweep()
