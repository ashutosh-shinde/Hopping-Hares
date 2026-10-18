# ======================================================================
# Plot_E1.R   figures of E1 (start: how every run begins)
# ======================================================================
# Plotting code written with AI assistance (Claude).
#
#   Outcome page, section 3. Data: results/E1/ (run/Run_E1.R) and,
#   for the 200-year drift, E2's runs at rain 1.0.
#
#   Plot_E1_runs()        every run at rain 1.0: 4 starts x 7 boards
#   Plot_E1_starts()      mean of each start per board (rain 1.0 or
#                         0.8)
#   Plot_E1_drift_mean()  mean hares for 200 more years at rain 1.0
#   Plot_E1_drift_quad()  poor quadrants of the 10x10 boards
#   nicknames at the end of the file
# ======================================================================

E1_starts = c("low", "mid", "high", "random", "one")
E1_start_col = c(low = "#1baf7a", mid = "#eda100", high = "#2a78d6",
                 random = "#6250d6", one = "#eb6834")

# ----------------------------------------------------------------------
# E1_start_label: legend label of a start (x_minus and x_plus are the
# tipping and settled counts of each cell at the run's rain)
#
#   k   start name
#
#   returns  an expression
# ----------------------------------------------------------------------
E1_start_label = function(k) {
  lab = switch(k,
               low = expression("low   " * x[-""] + 2),
               mid = expression("mid   " * (x[-""] + x[+""]) / 2 + 2),
               high = expression("high   " * x[+""] + 2),
               random = expression("random   0" - 20),
               one = expression("1 per cell"))
  return(lab)
}

# ----------------------------------------------------------------------
# E1_runs: the saved runs of one board, rain and start
#
#   returns  list(N = runs x 51 matrix of hares at years 0-50, t_ext)
# ----------------------------------------------------------------------
E1_runs = function(b, rain, start) {
  d = Res("E1/runs")
  d = d[d$board == b & abs(d$rain - rain) < 1e-9 & d$start == start, ]
  return(list(N = Wide_years(d, "N_"), t_ext = d$t_ext))
}

# ----------------------------------------------------------------------
# E1_det: the deterministic path from the same start (years 0-50)
#
#   returns  length 51
# ----------------------------------------------------------------------
E1_det = function(b, rain, start) {
  d = Res("E1/det_paths")
  d = d[d$board == b & abs(d$rain - rain) < 1e-9 & d$start == start, ]
  return(as.vector(Wide_years(d, "N_")[1, ]))
}

# ----------------------------------------------------------------------
# E1_band: the settled band of board b at rain 1.0
#
#   level = mean N over years 30-50 of all 80 runs (four starts);
#   band = level +/- 2 SD of those yearly values (one run's wobble)
#
#   returns  c(level, sd)
# ----------------------------------------------------------------------
E1_band = function(b) {
  x = unlist(lapply(E1_starts[1:4], function(k) {
    return(E1_runs(b, 1.0, k)$N[, as.character(30:50)])
  }))
  return(c(level = mean(x), sd = sqrt(mean((x - mean(x))^2))))
}

# ----------------------------------------------------------------------
# Plot_E1_runs: every run at rain 1.0, rows = starts, columns = boards
#
#   Thin coloured lines: the 20 runs; black: their mean; dashed: the
#   deterministic path from the same start; pink band: settled level
#   +/- 2 SD of one run; grey line: year 20, where E2-E4 begin.
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_E1_runs = function() {
  op = par(no.readonly = TRUE)
  on.exit(par(op))
  layout(rbind(rep(1, 7), matrix(2:29, 4, 7, byrow = TRUE)),
         heights = c(0.35, 1, 1, 1, 1))
  par(oma = c(2.6, 4.6, 0, 0.5), cex = 0.8)
  par(mar = c(0, 0, 0, 0))
  n_runs = nrow(E1_runs(1, 1.0, "low")$N)    # 20 in the full design
  Legend_panel(legend = c("one run (coloured by start)",
                          sprintf("mean of the %d runs", n_runs),
                          "deterministic path",
                          "settled band (mean \u00b1 2 SD of one run)",
                          "year 20, where E2-E4 begin"),
               col = c(Col$muted, Col$ink, Col$ink2, Col$band,
                       Col$muted),
               lwd = c(0.8, 1.8, 1.2, 10, 1), lty = c(1, 1, 2, 1, 1),
               where = "center", horiz = TRUE)
  yrs = 0:50
  for (j in 1:4) {
    k = E1_starts[j]
    for (b in 1:7) {
      par(mar = c(0.6, 0.6, if (j == 1) 2.0 else 0.6, 0.6))
      bd = E1_band(b)
      Panel(c(0, 50), c(1500, 5000),
            main = if (j == 1) Board_label[b] else "",
            x_lab = j == 4, y_lab = b == 1)
      rect(-5, bd[1] - 2 * bd[2], 55, bd[1] + 2 * bd[2], col = Col$band,
           border = NA)
      abline(v = 20, col = Col$muted, lwd = 0.9)
      r = E1_runs(b, 1.0, k)
      for (i in seq_len(nrow(r$N))) {
        lines(yrs, r$N[i, ], col = adjustcolor(E1_start_col[k], 0.45),
              lwd = 0.7)
      }
      lines(yrs, colMeans(r$N), col = Col$ink, lwd = 1.6)
      lines(yrs, E1_det(b, 1.0, k), col = Col$ink2, lwd = 1.1, lty = 2)
      if (b == 1) {
        mtext(E1_start_label(k), side = 2, line = 2.6, cex = 0.7,
              col = Col$ink)
      }
      if (j == 4) mtext("year", side = 1, line = 1.6, cex = 0.7,
                        col = Col$ink2)
    }
  }
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Plot_E1_starts: the mean of each start on every board
#
#   rain 1.0: four starts, the settled band and the year-20 line
#   rain 0.8: five starts (also 1 hare per cell) and, in each panel,
#             the runs alive at year 50
#   Solid: mean of 20 runs (collapsed runs counted as 0); dashed: the
#   deterministic path from the same start.
#
#   rain   1.0 or 0.8
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_E1_starts = function(rain = 1.0) {
  wet = rain > 0.9
  starts = if (wet) E1_starts[1:4] else E1_starts
  y_top = if (wet) 6200 else 4500
  op = par(mfrow = c(2, 4), mar = c(3.2, 4.4, 2.2, 0.8))
  on.exit(par(op))
  par(cex = 0.85)
  yrs = 0:50
  for (b in 1:7) {
    Panel(c(0, 50), c(0, y_top), main = Board_label[b],
          xlab = if (b >= 4) "year" else "",
          ylab = if (b %in% c(1, 5)) "hares on the board, N" else "")
    if (wet) {
      bd = E1_band(b)
      rect(-5, bd[1] - 2 * bd[2], 55, bd[1] + 2 * bd[2],
           col = Col$band, border = NA)
      abline(v = 20, col = Col$muted, lwd = 0.9)
    }
    alive = c()
    for (k in starts) {
      r = E1_runs(b, rain, k)
      lines(yrs, colMeans(r$N), col = E1_start_col[k], lwd = 2)
      lines(yrs, E1_det(b, rain, k), col = E1_start_col[k], lwd = 1,
            lty = 2)
      alive = c(alive, sprintf("%-6s %2d/20", k, sum(is.na(r$t_ext))))
    }
    if (!wet) {
      text(50, y_top * 0.98, paste(c("alive at 50:", alive),
                                   collapse = "\n"),
           adj = c(1, 1), family = "mono", cex = 0.72, col = Col$ink2)
    }
  }
  labs = do.call(c, lapply(starts, E1_start_label))
  labs = c(labs, expression("deterministic path"))
  cols = c(E1_start_col[starts], Col$ink2)
  ltys = c(rep(1, length(starts)), 2)
  lwds = c(rep(2, length(starts)), 1.2)
  if (wet) {
    band = expression("settled band (mean " %+-% " 2 SD of one run)")
    labs = c(labs, band)
    cols = c(cols, Col$band)
    ltys = c(ltys, 1)
    lwds = c(lwds, 10)
  }
  n_runs = nrow(E1_runs(1, 1.0, "low")$N)    # 20 in the full design
  note = if (wet) {
    paste0(sprintf("Rain 1.0. Solid: mean of %d runs.\n", n_runs),
           "Grey line: year 20, where E2-E4 begin.")
  } else {
    sprintf(paste0("Rain 0.8. Solid: mean of %d runs,\n",
                   "collapsed runs counted as 0."), n_runs)
  }
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  Legend_panel(legend = labs, col = cols, lty = ltys, lwd = lwds,
               note = note)
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Plot_E1_drift_mean: E2's runs at rain 1.0 continue E1's boards for
# 200 more years: mean hares per board (A)
#
#   Five boards stay flat; the 10x10 boards keep draining.
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_E1_drift_mean = function() {
  d = Res("E2/runs")
  op = par(mar = c(3.6, 4.6, 2.2, 1))
  on.exit(par(op))
  Panel(c(20, 220), c(2700, 3500),
        main = "A. Rain 1.0 held for 200 years",
        xlab = "year (E1 clock: E2 starts at year 20)",
        ylab = "mean hares, N")
  for (b in 1:7) {
    s = d[d$board == b & abs(d$rain - 1.0) < 1e-9, ]
    N = Wide_years(s, "N_")[, as.character(0:200), drop = FALSE]
    lines(20 + 0:200, colMeans(N), col = Board_col[b], lwd = 1.6)
  }
  legend("bottomleft", Board_label, col = Board_col, lwd = 2, ncol = 2,
         bty = "n", cex = 0.8, seg.len = 2)
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Plot_E1_drift_quad: hares in the two poor quadrants of the 10x10
# boards, at year 20 (the start boards, solid) and year 220 (runs
# alive after 200 years at rain 1.0, hatched) (B)
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_E1_drift_quad = function() {
  sb = Res("E1/start_boards")
  d = Res("E2/runs")
  poor = (sb$row <= 10 & sb$col > 10) | (sb$row > 10 & sb$col <= 10)
  op = par(mar = c(3.6, 4.6, 2.2, 1))
  on.exit(par(op))
  Panel(c(-0.5, 1.5), c(0, 1300),
        main = "B. Poor quadrants: year 20 (solid) vs 220 (hatched)",
        ylab = "hares in the two poor quadrants", x_at = c(0, 1),
        x_lab = Board_label[6:7], grid_x = FALSE)
  for (b in 6:7) {
    x = b - 6
    s = sb[sb$board == b & poor, ]
    start = mean(tapply(s$n, s$run, sum))           # mean of 20 boards
    e = d[d$board == b & abs(d$rain - 1.0) < 1e-9 & is.na(d$t_ext), ]
    end = mean(e$quad200)
    rect(x - 0.35, 0, x - 0.01, start, col = Board_col[b],
         border = NA)
    rect(x + 0.01, 0, x + 0.35, end, col = "white",
         border = Board_col[b])
    rect(x + 0.01, 0, x + 0.35, end, col = Board_col[b], density = 14,
         angle = 45, border = Board_col[b])
  }
  return(invisible(NULL))
}

# ---- nicknames: one line each, no input (README.Rmd calls these) -----
E1_runs_plot = function() Plot_E1_runs()
E1_start_plot = function() Plot_E1_starts(1.0)
E1_drift_plot_1 = function() Plot_E1_drift_mean()
E1_drift_plot_2 = function() Plot_E1_drift_quad()
E1_dry_plot = function() Plot_E1_starts(0.8)

