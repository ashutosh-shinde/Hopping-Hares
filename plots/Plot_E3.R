# ======================================================================
# Plot_E3.R   figures of E3 (falling rain: does the speed matter?)
# ======================================================================
# Plotting code written with AI assistance (Claude).
#
#   Outcome page, section 5. Data: results/E3/ (run/Run_E3.R), the
#   deterministic ramps and folds (results/deterministic/) and E2.
#
#   Rain falls along a straight line from 1.0 to 0.4 over D years
#   (D = 50, 150, 450, 1,350 or 4,050), in monthly steps, then holds.
#   Collapse rain: the rain when N first falls below 10% of its start
#   (rho_10 = 1 - (0.6 / D) t_frac), or when the last hare dies.
#
#   Plot_E3_collapse()     collapse rain against ramp length
#   Plot_E3_threshold()    the 10% rule against the last hare
#   Plot_E3_N_vs_rain()    hares against the current rain
#   Plot_E3_empty_first()  rain at which each cell empties
#   nicknames at the end of the file
# ======================================================================

E3_durs = c(50, 150, 450, 1350, 4050)
E3_speed_col = c("#c9d8ee", "#93b5e0", "#5b8fd0", "#2a63b5", "#123d7a")

# ----------------------------------------------------------------------
# E3_rain_at: the rain on a ramp of D years at times t (held at 0.4
# after the ramp; NA stays NA)
# ----------------------------------------------------------------------
E3_rain_at = function(t, D) {
  return(pmax(0.4, 1.0 - 0.6 / D * t))
}

# ----------------------------------------------------------------------
# E3_collapse_rain: collapse rain of every run of board b
#
#   rule   "frac" (N below 10% of its start) or "last" (last hare)
#
#   returns  data frame: dur, run, rain (NA = never collapsed)
# ----------------------------------------------------------------------
E3_collapse_rain = function(b, rule = "frac") {
  d = Res("E3/runs")
  d = d[d$board == b, ]
  t = if (rule == "frac") d$t_frac else d$t_ext
  return(data.frame(dur = d$dur, run = d$run,
                    rain = E3_rain_at(t, d$dur)))
}

# fixed jitter on the log axis, so plotting draws no random numbers
E3_jitter = function(run) {
  return(exp(((run * 7) %% 11 / 10 - 0.5) * 0.18))
}

# ----------------------------------------------------------------------
# Plot_E3_collapse: collapse rain (10% rule) against ramp length, one
# panel per board, and all boards' means in the last panel
#
#   Dots: single runs; coloured line: their mean. Dashed black: the
#   deterministic board on the same ramp (lag alone, Det_ramp()).
#   Dotted black: the deterministic fold. Dash-dot: E2's rho_50(200).
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_E3_collapse = function() {
  ramps = Res("deterministic/ramps")
  op = par(mfrow = c(2, 4), mar = c(3.6, 4.4, 2.2, 0.8))
  on.exit(par(op))
  par(cex = 0.85)
  means = list()
  for (b in 1:7) {
    cr = E3_collapse_rain(b, "frac")
    Panel(c(40, 5000), c(0.53, 0.86), log = "x", x_at = E3_durs,
          main = Board_label[b],
          xlab = if (b >= 4) "ramp length, years (1.0 to 0.4)" else "",
          ylab = if (b %in% c(1, 5)) {
            "rain at collapse (N < 10% of start)"
          } else "")
    abline(h = Fold_final(b), col = Col$ink, lwd = 1, lty = 3)
    abline(h = Rho50(b, 200)$rho50, col = Board_col[b], lwd = 1,
           lty = 4)
    ok = !is.na(cr$rain)
    points(cr$dur[ok] * E3_jitter(cr$run[ok]), cr$rain[ok], pch = 16,
           cex = 0.55, col = adjustcolor(Board_col[b], 0.45))
    m = tapply(cr$rain, cr$dur, mean, na.rm = TRUE)[as.character(
      E3_durs)]
    means[[b]] = m
    lines(E3_durs, m, col = Board_col[b], lwd = 1.8)
    Board_points(E3_durs, m, b, cex = 1.15)
    r = ramps[ramps$board == b, ]
    lines(r$dur, r$rho10, col = Col$ink, lwd = 1.2, lty = 2)
  }
  Panel(c(40, 5000), c(0.53, 0.86), log = "x", x_at = E3_durs,
        main = "All boards (means)",
        xlab = "ramp length, years (1.0 to 0.4)")
  for (b in 1:7) {
    lines(E3_durs, means[[b]], col = Board_col[b], lwd = 1.4)
    Board_points(E3_durs, means[[b]], b, cex = 0.9)
  }
  legend("bottomright", Board_label, col = Board_col, pch = Board_pch,
         pt.bg = Board_col, pt.lwd = ifelse(Board_pch >= 21, 0.8, 2),
         lwd = 1.4, bty = "n", cex = 0.75, seg.len = 1.6)
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Plot_E3_threshold: the same 700 runs measured two ways
#
#   rule "frac": the rain when fewer than 10% of the starting hares
#   remain (A); rule "last": the rain when the last hare dies (B).
#   Triangles at the right edge: E2's rho_50(200) for each board.
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_E3_threshold = function(rule = "frac") {
  op = par(mar = c(3.6, 4.4, 2.2, 1.6))
  on.exit(par(op))
  main = if (rule == "frac") {
    "A \u00b7 Collapse = fewer than 10% of the starting hares"
  } else {
    "B \u00b7 Collapse = the last hare dies"
  }
  Panel(c(38, 6000), c(0.48, 0.86), log = "x", x_at = E3_durs,
        main = main, xaxs = "i",
        xlab = "ramp length, years (rain 1.0 to 0.4)",
        ylab = if (rule == "frac") "rain at collapse" else "")
  for (b in 1:7) {
    cr = E3_collapse_rain(b, rule)
    ok = !is.na(cr$rain)
    points(cr$dur[ok] * E3_jitter(cr$run[ok]), cr$rain[ok], pch = 16,
           cex = 0.45, col = adjustcolor(Board_col[b], 0.35))
    m = tapply(cr$rain, cr$dur, mean, na.rm = TRUE)[as.character(
      E3_durs)]
    lines(E3_durs, m, col = Board_col[b], lwd = 1.8)
    Board_points(E3_durs, m, b, cex = 1.15)
    y = Rho50(b, 200)$rho50                # an arrow at the right edge
    polygon(c(5600, 6400, 6400), c(y, y + 0.006, y - 0.006),
            col = Board_col[b], border = NA, xpd = NA)
  }
  if (rule == "last") {
    legend("bottomright", Board_label, col = Board_col, pch = Board_pch,
           pt.bg = Board_col, pt.lwd = ifelse(Board_pch >= 21, 0.8, 2),
           lwd = 1.8, bty = "n", cex = 0.85, seg.len = 2)
  }
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# E3_fixed_rain_N: E2's mean hares at fixed rain (years 100-200, runs
# alive at year 200; rains with at least 5 such runs)
#
#   returns  data frame: rain, N
# ----------------------------------------------------------------------
E3_fixed_rain_N = function(b) {
  d = Res("E2/runs")
  d = d[d$board == b & is.na(d$t_ext), ]
  N = Wide_years(d, "N_")[, as.character(100:200), drop = FALSE]
  out = data.frame(rain = numeric(0), N = numeric(0))
  for (v in sort(unique(d$rain))) {
    sel = abs(d$rain - v) < 1e-9
    if (sum(sel) >= 5) out = rbind(out, data.frame(rain = v,
                                                   N = mean(N[sel, ])))
  }
  return(out)
}

# ----------------------------------------------------------------------
# Plot_E3_N_vs_rain: mean hares (collapsed runs = 0) against the
# current rain for each ramp (light: fast, dark: slow). Black: the
# deterministic settled board; grey dots: E2's mean hares at fixed
# rain. Slow ramps follow E2; fast ramps stay above it (lag).
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_E3_N_vs_rain = function() {
  mn = Res("E3/mean_N")
  op = par(mfrow = c(2, 4), mar = c(3.6, 4.4, 2.2, 0.8))
  on.exit(par(op))
  par(cex = 0.85)
  for (b in 1:7) {
    Panel(c(1.0, 0.55), c(0, 4200), main = Board_label[b],
          xaxs = "i",
          xlab = if (b >= 4) "current rain (falling)" else "",
          ylab = if (b %in% c(1, 5)) "hares on the board, N" else "")
    Det_branch_lines(b, Col$ink, lwd = 1.2)
    e2 = E3_fixed_rain_N(b)
    points(e2$rain, e2$N, pch = 16, cex = 0.8, col = Col$muted)
    for (k in seq_along(E3_durs)) {
      s = mn[mn$board == b & mn$dur == E3_durs[k], ]
      lines(s$rain, s$N_mean, col = E3_speed_col[k], lwd = 1.5)
    }
  }
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  Legend_panel(legend = c(sprintf("ramp of %s years",
                                  format(E3_durs, big.mark = ",",
                                         trim = TRUE)),
                          "deterministic settled board",
                          "E2: mean N at fixed rain (years 100-200)"),
               col = c(E3_speed_col, Col$ink, Col$muted),
               lwd = c(rep(1.8, 5), 1.2, NA), pch = c(rep(NA, 6), 16),
               note = sprintf(paste0("Ramp lines: mean of %d runs,\n",
                                     "collapsed runs counted as 0."),
                              E3_n_runs()))
  return(invisible(NULL))
}

# E3_n_runs: runs per board and ramp in the saved tables (20 in the
# full design, 3 in a quick run)
E3_n_runs = function() {
  d = Res("E3/runs")
  return(sum(d$board == 1 & d$dur == d$dur[1]))
}

# ----------------------------------------------------------------------
# Plot_E3_empty_first: the rain at which each cell is empty in more
# than half the runs (boards saved every 0.005 of rain), on the
# fastest (50-year) and slowest (4,050-year) ramps. Brown: empties
# early; blue: holds out. Top: the grass.
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_E3_empty_first = function() {
  B = Read_boards()
  er = Res("E3/empty_rain")
  op = par(no.readonly = TRUE)
  on.exit(par(op))
  layout(matrix(c(1:7, 22, 8:14, 23, 15:21, 23), 3, 8, byrow = TRUE),
         widths = c(rep(1, 7), 0.5))
  par(oma = c(0, 4.6, 0, 0.2), cex = 0.9)
  row_lab = c("grass q", "50-year\nramp", "4,050-year\nramp")
  for (r in 1:3) {
    for (b in 1:7) {
      par(mar = c(0.4, 0.4, if (r == 1) 2.0 else 0.6, 0.4))
      if (r == 1) {
        Map_panel(B[[b]], c(0.3, 0.96), Pal$q, main = Board_label[b],
                  frame = TRUE)
      } else {
        dur = if (r == 2) 50 else 4050
        s = er[er$board == b & er$dur == dur, ]
        g = matrix(NA_real_, 20, 20)
        g[cbind(s$row, s$col)] = s$rain_empty
        Map_panel(g, c(0.6, 1.0), Pal$empty, frame = TRUE)
      }
      if (b == 1) {
        mtext(row_lab[r], side = 2, line = 0.8, las = 1, adj = 1,
              cex = 0.9, col = Col$ink)
      }
    }
  }
  plot.new()
  par(mar = c(1, 0.2, 1, 3.2))
  Colour_bar(Pal$empty, c(0.6, 1.0),
             "rain at which the cell is empty in most runs",
             at = seq(0.6, 1.0, 0.05))
  return(invisible(NULL))
}

# ---- nicknames: one line each, no input (README.Rmd calls these) -----
E3_collapse_plot = function() Plot_E3_collapse()
E3_threshold_plot_1 = function() Plot_E3_threshold("frac")
E3_threshold_plot_2 = function() Plot_E3_threshold("last")
E3_N_vs_rain_plot = function() Plot_E3_N_vs_rain()
E3_empty_first_plot = function() Plot_E3_empty_first()

