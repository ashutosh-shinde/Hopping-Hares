# ======================================================================
# Plot_E4.R   figures of E4 (coming back: how big must a colony be?)
# ======================================================================
# Plotting code written with AI assistance (Claude).
#
#   Outcome page, section 6. Data: results/E4/ (run/Run_E4.R).
#
#   The Real board starts empty apart from the released colony: k x k
#   cells with 14 hares each, where the mean grass over the square is
#   highest. Arm A: k = 1, 2, 3, 4, 5, 6, 8. Arm B: 224 hares as one
#   4x4 colony, four 2x2 colonies or sixteen single cells. Rains 1.0,
#   0.9, 0.85 and 0.785; 50 runs of 50 years each. Success: hares
#   alive at year 50.
#
#   Plot_E4_sites()     where the colonies go
#   Plot_E4_size()      share of runs alive against colony size
#   Plot_E4_designs()   one big colony or several small ones
#   Plot_E4_growth()    hares after release (surviving runs)
#   nicknames at the end of the file
# ======================================================================

E4_rains = c(1.0, 0.9, 0.85, 0.785)
E4_rain_col = c("#2a78d6", "#1baf7a", "#eda100", "#eb6834")
E4_ks = c(1, 2, 3, 4, 5, 6, 8)
E4_designs = c("1 x 4x4", "4 x 2x2", "16 x 1x1")

# ----------------------------------------------------------------------
# E4_runs: the 50 runs of one design at one rain
#
#   arm      "A" or "B"
#   design   "2x2" (arm A) or "4 x 2x2" (arm B), ...
#   rain     one of E4_rains
#
#   returns  data frame (one row per run)
# ----------------------------------------------------------------------
E4_runs = function(arm, design, rain) {
  d = Res("E4/runs")
  return(d[d$arm == arm & d$design == design &
             abs(d$rain - rain) < 1e-9, ])
}

E4_label = function(design) {
  return(gsub("x", "\u00d7", design))
}

E4_rain_text = function(rain) {
  return(sub("[.]?0+$", "", sprintf("%.3f", rain)))
}

# ----------------------------------------------------------------------
# E4_block: draw a k x k colony with its top-left cell at row r0,
# column c0 on a Map_panel of a 20 x 20 grid
#
#   ...   rect() settings: border, col, lwd
# ----------------------------------------------------------------------
E4_block = function(r0, c0, k, ...) {
  rect(c0 - 0.5, 20.5 - (r0 + k - 1), c0 + k - 0.5, 21.5 - r0, ...)
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Plot_E4_sites: where the colonies go on the Real board
#
#   A: the best k x k square for each size (all in the grassy north).
#   B: the three designs of 224 hares (the best spot inside each
#   quadrant, or inside each 5 x 5 part).
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_E4_sites = function() {
  q = Read_boards()[[1]]
  s = Res("E4/sites")
  kc = c("1" = "#eb6834", "2" = "#eda100", "3" = "#1baf7a",
         "4" = "#2a78d6", "5" = "#6250d6", "6" = "#e87ba4",
         "8" = "#0b0b0b")
  op = par(no.readonly = TRUE)
  on.exit(par(op))
  layout(matrix(1:5, 1, 5), widths = c(0.42, 1, 1, 1, 1))
  par(cex = 0.9, mar = c(0.4, 0.2, 2.2, 0.4))
  plot.new()
  legend("topleft", sprintf("%d\u00d7%d", E4_ks, E4_ks), fill = NA,
         border = kc, bty = "n", cex = 0.85, inset = c(0, 0.04))
  Map_panel(q, c(0.3, 0.96), Pal$q, main = "A. One colony, side 1-8")
  for (k in E4_ks) {
    a = s[s$arm == "A" & s$design == sprintf("%dx%d", k, k), ]
    E4_block(a$row, a$col, k, border = kc[as.character(k)], lwd = 2)
  }
  for (dsn in E4_designs) {
    Map_panel(q, c(0.3, 0.96), Pal$q,
              main = paste("B. 224 hares as", E4_label(dsn)))
    a = s[s$arm == "B" & s$design == dsn, ]
    for (j in seq_len(nrow(a))) {
      E4_block(a$row[j], a$col[j], a$k[j],
               col = adjustcolor("#eb6834", 0.85), border = "#8a3b0c")
    }
  }
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Plot_E4_size: share of runs alive at year 50 against colony size,
# one panel per rain; shaded: 95% Wilson interval over 50 runs
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_E4_size = function() {
  op = par(mfrow = c(1, 4), mar = c(3.6, 1.2, 2.2, 0.8),
           oma = c(1.4, 3.4, 0, 0))
  on.exit(par(op))
  par(cex = 0.9)
  for (j in seq_along(E4_rains)) {
    v = E4_rains[j]
    alive = sapply(E4_ks, function(k) {
      return(sum(is.na(E4_runs("A", sprintf("%dx%d", k, k), v)$t_ext)))
    })
    n = sapply(E4_ks, function(k) {
      return(nrow(E4_runs("A", sprintf("%dx%d", k, k), v)))
    })
    ci = sapply(seq_along(E4_ks), function(i) Wilson(alive[i], n[i]))
    main = paste("rain", E4_rain_text(v))
    if (abs(v - 0.785) < 1e-9) {
      main = bquote(bold("rain 0.785  (E2's " * rho[50] * " for Real)"))
    }
    Panel(c(0.5, 8.5), c(-0.03, 1.05), main = main, x_at = E4_ks,
          x_lab = sprintf("%d\u00d7%d\n%d", E4_ks, E4_ks,
                          14 * E4_ks^2),
          y_lab = j == 1, ylab = "")
    polygon(c(E4_ks, rev(E4_ks)), c(ci[1, ], rev(ci[2, ])),
            col = adjustcolor(E4_rain_col[j], 0.15), border = NA)
    lines(E4_ks, alive / n, col = E4_rain_col[j], lwd = 2)
    points(E4_ks, alive / n, pch = 21, bg = E4_rain_col[j],
           col = "white", cex = 1.2)
  }
  mtext("share of runs alive at year 50", side = 2, outer = TRUE,
        line = 2.2, col = Col$ink2, cex = 0.9)
  mtext(paste0("colony side and hares released (14 per cell). ",
               "Shaded band: 95% interval over 50 runs."),
        side = 1, outer = TRUE, line = 0.4, col = Col$ink2, cex = 0.85)
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Plot_E4_designs: 224 hares as one 4x4 colony, four 2x2 colonies or
# sixteen single cells: share of runs alive at year 50
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_E4_designs = function() {
  dcol = c("#2a78d6", "#eda100", "#eb6834")
  op = par(mar = c(3.0, 4.4, 2.2, 1))
  on.exit(par(op))
  Panel(c(0.5, 4.5), c(0, 1.05), x_at = 1:4,
        x_lab = paste("rain", sapply(E4_rains, E4_rain_text)),
        main = "224 hares: one big colony or several small ones",
        ylab = "share of runs alive at year 50", grid_x = FALSE,
        yaxs = "i")
  w = 0.26
  for (j in seq_along(E4_designs)) {
    share = sapply(E4_rains, function(v) {
      r = E4_runs("B", E4_designs[j], v)
      return(mean(is.na(r$t_ext)))
    })
    x = 1:4 + (j - 2) * w
    rect(x - w * 0.46, 0, x + w * 0.46, share, col = dcol[j],
         border = NA)
  }
  legend("topright", E4_label(E4_designs), fill = dcol, border = NA,
         bty = "n", cex = 0.9)
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Plot_E4_growth: mean hares of the runs alive at year 50, after the
# release, for several colony sizes (sizes with at least 3 such runs)
#
#   rain   1.0 or 0.85
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_E4_growth = function(rain = 1.0) {
  kcol = c("2" = "#eda100", "3" = "#1baf7a", "4" = "#2a78d6",
           "6" = "#6250d6", "8" = "#0b0b0b")
  op = par(mar = c(3.6, 4.6, 2.2, 1))
  on.exit(par(op))
  Panel(c(0, 50), c(0, 3400), main = paste("rain", E4_rain_text(rain)),
        xlab = "years after release",
        ylab = if (rain > 0.9) "hares on the board (surviving runs)"
        else "")
  labs = c()
  cols = c()
  for (k in c(2, 3, 4, 6, 8)) {
    r = E4_runs("A", sprintf("%dx%d", k, k), rain)
    alive = is.na(r$t_ext)
    if (sum(alive) < 3) next
    N = Wide_years(r[alive, ], "N_")
    lines(0:50, colMeans(N), col = kcol[as.character(k)], lwd = 2)
    labs = c(labs, sprintf("%d\u00d7%d (%d/%d alive)", k, k,
                           sum(alive), nrow(r)))
    cols = c(cols, kcol[as.character(k)])
  }
  legend("bottomright", labs, col = cols, lwd = 2, bty = "n",
         cex = 0.85, seg.len = 2)
  return(invisible(NULL))
}

# ---- nicknames: one line each, no input (README.Rmd calls these) -----
E4_sites_plot = function() Plot_E4_sites()
E4_size_plot = function() Plot_E4_size()
E4_designs_plot = function() Plot_E4_designs()
E4_growth_plot_1 = function() Plot_E4_growth(1.0)
E4_growth_plot_2 = function() Plot_E4_growth(0.85)
