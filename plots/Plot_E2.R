# ======================================================================
# Plot_E2.R   figures of E2 (fixed rain: where does each board tip?)
# ======================================================================
# Plotting code written with AI assistance (Claude).
#
#   Outcome page, section 4. Data: results/E2/ (run/Run_E2.R) and
#   the deterministic folds (results/deterministic/).
#
#   Tipping rain rho_50(T): the rain at which half the runs lose every
#   hare within T years (Fit_logistic() in Measure.R).
#
#   Plot_E2_tipping()        tipping curves at 200 years (A)
#   Plot_E2_rho50()          rho_50(200) by board, bootstrap bars (B)
#   Plot_E2_fold_gap()       deterministic fold against rho_50
#   Plot_E2_hares_vs_det()   hares at year 200 over the settled board
#   Plot_E2_mean_N()         mean hares over time at six rains
#   Plot_E2_story()          every run at six rains per board
#   Plot_E2_pullback()       poor cells against good cells over time
#   Plot_E2_refuges()        where the last hares live
#   Plot_E2_survival()       Kaplan-Meier curves near rho_50
#   Plot_E2_window()         rho_50 against the waiting window
#   Plot_E2_long()           tipping curves up to 2,000 years
#   Plot_E2_theory_*()       exact one-cell curves (windows 10-10^8)
#   nicknames at the end of the file
# ======================================================================

# ----------------------------------------------------------------------
# E2 data helpers
#
#   E2_runs(b)        E2's runs on board b (one row per run)
#   E2_dead(d, T)     list(dead, usable): which runs lost every hare
#                     by year T. Past year 200 only runs at rains up
#                     to top300 were followed to year 300, so only
#                     those are usable.
#   E2_traj(d, rain)  runs x years matrix of hares at one rain; years
#                     201-300 where runs were followed (dead = 0)
#   Rho50(b, T)       saved rho_50, width and bootstrap bounds
#   Fold_final(b)     the board's last deterministic fold
# ----------------------------------------------------------------------
E2_runs = function(b) {
  d = Res("E2/runs")
  return(d[d$board == b, ])
}

E2_dead = function(d, T) {
  if (T <= 200) {
    dead = !is.na(d$t_ext) & d$t_ext <= T
    usable = rep(TRUE, nrow(d))
  } else {
    usable = d$rain <= d$top300 + 1e-9
    dead = !is.na(d$t_ext) | (!is.na(d$t_ext300) & d$t_ext300 <= T)
  }
  return(list(dead = dead, usable = usable))
}

E2_traj = function(d, rain) {
  s = d[abs(d$rain - rain) < 1e-9, ]
  N = Wide_years(s, "N_")
  if (rain <= s$top300[1] + 1e-9) {
    N[is.na(N)] = 0                      # dead by 200: 0 afterwards
  } else {
    N = N[, as.character(0:200), drop = FALSE]
  }
  return(N)
}

Rho50 = function(b, T = 200) {
  r = Res("E2/rho50")
  return(r[r$board == b & r$window == T, ])
}

Fold_final = function(b) {
  f = Res("deterministic/folds")
  f = f[f$board == b, ]
  return(f$rho[nrow(f)])
}

# the plan's two rough guesses: f_c divided by the mean grass of all
# cells (well-mixed) and of the better half (refuge)
E2_lines = function() {
  f_c = F_crit(Set_params())
  return(c(mixed = f_c / 0.662, refuge = f_c / 0.766))
}

# ----------------------------------------------------------------------
# Plot_E2_tipping: share of runs collapsed by year 200 against rain,
# with logistic fits (A)
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_E2_tipping = function() {
  op = par(mar = c(3.6, 4.4, 3.6, 1))
  on.exit(par(op))
  Panel(c(0.45, 1.05), c(-0.03, 1.03),
        xlab = expression("rain " * rho *
                            " (held fixed for 200 years)"),
        ylab = "share of runs collapsed by year 200")
  Panel_title("A. Tipping curves", line = 2.4)
  ln = E2_lines()
  abline(v = ln, col = Col$ink2, lwd = 1, lty = 2)
  text(ln, 1.06, c("well-mixed\nline 0.75", "refuge\nline 0.65"),
       adj = c(0.5, 0), cex = 0.8, col = Col$ink2, xpd = NA)
  xx = seq(0.3, 1.25, length.out = 400)
  for (b in 1:7) {
    d = E2_runs(b)
    r = Rho50(b, 200)
    lines(xx, 1 / (1 + exp((xx - r$rho50) / r$width)),
          col = Board_col[b], lwd = 2)
    xs = sort(unique(d$rain))
    sh = sapply(xs, function(v) mean(!is.na(d$t_ext[d$rain == v])))
    Board_points(xs, sh, b, cex = 0.75)
  }
  legend("topright", Board_label, col = Board_col, pch = Board_pch,
         pt.bg = Board_col, pt.lwd = ifelse(Board_pch >= 21, 0.8, 2),
         bty = "n", cex = 0.85, inset = c(0.01, 0.02))
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Plot_E2_rho50: rho_50(200) by board with 95% bootstrap intervals
# (narrower than the markers), sorted (B)
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_E2_rho50 = function() {
  r = Res("E2/rho50")
  r = r[r$window == 200, ]
  r = r[order(r$rho50), ]
  op = par(mar = c(3.6, 8.4, 3.6, 1))
  on.exit(par(op))
  Panel(c(0.64, 0.87), c(-0.5, 7.1), xlab = expression(
    rho[50] * " (200 years), 95% bootstrap interval"),
    y_at = 0:6, y_lab = Board_label[r$board], grid_y = FALSE)
  Panel_title(expression(bold("B. Tipping rain " * rho[50] *
                                " by board")), line = 2.4)
  ln = E2_lines()
  abline(v = ln, col = Col$ink2, lwd = 1, lty = 2)
  text(ln, 6.75, c("well-mixed 0.75", "refuge 0.65"), adj = c(0.5, 0),
       cex = 0.8, col = Col$ink2, xpd = NA)
  for (k in seq_len(nrow(r))) {
    b = r$board[k]
    segments(r$lo[k], k - 1, r$hi[k], k - 1, col = Board_col[b],
             lwd = 3, lend = 0)
    Board_points(r$rho50[k], k - 1, b, cex = 1.3)
    text(r$hi[k] + 0.003, k - 1, sprintf("%.3f", r$rho50[k]),
         adj = c(0, 0.5), cex = 0.85, col = Col$ink2)
  }
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Plot_E2_fold_gap: the deterministic final fold (filled) against
# rho_50(200) (hollow). The number is the gap, the part of tipping
# that comes from luck; ticks mark the first fold of the clumped
# boards.
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_E2_fold_gap = function() {
  folds = Res("deterministic/folds")
  r = Res("E2/rho50")
  r = r[r$window == 200, ]
  r = r[order(r$rho50), ]
  op = par(mar = c(3.6, 8.4, 2.4, 1))
  on.exit(par(op))
  Panel(c(0.69, 0.93), c(-0.6, 7.0), y_at = 0:6,
        y_lab = Board_label[r$board], grid_y = FALSE,
        xlab = expression("rain " * rho))
  Panel_title(expression(bold("Final fold (filled) vs E2 tipping " *
                                "rain " * rho[50] *
                                " at 200 years (hollow)")))
  mixed = E2_lines()[["mixed"]]
  abline(v = mixed, col = Col$ink2, lwd = 1, lty = 2)
  text(mixed, 6.75, "f_c / mean grass", cex = 0.8, col = Col$ink2)
  for (k in seq_len(nrow(r))) {
    b = r$board[k]
    f = folds[folds$board == b, ]
    fin = f$rho[nrow(f)]
    segments(fin, k - 1, r$rho50[k], k - 1,
             col = adjustcolor(Board_col[b], 0.5), lwd = 4)
    Board_points(fin, k - 1, b, cex = 1.5)
    Board_points(r$rho50[k], k - 1, b, cex = 1.5, hollow = TRUE)
    if (nrow(f) > 1) {
      segments(f$rho[1], k - 1.22, f$rho[1], k - 0.78,
               col = Board_col[b], lwd = 2)
    }
    text(r$rho50[k] + 0.006, k - 1, sprintf("+%.3f", r$rho50[k] - fin),
         adj = c(0, 0.5), cex = 0.85, col = Col$ink2)
  }
  text(0.93, -0.5, "ticks: first fold (first region lost)",
       adj = c(1, 0.5), cex = 0.8, col = Col$ink2)
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Plot_E2_hares_vs_det: hares at year 200 in every run against rain,
# over the deterministic settled board
#
#   Dots: runs alive at year 200; grey crosses: runs with no hare
#   left. Black: the settled board, dotted at each fold's drop.
#   Dashed: rho_50(200). Runs sit below the black line at every rain
#   and reach 0 at rho_50, where the settled board still holds
#   1,200-3,000 hares.
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_E2_hares_vs_det = function() {
  op = par(mfrow = c(2, 4), mar = c(3.4, 4.4, 2.2, 0.8))
  on.exit(par(op))
  par(cex = 0.85)
  for (b in 1:7) {
    d = E2_runs(b)
    Panel(c(0.6, 1.2), c(-150, 5000), main = Board_label[b],
          xlab = if (b >= 4) expression("rain " * rho) else "",
          ylab = if (b %in% c(1, 5)) "hares at year 200" else "",
          x_at = c(0.6, 0.8, 1.0, 1.2))
    dead = !is.na(d$t_ext)
    jit = ((d$run * 7) %% 11 / 10 - 0.5) * 0.004  # fixed small jitter
    points(d$rain[!dead] + jit[!dead], d$N_200[!dead], pch = 16,
           cex = 0.45, col = adjustcolor(Board_col[b], 0.55))
    points(d$rain[dead] + jit[dead], rep(0, sum(dead)), pch = 4,
           cex = 0.5, col = Col$muted)
    Det_branch_lines(b, Col$ink, lwd = 1.5)
    r50 = Rho50(b, 200)$rho50
    fin = Fold_final(b)
    abline(v = r50, col = Board_col[b], lwd = 1.4, lty = 2)
    text(r50 + 0.004, 4800, bquote(rho[50] ~ .(sprintf("%.3f", r50))),
         adj = c(0, 0.5), cex = 0.85, col = Board_col[b])
    text(fin - 0.004, 4800, sprintf("fold %.3f", fin), adj = c(1, 0.5),
         cex = 0.85, col = Col$ink)
  }
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  Legend_panel(legend = c("deterministic settled board",
                          "one run, alive at year 200",
                          "one run, no hare left",
                          expression("E2 tipping rain " * rho[50] *
                                       "(200)")),
               col = c(Col$ink, Board_col[2], Col$muted, Col$ink2),
               lty = c(1, NA, NA, 2), lwd = c(1.5, NA, NA, 1.4),
               pch = c(NA, 16, 4, NA),
               note = paste0("Dotted lines: the jump at each\n",
                             "deterministic fold."))
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Plot_E2_mean_N: mean hares over time for every board at six rains
# (collapsed runs counted as 0; dotted: year 200)
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_E2_mean_N = function(rains = c(0.75, 0.775, 0.8, 0.825, 0.85,
                                     0.9)) {
  op = par(mfrow = c(2, 3), mar = c(3.4, 4.6, 2.2, 0.8))
  on.exit(par(op))
  par(cex = 0.9)
  for (k in seq_along(rains)) {
    Panel(c(0, 300), c(0, 3000), main = sprintf("rain %.3f", rains[k]),
          xlab = if (k > 3) "year" else "",
          ylab = if (k %in% c(1, 4)) "mean hares N (collapsed runs = 0)"
          else "", xaxs = "i", yaxs = "i")
    abline(v = 200, col = Col$muted, lwd = 0.8, lty = 3)
    for (b in 1:7) {
      d = E2_runs(b)
      if (!any(abs(d$rain - rains[k]) < 1e-9)) next
      N = E2_traj(d, rains[k])
      lines(0:(ncol(N) - 1), colMeans(N), col = Board_col[b], lwd = 1.8)
    }
    if (k == 1) {
      legend("topright", Board_label, col = Board_col, lwd = 2,
             bty = "n", cex = 0.8, seg.len = 2)
    }
  }
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# E2_story_label: what the runs at one rain do (the panel titles of
# Plot_E2_story), from the rules of the outcome page
#
#   below the fold                      "below fold: all slide to 0"
#   every run dies, by the poor cells'  "all die: poor half empties
#   share at year 50 (pullback runs)     first" (< 20%), "... poor
#                                        cells thin faster" (< 30%),
#                                        "... thin out evenly"
#   some die                            "k/n die: luck decides"
#   rain 1.0                            "full and steady"
#   poor cells hold < 15% at year 200   "holds in refuges, poor half
#                                        empty"
#   otherwise                           "persists, thinned to x%"
#
#   b, rain   board and rain
#   N         runs x years from E2_traj()
#   s         the runs' rows of E2_runs(b) at this rain
#
#   returns  one line of text
# ----------------------------------------------------------------------
E2_story_label = function(b, rain, N, s) {
  dead = sum(N[, ncol(N)] == 0)
  n = nrow(N)
  if (rain < Fold_final(b)) return("below fold: all slide to 0")
  if (dead == n) {
    pb = Res("E2/pullback")
    pb = pb[pb$board == b & pb$edge == "wet", ]
    ps = sum(pb$Y_50[pb$group == "poor"]) / sum(pb$Y_50)
    if (ps < 0.2) return("all die: poor half empties first")
    if (ps < 0.3) return("all die: poor cells thin faster")
    return("all die: thin out evenly")
  }
  if (dead > 0) return(sprintf("%d/%d die: luck decides", dead, n))
  if (abs(rain - 1.0) < 1e-9) return("full and steady")
  alive = s$N_200 > 0
  ps = mean(s$poor200[alive] / s$N_200[alive])
  if (ps < 0.15) return("holds in refuges, poor half empty")
  keep = mean(N[, "200"]) / mean(N[, "0"])
  return(sprintf("persists, thinned to %d%%", round(100 * keep)))
}

# ----------------------------------------------------------------------
# Plot_E2_story: every run (thin) and the mean (black, collapsed
# runs = 0) at six rains per board: 0.10, 0.025 and 0.005 below its
# rho_50(200), 0.005 and 0.025 above it, and rain 1.0. Each title
# says what the hares do; dotted: year 200.
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_E2_story = function() {
  offs = c(-0.10, -0.025, -0.005, 0.005, 0.025, NA)  # NA = rain 1.0
  op = par(no.readonly = TRUE)
  on.exit(par(op))
  layout(rbind(rep(1, 6), matrix(2:43, 7, 6, byrow = TRUE)),
         heights = c(0.3, rep(1, 7)))
  par(cex = 0.75, oma = c(2.4, 3.8, 0, 0))
  par(mar = c(0, 0, 0, 0))
  Legend_panel(legend = c("one run",
                          "mean of all runs (collapsed = 0)"),
               col = c(Col$muted, Col$ink), lwd = c(0.8, 2),
               where = "center", horiz = TRUE)
  for (b in 1:7) {
    d = E2_runs(b)
    xs = sort(unique(d$rain))
    r50 = Rho50(b, 200)$rho50
    for (j in seq_along(offs)) {
      v = 1.0                            # NA offset: rain 1.0
      if (!is.na(offs[j])) v = xs[which.min(abs(xs - (r50 + offs[j])))]
      N = E2_traj(d, v)
      s = d[abs(d$rain - v) < 1e-9, ]
      par(mar = c(1.6, 1.6, 2.6, 0.6))
      Panel(c(0, 300), c(0, 4400), x_lab = b == 7, y_lab = j == 1,
            xaxs = "i", yaxs = "i")
      if (j == 1) {
        mtext(Board_label[b], side = 2, line = 3.6, font = 2,
              cex = 0.8, col = Col$ink)
        mtext("hares N", side = 2, line = 2.6, cex = 0.75,
              col = Col$ink2)
      }
      if (b == 7) mtext("year", side = 1, line = 2.0, cex = 0.8,
                        col = Col$ink2)
      lab = E2_story_label(b, v, N, s)
      title(main = lab, adj = 0, line = 1.25, cex.main = 0.95,
            font.main = 2, col.main = Col$ink)
      title(main = bquote(rho == .(sprintf("%.3f", v))), adj = 0,
            line = 0.3, cex.main = 0.95, font.main = 1,
            col.main = Col$ink)
      abline(v = 200, col = Col$muted, lwd = 0.8, lty = 3)
      yrs = 0:(ncol(N) - 1)
      for (i in seq_len(nrow(N))) {
        lines(yrs, N[i, ], col = adjustcolor(Board_col[b], 0.35),
              lwd = 0.6)
      }
      lines(yrs, colMeans(N), col = Col$ink, lwd = 1.8)
    }
  }
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Plot_E2_pullback: hares in the 200 poorest cells against hares in
# the 200 best, year 0 to 200, at the wet edge rho_w(200) (blue: 10%
# of runs die within 200 years) and the dry edge rho_d(200) (orange:
# 90% die)
#
#   Bold: mean path; circle: mean at year 0; square: mean at year 200
#   (dead runs = 0). Dotted: poor = good. Dashed: even thinning (the
#   start ratio kept). Fine boards slide down the dashed line; clumped
#   boards empty the poor half first.
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_E2_pullback = function() {
  pb = Res("E2/pullback")
  cw = "#2a78d6"
  cd = "#eb6834"
  op = par(mfrow = c(2, 4), mar = c(3.6, 4.4, 2.2, 0.8))
  on.exit(par(op))
  par(cex = 0.85)
  for (b in 1:7) {
    s = pb[pb$board == b, ]
    hi = 1.05 * max(as.matrix(s[, grep("^Y_", names(s))]))
    Panel(c(0, hi), c(0, hi), main = Board_label[b],
          xlab = "hares in good cells (top 200 q)",
          ylab = "hares in poor cells (bottom 200 q)")
    lines(c(0, hi), c(0, hi), lty = 3, col = Col$ink2, lwd = 1.2)
    labs = c()
    for (edge in c("wet", "dry")) {
      e = s[s$edge == edge, ]
      poor = Wide_years(e[e$group == "poor", ], "Y_")
      good = Wide_years(e[e$group == "good", ], "Y_")
      col = if (edge == "wet") cw else cd
      for (i in seq_len(nrow(poor))) {
        lines(good[i, ], poor[i, ], col = adjustcolor(col, 0.28),
              lwd = 0.5)
      }
      mp = colMeans(poor)
      mg = colMeans(good)
      lines(mg, mp, col = col, lwd = 2.6)
      if (edge == "wet") {
        lines(c(0, mg[1]), c(0, mp[1]), lty = 2, col = Col$muted)
      }
      points(mg[1], mp[1], pch = 21, bg = "white", col = Col$ink,
             cex = 1.5, lwd = 1.6)
      points(mg[length(mg)], mp[length(mp)], pch = 22, bg = col,
             col = Col$ink, cex = 1.4, lwd = 1.2)
      dead = sum(!is.na(e$t_ext[e$group == "poor"]))
      labs = c(labs, as.expression(bquote(
        rho[.(if (edge == "wet") "w" else "d")] ==
          .(sprintf("%.3f (%d/20 die)", e$rain[1], dead)))))
    }
    legend("topleft", legend = labs, col = c(cw, cd), lwd = 2.6,
           bty = "n", cex = 0.85, seg.len = 1.6)
  }
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  Legend_panel(legend = c(expression("wet edge " * rho[w] * "(200): " *
                                       "10% of runs die in 200 yr"),
                          expression("dry edge " * rho[d] * "(200): " *
                                       "90% die"),
                          "one run (year 0-200)", "mean at year 0",
                          "mean at year 200 (dead runs = 0)",
                          "diagonal: poor = good",
                          "even thinning: start ratio kept"),
               col = c(cw, cd, Col$muted, Col$ink, Col$ink, Col$ink2,
                       Col$muted),
               lwd = c(2.6, 2.6, 0.8, NA, NA, 1.2, 1),
               lty = c(1, 1, 1, NA, NA, 3, 2),
               pch = c(NA, NA, NA, 21, 22, NA, NA),
               pt.bg = c(NA, NA, NA, "white", Col$muted, NA, NA),
               pt.cex = 1.4)
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Plot_E2_refuges: where the last hares live. Top: grass. Bottom:
# share of years 180-200 each cell held hares, averaged over the runs
# alive at year 200 with rain from rho_50(200) to rho_50 + 0.05.
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_E2_refuges = function() {
  B = Read_boards()
  occ = Res("E2/occupancy")
  op = par(no.readonly = TRUE)
  on.exit(par(op))
  layout(matrix(c(1:7, 15, 8:14, 16), 2, 8, byrow = TRUE),
         widths = c(rep(1, 7), 0.45))
  par(oma = c(0, 3.2, 0, 0.4), cex = 0.9)
  for (b in 1:7) {
    par(mar = c(0.3, 0.4, 2.2, 0.4))
    Map_panel(B[[b]], c(0.3, 0.96), Pal$q, main = Board_label[b])
    if (b == 1) mtext("grass q", side = 2, line = 0.6, cex = 0.9)
  }
  for (b in 1:7) {
    r50 = Rho50(b, 200)$rho50
    s = occ[occ$board == b & occ$rain >= r50 - 1e-9 &
              occ$rain <= r50 + 0.05 + 1e-9, ]
    g = matrix(NA_real_, 20, 20)
    lab = "no survivors"
    if (nrow(s) > 0) {
      w = matrix(0, 20, 20)              # share x runs, summed
      for (k in seq_len(nrow(s))) {
        w[s$row[k], s$col[k]] = w[s$row[k], s$col[k]] +
          s$occ[k] * s$n_alive[k]
      }
      n_runs = sum(tapply(s$n_alive, s$rain, function(x) x[1]))
      g = w / n_runs
      lab = bquote(.(n_runs) ~ "runs," ~ rho ~
                     .(sprintf("%.3f-%.3f", min(s$rain), max(s$rain))))
    }
    par(mar = c(0.3, 0.4, 2.2, 0.4))
    Map_panel(g, c(0, 1), Pal$occ)
    title(main = lab, adj = 0, line = 0.4, cex.main = 0.85,
          font.main = 1, col.main = Col$ink2)
    if (b == 1) {
      mtext("share of years\n180-200 occupied", side = 2, line = 0.6,
            cex = 0.9)
    }
  }
  par(mar = c(1.5, 0.2, 2.5, 2.6))
  Colour_bar(Pal$q, c(0.3, 0.96), "q", at = seq(0.3, 0.9, 0.1))
  Colour_bar(Pal$occ, c(0, 1), "occupied share")
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Survival colours: orange for rains below the board's rho_50, blue
# above, darker further away (v = rain - rho_50, from -0.1 to 0.1)
# ----------------------------------------------------------------------
E2_surv_col = function(v) {
  v = pmin(pmax(v, -0.1), 0.1)
  dry = colorRamp(c("#a33a00", "#f08a4b"))
  wet = colorRamp(c("#3a7bd5", "#08306b"))
  rgbm = if (v < 0) dry((v + 0.1) / 0.1) else wet(v / 0.1)
  return(rgb(rgbm[1], rgbm[2], rgbm[3], maxColorValue = 255))
}

# ----------------------------------------------------------------------
# Plot_E2_survival: share of runs alive against year, one curve per
# rain within 0.1 of rho_50(200) (Kaplan-Meier; years 200-300 only
# where runs were followed)
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_E2_survival = function() {
  op = par(mfrow = c(2, 4), mar = c(3.4, 4.4, 2.2, 0.8))
  on.exit(par(op))
  par(cex = 0.85)
  y_at = c(0.1, 1)
  for (b in 1:7) {
    d = E2_runs(b)
    r50 = Rho50(b, 200)$rho50
    Panel(c(0, 300), c(0.03, 1.1), log = "y", main = Board_label[b],
          xlab = if (b >= 4) "year" else "",
          ylab = if (b %in% c(1, 5)) "share alive (log)" else "",
          y_at = y_at, y_lab = Log_labels(y_at), xaxs = "i")
    abline(v = 200, col = Col$muted, lwd = 0.8, lty = 3)
    xs = sort(unique(d$rain))
    xs = xs[abs(xs - r50) <= 0.1 + 1e-9]
    for (v in xs) {
      s = d[abs(d$rain - v) < 1e-9, ]
      te = ifelse(s$cont == 1 & !is.na(s$t_ext300), s$t_ext300, s$t_ext)
      cens = if (any(s$cont == 1) || v <= s$top300[1]) 300 else 200
      km = Km(te, cens)
      tt = c(km$t, cens)
      ss = pmax(c(km$S, km$S[nrow(km)]), 1e-3)
      lines(tt, ss, type = "s", col = E2_surv_col(v - r50), lwd = 1.3)
    }
  }
  par(mar = c(3.6, 1.5, 1.2, 1.5))
  plot.new()
  plot.window(c(0, 1), c(0, 1))
  text(0.5, c(0.86, 0.79, 0.72),
       c("S(t): share of runs", "still alive", "(years 200-300 only"),
       col = Col$ink2, cex = 1)
  text(0.5, 0.65, expression("near " * rho[50] * ")"), col = Col$ink2,
       cex = 1)
  vv = seq(-0.1, 0.1, length.out = 201)
  cols = sapply(vv[-201] + 0.0005, E2_surv_col)
  rect(vv[-201] * 4 + 0.5 - 0.0, 0.18, vv[-1] * 4 + 0.5, 0.26,
       col = cols, border = NA)
  rect(0.1, 0.18, 0.9, 0.26, border = Col$muted)
  at = c(-0.1, -0.05, 0, 0.05, 0.1)
  segments(at * 4 + 0.5, 0.18, at * 4 + 0.5, 0.15, col = Col$muted)
  text(at * 4 + 0.5, 0.11, format(at), cex = 0.85, col = Col$muted)
  text(0.5, 0.03, expression("rain minus the board's " * rho[50] *
                               "(200)"), cex = 0.9, col = Col$ink2)
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Plot_E2_window: rho_50 at windows of 50 to 300 years, with
# bootstrap bands
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_E2_window = function() {
  r = Res("E2/rho50")
  W = c(50, 100, 150, 200, 300)
  op = par(mar = c(3.6, 4.6, 2.2, 10.5))
  on.exit(par(op))
  Panel(c(46, 320), c(0.722, 0.849), log = "x", x_at = W,
        main = expression(bold(rho[50] *
                                 " against the waiting window")),
        xlab = "window: years waited (log scale)",
        ylab = expression(rho[50] * ": rain at which half the runs " *
                            "collapse"))
  for (b in 1:7) {
    s = r[r$board == b, ]
    s = s[order(s$window), ]
    polygon(c(s$window, rev(s$window)), c(s$lo, rev(s$hi)),
            col = adjustcolor(Board_col[b], 0.1), border = NA)
    lines(s$window, s$rho50, col = Board_col[b], lwd = 2)
    Board_points(s$window, s$rho50, b, cex = 1.1)
  }
  legend(10^par("usr")[2], mean(par("usr")[3:4]), Board_label,
         col = Board_col, pch = Board_pch, pt.bg = Board_col, lwd = 2,
         pt.lwd = ifelse(Board_pch >= 21, 0.8, 2), bty = "n",
         cex = 0.9, xpd = NA, yjust = 0.5, xjust = 0,
         inset = 0)
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Plot_E2_long: tipping curves of one board from 400 runs of 2,000
# years, at windows of 200 to 2,000 years. Longer windows move the
# switch up slowly and keep it sharp.
#
#   b   1 = Real, 2 = Shuffled
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_E2_long = function(b = 1) {
  d = Res("E2/long_runs")
  d = d[d$board == b, ]
  W = c(200, 300, 500, 1000, 2000)
  ramp = c("#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#0d366b")
  xs = sort(unique(d$rain))
  xx = seq(min(xs) - 0.005, max(xs) + 0.005, length.out = 400)
  op = par(mar = c(3.6, 4.4, 2.2, 1))
  on.exit(par(op))
  Panel(range(xx), c(-0.03, 1.03), xaxs = "i",
        main = sprintf("%s board, %d runs per rain, 2,000 years",
                       Board_label[b],
                       round(nrow(d) / length(unique(d$rain)))),
        xlab = expression("rain " * rho),
        ylab = if (b == 1) "share of runs collapsed" else "")
  labs = c()
  for (k in seq_along(W)) {
    y = !is.na(d$t_ext) & d$t_ext <= W[k]
    sh = sapply(xs, function(v) mean(y[d$rain == v]))
    fit = Fit_logistic(d$rain, y)
    w = max(fit[["width"]], 1e-4)
    lines(xx, 1 / (1 + exp((xx - fit[["rho50"]]) / w)), col = ramp[k],
          lwd = 1.8)
    points(xs, sh, pch = 21, bg = ramp[k], col = "white", cex = 0.9)
    labs = c(labs, as.expression(bquote(
      .(sprintf("%d years: ", W[k])) * rho[50] *
        .(sprintf(" %.3f, 90 to 10%% over %.3f", fit[["rho50"]],
                  2 * log(9) * fit[["width"]])))))
  }
  legend("topright", legend = labs, col = ramp, lwd = 2, bty = "n",
         cex = 0.8, title = "window", title.col = Col$ink)
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# E2 window theory: one sealed cell, solved exactly (Exact_P_empty()):
# the chance of being empty after windows of 10 to 10^8 years, for
# cells of about 14, 28, 71 and 142 hares (cap K = 20, 40, 100, 200)
#
#   Plot_E2_theory_curves()    the four sets of curves
#   Plot_E2_theory_position()  food at which half the cells are empty
#   Plot_E2_theory_width()     width of the switch (90% to 10% empty)
# ----------------------------------------------------------------------
E2_theory_blues = c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5",
                    "#256abf", "#184f95", "#104281", "#0d366b")

E2_theory_label = function(T) {
  if (T < 1e4) return(as.expression(sprintf("%g yr", T)))
  return(as.expression(bquote(10^.(round(log10(T))) ~ "yr")))
}

E2_theory_switch = function() {
  wt = Res("E2/window_theory")
  out = data.frame()
  for (K in sort(unique(wt$K))) {
    for (T in sort(unique(wt$window))) {
      s = wt[wt$K == K & wt$window == T, ]
      s = s[order(s$f), ]
      if (min(s$P_empty) < 0.1 && max(s$P_empty) > 0.9) {
        f_at = function(P) approx(rev(s$P_empty), rev(s$f), P,
                                  ties = mean)$y
        out = rbind(out, data.frame(K = K, window = T, f50 = f_at(0.5),
                                    span = f_at(0.1) - f_at(0.9)))
      }
    }
  }
  return(out)
}

Plot_E2_theory_curves = function() {
  wt = Res("E2/window_theory")
  Ks = sort(unique(wt$K))
  Ts = sort(unique(wt$window))
  op = par(mfrow = c(1, 4), mar = c(3.6, 1.2, 2.4, 0.8),
           oma = c(0, 3.4, 0, 0))
  on.exit(par(op))
  par(cex = 0.9)
  for (j in seq_along(Ks)) {
    s = wt[wt$K == Ks[j], ]
    Panel(c(0.4, 1.6), c(-0.03, 1.03), xaxs = "i",
          main = sprintf("cap K = %d (about %d hares)", Ks[j],
                         round(14.2 * Ks[j] / 20)),
          xlab = "food f", y_lab = j == 1,
          ylab = "")
    abline(v = s$fold[1], col = Col$ink2, lwd = 1, lty = 2)
    for (k in seq_along(Ts)) {
      r = s[s$window == Ts[k], ]
      r = r[order(r$f), ]
      lines(r$f, r$P_empty, col = E2_theory_blues[k], lwd = 1.8)
    }
    if (j == 1) {
      mtext("P(cell empty by the window's end)", side = 2, line = 2.6,
            col = Col$ink2, cex = 0.9)
      text(s$fold[1] + 0.01, 0.98, "fold", adj = c(0, 1), cex = 0.8,
           col = Col$ink2)
      legend("bottomleft", legend = do.call(c, lapply(Ts,
                                                      E2_theory_label)),
             col = E2_theory_blues, lwd = 2, bty = "n", cex = 0.75,
             title = "window", title.col = Col$ink)
    }
  }
  return(invisible(NULL))
}

Plot_E2_theory_position = function() {
  sw = E2_theory_switch()
  op = par(mar = c(3.6, 4.6, 2.2, 1))
  on.exit(par(op))
  x_at = 10^(2:8)
  Panel(c(30, 2e8), range(sw$f50) + c(-0.04, 0.04), log = "x",
        main = "Position: moves with the window", x_at = x_at,
        x_lab = Log_labels(x_at), xlab = "window (years, log scale)",
        ylab = "food at which half the cells are empty")
  Ks = sort(unique(sw$K))
  for (j in seq_along(Ks)) {
    s = sw[sw$K == Ks[j], ]
    lines(s$window, s$f50, col = Board_col[j], lwd = 2)
    Board_points(s$window, s$f50, j, cex = 1.1)
  }
  legend("topleft", paste("K =", Ks), col = Board_col[seq_along(Ks)],
         pch = Board_pch[seq_along(Ks)], pt.bg = Board_col[
           seq_along(Ks)], pt.lwd = 0.8, lwd = 2, bty = "n", cex = 0.9)
  return(invisible(NULL))
}

Plot_E2_theory_width = function() {
  sw = E2_theory_switch()
  op = par(mar = c(3.6, 4.6, 2.2, 1))
  on.exit(par(op))
  x_at = 10^(2:8)
  Panel(c(30, 2e8), c(0, 0.36), log = "x",
        main = "Width: set by population size, not the window",
        x_at = x_at, x_lab = Log_labels(x_at),
        xlab = "window (years, log scale)",
        ylab = "width of the switch (food, 90% to 10%)")
  Ks = sort(unique(sw$K))
  for (j in seq_along(Ks)) {
    s = sw[sw$K == Ks[j], ]
    lines(s$window, s$span, col = Board_col[j], lwd = 2)
    Board_points(s$window, s$span, j, cex = 1.1)
  }
  legend("topright", paste("K =", Ks), col = Board_col[seq_along(Ks)],
         pch = Board_pch[seq_along(Ks)], pt.bg = Board_col[
           seq_along(Ks)], pt.lwd = 0.8, lwd = 2, bty = "n", cex = 0.9)
  return(invisible(NULL))
}

# ---- nicknames: one line each, no input (README.Rmd calls these) -----
E2_tipping_plot_1 = function() Plot_E2_tipping()
E2_tipping_plot_2 = function() Plot_E2_rho50()
E2_fold_gap_plot = function() Plot_E2_fold_gap()
E2_hares_vs_det_plot = function() Plot_E2_hares_vs_det()
E2_mean_N_plot = function() Plot_E2_mean_N()
E2_story_plot = function() Plot_E2_story()
E2_pullback_plot = function() Plot_E2_pullback()
E2_refuge_plot = function() Plot_E2_refuges()
E2_survival_plot = function() Plot_E2_survival()
E2_window_plot = function() Plot_E2_window()
E2_long_plot_1 = function() Plot_E2_long(1)
E2_long_plot_2 = function() Plot_E2_long(2)
E2_theory_plot_1 = function() Plot_E2_theory_curves()
E2_theory_plot_2 = function() Plot_E2_theory_position()
E2_theory_plot_3 = function() Plot_E2_theory_width()


