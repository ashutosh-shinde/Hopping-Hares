# ======================================================================
# Checks.R   tests of the code against answers found another way
# ======================================================================
#
#   Each Check_ function prints one line, PASS or FAIL, and returns
#   TRUE or FALSE. Run_checks.R runs them all (a few minutes).
#
#   Structure, on a 5 x 5 test map (grass 1.0 down to 0.3, a wall in
#   the centre): links both ways, a_ij + a_ji = 1, conservation, the
#   q^alpha split, alpha set at run time, seeds, rate bookkeeping.
#   Exact answers (outcome page, section 1.6): mean years to
#   extinction against Exact_T_one() (1 cell) and Exact_T_two() (2
#   cells), with both engines; a 5 x 5 map has no exact answer.
#   Pass: simulated mean within 3 standard errors.
#   Engines: with the same seed the C and R engines give the same run.
#   Rain: lookups, monthly steps, the switch at a step boundary.
#   Deterministic board: Jacobian against finite differences, the
#   symmetric eigenvalue shortcut, the one-cell fold, the saved fold
#   of the Shuffled board.
#   Data: the grid, the seven boards, the saved results.
#   Guards: bad inputs must stop with an error.
#
#   Report()  Fails()  Test_map_5x5()  Tv_distance()
#   Check_map()  Check_rho()  Check_switch()  Check_conservation()
#   Check_split()  Check_alpha()  Check_seed()  Check_book()
#   Check_onecell()  Check_exact()  Check_engines()  Check_csv()
#   Check_place()  Check_fc()  Check_states()  Check_boards()
#   Check_groups()  Check_det_jacobian()  Check_det_eigen()
#   Check_det_onecell()  Check_det_fold()  Check_window_theory()
#   Check_results()  Check_guard_params()  Check_guard_grid()
#   Check_guard_rain()  Check_guard_n0()  Check_guard_tmax()
# ======================================================================

# ----------------------------------------------------------------------
# Report: print one PASS or FAIL line
#
#   name     check name
#   ok       TRUE or FALSE
#   detail   short text with the numbers behind the verdict
#
#   returns  ok, invisibly
# ----------------------------------------------------------------------
Report = function(name, ok, detail) {
  cat(sprintf("%s  %-20s %s\n", if (ok) "PASS" else "FAIL", name,
              detail))
  return(invisible(ok))
}

# ----------------------------------------------------------------------
# Fails: does an expression stop with an error?
#
#   expr   any R expression; it runs here, inside tryCatch()
#
#   returns  TRUE if it stopped with an error, FALSE if it ran
# ----------------------------------------------------------------------
Fails = function(expr) {
  failed = tryCatch({
    force(expr)
    FALSE
  }, error = function(e) {
    return(TRUE)
  })
  return(failed)
}

# ----------------------------------------------------------------------
# Test_map_5x5: the fixed 5 x 5 test map
#
#   grass from 1.0 (north-west corner) down to 0.3 (south-east), one
#   wall in the centre; 24 cells with 2, 3 or 4 neighbours
#
#   returns  map from Build_map()
# ----------------------------------------------------------------------
Test_map_5x5 = function() {
  g = matrix(c(1.0, 0.9, 0.8, 0.7, 0.6,
               0.9, 0.8, 0.7, 0.6, 0.5,
               0.8, 0.7,  NA, 0.5, 0.4,
               0.7, 0.6, 0.5, 0.4, 0.3,
               0.6, 0.5, 0.4, 0.3, 0.3), 5, 5, byrow = TRUE)
  map = Build_map(g)
  return(map)
}

# ----------------------------------------------------------------------
# Tv_distance: total-variation distance between two sets of shares
#
#   d = 0.5 * sum |x - y|   (0 = the same, 1 = no overlap)
#
#   x, y   shares that each add up to 1, same length
#
#   returns  d, one number from 0 to 1
# ----------------------------------------------------------------------
Tv_distance = function(x, y) {
  d = 0.5 * sum(abs(x - y))
  return(d)
}

# ----------------------------------------------------------------------
# Check_map: the map tables
#
#   (a) on the 5 x 5 map: if j is i's neighbour on one side, i is j's
#       neighbour on the opposite side, and a_ij + a_ji = 1
#   (b) the hand-worked 3 x 3 case: a cell of grass 0.4 with 0.8
#       above, 0.4 below, 0.2 on the left and a wall on the right
#       has a = 0.8, 0.5, 0.2, 0 at alpha = 2
#
#   p   parameters             uses alpha
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_map = function(p) {
  map = Test_map_5x5()
  a = Hop_table(map, p$alpha)
  back = c(2, 1, 4, 3)                   # opposite of up, down, l, r
  ok_link = TRUE
  worst = 0
  n_link = 0
  for (i in seq_len(map$C)) {
    for (d in 1:4) {
      j = map$nbr[i, d]
      if (j == 0) next
      n_link = n_link + 1
      if (map$nbr[j, back[d]] != i) ok_link = FALSE
      worst = max(worst, abs(a[i, d] + a[j, back[d]] - 1))
    }
  }
  g3 = matrix(c(0.5, 0.8, 0.5,
                0.2, 0.4,  NA,
                0.5, 0.4, 0.5), 3, 3, byrow = TRUE)
  m3 = Build_map(g3)
  a3 = Hop_table(m3, 2)
  ok_hand = isTRUE(all.equal(unname(a3[m3$id[2, 2], ]),
                             c(0.8, 0.5, 0.2, 0)))
  ok = ok_link && worst < 1e-12 && ok_hand
  Report("Check_map", ok,
         sprintf("%d links both ways, max |a_ij + a_ji - 1| %.1e; %s",
                 n_link, worst,
                 if (ok_hand) "3 x 3 case matches" else
                   "3 x 3 case WRONG"))
  return(ok)
}

# ----------------------------------------------------------------------
# Check_rho: the rain helpers and the rain lookup
#
#   lengths (yearly and monthly steps); Rho_at() at step boundaries
#   and after the last step; the E3 ramp takes the line's value at
#   the middle of each month; drought placement
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_rho = function() {
  ok_len = length(Rho_const(1, 10)) == 10 &&
    length(Rho_ramp(1.2, 0.3, 150)) == 150 &&
    length(Rho_ramp(1.0, 0.4, 50, 1 / 12)) == 600 &&
    length(Rho_shock(1, 0.5, 5, 3, 20)) == 20
  x = Rho_at(c(0, 0.99, 1, 1.5, 2, 9), c(1, 2, 3))
  ok_at = identical(x, c(1, 1, 2, 2, 3, 3))
  y = Rho_at(c(0, 1 / 12, 0.99 / 12), c(5, 6), dt_rain = 1 / 12)
  ok_month = identical(y, c(5, 6, 5))
  r = Rho_ramp(1.0, 0.4, 50, 1 / 12)
  ok_ramp = abs(r[1] - (1 - 0.6 / 50 / 24)) < 1e-12 &&
    abs(r[600] - (0.4 + 0.6 / 50 / 24)) < 1e-12 &&
    max(abs(diff(r) + 0.001)) < 1e-12    # 0.001 of rain per month
  s1 = Rho_shock(1, 0.4, start = 2, dry_years = 3, years = 8)
  ok_one = identical(s1, c(1, 1, 0.4, 0.4, 0.4, 1, 1, 1))
  ok = ok_len && ok_at && ok_month && ok_ramp && ok_one
  Report("Check_rho", ok,
         "lengths, lookup at step edges, monthly ramp, drought")
  return(ok)
}

# ----------------------------------------------------------------------
# Check_switch: the rain changes exactly at the step boundary, with
# both engines
#
#   One cell, no deaths, no hops, no cap: hares can only be born.
#   Rain 1, 1, then 0 from step 3: the count must still grow in step
#   2 and freeze from t = 2 steps on. Yearly and monthly steps.
#
#   p   parameters
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_switch = function(p) {
  p0 = p
  p0$mu = 0
  p0$m = 0
  p0$K = 1e6                             # no cap: growth never stalls
  map = Build_map(matrix(1))
  ok = TRUE
  msg = c()
  for (eng in c("R", "C")) {
    for (dt in c(1, 1 / 12)) {
      res = Run_sim(p0, map, rho = c(1, 1, 0) * (if (dt < 1) 12 else 1),
                    n0 = 3, t_max = 5 * dt, dt_rain = dt,
                    dt_rec = dt / 100, seed = 9, engine = eng,
                    record = FALSE)
      k1 = which.min(abs(res$t - dt))      # record nearest step 1
      k2 = which.min(abs(res$t - 2 * dt))  # record nearest step 2
      after = res$N[k2:length(res$N)]
      ok = ok && res$N[k2] > res$N[k1] && all(after == res$N[k2])
      msg = c(msg, sprintf("%s %s: %d -> %d", eng,
                           if (dt < 1) "monthly" else "yearly",
                           res$N[k1], res$N[k2]))
    }
  }
  Report("Check_switch", ok, paste(msg, collapse = "; "))
  return(ok)
}

# ----------------------------------------------------------------------
# Check_conservation: hops alone never change the total
#
#   5 x 5 map, b0 = mu = 0, 5 hares per cell, 20 years, both engines
#
#   p   parameters
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_conservation = function(p) {
  p0 = p
  p0$b0 = 0
  p0$mu = 0
  map = Test_map_5x5()
  ok = TRUE
  hops = c()
  for (eng in c("R", "C")) {
    res = Run_sim(p0, map, rho = 1, n0 = 5, t_max = 20, seed = 10,
                  engine = eng)
    ev = res$events
    ok = ok && all(res$N == 5 * map$C) && ev[["birth"]] == 0 &&
      ev[["death"]] == 0 && ev[["hop"]] > 0
    hops = c(hops, ev[["hop"]])
  }
  Report("Check_conservation", ok,
         sprintf("total stayed at %d through %d (R) and %d (C) hops",
                 5 * map$C, hops[1], hops[2]))
  return(ok)
}

# ----------------------------------------------------------------------
# Check_split: hops alone spread hares as q^alpha
#
#   share of cell i = q_i^alpha / sum(q^alpha)   (detailed balance:
#   share_i a_ij = share_j a_ji for every link)
#   5 x 5 map, b0 = mu = 0, 10 hares per cell, 200 years, first 20
#   skipped. Pass: total-variation distance below 0.03.
#
#   p   parameters             uses alpha
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_split = function(p) {
  p0 = p
  p0$b0 = 0
  p0$mu = 0
  map = Test_map_5x5()
  res = Run_sim(p0, map, rho = 1, n0 = 10, t_max = 200, seed = 12,
                engine = "C")
  keep = res$t >= 20
  share = colMeans(res$n[keep, ]) / sum(res$n0)
  expect = map$q^p0$alpha / sum(map$q^p0$alpha)
  tv = Tv_distance(share, expect)
  ok = tv < 0.03
  Report("Check_split", ok,
         sprintf("shares vs q^alpha / sum(q^alpha): distance %.4f",
                 tv))
  return(ok)
}

# ----------------------------------------------------------------------
# Check_alpha: the run uses the alpha in p, not one fixed in the map
#
#   One map, built once, run with alpha = 0 (shares should be equal)
#   and alpha = 2 (shares q^2 / sum(q^2)); hops only, 100 years,
#   first 10 skipped. Pass: both distances below 0.03.
#
#   p   parameters
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_alpha = function(p) {
  map = Test_map_5x5()
  tv = c(0, 0)
  for (k in 1:2) {
    p0 = p
    p0$b0 = 0
    p0$mu = 0
    p0$alpha = c(0, 2)[k]
    res = Run_sim(p0, map, rho = 1, n0 = 10, t_max = 100, seed = 13,
                  engine = "C")
    keep = res$t >= 10
    share = colMeans(res$n[keep, ]) / sum(res$n0)
    expect = map$q^p0$alpha / sum(map$q^p0$alpha)
    tv[k] = Tv_distance(share, expect)
  }
  ok = all(tv < 0.03)
  Report("Check_alpha", ok,
         sprintf("same map at alpha 0 and 2: distances %.4f, %.4f",
                 tv[1], tv[2]))
  return(ok)
}

# ----------------------------------------------------------------------
# Check_seed: the same seed gives the same run; another seed does not
#
#   p   parameters
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_seed = function(p) {
  map = Test_map_5x5()
  r1 = Run_sim(p, map, Rho_const(1, 10), n0 = 8, seed = 8)
  r2 = Run_sim(p, map, Rho_const(1, 10), n0 = 8, seed = 8)
  r3 = Run_sim(p, map, Rho_const(1, 10), n0 = 8, seed = 9)
  ok = identical(r1$n, r2$n) && !identical(r1$n, r3$n)
  Report("Check_seed", ok,
         sprintf("seed 8 twice: same %d events; seed 9 differs",
                 sum(r1$events)))
  return(ok)
}

# ----------------------------------------------------------------------
# Check_engines: the C and R engines run the same algorithm on the
# same random numbers, so one seed gives one run in both
#
#   5 x 5 map, a falling monthly ramp, 3 years, with group totals,
#   a snapshot, an occupancy window and the 10% rule all switched on
#
#   p   parameters
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_engines = function(p) {
  map = Test_map_5x5()
  rho = Rho_ramp(1.1, 0.5, 3, dt_rain = 1 / 12)
  groups = list(west = which(map$col <= 2), east = which(map$col > 2))
  run = function(eng) {
    return(Run_sim(p, map, rho, n0 = 8, t_max = 3, dt_rain = 1 / 12,
                   dt_rec = 1 / 12, record = FALSE, snap_t = c(1, 2),
                   occ_window = c(0.5, 2.5), frac = 0.9,
                   groups = groups, seed = 21, engine = eng))
  }
  a = run("R")
  b = run("C")
  ok = identical(as.integer(a$N), as.integer(b$N)) &&
    identical(as.integer(a$G), as.integer(b$G)) &&
    identical(as.integer(a$snap), as.integer(b$snap)) &&
    identical(as.integer(a$n_final), as.integer(b$n_final)) &&
    isTRUE(all.equal(a$occ_time, b$occ_time)) &&
    isTRUE(all.equal(a$t_frac, b$t_frac)) &&
    all(a$events == b$events)
  Report("Check_engines", ok,
         sprintf("seed 21: both engines give the same %d events",
                 sum(a$events)))
  return(ok)
}

# ----------------------------------------------------------------------
# Check_book: the rate bookkeeping (R engine)
#
#   5 x 5 map; births, deaths, hops and a falling rain. After the
#   run, every rate rebuilt from the final board must equal the table
#   the simulator kept updating; births - deaths must equal the
#   change in total; the last food must use the last rain.
#
#   p   parameters
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_book = function(p) {
  map = Test_map_5x5()
  rho = c(Rho_const(1.2, 5), Rho_ramp(1.2, 0.8, 10))
  res = Run_sim(p, map, rho, n0 = 8, seed = 11, engine = "R")
  gap = 0
  for (i in seq_len(map$C)) {
    fresh = Cell_rates(i, res$n_final[i], res$f_final[i], res$map, p)
    gap = max(gap, abs(fresh - res$rates_final[i, ]))
  }
  ev = res$events
  ok_count = sum(res$n_final) - sum(res$n0) ==
    ev[["birth"]] - ev[["death"]]
  ok_rain = res$extinct ||
    isTRUE(all.equal(res$f_final, rho[length(rho)] * map$q))
  ok = gap < 1e-9 && ok_count && ok_rain
  Report("Check_book", ok,
         sprintf("rate gap %.1e; births - deaths = change; %s", gap,
                 "last rain in use"))
  return(ok)
}

# ----------------------------------------------------------------------
# Check_onecell: one isolated cell against the exact mean years to
# extinction (Exact_T_one); food 0.8, 13 hares, both engines
#
#   p        parameters
#   n_runs   number of runs per engine
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_onecell = function(p, n_runs = 300) {
  map = Build_map(matrix(1))             # one cell, no neighbours
  n0 = 13
  exact = Exact_T_one(0.8, p)[n0]
  ok = TRUE
  msg = c()
  for (eng in c("R", "C")) {
    many = Run_many(p, map, rho = 0.8, n0 = n0, n_runs = n_runs,
                    t_max = 2000, dt_rec = NULL, seed = 5,
                    engine = eng)
    te = many$t_ext
    sim = mean(te)
    se = sd(te) / sqrt(n_runs)
    ok = ok && !any(is.na(te)) && abs(sim - exact) / se < 3
    msg = c(msg, sprintf("%s %.1f +/- %.1f", eng, sim, se))
  }
  Report("Check_onecell", ok,
         sprintf("%s; exact %.1f years", paste(msg, collapse = ", "),
                 exact))
  return(ok)
}

# ----------------------------------------------------------------------
# Check_exact: two linked cells (grass 1.0 and 0.6) against the exact
# mean years until both are empty (Exact_T_two); tests births,
# deaths and hops together, both engines
#
#   p        parameters
#   rho      rain
#   n0       starting hares, c(good cell, poor cell)
#   n_runs   number of runs per engine
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_exact = function(p, rho = 0.75, n0 = c(13, 5), n_runs = 300) {
  map = Build_map(matrix(c(1.0, 0.6), nrow = 1))
  exact = Exact_T_two(map, rho, p)[n0[1] + 1, n0[2] + 1]
  ok = TRUE
  msg = c()
  for (eng in c("R", "C")) {
    many = Run_many(p, map, rho = rho, n0 = n0, n_runs = n_runs,
                    t_max = 2000, dt_rec = NULL, seed = 6,
                    engine = eng)
    te = many$t_ext
    sim = mean(te)
    se = sd(te) / sqrt(n_runs)
    ok = ok && !any(is.na(te)) && abs(sim - exact) / se < 3
    msg = c(msg, sprintf("%s %.1f +/- %.1f", eng, sim, se))
  }
  Report("Check_exact", ok,
         sprintf("%s; exact %.1f years", paste(msg, collapse = ", "),
                 exact))
  return(ok)
}

# ----------------------------------------------------------------------
# Check_csv: the Mayureshwar grid
#
#   20 x 20: 399 cells, one wall at row 10, column 12, q from 0.306
#   to 0.962, mean 0.662
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_csv = function() {
  g = Read_grid("data/mayureshwar_q_grid.csv")
  m = Build_map(g)
  w = which(is.na(g), arr.ind = TRUE)
  ok = all(dim(g) == 20) && m$C == 399 && nrow(w) == 1 &&
    all(w[1, ] == c(10, 12)) &&
    isTRUE(all.equal(range(m$q), c(0.306, 0.962))) &&
    round(mean(m$q), 3) == 0.662
  Report("Check_csv", ok,
         sprintf(paste0("%d cells, wall at (10, 12), q %.3f to %.3f, ",
                        "mean %.3f"), m$C, min(m$q), max(m$q),
                 mean(m$q)))
  return(ok)
}

# ----------------------------------------------------------------------
# Check_place: Place_hares() puts hares on the right squares
#
#   20 x 20 map: squares before and after the wall in cell order,
#   and the last square; a wall square must stop with an error.
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_place = function() {
  map = Build_map("data/mayureshwar_q_grid.csv")
  rows = c(5, 15, 20)
  cols = c(5, 15, 20)
  counts = c(3, 10, 7)
  n0 = Place_hares(map, rows, cols, counts)
  ok = sum(n0) == sum(counts) && length(n0) == map$C
  for (s in seq_along(rows)) {
    k = which(map$row == rows[s] & map$col == cols[s])
    ok = ok && n0[k] == counts[s]
  }
  ok = ok && Fails(Place_hares(map, 10, 12, 5))
  Report("Check_place", ok,
         "rows 5, 15, 20 land right (cells after the wall too)")
  return(ok)
}

# ----------------------------------------------------------------------
# Check_fc: F_crit() agrees with a numerical maximum of phi(n) s(n)
# and gives 0.498 at the defaults
#
#   p   parameters
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_fc = function(p) {
  f_c = F_crit(p)
  o = optimize(function(n) {
    return(Partner_chance(n, p) * Free_space(n, p))
  }, c(1, p$K), maximum = TRUE, tol = 1e-10)
  f_opt = p$mu / (p$b0 * o$objective)
  f_def = F_crit(Set_params())
  ok = abs(f_c - f_opt) < 1e-6 && round(f_def, 3) == 0.498
  Report("Check_fc", ok,
         sprintf("f_c = %.5f (n* = %.2f), numerical max %.5f",
                 f_c, o$maximum, f_opt))
  return(ok)
}

# ----------------------------------------------------------------------
# Check_states: one cell's settled and tipping counts
#
#   At food 1.0: x_plus = 14.2 and x_minus = 1.8 (outcome page,
#   section 2.1); both solve b0 f phi(x) s(x) = mu; none exist below
#   the fold; at the fold they meet at n* = 5.5.
#
#   p   parameters
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_states = function(p) {
  s = Cell_states(1.0, p)
  g = function(x) p$b0 * (x - 1) / (x - 1 + p$h) * (1 - x / p$K)
  f_c = F_crit(p)
  near = Cell_states(f_c * (1 + 1e-6), p)
  ok = round(s[2], 1) == 14.2 && round(s[1], 1) == 1.8 &&
    max(abs(g(s) - p$mu)) < 1e-8 &&
    all(Cell_states(f_c * 0.99, p) == 0) && abs(diff(near)) < 0.05 &&
    abs(mean(near) - 5.48) < 0.03
  Report("Check_states", ok,
         sprintf("food 1.0: x- = %.2f, x+ = %.2f; they meet at %.2f",
                 s[1], s[2], mean(near)))
  return(ok)
}

# ----------------------------------------------------------------------
# Check_boards: the seven boards
#
#   data/boards.csv: every board holds the same 400 values, the Real
#   board is the grid with the village cell at the mean (0.662), the
#   good blocks of Block_grid() hold the best values, and
#   Make_boards() follows the same recipe (the same multiset of
#   values, good blocks on top, the gradient's best value at the
#   centre of its first good quadrant).
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_boards = function() {
  B = Read_boards()
  g = Read_grid("data/mayureshwar_q_grid.csv")
  g[is.na(g)] = mean(g, na.rm = TRUE)
  vals = sort(round(g, 4))
  same = all(sapply(B, function(b) identical(sort(round(b, 4)), vals)))
  real = max(abs(B[[1]] - g)) < 1e-4
  good = function(b, w) {
    m = (ceiling(row(b) / w) + ceiling(col(b) / w)) %% 2 == 0
    return(min(b[m]) >= max(b[!m]) - 1e-9)
  }
  blocks = good(B[[3]], 1) && good(B[[4]], 2) && good(B[[5]], 4) &&
    good(B[[6]], 10)
  M = Make_boards(seed = 1)
  made = all(sapply(M, function(b) {
    return(identical(sort(round(b, 4)), vals))
  })) && good(M[[3]], 1) && good(M[[6]], 10) &&
    max(M[[7]][5:6, 5:6]) == max(M[[7]])   # best value at the centre
  ok = same && real && blocks && made
  Report("Check_boards", ok,
         "same 400 values on all 7 boards; blocks; Make_boards recipe")
  return(ok)
}

# ----------------------------------------------------------------------
# Check_groups: Run_many() keeps group totals and the end boards, and
# takes one starting board per run
#
#   Two groups that split the 5 x 5 map must add up to the total in
#   every run at every time. With n0 as a C x n_runs matrix, run r
#   must start from column r; per-run seeds must repeat one run
#   whatever runs next to it.
#
#   p   parameters
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_groups = function(p) {
  map = Test_map_5x5()
  groups = list(west = which(map$col <= 2), east = which(map$col > 2))
  n0 = cbind(rep(4, map$C), rep(0:1, length.out = map$C))  # C x 2
  many = Run_many(p, map, Rho_const(1, 3), n0 = n0, n_runs = 2,
                  seeds = c(15, 16), groups = groups)
  sums = many$G[, "west", ] + many$G[, "east", ]   # n_rec x n_runs
  one = Run_many(p, map, Rho_const(1, 3), n0 = n0[, 2, drop = FALSE],
                 n_runs = 1, seeds = 16, groups = groups)
  ok = all(sums == many$N) && all(many$N[1, ] == colSums(n0)) &&
    all(colSums(many$end) == many$N[nrow(many$N), ]) &&
    identical(one$N[, 1], many$N[, 2])
  Report("Check_groups", ok,
         "groups add up; one board and one seed per run")
  return(ok)
}

# ----------------------------------------------------------------------
# Check_det_jacobian: Det_jacobian() against finite differences of
# Det_change() on the 5 x 5 map (a random board of 0-18 hares)
#
#   p   parameters
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_det_jacobian = function(p) {
  db = Det_board(Test_map_5x5(), p)
  x = c(0.5, 1.5, 3, 6, 9, 12, 15, 18, 2.5, 4, 7, 11, 13, 16, 1.2,
        5, 8, 10, 14, 17, 2, 3.5, 6.5, 9.5)
  J = Det_jacobian(x, 0.9, db)
  num = matrix(0, db$C, db$C)
  for (j in seq_len(db$C)) {
    e = rep(0, db$C)
    e[j] = 1e-6
    num[, j] = (Det_change(x + e, 0.9, db) - Det_change(x - e, 0.9,
                                                        db)) / 2e-6
  }
  gap = max(abs(J - num))
  ok = gap < 1e-5
  Report("Check_det_jacobian", ok,
         sprintf("largest gap to finite differences %.1e", gap))
  return(ok)
}

# ----------------------------------------------------------------------
# Check_det_eigen: the symmetric shortcut gives the same slowest rate
# and pattern as R's general eigen() on J itself
#
#   p   parameters
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_det_eigen = function(p) {
  db = Det_board(Test_map_5x5(), p)
  x = Det_start(db, 1.0)
  fast = Det_slowest(db, x, 1.0, vector = TRUE)
  e = eigen(Det_jacobian(x, 1.0, db))
  k = which.max(Re(e$values))
  v = abs(Re(e$vectors[, k]))
  ok = abs(fast$lambda - Re(e$values[k])) < 1e-9 &&
    max(abs(fast$pattern - v / max(v))) < 1e-6
  Report("Check_det_eigen", ok,
         sprintf("slowest rate %.6f both ways", fast$lambda))
  return(ok)
}

# ----------------------------------------------------------------------
# Check_det_onecell: one sealed cell of grass q folds at rain f_c / q
#
#   p   parameters
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_det_onecell = function(p) {
  q = 0.8
  db = Det_board(Build_map(matrix(q)), p)
  fol = Det_follow(db, rho_hi = 0.8, rho_lo = 0.5, step = 0.0025,
                   tol = 1e-7)
  expect = F_crit(p) / q
  ok = nrow(fol$folds) == 1 && abs(fol$folds$rho[1] - expect) < 1e-5
  Report("Check_det_onecell", ok,
         sprintf("fold at rain %.5f; f_c / q = %.5f", fol$folds$rho[1],
                 expect))
  return(ok)
}

# ----------------------------------------------------------------------
# Check_det_fold: the Shuffled board's fold, followed from rain 0.80,
# matches the saved one (results/deterministic/folds.csv)
#
#   p   parameters
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_det_fold = function(p) {
  db = Det_board(Board_map(2), p)
  fol = Det_follow(db, rho_hi = 0.80, rho_lo = 0.70)
  saved = Res("deterministic/folds")
  saved = saved[saved$board == 2, ]
  ok = nrow(fol$folds) == 1 &&
    abs(fol$folds$rho[1] - saved$rho[1]) < 1e-4 &&
    abs(fol$folds$N_before[1] - saved$N_before[1]) < 5
  Report("Check_det_fold", ok,
         sprintf("Shuffled folds at %.4f (saved %.4f), %d hares above",
                 fol$folds$rho[1], saved$rho[1],
                 round(fol$folds$N_before[1])))
  return(ok)
}

# ----------------------------------------------------------------------
# Check_window_theory: Exact_P_empty() matches the saved exact
# one-cell curves (results/E2/window_theory.csv) at a few points
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_window_theory = function() {
  wt = Res("E2/window_theory")
  f = c(0.6, 0.9, 1.2)
  T = c(50, 1e4)
  gap = 0
  for (K in c(20, 40, 100, 200)) {
    P = Exact_P_empty(f, T, Omega = K / 20)
    for (i in seq_along(T)) {
      for (j in seq_along(f)) {
        s = wt[wt$K == K & wt$window == T[i] &
                 abs(wt$f - f[j]) < 1e-9, ]
        gap = max(gap, abs(P[i, j] - s$P_empty))
      }
    }
  }
  ok = gap < 1e-6
  Report("Check_window_theory", ok,
         sprintf("largest gap to the saved curves %.1e", gap))
  return(ok)
}

# ----------------------------------------------------------------------
# Check_results: the saved results hold the runs the outcome page
# reports, and a refit of rho_50 matches the saved one
#
#   E1 1,260 runs; E2 2,890 runs (plus 400 long ones); E3 700 runs;
#   E4 2,000 runs; Real board rho_50(200) = 0.785
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_results = function() {
  n = c(E1 = nrow(Res("E1/runs")), E2 = nrow(Res("E2/runs")),
        E2_long = nrow(Res("E2/long_runs")), E3 = nrow(Res("E3/runs")),
        E4 = nrow(Res("E4/runs")))
  d = Res("E2/runs")
  d = d[d$board == 1, ]
  fit = Fit_logistic(d$rain, !is.na(d$t_ext) & d$t_ext <= 200)
  saved = Res("E2/rho50")
  saved = saved$rho50[saved$board == 1 & saved$window == 200]
  ok = all(n == c(1260, 2890, 400, 700, 2000)) &&
    abs(fit[["rho50"]] - saved) < 1e-4
  Report("Check_results", ok,
         sprintf("runs %s; Real rho50(200) %.4f (saved %.4f)",
                 paste(n, collapse = "/"), fit[["rho50"]], saved))
  return(ok)
}

# ----------------------------------------------------------------------
# Check_guard_params: Set_params() stops on bad values
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_guard_params = function() {
  bad = c(Fails(Set_params(K = 1)), Fails(Set_params(K = 20.5)),
          Fails(Set_params(mu = 0)), Fails(Set_params(b0 = -1)),
          Fails(Set_params(h = 0)), Fails(Set_params(m = -1)),
          Fails(Set_params(alpha = -1)), Fails(Set_params(k = 20)),
          Fails(Set_params(5)), Fails(Set_params(mu = NA)))
  good = !Fails(Set_params(alpha = 0, m = 0, b0 = 0))
  ok = all(bad) && good
  Report("Check_guard_params", ok,
         sprintf("%d of %d bad settings stopped", sum(bad),
                 length(bad)))
  return(ok)
}

# ----------------------------------------------------------------------
# Check_guard_grid: Build_map() and Read_grid() stop on bad grids
#
#   q above 1 or below 0, no usable square, text in the grid
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_guard_grid = function() {
  f = tempfile(fileext = ".csv")
  writeLines(c("0.5,wall", "0.4,0.3"), f)
  bad = c(Fails(Build_map(matrix(c(1.2, 0.5), 1))),
          Fails(Build_map(matrix(c(-0.1, 0.5), 1))),
          Fails(Build_map(matrix(NA, 2, 2))),
          Fails(Build_map(matrix(c("a", "b"), 1))),
          Fails(Read_grid(f)))
  unlink(f)
  ok = all(bad)
  Report("Check_guard_grid", ok,
         sprintf("%d of %d bad grids stopped", sum(bad), length(bad)))
  return(ok)
}

# ----------------------------------------------------------------------
# Check_guard_rain: Run_sim() stops on negative or missing rain, or a
# bad step length
#
#   p   parameters
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_guard_rain = function(p) {
  map = Test_map_5x5()
  bad = c(Fails(Run_sim(p, map, rho = c(1, -0.2), n0 = 5)),
          Fails(Run_sim(p, map, rho = c(1, NA), n0 = 5)),
          Fails(Run_sim(p, map, rho = -0.5, n0 = 5, t_max = 2)),
          Fails(Run_sim(p, map, rho = c(1, 1), n0 = 5, dt_rain = 0)))
  ok = all(bad)
  Report("Check_guard_rain", ok,
         sprintf("%d of %d bad rain inputs stopped", sum(bad),
                 length(bad)))
  return(ok)
}

# ----------------------------------------------------------------------
# Check_guard_n0: Run_sim() stops on a starting board of the wrong
# length or with bad counts
#
#   A 5 x 5 matrix (25 values) on the 24-cell map, two numbers on 24
#   cells, a negative count, half a hare, NA; a correct n0 must run.
#
#   p   parameters
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_guard_n0 = function(p) {
  map = Test_map_5x5()
  rho = Rho_const(1, 2)
  bad = c(Fails(Run_sim(p, map, rho, n0 = matrix(1, 5, 5))),
          Fails(Run_sim(p, map, rho, n0 = c(13, 5))),
          Fails(Run_sim(p, map, rho, n0 = -1)),
          Fails(Run_sim(p, map, rho, n0 = 2.5)),
          Fails(Run_sim(p, map, rho, n0 = NA)))
  good = !Fails(Run_sim(p, map, rho, n0 = rep(1, map$C), t_max = 0.1))
  ok = all(bad) && good
  Report("Check_guard_n0", ok,
         sprintf("%d of %d bad starting boards stopped", sum(bad),
                 length(bad)))
  return(ok)
}

# ----------------------------------------------------------------------
# Check_guard_tmax: runs need a finite end time
#
#   one rain number without t_max, t_max = Inf, a negative t_max,
#   and Run_many() with one rain number must all stop
#
#   p   parameters
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Check_guard_tmax = function(p) {
  map = Test_map_5x5()
  bad = c(Fails(Run_sim(p, map, rho = 1, n0 = 5)),
          Fails(Run_sim(p, map, rho = 1, n0 = 5, t_max = Inf)),
          Fails(Run_sim(p, map, rho = c(1, 1), n0 = 5, t_max = -1)),
          Fails(Run_many(p, map, rho = 1, n0 = 5, n_runs = 2)))
  ok = all(bad)
  Report("Check_guard_tmax", ok,
         sprintf("%d of %d bad end times stopped", sum(bad),
                 length(bad)))
  return(ok)
}
