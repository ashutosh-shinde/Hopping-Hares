# ======================================================================
# Rates.R   event rates of one cell, per year
# ======================================================================
#
#   Outcome page, section 1.2 (step 2 of the build chart, 1.1): the
#   four events as rates. The Gillespie simulator (Simulate.R and
#   src/hh_engine.c) fires them; the deterministic board
#   (Deterministic.R) averages them.
#
#   birth   B(n) = b0 * f * n * phi(n) * s(n)
#   death   D(n) = mu * n
#   hop     H_d  = (m / 4) * a_d * n      d = up, down, left, right
#
#   phi(n) = (n - 1) / (n - 1 + h)    partner chance (0 if n <= 1)
#   s(n)   = max(0, 1 - n / K)        free space (0 if n >= K)
#   f      = rho * q                  food in the cell
#
#   Example: a cell of grass 0.80 at rain 1.0 with 10 hares has
#   B = 4 * 0.8 * 10 * (9 / 11) * 0.5 = 13.1 births and D = 10 deaths
#   a year (the centre cell of the picture in section 1.1).
#
#   Partner_chance()  Free_space()  Birth_rate()  Death_rate()
#   Hop_rates()       Cell_rates()
# ======================================================================

# ----------------------------------------------------------------------
# Partner_chance: chance that a hare in a cell of n hares has a partner
#
#   phi(n) = (n - 1) / (n - 1 + h),   phi(0) = 0
#
#   n   hares in the cell      0, 1, 2, ...
#   p   parameters             uses h
#
#   returns  phi, one number from 0 to 1
# ----------------------------------------------------------------------
Partner_chance = function(n, p) {
  if (n < 1) return(0)
  phi = (n - 1) / (n - 1 + p$h)
  return(phi)
}

# ----------------------------------------------------------------------
# Free_space: share of the cell's room still free
#
#   s(n) = max(0, 1 - n / K)
#
#   n   hares in the cell      0, 1, 2, ... (hops can push n past K)
#   p   parameters             uses K
#
#   returns  s, one number from 0 to 1
# ----------------------------------------------------------------------
Free_space = function(n, p) {
  s = max(0, 1 - n / p$K)
  return(s)
}

# ----------------------------------------------------------------------
# Birth_rate: births per year in one cell
#
#   B = b0 * f * n * phi(n) * s(n)
#
#   n   hares in the cell      0, 1, 2, ...
#   f   food, rho * q          0 to about 1.3
#   p   parameters             uses b0, h, K
#
#   returns  B, one number (births per year)
# ----------------------------------------------------------------------
Birth_rate = function(n, f, p) {
  B = p$b0 * f * n * Partner_chance(n, p) * Free_space(n, p)
  return(B)
}

# ----------------------------------------------------------------------
# Death_rate: deaths per year in one cell
#
#   D = mu * n
#
#   n   hares in the cell      0, 1, 2, ...
#   p   parameters             uses mu
#
#   returns  D, one number (deaths per year)
# ----------------------------------------------------------------------
Death_rate = function(n, p) {
  D = p$mu * n
  return(D)
}

# ----------------------------------------------------------------------
# Hop_rates: hops per year from cell i toward each of its 4 sides
#
#   H_d = (m / 4) * a_d * n
#
#   i     cell number            1 to C
#   n     hares in the cell      0, 1, 2, ...
#   map   map with a, the C x 4 crossing chances (Run_sim adds it)
#   p     parameters             uses m
#
#   returns  H, 4 numbers: up, down, left, right (hops per year)
# ----------------------------------------------------------------------
Hop_rates = function(i, n, map, p) {
  H = (p$m / 4) * map$a[i, ] * n
  return(H)
}

# ----------------------------------------------------------------------
# Cell_rates: all six event rates of cell i
#
#   c(B, D, H_up, H_down, H_left, H_right)
#
#   i     cell number            1 to C
#   n     hares in the cell      0, 1, 2, ...
#   f     food in the cell, rho * q
#   map   map with a
#   p     parameters
#
#   returns  6 numbers (events per year)
# ----------------------------------------------------------------------
Cell_rates = function(i, n, f, map, p) {
  r = c(Birth_rate(n, f, p), Death_rate(n, p), Hop_rates(i, n, map, p))
  return(r)
}
