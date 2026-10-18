# ======================================================================
# Rho.R   rain, one value per step of dt_rain years
# ======================================================================
#
#   rho[k] is the rain in step k, from time (k - 1) * dt_rain to
#   k * dt_rain. Steps are a year long by default (dt_rain = 1); E3
#   uses monthly steps (dt_rain = 1/12). After the last step the last
#   value holds. Rain scales food in every cell: f_i = rho * q_i, so
#   rho = 1 is a normal year and rho = 0.6 means 60% of normal. Rain
#   must be 0 or more.
#
#   Experiments: E1, E2 and E4 hold the rain fixed (one number and a
#   t_max); E3 lets it fall along a straight line in monthly steps,
#   Rho_ramp(1.0, 0.4, D, dt_rain = 1/12), then holds 0.4.
#
#   Join pieces with c() (same dt_rain throughout):
#     rho = c(Rho_const(1.2, 20), Rho_ramp(1.2, 0.3, 150))
#
#   Rho_const()  Rho_ramp()  Rho_shock()  Rho_at()
# ======================================================================

# ----------------------------------------------------------------------
# Rho_const: the same rain in every step
#
#   value     rain, 0 or more
#   years     number of years
#   dt_rain   step length in years
#
#   returns  rain vector, one value per step (years / dt_rain)
# ----------------------------------------------------------------------
Rho_const = function(value, years, dt_rain = 1) {
  rho = rep(value, round(years / dt_rain))
  return(rho)
}

# ----------------------------------------------------------------------
# Rho_ramp: rain falling (or rising) along a straight line
#
#   Each step takes the line's value at its middle:
#     rho_k = from + (to - from) * (k - 1/2) / n_steps
#   so monthly steps never stray more than half a step from the line
#   (E3: at most 0.0005 of rain on the 50-year ramp).
#
#   from      rain at the start of the line, 0 or more
#   to        rain at its end, 0 or more
#   years     length of the line in years
#   dt_rain   step length in years (1/12 = monthly)
#
#   returns  rain vector, one value per step (years / dt_rain)
# ----------------------------------------------------------------------
Rho_ramp = function(from, to, years, dt_rain = 1) {
  n_steps = round(years / dt_rain)
  rho = from + (to - from) * (seq_len(n_steps) - 0.5) / n_steps
  return(rho)
}

# ----------------------------------------------------------------------
# Rho_shock: steady rain with one or more droughts
#
#   rain = low in years start + 1 to start + dry_years, for each start
#
#   base        rain in normal years
#   low         rain in drought years
#   start       years before each drought: one number or several
#   dry_years   length of each drought in years
#   years       total number of years
#
#   returns  rain vector, length years (yearly steps)
# ----------------------------------------------------------------------
Rho_shock = function(base, low, start, dry_years, years) {
  rho = Rho_const(base, years)
  for (s in start) {
    k = s + seq_len(dry_years)
    rho[k[k <= years]] = low
  }
  return(rho)
}

# ----------------------------------------------------------------------
# Rho_at: the rain in force at times t
#
#   step k = floor(t / dt_rain + 1e-9) + 1; the 1e-9 absorbs rounding
#   at step boundaries. After the last step the last value holds.
#
#   t         times in years, 0 or more
#   rho       rain vector
#   dt_rain   step length in years
#
#   returns  rain at each time, same length as t
# ----------------------------------------------------------------------
Rho_at = function(t, rho, dt_rain = 1) {
  k = floor(t / dt_rain + 1e-9) + 1
  r = rho[pmin(k, length(rho))]
  return(r)
}
