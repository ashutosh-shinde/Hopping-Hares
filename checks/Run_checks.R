# ======================================================================
# Run_checks.R   run every check (a few minutes)
# ======================================================================
#
#   Run from the Hopping_hare_final folder after any change to the
#   code:   source("checks/Run_checks.R")
#   Every line should say PASS. The exact checks compare a simulated
#   mean with the exact one in units of its standard error (pass:
#   below 3); every check uses fixed seeds, so the verdicts never
#   change. The engine checks need the C engine (src/hh_engine.c);
#   without a C compiler they report FAIL and the R engine still
#   works.
# ======================================================================

source("Load.R")
source(file.path("checks", "Checks.R"))
p = Set_params()
t0 = Sys.time()

ok = c(Check_map(p),
       Check_rho(),
       Check_switch(p),
       Check_conservation(p),
       Check_split(p),
       Check_alpha(p),
       Check_seed(p),
       Check_engines(p),
       Check_book(p),
       Check_onecell(p),
       Check_exact(p),
       Check_csv(),
       Check_place(),
       Check_fc(p),
       Check_states(p),
       Check_boards(),
       Check_groups(p),
       Check_det_jacobian(p),
       Check_det_eigen(p),
       Check_det_onecell(p),
       Check_det_fold(p),
       Check_window_theory(),
       Check_results(),
       Check_guard_params(),
       Check_guard_grid(),
       Check_guard_rain(p),
       Check_guard_n0(p),
       Check_guard_tmax(p))

cat(sprintf("\n%d of %d checks passed in %.0f seconds\n", sum(ok),
            length(ok), as.numeric(Sys.time() - t0, units = "secs")))
