# ======================================================================
# Engine.R   the optional C engine: compile, load, choose, call
# ======================================================================
#
#   src/hh_engine.c runs the same Gillespie steps as Gillespie_R()
#   (Simulate.R), about 100 times faster. The full experiments need
#   it (E2 alone is 478,000 run-years); the plots do not, because
#   they read the saved results in results/.
#
#   Engine_load() compiles the file once with R CMD SHLIB (it needs a
#   C compiler: Rtools on Windows, the Xcode command line tools on a
#   Mac, gcc on Linux) and loads it. If that fails, every run uses
#   plain R and gives the same results in distribution.
#
#   options(hh.engine = "auto")   C when it loads, else R (default)
#   options(hh.engine = "R")      always plain R
#   options(hh.engine = "C")      C, or stop with an error
#
#   Engine_load()  Engine_ready()  Engine_choose()  Gillespie_C()
# ======================================================================

# ----------------------------------------------------------------------
# Engine_load: compile src/hh_engine.c (if needed) and load it
#
#   src     path of the C file
#   quiet   FALSE = print the compiler's messages when it fails
#
#   returns  TRUE if the engine is loaded, FALSE if not (invisibly)
# ----------------------------------------------------------------------
Engine_load = function(src = file.path("src", "hh_engine.c"),
                       quiet = TRUE) {
  if (Engine_ready()) return(invisible(TRUE))
  if (!file.exists(src)) return(invisible(FALSE))
  lib = sub("[.]c$", .Platform$dynlib.ext, src)
  if (!file.exists(lib) || file.mtime(lib) < file.mtime(src)) {
    r_bin = file.path(R.home("bin"), "R")
    msg = tryCatch(
      suppressWarnings(system2(r_bin, c("CMD", "SHLIB", "-o",
                                        shQuote(lib), shQuote(src)),
                               stdout = TRUE, stderr = TRUE)),
      error = function(e) conditionMessage(e))
    if (!file.exists(lib)) {
      if (!quiet) cat(msg, sep = "\n")
      return(invisible(FALSE))
    }
  }
  ok = tryCatch({
    dyn.load(lib)
    TRUE
  }, error = function(e) FALSE)
  return(invisible(ok && Engine_ready()))
}

# ----------------------------------------------------------------------
# Engine_ready: is the C engine loaded?
#
#   returns  TRUE or FALSE
# ----------------------------------------------------------------------
Engine_ready = function() {
  return(is.loaded("hh_run"))
}

# ----------------------------------------------------------------------
# Engine_choose: which engine a run uses
#
#   engine   "R", "C", "auto" or NULL (= getOption("hh.engine",
#            "auto"))
#
#   returns  "R" or "C"
# ----------------------------------------------------------------------
Engine_choose = function(engine = NULL) {
  if (is.null(engine)) engine = getOption("hh.engine", "auto")
  if (!(engine %in% c("auto", "R", "C"))) {
    stop("Engine_choose: engine must be \"auto\", \"R\" or \"C\"")
  }
  if (engine == "R") return("R")
  if (Engine_ready() || Engine_load()) return("C")
  if (engine == "C") {
    stop("Engine_choose: the C engine did not compile; see ",
         "Engine_load(quiet = FALSE), or use ",
         "options(hh.engine = \"R\")")
  }
  return("R")
}

# ----------------------------------------------------------------------
# Gillespie_C: the C engine with the same inputs and outputs as
# Gillespie_R() (Simulate.R); called by Run_sim()
#
#   returns  a list: N, G, snap, occ, t_frac, end, t_ext, t_end, events
# ----------------------------------------------------------------------
Gillespie_C = function(p, map, rho, dt_rain, n0, t_max, rec_t, snap_t,
                       occ_window, frac, gid, n_g) {
  C = map$C
  n_rec = length(rec_t)
  n_snap = length(snap_t)
  win = if (is.null(occ_window)) c(0, 0) else as.double(occ_window)
  z = .C("hh_run",
         C = as.integer(C), q = as.double(map$q),
         nbr = as.integer(t(map$nbr) - 1L),      # by rows, 0-based
         a = as.double(t(map$a)),                 # by rows
         n0 = as.integer(n0), rho = as.double(rho),
         n_rho = as.integer(length(rho)), dt_rain = as.double(dt_rain),
         t_max = as.double(t_max), rec_t = as.double(rec_t),
         n_rec = as.integer(n_rec), out_N = integer(n_rec),
         group = as.integer(gid), n_group = as.integer(n_g),
         out_G = integer(n_rec * n_g), snap_t = as.double(snap_t),
         n_snap = as.integer(n_snap), out_snap = integer(n_snap * C),
         occ_win = win, out_occ = double(C),
         frac = as.double(if (is.na(frac)) -1 else frac),
         out_tfrac = NA_real_,
         par = as.double(c(p$K, p$mu, p$b0, p$h, p$m)),
         out_end = integer(C), out_text = NA_real_, out_tend = 0,
         out_events = double(3), NAOK = TRUE)
  out = list(N = z$out_N, G = matrix(z$out_G, n_rec, n_g),
             snap = matrix(z$out_snap, n_snap, C), occ = z$out_occ,
             t_frac = z$out_tfrac, end = z$out_end, t_ext = z$out_text,
             t_end = z$out_tend,
             events = c(birth = z$out_events[1],
                        death = z$out_events[2],
                        hop = z$out_events[3]))
  return(out)
}
