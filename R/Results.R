# ======================================================================
# Results.R   read the saved results (the plots draw from these)
# ======================================================================
#
#   The run files (run/) write CSV tables into results/; the plot
#   files (plots/) read them back through Res(). Each table is read
#   once per session and kept in memory. Big tables are gzipped
#   (.csv.gz); read.csv() opens them directly.
#
#   Where results live: options(hh.results = "results") by default.
#   A quick run writes to "results_quick"; set the option to that
#   folder to draw its figures instead.
#
#   Every table names cells by grid row and column (row 1 = north),
#   never by a cell number, so no cell order can be mixed up.
#
#   Res()  Res_dir()  Res_write()  Res_clear()  Wide_years()
# ======================================================================

Res_cache = new.env()

# ----------------------------------------------------------------------
# Res_dir: the folder the results are read from and written to
#
#   returns  a folder path
# ----------------------------------------------------------------------
Res_dir = function() {
  return(getOption("hh.results", "results"))
}

# ----------------------------------------------------------------------
# Res: one saved table, read once and kept
#
#   name   path inside the results folder without ".csv" or
#          ".csv.gz", e.g. "E2/runs"
#
#   returns  data frame
# ----------------------------------------------------------------------
Res = function(name) {
  key = file.path(Res_dir(), name)
  if (!is.null(Res_cache[[key]])) return(Res_cache[[key]])
  file = paste0(key, ".csv")
  if (!file.exists(file)) file = paste0(key, ".csv.gz")
  if (!file.exists(file)) {
    stop("Res: no table '", name, "' in ", Res_dir(), "/; run the ",
         "matching file in run/ first")
  }
  d = read.csv(file, stringsAsFactors = FALSE)
  assign(key, d, envir = Res_cache)
  return(d)
}

# ----------------------------------------------------------------------
# Res_write: save one table into the results folder
#
#   d      data frame
#   name   as in Res(); gzip = TRUE adds ".csv.gz"
#
#   returns  the file path (invisibly)
# ----------------------------------------------------------------------
Res_write = function(d, name, gzip = FALSE) {
  file = file.path(Res_dir(), paste0(name, if (gzip) ".csv.gz" else
    ".csv"))
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  con = if (gzip) gzfile(file, "w") else file(file, "w")
  write.csv(d, con, row.names = FALSE, na = "NA")
  close(con)
  key = file.path(Res_dir(), name)
  if (!is.null(Res_cache[[key]])) rm(list = key, envir = Res_cache)
  return(invisible(file))
}

# ----------------------------------------------------------------------
# Res_clear: forget every table read so far (after new runs)
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Res_clear = function() {
  rm(list = ls(Res_cache), envir = Res_cache)
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Wide_years: the yearly columns of a table as a matrix
#
#   d        data frame with columns N_0, N_1, ... (or Y_0, ...)
#   prefix   "N_" or "Y_"
#
#   returns  matrix, one row per table row, one column per year
# ----------------------------------------------------------------------
Wide_years = function(d, prefix = "N_") {
  cols = grep(paste0("^", prefix, "[0-9]+$"), names(d), value = TRUE)
  yrs = as.integer(sub(prefix, "", cols))
  m = as.matrix(d[, cols[order(yrs)], drop = FALSE])
  colnames(m) = sort(yrs)
  return(m)
}
