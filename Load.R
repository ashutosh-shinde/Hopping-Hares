# ======================================================================
# Load.R   read every function file
# ======================================================================
#
#   Run from this folder (Hopping_hare_final). In RStudio: Session >
#   Set Working Directory > To Source File Location. Base R only.
#
#     source("Load.R")
#     deterministic_fold_plot_2()      # any nickname draws a figure
#
#   R/       the model, the deterministic board and the measures
#   plots/   figure functions and their no-input nicknames
#   The run files (run/) and the checks (checks/) are not loaded
#   here; source them on their own.
# ======================================================================

if (!file.exists(file.path("R", "Params.R"))) {
  stop("Load.R: set the working directory to the Hopping_hare_final ",
       "folder first (it holds R/, plots/ and results/)")
}
for (load_file in c("Params.R", "Map.R", "Rates.R", "Rho.R",
                    "Simulate.R", "Engine.R", "Deterministic.R",
                    "Measure.R", "Results.R")) {
  source(file.path("R", load_file))
}
for (load_file in c("Plots.R", "Plot_deterministic.R", "Plot_E1.R",
                    "Plot_E2.R", "Plot_E3.R", "Plot_E4.R")) {
  source(file.path("plots", load_file))
}
rm(load_file)
