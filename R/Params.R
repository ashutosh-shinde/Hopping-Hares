# ======================================================================
# Params.R   the model parameters and their defaults
# ======================================================================
#
#   Outcome page, section 1.2 (step 2 of the build chart, 1.1).
#
#   name    meaning                                         default
#   K       cap: births stop when a cell holds K hares           20
#   mu      deaths per hare per year                              1
#   b0      top births per hare per year (partner, food 1)        4
#   h       partner constant in phi(n) = (n - 1) / (n - 1 + h)    2
#   m       hop tries per hare per year                          10
#   alpha   pickiness in a_ij = q_j^alpha / (q_i^alpha +          2
#           q_j^alpha); 0 = hares ignore grass when hopping
#
#   Time is in years. A hare lives 1 / mu = 1 year on average.
#   Every experiment (E1-E4) and the deterministic board use the
#   defaults; only the hop sweep (outcome page, section 2.5) changes m.
#
#   Set_params()
# ======================================================================

# ----------------------------------------------------------------------
# Set_params: the default parameters, with any of them changed by name
#
#   p = Set_params()               # defaults
#   p = Set_params(alpha = 0)      # hares ignore grass when hopping
#
#   ...   new values by name: K, mu, b0, h, m, alpha
#
#   Stops on an unknown name or a value out of range:
#   K a whole number, at least 2; mu > 0; b0 >= 0; h > 0; m >= 0;
#   alpha >= 0.
#
#   returns  p, a list with K, mu, b0, h, m, alpha
# ----------------------------------------------------------------------
Set_params = function(...) {
  p = list(K = 20, mu = 1, b0 = 4, h = 2, m = 10, alpha = 2)
  new = list(...)
  if (length(new) > 0 && (is.null(names(new)) ||
                          any(names(new) == ""))) {
    stop("Set_params: give each value a name, e.g. alpha = 0")
  }
  for (nm in names(new)) {
    if (!(nm %in% names(p))) {
      stop("Set_params: unknown parameter '", nm, "'")
    }
    p[[nm]] = new[[nm]]
  }
  for (nm in names(p)) {
    x = p[[nm]]
    if (!is.numeric(x) || length(x) != 1 || !is.finite(x)) {
      stop("Set_params: ", nm, " must be one finite number")
    }
  }
  if (p$K < 2 || p$K != round(p$K)) {
    stop("Set_params: K must be a whole number, at least 2")
  }
  if (p$mu <= 0) stop("Set_params: mu must be above 0")
  if (p$b0 < 0) stop("Set_params: b0 must be 0 or more")
  if (p$h <= 0) stop("Set_params: h must be above 0")
  if (p$m < 0) stop("Set_params: m must be 0 or more")
  if (p$alpha < 0) stop("Set_params: alpha must be 0 or more")
  return(p)
}
