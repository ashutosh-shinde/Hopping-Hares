/* =====================================================================
   hh_engine.c   the Gillespie simulator in C (optional speed-up)

   The same algorithm as Gillespie_R() in R/Simulate.R, step for step,
   about 100 times faster. R/Engine.R compiles this file once with
   R CMD SHLIB and calls hh_run() through .C(). Without a C compiler
   the R version runs instead, with the same results in distribution.

   Random numbers come from R's own generator (unif_rand), so
   set.seed() makes a run repeatable with either engine.

   Outcome page, section 1.6 (steps 10-13 of the build chart, 1.1):
     10. add up the rates of every event: the total rate R
     11. wait an exponential time with mean 1 / R
     12. pick one event with chance proportional to its rate
     13. apply it, update the rates of the cells that changed, repeat
   Rain changes in steps of dt_rain years (monthly in E3): the clock
   stops at each step, the birth rates take the new rain and a fresh
   wait is drawn; waits have no memory, so this stays exact.
   ================================================================== */
#include <R.h>
#include <R_ext/Random.h>
#include <math.h>

/* births per year in one cell: B = b0 f n phi(n) s(n) (Rates.R) */
static double birth_rate(int n, double f, double b0, double h,
                         double K) {
  if (n <= 1) return 0.0;                  /* phi(1) = 0: no partner */
  double phi = (n - 1.0) / (n - 1.0 + h);
  double s = 1.0 - n / K;
  if (s < 0) s = 0;                           /* no room past K */
  return b0 * f * n * phi * s;
}

/* the six rates of one cell: birth, death, hop up, down, left, right */
static void cell_rates(double *r, int n, double f, const double *a4,
                       double b0, double h, double K, double mu,
                       double m) {
  r[0] = birth_rate(n, f, b0, h, K);
  r[1] = mu * n;
  for (int d = 0; d < 4; d++) r[2 + d] = 0.25 * m * a4[d] * n;
}

/* time inside the window [w0, w1) between times a and b */
static double overlap(double a, double b, double w0, double w1) {
  double lo = a > w0 ? a : w0, hi = b < w1 ? b : w1;
  return hi > lo ? hi - lo : 0.0;
}

/* ---------------------------------------------------------------------
   hh_run: one run. Every argument is a pointer (the .C interface).

   C          number of cells
   q          grass per cell (C)
   nbr        neighbour cell, 0-based, C x 4 by rows (up, down, left,
              right); -1 = none
   a          crossing chances, C x 4 by rows
   n0         starting hares (C)
   rho        rain steps (n_rho); step k holds over [k dt, (k+1) dt),
              the last one holds to the end
   t_max      end of the run, years
   rec_t      record times of N (n_rec, increasing)
   group      group of each cell, 1..n_group (0 = none) (C)
   snap_t     times to save the whole board (n_snap, increasing)
   occ_win    window [w0, w1): years each cell holds hares inside it;
              w1 <= w0 switches it off
   frac       first time N < frac N(0) goes to out_tfrac (frac <= 0:
              off)
   par        K, mu, b0, h, m

   outputs: out_N (n_rec), out_G (n_rec x n_group, by records),
   out_snap (n_snap x C, by snapshots), out_occ (C), out_tfrac,
   out_end (C), out_text (extinction time, NA if alive), out_tend,
   out_events (births, deaths, hops)
   ------------------------------------------------------------------ */
void hh_run(int *C_, double *q, int *nbr, double *a, int *n0,
            double *rho, int *n_rho_, double *dt_rain_, double *t_max_,
            double *rec_t, int *n_rec_, int *out_N,
            int *group, int *n_group_, int *out_G,
            double *snap_t, int *n_snap_, int *out_snap,
            double *occ_win, double *out_occ,
            double *frac_, double *out_tfrac,
            double *par,
            int *out_end, double *out_text, double *out_tend,
            double *out_events) {
  const int C = *C_, n_rho = *n_rho_, n_rec = *n_rec_;
  const int n_group = *n_group_, n_snap = *n_snap_;
  const double dt_rain = *dt_rain_, t_max = *t_max_, frac = *frac_;
  const double K = par[0], mu = par[1], b0 = par[2], h = par[3];
  const double m = par[4];
  const double w0 = occ_win[0], w1 = occ_win[1];
  const int do_occ = w1 > w0;

  /* R_alloc memory is freed by R when .C() returns */
  int *n = (int *) R_alloc(C, sizeof(int));
  double *rates = (double *) R_alloc((size_t) C * 6, sizeof(double));
  double *tot = (double *) R_alloc(C, sizeof(double));
  double *f = (double *) R_alloc(C, sizeof(double));
  double *last = (double *) R_alloc(C, sizeof(double));
  long *G_tot = (long *) R_alloc(n_group > 0 ? n_group : 1,
                                 sizeof(long));

  GetRNGstate();
  long N = 0;
  for (int g = 0; g < n_group; g++) G_tot[g] = 0;
  for (int i = 0; i < C; i++) {
    n[i] = n0[i];
    N += n[i];
    if (group[i] > 0) G_tot[group[i] - 1] += n[i];
    last[i] = 0.0;
    out_occ[i] = 0.0;
  }
  const double thresh = frac > 0 ? frac * (double) N : -1.0;
  *out_tfrac = NA_REAL;

  /* rain step in force and the time it ends */
  int k_rho = 0;
  double t_change = n_rho > 1 ? dt_rain : R_PosInf;
  double R = 0.0;
  for (int i = 0; i < C; i++) {
    f[i] = rho[0] * q[i];
    cell_rates(rates + 6 * i, n[i], f[i], a + 4 * i, b0, h, K, mu, m);
    double s = 0;
    for (int e = 0; e < 6; e++) s += rates[6 * i + e];
    tot[i] = s;
    R += s;
  }

  double t = 0.0, t_ext = NA_REAL;
  double ev_b = 0, ev_d = 0, ev_h = 0;
  int k = 0, s_k = 0;
  long ev = 0;
  for (;;) {
    if (N == 0) { t_ext = t; break; }      /* no hare left: the trap */
    if ((ev & 2047) == 0) {                   /* refresh R: no drift */
      R = 0;
      for (int i = 0; i < C; i++) R += tot[i];
    }
    /* 10-11. total rate and the wait */
    double tau = R > 0 ? -log(unif_rand()) / R : R_PosInf;
    double t_next = t + tau;
    double t_stop = t_change < t_max ? t_change : t_max;
    double t_fill = t_next < t_stop ? t_next : t_stop;
    /* records: the board holds until t_fill */
    while (k < n_rec && rec_t[k] < t_fill) {
      out_N[k] = (int) N;
      for (int g = 0; g < n_group; g++)
        out_G[(size_t) g * n_rec + k] = (int) G_tot[g];
      k++;
    }
    while (s_k < n_snap && snap_t[s_k] < t_fill) {
      for (int i = 0; i < C; i++)
        out_snap[(size_t) i * n_snap + s_k] = n[i];
      s_k++;
    }
    /* 11. no event before the rain changes or the run ends: move the
       clock there, take the new rain, draw a fresh wait */
    if (t_next >= t_stop) {
      t = t_stop;
      if (t_stop >= t_max) break;
      k_rho++;
      t_change = k_rho < n_rho - 1 ? (k_rho + 1) * dt_rain : R_PosInf;
      R = 0;
      for (int i = 0; i < C; i++) {
        f[i] = rho[k_rho] * q[i];
        double old_b = rates[6 * i];
        rates[6 * i] = birth_rate(n[i], f[i], b0, h, K);
        tot[i] += rates[6 * i] - old_b;
        R += tot[i];
      }
      continue;
    }
    t = t_next;

    /* 12. pick the cell, then one of its six events */
    double u = unif_rand() * R, cs = 0;
    int ci = -1, last_c = -1;
    for (int i = 0; i < C; i++) {
      if (tot[i] > 0) last_c = i;
      cs += tot[i];
      if (cs > u) { ci = i; break; }
    }
    if (ci < 0) ci = last_c;                  /* rounding at the end */
    double *r = rates + 6 * ci;
    double ue = unif_rand() * tot[ci], ce = 0;
    int e = -1, last_e = -1;
    for (int j = 0; j < 6; j++) {
      if (r[j] > 0) last_e = j;
      ce += r[j];
      if (ce > ue) { e = j; break; }
    }
    if (e < 0) e = last_e;
    int cj = e >= 2 ? nbr[4 * ci + (e - 2)] : -1;

    /* close the occupancy spells of the cells about to change */
    int cells[2] = {ci, cj};
    if (do_occ) {
      for (int z = 0; z < 2; z++) {
        int c = cells[z];
        if (c < 0) continue;
        if (n[c] > 0) out_occ[c] += overlap(last[c], t, w0, w1);
        last[c] = t;
      }
    }

    /* 13. apply the event */
    if (e == 0) {                             /* birth */
      n[ci]++; N++; ev_b++;
      if (group[ci] > 0) G_tot[group[ci] - 1]++;
    } else if (e == 1) {                      /* death */
      n[ci]--; N--; ev_d++;
      if (group[ci] > 0) G_tot[group[ci] - 1]--;
    } else {                                  /* hop ci -> cj */
      n[ci]--; n[cj]++; ev_h++;
      if (group[ci] > 0) G_tot[group[ci] - 1]--;
      if (group[cj] > 0) G_tot[group[cj] - 1]++;
    }
    if (thresh >= 0 && ISNA(*out_tfrac) && (double) N < thresh)
      *out_tfrac = t;

    /* update only the cells that changed */
    for (int z = 0; z < 2; z++) {
      int c = cells[z];
      if (c < 0) continue;
      double old = tot[c];
      cell_rates(rates + 6 * c, n[c], f[c], a + 4 * c, b0, h, K, mu, m);
      double sum = 0;
      for (int j = 0; j < 6; j++) sum += rates[6 * c + j];
      tot[c] = sum;
      R += sum - old;
    }
    ev++;
    if ((ev & 1048575) == 0) R_CheckUserInterrupt();
  }

  /* the last board holds to t_max (all zeros after extinction) */
  while (k < n_rec) {
    out_N[k] = (int) N;
    for (int g = 0; g < n_group; g++)
      out_G[(size_t) g * n_rec + k] = (int) G_tot[g];
    k++;
  }
  while (s_k < n_snap) {
    for (int i = 0; i < C; i++)
      out_snap[(size_t) i * n_snap + s_k] = n[i];
    s_k++;
  }
  for (int i = 0; i < C; i++) {
    if (do_occ && n[i] > 0) out_occ[i] += overlap(last[i], t, w0, w1);
    out_end[i] = n[i];
  }
  *out_text = t_ext;
  *out_tend = t;
  out_events[0] = ev_b;
  out_events[1] = ev_d;
  out_events[2] = ev_h;
  PutRNGstate();
}
