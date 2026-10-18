# ======================================================================
# Plots.R   shared style and drawing helpers for every figure
# ======================================================================
# Plotting code written with AI assistance (Claude).
#
#   These functions only draw; they never run the model. Each figure
#   file (Plot_deterministic.R, Plot_E1.R ... Plot_E4.R) builds on
#   them and ends with one-line nicknames that draw a figure with no
#   input, e.g. deterministic_fold_plot_1(). README.Rmd calls those.
#
#   To save a figure, open a file first and close it after:
#     png("my_figure.png", width = 1400, height = 600, res = 120)
#     deterministic_fold_plot_1()
#     dev.off()
#
#   Style: white page, light grid, quiet grey axes on the left and
#   bottom, bold titles on the left, no boxes around legends. The
#   seven boards keep one colour and one symbol each everywhere
#   (Board_col, Board_pch), so a board reads the same in every
#   figure.
#
#   Col  Board_col  Board_pch  Board_label  Pal
#   Ramp()  Ramp_col()  Count_col()  Text_on()  Panel()  Log_labels()
#   Board_points()  Map_panel()  Colour_bar()  Legend_panel()
#   Panel_title()  Plot_board()  Save_frames()  Animate_board()
# ======================================================================

Col = list(
  ink   = "#0b0b0b",               # main text, mean lines
  ink2  = "#52514e",               # secondary text, references
  muted = "#898781",               # axis numbers, collapsed runs
  grid  = "#e1e0d9",               # gridlines
  faint = "#c3c2b7",               # quiet lines and ticks
  empty = "#f0efec",               # a cell with no hares
  wall  = "#52514e",               # a village or water square
  band  = "#f9dce6",               # settled band (E1)
  cat   = c("#2a78d6", "#eb6834",  # series 1 to 4, fixed order:
            "#1baf7a", "#eda100"), #   blue, orange, aqua, yellow
  ramp  = c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5",
            "#256abf", "#184f95", "#0d366b"),  # hares: few -> many
  grass = c("#f2f6ea", "#d4e8c4", "#a9d394", "#78b765",
            "#4b9442", "#2b7130", "#124d1c")   # grass: poor -> good
)

# one colour, one symbol and one label per board, in board order
Board_col = c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4",
              "#008300", "#6250d6")
Board_pch = c(21, 22, 24, 23, 25, 3, 4)  # circle, square, triangles,
                                         #   diamond, plus, cross
Board_label = c("Real", "Shuffled", "Checker 1\u00d71",
                "Blocks 2\u00d72", "Blocks 4\u00d74",
                "10\u00d710 random", "10\u00d710 gradient")

# colour ramps for maps (low -> high)
Pal = list(
  q       = c("#f2f6ea", "#a9d394", "#4b9442", "#124d1c"),
  n       = c("#ffffff", "#f6d9b8", "#e08a3c", "#8a3b0c"),
  pattern = c("#ffffff", "#c9b8f0", "#6250d6", "#2a1a7a"),
  occ     = c("#f0efec", "#9ec5f4", "#3987e5", "#184f95", "#0d366b"),
  empty   = c("#123d7a", "#5b8fd0", "#e9f1f8", "#f6d9b8", "#e08a3c",
              "#8a3b0c")
)

# ----------------------------------------------------------------------
# Ramp: n colours along a ramp
#
#   pal   ramp colours, low -> high
#   n     number of colours
#
#   returns  n colours
# ----------------------------------------------------------------------
Ramp = function(pal, n = 100) {
  return(colorRampPalette(pal)(n))
}

# ----------------------------------------------------------------------
# Ramp_col: colours for values on a ramp
#
#   x      values
#   zlim   value range mapped onto the ramp, c(low, high)
#   pal    ramp colours, low -> high
#
#   returns  one colour per value (values outside zlim are clipped)
# ----------------------------------------------------------------------
Ramp_col = function(x, zlim, pal) {
  cols = Ramp(pal, 100)
  z = (x - zlim[1]) / max(zlim[2] - zlim[1], 1e-12)
  k = pmin(100, pmax(1, round(1 + 99 * z)))
  out = cols[k]
  return(out)
}

# ----------------------------------------------------------------------
# Count_col: colour of a cell holding x hares
#
#   grey when empty, then light to dark blue up to K hares (more than
#   K gets the darkest blue)
#
#   x   hare counts
#   K   the cap
#
#   returns  one colour per count
# ----------------------------------------------------------------------
Count_col = function(x, K) {
  pal = colorRampPalette(Col$ramp)(K)
  out = pal[pmin(pmax(x, 1), K)]
  out[x == 0] = Col$empty
  return(out)
}

# ----------------------------------------------------------------------
# Text_on: white text on dark fills, ink on light fills
#
#   fill   fill colours
#
#   returns  one text colour per fill
# ----------------------------------------------------------------------
Text_on = function(fill) {
  v = col2rgb(fill) / 255
  lum = 0.2126 * v[1, ] + 0.7152 * v[2, ] + 0.0722 * v[3, ]
  out = ifelse(lum < 0.45, "white", Col$ink)
  return(out)
}

# ----------------------------------------------------------------------
# Log_labels: tick labels 10^k for a log axis
#
#   at   tick positions (powers of ten)
#
#   returns  expressions for axis(labels = ...)
# ----------------------------------------------------------------------
Log_labels = function(at) {
  k = round(log10(at))
  return(parse(text = paste0("10^", k)))
}

# ----------------------------------------------------------------------
# Panel: an empty panel with light gridlines and quiet axes
#
#   xlim, ylim   axis ranges (xlim may run backwards, e.g. falling
#                rain)
#   main         panel title (bold, left)
#   xlab, ylab   axis titles (text or expressions)
#   log          "", "x", "y" or "xy"
#   x_at, y_at   tick positions (default: pretty ones)
#   x_lab, y_lab tick labels (default: the numbers; FALSE = none)
#   grid_x, grid_y   draw gridlines at the ticks
#   xaxs, yaxs   "r" adds 4% room at each end, "i" none
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Panel = function(xlim, ylim, main = "", xlab = "", ylab = "", log = "",
                 x_at = NULL, y_at = NULL, x_lab = TRUE, y_lab = TRUE,
                 grid_x = TRUE, grid_y = TRUE, xaxs = "r", yaxs = "r") {
  plot.new()
  plot.window(xlim = xlim, ylim = ylim, log = log, xaxs = xaxs,
              yaxs = yaxs)
  usr = par("usr")
  if (is.null(x_at)) x_at = axTicks(1)
  if (is.null(y_at)) y_at = axTicks(2)
  if (grepl("x", log)) {
    in_x = log10(x_at) >= min(usr[1:2]) - 1e-9 &
      log10(x_at) <= max(usr[1:2]) + 1e-9
  } else {
    in_x = x_at >= min(usr[1:2]) - 1e-9 & x_at <= max(usr[1:2]) + 1e-9
  }
  if (grepl("y", log)) {
    in_y = log10(y_at) >= usr[3] - 1e-9 & log10(y_at) <= usr[4] + 1e-9
  } else {
    in_y = y_at >= usr[3] - 1e-9 & y_at <= usr[4] + 1e-9
  }
  if (!isTRUE(x_lab) && !isFALSE(x_lab)) x_lab = x_lab[in_x]
  if (!isTRUE(y_lab) && !isFALSE(y_lab)) y_lab = y_lab[in_y]
  x_at = x_at[in_x]
  y_at = y_at[in_y]
  if (grid_x) abline(v = x_at, col = Col$grid, lwd = 0.8)
  if (grid_y) abline(h = y_at, col = Col$grid, lwd = 0.8)
  box(bty = "l", col = Col$muted, lwd = 0.9)
  axis(1, at = x_at, labels = x_lab, col = NA, col.ticks = Col$muted,
       col.axis = Col$muted, tcl = -0.3, mgp = c(3, 0.45, 0),
       cex.axis = 0.85)
  axis(2, at = y_at, labels = y_lab, col = NA, col.ticks = Col$muted,
       col.axis = Col$muted, tcl = -0.3, mgp = c(3, 0.6, 0), las = 1,
       cex.axis = 0.85)
  Panel_title(main)
  title(xlab = xlab, col.lab = Col$ink2, line = 2.0, cex.lab = 0.95)
  title(ylab = ylab, col.lab = Col$ink2, line = 2.9, cex.lab = 0.95)
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Panel_title: a bold title on the left above the panel
#
#   main   title (text or expression); "" = none
#   line   distance above the panel, in text lines
#   cex    size
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Panel_title = function(main, line = 0.6, cex = 1.0) {
  if (is.language(main) || is.expression(main) ||
      (is.character(main) && nzchar(main[1]))) {
    title(main = main, adj = 0, font.main = 2, col.main = Col$ink,
          line = line, cex.main = cex)
  }
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Board_points: the symbols of board b at (x, y)
#
#   x, y     positions
#   b        board number (1-7)
#   cex      symbol size
#   hollow   TRUE = white inside, board colour outside
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Board_points = function(x, y, b, cex = 1.1, hollow = FALSE) {
  pch = Board_pch[b]
  if (pch >= 21) {
    if (hollow) {
      points(x, y, pch = pch, bg = "white", col = Board_col[b],
             cex = cex, lwd = 1.6)
    } else {
      points(x, y, pch = pch, bg = Board_col[b], col = "white",
             cex = cex, lwd = 0.8)
    }
  } else {
    points(x, y, pch = pch, col = Board_col[b], cex = cex * 0.95,
           lwd = if (hollow) 1.6 else 2.4)
  }
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Map_panel: a grid of values as coloured squares, row 1 at the top
#
#   z        n_row x n_col matrix (NA = blank)
#   zlim     value range of the colours (values outside are clipped)
#   pal      ramp colours, low -> high
#   main     title above the map
#   frame    draw a thin grey frame
#   cex_main title size
#
#   returns  invisible(NULL); the map spans x, y = 0.5 .. n + 0.5,
#            so contour() and rect() can draw on it afterwards
# ----------------------------------------------------------------------
Map_panel = function(z, zlim, pal, main = "", frame = FALSE,
                     cex_main = 0.95) {
  nr = nrow(z)
  nc = ncol(z)
  zz = pmin(pmax(z, zlim[1]), zlim[2])
  image(seq_len(nc), seq_len(nr), t(zz[nr:1, , drop = FALSE]),
        zlim = zlim, col = Ramp(pal, 200), axes = FALSE, xlab = "",
        ylab = "", asp = 1, useRaster = TRUE)
  if (frame) rect(0.5, 0.5, nc + 0.5, nr + 0.5, border = Col$muted,
                  lwd = 0.8)
  Panel_title(main, line = 0.4, cex = cex_main)
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Colour_bar: a colour key drawn in its own panel
#
#   pal     ramp colours, low -> high
#   zlim    value range
#   label   key title
#   at      tick positions (default: pretty ones)
#   side    4 = vertical bar, ticks on the right; 1 = horizontal bar,
#           ticks below
#   width   bar width as a share of the panel
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Colour_bar = function(pal, zlim, label = "", at = NULL, side = 4,
                      width = 0.35) {
  if (is.null(at)) at = pretty(zlim, 5)
  at = at[at >= zlim[1] - 1e-9 & at <= zlim[2] + 1e-9]
  cols = Ramp(pal, 200)
  v = seq(zlim[1], zlim[2], length.out = 201)
  plot.new()
  if (side == 4) {
    plot.window(xlim = c(0, 1), ylim = zlim, xaxs = "i", yaxs = "i")
    rect(0, v[-201], width, v[-1], col = cols, border = NA)
    rect(0, zlim[1], width, zlim[2], border = Col$muted, lwd = 0.8)
    axis(4, at = at, pos = width, las = 1, col = NA,
         col.ticks = Col$muted, col.axis = Col$muted, tcl = -0.3,
         cex.axis = 0.85, mgp = c(3, 0.4, 0))
    mtext(label, side = 4, line = 0.4, col = Col$ink2, cex = 0.85)
  } else {
    plot.window(xlim = zlim, ylim = c(0, 1), xaxs = "i", yaxs = "i")
    rect(v[-201], 1 - width, v[-1], 1, col = cols, border = NA)
    rect(zlim[1], 1 - width, zlim[2], 1, border = Col$muted, lwd = 0.8)
    axis(1, at = at, pos = 1 - width, col = NA, col.ticks = Col$muted,
         col.axis = Col$muted, tcl = -0.3, cex.axis = 0.85,
         mgp = c(3, 0.4, 0))
    text(mean(zlim), 1 - width - 0.42, label, col = Col$ink2,
         cex = 0.9)
  }
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Legend_panel: an empty panel holding a legend (and a note)
#
#   ...    arguments for legend(), e.g. legend = , col = , lwd =
#   note   text under the legend ("" = none)
#   where  legend position inside the panel
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Legend_panel = function(..., note = "", where = "left") {
  plot.new()
  plot.window(c(0, 1), c(0, 1))
  legend(where, ..., bty = "n", text.col = Col$ink, cex = 0.95,
         inset = 0.02, seg.len = 2.2, y.intersp = 1.15)
  if (nzchar(note)) {
    text(0.02, 0.08, note, adj = c(0, 0), col = Col$ink2, cex = 0.85)
  }
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Plot_board: the map of one run at recording number k, one frame
#
#   Each square is a cell, coloured by its number of hares; a strip
#   below shows total hares, with a dot at the current time. Needs a
#   run made with Run_sim(..., record = TRUE).
#
#   res        one run from Run_sim()
#   k          recording number, 1 to length(res$t)
#   numbers    write the count in each cell (default: small maps)
#   show_q     write each cell's q in its corner (default: numbers)
#   timeline   draw the strip of total hares
#
#   returns  invisible(NULL)
# ----------------------------------------------------------------------
Plot_board = function(res, k, numbers = NULL, show_q = NULL,
                      timeline = TRUE) {
  if (res$type != "run" || is.null(res$n)) {
    stop("Plot_board needs one run made with record = TRUE")
  }
  map = res$map
  K = res$p$K
  if (is.null(numbers)) numbers = map$C <= 60
  if (is.null(show_q)) show_q = numbers
  op = par(no.readonly = TRUE)
  on.exit(par(op))
  if (timeline) {
    layout(matrix(c(1, 2, 3, 3), 2, 2, byrow = TRUE),
           widths = c(6, 1), heights = c(3, 1.4))
  } else {
    layout(matrix(1:2, 1, 2), widths = c(6, 1))
  }

  # 1. the board: row 1 of the grid at the top
  par(mar = c(0.5, 1, 3, 0.5))
  n_row = nrow(map$q_grid)
  n_col = ncol(map$q_grid)
  plot.new()
  plot.window(xlim = c(0, n_col), ylim = c(0, n_row), asp = 1)
  counts = res$n[k, ]
  fills = Count_col(counts, K)
  x0 = map$col - 1
  y0 = n_row - map$row
  lw = 30 / max(n_row, n_col)
  rect(x0, y0, x0 + 1, y0 + 1, col = fills, border = "white", lwd = lw)
  walls = which(is.na(map$q_grid), arr.ind = TRUE)
  if (nrow(walls) > 0) {
    rect(walls[, 2] - 1, n_row - walls[, 1], walls[, 2],
         n_row - walls[, 1] + 1, col = Col$wall, border = "white",
         lwd = lw)
  }
  size = min(2.5, 18 / max(n_row, n_col))
  if (numbers) {
    y_mid = if (show_q) 0.42 else 0.5    # room for the q label
    text(x0 + 0.5, y0 + y_mid, counts, col = Text_on(fills),
         cex = size)
  }
  if (show_q) {                          # small grass label
    text(x0 + 0.06, y0 + 0.94, sprintf("q %.2f", map$q),
         adj = c(0, 1), col = Text_on(fills), cex = 0.4 * size)
  }
  title(main = sprintf("Year %.0f     rain %.2f     %d hares",
                       res$t[k], res$rho_t[k], res$N[k]),
        adj = 0, font.main = 1, col.main = Col$ink)

  # 2. the colour key
  par(mar = c(0.5, 0.5, 3, 3.5))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, K + 1))
  lev = 0:K
  rect(0, lev, 0.7, lev + 1, col = Count_col(lev, K), border = NA)
  ticks = pretty(c(0, K))
  ticks = ticks[ticks <= K]
  labs = as.character(ticks)
  labs[ticks == K] = paste0(K, "+")
  axis(4, at = ticks + 0.5, labels = labs, las = 1, col = NA,
       col.ticks = Col$faint, col.axis = Col$muted)
  mtext("hares", side = 3, line = 0.5, col = Col$ink2, cex = 0.8)

  # 3. the timeline: total hares, with a dot at the current time
  if (timeline) {
    par(mar = c(4.2, 4.5, 0.5, 1))
    Panel(range(res$t), c(0, max(1, res$N) * 1.1), xlab = "year",
          ylab = "all hares")
    lines(res$t, res$N, type = "s", col = Col$grid, lwd = 1.5)
    lines(res$t[1:k], res$N[1:k], type = "s", col = Col$ink2,
          lwd = 1.5)
    points(res$t[k], res$N[k], pch = 21, bg = Col$cat[1],
           col = "white", cex = 1.8, lwd = 2)
  }
  return(invisible(NULL))
}

# ----------------------------------------------------------------------
# Save_frames: write one PNG per frame, folder/frame_0001.png, ...
#
#   res             one run from Run_sim() with record = TRUE
#   frames          recording numbers to draw
#   folder          output folder (made if missing)
#   width, height   frame size in pixels
#
#   returns  invisible(the file names)
# ----------------------------------------------------------------------
Save_frames = function(res, frames, folder, width, height) {
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  png(file.path(folder, "frame_%04d.png"), width = width,
      height = height)
  for (k in frames) Plot_board(res, k)
  dev.off()
  files = sort(list.files(folder, pattern = "^frame_[0-9]+[.]png$",
                          full.names = TRUE))
  return(invisible(files))
}

# ----------------------------------------------------------------------
# Animate_board: play Plot_board() frame after frame
#
#   res     one run from Run_sim() with record = TRUE
#   every   use every k-th recording (monthly records: every = 12
#           gives one frame per year)
#   file    NULL      play on screen
#           "x.gif"   save a GIF (needs install.packages("gifski"))
#           "folder"  save numbered PNG frames into that folder
#   delay   seconds per frame
#   width, height   frame size in pixels
#
#   returns  invisible(the recording numbers drawn)
# ----------------------------------------------------------------------
Animate_board = function(res, every = 12, file = NULL, delay = 0.12,
                         width = 800, height = 560) {
  if (res$type != "run" || is.null(res$n)) {
    stop("Animate_board needs one run made with record = TRUE")
  }
  frames = seq(1, length(res$t), by = every)
  if (is.null(file)) {                   # on screen
    for (k in frames) {
      dev.hold()
      Plot_board(res, k)
      dev.flush()
      Sys.sleep(delay)
    }
  } else if (grepl("[.]gif$", file)) {   # a GIF
    if (!requireNamespace("gifski", quietly = TRUE)) {
      stop("Animate_board: a GIF needs install.packages(\"gifski\"); ",
           "or give a folder name to save PNG frames instead")
    }
    tmp = file.path(tempdir(), "hare_frames")
    unlink(tmp, recursive = TRUE)
    pngs = Save_frames(res, frames, tmp, width, height)
    gifski::gifski(pngs, gif_file = file, width = width,
                   height = height, delay = delay, progress = FALSE)
  } else {                               # PNG frames
    Save_frames(res, frames, file, width, height)
  }
  return(invisible(frames))
}
