# ======================================================================
# Map.R   the grid: grass quality, walls, neighbours, the seven boards
# ======================================================================
#
#   Outcome page, section 1.2 and step 1 of the build chart (1.1).
#
#   A grid is a matrix of grass quality q between 0 and 1. NA marks a
#   wall (village core, water) that no hare enters. Row 1 is the
#   north edge, column 1 the west edge. Usable squares become cells
#   1..C, numbered column by column, skipping walls.
#
#   crossing chance from cell i into neighbour j:
#     a_ij = q_j^alpha / (q_i^alpha + q_j^alpha)
#     0 toward an edge or a wall; 0.5 when q_i = q_j = 0
#   Edges are closed: a try across the map's edge fails and the hare
#   stays.
#
#   data/mayureshwar_q_grid.csv comes from ESA WorldCover 2021
#   (10 m pixels), about 10 x 10 km around Mayureshwar Wildlife
#   Sanctuary (74.32 to 74.42 E, 18.29 to 18.38 N), averaged into
#   500 m cells with weights grassland 1.0, shrubland 0.6, cropland
#   0.5, trees 0.3, bare 0.2. A cell more than half built-up or water
#   is a wall: only the village core (row 10, column 12).
#
#   The seven boards (outcome page, section 2.2): the same 400 grass
#   values in seven arrangements. The Real board fills the village
#   square with the mean of the 399 real cells (0.662), so every board
#   has 400 cells and no wall. data/boards.csv holds the boards every
#   result used (made with seed 20261010 in the first, Python version
#   of the code); Make_boards() follows the same recipe with R's random
#   numbers, so its random layouts differ in detail.
#
#   Read_grid()     Build_map()     Crossing_chance()  Hop_table()
#   Cell_id()       Place_hares()   Shuffle_grid()     Block_grid()
#   Gradient_grid() Make_boards()   Read_boards()      Board_map()
#   Cells_to_grid() Grid_to_cells()
# ======================================================================

Board_names = c("Real", "Shuffled", "Checker 1x1", "Blocks 2x2",
                "Blocks 4x4", "10x10 random", "10x10 gradient")

# ----------------------------------------------------------------------
# Read_grid: read a grid of grass quality from a CSV file
#
#   file   CSV path: no header, one line per grid row (first line =
#          north edge), NA = wall
#
#   returns  q matrix, n_row x n_col (NA = wall)
# ----------------------------------------------------------------------
Read_grid = function(file) {
  g = as.matrix(read.csv(file, header = FALSE))
  dimnames(g) = NULL
  if (!is.numeric(g)) {
    stop("Read_grid: ", file, " holds text; use numbers and NA")
  }
  return(g)
}

# ----------------------------------------------------------------------
# Build_map: number the usable squares and find their neighbours
#
#   grid   a q matrix (NA = wall) or the path of a CSV file
#
#   returns  map, a list:
#     C       number of usable cells
#     q       grass quality of each cell, length C
#     row     grid row of each cell, length C
#     col     grid column of each cell, length C
#     id      n_row x n_col: cell number of each square (0 = wall)
#     nbr     C x 4: neighbour up, down, left, right (0 = edge or
#             wall)
#     q_grid  the grid, kept for drawing
#   Run_sim() adds a, the C x 4 crossing chances, from p$alpha.
# ----------------------------------------------------------------------
Build_map = function(grid) {
  if (is.character(grid)) grid = Read_grid(grid)
  q_grid = as.matrix(grid)
  use = !is.na(q_grid)
  if (!any(use)) stop("Build_map: the grid has no usable square")
  if (!is.numeric(q_grid)) stop("Build_map: the grid must hold numbers")
  if (any(q_grid[use] < 0 | q_grid[use] > 1)) {
    stop("Build_map: grass quality q must lie between 0 and 1")
  }
  n_row = nrow(q_grid)
  n_col = ncol(q_grid)
  C = sum(use)
  id = matrix(0L, n_row, n_col)          # n_row x n_col: cell numbers,
  id[use] = seq_len(C)                   #   column by column
  pos = which(use, arr.ind = TRUE)       # C x 2: row, col of each cell
  q = q_grid[use]                        # length C

  d_row = c(-1, 1, 0, 0)                 # up, down, left, right
  d_col = c(0, 0, -1, 1)
  nbr = matrix(0L, C, 4)                 # C x 4: neighbour up, down,
  for (i in seq_len(C)) {                #   left, right (0 = none)
    for (d in 1:4) {
      rr = pos[i, 1] + d_row[d]
      cc = pos[i, 2] + d_col[d]
      if (rr < 1 || rr > n_row || cc < 1 || cc > n_col) next
      nbr[i, d] = id[rr, cc]             # 0 if that square is a wall
    }
  }
  colnames(nbr) = c("up", "down", "left", "right")

  map = list(C = C, q = q, row = unname(pos[, 1]),
             col = unname(pos[, 2]), id = id, nbr = nbr,
             q_grid = q_grid)
  return(map)
}

# ----------------------------------------------------------------------
# Crossing_chance: chance that a hare in cell i, peeking into
# neighbour j, crosses
#
#   a_ij = q_j^alpha / (q_i^alpha + q_j^alpha)
#
#   q_i, q_j   grass quality of the two cells, 0 to 1
#   alpha      pickiness, 0 or more (0 = a coin flip)
#
#   Example: from q 0.80, a try toward q 0.95 crosses with chance
#   0.59, toward q 0.30 with chance 0.12 (alpha = 2).
#
#   returns  a_ij, one number from 0 to 1 (0.5 if both are bare)
# ----------------------------------------------------------------------
Crossing_chance = function(q_i, q_j, alpha) {
  bottom = q_i^alpha + q_j^alpha
  if (bottom == 0) return(0.5)
  a = q_j^alpha / bottom
  return(a)
}

# ----------------------------------------------------------------------
# Hop_table: crossing chances of every cell toward its 4 sides
#
#   a[i, d] = a_ij for the neighbour j = nbr[i, d]; 0 if nbr[i, d] = 0
#
#   map     map from Build_map()
#   alpha   pickiness, 0 or more
#
#   returns  a, C x 4 matrix: up, down, left, right
# ----------------------------------------------------------------------
Hop_table = function(map, alpha) {
  a = matrix(0, map$C, 4)                # C x 4: crossing chances
  for (i in seq_len(map$C)) {
    for (d in 1:4) {
      j = map$nbr[i, d]
      if (j > 0) a[i, d] = Crossing_chance(map$q[i], map$q[j], alpha)
    }
  }
  colnames(a) = c("up", "down", "left", "right")
  return(a)
}

# ----------------------------------------------------------------------
# Cell_id: cell numbers of grid squares given by row and column
#
#   map          map from Build_map()
#   rows, cols   grid rows and columns, same length (row 1 = north)
#
#   returns  cell numbers, one per square; stops on a wall or on a
#            square off the grid
# ----------------------------------------------------------------------
Cell_id = function(map, rows, cols) {
  if (length(rows) != length(cols)) {
    stop("Cell_id: rows and cols need the same length")
  }
  off = rows < 1 | rows > nrow(map$id) | cols < 1 |
    cols > ncol(map$id) | rows != round(rows) | cols != round(cols)
  if (any(off)) stop("Cell_id: a square lies off the grid")
  ids = map$id[cbind(rows, cols)]
  if (any(ids == 0)) {
    k = which(ids == 0)[1]
    stop("Cell_id: row ", rows[k], ", column ", cols[k], " is a wall")
  }
  return(ids)
}

# ----------------------------------------------------------------------
# Place_hares: a starting board with hares on chosen grid squares
#
#   map          map from Build_map()
#   rows, cols   grid rows and columns of the squares
#   counts       hares on each square: one number, or one per square
#
#   Example (E4): a 2 x 2 colony of 14 hares per cell with its top
#   left corner at row 3, column 9:
#     Place_hares(map, c(3, 4, 3, 4), c(9, 9, 10, 10), 14)
#
#   returns  n0, length C: hares in each cell (0 elsewhere)
# ----------------------------------------------------------------------
Place_hares = function(map, rows, cols, counts) {
  ids = Cell_id(map, rows, cols)
  if (any(duplicated(ids))) stop("Place_hares: a square appears twice")
  if (!(length(counts) %in% c(1, length(ids)))) {
    stop("Place_hares: give one count, or one per square")
  }
  if (any(counts < 0 | counts != round(counts))) {
    stop("Place_hares: counts must be whole numbers, 0 or more")
  }
  n0 = rep(0, map$C)                     # length C
  n0[ids] = counts
  return(n0)
}

# ----------------------------------------------------------------------
# Shuffle_grid: the same q values, moved to random usable squares
#
#   grid   a q matrix (NA = wall); walls stay in place
#   seed   random seed, so the shuffle is repeatable (it resets R's
#          random numbers); NULL = no reset
#
#   returns  q matrix, same shape and walls as grid
# ----------------------------------------------------------------------
Shuffle_grid = function(grid, seed = NULL) {
  g = as.matrix(grid)
  use = !is.na(g)
  if (!is.null(seed)) set.seed(seed)
  g[use] = g[use][sample.int(sum(use))]
  return(g)
}

# ----------------------------------------------------------------------
# Block_grid: w x w blocks of good and poor grass, like a chessboard,
# made from the grid's own values
#
#   A block is good when (ceiling(row / w) + ceiling(col / w)) is
#   even. Good blocks take the best values of the grid (as many as
#   there are good cells), poor blocks the rest, each in random order.
#   So every board keeps the same 400 values: w = 1 is the Checker
#   1x1 board, w = 2 and 4 the Blocks boards, w = 10 the 10x10 random
#   board (outcome page, section 2.2).
#
#   grid   a q matrix with no walls
#   w      block width in cells
#   seed   random seed (NULL = no reset)
#
#   returns  q matrix, same shape as grid
# ----------------------------------------------------------------------
Block_grid = function(grid, w, seed = NULL) {
  g = as.matrix(grid)
  if (anyNA(g)) stop("Block_grid: fill the walls first")
  if (!is.null(seed)) set.seed(seed)
  good = (ceiling(row(g) / w) + ceiling(col(g) / w)) %% 2 == 0
  vals = sort(g, decreasing = TRUE)
  n_good = sum(good)
  out = g
  out[good] = sample(vals[seq_len(n_good)])
  out[!good] = sample(vals[-seq_len(n_good)])
  return(out)
}

# ----------------------------------------------------------------------
# Gradient_grid: the 10x10 gradient board
#
#   The same chessboard of 10 x 10 quadrants as Block_grid(grid, 10):
#   good quadrants top left and bottom right. Inside each quadrant
#   the values fall (good) or rise (poor) with distance from its
#   centre: the best value sits at the centre of a good quadrant, the
#   worst at the centre of a poor one. The two quadrants of a kind
#   take the sorted values in turn; equal distances break at random.
#
#   grid   a 20 x 20 q matrix with no walls
#   seed   random seed (NULL = no reset)
#
#   returns  20 x 20 q matrix
# ----------------------------------------------------------------------
Gradient_grid = function(grid, seed = NULL) {
  g = as.matrix(grid)
  if (!all(dim(g) == 20) || anyNA(g)) {
    stop("Gradient_grid: needs a 20 x 20 grid with no walls")
  }
  if (!is.null(seed)) set.seed(seed)
  vals = sort(g, decreasing = TRUE)
  dist = sqrt((row(matrix(0, 10, 10)) - 5.5)^2 +
                (col(matrix(0, 10, 10)) - 5.5)^2)   # 10 x 10
  corner = list(good = list(c(0, 0), c(10, 10)),   # top-left offsets
                poor = list(c(0, 10), c(10, 0)))
  out = matrix(0, 20, 20)
  for (kind in c("good", "poor")) {
    v = if (kind == "good") vals[1:200] else rev(vals[201:400])
    for (k in 1:2) {
      mine = v[seq(k, 200, by = 2)]      # dealt in turn
      ord = order(dist, runif(100))      # by distance, random ties
      blk = matrix(0, 10, 10)
      blk[ord] = mine
      rr = corner[[kind]][[k]][1] + 1:10
      cc = corner[[kind]][[k]][2] + 1:10
      out[rr, cc] = blk
    }
  }
  return(out)
}

# ----------------------------------------------------------------------
# Make_boards: the seven boards from one grid, by the recipe above
#
#   grid   a 20 x 20 q matrix; a wall takes the mean of the others
#   seed   random seed for the random layouts
#
#   returns  named list of seven 20 x 20 q matrices (Board_names)
# ----------------------------------------------------------------------
Make_boards = function(grid = Read_grid("data/mayureshwar_q_grid.csv"),
                       seed = 20261010) {
  g = as.matrix(grid)
  g[is.na(g)] = mean(g, na.rm = TRUE)    # village cell: the mean, 0.662
  set.seed(seed)
  b = list(g, Shuffle_grid(g), Block_grid(g, 1), Block_grid(g, 2),
           Block_grid(g, 4), Block_grid(g, 10), Gradient_grid(g))
  names(b) = Board_names
  return(b)
}

# ----------------------------------------------------------------------
# Read_boards: the seven boards every result used
#
#   file   data/boards.csv: board (1-7), name, row, col, q
#
#   returns  named list of seven 20 x 20 q matrices
# ----------------------------------------------------------------------
Read_boards = function(file = "data/boards.csv") {
  d = read.csv(file)
  b = lapply(1:7, function(k) {
    s = d[d$board == k, ]
    g = matrix(NA_real_, max(s$row), max(s$col))
    g[cbind(s$row, s$col)] = s$q
    return(g)
  })
  names(b) = Board_names
  return(b)
}

# ----------------------------------------------------------------------
# Board_map: the map of board k (1 = Real, ..., 7 = 10x10 gradient)
#
#   k       board number
#   boards  boards from Read_boards()
#
#   returns  map from Build_map()
# ----------------------------------------------------------------------
Board_map = function(k, boards = Read_boards()) {
  map = Build_map(boards[[k]])
  return(map)
}

# ----------------------------------------------------------------------
# Cells_to_grid: per-cell values laid out on the grid (NA on walls)
#
#   x     one value per cell, length C
#   map   map from Build_map()
#
#   returns  n_row x n_col matrix
# ----------------------------------------------------------------------
Cells_to_grid = function(x, map) {
  g = matrix(NA_real_, nrow(map$id), ncol(map$id))
  g[cbind(map$row, map$col)] = x
  return(g)
}

# ----------------------------------------------------------------------
# Grid_to_cells: grid values read out cell by cell
#
#   g     n_row x n_col matrix
#   map   map from Build_map()
#
#   returns  length C, in the map's cell order
# ----------------------------------------------------------------------
Grid_to_cells = function(g, map) {
  x = g[cbind(map$row, map$col)]
  return(x)
}
