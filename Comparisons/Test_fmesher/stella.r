cat("SCRIPT R PARTITO\n")
library(fmesher)


# 1. Crea la stella
x <- c(0.0, -22.45, -95.11, -36.33, -58.78, 
       0.0, 58.78, 36.33, 95.11, 22.45)
y <- c(100.0, 30.90, 30.90, -11.80, -80.90,
       -38.20, -80.90, -11.80, 30.90, 30.90)
loc <- cbind(x, y)

# 2. Segmenti come matrice
segments <- cbind(1:nrow(loc), c(2:nrow(loc), 1))

# 3. Converti in oggetto segmento
segm_obj <- fm_segm(loc = loc, idx = segments, closed=TRUE)
boundary <- fm_as_segm(segm_obj)

# 4. Triangolazione
mesh <- fmesher::fm_mesh_2d(boundary = boundary)
max_area <- 200
max_edge <- sqrt(4 * max_area / sqrt(3))   # per triangoli quasi equilateri

mesh <- fm_mesh_2d(
  boundary = boundary,
  max.edge = 24.389,             # (opzionale) controllo sulla dimensione massima degli spigoli
  min.angle = 25,              # <-- vincolo sugli angoli minimi
  delaunay = TRUE
)



# 6. Statistiche
mesh_statistics <- function(mesh) {
  triangles <- mesh$graph$tv
  loc <- mesh$loc

  n_tri <- nrow(triangles)

  min_area <- Inf
  max_area <- 0
  min_edge <- Inf
  max_edge <- 0
  min_altitude <- Inf
  max_aspect_ratio <- 0
  min_angle <- 180
  max_angle <- 0

  triangle_area <- function(A, B, C) {
    0.5 * abs((B[1] - A[1]) * (C[2] - A[2]) - (B[2] - A[2]) * (C[1] - A[1]))
  }

  edge_length <- function(P, Q) sqrt(sum((P - Q)^2))

  triangle_angles <- function(A, B, C) {
    a <- sqrt(sum((B - C)^2))
    b <- sqrt(sum((A - C)^2))
    c <- sqrt(sum((A - B)^2))
    angleA <- acos(pmin(pmax((b^2 + c^2 - a^2) / (2 * b * c), -1), 1)) * 180 / pi
    angleB <- acos(pmin(pmax((a^2 + c^2 - b^2) / (2 * a * c), -1), 1)) * 180 / pi
    angleC <- acos(pmin(pmax((a^2 + b^2 - c^2) / (2 * a * b), -1), 1)) * 180 / pi
    c(angleA, angleB, angleC)
  }

  for (i in 1:n_tri) {
    idx <- triangles[i, ]
    A <- loc[idx[1], ]
    B <- loc[idx[2], ]
    C <- loc[idx[3], ]

    area <- triangle_area(A, B, C)
    lengths <- c(edge_length(A, B), edge_length(B, C), edge_length(C, A))
    longest <- max(lengths)
    altitude <- 2 * area / longest
    aspect <- longest / altitude
    angles <- triangle_angles(A, B, C)

    min_area <- min(min_area, area)
    max_area <- max(max_area, area)
    min_edge <- min(min_edge, lengths)
    max_edge <- max(max_edge, lengths)
    min_altitude <- min(min_altitude, altitude)
    max_aspect_ratio <- max(max_aspect_ratio, aspect)
    min_angle <- min(min_angle, angles)
    max_angle <- max(max_angle, angles)
  }

  list(
    mesh_vertices = nrow(loc),
    mesh_triangles = n_tri,
    min_area = min_area,
    max_area = max_area,
    min_edge = min_edge,
    max_edge = max_edge,
    min_altitude = min_altitude,
    max_aspect_ratio = max_aspect_ratio,
    min_angle = min_angle,
    max_angle = max_angle
  )
}

stats <- mesh_statistics(mesh)
cat("\nMesh Statistics:\n")
cat("  Vertices:            ", stats$mesh_vertices, "\n")
cat("  Triangles:           ", stats$mesh_triangles, "\n")
cat("  Min area:            ", round(stats$min_area, 5), "\n")
cat("  Max area:            ", round(stats$max_area, 5), "\n")
cat("  Shortest edge:       ", round(stats$min_edge, 5), "\n")
cat("  Longest edge:        ", round(stats$max_edge, 5), "\n")
cat("  Min altitude:        ", round(stats$min_altitude, 5), "\n")
cat("  Max aspect ratio:    ", round(stats$max_aspect_ratio, 5), "\n")
cat("  Min angle (degrees): ", round(stats$min_angle, 2), "\n")
cat("  Max angle (degrees): ", round(stats$max_angle, 2), "\n")

# 5. Plot
plot(mesh, asp = 1)


