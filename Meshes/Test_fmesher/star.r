library(fmesher)

analyze_mesh_quality <- function(mesh) {
  triangles <- mesh$graph$tv
  coords <- mesh$loc

  triangle_areas <- numeric(nrow(triangles))
  edge_lengths <- matrix(0, nrow = nrow(triangles), ncol = 3)
  aspect_ratios <- numeric(nrow(triangles))
  min_altitudes <- numeric(nrow(triangles))
  angles_deg <- matrix(0, nrow = nrow(triangles), ncol = 3)

  for (i in seq_len(nrow(triangles))) {
    idx <- triangles[i, ]
    p <- coords[idx, ]

    a <- sqrt(sum((p[2, ] - p[3, ])^2))
    b <- sqrt(sum((p[1, ] - p[3, ])^2))
    c <- sqrt(sum((p[1, ] - p[2, ])^2))
    edge_lengths[i, ] <- c(a, b, c)

    s <- (a + b + c) / 2
    area <- sqrt(s * (s - a) * (s - b) * (s - c))
    triangle_areas[i] <- area

    h1 <- 2 * area / a
    h2 <- 2 * area / b
    h3 <- 2 * area / c
    min_altitudes[i] <- min(c(h1, h2, h3))

    aspect_ratios[i] <- max(c(a, b, c)) / min_altitudes[i]

    angleA <- acos((b^2 + c^2 - a^2) / (2 * b * c)) * 180 / pi
    angleB <- acos((a^2 + c^2 - b^2) / (2 * a * c)) * 180 / pi
    angleC <- 180 - angleA - angleB
    angles_deg[i, ] <- c(angleA, angleB, angleC)
  }

  n_vertices <- nrow(coords)
  n_triangles <- nrow(triangles)

  edge_list <- lapply(1:nrow(triangles), function(i) {
    tri <- triangles[i, ]
    list(
      sort(c(tri[1], tri[2])),
      sort(c(tri[2], tri[3])),
      sort(c(tri[3], tri[1]))
    )
  })
  all_edges <- do.call(rbind, unlist(edge_list, recursive = FALSE))
  unique_edges <- unique(all_edges)
  edge_counts <- table(apply(all_edges, 1, paste, collapse = "-"))
  n_edges <- nrow(unique_edges)
  n_hull_edges <- sum(edge_counts == 1)

  summary_metrics <- data.frame(
    Metric = c("Vertices", "Triangles", "Edges", "MinArea", "MaxArea", 
               "MinEdge", "MaxEdge", "MinAltitude", "MaxAspectRatio", 
               "MinAngle", "MaxAngle"),
    Value = c(
      n_vertices,
      n_triangles,
      n_edges,
      min(triangle_areas),
      max(triangle_areas),
      min(edge_lengths),
      max(edge_lengths),
      min(min_altitudes),
      max(aspect_ratios),
      min(angles_deg),
      max(angles_deg)
    )
  )

  cat("\nMesh structural info:\n\n")
  print(summary_metrics)

  ar_bins <- c(1.1547, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 10.0, 15.0, 25.0,
               50.0, 100.0, 300.0, 1000.0, 10000.0, 100000.0, Inf)
  ar_labels <- paste0("\"", format(ar_bins[-length(ar_bins)]), " - ", format(ar_bins[-1]), "\"")
  ar_hist <- table(cut(aspect_ratios, breaks = ar_bins, labels = ar_labels, right = FALSE))
  df_ar <- data.frame(AspectRatioRange = names(ar_hist), Count = as.integer(ar_hist))

  cat("\nAspect Ratio Histogram:\n")
  print(df_ar)

  angle_bins <- seq(0, 180, by = 10)
  angle_labels <- paste0("\"", angle_bins[-length(angle_bins)], " - ", angle_bins[-1], " degrees\"")
  angle_hist <- table(cut(as.vector(angles_deg), breaks = angle_bins, labels = angle_labels, right = FALSE))
  df_ang <- data.frame(AngleRange = names(angle_hist), Count = as.integer(angle_hist))

  cat("\nAngle Histogram:\n")
  print(df_ang)
}


