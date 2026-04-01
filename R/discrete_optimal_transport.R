# Discrete Optimal Transport (R translation of the MATLAB version)
# cp     : convex polygon (n x 2)
# face   : (nf x 3) triangle connectivity
# uv     : (k x 2) sites
# sigma  : function(xy) -> scalar or vector (density over R^2)
# delta  : length-k target masses, sum(delta) = ∫_cp sigma
# h      : optional initial weights (length k); if NULL, zeros
# option : list with fields $eps (default 1e-5) and $max_iter (default 20)
#
# Returns: pd (power-diagram list; includes updated h inside pd$h if your implementation sets it)
discrete_optimal_transport <- function(cp, face, uv, sigma, delta, h = NULL, option = list()) {
  # validate
  stopifnot(is.matrix(cp), ncol(cp) == 2)
  stopifnot(is.matrix(uv), ncol(uv) == 2)
  stopifnot(is.function(sigma))
  stopifnot(is.numeric(delta), length(delta) == nrow(uv))
  
  np <- nrow(uv)
  if (is.null(h) || length(h) == 0) h <- rep(0, np)
  
  # options
  eps <- if (!is.null(option$eps)) option$eps else 1e-5
  max_iter <- if (!is.null(option$max_iter)) option$max_iter else 20
  
  # initial PD
  pd <- power_diagram(face, uv, h)
  
  k <- 1L
  t0 <- proc.time()[["elapsed"]]
  
  while (k < max_iter) {
    # Gradient (cell masses under sigma, clipped to cp)
    G <- calculate_gradient(cp, pd, sigma)
    G <- G / sum(G) * sum(delta)  # normalize to match total mass of delta
    D <- G - delta
    
    # Hessian
    H <- calculate_hessian(cp, pd, sigma)
    
    # stabilize (fix gauge freedom)
    H[1, 1] <- H[1, 1] + 1
    
    # Solve H * dh = D  (MATLAB: dh = H\D)
    dh <- tryCatch(Matrix::solve(H, D), error = function(e) NA_real_)
    
    if (!all(is.finite(dh))) {
      stop(paste0(
        "ERROR: |dh| went infinite/NA, likely because the convex hull step failed ",
        "(a cell disappeared). Underlying cause is often poor mesh/measure."
      ))
    }
    
    # remove mean shift (twice, matching MATLAB code)
    dh <- dh - mean(dh)
    dh <- dh - mean(dh)
    
    # cat(sprintf("#%02d: max|dh| = %.10f\n", k, max(abs(dh))))
    
    if (max(abs(dh)) < eps) break
    
    # Backtracking & update happen inside power_diagram(face, uv, h, dh)
    pd <- power_diagram(face, uv, h, dh)
    
    # If your power_diagram returns updated weights in pd$h, adopt them for the next step:
    if (!is.null(pd$h)) h <- pd$h
    
    k <- k + 1L
  }
  
  # cat(sprintf("Elapsed: %.3fs\n", proc.time()[['elapsed']] - t0))
  pd
}

