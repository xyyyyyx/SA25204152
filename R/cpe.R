#' Circular Projection Estimator (CPE) for Comparison
#'
#' Implements the iterative projection method (Chen, 2022) to serve as a benchmark.
#'
#' @param Y_list A list of matrices (T x N_m).
#' @param r0 The number of global factors to estimate.
#' @param r_max The dimension of the subspace to project onto in each group.
#' @param max_iter Maximum iterations.
#'
#' @return A list containing the estimated global factors `G_hat`.
#' @export
#' @importFrom stats rnorm
cpe <- function(Y_list, r0 = 2, r_max = 8, max_iter = 100) {
  if (r0 == 0) return(list(G_hat = matrix(0, nrow(Y_list[[1]]), 0)))

  T_obs <- nrow(Y_list[[1]])
  M <- length(Y_list)

  # 1. Prepare Subspaces K_m
  K_list <- list()
  for (m in 1:M) {
    Y_c <- scale(Y_list[[m]], center = TRUE, scale = FALSE)
    factors <- svd(Y_c)$u[, 1:r_max, drop = FALSE] * sqrt(T_obs)
    K_list[[m]] <- factors
  }

  # 2. Initial Guess for G (Random)
  G_init <- matrix(rnorm(T_obs * r0), T_obs, r0)
  G_init <- svd(G_init)$u * sqrt(T_obs) # Orthogonalize

  # 3. Call C++ Implementation
  G_hat <- fit_cpe_cpp(K_list, G_init, max_iter)

  return(list(G_hat = G_hat))
}
