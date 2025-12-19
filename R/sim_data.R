#' Generate Simulation Data for Group Factor Model
#'
#' Generates data based on the settings in Li et al. (2024), Slide 35.
#'
#' @param M Number of groups.
#' @param N Number of individuals per group (can be a vector or scalar).
#' @param T_obs Number of time observations.
#' @param r0 Number of global factors.
#' @param r_local Number of local factors per group (scalar).
#' @param phi_G AR(1) coefficient for global factors.
#' @param phi_F AR(1) coefficient for local factors.
#' @param phi_e AR(1) coefficient for idiosyncratic errors.
#' @return A list containing `Y_list` (data) and true factors/loadings.
#' @export
#' @importFrom stats rnorm
sim_group_data <- function(M = 3, N = 50, T_obs = 100, r0 = 2, r_local = 2,
                           phi_G = 0.5, phi_F = 0.5, phi_e = 0.5) {

  if (length(N) == 1) N <- rep(N, M)

  # 1. Generate Global Factors G_t (AR(1))
  G <- matrix(0, T_obs, r0)
  if (r0 > 0) {
    G[1, ] <- rnorm(r0)
    for (t in 2:T_obs) {
      G[t, ] <- phi_G * G[t-1, ] + rnorm(r0)
    }
    # Standardize to satisfy assumption G'G/T = I asymptotically
    G <- scale(G, center = TRUE, scale = FALSE)
  }

  Y_list <- list()
  F_list <- list()
  Gamma_list <- list()
  Lambda_list <- list()

  for (m in 1:M) {
    Nm <- N[m]

    # 2. Generate Local Factors F_m,t (Case 1: Independent, AR(1))
    F_m <- matrix(0, T_obs, r_local)
    if (r_local > 0) {
      F_m[1, ] <- rnorm(r_local)
      for (t in 2:T_obs) {
        F_m[t, ] <- phi_F * F_m[t-1, ] + rnorm(r_local)
      }
      F_m <- scale(F_m, center = TRUE, scale = FALSE)
    }
    F_list[[m]] <- F_m

    # 3. Generate Loadings (Normal(0,1))
    if (r0 > 0) Gamma_m <- matrix(rnorm(Nm * r0), Nm, r0) else Gamma_m <- matrix(0, Nm, 0)
    if (r_local > 0) Lambda_m <- matrix(rnorm(Nm * r_local), Nm, r_local) else Lambda_m <- matrix(0, Nm, 0)

    Gamma_list[[m]] <- Gamma_m
    Lambda_list[[m]] <- Lambda_m

    # 4. Generate Errors e_m,i,t (Weakly correlated)
    E_m <- matrix(0, T_obs, Nm)
    V <- matrix(rnorm(T_obs * Nm), T_obs, Nm)
    E_m[1, ] <- V[1, ]
    for (t in 2:T_obs) {
      E_m[t, ] <- phi_e * E_m[t-1, ] + V[t, ]
    }

    # 5. Combine: Y = G * Gamma' + F * Lambda' + E
    Signal <- matrix(0, T_obs, Nm)
    if (r0 > 0) Signal <- Signal + G %*% t(Gamma_m)
    if (r_local > 0) Signal <- Signal + F_m %*% t(Lambda_m)

    Y_m <- Signal + E_m
    Y_list[[m]] <- Y_m
  }

  list(Y_list = Y_list, G = G, F_list = F_list,
       Gamma_list = Gamma_list, Lambda_list = Lambda_list)
}
