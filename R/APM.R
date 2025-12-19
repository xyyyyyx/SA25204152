#' Aggregated Projection Method (APM) for Group Factor Model
#'
#' This function implements the Algorithm 1 proposed by Li et al. (2024).
#' It separates global and local factors from grouped panel data.
#'
#' @param Y_list A list of matrices, where each element is a (T x N_m) matrix for group m.
#' @param r_max An integer. The maximum number of factors to consider. Default is 8.
#' @param method A character string. The method to determine the number of global factors ("Gap" or "IC"). Default is "Gap".
#'
#' @return A list containing:
#' \item{G_hat}{The estimated global factors (T x r0).}
#' \item{F_local}{A list of estimated local factors for each group.}
#' \item{Loadings_Global}{A list of global factor loadings for each group.}
#' \item{Loadings_Local}{A list of local factor loadings for each group.}
#' \item{r0}{The estimated number of global factors.}
#' \item{eigen_values}{The eigenvalues of the aggregated projection matrix.}
#'
#' @examples
#' \dontrun{
#' # Generate simulation data
#' data <- sim_group_data(M=3, N=30, T_obs=100)
#' # Run APM
#' result <- apm(data, r_max=5)
#' print(result$r0)
#' }
#' @export
apm <- function(Y_list, r_max = 8, method = "Gap") {

  M <- length(Y_list)
  if (M == 0) stop("Y_list must containing at least one group.")
  T_obs <- nrow(Y_list[[1]])

  K_hat_list <- list()

  # Step 2: Initial PCA for each group (Algorithm 1, Line 2)
  # 这里我们简化处理，假设每个组先提取固定的 r_max 个因子作为初始空间
  for (m in 1:M) {
    # 对协方差矩阵进行特征分解
    Y_c <- scale(Y_list[[m]], center = TRUE, scale = FALSE) # 中心化

    # 提取前 r_max 个因子 (使用 SVD 是更稳健的计算方式)
    # Factor F = sqrt(T) * eigenvectors
    factors <- svd(Y_c)$u[, 1:r_max, drop = FALSE] * sqrt(T_obs)
    K_hat_list[[m]] <- factors
  }

  # Step 3: Construct Aggregated Projection Matrix (Algorithm 1, Line 3)
  # 调用我们在 src/apm_core.cpp 里写的 C++ 函数
  Psi_hat <- calc_psi_cpp(K_hat_list)

  # Step 4: Eigen-decomposition (Algorithm 1, Line 4)
  eigen_res <- eigen(Psi_hat, symmetric = TRUE)
  eigen_vals <- eigen_res$values
  eigen_vecs <- eigen_res$vectors

  # Step 5: Determine r0 (Algorithm 1, Line 5)
  r0_hat <- 0
  if (method == "Gap") {
    # Eigenvalue Gap Method
    # 寻找特征值跌落最大的位置
    diffs <- -diff(eigen_vals[1:(r_max + 1)])
    r0_hat <- which.max(diffs)
  } else {
    # 默认使用 Gap
    diffs <- -diff(eigen_vals[1:(r_max + 1)])
    r0_hat <- which.max(diffs)
  }

  # Step 6: Estimate Global Factors (Algorithm 1, Line 6)
  if (r0_hat > 0) {
    G_hat <- eigen_vecs[, 1:r0_hat, drop = FALSE] * sqrt(T_obs)
  } else {
    G_hat <- matrix(0, nrow = T_obs, ncol = 0)
  }

  # Step 7: Estimate Local Factors (Algorithm 1, Line 7)
  F_local_list <- list()
  Gamma_list <- list()
  Lambda_list <- list()

  # 构造投影矩阵 P(G)
  if (r0_hat > 0) {
    # P_G = G(G'G)^-1G'
    P_G <- G_hat %*% solve(crossprod(G_hat), t(G_hat))
  } else {
    P_G <- matrix(0, T_obs, T_obs)
  }
  I_T <- diag(T_obs)
  M_G <- I_T - P_G # 这里的 M_G 是把 Global 投影掉剩下的部分

  for (m in 1:M) {
    Y_m <- Y_list[[m]]
    Y_m_c <- scale(Y_m, center = TRUE, scale = FALSE)

    # 估计 Global Loadings: Gamma = Y'G / T
    if (r0_hat > 0) {
      Gamma_m <- t(Y_m_c) %*% G_hat / T_obs
    } else {
      Gamma_m <- matrix(0, ncol(Y_m), 0)
    }
    Gamma_list[[m]] <- Gamma_m

    # 投影掉 Global 因子: Y_hat = (I - P_G)Y
    Y_m_tilde <- M_G %*% Y_m_c

    # 对残差做 PCA 得到 Local Factors
    # 因子数简单取 max(1, r_max - r0)
    r_local <- max(1, r_max - r0_hat)
    svd_local <- svd(Y_m_tilde)
    F_m <- svd_local$u[, 1:r_local, drop = FALSE] * sqrt(T_obs)
    Lambda_m <- t(Y_m_tilde) %*% F_m / T_obs

    F_local_list[[m]] <- F_m
    Lambda_list[[m]] <- Lambda_m
  }

  # Step 8: Output (Algorithm 1, Line 8)
  return(list(
    G_hat = G_hat,
    F_local = F_local_list,
    Loadings_Global = Gamma_list,
    Loadings_Local = Lambda_list,
    r0 = r0_hat,
    eigen_values = eigen_vals[1:(r_max + 1)]
  ))
}
