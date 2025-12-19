#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

//' Calculate Aggregated Projection Matrix (Core for APM)
//'
//' @param K_list A list of matrices (T x r) representing factor spaces for each group.
//' @return A matrix (T x T) representing the average projection matrix.
// [[Rcpp::export]]
 arma::mat calc_psi_cpp(const Rcpp::List& K_list) {
   int M = K_list.size();
   if (M == 0) Rcpp::stop("K_list is empty");

   arma::mat K1 = Rcpp::as<arma::mat>(K_list[0]);
   int T = K1.n_rows;

   arma::mat Psi = arma::zeros(T, T);

   for (int m = 0; m < M; m++) {
     arma::mat Km = Rcpp::as<arma::mat>(K_list[m]);

     // Calculate Projection Matrix P = K(K'K)^-1K'
     // Using QR decomposition for numerical stability: K = QR => P = QQ'
     arma::mat Q, R;
     arma::qr_econ(Q, R, Km);

     Psi += Q * Q.t();
   }

   return Psi / M;
 }

 //' Circular Projection Estimator (Iterative Method)
 //'
 //' @param K_list A list of matrices (T x r) representing factor spaces.
 //' @param G_init Initial guess for Global Factors (T x r0).
 //' @param max_iter Maximum iterations.
 //' @param tol Convergence tolerance.
 //' @return A matrix (T x r0) of estimated global factors.
 // [[Rcpp::export]]
 arma::mat fit_cpe_cpp(const Rcpp::List& K_list, arma::mat G_init, int max_iter = 100, double tol = 1e-6) {
   int M = K_list.size();
   arma::mat G = G_init;
   arma::mat G_old = G;

   for (int iter = 0; iter < max_iter; iter++) {
     // Circular Projection: Project G onto K_1, then K_2, ..., then K_M
     for (int m = 0; m < M; m++) {
       arma::mat Km = Rcpp::as<arma::mat>(K_list[m]);

       // Projection P_m = K_m(K_m'K_m)^-1 K_m'
       // Efficiently: G_new = Q_m Q_m' G
       arma::mat Q, R;
       arma::qr_econ(Q, R, Km);

       G = Q * (Q.t() * G);
     }

     // Orthogonalize G to prevent collapse (Gram-Schmidt / QR)
     arma::mat Q_g, R_g;
     arma::qr_econ(Q_g, R_g, G);
     G = Q_g * sqrt(G.n_rows); // Rescale to magnitude similar to sqrt(T)

     // Check convergence (norm of difference)
     double diff = arma::norm(G - G_old, "fro");
     if (diff < tol) {
       break;
     }
     G_old = G;
   }

   return G;
 }
