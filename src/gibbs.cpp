#include <Rcpp.h>
using namespace Rcpp;

//' A Gibbs Sampler using Rcpp
//'
//' @param n_iter Number of iterations
//' @param y Input vector
//' @param start_theta Initial value
//' @return A vector of samples
//' @export
// [[Rcpp::export]]
 NumericVector gibbs_sampler_cpp(int n_iter, NumericVector y, double start_theta) {

   // Safety check: ensure y has 4 elements
   if (y.size() != 4) {
     stop("Input vector y must have exactly 4 elements.");
   }

   // y should be a vector of length 4: (y1, y2, y3, y4)
   // observed data
   int y1 = y[0]; // 125
   int y2 = y[1]; // 18
   int y3 = y[2]; // 20
   int y4 = y[3]; // 34

   // Storage for the chain
   NumericVector theta_chain(n_iter);

   // Initialization
   double theta = start_theta;

   // Prior parameters for Beta(a, b)
   // Assuming Uniform prior: a=1, b=1
   double a = 1.0;
   double b = 1.0;

   for(int i = 0; i < n_iter; i++) {

     // --- Step 1: Data Augmentation (Update latent variable x2) ---
     // x2 | y, theta ~ Binomial(y1, theta / (2 + theta))
     double prob_x2 = theta / (2.0 + theta);

     // R::rbinom returns a double, but x2 is a count (integer concept)
     double x2 = R::rbinom(y1, prob_x2);

     // --- Step 2: Parameter Update (Update theta) ---
     // theta | x2, y ~ Beta(x2 + y4 + a, y2 + y3 + b)
     double shape1 = x2 + y4 + a;
     double shape2 = y2 + y3 + b;

     theta = R::rbeta(shape1, shape2);

     // Store result
     theta_chain[i] = theta;
   }

   return theta_chain;
 }
