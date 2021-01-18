#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double sigmoid(double x) {
  return 1 / (1 + exp(-x));
}
