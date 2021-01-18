#include <Rcpp.h>
using namespace Rcpp;
#include "utilities.h"

// [[Rcpp::export]]
NumericMatrix calculateSijMatrix(NumericVector Fi, NumericMatrix Pij, NumericVector lambda_sig, CharacterVector region_IDs) {
  int numpatch = Fi.length();
  NumericMatrix Sij(numpatch, numpatch);
  
  NumericVector lambda(numpatch);
  for(int i=0; i < numpatch; i++){
    lambda[i] = sigmoid(lambda_sig[i]);
  }
  
  for(int i=0; i < numpatch; i++){
    for(int j=0; j < numpatch; j++){
      if(i != j) 
        Sij(i,j) = (Fi[i] * Pij(i,j)) / (lambda[i] + Fi[i]);
      else 
        Sij(i,i) = lambda[i] / (lambda[i] + Fi[i]);
    }
  }
  
  rownames(Sij) = region_IDs;
  colnames(Sij) = region_IDs;
  
  return Sij;
}