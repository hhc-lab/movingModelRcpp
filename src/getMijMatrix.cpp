#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix getMijMatrix(NumericVector Fi, NumericMatrix Sij, NumericMatrix Pij, 
                              NumericVector Ni, NumericVector lambda_sig, 
                              CharacterVector region_IDs) {
  int numpatch = Fi.length();
  NumericMatrix Mij_simu(numpatch, numpatch);
  
  NumericVector lambda(numpatch);
  for(int i=0; i < numpatch; i++){
    lambda[i] = 1 / (1 + exp(-lambda_sig[i]));
  }
  
  for(int i=0; i < numpatch; i++){
    for(int j=0; j < numpatch; j++){
      if(i != j)
        Mij_simu(i,j) = (Fi[i] * Pij(i,j) * Ni[i] * lambda[i]) / (lambda[i] + Fi[i]) + 
                        (Fi[j] * Pij(j,i) * Ni[j] * lambda[j]) / (lambda[j] + Fi[j]);
      else{
        Mij_simu(i,i) = ( (1 - Fi[i]) * Ni[i] * lambda[i] ) / (lambda[i] + Fi[i]);
        for(int k = 0; k < numpatch; k++){
          if(k == i) continue;
          Mij_simu(i,i) += ( (1 - lambda[k]) * Ni[k] * Fi[k] * Pij(k,i) ) / (lambda[k] + Fi[k]);
        }
      }
        
    }
  }
  
  rownames(Mij_simu) = region_IDs;
  colnames(Mij_simu) = region_IDs;
  
  return Mij_simu;
}


