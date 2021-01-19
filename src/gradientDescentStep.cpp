#include <Rcpp.h>
#include<cmath>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List gradientDescentStep(NumericVector Fi, NumericMatrix Sij, NumericMatrix Pij, NumericVector Ni, 
                               NumericMatrix Mij_simu, NumericMatrix Mij_real, NumericVector lambda_sig, 
                               CharacterVector region_IDs, double alpha_F, double alpha_lambda_sig) {
  int numpatch = Fi.length();
  NumericVector J_Fi(numpatch), J_lambda_sig(numpatch), Fi_next(numpatch), lambda_sig_next(lambda_sig);
  
  for(int i=0; i < numpatch; i++){
    
    J_Fi[i] = ( Mij_simu(i,i) - Mij_real(i,i) ) * (-Sij(i,i) * Ni[i]);
    J_lambda_sig[i] = ( Mij_simu(i,i) - Mij_real(i,i) ) * ((1-Fi[i]) * Fi[i] * Ni[i] * exp(-lambda_sig[i])) / 
                      pow((1 + exp(-lambda_sig[i])) * Fi[i] + 1, 2);
    
    for(int j=0; j < numpatch; j++){
      if(j == i) continue;
      else{
        J_Fi[i] += ( Mij_simu(i,j) - Mij_real(i,j) ) * Pij(i,j) * Sij(i,i) * Ni[i];
        J_lambda_sig[i] += ( Mij_simu(i,i) - Mij_real(i,i) ) * (pow(Fi[i], 2) * Pij(i,j) * Ni[i] * exp(-lambda_sig[i])) / 
                           pow((1 + exp(-lambda_sig[i])) * Fi[i] + 1, 2);
      }
    }
    
    Fi_next[i] = Fi[i] - ( alpha_F * J_Fi[i] / pow(Ni[i], 2) );
    if(Fi_next[i] <= 0) Fi_next[i] = Fi[i];
    lambda_sig_next[i] = lambda_sig[i] - ( alpha_lambda_sig * J_lambda_sig[i] / pow(Ni[i], 2) );
  }
  
  Fi_next.names() = region_IDs;
  lambda_sig_next.names() = region_IDs;
  
  return List::create(Named("Fi") = Fi_next, Named("lambda_sig") = lambda_sig_next);
}