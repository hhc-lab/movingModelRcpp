#include <Rcpp.h>
using namespace Rcpp;

//' Calculate the event rate of contact model based on SEIR vector.
//' 
//' @param SEIR_vector The numeric vector of current SEIR status.
//' @param numpatch length of location list.
//' @param Di The duration of infectious period (days).
//' @param De The duration of exposed period (days).
//' @param R0_matrix The matrix of R0 between two locations.
//' @param N The vector of population at each location.
//' @return A vector of event rate.
// [[Rcpp::export]]
NumericVector event_rate_coloc_cpp(NumericVector SEIR_vector, int numpatch, float Di, float De,
                                   NumericMatrix R0_matrix, IntegerVector N){
  double *S = &SEIR_vector[0], *E = &SEIR_vector[numpatch], *I = &SEIR_vector[numpatch*2];
  
  // things that can happen, S->E, E->I, I->R
  NumericVector event (numpatch*3);
  for(int i=0; i < numpatch; i++){
    double eventSE = 0;
    for (int j=0; j < numpatch; j++) {
      //if(j == i) continue;
      eventSE += S[j] * I[j] * R0_matrix(j,i) / N[j] / Di;
    }
    
    event[i] = eventSE;
    event[i+numpatch] = E[i] / De;
    event[i+2*numpatch] = I[i] / Di;
  }
  
  return event;
}
