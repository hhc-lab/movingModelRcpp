#include <Rcpp.h>
using namespace Rcpp;

//' Calculate the event rate of residence model based on SEIR vector.
//' 
//' @param SEIR_vector The numeric vector of current SEIR status.
//' @param numpatch length of location list.
//' @param Di The duration of infectious period (days).
//' @param De The duration of exposed period (days).
//' @param R0 The vector of R0 in each location.
//' @param N The vector of population at each location.
//' @param P The matrix of possibility that a resident in location i move to locaiton j in a unit of time.
//' @return A vector of event rate.
// [[Rcpp::export]]
NumericVector event_rate_move_cpp(NumericVector SEIR_vector, int numpatch, float Di, float De,
                                  NumericVector R0, IntegerVector N, NumericMatrix P){
  /*float *S = new float[numpatch], *E = new float[numpatch], *I = new float[numpatch], *R = new float[numpatch];
  for(int i=0; i < numpatch; i++){
    S[i] = SEIR_vector[i];
    E[i] = SEIR_vector[i + numpatch];
    I[i] = SEIR_vector[i + 2*numpatch];
    R[i] = SEIR_vector[i + 3*numpatch];
  }
  */
  double *S = &SEIR_vector[0], *E = &SEIR_vector[numpatch], *I = &SEIR_vector[numpatch*2];
  
  // things that can happen, S->E, E->I, I->R
  NumericVector event (numpatch*3);
  for(int i=0; i < numpatch; i++){
    double IP = 0,  NP = 0, PIR = 0;
    for (int j=0; j < numpatch; j++) {
      if(j == i) continue;
      IP = IP + I[j] * P(j, i);
      NP = NP + N[j] * P(j, i);
      PIR = PIR + P(i, j) * I[j] / N[j] * R0[j] / Di;
    }

    event[i] = (S[i] * P(i, i) * R0[i] / Di * (I[i] * P(i, i) + IP) /
               (N[i] * P(i, i) + NP)) + (S[i] * PIR);
    event[i+numpatch] = E[i] / De;
    event[i+2*numpatch] = I[i] / Di;
  }
  // delete [] S; delete [] E; delete [] I; delete [] R;
  return event;
}