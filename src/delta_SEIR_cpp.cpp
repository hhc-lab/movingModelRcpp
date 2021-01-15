#include <Rcpp.h>
using namespace Rcpp;

//' Calculate the SEIR change based on the event rate.
//' 
//' @param event The numeric vector from event_rate function.
//' @return A vector of SEIR change.
// [[Rcpp::export]]
NumericVector delta_SEIR_cpp(NumericVector event){
  int numpatch= event.length()/3;
  NumericVector change (numpatch*4);
  
  for(int i=0; i < numpatch; i++){
    change[i] = -event[i];
    change[i+numpatch] = event[i] - event[i+numpatch];
    change[i+2*numpatch] = event[i+numpatch] - event[i+2*numpatch];
    change[i+3*numpatch] = event[i+2*numpatch];
  }

  return change;
}