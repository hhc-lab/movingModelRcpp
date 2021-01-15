calculateSijMatrix = function(Fi, Pij, lambda_sig, region_IDs){
  Sij = (Fi * Pij) / ((1/(1+exp(-lambda_sig))) + Fi) + 
        diag((1/(1+exp(-lambda_sig))) / ((1/(1+exp(-lambda_sig))) + Fi))
  
  return(Sij)
}