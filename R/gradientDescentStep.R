gradientDescentStep = function(Fi, Sij, Pij, Ni, Mij_simu, Mij_real, lambda_sig, region_IDs, moving_date=1, alpha_F, alpha_lambda_sig){
  if(length(dim(Mij_real)) == 2){
    temp = array(NA, dim = c(1, nrow(Mij_real), ncol(Mij_real)), 
                 dimnames = list(moving_date, rownames(Mij_real), colnames(Mij_real)))
    for (name in colnames(Mij_real)) temp[1, , name] = Mij_real[, name]
    Mij_real = temp
    remove(temp, name)
  }
  
  J_Fi = array(NA, dim = length(region_IDs), dimnames = list(region_IDs))
  J_lambda_sig = array(NA, dim = length(region_IDs), dimnames = list(region_IDs))
  
  for (rName1 in region_IDs) {
    jNoti = region_IDs[-match(rName1, region_IDs)]
    
    J_Fi[rName1] = mean(
      rowSums(t((Mij_simu[rName1, jNoti] - t(Mij_real[moving_date, rName1, jNoti])) * as.numeric(Pij[rName1, jNoti]) * 
                  Sij[rName1, rName1] * Ni[rName1])) + 
        (Mij_simu[rName1, rName1] - Mij_real[moving_date, rName1, rName1]) * (-Sij[rName1, rName1] * Ni[rName1]) 
    ) 
    J_lambda_sig[rName1] = mean( 
      rowSums(t((Mij_simu[rName1, jNoti] - t(Mij_real[moving_date, rName1, jNoti])) * as.numeric(Pij[rName1, jNoti]))) * 
        Fi[rName1]^2 * Ni[rName1] * exp(-lambda_sig[rName1]) / 
        ((1 + exp(-lambda_sig[rName1])) * Fi[rName1] + 1)^2 + 
        (Mij_simu[rName1, rName1] - Mij_real[moving_date, rName1, rName1]) * 
        ((1 - Fi[rName1]) * Ni[rName1] * Fi[rName1] * exp(-lambda_sig[rName1])) / 
        ((1 + exp(-lambda_sig[rName1])) * Fi[rName1] + 1)^2
    )
    
    #message(rName1, ": \nFi: ", as.numeric(Fi[rName1]), "\nlambda_sig_i: ", as.numeric(lambda_sig[rName1]), "\n")
  }
  #message("\n===============================\n")
  Fi = Fi - alpha_F * J_Fi / Ni^2
  lambda_sig = lambda_sig - alpha_lambda_sig * J_lambda_sig / Ni^2
  # message((alpha_lambda_sig * J_lambda_sig / Ni^2)[2])
  
  return(list(Fi = Fi,lambda_sig =  lambda_sig))
}