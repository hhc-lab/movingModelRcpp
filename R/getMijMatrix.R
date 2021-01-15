getMijMatrix = function(Fi, Sij, Pij, Ni, lambda_sig, region_IDs){
  Mij_simu = diag(0, length(region_IDs), length(region_IDs))
  rownames(Mij_simu) = colnames(Mij_simu) = region_IDs
  
  for (rName1 in region_IDs) {
    jNoti = region_IDs[-match(rName1, region_IDs)]
    
    for (rName2 in region_IDs) {
      if(rName1 != rName2)
        Mij_simu[rName1, rName2] = Pij[rName1, rName2] * Fi[rName1] * Sij[rName1, rName1] * Ni[rName1] + 
                                   ((1/(1+exp(-lambda_sig[rName2]))) * Sij[rName2, rName1] * Ni[rName2])
      else if (rName1 == rName2)
        Mij_simu[rName1, rName1] = (1-Fi[rName1]) * Sij[rName1, rName1] * Ni[rName1] + 
                                   sum(Sij[jNoti, rName1] * (1-(1/(1+exp(-lambda_sig[jNoti])))) * Ni[jNoti])
    }
  }
  return(Mij_simu)
}