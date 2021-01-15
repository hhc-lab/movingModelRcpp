drawHeatMap = function(Mij_sim, Mij_real, lambda_ini, runCounter){
  matrix.heatmap((Mij_sim - Mij_real) * 100, 
                 main = paste0("Mij_simu_prop (%) - Mij_real_prop (%) \n", 
                               "initial lambda: ", round(1/(1+exp(-lambda_ini)), 2), "\n", 
                               "Run Number ", runCounter))
}