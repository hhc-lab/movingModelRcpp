#' Simulation of resident moving matrix
#'
#' A tool that uses movement data to find a local optimum residents moving matrix via gradient descent method.
#' 
#' @param Mij_real Real movement data matrix.
#' @param region_IDs A character vector that contains labels or names of each region.
#' @param population A numeric vector that contains population of each region.
#' @param lambda_sig_ini Initial sigmoid lambda value. Default is 10.
#' @param alpha_F The descent speed parameter of leaving probability. Default is 1.
#' @param alpha_lambda_sig The descent speed parameter of sigmoid lambda. Default is 1.
#' @param fit_fi Whether the leaving probability should be fitted. Default is TRUE.
#' @param fit_lambda Whether lambda should be fitted. Default is FALSE.
#' @param draw_plot Whether plots will be made at the beginning and end of fitting. Default is TRUE.
#' @param threshold The threshold that determines the stop timing of descent process. 
#'                  The process will be stopped when the differences of RMS between two descents are smaller than this threshold.
#'                  Default is 1e-6. 
#' 
#' @keywords gradient descent model
#' @details test detail
#' @return 
#' A list that contains four components: \itemize{  
#' \item{Fi}{: the vector that contains the probability of a resident leaves his/her residence region.} 
#' \item{Sij}{: the matrix that contains the probability of each region a resident will go when he/she leaves his/her residence region.} 
#' \item{Mij_simu}{: the simulated movement matrix.}
#' \item{lambda_sig}{: the sigmoid lambda value.}
#' }
#' 
#' @examples
#' location_list = c("loc1", "loc2", "loc3", "loc4", "loc5")
#' loc_N = seq(1000, 9000, 2000)
#' iniLoc = "loc3"
#' movMat = t(matrix(c(990, 3, 4, 1, 2, 
#'                     23, 2911, 41, 17, 8, 
#'                     44, 102, 4201, 650, 3, 
#'                     2, 387, 1446, 4683, 482, 
#'                     79, 521, 482, 901, 7017), ncol = 5))
#' result = fitMoveMatrix(movMat, location_list, loc_N, alpha_F = 0.1)

fitMoveMatrix = function(Mij_real, region_IDs, population, lambda_sig_ini = 10, 
                         alpha_F = 1, alpha_lambda_sig = 1, fit_fi = TRUE, fit_lambda = FALSE, 
                         draw_plot = TRUE, threshold = 1e-6){
  
  rownames(Mij_real) = colnames(Mij_real) = names(population) = region_IDs
  
  if(!all(rowSums(Mij_real) == population))
    stop("The sums of each row in movement matrix must match the population vector.")
  
  ## lambda_sigmoid
  lambda_sig = rep(lambda_sig_ini, length(region_IDs))
  names(lambda_sig) = region_IDs
  
  ## The vector that contains the probability of a resident leaves his/her residence region.
  Fi = apply(Mij_real, MARGIN = 1, function(x) (1 - max(x)/sum(x)))
  names(Fi) = region_IDs
  
  ## The matrix that contains the probability of each region a resident will go when he leaves his/her residence region.
  Pij = getPijMatrix(Mij_real, region_IDs)
  
  ## The matrix that contains the probability of the place a resident will gp 
  Sij = calculateSijMatrix(Fi, Pij, lambda_sig, region_IDs)
  
  ## Get the movement matrix
  Mij_simu = getMijMatrix(Fi, Pij, population, lambda_sig, region_IDs)
  Mij_simu_prop = t( apply(Mij_simu, MARGIN = 1, function(x) x / sum(x)) )
  Mij_real_prop = t( apply(Mij_real, MARGIN = 1, function(x) x / sum(x)) )
  
  runCounter = 0
  rms_previous = Inf
  
  ## Draw heatmap
  if(draw_plot) drawHeatMap(Mij_simu_prop, Mij_real_prop, lambda_sig_ini, runCounter)
  
  runCounter = runCounter + 1
  
  ## grediant descent
  repeat {
    ## Calculate the cost function and implement gradient descent
    listTemp = gradientDescentStep(Fi, Sij, Pij, population, Mij_simu_prop, Mij_real_prop, lambda_sig, 
                                   region_IDs, alpha_F = alpha_F, alpha_lambda_sig = alpha_lambda_sig)
    if(fit_fi) Fi = listTemp$Fi
    if(fit_lambda) lambda_sig = listTemp$lambda_sig
    remove(listTemp)
    
    ## Recalculate Sij and Mij matrix
    Sij = calculateSijMatrix(Fi, Pij, lambda_sig, region_IDs)
    Mij_simu = getMijMatrix(Fi, Pij, population, lambda_sig, region_IDs)
    Mij_simu_prop = t( apply(Mij_simu, MARGIN = 1, function(x) x / sum(x)) )
    
    # rms_now = (sum((as.vector(Mij_simu_prop) - apply(Mij_real_prop, MARGIN = 1, function(x) x))^2) / length(1))^0.5
    rms_now = ( sum((Mij_simu_prop - Mij_real_prop)^2) / length(Mij_real_prop) )^0.5
    
    if(runCounter == 1) message("Round ", runCounter, " RMS: ", rms_now)
    # drawHeatMap(Mij_simu_prop, Mij_real_prop, lambda_sig_ini, runCounter)
    
    ## Set the stop poiont of descent gradient
    if((abs(rms_now - rms_previous) <= threshold) | ((rms_now - rms_previous) > 0)) break
    # plot(1/(1+exp(-lambda_sig)), type = 'l')
    rms_previous = rms_now
    runCounter = runCounter + 1
  }
  
  if(draw_plot) drawHeatMap(Mij_simu_prop, Mij_real_prop, lambda_sig_ini, runCounter)
  message("Round ", runCounter, " RMS: ", rms_now)
  
  return(list(Fi = Fi, Sij = Sij, Mij_simu = Mij_simu, lambda_sig = lambda_sig))
}
