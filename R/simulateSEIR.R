#' Simple stochastic SEIR simulation model with travel reduction
#'
#' A SEIR simulation tool using residence or contact model.
#' @param dataMatrix Residence or contact data matrix for simulation.
#' @param locationList The character vector of locations which involved in the spreading progress.
#' @param populationList The numeric vector of populations of each location in the locationList vector.
#' @param initialLoc The location of initial infectious starts (character).
#' @param model The simulation model. Must be "residence" or "contact".
#' @param reduceMethod The method which will be used in travel reduction. Must be "intercity", "intracity", or "all".
#' @param Di The duration of infectious period (days).
#' @param De The duration of exposed period (days).
#' @param nsim Number of simulations. Default is 1000.
#' @param R0_ref The base reference R0.
#' @param R0_max The maximum value of R0. Default is 3.5.
#' @param caseThreshold The threshold that stops the simulation process when the number of infectious exceed it. Default is 1000.
#' @param susMonth The time period of travel reduction (month). Default is 3.
#' @param decline A 0 to 1 numeric variable which represents the ratio of reduction. 0 means completely reduction and 1 means no reduction. Default is 0.5.
#' @param initial.I The number of initial infectious number. Default is 5.
#' @param threads A parallel computation option, setting the number of threads this function will make. Default is 1 (no parallel computing).
#' @param seed Random seed. Default is 1.
#' 
#' @keywords SEIR
#' @details
#' This tool mainly use two models: residence model and contact model. The data input of the former model is a matrix that 
#' contains the directional movement of residents in each residental regions. The data matrix of contact model contains the 
#' value of contact possibilities between any two residents from their residental regions.
#' 
#' @return 
#' A list that contains three components: \itemize{  
#' \item{probability}{: the probability that the number of infectious exceeds the threshold.} 
#' \item{accumulate_case}{: the accumulated number of infectious in each location at the end of each simulation.} 
#' \item{spend_day}{: the days costs at the end of each simulation.}
#' }
#' 
#' @examples
#' location_list = c("loc1", "loc2", "loc3", "loc4", "loc5")
#' loc_N = c(5382, 6609, 7193, 4572, 7760)
#' movMat = t(matrix(c(4445, 32, 411, 132, 362, 
#'                     223, 5441, 541, 117, 287, 
#'                     443, 729, 4969, 650, 402, 
#'                     201, 387, 346, 3456, 182, 
#'                     398, 421, 1282, 571, 5088), ncol = 5))
#' iniLoc = "loc3"
#' Rmatrix = fitMoveMatrix(movMat, location_list, loc_N, alpha_F = 1000)[["Sij"]]
#' result = simulateSEIR(Rmatrix, locationList = location_list, populationList = loc_N, initialLoc = iniLoc, 
#'                       model = "residence", Di = 3, De = 3.5, R0_ref = 2.4, threads = 4)

simulateSEIR = function(dataMatrix, locationList, populationList, initialLoc, model = c("residence", "contact"), 
                        reduceMethod = c("all", "intercity", "intracity"), Di, De, nsim = 1000, R0_ref, R0_max = 3.5, 
                        caseThreshold = 1000, susMonth = 3, decline = 0.5, initial.I = 5, threads = 1, seed = 1){
  
  ##### Check illegal conditions #####
  model = match.arg(model)
  reduceMethod = match.arg(reduceMethod)
  
  if(any(c(length(locationList), length(populationList), Di, De, nsim, R0_ref, R0_max, caseThreshold, initial.I, threads) <= 0)){
    stop("Unexpected condition occurred.")
  }
  if(length(locationList) != length(populationList)){
    stop("The length of locationList is not equal to the length of populationList.")
  }
  if(!(reduceMethod %in% c("intercity", "intracity", "all"))){
    stop("ReduceMethod must be \"intercity\", \"intracity\", or \"all\".")
  }
  if(!(initialLoc %in% locationList)){
    stop("InitialLoc must be one of the names in locationList.")
  }
  
  
  ##### Set variables #####
  numpatch = length(locationList)
  set.seed(seed)
  
  if(model == "residence"){
    R0_vec = rep(R0_ref, numpatch)
    Rmatrix = dataMatrix
    remove(R0_ref, R0_max)
  } 
  else if(model == "contact"){
    adjust= 1
    loc_N_matrix = matrix(rep(populationList, length(populationList)), byrow=TRUE, nrow=length(populationList))
    diag(loc_N_matrix) = diag(loc_N_matrix) - 1
    
    Cmatrix = dataMatrix
    
    R0_matrix = R0_ref * adjust * t(Cmatrix) * loc_N_matrix / 
                mean( sapply(locationList, function(x) Cmatrix[x, x] * (populationList[x] - 1)) ) 
    R0_matrix = apply(t(R0_matrix), MARGIN =1 , function(x) {
      if(max(x) > R0_max) x = x / max(x) * R0_max
      else x = x
    })
    
    remove(R0_ref, R0_max, adjust, loc_N_matrix, Cmatrix)
  }
  else{
    stop("Error: model must be \"residence\" or \"contact\".")
    return(NULL)
  }
  
  ##### Start simulation #####
  
  ## initial conditions 
  E1 = I1 = R1 = rep(0, numpatch)
  S1 = (populationList - I1 - E1 - R1)
  names(E1) = names(I1) = names(R1) = names(S1) = locationList
  
  ## parallel computing
  cl = makeCluster(threads)
  registerDoParallel(cl)
  
  ## Implementing travel reduction
  if(model == "residence"){
    Rmatrix_original = Rmatrix
    R0_vec_original = R0_vec
    
    ## Reduce all Rij matrix except diagonal
    if (reduceMethod == "intercity"){
      Rmatrix = Rmatrix * decline * (diag(-1, numpatch, numpatch) + 1) + diag(diag(Rmatrix))
      diag(Rmatrix) = diag(Rmatrix) + 1-rowSums(Rmatrix)
    } 
    ## Reduce R0
    else if (reduceMethod == "intracity"){
      R0_vec = R0_vec * decline
    } 
    ## Reduce all Rij matrix except diagonal and R0
    else if (reduceMethod == "all"){
      Rmatrix = Rmatrix * decline * (diag(-1, numpatch, numpatch) + 1) + diag(diag(Rmatrix))
      diag(Rmatrix) = diag(Rmatrix) + 1-rowSums(Rmatrix)
      R0_vec = R0_vec * decline
    }
  } 
  else if (model == "contact"){
    R0_matrix_original = R0_matrix
    
    ## Reduce all R0 matrix except diagonal
    if(reduceMethod == "intercity"){       
      R0_matrix = R0_matrix * decline * (diag(-1, numpatch, numpatch) + 1) + diag(diag(R0_matrix))
    } 
    ## Reduce just diag of R0 matrix
    else if(reduceMethod == "intracity"){  
      R0_matrix = R0_matrix * (diag(-1, numpatch, numpatch) + 1) + diag(diag(R0_matrix)) * decline
    } 
    ## Reduce all of R0 matrix
    else if (reduceMethod == "all"){      
      R0_matrix = R0_matrix * decline
    }
  }
  
  prob_loc = accum_loc = time_loc = NULL
  
  I1 = rep(0, numpatch)
  names(I1) = locationList
  I1[initialLoc] = initial.I
  
  SEIR_matrix = foreach (sim = 1:as.integer(nsim), .combine = "rbind", .packages = c("deSolve", "movingModelRcpp")) %dopar% {
    if(model == "residence"){
      Rmatrix_temp = Rmatrix
      R0_vec_temp = R0_vec
    }
    else if(model == "contact"){
      R0_matrix_temp = R0_matrix
    }
    
    switchFlag = TRUE
    
    S1_alter = S1
    S1_alter[initialLoc] = S1_alter[initialLoc] - initial.I
    SEIR_vector = c(S1_alter, E1, I1, R1, I1) #C1= I1
    
    day_previous = 0
    step = 2
    while (TRUE) {
      ## update event rate based on SEIR vector
      if(model == "residence")
        e_rate = event_rate_move_cpp(SEIR_vector, numpatch, Di, De, R0_vec_temp, populationList, Rmatrix_temp) 
      else if(model == "contact")
        e_rate = event_rate_coloc_cpp(SEIR_vector, numpatch, Di, De, R0_matrix_temp, populationList)
      
      delta_t = rexp(1, rate= sum(e_rate))                 # sample time, exponential distribution
      event = as.vector(rmultinom(1, 1, prob= e_rate))     # sample event, multinomial distribution
      
      change = delta_SEIR_cpp(event)
      SEIR_vector[1: (numpatch*4)] = SEIR_vector[1: (numpatch*4)] + change     # update SEIR vector
      
      ## if E + I = 0, then stop this simulation (no one can infect others)
      if((sum(SEIR_vector[(numpatch+1): (numpatch*3)]) <= 0)) return(NULL)
      
      if (sum(change[(numpatch*2+1): (numpatch*3)]) >0) {
        SEIR_vector[(numpatch*4+1): (numpatch*5)] = SEIR_vector[(numpatch*4+1): (numpatch*5)] + 
                                                    change[(numpatch*2+1): (numpatch*3)]
      }
      
      ## if the accumulated number of people who had been infected is more than set number, stop the simulation
      if(sum(SEIR_vector[(numpatch*4+1): (numpatch*5)]) >= caseThreshold)
        return(c(sim, (day_previous + delta_t), SEIR_vector))
      
      if(((day_previous + delta_t) >= (susMonth * 30.5)) & switchFlag){
        if(model == "residence"){
          Rmatrix_temp = Rmatrix_original
          R0_vec_temp = R0_vec_original
        }
        else if(model == "contact"){
          R0_matrix_temp = R0_matrix_original
        }
        
        switchFlag = FALSE
      }
      
      day_previous = day_previous + delta_t
      step = step + 1
    }  
    return(NULL)
  }
  
  if(is.null(nrow(SEIR_matrix))) 
    prob_loc = 0
  else{
    prob_loc = nrow(SEIR_matrix) / nsim
    accum_loc = SEIR_matrix[,(numpatch*4+1+2): (numpatch*5+2)]
    time_loc = as.numeric(SEIR_matrix[,2])
  }
  
  stopCluster(cl)
  return(list(probability = prob_loc, accumulate_case = accum_loc, spend_day = time_loc))
}