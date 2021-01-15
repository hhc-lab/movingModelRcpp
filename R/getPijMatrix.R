getPijMatrix = function(input_moving_matrix, region_IDs){
  P_matrix = input_moving_matrix[region_IDs, region_IDs]
  diag(P_matrix) = 0
  P_matrix = t( apply(P_matrix, MARGIN = 1, function(x) (x/sum(x, na.rm = TRUE))) )
  P_matrix = data.matrix(P_matrix)
  P_matrix[is.na(P_matrix)] = 0
  
  dimnames(P_matrix)[[1]] = region_IDs
  dimnames(P_matrix)[[2]] = region_IDs
  
  return(P_matrix)
}