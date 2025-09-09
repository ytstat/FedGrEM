
estimation_error <- function(estimated_value, true_value, parameter = c("w", "theta"), outlier_index = NULL) {
  
  if (parameter == "w") {
    R <- ncol(estimated_value)
    K <- nrow(estimated_value)
    if (is.null(outlier_index)) {
      outlier_index <- K+1
    }
    
    pi_all <- permutations(R, R)
    v <- max(sapply(setdiff(1:K, outlier_index), function(k){
      min(sapply(1:nrow(pi_all), function(i){
        max(abs(estimated_value[k, pi_all[i, ]] - true_value[k,]))
      }))
    }))
    return(v)
  } else if (parameter == "theta") {
    K <- length(estimated_value)
    R <- ncol(estimated_value[[1]])
    if (is.null(outlier_index)) {
      outlier_index <- K+1
    }
    pi_all <- permutations(R, R)
    v <- max(sapply(setdiff(1:K, outlier_index), function(k){
      min(sapply(1:nrow(pi_all), function(i){
        max(sapply(1:R, function(r){
          vec_norm(estimated_value[[k]][, pi_all[i, r]] - true_value[[k]][, r])
        }))
      }))
    }))
    return(v)
  } 
}
