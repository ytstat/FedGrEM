library(gtools)

score <- function(pihat, theta_r) {
  K <- ncol(pihat)
  R <- nrow(pihat)

  sum(sapply(1:R, function(r){
    sum(sapply(1:K, function(i){
      sapply(1:K, function(j){
        vec_norm(theta_r[[i]][, pihat[r, i]] - theta_r[[j]][, pihat[r, j]])
      })
    }))
  }))
}


score_fast <- function(pihat, theta_r) {
  K <- ncol(pihat)
  R <- nrow(pihat)
  
  sum(sapply(1:R, function(r){
    sum(sapply(1:K, function(j){
        vec_norm(theta_r[[K]][, pihat[r, K]] - theta_r[[j]][, pihat[r, j]])
      }))
  }))
}



score_single <- function(theta_past, theta_kr) {
  R <- ncol(theta_past[[1]])
  
  sapply(1:R, function(r){
    sum(sapply(1:length(theta_past), function(j){
      vec_norm(theta_kr - theta_past[[j]][, r])
    }))
  })
}


alignment <- function(theta, method = c("exhaustive", "greedy"), B = 50, ncores = 1, fast = FALSE) {
  
  R <- ncol(theta[[1]])
  K <- length(theta)
  
  if (method == "exhaustive") {
    pihat.index.table <- permutations(R, R)
    task.index.table <- expand.grid(rep(list(1:nrow(pihat.index.table)), K))
    if (ncores == 1) {
      score.list <- sapply(1:nrow(task.index.table), function(i){
        pihat.table <- sapply(1:K, function(k){
          pihat.index.table[task.index.table[i, k], ]
        })
        score(pihat.table, theta)
      })
    } else {
      registerDoParallel(ncores)
      score.list <- foreach(i = 1:nrow(task.index.table), .combine = "c") %dopar% {
        pihat.table <- sapply(1:K, function(k){
          pihat.index.table[task.index.table[i, k], ]
        })
        score(pihat.table, theta)
      }
      stopImplicitCluster()
    }
    i0 <- which.min(score.list)
    pihat.table <- sapply(1:K, function(k){
      pihat.index.table[task.index.table[i0, k], ]
    })
    return(pihat.table)
  } else if (method == "greedy") { # not changed yet
    pihat.index.table <- permutations(R, R)
    if (ncores == 1) {
      pihat.list <- sapply(1:B, function(b){
        task_ind <- sample(K)
        pihat.table <- matrix(1:R)
        for (k in 2:K) {
          theta.cur <- sapply(1:k, function(j){
            theta[[task_ind[j]]]
          }, simplify = FALSE)
          
          if(!fast) {
            score.list <- sapply(1:nrow(pihat.index.table), function(i){
              pihat.table.cur <- cbind(pihat.table, pihat.index.table[i, ])
              score(pihat.table.cur, theta.cur)
            })
          } else {
            score.list <- sapply(1:nrow(pihat.index.table), function(i){
              pihat.table.cur <- cbind(pihat.table, pihat.index.table[i, ])
              score_fast(pihat.table.cur, theta.cur)
            })
          }
          
          i0 <- which.min(score.list)
          pihat.table <- cbind(pihat.table, pihat.index.table[i0, ])
        }
        pihat.reordered <- sapply(1:K, function(k){
          pihat.table[, task_ind == k]
        })
        pihat.reordered
      }, simplify = FALSE)
    } else {
      registerDoParallel(ncores)
      pihat.list <- foreach(b = 1:B) %dopar% {
        task_ind <- sample(K)
        pihat.table <- matrix(1:R)
        for (k in 2:K) {
          theta.cur <- sapply(1:k, function(j){
            theta[[task_ind[j]]]
          }, simplify = FALSE)
          
          if(!fast) {
            score.list <- sapply(1:nrow(pihat.index.table), function(i){
              pihat.table.cur <- cbind(pihat.table, pihat.index.table[i, ])
              score(pihat.table.cur, theta.cur)
            })
          } else {
            score.list <- sapply(1:nrow(pihat.index.table), function(i){
              pihat.table.cur <- cbind(pihat.table, pihat.index.table[i, ])
              score_fast(pihat.table.cur, theta.cur)
            })
          }
          # 
          # score.list <- sapply(1:nrow(pihat.index.table), function(i){
          #   pihat.table.cur <- cbind(pihat.table, pihat.index.table[i, ])
          #   score(pihat.table.cur, theta.cur)
          # })
          
          i0 <- which.min(score.list)
          pihat.table <- cbind(pihat.table, pihat.index.table[i0, ])
        }
        pihat.reordered <- sapply(1:K, function(k){
          pihat.table[, task_ind == k]
        })
        pihat.reordered
      }
      stopImplicitCluster()
    }
    score.list <- sapply(1:B, function(b){
      score(pihat.list[[b]], theta)
    })
    pihat.table <- pihat.list[[which.min(score.list)]]
    return(pihat.table)
    
  } else if (method == "greedy-single") {
    if (ncores == 1) {
      pihat.list <- sapply(1:B, function(b){
        task_ind <- sample(K)
        pihat.table <- matrix(1:R)
        for (k in 2:K) {
          theta.cur <- sapply(1:k, function(j){
            theta[[task_ind[j]]]
          }, simplify = FALSE)
          cluster_num_k <- 1:R
          if (k == 2) {
            theta_past <- list(theta.cur[[1]])
          }
          for (r in 1:R) {
            score_r <- score_single(theta_past, theta.cur[[k]][, r])
            cluster_num_k[r] <- setdiff(order(score_r), cluster_num_k[0:(r-1)])[1]
          }
          # update theta_past and pihat.table
          theta_past[[k]] <- theta_past[[1]]
          pihat.table <- cbind(pihat.table, 1:R)
          for (r in 1:R) {
            theta_past[[k]][, r] <- theta.cur[[k]][, cluster_num_k[r]]
            pihat.table[r, k] <- which(cluster_num_k == r)
          }
        }
        pihat.reordered <- sapply(1:K, function(k){
          pihat.table[, task_ind == k]
        })
        pihat.reordered
      }, simplify = FALSE)
      
      score.list <- sapply(1:B, function(b){
        score(pihat.list[[b]], theta)
      })
      pihat.table <- pihat.list[[which.min(score.list)]]
      return(pihat.table)
    }
    
  }
  
}




alignment_swap <- function(pihat, initial_value_list) {
  K <- nrow(initial_value_list$w)
  R <- ncol(initial_value_list$w)
  p <- nrow(initial_value_list$theta[[1]])
  initial_value_list_new <- initial_value_list
  for (k in 1:K) {
    for (r in 1:R) {
      initial_value_list_new$w[k, r] <- initial_value_list$w[k, pihat[r, k]]
      initial_value_list_new$theta[[k]][, r] <- initial_value_list$theta[[k]][, pihat[r, k]]
    }
  }
  return(initial_value_list_new)
}

