
vec_norm <- function(x) {
  sqrt(sum(x^2))
}


# Local gradient EM
grem <- function(x, y = NULL, R = 5, w_0, theta_0, eta = 0.1, T_num = 1000, ncores = 1) {
  K <- length(x)
  registerDoParallel(ncores)
  
  if(!is.null(y)){
    result <- foreach(k = 1:K) %dopar% {
      w_t <- w_0[k, ]
      theta_t <- theta_0[[k]]
      for(t in 1:T_num){
        
        # local update
        
        W <- sapply(1:R, function(r){
          w_t[r]*exp(-0.5*(as.vector((y[[k]]-x[[k]]%*%theta_t[,r])^2)))
        })
        if (!any(is.infinite(W)) && !any(is.na(W))) {
          W[W <= 1e-8] <- 1e-8
          gamma_t <- W/rowSums(W)
        } else {
          log_W <- sapply(1:R, function(r){
            sapply(1:n, function(i){
              -0.5*(as.vector((y[[k]]-x[[k]]%*%theta_t[,r])^2))
            }) + log(w_t[r])
          })
          for (i in 1:n) {
            if (any(is.infinite(W[i, ])) || any(is.na(W[i, ]))) {
              W_ind_sorted <- sort(log_W[i, ], decreasing = TRUE, index.return = TRUE)
              if (W_ind_sorted$x[1] - W_ind_sorted$x[2] > 8) {
                W[i, ] <- rep(0, R)
                W[i, W_ind_sorted$ix[1]] <- 1
              } else {
                log_W[i, ] <- log_W[i, ] - W_ind_sorted$x[2]
                W[i, ] <- exp(log_W[i, ])
              }
            }
          }
          W[W <= 1e-8] <- 1e-8
          gamma_t <- W/rowSums(W)
        }
        
        w_t <- colMeans(gamma_t)
        
        theta_t <- sapply(1:R, function(r){
          theta_t[, r] - eta*colMeans(x[[k]]*as.vector(gamma_t[, r]*(x[[k]]%*%theta_t[, r]-y[[k]])))
        })
      }
      list(w = w_t, theta = theta_t)
    }
  } else {
    result <- foreach(k = 1:K) %dopar% {
      w_t <- w_0[k, ]
      theta_t <- theta_0[[k]]
      for(t in 1:T_num){
        
        # local update
        
        W <- sapply(1:R, function(r){
          w_t[r]*exp(-0.5*(colSums((t(x[[k]])-theta_t[,r])^2)))
        })
        if (!any(is.infinite(W)) && !any(is.na(W))) {
          W[W == 0] <- 1e-8
          gamma_t <- W/rowSums(W)
        } else {
          log_W <- sapply(1:R, function(r){
            sapply(1:n, function(i){
              -0.5*colSums((t(x[[k]])-theta_t[,r])^2)
            }) + log(w_t[r])
          })
          for (i in 1:n) {
            if (any(is.infinite(W[i, ])) || any(is.na(W[i, ]))) {
              W_ind_sorted <- sort(log_W[i, ], decreasing = TRUE, index.return = TRUE)
              if (W_ind_sorted$x[1] - W_ind_sorted$x[2] > 8) {
                W[i, ] <- rep(0, R)
                W[i, W_ind_sorted$ix[1]] <- 1
              } else {
                log_W[i, ] <- log_W[i, ] - W_ind_sorted$x[2]
                W[i, ] <- exp(log_W[i, ])
              }
            }
          }
          W[W == 0] <- 1e-8
          gamma_t <- W/rowSums(W)
        }
        
        w_t <- colMeans(gamma_t)
        
        theta_t <- sapply(1:R, function(r){
          theta_t[, r] - eta*colMeans(gamma_t[, r]*t(theta_t[, r] - t(x[[k]])))
        })
      }
      list(w = w_t, theta = theta_t)
    }
  }
  
  
  w <- t(sapply(1:K, function(k){
    result[[k]]$w
  }))
  
  theta <- sapply(1:K, function(k){
    result[[k]]$theta
  }, simplify = FALSE)
  
  stopImplicitCluster()
  
  return(list(w = w, theta = theta))
}



# Local EM
em <- function(x, y = NULL, R = 5, w_0, theta_0, eta = 0.1, T_num = 1000, ncores = 1) {
  K <- length(x)
  registerDoParallel(ncores)
  
  if(!is.null(y)){
    result <- foreach(k = 1:K) %dopar% {
      w_t <- w_0[k, ]
      theta_t <- theta_0[[k]]
      for(t in 1:T_num){
        
        # local update
        
        W <- sapply(1:R, function(r){
          w_t[r]*exp(-0.5*(as.vector((y[[k]]-x[[k]]%*%theta_t[,r])^2)))
        })
        if (!any(is.infinite(W)) && !any(is.na(W))) {
          W[W <= 1e-8] <- 1e-8
          gamma_t <- W/rowSums(W)
        } else {
          log_W <- sapply(1:R, function(r){
            sapply(1:n, function(i){
              -0.5*(as.vector((y[[k]]-x[[k]]%*%theta_t[,r])^2))
            }) + log(w_t[r])
          })
          for (i in 1:n) {
            if (any(is.infinite(W[i, ])) || any(is.na(W[i, ]))) {
              W_ind_sorted <- sort(log_W[i, ], decreasing = TRUE, index.return = TRUE)
              if (W_ind_sorted$x[1] - W_ind_sorted$x[2] > 8) {
                W[i, ] <- rep(0, R)
                W[i, W_ind_sorted$ix[1]] <- 1
              } else {
                log_W[i, ] <- log_W[i, ] - W_ind_sorted$x[2]
                W[i, ] <- exp(log_W[i, ])
              }
            }
          }
          W[W <= 1e-8] <- 1e-8
          gamma_t <- W/rowSums(W)
        }
        
        w_t <- colMeans(gamma_t)
        
        theta_t <- sapply(1:R, function(r){
          theta_t[, r] - eta*colMeans(x[[k]]*as.vector(gamma_t[, r]*(x[[k]]%*%theta_t[, r]-y[[k]])))
        })
      }
      list(w = w_t, theta = theta_t)
    }
  } else {
    result <- foreach(k = 1:K) %dopar% {
      w_t <- w_0[k, ]
      theta_t <- theta_0[[k]]
      for(t in 1:T_num){
        
        # local update
        
        W <- sapply(1:R, function(r){
          w_t[r]*exp(-0.5*(colSums((t(x[[k]])-theta_t[,r])^2)))
        })
        if (!any(is.infinite(W)) && !any(is.na(W))) {
          W[W == 0] <- 1e-8
          gamma_t <- W/rowSums(W)
        } else {
          log_W <- sapply(1:R, function(r){
            sapply(1:n, function(i){
              -0.5*colSums((t(x[[k]])-theta_t[,r])^2)
            }) + log(w_t[r])
          })
          for (i in 1:n) {
            if (any(is.infinite(W[i, ])) || any(is.na(W[i, ]))) {
              W_ind_sorted <- sort(log_W[i, ], decreasing = TRUE, index.return = TRUE)
              if (W_ind_sorted$x[1] - W_ind_sorted$x[2] > 8) {
                W[i, ] <- rep(0, R)
                W[i, W_ind_sorted$ix[1]] <- 1
              } else {
                log_W[i, ] <- log_W[i, ] - W_ind_sorted$x[2]
                W[i, ] <- exp(log_W[i, ])
              }
            }
          }
          W[W == 0] <- 1e-8
          gamma_t <- W/rowSums(W)
        }
        
        w_t <- colMeans(gamma_t)
        
        theta_t <- sapply(1:R, function(r){
          colMeans(gamma_t[, r]*x[[k]])/max(w_t[r], 1e-10)
        })

      }
      list(w = w_t, theta = theta_t)
    }
  }
  
  
  w <- t(sapply(1:K, function(k){
    result[[k]]$w
  }))
  
  theta <- sapply(1:K, function(k){
    result[[k]]$theta
  }, simplify = FALSE)
  
  stopImplicitCluster()
  
  return(list(w = w, theta = theta))
}


# FedGrEM
fedgrem <- function(x, y = NULL, R = 5, w_0, theta_0, C_eta = 0.1, lambda = 0.5, T_num = 500, iter_num = 100, C = 0.2, kappa_0 = 0.5) {

  w_t <- w_0
  theta_t <- theta_0
  n <- nrow(x[[1]])
  p <- ncol(x[[1]])
  K <- length(x)

  lambda <- 1
  if(!is.null(y)) {
    for(t in 1:T_num){
      # local update
      
      gamma_t <- sapply(1:K, function(k){
        W <- sapply(1:R, function(r){
          w_t[k,r]*exp(-0.5*(as.vector((y[[k]]-x[[k]]%*%theta_t[[k]][,r])^2)))
        })
        if (!any(is.infinite(W)) && !any(is.na(W))) {
          W[W <= 1e-8] <- 1e-8
          W/rowSums(W)
        } else {
          log_W <- sapply(1:R, function(r){
            sapply(1:n, function(i){
              -0.5*(as.vector((y[[k]]-x[[k]]%*%theta_t[[k]][,r])^2))
            }) + log(w_t[k,r])
          })
          for (i in 1:n) {
            if (any(is.infinite(W[i, ])) || any(is.na(W[i, ]))) {
              W_ind_sorted <- sort(log_W[i, ], decreasing = TRUE, index.return = TRUE)
              if (W_ind_sorted$x[1] - W_ind_sorted$x[2] > 8) {
                W[i, ] <- rep(0, R)
                W[i, W_ind_sorted$ix[1]] <- 1
              } else {
                log_W[i, ] <- log_W[i, ] - W_ind_sorted$x[2]
                W[i, ] <- exp(log_W[i, ])
              }
            }
          }
          W[W <= 1e-8] <- 1e-8
          W/rowSums(W)
        }
      }, simplify = FALSE)
      
      
      w_t <- t(sapply(1:K, function(k){
        colMeans(gamma_t[[k]])
      }))
      
      ttheta_t <- sapply(1:K, function(k){
        sapply(1:R, function(r){
          # eta <- C_eta/w_t[k, r]
          eta <- C_eta/w_0[k, r]
          # eta = C_eta
          theta_t[[k]][, r] - eta*colMeans(x[[k]]*as.vector(gamma_t[[k]][, r]*(x[[k]]%*%theta_t[[k]][, r]-y[[k]])))
        })
      }, simplify = FALSE)
      
      
      # central update
      lambda <- kappa_0*lambda + C*sqrt(p+log(K))
      
      for(r in 1:R){
        ttheta_r <- sapply(1:K, function(k){
          ttheta_t[[k]][, r]
        })
        
        vbar = rowMeans(ttheta_r)
        
        delta = sapply(1:K, function(k){
          v <- max(1-(lambda/sqrt(n))/vec_norm(ttheta_r[, k]-vbar), 0)
          v*(ttheta_r[, k]-vbar)
        })
        
        for(i in iter_num) {
          vbar <- rowMeans(ttheta_r - delta)
          
          delta = sapply(1:K, function(k){
            v <- max(1-(lambda/sqrt(n))/vec_norm(ttheta_r[, k]-vbar), 0)
            v*(ttheta_r[, k]-vbar)
          })
        }
        for(k in 1:K){
          theta_t[[k]][, r] <- vbar + delta[, k]
          theta_t[[k]][theta_t[[k]][, r] > 1e4, r] <- 1e4
          theta_t[[k]][theta_t[[k]][, r] < -1e4, r] <- -1e4
        }
      }
      
    }
    return(list(w = w_t, theta = theta_t, timing_total = timing_total, timing_central = timing_central))
    
  } else{
    timing_total <- double(T_num)
    timing_central <- double(T_num)
    for(t in 1:T_num){
      t1 <- Sys.time()
      # local update
      gamma_t <- sapply(1:K, function(k){
        W <- sapply(1:R, function(r){
          w_t[k,r]*exp(-0.5*(colSums((t(x[[k]])-theta_t[[k]][,r])^2)))
        })
        if (!any(is.infinite(W)) && !any(is.na(W))) {
          W[W == 0] <- 1e-8
          W/rowSums(W)
        } else {
          log_W <- sapply(1:R, function(r){
            sapply(1:n, function(i){
              -0.5*colSums((t(x[[k]])-theta_t[[k]][,r])^2)
            }) + log(w_t[r, k])
          })
          for (i in 1:n) {
            if (any(is.infinite(W[i, ])) || any(is.na(W[i, ]))) {
              W_ind_sorted <- sort(log_W[i, ], decreasing = TRUE, index.return = TRUE)
              if (W_ind_sorted$x[1] - W_ind_sorted$x[2] > 8) {
                W[i, ] <- rep(0, R)
                W[i, W_ind_sorted$ix[1]] <- 1
              } else {
                log_W[i, ] <- log_W[i, ] - W_ind_sorted$x[2]
                W[i, ] <- exp(log_W[i, ])
              }
            }
          }
          W[W == 0] <- 1e-8
          W/rowSums(W)
        }
      }, simplify = FALSE)
      
      
      w_t <- t(sapply(1:K, function(k){
        colMeans(gamma_t[[k]])
      }))
      
      ttheta_t <- sapply(1:K, function(k){
        sapply(1:R, function(r){
          eta <- C_eta/w_0[k, r]
          theta_t[[k]][, r] - eta*colMeans(gamma_t[[k]][, r]*t(theta_t[[k]][, r] - t(x[[k]])))
        })
      }, simplify = FALSE)
      
      
      # central update
      lambda <- kappa_0*lambda + C*sqrt(p+log(K))
      t2 <- Sys.time()
      for(r in 1:R){
        ttheta_r <- sapply(1:K, function(k){
          ttheta_t[[k]][, r]
        })
        
        vbar = rowMeans(ttheta_r)
        
        delta = sapply(1:K, function(k){
          v <- max(1-(lambda/sqrt(n))/vec_norm(ttheta_r[, k]-vbar), 0)
          v*(ttheta_r[, k]-vbar)
        })

        for(i in iter_num) {
          vbar <- rowMeans(ttheta_r - delta)
          
          delta = sapply(1:K, function(k){
            v <- max(1-(lambda/sqrt(n))/vec_norm(ttheta_r[, k]-vbar), 0)
            v*(ttheta_r[, k]-vbar)
          })
        }
        
        for(k in 1:K){
          theta_t[[k]][, r] <- vbar + delta[, k]
        }
      }
      t3 <- Sys.time()
      timing_total[t] <- as.double(t3-t1)
      timing_central[t] <- as.double(t3-t2)
    }
    return(list(w = w_t, theta = theta_t, timing_total = timing_total, timing_central = timing_central))
    
  }
  
}


# FedEM
fedem <- function(x, y = NULL, R = 5, w_0, theta_0, T_num = 500, C = 0.2, kappa_0 = 0.5) {
  
  w_t <- w_0
  theta_t <- theta_0
  n <- nrow(x[[1]])
  p <- ncol(x[[1]])
  K <- length(x)
  
  if(!is.null(y)) {
    for(t in 1:T_num){
      # local update
      gamma_t <- sapply(1:K, function(k){
        W <- sapply(1:R, function(r){
          w_t[k,r]*exp(-0.5*(as.vector((y[[k]]-x[[k]]%*%theta_t[[k]][,r])^2)))
        })
        if (!any(is.infinite(W)) && !any(is.na(W))) {
          W[W == 0] <- 1e-8
          W/rowSums(W)
        } else {
          log_W <- sapply(1:R, function(r){
            sapply(1:n, function(i){
              -0.5*(as.vector((y[[k]]-x[[k]]%*%theta_t[[k]][,r])^2))
            }) + log(w_t[k,r])
          })
          for (i in 1:n) {
            if (any(is.infinite(W[i, ])) || any(is.na(W[i, ]))) {
              W_ind_sorted <- sort(log_W[i, ], decreasing = TRUE, index.return = TRUE)
              if (W_ind_sorted$x[1] - W_ind_sorted$x[2] > 8) {
                W[i, ] <- rep(0, R)
                W[i, W_ind_sorted$ix[1]] <- 1
              } else {
                log_W[i, ] <- log_W[i, ] - W_ind_sorted$x[2]
                W[i, ] <- exp(log_W[i, ])
              }
            }
          }
          W[W == 0] <- 1e-8
          W/rowSums(W)
        }
      }, simplify = FALSE)
      
      w_t <- t(sapply(1:K, function(k){
        colMeans(gamma_t[[k]])
      }))
      
      
      ttheta_t <- sapply(1:K, function(k){ 
        sapply(1:R, function(r){
          solve(t(x[[k]])%*%diag(gamma_t[[k]][, r])%*%x[[k]]+0.00001*diag(p)) %*% t(x[[k]])%*%diag(gamma_t[[k]][, r])%*%y[[k]]
        })
      }, simplify = FALSE)
      
      # central update
      ttheta_r <- sapply(1:R, function(r){
        sapply(1:K, function(k){
          ttheta_t[[k]][, r]
        })
      }, simplify = FALSE)
      
      ttheta_bar <- sapply(1:R, function(r){
        rowMeans(ttheta_r[[r]])
      })
      
      theta_t <- rep(list(ttheta_bar), K)
    }
  } else {
    for(t in 1:T_num){
      # local update
      gamma_t <- sapply(1:K, function(k){
        W <- sapply(1:R, function(r){
          w_t[k,r]*exp(-0.5*(colSums((t(x[[k]])-theta_t[[k]][,r])^2)))
        })
        if (!any(is.infinite(W)) && !any(is.na(W))) {
          W[W == 0] <- 1e-8
          W/rowSums(W)
        } else {
          log_W <- sapply(1:R, function(r){
            sapply(1:n, function(i){
              -0.5*colSums((t(x[[k]])-theta_t[[k]][,r])^2)
            }) + log(w_t[k,r])
          })
          for (i in 1:n) {
            if (any(is.infinite(W[i, ])) || any(is.na(W[i, ]))) {
              W_ind_sorted <- sort(log_W[i, ], decreasing = TRUE, index.return = TRUE)
              if (W_ind_sorted$x[1] - W_ind_sorted$x[2] > 8) {
                W[i, ] <- rep(0, R)
                W[i, W_ind_sorted$ix[1]] <- 1
              } else {
                log_W[i, ] <- log_W[i, ] - W_ind_sorted$x[2]
                W[i, ] <- exp(log_W[i, ])
              }
            }
          }
          W[W == 0] <- 1e-8
          W/rowSums(W)
        }
      }, simplify = FALSE)
      
      w_t <- t(sapply(1:K, function(k){
        colMeans(gamma_t[[k]])
      }))
      
      
      ttheta_t <- sapply(1:K, function(k){
        sapply(1:R, function(r){
          colMeans(gamma_t[[k]][, r]*x[[k]])/max(w_t[k,r], 1e-10)
        })
      }, simplify = FALSE)
      
      # central update
      ttheta_r <- sapply(1:R, function(r){
        sapply(1:K, function(k){
          ttheta_t[[k]][, r]
        })
      }, simplify = FALSE)
      
      ttheta_bar <- sapply(1:R, function(r){
        rowMeans(ttheta_r[[r]])
      })
      
      theta_t <- rep(list(ttheta_bar), K)
    }
  }
  
  return(list(w = w_t, theta = theta_t))
}



# TGMM
tem <- function(x, y = NULL, R = 5, w_0, theta_0, eta = 0.1, T_num = 1000, ncores = 1) {
  
  registerDoParallel(ncores)
  
  if(!is.null(y)){
    result <- foreach(k = 1:K) %dopar% {
      w_t <- w_0[k, ]
      theta_t <- theta_0[[k]]
      for(t in 1:T_num){
        
        # local update
        
        W <- sapply(1:R, function(r){
          w_t[r]*exp(-0.5*(as.vector((y[[k]]-x[[k]]%*%theta_t[,r])^2)))
        })
        if (!any(is.infinite(W)) && !any(is.na(W))) {
          W[W <= 1e-8] <- 1e-8
          gamma_t <- W/rowSums(W)
        } else {
          log_W <- sapply(1:R, function(r){
            sapply(1:n, function(i){
              -0.5*(as.vector((y[[k]]-x[[k]]%*%theta_t[,r])^2))
            }) + log(w_t[r])
          })
          for (i in 1:n) {
            if (any(is.infinite(W[i, ])) || any(is.na(W[i, ]))) {
              W_ind_sorted <- sort(log_W[i, ], decreasing = TRUE, index.return = TRUE)
              if (W_ind_sorted$x[1] - W_ind_sorted$x[2] > 8) {
                W[i, ] <- rep(0, R)
                W[i, W_ind_sorted$ix[1]] <- 1
              } else {
                log_W[i, ] <- log_W[i, ] - W_ind_sorted$x[2]
                W[i, ] <- exp(log_W[i, ])
              }
            }
          }
          W[W <= 1e-8] <- 1e-8
          gamma_t <- W/rowSums(W)
        }
        
        w_t <- colMeans(gamma_t)
        
        theta_t <- sapply(1:R, function(r){
          theta_t[, r] - eta*colMeans(x[[k]]*as.vector(gamma_t[, r]*(x[[k]]%*%theta_t[, r]-y[[k]])))
        })
        
      }
      list(w = w_t, theta = theta_t)
    }
  } else {
    result <- foreach(k = 1:K) %dopar% {
      w_t <- w_0[k, ]
      theta_t <- theta_0[[k]]
      for(t in 1:T_num){
        
        # local update
        
        W <- sapply(1:R, function(r){
          w_t[r]*exp(-0.5*(colSums((t(x[[k]])-theta_t[,r])^2)))
        })
        if (!any(is.infinite(W)) && !any(is.na(W))) {
          W[W == 0] <- 1e-8
          gamma_t <- W/rowSums(W)
        } else {
          log_W <- sapply(1:R, function(r){
            sapply(1:n, function(i){
              -0.5*colSums((t(x[[k]])-theta_t[,r])^2)
            }) + log(w_t[r])
          })
          for (i in 1:n) {
            if (any(is.infinite(W[i, ])) || any(is.na(W[i, ]))) {
              W_ind_sorted <- sort(log_W[i, ], decreasing = TRUE, index.return = TRUE)
              if (W_ind_sorted$x[1] - W_ind_sorted$x[2] > 8) {
                W[i, ] <- rep(0, R)
                W[i, W_ind_sorted$ix[1]] <- 1
              } else {
                log_W[i, ] <- log_W[i, ] - W_ind_sorted$x[2]
                W[i, ] <- exp(log_W[i, ])
              }
            }
          }
          W[W == 0] <- 1e-8
          gamma_t <- W/rowSums(W)
        }
        
        w_t <- colMeans(gamma_t)
        
        theta_t <- sapply(1:R, function(r){
          gamma_sum_t <- sum(gamma_t[, r])
          gamma_theta_sum_t <- rowSums(sapply(setdiff(1:K, k), function(f){
            gamma_sum_t*theta_0[[f]][, r]
          }))
          
          

          (colSums(gamma_t[, r]*x[[k]]) + 0.4*gamma_theta_sum_t)/(gamma_sum_t*(1+0.4*(K-1)))
        })
        
        
      }
      list(w = w_t, theta = theta_t)
    }
  }
  
  
  w <- t(sapply(1:K, function(k){
    result[[k]]$w
  }))
  
  theta <- sapply(1:K, function(k){
    result[[k]]$theta
  }, simplify = FALSE)
  
  stopImplicitCluster()
  
  return(list(w = w, theta = theta))
}


# FedGMM
fedGMM <- function(x, y = NULL, R = 5, w_0, theta_0, T_num = 500, R1 = 3) {
  
  w_t <- sapply(1:K, function(k){
    matrix(rep(w_0[k, ], each = R1), nrow = R1)/R1
  }, simplify = FALSE)
    
    
  theta_t <- theta_0
  n <- nrow(x[[1]])
  p <- ncol(x[[1]])
  K <- length(x)

  
  
  gamma_t <- matrix(list(matrix(nrow = n, ncol = K)), nrow = R1, ncol = R)
  
  fit_mcluster <- sapply(1:K, function(k){
    Mclust(x[[k]], G = R1, modelNames = "EEE")
  }, simplify = FALSE)
  
  ww_0 <- sapply(1:R1, function(r){
    sapply(1:K, function(k){
      fit_mcluster[[k]]$parameters$pro[r]
    })
  })

  mu_0 <- sapply(1:K, function(k){
    fit_mcluster[[k]]$parameters$mean
  }, simplify = FALSE)
  
  if (R1 > 1){
    pihat <- alignment(mu_0, method = "greedy", B = 50, ncores = ncores) 
    initial_value <- alignment_swap(pihat, list(w = ww_0, theta = mu_0))
    mu_0 <- initial_value$theta
  }
  
  
  Sigma_0 <- matrix(list(), nrow = R1, ncol = K)
  for (r1 in 1:R1) {
    for (k in 1:K) {
      Sigma_0[[r1, k]] <- fit_mcluster[[k]]$parameters$variance$Sigma
    }
  }
  
  mu_t = mu_0
  Sigma_t = Sigma_0
  

  for(t in 1:T_num){
   
    # local update
    
    for (r1 in 1:R1) {
      for (r in 1:R) {
        for (k in 1:K) {
          gamma_t[[r1, r]][, k] <- w_t[[k]][r1, r]*dmvnorm(x[[k]], mean = mu_t[[k]][, r1], sigma = Sigma_t[[r1, k]])*exp(-0.5*(as.vector((y[[k]]-x[[k]]%*%theta_t[[k]][,r])^2)))
        }
      }
    }
    
    # normalize gamma_t
    normalize_factor <- matrix(nrow = n, ncol = K)
    for (k in 1:K) {
      for (i in 1:n) {
        normalize_factor[i, k] <- sum(sapply(1:R1, function(r1){
          sapply(1:R, function(r){
            gamma_t[[r1, r]][i, k]
          })
        }))
      }
    }
    
    for (r1 in 1:R1) {
      for (r in 1:R) {
        gamma_t[[r1, r]] <- gamma_t[[r1, r]]/normalize_factor
        gamma_t[[r1, r]][gamma_t[[r1, r]] == 0] <- 1e-8
      }
    }
    
  
    w_t <- sapply(1:K, function(k){
      W <- matrix(nrow = R1, ncol = R)
      for (r1 in 1:R1) {
        for (r in 1:R) {
          W[r1, r] <- mean(gamma_t[[r1, r]][, k])
        }
      }
      W
    }, simplify = FALSE)
    
    mu_t <- sapply(1:K, function(k){
      gamma_k <- sapply(1:R1, function(r1){
        rowSums(sapply(1:R, function(r){
          gamma_t[[r1,r]][, k]
        }))
      })
      
      sapply(1:R1, function(r1){
        colSums(gamma_k[, r1]*x[[k]])/sum(gamma_k[,r1])
      })
    }, simplify = FALSE)
    
    
    for (k in 1:K) {
      gamma_k <- sapply(1:R1, function(r1){
        rowSums(sapply(1:R, function(r){
          gamma_t[[r1,r]][, k]
        }))
      })
      for (r1 in 1:R1) {
        Sigma_t[[r1, k]] <- (t(x[[k]]) - mu_t[[k]][, r1]) %*% diag(gamma_k[, r1]) %*% t(t(x[[k]]) - mu_t[[k]][, r1])/sum(gamma_k[,r1])
      }
    }
    

    ttheta_t <- sapply(1:K, function(k){ 
      gamma_k <- sapply(1:R, function(r){
        rowSums(sapply(1:R1, function(r1){
          gamma_t[[r1,r]][, k]
        }))
      })
      
      sapply(1:R, function(r){
        solve(t(x[[k]])%*%diag(gamma_k[, r])%*%x[[k]]+0.00001*diag(p)) %*% t(x[[k]])%*%diag(gamma_k[, r])%*%y[[k]]
      })
    }, simplify = FALSE)
    
    # central update
    w_matrix <- matrix(sapply(1:K, function(k){
      sapply(1:R1, function(r1){
        sum(w_t[[k]][r1,])
      })
    }), nrow = R1, ncol = K)
    
    mu_t_bar <- matrix(sapply(1:R1, function(r1){
      rowSums(sapply(1:K, function(k){
        mu_t[[k]][, r1]*w_matrix[r1, k]
      }))/sum(w_matrix[r1, ])
    }), nrow = p, ncol = R1)
    
    Sigma_t_bar <- Reduce("+", sapply(1:K, function(k){
      Sigma_t[[1, k]]
    }, simplify = FALSE))/K
    
    
    w_matrix <- matrix(sapply(1:K, function(k){
      sapply(1:R, function(r){
        sum(w_t[[k]][,r])
      })
    }), nrow = R, ncol = K)
    
   
    
    theta_bar <- Reduce("+", ttheta_t)/K
    
    mu_t <- rep(list(mu_t_bar), K)
    for (r1 in 1:R1) {
      for (k in 1:K) {
        Sigma_t[[r1, k]] <- Sigma_t_bar
      }
    }
    theta_t <- rep(list(theta_bar), K)
    
  }

  w_t <- sapply(1:R, function(r){
    sapply(1:K, function(k){
      sum(w_t[[k]][,r])
    })
  })
  
  return(list(w = w_t, theta = theta_t, mu_t = mu_t, Sigma_t_bar = Sigma_t_bar, ttheta_t = ttheta_t))
}
