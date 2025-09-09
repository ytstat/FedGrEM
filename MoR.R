.libPaths("/burg/home/yt2661/R/x86_64-pc-linux-gnu-library/4.1/")


library(gtools)
library(DirichletReg)
library(doParallel)
library(mclust)
library(mixreg)

source("FedGrEM.R")
source("estimation_error.R")
source("alignment.R")


Sys.setenv(LANG = "en_US.UTF-8")
seed <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
cat("seed=", seed, "\n")

filename <- paste("/MoR/result/", seed, ".RData", sep = "")


if (file.exists(filename)) {
  stop("Done!")
}

set.seed(seed, kind = "L'Ecuyer-CMRG")
if(Sys.getenv("SLURM_CPUS_PER_TASK") != "") {
  ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
} else {
  ncores <- detectCores()
}

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

# -------------------------------------------------------

K = 10
n = 200
p = 10
R = 5
R1 = 3
theta_bar = t(matrix(c(1,0,3,-1,1,-1,0,1,1,-1,
                       0,1,-1,-3,2,-1,2,-1,1,-1,
                       -3,-1,2,-1,2,-1,1,-3,-1,-2,
                       1,-2,0,-1,-2,2,1,3,1,-1,
                       3,1,2,-1,-2,1,2,-1,-1,2), byrow = TRUE, nrow = R))
outlier_num = 1


h_list = seq(0, 3, 0.25)

error_em <- matrix(nrow = 2, ncol = length(h_list))
error_grem <- matrix(nrow = 2, ncol = length(h_list))
error_fedem <- matrix(nrow = 2, ncol = length(h_list))
error_fedgrem <- matrix(nrow = 2, ncol = length(h_list))
error_fedgmm <- matrix(nrow = 2, ncol = length(h_list))
error_pooledem <- matrix(nrow = 2, ncol = length(h_list))
error_pooledgrem <- matrix(nrow = 2, ncol = length(h_list))

for(i in 1:length(h_list)) {
  print(i)
  
  theta <- sapply(1:(K-outlier_num), function(k){
    sapply(1:R, function(r){
      z <- rnorm(p)
      theta_bar[, r] + h_list[i]*z/sqrt(sum(z^2))
    })
  }, simplify = FALSE)
  
  x <- sapply(1:K, function(k){
    matrix(rnorm(n*p, sd = 1), nrow = n)
  }, simplify = FALSE)
  
  w <- rdirichlet(K, alpha = rep(5, R)) # task by row
  
  true_value <- list(w = w, theta = theta)
  
  y <- sapply(1:K, function(k){
    if (k <= K-outlier_num) {
      z <- sample(1:R, size = n, replace = TRUE, prob = w[k, ])
      epsilon <- rnorm(n)
      Y <- numeric(n)
      for (r in 1:R) {
        Y[z == r] <- as.vector(x[[k]][z == r, ]%*%theta[[k]][, r] + epsilon[z == r])
      }
      Y
    } else { # outlier
      as.vector(x[[k]]%*%rep(-3,p) + rnorm(n, sd = 1))
    }
  }, simplify = FALSE)
  
  
  # initilization by mixreg
  w_0 <- matrix(nrow = K, ncol = R)
  theta_0 <- rep(list(matrix(nrow = p, ncol = R)))
  for (k in 1:K) {
    fit_mixreg <- mixreg(y~x.1+x.2+x.3+x.4+x.5+x.6+x.7+x.8+x.9+x.10-1, ncomp = R, 
                               data = data.frame(x = x[[k]], y = y[[k]]))
    w_0[k, ] <- as.vector(tail(t(fit_mixreg$parmat),1))
    M <- as.matrix(head(t(fit_mixreg$parmat),-2))
    colnames(M) <- rownames(M) <- NULL
    theta_0[[k]] <- M
  }
  
  # swap
  pihat <- alignment(theta_0, method = "greedy", B = 50, ncores = ncores) 
  initial_value <- alignment_swap(pihat, list(w = w_0, theta = theta_0))
  
  error_em[1, i] <- estimation_error(initial_value$w, true_value$w, parameter = "w", outlier_index = (K-outlier_num+1):K)
  error_em[2, i] <- estimation_error(initial_value$theta, true_value$theta, parameter = "theta", outlier_index = (K-outlier_num+1):K)
  
  
  # GrEM
  grem_value <- grem(x = x, y = y, R = 5, w_0 = initial_value$w, eta = 0.1, theta_0 = initial_value$theta, ncores = ncores)
  
  error_grem[1, i] <- estimation_error(grem_value$w, true_value$w, parameter = "w", outlier_index = (K-outlier_num+1):K)
  error_grem[2, i] <- estimation_error(grem_value$theta, true_value$theta, parameter = "theta", outlier_index = (K-outlier_num+1):K)
  
  
 
  # FedEM
  fedem_value <- fedem(x = x, y = y, R = 5, w_0 = initial_value$w, theta_0 = initial_value$theta, T_num = 1000)
  
  error_fedem[1, i] <- estimation_error(fedem_value$w, true_value$w, parameter = "w", outlier_index = (K-outlier_num+1):K)
  error_fedem[2, i] <- estimation_error(fedem_value$theta, true_value$theta, parameter = "theta", outlier_index = (K-outlier_num+1):K)
  
  
  # FedGrEM
  C <- 2
  kappa_0 <- 0.1
  
  fedgrem_value <- fedgrem(x = x, y = y, w_0 = initial_value$w, theta_0 = initial_value$theta, C_eta = 0.95, C = C, kappa_0 = kappa_0,
                           T_num = 1000, iter_num = 20, R = 5)
  
  error_fedgrem[1, i] <- estimation_error(fedgrem_value$w, true_value$w, parameter = "w", outlier_index = (K-outlier_num+1):K)
  error_fedgrem[2, i] <- estimation_error(fedgrem_value$theta, true_value$theta, parameter = "theta", outlier_index = (K-outlier_num+1):K)
  
  # FedGMM
  fedgmm_value <- fedGMM(x = x, y = y, w_0 = initial_value$w, theta_0 = initial_value$theta, T_num = 1000, R = 5, R1 = 3)
  error_fedgmm[1, i] <- estimation_error(fedgmm_value$w, true_value$w, parameter = "w", outlier_index = (K-outlier_num+1):K)
  error_fedgmm[2, i] <- estimation_error(fedgmm_value$theta, true_value$theta, parameter = "theta", outlier_index = (K-outlier_num+1):K)
  
 
  
  # Pooled-EM
  x_pooled <- Reduce("rbind", x)
  y_pooled <- Reduce("c", y)
  pooledem_value <- list(w = matrix(nrow = K, ncol = R), theta = rep(list(matrix(nrow = p, ncol = R)), K))
  
  fit_mixreg_pooled <- mixreg(y~x.1+x.2+x.3+x.4+x.5+x.6+x.7+x.8+x.9+x.10-1, ncomp = R, 
                              data = data.frame(x = x_pooled, y = y_pooled))
  M <- as.matrix(head(t(fit_mixreg_pooled$parmat),-2))
  colnames(M) <- rownames(M) <- NULL
  for (k in 1:K) {
    pooledem_value$w[k, ] <- as.vector(tail(t(fit_mixreg_pooled$parmat),1))
    pooledem_value$theta[[k]] <- M
  }
  error_pooledem[1, i] <- estimation_error(pooledem_value$w, true_value$w, parameter = "w", outlier_index = (K-outlier_num+1):K)
  error_pooledem[2, i] <- estimation_error(pooledem_value$theta, true_value$theta, parameter = "theta", outlier_index = (K-outlier_num+1):K)
  
  
  # Pooled-GrEM
  pooledgrem_value_single <- grem(x = list(x_pooled), y = list(y_pooled), R = 5, w_0 = initial_value$w, eta = 0.1, theta_0 = initial_value$theta, ncores = ncores)
  pooledgrem_value <- list(w = matrix(nrow = K, ncol = R), theta = rep(list(matrix(nrow = p, ncol = R)), K))
  for (k in 1:K) {
    pooledgrem_value$w[k, ] <- pooledgrem_value_single$w
    pooledgrem_value$theta[[k]] <- pooledgrem_value_single$theta[[1]]
  }
  
  error_pooledgrem[1, i] <- estimation_error(pooledgrem_value$w, true_value$w, parameter = "w", outlier_index = (K-outlier_num+1):K)
  error_pooledgrem[2, i] <- estimation_error(pooledgrem_value$theta, true_value$theta, parameter = "theta", outlier_index = (K-outlier_num+1):K)
  
  
}


save(error_em, error_grem, error_fedem, error_fedgrem, error_fedgmm, error_pooledem, error_pooledgrem, file = filename)





