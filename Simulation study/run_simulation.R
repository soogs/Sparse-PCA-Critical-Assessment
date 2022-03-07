# run simulation #

# we have more evaluation criteria therefore!

# setwd("C:\\Users\\park\\Desktop\\project_spca_models\\industry\\sim_28_jan_2020\\") 
# (to be changed by the user)

source("../Functions/Psparsedata.R")
source("../Functions/WPsparsedata.R")
source("../Functions/Wsparsedata.R")

source("./conditions.R")

source("../Functions/spca_adj_lasso.R")
source("../Functions/paper1_uslpca.R")

zero <- function(A, B, nzeros){
  tucker <- RegularizedSCA::TuckerCoef(A, B)
  ratio <- sum((abs(A[,tucker$perm]) < 1e-7) + 
                 (abs(B) < 1e-7) == 2) / (nzeros)
  return(ratio)
}

nonzero <- function(A, B, nonzeros){
  tucker <- RegularizedSCA::TuckerCoef(A, B)
  ratio <- sum((abs(A[,tucker$perm]) > 1e-7) + 
                 (abs(B) > 1e-7) == 2) / (nonzeros)
  return(ratio)
}

# correct classification rate:
# (number of correct zeros + number of correct nonzeros) / total number of coefficients
corrects <- function(A, B, total){
  tucker <- RegularizedSCA::TuckerCoef(A, B)
  ratio <- (sum((abs(A[,tucker$perm]) > 1e-7) + (abs(B) > 1e-7) == 2) + 
              sum((abs(A[,tucker$perm]) < 1e-7) + (abs(B) < 1e-7) == 2)) / (total)
  return(ratio)
}

perf <- function(estimate, defined, zero_total, nonzero_total, total){
  tucker <- RegularizedSCA::TuckerCoef(defined, estimate)
  
  zerohits <- zero(estimate[, tucker$perm], defined, nzeros = zero_total)
  nonzerohits <- nonzero(estimate, defined, nonzeros = nonzero_total)
  correcthits <- corrects(estimate, defined, total = total)
  tucker <- tucker$tucker_value
  
  return(c(zerohits, nonzerohits, correcthits, tucker))
}


nrow(condition_df)

# replicating starts ####
# making the dataframe where the results are stored
results <- data.frame(matrix(NA, nrow = 3600, ncol = 51))

results[,1] <- factor(results[,1], levels = c("both", "weights", "loadings", "both_new"))
results[,2] <- factor(results[,2], levels = c("low", "high"))

colnames(results) <- c("setup", "dimensions", "ones", "vafx", "sparsity", "overlap",
                       "reps", "model_seed", "noise_seed",
                       
                       "w_zero_hits","w_nonzero_hits", "w_corrects", "w_tuckers",
                       "w_tmat_tucker", "w_recons", "w_var",
                       
                       "wm_zero_hits","wm_nonzero_hits", "wm_corrects", "wm_tuckers",
                       "wm_tmat_tucker", "wm_recons", "wm_var",
                       
                       "wo_zero_hits","wo_nonzero_hits", "wo_corrects", "wo_tuckers",
                       "wo_tmat_tucker", "wo_recons", "wo_var",
                       
                       "p_zero_hits","p_nonzero_hits", "p_corrects", "p_tuckers",
                       "p_tmat_tucker", "p_recons", "p_var",
                       
                       "pm_zero_hits","pm_nonzero_hits", "pm_corrects", "pm_tuckers",
                       "pm_tmat_tucker", "pm_recons", "pm_var",
                       
                       "po_zero_hits","po_nonzero_hits", "po_corrects", "po_tuckers",
                       "po_tmat_tucker", "po_recons", "po_var")

# we have 7200 replicates

reps_to_do <- 1:2400

set.seed(212)
model_seed <- sample(x = 1:100000, size = 10800)

noise_seed <- sample(x = 1:100000, size = 10800)

for (rrr in reps_to_do){
  
  modelseeding <- model_seed[rrr]
  noiseseeding <- noise_seed[rrr]
  
  load("./condition_df.Rdata")
  
  cond <- condition_df[rrr,]
  
  rm(condition_df)
  
  if (cond$dimensions == "high"){
    I <- 100
    J <- 500 
  }
  
  if (cond$dimensions == "low"){
    I <- 100
    J <- 50
  }
  
  cond$sparsity <- as.numeric(cond$sparsity)
  cond$vafx <- as.numeric(cond$vafx)
  cond$ones <- as.numeric(cond$ones)
  cond$reps <- as.numeric(cond$reps)
  cond$overlap <- as.numeric(cond$overlap)
  
  if (cond$setup == "both"){
    dat <- WPsparsedata(n = I, p = J, k = 2, n_zeros = (cond$sparsity * J), 
                        empirical = FALSE, 
                        modelseed = modelseeding, noiseseed = noiseseeding, 
                        VAFx = cond$vafx, nonzero_ones = cond$ones,
                        overlap = cond$overlap)
    
    if (sum(abs(dat$Portho[,1:2]) < 1e-10) != cond$sparsity * J * 2){
      stop("WPsparse - number of zeros in generated data not the same as the zeros specified by the conditions")
    }
    
    if (sum((diag(2) - t(dat$Portho) %*% dat$Portho)^2) > 1e-10){
      stop("WPsparse - the loadings matrix is not orthogonal")
    }
    
    X <- dat$X
    
    coefmat <- dat$P[,1:2]
    
    tmat <- dat$Tmat
  }
  
  if (cond$setup == "weights"){
    
    dat <- Wsparsedata(n = I, p = J, k = 2,
                       n_zeros = (cond$sparsity * J), Winit = NULL, 
                       empirical = FALSE, modelseed = modelseeding, 
                       noiseseed = noiseseeding, 
                       VAFx = cond$vafx,
                       Xnewscale = FALSE, nonzero_ones = cond$ones, overlap = cond$overlap)
    
    if (sum(abs(dat$W[,1:2]) < 1e-9) != cond$sparsity * J * 2){
      stop("wsparse - number of zeros in generated data not the same as the zeros specified by the conditions")
    }
    
    X <- dat$X
    
    coefmat <- dat$W[,1:2]
    
    tmat <- dat$Xtrue %*% dat$W
  }
  
  if (cond$setup == "loadings"){
    dat <- Psparsedata(n = I, p = J, k = 2, n_zeros = (cond$sparsity * J), 
                       empirical = FALSE, modelseed = modelseeding, 
                       noiseseed = noiseseeding,
                       VAFx = cond$vafx, nonzero_ones = cond$ones, overlap = cond$overlap)
    
    if (sum(abs(dat$P[,1:2]) < 1e-9) != cond$sparsity * J * 2){
      stop("psparse - number of zeros in generated data not the same as the zeros specified by the conditions")
    }
    
    X <- dat$X
    
    coefmat <- dat$P[,1:2]
    
    tmat <- dat$Tmat
  }
  
  zerototal <-  cond$sparsity * J * 2
  nonzerototal <- (J*2) - zerototal
  
  zero_per_component <- cond$sparsity * J
  nonzero_per_component <- J - (cond$sparsity * J)
  
  if (cond$vafx == 1){ zerocolumns <- T} else {zerocolumns <- F}
  
  # estimation ####
  
  # spca, ridge = 1e-6, default #
  w_result <- spca_adj_lasso(x = X, 
                             K = 2, 
                             para = c(nonzero_per_component, nonzero_per_component), 
                             type = "predictor", sparse = "varnum", lambda = 1e-6,
                             inits = "SVD", zerocolumns = zerocolumns)
  
  w_perf <- perf(estimate = w_result$Wraw, defined = coefmat, 
                 zero_total = zerototal, nonzero_total = nonzerototal, total = J*2)
  
  w_tmat <- X %*% w_result$Wraw
  
  w_tmat_tucker <- RegularizedSCA::TuckerCoef(w_tmat, tmat)$tucker_value
  
  w_recons <- sum((X %*% w_result$Wraw %*% t(w_result$Pmat) - X)^2) / sum((X^2))
  
  w_var <- sum(diag(t(w_tmat) %*% w_tmat)) / sum((X^2))
  
  w_perf <- c(w_perf, w_tmat_tucker, w_recons, w_var)
  
  sum((X %*% w_result$Wraw %*% t(w_result$Pmat) - X)^2)
  
  # spca, ridge = 1e-6, multistart #
  Winit <- coefmat 
  
  iter <- 0 
  
  lasso_fail <- T
  
  while (lasso_fail && iter < 5){
    
    wm_result <- tryCatch(spca_adj_lasso(x = X, 
                                         K = 2, 
                                         para = c(nonzero_per_component, nonzero_per_component), 
                                         type = "predictor", sparse = "varnum", 
                                         lambda = 1e-6, inits = "multistart", 
                                         nrstart = 20, zerocolumns = zerocolumns),
                          error = function(e) NA)
    
    iter <- iter + 1
    
    lasso_fail <- sum(is.na(wm_result)) != 0
    
    if(!lasso_fail){"multistart lasso failed, trying again"}
  }
  
  
  
  wm_perf <- perf(estimate = wm_result$bunch[[which.min(wm_result$LOSS)]]$Wraw, 
                  defined = coefmat, 
                  zero_total = zerototal, nonzero_total = nonzerototal, total = J*2)
  
  wm_tmat <- X %*% wm_result$bunch[[which.min(wm_result$LOSS)]]$Wraw
  
  wm_tmat_tucker <- RegularizedSCA::TuckerCoef(wm_tmat, tmat)$tucker_value
  
  wm_recons <- sum((X %*% wm_result$bunch[[which.min(wm_result$LOSS)]]$Wraw %*% t(wm_result$bunch[[which.min(wm_result$LOSS)]]$Pmat) - X)^2) / sum((X^2))
  
  wm_var <- sum(diag(t(wm_tmat) %*% wm_tmat)) / sum((X^2))
  
  wm_perf <- c(wm_perf, wm_tmat_tucker, wm_recons, wm_var)
  
  sum((X %*% wm_result$bunch[[which.min(wm_result$LOSS)]]$Wraw %*% t(wm_result$bunch[[which.min(wm_result$LOSS)]]$Pmat) - X)^2)
  
  # spca, ridge = 1e-6, oracle info #
  Winit <- coefmat 
  # + matrix(runif(n = (J*R), min = -0.05, max = 0.05), nrow = J)
  
  wo_result <- spca_adj_lasso(x = X, 
                              K = 2, 
                              para = c(nonzero_per_component, nonzero_per_component), 
                              type = "predictor", sparse = "varnum", 
                              lambda = 1e-6, inits = "oracle", 
                              nrstart = NULL, oracle = Winit, zerocolumns = zerocolumns)
  
  wo_perf <- perf(estimate = wo_result$Wraw, defined = coefmat, 
                  zero_total = zerototal, nonzero_total = nonzerototal, total = J*2)
  
  wo_tmat <- X %*% wo_result$Wraw
  
  wo_tmat_tucker <- RegularizedSCA::TuckerCoef(wo_tmat, tmat)$tucker_value
  
  wo_recons <- sum((X %*% wo_result$Wraw %*% t(wo_result$Pmat) - X)^2) / sum((X^2))
  
  wo_var <- sum(diag(t(wo_tmat) %*% wo_tmat)) / sum((X^2))
  
  wo_perf <- c(wo_perf, wo_tmat_tucker, wo_recons, wo_var)
  
  
  # uslpca, default #
  p_result <- spca_P_cardinality(X = X, R = 2, P = NULL, 
                                 n_zeros = zero_per_component, 
                                 MAXITER = 100000, inits = "SVD")
  
  p_perf <- perf(estimate = p_result$P, defined = coefmat, 
                 zero_total = zerototal, nonzero_total = nonzerototal, total = J*2)
  
  p_tmat <- p_result$Tmat
  
  p_tmat_tucker <- RegularizedSCA::TuckerCoef(p_tmat, tmat)$tucker_value
  
  p_recons <- sum((p_result$Tmat %*% t(p_result$P) - X)^2) / sum((X^2))
  
  p_var <- sum(diag(t(p_result$P) %*% p_result$P)) / sum((X^2))
  
  p_perf <- c(p_perf, p_tmat_tucker, p_recons, p_var)
  
  
  # uslpca, mulitstart #
  pm_result <- spca_P_cardinality(X = X, R = 2, P = coefmat, 
                                  n_zeros = zero_per_component, 
                                  MAXITER = 100000, inits = "multistart", nrstart = 20)
  
  pm_perf <- perf(estimate = pm_result$bunch[[which.min(pm_result$LOSS)]]$P, 
                  defined = coefmat, zero_total = zerototal, nonzero_total = nonzerototal, total = J*2)
  
  pm_tmat <- pm_result$bunch[[which.min(pm_result$LOSS)]]$Tmat
  
  pm_tmat_tucker <- RegularizedSCA::TuckerCoef(pm_tmat, tmat)$tucker_value
  
  pm_recons <- sum((pm_tmat %*% t(pm_result$bunch[[which.min(pm_result$LOSS)]]$P) - X)^2) / sum((X^2))
  
  pm_var <- sum(diag(t(pm_result$bunch[[which.min(pm_result$LOSS)]]$P) %*% pm_result$bunch[[which.min(pm_result$LOSS)]]$P)) / sum((X^2))
  
  pm_perf <- c(pm_perf, pm_tmat_tucker, pm_recons, pm_var)
  
  
  
  # uslpca, oracle #
  po_result <- spca_P_cardinality(X = X, R = 2, P = coefmat, 
                                  n_zeros = zero_per_component, 
                                  MAXITER = 100000, inits = "oracle")
  
  po_perf <- perf(estimate = po_result$P, defined = coefmat, 
                  zero_total = zerototal, nonzero_total = nonzerototal, total = J*2)
  
  po_tmat <- po_result$Tmat
  
  po_tmat_tucker <- RegularizedSCA::TuckerCoef(po_tmat, tmat)$tucker_value
  
  po_recons <- sum((po_result$Tmat %*% t(po_result$P) - X)^2) / sum((X^2))
  
  po_var <- sum(diag(t(po_result$P) %*% po_result$P)) / sum((X^2))
  
  po_perf <- c(po_perf, po_tmat_tucker, po_recons, po_var)
  
  
  performances <- t(matrix(c(w_perf, wm_perf, wo_perf, 
                             p_perf, pm_perf, po_perf)))
  
  resulting <- data.frame(cond$setup, cond$dimensions, cond$ones,
                          cond$vafx, cond$sparsity, cond$overlap, cond$reps,
                          modelseeding, noiseseeding)
  
  resulting <- cbind(resulting, performances)
  
  results[rrr,] <- resulting
  
  if(anyNA(results[rrr,])){
    stop("NA resulted")
  }
  
  print(rep(rrr, 10))
  
  save("results", file = "./sim_1_2400_29_july_2020.Rdata")
  
  rm(X)
  rm(coefmat)
  
  rm(w_result)
  rm(wm_result)
  rm(p_result)
  rm(pm_result)
  rm(po_result)
  
  flush.console()
  
}




