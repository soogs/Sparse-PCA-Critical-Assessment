# IPIP Big 5 personality data analysis #

# data load #
library(haven)

# dat <- read_spss(file = "./IPIP120.por") 
# (publicly available from https://osf.io/wxvth/)

n_index <- c(1,1 + 5 * 1:23)
e_index <- c(2, 2+ 5 * 1:23)
o_index <- c(3, 3+ 5 * 1:23)
a_index <- c(4, 4+ 5 * 1:23)
c_index <- c(5, 5+ 5 * 1:23)


source("../../../functions/spca_adj_lasso.R")
source("../../../functions/paper1_uslpca.R")


nonzeroOCEAN <- function(p){
  return(c(sum(p[o_index] != 0), sum(p[c_index] != 0), 
           sum(p[e_index] != 0), sum(p[a_index] != 0), sum(p[n_index] != 0)))
}


results <- matrix(NA, nrow = length(samplesizes), ncol = 30)

dat_i <- dat[1:1000,]

# removing the auxiliary variables #
dat_i <- dat_i[,-c(1:10)]

dat_i <- as.matrix(dat_i)

# only centered! #
dat_i <- scale(dat_i, T, F)

# p-svd #
p_svd <- spca_P_cardinality(X = dat_i, R = 5, P=NULL, n_zeros = 120 - 24, 
                            MAXITER = 10000, inits = "SVD", seed = 11)


p_svd_P <- apply(p_svd$P,2,nonzeroOCEAN)

rownames(p_svd_P) <- c("O", "C", "E", "A", "N")

p_svd_result <- apply(p_svd_P, 2, max)

p_svd_sum <- sum(p_svd_result)

# p-multi #
p_multi <- spca_P_cardinality(X = dat_i, R = 5, P=NULL, n_zeros = 120 - 24, 
                              MAXITER = 10000, inits = "multistart",nrstart = 50, seed = 11)

p_multi <- p_multi$bunch[[which.min(p_multi$LOSS)]]

p_multi$P

p_multi_P <- apply(p_multi$P,2,nonzeroOCEAN)

p_multi_result <- apply(p_multi_P, 2, max)

p_multi_sum <- sum(p_multi_result)


# W-svd #
w_svd <- spca_adj_lasso(x = dat_i, K = 5, para = rep(24, 5), 
                        type = "predictor", sparse = "varnum", 
                        lambda = 1e-6, inits = "SVD", seed = 11, zerocolumns = F)

w_svd$Wraw

w_svd_W <- apply(w_svd$Wraw,2,nonzeroOCEAN)


w_svd_result <- apply(w_svd_W, 2, max)

w_svd_sum <- sum(w_svd_result)


# w-multi #
w_multi <- spca_adj_lasso(x = dat_i, K = 5, para = rep(24, 5), 
                          type = "predictor", sparse = "varnum", 
                          lambda = 1e-6, inits = "multistart", nrstart = 50,  seed = 2424, zerocolumns = F)

w_multi_chosen <- w_multi$bunch[[which.min(w_multi$LOSS)]]

w_multi_chosen$Wraw

w_multi_W <- apply(w_multi_chosen$Wraw,2,nonzeroOCEAN)

w_multi_result <- apply(w_multi_W, 2, max)

w_multi_sum <- sum(w_multi_result)



nonzero <- function(A, B, nonzeros){
  tucker <- RegularizedSCA::TuckerCoef(A, B)
  ratio <- sum((abs(A) > 1e-7) + 
                 (abs(B[,tucker$perm]) > 1e-7) == 2) / (nonzeros)
  return(ratio)
}

nonzero(w_multi_chosen$Wraw, w_svd$Wraw, 24*5)

# explained variance
w_svd_vaf <- 1 - sum((dat_i - dat_i %*% w_svd$Wraw %*% t(w_svd$Pmat))^2) / sum(dat_i^2)

w_multi_vaf <- 1 - sum((dat_i - dat_i %*% w_multi_chosen$Wraw %*% t(w_multi_chosen$Pmat))^2) / sum(dat_i^2)


p_svd_vaf <- 1 - sum((dat_i - p_svd$Tmat %*% t(p_svd$P))^2) / sum(dat_i^2)

p_multi_vaf <- 1 - sum((dat_i - p_multi$Tmat %*% t(p_multi$P))^2) / sum(dat_i^2)


