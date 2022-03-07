# autism_anova: 36 nonzero coefficients #

# the data is availble on the NCBI GEO database (accession code GSE7329)

# these results are provided in the appendix

# 0. load package, data, functions ####

library(RegularizedSCA)
library(rgl)

setwd("C:\\Users\\park\\Desktop\\research\\project_spca_models\\industry\\sim_28_jan_2020\\")

source("../../functions/spca_adj_lasso.R")
source("../../functions/paper1_uslpca.R")


# load data #
dat <- read.table("C:\\Users\\park\\Desktop\\research\\project_spca_models\\industry\\empirical\\autism\\GSE7329_series_matrix_processed.txt",
                  header = T)

# 1. data pre-processing ####

head(dat)

dat <- as.data.frame(dat)

rownames(dat) <- dat$ID_REF

dat <- t(dat)

dat <- as.data.frame(dat)

dat <- dat[-1,] # removing the ID_REF column

type <- c(rep("dup", 7), rep("FM", 8), rep("control", 15))

type <- as.factor(type)

# dat <- cbind(type, dat)

# remove some individuals:
# GSM176586, GSM176589, GSM176615

which(rownames(dat) == "GSM176586") # 12
which(rownames(dat) == "GSM176589") # 13
which(rownames(dat) == "GSM176615") # 27

dat <- dat[-c(12, 13, 27),]

type <- type[-c(12, 13, 27)]

# center and scale the variables #
dat <- scale(dat, T, T)

datmeanz <- colMeans(dat)

datvarz <- apply(dat,2,sd)

which(abs(apply(dat,2,sd) - 1) > 1e-15)

sum(abs(apply(dat,2,sd) - 1) > 1e-15) # the variables are correctly standardized

# remove the NA columns (these columsn were zero columns, but became NA colums due to scaling)
NAindex <- which(is.na(as.vector(datmeanz)))

dat <- dat[,-NAindex]

anyNA(dat)

# 2. anova session #####

hi <- aov(dat[,2] ~ type)

summary(hi)[[1]][["Pr(>F)"]]

# function that performs anova and extracts the P value #
anova_pvalue <- function(x){
  
  anovaz <- aov(x ~ type)
  
  (summary(anovaz)[[1]][["Pr(>F)"]])[1]
  
}

anova_pvalue(x = dat[,2])

anova_go <- apply(dat, 2, function(x){anova_pvalue(x)}) # this takes about 2 minutes

length(anova_go)

is.vector(anova_go)

anyNA(anova_go)

sum(anova_go < 0.05, na.rm = T) # 9682 variables.. much more than suggested by the paper

pvalues_adj <- p.adjust(as.vector(anova_go), "bonferroni")

sum(pvalues_adj < 0.05, na.rm = T) # 107 variables

important <- which(pvalues_adj < 0.05)

which(pvalues_adj > 0.5) 
# most of the bonferroni corrected pvalues are 1

length(which(pvalues_adj > 0.5) )

largerthan050 <- which(as.vector(pvalues_adj) > 0.5)

# 3. setting seed ####
set.seed(4360) # seed 4360 here

# let's take 1000 more probes to see if that is in any way feasible

redundant <- sample(x = largerthan050, size = 1000)

dat_subset <- dat[,c(important, redundant)]

colnames(dat_subset)[1]

colnames(dat)[322]

# 4. applying the methods (W-svd, W-multi, P-svd, P-multi) with 36 nonzero coefficients ####

# P svd #
P_svd <- spca_P_cardinality(X = dat_subset, R = 3, 
                            P = NULL, n_zeros = ncol(dat_subset) - 36, 
                            MAXITER = 100000, stop_value = 1e-13, 
                            inits = "SVD", seed = 1)

P_svd_vaf <- 1 - sum((dat_subset - P_svd$Tmat %*% t(P_svd$P))^2) / sum(dat_subset^2) 
# 6.344% explained

P_svd_nonzero <- colSums(P_svd$P != 0)[1]

# P multistart #

P_multi <- list()

for (MULTI in 1:50){
  
  P_multi[[MULTI]] <- spca_P_cardinality(X = dat_subset, R = 3, 
                                         P = NULL, 
                                         n_zeros = ncol(dat_subset) - 36, 
                                         MAXITER = 100000, stop_value = 1e-13, 
                                         inits = "multistart", 
                                         nrstart = 1,  seed = MULTI)
  
  print(MULTI)  
}

losses <- unlist(lapply(P_multi, function(x){x$Loss}))

minloss_index <- which.min(losses)

P_multi_final <- P_multi[[minloss_index]]

P_multi_nonzero <- colSums(P_multi_final$P != 0)[1]

P_multi_vaf <- 1 - sum((dat_subset - P_multi_final$Tmat %*% t(P_multi_final$P))^2) / sum(dat_subset^2) 
# 7.234% explained

# W svd #

W_svd <- spca_adj_lasso(x = dat_subset, K = 3, 
                        para = rep(36, 3), 
                        type = "predictor", 
                        sparse = "varnum", 
                        inits = "SVD", seed = 1, 
                        zerocolumns = F)

W_svd_zero <- colSums(W_svd$Wraw != 0)[1]

w_svd_vaf <- 1 - sum((dat_subset - dat_subset %*% W_svd$Wraw %*% t(W_svd$Pmat))^2) / sum(dat_subset^2) 
# 36.628% variance explained

# W multi #

W_multi <- list()

for (MULTI in 1:50){
  
  W_multi[[MULTI]] <- spca_adj_lasso(x = dat_subset, K = 3, 
                                     para = rep(36, 3), 
                                     type = "predictor", 
                                     sparse = "varnum", 
                                     inits = "multistart", 
                                     seed = MULTI, nrstart = 1, 
                                     zerocolumns = F)
  
  print(MULTI)
  
}

losses <- unlist(lapply(W_multi, function(x){x$loss}))

minloss_index <- which.min(losses) # 25th start led to the smallest loss

W_multi_final <- W_multi[[minloss_index]]


# variance accounted for: 36.628 (same as SVD-start)
W_multi_vaf <- 1 - sum((dat_subset - dat_subset %*% W_multi_final$Wraw %*% t(W_multi_final$Pmat))^2) / sum(dat_subset^2) 


# comparing the results between the methods #
corrects <- function(A, B, total){
  tucker <- RegularizedSCA::TuckerCoef(A, B)
  ratio <- (sum((abs(A) > 1e-7) + (abs(B[,tucker$perm]) > 1e-7) == 2) + 
              sum((abs(A) < 1e-7) + (abs(B[,tucker$perm]) < 1e-7) == 2)) / (total)
  return(ratio)
}


nonzero <- function(A, B, nonzeros){
  tucker <- RegularizedSCA::TuckerCoef(A, B)
  ratio <- sum((abs(A) > 1e-7) + 
                 (abs(B[,tucker$perm]) > 1e-7) == 2) / (nonzeros)
  return(ratio)
}

W_tucker <- TuckerCoef(W_svd$Wraw, W_multi_final$Wraw)$tucker_value

WP_svd_tucker <- TuckerCoef(W_svd$Wraw, P_svd$P)$tucker_value

WP_multi_tucker <- TuckerCoef(W_multi_final$Wraw, P_multi_final$P)$tucker_value

P_tucker <-  TuckerCoef(P_svd$P, P_multi_final$P)$tucker_value

W_corrects <- corrects(W_multi_final$Wraw, W_svd$Wraw, total = ncol(dat_subset)*3)

WP_svd_corrects <- corrects(W_svd$Wraw, P_svd$P, total = ncol(dat_subset)*3)

WP_multi_corrects <- corrects(W_multi_final$Wraw, P_multi_final$P, total = ncol(dat_subset)*3)

P_corrects <- corrects(P_svd$P, P_multi_final$P, total = ncol(dat_subset)*3)

W_nonzero <- nonzero(W_multi_final$Wraw, W_svd$Wraw, nonzeros = 36*3)

WP_svd_nonzero <- nonzero(W_svd$Wraw, P_svd$P, nonzeros = 36*3)

WP_multi_nonzero <- nonzero(W_multi_final$Wraw, P_multi_final$P, nonzeros = 36*3)

P_nonzero <- nonzero(P_svd$P, P_multi_final$P, nonzeros = 36*3)


