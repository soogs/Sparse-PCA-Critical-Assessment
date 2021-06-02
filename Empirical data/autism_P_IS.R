# autism data: model selection iva Index of Sparseness for sparse loadings #
# 31-may-2021 #
# soogeun park #

# seed 4360 #

# 0. load package, data, functions ####

library(RegularizedSCA)
library(rgl)

setwd("E:\\Users\\park\\Desktop\\spca_optimism_empirical\\")

source("./spca_adj_lasso.R")
source("./paper1_uslpca.R")
source("./eigenvectorCV_adapted.R")
source("./Index_of_Sparseness.R")

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




# load data #
dat <- read.table("./GSE7329_series_matrix_processed.txt",
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

# 3. seeding and sampling ####
set.seed(4360)
  
redundant <- sample(x = largerthan050, size = 1000)

dat_subset <- dat[,c(important, redundant)]

colnames(dat_subset)[1]

colnames(dat)[322]

# 4. Sparse P model selection: index of sparseness ####

nonzero_range <- c(seq(from = 5, to = 275, by = 20),
                   c(285, 295, 305, 315, 325, 335, 345),
                   seq(from = 355, to = 1015, by = 50), 
                   1017)

P_svd_IS <- data.frame(nonzero = NA, IS = NA)

for (i in 1:length(nonzero_range)){
  
  P_svd_IS[i,]$IS  <- index_sparseness(dat = dat_subset, ncomp = 3, nonzero = nonzero_range[i], 
                                       inits = "SVD", nrstart = NULL, seed_method = 1)
  
  P_svd_IS[i,]$nonzero <- nonzero_range[i]
  
  print(i)
}

IS_chosen <- P_svd_IS[which.max(P_svd_IS$IS),]


# 5. applying the chosen number of nonzero coefficients on the methods ####
# P svd #
P_svd <- spca_P_cardinality(X = dat_subset, R = 3, 
                            P = NULL, n_zeros = ncol(dat_subset) - IS_chosen$nonzero, 
                            MAXITER = 100000, stop_value = 1e-13, 
                            inits = "SVD", seed = 1)

P_svd_vaf <- 1 - sum((dat_subset - P_svd$Tmat %*% t(P_svd$P))^2) / sum(dat_subset^2) 
# 27.85% variance explained

P_svd_nonzero <- colSums(P_svd$P != 0)[1]

# P multistart #

P_multi <- list()

for (MULTI in 1:30){
  
  P_multi[[MULTI]] <- spca_P_cardinality(X = dat_subset, R = 3, 
                                         P = NULL, 
                                         n_zeros = ncol(dat_subset) - IS_chosen$nonzero, 
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
# 27.91% variance explained

# W svd #

W_svd <- spca_adj_lasso(x = dat_subset, K = 3, 
                        para = rep(IS_chosen$nonzero, 3), 
                        type = "predictor", 
                        sparse = "varnum", 
                        inits = "SVD", seed = 1, 
                        zerocolumns = F)

W_svd_zero <- colSums(W_svd$Wraw != 0)[1]

w_svd_vaf <- 1 - sum((dat_subset - dat_subset %*% W_svd$Wraw %*% t(W_svd$Pmat))^2) / sum(dat_subset^2) 
# 35.33% variance explained

# W multi #

W_multi <- list()

for (MULTI in 1:30){
  
  W_multi[[MULTI]] <- spca_adj_lasso(x = dat_subset, K = 3, 
                                     para = rep(IS_chosen$nonzero, 3), 
                                     type = "predictor", 
                                     sparse = "varnum", 
                                     inits = "multistart", 
                                     seed = MULTI, nrstart = 1, 
                                     zerocolumns = F)
  
  print(MULTI)
  
}

losses <- unlist(lapply(W_multi, function(x){x$loss}))

minloss_index <- which.min(losses) # 13th start led to the smallest loss

W_multi_final <- W_multi[[minloss_index]]


# variance accounted for: 35.7% (same as SVD-start)
W_multi_vaf <- 1 - sum((dat_subset - dat_subset %*% W_multi_final$Wraw %*% t(W_multi_final$Pmat))^2) / sum(dat_subset^2) 


# 6. comparing the results ####

W_tucker <- TuckerCoef(W_svd$Wraw, W_multi_final$Wraw)$tucker_value

WP_svd_tucker <- TuckerCoef(W_svd$Wraw, P_svd$P)$tucker_value

WP_multi_tucker <- TuckerCoef(W_multi_final$Wraw, P_multi_final$P)$tucker_value

P_tucker <-  TuckerCoef(P_svd$P, P_multi_final$P)$tucker_value

W_corrects <- corrects(W_multi_final$Wraw, W_svd$Wraw, total = ncol(dat_subset)*3)

WP_svd_corrects <- corrects(W_svd$Wraw, P_svd$P, total = ncol(dat_subset)*3)

WP_multi_corrects <- corrects(W_multi_final$Wraw, P_multi_final$P, total = ncol(dat_subset)*3)

P_corrects <- corrects(P_svd$P, P_multi_final$P, total = ncol(dat_subset)*3)

W_nonzero <- nonzero(W_multi_final$Wraw, W_svd$Wraw, nonzeros = IS_chosen$nonzero*3)

WP_svd_nonzero <- nonzero(W_svd$Wraw, P_svd$P, nonzeros = IS_chosen$nonzero*3)

WP_multi_nonzero <- nonzero(W_multi_final$Wraw, P_multi_final$P, nonzeros = IS_chosen$nonzero*3)

P_nonzero <- nonzero(P_svd$P, P_multi_final$P, nonzeros = IS_chosen$nonzero*3)



# save(seed, IS_chosen,
#      nonzero_range, P_svd_IS,
#      W_svd, W_multi, W_multi_final, P_svd, P_multi, P_multi_final,
#      file = "./autism_seed4360_P_IS_results.Rdata")

P_nonzero

W_nonzero

WP_svd_nonzero

WP_multi_nonzero


apply(W_svd$Wraw[1:107,] != 0,2,sum)
apply(W_svd$Wraw[108:1107,] != 0,2,sum)

apply(W_multi_final$Wraw[1:107,] != 0,2,sum)
apply(W_multi_final$Wraw[108:1107,] != 0,2,sum)



apply(P_svd$P[1:107,] != 0,2,sum)
apply(P_svd$P[108:1107,] != 0,2,sum)

apply(P_multi_final$P[1:107,] != 0,2,sum)
apply(P_multi_final$P[108:1107,] != 0,2,sum)



P_svd_imp <- intersect(names(which(apply(P_svd$P, 1, function(x){sum(x) != 0}))), colnames(dat_subset[,1:107]))

P_multi_imp <- intersect(names(which(apply(P_multi_final$P, 1, function(x){sum(x) != 0}))), colnames(dat_subset[,1:107]))

intersect(P_svd_imp, P_multi_imp)

length(intersect(P_svd_imp, P_multi_imp))

rownames(W_svd$Wraw) <- colnames(dat_subset)
rownames(W_multi_final$Wraw) <- colnames(dat_subset)

W_svd_imp <- intersect(names(which(apply(W_svd$Wraw, 1, function(x){sum(x) != 0}))), colnames(dat_subset[,1:107]))

W_multi_imp <- intersect(names(which(apply(W_multi_final$Wraw, 1, function(x){sum(x) != 0}))), colnames(dat_subset[,1:107]))

length(intersect(W_svd_imp, W_multi_imp))



plot(x = P_svd_IS$nonzero, y = P_svd_IS$IS)

W_multi_vaf
w_svd_vaf

P_svd_vaf
P_multi_vaf



# plotting 3d pca #
library(rgl)

# W-svd #
type_color <- as.numeric(type)

w_svd_tmat <- dat_subset %*% W_svd$Wraw

plot3d(w_svd_tmat[,c(1,3,2)], col=type_color, size=10, type='p')




# W-multi #
type_color <- as.numeric(type)

w_multi_tmat <- dat_subset %*% W_multi_final$Wraw

plot3d(w_multi_tmat[,c(1,3,2)], col=type_color, size=10, type='p')



# P-svd #
type_color <- as.numeric(type)

P_svd_tmat <- P_svd$Tmat

plot3d(P_svd_tmat[,c(1,3,2)], col=type_color, size=10, type='p')



# P-multi #
type_color <- as.numeric(type)

P_multi_tmat <- P_multi_final$Tmat

plot3d(P_multi_tmat[,c(1,3,2)], col=type_color, size=10, type='p')
