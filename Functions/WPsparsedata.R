# WPsparsedata - 28jan #
# 5-feb-2020 #
# soogeun #

# please refer to the document "28-jan-simulation.docx"
# the 28-jan-simulation has the goal to incorporate overlap into DGM

# a DGM with both loadings and weights sparse. 
# it is actually exactly same as the loadings sparse,
# except now the loadings are made column-orthogonal

WPsparsedata <- function(n, p, k, n_zeros = NULL, 
                         empirical = FALSE, 
                         modelseed, noiseseed, VAFx, 
                         nonzero_ones = c(TRUE, FALSE), overlap){
  
  scaleData <- function(X, value = 0){
    
    X <- scale(X, scale = FALSE)
    attr(X, "scaled:center") <- NULL
    sdX <-  apply(X, 2, function(x) sqrt( sum( x^2 ) / (length(x) - value )   ))  #compute the sd for each column
    
    sdX[sdX == 0] <- 1
    # to account for a column that is completely 0,
    # i make the sd into 1.
    
    sdX <- matrix(sdX, nrow(X), ncol(X), byrow = T)                     #put all the sd's in a matrix
    
    sdX
    
    X <- X * (1 / sdX)      #divide each entry in X by its sd
    return(X)
  }
  
  n_nonzeros <- p - n_zeros
  
  mu <- rep(0,p)
  Sigma <- diag(1,p)
  
  # first generating the T matrix
  set.seed(modelseed)
  
  # generate the X initial dataset
  Xinit <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma, empirical = FALSE)
  
  # center and normalize the columns
  Xinit <- scale(Xinit, T, F)
  
  svd1 <- svd(Xinit)
  
  tmat <- svd1$u[,1:k]
  
  tmat <- scale(tmat, T, F)
  
  # tmat orthogonalized
  # now tmat is column-orthogonal, and each column is a unit vector
  tmat <- qr.Q(qr(tmat))
  
  # P matrix is the right singular vectors
  P <- svd1$v[,1:k]
  
  # first find the biggest nonzeros
  Psize <- order(abs(P), decreasing = T)
  
  nons <- Psize[1:(n_nonzeros*2)]
  
  nons <- P[nons]
  
  Pnew <- P
  Pnew[,] <- 0
  
  Pnew1 <- sample(nons, size = n_nonzeros)
  Pnew2 <- sample(setdiff(x = nons, y = Pnew1))
  
  spot1 <- sample(1:p, n_nonzeros)
  
  if (overlap == 0){
    spot2 <- sample(setdiff(1:p, spot1), n_nonzeros)
  }
  
  if (overlap != 0){
    n_overlap <- n_nonzeros * overlap
    
    spot2 <- sample(setdiff(1:p, spot1), (n_nonzeros - n_overlap))
    
    spot2 <- append(spot2, spot1[1:n_overlap])
  }
  
  Pnew[spot1, 1] <- Pnew1
  Pnew[spot2, 2] <- Pnew2
  
  P <- Pnew
  
  if (nonzero_ones){
    P[P != 0] <- 1
  }
  
  nonzerorows <- rowSums(abs(P) > 1e-10) == 2
  
  # orthogonalize the parts that are overlapping, and insert them back
  nonzerorows_ortho <- qr.Q(qr(P[nonzerorows,]))
  
  P[nonzerorows,] <- nonzerorows_ortho
  
  # P normalized
  # normalize the columns of P
  P  <- P %*% diag(1/sqrt(diag(t(P) %*% P)))
  
  # check number of overlapping items
  overlaps <- sum(rowSums(abs(P) > 1e-10) == 2)
  
  overlaps <- overlaps / n_nonzeros
  
  # if the overlapping items do not match the number of overlaps defined,
  # stop the function
  if (overlap != overlaps){
    stop()
  }
  
  
  # variance of the components defined
  D <- diag(svd1$d[1:k])
  
  # so it's as if we have 2 big components 
  # and the remaining other eigenvalues are 0
  tdp <- tmat %*% D %*% t(P)
  
  VD <- P %*% D
  
  Ptrue <- VD
  
  Xtrue <- tmat %*% t(Ptrue)
  
  # VAFx to control the amount of noise #
  # sum of squares of the X2 dataset
  ssqXtrue <- sum(Xtrue^2)
  
  set.seed(noiseseed)
  # sample from normal distribution (Ex = Error of X)
  Ex <- matrix(MASS::mvrnorm(n = n, mu = rep(1,p), Sigma = diag(p)),nrow = n, ncol = p)
  
  # center and scale the noise matrix
  Ex <- scaleData(Ex)
  
  # sum of squares of EX
  ssqEx <- sum(Ex^2)
  
  # Rescale noise to desired level
  fx <- sqrt(ssqXtrue*(1-VAFx)/(VAFx * ssqEx))
  # fx=sqrt(ssqXtrue*(1-VAFx(vx))/(VAFx(vx)*ssqEx));
  # (matlab code - here vx is an index because VAFx in the matlab code is a vector)
  
  # 1. the VAFx and SSQ Ex are multiplied
  # 2. (1-VAFx) is multiplied with SSQ Xtrue
  # 3. ratio of these two are calculated
  # 4. square-rooted
  
  Xnew <- Xtrue + fx*Ex
  
  returnobj <- list(X = Xnew, Xtrue = Xtrue, Tmat = tmat, P = Ptrue, Portho = P,
                    D = D)
  return(returnobj)
}

# testing ####
hi <- WPsparsedata(n = 100, p = 50, k = 2, n_zeros = 50*0.5,
                   empirical = FALSE, noiseseed = sample(1:100, 1),
                   modelseed = sample(1:100, 1),
                   nonzero_ones = T, VAFx = 0.7, overlap = 0.4)
 
dim(hi$X)
sum(colMeans(hi$X)^2)
hi$P
apply(hi$P == 0,2,sum)
sum(rowSums(hi$P != 0)== 2)

colMeans(hi$Tmat)
hi$D
t(hi$Portho) %*% hi$Portho
sum(hi$Xtrue^2) / sum(hi$X^2)
sum(hi$Tmat %*% t(hi$P) - hi$Xtrue)


# hi$P
# 
# t(hi$P) %*% hi$P
# 
# hi$Portho
# 
# pca <- prcomp(hi$Xtrue)
# round(pca$sdev,5)
# sum(abs(pca$rotation[,1:2]) - abs(hi$Portho))
# 
# 
# colMeans(hi$Xtrue)
# colMeans(hi$X)
# 
# hi$X <- scaleData(hi$X)
# 
# RegularizedSCA::TuckerCoef(pca$rotation[,1:2], hi$P)
# 
# spcap <- sparks::spca_P_cardinality(X = hi$X, R = 2, P = NULL, n_zeros = 6, MAXITER = 10000)
# RegularizedSCA::TuckerCoef(spcap$P, hi$P)