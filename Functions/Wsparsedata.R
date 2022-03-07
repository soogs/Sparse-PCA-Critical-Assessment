# Wsparsedata - 28jan #
# 4-feb-2020 #
# soogeun #

# please refer to the document "28-jan-simulation.docx"

Wsparsedata <- function(n, p, k, n_zeros = NULL, Winit = NULL, empirical = FALSE, 
                        modelseed, noiseseed, VAFx, 
                        Xnewscale = FALSE, nonzero_ones = c(TRUE, FALSE), overlap){
  
  
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
  
  
  # first define the makesparse function
  # which puts the smallest elements to zero.
  makesparse <- function(eachrow) {
    eachrow[sort(abs(eachrow), index.return= T)$ix[1:n_zeros]] <- 0
    return(eachrow)
  }
  
  set.seed(modelseed)
  
  mu <- rep(0,p)
  Sigma <- diag(1,p)
  
  n_nonzeros <- p - n_zeros
  
  # generate the X initial dataset
  dat <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma, empirical = empirical)
  
  # center and normalize the columns
  dat <- scale(dat, T, F)
  
  # if you want to also normalize the columns, because you
  # might want to manipulate the eigenvalues
  # dat <- sparks::normalize(dat)
  
  # SVD on X initial
  svd1 <- svd(dat)
  
  v1 <- svd1$v[,1:k]
  
  # first find the biggest nonzeros
  v1size <- order(abs(v1), decreasing = T)
  
  nons <- v1size[1:(n_nonzeros*2)]
  
  nons <- v1[nons]
  
  v1new <- v1
  v1new[,] <- 0
  
  v1new1 <- sample(nons, size = n_nonzeros)
  v1new2 <- sample(setdiff(x = nons, y = v1new1))
  
  spot1 <- sample(1:p, n_nonzeros)
  
  if (overlap == 0){
    spot2 <- sample(setdiff(1:p, spot1), n_nonzeros)
  }
  
  if (overlap != 0){
    n_overlap <- n_nonzeros * overlap
    
    spot2 <- sample(setdiff(1:p, spot1), (n_nonzeros - n_overlap))
    
    spot2 <- append(spot2, spot1[1:n_overlap])
  }
  
  v1new[spot1, 1] <- v1new1
  v1new[spot2, 2] <- v1new2
  
  v1 <- v1new
  
  # assigning all of the nonzero elements to zeros, if that is specified
  if (nonzero_ones){
    v1[v1 != 0] <- 1
  }
  
  # normalize the columns of W
  v1  <- v1 %*% diag(1/sqrt(diag(t(v1) %*% v1)))
  
  # check overlap #
  overlaps <- sum(rowSums(abs(v1) > 1e-10) == 2)
  
  overlaps <- overlaps / n_nonzeros
  
  # if the overlapping items do not match the number of overlaps defined,
  # stop the function
  if (overlap != overlaps){
    stop()
  }
  
  
  
  # if you want to put some eigenvalues in this setup as well:
  # v1 <- v1 %*% diag(eigenvalue)
  
  # calculating the P matrix, with the least squares orthogonality constraint
  svdp1 <- svd(t(dat %*% v1) %*% dat)
  pnew1 <- t(svdp1$u %*% t(svdp1$v))
  
  # X1 <- X W P'
  # this is our Xtrue
  newdat1 <- dat %*% v1 %*% t(pnew1)
  
  # VAFx to control the amount of noise #
  # sum of squares of the X2 dataset
  ssqXtrue <- sum(newdat1^2)
  
  set.seed(noiseseed)
  # sample from normal distribution (Ex = Error of X)
  Ex <- matrix(MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma), nrow = n, ncol = p)
  
  # centering and scaling the noise matrix
  
  Ex <- scaleData(Ex)
  
  # sum of squares of EX
  ssqEx <- sum(Ex^2)
  
  # Rescale noise to desired level
  fx <- sqrt(ssqXtrue*(1-VAFx)/(VAFx * ssqEx))
  
  Xnew <- newdat1 + fx*Ex
  
  if (Xnewscale){
    Xnew <- scaleData(newdat2)
  }
  
  returnobj <- list(X = Xnew, Xtrue  = newdat1, P = pnew1, W = v1)
  
  return(returnobj)
}

# tests ####
# hi <- Wsparsedata(n = 100, p = 500, k = 2, n_zeros = 500*0.5, empirical = FALSE,
#                   modelseed = sample(1:100,1),
#                   noiseseed = sample(1:100,1), VAFx = 1, Winit = NULL,
#                   nonzero_ones = F, overlap = 0.4)
# 
# 
# dim(hi$X)
# sum(abs(colMeans(hi$X))) # checking if all of the columns are centered
# hi$W
# hi$P
# apply(hi$W != 0,2,sum) # number of non-zeros
# 
# sum(rowSums(hi$W != 0) == 2) # number of overlaps
# 
# t(hi$W) %*% hi$W # checking if columns orthogonal
# 
# t(hi$P) %*% hi$P # checking if P columns are also orthogonal
# 
# colMeans(hi$Xtrue %*% hi$W)
# sqrt(t(hi$Xtrue %*% hi$W) %*% (hi$Xtrue %*% hi$W))
# sum(hi$Xtrue^2) / sum(hi$X^2)
# sum(hi$Xtrue %*% hi$W %*% t(hi$P) - hi$Xtrue)
