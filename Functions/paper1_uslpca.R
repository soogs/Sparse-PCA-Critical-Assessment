# USLPCA #
# implementation for the USLPCA method from Adachi and Trendafilov, 2016. 
# This is a sparse PCA method that makes loadings matrix sparse, via cardinality constraints.

spca_P_cardinality <- function(X, R, P = NULL, n_zeros, MAXITER, 
                               stop_value = 1e-15, avoid0rows = FALSE,
                               inits = c("SVD", "oracle", "multistart"), 
                               nrstart = NULL, 
                               seed) {

  # INPUT:
  # X:	data. centered.
  # R:	number of components
  # P: the initial rational initialization of loading matrices
  # n_zeros: number of zeros which corresponds to the level of sparseness (for each PC component)
  # MAXITER: maximum number of iteration
  # stop_value: the convergence criteria
  # note that variable JkSize and all_components have not been used in the algorithm

  
  if((inits == "SVD" || inits == "oracle") && !is.null(nrstart)){
    stop("Cannot define nrstart while inits = SVD or oracle")
  }
  
  if(inits == "oracle" && is.null(P)){
    stop("Please specify the oracle P matrix initial values")
  }
  
  ## makesparse function
  # this induces sparsity for each row of P matrix
  makesparse <- function(eachrow) {
    eachrow[sort(abs(eachrow), index.return= T)$ix[1:n_zeros]] <- 0
    return(eachrow)
  }
  
  oracleP <- P
  
  ## number of variables
  J<-ncol(X)
  
  ## number of records
  I<-nrow(X)
  
  # if inits = SVD #
  if (inits == "SVD"){
    svdP <- svd(X)

    if (R == 1){
      P <- svdP$v[,1:R] %*% as.matrix(svdP$d[1:R])
    } else {
      P <- svdP$v[,1:R] %*% diag(svdP$d[1:R])
    }
    
    nrstart <- 1
  }
  
  if (inits == "oracle"){
    nrstart <- 1
  }
  
  resultbunch <- list()
  LOSS <- c()
  
  set.seed(seed)
  
  for (nr in 1:nrstart){
    
    # if inits = multistart #
    if (inits == "multistart"){
      P <- matrix(stats::runif(n = J*R, min = -1, max = 1), nrow = J, ncol = R)
    } 
    
    # include SVD starting values for the last multistart run
    if (inits == "multistart" & nr == nrstart & nrstart > 1){
      svdP <- svd(X)
      P <- svdP$v[,1:R] %*% diag(svdP$d[1:R])
    }
    
    ## initialization of the loss function
    Loss_P <- 1e9
    
    # alternating routine
    conv <- 0
    iter <- 0
    
    Loss_P.history <- Loss_P
    
    P_old <- P
    
    # iteration starts #
    while (conv == 0)
    {
      
      # conditional estimation of T given P
      svdT <- svd(t(X %*% P_old))
      Tmat_new <- svdT$v %*% t(svdT$u)
      
      ## update corresponding to each component (described in Gu & Van Deun, 2016)
      # this P calculation makes sense because our T matrix has orthogonal columns.. meaning that
      # T'T = I
      # P <- t(solve(t(U) %*% U) %*% t(U) %*% Xc)
      # the above code is to find the regression coefficients
      # and this becomes P = t(X) %*% Tmat
      
      # conditional estimation of P given T #
      P_new <- t(X) %*% Tmat_new
      # new P is calculated
      
      ## Comparing the absolute values of all non-zero elements in P, impose the smallest n_zeros elements to zero
      # this bit is done differently from Shuai's code
      # because i would like to impose sparsity per component
      # first col of P contains the regression coefficient that connects the components that you have with the first variable
      
      if (n_zeros > 0){
        
        if (avoid0rows){
          Pdat <- as.data.frame(P_new)
          
          sparseindex <- matrix(FALSE, nrow = J, ncol = R)
          
          for(i in 1:R){
            sparseindex[,i] <- abs(Pdat[,i]) <= max(sort(abs(Pdat[,i]))[1:n_zeros])
          }
          
          if (sum(rowSums(as.matrix(sparseindex[,-R])) == (R-1)) > 0){
            
            rowtofix <- which(rowSums(as.matrix(sparseindex[,-R])) == (R-1))
            
            rnames <- rownames(abs(Pdat[-rowtofix,]))
            
            o_vec <- order(abs(Pdat[-rowtofix, R]))[1:n_zeros]
            
            sparseindex[,R] <- FALSE
            sparseindex[as.numeric(rnames[o_vec]),R] <- TRUE
            
          }
          
          P[sparseindex] <- 0
          
        } else {
          
          P_new <- apply(P_new,2,makesparse)
          
        }
      }
      
      # here, there is a transpose because the reported matrix is kind of transposed..
      # i know it is not sufficient explanation but you'll understand, once you think carefully
      
      #update the current loss values
      DEV <- X - (Tmat_new%*%t(P_new))
      DEVsq <- DEV^2
      fit <- sum(DEVsq)
      update_loss <- (abs(Loss_P - fit))
      # here, Shuai doesn't take an absolute value..
      
      # if the loss increases:
      if (fit > Loss_P){
        print ("loss increased, estimates from previous iteration taken")
        P_new <- P_old
        Tmat_new <- Tmat_old
        conv <- 1
      }
      
      Loss_P <- fit
      
      Loss_P.history <- append(Loss_P.history, fit)
      
      
      #cat("Update P Loss: ",Loss_P, sep="\n")
      
      iter <- iter + 1
      
      Tmat_old <- Tmat_new
      P_old <- P_new
      
      ## stopping criterion
      if(MAXITER==iter) conv <- 1
      if(update_loss < stop_value) conv <- 1
      
      
      
    }
    
    obj <- list(Tmat = Tmat_new, P = P_new, Loss = Loss_P, Loss_hist = Loss_P.history)
    
    resultbunch[[nr]] <- obj
    LOSS[nr] <- Loss_P
  }
  
  
  if(nr == 1){
    results <- resultbunch[[1]]
  } else {
    results <- list(bunch = resultbunch, LOSS = LOSS)
  }
  
  return (results)
}


# testing ####

# set.seed(11)
# haro <- matrix(rnorm(n = 100), nrow = 5, ncol = 20)
#
# haro.c <- scaleData(haro)
#
# svdd <- svd(haro.c)
#
# haro.P <- (svdd$v %*% diag(svdd$d))[,1:5]
#
# ## CAREFUL ##
# # P is equal to (DV')'
#
# SPCA_P(X = haro.c, P = haro.P, n_zeros = 1, MAXITER =1000)
