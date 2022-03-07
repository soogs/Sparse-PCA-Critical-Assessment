# spca_adj_lasso #

# This script adjusts the spca function within the package elasticnet
# in such a way that it calculates the loss value, even if we specify the number of non-zero weights
# uses most of the code of the spca() function in elasticnet
# also uses the solvebeta function from elasticnet package.


spca_adj_lasso<-
  function(x,K,para,type=c("predictor","Gram"),sparse=c("penalty","varnum"),use.corr=FALSE,
           lambda=1e-6,max.iter=200,trace=FALSE,eps.conv=1e-3, inits = c("SVD", "oracle", "multistart"), nrstart = NULL, oracle = NULL,
           zerocolumns, seed)
  {
    
    if((inits == "SVD" || inits == "oracle") && !is.null(nrstart)){
      stop("Cannot define nrstart while inits = SVD or oracle")
    }
    
    if(inits == "oracle" && is.null(oracle)){
      stop("Please specify the oracle W matrix initial values")
    }
    
    call <- match.call()
    type <- match.arg(type)
    sparse <- match.arg(sparse)
    vn <- dimnames(x)[[2]]
    x<-switch(type,
              predictor = {
                n<-dim(x)[1]
                p<-dim(x)[2]
                if (n/p>=100){
                  cat("You may wish to restart and use a more efficient way \n")
                  cat("let the argument x be the sample covariance/correlation matrix and set type=Gram \n")
                }
                x<-scale(x,center=TRUE,scale=use.corr)
              },
              Gram = {x<-elasticnet::rootmatrix(x)}
    )
    
    svdobj<-svd(x)
    
    if (inits == "SVD"){
      v <- svdobj$v
      nrstart <- 1
    } else if (inits == "oracle"){
      v <- oracle
      
      nrstart <- 1
    }
    
    resultbunch <- list()
    LOSS <- c()
    
    set.seed(seed)
    
    for (nr in 1:nrstart){
      
      if (inits == "multistart"){
        p <- ncol(x)
        v <- matrix(stats::runif(n = p*K, min = -1, max = 1), nrow = p, ncol = K)
      }
      
      p <- ncol(x)
      
      # include the SVD starting values in the multistart setting, 
      # for the last multistart run,
      # only if more than 1 starting values is used
      if (inits == "multistart" & nr == nrstart & nrstart > 1){
        v <- svdobj$v
      }
      
      totalvariance<-sum((svdobj$d)^2)
      alpha<-as.matrix(v[,1:K,drop=FALSE])
      beta<-alpha
      
      for ( i in 1:K) {
        y<-drop(x%*%alpha[,i])
        beta[,i] <- elasticnet::solvebeta(x,y,paras=c(lambda,para[i]),sparse=sparse)
        
      }
      # the w matrix is updated here by the solvebeta function
      
      xtx<-t(x)%*%x
      temp<-beta
      normtemp<-sqrt(apply(temp^2,2,sum))
      # this calculates the norm of the column of my "beta"s = W matrix at the given timepoint
      
      normtemp[normtemp==0]<-1
      
      temp<-t(t(temp)/normtemp)
      
      k<-0
      
      diff<-1
      
      # i added these: history of the P matrix and the i index ####
      Pmat.history <- list()
      
      Wmat.history <- list()
      
      loop.index <- 0
      #####
      
      while((k<max.iter) & (diff>eps.conv)){
        # this is my i index to use for the list for the Pmat.history
        loop.index <- loop.index+1
        
        k<-k+1
        alpha<-xtx%*%beta
        z<-svd(alpha)
        
        alpha<-(z$u)%*%t(z$v)
        # this alpha is our P matrix that we are trying to get
        
        Pmat.history[[loop.index]] <- alpha
        # Pmat.history is being recorded per each iteration step ####
        
        # this is the exact same step to calculate the W matrix above
        # so, bascially, with the updated P matrix above at this timepoint,
        # the W is re-calculated
      
        for ( i in 1:K) {
          y<-drop(x%*%alpha[,i])
          beta[,i] <- elasticnet::solvebeta(x,y,paras=c(lambda,para[i]),sparse=sparse)
          
        }
        # the w matrix is updated here by the solvebeta function
      
        
        Wmat.history[[loop.index]] <- beta
        
        normbeta<-sqrt(apply(beta^2,2,sum))
        normbeta[normbeta==0]<-1
        beta2<-t(t(beta)/normbeta)
        diff<-max(abs(beta2-temp))
        
        
        temp<-beta2
        if(trace){
          if (k%%10==0){
            cat("iterations",k,fill=TRUE)
          }
        }
      }
      # when the while loop is done and the solution converges,
      # we calculate the lasso value (lambda1 penalty for the loss)
      
      xstar <- rbind(x, diag(p) * sqrt(lambda)) * (1 + lambda)^(-1/2)
      
      gamma <- rep(0, K)
      
      lasso <- rep(0, K)
      
      for ( i in 1:K) {
        
        ystar <- drop(append(x%*%alpha[,i], rep(0,p)))
        
        y <- drop(x %*% alpha[,i])
        
        res <- ystar - xstar %*% beta[,i]
        
        if (p > nrow(x)){
          gamma[i] <- max(abs(t(xstar) %*% res)) * 2
        } else {
          gamma[i] <- median(abs(t(xstar) %*% res)) * 2
        }

        lasso[i] <- gamma[i] * sqrt(1+lambda)
      }
      
      if (!zerocolumns){
        
        beta_pen <- beta
        
        beta_pen[,] <- 0
        
        beta_pen2 <- beta
        
        beta_pen2[,] <- 0
        
        lasso <- rep(0, K)
        
        gamma <- rep(0, K)
        
        for ( i in 1:K) {
          
          ystar <- drop(append(x%*%alpha[,i], rep(0,p)))
          
          y <- drop(x %*% alpha[,i])
          
          res <- ystar - xstar %*% beta[,i]
          
          if (p > nrow(x)){
            gamma[i] <- max(abs(t(xstar) %*% res)) * 2
          } else {
            gamma[i] <- median(abs(t(xstar) %*% res)) * 2
          }
          
          lasso[i] <- gamma[i] * sqrt(1+lambda)
          
          beta_pen[,i] <- elasticnet::solvebeta(xstar,ystar,paras=c(0,gamma[i]),sparse="penalty")
          
          beta_pen2[,i] <- elasticnet::solvebeta(x,y,paras=c(lambda,lasso[i]),sparse="penalty")
          
        }
        
        beta_pen2 <- beta_pen2 * sqrt(1+lambda)
        
        beta_diff <- beta_pen - beta
        beta_diff2 <- beta_pen2 - beta
        
        if (max(abs(beta_pen2 - beta)) > 1e-6){
          print("lasso value wrong; trying different method")
          
          beta_pen <- beta
          
          beta_pen[,] <- 0
          
          beta_pen2 <- beta
          
          beta_pen2[,] <- 0
          
          lasso <- rep(0,K)
          
          gamma <- rep(0,K)
          
          for ( i in 1:K) {
            
            ystar <- drop(append(x%*%alpha[,i], rep(0,p)))
            
            y <- drop(x %*% alpha[,i])
            
            res <- ystar - xstar %*% beta[,i]
            
            # opposite from the code above #
            if (p > nrow(x)){
              gamma[i] <- median(abs(t(xstar) %*% res)) * 2
            } else {
              gamma[i] <- max(abs(t(xstar) %*% res)) * 2
            }
            
            lasso[i] <- gamma[i] * sqrt(1+lambda)
            
            beta_pen[,i] <- elasticnet::solvebeta(xstar,ystar,paras=c(0,gamma[i]),sparse="penalty")
            
            beta_pen2[,i] <- elasticnet::solvebeta(x,y,paras=c(lambda,lasso[i]),sparse="penalty")
            
          }
          
          beta_pen2 <- beta_pen2 * sqrt(1+lambda)
          
          beta_diff <- beta_pen - beta
          beta_diff2 <- beta_pen2 - beta
          
        }
      
        
        if (max(abs(beta_pen2 - beta)) > 1e-6){
          stop("lasso value wrong in the end")
        }
      }
      
      # we can now calculate the exact loss #
      
      if (zerocolumns){
        loss <- sum((rbind(x %*% alpha, matrix(0, ncol = K, nrow = p)) - xstar %*% beta)^2) + sum(gamma * colSums(abs(beta)))
        
        loss <- sum((x %*% alpha - x %*% beta)^2) + lambda * sum(beta^2) + sum(lasso * colSums(abs(beta)))
        
      } else {
        loss <- sum((rbind(x %*% alpha, matrix(0, ncol = K, nrow = p)) - xstar %*% beta_pen2)^2) + sum(gamma * colSums(abs(beta_pen2)))
        
        loss <- sum((x %*% alpha - x %*% beta_pen2)^2) + lambda * sum(beta_pen^2) + sum(lasso * colSums(abs(beta_pen2)))
        
      }
      
      Wraw <- beta
      
      reconstruct <- sum((x - x %*% Wraw %*% t(alpha))^2)
      
      normbeta<-sqrt(apply(beta^2,2,sum))
      normbeta[normbeta==0]<-1
      beta<-t(t(beta)/normbeta)
      # at this step, beta <- normalized weights
      # so beta = beta2 at this step
      
      dimnames(beta)<-list(vn,paste("PC",1:K,sep=""))
      u<-x%*%beta
      R<-qr.R(qr(u))
      pev<-diag(R^2)/totalvariance
      
      # final Pmat calculation
      # alpha<-xtx%*%beta
      # z<-svd(alpha)
      # alpha<-(z$u)%*%t(z$v)
      
      obj<-list(call = call, type=type, K=K,loadings=beta,pev=pev,var.all=totalvariance,
                vn=vn,para=para,lambda=lambda, Pmat.history = Pmat.history,
                loop.index = loop.index, Pmat = alpha, Wraw = Wraw,
                Wmat.history = Wmat.history, loss = loss, reconstruct = reconstruct,
                gamma = gamma)
      class(obj) <- "spca"
      
      resultbunch[[nr]] <- obj
      LOSS[nr] <- loss
    }
    
    # beststart <- which(LOSS == min(LOSS))
    # print(beststart)
    
    if(nr == 1){
      results <- resultbunch[[1]]
    } else {
      results <- list(bunch = resultbunch, LOSS = LOSS)
    }
    
    return (results)
  }

