

spca_adj<-
  function(x,K,para,type=c("predictor","Gram"),sparse=c("penalty","varnum"),use.corr=FALSE,
           lambda=1e-6,max.iter=200,trace=FALSE,eps.conv=1e-3, inits = c("SVD", "oracle", "multistart"), nrstart = NULL, oracle = NULL)
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
    
    for (nr in 1:nrstart){
      
      if (inits == "multistart"){
        p <- ncol(x)
        v <- matrix(stats::runif(n = p*K, min = -1, max = 1), nrow = p, ncol = K)
      }
      
      # include the SVD starting values in the multistart setting, 
      # for the last multistart run
      if (inits == "multistart" & nr == nrstart){
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
      # i think this is where the sparsity is imposed..
      
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
      
      loss.history <- list()
      
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
        
        for ( i in 1:K) {
          y<-drop(x%*%alpha[,i])
          beta[,i] <- elasticnet::solvebeta(x,y,paras=c(lambda,para[i]),sparse=sparse)
        }
        # this is the exact same step to calculate the W matrix above
        # so, bascially, with the updated P matrix above at this timepoint,
        # the W is re-calculated
        
        Wmat.history[[loop.index]] <- beta
        
        normbeta<-sqrt(apply(beta^2,2,sum))
        normbeta[normbeta==0]<-1
        beta2<-t(t(beta)/normbeta)
        diff<-max(abs(beta2-temp))
        
        loss <- sum((x - x %*% beta %*% t(alpha))^2) + 
          sum(colSums(abs(beta)) * para) +
          sum(beta^2) * lambda
        # i've excluded the ridge and the lambda penalty parts from the loss
        
        reconstruct <- sum((x - x %*% beta %*% t(alpha))^2)
        
        loss.history[[loop.index]] <- loss
        
        # loss2 <- sum((x - x %*% beta2 %*% t(alpha))^2) + (lambda * sum(beta2^2))
        #
        # loss.history2[[loop.index]] <- loss2
        # # we should NOT use loss2
        # # loss2 uses normalized weights
        # # and this can be way off..
        # # 'loss' uses the raw weights and these are better
        
        temp<-beta2
        if(trace){
          if (k%%10==0){
            cat("iterations",k,fill=TRUE)
          }
        }
      }
      Wraw <- beta
      
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
                Wmat.history = Wmat.history, loss = loss, loss.history = loss.history, reconstruct = reconstruct)
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

