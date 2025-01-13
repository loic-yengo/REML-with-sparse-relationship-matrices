## This script implements an Maximum Likelihood
library(Matrix)
library(MASS)

sparseREML_worker_AIREML <- function(y,X,G_list,niter=50,tol=1e-6,verbose=TRUE,
                                     niterEM=1,startValues=NULL,indexVg=NULL,minVC=0){ 
  N  <- length(y)
  Id <- Diagonal(n = N)
  Vy <- var(y)
  ## Initialize
  K  <- length(G_list)
  if(is.null(startValues)){
    VC <- c(rep(Vy/N,K),Vy*(1-K/N))
  }else{
    eVC <- K + 1
    nVC <- length(startValues)
    if(nVC != eVC){
      cat(paste("[Error] The number of input parameters ([`startValues=` option]) is ",
                nVC," and does not match the number of variance components expected (i.e.,",
                eVC,").\n"))
      cat("[Error] Returning NULL")
      return(NULL)
    }else{
      VC <- startValues
    }
  }  
  Vg <- VC[1:K]
  Ve <- VC[1+K]
  V  <- Ve * Id
  for(k in 1:K){
    V <- V + Vg[k] * G_list[[k]]
  }
  if(is.null(indexVg)){
    indexVg <- 1:K
  }else{
    indexVg <- indexVg
  }
  AI  <- matrix(0,nrow=1+K,ncol=1+K)
  dL  <- rep(0,1+K)
  eps <- 1
  
  ## History
  logLikHistory <- NULL
  ParamHistory  <- NULL
  
  ## Running AI-REML algorithm
  if(verbose){
    cat("Running AI-REML algorithm...\n")
  }
  for(it in 1:niter){
    ## P = V_ - V_X(X'V_X)_X'V_
    V_     <- solve(V)
    V_y    <- crossprod(V_,y)#solve(V,y)
    V_X    <- crossprod(V_,X)#solve(V,X)
    xTV_x  <- crossprod(X,V_X)
    xTV_x_ <- solve(xTV_x)
    xTV_y  <- crossprod(X,V_y)
    V_X.xTV_x_ <- V_X %*% xTV_x_
    Py     <- V_y - V_X.xTV_x_ %*% xTV_y
    
    V_Py    <- crossprod(V_,Py)
    xTV_Py  <- crossprod(X,V_Py)
    PPy     <- V_Py - V_X.xTV_x_ %*% xTV_Py
    
    APy     <- lapply(1:K, function(k) crossprod(G_list[[k]],Py))
    V_APy   <- lapply(1:K, function(k) crossprod(V_,APy[[k]]) ) # solve(V,APy[[k]]))
    xTV_APy <- lapply(1:K, function(k) crossprod(X,V_APy[[k]]) ) 
    PAPy    <- lapply(1:K, function(k) V_APy[[k]] - V_X.xTV_x_ %*% xTV_APy[[k]] )
    yPAPy   <- sapply(1:K, function(k) sum(y * PAPy[[k]]))
    yPPy    <- sum(Py^2)
    
    ## Fixed effects
    fixed_effects <- crossprod(xTV_x_,xTV_y)[,1]
    
    ## Trace calculations
    ## Tr(PA) = Tr(V_A) - Tr[ V_X(X'V_X)_X'V_A ]
    ## Tr(PA) = Tr(V_A) - Tr[(X'V_A).(V_X).(X'V_X)_]
    V_A        <- lapply(1:K, function(k) crossprod(V_,G_list[[k]]))
    xTV_A      <- lapply(1:K, function(k) crossprod(X,V_A[[k]]) )
    subTr1     <- lapply(1:K, function(k) xTV_A[[k]] %*% V_X.xTV_x_ )
    trPA       <- sapply(1:K, function(k){ sum(diag(V_A[[k]])) - sum(diag(subTr1[[k]]))})
    
    ## Tr(P) = Tr(V_) - Tr[ (X'V_V_X)(X'V_X)_ ]
    subTr2 <- crossprod(V_X) %*% xTV_x_
    trP    <- sum(diag(V_)) - sum(diag(subTr2))
    
    ## Likelihood
    ## Calculate log-likelihood
    logDet_V     <- as.numeric( determinant(V,logarithm = TRUE)$modulus )
    yPy          <- crossprod(y,Py)[1,1]
    logDet_xTV_x <- as.numeric( determinant(xTV_x,logarithm = TRUE)$modulus )
    logLik       <- -0.5 * ( logDet_V + logDet_xTV_x + yPy )
    
    ## Record history
    logLikHistory  <- c(logLikHistory,logLik)
    ParamHistory   <- rbind(ParamHistory,VC)
    
    if(it>1){
      eps <- abs( logLik - prev_logLik )
    }
    prev_logLik <- logLik
    pf <- paste0("[it=",it,"]")
    for(k in 1:K){
      pf <- c(pf,paste0(" - Vg(",k,") = ",round(Vg[k],5)))
    }
    pf <- c(pf,paste0(" - Ve = ",round(Ve,5),
                      " - logLik = ",round(logLik,5),
                      " - eps = ",format(eps,scientific=T,digits=5)))
    if(verbose){
      cat(paste0(paste(pf,collapse = ""),".\n"))
    }
    
    ## Calculate AI matrix
    for(k in 1:K){
      for(l in k:K){
        AI[k,l] <- AI[l,k] <- sum(APy[[k]] * PAPy[[l]])
      }
      AI[k,K+1] <- AI[K+1,k] <- sum(APy[[k]] * PPy)
    }
    AI[K+1,K+1] <- sum( Py * PPy )
    AI <- 0.5 * AI
    
    ## Calculate gradient
    for(k in 1:K){
      dL[k] <- -0.5 * (trPA[k] - yPAPy[k])
    }
    dL[1+K] <- -0.5 * (trP - yPPy)
    
    ## Update parameters
    if(it<=niterEM){ # EM iteration
      dx <- VC*VC*(2*dL)/N
      VC <- VC + dx
    }else{
      dx <- solve(AI,dL)
      VC <- VC + dx
    }
    if(minVC>0){
      VC <- pmax(minVC,VC)
    }
    Vg <- VC[1:K]
    Ve <- VC[1+K]

    ## Update V
    V <- Ve * Id
    for(k in 1:K){
      V <- V + Vg[k] * G_list[[k]]
    }
    if(eps<tol){
      break
    }
  }
  
  ## Prepare output  
  COV_VC <- solve(AI)
  SE_VC  <- sqrt(diag(COV_VC))
  Vp     <- sum(Vg[indexVg]) + Ve
  
  var_Vg          <- diag(COV_VC)[indexVg]
  var_Vp          <- sum(COV_VC[c(indexVg,K+1),c(indexVg,K+1)])
  var_Vg_total    <- sum(COV_VC[indexVg,indexVg])
  
  if(length(indexVg)>1){
    cov_Vp_Vg       <- rowSums(COV_VC[indexVg,c(indexVg,K+1)])
  }else{
    cov_Vp_Vg       <- sum(COV_VC[indexVg,c(indexVg,K+1)])
  }
  
  cov_Vp_Vg_total <- sum(COV_VC[indexVg,c(indexVg,K+1)])
  
  ## Heritability and SE for each component
  h2_g     <- Vg[indexVg] / Vp
  #var_h2_g <- (h2_g^2)*( var_Vg/(Vg[indexVg]^2) - 2*cov_Vp_Vg[indexVg]/(Vg[indexVg]*Vp) + var_Vp/(Vp^2) )
  var_h2_g <- (h2_g^2)*( var_Vg/(Vg[indexVg]^2) - 2*cov_Vp_Vg/(Vg[indexVg]*Vp) + var_Vp/(Vp^2) )
  
  ## If intercept included
  colnames(ParamHistory) <- c(paste0("Varcomp",1:K),"Ve")
  ParamHistory <- cbind(Iter=1:nrow(ParamHistory),ParamHistory)
  
  logLikHistory <- cbind(Iter=1:length(logLikHistory),LogLik=logLikHistory)
  results  <- list(h2_g=rbind(Estimate=h2_g,SE=sqrt(var_h2_g)),
                   #h2_total=c(Estimate=h2_total,SE=sqrt(var_h2_total)),
                   FixedEffect=fixed_effects,
                   Vg=rbind(Estimate=Vg,SE=SE_VC[1:K]),
                   Ve=c(Estimate=Ve,SE=SE_VC[1+K]),
                   Vp=c(Estimate=Vp,SE=sqrt(var_Vp)),
                   VC=rbind(Estimate=VC,SE=SE_VC),
                   COV_VC=COV_VC,
                   N=N,
                   nparam=(length(fixed_effects) + length(VC)),
                   logLik=logLik,
                   BIC=-2*logLik + log(N) * (length(fixed_effects) + length(VC)),
                   AIC=-2*logLik + 2 * (length(fixed_effects) + length(VC)),
                   logLik_history=logLikHistory,
                   Param_history=ParamHistory)
  return(results)
}


sparseREML_worker_ML <- function(y,X,G_list,niter=50,tol=1e-6,verbose=TRUE,
                                 niterEM=1,startValues=NULL,indexVg=NULL,minVC=0){ 
  N  <- length(y)
  Id <- Diagonal(n = N)
  Vy <- var(y)
  ## Initialize
  K  <- length(G_list)
  if(is.null(startValues)){
    VC <- c(rep(Vy/N,K),Vy*(1-K/N))
  }else{
    eVC <- K + 1
    nVC <- length(startValues)
    if(nVC != eVC){
      cat(paste("[Error] The number of input parameters ([`startValues=` option]) is ",
                nVC," and does not match the number of variance components expected (i.e.,",
                eVC,").\n"))
      cat("[Error] Returning NULL")
      return(NULL)
    }else{
      VC <- startValues
    }
  }  
  Vg <- VC[1:K]
  Ve <- VC[1+K]
  V  <- Ve * Id
  for(k in 1:K){
    V <- V + Vg[k] * G_list[[k]]
  }
  if(is.null(indexVg)){
    indexVg <- 1:K
  }else{
    indexVg <- indexVg
  }
  AI  <- matrix(0,nrow=1+K,ncol=1+K)
  FIM <- matrix(0,nrow=1+K,ncol=1+K)
  dL  <- rep(0,1+K)
  eps <- 1
  
  ## History
  logLikHistory <- NULL
  ParamHistory  <- NULL
  
  ## Running AI-REML algorithm
  if(verbose){
    cat("Running ML algorithm...\n")
  }
  for(it in 1:niter){
    ## Step 1 - Optimize wrt fixed effects
    V_     <- solve(V)
    V_y    <- crossprod(V_,y)
    V_X    <- crossprod(V_,X)
    xTV_x  <- crossprod(X,V_X)
    xTV_x_ <- solve(xTV_x)
    xTV_y  <- crossprod(X,V_y)
    fixed_effects <- crossprod(xTV_x_,xTV_y)[,1]
    
    ## Step 2 - Optimize wrt random effects
    e      <- y - X%*%fixed_effects
    V_e    <- crossprod(V_,e)
    V2_e   <- crossprod(V_,V_e)

    AV_e    <- lapply(1:K, function(k) crossprod(G_list[[k]],V_e))
    V_AV_e  <- lapply(1:K, function(k) crossprod(V_,AV_e[[k]]) ) # solve(V,APy[[k]]))
    eV_AV_e <- sapply(1:K, function(k) sum(e * V_AV_e[[k]]))
    eV2_e   <- sum(V_e^2)

    ## Trace calculations
    V_A    <- lapply(1:K, function(k) crossprod(V_,G_list[[k]]))
    trV_A  <- sapply(1:K, function(k){ sum(diag(V_A[[k]])) })
    trV_   <- sum(diag(V_))
    
    ## Likelihood
    ## Calculate log-likelihood
    logDet_V     <- as.numeric( determinant(V,logarithm = TRUE)$modulus )
    eV_e         <- crossprod(e,V_e)[1,1]
    logLik       <- -0.5 * ( logDet_V + eV_e )
    
    ## Record history
    logLikHistory  <- c(logLikHistory,logLik)
    ParamHistory   <- rbind(ParamHistory,VC)
    
    if(it>1){
      eps <- abs( logLik - prev_logLik )
    }
    prev_logLik <- logLik
    pf <- paste0("[it=",it,"]")
    for(k in 1:K){
      pf <- c(pf,paste0(" - Vg(",k,") = ",round(Vg[k],5)))
    }
    pf <- c(pf,paste0(" - Ve = ",round(Ve,5),
                      " - logLik = ",round(logLik,5),
                      " - eps = ",format(eps,scientific=T,digits=5)))
    if(verbose){
      cat(paste0(paste(pf,collapse = ""),".\n"))
    }
    
    ## Calculate AI matrix
    # for(k in 1:K){
    #   for(l in k:K){
    #     AI[k,l] <- AI[l,k] <- sum(AV_e[[k]] * V_AV_e[[l]])
    #   }
    #   AI[k,K+1] <- AI[K+1,k] <- sum(AV_e[[k]] * V2_e)
    # }
    # AI[K+1,K+1] <- sum( V_e * V2_e )
    # AI <- 0.5 * AI
    ## Fisher-scoring
    for(k in 1:K){
      for(l in k:K){
        FIM[k,l] <- FIM[l,k] <- sum(V_A[[k]] * V_A[[l]])
      }
      FIM[k,K+1] <- FIM[K+1,k] <- sum(V_A[[k]] * V_)
    }
    FIM[K+1,K+1] <- sum( V_ * V_ )
    FIM <- 0.5 * FIM
    
    #H <- 2 * AI - FIM
    H <- FIM ## Fisher scoring
    
    ## Calculate gradient
    for(k in 1:K){
      dL[k] <- -0.5 * (trV_A[k] - eV_AV_e[k])
    }
    dL[1+K] <- -0.5 * (trV_ - eV2_e)
    
    ## Update parameters
    if(it<=niterEM){ # EM iteration
      dx <- VC*VC*(2*dL)/N
    }else{
      dx <- solve(H,dL)
    }
    VC <- VC + dx
    
    if(minVC>0){
      VC <- pmax(minVC,VC)
    }
    Vg <- VC[1:K]
    Ve <- VC[1+K]
    
    ## Update V
    V <- Ve * Id
    for(k in 1:K){
      V <- V + Vg[k] * G_list[[k]]
    }
    if(eps<tol){
      break
    }
  }
  
  ## Prepare output  
  COV_VC <- solve(FIM)
  SE_VC  <- sqrt(diag(COV_VC))
  Vp     <- sum(Vg[indexVg]) + Ve
  
  var_Vg          <- diag(COV_VC)[indexVg]
  var_Vp          <- sum(COV_VC[c(indexVg,K+1),c(indexVg,K+1)])
  var_Vg_total    <- sum(COV_VC[indexVg,indexVg])
  
  if(length(indexVg)>1){
    cov_Vp_Vg       <- rowSums(COV_VC[indexVg,c(indexVg,K+1)])
  }else{
    cov_Vp_Vg       <- sum(COV_VC[indexVg,c(indexVg,K+1)])
  }
  
  cov_Vp_Vg_total <- sum(COV_VC[indexVg,c(indexVg,K+1)])
  
  ## Heritability and SE for each component
  h2_g     <- Vg[indexVg] / Vp
  #var_h2_g <- (h2_g^2)*( var_Vg/(Vg[indexVg]^2) - 2*cov_Vp_Vg[indexVg]/(Vg[indexVg]*Vp) + var_Vp/(Vp^2) )
  var_h2_g <- (h2_g^2)*( var_Vg/(Vg[indexVg]^2) - 2*cov_Vp_Vg/(Vg[indexVg]*Vp) + var_Vp/(Vp^2) )
  
  ## If intercept included
  colnames(ParamHistory) <- c(paste0("Varcomp",1:K),"Ve")
  ParamHistory <- cbind(Iter=1:nrow(ParamHistory),ParamHistory)
  
  logLikHistory <- cbind(Iter=1:length(logLikHistory),LogLik=logLikHistory)
  results  <- list(h2_g=rbind(Estimate=h2_g,SE=sqrt(var_h2_g)),
                   FixedEffect=fixed_effects,
                   Vg=rbind(Estimate=Vg,SE=SE_VC[1:K]),
                   Ve=c(Estimate=Ve,SE=SE_VC[1+K]),
                   Vp=c(Estimate=Vp,SE=sqrt(var_Vp)),
                   VC=rbind(Estimate=VC,SE=SE_VC),
                   COV_VC=COV_VC,
                   N=N,
                   nparam=(length(fixed_effects) + length(VC)),
                   logLik=logLik,
                   BIC=-2*logLik + log(N) * (length(fixed_effects) + length(VC)),
                   AIC=-2*logLik + 2 * (length(fixed_effects) + length(VC)),
                   logLik_history=logLikHistory,
                   Param_history=ParamHistory)
  return(results)
}

sparseREML_worker_MLC <- function(y,X,G_list,niter=50,tol=1e-6,verbose=TRUE,startValues=NULL,indexVg=NULL){ 
  N  <- length(y)
  Id <- Diagonal(n = N)
  Vy <- var(y)
  
  ## Initialize
  K  <- length(G_list)
  if(is.null(startValues)){
    VC    <- c(rep(1/N,K),Vy*(1-K/N))
    VC    <- rep(Vy/(K+1),K+1)
  }else{
    eVC <- K + 1
    nVC <- length(startValues)
    if(nVC != eVC){
      cat(paste("[Error] The number of input parameters ([`startValues=` option]) is ",
                nVC," and does not match the number of variance components expected (i.e.,",
                eVC,").\n"))
      cat("[Error] Returning NULL")
      return(NULL)
    }else{
      VC <- startValues
    }
  }
  
  #theta <- sqrt(VC) #-log(Vy/VC-1)
  #theta <- -log(Vy/VC-1)
  theta <- acos(sqrt(VC/Vy))
  
  AI  <- matrix(0,nrow=1+K,ncol=1+K)
  FIM <- matrix(0,nrow=1+K,ncol=1+K)
  dL  <- rep(0,1+K)
  eps <- 1

  ## Running AI-REML algorithm
  if(verbose){
    cat("Running Constrained ML algorithm...\n")
  }
  for(it in 1:niter){
    ## Initialize/Update
    Vg <- VC[1:K]
    Ve <- VC[1+K]
    V  <- Ve * Id
    for(k in 1:K){
      V <- V + Vg[k] * G_list[[k]]
    }
    
    ## Step 1 - Optimize wrt fixed effects
    V_     <- solve(V)
    V_y    <- crossprod(V_,y)
    V_X    <- crossprod(V_,X)
    xTV_x  <- crossprod(X,V_X)
    xTV_x_ <- solve(xTV_x)
    xTV_y  <- crossprod(X,V_y)
    fixed_effects <- crossprod(xTV_x_,xTV_y)[,1]
    
    ## Step 2 - Optimize wrt random effects
    e      <- y - X%*%fixed_effects
    V_e    <- crossprod(V_,e)
    V2_e   <- crossprod(V_,V_e)
    
    AV_e    <- lapply(1:K, function(k) crossprod(G_list[[k]],V_e))
    V_AV_e  <- lapply(1:K, function(k) crossprod(V_,AV_e[[k]]) ) # solve(V,APy[[k]]))
    eV_AV_e <- sapply(1:K, function(k) sum(e * V_AV_e[[k]]))
    eV2_e   <- sum(V_e^2)
    
    ## Trace calculations
    V_A    <- lapply(1:K, function(k) crossprod(V_,G_list[[k]]))
    trV_A  <- sapply(1:K, function(k){ sum(diag(V_A[[k]])) })
    trV_   <- sum(diag(V_))
    
    ## Likelihood
    ## Calculate log-likelihood
    logDet_V     <- as.numeric( determinant(V,logarithm = TRUE)$modulus )
    eV_e         <- crossprod(e,V_e)[1,1]
    logLik       <- -0.5 * ( logDet_V + eV_e )

    if(it>1){
      eps <- abs( logLik - prev_logLik )
    }
    prev_logLik <- logLik
    pf <- paste0("[it=",it,"]")
    for(k in 1:K){
      pf <- c(pf,paste0(" - Vg(",k,") = ",round(Vg[k],5)))
    }
    pf <- c(pf,paste0(" - Ve = ",round(Ve,5)," - logLik = ",round(logLik,5)," - eps = ",format(eps,scientific=T,digits=5)))
    if(verbose){
      cat(paste0(paste(pf,collapse = ""),".\n"))
    }
    
    ## Calculate AI matrix
    for(k in 1:K){
      for(l in k:K){
        AI[k,l] <- AI[l,k] <- sum(AV_e[[k]] * V_AV_e[[l]])
      }
      AI[k,K+1] <- AI[K+1,k] <- sum(AV_e[[k]] * V2_e)
    }
    AI[K+1,K+1] <- sum( V_e * V2_e )
    AI          <- 0.5 * AI
    
    ## Fisher-scoring
    for(k in 1:K){
      for(l in k:K){
        FIM[k,l] <- FIM[l,k] <- sum(V_A[[k]] * V_A[[l]])
      }
      FIM[k,K+1] <- FIM[K+1,k] <- sum(V_A[[k]] * V_)
    }
    FIM[K+1,K+1] <- sum( V_ * V_ )
    FIM          <- 0.5 * FIM
    
    ## Hessian Matrix
    H <- 2 * AI - FIM
    
    ## Calculate gradient
    for(k in 1:K){
      dL[k] <- -0.5 * (trV_A[k] - eV_AV_e[k])
    }
    dL[1+K] <- -0.5 * (trV_ - eV2_e)
    
    ## --> Chain rules
    ## VC = Vy * cos(theta)^2
    dt <- -Vy * sin(2*theta)
    dL <- dL * dt 
    H  <- H * (tcrossprod(dt) )
    
    ## Update parameters
    H_    <- ginv(H)
    dx    <- crossprod(H_,dL)
    theta <- theta + dx
    VC    <- Vy * cos(theta) * cos(theta)
    
    if(eps<tol){
      break
    }
  }
  
  ## Prepare Output
  if(is.null(indexVg)){
    indexVg <- 1:K
  }else{
    indexVg <- indexVg
  }
  
  COV_VC <- solve(FIM)
  SE_VC  <- sqrt(diag(COV_VC))
  Vp     <- sum(Vg[indexVg]) + Ve
  
  cat(paste0("length(VC) = ",length(VC)," vs. length(SE_VC) = ",length(SE_VC),".\n"))
  
  var_Vg          <- diag(COV_VC)[indexVg]
  var_Vp          <- sum(COV_VC[c(indexVg,K+1),c(indexVg,K+1)])
  var_Vg_total    <- sum(COV_VC[indexVg,indexVg])
  
  if(length(indexVg)>1){
    cov_Vp_Vg       <- rowSums(COV_VC[indexVg,c(indexVg,K+1)])
  }else{
    cov_Vp_Vg       <- sum(COV_VC[indexVg,c(indexVg,K+1)])
  }
  
  cov_Vp_Vg_total <- sum(COV_VC[indexVg,c(indexVg,K+1)])
  
  ## Heritability and SE for each component
  h2_g     <- Vg[indexVg] / Vp
  #var_h2_g <- (h2_g^2)*( var_Vg/(Vg[indexVg]^2) - 2*cov_Vp_Vg[indexVg]/(Vg[indexVg]*Vp) + var_Vp/(Vp^2) )
  var_h2_g <- (h2_g^2)*( var_Vg/(Vg[indexVg]^2) - 2*cov_Vp_Vg/(Vg[indexVg]*Vp) + var_Vp/(Vp^2) )
  
  ## If intercept included
  results  <- list(h2_g=rbind(Estimate=h2_g,SE=sqrt(var_h2_g)),
                   FixedEffect=fixed_effects,
                   Vg=rbind(Estimate=Vg,SE=SE_VC[1:K]),
                   Ve=c(Estimate=Ve,SE=SE_VC[1+K]),
                   Vp=c(Estimate=Vp,SE=sqrt(var_Vp)),
                   VC=rbind(Estimate=as.numeric(VC),SE=as.numeric(SE_VC)),
                   COV_VC=COV_VC,
                   N=N,
                   nparam=(length(fixed_effects) + length(VC)),
                   logLik=logLik,
                   BIC=-2*logLik + log(N) * (length(fixed_effects) + length(VC)),
                   AIC=-2*logLik + 2 * (length(fixed_effects) + length(VC)))
  return(results)
}


sparseREML <- function(
  prefixGRM="UKB_HM3_m01_sp", ## GCTA fastGWA format
  phenoFile="phenotypes_preAdj_all/HT.phen",
  covarFile="20_pcs.txt",
  model="AE", ## Models are AE, ACE, ACE+, AAE, ACAE
  mpheno=1,
  GRM_range=c(0.05,2.0),
  niter=50,tol=1e-6,
  verbose=TRUE,niterEM=1,
  startValues=NULL,
  RINT=FALSE,
  minVC=0,
  algorithm="AI-REML"){
  
  ## Reading GRM
  if(verbose) cat("Reading GRM...")
  GRM           <- read.table(paste0(prefixGRM,".grm.sp"),stringsAsFactors = F)
  iidsall       <- read.table(paste0(prefixGRM,".grm.id"),stringsAsFactors = F)[,2]
  GRM[,1]       <- iidsall[GRM[,1]+1]
  GRM[,2]       <- iidsall[GRM[,2]+1]
  colnames(GRM) <- c("IID1","IID2","GRM")

  if(TRUE){#model%in%c("AE","ACE","ACE+","AME","AMCE","AMCE+","AM2E","AM2CE")){
    GRM_d <- GRM[+which(GRM[,1]==GRM[,2]),]
    GRM_o <- GRM[-which(GRM[,1]==GRM[,2]),]
    ikeep <- which(GRM_o[,3]>=GRM_range[1] & GRM_o[,3]<GRM_range[2])
    if(length(ikeep)==0){
      cat("Issues with GRM thresholds.\n")
      return(NULL)
    }else{
      GRM <- rbind(GRM_d,GRM_o[ikeep,])
    }
  }else{
    cat("[NOTE: 'GRM_range' can only be used with AE and ACE models.]\n")
  }
    
  if(verbose) cat(" -> [done].\n")
  
  ## Reading phenotype
  if(verbose) cat("Reading phenotype...")
  px <- na.omit(read.table(phenoFile)[,c(2,mpheno+2)]); colnames(px) <- c("IID","Y")
  if(verbose) cat(" -> [done].\n")
  
  if(!is.null(covarFile)){
    if(verbose) cat("Reading phenotype...")
    cx <- na.omit( read.table(covarFile,colClasses = "numeric")[,-1] ); colnames(cx) <- c("IID",paste0("X",1:(ncol(cx)-1)))
    if(verbose) cat(" -> [done].\n")
    mx <- merge(px,cx,by="IID")
  }else{
    mx <- px
  }
  rownames(mx) <- mx[,"IID"]

  if(verbose) cat("Merging data and preparing GRM(s) for analysis...")
  ## Build sparse GRM
  iidnomiss   <- mx[,"IID"]
  n_pheno     <- nrow(mx)
  grm         <- GRM[which(GRM[,1]%in%iidnomiss & GRM[,2]%in%iidnomiss),]
  idsgrm      <- unique(c(as.matrix(grm[,1:2])))
  n_grm       <- length(idsgrm)
  pheno       <- mx[which(mx[,"IID"]%in%idsgrm),]
  N           <- length(idsgrm)
  inds        <- 1:N
  names(inds) <- idsgrm

  ## Main GRM
  i  <- c( inds[as.character(grm[,1])], inds[as.character(grm[,2])] )
  j  <- c( inds[as.character(grm[,2])], inds[as.character(grm[,1])] )
  g  <- c( grm[,3],grm[,3])
  G  <- sparseMatrix(i=i,j=j,x=g) ## Check here if we can do better?
  diag(G) <- diag(G)/2
  colnames(G) <- rownames(G) <- idsgrm
  
  ## Additive by Additive
  Gsq <- G * G
  
  ## AM - ala Kemper et al. (2021)
  gam <- g * log(g) / log(0.5)
  Gam <- sparseMatrix(i=i,j=j,x=gam) ## Check here if we can do better?
  diag(Gam) <- diag(Gam)/2
  colnames(Gam) <- rownames(Gam) <- idsgrm
  # diag(Gam) <- 0
  
  gam_sq <- g * ( log(g) / log(0.5) ) * ( log(g) / log(0.5) - 1 ) * 0.5
  Gam_sq <- sparseMatrix(i=i,j=j,x=gam_sq) ## Check here if we can do better?
  diag(Gam_sq) <- diag(Gam_sq)/2
  colnames(Gam_sq) <- rownames(Gam_sq) <- idsgrm
  # diag(Gam_sq) <- 0
  
  ## Shared Environment
  C <- sparseMatrix(i=i,j=j,x=as.numeric(g>0)) ## Check here if we can do better?
  diag(C) <- diag(C)/2
  colnames(C) <- rownames(C) <- idsgrm
  
  ## Shared Environment for nuclear families only
  CN <- sparseMatrix(i=i,j=j,x=as.numeric(g>=0.35)) ## Check here if we can do better?
  diag(CN) <- diag(CN)/2
  colnames(CN) <- rownames(CN) <- idsgrm
  
  ## Interaction between G and NUC
  GI <- G * CN
  
  ## Phenotype and covariates
  mx <- mx[as.character(idsgrm),]

  if(!sum(table( mx$IID == idsgrm )) == nrow(mx)){
    cat("Merging failed!\n")
    return(NULL)
  }else{
    y <- mx[,"Y"]
    if(RINT){
      y <- qnorm( (rank(scale(y))-0.5)/length(y) )
    }
    if(!is.null(covarFile)){
      X <- as.matrix(cbind(Intercept=1,mx[,colnames(cx)[-1]]))
    }else{
      X <- matrix(1,nrow=N,ncol=1); colnames(X) <- "Intercept"
    }
    if(verbose) cat(" -> [done].\n")
    
    ## Running analysis
    indexVg=NULL
    if(model == "AE")   G_list = list(G)
    if(model == "AEI")  G_list = list(G,GI)
    if(model == "AME"){
      G_list = list(G,Gam)
      indexVg = 1
    }
    if(model == "AMEI"){
      G_list = list(G,Gam,GI)
      indexVg = c(1,3)
    }
    if(model == "ACE")  G_list = list(G,C)
    if(model == "ACE+") G_list = list(G,CN)
    if(model == "ACEI") G_list = list(G,GI,C)
    if(model == "ACEI+") G_list = list(G,GI,CN)
    if(model == "AMCEI+"){
      G_list = list(G,Gam,GI,CN)
      indexVg = c(1,3,4)
    } 
    
    if(model == "AMCE"){
      G_list = list(G,Gam,C)
      indexVg = c(1,3)
    }
    if(model == "AMCE+"){
      G_list = list(G,Gam,CN)
      indexVg = c(1,3)
    }
    if(model == "AMCEI"){
      G_list = list(G,Gam,GI,CN,C)
      indexVg = c(1,3,4)
    }
    if(model == "AM2E"){
      G_list = list(G,Gam,Gam_sq)
      indexVg = 1
    }  
    if(model == "AM2CE"){
      G_list = list(G,Gam,Gam_sq,C)
      indexVg = c(1,4)
    }
    if(model == "AAE")  G_list = list(G,Gsq)
    if(model == "AACE") G_list = list(G,Gsq,C)
    if(model == "AACE+") G_list = list(G,Gsq,CN)
    if(model == "AAMCE"){
      G_list = list(G,Gsq,Gam,C)
      indexVg = c(1,2,4)
    }
    if(model == "AAMCE+"){
      G_list = list(G,Gsq,Gam,CN)
      indexVg = c(1,2,4)
    }
    if(model == "AAE+"){
      G_list = list(G,Gsq,CN)
    }
    
    ## List of model for main analyes
    ## ACEI+ and AMCEI+
    ## AACE+ and AAMCE+
    if(model == "ACEI+") G_list = list(G,GI,CN)
    if(model == "AMCEI+"){
      G_list = list(G,Gam,GI,CN)
      indexVg = c(1,3,4)
    }
    if(model == "AACE+") G_list = list(G,Gsq,CN)
    if(model == "AAMCE"){
      G_list = list(G,Gsq,Gam,C)
      indexVg = c(1,2,4)
    }
    if(model == "AAME"){
      G_list = list(G,Gsq,Gam)
      indexVg = c(1,2)
    }
    if(algorithm=="AI-REML"){
      runModel <- sparseREML_worker_AIREML(y,X,G_list,niter=niter,tol=tol,verbose=verbose,
                                           niterEM=niterEM,startValues=startValues,
                                           indexVg=indexVg,minVC=minVC)
    }
    
    if(algorithm=="ML"){
      runModel <- sparseREML_worker_ML(y,X,G_list,niter=niter,tol=tol,verbose=verbose,
                                       niterEM=niterEM,startValues=startValues,
                                       indexVg=indexVg,minVC=minVC)
    }
    
    if(algorithm=="ML-C"){
      runModel <- sparseREML_worker_MLC(y,X,G_list,niter=niter,tol=tol,verbose=verbose,
                                        startValues=startValues,indexVg=indexVg)
    }
    
    return(runModel)
  }
}

