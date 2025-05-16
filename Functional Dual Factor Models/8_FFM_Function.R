#### Function for high-dimensional functional time series estimation####
####FFM Functions
# Loading the required packages
library(methods)
#library(profvis)
library(LaplacesDemon)
require('fda')
require('pracma')
require('MASS')
require('MTS')
require('PANICr')
require('xts')
require(gdata)
require('doParallel')
require('foreach')
require('RSpectra')
require('tictoc')
require('VAR.etp')
eigs<-RSpectra::eigs
#Definition of the class "MITS" that deals with mixed-nature FTS
mits <- setRefClass('mits',
                    fields = list(
                      indexFunc="vector",
                      indexScalar="vector",
                      indexType="vector",
                      dimFunc="vector",
                      dimCum="vector",
                      nFunc="numeric",
                      nScalar="numeric",
                      nTotal="numeric",
                      nDim="numeric",
                      sampleSize="numeric",
                      func="list",
                      scalar="list",
                      coefs="matrix",
                      inProd="matrix",
                      V="matrix"),
                    methods= list(
                      initialize = function() {
                        indexFunc <<- vector()
                        indexScalar <<- vector()
                        indexType <<- vector()
                        dimFunc <<- vector()
                        dimCum <<- vector()
                        nFunc <<- 0
                        nScalar <<- 0
                        nTotal <<- 0
                        nDim <<- 0 
                        sampleSize <<- 0 
                        func <<- list()
                        scalar <<- list()
                        coefs <<- matrix()
                        inProd <<- matrix()
                        V <<- matrix()
                      },
                      loadFunc = function(funcList,position=(nTotal+1)) {
                        if(is.fd(funcList)) { funcList <- list(funcList)}
                        if(!is.list(funcList)) stop('The first argument is not a list nor a fd')
                        if(position<1 || position > nTotal +1) stop('The position that the user has specified is not in the correct range of pre-existing indexes')
                        size <- length(funcList)
                        saveType <- indexType[position:nTotal]
                        saveFunc <- indexFunc[position:nTotal]
                        saveScalar <- indexScalar[position:nTotal]
                        saveN <- nTotal+1
                        for(i in 1:size){
                          pushFunc(funcList[[i]])
                          indexType[position+i-1] <<- TRUE
                          indexFunc[position+i-1] <<- nFunc
                          indexScalar[position+i-1] <<- 0
                        }
                        if(position < saveN) {
                          indexType <<- c(indexType[1:(position+size-1)],saveType)
                          indexFunc <<- c(indexFunc[1:(position+size-1)],saveFunc)
                          indexScalar <<- c(indexScalar[1:(position+size-1)],saveScalar)
                        }
                        #computeV();computeCoefs();
                      },
                      pushFunc = function(funcA) {
                        if(!is.fd(funcA)) stop('The argument is not of class fd')
                        if(nTotal==0) sampleSize <<- dim(funcA$coefs)[2]
                        if (dim(funcA$coefs)[2] != sampleSize ) stop('The length of the functional time series introduced does not fit with previously introduced data')
                        nFunc <<- nFunc+1
                        nTotal <<- nTotal+1
                        nDim <<- nDim + dim(funcA$coefs)[1]
                        dimFunc[nFunc] <<- dim(funcA$coefs)[1]
                        func[[nFunc]] <<- funcA			
                      },
                      loadScalar = function(scalarList,position=(nTotal+1)) {
                        if(is.vector(scalarList)) { scalarList <- t(as.matrix(scalarList))}
                        if(!is.list(scalarList) && !is.matrix(scalarList)) stop('The first argument is not a matrix nor a ts')
                        if(position<1 || position > nTotal +1) stop('The position that the user has specified is not in the correct range of pre-existing indexes')
                        if(is.list(scalarList)) {
                          size <- length(scalarList)	
                        } else if(is.matrix(scalarList)) {
                          size <- dim(scalarList)[1]
                        }		
                        saveType <- indexType[position:nTotal]
                        saveScalar <- indexScalar[position:nTotal]
                        saveFunc <- indexFunc[position:nTotal]
                        saveN <- nTotal+1
                        for(i in 1:size){
                          if(is.list(scalarList)){
                            pushScalar(scalarList[[i]])
                          } else if(is.matrix(scalarList)) {
                            pushScalar(scalarList[i,])
                          }				
                          indexType[position+i-1] <<- FALSE
                          indexScalar[position+i-1] <<- nScalar
                          indexFunc[position+i-1] <<- 0
                        }
                        if(position < saveN) {
                          indexType <<- c(indexType[1:(position+size-1)],saveType)
                          indexScalar <<- c(indexScalar[1:(position+size-1)],saveScalar)
                          indexFunc <<- c(indexFunc[1:(position+size-1)],saveFunc)
                        }
                        #computeV();computeCoefs();
                      },
                      pushScalar = function(scalarA) {
                        if(!is.vector(scalarA)) stop('The argument is not a vector')
                        if(nTotal==0) sampleSize <<- length(scalarA)
                        if (length(scalarA) != sampleSize ) stop('The length of the time series introduced does not fit with previously introduced data')
                        nScalar <<- nScalar+1
                        nTotal <<- nTotal+1
                        nDim <<- nDim + 1
                        scalar[[nScalar]] <<- scalarA			
                      },
                      computeV = function() {
                        computeDimCum()
                        V <<- matrix(0,nDim,nDim)
                        for (i in 1:nTotal) {
                          if(!indexType[i]) {
                            V[(dimCum[i]+1),(dimCum[i]+1)] <<- 1
                          } else {
                            V[(dimCum[i]+1):(dimCum[i]+dimFunc[indexFunc[i]]),(dimCum[i]+1):(dimCum[i]+dimFunc[indexFunc[i]])] <<- inprod(func[[indexFunc[i]]]$basis,func[[indexFunc[i]]]$basis)
                          }
                        }
                      },
                      computeDimCum = function() {
                        dimCum[1] <<- 0
                        for(i in 2:(nTotal+1)) {
                          if(!indexType[i-1]) {
                            dimCum[i] <<- dimCum[i-1]+1
                          } else {
                            dimCum[i] <<- dimCum[i-1]+ dimFunc[indexFunc[i-1]]
                          }
                        }
                        
                      },
                      computeCoefs = function() {
                        computeDimCum()
                        coefs <<- matrix(0, nDim, sampleSize)
                        for (j in 1:nTotal) {
                          if(!indexType[j]) {
                            coefs[(dimCum[j]+1),] <<- scalar[[indexScalar[j]]]
                          } else {
                            coefs[(dimCum[j]+1):(dimCum[j]+dimFunc[indexFunc[j]]),] <<- func[[indexFunc[j]]]$coefs
                          }
                        }
                      },
                      computeInProd = function() {
                        computeDimCum()
                        inProd <<- matrix(0, nDim, sampleSize)
                        for (j in 1:nTotal) {
                          if(!indexType[j]) {
                            inProd[(dimCum[j]+1),] <<- scalar[[indexScalar[j]]]
                          } else {
                            inProd[(dimCum[j]+1):(dimCum[j]+dimFunc[indexFunc[j]]),] <<- inprod(func[[indexFunc[j]]]$basis,func[[indexFunc[j]]])
                          }
                        }
                      },
                      actualize = function() {
                        computeV()
                        computeInProd()
                        computeCoefs()
                      },
                      trimSampleSize = function(N) {
                        if(nFunc>= 1) {
                          for(i in 1:nFunc) {
                            func[[i]]$coefs <<- func[[i]]$coefs[,1:N]
                            
                          }
                        }
                        if(nScalar>=1) {
                          for(i in 1:nScalar) {
                            scalar[[i]] <<-scalar[[i]][1:N] 
                          }
                        }
                        sampleSize <<- N
                        coefs <<- coefs[,1:N]
                        inProd <<- inProd[,1:N]
                      },
                      shuffleAll =function() {
                        perm <- sample.int(nFunc+nScalar)
                        indexAll <- pmax(indexScalar,indexFunc)
                        indexAllPerm <- indexAll[perm]
                        indexType <<- indexType[perm]
                        indexScalar <<- rep(0,nTotal)
                        indexScalar[indexType==0] <<- indexAllPerm[indexType==0]
                        indexFunc <<- rep(0,nTotal)
                        indexFunc[indexType==1] <<- indexAllPerm[indexType==1]
                        actualize()
                      },
                      reshapeFunc = function(index,trim) {
                        new_dim <- dimFunc[index]-trim
                        if(dimFunc[index]<=trim){
                          stop('number of dimension to remove smaller than number of basis functions of the function of the MTS to reshape')
                        }
                        else {
                          func[[index]]$basis$nbasis <<- new_dim
                          func[[index]]$basis$names <<- func[[1]]$basis$names[1:new_dim]
                          tmp <- func[[index]]$coefs[1:new_dim,]
                          func[[index]]$coefs <<- tmp
                        }
                      },
                      cut = function(n) {
                        indexType <<- indexType[1:n]
                        f <- rep(FALSE,nTotal-n)
                        nTotal <<- n
                        keep_func <- indexFunc[c(indexType,f)]
                        nFunc <<- length(keep_func)
                        indexFunc <<- indexFunc[1:nFunc]
                        func <<- func[sort(keep_func)]
                        indexFunc <<- rank(indexFunc)
                        dimFunc <<-dimFunc[keep_func]
                        if(nScalar>0){
                          keep_scalar <- indexScalar[c(!indexType,f)]
                          nScalar <<- length(keep_scalar)
                          indexScalar <<-indexScalar[1:nScalar]
                          scalar <<-scalar[sort(keep_scalar)]
                          indexScalar <<- rank(indexScalar)
                        }
                        
                        
                        dimCum <<- dimCum[1:(nFunc+1)]
                        nDim <<- dimCum[nFunc+1]
                        
                        coefs <<- coefs[1:nDim,]
                        inProd <<- inProd[1:nDim,]
                        V <<- V[1:nDim,1:nDim]
                      },
                      removeMean = function() {
                        mu <- NULL
                        if(nScalar>=1) {
                          for(i in 1:nScalar) {
                            mum <- mean(scalar[[i]])
                            mu <- c(mu,mum)
                            scalar[[i]] <<- scalar[[i]]-mum
                          }
                        }
                        muf <- mu
                        if(nFunc>= 1) {
                          for(i in 1:nFunc) {
                            mum <- rowMeans(func[[i]]$coefs)
                            mu <- c(mu,mum)
                            func[[i]]$coefs <<- func[[i]]$coefs - mum
                            muf <- c(muf,inprod(func[[i]],func[[i]]$basis))
                          }
                        }
                        #The next two lines assume that the cross-sectional order is such that the scalar have been loaded first then the functions. If it is not the case, the next two lines must be commented and replaced by "actualize()".
                        coefs <<- coefs-mu
                        inProd <<- inProd-muf
                        return(mu)
                      }
                    )
)


normVec <- function (x) {
  # Next function : Euclidian norm.
  # Args:
  #   x : vector
  # Returns:
  #  The Euclian norm of the vector x 
  return(sqrt(sum(x^2)))
} 



is.single_numeric <- function(x){
  # This function returns TRUE when p is a single numeric value
  # Input :
  #x : the value to be checked
  # Output : boolean
  if(is.numeric(x) & length(x)==1){
    return(TRUE)
  } else{
    return(FALSE)
  }
}

findCoefsReduced <- function(ts,p) {
  # This function computes the projection onto the first $p$ principal components of the mixed-nature time series ts 
  # Input :
  #ts : the mixed-nature functional time series (must be of class mits)
  # p : the number of principal components on which to project. 
  #If single value : all FTS should be projected on p principal components
  #Else : p represents the vector of p_i 
  # Output : a list containing
  #ei (length = n_fu, dim(ei[[i]]) = p_i x P) :  a list containing the set of eigenfunction (a vector containing the components of the eigenfunction in the basis of ts$func[[i]])
  #coef (dim = n_p x T) : a matrix containing the stacking of the p scores linked to the first p principal component of every function. 
  
  if(class(ts)!='mits'){
    stop('ts is not of class mits')
  }
  if(!is.numeric(p)){
    stop('p is not numeric')
  }
  tsc <- ts$copy()
  coef <- NULL
  ei <- list()
  n_fu <- tsc$nFunc
  
  if(is.single_numeric(p)){
    p <- rep(p,n_fu)
  }
  
  
  for(i in 1:n_fu) {
    ei[[i]] <- staticComp(tsc$func[[ts$indexFunc[i]]]$coefs,p[i])
    coef <- rbind(coef,crossprod(ei[[i]],tsc$func[[ts$indexFunc[i]]]$coefs))
  }
  return(list(ei,coef))
}

recover <- function(pred,ei) {
  #This function recovers the P-dimensional representation based on its lower-dimensional representation (i.e. the -p-dimension vector of scores) and on the set of eigenfunctions
  #Input:
  #pred (vector of length=n_p or matrix of length n_p x T) : A vector containing the stacking of the coefficients of every function
  #ei : voir output of findCoefsReduced
  #Output : 
  #coef (lenght= n_P) : a vector containing the full-dimensional representation of pred
  if(!is.vector(pred) && !is.matrix(pred)){
    stop('pred is not a vector nor a matrix')
  }
  if(!is.list(ei)){
    stop('ei is not a list')
  }
  
  n_fu <- length(ei)
  P_i <- rep(0,n_fu)
  p_i <- rep(0,n_fu)
  for(i in 1:n_fu) {
    P_i[i] <- dim(ei[[i]])[1]
    p_i[i] <- dim(ei[[i]])[2]
    
  }
  
  if(is.vector(pred)){
    coef <- rep(0,sum(P_i))
    n_i <- 1
    N_i <- 1
    for(i in 1:n_fu) {
      coef[N_i:(N_i+P_i[i]-1)] <- as.vector(crossprod(t(ei[[i]][,,drop=FALSE]),pred[n_i:(n_i+p_i[i]-1)]))
      n_i <- n_i+p_i[i]
      N_i <- N_i+P_i[i]
      
    }
  }
  else if(is.matrix(pred)){
    coef <- matrix(0,sum(P_i),dim(pred)[2])
    n_i <- 1
    N_i <- 1
    for(i in 1:n_fu) {
      coef[N_i:(N_i+P_i[i]-1),] <- crossprod(t(ei[[i]][,,drop=FALSE]),pred[n_i:(n_i+p_i[i]-1),])
      n_i <- n_i+p_i[i]
      N_i <- N_i+P_i[i]
      
    }
  }
  
  return(coef)
}


C_h <- function (X, h) {
  # Estimate the auto-covariance matrix of 'X' at lag 'h' by the usual emperical covariance matrix
  # Input:
  #   X (matrix (nxT) or vector) : the (multiple) time series 
  #   h (single numeric) : the lag
  # Output:
  #   matrix (n x n) the estimator of the auto-covariance matrix
  if(!is.single_numeric(h)){
    stop('h is not a single numeric')
  }
  if(is.vector(X)){
    X <- matrix(X,nrow=1,ncol=length(X))
  }
  p <- dim(X)[1]
  n <- dim(X)[2]
  X <- X - rowMeans(X)
  C <- matrix(0, p, p)
  for (i in 1:(n - h)) {
    #C <- C + X[, i] %*% t(X[, i + h]) #utiliser cross p
    C <- C + tcrossprod(X[,i+h],X[,i])
  }
  return((C) / n)
}

staticComp <- function (ts,r) {
  # This function computes the r first principal components of ts (that is the first r eigenvectors of the covariance matrix of ts)
  # Input : 
  #ts (class mits or matrix) : the function for which we want to compute the principal components. Dimension of the scalar-equivalent representation of the cross-section : N (can be reduced or complete)
  #r (single numeric) : the number of principal components required
  # Output : 
  # (matrix of dimension (N x r) ) : the r first principal components
  if(class(ts)=="mits"){
    coefs <- ts$coefs
  } else if(class(ts)=="matrix"){
    coefs <- ts
  } else {
    stop('ts is neither a MITS nor a matrix')
  }
  if(!is.single_numeric(r)){
    stop('r is not a single numeric')
  }
  
  C <- C_h(coefs,0)
  return(suppressWarnings(eigs(C,r)$vectors))	
}
normi <- function(coefs,V){
  # This function computes the Frobenius norm of coefs
  # Input : coefs(matrix) : contains the coefficients of the mixed-nature FTS
  #-V(matrix) : change of basis matrix, field (attribute) of mits class.
  # Output : Frobenius norm of coefs
  T <- dim(coefs)[2]
  norm2 <- 0
  for(t in 1:T){
    norm2 <- norm2+ t(coefs[,t,drop=FALSE])%*%V%*%coefs[,t,drop=FALSE]
  }
  return(sqrt(norm2))
}

autocov_est <- function(X,S){
  #Compute the matrix containing the auto-covariances of the whole cross-section
  # ts is the time series, S is the maximum lag of the auto-covariance considered, N_theta is the number of frequencies considered in the spectral density estimation, q is the number of dynamic factor and window_size is the length of the window used to estimate the spectral density
  Gamma <- array(dim=c(2*S+1,dim(X)[1],dim(X)[1]))
  
  for(i in 0:S){		
    C <- C_h(X,i)
    Gamma[i+S+1,,] <- C
    Gamma[S+1-i,,] <- t(C)
  }
  return(Gamma)
}

form_BC <- function(Gamma,order,S){
  C <- NULL
  B <- NULL
  index <- 0
  for(k in 1:order){
    B <- cbind(B,Gamma[k+S+1,,])
    CC <- NULL		
    for(j in 1:order){
      CC <- cbind(CC,Gamma[index+S+j,,])
    }
    C <- rbind(C,CC)
    index <- index-1
    
  }
  return(list(B,C))
}
form_AL_i <- function(B,C){
  # Compute the coefficients of the VAR processes  based on the autocovariance matrix at different lags (Gamma) using Yule-Walker equations
  
  AL <- try(B%*%solve(C),silent=TRUE)
  if(!is.matrix(AL)) {
    print('matrix C could not be inverted, replaced by zeroes')
    #AL <- NULL
  } else{ 
    return(AL)
  }
}
form_AL_trimming_i2 <- function(B,C,order,r,sigma_e,T,objectif=1)
{
  # Compute the coefficients of the VAR processes (there are N_ts time series => N_ts VAR processes) based on the autocovariance matrix at different lags (Gamma) using Yule-Walker equations
  tmp <- eigen(C)
  val <- tmp$values
  ev <- tmp$vectors
  
  
  #A_base <- (form_AL_i(B,C))%*%ev[,1,drop=FALSE]%*%t(ev[,1,drop=FALSE])
  A_base <- (form_AL_i(B,C))
  m <- (order)*r
  if(!is.matrix(A_base)){
    return(list(matrix(0,r,m),m))
  } else {
    
    
    
    if(sigma_e==Inf){
      
      return(list(A_base,m))
    } else {
      
      sigma <- rep(0,m)
      sigma1 <- rep(0,m)
      sigma2 <- rep(0,m)
      sigma3 <- rep(0,m)
      
      
      sigma_X <- sum(eigen(C[1:r,1:r]-A_base%*%C%*%t(A_base))$values)
      # cat('sigma_X',sigma_X,'\n')
      for(i in 1:m){
        # P_k <- ev[,1:i]%*%diag(val[1:i])%*%t(ev[,1:i])
        P_k <- ev[,1:i,drop=FALSE]%*%t(ev[,1:i,drop=FALSE])
        sigma2[i] <- sigma_e*sigma_X/(val[i]*T)
        if(objectif==1){
          sigma1[i] <- sigma_e*tr(A_base%*%P_k%*%t(A_base))
        } else {
          sigma1[i] <- sigma_e*tr(A_base%*%(diag(m)-P_k)%*%t(A_base))
        }
        if(i<m){
          if(i==m-1){
            diagmat <- val[m]
          }
          else{
            diagmat <- diag(val[(i+1):m])
          }
          P_k_superscript <- (ev[,(i+1):m,drop=FALSE])%*%diagmat%*%t(ev[,(i+1):m,drop=FALSE])
          sigma3[i] <- tr(A_base%*%P_k_superscript%*%t(A_base))
          sigma[i] <- sigma_X*i/T+sigma1[i]+sum(sigma2[1:i])+sigma3[i]
        }
        else {
          sigma[i] <- sigma_X*i/T+sigma1[i]+sum(sigma2[1:i])
        }
      }
      sigma <- abs(sigma)
      trim <- which(sigma==min(sigma))[1]
      if(is.na(trim)){
        trim <- 1
      }
      ev_final <- ev[,1:trim]
      AL <- try(A_base%*%ev_final%*%t(ev_final),silent=TRUE)
      
      return(list(AL,trim))
    }
  }
} 

stack_ts <- function(F,order) {
  r <- dim(F)[1]
  T <- dim(F)[2]
  Fst <- matrix(0,r*order,T+1-order)
  for(i in 1:order){
    Fst[(1+(i-1)*r):(i*r),] <- F[,(order+1-i):(T+1-i)]
  }
  return(Fst)
}

res_VAR <- function(F,ALi,order,r,T){
  Fst <- stack_ts(F,order)[,1:(T-order),drop=FALSE]
  Fhat <- ALi%*%Fst
  F_lim <- F[,(order+1):T,drop=FALSE]
  
  
  return(F_lim-Fhat)
}

findq_baing <- function(eigs,n,T,delta=0.1,m=0.1){
  l2 <- sum(eigs^2)
  s <- length(eigs)
  D <- rep(0,s)
  lims <-  m/(min(n^(0.5-delta),T^(0.5-delta)))
  
  for(i in 1:(s-1)){
    D[i] <- sqrt(eigs[i+1]^2/l2)
  }
  q <- which(D<lims)[1]
  return(q)
}

est_u <- function(ts,r,type="tilde"){
  Omega <- crossprod(ts$coefs,crossprod(t(ts$V),ts$coefs))
  T <- dim(Omega)[1]
  N <- dim(ts$coefs)[1]
  ei <- eigs(Omega,r)
  if(type=='hat'){
    return(t(tcrossprod(ei$vectors,diag(1/sqrt(N)*ei$values))))
  } else if(type=="tilde"){
    return(t(tcrossprod(ei$vectors,diag(sqrt(T)*rep(1,r)))))
  }
}


est_loadings <- function(ts,u) {
  # This function assumes that the u is the utilde (otherwise problem of scaling)
  T <- dim(ts$coefs)[2]
  # print(str(ts$coefs))
  #print(str(u))
  return(ts$coefs%*%t(u)/T)
}

computeTotalVar <- function(data,V){
  # This function computes the total variance (the sum of the eigenvalues of the empirical covariance operator) contained in a mixed-nature FTS
  # Input : data(matrix) : coefficients of the mixed-nature panel
  #-V(matrix) : change of basis matrix, field (attribute) of mits class.
  # Output : (scalar) the total variance
  T <- dim(data$coefs)[2]
  p <- dim(data$coefs)[1]
  C <- 0
  coefs <- data$coefs-rowMeans(data$coefs)
  for(t in 1:T){
    C <- C+t(coefs[,t,drop=FALSE])%*%V%*%coefs[,t,drop=FALSE] #Estimator of the covariance matrix. Division by T is done at the return to avoid too many calculations
  }
  return(as.numeric(C)/T)
}





###Convert the data into MITS###
convert.dat<- function (dat,standardize=TRUE) {
  # This function convert in a usable format. More details in the code
  # Input : dat: a list of functional time series
  # p (integer): number of basis functions
  # 
  #standardize (boolean) : for every time series : should we remove the mean and divide by the average standard deviation in the same category of data (cf. paper)
  
  #Initializing list that will be used later
  output <- list() 
  input <- list()
  nb_series <- length(dat)
  #ts will contain all the time series together using a mits class
  ts <- mits$new()
  N <- dim(dat[[1]])[2]
  # It also computes the total variance for every FTS. Finally it loads these data info ts, potentially after standardizing every FTS by removing the mean and dividing by the square root of the average variance of every SP100 FTS.
  args <- seq(0, 1, length=N)
  basis <- create.bspline.basis(c(0, 1), nbasis = 101, norder = 4)
  V <- inprod(basis,basis)
  totalvar <- rep(0,nb_series)
  func <- list()
  for(i in 1:nb_series) 
  {
    func[[i]] <- Data2fd(args, t(dat[[i]]), basis)
    totalvar[i] <- computeTotalVar(func[[i]],V)
  }
  
  scale <- sqrt(mean(totalvar))
  for(i in 1:nb_series){
    if(standardize){
      func[[i]]$coefs <- (func[[i]]$coefs-rowMeans(func[[i]]$coefs))/scale
    }
    ts$loadFunc(func[[i]])
  }
  ts$actualize()
  #Important information (ts and others) are put into the list output.
  output$ts <- ts
  output$basis <- basis
  output$N <- N
  output$nb_series <- nb_series
  output$scale<-scale
  return(output)
}

###Input is a list of FTS
ffm<- function(X)
{
  dat.mean<-list()
  for (ik in 1:length(X)){
    dat.mean[[ik]]<-matrix(0, nrow = nrow(X[[1]]), ncol = ncol(X[[1]]))
    for (i in 1:nrow(X[[1]])){
      dat.mean[[ik]][i,]<-colMeans(X[[ik]])
    }
  }
  dat.center<-list()
  for (ik in 1:length(X)){
    dat.center[[ik]]<-X[[ik]]-dat.mean[[ik]]
  }
  HDFFM.mits <- convert.dat(X)
  HDFFM.ts <- HDFFM.mits$ts
  dat.scale<-HDFFM.mits$scale
  #Important values stored in an easily accessible way
  T <- HDFFM.ts$sampleSize
  nDim <- HDFFM.ts$nDim
  n_ts <- HDFFM.ts$nTotal
  
  r <- 20 # number of factors considered. Big r case include small r cases by selectiong appropriate row/columns on B and u.
  u <- est_u(HDFFM.ts,r) #Estimating the factors
  B <- est_loadings(HDFFM.ts,u) #Estimating the factors loadings
  
  #Initializing important arrays, matrix, and vectors.
  chi_est <- array(0,dim=c(nDim,T,r))
  chi_est_cum <- array(0,dim=c(nDim,T,r))
  approx <- matrix(0,n_ts,r)
  sum_approx <- rep(0,n_ts)
  av_approx <- rep(0,r)
  
  
  for(i in 1:r){
    chi_est[,,i] <- B[,i,drop=FALSE]%*%u[i,,drop=FALSE] #Estimation of the common component part due to the i-th factor
    chi_est_cum[,,i] <- B[,1:i,drop=FALSE]%*%u[1:i,,drop=FALSE]#Estimation of the common component if the model contain $i$ factors
    for(n in 1:n_ts){
      dim <- (HDFFM.ts$dimCum[n]+1):HDFFM.ts$dimCum[n+1] #Dimension of the $n$-th time series (scalar or functional) in the cross-section
      chi_est_n <- matrix(chi_est[dim,,i,drop=FALSE],length(dim),T) #Selection of the estimated common component of the $n$-th time series
      ts_n <-HDFFM.ts$coefs[dim,,drop=FALSE] #Selection of the $n$-th functional time series of the true TS
      
      approx[n,i] <- 1-normi(chi_est_n-ts_n,HDFFM.ts$V[dim,dim])^2/normi(ts_n,HDFFM.ts$V[dim,dim])^2 #Estimating the variability explained by the $i$-th factor for the $n$-th time series
    }
    av_approx[i] <- 1-normi(chi_est_cum[,,i]-HDFFM.ts$coefs,HDFFM.ts$V)^2/normi(HDFFM.ts$coefs,HDFFM.ts$V)^2 # Estimating the average explained variance by the $i$-th compnent
  }
  if(max(av_approx>=0.99)){
    r.sel<-which(av_approx>=0.99)[1]
  }else{
    r.sel<-r
  }
  u.sel<-u[1:r.sel,]
  B.sel<-list()
  for (ik in 1:length(X)){
    B.sel[[ik]]<-B[((ik-1)*(ncol(X[[1]]))+1):(ik*ncol(X[[1]])),1:r.sel]
  }
  return( list(mu=dat.mean, u=u.sel, B=B.sel, dat.scale = dat.scale))
}

### Forecast Function for HDFTS
ffm.forc <- function(dat, year_horizon, years, ages, n_pop)
{ 
  ffm.dat<-ffm(dat)
  u.dat<-ffm.dat$u
  B.dat<-ffm.dat$B
  u.forc<-matrix(0, year_horizon, nrow(u.dat))
  for (i in 1:nrow(u.dat)){
    u.forc[,i]<-forecast(auto.arima(u.dat[,i]), h=year_horizon)$mean
  }
  X.forc<-list()
  for (ik in 1:length(dat)){
    X.forc[[ik]]<-(u.forc%*%t(B.dat[[ik]]))*ffm.dat$dat.scale+ffm.dat$mu[[ik]][1:year_horizon,]
  }
  mort.forc<-list()
  for (ik in 1:n_pop){
    mort.forc[[ik]]<-matrix(0, year_horizon, length(ages))
    for (ij in 1:year_horizon)
      mort.forc[[ik]][ij,]<-exp(X.forc[[ik]][ij,])
  }
  
  return(mort.forc)
}














