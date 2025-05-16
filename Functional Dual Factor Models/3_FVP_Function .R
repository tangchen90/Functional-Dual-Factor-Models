#######Estimaiton Procedures
library(fda)
#Input 'X' is a list of FTS
fvp<-function(X){
  ###Step 1: Perform FPCA to the M matrix to obtain Phi
  dat.t<-list()
  for (j in 1:nrow(X[[1]])){
    dat.t[[j]]<-matrix(0, ncol(X[[1]]), length(X))
    for (k in 1:length(X)){
      dat.t[[j]][,k]<-X[[k]][j,]
    }
  }
  mean.list<-Reduce("+", dat.t) / length(dat.t)
  dat_center<-list()
  for (l in 1:nrow(X[[1]])){
    dat_center[[l]]<-dat.t[[l]]-mean.list
  }
  M1_hat<-m.mat(dat_center)
  R.factor<-eigen(M1_hat)
  lambda_1<-ifelse(R.factor$values>=0, R.factor$values, 0)
  percent1 <- (lambda_1)/sum(lambda_1)
  k1<-max((min(which(cumsum(percent1) > 0.95))),2)
  Phi<-R.factor$vectors[,1:k1]
  
  #######Step 2: Estimate Ft
  dat_center_trans<-list()
  for (m in 1:nrow(X[[1]])){
    dat_center_trans[[m]]<-t(dat.t[[m]]-mean.list)
  }
  M2_hat<-m.mat(dat_center_trans)
  C.factor<-eigen(M2_hat)
  lambda_2<-ifelse(C.factor$values>=0, C.factor$values, 0)
  percent2 <- (lambda_2)/sum(lambda_2)
  k2<-max((min(which(cumsum(percent2) > 0.95))),2)
  Q2_hat<-C.factor$vectors[,1:k2]
  Ft<-list()
  for (n in 1:nrow(X[[1]])){
    Ft[[n]]<-t(Phi)%*%dat_center[[n]]%*%Q2_hat
  }
  #### Step 3: Estimate Lambda
  phi.u<-list()
  for (u in 1: ncol(X[[1]])){
    phi.u[[u]]<-rep(0, k1)
    for (i in 1:ncol(Phi)){
      phi.u[[u]][i]<-Phi[,i][u]
    }
  }
  
  Ft.u<-list()
  for (t in 1:length(Ft)){
    Ft.u[[t]]<-matrix(0, ncol(X[[1]]), ncol(Ft[[t]]))
    for (u in 1: ncol(X[[1]])){
      Ft.u[[t]][u,]<-t(as.matrix(phi.u[[u]]))%*%Ft[[t]]
    }
  }
  xlist<-list()
  for (k in 1:ncol(Ft.u[[1]])){
    xlist[[k]]<-matrix(0, ncol(X[[1]]), length(Ft))
    for (t in 1:length(Ft)){
      xlist[[k]][,t]<-Ft.u[[t]][,k]
    }
  }
  lambda<-list()
  for (n in 1:length(X)){
    basis.sel <- create.bspline.basis(c(0, (ncol(X[[1]])-1)), nbasis = ncol(X[[1]]))
    Y.n<-Data2fd(argvals=c(0:(ncol(X[[1]])-1)), t(X[[n]]), basis.sel)
    Xfd<-list()
    for (k in 1:length(xlist)){
      Xfd[[k]]<-Data2fd(argvals=c(0:(ncol(X[[1]])-1)), xlist[[k]], basis.sel)
    }
    betalist<-list()
    for (k in 1:length(xlist)){
      betalist[[k]]<- fdPar(with(Xfd[[k]], fd(basisobj=basis, fdnames=fdnames)))
    }
    fRegressout <- fRegress(Y.n, Xfd, betalist)
    lambda[[n]]<-sapply(fRegressout$betaestlist, coef)
  }
  return( list(mu=mean.list, phi=Phi, Ft=Ft, lambda=lambda) )
}
#### Forecast function of FVP
fvpts <- function(dat, year_horizon, n_pop, ages)
{ 
  fvp.dat<-fvp(dat)
  phi.dat<-fvp.dat$phi
  Ft.dat<-fvp.dat$Ft
  lambda.dat<-fvp.dat$lambda
  Ft.forc<-list()
  for (i in 1:year_horizon){
    Ft.forc[[i]]<-matrix(0, nrow=nrow(Ft.dat[[1]]), ncol = ncol(Ft.dat[[1]]))
  }
  
  for(r in 1: nrow(Ft.dat[[1]])){
    for( k in 1: ncol(Ft.dat[[1]])){
      F.vec<-rep(0, length(Ft.dat))
      for( t in 1: length(Ft.dat)){
        F.vec[t]<-Ft.dat[[t]][r,k]
      }
      F.forc<-forecast(auto.arima(F.vec), h=year_horizon)$mean
      for ( i in 1:year_horizon){
        Ft.forc[[i]][r,k]<-F.forc[i]
      }
    }
  }
  
  X.forc<-list()
  for(iq in 1:n_pop){
    X.forc[[iq]]<-matrix(0, nrow=year_horizon, ncol=ncol(dat[[1]]))
    for (t in 1:year_horizon){
      X.pred<-rep(0,length(ages))
      for (u in 1:ncol(X.forc[[iq]])){
        sum.pred<-0
        for (k in 1:ncol(phi.dat)){
          sum.pred <- sum.pred + phi.dat[u,k]*(Ft.forc[[t]]%*%t(lambda.dat[[iq]]))[k,u]
        }
        X.pred[u]<-sum.pred
      }
      X.forc[[iq]][t,]<-X.pred + fvp.dat$mu[,iq]
    }
  }
  mort.forc<-list()
  for (iq in 1:n_pop){
    mort.forc[[iq]]<-matrix(0, year_horizon, ncol(dat[[1]]))
    for (ij in 1:year_horizon)
      mort.forc[[iq]][ij,]<-exp(X.forc[[iq]][ij,])
  }
  fore_res <-list()
  for (im in 1:n_pop){
    fore_res[[im]]<-t(mort.forc[[im]])
  }
  fun_forc<-do.call(rbind, fore_res)
  return(fun_forc)
}







