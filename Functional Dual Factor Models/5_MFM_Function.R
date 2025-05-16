#Input 'X' is a list of FTS
mfm<-function(X){
  ###Step 1: Estimate the front loading
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
  k1<-(min(which(cumsum(percent1) > 0.99)))
  Q1_hat<-R.factor$vectors[,1:k1]
  
  #######Step 2: Estimate the back loading
  dat_center_trans<-list()
  for (m in 1:nrow(X[[1]])){
    dat_center_trans[[m]]<-t(dat.t[[m]]-mean.list)
  }
  M2_hat<-m.mat(dat_center_trans)
  C.factor<-eigen(M2_hat)
  lambda_2<-ifelse(C.factor$values>=0, C.factor$values, 0)
  percent2 <- (lambda_2)/sum(lambda_2)
  k2<-(min(which(cumsum(percent2) > 0.99)))
  Q2_hat<-C.factor$vectors[,1:k2]
  
  #####Step 3: Estimate the Ft
  Ft<-list()
  for (n in 1:nrow(X[[1]])){
    Ft[[n]]<-t(Q1_hat)%*%dat_center[[n]]%*%Q2_hat
  }
  return( list(mu=mean.list, R=Q1_hat, Ft=Ft, C=Q2_hat) )
}

######Forecast function of MFM
mfm.forc <- function(dat, year_horizon, years, ages, n_pop)
{ 
  mfm.dat<-mfm(dat)
  R.dat<-mfm.dat$R
  Ft.dat<-mfm.dat$Ft
  C.dat<-mfm.dat$C
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
  for(ik in 1:n_pop){
    X.forc[[ik]]<-matrix(0, nrow=year_horizon, ncol=length(ages))
    for (t in 1:year_horizon){
      X.forc[[ik]][t,]<- R.dat%*%Ft.forc[[t]]%*%t(C.dat)[,ik] + mfm.dat$mu[,ik]
    }
  }
  mort.forc<-list()
  for (ik in 1:n_pop){
    mort.forc[[ik]]<-matrix(0, year_horizon, length(ages))
    for (ij in 1:year_horizon)
      mort.forc[[ik]][ij,]<-exp(X.forc[[ik]][ij,])
  }
  return(mort.forc)
}










