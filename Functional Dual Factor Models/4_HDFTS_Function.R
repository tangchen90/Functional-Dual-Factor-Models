#### Function for high-dimensional functional time series estimation####
###Input is a list of FTS
hdfts<- function(X, years)
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
  dfpca.vector<-list()
  num.component<-rep(0,length(X))
  for (ik in 1:length(X)){
    T=length(years)
    Gamma_l <- long_run_covariance_estimation(t(dat.center[[ik]]))
    FPCA.dym <- eigen(Gamma_l)
    dfpca.value <- FPCA.dym$values
    dfpca.value <- ifelse(dfpca.value>=0, dfpca.value, 0)
    dfpca.vector[[ik]]<-FPCA.dym$vectors
    percent_var <- (dfpca.value)/sum(dfpca.value)
    ratio <- dfpca.value[1]/dfpca.value
    num.component[ik] <-  max(min(which(cumsum(percent_var) > 0.99)), min(which(ratio>(sqrt(T)/log10(T)))-1))
  }
  
  K <- max(num.component)
  dfpca.score<-list()
  for(ik in 1:length(X)){
    dfpca.score[[ik]]<-dat.center[[ik]]%*%dfpca.vector[[ik]][,1:K]
  }
  phi<-list()
  for (ik in 1:length(X)){
    phi[[ik]]<-dfpca.vector[[ik]][,1:K]
  }
  
  
  beta_t<-array(0,dim=c(nrow(dfpca.score[[1]]), length(X), K))
  for (i in 1:K){
    for(j in 1:length(X)){
      beta_t[,,i][,j]<-dfpca.score[[j]][,K]
    }
  }
  
  A_p<-list()
  Ft<-list()
  for (l in 1:K){
    beta.hat<-beta_t[,,l]
    M_hat<-vfm(beta.hat)
    A.factor<-eigen(M_hat)
    lambda<-ifelse(A.factor$values>=0, A.factor$values, 0)
    percent <- lambda/sum(lambda)
    k1<-(min(which(cumsum(percent) > 0.99)))
    A_hat<-A.factor$vectors[,1:k1]
    A_p[[l]]<-A_hat
    F_t<-matrix(0, nrow = nrow(beta.hat), ncol = k1)
    for ( i in 1: nrow(beta.hat)){
      F_t[i,]<-t(A_hat)%*%beta.hat[i,]
    }
    Ft[[l]]<-F_t
  }
  return( list(mu=dat.mean, phi=phi, beta = beta_t, A=A_p, Ft=Ft) )
}

### Forecast Function for HDFTS
hdfts.forc <- function(dat, year_horizon, years, ages, n_pop)
{ 
  hdfts.dat<-hdfts(dat, years)
  Ft.dat<-hdfts.dat$Ft
  beta.dat<-hdfts.dat$beta
  A.dat<-hdfts.dat$A
  beta.forc<-list()
  for (l in 1:length(Ft.dat)){
    F_t <- Ft.dat[[l]]
    F.forc <- matrix(0, nrow = year_horizon, ncol=ncol(F_t))
    for ( j in 1:ncol(F_t)){
      F.forc[,j]<-forecast(auto.arima(F_t[,j]), h=year_horizon)$mean
    }
    beta.forc[[l]]<-matrix(0, nrow =year_horizon, ncol = length(dat) )
    for (k in 1:year_horizon){
      beta.forc[[l]][k,]<-A.dat[[l]]%*%F.forc[k,]
    }
  }
  X.forc<-list()
  for (t in 1:length(dat)){
    X.forc[[t]]<-hdfts.dat$mu[[t]][1:year_horizon,]
    for (s in 1:year_horizon){
      for (m in 1:length(Ft.dat)){
        X.forc[[t]][s,]<-X.forc[[t]][s,]+ beta.forc[[m]][s,t]*hdfts.dat$phi[[t]][,m]
      }
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














