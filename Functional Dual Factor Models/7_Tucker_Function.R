library(rTensor)
### MSPE ###
mse<- function(x1, x2){
  mean((x1 - x2) ^ 2)
}


#### Function for high-dimensional functional time series estimation####
###Input is a list of FTS
tucker.tnsr<- function(X, ranks)
{
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
  dat.arr<-array(0, c(nrow(dat_center[[1]]), ncol(dat_center[[1]]), length(dat_center)))
  for (ik in 1:length(dat_center)){
    dat.arr[,,ik]<-dat_center[[ik]]
  }
  tsr.dat<-as.tensor(dat.arr)
  tucker.est <- tucker(tsr.dat, ranks)
  return( list(mu=mean.list, Z = tucker.est$Z, Ulist = tucker.est$U, tsr.dat = tsr.dat) )
}

### Forecast Function for HDFTS
tucker.forc <- function(dat, year_horizon, years, ages, n_pop, ranks)
{ 
  tucker.dat<-tucker.tnsr(dat, ranks)
  U.dat<-tucker.dat$Ulist
  year.mat<-U.dat[[3]]
  year.forc<-matrix(0, year_horizon, ncol = ncol(year.mat))
  for (i in 1:ncol(year.mat)){
    year.forc[,i]<-forecast(Arima(year.mat[,i], order=c(0,1,0)), h=year_horizon)$mean
  }
  U.pred<-U.dat
  U.pred[[3]]<-year.forc
  tsr.pred <- ttl(tucker.dat$Z, U.pred, ms = 1:tucker.dat$tsr.dat@num_modes)
  
  X.forc<-list()
  for (t in 1:length(dat)){
    X.forc[[t]]<-t(tsr.pred@data[,t,] + tucker.dat$mu[,t])
  }
  mort.forc<-list()
  for (ik in 1:n_pop){
    mort.forc[[ik]]<-matrix(0, year_horizon, length(ages))
    for (ij in 1:year_horizon)
      mort.forc[[ik]][ij,]<-exp(X.forc[[ik]][ij,])
  }
  
  
  return(mort.forc)
}














