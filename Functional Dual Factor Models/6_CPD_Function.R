library(rTensor)

####Selecting number of rank-1 tensors
### MSPE ###
mse<- function(x1, x2){
  mean((x1 - x2) ^ 2)
}
########
.superdiagonal_tensor<-function (num_modes, len, elements=1L) 
{
  modes <- rep(len, num_modes)
  arr <- array(0, dim = modes)
  if (length(elements) == 1) 
    elements <- rep(elements, len)
  for (i in 1:len) {
    txt <- paste("arr[", paste(rep("i", num_modes), collapse = ","), 
                 "] <- ", elements[i], sep = "")
    eval(parse(text = txt))
  }
  as.tensor(arr)
}
########
dat.train<-list()
for ( i in 1:47){
  dat.train[[i]]<-extract.ages(extract.years(jpn.data[[i]], years = 1973:2012), 0:100, combine.upper=TRUE)$rate$female
}

dat.test<-list()
for ( i in 1:47){
  dat.test[[i]]<-extract.ages(extract.years(jpn.data[[i]], years = 2013:2022), 0:100, combine.upper=TRUE)$rate$female
}

train.omit<-list()
for (ik in 1:47){
  train.omit[[ik]]<-matrix(0, nrow(dat.train[[1]]),ncol(dat.train[[1]]) )
  for(i in 1:nrow(dat.train[[1]])){
    for (j in 1:ncol(dat.train[[1]])){
      if (dat.train[[ik]][i,j]==0){
        train.omit[[ik]][i,j]<-0.0000001
      }else{
        train.omit[[ik]][i,j]<-dat.train[[ik]][i,j]
      }
    }
  }
}

test.omit<-list()
for (ik in 1:47){
  test.omit[[ik]]<-matrix(0, nrow(dat.test[[1]]),ncol(dat.test[[1]]) )
  for(i in 1:nrow(dat.test[[1]])){
    for (j in 1:ncol(dat.test[[1]])){
      if (dat.test[[ik]][i,j]==0){
        test.omit[[ik]][i,j]<-0.0000001
      }else{
        test.omit[[ik]][i,j]<-dat.test[[ik]][i,j]
      }
    }
  }
}

data.train<-list()
for (ik in 1:47){
  data.train[[ik]]<-t(log(train.omit[[ik]]))
}
data.test<-list()
for (ik in 1:47){
  data.test[[ik]]<-t(log(test.omit[[ik]]))
}

dat.t<-list()
for (j in 1:nrow(data.train[[1]])){
  dat.t[[j]]<-matrix(0, ncol(data.train[[1]]), length(data.train))
  for (k in 1:length(data.train)){
    dat.t[[j]][,k]<-data.train[[k]][j,]
  }
}
mean.list<-Reduce("+", dat.t) / length(dat.t)
dat_center<-list()
for (l in 1:nrow(data.train[[1]])){
  dat_center[[l]]<-dat.t[[l]]-mean.list
}

dat.arr<-array(0, c(nrow(dat_center[[1]]), ncol(dat_center[[1]]), length(dat_center)))
for (ik in 1:length(dat_center)){
  dat.arr[,,ik]<-dat_center[[ik]]
}
test.dat<-list()
for (j in 1:nrow(data.test[[1]])){
  test.dat[[j]]<-matrix(0, ncol(data.test[[1]]), length(data.test))
  for (k in 1:length(data.test)){
    test.dat[[j]][,k]<-data.test[[k]][j,]
  }
}
mean.list<-Reduce("+", test.dat) / length(test.dat)
test_center<-list()
for (l in 1:nrow(data.test[[1]])){
  test_center[[l]]<-test.dat[[l]]-mean.list
}

test.arr<-array(0, c(nrow(test_center[[1]]), ncol(test_center[[1]]), length(test_center)))
for (ik in 1:length(test_center)){
  test.arr[,,ik]<-test_center[[ik]]
}
mse.pred<-rep(0,10)
tsr.dat<-as.tensor(dat.arr)
for (k in 1:10){
  cp.tnsr<-cp(tsr.dat, k)
  t.tnsr<-cp.tnsr$U[[3]]
  t.forc<-matrix(0, 10, ncol = ncol(t.tnsr))
  for (i in 1:ncol(t.tnsr)){
    t.forc[,i]<-forecast(Arima(t.tnsr[,i], order=c(0,1,0)), h=10)$mean
  }
  U.pred<-cp.tnsr$U
  U.pred[[3]]<-t.forc
  Z <-.superdiagonal_tensor(num_modes = tsr.dat@num_modes, len = k, elements = cp.tnsr$lambdas)
  tsr.pred <- ttl(Z, U.pred, ms = 1:tsr.dat@num_modes)
  mse.pred[k]<-mse(tsr.pred@data, test.arr)
}


#### Function for high-dimensional functional time series estimation####
###Input is a list of FTS
cp.tnsr<- function(X, K)
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
  cp.est <- cp(tsr.dat, K)
  return( list(mu=mean.list, lambda = cp.est$lambdas, Ulist = cp.est$U, tsr.dat = tsr.dat) )
}

### Forecast Function for HDFTS
cp.forc <- function(dat, year_horizon, years, ages, n_pop, K)
{ 
  cp.dat<-cp.tnsr(dat, K)
  U.dat<-cp.dat$Ulist
  year.mat<-U.dat[[3]]
  year.forc<-matrix(0, year_horizon, ncol = ncol(year.mat))
  for (i in 1:ncol(year.mat)){
    year.forc[,i]<-forecast(Arima(year.mat[,i], order=c(0,1,0)), h=year_horizon)$mean
  }
  U.pred<-U.dat
  U.pred[[3]]<-year.forc
  Z <-.superdiagonal_tensor(num_modes = cp.dat$tsr.dat@num_modes, len = K, elements = cp.dat$lambda)
  tsr.pred <- ttl(Z, U.pred, ms = 1:cp.dat$tsr.dat@num_modes)

  X.forc<-list()
  for (t in 1:length(dat)){
    X.forc[[t]]<-t(tsr.pred@data[,t,] + cp.dat$mu[,t])
  }
  mort.forc<-list()
  for (ik in 1:n_pop){
    mort.forc[[ik]]<-matrix(0, year_horizon, length(ages))
    for (ij in 1:year_horizon)
      mort.forc[[ik]][ij,]<-exp(X.forc[[ik]][ij,])
  }
  
  
  return(mort.forc)
}














