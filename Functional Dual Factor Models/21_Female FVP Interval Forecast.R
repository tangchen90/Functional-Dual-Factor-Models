library(fda)
#############Female###

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


###########################################################
# constructing multivariate pointwise prediction intervals
###########################################################

# data_series: specific data series
# fh: forecast horizon
# nboot: number of bootstrap replication
# alpha: nominal coverage probability

#################################################
############## bootstrap error function##########
#################################################
boot_female <-function(err, nboot)
{
  boots <-matrix(, nrow(err), nboot)
  for(ij in 1:nboot)
  {
    boots [,ij] = err[, sample(1:ncol(err), 1, replace = TRUE)]
  }
  return(boots)
} 

######################################
#### Interval estimation Function ####
######################################
find_enlarge_val_fvp <- function(data_series, fh , nboot = 1000, alpha = 0.8, n_pop)
{
  n_test <- 25
  
  # calculate  in-sample forecast curves
  
  fore_curve_female  = matrix( , length(data_series[[1]]$age)*n_pop, length(data_series[[1]]$year) - n_test - fh )
  for(ij in 1:(length(data_series[[1]]$year) - n_test - fh))
  {
    dat.rate.female_h_step<-list()
    for (ik in 1:n_pop){
      dat.rate.female_h_step[[ik]]<-t(log(extract.years(data_series[[ik]], years = 1973:(1973 + n_test + ij - 1))$rate$female))
    }
    fore_res = fvpts(dat=dat.rate.female_h_step, year_horizon = fh, n_pop = length(data_series), ages = 0:100)
    # fill in gender specific in-sample forecast curves
    fore_curve_female[,ij] = log(fore_res[,fh])
  }
  
  # holdout data samples
  true_dat = list()
  for (im in 1: n_pop){
    true_dat[[im]] = as.data.frame(log(data_series[[im]]$rate$female[, (n_test + fh + 1):length(data_series[[im]]$year)]))
  }
  
  
  
  # female sample
  holdout_val_female = do.call(rbind, true_dat)
  
  
  err_female = holdout_val_female - fore_curve_female
  
  female.err<-list()
  for (il in 1: n_pop){
    female.err[[il]]<-err_female[(1:101)+(il-1)*101 ,]
  }
  female.boot_err<-list()
  for (il in 1: n_pop){
    female.boot_err[[il]]<-boot_female(female.err[[il]], 1000)
  }
  
  # constructing PI
  dat.series<-list()
  for (ik in 1:n_pop){
    dat.series[[ik]]<-t(log(data_series[[ik]]$rate$female))
  }
  
  
  fore_res  = fvpts(dat=dat.series, year_horizon = fh, n_pop = length(dat.series), ages = 0:100)
  
  fore_fvp_female <- list()
  for (il in 1: n_pop){
    fore_fvp_female[[il]]<-matrix(rep(log(fore_res[(1:101)+(il-1)*101 ,fh]), nboot), nrow = 101, ncol = nboot)
  }
  
  
  
  boot_PI_female <- list()
  for (il in 1: n_pop){
    boot_PI_female[[il]]<- female.boot_err[[il]] + fore_fvp_female[[il]]
  }
  
  boot_PI_lb_ub_female <- list()
  for (il in 1: n_pop){
    boot_PI_lb_ub_female[[il]]<- apply(boot_PI_female[[il]], 1, quantile, c((1 - alpha)/2, (1 + alpha)/2))
  }
  
  
  return(list(boot_PI_lb_ub_female = boot_PI_lb_ub_female, boot_sample_female = boot_PI_female))
}    



###############################################################################
# compute interval scores  for a particular forecast horizon 
###############################################################################

# PI_val: one- to 10-step-ahead prediction intervals
# alpha: nominal coverage probability

interval_score_fvp_female <- function(PI_val, data_series, fh, alpha = 0.8)
{
  
  test_val_female<- extract.years(data_series, (2012+fh):2022)$rate$female
  
  # transform back to the original scale
  
  boot_sample_female = exp(PI_val)
  boot_index_female = which(boot_sample_female > 1)
  boot_index_below_female = which(boot_sample_female < 0)
  if(length(boot_index_female) > 0)
  {
    boot_sample_v1_female = replace(boot_sample_female, boot_index_female, 1)
  } else 
  {
    boot_sample_v1_female = boot_sample_female
  }
  if(length(boot_index_below_female) > 0)
  {
    boot_sample_v2_female = replace(boot_sample_v1_female, boot_index_below_female, 0)
  } else
  {
    boot_sample_v2_female = boot_sample_v1_female
  }
  
  # lower and upper bounds  
  # female series
  
  dummy_female = apply(boot_sample_female, c(1,3), quantile, c((1 - alpha)/2, (1 + alpha)/2), na.rm = TRUE)
  PI_lb_val_female = dummy_female[1,,]
  PI_ub_val_female  = dummy_female[2,,]
  
  lb_ind_female  = ifelse(test_val_female < PI_lb_val_female, 1, 0)
  ub_ind_female  = ifelse(test_val_female > PI_ub_val_female, 1, 0)
  score_female  = mean((PI_ub_val_female - PI_lb_val_female) + 2/(1 - alpha) * (PI_lb_val_female - test_val_female) * lb_ind_female + 
                         2/(1 - alpha) * (test_val_female - PI_ub_val_female) * ub_ind_female)
  
  return(list(score_female = score_female))
}



dat.female<-list()
for ( i in 1:length(jpn.data)){
  dat.female[[i]]<-smooth.demogdata(extract.ages(extract.years(jpn.data[[i]], years = 1973:2022), 0:100, combine.upper=TRUE))
}

####### Prediction Intervals Constructions #######
PI_fh_fvp<- function(data_series, fh, nboot = 1000)
{
  PI_boot_fvp_Hokkaido = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Aomori = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Iwate = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Miyagi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Akita  = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Yamagata = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Fukushima = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Ibaraki = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Tochigi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Gunma = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Saitama = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Chiba = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Tokyo = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Kanagawa = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Niigata = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Toyama = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Ishikawa = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Fukui = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Yamanashi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Nagano = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Gifu = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Shizuoka = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Aichi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Mie = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Shiga = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Kyoto = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Osaka = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Hyogo = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Nara = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Wakayama = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Tottori = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Shimane = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Okayama = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Hiroshima = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Yamaguchi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Tokushima = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Kagawa = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Ehime  = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Kochi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Fukuoka = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Saga = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Nagasaki = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Kumamoto = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Oita = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Miyazaki = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Kagoshima = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Okinawa = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  
  for(io in 1:(11-fh))
  {
    dat.train<-list()
    for ( i in 1:length(data_series)){
      dat.train[[i]]<-extract.years(data_series[[i]], years = 1973:(2012+io))
    }
    dum_pointwise = find_enlarge_val_fvp(data_series = dat.train, fh = fh, n_pop = length(data_series))
    PI_boot_fvp_Hokkaido[,,io] = dum_pointwise$boot_sample_female[[1]]
    PI_boot_fvp_Aomori[,,io] = dum_pointwise$boot_sample_female[[2]]
    PI_boot_fvp_Iwate[,,io] = dum_pointwise$boot_sample_female[[3]]
    PI_boot_fvp_Miyagi[,,io] = dum_pointwise$boot_sample_female[[4]]
    PI_boot_fvp_Akita [,,io] = dum_pointwise$boot_sample_female[[5]]
    PI_boot_fvp_Yamagata[,,io] = dum_pointwise$boot_sample_female[[6]]
    PI_boot_fvp_Fukushima[,,io] = dum_pointwise$boot_sample_female[[7]]
    PI_boot_fvp_Ibaraki[,,io] = dum_pointwise$boot_sample_female[[8]]
    PI_boot_fvp_Tochigi[,,io] = dum_pointwise$boot_sample_female[[9]]
    PI_boot_fvp_Gunma[,,io] = dum_pointwise$boot_sample_female[[10]]
    PI_boot_fvp_Saitama[,,io] = dum_pointwise$boot_sample_female[[11]]
    PI_boot_fvp_Chiba[,,io] = dum_pointwise$boot_sample_female[[12]]
    PI_boot_fvp_Tokyo[,,io] = dum_pointwise$boot_sample_female[[13]]
    PI_boot_fvp_Kanagawa[,,io] = dum_pointwise$boot_sample_female[[14]]
    PI_boot_fvp_Niigata[,,io] = dum_pointwise$boot_sample_female[[15]]
    PI_boot_fvp_Toyama[,,io] = dum_pointwise$boot_sample_female[[16]]
    PI_boot_fvp_Ishikawa[,,io] = dum_pointwise$boot_sample_female[[17]]
    PI_boot_fvp_Fukui[,,io] = dum_pointwise$boot_sample_female[[18]]
    PI_boot_fvp_Yamanashi[,,io] = dum_pointwise$boot_sample_female[[19]]
    PI_boot_fvp_Nagano[,,io] = dum_pointwise$boot_sample_female[[20]]
    PI_boot_fvp_Gifu[,,io] = dum_pointwise$boot_sample_female[[21]]
    PI_boot_fvp_Shizuoka[,,io] = dum_pointwise$boot_sample_female[[22]]
    PI_boot_fvp_Aichi[,,io] = dum_pointwise$boot_sample_female[[23]]
    PI_boot_fvp_Mie[,,io] = dum_pointwise$boot_sample_female[[24]]
    PI_boot_fvp_Shiga[,,io] = dum_pointwise$boot_sample_female[[25]]
    PI_boot_fvp_Kyoto[,,io] = dum_pointwise$boot_sample_female[[26]]
    PI_boot_fvp_Osaka[,,io] = dum_pointwise$boot_sample_female[[27]]
    PI_boot_fvp_Hyogo[,,io] = dum_pointwise$boot_sample_female[[28]]
    PI_boot_fvp_Nara[,,io] = dum_pointwise$boot_sample_female[[29]]
    PI_boot_fvp_Wakayama[,,io] = dum_pointwise$boot_sample_female[[30]]
    PI_boot_fvp_Tottori[,,io] = dum_pointwise$boot_sample_female[[31]]
    PI_boot_fvp_Shimane[,,io] = dum_pointwise$boot_sample_female[[32]]
    PI_boot_fvp_Okayama[,,io] = dum_pointwise$boot_sample_female[[33]]
    PI_boot_fvp_Hiroshima[,,io] = dum_pointwise$boot_sample_female[[34]]
    PI_boot_fvp_Yamaguchi[,,io] = dum_pointwise$boot_sample_female[[35]]
    PI_boot_fvp_Tokushima[,,io] = dum_pointwise$boot_sample_female[[36]]
    PI_boot_fvp_Kagawa[,,io] = dum_pointwise$boot_sample_female[[37]]
    PI_boot_fvp_Ehime [,,io] = dum_pointwise$boot_sample_female[[38]]
    PI_boot_fvp_Kochi[,,io] = dum_pointwise$boot_sample_female[[39]]
    PI_boot_fvp_Fukuoka[,,io] = dum_pointwise$boot_sample_female[[40]]
    PI_boot_fvp_Saga[,,io] = dum_pointwise$boot_sample_female[[41]]
    PI_boot_fvp_Nagasaki[,,io] = dum_pointwise$boot_sample_female[[42]]
    PI_boot_fvp_Kumamoto[,,io] = dum_pointwise$boot_sample_female[[43]]
    PI_boot_fvp_Oita[,,io] = dum_pointwise$boot_sample_female[[44]]
    PI_boot_fvp_Miyazaki[,,io] = dum_pointwise$boot_sample_female[[45]]
    PI_boot_fvp_Kagoshima[,,io] = dum_pointwise$boot_sample_female[[46]]
    PI_boot_fvp_Okinawa[,,io] = dum_pointwise$boot_sample_female[[47]]
  }
  return(list(PI_boot_fvp_Hokkaido = PI_boot_fvp_Hokkaido,PI_boot_fvp_Aomori = PI_boot_fvp_Aomori,
              PI_boot_fvp_Iwate = PI_boot_fvp_Iwate,PI_boot_fvp_Miyagi = PI_boot_fvp_Miyagi,
              PI_boot_fvp_Akita  = PI_boot_fvp_Akita ,PI_boot_fvp_Yamagata = PI_boot_fvp_Yamagata,
              PI_boot_fvp_Fukushima = PI_boot_fvp_Fukushima, PI_boot_fvp_Ibaraki = PI_boot_fvp_Ibaraki,
              PI_boot_fvp_Tochigi = PI_boot_fvp_Tochigi, PI_boot_fvp_Gunma = PI_boot_fvp_Gunma,
              PI_boot_fvp_Saitama = PI_boot_fvp_Saitama, PI_boot_fvp_Chiba = PI_boot_fvp_Chiba, 
              PI_boot_fvp_Tokyo = PI_boot_fvp_Tokyo, PI_boot_fvp_Kanagawa = PI_boot_fvp_Kanagawa,
              PI_boot_fvp_Niigata = PI_boot_fvp_Niigata, PI_boot_fvp_Toyama = PI_boot_fvp_Toyama,
              PI_boot_fvp_Ishikawa = PI_boot_fvp_Ishikawa,PI_boot_fvp_Fukui = PI_boot_fvp_Fukui, 
              PI_boot_fvp_Yamanashi = PI_boot_fvp_Yamanashi,PI_boot_fvp_Nagano = PI_boot_fvp_Nagano,
              PI_boot_fvp_Gifu = PI_boot_fvp_Gifu,PI_boot_fvp_Shizuoka = PI_boot_fvp_Shizuoka,
              PI_boot_fvp_Aichi = PI_boot_fvp_Aichi,PI_boot_fvp_Mie = PI_boot_fvp_Mie,
              PI_boot_fvp_Shiga = PI_boot_fvp_Shiga,PI_boot_fvp_Kyoto = PI_boot_fvp_Kyoto,
              PI_boot_fvp_Osaka = PI_boot_fvp_Osaka, PI_boot_fvp_Hyogo = PI_boot_fvp_Hyogo,
              PI_boot_fvp_Nara = PI_boot_fvp_Nara, PI_boot_fvp_Wakayama = PI_boot_fvp_Wakayama, 
              PI_boot_fvp_Tottori = PI_boot_fvp_Tottori,PI_boot_fvp_Shimane = PI_boot_fvp_Shimane,
              PI_boot_fvp_Okayama = PI_boot_fvp_Okayama,PI_boot_fvp_Hiroshima = PI_boot_fvp_Hiroshima,
              PI_boot_fvp_Yamaguchi = PI_boot_fvp_Yamaguchi, PI_boot_fvp_Tokushima = PI_boot_fvp_Tokushima,
              PI_boot_fvp_Kagawa = PI_boot_fvp_Kagawa, PI_boot_fvp_Ehime  = PI_boot_fvp_Ehime, 
              PI_boot_fvp_Kochi = PI_boot_fvp_Kochi, PI_boot_fvp_Fukuoka = PI_boot_fvp_Fukuoka,
              PI_boot_fvp_Saga = PI_boot_fvp_Saga,PI_boot_fvp_Nagasaki = PI_boot_fvp_Nagasaki,
              PI_boot_fvp_Kumamoto = PI_boot_fvp_Kumamoto,PI_boot_fvp_Oita = PI_boot_fvp_Oita,
              PI_boot_fvp_Miyazaki = PI_boot_fvp_Miyazaki,PI_boot_fvp_Kagoshima = PI_boot_fvp_Kagoshima,
              PI_boot_fvp_Okinawa = PI_boot_fvp_Okinawa))
}
library(doParallel)
cl <- makeCluster(10) 
registerDoParallel(cl)

FVP_fh_PI = foreach(mm = 1:10, .packages = c("demography", "fda")) %dopar% PI_fh_fvp(data_series = dat.female,  fh = mm)

save.image("14-3-2024.RData")

female_score.fvp.Hokkaido = female_score.fvp.Aomori = female_score.fvp.Iwate = female_score.fvp.Miyagi = 
  female_score.fvp.Akita = female_score.fvp.Yamagata = female_score.fvp.Fukushima = female_score.fvp.Ibaraki = 
  female_score.fvp.Tochigi = female_score.fvp.Gunma = female_score.fvp.Saitama= female_score.fvp.Chiba = 
  female_score.fvp.Tokyo= female_score.fvp.Kanagawa = female_score.fvp.Niigata = female_score.fvp.Toyama = 
  female_score.fvp.Ishikawa = female_score.fvp.Fukui = female_score.fvp.Yamanashi = female_score.fvp.Nagano = 
  female_score.fvp.Gifu = female_score.fvp.Shizuoka = female_score.fvp.Aichi = female_score.fvp.Mie = 
  female_score.fvp.Shiga = female_score.fvp.Kyoto = female_score.fvp.Osaka = female_score.fvp.Hyogo = 
  female_score.fvp.Nara = female_score.fvp.Wakayama = female_score.fvp.Tottori = female_score.fvp.Shimane = 
  female_score.fvp.Okayama = female_score.fvp.Hiroshima = female_score.fvp.Yamaguchi = female_score.fvp.Tokushima = 
  female_score.fvp.Kagawa = female_score.fvp.Ehime = female_score.fvp.Kochi = female_score.fvp.Fukuoka = 
  female_score.fvp.Saga = female_score.fvp.Nagasaki = female_score.fvp.Kumamoto = female_score.fvp.Oita =   
  female_score.fvp.Miyazaki = female_score.fvp.Kagoshima = female_score.fvp.Okinawa =vector( ,year_horizon)

year_horizon = 10

for(ik in 1:year_horizon)
{
  # me 
  #### Hokkaido###
  female_score.fvp.Hokkaido[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Hokkaido,
                                                     data_series = dat.female[[1]], fh = ik)$score_female
  
  
  #### Aomori###
  female_score.fvp.Aomori[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Aomori,
                                                   data_series = dat.female[[2]], fh = ik)$score_female
  #### Iwate###
  female_score.fvp.Iwate[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Iwate,
                                                  data_series = dat.female[[3]], fh = ik)$score_female
  #### Miyagi###
  female_score.fvp.Miyagi[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Miyagi,
                                                   data_series = dat.female[[4]], fh = ik)$score_female
  #### Akita###
  female_score.fvp.Akita[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Akita,
                                                  data_series = dat.female[[5]], fh = ik)$score_female
  #### Yamagata###
  female_score.fvp.Yamagata[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Yamagata,
                                                     data_series = dat.female[[6]], fh = ik)$score_female
  #### Fukushima###
  female_score.fvp.Fukushima[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Fukushima,
                                                      data_series = dat.female[[7]], fh = ik)$score_female
  #### Ibaraki###
  female_score.fvp.Ibaraki[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Ibaraki,
                                                    data_series = dat.female[[8]], fh = ik)$score_female
  #### Tochigi###
  female_score.fvp.Tochigi[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Tochigi,
                                                    data_series = dat.female[[9]], fh = ik)$score_female
  #### Gunma###
  female_score.fvp.Gunma[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Gunma,
                                                  data_series = dat.female[[11]], fh = ik)$score_female
  #### Saitama###
  female_score.fvp.Saitama[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Saitama,
                                                    data_series = dat.female[[12]], fh = ik)$score_female
  #### Chiba###
  female_score.fvp.Chiba[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Chiba,
                                                  data_series = dat.female[[13]], fh = ik)$score_female
  #### Tokyo###
  female_score.fvp.Tokyo[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Tokyo,
                                                  data_series = dat.female[[14]], fh = ik)$score_female
  #### Kanagawa###
  female_score.fvp.Kanagawa[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Kanagawa,
                                                     data_series = dat.female[[1]], fh = ik)$score_female
  #### Niigata###
  female_score.fvp.Niigata[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Niigata,
                                                    data_series = dat.female[[15]], fh = ik)$score_female
  #### Toyama###
  female_score.fvp.Toyama[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Toyama,
                                                   data_series = dat.female[[16]], fh = ik)$score_female
  #### Ishikawa###
  female_score.fvp.Ishikawa[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Ishikawa,
                                                     data_series = dat.female[[17]], fh = ik)$score_female
  #### Fukui###
  female_score.fvp.Fukui[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Fukui,
                                                  data_series = dat.female[[18]], fh = ik)$score_female
  #### Yamanashi###
  female_score.fvp.Yamanashi[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Yamanashi,
                                                      data_series = dat.female[[19]], fh = ik)$score_female
  #### Nagano###
  female_score.fvp.Nagano[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Nagano,
                                                   data_series = dat.female[[20]], fh = ik)$score_female
  #### Gifu###
  female_score.fvp.Gifu[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Gifu,
                                                 data_series = dat.female[[21]], fh = ik)$score_female
  #### Shizuoka###
  female_score.fvp.Shizuoka[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Shizuoka,
                                                     data_series = dat.female[[22]], fh = ik)$score_female
  #### Aichi###
  female_score.fvp.Aichi[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Aichi,
                                                  data_series = dat.female[[23]], fh = ik)$score_female
  #### Mie###
  female_score.fvp.Mie[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Mie,
                                                data_series = dat.female[[24]], fh = ik)$score_female
  #### Shiga###
  female_score.fvp.Shiga[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Shiga,
                                                  data_series = dat.female[[25]], fh = ik)$score_female
  #### Kyoto###
  female_score.fvp.Kyoto[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Kyoto,
                                                  data_series = dat.female[[26]], fh = ik)$score_female
  #### Osaka###
  female_score.fvp.Osaka[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Osaka,
                                                  data_series = dat.female[[27]], fh = ik)$score_female
  #### Hyogo###
  female_score.fvp.Hyogo[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Hyogo,
                                                  data_series = dat.female[[28]], fh = ik)$score_female
  #### Nara###
  female_score.fvp.Nara[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Nara,
                                                 data_series = dat.female[[29]], fh = ik)$score_female
  #### Wakayama###
  female_score.fvp.Wakayama[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Wakayama,
                                                     data_series = dat.female[[30]], fh = ik)$score_female
  #### Tottori###
  female_score.fvp.Tottori[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Tottori,
                                                    data_series = dat.female[[31]], fh = ik)$score_female
  #### Shimane###
  female_score.fvp.Shimane[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Shimane,
                                                    data_series = dat.female[[32]], fh = ik)$score_female
  #### Okayama###
  female_score.fvp.Okayama[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Okayama,
                                                    data_series = dat.female[[33]], fh = ik)$score_female
  #### Hiroshima###
  female_score.fvp.Hiroshima[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Hiroshima,
                                                      data_series = dat.female[[34]], fh = ik)$score_female
  #### Yamaguchi###
  female_score.fvp.Yamaguchi[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Yamaguchi,
                                                      data_series = dat.female[[35]], fh = ik)$score_female
  #### Tokushima###
  female_score.fvp.Tokushima[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Tokushima,
                                                      data_series = dat.female[[36]], fh = ik)$score_female
  #### Kagawa###
  female_score.fvp.Kagawa[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Kagawa,
                                                   data_series = dat.female[[37]], fh = ik)$score_female
  #### Ehime###
  female_score.fvp.Ehime[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Ehime,
                                                  data_series = dat.female[[38]], fh = ik)$score_female
  #### Kochi###
  female_score.fvp.Kochi[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Kochi,
                                                  data_series = dat.female[[39]], fh = ik)$score_female
  #### Fukuoka###
  female_score.fvp.Fukuoka[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Fukuoka,
                                                    data_series = dat.female[[40]], fh = ik)$score_female
  #### Saga###
  female_score.fvp.Saga[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Saga,
                                                 data_series = dat.female[[41]], fh = ik)$score_female
  #### Nagasaki###
  female_score.fvp.Nagasaki[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Nagasaki,
                                                     data_series = dat.female[[42]], fh = ik)$score_female
  #### Kumamoto###
  female_score.fvp.Kumamoto[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Kumamoto,
                                                     data_series = dat.female[[43]], fh = ik)$score_female
  #### Oita###
  female_score.fvp.Oita[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Oita,
                                                 data_series = dat.female[[44]], fh = ik)$score_female
  #### Miyazaki###
  female_score.fvp.Miyazaki[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Miyazaki,
                                                     data_series = dat.female[[45]], fh = ik)$score_female
  #### Kagoshima###
  female_score.fvp.Kagoshima[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Kagoshima,
                                                      data_series = dat.female[[46]], fh = ik)$score_female
  #### Okinawa###
  female_score.fvp.Okinawa[ik] = interval_score_fvp_female(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Okinawa,
                                                    data_series = dat.female[[47]], fh = ik)$score_female
}


fvp.jpfemale.score<-cbind(female_score.fvp.Hokkaido , female_score.fvp.Aomori , female_score.fvp.Iwate ,
                          female_score.fvp.Miyagi , female_score.fvp.Akita , female_score.fvp.Yamagata ,
                          female_score.fvp.Fukushima , female_score.fvp.Ibaraki , 
                          female_score.fvp.Tochigi , female_score.fvp.Gunma , female_score.fvp.Saitama,
                          female_score.fvp.Chiba ,  female_score.fvp.Tokyo, female_score.fvp.Kanagawa ,
                          female_score.fvp.Niigata , female_score.fvp.Toyama , female_score.fvp.Ishikawa ,
                          female_score.fvp.Fukui , female_score.fvp.Yamanashi , female_score.fvp.Nagano , 
                          female_score.fvp.Gifu , female_score.fvp.Shizuoka , female_score.fvp.Aichi ,
                          female_score.fvp.Mie ,  female_score.fvp.Shiga , female_score.fvp.Kyoto ,
                          female_score.fvp.Osaka , female_score.fvp.Hyogo ,  female_score.fvp.Nara ,
                          female_score.fvp.Wakayama , female_score.fvp.Tottori , female_score.fvp.Shimane ,
                          female_score.fvp.Okayama , female_score.fvp.Hiroshima , 
                          female_score.fvp.Yamaguchi ,  female_score.fvp.Tokushima , 
                          female_score.fvp.Kagawa , female_score.fvp.Ehime , female_score.fvp.Kochi ,
                          female_score.fvp.Fukuoka , female_score.fvp.Saga , female_score.fvp.Nagasaki ,
                          female_score.fvp.Kumamoto , female_score.fvp.Oita , female_score.fvp.Miyazaki ,
                          female_score.fvp.Kagoshima , female_score.fvp.Okinawa )

colnames(fvp.jpfemale.score)<-c("Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima", "Ibaraki",  "Tochigi", "Gunma", "Saitama",  
                             "Chiba", "Tokyo", "Kanagawa", "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",    
                             "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi",
                             "Tokushima", "Kagawa", "Ehime", "Kochi", "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")
fvp.jpfemale.score = rbind(fvp.jpfemale.score, apply(fvp.jpfemale.score, 2, mean))
write.csv(fvp.jpfemale.score, "fvp.jpfemale.score.csv")





