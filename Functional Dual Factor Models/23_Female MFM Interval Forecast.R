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
mfmfts <- function(dat, year_horizon, years, ages, n_pop)
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
find_enlarge_val_mfm <- function(data_series, fh , nboot = 1000, alpha = 0.8, n_pop)
{
  n_test<-17
  
  # calculate  in-sample forecast curves
  fore_curve_female  = matrix( , length(data_series[[1]]$age)*n_pop, length(data_series[[1]]$year) - n_test - fh )
  for(ij in 1:(length(data_series[[1]]$year) - n_test - fh))
  {
    dat.rate.female_h_step<-list()
    for (ik in 1:n_pop){
      dat.rate.female_h_step[[ik]]<-t(log(extract.years(data_series[[ik]], years = 1973:(1973 + n_test + ij - 1))$rate$female))
    }
    fore_res = mfmfts(dat=dat.rate.female_h_step, year_horizon = fh, n_pop = length(data_series),ages =data_series[[1]]$age)
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
  
  
  fore_res  = mfmfts(dat=dat.series, year_horizon = fh, n_pop = length(dat.series), ages =data_series[[1]]$age)
  
  fore_mfm_female <- list()
  for (il in 1: n_pop){
    fore_mfm_female[[il]]<-matrix(rep(log(fore_res[(1:101)+(il-1)*101 ,fh]), nboot), nrow = 101, ncol = nboot)
  }
  
  
  
  boot_PI_female <- list()
  for (il in 1: n_pop){
    boot_PI_female[[il]]<- female.boot_err[[il]] + fore_mfm_female[[il]]
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

interval_score_mfm <- function(PI_val, data_series, fh, alpha = 0.8)
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
  dat.female[[i]]<-extract.ages(extract.years(jpn.data[[i]], years = 1973:2022), 0:100, combine.upper=TRUE)
}

for (ik in 1:length(dat.female)){
  for(i in 1:length(dat.female[[1]]$age)){
    for (j in 1:length(dat.female[[1]]$year)){
      if (dat.female[[ik]]$rate$female[i,j]==0){
        dat.female[[ik]]$rate$female[i,j]<-0.0000001
      }else{
        dat.female[[ik]]$rate$female[i,j]<-dat.female[[ik]]$rate$female[i,j]
      }
    }
  }
}

####### Prediction Intervals Constructions #######
PI_fh_mfm<- function(data_series, fh, nboot = 1000)
{
  PI_boot_mfm_Hokkaido = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Aomori = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Iwate = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Miyagi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Akita  = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Yamagata = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Fukushima = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Ibaraki = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Tochigi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Gunma = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Saitama = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Chiba = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Tokyo = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Kanagawa = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Niigata = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Toyama = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Ishikawa = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Fukui = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Yamanashi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Nagano = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Gifu = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Shizuoka = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Aichi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Mie = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Shiga = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Kyoto = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Osaka = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Hyogo = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Nara = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Wakayama = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Tottori = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Shimane = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Okayama = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Hiroshima = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Yamaguchi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Tokushima = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Kagawa = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Ehime  = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Kochi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Fukuoka = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Saga = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Nagasaki = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Kumamoto = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Oita = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Miyazaki = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Kagoshima = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_mfm_Okinawa = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  
  for(io in 1:(11-fh))
  {
    dat.train<-list()
    for ( i in 1:length(jpn.data)){
      dat.train[[i]]<-extract.years(data_series[[i]], years = 1973:(2012+io))
    }
    dum_pointwise = find_enlarge_val_mfm(data_series = dat.train, fh = fh, n_pop = length(data_series))
    PI_boot_mfm_Hokkaido[,,io] = dum_pointwise$boot_sample_female[[1]]
    PI_boot_mfm_Aomori[,,io] = dum_pointwise$boot_sample_female[[2]]
    PI_boot_mfm_Iwate[,,io] = dum_pointwise$boot_sample_female[[3]]
    PI_boot_mfm_Miyagi[,,io] = dum_pointwise$boot_sample_female[[4]]
    PI_boot_mfm_Akita [,,io] = dum_pointwise$boot_sample_female[[5]]
    PI_boot_mfm_Yamagata[,,io] = dum_pointwise$boot_sample_female[[6]]
    PI_boot_mfm_Fukushima[,,io] = dum_pointwise$boot_sample_female[[7]]
    PI_boot_mfm_Ibaraki[,,io] = dum_pointwise$boot_sample_female[[8]]
    PI_boot_mfm_Tochigi[,,io] = dum_pointwise$boot_sample_female[[9]]
    PI_boot_mfm_Gunma[,,io] = dum_pointwise$boot_sample_female[[10]]
    PI_boot_mfm_Saitama[,,io] = dum_pointwise$boot_sample_female[[11]]
    PI_boot_mfm_Chiba[,,io] = dum_pointwise$boot_sample_female[[12]]
    PI_boot_mfm_Tokyo[,,io] = dum_pointwise$boot_sample_female[[13]]
    PI_boot_mfm_Kanagawa[,,io] = dum_pointwise$boot_sample_female[[14]]
    PI_boot_mfm_Niigata[,,io] = dum_pointwise$boot_sample_female[[15]]
    PI_boot_mfm_Toyama[,,io] = dum_pointwise$boot_sample_female[[16]]
    PI_boot_mfm_Ishikawa[,,io] = dum_pointwise$boot_sample_female[[17]]
    PI_boot_mfm_Fukui[,,io] = dum_pointwise$boot_sample_female[[18]]
    PI_boot_mfm_Yamanashi[,,io] = dum_pointwise$boot_sample_female[[19]]
    PI_boot_mfm_Nagano[,,io] = dum_pointwise$boot_sample_female[[20]]
    PI_boot_mfm_Gifu[,,io] = dum_pointwise$boot_sample_female[[21]]
    PI_boot_mfm_Shizuoka[,,io] = dum_pointwise$boot_sample_female[[22]]
    PI_boot_mfm_Aichi[,,io] = dum_pointwise$boot_sample_female[[23]]
    PI_boot_mfm_Mie[,,io] = dum_pointwise$boot_sample_female[[24]]
    PI_boot_mfm_Shiga[,,io] = dum_pointwise$boot_sample_female[[25]]
    PI_boot_mfm_Kyoto[,,io] = dum_pointwise$boot_sample_female[[26]]
    PI_boot_mfm_Osaka[,,io] = dum_pointwise$boot_sample_female[[27]]
    PI_boot_mfm_Hyogo[,,io] = dum_pointwise$boot_sample_female[[28]]
    PI_boot_mfm_Nara[,,io] = dum_pointwise$boot_sample_female[[29]]
    PI_boot_mfm_Wakayama[,,io] = dum_pointwise$boot_sample_female[[30]]
    PI_boot_mfm_Tottori[,,io] = dum_pointwise$boot_sample_female[[31]]
    PI_boot_mfm_Shimane[,,io] = dum_pointwise$boot_sample_female[[32]]
    PI_boot_mfm_Okayama[,,io] = dum_pointwise$boot_sample_female[[33]]
    PI_boot_mfm_Hiroshima[,,io] = dum_pointwise$boot_sample_female[[34]]
    PI_boot_mfm_Yamaguchi[,,io] = dum_pointwise$boot_sample_female[[35]]
    PI_boot_mfm_Tokushima[,,io] = dum_pointwise$boot_sample_female[[36]]
    PI_boot_mfm_Kagawa[,,io] = dum_pointwise$boot_sample_female[[37]]
    PI_boot_mfm_Ehime [,,io] = dum_pointwise$boot_sample_female[[38]]
    PI_boot_mfm_Kochi[,,io] = dum_pointwise$boot_sample_female[[39]]
    PI_boot_mfm_Fukuoka[,,io] = dum_pointwise$boot_sample_female[[40]]
    PI_boot_mfm_Saga[,,io] = dum_pointwise$boot_sample_female[[41]]
    PI_boot_mfm_Nagasaki[,,io] = dum_pointwise$boot_sample_female[[42]]
    PI_boot_mfm_Kumamoto[,,io] = dum_pointwise$boot_sample_female[[43]]
    PI_boot_mfm_Oita[,,io] = dum_pointwise$boot_sample_female[[44]]
    PI_boot_mfm_Miyazaki[,,io] = dum_pointwise$boot_sample_female[[45]]
    PI_boot_mfm_Kagoshima[,,io] = dum_pointwise$boot_sample_female[[46]]
    PI_boot_mfm_Okinawa[,,io] = dum_pointwise$boot_sample_female[[47]]
  }
  return(list(PI_boot_mfm_Hokkaido = PI_boot_mfm_Hokkaido,PI_boot_mfm_Aomori = PI_boot_mfm_Aomori,
              PI_boot_mfm_Iwate = PI_boot_mfm_Iwate,PI_boot_mfm_Miyagi = PI_boot_mfm_Miyagi,
              PI_boot_mfm_Akita  = PI_boot_mfm_Akita ,PI_boot_mfm_Yamagata = PI_boot_mfm_Yamagata,
              PI_boot_mfm_Fukushima = PI_boot_mfm_Fukushima, PI_boot_mfm_Ibaraki = PI_boot_mfm_Ibaraki,
              PI_boot_mfm_Tochigi = PI_boot_mfm_Tochigi, PI_boot_mfm_Gunma = PI_boot_mfm_Gunma,
              PI_boot_mfm_Saitama = PI_boot_mfm_Saitama, PI_boot_mfm_Chiba = PI_boot_mfm_Chiba, 
              PI_boot_mfm_Tokyo = PI_boot_mfm_Tokyo, PI_boot_mfm_Kanagawa = PI_boot_mfm_Kanagawa,
              PI_boot_mfm_Niigata = PI_boot_mfm_Niigata, PI_boot_mfm_Toyama = PI_boot_mfm_Toyama,
              PI_boot_mfm_Ishikawa = PI_boot_mfm_Ishikawa,PI_boot_mfm_Fukui = PI_boot_mfm_Fukui, 
              PI_boot_mfm_Yamanashi = PI_boot_mfm_Yamanashi,PI_boot_mfm_Nagano = PI_boot_mfm_Nagano,
              PI_boot_mfm_Gifu = PI_boot_mfm_Gifu,PI_boot_mfm_Shizuoka = PI_boot_mfm_Shizuoka,
              PI_boot_mfm_Aichi = PI_boot_mfm_Aichi,PI_boot_mfm_Mie = PI_boot_mfm_Mie,
              PI_boot_mfm_Shiga = PI_boot_mfm_Shiga,PI_boot_mfm_Kyoto = PI_boot_mfm_Kyoto,
              PI_boot_mfm_Osaka = PI_boot_mfm_Osaka, PI_boot_mfm_Hyogo = PI_boot_mfm_Hyogo,
              PI_boot_mfm_Nara = PI_boot_mfm_Nara, PI_boot_mfm_Wakayama = PI_boot_mfm_Wakayama, 
              PI_boot_mfm_Tottori = PI_boot_mfm_Tottori,PI_boot_mfm_Shimane = PI_boot_mfm_Shimane,
              PI_boot_mfm_Okayama = PI_boot_mfm_Okayama,PI_boot_mfm_Hiroshima = PI_boot_mfm_Hiroshima,
              PI_boot_mfm_Yamaguchi = PI_boot_mfm_Yamaguchi, PI_boot_mfm_Tokushima = PI_boot_mfm_Tokushima,
              PI_boot_mfm_Kagawa = PI_boot_mfm_Kagawa, PI_boot_mfm_Ehime  = PI_boot_mfm_Ehime, 
              PI_boot_mfm_Kochi = PI_boot_mfm_Kochi, PI_boot_mfm_Fukuoka = PI_boot_mfm_Fukuoka,
              PI_boot_mfm_Saga = PI_boot_mfm_Saga,PI_boot_mfm_Nagasaki = PI_boot_mfm_Nagasaki,
              PI_boot_mfm_Kumamoto = PI_boot_mfm_Kumamoto,PI_boot_mfm_Oita = PI_boot_mfm_Oita,
              PI_boot_mfm_Miyazaki = PI_boot_mfm_Miyazaki,PI_boot_mfm_Kagoshima = PI_boot_mfm_Kagoshima,
              PI_boot_mfm_Okinawa = PI_boot_mfm_Okinawa))
}
library(doParallel)
cl <- makeCluster(4) 
registerDoParallel(cl)

mfm_fh_PI = foreach(mm = 1:10, .packages = "demography") %dopar% PI_fh_mfm(data_series = dat.female,  fh = mm)

female_score.mfm.Hokkaido = female_score.mfm.Aomori = female_score.mfm.Iwate = female_score.mfm.Miyagi = 
  female_score.mfm.Akita = female_score.mfm.Yamagata = female_score.mfm.Fukushima = female_score.mfm.Ibaraki = 
  female_score.mfm.Tochigi = female_score.mfm.Gunma = female_score.mfm.Saitama= female_score.mfm.Chiba = 
  female_score.mfm.Tokyo= female_score.mfm.Kanagawa = female_score.mfm.Niigata = female_score.mfm.Toyama = 
  female_score.mfm.Ishikawa = female_score.mfm.Fukui = female_score.mfm.Yamanashi = female_score.mfm.Nagano = 
  female_score.mfm.Gifu = female_score.mfm.Shizuoka = female_score.mfm.Aichi = female_score.mfm.Mie = 
  female_score.mfm.Shiga = female_score.mfm.Kyoto = female_score.mfm.Osaka = female_score.mfm.Hyogo = 
  female_score.mfm.Nara = female_score.mfm.Wakayama = female_score.mfm.Tottori = female_score.mfm.Shimane = 
  female_score.mfm.Okayama = female_score.mfm.Hiroshima = female_score.mfm.Yamaguchi = female_score.mfm.Tokushima = 
  female_score.mfm.Kagawa = female_score.mfm.Ehime = female_score.mfm.Kochi = female_score.mfm.Fukuoka = 
  female_score.mfm.Saga = female_score.mfm.Nagasaki = female_score.mfm.Kumamoto = female_score.mfm.Oita =   
  female_score.mfm.Miyazaki = female_score.mfm.Kagoshima = female_score.mfm.Okinawa =vector( ,year_horizon)

for(ik in 1:year_horizon)
{
  # me 
  #### Hokkaido###
  female_score.mfm.Hokkaido[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Hokkaido,
                                                   data_series = dat.female[[1]], fh = ik)$score_female
  
  
  #### Aomori###
  female_score.mfm.Aomori[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Aomori,
                                                 data_series = dat.female[[2]], fh = ik)$score_female
  #### Iwate###
  female_score.mfm.Iwate[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Iwate,
                                                data_series = dat.female[[3]], fh = ik)$score_female
  #### Miyagi###
  female_score.mfm.Miyagi[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Miyagi,
                                                 data_series = dat.female[[4]], fh = ik)$score_female
  #### Akita###
  female_score.mfm.Akita[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Akita,
                                                data_series = dat.female[[5]], fh = ik)$score_female
  #### Yamagata###
  female_score.mfm.Yamagata[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Yamagata,
                                                   data_series = dat.female[[6]], fh = ik)$score_female
  #### Fukushima###
  female_score.mfm.Fukushima[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Fukushima,
                                                    data_series = dat.female[[7]], fh = ik)$score_female
  #### Ibaraki###
  female_score.mfm.Ibaraki[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Ibaraki,
                                                  data_series = dat.female[[8]], fh = ik)$score_female
  #### Tochigi###
  female_score.mfm.Tochigi[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Tochigi,
                                                  data_series = dat.female[[9]], fh = ik)$score_female
  #### Gunma###
  female_score.mfm.Gunma[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Gunma,
                                                data_series = dat.female[[11]], fh = ik)$score_female
  #### Saitama###
  female_score.mfm.Saitama[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Saitama,
                                                  data_series = dat.female[[12]], fh = ik)$score_female
  #### Chiba###
  female_score.mfm.Chiba[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Chiba,
                                                data_series = dat.female[[13]], fh = ik)$score_female
  #### Tokyo###
  female_score.mfm.Tokyo[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Tokyo,
                                                data_series = dat.female[[14]], fh = ik)$score_female
  #### Kanagawa###
  female_score.mfm.Kanagawa[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Kanagawa,
                                                   data_series = dat.female[[1]], fh = ik)$score_female
  #### Niigata###
  female_score.mfm.Niigata[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Niigata,
                                                  data_series = dat.female[[15]], fh = ik)$score_female
  #### Toyama###
  female_score.mfm.Toyama[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Toyama,
                                                 data_series = dat.female[[16]], fh = ik)$score_female
  #### Ishikawa###
  female_score.mfm.Ishikawa[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Ishikawa,
                                                   data_series = dat.female[[17]], fh = ik)$score_female
  #### Fukui###
  female_score.mfm.Fukui[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Fukui,
                                                data_series = dat.female[[18]], fh = ik)$score_female
  #### Yamanashi###
  female_score.mfm.Yamanashi[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Yamanashi,
                                                    data_series = dat.female[[19]], fh = ik)$score_female
  #### Nagano###
  female_score.mfm.Nagano[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Nagano,
                                                 data_series = dat.female[[20]], fh = ik)$score_female
  #### Gifu###
  female_score.mfm.Gifu[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Gifu,
                                               data_series = dat.female[[21]], fh = ik)$score_female
  #### Shizuoka###
  female_score.mfm.Shizuoka[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Shizuoka,
                                                   data_series = dat.female[[22]], fh = ik)$score_female
  #### Aichi###
  female_score.mfm.Aichi[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Aichi,
                                                data_series = dat.female[[23]], fh = ik)$score_female
  #### Mie###
  female_score.mfm.Mie[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Mie,
                                              data_series = dat.female[[24]], fh = ik)$score_female
  #### Shiga###
  female_score.mfm.Shiga[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Shiga,
                                                data_series = dat.female[[25]], fh = ik)$score_female
  #### Kyoto###
  female_score.mfm.Kyoto[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Kyoto,
                                                data_series = dat.female[[26]], fh = ik)$score_female
  #### Osaka###
  female_score.mfm.Osaka[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Osaka,
                                                data_series = dat.female[[27]], fh = ik)$score_female
  #### Hyogo###
  female_score.mfm.Hyogo[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Hyogo,
                                                data_series = dat.female[[28]], fh = ik)$score_female
  #### Nara###
  female_score.mfm.Nara[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Nara,
                                               data_series = dat.female[[29]], fh = ik)$score_female
  #### Wakayama###
  female_score.mfm.Wakayama[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Wakayama,
                                                   data_series = dat.female[[30]], fh = ik)$score_female
  #### Tottori###
  female_score.mfm.Tottori[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Tottori,
                                                  data_series = dat.female[[31]], fh = ik)$score_female
  #### Shimane###
  female_score.mfm.Shimane[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Shimane,
                                                  data_series = dat.female[[32]], fh = ik)$score_female
  #### Okayama###
  female_score.mfm.Okayama[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Okayama,
                                                  data_series = dat.female[[33]], fh = ik)$score_female
  #### Hiroshima###
  female_score.mfm.Hiroshima[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Hiroshima,
                                                    data_series = dat.female[[34]], fh = ik)$score_female
  #### Yamaguchi###
  female_score.mfm.Yamaguchi[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Yamaguchi,
                                                    data_series = dat.female[[35]], fh = ik)$score_female
  #### Tokushima###
  female_score.mfm.Tokushima[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Tokushima,
                                                    data_series = dat.female[[36]], fh = ik)$score_female
  #### Kagawa###
  female_score.mfm.Kagawa[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Kagawa,
                                                 data_series = dat.female[[37]], fh = ik)$score_female
  #### Ehime###
  female_score.mfm.Ehime[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Ehime,
                                                data_series = dat.female[[38]], fh = ik)$score_female
  #### Kochi###
  female_score.mfm.Kochi[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Kochi,
                                                data_series = dat.female[[39]], fh = ik)$score_female
  #### Fukuoka###
  female_score.mfm.Fukuoka[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Fukuoka,
                                                  data_series = dat.female[[40]], fh = ik)$score_female
  #### Saga###
  female_score.mfm.Saga[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Saga,
                                               data_series = dat.female[[41]], fh = ik)$score_female
  #### Nagasaki###
  female_score.mfm.Nagasaki[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Nagasaki,
                                                   data_series = dat.female[[42]], fh = ik)$score_female
  #### Kumamoto###
  female_score.mfm.Kumamoto[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Kumamoto,
                                                   data_series = dat.female[[43]], fh = ik)$score_female
  #### Oita###
  female_score.mfm.Oita[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Oita,
                                               data_series = dat.female[[44]], fh = ik)$score_female
  #### Miyazaki###
  female_score.mfm.Miyazaki[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Miyazaki,
                                                   data_series = dat.female[[45]], fh = ik)$score_female
  #### Kagoshima###
  female_score.mfm.Kagoshima[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Kagoshima,
                                                    data_series = dat.female[[46]], fh = ik)$score_female
  #### Okinawa###
  female_score.mfm.Okinawa[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Okinawa,
                                                  data_series = dat.female[[47]], fh = ik)$score_female
}


mfm.jpfemale.score<-cbind(female_score.mfm.Hokkaido , female_score.mfm.Aomori , female_score.mfm.Iwate ,
         female_score.mfm.Miyagi , female_score.mfm.Akita , female_score.mfm.Yamagata ,
         female_score.mfm.Fukushima , female_score.mfm.Ibaraki , 
         female_score.mfm.Tochigi , female_score.mfm.Gunma , female_score.mfm.Saitama,
         female_score.mfm.Chiba ,  female_score.mfm.Tokyo, female_score.mfm.Kanagawa ,
         female_score.mfm.Niigata , female_score.mfm.Toyama , female_score.mfm.Ishikawa ,
         female_score.mfm.Fukui , female_score.mfm.Yamanashi , female_score.mfm.Nagano , 
         female_score.mfm.Gifu , female_score.mfm.Shizuoka , female_score.mfm.Aichi ,
         female_score.mfm.Mie ,  female_score.mfm.Shiga , female_score.mfm.Kyoto ,
         female_score.mfm.Osaka , female_score.mfm.Hyogo ,  female_score.mfm.Nara ,
         female_score.mfm.Wakayama , female_score.mfm.Tottori , female_score.mfm.Shimane ,
         female_score.mfm.Okayama , female_score.mfm.Hiroshima , 
         female_score.mfm.Yamaguchi ,  female_score.mfm.Tokushima , 
         female_score.mfm.Kagawa , female_score.mfm.Ehime , female_score.mfm.Kochi ,
         female_score.mfm.Fukuoka , female_score.mfm.Saga , female_score.mfm.Nagasaki ,
         female_score.mfm.Kumamoto , female_score.mfm.Oita , female_score.mfm.Miyazaki ,
         female_score.mfm.Kagoshima , female_score.mfm.Okinawa )

colnames(mfm.jpfemale.score)<-c("Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima", "Ibaraki",  "Tochigi", "Gunma", "Saitama",  
                               "Chiba", "Tokyo", "Kanagawa", "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",    
                               "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi",
                               "Tokushima", "Kagawa", "Ehime", "Kochi", "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")



write.csv(mfm.jpfemale.score, "mfm.jpfemale.score.csv")
