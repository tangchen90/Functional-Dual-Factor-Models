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
boot_male <-function(err, nboot)
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
  fore_curve_male  = matrix( , length(data_series[[1]]$age)*n_pop, length(data_series[[1]]$year) - n_test - fh )
  for(ij in 1:(length(data_series[[1]]$year) - n_test - fh))
  {
    dat.rate.male_h_step<-list()
    for (ik in 1:n_pop){
      dat.rate.male_h_step[[ik]]<-t(log(extract.years(data_series[[ik]], years = 1973:(1973 + n_test + ij - 1))$rate$male))
    }
    fore_res = mfmfts(dat=dat.rate.male_h_step, year_horizon = fh, n_pop = length(data_series),ages =data_series[[1]]$age)
    # fill in gender specific in-sample forecast curves
    fore_curve_male[,ij] = log(fore_res[,fh])
  }
  
  # holdout data samples
  true_dat = list()
  for (im in 1: n_pop){
    true_dat[[im]] = as.data.frame(log(data_series[[im]]$rate$male[, (n_test + fh + 1):length(data_series[[im]]$year)]))
  }
  
  
  
  # male sample
  holdout_val_male = do.call(rbind, true_dat)
  
  
  err_male = holdout_val_male - fore_curve_male
  
  male.err<-list()
  for (il in 1: n_pop){
    male.err[[il]]<-err_male[(1:101)+(il-1)*101 ,]
  }
  male.boot_err<-list()
  for (il in 1: n_pop){
    male.boot_err[[il]]<-boot_male(male.err[[il]], 1000)
  }
  
  # constructing PI
  dat.series<-list()
  for (ik in 1:n_pop){
    dat.series[[ik]]<-t(log(data_series[[ik]]$rate$male))
  }
  
  
  fore_res  = mfmfts(dat=dat.series, year_horizon = fh, n_pop = length(dat.series), ages =data_series[[1]]$age)
  
  fore_mfm_male <- list()
  for (il in 1: n_pop){
    fore_mfm_male[[il]]<-matrix(rep(log(fore_res[(1:101)+(il-1)*101 ,fh]), nboot), nrow = 101, ncol = nboot)
  }
  
  
  
  boot_PI_male <- list()
  for (il in 1: n_pop){
    boot_PI_male[[il]]<- male.boot_err[[il]] + fore_mfm_male[[il]]
  }
  
  boot_PI_lb_ub_male <- list()
  for (il in 1: n_pop){
    boot_PI_lb_ub_male[[il]]<- apply(boot_PI_male[[il]], 1, quantile, c((1 - alpha)/2, (1 + alpha)/2))
  }
  
  
  return(list(boot_PI_lb_ub_male = boot_PI_lb_ub_male, boot_sample_male = boot_PI_male))
}    



###############################################################################
# compute interval scores  for a particular forecast horizon 
###############################################################################

# PI_val: one- to 10-step-ahead prediction intervals
# alpha: nominal coverage probability

interval_score_mfm <- function(PI_val, data_series, fh, alpha = 0.8)
{
  
  test_val_male<- extract.years(data_series, (2012+fh):2022)$rate$male
  
  # transform back to the original scale
  
  boot_sample_male = exp(PI_val)
  boot_index_male = which(boot_sample_male > 1)
  boot_index_below_male = which(boot_sample_male < 0)
  if(length(boot_index_male) > 0)
  {
    boot_sample_v1_male = replace(boot_sample_male, boot_index_male, 1)
  } else 
  {
    boot_sample_v1_male = boot_sample_male
  }
  if(length(boot_index_below_male) > 0)
  {
    boot_sample_v2_male = replace(boot_sample_v1_male, boot_index_below_male, 0)
  } else
  {
    boot_sample_v2_male = boot_sample_v1_male
  }
  
  # lower and upper bounds  
  # male series
  
  dummy_male = apply(boot_sample_male, c(1,3), quantile, c((1 - alpha)/2, (1 + alpha)/2), na.rm = TRUE)
  PI_lb_val_male = dummy_male[1,,]
  PI_ub_val_male  = dummy_male[2,,]
  
  lb_ind_male  = ifelse(test_val_male < PI_lb_val_male, 1, 0)
  ub_ind_male  = ifelse(test_val_male > PI_ub_val_male, 1, 0)
  score_male  = mean((PI_ub_val_male - PI_lb_val_male) + 2/(1 - alpha) * (PI_lb_val_male - test_val_male) * lb_ind_male + 
                         2/(1 - alpha) * (test_val_male - PI_ub_val_male) * ub_ind_male)
  
  return(list(score_male = score_male))
}



dat.male<-list()
for ( i in 1:length(jpn.data)){
  dat.male[[i]]<-extract.ages(extract.years(jpn.data[[i]], years = 1973:2022), 0:100, combine.upper=TRUE)
}

for (ik in 1:length(dat.male)){
  for(i in 1:length(dat.male[[1]]$age)){
    for (j in 1:length(dat.male[[1]]$year)){
      if (dat.male[[ik]]$rate$male[i,j]==0){
        dat.male[[ik]]$rate$male[i,j]<-0.0000001
      }else{
        dat.male[[ik]]$rate$male[i,j]<-dat.male[[ik]]$rate$male[i,j]
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
    PI_boot_mfm_Hokkaido[,,io] = dum_pointwise$boot_sample_male[[1]]
    PI_boot_mfm_Aomori[,,io] = dum_pointwise$boot_sample_male[[2]]
    PI_boot_mfm_Iwate[,,io] = dum_pointwise$boot_sample_male[[3]]
    PI_boot_mfm_Miyagi[,,io] = dum_pointwise$boot_sample_male[[4]]
    PI_boot_mfm_Akita [,,io] = dum_pointwise$boot_sample_male[[5]]
    PI_boot_mfm_Yamagata[,,io] = dum_pointwise$boot_sample_male[[6]]
    PI_boot_mfm_Fukushima[,,io] = dum_pointwise$boot_sample_male[[7]]
    PI_boot_mfm_Ibaraki[,,io] = dum_pointwise$boot_sample_male[[8]]
    PI_boot_mfm_Tochigi[,,io] = dum_pointwise$boot_sample_male[[9]]
    PI_boot_mfm_Gunma[,,io] = dum_pointwise$boot_sample_male[[10]]
    PI_boot_mfm_Saitama[,,io] = dum_pointwise$boot_sample_male[[11]]
    PI_boot_mfm_Chiba[,,io] = dum_pointwise$boot_sample_male[[12]]
    PI_boot_mfm_Tokyo[,,io] = dum_pointwise$boot_sample_male[[13]]
    PI_boot_mfm_Kanagawa[,,io] = dum_pointwise$boot_sample_male[[14]]
    PI_boot_mfm_Niigata[,,io] = dum_pointwise$boot_sample_male[[15]]
    PI_boot_mfm_Toyama[,,io] = dum_pointwise$boot_sample_male[[16]]
    PI_boot_mfm_Ishikawa[,,io] = dum_pointwise$boot_sample_male[[17]]
    PI_boot_mfm_Fukui[,,io] = dum_pointwise$boot_sample_male[[18]]
    PI_boot_mfm_Yamanashi[,,io] = dum_pointwise$boot_sample_male[[19]]
    PI_boot_mfm_Nagano[,,io] = dum_pointwise$boot_sample_male[[20]]
    PI_boot_mfm_Gifu[,,io] = dum_pointwise$boot_sample_male[[21]]
    PI_boot_mfm_Shizuoka[,,io] = dum_pointwise$boot_sample_male[[22]]
    PI_boot_mfm_Aichi[,,io] = dum_pointwise$boot_sample_male[[23]]
    PI_boot_mfm_Mie[,,io] = dum_pointwise$boot_sample_male[[24]]
    PI_boot_mfm_Shiga[,,io] = dum_pointwise$boot_sample_male[[25]]
    PI_boot_mfm_Kyoto[,,io] = dum_pointwise$boot_sample_male[[26]]
    PI_boot_mfm_Osaka[,,io] = dum_pointwise$boot_sample_male[[27]]
    PI_boot_mfm_Hyogo[,,io] = dum_pointwise$boot_sample_male[[28]]
    PI_boot_mfm_Nara[,,io] = dum_pointwise$boot_sample_male[[29]]
    PI_boot_mfm_Wakayama[,,io] = dum_pointwise$boot_sample_male[[30]]
    PI_boot_mfm_Tottori[,,io] = dum_pointwise$boot_sample_male[[31]]
    PI_boot_mfm_Shimane[,,io] = dum_pointwise$boot_sample_male[[32]]
    PI_boot_mfm_Okayama[,,io] = dum_pointwise$boot_sample_male[[33]]
    PI_boot_mfm_Hiroshima[,,io] = dum_pointwise$boot_sample_male[[34]]
    PI_boot_mfm_Yamaguchi[,,io] = dum_pointwise$boot_sample_male[[35]]
    PI_boot_mfm_Tokushima[,,io] = dum_pointwise$boot_sample_male[[36]]
    PI_boot_mfm_Kagawa[,,io] = dum_pointwise$boot_sample_male[[37]]
    PI_boot_mfm_Ehime [,,io] = dum_pointwise$boot_sample_male[[38]]
    PI_boot_mfm_Kochi[,,io] = dum_pointwise$boot_sample_male[[39]]
    PI_boot_mfm_Fukuoka[,,io] = dum_pointwise$boot_sample_male[[40]]
    PI_boot_mfm_Saga[,,io] = dum_pointwise$boot_sample_male[[41]]
    PI_boot_mfm_Nagasaki[,,io] = dum_pointwise$boot_sample_male[[42]]
    PI_boot_mfm_Kumamoto[,,io] = dum_pointwise$boot_sample_male[[43]]
    PI_boot_mfm_Oita[,,io] = dum_pointwise$boot_sample_male[[44]]
    PI_boot_mfm_Miyazaki[,,io] = dum_pointwise$boot_sample_male[[45]]
    PI_boot_mfm_Kagoshima[,,io] = dum_pointwise$boot_sample_male[[46]]
    PI_boot_mfm_Okinawa[,,io] = dum_pointwise$boot_sample_male[[47]]
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

mfm_fh_PI = foreach(mm = 1:10, .packages = "demography") %dopar% PI_fh_mfm(data_series = dat.male,  fh = mm)

male_score.mfm.Hokkaido = male_score.mfm.Aomori = male_score.mfm.Iwate = male_score.mfm.Miyagi = 
  male_score.mfm.Akita = male_score.mfm.Yamagata = male_score.mfm.Fukushima = male_score.mfm.Ibaraki = 
  male_score.mfm.Tochigi = male_score.mfm.Gunma = male_score.mfm.Saitama= male_score.mfm.Chiba = 
  male_score.mfm.Tokyo= male_score.mfm.Kanagawa = male_score.mfm.Niigata = male_score.mfm.Toyama = 
  male_score.mfm.Ishikawa = male_score.mfm.Fukui = male_score.mfm.Yamanashi = male_score.mfm.Nagano = 
  male_score.mfm.Gifu = male_score.mfm.Shizuoka = male_score.mfm.Aichi = male_score.mfm.Mie = 
  male_score.mfm.Shiga = male_score.mfm.Kyoto = male_score.mfm.Osaka = male_score.mfm.Hyogo = 
  male_score.mfm.Nara = male_score.mfm.Wakayama = male_score.mfm.Tottori = male_score.mfm.Shimane = 
  male_score.mfm.Okayama = male_score.mfm.Hiroshima = male_score.mfm.Yamaguchi = male_score.mfm.Tokushima = 
  male_score.mfm.Kagawa = male_score.mfm.Ehime = male_score.mfm.Kochi = male_score.mfm.Fukuoka = 
  male_score.mfm.Saga = male_score.mfm.Nagasaki = male_score.mfm.Kumamoto = male_score.mfm.Oita =   
  male_score.mfm.Miyazaki = male_score.mfm.Kagoshima = male_score.mfm.Okinawa =vector( ,year_horizon)

for(ik in 1:year_horizon)
{
  # me 
  #### Hokkaido###
  male_score.mfm.Hokkaido[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Hokkaido,
                                                   data_series = dat.male[[1]], fh = ik)$score_male
  
  
  #### Aomori###
  male_score.mfm.Aomori[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Aomori,
                                                 data_series = dat.male[[2]], fh = ik)$score_male
  #### Iwate###
  male_score.mfm.Iwate[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Iwate,
                                                data_series = dat.male[[3]], fh = ik)$score_male
  #### Miyagi###
  male_score.mfm.Miyagi[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Miyagi,
                                                 data_series = dat.male[[4]], fh = ik)$score_male
  #### Akita###
  male_score.mfm.Akita[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Akita,
                                                data_series = dat.male[[5]], fh = ik)$score_male
  #### Yamagata###
  male_score.mfm.Yamagata[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Yamagata,
                                                   data_series = dat.male[[6]], fh = ik)$score_male
  #### Fukushima###
  male_score.mfm.Fukushima[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Fukushima,
                                                    data_series = dat.male[[7]], fh = ik)$score_male
  #### Ibaraki###
  male_score.mfm.Ibaraki[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Ibaraki,
                                                  data_series = dat.male[[8]], fh = ik)$score_male
  #### Tochigi###
  male_score.mfm.Tochigi[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Tochigi,
                                                  data_series = dat.male[[9]], fh = ik)$score_male
  #### Gunma###
  male_score.mfm.Gunma[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Gunma,
                                                data_series = dat.male[[11]], fh = ik)$score_male
  #### Saitama###
  male_score.mfm.Saitama[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Saitama,
                                                  data_series = dat.male[[12]], fh = ik)$score_male
  #### Chiba###
  male_score.mfm.Chiba[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Chiba,
                                                data_series = dat.male[[13]], fh = ik)$score_male
  #### Tokyo###
  male_score.mfm.Tokyo[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Tokyo,
                                                data_series = dat.male[[14]], fh = ik)$score_male
  #### Kanagawa###
  male_score.mfm.Kanagawa[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Kanagawa,
                                                   data_series = dat.male[[1]], fh = ik)$score_male
  #### Niigata###
  male_score.mfm.Niigata[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Niigata,
                                                  data_series = dat.male[[15]], fh = ik)$score_male
  #### Toyama###
  male_score.mfm.Toyama[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Toyama,
                                                 data_series = dat.male[[16]], fh = ik)$score_male
  #### Ishikawa###
  male_score.mfm.Ishikawa[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Ishikawa,
                                                   data_series = dat.male[[17]], fh = ik)$score_male
  #### Fukui###
  male_score.mfm.Fukui[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Fukui,
                                                data_series = dat.male[[18]], fh = ik)$score_male
  #### Yamanashi###
  male_score.mfm.Yamanashi[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Yamanashi,
                                                    data_series = dat.male[[19]], fh = ik)$score_male
  #### Nagano###
  male_score.mfm.Nagano[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Nagano,
                                                 data_series = dat.male[[20]], fh = ik)$score_male
  #### Gifu###
  male_score.mfm.Gifu[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Gifu,
                                               data_series = dat.male[[21]], fh = ik)$score_male
  #### Shizuoka###
  male_score.mfm.Shizuoka[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Shizuoka,
                                                   data_series = dat.male[[22]], fh = ik)$score_male
  #### Aichi###
  male_score.mfm.Aichi[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Aichi,
                                                data_series = dat.male[[23]], fh = ik)$score_male
  #### Mie###
  male_score.mfm.Mie[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Mie,
                                              data_series = dat.male[[24]], fh = ik)$score_male
  #### Shiga###
  male_score.mfm.Shiga[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Shiga,
                                                data_series = dat.male[[25]], fh = ik)$score_male
  #### Kyoto###
  male_score.mfm.Kyoto[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Kyoto,
                                                data_series = dat.male[[26]], fh = ik)$score_male
  #### Osaka###
  male_score.mfm.Osaka[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Osaka,
                                                data_series = dat.male[[27]], fh = ik)$score_male
  #### Hyogo###
  male_score.mfm.Hyogo[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Hyogo,
                                                data_series = dat.male[[28]], fh = ik)$score_male
  #### Nara###
  male_score.mfm.Nara[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Nara,
                                               data_series = dat.male[[29]], fh = ik)$score_male
  #### Wakayama###
  male_score.mfm.Wakayama[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Wakayama,
                                                   data_series = dat.male[[30]], fh = ik)$score_male
  #### Tottori###
  male_score.mfm.Tottori[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Tottori,
                                                  data_series = dat.male[[31]], fh = ik)$score_male
  #### Shimane###
  male_score.mfm.Shimane[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Shimane,
                                                  data_series = dat.male[[32]], fh = ik)$score_male
  #### Okayama###
  male_score.mfm.Okayama[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Okayama,
                                                  data_series = dat.male[[33]], fh = ik)$score_male
  #### Hiroshima###
  male_score.mfm.Hiroshima[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Hiroshima,
                                                    data_series = dat.male[[34]], fh = ik)$score_male
  #### Yamaguchi###
  male_score.mfm.Yamaguchi[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Yamaguchi,
                                                    data_series = dat.male[[35]], fh = ik)$score_male
  #### Tokushima###
  male_score.mfm.Tokushima[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Tokushima,
                                                    data_series = dat.male[[36]], fh = ik)$score_male
  #### Kagawa###
  male_score.mfm.Kagawa[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Kagawa,
                                                 data_series = dat.male[[37]], fh = ik)$score_male
  #### Ehime###
  male_score.mfm.Ehime[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Ehime,
                                                data_series = dat.male[[38]], fh = ik)$score_male
  #### Kochi###
  male_score.mfm.Kochi[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Kochi,
                                                data_series = dat.male[[39]], fh = ik)$score_male
  #### Fukuoka###
  male_score.mfm.Fukuoka[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Fukuoka,
                                                  data_series = dat.male[[40]], fh = ik)$score_male
  #### Saga###
  male_score.mfm.Saga[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Saga,
                                               data_series = dat.male[[41]], fh = ik)$score_male
  #### Nagasaki###
  male_score.mfm.Nagasaki[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Nagasaki,
                                                   data_series = dat.male[[42]], fh = ik)$score_male
  #### Kumamoto###
  male_score.mfm.Kumamoto[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Kumamoto,
                                                   data_series = dat.male[[43]], fh = ik)$score_male
  #### Oita###
  male_score.mfm.Oita[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Oita,
                                               data_series = dat.male[[44]], fh = ik)$score_male
  #### Miyazaki###
  male_score.mfm.Miyazaki[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Miyazaki,
                                                   data_series = dat.male[[45]], fh = ik)$score_male
  #### Kagoshima###
  male_score.mfm.Kagoshima[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Kagoshima,
                                                    data_series = dat.male[[46]], fh = ik)$score_male
  #### Okinawa###
  male_score.mfm.Okinawa[ik] = interval_score_mfm(PI_val = mfm_fh_PI[[ik]]$PI_boot_mfm_Okinawa,
                                                  data_series = dat.male[[47]], fh = ik)$score_male
}


mfm.jpmale.score<-cbind(male_score.mfm.Hokkaido , male_score.mfm.Aomori , male_score.mfm.Iwate ,
         male_score.mfm.Miyagi , male_score.mfm.Akita , male_score.mfm.Yamagata ,
         male_score.mfm.Fukushima , male_score.mfm.Ibaraki , 
         male_score.mfm.Tochigi , male_score.mfm.Gunma , male_score.mfm.Saitama,
         male_score.mfm.Chiba ,  male_score.mfm.Tokyo, male_score.mfm.Kanagawa ,
         male_score.mfm.Niigata , male_score.mfm.Toyama , male_score.mfm.Ishikawa ,
         male_score.mfm.Fukui , male_score.mfm.Yamanashi , male_score.mfm.Nagano , 
         male_score.mfm.Gifu , male_score.mfm.Shizuoka , male_score.mfm.Aichi ,
         male_score.mfm.Mie ,  male_score.mfm.Shiga , male_score.mfm.Kyoto ,
         male_score.mfm.Osaka , male_score.mfm.Hyogo ,  male_score.mfm.Nara ,
         male_score.mfm.Wakayama , male_score.mfm.Tottori , male_score.mfm.Shimane ,
         male_score.mfm.Okayama , male_score.mfm.Hiroshima , 
         male_score.mfm.Yamaguchi ,  male_score.mfm.Tokushima , 
         male_score.mfm.Kagawa , male_score.mfm.Ehime , male_score.mfm.Kochi ,
         male_score.mfm.Fukuoka , male_score.mfm.Saga , male_score.mfm.Nagasaki ,
         male_score.mfm.Kumamoto , male_score.mfm.Oita , male_score.mfm.Miyazaki ,
         male_score.mfm.Kagoshima , male_score.mfm.Okinawa )

colnames(mfm.jpmale.score)<-c("Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima", "Ibaraki",  "Tochigi", "Gunma", "Saitama",  
                               "Chiba", "Tokyo", "Kanagawa", "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",    
                               "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi",
                               "Tokushima", "Kagawa", "Ehime", "Kochi", "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")



write.csv(mfm.jpmale.score, "mfm.jpmale.score.csv")
