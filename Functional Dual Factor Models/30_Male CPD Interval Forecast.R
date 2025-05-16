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
  tsr.dat<- rTensor::as.tensor(dat.arr)
  cp.est <- rTensor::cp(tsr.dat, K)
  return( list(mu=mean.list, lambda = cp.est$lambdas, Ulist = cp.est$U, tsr.dat = tsr.dat) )
}

### Forecast Function for HDFTS
cpfts <- function(dat, year_horizon, years, ages, n_pop, K)
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
find_enlarge_val_cp <- function(data_series, fh , nboot = 1000, alpha = 0.8, n_pop)
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
    fore_res = cpfts(dat=dat.rate.male_h_step, year_horizon = fh, n_pop = length(data_series), ages =data_series[[1]]$age,K=1)
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
  
  
  fore_res  = cpfts(dat=dat.series, year_horizon = fh, ages =data_series[[1]]$age, n_pop = length(dat.series), K=1)
  
  fore_cp_male <- list()
  for (il in 1: n_pop){
    fore_cp_male[[il]]<-matrix(rep(log(fore_res[(1:101)+(il-1)*101 ,fh]), nboot), nrow = 101, ncol = nboot)
  }
  
  
  
  boot_PI_male <- list()
  for (il in 1: n_pop){
    boot_PI_male[[il]]<- male.boot_err[[il]] + fore_cp_male[[il]]
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

interval_score_cp <- function(PI_val, data_series, fh, alpha = 0.8)
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
PI_fh_cp<- function(data_series, fh, nboot = 1000)
{
  PI_boot_cp_Hokkaido = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Aomori = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Iwate = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Miyagi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Akita  = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Yamagata = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Fukushima = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Ibaraki = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Tochigi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Gunma = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Saitama = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Chiba = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Tokyo = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Kanagawa = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Niigata = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Toyama = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Ishikawa = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Fukui = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Yamanashi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Nagano = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Gifu = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Shizuoka = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Aichi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Mie = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Shiga = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Kyoto = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Osaka = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Hyogo = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Nara = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Wakayama = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Tottori = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Shimane = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Okayama = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Hiroshima = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Yamaguchi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Tokushima = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Kagawa = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Ehime  = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Kochi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Fukuoka = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Saga = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Nagasaki = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Kumamoto = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Oita = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Miyazaki = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Kagoshima = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_cp_Okinawa = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  
  for(io in 1:(11-fh))
  {
    dat.train<-list()
    for ( i in 1:length(jpn.data)){
      dat.train[[i]]<-extract.years(data_series[[i]], years = 1973:(2012+io))
    }
    dum_pointwise = find_enlarge_val_cp(data_series = dat.train, fh = fh, n_pop = length(data_series))
    PI_boot_cp_Hokkaido[,,io] = dum_pointwise$boot_sample_male[[1]]
    PI_boot_cp_Aomori[,,io] = dum_pointwise$boot_sample_male[[2]]
    PI_boot_cp_Iwate[,,io] = dum_pointwise$boot_sample_male[[3]]
    PI_boot_cp_Miyagi[,,io] = dum_pointwise$boot_sample_male[[4]]
    PI_boot_cp_Akita [,,io] = dum_pointwise$boot_sample_male[[5]]
    PI_boot_cp_Yamagata[,,io] = dum_pointwise$boot_sample_male[[6]]
    PI_boot_cp_Fukushima[,,io] = dum_pointwise$boot_sample_male[[7]]
    PI_boot_cp_Ibaraki[,,io] = dum_pointwise$boot_sample_male[[8]]
    PI_boot_cp_Tochigi[,,io] = dum_pointwise$boot_sample_male[[9]]
    PI_boot_cp_Gunma[,,io] = dum_pointwise$boot_sample_male[[10]]
    PI_boot_cp_Saitama[,,io] = dum_pointwise$boot_sample_male[[11]]
    PI_boot_cp_Chiba[,,io] = dum_pointwise$boot_sample_male[[12]]
    PI_boot_cp_Tokyo[,,io] = dum_pointwise$boot_sample_male[[13]]
    PI_boot_cp_Kanagawa[,,io] = dum_pointwise$boot_sample_male[[14]]
    PI_boot_cp_Niigata[,,io] = dum_pointwise$boot_sample_male[[15]]
    PI_boot_cp_Toyama[,,io] = dum_pointwise$boot_sample_male[[16]]
    PI_boot_cp_Ishikawa[,,io] = dum_pointwise$boot_sample_male[[17]]
    PI_boot_cp_Fukui[,,io] = dum_pointwise$boot_sample_male[[18]]
    PI_boot_cp_Yamanashi[,,io] = dum_pointwise$boot_sample_male[[19]]
    PI_boot_cp_Nagano[,,io] = dum_pointwise$boot_sample_male[[20]]
    PI_boot_cp_Gifu[,,io] = dum_pointwise$boot_sample_male[[21]]
    PI_boot_cp_Shizuoka[,,io] = dum_pointwise$boot_sample_male[[22]]
    PI_boot_cp_Aichi[,,io] = dum_pointwise$boot_sample_male[[23]]
    PI_boot_cp_Mie[,,io] = dum_pointwise$boot_sample_male[[24]]
    PI_boot_cp_Shiga[,,io] = dum_pointwise$boot_sample_male[[25]]
    PI_boot_cp_Kyoto[,,io] = dum_pointwise$boot_sample_male[[26]]
    PI_boot_cp_Osaka[,,io] = dum_pointwise$boot_sample_male[[27]]
    PI_boot_cp_Hyogo[,,io] = dum_pointwise$boot_sample_male[[28]]
    PI_boot_cp_Nara[,,io] = dum_pointwise$boot_sample_male[[29]]
    PI_boot_cp_Wakayama[,,io] = dum_pointwise$boot_sample_male[[30]]
    PI_boot_cp_Tottori[,,io] = dum_pointwise$boot_sample_male[[31]]
    PI_boot_cp_Shimane[,,io] = dum_pointwise$boot_sample_male[[32]]
    PI_boot_cp_Okayama[,,io] = dum_pointwise$boot_sample_male[[33]]
    PI_boot_cp_Hiroshima[,,io] = dum_pointwise$boot_sample_male[[34]]
    PI_boot_cp_Yamaguchi[,,io] = dum_pointwise$boot_sample_male[[35]]
    PI_boot_cp_Tokushima[,,io] = dum_pointwise$boot_sample_male[[36]]
    PI_boot_cp_Kagawa[,,io] = dum_pointwise$boot_sample_male[[37]]
    PI_boot_cp_Ehime [,,io] = dum_pointwise$boot_sample_male[[38]]
    PI_boot_cp_Kochi[,,io] = dum_pointwise$boot_sample_male[[39]]
    PI_boot_cp_Fukuoka[,,io] = dum_pointwise$boot_sample_male[[40]]
    PI_boot_cp_Saga[,,io] = dum_pointwise$boot_sample_male[[41]]
    PI_boot_cp_Nagasaki[,,io] = dum_pointwise$boot_sample_male[[42]]
    PI_boot_cp_Kumamoto[,,io] = dum_pointwise$boot_sample_male[[43]]
    PI_boot_cp_Oita[,,io] = dum_pointwise$boot_sample_male[[44]]
    PI_boot_cp_Miyazaki[,,io] = dum_pointwise$boot_sample_male[[45]]
    PI_boot_cp_Kagoshima[,,io] = dum_pointwise$boot_sample_male[[46]]
    PI_boot_cp_Okinawa[,,io] = dum_pointwise$boot_sample_male[[47]]
  }
  return(list(PI_boot_cp_Hokkaido = PI_boot_cp_Hokkaido,PI_boot_cp_Aomori = PI_boot_cp_Aomori,
              PI_boot_cp_Iwate = PI_boot_cp_Iwate,PI_boot_cp_Miyagi = PI_boot_cp_Miyagi,
              PI_boot_cp_Akita  = PI_boot_cp_Akita ,PI_boot_cp_Yamagata = PI_boot_cp_Yamagata,
              PI_boot_cp_Fukushima = PI_boot_cp_Fukushima, PI_boot_cp_Ibaraki = PI_boot_cp_Ibaraki,
              PI_boot_cp_Tochigi = PI_boot_cp_Tochigi, PI_boot_cp_Gunma = PI_boot_cp_Gunma,
              PI_boot_cp_Saitama = PI_boot_cp_Saitama, PI_boot_cp_Chiba = PI_boot_cp_Chiba, 
              PI_boot_cp_Tokyo = PI_boot_cp_Tokyo, PI_boot_cp_Kanagawa = PI_boot_cp_Kanagawa,
              PI_boot_cp_Niigata = PI_boot_cp_Niigata, PI_boot_cp_Toyama = PI_boot_cp_Toyama,
              PI_boot_cp_Ishikawa = PI_boot_cp_Ishikawa,PI_boot_cp_Fukui = PI_boot_cp_Fukui, 
              PI_boot_cp_Yamanashi = PI_boot_cp_Yamanashi,PI_boot_cp_Nagano = PI_boot_cp_Nagano,
              PI_boot_cp_Gifu = PI_boot_cp_Gifu,PI_boot_cp_Shizuoka = PI_boot_cp_Shizuoka,
              PI_boot_cp_Aichi = PI_boot_cp_Aichi,PI_boot_cp_Mie = PI_boot_cp_Mie,
              PI_boot_cp_Shiga = PI_boot_cp_Shiga,PI_boot_cp_Kyoto = PI_boot_cp_Kyoto,
              PI_boot_cp_Osaka = PI_boot_cp_Osaka, PI_boot_cp_Hyogo = PI_boot_cp_Hyogo,
              PI_boot_cp_Nara = PI_boot_cp_Nara, PI_boot_cp_Wakayama = PI_boot_cp_Wakayama, 
              PI_boot_cp_Tottori = PI_boot_cp_Tottori,PI_boot_cp_Shimane = PI_boot_cp_Shimane,
              PI_boot_cp_Okayama = PI_boot_cp_Okayama,PI_boot_cp_Hiroshima = PI_boot_cp_Hiroshima,
              PI_boot_cp_Yamaguchi = PI_boot_cp_Yamaguchi, PI_boot_cp_Tokushima = PI_boot_cp_Tokushima,
              PI_boot_cp_Kagawa = PI_boot_cp_Kagawa, PI_boot_cp_Ehime  = PI_boot_cp_Ehime, 
              PI_boot_cp_Kochi = PI_boot_cp_Kochi, PI_boot_cp_Fukuoka = PI_boot_cp_Fukuoka,
              PI_boot_cp_Saga = PI_boot_cp_Saga,PI_boot_cp_Nagasaki = PI_boot_cp_Nagasaki,
              PI_boot_cp_Kumamoto = PI_boot_cp_Kumamoto,PI_boot_cp_Oita = PI_boot_cp_Oita,
              PI_boot_cp_Miyazaki = PI_boot_cp_Miyazaki,PI_boot_cp_Kagoshima = PI_boot_cp_Kagoshima,
              PI_boot_cp_Okinawa = PI_boot_cp_Okinawa))
}
library(doParallel)
cl <- makeCluster(4) 
registerDoParallel(cl)

cp_fh_PI = foreach(mm = 1:10, .packages = c("demography","rTensor")) %dopar% PI_fh_cp(data_series = dat.male,  fh = mm)

male_score.cp.Hokkaido = male_score.cp.Aomori = male_score.cp.Iwate = male_score.cp.Miyagi = 
  male_score.cp.Akita = male_score.cp.Yamagata = male_score.cp.Fukushima = male_score.cp.Ibaraki = 
  male_score.cp.Tochigi = male_score.cp.Gunma = male_score.cp.Saitama= male_score.cp.Chiba = 
  male_score.cp.Tokyo= male_score.cp.Kanagawa = male_score.cp.Niigata = male_score.cp.Toyama = 
  male_score.cp.Ishikawa = male_score.cp.Fukui = male_score.cp.Yamanashi = male_score.cp.Nagano = 
  male_score.cp.Gifu = male_score.cp.Shizuoka = male_score.cp.Aichi = male_score.cp.Mie = 
  male_score.cp.Shiga = male_score.cp.Kyoto = male_score.cp.Osaka = male_score.cp.Hyogo = 
  male_score.cp.Nara = male_score.cp.Wakayama = male_score.cp.Tottori = male_score.cp.Shimane = 
  male_score.cp.Okayama = male_score.cp.Hiroshima = male_score.cp.Yamaguchi = male_score.cp.Tokushima = 
  male_score.cp.Kagawa = male_score.cp.Ehime = male_score.cp.Kochi = male_score.cp.Fukuoka = 
  male_score.cp.Saga = male_score.cp.Nagasaki = male_score.cp.Kumamoto = male_score.cp.Oita =   
  male_score.cp.Miyazaki = male_score.cp.Kagoshima = male_score.cp.Okinawa =vector( ,year_horizon)

for(ik in 1:year_horizon)
{
  # me 
  #### Hokkaido###
  male_score.cp.Hokkaido[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Hokkaido,
                                                     data_series = dat.male[[1]], fh = ik)$score_male
  
  
  #### Aomori###
  male_score.cp.Aomori[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Aomori,
                                                   data_series = dat.male[[2]], fh = ik)$score_male
  #### Iwate###
  male_score.cp.Iwate[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Iwate,
                                                  data_series = dat.male[[3]], fh = ik)$score_male
  #### Miyagi###
  male_score.cp.Miyagi[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Miyagi,
                                                   data_series = dat.male[[4]], fh = ik)$score_male
  #### Akita###
  male_score.cp.Akita[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Akita,
                                                  data_series = dat.male[[5]], fh = ik)$score_male
  #### Yamagata###
  male_score.cp.Yamagata[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Yamagata,
                                                     data_series = dat.male[[6]], fh = ik)$score_male
  #### Fukushima###
  male_score.cp.Fukushima[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Fukushima,
                                                      data_series = dat.male[[7]], fh = ik)$score_male
  #### Ibaraki###
  male_score.cp.Ibaraki[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Ibaraki,
                                                    data_series = dat.male[[8]], fh = ik)$score_male
  #### Tochigi###
  male_score.cp.Tochigi[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Tochigi,
                                                    data_series = dat.male[[9]], fh = ik)$score_male
  #### Gunma###
  male_score.cp.Gunma[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Gunma,
                                                  data_series = dat.male[[11]], fh = ik)$score_male
  #### Saitama###
  male_score.cp.Saitama[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Saitama,
                                                    data_series = dat.male[[12]], fh = ik)$score_male
  #### Chiba###
  male_score.cp.Chiba[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Chiba,
                                                  data_series = dat.male[[13]], fh = ik)$score_male
  #### Tokyo###
  male_score.cp.Tokyo[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Tokyo,
                                                  data_series = dat.male[[14]], fh = ik)$score_male
  #### Kanagawa###
  male_score.cp.Kanagawa[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Kanagawa,
                                                     data_series = dat.male[[1]], fh = ik)$score_male
  #### Niigata###
  male_score.cp.Niigata[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Niigata,
                                                    data_series = dat.male[[15]], fh = ik)$score_male
  #### Toyama###
  male_score.cp.Toyama[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Toyama,
                                                   data_series = dat.male[[16]], fh = ik)$score_male
  #### Ishikawa###
  male_score.cp.Ishikawa[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Ishikawa,
                                                     data_series = dat.male[[17]], fh = ik)$score_male
  #### Fukui###
  male_score.cp.Fukui[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Fukui,
                                                  data_series = dat.male[[18]], fh = ik)$score_male
  #### Yamanashi###
  male_score.cp.Yamanashi[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Yamanashi,
                                                      data_series = dat.male[[19]], fh = ik)$score_male
  #### Nagano###
  male_score.cp.Nagano[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Nagano,
                                                   data_series = dat.male[[20]], fh = ik)$score_male
  #### Gifu###
  male_score.cp.Gifu[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Gifu,
                                                 data_series = dat.male[[21]], fh = ik)$score_male
  #### Shizuoka###
  male_score.cp.Shizuoka[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Shizuoka,
                                                     data_series = dat.male[[22]], fh = ik)$score_male
  #### Aichi###
  male_score.cp.Aichi[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Aichi,
                                                  data_series = dat.male[[23]], fh = ik)$score_male
  #### Mie###
  male_score.cp.Mie[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Mie,
                                                data_series = dat.male[[24]], fh = ik)$score_male
  #### Shiga###
  male_score.cp.Shiga[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Shiga,
                                                  data_series = dat.male[[25]], fh = ik)$score_male
  #### Kyoto###
  male_score.cp.Kyoto[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Kyoto,
                                                  data_series = dat.male[[26]], fh = ik)$score_male
  #### Osaka###
  male_score.cp.Osaka[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Osaka,
                                                  data_series = dat.male[[27]], fh = ik)$score_male
  #### Hyogo###
  male_score.cp.Hyogo[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Hyogo,
                                                  data_series = dat.male[[28]], fh = ik)$score_male
  #### Nara###
  male_score.cp.Nara[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Nara,
                                                 data_series = dat.male[[29]], fh = ik)$score_male
  #### Wakayama###
  male_score.cp.Wakayama[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Wakayama,
                                                     data_series = dat.male[[30]], fh = ik)$score_male
  #### Tottori###
  male_score.cp.Tottori[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Tottori,
                                                    data_series = dat.male[[31]], fh = ik)$score_male
  #### Shimane###
  male_score.cp.Shimane[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Shimane,
                                                    data_series = dat.male[[32]], fh = ik)$score_male
  #### Okayama###
  male_score.cp.Okayama[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Okayama,
                                                    data_series = dat.male[[33]], fh = ik)$score_male
  #### Hiroshima###
  male_score.cp.Hiroshima[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Hiroshima,
                                                      data_series = dat.male[[34]], fh = ik)$score_male
  #### Yamaguchi###
  male_score.cp.Yamaguchi[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Yamaguchi,
                                                      data_series = dat.male[[35]], fh = ik)$score_male
  #### Tokushima###
  male_score.cp.Tokushima[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Tokushima,
                                                      data_series = dat.male[[36]], fh = ik)$score_male
  #### Kagawa###
  male_score.cp.Kagawa[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Kagawa,
                                                   data_series = dat.male[[37]], fh = ik)$score_male
  #### Ehime###
  male_score.cp.Ehime[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Ehime,
                                                  data_series = dat.male[[38]], fh = ik)$score_male
  #### Kochi###
  male_score.cp.Kochi[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Kochi,
                                                  data_series = dat.male[[39]], fh = ik)$score_male
  #### Fukuoka###
  male_score.cp.Fukuoka[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Fukuoka,
                                                    data_series = dat.male[[40]], fh = ik)$score_male
  #### Saga###
  male_score.cp.Saga[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Saga,
                                                 data_series = dat.male[[41]], fh = ik)$score_male
  #### Nagasaki###
  male_score.cp.Nagasaki[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Nagasaki,
                                                     data_series = dat.male[[42]], fh = ik)$score_male
  #### Kumamoto###
  male_score.cp.Kumamoto[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Kumamoto,
                                                     data_series = dat.male[[43]], fh = ik)$score_male
  #### Oita###
  male_score.cp.Oita[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Oita,
                                                 data_series = dat.male[[44]], fh = ik)$score_male
  #### Miyazaki###
  male_score.cp.Miyazaki[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Miyazaki,
                                                     data_series = dat.male[[45]], fh = ik)$score_male
  #### Kagoshima###
  male_score.cp.Kagoshima[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Kagoshima,
                                                      data_series = dat.male[[46]], fh = ik)$score_male
  #### Okinawa###
  male_score.cp.Okinawa[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Okinawa,
                                                    data_series = dat.male[[47]], fh = ik)$score_male
}


cp.jpmale.score<-cbind(male_score.cp.Hokkaido , male_score.cp.Aomori , male_score.cp.Iwate ,
                          male_score.cp.Miyagi , male_score.cp.Akita , male_score.cp.Yamagata ,
                          male_score.cp.Fukushima , male_score.cp.Ibaraki , 
                          male_score.cp.Tochigi , male_score.cp.Gunma , male_score.cp.Saitama,
                          male_score.cp.Chiba ,  male_score.cp.Tokyo, male_score.cp.Kanagawa ,
                          male_score.cp.Niigata , male_score.cp.Toyama , male_score.cp.Ishikawa ,
                          male_score.cp.Fukui , male_score.cp.Yamanashi , male_score.cp.Nagano , 
                          male_score.cp.Gifu , male_score.cp.Shizuoka , male_score.cp.Aichi ,
                          male_score.cp.Mie ,  male_score.cp.Shiga , male_score.cp.Kyoto ,
                          male_score.cp.Osaka , male_score.cp.Hyogo ,  male_score.cp.Nara ,
                          male_score.cp.Wakayama , male_score.cp.Tottori , male_score.cp.Shimane ,
                          male_score.cp.Okayama , male_score.cp.Hiroshima , 
                          male_score.cp.Yamaguchi ,  male_score.cp.Tokushima , 
                          male_score.cp.Kagawa , male_score.cp.Ehime , male_score.cp.Kochi ,
                          male_score.cp.Fukuoka , male_score.cp.Saga , male_score.cp.Nagasaki ,
                          male_score.cp.Kumamoto , male_score.cp.Oita , male_score.cp.Miyazaki ,
                          male_score.cp.Kagoshima , male_score.cp.Okinawa )

colnames(cp.jpmale.score)<-c("Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima", "Ibaraki",  "Tochigi", "Gunma", "Saitama",  
                                "Chiba", "Tokyo", "Kanagawa", "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",    
                                "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi",
                                "Tokushima", "Kagawa", "Ehime", "Kochi", "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")


write.csv(cp.jpmale.score, "cp.jpmale.score.csv")





