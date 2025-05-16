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
find_enlarge_val_cp <- function(data_series, fh , nboot = 1000, alpha = 0.8, n_pop)
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
    fore_res = cpfts(dat=dat.rate.female_h_step, year_horizon = fh, n_pop = length(data_series), ages =data_series[[1]]$age,K=1)
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
  
  
  fore_res  = cpfts(dat=dat.series, year_horizon = fh, ages =data_series[[1]]$age, n_pop = length(dat.series), K=1)
  
  fore_cp_female <- list()
  for (il in 1: n_pop){
    fore_cp_female[[il]]<-matrix(rep(log(fore_res[(1:101)+(il-1)*101 ,fh]), nboot), nrow = 101, ncol = nboot)
  }
  
  
  
  boot_PI_female <- list()
  for (il in 1: n_pop){
    boot_PI_female[[il]]<- female.boot_err[[il]] + fore_cp_female[[il]]
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

interval_score_cp <- function(PI_val, data_series, fh, alpha = 0.8)
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
      dat.train[[i]]<-extract.years(data_series[[i]], years = 1973:(2022+io))
    }
    dum_pointwise = find_enlarge_val_cp(data_series = dat.train, fh = fh, n_pop = length(data_series))
    PI_boot_cp_Hokkaido[,,io] = dum_pointwise$boot_sample_female[[1]]
    PI_boot_cp_Aomori[,,io] = dum_pointwise$boot_sample_female[[2]]
    PI_boot_cp_Iwate[,,io] = dum_pointwise$boot_sample_female[[3]]
    PI_boot_cp_Miyagi[,,io] = dum_pointwise$boot_sample_female[[4]]
    PI_boot_cp_Akita [,,io] = dum_pointwise$boot_sample_female[[5]]
    PI_boot_cp_Yamagata[,,io] = dum_pointwise$boot_sample_female[[6]]
    PI_boot_cp_Fukushima[,,io] = dum_pointwise$boot_sample_female[[7]]
    PI_boot_cp_Ibaraki[,,io] = dum_pointwise$boot_sample_female[[8]]
    PI_boot_cp_Tochigi[,,io] = dum_pointwise$boot_sample_female[[9]]
    PI_boot_cp_Gunma[,,io] = dum_pointwise$boot_sample_female[[10]]
    PI_boot_cp_Saitama[,,io] = dum_pointwise$boot_sample_female[[11]]
    PI_boot_cp_Chiba[,,io] = dum_pointwise$boot_sample_female[[12]]
    PI_boot_cp_Tokyo[,,io] = dum_pointwise$boot_sample_female[[13]]
    PI_boot_cp_Kanagawa[,,io] = dum_pointwise$boot_sample_female[[14]]
    PI_boot_cp_Niigata[,,io] = dum_pointwise$boot_sample_female[[15]]
    PI_boot_cp_Toyama[,,io] = dum_pointwise$boot_sample_female[[16]]
    PI_boot_cp_Ishikawa[,,io] = dum_pointwise$boot_sample_female[[17]]
    PI_boot_cp_Fukui[,,io] = dum_pointwise$boot_sample_female[[18]]
    PI_boot_cp_Yamanashi[,,io] = dum_pointwise$boot_sample_female[[19]]
    PI_boot_cp_Nagano[,,io] = dum_pointwise$boot_sample_female[[20]]
    PI_boot_cp_Gifu[,,io] = dum_pointwise$boot_sample_female[[21]]
    PI_boot_cp_Shizuoka[,,io] = dum_pointwise$boot_sample_female[[22]]
    PI_boot_cp_Aichi[,,io] = dum_pointwise$boot_sample_female[[23]]
    PI_boot_cp_Mie[,,io] = dum_pointwise$boot_sample_female[[24]]
    PI_boot_cp_Shiga[,,io] = dum_pointwise$boot_sample_female[[25]]
    PI_boot_cp_Kyoto[,,io] = dum_pointwise$boot_sample_female[[26]]
    PI_boot_cp_Osaka[,,io] = dum_pointwise$boot_sample_female[[27]]
    PI_boot_cp_Hyogo[,,io] = dum_pointwise$boot_sample_female[[28]]
    PI_boot_cp_Nara[,,io] = dum_pointwise$boot_sample_female[[29]]
    PI_boot_cp_Wakayama[,,io] = dum_pointwise$boot_sample_female[[30]]
    PI_boot_cp_Tottori[,,io] = dum_pointwise$boot_sample_female[[31]]
    PI_boot_cp_Shimane[,,io] = dum_pointwise$boot_sample_female[[32]]
    PI_boot_cp_Okayama[,,io] = dum_pointwise$boot_sample_female[[33]]
    PI_boot_cp_Hiroshima[,,io] = dum_pointwise$boot_sample_female[[34]]
    PI_boot_cp_Yamaguchi[,,io] = dum_pointwise$boot_sample_female[[35]]
    PI_boot_cp_Tokushima[,,io] = dum_pointwise$boot_sample_female[[36]]
    PI_boot_cp_Kagawa[,,io] = dum_pointwise$boot_sample_female[[37]]
    PI_boot_cp_Ehime [,,io] = dum_pointwise$boot_sample_female[[38]]
    PI_boot_cp_Kochi[,,io] = dum_pointwise$boot_sample_female[[39]]
    PI_boot_cp_Fukuoka[,,io] = dum_pointwise$boot_sample_female[[40]]
    PI_boot_cp_Saga[,,io] = dum_pointwise$boot_sample_female[[41]]
    PI_boot_cp_Nagasaki[,,io] = dum_pointwise$boot_sample_female[[42]]
    PI_boot_cp_Kumamoto[,,io] = dum_pointwise$boot_sample_female[[43]]
    PI_boot_cp_Oita[,,io] = dum_pointwise$boot_sample_female[[44]]
    PI_boot_cp_Miyazaki[,,io] = dum_pointwise$boot_sample_female[[45]]
    PI_boot_cp_Kagoshima[,,io] = dum_pointwise$boot_sample_female[[46]]
    PI_boot_cp_Okinawa[,,io] = dum_pointwise$boot_sample_female[[47]]
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

cp_fh_PI = foreach(mm = 1:10, .packages = c("demography","rTensor")) %dopar% PI_fh_cp(data_series = dat.female,  fh = mm)

female_score.cp.Hokkaido = female_score.cp.Aomori = female_score.cp.Iwate = female_score.cp.Miyagi = 
  female_score.cp.Akita = female_score.cp.Yamagata = female_score.cp.Fukushima = female_score.cp.Ibaraki = 
  female_score.cp.Tochigi = female_score.cp.Gunma = female_score.cp.Saitama= female_score.cp.Chiba = 
  female_score.cp.Tokyo= female_score.cp.Kanagawa = female_score.cp.Niigata = female_score.cp.Toyama = 
  female_score.cp.Ishikawa = female_score.cp.Fukui = female_score.cp.Yamanashi = female_score.cp.Nagano = 
  female_score.cp.Gifu = female_score.cp.Shizuoka = female_score.cp.Aichi = female_score.cp.Mie = 
  female_score.cp.Shiga = female_score.cp.Kyoto = female_score.cp.Osaka = female_score.cp.Hyogo = 
  female_score.cp.Nara = female_score.cp.Wakayama = female_score.cp.Tottori = female_score.cp.Shimane = 
  female_score.cp.Okayama = female_score.cp.Hiroshima = female_score.cp.Yamaguchi = female_score.cp.Tokushima = 
  female_score.cp.Kagawa = female_score.cp.Ehime = female_score.cp.Kochi = female_score.cp.Fukuoka = 
  female_score.cp.Saga = female_score.cp.Nagasaki = female_score.cp.Kumamoto = female_score.cp.Oita =   
  female_score.cp.Miyazaki = female_score.cp.Kagoshima = female_score.cp.Okinawa =vector( ,year_horizon)

for(ik in 1:year_horizon)
{
  # me 
  #### Hokkaido###
  female_score.cp.Hokkaido[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Hokkaido,
                                                     data_series = dat.female[[1]], fh = ik)$score_female
  
  
  #### Aomori###
  female_score.cp.Aomori[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Aomori,
                                                   data_series = dat.female[[2]], fh = ik)$score_female
  #### Iwate###
  female_score.cp.Iwate[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Iwate,
                                                  data_series = dat.female[[3]], fh = ik)$score_female
  #### Miyagi###
  female_score.cp.Miyagi[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Miyagi,
                                                   data_series = dat.female[[4]], fh = ik)$score_female
  #### Akita###
  female_score.cp.Akita[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Akita,
                                                  data_series = dat.female[[5]], fh = ik)$score_female
  #### Yamagata###
  female_score.cp.Yamagata[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Yamagata,
                                                     data_series = dat.female[[6]], fh = ik)$score_female
  #### Fukushima###
  female_score.cp.Fukushima[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Fukushima,
                                                      data_series = dat.female[[7]], fh = ik)$score_female
  #### Ibaraki###
  female_score.cp.Ibaraki[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Ibaraki,
                                                    data_series = dat.female[[8]], fh = ik)$score_female
  #### Tochigi###
  female_score.cp.Tochigi[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Tochigi,
                                                    data_series = dat.female[[9]], fh = ik)$score_female
  #### Gunma###
  female_score.cp.Gunma[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Gunma,
                                                  data_series = dat.female[[11]], fh = ik)$score_female
  #### Saitama###
  female_score.cp.Saitama[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Saitama,
                                                    data_series = dat.female[[12]], fh = ik)$score_female
  #### Chiba###
  female_score.cp.Chiba[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Chiba,
                                                  data_series = dat.female[[13]], fh = ik)$score_female
  #### Tokyo###
  female_score.cp.Tokyo[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Tokyo,
                                                  data_series = dat.female[[14]], fh = ik)$score_female
  #### Kanagawa###
  female_score.cp.Kanagawa[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Kanagawa,
                                                     data_series = dat.female[[1]], fh = ik)$score_female
  #### Niigata###
  female_score.cp.Niigata[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Niigata,
                                                    data_series = dat.female[[15]], fh = ik)$score_female
  #### Toyama###
  female_score.cp.Toyama[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Toyama,
                                                   data_series = dat.female[[16]], fh = ik)$score_female
  #### Ishikawa###
  female_score.cp.Ishikawa[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Ishikawa,
                                                     data_series = dat.female[[17]], fh = ik)$score_female
  #### Fukui###
  female_score.cp.Fukui[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Fukui,
                                                  data_series = dat.female[[18]], fh = ik)$score_female
  #### Yamanashi###
  female_score.cp.Yamanashi[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Yamanashi,
                                                      data_series = dat.female[[19]], fh = ik)$score_female
  #### Nagano###
  female_score.cp.Nagano[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Nagano,
                                                   data_series = dat.female[[20]], fh = ik)$score_female
  #### Gifu###
  female_score.cp.Gifu[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Gifu,
                                                 data_series = dat.female[[21]], fh = ik)$score_female
  #### Shizuoka###
  female_score.cp.Shizuoka[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Shizuoka,
                                                     data_series = dat.female[[22]], fh = ik)$score_female
  #### Aichi###
  female_score.cp.Aichi[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Aichi,
                                                  data_series = dat.female[[23]], fh = ik)$score_female
  #### Mie###
  female_score.cp.Mie[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Mie,
                                                data_series = dat.female[[24]], fh = ik)$score_female
  #### Shiga###
  female_score.cp.Shiga[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Shiga,
                                                  data_series = dat.female[[25]], fh = ik)$score_female
  #### Kyoto###
  female_score.cp.Kyoto[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Kyoto,
                                                  data_series = dat.female[[26]], fh = ik)$score_female
  #### Osaka###
  female_score.cp.Osaka[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Osaka,
                                                  data_series = dat.female[[27]], fh = ik)$score_female
  #### Hyogo###
  female_score.cp.Hyogo[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Hyogo,
                                                  data_series = dat.female[[28]], fh = ik)$score_female
  #### Nara###
  female_score.cp.Nara[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Nara,
                                                 data_series = dat.female[[29]], fh = ik)$score_female
  #### Wakayama###
  female_score.cp.Wakayama[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Wakayama,
                                                     data_series = dat.female[[30]], fh = ik)$score_female
  #### Tottori###
  female_score.cp.Tottori[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Tottori,
                                                    data_series = dat.female[[31]], fh = ik)$score_female
  #### Shimane###
  female_score.cp.Shimane[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Shimane,
                                                    data_series = dat.female[[32]], fh = ik)$score_female
  #### Okayama###
  female_score.cp.Okayama[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Okayama,
                                                    data_series = dat.female[[33]], fh = ik)$score_female
  #### Hiroshima###
  female_score.cp.Hiroshima[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Hiroshima,
                                                      data_series = dat.female[[34]], fh = ik)$score_female
  #### Yamaguchi###
  female_score.cp.Yamaguchi[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Yamaguchi,
                                                      data_series = dat.female[[35]], fh = ik)$score_female
  #### Tokushima###
  female_score.cp.Tokushima[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Tokushima,
                                                      data_series = dat.female[[36]], fh = ik)$score_female
  #### Kagawa###
  female_score.cp.Kagawa[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Kagawa,
                                                   data_series = dat.female[[37]], fh = ik)$score_female
  #### Ehime###
  female_score.cp.Ehime[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Ehime,
                                                  data_series = dat.female[[38]], fh = ik)$score_female
  #### Kochi###
  female_score.cp.Kochi[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Kochi,
                                                  data_series = dat.female[[39]], fh = ik)$score_female
  #### Fukuoka###
  female_score.cp.Fukuoka[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Fukuoka,
                                                    data_series = dat.female[[40]], fh = ik)$score_female
  #### Saga###
  female_score.cp.Saga[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Saga,
                                                 data_series = dat.female[[41]], fh = ik)$score_female
  #### Nagasaki###
  female_score.cp.Nagasaki[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Nagasaki,
                                                     data_series = dat.female[[42]], fh = ik)$score_female
  #### Kumamoto###
  female_score.cp.Kumamoto[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Kumamoto,
                                                     data_series = dat.female[[43]], fh = ik)$score_female
  #### Oita###
  female_score.cp.Oita[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Oita,
                                                 data_series = dat.female[[44]], fh = ik)$score_female
  #### Miyazaki###
  female_score.cp.Miyazaki[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Miyazaki,
                                                     data_series = dat.female[[45]], fh = ik)$score_female
  #### Kagoshima###
  female_score.cp.Kagoshima[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Kagoshima,
                                                      data_series = dat.female[[46]], fh = ik)$score_female
  #### Okinawa###
  female_score.cp.Okinawa[ik] = interval_score_cp(PI_val = cp_fh_PI[[ik]]$PI_boot_cp_Okinawa,
                                                    data_series = dat.female[[47]], fh = ik)$score_female
}


cp.jpfemale.score<-cbind(female_score.cp.Hokkaido , female_score.cp.Aomori , female_score.cp.Iwate ,
                          female_score.cp.Miyagi , female_score.cp.Akita , female_score.cp.Yamagata ,
                          female_score.cp.Fukushima , female_score.cp.Ibaraki , 
                          female_score.cp.Tochigi , female_score.cp.Gunma , female_score.cp.Saitama,
                          female_score.cp.Chiba ,  female_score.cp.Tokyo, female_score.cp.Kanagawa ,
                          female_score.cp.Niigata , female_score.cp.Toyama , female_score.cp.Ishikawa ,
                          female_score.cp.Fukui , female_score.cp.Yamanashi , female_score.cp.Nagano , 
                          female_score.cp.Gifu , female_score.cp.Shizuoka , female_score.cp.Aichi ,
                          female_score.cp.Mie ,  female_score.cp.Shiga , female_score.cp.Kyoto ,
                          female_score.cp.Osaka , female_score.cp.Hyogo ,  female_score.cp.Nara ,
                          female_score.cp.Wakayama , female_score.cp.Tottori , female_score.cp.Shimane ,
                          female_score.cp.Okayama , female_score.cp.Hiroshima , 
                          female_score.cp.Yamaguchi ,  female_score.cp.Tokushima , 
                          female_score.cp.Kagawa , female_score.cp.Ehime , female_score.cp.Kochi ,
                          female_score.cp.Fukuoka , female_score.cp.Saga , female_score.cp.Nagasaki ,
                          female_score.cp.Kumamoto , female_score.cp.Oita , female_score.cp.Miyazaki ,
                          female_score.cp.Kagoshima , female_score.cp.Okinawa )

colnames(cp.jpfemale.score)<-c("Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima", "Ibaraki",  "Tochigi", "Gunma", "Saitama",  
                                "Chiba", "Tokyo", "Kanagawa", "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",    
                                "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi",
                                "Tokushima", "Kagawa", "Ehime", "Kochi", "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")
write.csv(cp.jpfemale.score, "cp.jpfemale.score.csv")





