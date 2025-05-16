library(doParallel)
library(rTensor)
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
tuckerfts<- function(dat, year_horizon, years, ages, n_pop, ranks)
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
find_enlarge_val_tucker <- function(data_series, fh , nboot = 1000, alpha = 0.8, n_pop)
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
    fore_res = tuckerfts(dat=dat.rate.male_h_step, year_horizon = fh, n_pop = length(data_series), ages =data_series[[1]]$age,ranks=c(2,2,2))
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
  
  
  fore_res  = tuckerfts(dat=dat.series, year_horizon = fh, ages =data_series[[1]]$age, n_pop = length(dat.series), ranks=c(2,2,2))
  
  fore_tucker_male <- list()
  for (il in 1: n_pop){
    fore_tucker_male[[il]]<-matrix(rep(log(fore_res[(1:101)+(il-1)*101 ,fh]), nboot), nrow = 101, ncol = nboot)
  }
  
  
  
  boot_PI_male <- list()
  for (il in 1: n_pop){
    boot_PI_male[[il]]<- male.boot_err[[il]] + fore_tucker_male[[il]]
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

interval_score_tucker <- function(PI_val, data_series, fh, alpha = 0.8)
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
PI_fh_tucker<- function(data_series, fh, nboot = 1000)
{
  PI_boot_tucker_Hokkaido = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Aomori = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Iwate = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Miyagi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Akita  = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Yamagata = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Fukushima = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Ibaraki = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Tochigi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Gunma = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Saitama = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Chiba = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Tokyo = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Kanagawa = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Niigata = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Toyama = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Ishikawa = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Fukui = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Yamanashi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Nagano = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Gifu = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Shizuoka = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Aichi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Mie = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Shiga = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Kyoto = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Osaka = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Hyogo = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Nara = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Wakayama = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Tottori = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Shimane = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Okayama = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Hiroshima = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Yamaguchi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Tokushima = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Kagawa = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Ehime  = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Kochi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Fukuoka = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Saga = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Nagasaki = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Kumamoto = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Oita = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Miyazaki = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Kagoshima = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_tucker_Okinawa = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  
  for(io in 1:(11-fh))
  {
    dat.train<-list()
    for ( i in 1:length(jpn.data)){
      dat.train[[i]]<-extract.years(data_series[[i]], years = 1973:(2012+io))
    }
    dum_pointwise = find_enlarge_val_tucker(data_series = dat.train, fh = fh, n_pop = length(data_series))
    PI_boot_tucker_Hokkaido[,,io] = dum_pointwise$boot_sample_male[[1]]
    PI_boot_tucker_Aomori[,,io] = dum_pointwise$boot_sample_male[[2]]
    PI_boot_tucker_Iwate[,,io] = dum_pointwise$boot_sample_male[[3]]
    PI_boot_tucker_Miyagi[,,io] = dum_pointwise$boot_sample_male[[4]]
    PI_boot_tucker_Akita [,,io] = dum_pointwise$boot_sample_male[[5]]
    PI_boot_tucker_Yamagata[,,io] = dum_pointwise$boot_sample_male[[6]]
    PI_boot_tucker_Fukushima[,,io] = dum_pointwise$boot_sample_male[[7]]
    PI_boot_tucker_Ibaraki[,,io] = dum_pointwise$boot_sample_male[[8]]
    PI_boot_tucker_Tochigi[,,io] = dum_pointwise$boot_sample_male[[9]]
    PI_boot_tucker_Gunma[,,io] = dum_pointwise$boot_sample_male[[10]]
    PI_boot_tucker_Saitama[,,io] = dum_pointwise$boot_sample_male[[11]]
    PI_boot_tucker_Chiba[,,io] = dum_pointwise$boot_sample_male[[12]]
    PI_boot_tucker_Tokyo[,,io] = dum_pointwise$boot_sample_male[[13]]
    PI_boot_tucker_Kanagawa[,,io] = dum_pointwise$boot_sample_male[[14]]
    PI_boot_tucker_Niigata[,,io] = dum_pointwise$boot_sample_male[[15]]
    PI_boot_tucker_Toyama[,,io] = dum_pointwise$boot_sample_male[[16]]
    PI_boot_tucker_Ishikawa[,,io] = dum_pointwise$boot_sample_male[[17]]
    PI_boot_tucker_Fukui[,,io] = dum_pointwise$boot_sample_male[[18]]
    PI_boot_tucker_Yamanashi[,,io] = dum_pointwise$boot_sample_male[[19]]
    PI_boot_tucker_Nagano[,,io] = dum_pointwise$boot_sample_male[[20]]
    PI_boot_tucker_Gifu[,,io] = dum_pointwise$boot_sample_male[[21]]
    PI_boot_tucker_Shizuoka[,,io] = dum_pointwise$boot_sample_male[[22]]
    PI_boot_tucker_Aichi[,,io] = dum_pointwise$boot_sample_male[[23]]
    PI_boot_tucker_Mie[,,io] = dum_pointwise$boot_sample_male[[24]]
    PI_boot_tucker_Shiga[,,io] = dum_pointwise$boot_sample_male[[25]]
    PI_boot_tucker_Kyoto[,,io] = dum_pointwise$boot_sample_male[[26]]
    PI_boot_tucker_Osaka[,,io] = dum_pointwise$boot_sample_male[[27]]
    PI_boot_tucker_Hyogo[,,io] = dum_pointwise$boot_sample_male[[28]]
    PI_boot_tucker_Nara[,,io] = dum_pointwise$boot_sample_male[[29]]
    PI_boot_tucker_Wakayama[,,io] = dum_pointwise$boot_sample_male[[30]]
    PI_boot_tucker_Tottori[,,io] = dum_pointwise$boot_sample_male[[31]]
    PI_boot_tucker_Shimane[,,io] = dum_pointwise$boot_sample_male[[32]]
    PI_boot_tucker_Okayama[,,io] = dum_pointwise$boot_sample_male[[33]]
    PI_boot_tucker_Hiroshima[,,io] = dum_pointwise$boot_sample_male[[34]]
    PI_boot_tucker_Yamaguchi[,,io] = dum_pointwise$boot_sample_male[[35]]
    PI_boot_tucker_Tokushima[,,io] = dum_pointwise$boot_sample_male[[36]]
    PI_boot_tucker_Kagawa[,,io] = dum_pointwise$boot_sample_male[[37]]
    PI_boot_tucker_Ehime [,,io] = dum_pointwise$boot_sample_male[[38]]
    PI_boot_tucker_Kochi[,,io] = dum_pointwise$boot_sample_male[[39]]
    PI_boot_tucker_Fukuoka[,,io] = dum_pointwise$boot_sample_male[[40]]
    PI_boot_tucker_Saga[,,io] = dum_pointwise$boot_sample_male[[41]]
    PI_boot_tucker_Nagasaki[,,io] = dum_pointwise$boot_sample_male[[42]]
    PI_boot_tucker_Kumamoto[,,io] = dum_pointwise$boot_sample_male[[43]]
    PI_boot_tucker_Oita[,,io] = dum_pointwise$boot_sample_male[[44]]
    PI_boot_tucker_Miyazaki[,,io] = dum_pointwise$boot_sample_male[[45]]
    PI_boot_tucker_Kagoshima[,,io] = dum_pointwise$boot_sample_male[[46]]
    PI_boot_tucker_Okinawa[,,io] = dum_pointwise$boot_sample_male[[47]]
  }
  return(list(PI_boot_tucker_Hokkaido = PI_boot_tucker_Hokkaido,PI_boot_tucker_Aomori = PI_boot_tucker_Aomori,
              PI_boot_tucker_Iwate = PI_boot_tucker_Iwate,PI_boot_tucker_Miyagi = PI_boot_tucker_Miyagi,
              PI_boot_tucker_Akita  = PI_boot_tucker_Akita ,PI_boot_tucker_Yamagata = PI_boot_tucker_Yamagata,
              PI_boot_tucker_Fukushima = PI_boot_tucker_Fukushima, PI_boot_tucker_Ibaraki = PI_boot_tucker_Ibaraki,
              PI_boot_tucker_Tochigi = PI_boot_tucker_Tochigi, PI_boot_tucker_Gunma = PI_boot_tucker_Gunma,
              PI_boot_tucker_Saitama = PI_boot_tucker_Saitama, PI_boot_tucker_Chiba = PI_boot_tucker_Chiba, 
              PI_boot_tucker_Tokyo = PI_boot_tucker_Tokyo, PI_boot_tucker_Kanagawa = PI_boot_tucker_Kanagawa,
              PI_boot_tucker_Niigata = PI_boot_tucker_Niigata, PI_boot_tucker_Toyama = PI_boot_tucker_Toyama,
              PI_boot_tucker_Ishikawa = PI_boot_tucker_Ishikawa,PI_boot_tucker_Fukui = PI_boot_tucker_Fukui, 
              PI_boot_tucker_Yamanashi = PI_boot_tucker_Yamanashi,PI_boot_tucker_Nagano = PI_boot_tucker_Nagano,
              PI_boot_tucker_Gifu = PI_boot_tucker_Gifu,PI_boot_tucker_Shizuoka = PI_boot_tucker_Shizuoka,
              PI_boot_tucker_Aichi = PI_boot_tucker_Aichi,PI_boot_tucker_Mie = PI_boot_tucker_Mie,
              PI_boot_tucker_Shiga = PI_boot_tucker_Shiga,PI_boot_tucker_Kyoto = PI_boot_tucker_Kyoto,
              PI_boot_tucker_Osaka = PI_boot_tucker_Osaka, PI_boot_tucker_Hyogo = PI_boot_tucker_Hyogo,
              PI_boot_tucker_Nara = PI_boot_tucker_Nara, PI_boot_tucker_Wakayama = PI_boot_tucker_Wakayama, 
              PI_boot_tucker_Tottori = PI_boot_tucker_Tottori,PI_boot_tucker_Shimane = PI_boot_tucker_Shimane,
              PI_boot_tucker_Okayama = PI_boot_tucker_Okayama,PI_boot_tucker_Hiroshima = PI_boot_tucker_Hiroshima,
              PI_boot_tucker_Yamaguchi = PI_boot_tucker_Yamaguchi, PI_boot_tucker_Tokushima = PI_boot_tucker_Tokushima,
              PI_boot_tucker_Kagawa = PI_boot_tucker_Kagawa, PI_boot_tucker_Ehime  = PI_boot_tucker_Ehime, 
              PI_boot_tucker_Kochi = PI_boot_tucker_Kochi, PI_boot_tucker_Fukuoka = PI_boot_tucker_Fukuoka,
              PI_boot_tucker_Saga = PI_boot_tucker_Saga,PI_boot_tucker_Nagasaki = PI_boot_tucker_Nagasaki,
              PI_boot_tucker_Kumamoto = PI_boot_tucker_Kumamoto,PI_boot_tucker_Oita = PI_boot_tucker_Oita,
              PI_boot_tucker_Miyazaki = PI_boot_tucker_Miyazaki,PI_boot_tucker_Kagoshima = PI_boot_tucker_Kagoshima,
              PI_boot_tucker_Okinawa = PI_boot_tucker_Okinawa))
}
library(doParallel)
cl <- makeCluster(4) 
registerDoParallel(cl)

tucker_fh_PI = foreach(mm = 1:10, .packages = c("demography","rTensor")) %dopar% PI_fh_tucker(data_series = dat.male,  fh = mm)

male_score.tucker.Hokkaido = male_score.tucker.Aomori = male_score.tucker.Iwate = male_score.tucker.Miyagi = 
  male_score.tucker.Akita = male_score.tucker.Yamagata = male_score.tucker.Fukushima = male_score.tucker.Ibaraki = 
  male_score.tucker.Tochigi = male_score.tucker.Gunma = male_score.tucker.Saitama= male_score.tucker.Chiba = 
  male_score.tucker.Tokyo= male_score.tucker.Kanagawa = male_score.tucker.Niigata = male_score.tucker.Toyama = 
  male_score.tucker.Ishikawa = male_score.tucker.Fukui = male_score.tucker.Yamanashi = male_score.tucker.Nagano = 
  male_score.tucker.Gifu = male_score.tucker.Shizuoka = male_score.tucker.Aichi = male_score.tucker.Mie = 
  male_score.tucker.Shiga = male_score.tucker.Kyoto = male_score.tucker.Osaka = male_score.tucker.Hyogo = 
  male_score.tucker.Nara = male_score.tucker.Wakayama = male_score.tucker.Tottori = male_score.tucker.Shimane = 
  male_score.tucker.Okayama = male_score.tucker.Hiroshima = male_score.tucker.Yamaguchi = male_score.tucker.Tokushima = 
  male_score.tucker.Kagawa = male_score.tucker.Ehime = male_score.tucker.Kochi = male_score.tucker.Fukuoka = 
  male_score.tucker.Saga = male_score.tucker.Nagasaki = male_score.tucker.Kumamoto = male_score.tucker.Oita =   
  male_score.tucker.Miyazaki = male_score.tucker.Kagoshima = male_score.tucker.Okinawa =vector( ,year_horizon)

for(ik in 1:year_horizon)
{
  # me 
  #### Hokkaido###
  male_score.tucker.Hokkaido[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Hokkaido,
                                                   data_series = dat.male[[1]], fh = ik)$score_male
  
  
  #### Aomori###
  male_score.tucker.Aomori[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Aomori,
                                                 data_series = dat.male[[2]], fh = ik)$score_male
  #### Iwate###
  male_score.tucker.Iwate[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Iwate,
                                                data_series = dat.male[[3]], fh = ik)$score_male
  #### Miyagi###
  male_score.tucker.Miyagi[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Miyagi,
                                                 data_series = dat.male[[4]], fh = ik)$score_male
  #### Akita###
  male_score.tucker.Akita[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Akita,
                                                data_series = dat.male[[5]], fh = ik)$score_male
  #### Yamagata###
  male_score.tucker.Yamagata[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Yamagata,
                                                   data_series = dat.male[[6]], fh = ik)$score_male
  #### Fukushima###
  male_score.tucker.Fukushima[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Fukushima,
                                                    data_series = dat.male[[7]], fh = ik)$score_male
  #### Ibaraki###
  male_score.tucker.Ibaraki[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Ibaraki,
                                                  data_series = dat.male[[8]], fh = ik)$score_male
  #### Tochigi###
  male_score.tucker.Tochigi[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Tochigi,
                                                  data_series = dat.male[[9]], fh = ik)$score_male
  #### Gunma###
  male_score.tucker.Gunma[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Gunma,
                                                data_series = dat.male[[11]], fh = ik)$score_male
  #### Saitama###
  male_score.tucker.Saitama[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Saitama,
                                                  data_series = dat.male[[12]], fh = ik)$score_male
  #### Chiba###
  male_score.tucker.Chiba[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Chiba,
                                                data_series = dat.male[[13]], fh = ik)$score_male
  #### Tokyo###
  male_score.tucker.Tokyo[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Tokyo,
                                                data_series = dat.male[[14]], fh = ik)$score_male
  #### Kanagawa###
  male_score.tucker.Kanagawa[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Kanagawa,
                                                   data_series = dat.male[[1]], fh = ik)$score_male
  #### Niigata###
  male_score.tucker.Niigata[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Niigata,
                                                  data_series = dat.male[[15]], fh = ik)$score_male
  #### Toyama###
  male_score.tucker.Toyama[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Toyama,
                                                 data_series = dat.male[[16]], fh = ik)$score_male
  #### Ishikawa###
  male_score.tucker.Ishikawa[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Ishikawa,
                                                   data_series = dat.male[[17]], fh = ik)$score_male
  #### Fukui###
  male_score.tucker.Fukui[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Fukui,
                                                data_series = dat.male[[18]], fh = ik)$score_male
  #### Yamanashi###
  male_score.tucker.Yamanashi[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Yamanashi,
                                                    data_series = dat.male[[19]], fh = ik)$score_male
  #### Nagano###
  male_score.tucker.Nagano[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Nagano,
                                                 data_series = dat.male[[20]], fh = ik)$score_male
  #### Gifu###
  male_score.tucker.Gifu[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Gifu,
                                               data_series = dat.male[[21]], fh = ik)$score_male
  #### Shizuoka###
  male_score.tucker.Shizuoka[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Shizuoka,
                                                   data_series = dat.male[[22]], fh = ik)$score_male
  #### Aichi###
  male_score.tucker.Aichi[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Aichi,
                                                data_series = dat.male[[23]], fh = ik)$score_male
  #### Mie###
  male_score.tucker.Mie[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Mie,
                                              data_series = dat.male[[24]], fh = ik)$score_male
  #### Shiga###
  male_score.tucker.Shiga[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Shiga,
                                                data_series = dat.male[[25]], fh = ik)$score_male
  #### Kyoto###
  male_score.tucker.Kyoto[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Kyoto,
                                                data_series = dat.male[[26]], fh = ik)$score_male
  #### Osaka###
  male_score.tucker.Osaka[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Osaka,
                                                data_series = dat.male[[27]], fh = ik)$score_male
  #### Hyogo###
  male_score.tucker.Hyogo[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Hyogo,
                                                data_series = dat.male[[28]], fh = ik)$score_male
  #### Nara###
  male_score.tucker.Nara[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Nara,
                                               data_series = dat.male[[29]], fh = ik)$score_male
  #### Wakayama###
  male_score.tucker.Wakayama[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Wakayama,
                                                   data_series = dat.male[[30]], fh = ik)$score_male
  #### Tottori###
  male_score.tucker.Tottori[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Tottori,
                                                  data_series = dat.male[[31]], fh = ik)$score_male
  #### Shimane###
  male_score.tucker.Shimane[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Shimane,
                                                  data_series = dat.male[[32]], fh = ik)$score_male
  #### Okayama###
  male_score.tucker.Okayama[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Okayama,
                                                  data_series = dat.male[[33]], fh = ik)$score_male
  #### Hiroshima###
  male_score.tucker.Hiroshima[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Hiroshima,
                                                    data_series = dat.male[[34]], fh = ik)$score_male
  #### Yamaguchi###
  male_score.tucker.Yamaguchi[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Yamaguchi,
                                                    data_series = dat.male[[35]], fh = ik)$score_male
  #### Tokushima###
  male_score.tucker.Tokushima[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Tokushima,
                                                    data_series = dat.male[[36]], fh = ik)$score_male
  #### Kagawa###
  male_score.tucker.Kagawa[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Kagawa,
                                                 data_series = dat.male[[37]], fh = ik)$score_male
  #### Ehime###
  male_score.tucker.Ehime[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Ehime,
                                                data_series = dat.male[[38]], fh = ik)$score_male
  #### Kochi###
  male_score.tucker.Kochi[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Kochi,
                                                data_series = dat.male[[39]], fh = ik)$score_male
  #### Fukuoka###
  male_score.tucker.Fukuoka[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Fukuoka,
                                                  data_series = dat.male[[40]], fh = ik)$score_male
  #### Saga###
  male_score.tucker.Saga[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Saga,
                                               data_series = dat.male[[41]], fh = ik)$score_male
  #### Nagasaki###
  male_score.tucker.Nagasaki[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Nagasaki,
                                                   data_series = dat.male[[42]], fh = ik)$score_male
  #### Kumamoto###
  male_score.tucker.Kumamoto[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Kumamoto,
                                                   data_series = dat.male[[43]], fh = ik)$score_male
  #### Oita###
  male_score.tucker.Oita[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Oita,
                                               data_series = dat.male[[44]], fh = ik)$score_male
  #### Miyazaki###
  male_score.tucker.Miyazaki[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Miyazaki,
                                                   data_series = dat.male[[45]], fh = ik)$score_male
  #### Kagoshima###
  male_score.tucker.Kagoshima[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Kagoshima,
                                                    data_series = dat.male[[46]], fh = ik)$score_male
  #### Okinawa###
  male_score.tucker.Okinawa[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Okinawa,
                                                  data_series = dat.male[[47]], fh = ik)$score_male
}


tucker.jpmale.score<-cbind(male_score.tucker.Hokkaido , male_score.tucker.Aomori , male_score.tucker.Iwate ,
                         male_score.tucker.Miyagi , male_score.tucker.Akita , male_score.tucker.Yamagata ,
                         male_score.tucker.Fukushima , male_score.tucker.Ibaraki , 
                         male_score.tucker.Tochigi , male_score.tucker.Gunma , male_score.tucker.Saitama,
                         male_score.tucker.Chiba ,  male_score.tucker.Tokyo, male_score.tucker.Kanagawa ,
                         male_score.tucker.Niigata , male_score.tucker.Toyama , male_score.tucker.Ishikawa ,
                         male_score.tucker.Fukui , male_score.tucker.Yamanashi , male_score.tucker.Nagano , 
                         male_score.tucker.Gifu , male_score.tucker.Shizuoka , male_score.tucker.Aichi ,
                         male_score.tucker.Mie ,  male_score.tucker.Shiga , male_score.tucker.Kyoto ,
                         male_score.tucker.Osaka , male_score.tucker.Hyogo ,  male_score.tucker.Nara ,
                         male_score.tucker.Wakayama , male_score.tucker.Tottori , male_score.tucker.Shimane ,
                         male_score.tucker.Okayama , male_score.tucker.Hiroshima , 
                         male_score.tucker.Yamaguchi ,  male_score.tucker.Tokushima , 
                         male_score.tucker.Kagawa , male_score.tucker.Ehime , male_score.tucker.Kochi ,
                         male_score.tucker.Fukuoka , male_score.tucker.Saga , male_score.tucker.Nagasaki ,
                         male_score.tucker.Kumamoto , male_score.tucker.Oita , male_score.tucker.Miyazaki ,
                         male_score.tucker.Kagoshima , male_score.tucker.Okinawa )

colnames(tucker.jpmale.score)<-c("Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima", "Ibaraki",  "Tochigi", "Gunma", "Saitama",  
                               "Chiba", "Tokyo", "Kanagawa", "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",    
                               "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi",
                               "Tokushima", "Kagawa", "Ehime", "Kochi", "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")


write.csv(tucker.jpmale.score, "tucker.jpmale.score.csv")





