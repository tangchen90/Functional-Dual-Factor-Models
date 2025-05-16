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
find_enlarge_val_tucker <- function(data_series, fh , nboot = 1000, alpha = 0.8, n_pop)
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
    fore_res = tuckerfts(dat=dat.rate.female_h_step, year_horizon = fh, n_pop = length(data_series), ages =data_series[[1]]$age,ranks=c(2,2,2))
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
  
  
  fore_res  = tuckerfts(dat=dat.series, year_horizon = fh, ages =data_series[[1]]$age, n_pop = length(dat.series), ranks=c(2,2,2))
  
  fore_tucker_female <- list()
  for (il in 1: n_pop){
    fore_tucker_female[[il]]<-matrix(rep(log(fore_res[(1:101)+(il-1)*101 ,fh]), nboot), nrow = 101, ncol = nboot)
  }
  
  
  
  boot_PI_female <- list()
  for (il in 1: n_pop){
    boot_PI_female[[il]]<- female.boot_err[[il]] + fore_tucker_female[[il]]
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

interval_score_tucker <- function(PI_val, data_series, fh, alpha = 0.8)
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
    PI_boot_tucker_Hokkaido[,,io] = dum_pointwise$boot_sample_female[[1]]
    PI_boot_tucker_Aomori[,,io] = dum_pointwise$boot_sample_female[[2]]
    PI_boot_tucker_Iwate[,,io] = dum_pointwise$boot_sample_female[[3]]
    PI_boot_tucker_Miyagi[,,io] = dum_pointwise$boot_sample_female[[4]]
    PI_boot_tucker_Akita [,,io] = dum_pointwise$boot_sample_female[[5]]
    PI_boot_tucker_Yamagata[,,io] = dum_pointwise$boot_sample_female[[6]]
    PI_boot_tucker_Fukushima[,,io] = dum_pointwise$boot_sample_female[[7]]
    PI_boot_tucker_Ibaraki[,,io] = dum_pointwise$boot_sample_female[[8]]
    PI_boot_tucker_Tochigi[,,io] = dum_pointwise$boot_sample_female[[9]]
    PI_boot_tucker_Gunma[,,io] = dum_pointwise$boot_sample_female[[10]]
    PI_boot_tucker_Saitama[,,io] = dum_pointwise$boot_sample_female[[11]]
    PI_boot_tucker_Chiba[,,io] = dum_pointwise$boot_sample_female[[12]]
    PI_boot_tucker_Tokyo[,,io] = dum_pointwise$boot_sample_female[[13]]
    PI_boot_tucker_Kanagawa[,,io] = dum_pointwise$boot_sample_female[[14]]
    PI_boot_tucker_Niigata[,,io] = dum_pointwise$boot_sample_female[[15]]
    PI_boot_tucker_Toyama[,,io] = dum_pointwise$boot_sample_female[[16]]
    PI_boot_tucker_Ishikawa[,,io] = dum_pointwise$boot_sample_female[[17]]
    PI_boot_tucker_Fukui[,,io] = dum_pointwise$boot_sample_female[[18]]
    PI_boot_tucker_Yamanashi[,,io] = dum_pointwise$boot_sample_female[[19]]
    PI_boot_tucker_Nagano[,,io] = dum_pointwise$boot_sample_female[[20]]
    PI_boot_tucker_Gifu[,,io] = dum_pointwise$boot_sample_female[[21]]
    PI_boot_tucker_Shizuoka[,,io] = dum_pointwise$boot_sample_female[[22]]
    PI_boot_tucker_Aichi[,,io] = dum_pointwise$boot_sample_female[[23]]
    PI_boot_tucker_Mie[,,io] = dum_pointwise$boot_sample_female[[24]]
    PI_boot_tucker_Shiga[,,io] = dum_pointwise$boot_sample_female[[25]]
    PI_boot_tucker_Kyoto[,,io] = dum_pointwise$boot_sample_female[[26]]
    PI_boot_tucker_Osaka[,,io] = dum_pointwise$boot_sample_female[[27]]
    PI_boot_tucker_Hyogo[,,io] = dum_pointwise$boot_sample_female[[28]]
    PI_boot_tucker_Nara[,,io] = dum_pointwise$boot_sample_female[[29]]
    PI_boot_tucker_Wakayama[,,io] = dum_pointwise$boot_sample_female[[30]]
    PI_boot_tucker_Tottori[,,io] = dum_pointwise$boot_sample_female[[31]]
    PI_boot_tucker_Shimane[,,io] = dum_pointwise$boot_sample_female[[32]]
    PI_boot_tucker_Okayama[,,io] = dum_pointwise$boot_sample_female[[33]]
    PI_boot_tucker_Hiroshima[,,io] = dum_pointwise$boot_sample_female[[34]]
    PI_boot_tucker_Yamaguchi[,,io] = dum_pointwise$boot_sample_female[[35]]
    PI_boot_tucker_Tokushima[,,io] = dum_pointwise$boot_sample_female[[36]]
    PI_boot_tucker_Kagawa[,,io] = dum_pointwise$boot_sample_female[[37]]
    PI_boot_tucker_Ehime [,,io] = dum_pointwise$boot_sample_female[[38]]
    PI_boot_tucker_Kochi[,,io] = dum_pointwise$boot_sample_female[[39]]
    PI_boot_tucker_Fukuoka[,,io] = dum_pointwise$boot_sample_female[[40]]
    PI_boot_tucker_Saga[,,io] = dum_pointwise$boot_sample_female[[41]]
    PI_boot_tucker_Nagasaki[,,io] = dum_pointwise$boot_sample_female[[42]]
    PI_boot_tucker_Kumamoto[,,io] = dum_pointwise$boot_sample_female[[43]]
    PI_boot_tucker_Oita[,,io] = dum_pointwise$boot_sample_female[[44]]
    PI_boot_tucker_Miyazaki[,,io] = dum_pointwise$boot_sample_female[[45]]
    PI_boot_tucker_Kagoshima[,,io] = dum_pointwise$boot_sample_female[[46]]
    PI_boot_tucker_Okinawa[,,io] = dum_pointwise$boot_sample_female[[47]]
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

tucker_fh_PI = foreach(mm = 1:10, .packages = c("demography","rTensor")) %dopar% PI_fh_tucker(data_series = dat.female,  fh = mm)

female_score.tucker.Hokkaido = female_score.tucker.Aomori = female_score.tucker.Iwate = female_score.tucker.Miyagi = 
  female_score.tucker.Akita = female_score.tucker.Yamagata = female_score.tucker.Fukushima = female_score.tucker.Ibaraki = 
  female_score.tucker.Tochigi = female_score.tucker.Gunma = female_score.tucker.Saitama= female_score.tucker.Chiba = 
  female_score.tucker.Tokyo= female_score.tucker.Kanagawa = female_score.tucker.Niigata = female_score.tucker.Toyama = 
  female_score.tucker.Ishikawa = female_score.tucker.Fukui = female_score.tucker.Yamanashi = female_score.tucker.Nagano = 
  female_score.tucker.Gifu = female_score.tucker.Shizuoka = female_score.tucker.Aichi = female_score.tucker.Mie = 
  female_score.tucker.Shiga = female_score.tucker.Kyoto = female_score.tucker.Osaka = female_score.tucker.Hyogo = 
  female_score.tucker.Nara = female_score.tucker.Wakayama = female_score.tucker.Tottori = female_score.tucker.Shimane = 
  female_score.tucker.Okayama = female_score.tucker.Hiroshima = female_score.tucker.Yamaguchi = female_score.tucker.Tokushima = 
  female_score.tucker.Kagawa = female_score.tucker.Ehime = female_score.tucker.Kochi = female_score.tucker.Fukuoka = 
  female_score.tucker.Saga = female_score.tucker.Nagasaki = female_score.tucker.Kumamoto = female_score.tucker.Oita =   
  female_score.tucker.Miyazaki = female_score.tucker.Kagoshima = female_score.tucker.Okinawa =vector( ,year_horizon)

for(ik in 1:year_horizon)
{
  # me 
  #### Hokkaido###
  female_score.tucker.Hokkaido[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Hokkaido,
                                                   data_series = dat.female[[1]], fh = ik)$score_female
  
  
  #### Aomori###
  female_score.tucker.Aomori[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Aomori,
                                                 data_series = dat.female[[2]], fh = ik)$score_female
  #### Iwate###
  female_score.tucker.Iwate[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Iwate,
                                                data_series = dat.female[[3]], fh = ik)$score_female
  #### Miyagi###
  female_score.tucker.Miyagi[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Miyagi,
                                                 data_series = dat.female[[4]], fh = ik)$score_female
  #### Akita###
  female_score.tucker.Akita[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Akita,
                                                data_series = dat.female[[5]], fh = ik)$score_female
  #### Yamagata###
  female_score.tucker.Yamagata[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Yamagata,
                                                   data_series = dat.female[[6]], fh = ik)$score_female
  #### Fukushima###
  female_score.tucker.Fukushima[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Fukushima,
                                                    data_series = dat.female[[7]], fh = ik)$score_female
  #### Ibaraki###
  female_score.tucker.Ibaraki[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Ibaraki,
                                                  data_series = dat.female[[8]], fh = ik)$score_female
  #### Tochigi###
  female_score.tucker.Tochigi[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Tochigi,
                                                  data_series = dat.female[[9]], fh = ik)$score_female
  #### Gunma###
  female_score.tucker.Gunma[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Gunma,
                                                data_series = dat.female[[11]], fh = ik)$score_female
  #### Saitama###
  female_score.tucker.Saitama[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Saitama,
                                                  data_series = dat.female[[12]], fh = ik)$score_female
  #### Chiba###
  female_score.tucker.Chiba[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Chiba,
                                                data_series = dat.female[[13]], fh = ik)$score_female
  #### Tokyo###
  female_score.tucker.Tokyo[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Tokyo,
                                                data_series = dat.female[[14]], fh = ik)$score_female
  #### Kanagawa###
  female_score.tucker.Kanagawa[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Kanagawa,
                                                   data_series = dat.female[[1]], fh = ik)$score_female
  #### Niigata###
  female_score.tucker.Niigata[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Niigata,
                                                  data_series = dat.female[[15]], fh = ik)$score_female
  #### Toyama###
  female_score.tucker.Toyama[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Toyama,
                                                 data_series = dat.female[[16]], fh = ik)$score_female
  #### Ishikawa###
  female_score.tucker.Ishikawa[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Ishikawa,
                                                   data_series = dat.female[[17]], fh = ik)$score_female
  #### Fukui###
  female_score.tucker.Fukui[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Fukui,
                                                data_series = dat.female[[18]], fh = ik)$score_female
  #### Yamanashi###
  female_score.tucker.Yamanashi[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Yamanashi,
                                                    data_series = dat.female[[19]], fh = ik)$score_female
  #### Nagano###
  female_score.tucker.Nagano[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Nagano,
                                                 data_series = dat.female[[20]], fh = ik)$score_female
  #### Gifu###
  female_score.tucker.Gifu[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Gifu,
                                               data_series = dat.female[[21]], fh = ik)$score_female
  #### Shizuoka###
  female_score.tucker.Shizuoka[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Shizuoka,
                                                   data_series = dat.female[[22]], fh = ik)$score_female
  #### Aichi###
  female_score.tucker.Aichi[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Aichi,
                                                data_series = dat.female[[23]], fh = ik)$score_female
  #### Mie###
  female_score.tucker.Mie[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Mie,
                                              data_series = dat.female[[24]], fh = ik)$score_female
  #### Shiga###
  female_score.tucker.Shiga[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Shiga,
                                                data_series = dat.female[[25]], fh = ik)$score_female
  #### Kyoto###
  female_score.tucker.Kyoto[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Kyoto,
                                                data_series = dat.female[[26]], fh = ik)$score_female
  #### Osaka###
  female_score.tucker.Osaka[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Osaka,
                                                data_series = dat.female[[27]], fh = ik)$score_female
  #### Hyogo###
  female_score.tucker.Hyogo[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Hyogo,
                                                data_series = dat.female[[28]], fh = ik)$score_female
  #### Nara###
  female_score.tucker.Nara[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Nara,
                                               data_series = dat.female[[29]], fh = ik)$score_female
  #### Wakayama###
  female_score.tucker.Wakayama[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Wakayama,
                                                   data_series = dat.female[[30]], fh = ik)$score_female
  #### Tottori###
  female_score.tucker.Tottori[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Tottori,
                                                  data_series = dat.female[[31]], fh = ik)$score_female
  #### Shimane###
  female_score.tucker.Shimane[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Shimane,
                                                  data_series = dat.female[[32]], fh = ik)$score_female
  #### Okayama###
  female_score.tucker.Okayama[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Okayama,
                                                  data_series = dat.female[[33]], fh = ik)$score_female
  #### Hiroshima###
  female_score.tucker.Hiroshima[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Hiroshima,
                                                    data_series = dat.female[[34]], fh = ik)$score_female
  #### Yamaguchi###
  female_score.tucker.Yamaguchi[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Yamaguchi,
                                                    data_series = dat.female[[35]], fh = ik)$score_female
  #### Tokushima###
  female_score.tucker.Tokushima[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Tokushima,
                                                    data_series = dat.female[[36]], fh = ik)$score_female
  #### Kagawa###
  female_score.tucker.Kagawa[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Kagawa,
                                                 data_series = dat.female[[37]], fh = ik)$score_female
  #### Ehime###
  female_score.tucker.Ehime[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Ehime,
                                                data_series = dat.female[[38]], fh = ik)$score_female
  #### Kochi###
  female_score.tucker.Kochi[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Kochi,
                                                data_series = dat.female[[39]], fh = ik)$score_female
  #### Fukuoka###
  female_score.tucker.Fukuoka[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Fukuoka,
                                                  data_series = dat.female[[40]], fh = ik)$score_female
  #### Saga###
  female_score.tucker.Saga[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Saga,
                                               data_series = dat.female[[41]], fh = ik)$score_female
  #### Nagasaki###
  female_score.tucker.Nagasaki[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Nagasaki,
                                                   data_series = dat.female[[42]], fh = ik)$score_female
  #### Kumamoto###
  female_score.tucker.Kumamoto[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Kumamoto,
                                                   data_series = dat.female[[43]], fh = ik)$score_female
  #### Oita###
  female_score.tucker.Oita[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Oita,
                                               data_series = dat.female[[44]], fh = ik)$score_female
  #### Miyazaki###
  female_score.tucker.Miyazaki[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Miyazaki,
                                                   data_series = dat.female[[45]], fh = ik)$score_female
  #### Kagoshima###
  female_score.tucker.Kagoshima[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Kagoshima,
                                                    data_series = dat.female[[46]], fh = ik)$score_female
  #### Okinawa###
  female_score.tucker.Okinawa[ik] = interval_score_tucker(PI_val = tucker_fh_PI[[ik]]$PI_boot_tucker_Okinawa,
                                                  data_series = dat.female[[47]], fh = ik)$score_female
}


tucker.jpfemale.score<-cbind(female_score.tucker.Hokkaido , female_score.tucker.Aomori , female_score.tucker.Iwate ,
                         female_score.tucker.Miyagi , female_score.tucker.Akita , female_score.tucker.Yamagata ,
                         female_score.tucker.Fukushima , female_score.tucker.Ibaraki , 
                         female_score.tucker.Tochigi , female_score.tucker.Gunma , female_score.tucker.Saitama,
                         female_score.tucker.Chiba ,  female_score.tucker.Tokyo, female_score.tucker.Kanagawa ,
                         female_score.tucker.Niigata , female_score.tucker.Toyama , female_score.tucker.Ishikawa ,
                         female_score.tucker.Fukui , female_score.tucker.Yamanashi , female_score.tucker.Nagano , 
                         female_score.tucker.Gifu , female_score.tucker.Shizuoka , female_score.tucker.Aichi ,
                         female_score.tucker.Mie ,  female_score.tucker.Shiga , female_score.tucker.Kyoto ,
                         female_score.tucker.Osaka , female_score.tucker.Hyogo ,  female_score.tucker.Nara ,
                         female_score.tucker.Wakayama , female_score.tucker.Tottori , female_score.tucker.Shimane ,
                         female_score.tucker.Okayama , female_score.tucker.Hiroshima , 
                         female_score.tucker.Yamaguchi ,  female_score.tucker.Tokushima , 
                         female_score.tucker.Kagawa , female_score.tucker.Ehime , female_score.tucker.Kochi ,
                         female_score.tucker.Fukuoka , female_score.tucker.Saga , female_score.tucker.Nagasaki ,
                         female_score.tucker.Kumamoto , female_score.tucker.Oita , female_score.tucker.Miyazaki ,
                         female_score.tucker.Kagoshima , female_score.tucker.Okinawa )

colnames(tucker.jpfemale.score)<-c("Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima", "Ibaraki",  "Tochigi", "Gunma", "Saitama",  
                               "Chiba", "Tokyo", "Kanagawa", "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",    
                               "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi",
                               "Tokushima", "Kagawa", "Ehime", "Kochi", "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")
write.csv(tucker.jpfemale.score, "tucker.jpfemale.score.csv")





