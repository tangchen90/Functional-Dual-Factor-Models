#### Function for high-dimensional functional time series estimation####
###Input is a list of FTS
hdfts<- function(X)
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
    Gamma_l <- long_run_covariance_estimation(t(dat.center[[ik]]))
    FPCA.dym <- eigen(Gamma_l)
    dfpca.value <- FPCA.dym$values
    dfpca.value <- ifelse(dfpca.value>=0, dfpca.value, 0)
    dfpca.vector[[ik]]<-FPCA.dym$vectors
    percent_var <- (dfpca.value)/sum(dfpca.value)
    ratio <- dfpca.value[1]/dfpca.value
    num.component[ik] <-  min(which(cumsum(percent_var) > 0.99))
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
hdftsfts<- function(dat, year_horizon, years, ages, n_pop)
{ 
  hdfts.dat<-hdfts(dat)
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
  if (year_horizon==1){
    X.forc<-list()
    for (t in 1:length(dat)){
      X.forc[[t]]<-t(as.matrix(hdfts.dat$mu[[t]][1:year_horizon,]))
      for (s in 1:year_horizon){
        for (m in 1:length(Ft.dat)){
          X.forc[[t]][s,]<-X.forc[[t]][s,]+ beta.forc[[m]][s,t]*hdfts.dat$phi[[t]][,m]
        }
      }
    }
  }else{
    X.forc<-list()
    for (t in 1:length(dat)){
      X.forc[[t]]<-hdfts.dat$mu[[t]][1:year_horizon,]
      for (s in 1:year_horizon){
        for (m in 1:length(Ft.dat)){
          X.forc[[t]][s,]<-X.forc[[t]][s,]+ beta.forc[[m]][s,t]*hdfts.dat$phi[[t]][,m]
        }
      }
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
find_enlarge_val_hdfts <- function(data_series, fh , nboot = 1000, alpha = 0.8, n_pop)
{
  n_test<-25
  
  # calculate  in-sample forecast curves
  
  fore_curve_female  = matrix( , length(data_series[[1]]$age)*n_pop, length(data_series[[1]]$year) - n_test - fh )
  for(ij in 1:(length(data_series[[1]]$year) - n_test - fh))
  {
    dat.rate.female_h_step<-list()
    for (ik in 1:n_pop){
      dat.rate.female_h_step[[ik]]<-t(log(extract.years(data_series[[ik]], years = 1973:(1973 + n_test + ij - 1))$rate$female))
    }
    fore_res = hdftsfts(dat=dat.rate.female_h_step, year_horizon = fh, n_pop = length(data_series), ages =data_series[[1]]$age)
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
  
  
  fore_res  = hdftsfts(dat=dat.series, year_horizon = fh, n_pop = length(dat.series), ages =data_series[[1]]$age)
  
  fore_hdfts_female <- list()
  for (il in 1: n_pop){
    fore_hdfts_female[[il]]<-matrix(rep(log(fore_res[(1:101)+(il-1)*101 ,fh]), nboot), nrow = 101, ncol = nboot)
  }
  
  
  
  boot_PI_female <- list()
  for (il in 1: n_pop){
    boot_PI_female[[il]]<- female.boot_err[[il]] + fore_hdfts_female[[il]]
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

interval_score_hdfts <- function(PI_val, data_series, fh, alpha = 0.8)
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
PI_fh_hdfts<- function(data_series, fh, nboot = 1000)
{
  PI_boot_hdfts_Hokkaido = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Aomori = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Iwate = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Miyagi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Akita  = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Yamagata = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Fukushima = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Ibaraki = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Tochigi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Gunma = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Saitama = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Chiba = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Tokyo = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Kanagawa = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Niigata = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Toyama = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Ishikawa = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Fukui = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Yamanashi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Nagano = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Gifu = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Shizuoka = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Aichi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Mie = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Shiga = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Kyoto = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Osaka = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Hyogo = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Nara = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Wakayama = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Tottori = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Shimane = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Okayama = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Hiroshima = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Yamaguchi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Tokushima = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Kagawa = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Ehime  = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Kochi = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Fukuoka = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Saga = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Nagasaki = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Kumamoto = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Oita = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Miyazaki = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Kagoshima = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_hdfts_Okinawa = array( , dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  
  for(io in 1:(11-fh))
  {
    dat.train<-list()
    for ( i in 1:length(jpn.data)){
      dat.train[[i]]<-extract.years(data_series[[i]], years = 1973:(2012+io))
    }
    dum_pointwise = find_enlarge_val_hdfts(data_series = dat.train, fh = fh, n_pop = length(data_series))
    PI_boot_hdfts_Hokkaido[,,io] = dum_pointwise$boot_sample_female[[1]]
    PI_boot_hdfts_Aomori[,,io] = dum_pointwise$boot_sample_female[[2]]
    PI_boot_hdfts_Iwate[,,io] = dum_pointwise$boot_sample_female[[3]]
    PI_boot_hdfts_Miyagi[,,io] = dum_pointwise$boot_sample_female[[4]]
    PI_boot_hdfts_Akita [,,io] = dum_pointwise$boot_sample_female[[5]]
    PI_boot_hdfts_Yamagata[,,io] = dum_pointwise$boot_sample_female[[6]]
    PI_boot_hdfts_Fukushima[,,io] = dum_pointwise$boot_sample_female[[7]]
    PI_boot_hdfts_Ibaraki[,,io] = dum_pointwise$boot_sample_female[[8]]
    PI_boot_hdfts_Tochigi[,,io] = dum_pointwise$boot_sample_female[[9]]
    PI_boot_hdfts_Gunma[,,io] = dum_pointwise$boot_sample_female[[10]]
    PI_boot_hdfts_Saitama[,,io] = dum_pointwise$boot_sample_female[[11]]
    PI_boot_hdfts_Chiba[,,io] = dum_pointwise$boot_sample_female[[12]]
    PI_boot_hdfts_Tokyo[,,io] = dum_pointwise$boot_sample_female[[13]]
    PI_boot_hdfts_Kanagawa[,,io] = dum_pointwise$boot_sample_female[[14]]
    PI_boot_hdfts_Niigata[,,io] = dum_pointwise$boot_sample_female[[15]]
    PI_boot_hdfts_Toyama[,,io] = dum_pointwise$boot_sample_female[[16]]
    PI_boot_hdfts_Ishikawa[,,io] = dum_pointwise$boot_sample_female[[17]]
    PI_boot_hdfts_Fukui[,,io] = dum_pointwise$boot_sample_female[[18]]
    PI_boot_hdfts_Yamanashi[,,io] = dum_pointwise$boot_sample_female[[19]]
    PI_boot_hdfts_Nagano[,,io] = dum_pointwise$boot_sample_female[[20]]
    PI_boot_hdfts_Gifu[,,io] = dum_pointwise$boot_sample_female[[21]]
    PI_boot_hdfts_Shizuoka[,,io] = dum_pointwise$boot_sample_female[[22]]
    PI_boot_hdfts_Aichi[,,io] = dum_pointwise$boot_sample_female[[23]]
    PI_boot_hdfts_Mie[,,io] = dum_pointwise$boot_sample_female[[24]]
    PI_boot_hdfts_Shiga[,,io] = dum_pointwise$boot_sample_female[[25]]
    PI_boot_hdfts_Kyoto[,,io] = dum_pointwise$boot_sample_female[[26]]
    PI_boot_hdfts_Osaka[,,io] = dum_pointwise$boot_sample_female[[27]]
    PI_boot_hdfts_Hyogo[,,io] = dum_pointwise$boot_sample_female[[28]]
    PI_boot_hdfts_Nara[,,io] = dum_pointwise$boot_sample_female[[29]]
    PI_boot_hdfts_Wakayama[,,io] = dum_pointwise$boot_sample_female[[30]]
    PI_boot_hdfts_Tottori[,,io] = dum_pointwise$boot_sample_female[[31]]
    PI_boot_hdfts_Shimane[,,io] = dum_pointwise$boot_sample_female[[32]]
    PI_boot_hdfts_Okayama[,,io] = dum_pointwise$boot_sample_female[[33]]
    PI_boot_hdfts_Hiroshima[,,io] = dum_pointwise$boot_sample_female[[34]]
    PI_boot_hdfts_Yamaguchi[,,io] = dum_pointwise$boot_sample_female[[35]]
    PI_boot_hdfts_Tokushima[,,io] = dum_pointwise$boot_sample_female[[36]]
    PI_boot_hdfts_Kagawa[,,io] = dum_pointwise$boot_sample_female[[37]]
    PI_boot_hdfts_Ehime [,,io] = dum_pointwise$boot_sample_female[[38]]
    PI_boot_hdfts_Kochi[,,io] = dum_pointwise$boot_sample_female[[39]]
    PI_boot_hdfts_Fukuoka[,,io] = dum_pointwise$boot_sample_female[[40]]
    PI_boot_hdfts_Saga[,,io] = dum_pointwise$boot_sample_female[[41]]
    PI_boot_hdfts_Nagasaki[,,io] = dum_pointwise$boot_sample_female[[42]]
    PI_boot_hdfts_Kumamoto[,,io] = dum_pointwise$boot_sample_female[[43]]
    PI_boot_hdfts_Oita[,,io] = dum_pointwise$boot_sample_female[[44]]
    PI_boot_hdfts_Miyazaki[,,io] = dum_pointwise$boot_sample_female[[45]]
    PI_boot_hdfts_Kagoshima[,,io] = dum_pointwise$boot_sample_female[[46]]
    PI_boot_hdfts_Okinawa[,,io] = dum_pointwise$boot_sample_female[[47]]
  }
  return(list(PI_boot_hdfts_Hokkaido = PI_boot_hdfts_Hokkaido,PI_boot_hdfts_Aomori = PI_boot_hdfts_Aomori,
              PI_boot_hdfts_Iwate = PI_boot_hdfts_Iwate,PI_boot_hdfts_Miyagi = PI_boot_hdfts_Miyagi,
              PI_boot_hdfts_Akita  = PI_boot_hdfts_Akita ,PI_boot_hdfts_Yamagata = PI_boot_hdfts_Yamagata,
              PI_boot_hdfts_Fukushima = PI_boot_hdfts_Fukushima, PI_boot_hdfts_Ibaraki = PI_boot_hdfts_Ibaraki,
              PI_boot_hdfts_Tochigi = PI_boot_hdfts_Tochigi, PI_boot_hdfts_Gunma = PI_boot_hdfts_Gunma,
              PI_boot_hdfts_Saitama = PI_boot_hdfts_Saitama, PI_boot_hdfts_Chiba = PI_boot_hdfts_Chiba, 
              PI_boot_hdfts_Tokyo = PI_boot_hdfts_Tokyo, PI_boot_hdfts_Kanagawa = PI_boot_hdfts_Kanagawa,
              PI_boot_hdfts_Niigata = PI_boot_hdfts_Niigata, PI_boot_hdfts_Toyama = PI_boot_hdfts_Toyama,
              PI_boot_hdfts_Ishikawa = PI_boot_hdfts_Ishikawa,PI_boot_hdfts_Fukui = PI_boot_hdfts_Fukui, 
              PI_boot_hdfts_Yamanashi = PI_boot_hdfts_Yamanashi,PI_boot_hdfts_Nagano = PI_boot_hdfts_Nagano,
              PI_boot_hdfts_Gifu = PI_boot_hdfts_Gifu,PI_boot_hdfts_Shizuoka = PI_boot_hdfts_Shizuoka,
              PI_boot_hdfts_Aichi = PI_boot_hdfts_Aichi,PI_boot_hdfts_Mie = PI_boot_hdfts_Mie,
              PI_boot_hdfts_Shiga = PI_boot_hdfts_Shiga,PI_boot_hdfts_Kyoto = PI_boot_hdfts_Kyoto,
              PI_boot_hdfts_Osaka = PI_boot_hdfts_Osaka, PI_boot_hdfts_Hyogo = PI_boot_hdfts_Hyogo,
              PI_boot_hdfts_Nara = PI_boot_hdfts_Nara, PI_boot_hdfts_Wakayama = PI_boot_hdfts_Wakayama, 
              PI_boot_hdfts_Tottori = PI_boot_hdfts_Tottori,PI_boot_hdfts_Shimane = PI_boot_hdfts_Shimane,
              PI_boot_hdfts_Okayama = PI_boot_hdfts_Okayama,PI_boot_hdfts_Hiroshima = PI_boot_hdfts_Hiroshima,
              PI_boot_hdfts_Yamaguchi = PI_boot_hdfts_Yamaguchi, PI_boot_hdfts_Tokushima = PI_boot_hdfts_Tokushima,
              PI_boot_hdfts_Kagawa = PI_boot_hdfts_Kagawa, PI_boot_hdfts_Ehime  = PI_boot_hdfts_Ehime, 
              PI_boot_hdfts_Kochi = PI_boot_hdfts_Kochi, PI_boot_hdfts_Fukuoka = PI_boot_hdfts_Fukuoka,
              PI_boot_hdfts_Saga = PI_boot_hdfts_Saga,PI_boot_hdfts_Nagasaki = PI_boot_hdfts_Nagasaki,
              PI_boot_hdfts_Kumamoto = PI_boot_hdfts_Kumamoto,PI_boot_hdfts_Oita = PI_boot_hdfts_Oita,
              PI_boot_hdfts_Miyazaki = PI_boot_hdfts_Miyazaki,PI_boot_hdfts_Kagoshima = PI_boot_hdfts_Kagoshima,
              PI_boot_hdfts_Okinawa = PI_boot_hdfts_Okinawa))
}
library(doParallel)
cl <- makeCluster(4) 
registerDoParallel(cl)

hdfts_fh_PI = foreach(mm = 1:10, .packages = "demography") %dopar% PI_fh_hdfts(data_series = dat.female,  fh = mm)

female_score.hdfts.Hokkaido = female_score.hdfts.Aomori = female_score.hdfts.Iwate = female_score.hdfts.Miyagi = 
  female_score.hdfts.Akita = female_score.hdfts.Yamagata = female_score.hdfts.Fukushima = female_score.hdfts.Ibaraki = 
  female_score.hdfts.Tochigi = female_score.hdfts.Gunma = female_score.hdfts.Saitama= female_score.hdfts.Chiba = 
  female_score.hdfts.Tokyo= female_score.hdfts.Kanagawa = female_score.hdfts.Niigata = female_score.hdfts.Toyama = 
  female_score.hdfts.Ishikawa = female_score.hdfts.Fukui = female_score.hdfts.Yamanashi = female_score.hdfts.Nagano = 
  female_score.hdfts.Gifu = female_score.hdfts.Shizuoka = female_score.hdfts.Aichi = female_score.hdfts.Mie = 
  female_score.hdfts.Shiga = female_score.hdfts.Kyoto = female_score.hdfts.Osaka = female_score.hdfts.Hyogo = 
  female_score.hdfts.Nara = female_score.hdfts.Wakayama = female_score.hdfts.Tottori = female_score.hdfts.Shimane = 
  female_score.hdfts.Okayama = female_score.hdfts.Hiroshima = female_score.hdfts.Yamaguchi = female_score.hdfts.Tokushima = 
  female_score.hdfts.Kagawa = female_score.hdfts.Ehime = female_score.hdfts.Kochi = female_score.hdfts.Fukuoka = 
  female_score.hdfts.Saga = female_score.hdfts.Nagasaki = female_score.hdfts.Kumamoto = female_score.hdfts.Oita =   
  female_score.hdfts.Miyazaki = female_score.hdfts.Kagoshima = female_score.hdfts.Okinawa =vector( ,year_horizon)

for(ik in 1:year_horizon)
{
  # me 
  #### Hokkaido###
  female_score.hdfts.Hokkaido[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Hokkaido,
                                                     data_series = dat.female[[1]], fh = ik)$score_female
  
  
  #### Aomori###
  female_score.hdfts.Aomori[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Aomori,
                                                   data_series = dat.female[[2]], fh = ik)$score_female
  #### Iwate###
  female_score.hdfts.Iwate[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Iwate,
                                                  data_series = dat.female[[3]], fh = ik)$score_female
  #### Miyagi###
  female_score.hdfts.Miyagi[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Miyagi,
                                                   data_series = dat.female[[4]], fh = ik)$score_female
  #### Akita###
  female_score.hdfts.Akita[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Akita,
                                                  data_series = dat.female[[5]], fh = ik)$score_female
  #### Yamagata###
  female_score.hdfts.Yamagata[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Yamagata,
                                                     data_series = dat.female[[6]], fh = ik)$score_female
  #### Fukushima###
  female_score.hdfts.Fukushima[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Fukushima,
                                                      data_series = dat.female[[7]], fh = ik)$score_female
  #### Ibaraki###
  female_score.hdfts.Ibaraki[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Ibaraki,
                                                    data_series = dat.female[[8]], fh = ik)$score_female
  #### Tochigi###
  female_score.hdfts.Tochigi[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Tochigi,
                                                    data_series = dat.female[[9]], fh = ik)$score_female
  #### Gunma###
  female_score.hdfts.Gunma[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Gunma,
                                                  data_series = dat.female[[11]], fh = ik)$score_female
  #### Saitama###
  female_score.hdfts.Saitama[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Saitama,
                                                    data_series = dat.female[[12]], fh = ik)$score_female
  #### Chiba###
  female_score.hdfts.Chiba[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Chiba,
                                                  data_series = dat.female[[13]], fh = ik)$score_female
  #### Tokyo###
  female_score.hdfts.Tokyo[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Tokyo,
                                                  data_series = dat.female[[14]], fh = ik)$score_female
  #### Kanagawa###
  female_score.hdfts.Kanagawa[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Kanagawa,
                                                     data_series = dat.female[[1]], fh = ik)$score_female
  #### Niigata###
  female_score.hdfts.Niigata[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Niigata,
                                                    data_series = dat.female[[15]], fh = ik)$score_female
  #### Toyama###
  female_score.hdfts.Toyama[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Toyama,
                                                   data_series = dat.female[[16]], fh = ik)$score_female
  #### Ishikawa###
  female_score.hdfts.Ishikawa[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Ishikawa,
                                                     data_series = dat.female[[17]], fh = ik)$score_female
  #### Fukui###
  female_score.hdfts.Fukui[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Fukui,
                                                  data_series = dat.female[[18]], fh = ik)$score_female
  #### Yamanashi###
  female_score.hdfts.Yamanashi[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Yamanashi,
                                                      data_series = dat.female[[19]], fh = ik)$score_female
  #### Nagano###
  female_score.hdfts.Nagano[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Nagano,
                                                   data_series = dat.female[[20]], fh = ik)$score_female
  #### Gifu###
  female_score.hdfts.Gifu[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Gifu,
                                                 data_series = dat.female[[21]], fh = ik)$score_female
  #### Shizuoka###
  female_score.hdfts.Shizuoka[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Shizuoka,
                                                     data_series = dat.female[[22]], fh = ik)$score_female
  #### Aichi###
  female_score.hdfts.Aichi[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Aichi,
                                                  data_series = dat.female[[23]], fh = ik)$score_female
  #### Mie###
  female_score.hdfts.Mie[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Mie,
                                                data_series = dat.female[[24]], fh = ik)$score_female
  #### Shiga###
  female_score.hdfts.Shiga[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Shiga,
                                                  data_series = dat.female[[25]], fh = ik)$score_female
  #### Kyoto###
  female_score.hdfts.Kyoto[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Kyoto,
                                                  data_series = dat.female[[26]], fh = ik)$score_female
  #### Osaka###
  female_score.hdfts.Osaka[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Osaka,
                                                  data_series = dat.female[[27]], fh = ik)$score_female
  #### Hyogo###
  female_score.hdfts.Hyogo[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Hyogo,
                                                  data_series = dat.female[[28]], fh = ik)$score_female
  #### Nara###
  female_score.hdfts.Nara[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Nara,
                                                 data_series = dat.female[[29]], fh = ik)$score_female
  #### Wakayama###
  female_score.hdfts.Wakayama[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Wakayama,
                                                     data_series = dat.female[[30]], fh = ik)$score_female
  #### Tottori###
  female_score.hdfts.Tottori[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Tottori,
                                                    data_series = dat.female[[31]], fh = ik)$score_female
  #### Shimane###
  female_score.hdfts.Shimane[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Shimane,
                                                    data_series = dat.female[[32]], fh = ik)$score_female
  #### Okayama###
  female_score.hdfts.Okayama[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Okayama,
                                                    data_series = dat.female[[33]], fh = ik)$score_female
  #### Hiroshima###
  female_score.hdfts.Hiroshima[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Hiroshima,
                                                      data_series = dat.female[[34]], fh = ik)$score_female
  #### Yamaguchi###
  female_score.hdfts.Yamaguchi[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Yamaguchi,
                                                      data_series = dat.female[[35]], fh = ik)$score_female
  #### Tokushima###
  female_score.hdfts.Tokushima[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Tokushima,
                                                      data_series = dat.female[[36]], fh = ik)$score_female
  #### Kagawa###
  female_score.hdfts.Kagawa[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Kagawa,
                                                   data_series = dat.female[[37]], fh = ik)$score_female
  #### Ehime###
  female_score.hdfts.Ehime[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Ehime,
                                                  data_series = dat.female[[38]], fh = ik)$score_female
  #### Kochi###
  female_score.hdfts.Kochi[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Kochi,
                                                  data_series = dat.female[[39]], fh = ik)$score_female
  #### Fukuoka###
  female_score.hdfts.Fukuoka[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Fukuoka,
                                                    data_series = dat.female[[40]], fh = ik)$score_female
  #### Saga###
  female_score.hdfts.Saga[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Saga,
                                                 data_series = dat.female[[41]], fh = ik)$score_female
  #### Nagasaki###
  female_score.hdfts.Nagasaki[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Nagasaki,
                                                     data_series = dat.female[[42]], fh = ik)$score_female
  #### Kumamoto###
  female_score.hdfts.Kumamoto[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Kumamoto,
                                                     data_series = dat.female[[43]], fh = ik)$score_female
  #### Oita###
  female_score.hdfts.Oita[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Oita,
                                                 data_series = dat.female[[44]], fh = ik)$score_female
  #### Miyazaki###
  female_score.hdfts.Miyazaki[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Miyazaki,
                                                     data_series = dat.female[[45]], fh = ik)$score_female
  #### Kagoshima###
  female_score.hdfts.Kagoshima[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Kagoshima,
                                                      data_series = dat.female[[46]], fh = ik)$score_female
  #### Okinawa###
  female_score.hdfts.Okinawa[ik] = interval_score_hdfts(PI_val = hdfts_fh_PI[[ik]]$PI_boot_hdfts_Okinawa,
                                                    data_series = dat.female[[47]], fh = ik)$score_female
}


hdfts.jpfemale.score<-cbind(female_score.hdfts.Hokkaido , female_score.hdfts.Aomori , female_score.hdfts.Iwate ,
                          female_score.hdfts.Miyagi , female_score.hdfts.Akita , female_score.hdfts.Yamagata ,
                          female_score.hdfts.Fukushima , female_score.hdfts.Ibaraki , 
                          female_score.hdfts.Tochigi , female_score.hdfts.Gunma , female_score.hdfts.Saitama,
                          female_score.hdfts.Chiba ,  female_score.hdfts.Tokyo, female_score.hdfts.Kanagawa ,
                          female_score.hdfts.Niigata , female_score.hdfts.Toyama , female_score.hdfts.Ishikawa ,
                          female_score.hdfts.Fukui , female_score.hdfts.Yamanashi , female_score.hdfts.Nagano , 
                          female_score.hdfts.Gifu , female_score.hdfts.Shizuoka , female_score.hdfts.Aichi ,
                          female_score.hdfts.Mie ,  female_score.hdfts.Shiga , female_score.hdfts.Kyoto ,
                          female_score.hdfts.Osaka , female_score.hdfts.Hyogo ,  female_score.hdfts.Nara ,
                          female_score.hdfts.Wakayama , female_score.hdfts.Tottori , female_score.hdfts.Shimane ,
                          female_score.hdfts.Okayama , female_score.hdfts.Hiroshima , 
                          female_score.hdfts.Yamaguchi ,  female_score.hdfts.Tokushima , 
                          female_score.hdfts.Kagawa , female_score.hdfts.Ehime , female_score.hdfts.Kochi ,
                          female_score.hdfts.Fukuoka , female_score.hdfts.Saga , female_score.hdfts.Nagasaki ,
                          female_score.hdfts.Kumamoto , female_score.hdfts.Oita , female_score.hdfts.Miyazaki ,
                          female_score.hdfts.Kagoshima , female_score.hdfts.Okinawa )

colnames(hdfts.jpfemale.score)<-c("Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima", "Ibaraki",  "Tochigi", "Gunma", "Saitama",  
                                "Chiba", "Tokyo", "Kanagawa", "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",    
                                "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi",
                                "Tokushima", "Kagawa", "Ehime", "Kochi", "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")
write.csv(hdfts.jpfemale.score, "hdfts.jpfemale.score.csv")













