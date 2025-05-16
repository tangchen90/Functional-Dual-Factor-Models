#################################################
############## bootstrap error function##########
#################################################
boot_male <-function(err, nboot)
{
    boots <-matrix(, nrow(err), nboot)
    for(ij in 1:nboot)
    {
        boots[,ij] = err[, sample(1:ncol(err), 1, replace = TRUE)]
    }
    return(boots)
}

###########################################################
# constructing multivariate pointwise prediction intervals
###########################################################

# data_series: specific data series
# fh: forecast horizon
# nboot: number of bootstrap replication
# alpha: nominal coverage probability

find_enlarge_val_fvp <- function(data_series, fh , nboot = 1000, alpha = 0.8, n_pop)
{
    n_test <- 25

    # calculate in-sample forecast curves

    fore_curve_male  = matrix(NA, length(data_series[[1]]$age)*n_pop, length(data_series[[1]]$year) - n_test - fh )
    for(ij in 1:(length(data_series[[1]]$year) - n_test - fh))
    {
        dat.rate.male_h_step <- list()
        for(ik in 1:n_pop)
        {
          dat.rate.male_h_step[[ik]]<-t(log(extract.years(data_series[[ik]], years = 1973:(1973 + n_test + ij - 1))$rate$male))
        }
        fore_res = fvpts(dat=dat.rate.male_h_step, year_horizon = fh, n_pop = length(data_series), ages = 0:100)
        # fill in gender specific in-sample forecast curves
        fore_curve_male[,ij] = log(fore_res[,fh])
    }

    # holdout data samples
    true_dat = list()
    for(im in 1: n_pop)
    {
        true_dat[[im]] = as.data.frame(log(data_series[[im]]$rate$male[, (n_test + fh + 1):length(data_series[[im]]$year)]))
    }

    # male sample
    holdout_val_male = do.call(rbind, true_dat)

    err_male = holdout_val_male - fore_curve_male

    male.err <- list()
    for(il in 1: n_pop)
    {
        male.err[[il]]<-err_male[(1:101)+(il-1)*101 ,]
    }
    male.boot_err <- list()
    for(il in 1: n_pop)
    {
       male.boot_err[[il]]<-boot_male(male.err[[il]], 1000)
    }

    # constructing PI
    dat.series<-list()
    for(ik in 1:n_pop)
    {
        dat.series[[ik]] <- t(log(data_series[[ik]]$rate$male))
    }

    fore_res = fvpts(dat = dat.series, year_horizon = fh, n_pop = length(dat.series), ages = 0:100)

    fore_fvp_male <- list()
    for(il in 1: n_pop)
    {
        fore_fvp_male[[il]]<-matrix(rep(log(fore_res[(1:101)+(il-1)*101 ,fh]), nboot), nrow = 101, ncol = nboot)
    }

    boot_PI_male <- list()
    for(il in 1: n_pop)
    {
        boot_PI_male[[il]] <- male.boot_err[[il]] + fore_fvp_male[[il]]
    }

    boot_PI_lb_ub_male <- list()
    for(il in 1: n_pop)
    {
        boot_PI_lb_ub_male[[il]]<- apply(boot_PI_male[[il]], 1, quantile, c((1 - alpha)/2, (1 + alpha)/2))
    }
    return(list(boot_PI_lb_ub_male = boot_PI_lb_ub_male, boot_sample_male = boot_PI_male))
}

###############################################################################
# compute interval scores  for a particular forecast horizon
###############################################################################

# PI_val: one- to 10-step-ahead prediction intervals
# alpha: nominal coverage probability

interval_score_fvp <- function(PI_val, data_series, fh, alpha = 0.8)
{
  test_val_male<- extract.years(data_series, (2012+fh):2022)$rate$male

  # transform back to the original scale

  boot_sample_male = exp(PI_val)
  boot_index_male = which(boot_sample_male > 1)
  boot_index_below_male = which(boot_sample_male < 0)
  if(length(boot_index_male) > 0)
  {
    boot_sample_v1_male = replace(boot_sample_male, boot_index_male, 1)
  }
  else
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
for(i in 1:length(jpn.data))
{
    dat.male[[i]]<-smooth.demogdata(extract.ages(extract.years(jpn.data[[i]], years = 1973:2022), 0:100, combine.upper=TRUE))
}

####### Prediction Intervals Constructions #######
PI_fh_fvp<- function(data_series, fh, nboot = 1000)
{
  PI_boot_fvp_Hokkaido = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Aomori = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Iwate = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Miyagi = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Akita  = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Yamagata = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Fukushima = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Ibaraki = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Tochigi = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Gunma = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Saitama = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Chiba = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Tokyo = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Kanagawa = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Niigata = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Toyama = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Ishikawa = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Fukui = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Yamanashi = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Nagano = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Gifu = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Shizuoka = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Aichi = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Mie = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Shiga = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Kyoto = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Osaka = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Hyogo = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Nara = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Wakayama = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Tottori = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Shimane = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Okayama = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Hiroshima = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Yamaguchi = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Tokushima = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Kagawa = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Ehime  = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Kochi = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Fukuoka = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Saga = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Nagasaki = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Kumamoto = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Oita = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Miyazaki = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Kagoshima = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  PI_boot_fvp_Okinawa = array(NA, dim = c(length(data_series[[1]]$age), nboot, (11-fh)))
  
  for(io in 1:(11-fh))
  {
    dat.train<-list()
    for ( i in 1:length(jpn.data)){
      dat.train[[i]]<-extract.years(data_series[[i]], years = 1973:(2012+io))
    }
    dum_pointwise = find_enlarge_val_fvp(data_series = dat.train, fh = fh, n_pop = length(data_series))
    PI_boot_fvp_Hokkaido[,,io] = dum_pointwise$boot_sample_male[[1]]
    PI_boot_fvp_Aomori[,,io] = dum_pointwise$boot_sample_male[[2]]
    PI_boot_fvp_Iwate[,,io] = dum_pointwise$boot_sample_male[[3]]
    PI_boot_fvp_Miyagi[,,io] = dum_pointwise$boot_sample_male[[4]]
    PI_boot_fvp_Akita [,,io] = dum_pointwise$boot_sample_male[[5]]
    PI_boot_fvp_Yamagata[,,io] = dum_pointwise$boot_sample_male[[6]]
    PI_boot_fvp_Fukushima[,,io] = dum_pointwise$boot_sample_male[[7]]
    PI_boot_fvp_Ibaraki[,,io] = dum_pointwise$boot_sample_male[[8]]
    PI_boot_fvp_Tochigi[,,io] = dum_pointwise$boot_sample_male[[9]]
    PI_boot_fvp_Gunma[,,io] = dum_pointwise$boot_sample_male[[10]]
    PI_boot_fvp_Saitama[,,io] = dum_pointwise$boot_sample_male[[11]]
    PI_boot_fvp_Chiba[,,io] = dum_pointwise$boot_sample_male[[12]]
    PI_boot_fvp_Tokyo[,,io] = dum_pointwise$boot_sample_male[[13]]
    PI_boot_fvp_Kanagawa[,,io] = dum_pointwise$boot_sample_male[[14]]
    PI_boot_fvp_Niigata[,,io] = dum_pointwise$boot_sample_male[[15]]
    PI_boot_fvp_Toyama[,,io] = dum_pointwise$boot_sample_male[[16]]
    PI_boot_fvp_Ishikawa[,,io] = dum_pointwise$boot_sample_male[[17]]
    PI_boot_fvp_Fukui[,,io] = dum_pointwise$boot_sample_male[[18]]
    PI_boot_fvp_Yamanashi[,,io] = dum_pointwise$boot_sample_male[[19]]
    PI_boot_fvp_Nagano[,,io] = dum_pointwise$boot_sample_male[[20]]
    PI_boot_fvp_Gifu[,,io] = dum_pointwise$boot_sample_male[[21]]
    PI_boot_fvp_Shizuoka[,,io] = dum_pointwise$boot_sample_male[[22]]
    PI_boot_fvp_Aichi[,,io] = dum_pointwise$boot_sample_male[[23]]
    PI_boot_fvp_Mie[,,io] = dum_pointwise$boot_sample_male[[24]]
    PI_boot_fvp_Shiga[,,io] = dum_pointwise$boot_sample_male[[25]]
    PI_boot_fvp_Kyoto[,,io] = dum_pointwise$boot_sample_male[[26]]
    PI_boot_fvp_Osaka[,,io] = dum_pointwise$boot_sample_male[[27]]
    PI_boot_fvp_Hyogo[,,io] = dum_pointwise$boot_sample_male[[28]]
    PI_boot_fvp_Nara[,,io] = dum_pointwise$boot_sample_male[[29]]
    PI_boot_fvp_Wakayama[,,io] = dum_pointwise$boot_sample_male[[30]]
    PI_boot_fvp_Tottori[,,io] = dum_pointwise$boot_sample_male[[31]]
    PI_boot_fvp_Shimane[,,io] = dum_pointwise$boot_sample_male[[32]]
    PI_boot_fvp_Okayama[,,io] = dum_pointwise$boot_sample_male[[33]]
    PI_boot_fvp_Hiroshima[,,io] = dum_pointwise$boot_sample_male[[34]]
    PI_boot_fvp_Yamaguchi[,,io] = dum_pointwise$boot_sample_male[[35]]
    PI_boot_fvp_Tokushima[,,io] = dum_pointwise$boot_sample_male[[36]]
    PI_boot_fvp_Kagawa[,,io] = dum_pointwise$boot_sample_male[[37]]
    PI_boot_fvp_Ehime [,,io] = dum_pointwise$boot_sample_male[[38]]
    PI_boot_fvp_Kochi[,,io] = dum_pointwise$boot_sample_male[[39]]
    PI_boot_fvp_Fukuoka[,,io] = dum_pointwise$boot_sample_male[[40]]
    PI_boot_fvp_Saga[,,io] = dum_pointwise$boot_sample_male[[41]]
    PI_boot_fvp_Nagasaki[,,io] = dum_pointwise$boot_sample_male[[42]]
    PI_boot_fvp_Kumamoto[,,io] = dum_pointwise$boot_sample_male[[43]]
    PI_boot_fvp_Oita[,,io] = dum_pointwise$boot_sample_male[[44]]
    PI_boot_fvp_Miyazaki[,,io] = dum_pointwise$boot_sample_male[[45]]
    PI_boot_fvp_Kagoshima[,,io] = dum_pointwise$boot_sample_male[[46]]
    PI_boot_fvp_Okinawa[,,io] = dum_pointwise$boot_sample_male[[47]]
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

FVP_fh_PI = foreach(mm = 1:10, .packages = c("demography","fda")) %dopar% PI_fh_fvp(data_series = dat.male,  fh = mm)

save.image("14-3-2024_2.RData")

interval_score_fvp_male <- function(PI_val, data_series, fh, alpha = 0.8)
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



male_score.fvp.Hokkaido = male_score.fvp.Aomori = male_score.fvp.Iwate = male_score.fvp.Miyagi =
  male_score.fvp.Akita = male_score.fvp.Yamagata = male_score.fvp.Fukushima = male_score.fvp.Ibaraki =
  male_score.fvp.Tochigi = male_score.fvp.Gunma = male_score.fvp.Saitama= male_score.fvp.Chiba =
  male_score.fvp.Tokyo= male_score.fvp.Kanagawa = male_score.fvp.Niigata = male_score.fvp.Toyama =
  male_score.fvp.Ishikawa = male_score.fvp.Fukui = male_score.fvp.Yamanashi = male_score.fvp.Nagano =
  male_score.fvp.Gifu = male_score.fvp.Shizuoka = male_score.fvp.Aichi = male_score.fvp.Mie =
  male_score.fvp.Shiga = male_score.fvp.Kyoto = male_score.fvp.Osaka = male_score.fvp.Hyogo =
  male_score.fvp.Nara = male_score.fvp.Wakayama = male_score.fvp.Tottori = male_score.fvp.Shimane =
  male_score.fvp.Okayama = male_score.fvp.Hiroshima = male_score.fvp.Yamaguchi = male_score.fvp.Tokushima =
  male_score.fvp.Kagawa = male_score.fvp.Ehime = male_score.fvp.Kochi = male_score.fvp.Fukuoka =
  male_score.fvp.Saga = male_score.fvp.Nagasaki = male_score.fvp.Kumamoto = male_score.fvp.Oita =
  male_score.fvp.Miyazaki = male_score.fvp.Kagoshima = male_score.fvp.Okinawa =vector( ,year_horizon)

for(ik in 1:year_horizon)
{
  # me
  #### Hokkaido###
  male_score.fvp.Hokkaido[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Hokkaido,
                                                     data_series = dat.male[[1]], fh = ik)$score_male


  #### Aomori###
  male_score.fvp.Aomori[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Aomori,
                                                   data_series = dat.male[[2]], fh = ik)$score_male
  #### Iwate###
  male_score.fvp.Iwate[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Iwate,
                                                  data_series = dat.male[[3]], fh = ik)$score_male
  #### Miyagi###
  male_score.fvp.Miyagi[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Miyagi,
                                                   data_series = dat.male[[4]], fh = ik)$score_male
  #### Akita###
  male_score.fvp.Akita[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Akita,
                                                  data_series = dat.male[[5]], fh = ik)$score_male
  #### Yamagata###
  male_score.fvp.Yamagata[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Yamagata,
                                                     data_series = dat.male[[6]], fh = ik)$score_male
  #### Fukushima###
  male_score.fvp.Fukushima[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Fukushima,
                                                      data_series = dat.male[[7]], fh = ik)$score_male
  #### Ibaraki###
  male_score.fvp.Ibaraki[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Ibaraki,
                                                    data_series = dat.male[[8]], fh = ik)$score_male
  #### Tochigi###
  male_score.fvp.Tochigi[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Tochigi,
                                                    data_series = dat.male[[9]], fh = ik)$score_male
  #### Gunma###
  male_score.fvp.Gunma[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Gunma,
                                                  data_series = dat.male[[11]], fh = ik)$score_male
  #### Saitama###
  male_score.fvp.Saitama[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Saitama,
                                                    data_series = dat.male[[12]], fh = ik)$score_male
  #### Chiba###
  male_score.fvp.Chiba[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Chiba,
                                                  data_series = dat.male[[13]], fh = ik)$score_male
  #### Tokyo###
  male_score.fvp.Tokyo[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Tokyo,
                                                  data_series = dat.male[[14]], fh = ik)$score_male
  #### Kanagawa###
  male_score.fvp.Kanagawa[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Kanagawa,
                                                     data_series = dat.male[[1]], fh = ik)$score_male
  #### Niigata###
  male_score.fvp.Niigata[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Niigata,
                                                    data_series = dat.male[[15]], fh = ik)$score_male
  #### Toyama###
  male_score.fvp.Toyama[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Toyama,
                                                   data_series = dat.male[[16]], fh = ik)$score_male
  #### Ishikawa###
  male_score.fvp.Ishikawa[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Ishikawa,
                                                     data_series = dat.male[[17]], fh = ik)$score_male
  #### Fukui###
  male_score.fvp.Fukui[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Fukui,
                                                  data_series = dat.male[[18]], fh = ik)$score_male
  #### Yamanashi###
  male_score.fvp.Yamanashi[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Yamanashi,
                                                      data_series = dat.male[[19]], fh = ik)$score_male
  #### Nagano###
  male_score.fvp.Nagano[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Nagano,
                                                   data_series = dat.male[[20]], fh = ik)$score_male
  #### Gifu###
  male_score.fvp.Gifu[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Gifu,
                                                 data_series = dat.male[[21]], fh = ik)$score_male
  #### Shizuoka###
  male_score.fvp.Shizuoka[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Shizuoka,
                                                     data_series = dat.male[[22]], fh = ik)$score_male
  #### Aichi###
  male_score.fvp.Aichi[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Aichi,
                                                  data_series = dat.male[[23]], fh = ik)$score_male
  #### Mie###
  male_score.fvp.Mie[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Mie,
                                                data_series = dat.male[[24]], fh = ik)$score_male
  #### Shiga###
  male_score.fvp.Shiga[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Shiga,
                                                  data_series = dat.male[[25]], fh = ik)$score_male
  #### Kyoto###
  male_score.fvp.Kyoto[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Kyoto,
                                                  data_series = dat.male[[26]], fh = ik)$score_male
  #### Osaka###
  male_score.fvp.Osaka[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Osaka,
                                                  data_series = dat.male[[27]], fh = ik)$score_male
  #### Hyogo###
  male_score.fvp.Hyogo[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Hyogo,
                                                  data_series = dat.male[[28]], fh = ik)$score_male
  #### Nara###
  male_score.fvp.Nara[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Nara,
                                                 data_series = dat.male[[29]], fh = ik)$score_male
  #### Wakayama###
  male_score.fvp.Wakayama[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Wakayama,
                                                     data_series = dat.male[[30]], fh = ik)$score_male
  #### Tottori###
  male_score.fvp.Tottori[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Tottori,
                                                    data_series = dat.male[[31]], fh = ik)$score_male
  #### Shimane###
  male_score.fvp.Shimane[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Shimane,
                                                    data_series = dat.male[[32]], fh = ik)$score_male
  #### Okayama###
  male_score.fvp.Okayama[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Okayama,
                                                    data_series = dat.male[[33]], fh = ik)$score_male
  #### Hiroshima###
  male_score.fvp.Hiroshima[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Hiroshima,
                                                      data_series = dat.male[[34]], fh = ik)$score_male
  #### Yamaguchi###
  male_score.fvp.Yamaguchi[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Yamaguchi,
                                                      data_series = dat.male[[35]], fh = ik)$score_male
  #### Tokushima###
  male_score.fvp.Tokushima[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Tokushima,
                                                      data_series = dat.male[[36]], fh = ik)$score_male
  #### Kagawa###
  male_score.fvp.Kagawa[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Kagawa,
                                                   data_series = dat.male[[37]], fh = ik)$score_male
  #### Ehime###
  male_score.fvp.Ehime[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Ehime,
                                                  data_series = dat.male[[38]], fh = ik)$score_male
  #### Kochi###
  male_score.fvp.Kochi[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Kochi,
                                                  data_series = dat.male[[39]], fh = ik)$score_male
  #### Fukuoka###
  male_score.fvp.Fukuoka[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Fukuoka,
                                                    data_series = dat.male[[40]], fh = ik)$score_male
  #### Saga###
  male_score.fvp.Saga[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Saga,
                                                 data_series = dat.male[[41]], fh = ik)$score_male
  #### Nagasaki###
  male_score.fvp.Nagasaki[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Nagasaki,
                                                     data_series = dat.male[[42]], fh = ik)$score_male
  #### Kumamoto###
  male_score.fvp.Kumamoto[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Kumamoto,
                                                     data_series = dat.male[[43]], fh = ik)$score_male
  #### Oita###
  male_score.fvp.Oita[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Oita,
                                                 data_series = dat.male[[44]], fh = ik)$score_male
  #### Miyazaki###
  male_score.fvp.Miyazaki[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Miyazaki,
                                                     data_series = dat.male[[45]], fh = ik)$score_male
  #### Kagoshima###
  male_score.fvp.Kagoshima[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Kagoshima,
                                                      data_series = dat.male[[46]], fh = ik)$score_male
  #### Okinawa###
  male_score.fvp.Okinawa[ik] = interval_score_fvp_male(PI_val = FVP_fh_PI[[ik]]$PI_boot_fvp_Okinawa,
                                                    data_series = dat.male[[47]], fh = ik)$score_male
}


fvp.jpmale.score<-cbind(male_score.fvp.Hokkaido , male_score.fvp.Aomori , male_score.fvp.Iwate ,
                          male_score.fvp.Miyagi , male_score.fvp.Akita , male_score.fvp.Yamagata ,
                          male_score.fvp.Fukushima , male_score.fvp.Ibaraki ,
                          male_score.fvp.Tochigi , male_score.fvp.Gunma , male_score.fvp.Saitama,
                          male_score.fvp.Chiba ,  male_score.fvp.Tokyo, male_score.fvp.Kanagawa ,
                          male_score.fvp.Niigata , male_score.fvp.Toyama , male_score.fvp.Ishikawa ,
                          male_score.fvp.Fukui , male_score.fvp.Yamanashi , male_score.fvp.Nagano ,
                          male_score.fvp.Gifu , male_score.fvp.Shizuoka , male_score.fvp.Aichi ,
                          male_score.fvp.Mie ,  male_score.fvp.Shiga , male_score.fvp.Kyoto ,
                          male_score.fvp.Osaka , male_score.fvp.Hyogo ,  male_score.fvp.Nara ,
                          male_score.fvp.Wakayama , male_score.fvp.Tottori , male_score.fvp.Shimane ,
                          male_score.fvp.Okayama , male_score.fvp.Hiroshima ,
                          male_score.fvp.Yamaguchi ,  male_score.fvp.Tokushima ,
                          male_score.fvp.Kagawa , male_score.fvp.Ehime , male_score.fvp.Kochi ,
                          male_score.fvp.Fukuoka , male_score.fvp.Saga , male_score.fvp.Nagasaki ,
                          male_score.fvp.Kumamoto , male_score.fvp.Oita , male_score.fvp.Miyazaki ,
                          male_score.fvp.Kagoshima , male_score.fvp.Okinawa )

colnames(fvp.jpmale.score)<-c("Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima", "Ibaraki",  "Tochigi", "Gunma", "Saitama",
                             "Chiba", "Tokyo", "Kanagawa", "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",
                             "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi",
                             "Tokushima", "Kagawa", "Ehime", "Kochi", "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")
fvp.jpmale.score = rbind(fvp.jpmale.score, apply(fvp.jpmale.score, 2, mean))
write.csv(fvp.jpmale.score, "fvp.jpmale.score.csv")
