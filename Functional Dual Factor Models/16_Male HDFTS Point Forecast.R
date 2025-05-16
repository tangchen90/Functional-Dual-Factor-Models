######Forecast 
year_horizon=10
forecast_year = (year_horizon + 1)
n_year = 2022 - forecast_year
res_male.Hokkaido = res_male.Aomori = res_male.Iwate = res_male.Miyagi = 
  res_male.Akita = res_male.Yamagata = res_male.Fukushima = 
  res_male.Ibaraki = res_male.Tochigi = res_male.Gunma = 
  res_male.Saitama = res_male.Chiba = res_male.Tokyo= res_male.Kanagawa = res_male.Niigata = res_male.Toyama = res_male.Ishikawa = res_male.Fukui = res_male.Yamanashi = res_male.Nagano = res_male.Gifu = 
  res_male.Shizuoka = res_male.Aichi = res_male.Mie = res_male.Shiga = res_male.Kyoto = res_male.Osaka = res_male.Hyogo = res_male.Nara = res_male.Wakayama = res_male.Tottori = res_male.Shimane = 
  res_male.Okayama = res_male.Hiroshima = res_male.Yamaguchi = 
  res_male.Tokushima = res_male.Kagawa = res_male.Ehime = 
  res_male.Kochi = res_male.Fukuoka = res_male.Saga = res_male.Nagasaki = res_male.Kumamoto = res_male.Oita = res_male.Miyazaki = 
  res_male.Kagoshima = res_male.Okinawa =array(NA, dim = c(year_horizon, 101, year_horizon))

for(ij in 1:year_horizon)
{
  jpn.trunc<-list()
  for ( i in 1:47){
    jpn.trunc[[i]]<-smooth.demogdata(extract.ages(extract.years(jpn.data[[i]], years = 1973:(n_year+ij)), 0:100, combine.upper=TRUE))
  }
  jpn.truc<-list()
  for (ik in 1:47){
    jpn.truc[[ik]]<-t(log(jpn.trunc[[ik]]$rate$male))
  }
  res_forc<-hdfts.forc(jpn.truc, year_horizon = year_horizon, years = 1973:(n_year+ij), ages = 0:100, n_pop =47)
  res_male.Hokkaido[,,ij] = res_forc[[1]]
  res_male.Aomori[,,ij] = res_forc[[2]]
  res_male.Iwate[,,ij] = res_forc[[3]]
  res_male.Miyagi[,,ij] = res_forc[[4]]
  res_male.Akita [,,ij] = res_forc[[5]]
  res_male.Yamagata[,,ij] = res_forc[[6]]
  res_male.Fukushima[,,ij] = res_forc[[7]]
  res_male.Ibaraki[,,ij] = res_forc[[8]]
  res_male.Tochigi[,,ij] = res_forc[[9]]
  res_male.Gunma[,,ij] = res_forc[[10]]
  res_male.Saitama[,,ij] = res_forc[[11]]
  res_male.Chiba[,,ij] = res_forc[[12]]
  res_male.Tokyo[,,ij] = res_forc[[13]]
  res_male.Kanagawa[,,ij] = res_forc[[14]]
  res_male.Niigata[,,ij] = res_forc[[15]]
  res_male.Toyama[,,ij] = res_forc[[16]]
  res_male.Ishikawa[,,ij] = res_forc[[17]]
  res_male.Fukui[,,ij] = res_forc[[18]]
  res_male.Yamanashi[,,ij] = res_forc[[19]]
  res_male.Nagano[,,ij] = res_forc[[20]]
  res_male.Gifu[,,ij] = res_forc[[21]]
  res_male.Shizuoka[,,ij] = res_forc[[22]]
  res_male.Aichi[,,ij] = res_forc[[23]]
  res_male.Mie[,,ij] = res_forc[[24]]
  res_male.Shiga[,,ij] = res_forc[[25]]
  res_male.Kyoto[,,ij] = res_forc[[26]]
  res_male.Osaka[,,ij] = res_forc[[27]]
  res_male.Hyogo[,,ij] = res_forc[[28]]
  res_male.Nara[,,ij] = res_forc[[29]]
  res_male.Wakayama[,,ij] = res_forc[[30]]
  res_male.Tottori[,,ij] = res_forc[[31]]
  res_male.Shimane[,,ij] = res_forc[[32]]
  res_male.Okayama[,,ij] = res_forc[[33]]
  res_male.Hiroshima[,,ij] = res_forc[[34]]
  res_male.Yamaguchi[,,ij] = res_forc[[35]]
  res_male.Tokushima[,,ij] = res_forc[[36]]
  res_male.Kagawa[,,ij] = res_forc[[37]]
  res_male.Ehime [,,ij] = res_forc[[38]]
  res_male.Kochi[,,ij] = res_forc[[39]]
  res_male.Fukuoka[,,ij] = res_forc[[40]]
  res_male.Saga[,,ij] = res_forc[[41]]
  res_male.Nagasaki[,,ij] = res_forc[[42]]
  res_male.Kumamoto[,,ij] = res_forc[[43]]
  res_male.Oita[,,ij] = res_forc[[44]]
  res_male.Miyazaki[,,ij] = res_forc[[45]]
  res_male.Kagoshima[,,ij] = res_forc[[46]]
  res_male.Okinawa[,,ij] = res_forc[[47]]
}


me   = ftsa:::me
mae  = ftsa:::mae
rmse = ftsa:::rmse
# MAE & RMSE
male_me.Hokkaido = male_me.Aomori = male_me.Iwate = male_me.Miyagi = 
  male_me.Akita = male_me.Yamagata = male_me.Fukushima = male_me.Ibaraki = 
  male_me.Tochigi = male_me.Gunma = male_me.Saitama= male_me.Chiba = 
  male_me.Tokyo= male_me.Kanagawa = male_me.Niigata = male_me.Toyama = 
  male_me.Ishikawa = male_me.Fukui = male_me.Yamanashi = male_me.Nagano = 
  male_me.Gifu = male_me.Shizuoka = male_me.Aichi = male_me.Mie = 
  male_me.Shiga = male_me.Kyoto = male_me.Osaka = male_me.Hyogo = 
  male_me.Nara = male_me.Wakayama = male_me.Tottori = male_me.Shimane = 
  male_me.Okayama = male_me.Hiroshima = male_me.Yamaguchi = male_me.Tokushima = 
  male_me.Kagawa = male_me.Ehime = male_me.Kochi = male_me.Fukuoka = 
  male_me.Saga = male_me.Nagasaki = male_me.Kumamoto = male_me.Oita =   
  male_me.Miyazaki = male_me.Kagoshima = male_me.Okinawa =
  male_mae.Hokkaido = male_mae.Aomori = male_mae.Iwate = male_mae.Miyagi = 
  male_mae.Akita = male_mae.Yamagata = male_mae.Fukushima = male_mae.Ibaraki = 
  male_mae.Tochigi = male_mae.Gunma = male_mae.Saitama= male_mae.Chiba = 
  male_mae.Tokyo= male_mae.Kanagawa = male_mae.Niigata = male_mae.Toyama = 
  male_mae.Ishikawa = male_mae.Fukui = male_mae.Yamanashi = male_mae.Nagano = 
  male_mae.Gifu = male_mae.Shizuoka = male_mae.Aichi = male_mae.Mie = 
  male_mae.Shiga = male_mae.Kyoto = male_mae.Osaka = male_mae.Hyogo = 
  male_mae.Nara = male_mae.Wakayama = male_mae.Tottori = male_mae.Shimane = 
  male_mae.Okayama = male_mae.Hiroshima = male_mae.Yamaguchi = male_mae.Tokushima = 
  male_mae.Kagawa = male_mae.Ehime = male_mae.Kochi = male_mae.Fukuoka = 
  male_mae.Saga = male_mae.Nagasaki = male_mae.Kumamoto = male_mae.Oita =   
  male_mae.Miyazaki = male_mae.Kagoshima = male_mae.Okinawa = 
  male_rmse.Hokkaido = male_rmse.Aomori = male_rmse.Iwate = male_rmse.Miyagi = 
  male_rmse.Akita = male_rmse.Yamagata = male_rmse.Fukushima = male_rmse.Ibaraki = 
  male_rmse.Tochigi = male_rmse.Gunma = male_rmse.Saitama= male_rmse.Chiba = 
  male_rmse.Tokyo= male_rmse.Kanagawa = male_rmse.Niigata = male_rmse.Toyama = 
  male_rmse.Ishikawa = male_rmse.Fukui = male_rmse.Yamanashi = male_rmse.Nagano = 
  male_rmse.Gifu = male_rmse.Shizuoka = male_rmse.Aichi = male_rmse.Mie = 
  male_rmse.Shiga = male_rmse.Kyoto = male_rmse.Osaka = male_rmse.Hyogo = 
  male_rmse.Nara = male_rmse.Wakayama = male_rmse.Tottori = male_rmse.Shimane = 
  male_rmse.Okayama = male_rmse.Hiroshima = male_rmse.Yamaguchi = 
  male_rmse.Tokushima = male_rmse.Kagawa = male_rmse.Ehime = male_rmse.Kochi = male_rmse.Fukuoka = male_rmse.Saga = male_rmse.Nagasaki = male_rmse.Kumamoto = male_rmse.Oita = male_rmse.Miyazaki = male_rmse.Kagoshima = male_rmse.Okinawa = 
  vector("numeric",year_horizon)

for(ik in 1:year_horizon)
{
  # me 
  #### Hokkaido###
  male_me.Hokkaido[ik] = me(res_male.Hokkaido[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[1]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  
  #### Aomori###
  male_me.Aomori[ik] = me(res_male.Aomori[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[2]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Iwate###
  male_me.Iwate[ik] = me(res_male.Iwate[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[3]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Miyagi###
  male_me.Miyagi[ik] = me(res_male.Miyagi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[4]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Akita###
  male_me.Akita[ik] = me(res_male.Akita[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[5]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Yamagata###
  male_me.Yamagata[ik] = me(res_male.Yamagata[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[6]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Fukushima###
  male_me.Fukushima[ik] = me(res_male.Fukushima[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[7]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Ibaraki###
  male_me.Ibaraki[ik] = me(res_male.Ibaraki[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[8]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Tochigi###
  male_me.Tochigi[ik] = me(res_male.Tochigi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[9]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Gunma###
  male_me.Gunma[ik] = me(res_male.Gunma[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[10]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Saitama###
  male_me.Saitama[ik] = me(res_male.Saitama[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[11]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Chiba###
  male_me.Chiba[ik] = me(res_male.Chiba[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[12]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Tokyo###
  male_me.Tokyo[ik] = me(res_male.Tokyo[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[13]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Kanagawa###
  male_me.Kanagawa[ik] = me(res_male.Kanagawa[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[14]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Niigata###
  male_me.Niigata[ik] = me(res_male.Niigata[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[15]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Toyama###
  male_me.Toyama[ik] = me(res_male.Toyama[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[16]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Ishikawa###
  male_me.Ishikawa[ik] = me(res_male.Ishikawa[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[17]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Fukui###
  male_me.Fukui[ik] = me(res_male.Fukui[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[18]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Yamanashi###
  male_me.Yamanashi[ik] = me(res_male.Yamanashi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[19]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Nagano###
  male_me.Nagano[ik] = me(res_male.Nagano[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[20]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Gifu###
  male_me.Gifu[ik] = me(res_male.Gifu[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[21]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Shizuoka###
  male_me.Shizuoka[ik] = me(res_male.Shizuoka[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[22]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Aichi###
  male_me.Aichi[ik] = me(res_male.Aichi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[23]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Mie###
  male_me.Mie[ik] = me(res_male.Mie[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[24]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Shiga###
  male_me.Shiga[ik] = me(res_male.Shiga[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[25]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Kyoto###
  male_me.Kyoto[ik] = me(res_male.Kyoto[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[26]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Osaka###
  male_me.Osaka[ik] = me(res_male.Osaka[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[27]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Hyogo###
  male_me.Hyogo[ik] = me(res_male.Hyogo[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[28]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Nara###
  male_me.Nara[ik] = me(res_male.Nara[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[29]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Wakayama###
  male_me.Wakayama[ik] = me(res_male.Wakayama[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[30]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Tottori###
  male_me.Tottori[ik] = me(res_male.Tottori[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[31]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Shimane###
  male_me.Shimane[ik] = me(res_male.Shimane[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[32]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Okayama###
  male_me.Okayama[ik] = me(res_male.Okayama[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[33]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Hiroshima###
  male_me.Hiroshima[ik] = me(res_male.Hiroshima[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[34]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Yamaguchi###
  male_me.Yamaguchi[ik] = me(res_male.Yamaguchi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[35]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Tokushima###
  male_me.Tokushima[ik] = me(res_male.Tokushima[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[36]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Kagawa###
  male_me.Kagawa[ik] = me(res_male.Kagawa[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[37]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Ehime###
  male_me.Ehime[ik] = me(res_male.Ehime[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[38]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Kochi###
  male_me.Kochi[ik] = me(res_male.Kochi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[39]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Fukuoka###
  male_me.Fukuoka[ik] = me(res_male.Fukuoka[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[40]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Saga###
  male_me.Saga[ik] = me(res_male.Saga[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[41]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Nagasaki###
  male_me.Nagasaki[ik] = me(res_male.Nagasaki[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[42]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Kumamoto###
  male_me.Kumamoto[ik] = me(res_male.Kumamoto[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[43]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Oita###
  male_me.Oita[ik] = me(res_male.Oita[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[44]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Miyazaki###
  male_me.Miyazaki[ik] = me(res_male.Miyazaki[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[45]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Kagoshima###
  male_me.Kagoshima[ik] = me(res_male.Kagoshima[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[46]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Okinawa###
  male_me.Okinawa[ik] = me(res_male.Okinawa[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[47]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  # mae 
  #### Hokkaido###
  male_mae.Hokkaido[ik] = mae(res_male.Hokkaido[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[1]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Aomori###
  male_mae.Aomori[ik] = mae(res_male.Aomori[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[2]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Iwate###
  male_mae.Iwate[ik] = mae(res_male.Iwate[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[3]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Miyagi###
  male_mae.Miyagi[ik] = mae(res_male.Miyagi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[4]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Akita###
  male_mae.Akita[ik] = mae(res_male.Akita[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[5]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Yamagata###
  male_mae.Yamagata[ik] = mae(res_male.Yamagata[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[6]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Fukushima###
  male_mae.Fukushima[ik] = mae(res_male.Fukushima[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[7]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Ibaraki###
  male_mae.Ibaraki[ik] = mae(res_male.Ibaraki[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[8]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Tochigi###
  male_mae.Tochigi[ik] = mae(res_male.Tochigi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[9]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Gunma###
  male_mae.Gunma[ik] = mae(res_male.Gunma[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[10]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Saitama###
  male_mae.Saitama[ik] = mae(res_male.Saitama[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[11]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Chiba###
  male_mae.Chiba[ik] = mae(res_male.Chiba[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[12]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Tokyo###
  male_mae.Tokyo[ik] = mae(res_male.Tokyo[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[13]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Kanagawa###
  male_mae.Kanagawa[ik] = mae(res_male.Kanagawa[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[14]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Niigata###
  male_mae.Niigata[ik] = mae(res_male.Niigata[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[15]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Toyama###
  male_mae.Toyama[ik] = mae(res_male.Toyama[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[16]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Ishikawa###
  male_mae.Ishikawa[ik] = mae(res_male.Ishikawa[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[17]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Fukui###
  male_mae.Fukui[ik] = mae(res_male.Fukui[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[18]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Yamanashi###
  male_mae.Yamanashi[ik] = mae(res_male.Yamanashi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[19]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Nagano###
  male_mae.Nagano[ik] = mae(res_male.Nagano[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[20]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Gifu###
  male_mae.Gifu[ik] = mae(res_male.Gifu[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[21]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Shizuoka###
  male_mae.Shizuoka[ik] = mae(res_male.Shizuoka[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[22]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Aichi###
  male_mae.Aichi[ik] = mae(res_male.Aichi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[23]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Mie###
  male_mae.Mie[ik] = mae(res_male.Mie[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[24]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Shiga###
  male_mae.Shiga[ik] = mae(res_male.Shiga[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[25]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Kyoto###
  male_mae.Kyoto[ik] = mae(res_male.Kyoto[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[26]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Osaka###
  male_mae.Osaka[ik] = mae(res_male.Osaka[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[27]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Hyogo###
  male_mae.Hyogo[ik] = mae(res_male.Hyogo[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[28]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Nara###
  male_mae.Nara[ik] = mae(res_male.Nara[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[29]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Wakayama###
  male_mae.Wakayama[ik] = mae(res_male.Wakayama[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[30]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Tottori###
  male_mae.Tottori[ik] = mae(res_male.Tottori[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[31]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Shimane###
  male_mae.Shimane[ik] = mae(res_male.Shimane[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[32]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Okayama###
  male_mae.Okayama[ik] = mae(res_male.Okayama[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[33]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Hiroshima###
  male_mae.Hiroshima[ik] = mae(res_male.Hiroshima[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[34]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Yamaguchi###
  male_mae.Yamaguchi[ik] = mae(res_male.Yamaguchi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[35]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Tokushima###
  male_mae.Tokushima[ik] = mae(res_male.Tokushima[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[36]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Kagawa###
  male_mae.Kagawa[ik] = mae(res_male.Kagawa[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[37]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Ehimae###
  male_mae.Ehime[ik] = mae(res_male.Ehime[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[38]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Kochi###
  male_mae.Kochi[ik] = mae(res_male.Kochi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[39]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Fukuoka###
  male_mae.Fukuoka[ik] = mae(res_male.Fukuoka[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[40]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Saga###
  male_mae.Saga[ik] = mae(res_male.Saga[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[41]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Nagasaki###
  male_mae.Nagasaki[ik] = mae(res_male.Nagasaki[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[42]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Kumamoto###
  male_mae.Kumamoto[ik] = mae(res_male.Kumamoto[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[43]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Oita###
  male_mae.Oita[ik] = mae(res_male.Oita[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[44]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Miyazaki###
  male_mae.Miyazaki[ik] = mae(res_male.Miyazaki[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[45]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Kagoshima###
  male_mae.Kagoshima[ik] = mae(res_male.Kagoshima[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[46]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Okinawa###
  male_mae.Okinawa[ik] = mae(res_male.Okinawa[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[47]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  # rmse 
  #### Hokkaido###
  male_rmse.Hokkaido[ik] = rmse(res_male.Hokkaido[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[1]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Aomori###
  male_rmse.Aomori[ik] = rmse(res_male.Aomori[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[2]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Iwate###
  male_rmse.Iwate[ik] = rmse(res_male.Iwate[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[3]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Miyagi###
  male_rmse.Miyagi[ik] = rmse(res_male.Miyagi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[4]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Akita###
  male_rmse.Akita[ik] = rmse(res_male.Akita[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[5]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Yamagata###
  male_rmse.Yamagata[ik] = rmse(res_male.Yamagata[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[6]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Fukushima###
  male_rmse.Fukushima[ik] = rmse(res_male.Fukushima[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[7]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Ibaraki###
  male_rmse.Ibaraki[ik] = rmse(res_male.Ibaraki[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[8]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Tochigi###
  male_rmse.Tochigi[ik] = rmse(res_male.Tochigi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[9]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Gunma###
  male_rmse.Gunma[ik] = rmse(res_male.Gunma[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[10]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Saitama###
  male_rmse.Saitama[ik] = rmse(res_male.Saitama[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[11]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Chiba###
  male_rmse.Chiba[ik] = rmse(res_male.Chiba[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[12]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Tokyo###
  male_rmse.Tokyo[ik] = rmse(res_male.Tokyo[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[13]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Kanagawa###
  male_rmse.Kanagawa[ik] = rmse(res_male.Kanagawa[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[14]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Niigata###
  male_rmse.Niigata[ik] = rmse(res_male.Niigata[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[15]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Toyama###
  male_rmse.Toyama[ik] = rmse(res_male.Toyama[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[16]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Ishikawa###
  male_rmse.Ishikawa[ik] = rmse(res_male.Ishikawa[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[17]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Fukui###
  male_rmse.Fukui[ik] = rmse(res_male.Fukui[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[18]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Yamanashi###
  male_rmse.Yamanashi[ik] = rmse(res_male.Yamanashi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[19]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Nagano###
  male_rmse.Nagano[ik] = rmse(res_male.Nagano[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[20]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Gifu###
  male_rmse.Gifu[ik] = rmse(res_male.Gifu[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[21]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Shizuoka###
  male_rmse.Shizuoka[ik] = rmse(res_male.Shizuoka[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[22]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Aichi###
  male_rmse.Aichi[ik] = rmse(res_male.Aichi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[23]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Mie###
  male_rmse.Mie[ik] = rmse(res_male.Mie[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[24]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Shiga###
  male_rmse.Shiga[ik] = rmse(res_male.Shiga[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[25]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Kyoto###
  male_rmse.Kyoto[ik] = rmse(res_male.Kyoto[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[26]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Osaka###
  male_rmse.Osaka[ik] = rmse(res_male.Osaka[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[27]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Hyogo###
  male_rmse.Hyogo[ik] = rmse(res_male.Hyogo[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[28]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Nara###
  male_rmse.Nara[ik] = rmse(res_male.Nara[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[29]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Wakayama###
  male_rmse.Wakayama[ik] = rmse(res_male.Wakayama[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[30]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Tottori###
  male_rmse.Tottori[ik] = rmse(res_male.Tottori[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[31]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Shimane###
  male_rmse.Shimane[ik] = rmse(res_male.Shimane[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[32]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Okayama###
  male_rmse.Okayama[ik] = rmse(res_male.Okayama[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[33]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Hiroshima###
  male_rmse.Hiroshima[ik] = rmse(res_male.Hiroshima[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[34]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Yamaguchi###
  male_rmse.Yamaguchi[ik] = rmse(res_male.Yamaguchi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[35]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Tokushima###
  male_rmse.Tokushima[ik] = rmse(res_male.Tokushima[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[36]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Kagawa###
  male_rmse.Kagawa[ik] = rmse(res_male.Kagawa[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[37]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Ehime###
  male_rmse.Ehime[ik] = rmse(res_male.Ehime[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[38]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Kochi###
  male_rmse.Kochi[ik] = rmse(res_male.Kochi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[39]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Fukuoka###
  male_rmse.Fukuoka[ik] = rmse(res_male.Fukuoka[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[40]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Saga###
  male_rmse.Saga[ik] = rmse(res_male.Saga[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[41]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Nagasaki###
  male_rmse.Nagasaki[ik] = rmse(res_male.Nagasaki[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[42]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Kumamoto###
  male_rmse.Kumamoto[ik] = rmse(res_male.Kumamoto[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[43]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Oita###
  male_rmse.Oita[ik] = rmse(res_male.Oita[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[44]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Miyazaki###
  male_rmse.Miyazaki[ik] = rmse(res_male.Miyazaki[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[45]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Kagoshima###
  male_rmse.Kagoshima[ik] = rmse(res_male.Kagoshima[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[46]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
  
  #### Okinawa###
  male_rmse.Okinawa[ik] = rmse(res_male.Okinawa[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[47]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$male)
}


hdfts.jpmale.me<-cbind(male_me.Hokkaido , male_me.Aomori , male_me.Iwate ,
                       male_me.Miyagi , male_me.Akita , male_me.Yamagata ,
                       male_me.Fukushima , male_me.Ibaraki , 
                       male_me.Tochigi , male_me.Gunma , male_me.Saitama,
                       male_me.Chiba ,  male_me.Tokyo, male_me.Kanagawa ,
                       male_me.Niigata , male_me.Toyama , male_me.Ishikawa ,
                       male_me.Fukui , male_me.Yamanashi , male_me.Nagano , 
                       male_me.Gifu , male_me.Shizuoka , male_me.Aichi ,
                       male_me.Mie ,  male_me.Shiga , male_me.Kyoto ,
                       male_me.Osaka , male_me.Hyogo ,  male_me.Nara ,
                       male_me.Wakayama , male_me.Tottori , male_me.Shimane ,
                       male_me.Okayama , male_me.Hiroshima , 
                       male_me.Yamaguchi ,  male_me.Tokushima , 
                       male_me.Kagawa , male_me.Ehime , male_me.Kochi ,
                       male_me.Fukuoka , male_me.Saga , male_me.Nagasaki ,
                       male_me.Kumamoto , male_me.Oita , male_me.Miyazaki ,
                       male_me.Kagoshima , male_me.Okinawa )

hdfts.jpmale.mae<-cbind(male_mae.Hokkaido , male_mae.Aomori , male_mae.Iwate ,
                        male_mae.Miyagi , male_mae.Akita , male_mae.Yamagata ,
                        male_mae.Fukushima , male_mae.Ibaraki , 
                        male_mae.Tochigi , male_mae.Gunma , male_mae.Saitama,
                        male_mae.Chiba ,  male_mae.Tokyo, male_mae.Kanagawa ,
                        male_mae.Niigata , male_mae.Toyama , male_mae.Ishikawa ,
                        male_mae.Fukui , male_mae.Yamanashi , male_mae.Nagano , 
                        male_mae.Gifu , male_mae.Shizuoka , male_mae.Aichi ,
                        male_mae.Mie ,  male_mae.Shiga , male_mae.Kyoto ,
                        male_mae.Osaka , male_mae.Hyogo ,  male_mae.Nara ,
                        male_mae.Wakayama , male_mae.Tottori , male_mae.Shimane ,
                        male_mae.Okayama , male_mae.Hiroshima , 
                        male_mae.Yamaguchi ,  male_mae.Tokushima , 
                        male_mae.Kagawa , male_mae.Ehime , male_mae.Kochi ,
                        male_mae.Fukuoka , male_mae.Saga , male_mae.Nagasaki ,
                        male_mae.Kumamoto , male_mae.Oita , male_mae.Miyazaki ,
                        male_mae.Kagoshima , male_mae.Okinawa )
hdfts.jpmale.rmse<-cbind(male_rmse.Hokkaido , male_rmse.Aomori , male_rmse.Iwate ,
                         male_rmse.Miyagi , male_rmse.Akita , male_rmse.Yamagata ,
                         male_rmse.Fukushima , male_rmse.Ibaraki , 
                         male_rmse.Tochigi , male_rmse.Gunma , male_rmse.Saitama,
                         male_rmse.Chiba ,  male_rmse.Tokyo, male_rmse.Kanagawa ,
                         male_rmse.Niigata , male_rmse.Toyama , male_rmse.Ishikawa ,
                         male_rmse.Fukui , male_rmse.Yamanashi , male_rmse.Nagano , 
                         male_rmse.Gifu , male_rmse.Shizuoka , male_rmse.Aichi ,
                         male_rmse.Mie ,  male_rmse.Shiga , male_rmse.Kyoto ,
                         male_rmse.Osaka , male_rmse.Hyogo ,  male_rmse.Nara ,
                         male_rmse.Wakayama , male_rmse.Tottori , male_rmse.Shimane ,
                         male_rmse.Okayama , male_rmse.Hiroshima , 
                         male_rmse.Yamaguchi ,  male_rmse.Tokushima , 
                         male_rmse.Kagawa , male_rmse.Ehime , male_rmse.Kochi ,
                         male_rmse.Fukuoka , male_rmse.Saga , male_rmse.Nagasaki ,
                         male_rmse.Kumamoto , male_rmse.Oita , male_rmse.Miyazaki ,
                         male_rmse.Kagoshima , male_rmse.Okinawa )
colnames(hdfts.jpmale.me)<-c("Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima", "Ibaraki",  "Tochigi", "Gunma", "Saitama",  
                             "Chiba", "Tokyo", "Kanagawa", "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",    
                             "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi",
                             "Tokushima", "Kagawa", "Ehime", "Kochi", "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")

colnames(hdfts.jpmale.mae)<-c("Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima", "Ibaraki",  "Tochigi", "Gunma", "Saitama",  
                              "Chiba", "Tokyo", "Kanagawa", "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",    
                              "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi",
                              "Tokushima", "Kagawa", "Ehime", "Kochi", "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")
colnames(hdfts.jpmale.rmse)<-c("Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima", "Ibaraki",  "Tochigi", "Gunma", "Saitama",  
                               "Chiba", "Tokyo", "Kanagawa", "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",    
                               "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi",
                               "Tokushima", "Kagawa", "Ehime", "Kochi", "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")


write.csv(hdfts.jpmale.me, "hdfts.jpmale.me.csv")
write.csv(hdfts.jpmale.mae, "hdfts.jpmale.mae.csv")
write.csv(hdfts.jpmale.rmse, "hdfts.jpmale.rmse.csv")



