######Forecast 
year_horizon=10
forecast_year = (year_horizon + 1)
n_year = 2022 - forecast_year
res_female.Hokkaido = res_female.Aomori = res_female.Iwate = res_female.Miyagi = 
  res_female.Akita = res_female.Yamagata = res_female.Fukushima = 
  res_female.Ibaraki = res_female.Tochigi = res_female.Gunma = 
  res_female.Saitama = res_female.Chiba = res_female.Tokyo= res_female.Kanagawa = res_female.Niigata = res_female.Toyama = res_female.Ishikawa = res_female.Fukui = res_female.Yamanashi = res_female.Nagano = res_female.Gifu = 
  res_female.Shizuoka = res_female.Aichi = res_female.Mie = res_female.Shiga = res_female.Kyoto = res_female.Osaka = res_female.Hyogo = res_female.Nara = res_female.Wakayama = res_female.Tottori = res_female.Shimane = 
  res_female.Okayama = res_female.Hiroshima = res_female.Yamaguchi = 
  res_female.Tokushima = res_female.Kagawa = res_female.Ehime = 
  res_female.Kochi = res_female.Fukuoka = res_female.Saga = res_female.Nagasaki = res_female.Kumamoto = res_female.Oita = res_female.Miyazaki = 
  res_female.Kagoshima = res_female.Okinawa =array(NA, dim = c(year_horizon, 101, year_horizon))

for(ij in 1:year_horizon)
{
  jpn.trunc<-list()
  for ( i in 1:47){
    jpn.trunc[[i]]<-smooth.demogdata(extract.ages(extract.years(jpn.data[[i]], years = 1973:(n_year+ij)), 0:100, combine.upper=TRUE))
  }
  jpn.truc<-list()
  for (ik in 1:47){
    jpn.truc[[ik]]<-t(log(jpn.trunc[[ik]]$rate$female))
  }
  res_forc<-ffm.forc(jpn.truc, year_horizon = year_horizon, years = 1973:(n_year+ij), ages = 0:100, n_pop =47)
  res_female.Hokkaido[,,ij] = res_forc[[1]]
  res_female.Aomori[,,ij] = res_forc[[2]]
  res_female.Iwate[,,ij] = res_forc[[3]]
  res_female.Miyagi[,,ij] = res_forc[[4]]
  res_female.Akita [,,ij] = res_forc[[5]]
  res_female.Yamagata[,,ij] = res_forc[[6]]
  res_female.Fukushima[,,ij] = res_forc[[7]]
  res_female.Ibaraki[,,ij] = res_forc[[8]]
  res_female.Tochigi[,,ij] = res_forc[[9]]
  res_female.Gunma[,,ij] = res_forc[[10]]
  res_female.Saitama[,,ij] = res_forc[[11]]
  res_female.Chiba[,,ij] = res_forc[[12]]
  res_female.Tokyo[,,ij] = res_forc[[13]]
  res_female.Kanagawa[,,ij] = res_forc[[14]]
  res_female.Niigata[,,ij] = res_forc[[15]]
  res_female.Toyama[,,ij] = res_forc[[16]]
  res_female.Ishikawa[,,ij] = res_forc[[17]]
  res_female.Fukui[,,ij] = res_forc[[18]]
  res_female.Yamanashi[,,ij] = res_forc[[19]]
  res_female.Nagano[,,ij] = res_forc[[20]]
  res_female.Gifu[,,ij] = res_forc[[21]]
  res_female.Shizuoka[,,ij] = res_forc[[22]]
  res_female.Aichi[,,ij] = res_forc[[23]]
  res_female.Mie[,,ij] = res_forc[[24]]
  res_female.Shiga[,,ij] = res_forc[[25]]
  res_female.Kyoto[,,ij] = res_forc[[26]]
  res_female.Osaka[,,ij] = res_forc[[27]]
  res_female.Hyogo[,,ij] = res_forc[[28]]
  res_female.Nara[,,ij] = res_forc[[29]]
  res_female.Wakayama[,,ij] = res_forc[[30]]
  res_female.Tottori[,,ij] = res_forc[[31]]
  res_female.Shimane[,,ij] = res_forc[[32]]
  res_female.Okayama[,,ij] = res_forc[[33]]
  res_female.Hiroshima[,,ij] = res_forc[[34]]
  res_female.Yamaguchi[,,ij] = res_forc[[35]]
  res_female.Tokushima[,,ij] = res_forc[[36]]
  res_female.Kagawa[,,ij] = res_forc[[37]]
  res_female.Ehime [,,ij] = res_forc[[38]]
  res_female.Kochi[,,ij] = res_forc[[39]]
  res_female.Fukuoka[,,ij] = res_forc[[40]]
  res_female.Saga[,,ij] = res_forc[[41]]
  res_female.Nagasaki[,,ij] = res_forc[[42]]
  res_female.Kumamoto[,,ij] = res_forc[[43]]
  res_female.Oita[,,ij] = res_forc[[44]]
  res_female.Miyazaki[,,ij] = res_forc[[45]]
  res_female.Kagoshima[,,ij] = res_forc[[46]]
  res_female.Okinawa[,,ij] = res_forc[[47]]
}


me   = ftsa:::me
mae  = ftsa:::mae
rmse = ftsa:::rmse
# MAE & RMSE
female_me.Hokkaido = female_me.Aomori = female_me.Iwate = female_me.Miyagi = 
  female_me.Akita = female_me.Yamagata = female_me.Fukushima = female_me.Ibaraki = 
  female_me.Tochigi = female_me.Gunma = female_me.Saitama= female_me.Chiba = 
  female_me.Tokyo= female_me.Kanagawa = female_me.Niigata = female_me.Toyama = 
  female_me.Ishikawa = female_me.Fukui = female_me.Yamanashi = female_me.Nagano = 
  female_me.Gifu = female_me.Shizuoka = female_me.Aichi = female_me.Mie = 
  female_me.Shiga = female_me.Kyoto = female_me.Osaka = female_me.Hyogo = 
  female_me.Nara = female_me.Wakayama = female_me.Tottori = female_me.Shimane = 
  female_me.Okayama = female_me.Hiroshima = female_me.Yamaguchi = female_me.Tokushima = 
  female_me.Kagawa = female_me.Ehime = female_me.Kochi = female_me.Fukuoka = 
  female_me.Saga = female_me.Nagasaki = female_me.Kumamoto = female_me.Oita =   
  female_me.Miyazaki = female_me.Kagoshima = female_me.Okinawa =
  female_mae.Hokkaido = female_mae.Aomori = female_mae.Iwate = female_mae.Miyagi = 
  female_mae.Akita = female_mae.Yamagata = female_mae.Fukushima = female_mae.Ibaraki = 
  female_mae.Tochigi = female_mae.Gunma = female_mae.Saitama= female_mae.Chiba = 
  female_mae.Tokyo= female_mae.Kanagawa = female_mae.Niigata = female_mae.Toyama = 
  female_mae.Ishikawa = female_mae.Fukui = female_mae.Yamanashi = female_mae.Nagano = 
  female_mae.Gifu = female_mae.Shizuoka = female_mae.Aichi = female_mae.Mie = 
  female_mae.Shiga = female_mae.Kyoto = female_mae.Osaka = female_mae.Hyogo = 
  female_mae.Nara = female_mae.Wakayama = female_mae.Tottori = female_mae.Shimane = 
  female_mae.Okayama = female_mae.Hiroshima = female_mae.Yamaguchi = female_mae.Tokushima = 
  female_mae.Kagawa = female_mae.Ehime = female_mae.Kochi = female_mae.Fukuoka = 
  female_mae.Saga = female_mae.Nagasaki = female_mae.Kumamoto = female_mae.Oita =   
  female_mae.Miyazaki = female_mae.Kagoshima = female_mae.Okinawa = 
  female_rmse.Hokkaido = female_rmse.Aomori = female_rmse.Iwate = female_rmse.Miyagi = 
  female_rmse.Akita = female_rmse.Yamagata = female_rmse.Fukushima = female_rmse.Ibaraki = 
  female_rmse.Tochigi = female_rmse.Gunma = female_rmse.Saitama= female_rmse.Chiba = 
  female_rmse.Tokyo= female_rmse.Kanagawa = female_rmse.Niigata = female_rmse.Toyama = 
  female_rmse.Ishikawa = female_rmse.Fukui = female_rmse.Yamanashi = female_rmse.Nagano = 
  female_rmse.Gifu = female_rmse.Shizuoka = female_rmse.Aichi = female_rmse.Mie = 
  female_rmse.Shiga = female_rmse.Kyoto = female_rmse.Osaka = female_rmse.Hyogo = 
  female_rmse.Nara = female_rmse.Wakayama = female_rmse.Tottori = female_rmse.Shimane = 
  female_rmse.Okayama = female_rmse.Hiroshima = female_rmse.Yamaguchi = 
  female_rmse.Tokushima = female_rmse.Kagawa = female_rmse.Ehime = female_rmse.Kochi = female_rmse.Fukuoka = female_rmse.Saga = female_rmse.Nagasaki = female_rmse.Kumamoto = female_rmse.Oita = female_rmse.Miyazaki = female_rmse.Kagoshima = female_rmse.Okinawa = 
  vector("numeric",year_horizon)

for(ik in 1:year_horizon)
{
  # me 
  #### Hokkaido###
  female_me.Hokkaido[ik] = me(res_female.Hokkaido[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[1]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  
  #### Aomori###
  female_me.Aomori[ik] = me(res_female.Aomori[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[2]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Iwate###
  female_me.Iwate[ik] = me(res_female.Iwate[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[3]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Miyagi###
  female_me.Miyagi[ik] = me(res_female.Miyagi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[4]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Akita###
  female_me.Akita[ik] = me(res_female.Akita[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[5]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Yamagata###
  female_me.Yamagata[ik] = me(res_female.Yamagata[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[6]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Fukushima###
  female_me.Fukushima[ik] = me(res_female.Fukushima[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[7]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Ibaraki###
  female_me.Ibaraki[ik] = me(res_female.Ibaraki[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[8]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Tochigi###
  female_me.Tochigi[ik] = me(res_female.Tochigi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[9]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Gunma###
  female_me.Gunma[ik] = me(res_female.Gunma[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[10]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Saitama###
  female_me.Saitama[ik] = me(res_female.Saitama[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[11]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Chiba###
  female_me.Chiba[ik] = me(res_female.Chiba[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[12]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Tokyo###
  female_me.Tokyo[ik] = me(res_female.Tokyo[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[13]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Kanagawa###
  female_me.Kanagawa[ik] = me(res_female.Kanagawa[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[14]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Niigata###
  female_me.Niigata[ik] = me(res_female.Niigata[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[15]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Toyama###
  female_me.Toyama[ik] = me(res_female.Toyama[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[16]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Ishikawa###
  female_me.Ishikawa[ik] = me(res_female.Ishikawa[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[17]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Fukui###
  female_me.Fukui[ik] = me(res_female.Fukui[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[18]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Yamanashi###
  female_me.Yamanashi[ik] = me(res_female.Yamanashi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[19]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Nagano###
  female_me.Nagano[ik] = me(res_female.Nagano[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[20]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Gifu###
  female_me.Gifu[ik] = me(res_female.Gifu[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[21]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Shizuoka###
  female_me.Shizuoka[ik] = me(res_female.Shizuoka[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[22]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Aichi###
  female_me.Aichi[ik] = me(res_female.Aichi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[23]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Mie###
  female_me.Mie[ik] = me(res_female.Mie[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[24]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Shiga###
  female_me.Shiga[ik] = me(res_female.Shiga[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[25]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Kyoto###
  female_me.Kyoto[ik] = me(res_female.Kyoto[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[26]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Osaka###
  female_me.Osaka[ik] = me(res_female.Osaka[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[27]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Hyogo###
  female_me.Hyogo[ik] = me(res_female.Hyogo[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[28]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Nara###
  female_me.Nara[ik] = me(res_female.Nara[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[29]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Wakayama###
  female_me.Wakayama[ik] = me(res_female.Wakayama[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[30]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Tottori###
  female_me.Tottori[ik] = me(res_female.Tottori[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[31]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Shimane###
  female_me.Shimane[ik] = me(res_female.Shimane[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[32]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Okayama###
  female_me.Okayama[ik] = me(res_female.Okayama[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[33]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Hiroshima###
  female_me.Hiroshima[ik] = me(res_female.Hiroshima[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[34]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Yamaguchi###
  female_me.Yamaguchi[ik] = me(res_female.Yamaguchi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[35]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Tokushima###
  female_me.Tokushima[ik] = me(res_female.Tokushima[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[36]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Kagawa###
  female_me.Kagawa[ik] = me(res_female.Kagawa[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[37]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Ehime###
  female_me.Ehime[ik] = me(res_female.Ehime[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[38]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Kochi###
  female_me.Kochi[ik] = me(res_female.Kochi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[39]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Fukuoka###
  female_me.Fukuoka[ik] = me(res_female.Fukuoka[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[40]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Saga###
  female_me.Saga[ik] = me(res_female.Saga[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[41]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Nagasaki###
  female_me.Nagasaki[ik] = me(res_female.Nagasaki[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[42]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Kumamoto###
  female_me.Kumamoto[ik] = me(res_female.Kumamoto[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[43]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Oita###
  female_me.Oita[ik] = me(res_female.Oita[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[44]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Miyazaki###
  female_me.Miyazaki[ik] = me(res_female.Miyazaki[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[45]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Kagoshima###
  female_me.Kagoshima[ik] = me(res_female.Kagoshima[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[46]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Okinawa###
  female_me.Okinawa[ik] = me(res_female.Okinawa[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[47]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  # mae 
  #### Hokkaido###
  female_mae.Hokkaido[ik] = mae(res_female.Hokkaido[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[1]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Aomori###
  female_mae.Aomori[ik] = mae(res_female.Aomori[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[2]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Iwate###
  female_mae.Iwate[ik] = mae(res_female.Iwate[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[3]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Miyagi###
  female_mae.Miyagi[ik] = mae(res_female.Miyagi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[4]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Akita###
  female_mae.Akita[ik] = mae(res_female.Akita[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[5]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Yamagata###
  female_mae.Yamagata[ik] = mae(res_female.Yamagata[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[6]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Fukushima###
  female_mae.Fukushima[ik] = mae(res_female.Fukushima[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[7]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Ibaraki###
  female_mae.Ibaraki[ik] = mae(res_female.Ibaraki[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[8]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Tochigi###
  female_mae.Tochigi[ik] = mae(res_female.Tochigi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[9]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Gunma###
  female_mae.Gunma[ik] = mae(res_female.Gunma[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[10]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Saitama###
  female_mae.Saitama[ik] = mae(res_female.Saitama[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[11]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Chiba###
  female_mae.Chiba[ik] = mae(res_female.Chiba[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[12]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Tokyo###
  female_mae.Tokyo[ik] = mae(res_female.Tokyo[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[13]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Kanagawa###
  female_mae.Kanagawa[ik] = mae(res_female.Kanagawa[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[14]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Niigata###
  female_mae.Niigata[ik] = mae(res_female.Niigata[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[15]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Toyama###
  female_mae.Toyama[ik] = mae(res_female.Toyama[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[16]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Ishikawa###
  female_mae.Ishikawa[ik] = mae(res_female.Ishikawa[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[17]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Fukui###
  female_mae.Fukui[ik] = mae(res_female.Fukui[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[18]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Yamanashi###
  female_mae.Yamanashi[ik] = mae(res_female.Yamanashi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[19]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Nagano###
  female_mae.Nagano[ik] = mae(res_female.Nagano[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[20]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Gifu###
  female_mae.Gifu[ik] = mae(res_female.Gifu[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[21]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Shizuoka###
  female_mae.Shizuoka[ik] = mae(res_female.Shizuoka[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[22]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Aichi###
  female_mae.Aichi[ik] = mae(res_female.Aichi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[23]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Mie###
  female_mae.Mie[ik] = mae(res_female.Mie[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[24]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Shiga###
  female_mae.Shiga[ik] = mae(res_female.Shiga[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[25]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Kyoto###
  female_mae.Kyoto[ik] = mae(res_female.Kyoto[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[26]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Osaka###
  female_mae.Osaka[ik] = mae(res_female.Osaka[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[27]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Hyogo###
  female_mae.Hyogo[ik] = mae(res_female.Hyogo[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[28]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Nara###
  female_mae.Nara[ik] = mae(res_female.Nara[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[29]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Wakayama###
  female_mae.Wakayama[ik] = mae(res_female.Wakayama[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[30]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Tottori###
  female_mae.Tottori[ik] = mae(res_female.Tottori[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[31]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Shimane###
  female_mae.Shimane[ik] = mae(res_female.Shimane[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[32]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Okayama###
  female_mae.Okayama[ik] = mae(res_female.Okayama[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[33]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Hiroshima###
  female_mae.Hiroshima[ik] = mae(res_female.Hiroshima[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[34]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Yamaguchi###
  female_mae.Yamaguchi[ik] = mae(res_female.Yamaguchi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[35]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Tokushima###
  female_mae.Tokushima[ik] = mae(res_female.Tokushima[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[36]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Kagawa###
  female_mae.Kagawa[ik] = mae(res_female.Kagawa[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[37]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Ehimae###
  female_mae.Ehime[ik] = mae(res_female.Ehime[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[38]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Kochi###
  female_mae.Kochi[ik] = mae(res_female.Kochi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[39]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Fukuoka###
  female_mae.Fukuoka[ik] = mae(res_female.Fukuoka[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[40]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Saga###
  female_mae.Saga[ik] = mae(res_female.Saga[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[41]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Nagasaki###
  female_mae.Nagasaki[ik] = mae(res_female.Nagasaki[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[42]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Kumamoto###
  female_mae.Kumamoto[ik] = mae(res_female.Kumamoto[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[43]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Oita###
  female_mae.Oita[ik] = mae(res_female.Oita[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[44]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Miyazaki###
  female_mae.Miyazaki[ik] = mae(res_female.Miyazaki[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[45]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Kagoshima###
  female_mae.Kagoshima[ik] = mae(res_female.Kagoshima[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[46]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Okinawa###
  female_mae.Okinawa[ik] = mae(res_female.Okinawa[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[47]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  # rmse 
  #### Hokkaido###
  female_rmse.Hokkaido[ik] = rmse(res_female.Hokkaido[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[1]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Aomori###
  female_rmse.Aomori[ik] = rmse(res_female.Aomori[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[2]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Iwate###
  female_rmse.Iwate[ik] = rmse(res_female.Iwate[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[3]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Miyagi###
  female_rmse.Miyagi[ik] = rmse(res_female.Miyagi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[4]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Akita###
  female_rmse.Akita[ik] = rmse(res_female.Akita[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[5]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Yamagata###
  female_rmse.Yamagata[ik] = rmse(res_female.Yamagata[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[6]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Fukushima###
  female_rmse.Fukushima[ik] = rmse(res_female.Fukushima[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[7]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Ibaraki###
  female_rmse.Ibaraki[ik] = rmse(res_female.Ibaraki[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[8]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Tochigi###
  female_rmse.Tochigi[ik] = rmse(res_female.Tochigi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[9]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Gunma###
  female_rmse.Gunma[ik] = rmse(res_female.Gunma[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[10]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Saitama###
  female_rmse.Saitama[ik] = rmse(res_female.Saitama[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[11]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Chiba###
  female_rmse.Chiba[ik] = rmse(res_female.Chiba[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[12]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Tokyo###
  female_rmse.Tokyo[ik] = rmse(res_female.Tokyo[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[13]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Kanagawa###
  female_rmse.Kanagawa[ik] = rmse(res_female.Kanagawa[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[14]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Niigata###
  female_rmse.Niigata[ik] = rmse(res_female.Niigata[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[15]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Toyama###
  female_rmse.Toyama[ik] = rmse(res_female.Toyama[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[16]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Ishikawa###
  female_rmse.Ishikawa[ik] = rmse(res_female.Ishikawa[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[17]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Fukui###
  female_rmse.Fukui[ik] = rmse(res_female.Fukui[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[18]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Yamanashi###
  female_rmse.Yamanashi[ik] = rmse(res_female.Yamanashi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[19]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Nagano###
  female_rmse.Nagano[ik] = rmse(res_female.Nagano[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[20]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Gifu###
  female_rmse.Gifu[ik] = rmse(res_female.Gifu[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[21]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Shizuoka###
  female_rmse.Shizuoka[ik] = rmse(res_female.Shizuoka[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[22]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Aichi###
  female_rmse.Aichi[ik] = rmse(res_female.Aichi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[23]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Mie###
  female_rmse.Mie[ik] = rmse(res_female.Mie[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[24]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Shiga###
  female_rmse.Shiga[ik] = rmse(res_female.Shiga[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[25]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Kyoto###
  female_rmse.Kyoto[ik] = rmse(res_female.Kyoto[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[26]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Osaka###
  female_rmse.Osaka[ik] = rmse(res_female.Osaka[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[27]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Hyogo###
  female_rmse.Hyogo[ik] = rmse(res_female.Hyogo[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[28]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Nara###
  female_rmse.Nara[ik] = rmse(res_female.Nara[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[29]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Wakayama###
  female_rmse.Wakayama[ik] = rmse(res_female.Wakayama[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[30]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Tottori###
  female_rmse.Tottori[ik] = rmse(res_female.Tottori[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[31]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Shimane###
  female_rmse.Shimane[ik] = rmse(res_female.Shimane[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[32]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Okayama###
  female_rmse.Okayama[ik] = rmse(res_female.Okayama[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[33]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Hiroshima###
  female_rmse.Hiroshima[ik] = rmse(res_female.Hiroshima[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[34]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Yamaguchi###
  female_rmse.Yamaguchi[ik] = rmse(res_female.Yamaguchi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[35]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Tokushima###
  female_rmse.Tokushima[ik] = rmse(res_female.Tokushima[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[36]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Kagawa###
  female_rmse.Kagawa[ik] = rmse(res_female.Kagawa[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[37]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Ehime###
  female_rmse.Ehime[ik] = rmse(res_female.Ehime[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[38]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Kochi###
  female_rmse.Kochi[ik] = rmse(res_female.Kochi[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[39]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Fukuoka###
  female_rmse.Fukuoka[ik] = rmse(res_female.Fukuoka[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[40]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Saga###
  female_rmse.Saga[ik] = rmse(res_female.Saga[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[41]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Nagasaki###
  female_rmse.Nagasaki[ik] = rmse(res_female.Nagasaki[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[42]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Kumamoto###
  female_rmse.Kumamoto[ik] = rmse(res_female.Kumamoto[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[43]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Oita###
  female_rmse.Oita[ik] = rmse(res_female.Oita[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[44]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Miyazaki###
  female_rmse.Miyazaki[ik] = rmse(res_female.Miyazaki[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[45]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Kagoshima###
  female_rmse.Kagoshima[ik] = rmse(res_female.Kagoshima[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[46]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
  
  #### Okinawa###
  female_rmse.Okinawa[ik] = rmse(res_female.Okinawa[ik,,1:(forecast_year-ik)], extract.ages(extract.years(jpn.data[[47]], years = (2012+ik):2022), 0:100, combine.upper=TRUE)$rate$female)
}


ffm.jpfemale.me<-cbind(female_me.Hokkaido , female_me.Aomori , female_me.Iwate ,
                       female_me.Miyagi , female_me.Akita , female_me.Yamagata ,
                       female_me.Fukushima , female_me.Ibaraki , 
                       female_me.Tochigi , female_me.Gunma , female_me.Saitama,
                       female_me.Chiba ,  female_me.Tokyo, female_me.Kanagawa ,
                       female_me.Niigata , female_me.Toyama , female_me.Ishikawa ,
                       female_me.Fukui , female_me.Yamanashi , female_me.Nagano , 
                       female_me.Gifu , female_me.Shizuoka , female_me.Aichi ,
                       female_me.Mie ,  female_me.Shiga , female_me.Kyoto ,
                       female_me.Osaka , female_me.Hyogo ,  female_me.Nara ,
                       female_me.Wakayama , female_me.Tottori , female_me.Shimane ,
                       female_me.Okayama , female_me.Hiroshima , 
                       female_me.Yamaguchi ,  female_me.Tokushima , 
                       female_me.Kagawa , female_me.Ehime , female_me.Kochi ,
                       female_me.Fukuoka , female_me.Saga , female_me.Nagasaki ,
                       female_me.Kumamoto , female_me.Oita , female_me.Miyazaki ,
                       female_me.Kagoshima , female_me.Okinawa )

ffm.jpfemale.mae<-cbind(female_mae.Hokkaido , female_mae.Aomori , female_mae.Iwate ,
                        female_mae.Miyagi , female_mae.Akita , female_mae.Yamagata ,
                        female_mae.Fukushima , female_mae.Ibaraki , 
                        female_mae.Tochigi , female_mae.Gunma , female_mae.Saitama,
                        female_mae.Chiba ,  female_mae.Tokyo, female_mae.Kanagawa ,
                        female_mae.Niigata , female_mae.Toyama , female_mae.Ishikawa ,
                        female_mae.Fukui , female_mae.Yamanashi , female_mae.Nagano , 
                        female_mae.Gifu , female_mae.Shizuoka , female_mae.Aichi ,
                        female_mae.Mie ,  female_mae.Shiga , female_mae.Kyoto ,
                        female_mae.Osaka , female_mae.Hyogo ,  female_mae.Nara ,
                        female_mae.Wakayama , female_mae.Tottori , female_mae.Shimane ,
                        female_mae.Okayama , female_mae.Hiroshima , 
                        female_mae.Yamaguchi ,  female_mae.Tokushima , 
                        female_mae.Kagawa , female_mae.Ehime , female_mae.Kochi ,
                        female_mae.Fukuoka , female_mae.Saga , female_mae.Nagasaki ,
                        female_mae.Kumamoto , female_mae.Oita , female_mae.Miyazaki ,
                        female_mae.Kagoshima , female_mae.Okinawa )
ffm.jpfemale.rmse<-cbind(female_rmse.Hokkaido , female_rmse.Aomori , female_rmse.Iwate ,
                         female_rmse.Miyagi , female_rmse.Akita , female_rmse.Yamagata ,
                         female_rmse.Fukushima , female_rmse.Ibaraki , 
                         female_rmse.Tochigi , female_rmse.Gunma , female_rmse.Saitama,
                         female_rmse.Chiba ,  female_rmse.Tokyo, female_rmse.Kanagawa ,
                         female_rmse.Niigata , female_rmse.Toyama , female_rmse.Ishikawa ,
                         female_rmse.Fukui , female_rmse.Yamanashi , female_rmse.Nagano , 
                         female_rmse.Gifu , female_rmse.Shizuoka , female_rmse.Aichi ,
                         female_rmse.Mie ,  female_rmse.Shiga , female_rmse.Kyoto ,
                         female_rmse.Osaka , female_rmse.Hyogo ,  female_rmse.Nara ,
                         female_rmse.Wakayama , female_rmse.Tottori , female_rmse.Shimane ,
                         female_rmse.Okayama , female_rmse.Hiroshima , 
                         female_rmse.Yamaguchi ,  female_rmse.Tokushima , 
                         female_rmse.Kagawa , female_rmse.Ehime , female_rmse.Kochi ,
                         female_rmse.Fukuoka , female_rmse.Saga , female_rmse.Nagasaki ,
                         female_rmse.Kumamoto , female_rmse.Oita , female_rmse.Miyazaki ,
                         female_rmse.Kagoshima , female_rmse.Okinawa )
colnames(ffm.jpfemale.me)<-c("Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima", "Ibaraki",  "Tochigi", "Gunma", "Saitama",  
                             "Chiba", "Tokyo", "Kanagawa", "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",    
                             "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi",
                             "Tokushima", "Kagawa", "Ehime", "Kochi", "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")

colnames(ffm.jpfemale.mae)<-c("Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima", "Ibaraki",  "Tochigi", "Gunma", "Saitama",  
                              "Chiba", "Tokyo", "Kanagawa", "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",    
                              "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi",
                              "Tokushima", "Kagawa", "Ehime", "Kochi", "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")
colnames(ffm.jpfemale.rmse)<-c("Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima", "Ibaraki",  "Tochigi", "Gunma", "Saitama",  
                               "Chiba", "Tokyo", "Kanagawa", "Niigata", "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",    
                               "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", "Tottori", "Shimane", "Okayama", "Hiroshima", "Yamaguchi",
                               "Tokushima", "Kagawa", "Ehime", "Kochi", "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")

write.csv(ffm.jpfemale.me, "ffm.jpfemale.me.csv")
write.csv(ffm.jpfemale.mae, "ffm.jpfemale.mae.csv")
write.csv(ffm.jpfemale.rmse, "ffm.jpfemale.rmse.csv")




