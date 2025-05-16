library(demography)

#############################################
# read data from Japanese Mortality Database
#############################################
read.jpn <- function (region,  label = region) 
{
  path_death <- paste("https://www.ipss.go.jp/p-toukei/JMD/", region, "/STATS/",   "Mx_1x1.txt", sep = "")
  mx <- read.csv(path_death, skip = 2, header = TRUE, sep = "")
  
  path_pop <- paste("https://www.ipss.go.jp/p-toukei/JMD/", region, "/STATS/",   "Exposures_1x1.txt", sep = "")
  pop <- read.csv(path_pop, skip = 2, header = TRUE, sep = "")

  obj <- list(type = "mortality", label = label, lambda = 0)
  obj$year <- sort(unique(mx[, 1]))
  n <- length(obj$year)
  m <- length(unique(mx[, 2]))
  obj$age <- mx[1:m, 2]
  mnames <- names(mx)[-c(1, 2)]
  n.mort <- length(mnames)
  obj$rate <- obj$pop <- list()
  for (i in 1:n.mort) 
  {
    obj$rate[[i]] <- matrix(as.numeric(mx[, i + 2]), nrow = m, ncol = n)
    obj$rate[[i]][obj$rate[[i]] < 0] <- NA
    obj$pop[[i]] <- matrix(as.numeric(pop[, i + 2]), nrow = m, ncol = n)
    obj$pop[[i]][obj$pop[[i]] < 0] <- NA
    dimnames(obj$rate[[i]]) <- dimnames(obj$pop[[i]]) <- list(obj$age, obj$year)
  }
  names(obj$pop) = (names(obj$rate) <- tolower(mnames))
  obj$age <- as.numeric(as.character(obj$age))
  if (is.na(obj$age[m])) 
    obj$age[m] <- 2 * obj$age[m - 1] - obj$age[m - 2]
  return(structure(obj, class = "demogdata"))
}



#Read Data#
jpn.data<-list()
jpn.data[[1]] = read.jpn("01", "Hokkaido")
jpn.data[[2]] = read.jpn("02", "Aomori")
jpn.data[[3]] = read.jpn("03", "Iwate")
jpn.data[[4]] = read.jpn("04", "Miyagi")
jpn.data[[5]] = read.jpn("05", "Akita")
jpn.data[[6]] = read.jpn("06", "Yamagata")
jpn.data[[7]] = read.jpn("07", "Fukushima")
jpn.data[[8]] = read.jpn("08", "Ibaraki")
jpn.data[[9]] = read.jpn("09", "Tochigi")
jpn.data[[10]] = read.jpn("10", "Gunma")
jpn.data[[11]] = read.jpn("11", "Saitama")
jpn.data[[12]] = read.jpn("12", "Chiba")
jpn.data[[13]] = read.jpn("13", "Tokyo")
jpn.data[[14]] = read.jpn("14", "Kanagawa")
jpn.data[[15]] = read.jpn("15", "Niigata")
jpn.data[[16]] = read.jpn("16", "Toyama")
jpn.data[[17]] = read.jpn("17", "Ishikawa")
jpn.data[[18]]  = read.jpn("18", "Fukui")
jpn.data[[19]] = read.jpn("19", "Yamanashi")
jpn.data[[20]] = read.jpn("20", "Nagano")
jpn.data[[21]] = read.jpn("21", "Gifu")
jpn.data[[22]] = read.jpn("22", "Shizuoka")
jpn.data[[23]] = read.jpn("23", "Aichi")
jpn.data[[24]] = read.jpn("24", "Mie")
jpn.data[[25]] = read.jpn("25", "Shiga")
jpn.data[[26]] = read.jpn("26", "Kyoto")
jpn.data[[27]] = read.jpn("27", "Osaka")
jpn.data[[28]] = read.jpn("28", "Hyogo")
jpn.data[[29]] = read.jpn("29", "Nara")
jpn.data[[30]] = read.jpn("30", "Wakayama")
jpn.data[[31]] = read.jpn("31", "Tottori")
jpn.data[[32]] = read.jpn("32", "Shimane")
jpn.data[[33]] = read.jpn("33", "Okayama")
jpn.data[[34]] = read.jpn("34", "Hiroshima")
jpn.data[[35]] = read.jpn("35", "Yamaguchi")
jpn.data[[36]] = read.jpn("36", "Tokushima")
jpn.data[[37]] = read.jpn("37", "Kagawa")
jpn.data[[38]]  = read.jpn("38", "Ehime")
jpn.data[[39]] = read.jpn("39", "Kochi")
jpn.data[[40]] = read.jpn("40", "Fukuoka")
jpn.data[[41]] = read.jpn("41", "Saga")
jpn.data[[42]] = read.jpn("42", "Nagasaki")
jpn.data[[43]] = read.jpn("43", "Kumamoto")
jpn.data[[44]] = read.jpn("44", "Oita")
jpn.data[[45]] = read.jpn("45", "Miyazaki")
jpn.data[[46]] = read.jpn("46", "Kagoshima")
jpn.data[[47]] = read.jpn("47", "Okinawa")




