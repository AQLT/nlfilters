library(rjd3filters)
library(patchwork)
library(robfilter)
hf <- rjd3filters::lp_filter()

y_obs <- readRDS("data/CE16OV.RDS")
y_obs <- y_obs[[length(y_obs)]]
y_simul <- readRDS("data/mediumvariability2.RDS")
y_simul <- y_simul[[length(y_simul)]]

estim_rob <- function(x, length = 13) {
  data_rob <- robreg.filter(x, length)
  data_rob_level <- ts(data_rob$level,
                       start = start(x),
                       frequency = frequency(x))
  h_f <- x * hf
  res <- ts.union(h_f, data_rob_level)
  colnames(res) <- c("Henderson", colnames(data_rob_level))
  res
}

data_obs <- estim_rob(y_obs)
data_simul <- estim_rob(y_simul)
AQLTools::hc_stocks(window(data_obs))
AQLTools::hc_stocks(window(data_obs, start = c(2020), end = c(2020, 12)))
AQLTools::graph_ts(window(estim_rob(y_obs, 5), start = 2020, end = c(2020, 12)))

AQLTools::hc_stocks(window(estim_rob(y_simul, 25), start = 1966, end = 1968))

h_f <- y_obs * hf

med_f <- pracma::hampel(y_obs, 6, 0)
med_f <- med_f$y
med_f <- ts(runmed(y_obs, 9), start = start(y_obs), frequency = frequency(y_obs))
hp_f <- pracma::hampel(y_obs, 6)

data <- ts.union(y_obs, h_f, med_f, hp_f$y)
colnames(data) <- c("y", "Henderson", "Median filter", "Hampel")

data_rob <- robreg.filter(y_obs, 13)
data_rob_level <- ts(data_rob$level,
                     start = start(y_obs),
                     frequency = frequency(y_obs))
data <- ts.union(data, data_rob_level)
colnames(data) <- gsub("(data.)|(data_rob_level.)", "", colnames(data))
AQLTools::hc_stocks(window(data, start = 2000, end = c(2001,12)))

AQLTools::graph_ts(window(data, start = 2000, end = c(2001,12)),size = 1) /
AQLTools::graph_ts(window(data, start = 2019, end = c(2021,12)),size = 1)

AQLTools::hc_stocks(window(data, start = 2000, end = c(2001,12)))
AQLTools::hc_stocks(data)


h_f <- y_simul * hf

hp_f <- pracma::hampel(y_simul, 6)

data <- ts.union(y_simul, h_f, hp_f$y)
colnames(data) <- c("y", "Henderson", "Hampel")

data_rob <- robreg.filter(y_obs, 13)
data_rob_level <- ts(data_rob$level,
                     start = start(y_obs),
                     frequency = frequency(y_obs))
data_rob_slope <- ts(data_rob$slope,
                     start = start(y_obs),
                     frequency = frequency(y_obs))
data <- ts.union(data, data_rob_level)
colnames(data) <- gsub("(data.)|(data_rob_level.)", "", colnames(data))
AQLTools::hc_stocks(data)
AQLTools::hc_stocks(data_rob_level)
AQLTools::hc_stocks(data_rob_slope)
lms.filter

robreg.filter



AQLTools::hc_stocks(data)