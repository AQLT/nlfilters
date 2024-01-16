library(rjd3filters)
library(patchwork)
hf <- rjd3filters::lp_filter()

y_obs <- readRDS("data/CE16OV.RDS")
y_obs <- y_obs[[length(y_obs)]]
y_simul <- readRDS("data/mediumvariability2.RDS")
y_simul <- y_simul[[length(y_simul)]]

h_f <- y_obs * hf

med_f <- pracma::hampel(y_obs, 6, 0)
med_f <- med_f$y
med_f <- ts(runmed(y_obs, 9), start = start(y_obs), frequency = frequency(y_obs))
hp_f <- pracma::hampel(y_obs, 6)

data <- ts.union(y_obs, h_f, med_f, hp_f$y)
colnames(data) <- c("y", "Henderson", "Median filter", "Hampel")

AQLTools::graph_ts(window(data, start = 2000, end = c(2001,12)),size = 1) /
AQLTools::graph_ts(window(data, start = 2019, end = c(2021,12)),size = 1)

AQLTools::hc_stocks(window(data, start = 2000, end = c(2001,12)))
plot(y_obs)
plot(window(y_obs, start = (2001+1/12) - 6/12, end = (2001+1/12) + 6/12))
abline(h = median(window(y_obs, start = (2001+1/12) - 6/12, end = (2001+1/12) + 6/12)))
