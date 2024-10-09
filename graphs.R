library(rjd3filters)
library(patchwork)
library(robfilter)
library(ggplot2)
library(AQLThesis)
export_dir <- "/Users/alainquartierlatente/Desktop/Carriere/Slides/2024 - 06 - CSI/img"
labels <- c("Henderson" = "Henderson", 
            "Hendout1"  = "Hend. out. m",
            "Hendout2"  = "Hend. out. m-(m+1)",
            "LQD" = "Least Quartile Difference" , 
            "RM" = "Repeated Median", 
            "LMS" = "Least Median of Squares", 
            "LTS" = "Least Trimmed Squares", 
            "DR" = "Deepest Regression", 
            "MED" = "Median")

set.seed(100)
start = 1960
frequency = 12
time = seq_along(seq(start, 2019+11/12, by = 1/12))
series_simul  = list(
  highvariability1 = simulated_tci(time,sigma_nu = 0.08,sigma_e = 0.40,lambda = 72,rho = 0.5),
  highvariability2 = simulated_tci(time,sigma_nu = 0.08,sigma_e = 0.40,lambda = 72,rho = 0.7),
  highvariability3 = simulated_tci(time,sigma_nu = 0.08,sigma_e = 0.40,lambda = 72,rho = 1),
  mediumvariability1 = simulated_tci(time,sigma_nu = 0.08,sigma_e = 0.30,lambda = 72,rho = 1.5),
  mediumvariability2 = simulated_tci(time,sigma_nu = 0.08,sigma_e = 0.30,lambda = 72,rho = 2),
  mediumvariability3 = simulated_tci(time,sigma_nu = 0.08,sigma_e = 0.30,lambda = 72,rho = 3),
  lowvariability1 = simulated_tci(time,sigma_nu = 0.08,sigma_e = 0.20,lambda = 72,rho = 3),
  lowvariability2 = simulated_tci(time,sigma_nu = 0.08,sigma_e = 0.20,lambda = 72,rho = 3.5),
  lowvariability3 = simulated_tci(time,sigma_nu = 0.08,sigma_e = 0.20,lambda = 72,rho = 4)
)
series_simul = lapply(series_simul,ts, start = start, frequency = frequency)
tc_simul <- do.call(cbind, lapply(series_simul,  function(x) x[, "trend"] +x[, "cycle"]))

dim <- c(5.5, 3.8)*1.5
gen_MM <- function(p=6, q=p, d=2, reg = NULL){
  X_gen <- function(d = 1, p = 6, q = p){
    sapply(0:d, function(exp) seq(-p, q)^exp)
  }
  k = rjd3filters::get_kernel("Henderson", h = p)
  k = c(rev(k$coef[-1]), k$coef[seq(0,q)+1])
  K = diag(k)
  X = cbind(X_gen(d=d, p = p, q = q), reg)
  e1 = matrix(0, ncol = 1, nrow = ncol(X))
  e1[1] = 1
  # Estimator of the constant
  M1 = K %*% X %*% solve(t(X) %*% K %*% X, e1)
  moving_average(M1, lags = -p)
}
hf <- rjd3filters::lp_filter()
hf_s <- hf@sfilter
y_obs <- readRDS("data/CE16OV.RDS")
y_obs <- y_obs[[length(y_obs)]]
y_simul_low <- readRDS("data/lowvariability2.RDS")
y_simul_low <- y_simul_low[[length(y_simul_low)]]
y_simul_med <- readRDS("data/mediumvariability2.RDS")
y_simul_med <- y_simul_med[[length(y_simul_med)]]
y_simul_high <- readRDS("data/highvariability2.RDS")
y_simul_high <- y_simul_high[[length(y_simul_high)]]

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

h_f_out1 <- 
  h_f_out2 <- 
  y_obs * hf_s
for (i in 1:13) {
  tmp_f <- hf_s
  tmp_f[13 - i + 1] <- 0
  tmp_f <- tmp_f / sum(tmp_f)
  current_date <- 2020+ (4-1)/12 + (i-7)/12
  new_value <- window(y_obs * tmp_f, start = current_date, end = current_date)
  window(h_f_out1, start = current_date, end = current_date) <- 
    new_value
  
  current_date <- 2001.16666666667 + (i-7)/12
  new_value <- window(y_obs * tmp_f, start = current_date, end = current_date)
  window(h_f_out1, start = current_date, end = current_date) <- 
    new_value
}
for (i in 1:14) {
  tmp_f <- hf_s
  tmp_f[max(13 - i + 1, 1):min(13 - i + 2, 13)] <- 0
  tmp_f <- tmp_f / sum(tmp_f)
  current_date <- 2020+ (4-1)/12 + (i-7)/12
  new_value <- window(y_obs * tmp_f, start = current_date, end = current_date)
  window(h_f_out2, start = current_date, end = current_date) <- 
    new_value
  
  current_date <- 2001.16666666667 + (i-7)/12
  new_value <- window(y_obs * tmp_f, start = current_date, end = current_date)
  window(h_f_out2, start = current_date, end = current_date) <- 
    new_value
}
data_obs_graph <- ts.union(data_obs[,1], h_f_out1, h_f_out2,
                           data_obs[, seq.int(ncol(data_obs), 2)])
colnames(data_obs_graph) <- c("Henderson", "Hendout1", "Hendout2",
                              colnames(data_obs)[seq.int(ncol(data_obs), 2)])

# data_obs_graph <- window(data_obs_graph, start = c(2019,10), end = (2021))
all_plot <- lapply(seq_len(ncol(data_obs_graph)), function(i) {
  data_plot <- data.frame(
    x = as.numeric(time(data_obs_graph)),
    y = as.numeric(window(y_obs, start = start(data_obs_graph), end = end(data_obs_graph))),
    mu = as.numeric(data_obs_graph[,i])
  )
  data_plot <- na.omit(data_plot)
  ggplot(data = data_plot, aes(x = x)) +
    geom_line(aes(y = y), color = "black", lty = 1) +
    geom_line(aes(y = mu), color = "orange") +
    labs(title = labels[colnames(data_obs_graph)[i]],
         x = NULL, y = NULL) +
    scale_x_continuous(labels = function(x) as.character(zoo::as.yearmon(x)),
                       breaks = scales::pretty_breaks(n = 3))
})
c_tp <- 2020+ (4-1)/12
p <- wrap_plots(all_plot, ncol = 3) & theme_bw() &
  theme(panel.grid.minor = element_blank()) &
  coord_cartesian(xlim = c(c_tp-.5, c_tp+.5),
                  ylim = range(window(y_obs, start = c_tp-4/12, end = c_tp+.5))) &
  geom_vline(xintercept = c_tp, lty = 2, color = "gray") 
p |> 
  ggsave(filename = file.path(export_dir, "obs_covid.pdf"), width = dim[1], height = dim[2])
c_tp <- 2001  + 1/12
(
  wrap_plots(all_plot, ncol = 3) & theme_bw() &
    theme(panel.grid.minor = element_blank()) &
    coord_cartesian(
      xlim = c(c_tp-.5, c_tp+.5),
      ylim = range(window(y_obs, start = c_tp-4/12, end = c_tp+.5))
      ) &
    geom_vline(xintercept = c_tp, lty = 2, color = "gray")
) |> 
  ggsave(filename = file.path(export_dir, "obs_2001.pdf"), width = dim[1], height = dim[2])



data_simul <- estim_rob(y_simul_med)
tp <- 1972.583
tp <- 1978.667
AQLThesis::turning_points(data_simul[,1])
data_simul_graph <- ts.union(data_simul[,1],
                             data_simul[, seq.int(ncol(data_simul), 2)])
colnames(data_simul_graph) <- c("Henderson", 
                                colnames(data_simul)[seq.int(ncol(data_simul), 2)])
data_simul_graph <- window(data_simul_graph, start = tp-1, end = tp+1)
all_plot <- lapply(seq_len(ncol(data_simul_graph)), function(i) {
  data_plot <- data.frame(
    x = as.numeric(time(data_simul_graph)),
    y = as.numeric(window(y_simul_med, start = start(data_simul_graph), end = end(data_simul_graph))),
    mu = as.numeric(data_simul_graph[,i])
  )
  ggplot(data = data_plot, aes(x = x)) +
    geom_vline(xintercept = tp, lty = 2, color = "gray") +
    geom_line(aes(y = y), color = "black", lty = 1) +
    geom_line(aes(y = mu), color = "orange") +
    labs(title = labels[colnames(data_simul_graph)[i]],
         x = NULL, y = NULL) +
    scale_x_continuous(labels = function(x) as.character(zoo::as.yearmon(x)))
})
p <- wrap_plots(all_plot, ncol = 3) & theme_bw() &
  theme(panel.grid.minor = element_blank())

p |> 
  ggsave(filename = file.path(export_dir, "simul_med.pdf"), width = dim[1], height = dim[2])
data_simul <- estim_rob(y_simul_high)
tp <- 1972.500
data_simul_graph <- ts.union(data_simul[,1],
                             data_simul[, seq.int(ncol(data_simul), 2)])
colnames(data_simul_graph) <- c("Henderson", 
                                colnames(data_simul)[seq.int(ncol(data_simul), 2)])
data_simul_graph <- window(data_simul_graph, start = tp-2, end = tp+2)
all_plot <- lapply(seq_len(ncol(data_simul_graph)), function(i) {
  data_plot <- data.frame(
    x = as.numeric(time(data_simul_graph)),
    y = as.numeric(window(y_simul_high, start = start(data_simul_graph), end = end(data_simul_graph))),
    mu = as.numeric(data_simul_graph[,i])
  )
  ggplot(data = data_plot, aes(x = x)) +
    geom_vline(xintercept = tp, lty = 2, color = "gray") +
    geom_line(aes(y = y), color = "black", lty = 1) +
    geom_line(aes(y = mu), color = "orange") +
    labs(title = labels[colnames(data_simul_graph)[i]],
         x = NULL, y = NULL) +
    scale_x_continuous(labels = function(x) as.character(zoo::as.yearmon(x)))
})
p <- wrap_plots(all_plot, ncol = 3) & theme_bw() &
  theme(panel.grid.minor = element_blank())
p |> 
  ggsave(filename = file.path(export_dir, "simul_high.pdf"), width = dim[1], height = dim[2])


data_simul <- estim_rob(y_simul_low)

AQLThesis::turning_points(data_simul[,1])
tp <- 1972.58333333333
data_simul_graph <- ts.union(data_simul[,1],
                             data_simul[, seq.int(ncol(data_simul), 2)])
colnames(data_simul_graph) <- c("Henderson", 
                                colnames(data_simul)[seq.int(ncol(data_simul), 2)])
data_simul_graph <- window(data_simul_graph, start = tp-7/12, end = tp+1)
all_plot <- lapply(seq_len(ncol(data_simul_graph)), function(i) {
  data_plot <- data.frame(
    x = as.numeric(time(data_simul_graph)),
    y = as.numeric(window(y_simul_low, start = start(data_simul_graph), end = end(data_simul_graph))),
    mu = as.numeric(data_simul_graph[,i])
  )
  ggplot(data = data_plot, aes(x = x)) +
    geom_vline(xintercept = tp, lty = 2, color = "gray") +
    geom_line(aes(y = y), color = "black", lty = 1) +
    geom_line(aes(y = mu), color = "orange") +
    labs(title = labels[colnames(data_simul_graph)[i]],
         x = NULL, y = NULL) +
    scale_x_continuous(labels = function(x) as.character(zoo::as.yearmon(x)))
})
p <- wrap_plots(all_plot, ncol = 3) & theme_bw() &
  theme(panel.grid.minor = element_blank())
p |> 
  ggsave(filename = file.path(export_dir, "simul_low.pdf"), width = dim[1], height = dim[2])


c(upturn1 = 1963.58333333333, upturn2 = 1969.58333333333, upturn3 = 1975.58333333333, 
  upturn4 = 1981.58333333333, upturn5 = 1987.58333333333, upturn6 = 1993.58333333333, 
  upturn7 = 1999.58333333333, upturn8 = 2005.58333333333, upturn9 = 2011.58333333333, 
  upturn10 = 2017.58333333333, downturn1 = 1960.58333333333, downturn2 = 1966.58333333333, 
  downturn3 = 1972.58333333333, downturn4 = 1978.58333333333, downturn5 = 1984.58333333333, 
  downturn6 = 1990.58333333333, downturn7 = 1996.58333333333, downturn8 = 2002.58333333333, 
  downturn9 = 2008.58333333333, downturn10 = 2014.58333333333)
