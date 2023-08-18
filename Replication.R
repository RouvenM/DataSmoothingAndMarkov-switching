
# Getting the real data ---------------------------------------------------

# install.packages(c("BasketballAnalyzeR", "tidyverse"))
library(BasketballAnalyzeR)
library(tidyverse)

# this function creates a data set, that contains the variables that are described by Sandri et al. (2020)
create_data = function(p){
  PbP.BDB = BasketballAnalyzeR::PbP.BDB
  PbP = BasketballAnalyzeR::PbPmanipulation(PbP.BDB)
  team = (as.character(PbP$team)[which(as.character(PbP$player) == p)][1])
  data = PbP %>% 
    filter(player == p) %>% 
    filter(event_type == "shot" | event_type == "miss" | event_type == "free throw") %>% 
    dplyr::select(game_id, a1, a2, a3, a4, a5, h1, h2, h3, h4, h5, totalTime, player, result, ShotType) %>% 
    distinct(.keep_all = T)
  # number of shots und time auf game level
  S = PbP %>% 
    filter(event_type == "shot" | event_type == "miss" | event_type == "free throw") %>%
    filter(team == team) %>% 
    group_by(game_id) %>% 
    summarise(S = n())
  Time = PbP %>% 
    filter(event_type == "shot" | event_type == "miss" | event_type == "free throw") %>%
    group_by(game_id) %>% 
    summarise(Time = max(totalTime))
  data = data %>% left_join(S, by = "game_id") %>% left_join(Time, by = "game_id")
  data$intensity_u = c(NA, 1/diff(data$totalTime))
  data = data %>% mutate(intensity = intensity_u/(S/Time))
  # calculate scoring probs
  scoring_probs = data %>% 
    group_by(game_id) %>% 
    summarise(p2 = sum(as.character(ShotType) == "2P" & as.character(result) == "made")/sum(as.character(ShotType) == "2P"),
              p3 = sum(as.character(ShotType) == "3P" & as.character(result) == "made")/sum(as.character(ShotType) == "3P"),
              ft = sum(as.character(ShotType) == "FT" & as.character(result) == "made")/sum(as.character(ShotType) == "FT")) %>% 
    ungroup()
  
  data = data %>% left_join(scoring_probs, by = "game_id")
  data$made = 0
  data$made[as.character(data$result) == "made"] = 1
  data = data %>% filter(intensity != Inf) %>% 
    mutate(ShotType = as.character(ShotType)) %>% 
    mutate(result = as.character(result))
  # compute efficiency
  data$efficiency = 0
  p2_ind = which(data$ShotType == "2P")
  p3_ind = which(data$ShotType == "3P")
  ft_ind = which(data$ShotType == "FT")
  data$efficiency[p2_ind] = data$made[p2_ind] - data$p2[p2_ind]
  data$efficiency[p3_ind] = data$made[p3_ind] - data$p3[p3_ind]
  data$efficiency[ft_ind] = data$made[ft_ind] - data$ft[ft_ind]
  # calculate shooting performance
  data$shooting_perf = data$efficiency*data$intensity
  data$color[data$ShotType == "2P"] = "#00b4c9"
  data$color[data$ShotType == "3P"] = "#e6c467"
  data$color[data$ShotType == "FT"] = "#e6896a"
  return(data)
}

data = create_data("Kevin Durant")
shots = data %>% group_by(game_id) %>% summarise(shots = n()) %>% pull(shots)
perc = quantile(shots, 0.15)

n = nrow(data)

data$shooting_perf_smooth = ksmooth(1:n, data$shooting_perf, kernel = "normal", bandwidth = perc, x.points = 1:n)$y



# Creating iid data -------------------------------------------------------

set.seed(666)

data_sim = rnorm(1400, 0, 0.18)
data_sim_smooth = ksmooth(1:n, data_sim, kernel = "normal", bandwidth = perc, x.points = 1:n)$y


# Analysis and Figures from the papaer ------------------------------------

## Figure 1

par(mfrow = c(1,2))
par(mar = c(5,4,4,1.2))
plot(data$shooting_perf, pch = 20, ylim = c(-1.5,1.5),
     main = "Original data", ylab = "Shooting performance", bty = "n")
lines(data$shooting_perf_smooth, lwd = 1.5, col = "orange")
plot(data_sim, pch = 20, ylim = c(-1.5,1.5), main = "Simulated data", ylab = "Shooting performance", bty = "n")
lines(data_sim_smooth, col = "orange", lwd = 1.5)

# it is evident that the real-data distribution has heavier tails, but this is ofter lesser importance here.


## Comparing ACF of real and simulated data before and after smoothing (this is not included in the paper)

par(mfrow = c(2,2))
acf(data$shooting_perf, main = "Shooting performance of Kevin Durant")
acf(data_sim, main = "White noise")
acf(data$shooting_perf_smooth, main = "Smoothed shooting performance")
acf(data_sim_smooth, main = "Smoothed white noise")



# Finding the optimal factor of k -----------------------------------------

# functions to amplify the data and compute the mean difference as described by Sandri et al. (2020)
amplification = function(data, k){
  sign(diff(data))*abs(diff(data))^k 
}
mean_diff = function(m1, m2, k){
 sign(m1)*abs(m1)^(1/k) - sign(m2)*abs(m2)^(1/k)
}

library(depmixS4)

L = 50 # grid of 50 values
l = 20 # each model is fitted 20 times to avoid local maxima
k = seq(0.01, 1.4, length.out = L)
m1 = m2 = numeric(L) # save the state-dependent means for all models
m_diff = numeric(L) # save mean differences for all models

for (i in 1:L){
  d = amplification(data$shooting_perf_smooth, k = k[i]) # amplify data
  data$response = c(NA, d)
  
  mods = list()
  llks = numeric(l)
  for (j in 1:l){
    m = depmix(response = response ~ 1, data = data, nstates = 2)
    mods[[j]] = tryCatch(fit(m), error = function(e) e)
    if(is(mods[[j]], "error")){
      llks[j] = NA
    } else{
      llks[j] = logLik(mods[[j]])
    }
  }
  
  pars = getpars(mods[[which.max(llks)]])
  mu = c(pars[7], pars[9])
  # order state-dependent means
  m1[i] = max(mu)
  m2[i] = min(mu)
  
  m_diff[i] = mean_diff(m1[i], m2[i], k = k[i])
}

par(mfrow = c(1,1))
plot(k, m_diff, type = "l", lwd = 2, bty = "n")
(kstar = k[which.max(m_diff)]) # optimal k following the criterion by Sandri et al. (2020)



# Fitting HMMs with optimal k ---------------------------------------------

# we use direct numerical maximum likelihood estimation here
mllk = function(theta.star, x){ # minus-log-likelihood
  Gamma = diag(plogis(theta.star[1:2]))
  Gamma[1, 2] = 1 - Gamma[1, 1]
  Gamma[2, 1] = 1 - Gamma[2, 2]
  delta = solve(t(diag(2) - Gamma + 1), c(1, 1), tol = 1e-30)
  mu = theta.star[3:4]
  sigma = exp(theta.star[5:6])
  allprobs = matrix(1, length(x), 2)
  ind = which(!is.na(x))
  allprobs[ind, ] = cbind(dnorm(x[ind], mu[1], sigma[1]),
                          dnorm(x[ind], mu[2], sigma[2]))
  foo = delta %*% diag(allprobs[1, ])
  l = log(sum(foo))
  phi = foo / sum(foo)
  for (t in 2:length(x)){
    foo = phi %*% Gamma %*% diag(allprobs[t, ])
    l = l + log(sum(foo))
    phi = foo / sum(foo)
  }
  return(-l)
}

theta.star = c(qlogis(c(0.7, 0.7)), -0.01, 0.01, log(c(0.005, 0.005)))

mod_real = nlm(mllk, theta.star, x = amplification(data$shooting_perf_smooth, kstar), print.level = 2)
mod_sim = nlm(mllk, theta.star, x = amplification(data_sim_smooth, kstar), print.level = 2)

# parameters of real-data model
theta.star = mod_real$estimate
Gamma_real = diag(plogis(theta.star[1:2]))
Gamma_real[1, 2] = 1 - Gamma_real[1, 1]
Gamma_real[2, 1] = 1 - Gamma_real[2, 2]
delta_real = solve(t(diag(2) - Gamma_real + 1), c(1, 1), tol = 1e-30)
mu_real = theta.star[3:4]
sigma_real = exp(theta.star[5:6])
# parameters of simulated-data model
theta.star = mod_sim$estimate
Gamma_sim = diag(plogis(theta.star[1:2]))
Gamma_sim[1, 2] = 1 - Gamma_sim[1, 1]
Gamma_sim[2, 1] = 1 - Gamma_sim[2, 2]
delta_sim = solve(t(diag(2) - Gamma_sim + 1), c(1, 1), tol = 1e-30)
mu_sim = theta.star[3:4]
sigma_sim = exp(theta.star[5:6])

# parameter estimates of real-data model
round(Gamma_real,2)
round(mu_real, 4)


## Figure 2
par(mfrow = c(1,2))
par(mar = c(5,4,4,1.2))

plot(amplification(data$shooting_perf_smooth, kstar), type = "l", lwd = 1.2, ylim = c(-0.04, 0.04), main = "Original data", bty = "n", ylab = "Amplified differences")
abline(h = mu_real[1], col = "deepskyblue", lwd = 2)
abline(h = mu_real[2], col = "orange", lwd = 2)

plot(amplification(data_sim_smooth, kstar), type = "l", lwd = 1.2, ylim = c(-0.04, 0.04), main = "Simulated data", bty = "n", ylab = "Amplified differences")
abline(h = mu_sim[1], col = "deepskyblue", lwd = 2)
abline(h = mu_sim[2], col = "orange", lwd = 2)


## Figure 3
par(mfrow = c(1,2))
par(mar = c(5,4,4,1.2))

hist(amplification(data$shooting_perf_smooth, kstar), prob = T, main = "Original data", xlab = "Amplified differences", breaks = 30, xlim = c(-0.04, 0.04), border = "white")
curve(delta_real[1]*dnorm(x, mu_real[1], sigma_real[1]), add = T, lwd = 3, col = "orange", n = 500)
curve(delta_real[2]*dnorm(x, mu_real[2], sigma_real[2]), add = T, lwd = 3, col = "deepskyblue", n = 500)

hist(amplification(data_sim_smooth, kstar), prob = T, main = "Simulated data", xlab = "Amplified differences", breaks = 30, xlim = c(-0.04, 0.04), border = "white")
curve(delta_sim[1]*dnorm(x, mu_sim[1], sigma_sim[1]), add = T, lwd = 3, col = "orange", n = 500)
curve(delta_sim[2]*dnorm(x, mu_sim[2], sigma_sim[2]), add = T, lwd = 3, col = "deepskyblue", n = 500)



# Model without preprocessing ---------------------------------------------

theta.star = c(qlogis(c(0.7, 0.7)), -0.1, 0.1, log(c(0.1, 0.1)))
mod_wo = nlm(mllk, theta.star, x = diff(data$shooting_perf), print.level = 2)

theta.star = mod_wo$estimate
Gamma = diag(plogis(theta.star[1:2]))
Gamma[1, 2] = 1 - Gamma[1, 1]
Gamma[2, 1] = 1 - Gamma[2, 2]
delta = solve(t(diag(2) - Gamma + 1), c(1, 1), tol = 1e-30)
mu = theta.star[3:4]
sigma = exp(theta.star[5:6])

par(mfrow = c(1,1))
hist(diff(data$shooting_perf), prob = T, main = "Analysis without pre-processing", xlab = "First-order differences", breaks = 30, border = "white")
curve(delta[1]*dnorm(x, mu[1], sigma[1]), add = T, lwd = 3, col = "orange", n = 500)
curve(delta[2]*dnorm(x, mu[2], sigma[2]), add = T, lwd = 3, col = "deepskyblue", n = 500)

# withouth preprocessing we do not find state separation regarging the means of the state dependent
# distributions, i.e. between increasing and decreasing performance, but only between phases of different performance variation.



