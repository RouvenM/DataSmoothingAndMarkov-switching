# Replication JRSSC Paper
library(dplyr)
library(ggplot2)
testiid <- data.frame(index = seq(1, 1500, 1), iid = rnorm(1500, 0, 0.2))
# testiid <- data.frame(index = seq(1, 1500, 1), iid = rep(NA, 1500), state = rep(NA, 1500))
# Gamma <- matrix(c(0.9, 0.1, 0.1, 0.9), 2, 2)
# testiid$state[1] <- sample(1:2, 1, 0.5)
# testiid$iid[1] <- rnorm(1, (0.5-testiid$state[1])/5, 0.3) 
# for (i in 2:nrow(testiid)) {
#   testiid$state[i] <- sample(1:2, size = 1, prob = Gamma[testiid$state[i - 1], ])
#   testiid$iid[i] <- rnorm(1, (0.5-testiid$state[1])/5, 0.3)
# }
# 
plot(testiid$iid, col = "grey")
testiid$kernel <- ksmooth(testiid$index, testiid$iid, kernel = c("normal"), bandwidth = 14,
                    range.x = range(x),
                    n.points = length(x), testiid$index)
ksmooth(testiid$index, testiid$iid, kernel = c("normal"), bandwidth = 14,
        range.x = range(x),
        n.points = length(x), testiid$index) %>% lines(., col = "red", add = T)


sign1 <- matrix(0, nrow = nrow(testiid), ncol = 10)
for (i in 1:10) {
  sign1[,i] <- c(0, sign(diff(testiid$kernel$y))* 
                  abs(diff(testiid$kernel$y))^(i/10))
}

sign1 <- as.data.frame(sign1)
library(tidyr)
sign2 <- gather(sign1)
i <- 1
hist(sign1[,i])
plot(sign1[,i])

mllk <- function(theta.star, x){
  theta <- c(plogis(theta.star[1]), plogis(theta.star[2]),
             (theta.star[3]), (theta.star[4]),
             exp(theta.star[5]), exp(theta.star[6]))
  Gamma <- diag(theta[1:2])
  Gamma[1, 2] <- 1 - Gamma[1, 1]
  Gamma[2, 1] <- 1 - Gamma[2, 2]
  delta <- solve(t(diag(2) - Gamma + 1), c(1, 1))
  mu <- theta[3:4]
  sigma <- (theta[5:6])
  allprobs <- matrix(1, length(x), 2)
  ind <- which(!is.na(x))
  allprobs[ind, ] <- cbind(
    dnorm(x[ind], mu[1], sigma[1]),
    dnorm(x[ind], mu[2], sigma[2]))
  foo <- delta %*% diag(allprobs[1, ])
  l <- log(sum(foo))
  phi <- foo / sum(foo)
  for (t in 2:length(x)){
    foo <- phi %*% Gamma %*% diag(allprobs[t, ])
    l <- l + log(sum(foo))
    phi <- foo / sum(foo)
  }
  return(-l)
}
theta.star <- c(qlogis(c(0.5, 0.5)), -0.1, 0.1, 0, 0)

amfactor <- rep(0, 10)
mod <- list()
for (i in 1:10) {
mod[[i]] <- nlm(mllk, theta.star, x = sign1[,i], print.level = 2, iterlim = 10000, stepmax = 10)
est <- mod[[i]]$estimate
amfactor[i] <- abs(sign(est[3])*abs(est[3])^(1/(i/10))-sign(est[4])*abs(est[4])^(1/(i/10)))
}

plot(x = (1:10)/10, y = amfactor, type = "l")

i = 10
theta.star <- mod[[i]]$estimate
theta <- c(plogis(theta.star[1]), plogis(theta.star[2]),
           (theta.star[3]), (theta.star[4]),
           exp(theta.star[5]), exp(theta.star[6]))
Gamma <- diag(theta[1:2])
Gamma[1, 2] <- 1 - Gamma[1, 1]
Gamma[2, 1] <- 1 - Gamma[2, 2]
delta <- solve(t(diag(2) - Gamma + 1), c(1, 1))
mu <- theta[3:4]
sigma <- (theta[5:6])
Gamma; mu; sigma

# it seems to be the case that the serial correlation is induced by chance

# Replication -------------------------------------------------------------
library(dplyr)
library(BasketballAnalyzeR)
data(package="BasketballAnalyzeR")
PbP <- PbPmanipulation(PbP.BDB)
PbP$madeormiss <- ifelse(PbP$result == "missed", 0, 1)


GS <- PbP %>% filter(event_type %in% c("free throw", "shot", "miss"), team == "GSW") %>% 
  group_by(game_id) %>% summarise(nr_shots = n(), oppTeam = unique(oppTeam))
KD <- PbP %>% filter(player == "Kevin Durant", event_type %in% c("free throw", "shot", "miss")) %>% 
  dplyr::select("game_id", "date", "elapsed", "remaining_time", 
         "event_type", "points", "result", "ShotType", "oppTeam", "totalTime")
KD <- KD %>% group_by(game_id) %>% mutate(elapsedshot = c(totalTime[1], diff(totalTime))) %>% 
              filter(elapsedshot != 0)

KD$elapsedshot[which(KD$elapsedshot == 0)] <- 1
KD$shoot_int <- 1/KD$elapsedshot

KD2 <- merge(KD, GS, by = c("game_id"))

KD2$adj_sh_int <- KD2$shoot_int / (KD2$nr_shots/3600)
KD2$madeormiss <- ifelse(KD2$result == "missed", 0, 1)

shot_perc <- KD2 %>% group_by(game_id, ShotType) %>% summarise(shot_perc = mean(madeormiss))
#shot_perc2 <- PbP %>% filter(event_type %in% c("free throw", "shot", "miss"), team == "GSW") %>% 
#  group_by(game_id, ShotType) %>% summarise(shot_perc = mean(madeormiss))
KD3 <- merge(KD2, shot_perc, by = c("game_id", "ShotType")) %>% arrange(date, totalTime) %>% 
  mutate(index = 1:n())
#KD3 <- merge(KD2, shot_perc2, by = c("game_id", "ShotType")) %>% arrange(date, totalTime) %>% 
#  mutate(index = 1:n())

KD3$eij <- KD3$madeormiss-KD3$shot_perc
KD3$shoot_perf <- KD3$eij*KD3$adj_sh_int

shots <- KD3 %>% group_by(game_id) %>% summarise(shots = n()) %>% pull(shots)
perc <- quantile(shots, 0.15)

#  ylim(c(-0.5, 0.5))
kernel <- ksmooth(KD3$index, KD3$shoot_perf, kernel = c("normal"), 
        bandwidth = perc,
        range.x = range(x),
        n.points = length(x), KD3$index)

ggplot(KD3, aes(x = index, y = shoot_perf)) + geom_point(aes(col = ShotType)) + 
  geom_line(y = kernel$y, size = 1.2) 

sign1 <- matrix(0, nrow = length(kernel$y), ncol = 20)
for (i in 1:20) {
  sign1[,i] <- c(0, sign(diff(kernel$y))* 
                   abs(diff(kernel$y))^(i/10))
}

sign1 <- as.data.frame(sign1)
nr <- 5
hist(sign1[,nr])
plot(sign1[,nr], type = "l")

mllk <- function(theta.star, x){
  theta <- c(plogis(theta.star[1]), plogis(theta.star[2]),
             (theta.star[3]), (theta.star[4]),
             (theta.star[5]), (theta.star[6]))
  Gamma <- diag(theta[1:2])
  Gamma[1, 2] <- 1 - Gamma[1, 1]
  Gamma[2, 1] <- 1 - Gamma[2, 2]
  delta <- solve(t(diag(2) - Gamma + 1), c(1, 1))
  mu <- theta[3:4]
  sigma <- exp(theta[5:6])
  allprobs <- matrix(1, length(x), 2)
  ind <- which(!is.na(x))
  allprobs[ind, ] <- cbind(
    dnorm(x[ind], mu[1], sigma[1]),
    dnorm(x[ind], mu[2], sigma[2]))
  foo <- delta %*% diag(allprobs[1, ])
  l <- log(sum(foo))
  phi <- foo / sum(foo)
  for (t in 2:length(x)){
    foo <- phi %*% Gamma %*% diag(allprobs[t, ])
    l <- l + log(sum(foo))
    phi <- foo / sum(foo)
  }
  return(-l)
}
theta.star <- c(qlogis(c(0.8, 0.8)), -0.02, 0.02, 0, 0)
i = 1
mod[[i]] <- nlm(mllk, theta.star, x = sign1[,i], print.level = 2, iterlim = 10000, stepmax = 150)


amfactor <- rep(0, 20)
mod <- list()
for (i in 1:20) {
  mod[[i]] <- nlm(mllk, theta.star, x = sign1[,i], print.level = 2, iterlim = 10000, stepmax = 150)
  est <- mod[[i]]$estimate
  amfactor[i] <- (sign(est[3])*abs(est[3])^(1/(i/10))-sign(est[4])*abs(est[4])^(1/(i/10)))
}

plot(x = (1:20)/10, y = amfactor, type = "l")
