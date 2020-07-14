## Data analysis ##
## V. 3 ##
## Andrea Guido ##
rm(list = ls())
# required libraries
library(haven)
library(psych)
library(ggplot2)
library(dplyr);library(xtable) # xtable is for latex tables
library(rethinking)

# define directories
root <- "C:/Users/andrea/Dropbox/Werk - Projects/Meta - Cooperation RT"
setwd(root)
dir_table <- "Output/Tables_R/"
dir_figure <- "Output/Figures_R/"

# IMPORT DATA
source("ancillary_scripts/data_import_bayes_analysis.R")

# PRODUCE DESCRIPTIVES
source("ancillary_scripts/descriptives.R")

# REGRESSIONS --------------------------------------------
## 1. MODELS USING SVO ANGLE --------------------------------
## these models do not consider all the studies for which the SVO angle is not provided in the original data
## as of now, Yamagishi et al. 2017 is the only one.
## M1: define basic model for prior predictive simulations  ----
m1 <- ulam(
  alist(
    coop ~ dnorm(mu, sigma),
    mu <- a + b*(rt_std),
    a ~ dnorm(50, 10),
    b ~ dnorm(0, 10),
    sigma ~ dexp(1)
  ),
  data = dat, chains = 4, iter = 2000, warmup = 1000, cores = 4, log_lik = T
)

## analyze Prior Predictive Simulations
prior <- extract.prior(m1)
mu <- link(m1, post=prior, data = list(rt_std=c(-3,3)))
plot( NULL, xlim=c(-3,3), ylim=c(0,100), main= "Prior Predictive Simulations")
for (i in 1:200) lines(c(-3,3), mu[i,] , col=col.alpha("blue", 0.4))

### analize posterior estimates
precis(m1,2)
p <- extract.samples(m1)
dens(p$b, show.HPDI = .95); HPDI(p$b, 0.95)
mu_p <- link(m1, post = p, data=list(rt_std=c(-3,3)))
plot( NULL, xlim=c(-3,3), ylim=c(0,100), main = "Posterior of b")
for (i in 1:500) lines(c(-3,3), mu_p[i,] , col=col.alpha("blue", 0.4))

## M1.1: define multilevel model with adaptive priors -----
m1.1 <- ulam(
  alist(
    coop ~ dnorm(mu, sigma),
    mu <- a[k] + b*(rt_std),
    a[k] ~ dnorm(mean_a, sigma_a),
    b ~ dnorm(0, 10),
    sigma ~ dexp(1),
    mean_a ~ dnorm(50, 10),
    sigma_a ~ dexp(1)
  ),
  data = dat, chains = 4, iter = 4000, warmup = 1000, cores = 4, log_lik=T
)

### analyze posterior
precis(m1.1,2)
p <- extract.samples(m1.1)
dens(p$b, show.HPDI = .95, show.zero = T)
mu_p <- link(m1.1, post = p, data=list(rt_std=c(-3,3), k=c(1:13)))
plot( NULL, xlim=c(-3,3), ylim=c(0,100))
for (i in 1:500) lines(c(-3,3), mu_p[i,] , col=col.alpha("blue", 0.4))

## M1.2: define multilevel model with adaptive priors and SVO angle ----
m1.2 <- ulam(
  alist(
    coop ~ dnorm(mu, sigma),
    mu <- a[k] + b*(rt_std) + c*svo + d*svo*rt_std,
    a[k] ~ dnorm(mean_a, sigma_a),
    b ~ dnorm(0, 10),
    c ~ dnorm(0, 10),
    d ~ dnorm(0, 10),
    sigma ~ dexp(1),
    mean_a ~ dnorm(50, 10),
    sigma_a ~ dexp(1)
  ),
  data = dat, chains = 4, iter = 4000, warmup = 1000, cores = 4, log_lik = T
)

### analyze posteriors
precis(m1.2,2)
p <- extract.samples(m1.2)
HPDI(p$b, 0.95);dens(p$b, show.HPDI = .95)
HPDI(p$c, 0.95);dens(p$c, show.HPDI = .95)
HPDI(p$d, 0.95);dens(p$d, show.HPDI = .95)
### save b for future tests
beta_1_1 <- p$b
beta_3_1 <- p$d

### plot the relation for some values of SVO
svo_values <- c(-16.260,7.815,21.417,20.717,34.875,61.390)
rt_values <- seq(-1,1, length.out = 20)
coop_values_min <- link(m1.2, data= data.frame(k=1, svo=-16.26, rt_std = seq(-1,1, length.out = 20)))
coop_values_med <- link(m1.2, data= data.frame(k=1, svo=21.41, rt_std = seq(-1,1, length.out = 20)))
coop_values_max <- link(m1.2, data= data.frame(k=1, svo=61.3, rt_std = seq(-1,1, length.out = 20)))
coop_values_min <- apply(coop_values_min, 2, mean)
coop_values_med <- apply(coop_values_med, 2, mean)
coop_values_max <- apply(coop_values_max, 2, mean)

plot(rt_values, coop_values_min, type = "l", ylim=c(0,100), ylab="Cooperation", xlab="RT_std", col="red")
lines(rt_values, coop_values_med, type = "l", col="black")
lines(rt_values, coop_values_max, type = "l", col="blue")
legend(0.38, 95, legend=c("SVO= - 14", "SVO = 22", "SVO = 62"),
       col=c("red", "black","blue"), lty=1, cex=0.6,
       box.lty=0)

## NOT USED - DELETE --------------------
mR1 <- ulam(
  alist(
    rt_std ~ dnorm(mu, sigma),
    mu <- a[k] + b*svo,
    a[k] ~ dnorm(mean_a, sigma_a),
    b ~ dnorm(0, 10),
    sigma ~ dexp(1),
    mean_a ~ dnorm(0, 10),
    sigma_a ~ dexp(1)
  ),
  data = dat, chains = 4, iter = 4000, warmup = 1000, cores = 4, log_lik = T
)

mR1.1 <- ulam(
  alist(
    rt_std ~ dnorm(mu, sigma),
    mu <- a[k] + b*svo + c*svo*svo,
    a[k] ~ dnorm(mean_a, sigma_a),
    b ~ dnorm(0, 10),
    c ~ dnorm(0, 10),
    sigma ~ dexp(1),
    mean_a ~ dnorm(0, 10),
    sigma_a ~ dexp(1)
  ),
  data = dat, chains = 4, iter = 4000, warmup = 1000, cores = 4, log_lik = T
)
###########################
## MODEL CHECKS AND COMPARISON ----------
traceplot(m1)
traceplot(m1.1)
traceplot(m1.2)
dev.off()
compare(m1, m1.1, m1.2, WAIC = T); plot(compare(m1, m1.1, m1.2, WAIC = T))



# M4 - ANALYSIS OF EXTREMITY (Model 3 in the manuscript) --------------------
dens(dat_2$extremity_coop)
ggplot(data=as.data.frame(dat_2), aes(x=rt_std, y=extremity_coop))+geom_smooth()+geom_point()

m4 <- ulam(
  alist(
    extremity_coop ~ dnorm(mu, sigma),
    mu <- a[k] + b*rt_std,
    a[k] ~ dnorm(mean_a, sigma_a),
    b ~ dnorm(0, 10),
    sigma ~ dexp(1),
    mean_a ~ dnorm(0.3, 5),
    sigma_a ~ dexp(1)
  ),
  data = dat_2, chains = 4, iter = 4000, warmup = 1000, cores = 4, log_lik = T
)

precis(m4,2)
p <- extract.samples(m4)
dens(p$b, show.HPDI = .95);HPDI(p$b, .95)

# plot relation extremity coop RT
temp.list_m4 <- list()
for (i in 1:length(levels(factor(df$study_id)))) {
  
  temp.k <- link(m4, data= data.frame(k=i, rt_std = seq(-1,1, length.out = 15)))
  temp.list_m4[length(temp.list_m4)+1] <- list(temp.k)
}

# compute average across studies (k)
temp.list_2 <- list()
temp.list_3 <- list()
for (j in 1:length(levels(factor(df$study_id)))){
  print(j)
  temp.k <- apply(temp.list_m4[[j]], 2, mean)
  temp.p <- apply(temp.list_m4[[j]], 2, PI, prob=.95)
  temp.list_2[length(temp.list_2)+1] <- list(temp.k)
  temp.list_3[length(temp.list_3)+1] <- list(temp.p)
}
temp.m <- matrix(unlist(temp.list_2), ncol=15, byrow=T)
extremecoop.mean <- apply(temp.m, 2, mean)
extremecoop.CI <- apply(simplify2array(temp.list_3), 1:2, mean)

#plot (Figure 2 in the paper)
png(filename = paste(dir_figure, "extremecoop.png",sep=""), width = 597)
plot( seq(-1,1, length.out = 15),extremecoop.mean, type="l", col="blue",
      ylim=c(0,.6), xlab="RT (centered)", ylab="Extremity in Cooperation (E)")
shade(extremecoop.CI, seq(-1,1, length.out = 15) )
dev.off()



# MEDIATION ANALYSIS ------------------
# M4a --------------
m4.1 <- ulam(
  alist(
    coop ~ dnorm(mu, sigma),
    mu <- a[k] + b*(rt_std) + c*svo + d*svo*(rt_std) + e*extremity_coop + f*extremity_coop*svo,
    a[k] ~ dnorm(mean_a, sigma_a),
    b ~ dnorm(0, 10),
    c ~ dnorm(0, 10),
    d ~ dnorm(0, 10),
    e ~ dnorm(0, 10),
    f ~ dnorm(0, 10),
    sigma ~ dexp(1),
    mean_a ~ dnorm(50, 10),
    sigma_a ~ dexp(1)
  ),
  data = dat, chains = 4, iter = 4000, warmup = 1000, cores = 4, log_lik = T
)

precis(m4.1,2)
p <- extract.samples(m4.1)
### beta 1
b <- precis(as.data.frame(p$b), prob=.95)[1:4]; b$HPDI_low <- HPDI(p$b)[1]; b$HPDI_high <- HPDI(p$b)[2]
### beta 3
d <- precis(as.data.frame(p$d), prob=.95)[1:4]; d$HPDI_low <- HPDI(p$d)[1]; d$HPDI_high <- HPDI(p$d)[2]
### beta 4
e <- precis(as.data.frame(p$e), prob=.95)[1:4]; e$HPDI_low <- HPDI(p$e)[1]; e$HPDI_high <- HPDI(p$e)[2]
### beta 5
f <- precis(as.data.frame(p$f), prob=.95)[1:4]; f$HPDI_low <- HPDI(p$f)[1]; f$HPDI_high <- HPDI(p$f)[2]
table_b_m4.1 <- rbind.data.frame(b, d, e, f)
print(xtable(table_b_m4.1, type="latex"), file = paste(dir_table,"model4a.tex", sep=""))

##HP2: ---------
### abs(beta_14) < abs(beta_11)
HP2b <- abs(p$b) - abs(beta_1_1); precis(as.data.frame(HP2b)); HPDI(HP2b, .95); dens(HP2b, show.HPDI = .95)
HPDI(HP2b, 0.518)
### abs(beta_34) < abs(beta_31)
HP2d <- abs(p$d) - abs(beta_3_1); precis(as.data.frame(HP2d)); HPDI(HP2d, .95); dens(HP2d, show.HPDI = .95)
HPDI(HP2d, 0.848)

# model 4b -------------
m4.2 <- ulam(
  alist(
    coop ~ dnorm(mu, sigma),
    mu <- a[k] + b*(rt_std) + c*extremity_coop,
    a[k] ~ dnorm(mean_a, sigma_a),
    b ~ dnorm(0, 10),
    c ~ dnorm(0, 10),
    sigma ~ dexp(1),
    mean_a ~ dnorm(50, 10),
    sigma_a ~ dexp(1)
  ),
  data = dat_2_4b, chains = 4, iter = 4000, warmup = 1000, cores = 4, log_lik = T
)

precis(m4.2,2)
p <- extract.samples(m4.2)
b <- precis(as.data.frame(p$b), prob=.95)[1:4]; b$HPDI_low <- HPDI(p$b)[1]; b$HPDI_high <- HPDI(p$b)[2]
c <- precis(as.data.frame(p$c), prob=.95)[1:4]; c$HPDI_low <- HPDI(p$c)[1]; c$HPDI_high <- HPDI(p$c)[2]
table_b_m4.2 <- rbind.data.frame(b, c); print(table_b_m4.2)
print(xtable(table_b_m4.2, type="latex"), file = paste(dir_table,"model4b.tex", sep=""))

##HP2b: abs(beta_14b) < abs(beta_11 + beta13)---------
HP2b <- abs(p$b) - abs(beta_1_1+beta_3_1); precis(as.data.frame(HP2b)); HPDI(HP2b, .95); dens(HP2b, show.HPDI = .95)
HPDI(HP2b, .88)
# model 5 -----------------
m4.3 <- ulam(
  alist(
    extremity_svo ~ dnorm(mu, sigma),
    mu <- a[k] + b*(rt_std),
    a[k] ~ dnorm(mean_a, sigma_a),
    b ~ dnorm(0, 10),
    sigma ~ dexp(1),
    mean_a ~ dnorm(0, 10),
    sigma_a ~ dexp(1)
  ),
  data = dat, chains = 4, iter = 4000, warmup = 1000, cores = 4, log_lik = T
)

precis(m4.3,2)
p <- extract.samples(m4.3)
b <- precis(as.data.frame(p$b), prob=.95)[1:4]; b$HPDI_low <- HPDI(p$b)[1]; b$HPDI_high <- HPDI(p$b)[2]
table_b_m4.3 <- rbind.data.frame(b); print(table_b_m4.3)
print(xtable(table_b_m4.3, type="latex"), file = paste(dir_table,"model5.tex", sep=""))
# save beta_1_5
beta_1_5 <- p$b

# plot relation extremity in SVO RT
temp.list_m4.3 <- list()
for (i in 1:length(levels(factor(df$study_id)))) {
  
  temp.k <- link(m4.3, data= data.frame(k=i, rt_std = seq(-1,1, length.out = 15)))
  temp.list_m4.3[length(temp.list_m4.3)+1] <- list(temp.k)
}


# compute average across studies (k)
temp.list_2 <- list()
temp.list_3 <- list()
for (j in 1:length(levels(factor(df$study_id)))){
  print(j)
  temp.k <- apply(temp.list_m4.3[[j]], 2, mean)
  temp.p <- apply(temp.list_m4.3[[j]], 2, PI, prob=.95)
  temp.list_2[length(temp.list_2)+1] <- list(temp.k)
  temp.list_3[length(temp.list_3)+1] <- list(temp.p)
}
temp.m <- matrix(unlist(temp.list_2), ncol=15, byrow=T)
extremesvo.mean <- apply(temp.m, 2, mean)
extremesvo.CI <- apply(simplify2array(temp.list_3), 1:2, mean)

#plot
png(filename = paste(dir_figure, "extremesvo.png",sep=""), width = 597)
plot( seq(-1,1, length.out = 15),extremesvo.mean, type="l", col="blue",
      ylim=c(7,20), xlab="RT (centered)", ylab="Extremity in SVO (P)")
shade(extremesvo.CI, seq(-1,1, length.out = 15) )
dev.off()
# make a panel fo model 3 and 5
png(filename = paste(dir_figure, "panel.png",sep=""), width = 597)
par(mfrow=c(1,2)) # make panel figure
# model 3
plot( seq(-1,1, length.out = 15),extremecoop.mean, type="l", col="blue",
      ylim=c(0,.6), xlab="RT (centered)", ylab="Extremity in Cooperation (E)")
shade(extremecoop.CI, seq(-1,1, length.out = 15) )
# model 5
plot( seq(-1,1, length.out = 15),extremesvo.mean, type="l", col="blue",
      ylim=c(7,20), xlab="RT (centered)", ylab="Extremity in SVO (P)")
shade(extremesvo.CI, seq(-1,1, length.out = 15) )

dev.off()
# model 6a --------------
m6a <- ulam(
  alist(
    coop ~ dnorm(mu, sigma),
    mu <- a[k] + b*(rt_std) + c*extremity_svo,
    a[k] ~ dnorm(mean_a, sigma_a),
    b ~ dnorm(0, 10),
    c ~ dnorm(0, 10),
    sigma ~ dexp(1),
    mean_a ~ dnorm(50, 10),
    sigma_a ~ dexp(1)
  ),
  data = dat_6a, chains = 4, iter = 4000, warmup = 1000, cores = 4, log_lik = T
)

precis(m6a,2)
p <- extract.samples(m6a)
b <- precis(as.data.frame(p$b), prob=.95)[1:4]; b$HPDI_low <- HPDI(p$b)[1]; b$HPDI_high <- HPDI(p$b)[2]
c <- precis(as.data.frame(p$c), prob=.95)[1:4]; c$HPDI_low <- HPDI(p$c)[1]; c$HPDI_high <- HPDI(p$c)[2]
table_b_m6a <- rbind.data.frame(b, c)
print(xtable(table_b_m6a, type="latex"), file = paste(dir_table,"model6a.tex", sep=""))

##HP3a: abs(beta_14b) < abs(beta_11 + beta13)------
HP3a <- abs(p$b) - abs(beta_1_1); precis(as.data.frame(HP3a)); HPDI(HP3a, .95); dens(HP3a, show.HPDI = .95)
HPDI(HP3a, .29)
# model 6b -------------
m6b <- ulam(
  alist(
    coop ~ dnorm(mu, sigma),
    mu <- a[k] + b*(rt_std) + c*extremity_svo,
    a[k] ~ dnorm(mean_a, sigma_a),
    b ~ dnorm(0, 10),
    c ~ dnorm(0, 10),
    sigma ~ dexp(1),
    mean_a ~ dnorm(50, 10),
    sigma_a ~ dexp(1)
  ),
  data = dat_6b, chains = 4, iter = 4000, warmup = 1000, cores = 4, log_lik = T
)

precis(m6b,2)
p <- extract.samples(m6b)
b <- precis(as.data.frame(p$b), prob=.95)[1:4]; b$HPDI_low <- HPDI(p$b)[1]; b$HPDI_high <- HPDI(p$b)[2]
c <- precis(as.data.frame(p$c), prob=.95)[1:4]; c$HPDI_low <- HPDI(p$c)[1]; c$HPDI_high <- HPDI(p$c)[2]
table_b_m6b <- rbind.data.frame(b, c); print(table_b_m6b)
print(xtable(table_b_m6b, type="latex"), file = paste(dir_table,"model6b.tex", sep=""))

##HP3b: abs(beta_14b) < abs(beta_11 + beta13)---------
HP3b <- abs(p$b) - abs(beta_1_1+beta_3_1); precis(as.data.frame(HP3b)); HPDI(HP3b, .95); dens(HP3b, show.HPDI = .95)
