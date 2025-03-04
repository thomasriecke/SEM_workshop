# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# packages
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(piecewiseSEM)
library(latex2exp)
library(semEff)
library(lavaan)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plotting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1), family = 'sans')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read data from GitHub
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dat <- read.csv("C:/Users/thomas.riecke/Desktop/SEM_workshop/data/wolf_rsf.csv")
x <- url("https://github.com/thomasriecke/SEM_workshop/blob/main/data/wolf_rsf.csv?raw=true")
dat <- read.csv(x)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# subset to Bow Valley pack
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
str(dat)
dat <- subset(dat, pack == 'Bow valley')
dat <- subset(dat, !is.na(distacc))
table(dat$used)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make a few quick maps
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ggplot2)
# map of elevation
ggplot(dat, aes(x = easting, y = northing, col = elevati)) + 
  geom_point(size = 2) + 
  coord_equal() +
  scale_colour_gradient(low = 'yellow', high = 'red')

# map of deer winter use
ggplot(dat, aes(easting, northing, col = deerwin)) + 
  geom_point(size = 2) + 
  coord_equal() +
  scale_colour_gradient(low = 'yellow', high = 'red')
# what's the relationships between deer winter use and elevation?
# how could we visualize or model that?


# map of used and random points
ggplot(dat, aes(easting, northing)) + 
  geom_point(size = dat$used*2+0.01, shape = dat$used) + 
  coord_equal()



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# scale covariates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat$z.elev <- as.numeric(scale(dat$elevati))
dat$z.deer <- as.numeric(scale(dat$deerwin))
dat$z.elk <- as.numeric(scale(dat$elkwin))
dat$z.moose <- as.numeric(scale(dat$moosewin))
dat$z.sheep <- as.numeric(scale(dat$sheepwin))
dat$z.goat <- as.numeric(scale(dat$goatwin))
dat$z.dist <- as.numeric(scale(dat$distacc))


cor(data.frame(dat$deerwin,dat$elkwin,dat$moosewin,dat$goatwin,dat$sheepwin))





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bow Valley wolves model 1
# structural equation model:
# latent low-elevation ungulate abundance is a function of both deer and elk abundance
# latent hi-elevation ungulate abundance is a function of both deer and elk abundance
# wolf rsf is a function of ungulate abundance
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plot(jitter(dat$z.deer) ~ jitter(dat$z.elk), las = 1, cex.lab = 1.25,
     ylab = 'Deer winter use (z-score)',
     xlab = 'Elk winter use (z-score)')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# note that we fix intercepts to 0 because we scaled 
# elk, moose, and wolf use, i.e.:
# z = (x - mean(x))/sd(x)
# this means that the intercepts are 0
#
# we also derive the total effect of elevation acting via deer
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model <- "ung =~ z.elk + z.deer + z.moose
          ung ~ b * z.elev
          used ~ a * ung + c * z.elev
          total := a * b + c
          z.deer ~ 0
          z.elk ~ 0
          z.moose ~ 0"


sem1 <- sem(model, data = dat, ordered = 'used')
out1 <- summary(sem1)
out1

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# note that lavaan is using ordinal regression! note a logit link
# (slightly different than slides)
# this is really interesting, but a bit complicated.
# https://en.wikipedia.org/wiki/Ordinal_regression#:~:text=In%20statistics%2C%20ordinal%20regression%2C%20also,between%20different%20values%20is%20significant.
# Fortunately here we're only considering two options (unused and used)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's plot the total effect of elevation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
str(out1)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we extract parameter values and generate uncertainty
# using Monte Carlo sampling
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
samps <- 1000
cut <- rnorm(samps, out1$pe$est[10], out1$pe$se[10]) # our cut-off point for ordinal (0 -> 1)
total <- rnorm(samps, out1$pe$est[18], out1$pe$se[18])
var <- out1$pe$est[14]
hist(total, breaks = 100) # Monte Carlo sampling our total effect

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# here we create matrices to store results and a range of covariate values
# res: the resolution for our line (number of dots to connect)
# pred.z: a range of possible covariate values
# pred.y: estimates on the real scale
# pred.p: estimates on the probability scale
# q.p: quantiles of estimates on the probability scale
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
res <- 100
pred.z <- seq(min(dat$z.elev), max(dat$z.elev), length.out = res)
pred.y <- matrix(NA, samps, res) # predicted on real scale
pred.p <- matrix(NA, samps, res) # predicted probability
q.p <- matrix(NA, res, 5)        # quantiles of predicted probability

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run a loop to generate estimates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (j in 1:res){
  pred.y[,j] <- pred.z[j] * total
  pred.p[,j] <- pnorm(cut, pred.y[,j], sqrt(var), lower.tail = F)
  q.p[j,] <- quantile(pred.p[,j], c(0.025,0.05,0.5,0.95,0.975))
}

plot(q.p[,3] ~ pred.z, type = 'l', ylab = 'P(Wolf Used)', las = 1, cex.lab = 1.5,
     xlab = 'Z-standardized elevation', ylim = c(0,0.75))
lines(q.p[,1] ~ pred.z, lty = 2)
lines(q.p[,5] ~ pred.z, lty = 2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make the same plot with real values of elevation (meters)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
real.elev <- seq(1400,2600,length.out = 6)
z.loc <- (real.elev - mean(dat$elevati))/sd(dat$elevati)

plot(q.p[,3] ~ pred.z, type = 'l', ylab = 'P(Wolf Used)', las = 1, cex.lab = 1.5,
     xlab = 'Z-standardized elevation', ylim = c(0,0.75), xaxt = 'n')
axis(side = 1, at = z.loc, labels = real.elev)
lines(q.p[,1] ~ pred.z, lty = 2)
lines(q.p[,5] ~ pred.z, lty = 2)






# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Homework
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) remove moose from the model. what changes and why?
# 2) try to make the same plot as above but just with the effect
#    of ungulates
# 3) Try to add sheep and goats to the model. Examine the latent variable
#    loadings. What is happening?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bow Valley wolves model 2
# structural equation model:
# latent low-elevation ungulate abundance is a function of both deer and elk abundance
# latent hi-elevation ungulate abundance is a function of both sheep and goat abundance
# wolf rsf is a function of ungulate abundance, do wolves select for lo-
# or hi- elevation ungulates? pretend this is another system where this
# is actually an interesting question :) 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model <- "lo.ung =~ z.elk + z.deer
          hi.ung =~ z.sheep + z.goat
          lo.ung ~ z.elev
          hi.ung ~ z.elev
          used ~ hi.ung + lo.ung
          z.deer ~ 0
          z.elk ~ 0
          z.sheep ~ 0
          z.goat ~ 0"

sem2 <- sem(model, data =dat, ordered = 'used')
out2 <- summary(sem2)


out2$pe

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Let's plot the effect on use given lo- and hi-elevation 
# ungualte abundance
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
latent.sem2 <- data.frame(predict(sem2))
plot(latent.sem2$lo.ung ~ latent.sem2$hi.ung)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we extract parameter values and generate uncertainty
# using Monte Carlo sampling
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
samps <- 1000
cut <- rnorm(samps, out2$pe$est[13], out2$pe$se[13]) # our cut-off point for ordinal (0 -> 1)
hi <- rnorm(samps, out2$pe$est[7], out2$pe$se[7])
lo <- rnorm(samps, out2$pe$est[8], out2$pe$se[8])
var <- out2$pe$est[18]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# here we create matrices to store results and a range of covariate values
# res: the resolution for our line (number of dots to connect)
# pred.z: a range of possible covariate values
# pred.y: estimates on the real scale
# pred.p: estimates on the probability scale
# q.p: quantiles of estimates on the probability scale
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
res <- 100
px.lo <- seq(min(latent.sem2$lo.ung), max(latent.sem2$lo.ung), length.out = res)
px.hi <- seq(min(latent.sem2$hi.ung), max(latent.sem2$hi.ung), length.out = res)

pred.lo <- matrix(NA, samps, res) # predicted on real scale
pred.hi <- matrix(NA, samps, res) # predicted on real scale
p.lo <- matrix(NA, samps, res) # predicted probability
p.hi <- matrix(NA, samps, res) # predicted probability
q.lo <- matrix(NA, res, 5)        # quantiles of predicted probability
q.hi <- matrix(NA, res, 5)        # quantiles of predicted probability
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run a loop to generate estimates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (j in 1:res){
  pred.lo[,j] <- px.lo[j] * lo
  pred.hi[,j] <- px.hi[j] * hi  
  p.lo[,j] <- pnorm(cut, pred.lo[,j], sqrt(var), lower.tail = F)
  p.hi[,j] <- pnorm(cut, pred.hi[,j], sqrt(var), lower.tail = F)  
  q.lo[j,] <- quantile(p.lo[,j], c(0.025,0.05,0.5,0.95,0.975))
  q.hi[j,] <- quantile(p.hi[,j], c(0.025,0.05,0.5,0.95,0.975))  
}

plot(q.lo[,3] ~ px.lo, type = 'l', ylab = 'P(Wolf Used)', las = 1, cex.lab = 1.5,
     xlab = 'Latent deer and elk abundance', ylim = c(0,1))
lines(q.lo[,1] ~ px.lo, lty = 2)
lines(q.lo[,5] ~ px.lo, lty = 2)

plot(q.hi[,3] ~ pred.z, type = 'l', ylab = 'P(Wolf Used)', las = 1, cex.lab = 1.5,
     xlab = 'Latent sheep and goat abundance', ylim = c(0,1))
lines(q.hi[,1] ~ pred.z, lty = 2)
lines(q.hi[,5] ~ pred.z, lty = 2)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Homework
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) use the model output to figure out the linkage between
#    elevation and sheep and goat abundance... make a new axis for elevation
#    (i.e., the hypothesized driver of hi-elev ungulate use)
# 2) what happens to inference if we use more limited model structures
#    (i.e., just deer and sheep?)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


