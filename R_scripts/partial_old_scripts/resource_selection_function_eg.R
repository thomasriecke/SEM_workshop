# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script will work use through our first 'path analysis,' a term
# often used to describe structural equation models without latent variables
#
# We'll use resource-selection data (used vs. available) from the Bow Valley
# wolf pack in Banff National Park and the R package piecewiseSEM
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# packages
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(piecewiseSEM)
library(latex2exp)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plotting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1), family = 'sans')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read and format data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- read.csv("C:/Users/thomas.riecke/Desktop/SEM_workshop/data/wolf_rsf.csv")
str(dat)
dat <- subset(dat, pack == 'Bow valley')
dat <- subset(dat, !is.na(distacc))
table(dat$used)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's run two models... one where we estimate wolf use as a function of 
# deer use and elevation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# glm
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
glm1 <- glm(used ~ deerwin + elevati, data = dat, family = 'binomial')
summary(glm1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot prediction of elevation effect
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
res <- 100 # resolution of line
px <- seq(min(dat$elevati), max(dat$elevati), length.out = res)
py <- predict(glm1, 
              newdata = data.frame(elevati = px, deerwin = rep(mean(dat$deerwin), res)), 
              se.fit = T)
lci <- with(py, plogis(fit + qnorm(0.025)*se.fit))
uci <- with(py, plogis(fit + qnorm(0.975)*se.fit))

plot(plogis(py$fit) ~ px, type = 'l', xlab = 'Elevation',
     ylab = '')

?predict()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# and another where we do the same thing, and also acknowledge 
# (technically, we assume) that elevation has an effect on deer use
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
path1 <- psem(
  glm(used ~ deerwin + elevati, data = dat, family = 'binomial'),
  glm(deerwin ~ elevati, data = dat),
  data = dat
)
summary(path1, intercepts = T)












# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# a bit of foreshadowing, it seems that ungulate use has interspecific correlations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(2,2), mar = c(5.1,5.1,2.21,2.1), family = 'sans')
plot(jitter(elkwin) ~ jitter(deerwin), data = dat, las = 1,
     xlab = 'Deer winter use', ylab = 'Elk winter use')
plot(jitter(moosewin) ~ jitter(deerwin), data = dat, las = 1,
     xlab = 'Deer winter use', ylab = 'Moose winter use')
plot(jitter(moosewin) ~ jitter(elkwin), data = dat, las = 1,
     ylab = 'Moose winter use', xlab = 'Elk winter use')

cor.test(dat$elkwin,dat$deerwin)
cor.test(dat$moosewin,dat$deerwin)
cor.test(dat$moosewin,dat$elkwin)

ungulate <- data.frame(dat$elkwin,dat$deerwin,dat$moosewin)
cov(ungulate)
cor(ungulate)
