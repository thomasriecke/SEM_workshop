# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Welcome to latent variables. This type of modelling is one of the most
# intuitive things I've ever done conceptually. There are also statistical components of
# this process that are incredibly non-intuitive :)
#
# Our goal this week is to familiarize ourselves with the concept (i.e.,
# we're observing multiple measurements of an underlying 'latent' process), and 
# to begin to familiarize ourselves with implementing these types of models
# in lavaan/blavaan and JAGS. We will take our time. This is a critical concept
# for understanding SEMs. Ask questions as they arise!
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# our first simulated dataset will be counts of warblers (y; Setophoga townsendi)
# in different age-classes of forests. Counts will be greater in older forests
# 
# We will (initially) measure two components of forest age, canopy height (c) and 
# sub-canopy height (s).
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load packages
library(jagsUI)
library(vioplot)
library(lavaan)
library(blavaan)

set.seed(1234)
n <- 200
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1), family = 'sans')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) Let's simulate forest age (i.e., maturity) as a z-standardized covariate
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
maturity <- rnorm(n, 0, 1)
hist(maturity, breaks = 50, main = NULL, xlab = 'Forest maturity', cex.lab =2)
# large positive values indicate 'old growth', very negative numbers 
# are doghair thickets (dense, short, regen following logging or stand-replacing
# fires)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2) let's simulate canopy height and sub-canopy height
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
canopy <- rlnorm(n, 3.65 + 0.25 * maturity, 0.05)
hist(canopy, breaks = 50, main = NULL, xlab = 'Canopy height (m)')
plot(canopy ~ maturity, ylab = 'Canopy height (m)',
     xlab = 'Forest maturity', cex.lab = 2)


# similarly, we'll simulate sub-canopy height
subcan <- rlnorm(n, 2 + 0.5 * maturity, 0.1) 
hist(subcan, breaks = 50, main = NULL, xlab = 'Sub-canopy height (s)')
plot(subcan ~ maturity, ylab = 'Sub-canopy height (m)',
     xlab = 'Forest maturity', cex.lab = 2)

# Note the extreme multicollinearity between the two covariates... this
# is because they're arising from the same underlying 'latent' process.
# As forest stand mature, trees get taller, similarly, as stands mature and 
# gaps develop, there is more room for regeneration and second-growth
plot(canopy ~ subcan, ylab = 'Canopy height (m)',
     xlab = 'Subcanopy height (m)', cex.lab = 2, las = 1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3) let's simulate warblers, more warblers in older forests
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
warblers <- rpois(n, exp(0.5 + 0.75 * maturity))

par(mfrow= c(1,1))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4) Let's visualize our data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,3), mar = c(5.1,5.1,2.1,2.1), family = 'sans')
plot(jitter(warblers) ~ canopy,
     ylab = 'Warbler abundance', xlab = 'Canopy height (m)', cex.lab = 2, las = 1)
plot(jitter(warblers) ~ subcan,
     ylab = 'Warbler abundance', xlab = 'Subcanopy height (m)', cex.lab = 2, las = 1)
plot(canopy ~ subcan, ylab = 'Canopy height (m)',
     xlab = 'Subcanopy height (m)', cex.lab = 2, las = 1)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# teaching figures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d <- data.frame(c = canopy, s = subcan, y = warblers)
sem1 <- sem('m =~ c + s
            y ~ m
            c ~ 1
            s ~ 1',
            data = d)
summary(sem1)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2x2 plot of all covariates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(2,2), mar = c(5.1,5.1,2.1,5.1), family = 'sans')
hist(maturity, breaks = 50, main = NULL, xlab = 'Forest maturity', cex.lab =2)
plot(subcan ~ maturity, ylab = 'Sub-canopy height (m)',las = 1,
     xlab = 'Forest maturity', cex.lab = 2)
axis(side = 4, labels = round(log(c(1,5,10,20)), digits = 1), at = c(1,5,10,20), las = 1)
mtext('ln(Sub-canopy height)', side = 4, line = 3, cex = 1.5)

plot(canopy ~ maturity, ylab = 'Canopy height (m)',las = 1,
     xlab = 'Forest maturity', cex.lab = 2)
axis(side = 4, labels = round(log(c(10,20,40,80)), digits = 1), at = c(10,20,40,80), las = 1)
mtext('ln(Canopy height)', side = 4, line = 3, cex = 1.5)

plot(c(warblers) ~ maturity, ylab = 'Warbler counts',las = 1,
     xlab = 'Forest maturity', cex.lab = 2)
axis(side = 4, labels = round(log(c(1,3,6,12)), digits = 1), at = c(1,3,6,12), las = 1)
mtext('ln(Warbler counts)', side = 4, line = 3, cex = 1.5)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# rough changes in y
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
max(canopy) - min(canopy)
max(subcan) - min(subcan)
max(warblers) - min(warblers)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2x2 plot of all covariates with no x-axis
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(2,2), mar = c(5.1,5.1,2.1,5.1), family = 'sans')
hist(maturity, breaks = 50, main = NULL, xlab = 'Forest maturity', cex.lab =2, xaxt = 'n')
plot(subcan ~ maturity, ylab = 'Sub-canopy height (m)',las = 1,
     xlab = 'Forest maturity', cex.lab = 2, xaxt = 'n')
axis(side = 4, labels = round(log(c(1,5,10,20)), digits = 1), at = c(1,5,10,20), las = 1)
mtext('ln(Sub-canopy height)', side = 4, line = 3, cex = 1.5)

plot(canopy ~ maturity, ylab = 'Canopy height (m)',las = 1,
     xlab = 'Forest maturity', cex.lab = 2, xaxt = 'n')
axis(side = 4, labels = round(log(c(10,20,40,80)), digits = 1), at = c(10,20,40,80), las = 1)
mtext('ln(Canopy height)', side = 4, line = 3, cex = 1.5)

plot(c(warblers) ~ maturity, ylab = 'Warbler counts',las = 1,
     xlab = 'Forest maturity', cex.lab = 2, xaxt = 'n')
axis(side = 4, labels = round(log(c(1,3,6,12)), digits = 1), at = c(1,3,6,12), las = 1)
mtext('ln(Warbler counts)', side = 4, line = 3, cex = 1.5)



d <- data.frame(c = canopy, s = subcan, y = warblers)
sem1 <- sem('m =~ c + s
            y ~ m
            c ~ 1
            s ~ 1',
            data = d)
summary(sem1)
m.pred <- as.numeric(predict(sem1))
sd(m.pred)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2x2 plot of all covariates with new x-axis
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(2,2), mar = c(5.1,5.1,2.1,5.1), family = 'sans')
hist(m.pred, breaks = 50, main = NULL, xlab = 'Forest maturity', cex.lab =2)
plot(subcan ~ m.pred, ylab = 'Sub-canopy height (m)',las = 1,
     xlab = 'Forest maturity', cex.lab = 2)
axis(side = 4, labels = round(log(c(1,5,10,20)), digits = 1), at = c(1,5,10,20), las = 1)
mtext('ln(Sub-canopy height)', side = 4, line = 3, cex = 1.5)

plot(canopy ~ m.pred, ylab = 'Canopy height (m)',las = 1,
     xlab = 'Forest maturity', cex.lab = 2)
axis(side = 4, labels = round(log(c(10,20,40,80)), digits = 1), at = c(10,20,40,80), las = 1)
mtext('ln(Canopy height)', side = 4, line = 3, cex = 1.5)

plot(c(warblers) ~ m.pred, ylab = 'Warbler counts',las = 1,
     xlab = 'Forest maturity', cex.lab = 2)
axis(side = 4, labels = round(log(c(1,3,6,12)), digits = 1), at = c(1,3,6,12), las = 1)
mtext('ln(Warbler counts)', side = 4, line = 3, cex = 1.5)





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5) models with different fixed slopes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d <- data.frame(c = canopy, s = subcan, y = warblers)
sem1 <- sem('m =~ c + s
            y ~ m
            c ~ 1
            s ~ 1',
            data = d)
summary(sem1)
m.pred1 <- as.numeric(predict(sem1))


sem2 <- sem('m =~ s + c
            y ~ m
            c ~ 1
            s ~ 1',
            data = d)
summary(sem2)
m.pred2 <- as.numeric(predict(sem2))


par(mfrow = c(1,1), mar = c(5.1,5.1,5.1,5.1))
plot(m.pred1 ~ m.pred2, 
     ylab = "Latent (m) if canopy beta = 1",
     xlab = "Latent (m) if subcanopy beta = 1", cex.lab = 1.5, las = 1)
axis(side = 3, labels = c(5,15,25), at = c(5,15,25) - 8.164)
mtext('Sub-canopy height (m)', side = 3, line = 2.5, cex = 1.5)

axis(side = 4, labels = c(25,52.5,80), at = c(25,52.5,80) - 39.451)
mtext('Canopy height (m)', side = 4, line = 2.5, cex = 1.5)

py1 <- 2.135 + m.pred1 * 0.196
py2 <- 2.135 + m.pred2 * 0.431

plot(py1 ~ py2, las = 1, cex.lab = 1.5,
     ylab = 'Predicted warblers if canopy beta = 1',
     xlab = 'Predicted warblers if subcanopy beta = 1')
abline(0,1)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5) blavaan model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getwd()
d <- data.frame(c = canopy, s = subcan, y = warblers)
sem1 <- bsem('m =~ c + s
            y ~ m
            c ~ 1
            s ~ 1',
            data = d,
            target = 'jags',
            mcmcfile = T)
summary(sem1)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# look in your working directory
# there will be a perfectly functioning JAGS model file
# that was used to generate the parameter estimates in a lavExport folder
# you may have to change some parameter names &| priors, but that's it!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~