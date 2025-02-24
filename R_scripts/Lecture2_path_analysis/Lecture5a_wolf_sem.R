# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# packages
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(piecewiseSEM)
library(latex2exp)
library(semEff)

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

model <- "ung =~ z.deer + z.elk
          ung ~ z.elev
          used ~ ung + z.elev
          z.deer ~ 0
          z.elk ~ 0"

sem1 <- sem(model, data =dat, ordered = 'used')
summary(sem1)


plot(dat$z.deer ~ dat$z.elev)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bow Valley wolves model 2
# structural equation model:
# latent low-elevation ungulate abundance is a function of both deer and elk abundance
# latent hi-elevation ungulate abundance is a function of both deer and elk abundance
# wolf rsf is a function of ungulate abundance
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model <- "lo.ung =~ z.deer + z.elk
          hi.ung =~ z.sheep + z.goat
          lo.ung ~ z.elev
          hi.ung ~ z.elev
          used ~ hi.ung + lo.ung
          z.deer ~ 0
          z.elk ~ 0
          z.sheep ~ 0
          z.goat ~ 0"

sem2 <- sem(model, data =dat, ordered = 'used')
summary(sem2)


