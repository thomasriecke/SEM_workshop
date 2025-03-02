# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# a dynamic SEM with environmental variation (x), variation in carrying 
# capactiy as a function of vegetation(K), and variation in abundance of an herbivore (n)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(jagsUI)
library(latex2exp)
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1), mfrow = c(1,1))


# length of time-series (nT)
nT <- 50

# number of simulations (nS)
nS <- 1000

beta <- rep(NA, nS)
pval <- rep(NA, nS)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# run the simulation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (ii in 1:nS){
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # variation in environmental conditions; 0 is average
  # note that here we use an auto-regressive approach to simulate variation in x
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  x <- rep(NA, nT)
  x[1] <- 0
  sigma <- 0.5
  rho <- 0.5
  
  for (t in 2:nT){
    x[t] <- rnorm(1, x[t-1] * rho, sigma)
  }
  
  # plot(x, ylab = 'Environmental conditions (x)', xlab = 'Time', type = 'b',
  #      las = 1, cex.lab = 1.5, cex.axis = 1.25, cex = 1.5)
  
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # now we simulate starting values for abundance and carrying capacity
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  n <- rep(NA, nT)
  k <- rep(NA, nT)
  
  k[1] <- 1000
  n[1] <- 1000
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # and simulate the population and carrying capacity forward through time
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  r <- 0.3
  delta <- c(0.05, 0.025)
  for (t in 2:nT){
    n[t] <- rpois(1, n[t-1] * exp(r * (1 - (n[t-1]/k[t-1]))))
    k[t] <- rpois(1, k[t-1] * exp(delta[1] * (1 - (n[t-1]/k[t-1])) + delta[2] * x[t]))
  }
  
  # derive population growth rate (lambda or 'l')
  l <- n[2:nT]/n[1:(nT-1)]
  
  
  # plot(k, type = 'b', cex.lab = 2, ylab = 'Abundance (n) and carrying capacity (k)', xlab = 'Time',
  #     ylim = c(0,max(n,k)), cex = 1.5)
  # points(n, pch = 19, type = 'b', cex = 1.5)
  # legend('bottomleft', legend = c('Abundance (n)','Carrying capacity (k)'), cex = 1.5,
  #        pch = c(19, 21), lty = c(1,1), pt.bg = c(NA, 'white'), bty = 'n')
  
  
  # plot(l ~ x[1:(nT-1)], cex = 1.5, cex.axis = 1.25, cex.lab = 1.5,
  #      ylab = TeX("Population growth rate ($\\lambda$)"), xlab = 'Environmental conditions (x)')
  # cor.test(l, x[1:(nT-1)])
  
  m <- lm(l ~ x[1:(nT-1)])
  
  # str(m)
  # str(summary(m))
  beta[ii] <- as.numeric(m$coefficients[2])
  pval[ii] <- as.numeric(summary(m)$coefficients[2,4])
}

par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1), mfrow = c(1,1))
boxplot(beta, outline = F, boxwex = 0.5,
        ylab = TeX("Regression coefficient ($\\beta$)"), cex.lab = 1.5,
        las = 1, cex.axix = 1.25)
abline(h = 0.0, lty = 1, col = 'grey', lwd = 2)
points(beta ~ jitter(rep(1, nS)))

boxplot(pval, outline = F, boxwex = 0.5,
        ylab = TeX("p-value"), cex.lab = 1.5,
        las = 1, cex.axix = 1.25)
abline(h = 0.05, lty = 1, col = 'grey', lwd = 2)
points(pval ~ jitter(rep(1, nS)))

length(which(beta > 0))/nS
length(which(pval > 0.05))/nS




sink("cross_lag.jags")
cat("
    model {

    delta[1] ~ dnorm(0,1)
    delta[2] ~ dnorm(0,1)
    
    r ~ dgamma(1,1)
    rho ~ dbeta(1,1)
    
    sigma ~ dgamma(1,1)
    tau = 1/pow(sigma,2)
    k[1] ~ dpois(1000)

    for (t in 2:nT){
      x[t] ~ dnorm(x[t-1] * rho, tau)
      n[t] ~ dpois(n[t-1] * exp(r * (1 - (n[t-1]/k[t-1]))))
      k[t] ~ dpois(k[t-1] * exp(delta[1] * (1 - (n[t-1]/k[t-1])) + delta[2] * x[t-1]))      
    }


    }
    ",fill = TRUE)
sink()



jags.data <- list(nT = nT, n = n, x = x)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we provide initial values. This is where JAGS will begin sampling
# to build a posterior distribution, in this case, I won't provide any 
# initial values
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inits <- function(){list()}  

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This is where we tell JAGS which parameters we wish to monitor
# i.e., which posterior distributions we want to save, 
# generally we'll save all of them, but 
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
parameters <- c('delta','rho','sigma','k')

# number of chains (nc), thinning rate (nt), number of iterations (ni), and number to burn-in
nc <- 4
nt <- 10
ni <- 20000
nb <- 10000


# Call JAGS from R 
# 20k iterations takes awhile on an i9
library(jagsUI)
Sys.time()
m <- jags(jags.data, inits, parameters, "cross_lag.jags", parallel = T, 
          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
Sys.time()

print(m)



plot(m$q50$k ~ k,
     ylab = 'Model estimates of carrying capacity (k)',
     xlab = 'True carrying capacity')
abline(0,1)

vioplot::vioplot(m$sims.list$k, outline = F, drawRect = F)
mtext('Estimated carrying capacity', side = 2, line = 2.5, cex = 1.5)
arrows(1:nT, m$q2.5$k, 1:nT, m$q97.5$k, length = 0, col = 'white')
points(m$q50$k, pch = 21, bg = 'white', cex = 1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Note, we supply zero information or data re: carrying capacity other
# than a prior for t = 1, abundance data (n), and environmental data (x)
#
# having some type of data would tighten CrI's substantially
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
