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

