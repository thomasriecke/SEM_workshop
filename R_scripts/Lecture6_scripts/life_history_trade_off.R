# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Let's imagine 1000 female ungulates we follow for up to 10 years (quite the study!).
# During each year we'll measure their reproductive output (rho; expected number of fawns)
# and monitor their survival as a function of their
# survival probability (phi). This will result in alive/dead data (z)
# and reproductive data (y).
#
# The 'catch' or 'trick' is that these individuals will be of different
# 'quality.' Some individuals will have high survival and reproductive potential
# while others are more likely to die and less likely to reproduce.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(latex2exp)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# number of individuals (nI) and number of years (nT)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nI <- 1000
nT <- 10

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# data matrices
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
z <- matrix(NA, nI, nT) # alive/dead (all alive in first year)
z[,1] <- 1
y <- matrix(NA, nI, nT) # number of fawns born

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# probability matrices
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
phi <- matrix(NA, nI, nT-1) # survival
rho <- matrix(NA, nI, nT) # recruitment

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate heterogeneity
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
eta <- rnorm(nI, 0, 1)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# cross-lag simulation
# E(phi[t]) is a function of rho[t]
# E(rho[t]) is a function of eta and rho[t-1]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha <- c(2, -0.25, 0.5)
beta <- c(1, -1)
q <- matrix(NA, nI, 3)
p <- matrix(NA, nI, 3)
c <- c(-2,2)





for (i in 1:nI){
  rho[i,1] = beta[1] * eta[i]
  
  # ordinal regression equations [3 categories, no reproduction, 1 fawn, 2 fawns]
  q[i,1] = plogis(c[1] - rho[i,1])
  p[i,1] = q[i,1]
  q[i,2] = plogis(c[2] - rho[i,1])
  p[i,2] = q[i,2] - q[i,1]
  p[i,3] = 1 - q[i,2]
  
  y[i,1] = which(rmultinom(1, 1, p[i,]) == 1) - 1
  
  for (t in 2:nT){
    
    if (z[i,t-1] == 1){
      phi[i,t-1] = plogis(alpha[1] + alpha[2] * y[i,t-1] + alpha[3] * eta[i])
      z[i,t] <- rbinom(1, z[i,t-1], phi[i,t-1])
      
        # reproduction value
        rho[i,t] = beta[1] * eta[i] + beta[2] * y[i,t-1]
        
        # ordinal regression equations [3 categories, no reproduction, 1 fawn, 2 fawns]
        q[i,1] = plogis(c[1] - rho[i,t])
        p[i,1] = q[i,1]
        q[i,2] = plogis(c[2] - rho[i,t])
        p[i,2] = q[i,2] - q[i,1]
        p[i,3] = 1 - q[i,2]

        y[i,t] <- ifelse(z[i,t] == 1, 
                         which(rmultinom(1, z[i,t], p[i,]) == 1) - 1, 
                         NA)
    }

    if (z[i,t-1] == 0){
      z[i,t] <- 0
      y[i,t] <- NA
    }

    
  }
  
}

par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1), family = 'sans')
boxplot(phi[,1] ~ y[,1], xlab = 'Number of fawns', boxwex = 0.25, outline = F,
        ylab = 'Survival probability', cex.lab = 1.5, las = 1)

plot(jitter(y[,2], 0.5) ~ jitter(y[,1], 0.5), xlab = TeX("$y_{t-1}$"),
        ylab = TeX("$y_{t}$"), cex.lab = 1.5, las = 1)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# We can plot the breakpoints as if they're integers (the lower limit)
# to retain a 'single line' type of approach rather than trying to deal 
# with probabilities of each outcome? Dealing with clutch size would just
# add breakpoints?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p.rho <- beta[2] * 0:2 # calculate expected number of fawns as a function of previous success


plot(p.rho ~ c(0:2), type = 'l', lwd = 2, yaxt = 'n', ylim = c(-2.5,2),
     xaxt = 'n', ylab = TeX("Expected number of fawns in t ($\\rho_{it}$)"),
     xlab = TeX("Number of fawns in t-1 ($y_{t-1}$)"), cex.lab = 1.25)
axis(side = 1, at = c(0:2), labels = 0:2)
axis(side = 2, at = c, labels = 1:2) # use breakpoints as integers


