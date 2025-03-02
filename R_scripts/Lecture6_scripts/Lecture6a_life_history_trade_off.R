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

hist(eta, xlab = TeX("Individual quality ($\\epsilon$)"), breaks = 50,
     main = NULL, cex.lab = 1.5)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# cross-lag simulation
# E(phi[t]) is a function of rho[t]
# E(rho[t]) is a function of eta and rho[t-1]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha <- c(2, -0.5, 0.75)
beta <- c(1.25, -1)
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




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# more ordinal regression fun
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p.rho <- beta[2] * 0:2 # calculate expected latent value

prob.tplus1 <- matrix(NA, 3,3)
for (j in 1:3){
  prob.tplus1[j,1] <- plogis(c[1] - p.rho[j])
  prob.tplus1[j,2] <- plogis(c[2] - p.rho[j]) - plogis(c[1] - p.rho[j])
  prob.tplus1[j,3] <- 1 - plogis(c[2] - p.rho[j])
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# expected carry-over effect on reproduction
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(prob.tplus1[1,] ~ c(0:2), ylim = c(0,1), type = 'b', pch = 1,
     ylab = "Probability of n. fawns",
     xlab = 'Number of fawns in t+1', cex = 2, cex.lab = 1.5)
points(prob.tplus1[2,] ~ c(0:2), type = 'b', pch = 19, cex = 2)
points(prob.tplus1[3,] ~ c(0:2), type = 'b', pch = 21, cex = 2, bg = 'grey')
legend('topleft', legend = c('0 fawns in t', '1 fawn in t', '2 fawns in t'),
       pch = c(1,19,21), pt.bg = c(NA, NA, 'grey'), cex = 1.5, bty = 'n')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# observed cost of reproduction
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(jitter(y[,2], 0.5) ~ jitter(y[,1], 0.5), xlab = TeX("$y_{t-1}$"),
     ylab = TeX("$y_{t}$"), cex.lab = 1.5, las = 1)
cor.test(jitter(y[,2], 0.5), jitter(y[,1], 0.5))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# expected carry-over effect on survival
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
e.phi <- plogis(alpha[1] + alpha[2] * c(0:2))
plot(e.phi ~ c(0:2), ylab = 'Expected survival from t-1 to t', type = 'b',
     xlab = 'Fawns in t-1', ylim = c(0.7,0.9), cex.lab = 1.5, cex = 2, las = 1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# observed survival
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1), family = 'sans')
boxplot(phi[,1] ~ y[,1], xlab = 'Number of fawns', boxwex = 0.25, outline = F,
        ylab = 'Survival probability', cex.lab = 1.5, las = 1)



get.l <- function(x){max(which(x == 1))}
l <- apply(z,1,get.l)

loop <- which(l > 1 & l < 10)
all <- which(l == 10)
one <- which(l == 1)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Our JAGS model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink("our_first_latent_variable.jags")
cat("
    model {

    alpha[1] ~ dlogis(0,1)
    alpha[2] ~ dnorm(0, 0.1)
    alpha[3] ~ dnorm(0, 0.1)
    
    beta[1] = 1
    beta[2] ~ dnorm(0, 0.1)

    sigma ~ dgamma(1,1)
    tau = 1/(sigma * sigma)

    # cut-off points for ordinal regression
    c0[1] ~ dnorm(0,0.1)
    c0[2] ~ dnorm(0,0.1)
    c[1:2] <- sort(c0[1:2])

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # individual random effects and first reproduction
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for (i in 1:nI){
      # non-centered trick, same as eta[i] ~ dnorm(0,tau)
      eta.star[i] ~ dnorm(0,1)
      eta[i] = eta.star[i] * sigma 
  
      rho[i,1] = beta[1] * eta[i]
      y[i,1] ~ dordered.logit(rho[i,1], c[1:2])

    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # immediate mortality
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for (i in one){
      logit(phi[i,1]) = alpha[1] + alpha[2] * y[i,1] + alpha[3] * eta[i]
      z[i,2] ~ dbern(phi[i,1])
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # loop to death for individuals that survived at least once
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for (i in loop){
      for (t in 2:l[i]){
        rho[i,t] = beta[1] * eta[i] + beta[2] * y[i,t-1]
        y[i,t] ~ dordered.logit(rho[i,t], c[1:2])
      }
      
      for (t in 2:(l[i] + 1)){
        logit(phi[i,t-1]) = alpha[1] + alpha[2] * y[i,t-1] + alpha[3] * eta[i]
        z[i,t] ~ dbern(z[i,t-1] * phi[i,t-1])
      }
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # loop to end of time series for individuals that survived all 10 years
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    for (i in all){
      for (t in 2:nT){
        rho[i,t] = beta[1] * eta[i] + beta[2] * y[i,t-1]
        y[i,t] ~ dordered.logit(rho[i,t], c[1:2])
      }
      
      for (t in 2:nT){
        logit(phi[i,t-1]) = alpha[1] + alpha[2] * y[i,t-1] + alpha[3] * eta[i]
        z[i,t] ~ dbern(z[i,t-1] * phi[i,t-1])
      }
    }


    }
    ",fill = TRUE)
sink()



jags.data <- list(y = y+1, z = z, nI = nI, nT = nT, 
                  loop = loop, one = one, all = all, l = l)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we provide initial values. This is where JAGS will begin sampling
# to build a posterior distribution, in this case, I won't provide any 
# initial values
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inits <- function(){list(c0 = c(-2,2))}  

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This is where we tell JAGS which parameters we wish to monitor
# i.e., which posterior distributions we want to save, 
# generally we'll save all of them, but 
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
parameters <- c('alpha','beta','sigma','m','y.new')

# number of chains (nc), thinning rate (nt), number of iterations (ni), and number to burn-in
nc <- 4
nt <- 10
ni <- 20000
nb <- 10000


# Call JAGS from R 
# 20k iterations takes awhile on an i9
library(jagsUI)
Sys.time()
m <- jags(jags.data, inits, parameters, "our_first_latent_variable.jags", parallel = T, 
          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
Sys.time()

print(m)








