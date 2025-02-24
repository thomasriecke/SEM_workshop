# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# a simple script to simulate correlated y and x data and calculate 
# variances, covariances, slopes and intercepts both by hand
# and using various simple functions in R...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Univariate regression example
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate some data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n <- 100
x <- rnorm(n, 0, 1)
y <- rnorm(n, x, 1)

plot(y ~ x, las = 1, cex.lab = 1.5, las = 1)
cov(data.frame(y,x))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate the variance of y (slide 5)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
var(y)
sum((y-mean(y))^2)/(n-1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate the variance of x (slide 6)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
var(x)
sum((x-mean(x))^2)/(n-1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate the covariance of y & x (slide 7)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cov(y,x)
sum((x-mean(x)) * (y - mean(y)))/(n-1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate the correlation between y & x (slide 8)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cor(y,x)
sum((x-mean(x)) * (y - mean(y)))/(n-1)/(sd(y) * sd(x))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate the correlation between y & x (slide 8)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lm(y ~ x)$coefficients[2]
cov(y,x)/var(x)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate the correlation between y & x (slide 9)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lm(y ~ x)$coefficients[2]
cov(y,x)/var(x)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate the intercept of y(slide 10)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lm(y ~ x)$coefficients[1]
mean(y) - cov(y,x)/var(x) * mean(x)





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# two covariate figure
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n <- 100
x <- matrix(NA, n, 2)
x[,1] <- rnorm(n, 0, 1)
x[,2] <- rnorm(n, 0, 1)
y <- rnorm(n, rowSums(x), 1)
par(mfrow = c(2,1), mar = c(4.1,4.1,1.1,2.1), oma = c(0,0,0,0))
plot(y ~ x[,1], xlab = TeX("$x_1$"))
plot(y ~ x[,2], xlab = TeX("$x_2$"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Path analysis example...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(piecewiseSEM)
n <- 100
x <- rnorm(n, 0, 1)
y <- rnorm(n, x, 1)
z <- rnorm(n, x + y, 1)
d <- data.frame(x,y,z)
sd <- data.frame(x = scale(x), y = scale(y), z = scale(z))

# correlation and variance-covariance matrices (unstandardized)
R <- cor(d)
V <- cov(d)

# correlation and variance-covariance matrices (standardized)
# R <- cor(sd)
# V <- cov(sd)

# correlation plot
plot(d, las = 1, cex.lab = 1.5, las = 1)

path1 <- psem(
  lm(y ~ x, data = d),
  lm(z ~ y + x, data = d),
  data = d
)
summary(path1)

# Calculating the relationship between y and x
V[1,2]/V[1,1]
R[1,2] * sqrt(V[2,2]/V[1,1])
summary(path1)$coefficients$Estimate[1]

# Calculating the relationship between z and y
(R[2,3] - (R[1,2]*R[1,3]))/(1 - R[1,2]^2) * sqrt(V[3,3]/V[2,2])
summary(path1)$coefficients$Estimate[2]

# Calculating the relationship between z and x
(R[1,3] - (R[1,2]*R[2,3]))/(1 - R[1,2]^2) * sqrt(V[3,3]/V[1,1])
summary(path1)$coefficients$Estimate[3]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# we can just run these models together
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary(lm(y ~ x, dat = d))

summary(lm(z ~ x + y, dat = d))

path1 <- psem(
  lm(y ~ x, data = d),
  lm(z ~ y + x, data = d),
  data = d
)
summary(path1)




