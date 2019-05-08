library(tidyverse)
library(invgamma)
# Height in cm
olympics <- read_csv("C:/Users/Cloe Wynn/Documents/BYU/Senior/Winter 2019/Stat 251/athlete_events.csv",
                     col_types = list(ID = col_character(),
                                      Year = col_character()))
str(olympics)

#females gymnastics vs female swimmers
gymnastics <- subset(olympics, Sport == "Gymnastics" & Sex == "F" 
                     & !is.na(olympics["Weight"]) & !is.na(olympics["Height"]))
swimming <- subset(olympics, Sport == "Swimming" & Sex == "F" 
                   & !is.na(olympics["Weight"]) & !is.na(olympics["Height"]))
#remove duplicate athletes
gymnastics <- gymnastics[!duplicated(gymnastics$Name), ]
swimming <- swimming[!duplicated(swimming$Name), ]

#5'0" 2"   152 cm   5 cm
#5'6" 2 "  168 cm   5 cm 

mcmc.func <- function(x, J, m, v, a, b) {
  
  out <- NULL
  
  n <- length(x)
  ybar <- mean(x)
  s2 <- var(x)
  
  # empty vectors to hold draws from posterior
  mu <- numeric()
  sig2 <- numeric()
  
  # starting values
  mu[1] <- ybar
  sig2[1] <- s2
  
  # make for loop to update everything
  for(j in 2:J) {
    
    # update mu
    mstar <- (n*ybar*v+m*sig2[j-1])/(n*v + sig2[j-1]) 
    vstar <- (v*sig2[j-1])/(n*v + sig2[j-1])
    
    mu[j] <- rnorm(1, mstar, sqrt(vstar))
    
    # update sig2
    astar <- 0.5*n + a
    bstar <- 0.5*sum((x - mu[j])^2) + b
    
    sig2[j] <- rinvgamma(1, astar, rate = bstar)
  }
  
  post.pred <- rnorm(J, mu, sqrt(sig2))
  prior.pred <- rnorm(J, rnorm(J, m, sqrt(v)), sqrt(rinvgamma(J, a, b)))
  
  out$mu <- mu
  out$sig2 <- sig2
  out$post.pred <- post.pred
  out$prior.pred <- prior.pred
  
  out
}

#IG(1,1)
sd(gymnastics$Height)
sd(swimming$Height)

# Prior for gymnasts
m1 <- 152
v1 <- 5^2

out.g <- mcmc.func(gymnastics$Height, 100000, m1, v1, a = 1, b = 1)
# Prior distribution
plot(density(rnorm(100000, m1, sqrt(v1))), lwd = 2, col = "red", 
     main = expression(paste("Prior for Olympic Gymnasts' Heights (", mu, ")")),
     xlab = "Support")
# Prior and posterior for gymnasts
plot(density(out.g$mu), lwd = 2, xlab = "Support", ylab = "Density",
     main = expression(paste("Prior & Posterior of Olympic Gymnasts' Heights (", mu, ")")))
lines(density(rnorm(100000, m1, sqrt(v1))), lwd = 2, col = "red")
legend("topright", legend = c("Prior", "Posterior"), col = c("red", "black"), 
       lty = 1:1, cex = 0.8)
# 95% credible interval
quantile(out.g$mu, c(0.025, 0.975))

# Prior for swimmers
m2 <- 168
v2 <- 5^2

out.s <- mcmc.func(swimming$Height, 100000, m2, v2, a = 1, b = 1)
# Prior distribution (mu)
plot(density(rnorm(100000, m2, sqrt(v2))), lwd = 2, col = "red",
     main = expression(paste("Prior for Olympic Swimmers' Heights (", mu, ")")),
     xlab = "Support")
# Prior and posterior for swimmers (mu)
plot(density(out.s$mu), lwd = 2, xlab = "Support", ylab = "Density",
     main = expression(paste("Prior & Posterior of Olympic Swimmers' Heights (", mu, ")")))
lines(density(rnorm(100000, m2, sqrt(v2))), lwd = 2, col = "red")
legend("topright", legend = c("Prior", "Posterior"), col = c("red", "black"), 
       lty = 1:1, cex = 0.8)
# 95% credible interval
quantile(out.s$mu, c(0.025, 0.975))

# Posterior of Gymnasts and Swimmers (mu)
plot(density(out.g$mu), lwd = 2, main = "Posterior of Olympic Gymnasts & Swimmers", 
     xlab = "Support", ylab = "Density", col = "orange", 
     xlim = c(155, 171), ylim = c(0, 3))
lines(density(out.s$mu), lwd = 2, col = "blue")
legend("topleft", legend = c("Gymnasts", "Swimmers"), col = c("orange", "blue"), 
       lty = 1:1, cex = 0.8)

# Difference of mus
diffs <- out.s$mu - out.g$mu
plot(density(diffs), main = "Difference of Means", xlab = "Difference", ylab = "Density",
     lwd = 2)
quantile(diffs, c(0.025, 0.975))