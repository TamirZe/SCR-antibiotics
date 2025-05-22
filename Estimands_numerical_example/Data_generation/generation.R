library(ggplot2)
library(dplyr)
# shape = k (alpha), scale = 1/lambda
lambda = 2; k=4
rate = lambda; scale = 1/lambda

exp = rexp(n=100000, rate = rate)
weib = rweibull(n=100000, shape=k, scale = scale)
mean(weib); var(weib); mean(exp); var(exp)

U = runif(100000)
X = (-log(1-U))^(1/k) * scale
mean(X); var(X)

dat <- data.frame(dens = c(exp, weib, X)
                  , lines = rep(c( "exp", "weib", "X"), each = length(exp)))
ggplot(dat %>% filter(lines != "X"), aes(x = dens, fill = lines)) + geom_density(alpha = 0.5)
ggplot(dat %>% filter(lines != "exp"), aes(x = dens, fill = lines)) + geom_density(alpha = 0.5)
