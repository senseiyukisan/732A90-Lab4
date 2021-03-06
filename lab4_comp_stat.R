### COMP STAT LAB 4 ####

# rm(list = ls())
library(ggplot2)
############### EXCERCISE 1 ###################


# target probability density function 
f_target  =  function(x){
  return((x^5)*exp(-x))
}


#########
### 1 ###
#########

f.MCMC.MH.lnorm = function(nstep, X0, props){
  vN = 1:nstep
  vX = rep(X0, nstep);
  for (i in 2:nstep){
    X = vX[i-1]
    Y = rlnorm(1, log(X), sdlog=props)
    u = runif(1)
    a = min(c(1, (f_target(Y)*dlnorm(X, meanlog=log(Y), sdlog=props)) / (f_target(X)*dlnorm(Y, meanlog=log(X), sdlog=props))))
    if (u <= a) {
      vX[i] = Y
    } else {
      vX[i] = X
    }    
  }
  return(vX)
}

set.seed(12345)
lnorm_vals = f.MCMC.MH.lnorm(nstep = 500, X0 = 1, props = f_target(1))

lnorm_plot = ggplot(data = data.frame(lnorm_vals), aes(x = 1:length(lnorm_vals), y = lnorm_vals)) +
  geom_line(color="darkgreen") +  
  xlab("t") + ylab("x(t)")

lnorm_plot

#########
### 2 ###
#########

f.MCMC.MH.chisq = function(nstep, X0){
  vN = 1:nstep
  vX = rep(X0,nstep)
  for (i in 2:nstep){
    X = vX[i-1]
    Y = rchisq(1, df=floor(X+1))
    u = runif(1)
    a = min(c(1, (f_target(Y)*dchisq(X, df=floor(X+1))) / (f_target(X)*dchisq(Y, df=floor(X+1)))))
    if (u <= a){
      vX[i] = Y
    } else {
      vX[i] = X
    }    
  }
  return(vX)
}

set.seed(12345)
chisq_vals = f.MCMC.MH.chisq(nstep = 10000, X0 = 1)


chisq_plot = ggplot(data = data.frame(chisq_vals), aes(x = 1:length(chisq_vals), y = chisq_vals)) +
  geom_line(color="darkgreen") +  
  xlab("t") + ylab("x(t)")

chisq_plot

#########
### 3 ###
#########

# Probably the second distribution works pretty well and we will see the burn in period and everything converging nicely.
library(gridExtra)

# Plotting only first 1000 values for a better view on burn-in period.
lnorm_vals_burn_in = lnorm_vals[1:1000]
lnorm_burn_in = ggplot(data = data.frame(lnorm_vals_burn_in), aes(x = 1:length(lnorm_vals_burn_in), y = lnorm_vals_burn_in)) +
  geom_line(color="darkgreen") +  
  xlab("t") + ylab("x(t)") + ggtitle("lnorm")

chisq_vals_burn_in = chisq_vals[1:1000]
chisq_burn_in = ggplot(data = data.frame(chisq_vals_burn_in), aes(x = 1:length(chisq_vals_burn_in), y = chisq_vals_burn_in)) +
  geom_line(color="darkgreen") +  
  xlab("t") + ylab("x(t)") + ggtitle("chi-square")

grid.arrange(lnorm_burn_in, chisq_burn_in, ncol=2)

#########
### 4 ###
#########
library(coda)

last_start_value = 10
mcmc_list = mcmc.list()

for (i in 1:last_start_value) {
  mcmc_list[[i]] = as.mcmc(f.MCMC.MH.chisq(nstep = 10000, X0 = i))
}

gelman.diag(mcmc_list)

#########
### 5 ###
#########

# We get the estimate by calculating the mean of our samples from 1 and 2
# lnorm
lnorm_estimate = mean(lnorm_vals[-(1:1000)])

# chisquare
chisq_estimate = mean(chisq_vals[-(1:1000)])

#########
### 6 ###
#########

# Question 5 is actually the calculation of the expectation of the gamma distrib that is given in the beginning. So the expectation of a gamma
# distrib is just alpha*beta so the result is 6

######################### Exercise 2 #####################


#########
### 1 ###
#########

load("chemical.RData")
data = data.frame("X" = X,"Y" = Y)

plot1  =  ggplot(data, aes(x = X, y = Y))+
  geom_point()+
  geom_line()
plot1


## Looks quadratic: positive a0 and negative a1

#########
### 2 ###
#########

## Need to write this down by hand on some paper first, but I suspect that this is a product of normal distributions that when you compute it out transforms into a 
## factor in the front to the power of how many steps you made and the product transforms into a sum inside the exponential term


#########
### 3 ###
#########

# THe posterior is the product of the previous two things that we calculated in the previous question.
# Porbably get some ugly expression
# THen we are asked to calculate the first, last and all the middle point posteriors. My guess is that the hints B and C are used in those

#########
### 4 ###
#########

num_values = 1000

f.MCMC.Gibbs<-function(nstep, X0, Y, variance){
  vN<-1:nstep
  mX<-matrix(0,nrow=nstep,ncol=length(Y))
  mX[1,]<-X0
  length_y = length(Y)
  for (i in 2:nstep){
    prev_sample<-mX[i-1,]
    curr_sample<-rep(0,length(Y))
    curr_sample[1]<-rnorm(1, mean=((prev_sample[2]+Y[1])/2), sqrt(variance)/2)
    for (j in 2:(length_y-1)){
      curr_sample[j]<-rnorm(1, mean=(Y[j] + prev_sample[j-1] + prev_sample[j+1])/3, sqrt(variance)/3)
    }
    curr_sample[length_y]<-rnorm(1, mean=((Y[length_y]+prev_sample[length_y-1])/2), sqrt(variance)/2)
    mX[i,]<-curr_sample
  }
  mX
}

gibbs_vals = f.MCMC.Gibbs(num_values, rep(0, length(Y)), Y, variance=0.2)
column_means = colMeans(gibbs_vals)

data = cbind(data, column_means)

plot1  =  ggplot(data) +
  geom_line(aes(x=X, y=Y), color="darkgreen") +
  geom_line(aes(x=X, y=column_means), color="red") + ggtitle("true ~ predicted")
plot1

#########
### 5 ###
#########

mu_n = gibbs_vals[,50]

chisq_burn_in = ggplot(data = mu_n, aes(x = 1:length(mu_n), y = mu_n)) +
  geom_line(color="darkgreen") +  
  xlab("t") + ylab("x(t)") + ggtitle("trace plot mu_n")
