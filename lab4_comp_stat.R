### COMP STAT LAB 4 ####


# rm(list = ls())
library(ggplot2)
############### EXCERCISE 1 ###################


# target probability density function 


f_target <- function(x){
  return((x^5)*exp(-x))
}



#########
### 1 ###
#########
# Metropolis Hastings
# Need to change the plotting because the horizontal lines are calculated from the original code (2sd away from mean in each direction)
#In principal, this works
f.MCMC.MH.lnorm<-function(nstep,X0,props){
  vN<-1:nstep
  vX<-rep(X0,nstep);
  for (i in 2:nstep){
    X<-vX[i-1]
    Y<-rlnorm(1,meanlog = X,sdlog=props)
    u<-runif(1)
    a<-min(c(1,(dlnorm(Y)*dlnorm(X,meanlog=Y,sdlog=props))/(dlnorm(X)*dlnorm(Y,meanlog=X,sdlog =props))))
    if (u <=a){vX[i]<-Y}else{vX[i]<-X}    
  }
  plot(vN,vX,pch=19,cex=0.3,col="black",xlab="t",ylab="X(t)",main="",ylim=c(min(X0-0.5,-5),max(5,X0+0.5)))
  abline(h=0)
  abline(h=1.96)
  abline(h=-1.96)
}

f.MCMC.MH.lnorm(nstep = 50000, X0 = 1, props = f_target(1))


#No burn-in perdiod, no convergences?

#########
### 2 ###
#########
# Doesn't work. Might be worth writing our own MH algorithm based on the slides.
# Not sure what I am supposed to floor here
f.MCMC.MH.chisq<-function(nstep,X0,props){
  vN<-1:nstep
  vX<-rep(X0,nstep);
  for (i in 2:nstep){
    X<-vX[i-1]
    Y<-rchisq(1,df = floor(X))
    u<-runif(1)
    a<-min(c(1,(dchisq(Y)*dchisq(X,df = floor(Y)))/(dchisq(X)*dchisq(Y,df=floor(X)))))
    if (u <=a){vX[i]<-Y}else{vX[i]<-X}    
  }
  plot(vN,vX,pch=19,cex=0.3,col="black",xlab="t",ylab="X(t)",main="",ylim=c(min(X0-0.5,-5),max(5,X0+0.5)))
  abline(h=0)
  abline(h=1.96)
  abline(h=-1.96)
}

f.MCMC.MH.chisq(nstep = 10000, X0 = 1, props = f_target(1))

#########
### 3 ###
#########

# Probably the second distribution works pretty well and we will see the burn in period and everything converging nicely.

#########
### 4 ###
#########
library(coda) #this library has a gelman.diag()function to calculate the gelman-rubin factor immediately
#see code for slide 22

## Bullshit to have an idea of what to do
for (i in 1:10){
  do_MCMC_chisq(nstep = 10000, X0 = i, props = f(target(i)))
}

# Store itemas a mcmc.list
mcmc_item = mcmc.list()
for(i in 1:10){
  mcmc_item[[i]] = as.mcmc(do_MCMC_chisq[[i]]) #proably we need to specify the sample here
}

gelman_factor = gelman.diag(mcmc_item)

## Result is probably close to 1 and shows convergence


#########
### 5 ###
#########


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

plot1 <- ggplot(data, aes(x = X, y = Y))+
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

# See code for slides 18 and 19

#########
### 5 ###
#########

# See code for slides 18 and 19