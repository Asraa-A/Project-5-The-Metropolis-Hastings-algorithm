library(MASS)
library(coda)
library(reshape2)
library(ggplot2)
library(ggmap)
library(osmar)
library(geosphere)
library(mvtnorm)
library(grid)
library(gridExtra)
library(cowplot)


install.packages("devtools")
devtools::install_github("dkahle/ggmap")

library(ggmap)
register_google(key = 'AIzaSyA5glCK4woGllJ6JZkeKgLnBfYh1gGpf4Q')

p0 <- c(-0.945724,51.439892)
p2 <- c(-0.9603548 ,51.5350965)
p3 <- c(-1.0906902,51.4816764)  
p1<- c(-0.898345,51.441257) 

B <- function(lambda1,phi1,lambda2,phi2){
  s <- (lambda2-lambda1)/(phi2-phi1)
  angle <- (360+atan(s)*180/pi) %% 360
}
set.seed(123)
sigma <- rexp(1,20)
error <- rnorm(1,0,sigma)
# lambda0 <- runif(1,-1,0)
# phi0 <- runif(1,51,52)
# p0=c(lambda0,phi0)
alpha0 <- B(p0[1],p0[2],p1[1],p1[2])-error 
#mu_alpha <- mean(alpha)
beta0 <- B(p0[1],p0[2],p2[1],p2[2]) +2* error
#mu_beta <- mean(beta)
gamma0 <- B(p0[1],p0[2],p3[1],p3[2]) -1.5*error
#mu_gamma <- mean(gamma)


#data <-data.frame(alpha=alpha0+c(-1,1)*error,beta=beta0+c(-1,1)*error,gamma=gamma0+c(-1,1)*error)
# mu_alpha <- mean(data$alpha)
# mu_beta <- mean(data$beta)
# mu_gamma <- mean(data$gamma)

##ploting the map

landmarks<-data.frame(lon=c(p1[1],p2[1],p3[1]),lat=c(p1[2],p2[2],p3[2]))
d <- seq(0,0.4,0.0001) # Length of the line
line1 <- data.frame(lon=landmarks[1,1] + d*sin(alpha0*pi/180+pi),
                    lat=landmarks[1,2] + d*cos(alpha0*pi/180+pi))
line2 <- data.frame(lon=landmarks[2,1] + d*sin(beta0*pi/180+pi),
                    lat=landmarks[2,2] + d*cos(beta0*pi/180+pi))
line3 <- data.frame(lon=landmarks[3,1] + d*sin(gamma0*pi/180+pi),
                    lat=landmarks[3,2] + d*cos(gamma0*pi/180+pi))
map<- get_map(c(-0.9456352,51.4398685),zoom=18,maptype = "road")

#map <- get_map(c(-0.9456352,51.4398685),zoom=11,maptype="watercolor")
mapPlot <- ggmap(map)+
  geom_point(aes(x = lon, y = lat), size = 1, data = landmarks, alpha = .5) +
  geom_line(aes(x=lon,y=lat),data=line1) +
  geom_line(aes(x=lon,y=lat),data=line2) +
  geom_line(aes(x=lon,y=lat),data=line3)
intersectBearings <- function(p1,b1,p2,b2) {
  x1 <- p1[1]
  x2 <- p1[1] + 0.1*sin(b1*pi/180)
  x3 <- p2[1]
  x4 <- p2[1] + 0.1*sin(b2*pi/180)
  y1 <- p1[2]
  y2 <- p1[2] + 0.1*cos(b1*pi/180)
  y3 <- p2[2]
  y4 <- p2[2] + 0.1*cos(b2*pi/180)
  x <- ((x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))
  y <- ((x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4))/((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4))
  return(as.numeric(c(x,y)))
}

intersection <- intersectBearings(landmarks[1,],alpha0,landmarks[2,],beta0)

################################################################################
loglikelihood <- function(theta,p1.=p1,p2.=p2,p3.=p3){
  lambda=theta[1]
  phi=theta[2]
  sigma= theta[3]
  
  if (sigma < 0){
    sigma <- 0 }
  
  alpha <- B(lambda,phi,p1[1],p1[2])
  beta <-B(lambda,phi,p2[1],p2[2])
  gamma <- B(lambda,phi,p3[1],p3[2])
  R <- sum(dnorm(alpha0,alpha,sigma,log=T)+
             dnorm(beta0,beta,sigma,log=T)+
             dnorm(gamma0,gamma,sigma,log=T))
  return(R)
}

logprior <- function(theta){
  lambda=theta[1]
  phi=theta[2]
  sigma= theta[3]
  if (sigma < 0){
    sigma <- 0 }
  lambda_prior <- sum(dunif(lambda,-1.5,-0.5,log=T))
  phi_prior <- sum(dunif(phi,51,52,log=T))
  sigma_prior <- sum(dexp(sigma,20,log=T))
  return (lambda_prior+phi_prior+sigma_prior)
}
##############################################################################################
#Metropolis-Hastings aribtrary variance
set.seed(12345)
draws <- array(0,dim=c(4000,3,3))
draws[4000,1,] <- runif(3,intersection[1] - 0.01, intersection[1] + 0.01)
draws[4000,2,] <- runif(3,intersection[2] - 0.01, intersection[2] + 0.01)
draws[4000,3,] <- rexp(3,20)

#covariance matrix
#prop.cov <- c(1e-6,1e-6,1e-2)*diag(3)
#prop.cov <- c(1e-8,1e-8,1e-4)*diag(3)
prop.cov <- c(1e-8,1e-8,0.001)*diag(3)
#prop.cov <- c(1e-10,1e-10,1e-4)*diag(3)
#prop.cov <- c(2.799102e-8,9.678784e-8, 0.015)*diag(3)
accepted <- rep(1,3)
converged <- FALSE

while (!converged) {
  draws[1,,] <- draws[4000,,]
  for (step in 2:4000) {
    for (chain in 1:3) {
      proposed <- rmvnorm(1,draws[step-1,,chain],prop.cov)
      r <- loglikelihood(proposed,p1,p2,p3)+logprior(proposed)-
            loglikelihood(draws[step-1,,chain],p1,p2,p3)-
            logprior(draws[step-1,,chain])
      logalpha <- min(0,r)
      u <- runif(1)
      if (log(u)< logalpha) {
        draws[step,,chain] <- proposed
        accepted[chain] <- accepted[chain]+1
      } else {
        draws[step,,chain] <- draws[step-1,,chain]
      }
    }
  }
  print(sprintf("Acceptance rate: %f",accepted/120))
  chainlist <- mcmc.list(Chain1=mcmc(draws[,,1]),
                         Chain2=mcmc(draws[,,2]),
                         Chain3=mcmc(draws[,,3]))
  converged <- all((gelman.diag(chainlist)$psrf[,2])<1.1)
  plot(chainlist) # This plots current state of the chains
  Sys.sleep(0.1)
}
############################
#Metropolis hasting,varinace from previous run
prop.cov1 <-cov(draws[,,1])
samples <- array(0,dim=c(4000,3,3))
samples[4000,1,] <- sample(draws[,1,1],3)
samples[4000,2,] <- sample(draws[,2,1],3)
samples[4000,3,] <- sample(draws[,3,1],3)
accepted <- rep(1,3)
converged <- FALSE

while (!converged) {
  samples[1,,] <- samples[4000,,]
  for (step in 2:4000) {
    for (chain in 1:3) {
      proposed <- rmvnorm(1,samples[step-1,,chain],prop.cov1)
      r <- loglikelihood(proposed,p1,p2,p3)+logprior(proposed)-
        loglikelihood(samples[step-1,,chain],p1,p2,p3)-
        logprior(samples[step-1,,chain])
      logalpha <- min(0,r)
      u <- runif(1)
      if (log(u)< logalpha) {
       samples[step,,chain] <- proposed
        accepted[chain] <- accepted[chain]+1
      } else {
        samples[step,,chain] <- samples[step-1,,chain]
      }
    }
  }
  print(sprintf("Acceptance rate: %f",accepted/120))
  chainlist <- mcmc.list(Chain1=mcmc(samples[,,1]),
                         Chain2=mcmc(samples[,,2]),
                         Chain3=mcmc(samples[,,3]))
  converged <- all((gelman.diag(chainlist)$psrf[,2])<1.2)
  plot(chainlist) # This plots current state of the chains
  Sys.sleep(0.1)
}
#########################
#one converged we start collecting samples from the posterior
#Once convergence is achieved,
# we begin collecting the final sample from the posterior
final_sample <- array(0,dim=c(100000,3,3))
final_sample[1,,] <- samples[4000,,]
for (step in 2:100000) {
  for (chain in 1:3) {
    proposed <- rmvnorm(1,final_sample[step-1,,chain],prop.cov1)
    r <- loglikelihood(proposed,p1,p2,p3)+logprior(proposed)-
      loglikelihood(final_sample[step-1,,chain],p1,p2,p3)-
      logprior(final_sample[step-1,,chain])
    logalpha <- min(0,r)
    u <- runif(1)
    if (log(u)< logalpha) {
      final_sample[step,,chain] <- proposed
      accepted[chain] <- accepted[chain]+1
    } else {
      final_sample[step,,chain] <- final_sample[step-1,,chain]
    }
  }
  
}
finalsample <- rbind(final_sample[,,1],final_sample[,,2],final_sample[,,3])


###########################################################################
#Plotting the Density Estimate on the Map
x11();dev.off()
D <- kde2d(as.vector(final_sample[,1,]),as.vector(final_sample[,2,]),
           h=c(sd(final_sample[,1,]),sd(final_sample[,2,])),
           n=1024,
           lims=c(-4.9055,-4.9035,55.7915,55.7935)) # Enough to cover map
z <- melt(D$z)
z$Var1<-D$x[z$Var1]
z$Var2<-D$y[z$Var2]
map <- get_map(c(mean(final_sample[,1,]),mean(final_sample[,2,])),zoom=18,maptype="road")
mapPoints <- ggmap(map)+
  geom_point(aes(x = lon, y = lat), size = 1, data = landmarks, alpha = .5) +
  geom_raster(data=z,aes(x=Var1,y=Var2,fill=value))+
  guides(fill=FALSE,alpha=FALSE)+
  scale_fill_gradientn(colours=c("#0000FF00","#0000FFFF"))+
  coord_cartesian() +
  geom_line(aes(x=lon,y=lat),data=line1) +
  geom_line(aes(x=lon,y=lat),data=line2) +
  geom_line(aes(x=lon,y=lat),data=line3)+
  geom_point(aes(x=lon,y=lat),
             data=data.frame(lon=mean(final_sample[,1,]),lat=mean(final_sample[,2,])),
             size=0.5,colour="#FF0000")
mapPoints
