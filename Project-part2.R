roads.box <- center_bbox(center_lon=intersection[1],
                         center_lat=intersection[2],
                         width=100,
                         height=100)
api <- osmsource_api(url="https://api.openstreetmap.org/api/0.6/")
roads <- get_osm(roads.box, source=api)
ways <- find(roads, way(tags(k=="highway")))
ways <- find_down(roads, way(ways))
ways <- subset(roads, ids=ways)
hw_lines <- as_sp(ways,"lines")

distanceToRoad <- function(par) {
  distance <- dist2Line(c(par[1],par[2]), hw_lines)
  proper.distance <- distance[1] * 3.28084
  return(proper.distance)
}


rho0 <-distanceToRoad(p0)


###################################################
#our prior
logprior2<- function(theta){
  lambda=theta[1]
  phi=theta[2]
  sigma=theta[3]
  
  if (sigma < 0){
    sigma <- 0 }
  
  p <- c(lambda,phi)
  rho <- distanceToRoad(p)
  lambda_phi_prior <- sum(dnorm(rho,0,6,log=T))
  sigma_prior <- sum(dexp(sigma,20,log=T))
  return (lambda_phi_prior+sigma_prior)
}

loglikelihood2 <- function(theta,p1.=p1,p2.=p2,p3.=p3){
  lambda=theta[1]
  phi=theta[2]
  p=c(lambda,phi)
  sigma=theta[3]
  
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

logpost <- logprior2(theta)+loglikelihood2(theta)
####################################################################
set.seed(54321)
draws2 <- array(0,dim=c(4000,3,3))
draws2[4000,1,] <- runif(3,intersection[1] - 0.01, intersection[1] + 0.01)
draws2[4000,2,] <- runif(3,intersection[2] - 0.01, intersection[2] + 0.01)
draws2[4000,3,] <- rexp(3,20)

#covariance matri
#prop.cov2<- c(1e-10,1e-10,1e-6)*diag(3)
prop.cov2 <- rbind(c(3.283647e-08, 3.850296e-09 ,2.620361e-06),c(3.850296e-09, 3.475803e-09 ,1.744348e-08),c(2.620361e-06 ,1.744348e-08 ,3.663148e-04))
accepted2 <- rep(1,3)
converged2 <- FALSE

while (!converged2) {
  draws2[1,,] <- draws2[4000,,]
  for (step in 2:4000) {
    for (chain in 1:3) {
      proposed <- rmvnorm(1,draws2[step-1,,chain],prop.cov2)
      r <- loglikelihood2(proposed,p1,p2,p3)+logprior2(proposed)-
        loglikelihood2(draws2[step-1,,chain],p1,p2,p3)-
        logprior2(draws2[step-1,,chain])
      logalpha <- min(0,r)
      u <- runif(1)
      if (log(u)< logalpha) {
        draws2[step,,chain] <- proposed
        accepted2[chain] <- accepted2[chain]+1
      } else {
        draws2[step,,chain] <- draws2[step-1,,chain]
      }
    }
  }
  print(sprintf("Acceptance rate: %f",accepted/120))
  chainlist <- mcmc.list(Chain1=mcmc(draws2[,,1]),
                         Chain2=mcmc(draws2[,,2]),
                         Chain3=mcmc(draws2[,,3]))
  converged2 <- all((gelman.diag(chainlist)$psrf[,2])<1.2)
  plot(chainlist) # This plots current state of the chains
  Sys.sleep(0.1)
}
############################
prop.cov22 <-cov(draws2[2000:4000,,1])
samples2 <- array(0,dim=c(4000,3,3))
samples2[4000,1,] <- sample(draws2[,1,1],3)
samples2[4000,2,] <- sample(draws2[,2,2],3)
samples2[4000,3,] <- sample(draws2[,3,3],3)

accepted2 <- rep(1,3)
converged2 <- FALSE

while (!converged2) {
  samples2[1,,] <- samples2[4000,,]
  for (step in 2:4000) {
    for (chain in 1:3) {
      proposed2 <- rmvnorm(1,samples2[step-1,,chain],prop.cov22)
      r <- loglikelihood2(proposed2,p1,p2,p3)+logprior2(proposed2)-
        loglikelihood2(samples2[step-1,,chain],p1,p2,p3)-
        logprior2(samples2[step-1,,chain])
      logalpha <- min(0,r)
      u <- runif(1)
      if (log(u)< logalpha) {
        samples2[step,,chain] <- proposed2
        accepted2[chain] <- accepted2[chain]+1
      } else {
        samples2[step,,chain] <- samples2[step-1,,chain]
      }
    }
  }
  print(sprintf("Acceptance rate: %f",accepted2/120))
  chainlist <- mcmc.list(Chain1=mcmc(samples2[,,1]),
                         Chain2=mcmc(samples2[,,2]),
                         Chain3=mcmc(samples2[,,3]))
  converged2 <- all((gelman.diag(chainlist)$psrf[,2])<1.2)
  plot(chainlist) # This plots current state of the chains
  Sys.sleep(0.1)
}
#########################
#one converged2 we start collecting samples2 from the posterior
#Once convergence is achieved,
# we begin collecting the final sample from the posterior
final_sample2 <- array(0,dim=c(20000,3,3))
final_sample2[1,,] <- samples2[4000,,]
for (step in 2:20000) {
  for (chain in 1:3) {
    proposed2 <- rmvnorm(1,final_sample2[step-1,,chain],prop.cov22)
    r <- loglikelihood2(proposed2,p1,p2,p3)+logprior2(proposed2)-
      loglikelihood2(final_sample2[step-1,,chain],p1,p2,p3)-
      logprior2(final_sample2[step-1,,chain])
    logalpha <- min(0,r)
    u <- runif(1)
    if (log(u)< logalpha) {
      final_sample2[step,,chain] <- proposed2
      accepted2[chain] <- accepted2[chain]+1
    } else {
      final_sample2[step,,chain] <- final_sample2[step-1,,chain]
    }
  }
  
}
finalsample2 <- rbind(final_sample2[,,1],final_sample2[,,2],final_sample2[,,3])
###############################################################################
###########################################################################
#Plotting the Density Estimate on the Map
x11();dev.off()
D <- kde2d(as.vector(final_sample2[,1,]),as.vector(final_sample2[,2,]),
           h=c(sd(final_sample2[,1,]),sd(final_sample2[,2,])),
           n=1024,
           lims=c(-4.9055,-4.9035,55.7915,55.7935)) # Enough to cover map
z <- melt(D$z)
z$Var1<-D$x[z$Var1]
z$Var2<-D$y[z$Var2]
map <- get_map(c(-0.945724,51.439892),zoom=18,maptype="road")
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
             data=data.frame(lon=mean(final_sample2[,1,]),lat=mean(final_sample2[,2,])),
             size=0.5,colour="#FF0000")
mapPoints

