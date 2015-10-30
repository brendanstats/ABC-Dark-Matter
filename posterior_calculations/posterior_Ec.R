###################
#Advanced Data Analysis
#Brendan McVeigh
#July 28, 2015
###################

#change ~line 4

setwd("/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/posterior_calculations")

#########################################################
#Define Functions
#########################################################

vesc <- function(x, r.lim = 1.5, vmax = 21, rmax = 1.5){
  rs = rmax/2.16
  Ps = (vmax/0.465)^2
  
  xlim = r.lim/rs  	 #turn rlim in unit of rs
  Pr   = Ps * ( 1 - (log(1+x))/x )
  Plim = Ps * ( 1 - (log(1+xlim))/xlim ) #lms 0.45*Ps
  
  vesc = (2 * (Plim - Pr))^0.5
  return(vesc)
}

hprob <- function(v, x, a = 2.0, d = -5.3, e = 2.5, Ec = .16, r.lim = 1.5,
                  b = -9.0, q = 6.9, Jb = 8.6*10^(-2), vmax = 21, rmax = 1.5){
  
  phi.s <- (vmax/.465)^2
  r.s <- rmax/2.16  
  phi.r <- phi.s * (1-log(1+x)/x)
  
  E <- v^2/2 + phi.r
  Ec <- Ec * phi.s
  xlim <- r.lim/r.s
  phi.lim <- phi.s * (1-log(1+xlim)/xlim)
  
  if(E<phi.lim){
    hE <- E^a*(E^q+Ec^q)^(d/q)*(phi.lim-E)^e
  } else hE = 0
  #return(hE*v^2*x^2)
  return(hE*v^2*x^2*r.s)
}

hprob.vec.v <- Vectorize(hprob, vectorize.args = "v")

intv.dens <- function(x, a = 2.0, d = -5.3, e = 2.5, Ec = .16, r.lim = 1.5,
                      b = -9.0, q = 6.9, Jb = 8.6*10^(-2), vmax = 21, rmax = 1.5){
  out <- integrate(hprob.vec.v, lower = 0, upper = vesc(x), x, a=a, d=d, e=e, Ec=Ec, r.lim=r.lim,
                   q=q, b=b, Jb=Jb, vmax=vmax, rmax=rmax)
  return(out$value)
}

intv.dens.vec.x <- Vectorize(intv.dens, vectorize.args = "x")

norm <- function(a = 2.0, d = -5.3, e = 2.5, Ec = .16, r.lim = 1.5,
                 b = -9.0, q = 6.9, Jb = 8.6*10^(-2), vmax = 21, rmax = 1.5){
  N <- integrate(intv.dens.vec.x, lower = 0, upper = r.lim/r.s, a=a, d=d, e=e, Ec=Ec, r.lim=r.lim,
          q=q, b=b, Jb=Jb, vmax=vmax, rmax=rmax)$value
  return(N)
}

hprob.norm <- function(xv, a = 2.0, d = -5.3, e = 2.5, Ec = .16, r.lim = 1.5,
                       b = -9.0, q = 6.9, Jb = 8.6*10^(-2), vmax = 21, rmax = 1.5){
  x <- xv[1]
  v <- xv[2]
  N <- norm(a, d, e, Ec, r.lim, b, q, Jb, vmax, rmax)
  prob <- hprob(v, x, a, d, e, Ec, r.lim, b, q, Jb, vmax, rmax)/N
  return(prob)
}

hprob.xv <- function(xv, a = 2.0, d = -5.3, e = 2.5, Ec = .16, r.lim = 1.5,
                  b = -9.0, q = 6.9, Jb = 8.6*10^(-2), vmax = 21, rmax = 1.5){
  x <- xv[1]
  v <- xv[2]
  phi.s <- (vmax/.465)^2
  r.s <- rmax/2.16  
  phi.r <- phi.s * (1-log(1+x)/x)
  
  E <- v^2/2 + phi.r
  Ec <- Ec * phi.s
  xlim <- r.lim/r.s
  phi.lim <- phi.s * (1-log(1+xlim)/xlim)
  
  if(E<phi.lim){
    hE <- E^a*(E^q+Ec^q)^(d/q)*(phi.lim-E)^e
  } else hE = 0
  return(hE*v^2*x^2)
}

nll.norm <- function(xv, a = 2.0, d = -5.3, e = 2.5, Ec = .16, r.lim = 1.5,
                     b = -9.0, q = 6.9, Jb = 8.6*10^(-2), vmax = 21, rmax = 1.5){
  
  liklihood <- apply(xv,1,hprob.xv, a=a, d=d, e=e, Ec=Ec, r.lim=r.lim,
                     q=q, b=b, Jb=Jb, vmax=vmax, rmax=rmax) / norm(a, d, e, Ec, r.lim, b, q, Jb, vmax, rmax)
  ll <- sum(log(liklihood))
  return(-ll)
}

nll.norm.Ec <- Vectorize(function(Ec,xv){nll.norm(xv,Ec=Ec)}, vectorize.args = "Ec")

density.adj <- function(xv, a = 2.0, d = -5.3, e = 2.5, Ec = .16, r.lim = 1.5,
                     b = -9.0, q = 6.9, Jb = 8.6*10^(-2), vmax = 21, rmax = 1.5, adj = 0){
  nll <- nll.norm(xv=xv, a=a, d=d, e=e, Ec=Ec, r.lim=r.lim,
                 q=q, b=b, Jb=Jb, vmax=vmax, rmax=rmax)
  
  return(exp((adj-nll)))
}

density.adj.Ec <- Vectorize(function(Ec,xv,adj){density.adj(xv,Ec=Ec,adj=adj)}, vectorize.args = "Ec")

######################################
#Set parameters
######################################

a <- 2.0
d <- -5.3
e <- 2.5
Ec <- .16
r.lim <- 1.5
b <- -9.0
q <- 6.9
Jb <- 8.6*10^(-2)
vmax <- 21
rmax <- 1.5

phi.s <- (vmax/.465)^2
r.s <- rmax/2.16
xlim <- r.lim / r.s

######################################
#Norm Calculation
######################################

#Ec.grid <- seq(.01,1,.001)
Ec.grid <- unique(c(seq(.01,.31,.1),seq(.11,.121,.01),seq(.133,.182,.001),seq(.148,.168,.0001)))
Ec.grid <- sort(Ec.grid)
# N.grid <- numeric(length=length(Ec.grid))
# 
# for(i in 1:length(Ec.grid)){
#   N.grid[i] <- norm(a=a,q=q,Ec=Ec.grid[i],d=d,r.lim=r.lim,e=e, Jb = 8.6*10^(-2), vmax = 21, rmax = 1.5)
# }
# 
# plot(Ec.grid,N.grid,"l", ylab = "1/N", xlab = "Ec")


#########################################################
#Calculate Posterior
#########################################################

data.E <- read.delim("../data/simulated_galaxy_reduced_model.txt")

r <- sqrt(data.E[,1]^2+data.E[,2]^2+data.E[,3]^2)
x <- r / r.s
v <- sqrt(data.E[,4]^2+data.E[,5]^2+data.E[,6]^2)
xv <- cbind(x,v)

#MLE Adjustment
Ec.mle.full <- optim(.16,nll.norm.Ec,xv=xv, method = 'L-BFGS-B', lower = .01, upper = 1)
denom.int <- integrate(density.adj.Ec,lower=.1, upper=.3,xv=xv, adj = Ec.mle.full$value)$value

Ec.posterior.full <- density.adj.Ec(Ec=Ec.grid, xv=xv, adj = Ec.mle.full$value)/denom.int

data.E <- read.delim("../data/simulated_galaxy_reduced_model_3.txt")

r <- sqrt(data.E[,1]^2+data.E[,2]^2+data.E[,3]^2)
x <- r / r.s
v <- sqrt(data.E[,4]^2+data.E[,5]^2+data.E[,6]^2)
xv <- cbind(x,v)

#MLE Adjustment
Ec.mle.reduced <- optim(.16,nll.norm.Ec,xv=xv, method = 'L-BFGS-B', lower = .01, upper = 1)
denom.int <- integrate(density.adj.Ec,lower=.1, upper=.3,xv=xv, adj = Ec.mle.reduced$value)$value

Ec.posterior.reduced <- density.adj.Ec(Ec=Ec.grid, xv=xv, adj = Ec.mle.reduced$value)/denom.int

save(Ec.grid, Ec.posterior.full, Ec.mle.full,
     Ec.posterior.reduced, Ec.mle.reduced, file="posterior_Ec_results.R")
