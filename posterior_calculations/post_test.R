setwd("/Users/Brendan/Google Drive/2015_S2_Fall/ADA/code/posterior_calculations")

a <- 2.0
d <- -5.3
e <- 2.5
Ec <- .16
rlim <- 1.5
b <- -9.0
q <- 6.9
Jb <- 8.6*10^(-2)
vmax <- 21
rmax <- 1.5

phi.s <- (vmax/.465)^2
r.s <- rmax/2.16

Phi <- function(r, vmax=21, rmax = 1.5){
  PhiS <- (vmax/.465)^2
  rs <- rmax / 2.16
  return(PhiS * (1-log(1+r/rs)/(r/rs)))
}

vesc <- function(r, rlim = 1.5, vmax=21, rmax = 1.5){
  return(sqrt(2*(Phi(rlim, vmax=vmax, rmax=rmax) - Phi(r, vmax=vmax, rmax=rmax))))
}

hE <- function(r,v,a = 2.0, d = -5.3, e = 2.5, Ec = .16, rlim = 1.5, b = -9.0,
               q = 6.9, Jb = 8.6*10^(-2), vmax = 21, rmax = 1.5){
  phiLim <- Phi(rlim, vmax=vmax, rmax=rmax)
  E <- v^2/2 + Phi(r)
  if(E < phiLim){
    return(E^a*(E^q+Ec^q)^(d/q)*(phiLim -E)^e)
  } else{
    return(0)
  }
}

dens.inner <- function(v, r, a = 2.0, d = -5.3, e = 2.5, Ec = .16, rlim = 1.5,
                       b = -9.0,q = 6.9, Jb = 8.6*10^(-2), vmax = 21, rmax = 1.5){

  out <- lapply(v, function(vel) hE(r,vel,a=a, d=d, e=e, Ec=Ec, rlim=rlim,q=q, b=b,
                                    Jb=Jb, vmax=vmax, rmax=rmax)*r^2*vel^2)
  return(unlist(out))
}

int.inner <- function(r, a = 2.0, d = -5.3, e = 2.5, Ec = .16, rlim = 1.5,
                      b = -9.0,q = 6.9, Jb = 8.6*10^(-2), vmax = 21, rmax = 1.5){
  low <- 0
  up <- vesc(r, rlim=rlim, vmax=vmax, rmax=rmax)
  out <- integrate(dens.inner, lower = low, upper = up, r=r, a=a, d=d, e=e, Ec=Ec, rlim=rlim,
            q=q, b=b, Jb=Jb, vmax=vmax, rmax=rmax)
  return(out$value)
}

dens.outer <- function(r, a = 2.0, d = -5.3, e = 2.5, Ec = .16, rlim = 1.5,
                       b = -9.0,q = 6.9, Jb = 8.6*10^(-2), vmax = 21, rmax = 1.5){
  out <- lapply(r, function(rad) int.inner(rad, a=a, d=d, e=e, Ec=Ec, rlim=rlim,q=q, b=b,
                                     Jb=Jb, vmax=vmax, rmax=rmax))
  return(unlist(out))
}

int.outer <- function(a = 2.0, d = -5.3, e = 2.5, Ec = .16, rlim = 1.5,
                      b = -9.0,q = 6.9, Jb = 8.6*10^(-2), vmax = 21, rmax = 1.5){
  low <- 0
  up <- rlim
  out <- integrate(dens.outer, lower = low, upper = up, a=a, d=d, e=e, Ec=Ec, rlim=rlim,
                   q=q, b=b, Jb=Jb, vmax=vmax, rmax=rmax)
  return(out$value)
}

Ec.grid <- seq(0,1,length=500)
y <- numeric(500)
for(ii in 1:500){
  y[ii] <- int.outer(Ec = Ec.grid[ii])
}

# a = 2.0, d = -5.3, e = 2.5, Ec = .16, rlim = 1.5,
# b = -9.0,q = 6.9, Jb = 8.6*10^(-2), vmax = 21, rmax = 1.5
# a=a, d=d, e=e, Ec=Ec, rlim=rlim,
# q=q, b=b, Jb=Jb, vmax=vmax, rmax=rmax
