#############################
#Advanced Data Analysis
#Brendan McVeigh
#May 25, 2015
#############################

###############################################################################
#Define Prior Functions
###############################################################################

#######################################
#Gaussian Example Priors
#######################################

###################
#Prior on Gaussian Mean
###################

prior.mean <- function(x,type){
  
  #random draws
  if(type=="s") {
    return(rnorm(x,0,3))
  }
  #compute density
  else if(type=="d"){
    return(dnorm(x,0,3))
  }
  else{
    cat("type specified incorrectly")
    return(NA)
  }
}

###################
#Prior Gaussian SD
###################

prior.sd <- function(x,type){
  
  #random draws
  if(type=="s") {
    return(rgamma(x,2,3))
  }
  #compute density
  else if(type=="d"){
    return(dgamma(x,2,3))
  }
  else{
    cat("type specified incorrectly")
    return(NA)
  }
}

#######################################
#King Model Priors
#######################################

prior.sigma <- function(x,type){
  
  #random draws
  if(type=="s") {
    return(runif(x, min = 5, max = 10))
  }
  #compute density
  else if(type=="d"){
    return(dunif(x, min = 5, max = 10))
  }
  else{
    cat("type specified incorrectly")
    return(NA)
  }
}

prior.ratio <- function(x,type){
  
  #random draws
  if(type=="s") {
    #return(runif(x, min = 3, max = 10))
    return(runif(x, min = 3, max = 6))
  }
  #compute density
  else if(type=="d"){
    #return(dunif(x, min = 3, max = 10))
    return(dunif(x, min = 3, max = 6))
  }  
  else{
    cat("type specified incorrectly")
    return(NA)
  }
}

#######################################
#Full NFW Priors
#######################################

###################
#Prior on a
###################

prior.a <- function(x,type){
  
  #random draws
  if(type=="s") {
    return(runif(x,min=0,max=10))
  }
  #compute density
  else if(type=="d"){
    return(dunif(x,min=0,max=10))
  }
  else{
    cat("type specified incorrectly")
    return(NA)
  }
}

###################
#Prior on d
###################

prior.d <- function(x,type){
  
  #random draws
  if(type=="s") {
    return(runif(x,min=-10,max=0))
  }
  #compute density
  else if(type=="d"){
    return(dunif(x,min=-10,max=0))
  }
  else{
    cat("type specified incorrectly")
    return(NA)
  }
}

###################
#Prior on e
###################

prior.e <- function(x,type){
  
  #random draws
  if(type=="s") {
    return(runif(x,min=0,max=10))
  }
  #compute density
  else if(type=="d"){
    return(dunif(x,min=0,max=10))
  }
  else{
    cat("type specified incorrectly")
    return(NA)
  }
}

###################
#Prior on Ec
###################

prior.Ec <- function(x,type){
  
  #random draws
  if(type=="s") {
    return(runif(x,min=.01,max=1))
  }
  #compute density
  else if(type=="d"){
    return(dunif(x,min=.01,max=1))
  }
  else{
    cat("type specified incorrectly")
    return(NA)
  }
}

###################
#Prior on rlim
###################

prior.rlim <- function(x,type){
  
  #random draws
  if(type=="s") {
    return(runif(x,min=.5,max=12))
  }
  #compute density
  else if(type=="d"){
    return(dunif(x,min=.5,max=12))
  }
  else{
    cat("type specified incorrectly")
    return(NA)
  }
}

prior.rlim2 <- function(x,type){
  
  #random draws
  if(type=="s") {
    return(rgamma(x,1.5,.05))
  }
  #compute density
  else if(type=="d"){
    return(rgamma(x,1.5,.05))
  }
  else{
    cat("type specified incorrectly")
    return(NA)
  }
}

###################
#Prior on b
###################

prior.b <- function(x,type){
  
  #random draws
  if(type=="s") {
    return(runif(x,min=-15,max=0))
  }
  #compute density
  else if(type=="d"){
    return(dunif(x,min=-15,max=0))
  }
  else{
    cat("type specified incorrectly")
    return(NA)
  }
}

###################
#Prior on q
###################

###################
#Prior on Jb
###################

###################
#Prior on vmax
###################

###################
#Prior on rmax
###################