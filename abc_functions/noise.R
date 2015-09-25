###################
#Advanced Data Analysis
#Brendan McVeigh
#September 22, 2015
###################

###################
#noise functions
###################

#gaussian noise
noise.normal <- function(mean=0,sd=1,type="s",obs=NA){
  
  #sample from function, will take vector of means
  if (type=="s"){
    n <- length(mean)
    return(rnorm(n,mean,sd))
  }
  
  #determine density of function
  else if(type=="d"){
    return(dnorm(obs,mean,sd))
  }
  else{
    cat("type specified incorrectly")
    return(NA)
  }
}

#sample/density of gamma based on mean and variance
noise.gamma <- function(mean=1,sd=1,type="s",obs=NA){
  
  #sample from function, will take vector of means
  if (type=="s"){
    n <- length(mean)
    return(rgamma(n,shape=mean^2/sd^2, scale=sd^2/mean))  
  }
  #determine density of function
  else if(type=="d"){
    return(dgamma(obs,shape=mean^2/sd^2, scale=sd^2/mean))
  }
  else{
    cat("type specified incorrectly")
    return(NA)
  }
}

#sample/density of gamma based on mean and variance, minium vale of .01
noise.gamma.01 <- function(mean=1,sd=1,type="s",obs=NA){
  
  #sample from function, will take vector of means
  if (type=="s"){
    n <- length(mean)
    mean <- mean - .01
    noise.value <- rgamma(n,shape=mean^2/sd^2, scale=sd^2/mean)
    noise.value <- noise.value + .01
    return(noise.value)  
  }
  #determine density of function
  else if(type=="d"){
    mean <- mean - .01
    obs <- obs - .01
    return(dgamma(obs,shape=mean^2/sd^2, scale=sd^2/mean))
  }
  else{
    cat("type specified incorrectly")
    return(NA)
  }
}

noise.beta.xy <- function(mean=1,sd=1,type="s",obs=NA){  
  width <- abs(x-y)
  mi <- min(x,y)
  ma <- max(x,y)
  #sample from function, will take vector of means
  if (type=="s"){
    n <- length(mean)
    mean <- (mean - mi)/width
    sd <- sd/width
    noise.value <- rbeta(n,shape1=(mean^2-(sd)^2*mean-mean^3)/sd^2,
                         shape2=(1-mean)*(mean-sd^2-mean^2)/sd^2)
    noise.value <- noise.value*width+x
    return(noise.value)  
  }
  #determine density of function
  else if(type=="d"){
    mean <- (mean - mi)/width
    obs <- (obs - mi)/width
    sd <- sd/width
    return(dgamma(obs,shape1=(mean^2-(sd)^2*mean-mean^3)/sd^2,
                  shape2=(1-mean)*(mean-sd^2-mean^2)/sd^2))
  }
  else{
    cat("type specified incorrectly")
    return(NA)
  }
}

noise.beta.3.10 <- function(mean=1,sd=1,type="s",obs=NA){  
  x <- 3
  y <- 10
  width <- abs(x-y)
  mi <- min(x,y)
  ma <- max(x,y)
  #sample from function, will take vector of means
  if (type=="s"){
    n <- length(mean)
    mean <- (mean - mi)/width
    sd <- sd/width
    noise.value <- rbeta(n,shape1=(mean^2-(sd)^2*mean-mean^3)/sd^2,
                         shape2=(1-mean)*(mean-sd^2-mean^2)/sd^2)
    noise.value <- noise.value*width+x
    return(noise.value)  
  }
  #determine density of function
  else if(type=="d"){
    mean <- (mean - mi)/width
    obs <- (obs - mi)/width
    sd <- sd/width
    return(dgamma(obs,shape1=(mean^2-(sd)^2*mean-mean^3)/sd^2,
                  shape2=(1-mean)*(mean-sd^2-mean^2)/sd^2))
  }
  else{
    cat("type specified incorrectly")
    return(NA)
  }
}

#Generates normal noise truncated between .01 and 1
noise.normal.01.1 <- function(mean=0,sd=1,type="s",obs=NA){
  
  noise.normal.01.1.single <- function(mean=0,sd=1,type="s",obs=NA){
    
    #sample from function, will take vector of means
    if (type=="s"){
      proposal <- rnorm(1,mean=mean, sd=sd)
      while(proposal<.01 | proposal > 1){
        proposal <- rnorm(1,mean=mean,sd=sd)
      }
      return(proposal)
    }
    
    #determine density of function
    else if(type=="d"){
      scale.factor <- 1/(pnorm(1,mean,sd)-pnorm(.01,mean,sd))
      return(scale.factor*dnorm(obs,mean,sd))
    }
    else{
      cat("type specified incorrectly")
      return(NA)
    }
  }
  
  n <- length(mean)
  result <- numeric(n)
  for(ii in 1:n){
    result[ii] <- noise.normal.01.1.single(mean[ii],sd,type,obs)
  }
  return(result)
}

#Generates normal noise truncated between 1 and 12
noise.normal.1.12 <- function(mean=0,sd=1,type="s",obs=NA){
  
  noise.normal.1.12.single <- function(mean=0,sd=1,type="s",obs=NA){
    
    #sample from function, will take vector of means
    if (type=="s"){
      proposal <- rnorm(1,mean=mean, sd=sd)
      while(proposal<1 | proposal > 12){
        proposal <- rnorm(1,mean=mean,sd=sd)
      }
      return(proposal)
    }
    
    #determine density of function
    else if(type=="d"){
      scale.factor <- 1/(pnorm(12,mean,sd)-pnorm(1,mean,sd))
      return(scale.factor*dnorm(obs,mean,sd))
    }
    else{
      cat("type specified incorrectly")
      return(NA)
    }
  }
  
  n <- length(mean)
  result <- numeric(n)
  for(ii in 1:n){
    result[ii] <- noise.normal.1.12.single(mean[ii],sd,type,obs)
  }
  return(result)
}
