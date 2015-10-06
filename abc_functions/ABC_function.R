###################
#Advanced Data Analysis
#Brendan McVeigh
#March 30, 2015
###################

library(parallel)
require(doParallel) 
require(foreach)

###################
#call a function with parameters
###################

#function to call functions based on character
call.funs <- function(fun,...){
  return(match.fun(fun)(...))
}

#calculated weight variance
weighted.var <- function(x,w){
  return(sum(w*(x-sum(x*w))^2)/sum(w))
}

#calculated magnitude of a vector
magnitude <- function(x){
  return(sqrt(sum(x^2)))
}

#scale errors based on minimum and sd
scale.errors <- function(x){
  scale.x <- (x-min(x))/sd(x)
  return(scale.x)
}

###################
#add noise
###################
add.noise <- function(thetas,sigma,noise.funs){
  results <- list()
  
  #if just one parameter
  if (dim(thetas)[2]==1){
    results <- as.numeric(match.fun(noise.funs)(mean = thetas, sd = sigma, type = "s"))
  }
  else { 
  	#for multiple parameters add noise to each parameter individually
    for (i in 1:dim(thetas)[2]){
      results[[i]]<- match.fun(noise.funs[i])(mean = thetas[,i], sd = sigma[i], type = "s")
    }
    results <- do.call(cbind,results)
  }
  return(results)
}

###################
#function to calculate weights
###################
update.weight <- function(theta.t, thetas, sigma, w, noise.funs, prior.funs) {
  
  param <- length(theta.t)
  denom.densities <- list()
  prior.densities <- list()
  
  for (i in 1:param){
    denom.densities[[i]] <- match.fun(noise.funs[i])(mean=thetas[,i], sd = sigma[i], type = "d",obs = theta.t[i])
    prior.densities[[i]] <- match.fun(prior.funs[i])(theta.t[i], type = "d")
  }
  denom.densities <- do.call(cbind,denom.densities)
  w.t <- do.call(prod,prior.densities)/sum(apply(denom.densities,1,prod)*w)
  return(w.t)
}

###################
#Function to sample until one acceptable parameter is found
###################
find.one <- function(thetas,w,sigma,epsilon,suff.stats,noise.funs,distance.fun){
  
  cat("find 1", "\n")
  if(class(thetas)=="numeric") thetas <- matrix(thetas,ncol=1)
  np <- dim(thetas)[1]
  
  #set counter values
  accept <- TRUE
  sampled <- 0
  
  #create empty lists
  theta.sampled <- list()
  theta.tested <- list()
  error.tested <- list()
  
  cat("find 2", "\n")
  while(accept){
    
    sampled <- sampled + 1
    
    #sample new value
    theta <- thetas[sample(np,size=1,prob=w),]
    theta <- matrix(theta,nrow=1)
    theta.sampled[[sampled]] <- theta
    cat("Test:", theta.sampled[[sampled]], "\n")
    
    #add noise
    theta <- add.noise(theta,sigma,noise.funs)
    theta.tested[[sampled]] <- theta
    cat("With Noise:", theta.sampled[[sampled]], "\n")
    
    #make into array
    error <- distance.fun(theta,suff.stats)
    error.tested[[sampled]] <- error
    cat("Distance:",error, "Cutoff:", epsilon, "\n")
    #check if parameter is accepted
    if(all(error<=epsilon)){
      accept <- FALSE
    }
  }
  accepted <- rep(FALSE,sampled)
  accepted[sampled] <- TRUE
  all.tested <- data.frame(theta.sampled = do.call(rbind,theta.sampled),
                           theta.tested = do.call(rbind,theta.tested), 
                           error.tested = do.call(rbind,error.tested),
                           accepted = accepted)
  #return parameter
  return(list(theta=theta,error=error,sampled=sampled, all.tested = all.tested))
}

###################
#function to iterate
###################
pmc.step <- function(thetas,w,epsilon,suff.stats,noise.funs,prior.funs,distance.fun,.python.load){
	#convert to matrix if necesary
  cat("pmc 1", "\n")
  if(class(thetas)=="numeric") thetas <- matrix(thetas,ncol=1)
  
  #get number of particles
  np <- dim(thetas)[1]
  stopifnot(np==length(w))
  
  #calculate sample variance for parameter posteriors
  cat("pmc 2", "\n")
  sigma <- sqrt(2*apply(thetas,2,weighted.var,w))
  
  #remove values with 0 weight
  thetas <- thetas[w>0,]
  w <- w[w>0]
  
  #sample np values
  cat("pmc 3", "\n")
  export <- c(noise.funs,"magnitude","calc.m2","find.one","add.noise")
  ncore <- detectCores()
  cl <- makeCluster(ncore)
  registerDoParallel(cl)
  foreach.out <- foreach(i=1:np, .packages=c("rPython","energy"),.export=export) %dopar% {

    for (j in 1:length(.python.load)){
      python.load(.python.load[j])
    }
    find.one(thetas,w,sigma,epsilon,suff.stats,noise.funs,distance.fun)
  }
    
  stopCluster(cl)
  cat("pmc 4", "\n")
  
  #Extract results from loop
  cat("pmc 4.1", "\n")
  thetas.t <- do.call(rbind,lapply(1:np,function(x) foreach.out[[x]][["theta"]]))
  cat("pmc 4.2", "\n")
  error.t <- do.call(rbind,lapply(1:np,function(x) foreach.out[[x]][["error"]]))
  cat("pmc 4.3", "\n")
  total <- do.call(rbind,lapply(1:np,function(x) foreach.out[[x]][["sampled"]]))
  cat("pmc 4.3", "\n")
  #save(foreach.out, file = "save_all.R")
  all.tested <- do.call(rbind,lapply(foreach.out, function(x) x[["all.tested"]]) )
	
  #update weights
  cat("pmc 4.4", "\n")
  if(class(thetas)=="numeric") thetas <- matrix(thetas,ncol=1)
  cat("pmc 4.5", "\n")
  if(class(thetas.t)=="numeric") thetas.t <- matrix(thetas.t,ncol=1) #changed 5/14/2015, might cause issues

  cat("pmc 4.6", "\n")
  wt <- apply(thetas.t,1,update.weight,thetas,sigma,w,noise.funs,prior.funs)
  cat("pmc 4.7", "\n")
  wt <- wt / sum(wt)
  cat("pmc 5", "\n")
  #return values
  return(list(parameters=thetas.t,weights=wt,
              distances=error.t,
              epsilon=epsilon,total.samples=total,finish.time=Sys.time(),
              all.tested=all.tested
  ))
}
######################################
#ABC function
######################################

abc_pmc <- function(suff.stats,n,np,q,steps,prior.funs,noise.funs,distance.fun,
                    parallel=TRUE,
                    progress=TRUE,progress.file="progress.txt",
                    write=FALSE,write.file="output.R",
                    .python.load="functions.py"){
  
  start.time <- Sys.time()
  
  ######################################
  #Output Session Info if recording progress
  ######################################
  if(progress){
    cat("--------------------","\n",file=progress.file,append=FALSE)
    cat("Session Info","\n",file=progress.file,append=TRUE)
    cat("--------------------","\n",file=progress.file,append=TRUE)
    cat(paste("Start Time:",as.character(start.time)),"\n",file=progress.file,append=TRUE)
    cat("Initial Samples:",n,"\n",file=progress.file,append=TRUE)
    cat("Required Samples:",np,"\n",file=progress.file,append=TRUE)
    cat("Shrinkage Rate:",q,"\n",file=progress.file,append=TRUE)
    cat("Time Steps:",steps,"\n",file=progress.file,append=TRUE)
    cat("Prior Functions:",prior.funs,"\n",file=progress.file,append=TRUE)
    cat("Noise Functions:",noise.funs,"\n",file=progress.file,append=TRUE)
    cat("Distance Function:",as.character(substitute(distance.fun)),"\n",file=progress.file,append=TRUE)
    cat("Progress File:",progress.file,"\n",file=progress.file,append=TRUE)
    cat("Output File:",write.file,"\n",file=progress.file,append=TRUE)
    cat("Python Load Call:",.python.load,"\n",file=progress.file,append=TRUE)
    cat("--------------------","\n","\n",file=progress.file,append=TRUE)
  }
  
  ######################################
  #Time Step 1
  ######################################
  cat("abc 1", "\n")
  
  stopifnot(length(prior.funs)==length(noise.funs))
  num.param <- length(prior.funs)
  
  prior.draw <- sapply(prior.funs,call.funs,n,type="s")
  
  #set up parallel nodes
  ncore <- detectCores()
  cl <- makeCluster(ncore)
  registerDoParallel(cl)
  
  error <- foreach(i=1:n,.combine=rbind,.packages=c("rPython","energy"),.export=c("magnitude","calc.m2")) %dopar% {
    for (i in 1:length(.python.load)){
      python.load(.python.load[i])
    }
    distance.fun(prior.draw[i,],suff.stats)
  }
  
  #close clusters
  stopCluster(cl)
  
  all.tested <- data.frame(theta.sampled = prior.draw, theta.tested = prior.draw, error.tested = error)
  
  cat("abc 2", "\n")
  
  if (class(error)=="numeric") error <- matrix(error,ncol=1)
  if(dim(error)[2]==n) error <- t(error)
  
  #rank using scaled errors
  #error.s <- apply(error,2,scale.errors)
  #ord <- order(apply(error.s,1,magnitude))[1:np]

  #simple ranking of errors  
  ranks <- apply(error,2,rank)
  ord <- order(apply(ranks,1,sum))[1:np]
  
  error <- error[ord,]
  if (class(error)=="numeric") error <- matrix(error,ncol=1)
  epsilon <- apply(error,2,max)
  
  cat("abc 3", "\n")
  
  #values at t=1
  thetas <- prior.draw[ord,]
  w <- rep(1/np,np)
  
  results <- list(list(parameters=thetas,
                       weights=w,
                       distances=error,
                       epsilon=epsilon,
                       total.samples=n,
                       finish.time=Sys.time(),
                       prior.draws=prior.draw,
                       all.tested))
  
  cat("abc 4", "\n")
  cat("error type", class(error), "\n")  
  cat("error dim", dim(error), "\n") 
  epsilon <- apply(error,2,quantile,q)
  stopifnot(epsilon>0)
  
  cat("abc 5", "\n")
  
  ######################################
  #Output Session Info if recording progress
  ######################################
  if(progress){
    step.time <- results[[1]][["finish.time"]]-start.time
    cat("--------------------","\n",file=progress.file,append=TRUE)
    cat("Time step", 1, "\n",file=progress.file,append=TRUE)
    cat("--------------------","\n",file=progress.file,append=TRUE)
    cat(paste("Time:",as.character(results[[1]][["finish.time"]])),"\n",file=progress.file,append=TRUE)
    cat(paste("Total Time:",as.character(step.time)),"\n",file=progress.file,append=TRUE)
    cat(paste("Step Time:",as.character(step.time)),"\n",file=progress.file,append=TRUE)
    cat("Epsilon:",results[[1]][["epsilon"]],"\n",file=progress.file,append=TRUE)
    cat("Drawn Samples:",results[[1]][["total.samples"]],"\n",file=progress.file,append=TRUE)
    cat("Accepted Samples:",np,"\n",file=progress.file,append=TRUE)
    cat("--------------------","\n","\n",file=progress.file,append=TRUE)
  }
  cat("Time Step",1,"\n")
  if(write){
    save(results,file=write.file)
  }
######################################
#Time Steps
######################################
  
  for(i in 2:steps){
    results[[i]] <- pmc.step(thetas,w,epsilon,suff.stats,noise.funs,prior.funs,distance.fun,.python.load)
    cat("abc 7", "\n")
    
    ######################################
    #Output Session Info if recording progress
    ######################################
    
    if(progress){
      total.time <- results[[i]][["finish.time"]]-start.time
      step.time <- results[[i]][["finish.time"]]-results[[i-1]][["finish.time"]]
      cat("--------------------","\n",file=progress.file,append=TRUE)
      cat("Time step", i, "\n",file=progress.file,append=TRUE)
      cat("--------------------","\n",file=progress.file,append=TRUE)
      cat(paste("Time:",as.character(results[[i]][["finish.time"]])),"\n",file=progress.file,append=TRUE)
      cat(paste("Total Time:",as.character(total.time)),"\n",file=progress.file,append=TRUE)
      cat(paste("Step Time:",as.character(step.time)),"\n",file=progress.file,append=TRUE)
      cat("Epsilon:",results[[i]][["epsilon"]],"\n",file=progress.file,append=TRUE)
      cat("Drawn Samples:",sum(results[[i]][["total.samples"]]),"\n",file=progress.file,append=TRUE)
      cat("Accepted Samples:",np,"\n",file=progress.file,append=TRUE)
      cat("--------------------","\n","\n",file=progress.file,append=TRUE)
    }
    if(write){
      save(results,file=write.file)
    }
    thetas <- results[[i]][["parameters"]]
    w <- results[[i]][["weights"]]
    error <- results[[i]][["distances"]]
        
    epsilon <- apply(error,2,quantile,q)
    stopifnot(epsilon>0)
    
  }
  return(results)
}

abc_pmc_add <- function(results.prev,
                        suff.stats,
                        q,
                        steps,
                        prior.funs,
                        noise.funs,
                        distance.fun,
                        parallel=TRUE,
                        progress=TRUE,
                        progress.file="progress.txt",
                        write=FALSE,
                        write.file="output.R",
                        .python.load){
  
  steps.complete <- length(results)
  stopifnot(steps.complete < steps)
  
  start.time <- Sys.time()
  
  thetas <- results[[steps.complete]]$parameters
  w <- results[[steps.complete]]$weights
  error <- results[[steps.complete]]$distances
  epsilon <- apply(error,2,quantile,q)
  np = length(thetas)
  
  ######################################
  #Output Session Info if recording progress
  ######################################
  if(progress){
    cat("--------------------","\n",file=progress.file,append=TRUE)
    cat("Session Info","\n",file=progress.file,append=TRUE)
    cat("--------------------","\n",file=progress.file,append=TRUE)
    cat(paste("Start Time:",as.character(start.time)),"\n",file=progress.file,append=TRUE)
    cat("Required Samples:",np,"\n",file=progress.file,append=TRUE)
    cat("Shrinkage Rate:",q,"\n",file=progress.file,append=TRUE)
    cat("Complete Time Steps:",steps.complete,"\n",file=progress.file,append=TRUE)
    cat("Time Steps:",steps,"\n",file=progress.file,append=TRUE)
    cat("Prior Functions:",prior.funs,"\n",file=progress.file,append=TRUE)
    cat("Noise Functions:",noise.funs,"\n",file=progress.file,append=TRUE)
    cat("Distance Function:",as.character(substitute(distance.fun)),"\n",file=progress.file,append=TRUE)
    cat("Progress File:",progress.file,"\n",file=progress.file,append=TRUE)
    cat("Output File:",write.file,"\n",file=progress.file,append=TRUE)
    cat("Python Load Call:",.python.load,"\n",file=progress.file,append=TRUE)
    cat("--------------------","\n","\n",file=progress.file,append=TRUE)
  }
  
  for(i in (steps.complete+1):steps){
    results[[i]] <- pmc.step(thetas,w,epsilon,suff.stats,noise.funs,prior.funs,distance.fun,.python.load)
    cat("abc 7", "\n")
    
    ######################################
    #Output Session Info if recording progress
    ######################################
    
    if(progress){
      total.time <- results[[i]][["finish.time"]]-start.time
      step.time <- results[[i]][["finish.time"]]-results[[i-1]][["finish.time"]]
      cat("--------------------","\n",file=progress.file,append=TRUE)
      cat("Time step", i, "\n",file=progress.file,append=TRUE)
      cat("--------------------","\n",file=progress.file,append=TRUE)
      cat(paste("Time:",as.character(results[[i]][["finish.time"]])),"\n",file=progress.file,append=TRUE)
      cat(paste("Total Time:",as.character(total.time)),"\n",file=progress.file,append=TRUE)
      cat(paste("Step Time:",as.character(step.time)),"\n",file=progress.file,append=TRUE)
      cat("Epsilon:",results[[i]][["epsilon"]],"\n",file=progress.file,append=TRUE)
      cat("Drawn Samples:",results[[i]][["total.samples"]],"\n",file=progress.file,append=TRUE)
      cat("Accepted Samples:",np,"\n",file=progress.file,append=TRUE)
      cat("--------------------","\n","\n",file=progress.file,append=TRUE)
    }
    if(write){
      save(results,file=write.file)
    }
    thetas <- results[[i]][["parameters"]]
    w <- results[[i]][["weights"]]
    error <- results[[i]][["distances"]]
    
    epsilon <- apply(error,2,quantile,q)
    stopifnot(epsilon>0)
    
  }
  return(results)
}
