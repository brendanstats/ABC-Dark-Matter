#####################
#Advanced Data Analysis
#Brendan McVeigh
#September 12, 2015
#####################

setwd("~/Google Drive/2015_S2_Fall/ADA/code/sim_galaxies")
library(rPython)
library(parallel)
library(doParallel) 
library(foreach)

python.load("../sampling_functions/functions_0721.py")
source("../abc_functions/priors.R")

##########################################
#Generate Datasets
##########################################

generate.sample <- function(python.file,sample.fun,param,ct,prior){
  
  python.load(python.file)
  
  proposal <- unlist(lapply(prior,function(x) match.fun(x)(1,"s")))
	param[4:5] <- proposal
  param <- as.numeric(param)
  
	sample <- do.call(cbind,python.call(sample.fun,param,25,ct,1.5))
	out.pair <- list(theta=proposal,sample=sample)
	return(out.pair)
}

file.name <- "Ec_rlim_samples3.R"
n <- 5000
ct <- 5000
param <- c(2.0,-5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5, 1, 3, 1)

ncore <- detectCores()
cl <- makeCluster(ncore)
registerDoParallel(cl)
sample <- foreach(i=1:n, .packages="rPython", .export = c("prior.Ec","prior.rlim")) %dopar% {
  source("../abc_functions/priors.R")
  generate.sample(python.file = "../sampling_functions/functions_0721.py", 
                  sample.fun = "sample3R", param = param,
                  ct = ct, prior = c("prior.Ec","prior.rlim"))
}
stopCluster(cl)

save(sample,file = file.name)

generate.matrix <- function(obs){
  n <- length(obs[["times"]])
  return(data.frame(rlim=rep(obs[["rlim"]],n),
                    Ec=rep(obs[["Ec"]],n),
                    steps=rep(obs[["steps"]],n),
                    time=obs[["times"]]))
}

sample.times <- do.call(rbind,lapply(out,generate.matrix))

summary(lm(time~steps+rlim+steps*rlim, data = sample.times))
summary(lm(time~steps+Ec+rlim+steps*Ec+steps*rlim+Ec*rlim+Ec*rlim*steps, data = sample.times))

anova(lm(time~steps+Ec+rlim+steps*Ec+steps*rlim+Ec*rlim+Ec*rlim*steps, data = sample.times))
anova(lm(time~steps+Ec+rlim+steps*Ec+steps*rlim+Ec*rlim+Ec*rlim*steps, data = sample.times),
      lm(time~steps+Ec+rlim+steps*Ec+steps*rlim+Ec*rlim, data = sample.times))

plot(sample.times$steps,sample.times$time)
plot(sample.times$Ec,sample.times$time)
plot(sample.times$rlim,sample.times$time)
