#####################
#Advanced Data Analysis
#Brendan McVeigh
#September 12, 2015
#####################

setwd("~/Google Drive/2015_S2_Fall/ADA/code")
library(rPython)
library(parallel)
library(doParallel) 
library(foreach)

python.load("functions_0624.py")
source("abc_functions/priors.R")

##########################################
#Generate Datasets
##########################################

generate.sample <- function(python.file,sample.fun,param,ct,prior){
  
  python.load(python.file)
  
	proposal <- prior(1,"s")
	param[5] <- proposal
  param <- as.numeric(param)
  
	sample <- do.call(cbind,python.call(sample.fun,param,ct))
	out.pair <- list(theta=proposal,sample=sample)
	return(out.pair)
}

file.name <- "samples_rlim.R"
n <- 2000
ct <- 5000
param <- c(2.0,-5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5)

ncore <- detectCores()
cl <- makeCluster(ncore)
registerDoParallel(cl)
sample <- foreach(i=1:n, .packages="rPython") %dopar%
  generate.sample(python.file = "functions_0624.py", 
  sample.fun = "sampleR",
  param = param, ct = ct, prior = prior.rlim)
stopCluster(cl)

save(sample,file = file.name)
