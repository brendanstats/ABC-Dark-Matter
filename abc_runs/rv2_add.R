###################
#Advanced Data Analysis
#Brendan McVeigh
#October 06, 2015
###################

###################
#Server Test Algorithm
###################

source("../abc_functions/ABC_function.R")
source("../abc_functions/priors.R")
source("../abc_functions/distances.R")
source("../abc_functions/noise.R")

stub <- "rv2_output/"
load(paste0(stub,"inputs.R"))
load(paste0(stub,"results.R"))

results.Ec <- abc_pmc_add(results.prev=results,
                          suff.stats=suff.stats,
                          q=.8,
                          steps=40,
                          prior.funs=c("prior.Ec"),
                          noise.funs=c("noise.normal.01.1"),
                          distance.fun=rv2.distance,
                          parallel=TRUE,
                          progress=TRUE,
                          progress.file=paste0(stub,"progress_add.txt"),
                          write=FALSE,
                          write.file=paste0(stub,"results_add.R"),
                          .python.load="../sampling_functions/functions_0721_E.py")