###################
#Advanced Data Analysis
#Brendan McVeigh
#October 02, 2015
###################

library(rPython)

python.load("../sampling_functions/functions_0721_E.py")
source("../abc_functions/ABC_function.R")
source("../abc_functions/priors.R")
source("../abc_functions/distances.R")
source("../abc_functions/noise.R")

###################
#Server Test Algorithm
###################

stub <- "dens2d_"
stub <- paste0(stub,"output")
if(!dir.exists(stub)){
  dir.create(stub)
}
stub <- paste0(stub,"/")

data <- read.delim("../data/simulated_galaxy_reduced_model.txt")
suff.stats <- list(python.fun="sample3R",
                   param=c(2.0,-5.3, 2.5, 0.16, 1.5, -9.0, 6.9, 0.086, 21.0, 1.5,1 ,3 ,1),
                   steps=25,
                   samplesize=5000,
                   resamplefactor=1.5,
                   data.true=data,
                   indicies=4)
save(data, suff.stats, file=paste0(stub,"inputs.R"))

results.Ec <- abc_pmc(suff.stats, n = 10000, np = 1000, q = .8, steps = 20,
                      prior.funs = c("prior.Ec"),
                      noise.funs = c("noise.normal.01.1"),
                      distance.fun = dens2d.distance,
                      progress=TRUE,
                      progress.file=paste0(stub,"prog.txt"),
                      write=TRUE,
                      write.file = paste0(stub,"results.R"),
                      .python.load="../sampling_functions/functions_0721_E.py")