# This file runs the SVI algorithm for simulation 2
rm(list = ls())

# read in arguments from the command line
args = commandArgs(trailingOnly=TRUE)
pttn = as.numeric(args[1])
noise = as.numeric(args[2])
size = as.numeric(args[3])


# for testing
# pttn = 1
# noise = 1
# size = 1


library(methods)
library(variationalMira)
library(microbenchmark)
#library(readr)

source("./code/corefun.R")
source("./code/finalInit.R")
source("./params.R")
# Parameters
# Range of frequency search
fmin = 1/1000
fmax = 1/100
nseq = 500
lev = c("I", "K")
nBand = length(lev)


# Data for the simulated dataset
allFiles = list.files(paste0(PARAMS_DATAPATH, 
                             "simu2/output/pttn", pttn,
                             "/noise", noise, 
                             "/size", size), 
                      full.names = TRUE)
nSample = length(allFiles)

if(!dir.exists("./result/simu2/SVI/")){
    dir.create("./result/simu2/SVI/", recursive = T)
}

# where to save the final result
finalsavepath = paste0("./result/simu2/SVI/result_pttn", pttn,
                       "_noise", noise, 
                       "_size", size, ".RData")



## Load the Priors
load("./data/prior2.RData")
set.seed(50)
deltaBar = 100
rbar = 1
nBar = 1

## This function uses the first 200 samples to initialize
## PLR intercept.
initList = initPLR(allFiles, 300, fmin, fmax,
                   nseq, lev,  alphaBarList)
gammaBar = initList$gammaBar
alphaBarList = initList$alphaBarList

kernelBetaList = getInitKernelPara(allFiles, 300, lev, 1e4,
                                   initList$freqAll,
                                   initList$keepS)
# kernelBetaList


print("Initialization finished")

## Construct a hierarchical model from the R package
modelObj = new(hierarchicalPE)

## The model will search over a dense sequency of frequency [fmin, fmax]
## The length of the sequence is `nseq`
modelObj$set_globalFreq(fmin, fmax, nseq)
## Set the number of bands
modelObj$set_nBand(nBand)
## Set hyperparameters
modelObj$set_alpha(alphaBarList, deltaBar)
modelObj$set_gamma(gammaBar, rbar)
modelObj$set_Omega(omegaBar, nBar)
## Set kernel parameter, `1` stands for squared exponential kernel.
modelObj$set_kernel(kernelBetaList, 1)

## Load all the samples
for(i in 1:nSample){
    filePath = allFiles[i]
    lcData = read.table(filePath)
    inputData = with(lcData, convertList(V1, V3, V4, V2, lev))
    modelObj$push_data(inputData[[1]], inputData[[2]], inputData[[3]], inputData[[4]])
}

## Start model iteration with 1000 rounds of iteration. 
## In Algorithm 2, the step size is set by (c_1+t)^{c2},
## for t = 1,2,..., maxIter
maxIter = 1000
c_1 = 2e3
c_2 = 0.8
print("Start SVI")
# microbenchmark::microbenchmark(
modelObj$compute(maxIter, c_1, -c_2)
# times = 1)
print("End SVI")


## 
## The first row: estimated period
## The Second row: sigma uncertainty of period
## The 3rd, 5th row: the mean magnitude for the I K bands respectively.
## The 4th, 6th row: the uncertainty for the mean magnitude for each band. 
allPL = modelObj$get_PLValues()

result = t(allPL[c(1,5),])
eta_alpha1 = modelObj$get_eta_alpha1()
eta_alpha2 = modelObj$get_eta_alpha2()
save(allPL, result, eta_alpha1, eta_alpha2, 
     file = finalsavepath)
#save(allPL, file = finalsavepath)


