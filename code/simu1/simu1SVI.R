## The SVI method for the first simulation dataset. 

rm(list = ls())

library(microbenchmark)
library(methods)
library(variationalMira)
library(readr)
source("./code/corefun.R")
source("./code/finalInit.R")
source("./params.R")

# Parameters
# Range of frequency search
fmin = 1/1000
fmax = 1/100
nseq = 400
lev = c("I", "J", "H", "K")
nBand = length(lev)

# directory to save computed results
if(!dir.exists("./result/simu1/")){
    dir.create("./result/simu1/", recursive = T)
}
# The path to save the final result
finalsavepath = paste0("./result/simu1/SVIresult.RData")



# Data for the first 5000 Miras
fpath = paste0(PARAMS_DATAPATH, "simu1")
allFiles = list.files(fpath, full.names = TRUE)
afname = list.files(fpath, full.names = FALSE)
nSample = length(allFiles)



## Load the Priors
load("./data/prior1.RData")
# Week prior
deltaBar = 100
nBar = 1
rbar = 1

## Robust estimation of the intercept
initList = initPLR(allFiles, 200, fmin, fmax,
                   nseq, lev,  alphaBarList)
gammaBar = initList$gammaBar # 1/sigma^2 for the PLR
alphaBarList = initList$alphaBarList # List of PLR coefficients

## Estimate kernel parameters
kernelBetaList = getInitKernelPara(allFiles, 200, # only the first 200 samples are used
                                   lev, 1e5,
                                   initList$freqAll, # estimated frequency for each sample
                                   initList$keepS) # outlier inidicator

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
    # Convert to required format.
    # Input for convertList: mjd, mag, sigma, band, lev
    inputData = with(lcData, convertList(V1, V3, V4, V2, lev))
    modelObj$push_data(inputData[[1]], inputData[[2]], inputData[[3]], inputData[[4]])
}


set.seed(50)
## Start model iteration with 1500 rounds of iteration. 
## In Algorithm 2, the step size is set by (c_1+t)^{c2}
## Here c_1=1e3 and c_2 = 0.8
res = microbenchmark::microbenchmark(
    modelObj$compute(1500, 2e3, -0.8),
    times = 1)
print(res)

# For allPL, each column corresponds to one observation. 
# The first 2 rows are the estimated period and its variance  
# The remaining rows are the estimated average mangnitude,
# and its *standard deviation* for each band. lev = c("I", "J", ..)
allPL = modelObj$get_PLValues()

SVIresult = t(allPL[c(1,3,5,7,9),]) # estimated period and mag for IJHK bands
eta_alpha1 = modelObj$get_eta_alpha1()
eta_alpha2 = modelObj$get_eta_alpha2()
#get_eta_Omega1 get_eta_Omega2 get_eta_gamma1 get_eta_gamma2

save(allPL, SVIresult, 
     eta_alpha1, eta_alpha2, 
     file = finalsavepath)


# Store all the posterior densities
qfDir = "./result/simu1/SVIqf/"
if(!dir.exists(qfDir)){
    dir.create(qfDir, recursive = T)
}

for(i in 1:nSample){
    # for the i-th sample, get the final q(f)
    # tmp has two columns. The first column is a dense grid of 
    # frequency over [fmin, fmax]. The second column is the 
    # corresponding value of q(f).
    # Note q(f) is a density function satisfying \int_fmin^fmax q(f) df = 1.
    tmp = modelObj$get_qf(i - 1)
    tmp = round(tmp, 7)
    write.table(tmp, file = paste0(qfDir, afname[i]), 
                row.names = F, col.names = F, quote = FALSE)
}
