## Apply SVI to the real dataset. 

rm(list = ls())


library(microbenchmark)
library(methods)
library(variationalMira)
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


# Data for the first 5000 Miras
dataPath = paste0(PARAMS_DATAPATH, "realData")
allFiles = list.files(dataPath, full.names = TRUE)
allFilesName = list.files(dataPath, full.names = FALSE)
nSample = length(allFiles)

# The path to save the final resutl
if(!dir.exists("./data/realData/")){
    dir.create("./data/realData/", recursive = TRUE)
}
finalsavepath = paste0("./result/realData/realData1265.RData")


## Load the Priors
load("./data/prior1.RData")

# Week prior
deltaBar = 100
nBar = 1
rbar = 1

## This function uses the first 200 samples to initialize
## PLR intercept.
initList = initPLR(allFiles, 200, fmin, fmax,
                   nseq, lev,  alphaBarList)
gammaBar = initList$gammaBar
alphaBarList = initList$alphaBarList



kernelBetaList = getInitKernelPara(allFiles, 200, lev, 1e4,
                                   initList$freqAll,
                                   initList$keepS)



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

set.seed(50)
## Start model iteration with 1000 rounds of iteration. 
## In Algorithm 2, the step size is set by (c_1+t)^{-c2}
## Here c_1=2e3 and c_2 = 0.8
res = microbenchmark::microbenchmark(
    modelObj$compute(1500, 2e3, -0.8),
    times = 1)
print(res)
allPL = modelObj$get_PLValues()
allPL = t(allPL)
colnames(allPL) = c("P", "PSigma", "I", "ISigma",
                    "J", "JSigma", "H", "HSigma", "K", "KSigma")
allTheta = modelObj$get_ThetaAll()

# plot(log(allPL[1,]), allPL[9,], pch = 20, cex = 0.2)
# plot(log(allPL[1,]), allPL[7,], pch = 20, cex = 0.2)
# plot(log(allPL[1,]), allPL[5,], pch = 20, cex = 0.2)
# plot(log(allPL[1,]), allPL[3,], pch = 20, cex = 0.2)

#SVIresult = t(allPL[c(1,3,5,7,9),])
eta_alpha1 = modelObj$get_eta_alpha1()
eta_alpha2 = modelObj$get_eta_alpha2()

eta_Omega1 = modelObj$get_eta_Omega1()
eta_Omega2 = modelObj$get_eta_Omega2()


eta_gamma1 = modelObj$get_eta_gamma1()
eta_gamma2 = modelObj$get_eta_gamma2()

## Compute the corrected magnitude in flux. 
## Compute the uncertainty of log P
fluxCorrectedPL = matrix(0, nSample, 10)
for(starI in 1:nSample){
    filePath = allFiles[starI]
    periodI = allPL[starI, 1]
    
    lcData = read.table(filePath)
    tmin = min(lcData$V1)
    tmax = max(lcData$V1)
    tmax = tmin + ceiling( (tmax - tmin) / periodI) * periodI +  2*periodI
    tmin = tmin - 2* periodI
    tseq = seq(tmin, tmax, by = 1)
    yseq = modelObj$get_fit(starI-1, tseq)
    tmp = modelObj$get_qf(starI - 1)
    # write.table(tmp, file = paste0("~/Data/m33Qf/", allFilesName[starI]),
    #             col.names = F, row.names = F, quote = F)
    
    logPSeq = log10(1/tmp[,1])
    probSeq = tmp[,2]
    logP2M = sum(logPSeq^2 * probSeq) / sum(probSeq)
    logP1M = sum(logPSeq * probSeq) / sum(probSeq)
    logPSigma = sqrt(logP2M - logP1M^2)
    fluxMean = colMeans(10^(-0.4 * yseq))
    magMean = -2.5 * log10(fluxMean)
    fluxCorrectedPL[starI, ] = c(log10(periodI), logPSigma,  
                               magMean[1], allPL[starI, 4],
                               magMean[2], allPL[starI, 6],
                               magMean[3], allPL[starI, 8],
                               magMean[4], allPL[starI, 10])   
}
fluxCorrectedPL = round(fluxCorrectedPL,5)
fluxCorrectedPL[fluxCorrectedPL==0] = NA
fluxCorrectedPL = data.frame(fluxCorrectedPL)
fluxCorrectedPL$ID = gsub(".txt", "", allFilesName)
fluxCorrectedPL = fluxCorrectedPL[,c(11,1:10)]

colnames(fluxCorrectedPL) = c("ID", "logP", "logPSigma", "I", "ISigma",
                            "J", "JSigma", "H", "HSigma", "K", "KSigma")


save(allPL, fluxCorrectedPL,
     allTheta,
     eta_alpha1, eta_alpha2, 
     eta_Omega1, eta_Omega2,
     eta_gamma1, eta_gamma2,
     file = finalsavepath)

