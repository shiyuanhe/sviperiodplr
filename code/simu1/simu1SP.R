# Period estimation with the Semi-parametric (SP) model for the first simuldation
rm(list = ls())
library(varStar)
library(snowfall)
source("./params.R")
lev = c("I", "J", "H", "K")

# Set the path to the simulated dataset
fpath = paste0(PARAMS_DATAPATH, "simu1/")
allFiles = list.files(fpath, full.names = TRUE)
nSample = length(allFiles)

# The path to save the final resutl
finalsavepath = paste0("./result/simu1/SPresult.RData")


fitSingle = function(ii){
    res = c(NA, NA)
    filePath = allFiles[ii]
    lcData = read.table(filePath)
    ## Only I-band Data is used
    lcData = subset(lcData, lcData$V2 == "I")
    priorV = c(mean( lcData$V3), 10,10)
    starObs = new(simple_gpModel, lcData$V1, lcData$V3, lcData$V4, priorV)
    try({
        # freq is searched over 1/1000 to 1/100
        fdelta = (1/100-1/1000)/500
        spc = starObs$freq_est(1/1000,1/100,fdelta)
        #plot(spc[,1], spc[,2], type = "l")
        optIndex = which.max(spc[,2])
        pHat = 1/spc[optIndex, 1]
        thetaHat = spc[optIndex, 3:4]
        starObs$set_theta(thetaHat)
        yHat = starObs$predict(1000)
        mHat = yHat$gamma_bar[1]
        res = c(pHat, mHat)
    })
    res = matrix(res, nrow = 1)
    return(res)    
}

sfInit(parallel = TRUE, cpus = PARAMS_NCPUS)
sfExportAll()
sfLibrary(varStar)
res = sfClusterApplyLB(1:nSample, fitSingle)
sfStop()

SPresult = do.call(rbind, res)
save(SPresult, file = finalsavepath)

