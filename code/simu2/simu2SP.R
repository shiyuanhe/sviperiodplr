# Period estimation with the Semi-parametric (SP) model
rm(list = ls())
library(varStar)
library(snowfall)
source("./params.R")

args = commandArgs(trailingOnly=TRUE)
pttn = as.numeric(args[1]) # pattern index
noise = as.numeric(args[2]) # noise index
size = as.numeric(args[3]) # size index

# for testing
# pttn = 3
# noise = 3
# size = 7


# Set the path to the simulated dataset
allFiles = list.files(paste0(PARAMS_DATAPATH, 
                             "/simu2/output/pttn", pttn,
                             "/noise", noise, 
                             "/size", size), 
                      full.names = TRUE)
nSample = length(allFiles)

if(!dir.exists("./result/simu2/SP/")){
    dir.create("./result/simu2/SP/", recursive = T)
}

# The path to save the final resutl
finalsavepath = paste0("./result/simu2/SP/result_pttn", pttn, 
                       "_noise", noise, 
                       "_size", size, ".RData")




fitSingle = function(ii){
    res = c(NA, NA)
    filePath = allFiles[ii]
    lcData = read.table(filePath)
    ## Only K-band Data is used
    lcData = subset(lcData, lcData$V2 == "K")
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

result = do.call(rbind, res)
save(result, file = finalsavepath)

