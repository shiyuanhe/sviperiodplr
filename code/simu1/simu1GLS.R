## Run the MGLS method on the first simulated dataset.
## Only the I-band data is employed for period estimation. 

rm(list = ls())
library(multiband)
library(snowfall)
source("./params.R")


# Parameters
# Range of frequency search
fmin = 1/1000
fmax = 1/100
nseq = 500
periodGrid = 1/seq(fmin, fmax, length.out = nseq)

# Set the path to the simulated dataset
fpath = paste0(PARAMS_DATAPATH, "simu1/")
allFiles = list.files(fpath, full.names = TRUE)
nSample = length(allFiles)

# The path to save the final resutl
finalsavepath = paste0("./result/simu1/GLSresult.RData")


## prepare data for the package
multibandData = function(lcData){
    bands = unique(lcData$V2)
    result = list()
    for(bb in bands){
        sel = (lcData$V2==bb)
        lcDataSub = subset(lcData, sel)
        result = c(result, list(lcDataSub[,-2]))
    }
    return(result)
}

fitSingle = function(ii){
    filePath = allFiles[ii]
    lcData = read.table(filePath)
    lcData = subset(lcData, lcData$V2=="I")
    inputData = multibandData(lcData)
    result = pgls(inputData, periods = periodGrid, gamma1 = 0, gamma2 = 0)
    #max_iter = 1000, fast_BCD_flag = FALSE)
    #plot(1/result$period_seq_all, result$rss_ls, type ="l")
    optIndex = which.min(result$rss_ls)
    pHat = result$period_seq_all[optIndex]
    mHat = result$best_fitLS[1]
    res = c(pHat, mHat)
    res = matrix(res, nrow = 1)
    return(res)    
}

sfInit(parallel = TRUE, cpus = PARAMS_NCPUS)
sfExportAll()
sfLibrary(multiband)
res = sfClusterApplyLB(1:nSample, fitSingle)
sfStop()

GLSresult = do.call(rbind, res)
save(GLSresult, file = finalsavepath)

