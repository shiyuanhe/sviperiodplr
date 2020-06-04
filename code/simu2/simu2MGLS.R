## Run the MGLS method on the second simulated dataset.
rm(list = ls())
library(multiband)
library(snowfall)
source("./params.R")

args = commandArgs(trailingOnly=TRUE)
pttn = as.numeric(args[1]) # pattern index
noise = as.numeric(args[2]) # noise index
size = as.numeric(args[3]) # size index

# for testing
# pttn = 1
# noise = 1
# size = 1


# Parameters
# Range of frequency search
fmin = 1/1000
fmax = 1/100
nseq = 500
periodGrid = 1/seq(fmin, fmax, length.out = nseq)

# Set the path to the simulated dataset
allFiles = list.files(paste0(PARAMS_DATAPATH, 
                             "simu2/output/pttn", pttn,
                             "/noise", noise, 
                             "/size", size), 
                      full.names = TRUE)
nSample = length(allFiles)

if(!dir.exists("./result/simu2/MGLS/")){
    dir.create("./result/simu2/MGLS/", recursive = T)
}

# The path to save the final resutl
finalsavepath = paste0("./result/simu2/MGLS/result_pttn", pttn, 
                       "_noise", noise, 
                       "_size", size, ".RData")



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
    inputData = multibandData(lcData)
    result = pgls(inputData, periods = periodGrid, gamma1 = 0, gamma2 = 0)
    #max_iter = 1000, fast_BCD_flag = FALSE)
    #plot(1/result$period_seq_all, result$rss_ls, type ="l")
    optIndex = which.min(result$rss_ls)
    pHat = result$period_seq_all[optIndex]
    mHat = result$best_fitLS[2,1]
    res = c(pHat, mHat)
    res = matrix(res, nrow = 1)
    return(res)    
}

sfInit(parallel = TRUE, cpus = PARAMS_NCPUS)
sfExportAll()
sfLibrary(multiband)
res = sfClusterApplyLB(1:nSample, fitSingle)
sfStop()

result = do.call(rbind, res)
save(result, file = finalsavepath)

