rm(list = ls())
source("./params.R")


args = commandArgs(trailingOnly=TRUE)
pttn = as.numeric(args[1]) # pattern index
noise = as.numeric(args[2]) # noise index
size = as.numeric(args[3]) # size index

# for testing
# pttn = 3
# noise = 3
# size = 3

# Set the path to the simulated dataset
allFiles = list.files(paste0(PARAMS_DATAPATH, 
                             "simu2/output/pttn", pttn,
                             "/noise", noise, 
                             "/size", size), 
                      full.names = TRUE)
nSample = length(allFiles)

# The path to save the final result
if(!dir.exists("./result/simu2/MSP/")){
    dir.create("./result/simu2/MSP/", recursive = T)
}
finalsavepath = paste0("./result/simu2/MSP/result_pttn", pttn, 
                       "_noise", noise, 
                       "_size", size, ".RData")


MSP_process = function(i){
    if(i %% 5 == 0){
        logPath = paste0("./log/msp_simu2.dat")
        write(i, file = logPath, append = TRUE)
    }
    lcData = allFiles[i]

    ## Run MSP algorithm
    spcPath = paste0("./MSP/gp_spectra/", Sys.getpid())
    cmd = paste("./MSP/mGP", lcData, spcPath, "1")
    system(cmd)
    spc = read.table(paste0(spcPath,".dat"), header = FALSE)
    pHat = 1/spc$V1[1]
    return(pHat)
}

library(snowfall)
sfInit(parallel = TRUE, cpus = PARAMS_NCPUS)
sfExportAll()
res = sfClusterApplyLB(1:nSample, MSP_process)
sfStop()

result = sapply(res, function(tmp) tmp)
result = matrix(result, nrow = nSample, ncol = 1)
save(result, file = finalsavepath)

