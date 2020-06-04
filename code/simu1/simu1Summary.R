# ------- simu1 ---------
source("./params.R")
# read in the true frequency
trueFreq = read.table(paste0(PARAMS_DATAPATH,"/id_map_O.dat"), 
                      header = TRUE, stringsAsFactors = FALSE)


computeLoss = function(pHat, p0){
    fHat = 1/pHat
    f0 = 1/p0
    diffF = abs(fHat - f0)
    l1 = mean( diffF < 0.00027, na.rm = TRUE)
    l2 = mean( diffF, na.rm = TRUE) 
    return(c(l1, l2))
}


load("./result/simu1/GLSresult.RData")
gls = computeLoss(GLSresult[,1], trueFreq$period)

load("./result/simu1/MGLSresult.RData")
mgls = computeLoss(MGLSresult[,1], trueFreq$period)

load("./result/simu1/MSPresult.RData")
msp = computeLoss(MSPresult, trueFreq$period)

load("./result/simu1/SVIresult.RData")
svi = computeLoss(SVIresult[,1], trueFreq$period)

load("./result/simu1/SPresult.RData")
sp = computeLoss(SPresult[,1], trueFreq$period)

rrformat = function(value){
    value = round(value * 100, 2)
    sprintf("%.2f",value)
}

maeformat = function(value){
    value = round(value * 1e4, 2)
    sprintf("%.2f",value)
}


qfList = "./result/simu1/SVIqf/"
starID = gsub(".txt","",list.files(qfList))
fList = list.files(qfList, full.names = T)

computeInterval = function(qf, confLevel){
    fSeqLong = qf[,1]
    qSeqLong = qf[,2]
    qD = sort(qSeqLong, decreasing = T)
    deltaF = fSeqLong[2] - fSeqLong[1]
    for(qValue in qD){
        sel = qSeqLong > qValue
        proV = sum(qSeqLong[sel]) * deltaF
        if(proV>confLevel) break
    }
    confInterval = fSeqLong[sel]
    return(confInterval)    
}

processSingle = function(fI,confLevel){
    freq0 = 1/trueFreq$period[fI]
    qf = read.table(fList[fI], header = F)
    cInt = computeInterval(qf,confLevel)
    res = FALSE
    if(any(freq0 > cInt) && any(freq0 < cInt)){
        c1 = max(which(freq0 > cInt))
        c2 = min(which(freq0 < cInt))
        if(c1+1==c2) res = TRUE
    }
    return(res)
}

# Check the actual coverage of the confidence sets
library(snowfall)
sfInit(parallel = TRUE, cpus = 8)
sfExportAll()
resList90 = sfClusterApplyLB(1:5000, processSingle,confLevel=0.9)
resList95 = sfClusterApplyLB(1:5000, processSingle,confLevel=0.95)
resList99 = sfClusterApplyLB(1:5000, processSingle,confLevel=0.99)
resList995 = sfClusterApplyLB(1:5000, processSingle,confLevel=0.995)
sfStop()

resList90 = Reduce("c",resList90)
cover90 = sum(resList90)/5000
resList95 = Reduce("c",resList95)
cover95 = sum(resList95)/5000
resList99 = Reduce("c",resList99)
cover99 = sum(resList99)/5000
resList995 = Reduce("c",resList995)
cover995 = sum(resList995)/5000

