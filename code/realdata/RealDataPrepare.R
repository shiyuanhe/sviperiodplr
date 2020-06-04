
## This file convert the real data to required data format.
source("./params.R")
inputFolder = "./data/elcs3/"
outputFolder = "./data/mira/realData/"
inputFiles = list.files(inputFolder)

for(ii in 1:length(inputFiles)){
    cf = paste0(inputFolder, inputFiles[ii])
    miraData = read.table(cf, stringsAsFactors = F)
    miraData$V5 = sqrt(miraData$V5^2 + miraData$V6^2)
    miraData = miraData[,c(1,2,4,5)]
    of = paste0(outputFolder, inputFiles[ii],".txt")
    write.table(miraData, file = of, quote = F, col.names = F, row.names = F)
}

