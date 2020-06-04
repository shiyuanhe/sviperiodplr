# convert the original format of Wenlong et al.(2018) to the required data format. 

iP = "./data/simu1Raw"
oP = "./data/simu1/"

filelist = list.files(iP, full.names = T)
starID = list.files(iP, full.names = F)

for(ri in 1:length(filelist)){
    fileI = filelist[ri]
    fileO = paste0(oP, starID[ri])
    tmp = read.table(fileI, skip = 1)
    tmp$V7 = sqrt(tmp$V5^2+tmp$V6^2)
    write.table(tmp[,c("V1","V2","V4","V7")], fileO, col.names = F, row.names = F, quote=F)
}

