## Summarize Simulation  2 Result 
# ----- prepare ----
options(stringsAsFactors = FALSE)
library(readr)
library(ggplot2)
library(scales)
source("./params.R")

lsfilePath = paste0(PARAMS_DATAPATH, "simu2/output")

computeAccuracy = function(pttn, noise, size, fhat){
    inputfile = paste0(lsfilePath,"/pttn", pttn,
                       "/noise", noise,"/pars/mag_",size,".txt")
    truePL = read_delim(inputfile, col_types = cols(),
                        col_names = TRUE, delim = " ")
    diffF = abs(fhat - truePL$freq)
    accuracy = mean(diffF<0.00027, na.rm = TRUE)
    #msre = sqrt(mean((fhat - truePL$freq)^2, na.rm = TRUE))
    msre = mean(diffF, na.rm = TRUE)
    return(c(accuracy, msre))
}


# computePLRLoss = function(freq, mag){
#     plresidual = mag - alphaTrue[2] * log10(freq) - alphaTrue[3] * log10(freq)^2
#     alpha0Hat = median(plresidual)
#     diff = alpha0Hat - alphaTrue[1]
#     return(abs(diff))
# }

# The following summarizes the accuracy and mean square root error 
# of the SVI, LombScarge and other methods
summarizeResult2 = function(method,pttn, noise){
    tmpRes = matrix(0, 10, 3)
    for(size in 1:10){ #1:10
        try({
            finalsavepath = paste0("./result/simu2/", method, 
                                   "/result_pttn", pttn,
                                   "_noise", noise, 
                                   "_size", size, ".RData")
            load(finalsavepath)
            tmpRes[size, 1:2]  = computeAccuracy(pttn, noise, size, 1/result[,1])
            #tmpRes[size, 3] = computePLRLoss(1/result[,1], result[,2])
        })
    }
    tmpRes = data.frame(
        size = 1:10,
        methods = method, 
        accuracy = tmpRes[,1],
        msre = tmpRes[,2],
        plr = tmpRes[,3])
    return(tmpRes)
}

summarizeLSresult = function(pttn, noise){
    result = matrix(0, 10, 2)
    
    for(size in 1:10){
        initEsti = paste0(lsfilePath, "results/periods_pttn", pttn, "_noise", 
                          noise, "_size", size, ".txt")
        res = read_delim(initEsti, col_names = TRUE, delim = " ", col_types = cols())
        colnames(res) = c("mag", "freq")
        result[size, ]  = computeAccuracy(pttn, noise, size, res$freq)
    }
    result = data.frame(
        size = 1:10,
        methods = "LS", 
        accuracy = result[,1],
        msre = result[,2])
    return(result)
}

summarizeOtherResult = function(pttn, noise, datasetname){
    result = matrix(0, 10, 2)
    ds = get(datasetname)
    for(size in 1:10){
        ename = paste0("p", pttn, "n", noise, "s", size)
        result[size, ]  = computeAccuracy(pttn, noise, size, ds[[ename]])
    }
    if(datasetname=="res_msp") mm = "MSP"
    if(datasetname=="res_sp") mm = "SP"
    if(datasetname=="res_pbmsp_sim2") mm="PBMSP2"
    result = data.frame(
        size = 1:10,
        methods = mm, 
        accuracy = result[,1],
        msre = result[,2])
    return(result)
    
}


# pack all result into one
summarizeAll = function(pttn, noise){
    # res1=  summarizeOtherResult(pttn, noise, "res_sp")
    # res3=  summarizeOtherResult(pttn, noise, "res_pbmsp_sim2")
    # res4=  summarizeLSresult(pttn,noise)
    res = summarizeResult2("GLS_I",pttn,noise)
    res1 = summarizeResult2("MGLS",pttn,noise)
    res2 = summarizeResult2("SP_I",pttn,noise)
    res3 = summarizeResult2("MSP",pttn,noise)
    #res3=  summarizeOtherResult(pttn, noise, "res_msp")
    res4 = summarizeResult2("SVI",pttn,noise)
    #resAll = rbind(res, res1, res2, res3, res4) #res3
    resAll = rbind(res, res1, res2,res3, res4) 
    resAll$pttn = paste0("C",pttn)
    resAll$noise = paste0("N", noise)
    return(resAll)
}


finalRes = data.frame()
for(pttn in 1:3){
    for(noise in 1:3){
        tmp = summarizeAll(pttn, noise)
        finalRes = rbind(finalRes, tmp)
    }
}

size10 = c("(5,5)","(5,10)","(5,20)","(5,30)","(10,10)","(10,20)","(10,30)","(20,20)","(20,30)","(30,30)")
size101 = factor(c("(5,5)","(5,10)","(10,10)","(5,20)","(10,20)","(20,20)","(5,30)","(10,30)","(20,30)","(30,30)"))

#size10 = c("5 5","5 10","5 20","5 30","10 10","10 20","10 30","20 20","20 30","30 30")
reorder = c(1,2,5,3,4,6,7,8,9,10)
finalRes$size = factor(finalRes$size, levels = 1:10, labels  = size10)
finalRes$methods = gsub("_I","", finalRes$methods)
finalRes$size = factor(finalRes$size, levels = size101)

# ---- accuracy ----
if(SAVEPLOT) pdf("./figures/simu2accuracy.pdf", height = 6, width = 6.5)
p = ggplot(finalRes, aes(x = size, y = accuracy, group = methods, color =methods, shape = methods)) + 
    geom_point(size = 2) + geom_line() + theme_bw() +
    facet_grid(pttn~noise) +#, labeller = labeller(.rows = label_both, .cols = label_both)) +
     theme(axis.text.x = element_text(angle = 90, hjust = 1))+ ylab("RR") + 
    xlab(expression(paste("Size (", n[K],", ", n[I], ")"))) + 
    scale_y_continuous(labels = scales::percent) +
    #theme(legend.position = c(0.32, 0.0), legend.direction="horizontal") #+ 
    theme(legend.position = "top", legend.direction="horizontal") #+ 
print(p)
if(SAVEPLOT) dev.off()

# ---- msre ----
if(SAVEPLOT) pdf("./figures/simu2msre.pdf", height = 6, width = 6.5)
p = ggplot(finalRes, aes(x = size, y = msre, group = methods, color =methods, shape = methods)) + 
    geom_point(size = 2) + geom_line() + theme_bw()+
    facet_grid(pttn~noise) +#, labeller = labeller(.rows = label_both, .cols = label_both)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+ ylab("ADE") + 
    xlab(expression(paste("Size (", n[K],", ", n[I], ")"))) + 
    scale_y_continuous(trans = "log", labels = function(x) paste0(sprintf("%.2e\n", x))) +
     theme(legend.position = "top", legend.direction="horizontal") #+ 
    #theme(legend.title = element_blank())
     #guides(shape = guide_legend(override.aes = list(size = 0.5)))
print(p)
if(SAVEPLOT) dev.off()


