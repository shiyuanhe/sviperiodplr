## Compute the hyper parameters for all simulations and real data 
## analysis

alphaBarList = list(
    c(13.7989511, -1.8358380, -0.6852287 ), # I-band PLR
    c(12.5734, -3.594, -1.54 ), # J-band PLR
    c(2.234, -12.05, -3.4), # H-band PLR
    c(8.4673, -6.488, -2.23 ) # K-band PLR
)

## Initial Estimation for Omega
templateFiles = list.files("./data/LMC/LMC_templates", full.names = TRUE)
catlog = read.table("./data/LMC/catalog.dat", header = TRUE)

# fit sinusoid curve and return coef
fitSinusoid = function(mjd, mag, period){
    tt = mjd / period * 2 * pi
    xMat = cbind(cos(tt), sin(tt))
    lmCoef = coef(lm(mag~xMat))
    
    tt = seq(1000,5000, length.out = 5000) / period * 2 * pi
    xMatNew = cbind(1,cos(tt), sin(tt))
    yy = xMatNew %*% lmCoef
    #lines(seq(1000,5000, length.out = 5000) , yy)
    return(lmCoef[2:3])
}

# get  coef for all template
betaAll = matrix(0, length(templateFiles), 8)
for(i in 1:length(templateFiles)){
    period = catlog$period[i]
    starObs = read.table(templateFiles[[i]])
    beta = rep(0,8)
    #plot(starObs$V1, starObs$V2, pch =20, col = "grey")
    beta[1:2] = fitSinusoid(starObs$V1, starObs$V2, period)
    #plot(starObs$V1, starObs$V3, pch =20, col = "grey")
    beta[3:4] = fitSinusoid(starObs$V1, starObs$V3, period)
    beta[5:6] = fitSinusoid(starObs$V1, starObs$V4, period)
    beta[7:8] = fitSinusoid(starObs$V1, starObs$V5, period)
    betaAll[i,] = beta
}
# pairs(betaAll)

SigmaMat = cov(betaAll)
omegaBar = solve(SigmaMat) / 1000 # Make the prior weak
save(alphaBarList, omegaBar, file = "./data/prior1.RData")

## Extract the I/J band information for the second simulation
alphaBarList = alphaBarList[c(1,4)]
sel = c(1,2,7,8)
omegaBar = solve(SigmaMat[sel, sel]) #omegaBar[sel, sel]
save(alphaBarList, omegaBar, file = "./data/prior2.RData")


#library(tidyverse)
#library(MASS)
# 
# 
# PLRRobustFit = function(period, aveMag, sigma){
#     datalm = data.frame(logp = log10(period), avemag = aveMag)
#     lmfit = rlm(avemag~logp + I(logp^2), datalm, weights = 1/sigma^2)
#     newData = data.frame(logp = log10(seq(100,1000,length.out = 100)))
#     ySeq = predict(lmfit, newData)
#     newData$ySeq = ySeq
#     return(newData)
# }
# 
# lmcData = read_csv("./data/LMC/LMCData.csv", col_names = FALSE)
# lmcDataSub = filter(lmcData, X2=="0" & X3=="O")
# ## J-band PLR
# plot(log10(lmcDataSub$X16), lmcDataSub$X4, pch = 20)
# lines(PLRRobustFit(lmcDataSub$X16, lmcDataSub$X4, lmcDataSub$X7), col = "blue")
# lines(log10(pSeq), LMC_PLR(pSeq, "J"), col = "red")
# ## H-band PLR
# plot(log10(lmcDataSub$X16), lmcDataSub$X8, pch = 20)
# ## K-band PLR
# plot(log10(lmcDataSub$X16), lmcDataSub$X12, pch = 20)


