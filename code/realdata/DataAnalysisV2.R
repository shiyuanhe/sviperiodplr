## Plot the final figures and Tables for the real data analysis
## This file only use 1265 O-rich Miras classified by Wenlong et al. (2018).



# ------- distInit -------
source("./params.R")
## After running SVI on M33 real data, 
source("./code/realdata/realdata_fitPLR.R")
load("./result/realData/realData1265.RData")
allPL = data.frame(allPL)
bandName = c("I", "J", "H", "K")
ylabVec = c(expression(paste(italic(I), "  [mag]")),
            expression(paste(italic(J), "  [mag]")),
            expression(paste(italic(H), "  [mag]")),
            expression(paste(italic(K)[s], "  [mag]")))

yuan2018Table = read.table(paste0(PARAMS_DATAPATH,"/yuan2018.dat"),
                         stringsAsFactors = F, header = TRUE)
sel = yuan2018Table$cl != "C"
selectedID = yuan2018Table$id[sel]
p1 = yuan2018Table$P[sel]


af = list.files(paste0(PARAMS_DATAPATH, "realData"))
af = gsub(".txt", "", af)
aPos = match(af, yuan2018Table$id)
yuan2018Table = yuan2018Table[aPos, ]
#yuan2018Table = yuan2018Table[-1500,]
msp_period = yuan2018Table$P
msp_Jmean = yuan2018Table$J
msp_Hmean = yuan2018Table$H
msp_Kmean = yuan2018Table$K
svi_period = allPL$P
diff = abs(svi_period -  msp_period)/ svi_period
IncP = diff>0.1
#plot(svi_period, msp_period)

## ---- betaCorr -----
allThetaNA = t(allTheta[c(2,3,5,6,8,9,11,12),])
allThetaNA[allThetaNA==0] = NA
allThetaNA = data.frame(allThetaNA)
lbAll = c(expression(beta[italic(I1)]),
          expression(beta[italic(I2)]),
          expression(beta[italic(J1)]),
          expression(beta[italic(J2)]),
          expression(beta[italic(H1)]),
          expression(beta[italic(H2)]),
          expression(beta[italic(K1)]),
          expression(beta[italic(K2)]))
if(SAVEPLOT) pdf(paste0("./figures/betaPair.pdf"), height = 6, width = 6)
par(mar = rep(0,4))
pairs(allThetaNA, pch = 20, cex = 0.4,
      labels = lbAll,
      gap = 0.2,ylim=c(-1.25,1.25),xlim=c(-1.25,1.25),
      upper.panel=NULL,  las=2,
      col = rgb(0, 0, 255, max = 255, alpha = 80, names = "blue50"))
if(SAVEPLOT) dev.off()

# library(GGally)
# ggpairs(allThetaNA, upper = NULL,
#         diag = NULL,
#         colour = rgb(0, 0, 255, max = 255, alpha = 80, names = "blue50"))
# OmegaMean = (2 * eta_Omega2 + 8 + 1) * solve(eta_Omega1) * (-0.5)
# SigmaMean = solve(OmegaMean)
# library(ggcorrplot)
# ggcorrplot(SigmaMean)
# plot(eigen(SigmaMean)$values )
# plot(eigen(SigmaMean)$vectors[,1])
# plot(eigen(SigmaMean)$vectors[,2])



# ------- compareMagnitude -------
## Compare the difference of the flux corrected average magnitude
## Compare the J-band result.
if(SAVEPLOT) pdf(paste0("./figures/realDeltaJ.pdf"), height = 5, width = 5)
par(mar  = c(4,4,1,1))
plotMagDiff(fluxCorrectedPL$logP,  allPL$J, msp_Jmean,
            ylab = expression(paste(Delta, italic(J)," [mag]")))
if(SAVEPLOT) dev.off()

## Compare the H-band result.
if(SAVEPLOT) pdf(paste0("./figures/realDeltaH.pdf"), height = 5, width = 5)
par(mar  = c(4,4,1,1))
plotMagDiff(fluxCorrectedPL$logP, allPL$H, msp_Hmean, 
            ylab = expression(paste(Delta,  italic(H), " [mag]")))
if(SAVEPLOT) dev.off()

## Compare the K-band result.
if(SAVEPLOT) pdf(paste0("./figures/realDeltaK.pdf"), height = 5, width = 5)
par(mar  = c(4,4,1,1))
plotMagDiff(fluxCorrectedPL$logP, allPL$K, msp_Kmean, 
            ylab = expression(paste(Delta, italic(K)[s], " [mag]")))
if(SAVEPLOT) dev.off()

# ------- comparePeriod -------
## Compare the SVI period and MSP period
if(SAVEPLOT) pdf(paste0("./figures/realPeriod.pdf"), height = 5, width = 5)
par(mar  = c(4,4,1,1))
plot(svi_period, msp_period, 
     pch = 20, col = "black", cex = 0.5,
     ylim = c(100,1000), xlim = c(100, 1000), 
     type = "n", xlab = "SVI Period", ylab = "MSP Period")
abline(a = 0, b = 1, col = "black", lty = 2)
arrows(svi_period - sqrt(allPL[,2]), msp_period,
       svi_period + sqrt(allPL[,2]), msp_period,
       code = 3, length = 0.02, angle = 90, col = "grey")
arrows(svi_period, msp_period + yuan2018Table$eP[aPos],
       svi_period, msp_period - yuan2018Table$eP[aPos],
       code = 3, length = 0.02, angle = 90, col = "grey")
points(svi_period, msp_period, cex = 0.5, pch = 20,
       col = rgb(0, 0, 255, max = 255, alpha = 125, names = "blue50"))
points(svi_period[IncP], msp_period[IncP], cex = 0.6, pch = 20, 
       col = "red")
if(SAVEPLOT) dev.off()


# ------- plotPLR -------
## Fit PLR for each band with the uncorrected magnitude.
## Get intercept and residual sigma.
qPLRmatrix = matrix(0, 3, 7)
colnames(qPLRmatrix) = c("a0", "a0 sigma", 
                        "a1","a1 sigma", 
                        "a2","a2 sigma", "residual sigma")
lPLRmatrix = matrix(0, 3, 5)
colnames(lPLRmatrix) = c("a0", "a0 sigma", 
                         "a1","a1 sigma", "residual sigma")

rownames(qPLRmatrix) = c("J band", "H band", "K band")
ConvertM = matrix(0, 3,3)
ConvertM[1,1] = 1
ConvertM[1,2] = -2.3
ConvertM[1,3] = -2.3^2 + 4.6*2.3
ConvertM[2,2] = -1
ConvertM[2,3] = 4.6
ConvertM[3,3] = 1

for(bI in 2:4){
  alphaMean = -0.5*solve(eta_alpha1[[bI]], eta_alpha2[[bI]])
  SigmaMat = -0.5 * solve(eta_alpha1[[bI]])
  qPLRmatrix[bI-1, c(1,3,5)] = ConvertM %*% alphaMean
  qPLRmatrix[bI-1, c(2,4,6)] = sqrt(diag(ConvertM %*% SigmaMat %*% t(ConvertM)))
  
}

# $19.01\pm 0.01$ & $-3.36\pm0.05$ &  $-1.52\pm0.02$
# $18.27\pm0.01$ & $-3.11\pm0.04$ & $-3.24\pm$0.02$
# $17.87\pm0.01$ & $-3.68\pm0.04$ & $-2.20\pm0.02$




for(bI in 2:4){
    ## remove sample without enough band points
    logP = fluxCorrectedPL$logP
    logPSigma = fluxCorrectedPL$logPSigma
    magSel = allPL[, bandName[bI]]
    magSigma = allPL[, paste0(bandName[bI], "Sigma") ]

    if(SAVEPLOT) pdf(paste0("./figures/realPLR-",bandName[bI], ".pdf"), height = 4.5, width = 4)
    par(mar = c(3,2,0.1,0.1))
    keppAll = rep(TRUE, length(logP))
    plotPLR(logP, logPSigma, magSel, magSigma,ylabVec[bI], keppAll)
    alphaMean = -0.5*solve(eta_alpha1[[bI]], eta_alpha2[[bI]])
    log10PSeq = log10(1/seq(100, 1000, length.out = 100))
    plrSeq = alphaMean[1] + alphaMean[2] * log10PSeq  + alphaMean[3] * (log10PSeq )^2
    lines((-log10PSeq), plrSeq, lwd = 2, col = "red")

    magHat = alphaMean[1] - alphaMean[2] * logP  + alphaMean[3] * (logP)^2
    qPLRmatrix[bI-1, 7] = sd(magSel - magHat)
    # fit linear PLR
    sel = logP < log10(400)
    y = magSel[sel]
    x1 = logP[sel] - 2.3
    x1Sigma = logPSigma[sel]
    lfit = summary(lm(y~x1))
    
    #lPLRmatrix[bI-1, c(1,3)] = lfit$coefficients[,1]
    lPLRmatrix[bI-1, c(2,4)] = lfit$coefficients[,2]
    lPLRmatrix[bI-1, 5] = sd(lfit$residuals)
    lPLRmatrix[bI-1,3] =     cov(x1,y) / (var(x1) - mean(x1Sigma^2))
    lPLRmatrix[bI-1, 1] = mean(y) - mean(x1) * lPLRmatrix[bI-1,3]
    
    log10PSeq = log10(seq(100,400,length.out = 100))
    mSeq = lPLRmatrix[bI-1, 1] + lPLRmatrix[bI-1,3] * (log10PSeq - 2.3)
    lines(log10PSeq, mSeq, lwd = 2, col = "black")
    
    if(SAVEPLOT) dev.off()
}
#PLRmatrix = data.frame(PLRmatrix)

## Compute again with the corrected magnitude
PLRmatrixCorrected = matrix(0, 3, 6)
colnames(PLRmatrixCorrected) = c("Linear a0", "Linear a0 sigma", 
                        "Linear residual sigma",
                        "Quadratic a0", "Quadratic a0 sigma", 
                        "Quadratic residual sigma")
rownames(PLRmatrixCorrected) = c("J band", "H band", "K band")
for(bI in 2:4){
    ## remove sample without enough band points
    sel =  !is.na(fluxCorrectedPL[, bandName[bI]])
    logP = fluxCorrectedPL$logP[sel]
    logPSigma = fluxCorrectedPL$logPSigma[sel]
    
    ## use the corrected mag
    magSel = fluxCorrectedPL[sel, bandName[bI]]
    magSigma = fluxCorrectedPL[sel, paste0(bandName[bI], "Sigma") ]
    res2 = fitQuadPLR(quadraticCoefLMC[[bI]], logP, logPSigma, magSel, magSigma)
    sel = logP < log10(400)
    res1 = fitLinearPLR(linearCoefLMC[[bI]], logP[sel], logPSigma[sel], magSel[sel], magSigma[sel])
    PLRmatrixCorrected[bI - 1, 1:2] = res1$interc
    PLRmatrixCorrected[bI - 1, 3] = res1$resSigma
    PLRmatrixCorrected[bI - 1, 4:5] = res2$interc
    PLRmatrixCorrected[bI - 1, 6] = res2$resSigma
}

# ------- distanceMod -------
## Compute the table for distance modulus
distModu = matrix(0,3,7)
distModuSigma = matrix(0,3,7)
bNames = c("I","J","H","K")
for(bI in 2:4){
    # Delta intercpet
    distModu[bI-1,1] = qPLRmatrix[bI-1,1] - quadraticCoefLMC[[bI]][1]
    # mean magnitude correction
    diffM = -allPL[,bNames[bI]] + fluxCorrectedPL[,bNames[bI]]
    distModu[bI-1,2] = mean(diffM, na.rm = T)
    distModuSigma[bI-1,2]  = sd(diffM, na.rm = T) / sqrt(sum(!is.na(diffM)))
}
# Delta A_lambda for extinction
distModu[,3] = c(0.029, 0.018, 0.012) 
# Delta ct for color term correction
distModu[,4] = c(0.016, 0.010, -0.007)
# Distance modulus difference
distModu[,5] = rowSums(distModu[,1:4])
# Distance to LMC
distModu[,6] = rep(18.493, 3)
# Distance to M33
distModu[,7] = distModu[,5] + distModu[,6]
distModu = round(distModu,3)

# Uncertainty propagation
distModuSigma[,1] =  2*0.01^2
distModuSigma[,3] = c(0.008, 0.005, 0.003)^2
distModuSigma[,4] = c(0.036, 0.040, 0.032)^2
distModuSigma[,5] = rowSums(distModuSigma[,1:4])
distModuSigma[,6] = rep(0.048^2, 3)
distModuSigma[,7] = distModuSigma[,5] + distModuSigma[,6]
distModuSigma = sqrt(distModuSigma)
distModuSigma = round(distModuSigma, 3)

distModu = data.frame(distModu)
distModuSigma = data.frame(distModuSigma)
colnames(distModu) = c("Da0", "Dmbar", "DAlambda",
                       "Dct", "Dmu", "muLMC","muM33")
colnames(distModuSigma) = c("Da0", "Dmbar", "DAlambda",
                       "Dct", "Dmu", "muLMC","muM33")
rownames(distModu) = c("J", "H", "Ks")
rownames(distModuSigma) = c("J", "H", "Ks")

