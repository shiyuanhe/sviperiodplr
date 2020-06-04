## Plot the final figures and Tables for the real data analysis
## This file works with 3sigma-clipping. The a_1 and a_2 are fixed.
## Only the PLR intercept are fitted for distance estimation.

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
sel = yuan2018Table$cl == "O"
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
      gap = 0.2,
      upper.panel=NULL,
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
PLRmatrix = matrix(0, 3, 6)
colnames(PLRmatrix) = c("Linear a0", "Linear a0 sigma", 
                        "Linear residual sigma",
                        "Quadratic a0", "Quadratic a0 sigma", 
                        "Quadratic residual sigma")
rownames(PLRmatrix) = c("J band", "H band", "K band")

keepMat = matrix(FALSE, nrow(fluxCorrectedPL),3)
for(bI in 2:4){
  ## remove sample without enough band points
  sel =  !is.na(fluxCorrectedPL[, bandName[bI]])
  keepIndex1 = which(sel)
  logP = fluxCorrectedPL$logP[sel]
  logPSigma = fluxCorrectedPL$logPSigma[sel]
  magSel = allPL[sel, bandName[bI]]
  magSigma = allPL[sel, paste0(bandName[bI], "Sigma") ]
  res2 = fitQuadPLR(quadraticCoefLMC[[bI]], logP, logPSigma, magSel, magSigma)
  keepMat[keepIndex1[res2$keepI], bI-1] = TRUE
}
orichIndex = rowSums(keepMat)==3




for(bI in 2:4){
    ## remove sample without enough band points
    logP = fluxCorrectedPL$logP
    logPSigma = fluxCorrectedPL$logPSigma
    magSel = allPL[, bandName[bI]]
    magSigma = allPL[, paste0(bandName[bI], "Sigma") ]

    if(SAVEPLOT) pdf(paste0("./figures/realPLR-",bandName[bI], ".pdf"), height = 4.5, width = 4)
    par(mar = c(3,2,0.1,0.1))

    plotPLR(logP, logPSigma, magSel, magSigma,ylabVec[bI], orichIndex)
    alphaMean = -0.5*solve(eta_alpha1[[bI]], eta_alpha2[[bI]])
    log10PSeq = log10(1/seq(100, 1000, length.out = 100))
    plrSeq = alphaMean[1] + alphaMean[2] * log10PSeq  + alphaMean[3] * (log10PSeq )^2
    lines((-log10PSeq), plrSeq, lwd = 2, col = "black")

        
    res2 = fitQuadPLR(quadraticCoefLMC[[bI]], 
                      logP[orichIndex], logPSigma[orichIndex],
                      magSel[orichIndex], magSigma[orichIndex])
    
    lines(res2$fittedCurve, lwd = 3, col = "red", lty = 2)
    sel = (logP < log10(400)) & orichIndex
    res1 = fitLinearPLR(linearCoefLMC[[bI]], 
                        logP[sel], logPSigma[sel], 
                        magSel[sel], magSigma[sel])
    lines(res1$fittedCurve, lwd = 3, col = "red")

    PLRmatrix[bI - 1, 1:2] = res1$interc
    PLRmatrix[bI - 1, 3] = res1$resSigma
    PLRmatrix[bI - 1, 4:5] = res2$interc
    PLRmatrix[bI - 1, 6] = res2$resSigma
    
    if(SAVEPLOT) dev.off()
}
PLRmatrix = data.frame(PLRmatrix)

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
for(bI in 2:4){
    # Delta intercpet
    distModu[bI-1,1] = PLRmatrix[bI-1,1] - linearCoefLMC[[bI]][1]
    # mean magnitude correction
    distModu[bI-1,2] = PLRmatrixCorrected[bI-1,1] - PLRmatrix[bI-1,1]
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
distModuSigma = matrix(0,3,7)
distModuSigma[,1] = PLRmatrix[,2]^2 + 0.01^2
distModuSigma[,2] = PLRmatrix[,2]^2 + PLRmatrixCorrected[,2]^2
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

# 
# cat("\n\n")
# 
# dMu = res1$interc[1] - linearCoef[[bI]][1] + dALambda[bI] + dct[bI]
# dMuSigma = sqrt(res1$interc[2]^2 + dALambdaSigma[bI]^2 + dctSigma[bI]^2)
# cat("Delta Mu and sigma", dMu, dMuSigma, "\n")
# 
# muM33 = dMu + 18.493
# muM33Sigma = sqrt(dMuSigma^2 + 0.048^2)
# cat("MuM33 and Sigma", muM33, muM33Sigma, "\n")
# 
# cat(bandName[bI], " band ", "\n")
# 
# 
# 
# ## Fit PLR for each band with the uncorrected magnitude.
# ## Get intercept and residual sigma.
# for(bI in 2:4){
#     logPAll = fluxCorrectedPL$logP
#     magAll = fluxCorrectedPL[, bandName[bI]]
#     sel =  !is.na(magAll)
#     logP = logPAll[sel]
#     logPSigma = fluxCorrectedPL$logPSigma[sel]
#     magSel = magAll[sel]
#     magSigma = fluxCorrectedPL[sel, paste0(bandName[bI], "Sigma") ]
#     IncP_tmp = IncP[sel]
#     mspSub = log10(msp_period[sel])
#     
#     pdf(paste0("./figures/realPLR-",bandName[bI], ".pdf"), height = 4, width = 7)
#     par(mar = c(4,4,0.5,0.5))
#     plotPLR(logP, magSel, ylabVec[bI])
#     
#     cat(bandName[bI], " band ", "\n")
#     res2 = fitQuadPLR(quadraticCoef[[bI]], logP, logPSigma, magSel, magSigma)
#     lines(res2$fittedCurve, lwd = 3, col = "blue", lty = 2)
#     cat("quadratic intercept and sigma: ", res2$interc, "\n")
#     
#     
#     sel = logP < log10(400)
#     res1 = fitLinearPLR(linearCoef[[bI]], logP[sel], logPSigma[sel], magSel[sel], magSigma[sel])
#     lines(res1$fittedCurve, lwd = 3, col = "red")
#     cat("linear intercept and sigma: ", res1$interc, "\n")
#     
#     cat("Delta a0 ", res1$interc[1] - linearCoef[[bI]][1], "\n")
#     
#     dMu = res1$interc[1] - linearCoef[[bI]][1] + dALambda[bI] + dct[bI]
#     dMuSigma = sqrt(res1$interc[2]^2 + dALambdaSigma[bI]^2 + dctSigma[bI]^2)
#     cat("Delta Mu and sigma", dMu, dMuSigma, "\n")
#     
#     muM33 = dMu + 18.493
#     muM33Sigma = sqrt(dMuSigma^2 + 0.048^2)
#     cat("MuM33 and Sigma", muM33, muM33Sigma, "\n")
#     
#     cat("\n\n")
#     dev.off()
#     
# }
# 
# 
# ## Fit PLR for each band
# for(bI in 2:4){
#     cat(bandName[bI], "\n")
#     logPAll = fluxCorrectedPL$logP
#     magAll = fluxCorrectedPL[, bandName[bI]]
#     sel =  !is.na(magAll)
#     logP = logPAll[sel]
#     logPSigma = fluxCorrectedPL$logPSigma[sel]
#     magSel = magAll[sel]
#     magSigma = fluxCorrectedPL[sel, paste0(bandName[bI], "Sigma") ]
#     IncP_tmp = IncP[sel]
#     
#     mspP = log10(yuan2018Table$P[sel])
#     mspPSigma = yuan2018Table$eP[sel] / mspP * log(10)
#     mspMag = yuan2018Table[sel, bandName[bI]]
#     mspMagSigma = yuan2018Table[sel, paste0("e", bandName[bI] ) ]
#     print(sum(sel))
#     res1 = fitQuadPLR(quadraticCoef[[bI]], logP, logPSigma, magSel, magSigma)
#     res2 = fitQuadPLR(quadraticCoef[[bI]], mspP, mspPSigma, mspMag, mspMagSigma)
# 
#     cat("quadratic sigma SVI ", res1$resSigma,
#         "MSP", res2$resSigma, 
#         "Reduction", 1 - res1$resSigma/ res2$resSigma,
#         "\n")
#     
#     sel = logP < log10(400)
#     print(sum(sel))
#     res3 = fitLinearPLR(linearCoef[[bI]], logP[sel], logPSigma[sel], magSel[sel], magSigma[sel])
#     sel = mspP < log10(400)
#     res4 = fitLinearPLR(linearCoef[[bI]], mspP[sel], mspPSigma[sel], mspMag[sel], mspMagSigma[sel])
# 
#     cat("linear sigma SVI ", res3$resSigma,
#         "MSP", res4$resSigma, 
#         "Reduction", 1 - res3$resSigma/ res4$resSigma,
#         "\n")
#     
# 
# }


#[IncP_tmp]
# res2 = fitQuadPLR(quadraticCoef[[bI]], logP, logPSigma, magSel, magSigma)
# res3 = fitQuadPLR(quadraticCoef[[bI]], mspSub, logPSigma, magSel, magSigma)
# plot(abs(res2$plresidual[IncP_tmp] - res2$interc[1])-
#      abs(res3$plresidual[IncP_tmp]- res3$interc[1]))
# PLR1 = data.frame(residual = res2$plresidual[IncP_tmp] - res2$interc[1], method = "SVI", stringsAsFactors = F)
# PLR2 = data.frame(residual = res3$plresidual[IncP_tmp]- res3$interc[1], method = "MSP", stringsAsFactors = F)
# PLR = rbind(PLR1, PLR2)
# colnames(PLR) = c("residual", "method")
# library(ggplot2)
# p = ggplot(PLR, aes(residual, ..count.., fill = method)) + geom_histogram()
# print(p)





# starI = 783
# tmp = modelObj$get_qf(starI -1)
# plot(tmp[,1], sqrt(tmp[,2]), type = "l")
# abline(v = 1 / allPL[1,starI])
# abline(v = 1 /yuan2018Table$P[aPos][starI])

# 
# 
# allPL[,1]
# inputData = with(lcData, convertList(V1, V3, V4, V2, lev))
# for(b in 1:nBand){
#     timeB = inputData[[1]][[b]]
#     magB = inputData[[2]][[b]]
#     tseq = as.matrix(seq(tmin, tmax, by = 1), ncol = 1)
#     plot(tseq, yseq[,b], type = "l")
#     points(timeB, magB, pch = 20)
#     allPL[,1]
#     plot(qf, type = "l")
#     1/qf[which.max(qf[,2]),]
# }
## Check convergance for gamma
## Each column is the parameter for each round of iteration. 
# x = modelObj$get_gammaTrace()
# plot(x[4,], pch = 20, type = "l")
# plot(x[3,], pch = 20, type = "l")
# plot(x[2,], pch = 20, type = "l")
# plot(x[1,], pch = 20, type = "l")
# # 
# # # Check convergence for alpha
# y = modelObj$get_alphaTrace()
# plot(y[4,], type = "l")
# plot(y[5,], type = "l")
# plot(y[6,], type = "l")
# 
# # Check convergence for kernel parameter
# z = modelObj$get_kernelTrace()
# plot(z[10,], type = "l")
# plot(z[11,], type = "l")
# plot(z[12,], type = "l")

## 
## The first row: estimated period
## The Second row: sigma uncertainty of period
## The 3rd, 5th, 7th, 9th row: the mean magnitude for the I J H K bands respectively.
## The 4th, 6th, 8th, 10th row: the uncertainty for the mean magnitude for each band. 

# 
# # H
# a0 = 18.27
# a1 = -3.59
# a2 = -3.4
# #  8.541 8.594147
# 
# # K Data
# a0 = 17.83
# a1 = -3.77
# a2 = -2.23
# # 14.7043 14.75115
# a0 - a1 * 2.3 + a2 * 2.3^2