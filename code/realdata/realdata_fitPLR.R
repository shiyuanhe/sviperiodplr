## fit the PLR to the real data
dALambda = c(0, 0.029, 0.018, 0.012)
dALambdaSigma = c(0, 8,5,3)/1000
dct = c(0, 0.016, 0.010, -0.007)
dctSigma = c(0,36,40,32)/1000

## This is the quadratic PLR coef from LMC
quadraticCoefLMC = list(c(0,0,0),
                     c(12.7, -3.49, -1.54),
                     c(11.96, -3.59, -3.40),
                     c(11.59, -3.77, -2.23))

## Fit the quadratic PLR to the dataset. 
## Only the intercept term is determined from the dataset. 
fitQuadPLR = function(initCoef, 
                      log10P, log10PSigma,
                      mag, magSigma){
    
    log10PC = log10P - 2.3
    plresidual = mag - initCoef[2] * log10PC  - initCoef[3] * log10PC^2
    plrSigma = sqrt(magSigma^2 + initCoef[2]^2 * log10PSigma^2 +  initCoef[3]^2 * log10PSigma^2 * 4)
    
    # weighted estimation of the mean
    tmp = iterativeClipping(plresidual, plrSigma)
    interc = tmp$interc
    intercSigma = tmp$intercSigma
    resSigma = tmp$resSigma
    keepI = tmp$keep
    
    # a dense curve of fitted PLR for drawing
    log10PSeq = log10(seq(100, 1000, length.out = 100))
    plrSeq = interc + initCoef[2] * (log10PSeq - 2.3) + initCoef[3] * (log10PSeq - 2.3)^2
    fittedCurve = cbind(log10PSeq, plrSeq)
    
    resultList = list(
        interc = c(interc, intercSigma),
        fittedCurve = fittedCurve,
        resSigma = resSigma,
        keepI = keepI,
        plresidual = plresidual)
    return( resultList )
}


## This is the linear PLR coef from LMC
linearCoefLMC = list(c(0,0,0),
                  c(12.67, -3.48),
                  c(11.91, -3.64),
                  c(11.56, -3.77))
# alphaTmp = linearCoef[[bI]]
# period = finalPLValues[,1]
# mag = finalPLValues[,bI+1]
# sigmaM = allPL[2*bI+2,]
# sel =  mag > 0 & period < 400
# period = period[sel]

fitLinearPLR = function(initCoef, 
                        log10P, log10PSigma,
                        mag, magSigma){
    
    log10PC = log10P - 2.3
    plresidual = mag - initCoef[2] * log10PC
    plrSigma = sqrt(magSigma^2 + initCoef[2]^2 * log10PSigma^2)

    # weighted estimation of the mean
    # interc = sum(plresidual / plrSigma^2) / sum(1/plrSigma^2)
    # # estimation uncertainty
    # intercSigma = sqrt(1/sum(1/plrSigma^2))
    # resSigma = sd(plresidual - interc)
    tmp = iterativeClipping(plresidual, plrSigma)
    interc = tmp$interc
    intercSigma = tmp$intercSigma
    resSigma = tmp$resSigma
    keepI = tmp$keep
    # a dense curve of fitted PLR for drawing
    log10PSeq = log10(seq(100, 400, length.out = 100))
    plrSeq = interc + initCoef[2] * (log10PSeq - 2.3)
    fittedCurve = cbind(log10PSeq, plrSeq)
    
    resultList = list(
        interc = c(interc, intercSigma),
        resSigma = resSigma,
        keepI = keepI,
        fittedCurve = fittedCurve)
    return( resultList )
}

iterativeClipping = function(plresidual, plrSigma){
    # weighted estimation of the mean
    keep = rep(TRUE, length(plrSigma))
    CTU = TRUE
    k = 0
     while(CTU){
        CTU = FALSE
        interc = sum(plresidual[keep] / plrSigma[keep]^2) / 
            sum(1/plrSigma[keep]^2)
        resSigma = median(abs(plresidual[keep] - interc))/  0.7978846
        sel = abs(plresidual - interc) / resSigma > 3
        selNew =  abs(plresidual[keep] - interc) / resSigma > 3
        keep[sel] = FALSE
        if(any(selNew)) CTU = TRUE
     }
    # estimation uncertainty
    interc = sum(plresidual[keep] / plrSigma[keep]^2) / sum(1/plrSigma[keep]^2)
    resSigma = sd(plresidual[keep] - interc)
    intercSigma = sqrt(1/sum(1/plrSigma[keep]^2))
    return(list(interc = interc,
                resSigma = resSigma,
                intercSigma = intercSigma,
                keep = keep))
}

# getLinearIntercept(2) - 12.67 + 0.029 + 0.016
# getLinearIntercept(3) - 11.91 + 0.018 + 0.01
# getLinearIntercept(4) - 11.56 + 0.012 - 0.007
plotPLR = function(logP, logPSigma, magSel, magSigma,ylab,keepI){
    plot(logP, magSel,
         ylim = rev(range(magSel)),
         xlab = "", ylab = "",
         pch = 20, type = "n")
    xlab = expression(paste(log[10], "(Period)"))
    mtext(xlab , side= 1, line = 2,cex = 0.9)
    mtext(ylab, side= 2, line = -1.2,cex = 1)
    
    arrows(logP - logPSigma, magSel, 
           logP + logPSigma, magSel, 
           length = 0.02, angle = 90, code = 3,
           col = "grey")
    arrows(logP, magSel - magSigma, 
           logP, magSel + magSigma, 
           length = 0.02, angle = 90, code = 3,
           col = "grey")
    points(logP[!keepI], magSel[!keepI], pch = 20, cex = 0.6, 
           col = rgb(255, 0, 0, max = 255, alpha = 125, names = "red50"))
    points(logP[keepI], magSel[keepI], pch = 20, cex = 0.6, 
           col = rgb(0, 0, 255, max = 255, alpha = 90, names = "blue50"))
    
}

plotMagDiff = function(logP, mag1, mag2, ylab){
    magDiff = mag1 - mag2
    plot(logP,  magDiff, pch = 19, cex = 0.5, 
         col = rgb(0, 0, 255, max = 255, alpha = 125, names = "blue50"),
         xlab = expression(paste(log[10], "(Period)") ),
         ylab = ylab)
    mDiff = median( magDiff, na.rm = TRUE)
    
    abline(h = mDiff, lty = 2)

}

