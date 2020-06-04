# This file contains the final proposed strategy for parameter initialization.


library(readr)
library(multiband)


## prepare data for the required format of the multiband package
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

##' @param allFiles: Full file path list
##' @param nSample: the number of samples to be used for initilization
##' @param fmin, fmax: the range of frequency search
##' @param nseq: number of trial frequency inside [fmin, fmax]
##' @param lev: vector of band names, e.g. c("I","K", "J", "H")
##' @param alphaBarList: list of PLR coefficients, each item in the list 
##'             is the coefficien for one band. The itercept will be updated. 
##' @return 
##'     gammaBar: the precision (1/sigma^2) for the PLR of each band. 
##'     alphaBarList: the updated list of PLR coefficients.
##'     keepS: a logical vector. keepS[i] = FALSE if the i-th sample is outlier.
##'     freqAll: estimated frequency for each observation.
initPLR = function(allFiles, nSample,
                   fmin, fmax, nseq,
                   lev,alphaBarList){
    nBands = length(lev)
    fSeq = seq(fmin, fmax, length.out = nseq)
    periodGrid = 1 / fSeq
    magfreq = matrix(NA, nSample, nBands + 1)
    magfreq = data.frame(magfreq)

    ## Load all the samples
    for(ii in 1:nSample){
        filePath = allFiles[ii]
        lcData = read.table(filePath)
        inputData = multibandData(lcData)
        result = pgls(inputData, periods = periodGrid, gamma1 = 0, gamma2 = 0)
        #max_iter = 1000, fast_BCD_flag = FALSE)
        #plot(1/result$period_seq_all, result$rss_ls, type ="l")
        optIndex = which.min(result$rss_ls)
        fHat = 1/result$period_seq_all[optIndex]
        mHat = result$best_fitLS[,1]
        
        mHatLong = rep(NA, nBands)
        bands = unique(lcData$V2)
        relPos = match(bands, lev)
        mHatLong[relPos] = mHat
        res = c(fHat, mHatLong)
        magfreq[ii, ] = res
    }
    

    gammaBar = rep(0, nBands)
    for(bI in 1:nBands){
        tmpData = magfreq[, c(1, bI+1)]
        colnames(tmpData) = c("freq", "mag")
        
        # where to remove missing values and outliers
        keepS = !is.na(tmpData$mag)
        errVec = abs( tmpData$mag - median(tmpData$mag, na.rm = T) )
        sdErr = median(errVec, na.rm = TRUE) / 0.8 * 10 
        sel2 = which(errVec > sdErr) # samples outside of 10*sigma
        keepS[sel2] = FALSE # marked as FALSE if the sample is outlier

        ## Robust fitting 
        alphaBar = alphaBarList[[bI]]
        tmp = initPLRPara(tmpData, alphaBar[2], alphaBar[3], keepS, FALSE)
        gammaBar[bI] = tmp$gammab
        alphaBar[1] = tmp$alpha1
        keepS = tmp$keepS
        alphaBarList[[bI]] = alphaBar
    }
    freqAll = magfreq[,1]
    
    return(list(gammaBar = gammaBar, 
                alphaBarList = alphaBarList, 
                keepS = keepS, # keepS[i] = FLASE if the i-th sample is outlier; otherwise,  keepS[i] = FLASE 
                freqAll = freqAll) # the estimated frequency for each star
           )
    
}

##' Get estimation for the PL relation intercept
##' Input:
##' @param  magfreq: dataframe with two columns, magnitude and frequency. Each row is a sample
##' @param  alpha2: fixed linear term of PLR.
##' @param alpha3: fixed quadratic term of PLR.
##' @param keepS: When keepS(i)==1, the i-th sample will be used for
##'                compute the gaussian process parameter. This is used to remove outlier in computation. 
##' @return 
##'      alpha1: the intercep for the PLR
##'      gammab: the precision for the fited PLR, 1/\hat{sigm}a^2
initPLRPara = function(magfreq, alpha2, alpha3, keepS, showplot = FALSE){
    
    ## Robust estimation of intercept with median residual
    plresidual = magfreq$mag - log10(magfreq$freq)^2 * alpha3 - log10(magfreq$freq) * alpha2
    alpha1 = median(plresidual[keepS])
    #  inverse of median absolute deviation
    plresidual = plresidual - alpha1
    gammab = 0.79/median( abs(plresidual)[keepS] ) 
    
    
    rmIndex = which(abs(plresidual) > 2/gammab)
    keepS[rmIndex] = FALSE
    
    if(showplot){
        xx = seq(0.0005, 0.01, length.out = 100)
        yy =  log10(xx)^2 * alpha3 + log10(xx) * alpha2 + alpha1
        rrg = rev(range(yy))
        rrg[1] = rrg[1] + 2
        rrg[2] = rrg[2] - 2
        plot(xx, yy, type = "l", ylim = rrg)
        points(magfreq$freq, magfreq$mag, pch = 20, 
               cex = 0.5, col = c("grey","black")[keepS+1]
        )
    }
    
    return(list(alpha1 = alpha1, # estimated intercept for PLR
                gammab = gammab^2, # inverse of 
                keepS = keepS))
}



##' This function get initial kernel parameter 
##' @param allFiles: full file path list
##' @param  nSample: the number of samples to be used for initilization
##' @param freqAll: the estimated frequency of each sample
##' @param keepS: When keepS(i)==1, the i-th sample will be used for
##'                compute the gaussian process parameter. This is used to remove outlier in computation. 
##' @param gradScaling: 
getInitKernelPara = function(allFiles, nSample,
                             lev, gradScaling, 
                             freqAll, keepS){
    nBands = length(lev)
    lsCollection = new(hierarchicalLS)
    ## Load all the samples
    for(i in 1:nSample){
        filePath = allFiles[i]
        lcData = read.table(filePath)
        inputData = with(lcData, convertList(V1, V3, V4, V2, lev))
        lsCollection$push_data(inputData[[1]], inputData[[2]],
                               inputData[[3]], inputData[[4]])
    }
    
    lsCollection$set_individualFreq(freqAll)

    
    kernelBetaList = list()
    for(bI in 1:nBands){
        ## Start to work on the bI-th band
        lsCollection$setGPBand(bI - 1)
        lsCollection$set_gpUse(keepS)
        
        newBeta = c(-3, -3, -3) # initial value
        ## theta_j = exp(beta_j) for j=1,2,3
        ## k(x, x') = theta_2 *exp(-||x-x'||^2/ beta_3) + beta_1*I(x=x')
        
        lsCollection$set_scaling(gradScaling)
        ## get optimal kernel parameter
        res = optim(newBeta, lsCollection$gpNegativeLogLik,
                    lsCollection$get_covDeriv, 
                    method = "BFGS", control = list(maxit = 1000))
        kernelBetaList = c(kernelBetaList, list(res$par))
    }
    
    return( kernelBetaList )
    
}









