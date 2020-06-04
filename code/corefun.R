
convertList = function(mjd, mag, sigma, band, lev){
    l1 = l2 = l3 = list()
    l4 = rep(0, length(lev))
    k = 0
    sigma = sigma^2
    for(ll in lev){
        sel = (band == ll)
        k = k + 1
        l4[k] = sum(sel)
        if(l4[k] ==0){
            l1 = c(l1, list(-1e5))
            l2 = c(l2, list(-1e5))
            l3 = c(l3, list(-1e5))
        }else{
            l1 = c(l1, list(mjd[sel]))
            l2 = c(l2, list(mag[sel]))
            l3 = c(l3, list(sigma[sel]^2))
        }
    }
    return(list(l1, l2, l3, l4))
}


# readInCompute = function(filePath){
#     lcData = read.table(filePath)
#     inputData = with(lcData, convertList(V1, V3, V4, V2, lev))
#     miraObj = new(miraClass)
#     miraObj$set_data(inputData[[1]], inputData[[2]], inputData[[3]], inputData[[4]])
#     miraObj$set_globalFreq(fmin, fmax, nseq)
#     Omega = diag(1, nBand *3) * 1e-3
#     miraObj$set_globalBeta(Omega)
#     
#     tmp = rep(0,3)
#     viU_alpha_ = list(tmp, tmp, tmp, tmp)
#     tmp = diag(1,3)
#     viS_alpha_ = list(tmp, tmp, tmp, tmp)
#     miraObj$set_globalPL(viU_alpha_, viS_alpha_)
#     miraObj$compute()
#     return(miraObj)
# }




readInCompute = function(filePath){
    lcData = read.table(filePath)
    inputData = with(lcData, convertList(V1, V3, V4, V2, lev))
    miraObj = new(miraClass)
    miraObj$set_data(inputData[[1]], inputData[[2]], inputData[[3]], inputData[[4]])
    miraObj$set_globalFreq(fmin, fmax, nseq)
    
    Iset = c(1,4, 7, 10)
    Jset = c(2,3,5,6,8,9,11,12)
    Omega = matrix(0, 12, 12)
    Omega[Iset, Iset] = diag(gammaExpt[,1])
    Omega[Jset, Jset] = OmegaExpt
    miraObj$set_globalBeta(Omega)
    
    miraObj$set_globalPL(viU_alpha, viS_alpha)
    miraObj$compute()
    return(miraObj)
}
