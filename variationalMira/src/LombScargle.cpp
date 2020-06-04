#include "LombScargle.hpp"

lombScargle::lombScargle(vector<vec> tt_, vector<vec> yy_,
                         vector<vec> ss_, vec nSample_):
    coreCov(zeros<vec>(3))
{
    tt = tt_;
    yy = yy_;
    ss = ss_;
    phase = tt;
    residual = yy;
    nSample = nSample_;
    nBand = nSample.size();
    aveMag = ones<vec>(nBand) * (-100000);
}



// find the optimal frequency by multi-band lomb-scargle
double lombScargle::compute()
{
    int bI, fI;
    double resid = 0.0;
    periodogram = zeros<vec>(nSeq);
    double cf = fmin;
    vec beta;
    for(fI = 0; fI < nSeq; fI++){
        resid = 0.0;
        for(bI = 0; bI < nBand; bI++){
            if (notEnoughBandSample(nSample.at(bI) - 2) ) continue;
            resid += fitBandI(cf, bI, beta);
        }
        cf += fdelta;
        periodogram(fI) = resid;
    }


    int optPos = periodogram.index_min();
    fhat = fmin + fdelta * optPos;
    set_freq(fhat);

    return fhat;
}

void lombScargle::set_freq(double freq_){
    int bI;
    vec beta;
    fhat = freq_;
    for(bI = 0; bI < nBand; bI++){
        if (notEnoughBandSample(nSample.at(bI) - 2) ) continue;
        phase.at(bI) = tt.at(bI) * fhat;
        fitBandI(fhat, bI, beta);
    }
}

double lombScargle::fitBandI(double freq, size_t bI, vec & beta){
    vec tmpy;
    mat basisMat = ones<mat>(nSample(bI), 3), tmpX;
    basisMat.col(1) = cos(tt.at(bI) * freq * 2 * 3.1415926);
    basisMat.col(2) = sin(tt.at(bI) * freq * 2 * 3.1415926);


    vec weights = sqrt(ss.at(bI)) + 1e-8;
    weights = 1 / weights;
    tmpy = yy.at(bI);// % weights;
    tmpX = basisMat;//.each_col() % weights;

    tmpy = tmpX.t() * tmpy;
    tmpX = tmpX.t() * tmpX;
    // tmpX(1,1) += 1e-3;
    // tmpX(2,2) += 1e-3;

    beta = solve(tmpX, tmpy);
    aveMag(bI) = beta(0);
    residual.at(bI) = yy.at(bI) - basisMat * beta;

    return sum(residual.at(bI) % residual.at(bI));
}


// Get covariance derivative for the negative log-likelihood
// Only to squared exponetial class
vec lombScargle::get_covDeriv(vec betaNew){


    mat covCoreMat;
    vec tmpY;
    vec result = zeros<vec>(3);

    if (notEnoughBandSample(nSample.at(currentB) - 2) ) return result;

    coreCov.set_beta(betaNew);
    coreCov.setupData(phase.at(currentB), ss.at(currentB));

    tmpY = coreCov.inverseOnVector(residual.at(currentB));

    covCoreMat = tmpY * tmpY.t();
    covCoreMat -=  coreCov.getInverse();

    result(0) = dot(covCoreMat, coreCov.covPartial_j(0));
    result(1) = dot(covCoreMat, coreCov.covPartial_j(1));
    result(2) = dot(covCoreMat, coreCov.covPartial_j(2));

    result *= -0.5;

    return result;
}

// This function serves to provide initial theta estimation.
// The expected negative loglikelihood under variational parameter
double lombScargle::gpNegativeLogLik(vec betaNew){
    double result = 0.0;

    if (notEnoughBandSample(nSample.at(currentB) - 2) ) return result;
    vec tmpy;
    coreCov.set_beta(betaNew);
    coreCov.setupData(phase.at(currentB), ss.at(currentB));

    tmpy = residual.at(currentB);
    tmpy = coreCov.halfActOnVector(tmpy);
    result = 0.5 * dot(tmpy, tmpy);
    result += 0.5 * coreCov.logdet();

    return result;
}




vec lombScargle::getBandIFreqMag(){
    vec result(2);
    result(0) = fhat;
    result(1) = aveMag(currentB);// mean(yy.at(currentB));
    if (notEnoughBandSample(nSample.at(currentB) - 2) ) result(1) = -10000;
    return result;
}
