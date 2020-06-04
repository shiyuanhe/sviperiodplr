#include "hierarchicalLS.hpp"


hierarchicalLS::hierarchicalLS()
{
    nObs = 0;
    grad_scaling = 10000;
}


void hierarchicalLS::compute()
{
    int oI;
    for(oI = 0; oI < nObs; oI++){
        allObs.at(oI).set_freqParam(fMin, fMax, nfSeq);
        allObs.at(oI).compute();
    }
}


void hierarchicalLS::push_data(vector<vec> tt, vector<vec> yy,
                               vector<vec> ss, vec nSample)
{

    lombScargle tmp(tt, yy, ss, nSample);
    //tmp.set_freqParam(fMin, fMax, nfSeq);
    allObs.push_back(tmp);
    nObs++;
}

// Compute the negative log-likelihood for optimization.
// Only for the currentB-th band, use the "setGPBand()" first.
double hierarchicalLS::gpNegativeLogLik(vec newBeta){
    int oI;
    double result = 0;
    try{
        for(oI = 0; oI < nObs; oI++){
            if(gpUse(oI)==0) continue;
            result += allObs.at(oI).gpNegativeLogLik(newBeta);
        }
    }catch(std::runtime_error e){
        result = 1e50;
    }
    return result / grad_scaling;
}


// Compute the first order gradient for optimization.
// Only for the currentB-th band, use the "setGPBand()" first.
vec hierarchicalLS::get_covDeriv(vec betaNew){
    vec result(3, fill::zeros);
    int oI;
    for(oI = 0; oI < nObs; oI++){
        if(gpUse(oI)==0) continue;
        result += allObs.at(oI).get_covDeriv(betaNew);
    }
    result /= grad_scaling;
    return result;

}

// Get a matrix for two columns: freq, average mag.
// Only for the currentB-th band, use the "setGPBand()" first.
mat hierarchicalLS::getFreqMag()
{
    int oI, cI;
    mat result(nObs, 2);
    cI = 0;
    vec tmp;
    for(oI = 0; oI < nObs; oI++){
        tmp = allObs.at(oI).getBandIFreqMag();
        result(oI, 0) = tmp(0);
        result(oI, 1) = tmp(1);
    }
    return result; //.rows(0, cI - 1);
}


void hierarchicalLS::set_individualFreq(vec freq_){
    int oI;
    for(oI = 0; oI < nObs; oI++){
        allObs.at(oI).set_freq(freq_(oI));
    }
}

