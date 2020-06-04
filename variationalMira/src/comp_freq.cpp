#include "comp_freq.hpp"

// #ifndef RSAMPLING
// #define RSAMPLING
// #include <RcppArmadilloExtensions/sample.h>
// #endif

// Sample according to VI q(f)
double freqClass::sample_freq(trng::yarn2 &r){
    if(!normalized) throw(runtime_error("density not normalized yet!"));
    trng::discrete_dist dist(qfSeqLong.begin(), qfSeqLong.end());
    size_t x = dist(r);
    // Rcpp::Rcout << "Sample " << size(fSeqLong) << " " <<
    //     size(qfSeqLong) << " " << x << std::endl;

    // return Rcpp::RcppArmadillo::sample(fSeqLong, 1, false, qfSeqLong)(0);
    return fSeqLong(x);
}

void freqClass::set_freq_param(double freqMin_, double freqMax_,
                    size_t nSeqf_)
{
    n_fSeq = nSeqf_;
    freqMin = freqMin_;
    freqMax = freqMax_;
    fSeq = linspace(freqMin, freqMax, n_fSeq);
    qfSeq = ones<vec>(n_fSeq);
    fSeqLong = linspace(freqMin, freqMax, 4 * n_fSeq);
    qfSeqLong = ones<vec>(4 * n_fSeq);
    deltaFreq = (freqMax - freqMin) / (4.0 * n_fSeq);

    //modelSel = qfSeq;
    //modelSelLong = qfSeqLong;
    normalized = false;
}


// Return the estimated period (posterior mode)
// and its uncertainty
vec freqClass::get_estimatedPeriod(){
    if(!normalized) throw(runtime_error("density not normalized yet!"));
        vec result(2, fill::zeros);
    vec pSeq = 1/fSeqLong;
    double pHat, p2mHat;
    //pHat = sum(fSeqLong % qfSeqLong) * deltaFreq;
    pHat = sum(pSeq % qfSeqLong) * deltaFreq;
    p2mHat = sum(pSeq % pSeq % qfSeqLong) * deltaFreq;

    uword imax = qfSeqLong.index_max();
    result(0) = 1/fSeqLong(imax);
    result(1) = p2mHat - pHat * pHat;
    return result;
}




// Linear interpolation for more accurate result.
void freqClass::normalize_freq(){
    qfSeq = qfSeq - max(qfSeq);
    interp1(fSeq, qfSeq, fSeqLong, qfSeqLong, "*linear");
    qfSeqLong = exp(qfSeqLong);
    qfSeqLong /= sum(qfSeqLong) * deltaFreq;
    normalized = true;
}



// double freqClass::getMiraModelEvidence(){
//     if(!normalized) throw(runtime_error("density not normalized yet!"));
//     interp1(fSeq, modelSel, fSeqLong, modelSelLong, "*linear");
//     return sum(modelSelLong % qfSeqLong) * deltaFreq;
// }






