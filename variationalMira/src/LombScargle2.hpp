#ifndef LOMBSCARGLECLASST2
#define LOMBSCARGLECLASST2
#include "include.hpp"
#include "covSquaredExp.hpp"


class lombScargleT2
{
public:
    // Data: list of time, mag, uncertainty, # obs for each band
    lombScargleT2(vector<vec> tt_, vector<vec> yy_,
              vector<vec> ss_, vec nSample_, vec prior_);

    // fmin, fmax and nSeq
    void set_freq(double fmin_, double fmax_, double nSeq_){
        fmin = fmin_;
        fmax = fmax_;
        nSeq = nSeq_;

        fdelta = (fmax - fmin) / nSeq;
    }



    void setGPBand(int cB){currentB = cB;}
    double gpNegativeLogLik(vec betaNew);
    vec get_covDeriv(vec betaNew);

    void computeCovMat(const vec &phase, const vec &beta);

    mat covPartial_j(int j, const vec &phase, const vec & beta);

protected:

    // Data
    size_t nBand;
    vec nSample;
    // mjd, mag, sigma^2
    vector<vec> tt, yy, ss;
    int currentB;
    vec PriorMean;
    double PriorVar;

    // min and max value of frequency.
    double fmin, fmax, fdelta;
    int nSeq;

    vec aveMag; // average magnitude for each band

    //
    mat Kmat, Sigma, SigmaChol, SigmaInv;
    squaredExp kernelObj;

};


#endif
