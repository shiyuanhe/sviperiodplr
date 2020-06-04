#ifndef LOMBSCARGLECLASS
#define LOMBSCARGLECLASS
#include "include.hpp"
#include "covSquaredExp.hpp"


class lombScargle
{
public:
    // Data: list of time, mag, uncertainty, # obs for each band
    lombScargle(vector<vec> tt_, vector<vec> yy_,
              vector<vec> ss_, vec nSample_);

    // fmin, fmax and nSeq
    void set_freqParam(double fmin_, double fmax_, double nSeq_){
        fmin = fmin_;
        fmax = fmax_;
        nSeq = nSeq_;

        fdelta = (fmax - fmin) / nSeq;
    }

    double fitBandI(double freq, size_t, vec &);
    vec getBandIFreqMag();


    // core of computation
    double compute();
    void set_freq(double freq);
    void setGPBand(int cB){currentB = cB;}
    double gpNegativeLogLik(vec betaNew);
    vec get_covDeriv(vec betaNew);


protected:

    // Data
    size_t nBand;
    vec nSample;
    // mjd, mag, sigma^2
    vector<vec> tt, yy, ss;
    vector<vec> phase, residual;
    covSquaredExp coreCov;
    int currentB;

    // min and max value of frequency.
    // fhat: estimated value.
    double fmin, fmax, fdelta, fhat;
    int nSeq;
    vec periodogram;

    vec aveMag; // average magnitude for each band
};


#endif
