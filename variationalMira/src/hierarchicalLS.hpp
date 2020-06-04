#ifndef HIERLOMBSCARGLECLASS
#define HIERLOMBSCARGLECLASS
#include "LombScargle.hpp"


class hierarchicalLS
{
public:
    hierarchicalLS();
    void push_data(vector<vec>, vector<vec>, vector<vec>, vec);

    void compute();

    void set_globalFreq(double fmin_, double fmax_, size_t nfseq_){
        fMin = fmin_;
        fMax = fmax_;
        nfSeq = nfseq_;
    }

    void setGPBand(int currentB){
        int oI;
        for(oI = 0; oI < nObs; oI++)
            allObs.at(oI).setGPBand(currentB);
        gpUse = ones<vec>(nObs);
    }

    void set_scaling(double ss_){
        grad_scaling = ss_;
    }
    double gpNegativeLogLik(vec newBeta);
    vec get_covDeriv(vec betaNew);

    mat getFreqMag();

    void set_gpUse(vec u_){
        // if gpUse(i)==1, the i-th sample will be used for
        // compute the gaussian process parameter.
        // This is used to remove outliers.
        gpUse = u_;
    }


    void set_individualFreq(vec freq_);
private:

    // Data
    vector<lombScargle> allObs;
    size_t nObs;
    double fMin, fMax;
    size_t nfSeq;
    // mat PLValues;
    double grad_scaling;
    vec gpUse;

};

#endif
