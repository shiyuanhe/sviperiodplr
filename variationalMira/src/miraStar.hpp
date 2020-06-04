#ifndef MIRA_CLASS_BASE
#define MIRA_CLASS_BASE
#include "include.hpp"
#include "comp_freq.hpp"
#include "comp_basis.hpp"
#include "plrelation.hpp"
#include "paramOmega.hpp"
#include "paramGamma.hpp"
#include "covPCA.hpp"
#include "covSquaredExp.hpp"


class miraClass
{
public:
    // Data: list of time, mag, uncertainty, # obs for each band
    miraClass(vector<vec> tt_, vector<vec> yy_,
              vector<vec> ss_, vec nSample_);

    // *** Set from C++ or R, both OK *****
    // fmin, fmax and nSeq
    void set_globalFreq(double, double, size_t);

    // *** Set from C++ *****
    // Mean of Theta under q( ), i.e.
    // mean of gamma and Omega
    void set_globalTheta(mat_ptr);
    // Set the mean and var of alpha under q(.)
    void set_globalPL(vector<plr_ptr> plrVec_);
    // Covariance class parameters from R
    void set_globalCov(vector<cov_ptr> covPtrVec_);

    // *** Set from R *****
    // Set PL relation from R
    void Rset_alpha(vector<vec> alphaBar_, double deltaBar_);
    // Set globalTheta from R
    void Rset_globalTheta(mat Omega_thetaR_){
        Omega_theta = make_shared<mat>(Omega_thetaR_);
        haveThetaPrior = true;
    }
    // Covariance class parameters from R
    void Rset_globalCov(vector<vec> beta, int typeI);
    void Rupdate_globalCov(vec betaNew);

    // core of computation
    void compute();
    void computeUpdateCore(double freq);
    void computeUpdateCore();

    //double getMiraModelEvidence();
    double sampleFreqComputeExpect(trng::yarn2 &r);

    vec get_viUBetaAll();

    double gpNegativeLogLik();
    vec get_covDeriv();
    //void lombScargle();


    // A group of functions for output
    mat get_fit(vec tSeq);
    vec get_theta();
    mat get_qfSeq(){ return component_freq.get_qfSeq(); }
    vec getPLValues();
    vec get_viUthetaVec(){ return component_basis.viU_theta;}
    mat get_viSthetaVec(){ return component_basis.viS_theta;}



protected:

    void init();
    // update the mean field family
    //void update_m(double stepsize);
    void update_theta();
    double update_Qf();
    // sample one f by q(f)
    //double sample_freq();
    // specify a value for f and compute basis
    void set_freq(double ff_);


    // Used to check if it has received global parameters
    bool haveFreq, havePLR,
    haveData, haveThetaPrior, haveCov;


    // Data
    size_t nBand;
    vec nSample;
    // mjd, mag, sigma^2
    vector<vec> tt, yy, ss;

    double freqValue;
    freqClass component_freq;
    basisClass component_basis;
    vec thetaPriorMean; // prior mean of theta_tilde

    // Global

    // Prior precision of theta_tilde.
    // It'ss the expected value given related q(.)
    mat_ptr Omega_theta;
    // the PL relation for all bands
    vector<plr_ptr> corePLR;
    vector<cov_ptr> coreCov;

};


#endif
