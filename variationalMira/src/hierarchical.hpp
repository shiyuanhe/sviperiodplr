#ifndef HIER_MODEL_CLASS
#define HIER_MODEL_CLASS
#include "miraStar.hpp"


class hierarchicalPE
{
public:
    hierarchicalPE();
    void push_data(vector<vec>, vector<vec>, vector<vec>, vec);

    void compute(size_t maxIter, double scale1, double scale2);

    // set hyperparameters
    void set_globalFreq(double, double, size_t);
    void set_nBand(int);
    void set_Omega(mat, int);
    void set_alpha(vector<vec>, double);
    void set_gamma(vec, double);

    void set_eta_Omega(mat, double);
    void set_eta_alpha(vector<mat>, vector<vec>);
    void set_eta_gamma(vec, vec);

    void set_kernel(vector<vec> beta, int typeI);
    void set_kernelScaling(double);
    mat get_qf(size_t miraI);
    mat get_fit(size_t miraI, vec tSeq);

    // set initial value
    // void set_eta_alpha(vector<mat> e1, vector<vec> e2){
    // }
    // void set_eta_Omega(mat e1, double e2){
    // }
    // void set_eta_gamma(vec e1, vec e2){
    // }

    // get fitted parameters
    Rcpp::List get_eta_alpha1();
    Rcpp::List get_eta_alpha2();
    mat get_eta_Omega1();
    double get_eta_Omega2();
    vec get_eta_gamma1();
    vec get_eta_gamma2();

    mat get_gammaTrace(){return gammaTrace;}
    mat get_alphaTrace(){return alphaTrace;}
    mat get_kernelTrace(){return kernelTrace;}
    vec get_negloglikeTrace() { return negloglikeTrace;}
    mat get_PLValues();
    mat get_ThetaAll();


    // Test for empirical GP
    double gpObjective(vec beta);
    vec gpGradient(vec beta);

private:
    void init();
    void check_ready();
    void computeMiniBatch(size_t startI);
    void updateOmega(double kappa);
    void updateGamma(double kappa);
    void updateAlpha(double kappa);
    void update_viU_Theta();
    void updateKernel(double kappa);

    void recordTrace(size_t );
    void initTrace(size_t );


    // Data
    vector<miraClass> allObs;
    size_t nObs, nBand;
    double fMin, fMax;
    size_t nfSeq;
    // mat PLValues;


    // bool set status
    bool haveFreq, haveOmega, haveNBand,
    havePLR, haveGamma, haveKernel;

    // Parameters. Xi: canonical param of full conditional
    // eta: canoncial parameters of Variational Inference
    // viU and viS are the variational mean and variational variance
    vector<plr_ptr> plrVec;
    paramOmegaClass paramOmega;
    vector<paramGammaClass> paramGamma;
    // Degree of freedom for PL relation basis and expansion basis

    // Parameters for GP kernel
    vec kernelBeta;
    vector<cov_ptr> coreCov;
    int kernelBetaDim, kernelType;

    // Theta for local update
    mat_ptr viU_Theta; // joint information of viU_Omega and viU_gamma
    // Iset: index for m in theta
    // Jset: index for beta in theta
    uvec  Jset;

    // Sampling
    //trng::mt19937_64 samplingMechanism;
    //trng::uniform_int_dist uniformInt;
    //Rcpp::NumericVector nObsVec;
    vector<trng::yarn2> rx;
    uvec samplingSequence; //0,1,...,sampleSize (nObs) - 1

    double kernelScaling;
    // Trace
    mat gammaTrace;
    mat alphaTrace;
    mat kernelTrace;
    vec negloglikeTrace;

    // Mini-batch Collection
    size_t miniBatchSize;
    mat viU_thetaVecCollect; // minibatch obs in columns
    mat kernelGradientCollect;
    mat plBasisCollect;
    vector<mat> viS_thetaVecCollect;
    mat viS_thetaVecSum;
    vec freqSample;
    vec negLogLike;


};

#endif
