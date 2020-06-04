#ifndef COVPCACLASS
#define COVPCACLASS

#include "covClass.hpp"

class covPCA:
    public covarianceClass
{
public:
    covPCA(vec covBeta_): covarianceClass(covBeta_){
        covCore.set_beta(covBeta_.rows(1,1000) );
    }

    virtual ~covPCA(){}


    void setupData(const vec & tSeq_, const vec & sigmaSqSeq_){
        Sigma = covCore.kernelVar(tSeq_);
        //Sigma = zeros<mat>(tSeq.n_rows, tSeq.n_rows);
        Sigma.diag() += sigmaSqSeq_ + exp( covBeta(0) );
        SigmaChol = chol(Sigma, "lower");
    }



private:

    pcaKernel covCore;
};


#endif
