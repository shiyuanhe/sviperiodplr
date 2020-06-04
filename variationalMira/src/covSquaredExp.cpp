#include "covSquaredExp.hpp"

covSquaredExp::covSquaredExp(vec covBeta_):
    covarianceClass(covBeta_){
    kernelObj.set_beta(covBeta.rows(1,2) );
    gradBetaCumu = vec(size(covBeta_), fill::zeros);

}

vec covSquaredExp::inverseOnVector(const vec & v){
    return SigmaInv * v;
}

mat covSquaredExp::inverseOnMat(const mat & m){
    return SigmaInv * m;
}

mat covSquaredExp::getInverse(){
    return SigmaInv;
}


void covSquaredExp::setupData(const vec & tSeq_,
                                     const vec & sigmaSqSeq_){
    Kmat = kernelObj.kernelVar(tSeq_, 0);

    Sigma = Kmat;
    Sigma.diag() += sigmaSqSeq_ + exp(covBeta(0));

    SigmaInv = inv(Sigma);
    SigmaChol = chol(Sigma, "lower"); // lower part??
    nSample = tSeq_.n_elem;
    tSeq = tSeq_;
}


// Compute the covariance matrix between the new time sequence tSeqNew
// and the original stored sequence tSeq.
mat covSquaredExp::computeCrossCov(const vec & tSeqNew){
    return kernelObj.kernelCov(tSeqNew, tSeq);
}

mat covSquaredExp::covPartial_j(int j){
    if(j==1)
        return Kmat;
    if(j==2)
        return kernelObj.kernelVar(tSeq, 2);

    mat result;
    result = eye<mat>(nSample, nSample);
    return exp(covBeta(0)) * result;

}


void covSquaredExp::update_beta(vec gradBeta, double kappa){
    gradBetaCumu = 0.8 * gradBetaCumu + gradBeta;
    covBeta -= kappa * gradBetaCumu;
    kernelObj.set_beta(covBeta.rows(1,2) );
}



void covSquaredExp::set_beta(vec betaNew){
    covBeta = betaNew;
    gradBetaCumu = vec(size(betaNew), fill::zeros);
    kernelObj.set_beta(covBeta.rows(1,2) );
}

