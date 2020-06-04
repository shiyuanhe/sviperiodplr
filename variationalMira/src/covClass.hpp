#ifndef COVARIANCE_CLASS
#define COVARIANCE_CLASS

#include "include.hpp"
#include "kernelPCA.hpp"
#include "kernelSquaredExp.hpp"

class covarianceClass{
public:
    covarianceClass(vec covBeta_){
        covBeta = covBeta_;
        betaDim = covBeta_.n_elem;
    }

    virtual ~covarianceClass(){}

    // Compute SigmaChol^{-1} v
    vec halfActOnVector(const vec & v){
        return solve(SigmaChol, v);
    }

    // Compute SigmaChol^{-1} m
    mat halfActOnMat(const mat & m){
        return solve(SigmaChol, m);
    }


    // Return the log-determinant of the covariance matrix
    double logdet() const{
        return 2 * sum(log(abs(SigmaChol.diag() )));
    }

    int getBetaDim() const { return betaDim; }
    vec getCovBeta() { return covBeta; }

    virtual vec inverseOnVector(const vec & v){
        throw(runtime_error("inverseOnVector not defined!"));
        return zeros<vec>(1);
    }

    virtual mat inverseOnMat(const mat & m){
        throw(runtime_error("inverseOnMat not defined!"));
        return zeros<mat>(1,1);
    }

    virtual mat getInverse(){
        throw(runtime_error("getInverse not defined!"));
        return zeros<mat>(1,1);
    }


    virtual void setupData(const vec & tSeq, const vec & sigmaSqSeq){
        throw(runtime_error("setupData not defined!"));
    }

    virtual mat computeCrossCov(const vec & tSeq){
        throw(runtime_error("computeCrossCov not defined!"));
    }


    virtual mat covPartial_j(int j){
        throw(runtime_error("covPartial_j not defined!"));
        return zeros<mat>(1,1);

    }

    virtual void update_beta(vec gradBeta, double kappa){
        throw(runtime_error("update_beta not defined!"));
    }

    virtual void set_beta(vec betaNew){
        throw(runtime_error("set_beta not defined!"));
    }


protected:
    // Sigma = SigmaChol * t(SigmaChol)
    mat Sigma, SigmaChol;
    vec covBeta;

    int betaDim; // the dimension of beta
};


typedef shared_ptr<covarianceClass> cov_ptr;

#endif
