#ifndef COVSQUAREDEXPCLASS
#define COVSQUAREDEXPCLASS

#include "covClass.hpp"

class covSquaredExp:
    public covarianceClass
{
public:
    covSquaredExp(vec covBeta_);
    virtual ~covSquaredExp(){}

    vec inverseOnVector(const vec & v);
    mat inverseOnMat(const mat & m);
    mat getInverse();

    void setupData(const vec & tSeq_, const vec & sigmaSqSeq_);
    mat computeCrossCov(const vec & tSeq);

    mat covPartial_j(int j);
    void update_beta(vec gradBeta, double kappa);
    void set_beta(vec betaNew);
private:
    int nSample;
    vec tSeq;
    mat Kmat, SigmaInv;
    squaredExp kernelObj;
    vec gradBetaCumu;
};


#endif


