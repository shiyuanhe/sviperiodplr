#include "hierarchical.hpp"

// This is a set to optimize GP kernel parameter.
// The objective function, gradient function is provided.
// A whole round of iteration is required beforehand.
double hierarchicalPE::gpObjective(vec beta){
    double result = 0;
    size_t mI, i, betaDim;
    vec tmpBeta;

    betaDim = coreCov.at(0)->getBetaDim();
    for(i = 0; i < nBand; i++){
        tmpBeta = beta.rows(i * betaDim,
                            (i+1) * betaDim - 1);
        coreCov.at(i)->set_beta(tmpBeta);
    }

    for (mI = 0; mI < nObs; mI++)
    {
        allObs.at(mI).computeUpdateCore();
        result += allObs.at(mI).gpNegativeLogLik();
    }

    return result;
}

vec hierarchicalPE::gpGradient(vec beta){
    vec result(size(beta), fill::zeros);

    size_t mI, i, betaDim;
    vec tmpBeta;

    betaDim = coreCov.at(0)->getBetaDim();
    for(i = 0; i < nBand; i++){
        tmpBeta = beta.rows(i * betaDim,
                            (i+1) * betaDim - 1);
        coreCov.at(i)->set_beta(tmpBeta);
    }

    for (mI = 0; mI < nObs; mI++)
    {
        allObs.at(mI).computeUpdateCore();
        result += allObs.at(mI).get_covDeriv();
    }

    return result;
}







