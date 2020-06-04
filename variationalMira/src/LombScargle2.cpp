#include "LombScargle2.hpp"

lombScargleT2::lombScargleT2(vector<vec> tt_, vector<vec> yy_,
                         vector<vec> ss_, vec nSample_, vec Prior_){
    tt = tt_;
    yy = yy_;
    ss = ss_;
    nSample = nSample_;
    nBand = nSample.size();
    PriorMean = Prior_.rows(1, nBand);
    PriorVar = Prior_(0);
}




// Get covariance derivative for the negative log-likelihood
// Only to squared exponetial class
vec lombScargleT2::get_covDeriv(vec betaNew){


    mat covCoreMat;
    vec tmpY;
    vec result = zeros<vec>(3);

    if (notEnoughBandSample(nSample.at(currentB) - 2) ) return result;

    vec phase;
    double cf = fmin;
    vec tmpy;
    int fI;

    for(fI = 0; fI < nSeq; fI++){
        phase = tt.at(currentB) * cf;

        //coreCov.set_beta(betaNew);
        //coreCov.setupData(phase, ss.at(currentB));
        // compute kernel
        computeCovMat(phase, betaNew);

        // subtract prior mean
        tmpy = yy.at(currentB) - PriorMean(currentB);

        tmpy = SigmaInv * tmpy;
        covCoreMat = tmpy * tmpy.t();
        covCoreMat -=  SigmaInv;

        result(0) += dot(covCoreMat, covPartial_j(0, phase, betaNew));
        result(1) += dot(covCoreMat, covPartial_j(1, phase, betaNew));
        result(2) += dot(covCoreMat, covPartial_j(2, phase, betaNew));

        cf += fdelta;
    }

    result *= -0.5;
    return result / static_cast<double>(nSeq);
}

// This function serves to provide initial theta estimation.
// The expected negative loglikelihood under variational parameter
double lombScargleT2::gpNegativeLogLik(vec betaNew){
    double result = 0.0;
    if (notEnoughBandSample(nSample.at(currentB) - 2) ) return result;

    vec phase;
    double cf = fmin;
    vec tmpy;
    int fI;

    for(fI = 0; fI < nSeq; fI++){
        phase = tt.at(currentB) * cf;

        // compute kernel
        computeCovMat(phase, betaNew);

        // subtract prior mean
        tmpy = yy.at(currentB) - PriorMean(currentB);

        // compute half
        tmpy = solve(SigmaChol, tmpy);
        result += 0.5 * dot(tmpy, tmpy);
        result += 0.5 * 2 * sum(log(abs(SigmaChol.diag() ))); // add the log-determinant
        cf += fdelta;
    }


    return result / static_cast<double>(nSeq);
}


void lombScargleT2::computeCovMat(const vec &phase,
                                  const vec &beta){
    kernelObj.set_beta(beta.rows(1,2) );
    Sigma = kernelObj.kernelVar(phase, 0);
    Kmat = Sigma;
    Sigma.diag() += ss.at(currentB) + exp(beta(0));

    mat basisMat = ones<mat>(nSample(currentB), 3);
    basisMat.col(1) = cos(phase * 2 * 3.1415926);
    basisMat.col(2) = sin(phase * 2 * 3.1415926);
    Sigma +=   basisMat * basisMat.t() * PriorVar;

    SigmaInv = inv(Sigma);
    SigmaChol = chol(Sigma, "lower"); // correct??
}


mat lombScargleT2::covPartial_j(int j,
                                const vec &phase,
                                const vec & beta){
    if(j==1)
        return Kmat;
    if(j==2)
        return kernelObj.kernelVar(phase, 2);

    mat result;
    result = eye<mat>(nSample(currentB), nSample(currentB));
    return exp(beta(0)) * result;
}

