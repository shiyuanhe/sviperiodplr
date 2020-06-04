#ifndef COMP_BASIS
#define COMP_BASIS

#include "include.hpp"

// The first is the constant 1 for average magnitude.
// The rest is abitrary basis.

const int BASIS_DOF = 3;


class basisClass{
public:
    // theta is the combination of (m, beta)
    // The m is the average magnitude.
    // The beta is the sinusoid basis coef.
    vector<mat> basisMats;
    size_t nBand;
    mat viS_theta; // Var under VI q(.)
    vec viU_theta; // Expectation under VI q(.)
    mat  eta1_theta; // canonical Param 1
    vec  eta2_theta; // canonical Param 2

    // two terms related to the computation of q(f)
    // qf1 = -0.25 (eta2)^T(eta1)^{-1}(eta2)
    // qf2 = -0.5 * log det (-eta1)
    double qf1, qf2;

    // input: the number of sample for each band
    basisClass(vec nSample_){
        nBand = nSample_.size();
        basisMats.resize(nBand);
        size_t i;
        for(i = 0; i < nBand; i++){
            basisMats.at(i) = ones<mat>(nSample_(i), BASIS_DOF);
        }
    }


    void generate_basis(const vec & tSeq, int bandI){
        basisMats.at(bandI).col(1) = cos(2*PI*tSeq);
        basisMats.at(bandI).col(2) = sin(2*PI*tSeq);
    }


    // get the fitted value for band bI
    vec get_fittedValue(size_t bI){
        auto colSel = get_thetaSpan(bI);
        return basisMats.at(bI) * viU_theta(colSel);
    }

    void update_theta(const mat & eta1_theta_, const vec & eta2_theta_){
        eta1_theta = eta1_theta_;
        eta2_theta = eta2_theta_;
        compute_viUS_theta();
    }

    void compute_viUS_theta(){
        // Compute Exptation and Variance of theta
        mat cholEta1;
        cholEta1 = chol( -2 *eta1_theta, "lower");

        // post mean, step 1
        // (-2 eta1)^{-1} eta2 = (LL^T)^{-1} eta_2
        viU_theta = solve(trimatl(cholEta1), eta2_theta);

        // use mid-value to compute qf1
        // -0.25 (eta2)^T(eta1)^{-1} eta2
        // = 0.5(eta2)^T(-2*eta1)^{-1} eta2
        // = 0.5(eta2)^T(LL^T)^{-1} eta2
        qf1  =  0.5 * sum(viU_theta % viU_theta);

        // post mean, step 2
        viU_theta = solve(trimatu(cholEta1.t()), viU_theta);

        // post variance
        viS_theta = -0.5 * inv(eta1_theta);

        // Compute the entropy term, ignoring additive constants
        // -0.5 log |-2eta1| = -0.5 log |L|^2 = -log|L|
        qf2 = - sum(log(cholEta1.diag()));
    }

    // The mean mag and its uncertainty for all bands.
    vec get_estimatedMeanMag(){
        vec result(2 * nBand, fill::zeros);
        size_t bI, dI;
        for(bI = 0; bI < nBand; bI++){
            dI = bI * BASIS_DOF;
            result(2 * bI) = viU_theta(dI);
            result(2 * bI + 1) = sqrt(viS_theta(dI, dI));
        }
        return result;
    }

    // The first two sinusoid coefs for all bands, beta
    vec get_viUBetaAll(){
        vec result(2 * nBand);
        size_t bI;
        for(bI = 0; bI < nBand; bI++){
            result(2*bI) = viU_theta(BASIS_DOF*bI + 1);
            result(2*bI + 1) = viU_theta(BASIS_DOF*bI + 2);
        }
        return result;
    }



    // // get the span for theta for band bI
    inline span get_thetaSpan(size_t bI){
        return span(BASIS_DOF * bI,
                    BASIS_DOF * (bI + 1)  - 1);
    }


};


// mat get_viSBeta(size_t bI){
//     span colSel = get_betaSpan(bI);
//     return viS_theta(colSel, colSel);
// }

// Get the span for beta for band bI.
// This excludes the first constant term.
// inline span get_betaSpan(size_t bI){
//     return span(BASIS_DOF * bI + 1,
//                 BASIS_DOF * (bI + 1) - 1);
// }

#endif



