#ifndef PLRELATION_CLASS
#define PLRELATION_CLASS

const int PLR_BASIS_DOF = 3;

#include "include.hpp"

class plRelationClass{
public:
    // Alpha (multivariate Gaussian)
    double deltaBar;
    vec alphaBar;
    mat eta_alpha1, eta_alpha1_gradient;
    vec eta_alpha2, eta_alpha2_gradient;
    vec viU_alpha;
    mat viS_alpha, viEvvT; // EvvT: noncentral second moment

    vec computePLBasis(double f){
        vec tmpbpl = ones<vec>(PLR_BASIS_DOF);
        tmpbpl(1) = log10(f);
        tmpbpl(2) = tmpbpl(1) * tmpbpl(1);
        //tmpbpl(3) = tmpbpl(2) * tmpbpl(1);
        return tmpbpl;
    }

    // prior mean for the average magnitude for one Mira
    double individual_priorMean(double f){
        vec plbasis = computePLBasis(f);
        return dot(plbasis, viU_alpha);
    }

    // prior var for the average magnitude for one Mira
    double individual_priorVar(double f){
        vec plbasis = computePLBasis(f);
        return dot(plbasis, viS_alpha * plbasis);
    }

    mat init_eta_alpha1(){
        return (eye<mat>(PLR_BASIS_DOF, PLR_BASIS_DOF) * deltaBar * (-0.5));
    }

    vec init_eta_alpha2(){
        return (deltaBar * alphaBar);
    }



    void set_alpha(vec alphaBar_, double deltaBar_){
        alphaBar = alphaBar_;
        deltaBar = deltaBar_;

        eta_alpha1 = init_eta_alpha1();
        eta_alpha2 = init_eta_alpha2();

        // viU_alpha = alphaBar;
        // viS_alpha = eye<mat>(PLR_BASIS_DOF, PLR_BASIS_DOF) / deltaBar;
        // viEvvT = viS_alpha + viU_alpha * viU_alpha.t();

        compute_viUS_alpha();

        eta_alpha1_gradient = zeros<mat>(PLR_BASIS_DOF, PLR_BASIS_DOF);
        eta_alpha2_gradient = zeros<vec>(PLR_BASIS_DOF);
    }

    // Set posterior (initialization)
    void  set_posterior(mat xi_eta1, vec xi_eta2){
        eta_alpha1 = xi_eta1;
        eta_alpha2 = xi_eta2;
        compute_viUS_alpha();

        eta_alpha1_gradient = zeros<mat>(PLR_BASIS_DOF, PLR_BASIS_DOF);
        eta_alpha2_gradient = zeros<vec>(PLR_BASIS_DOF);
    }


    // Input: posterior parameters and step size kappa
    void updateAlpha(const mat & xi_alpha1, const vec & xi_alpha2, double kappa)
    {
        // stochastic gradient descent with momentum.
        eta_alpha1_gradient = 0.8 * eta_alpha1_gradient +
            (xi_alpha1 - eta_alpha1);
        eta_alpha2_gradient = 0.8 * eta_alpha2_gradient +
            (xi_alpha2 - eta_alpha2);
        // eta_alpha1_gradient =  (xi_alpha1 - eta_alpha1);
        // eta_alpha2_gradient =  (xi_alpha2 - eta_alpha2);
        eta_alpha1 +=  kappa * eta_alpha1_gradient;
        eta_alpha2 +=  kappa * eta_alpha2_gradient;
        compute_viUS_alpha();
    }

    // compute the conditional expection from the canonical parameters,
    // eta_alpha1, eta_alpha2
    void compute_viUS_alpha()
    {
        // Compute Exptation and Variance of alphab
        // mat cholEta1;
        // cholEta1 = chol((-2.0) * eta_alpha1, "lower"); //
        // viU_alpha = solve(trimatl(cholEta1), eta_alpha2);
        // viU_alpha = solve(trimatu(cholEta1.t()), viU_alpha);
        //
        // viS_alpha = solve(trimatl(cholEta1),
        //                   eye<mat>(eta_alpha1.n_rows,
        //                            eta_alpha1.n_cols));
        // viS_alpha = solve(trimatu(cholEta1.t()), viS_alpha);

        viU_alpha = -0.5 * solve(eta_alpha1, eta_alpha2);
        viS_alpha = -0.5 * inv(eta_alpha1);
        viEvvT = viS_alpha + viU_alpha * viU_alpha.t();
    }
};

typedef shared_ptr<plRelationClass> plr_ptr;

#endif
