#ifndef GAMMA_CLASS
#define GAMMA_CLASS

#include "include.hpp"

class paramGammaClass{
public:
    // Gamma (Gamma distribution)
    double rbar;  // hyperparam
    double gammabar; // hyperparam
    double eta_gamma1, eta_gamma1_gradient;
    double eta_gamma2, eta_gamma2_gradient;
    double viU_gamma;// posterior expectation

    paramGammaClass(){
        // init p
    }

    double init_xi_gamma1(){
        return -rbar;  // store for computing saving
    }

    double init_xi_gamma2(){
        return (rbar * gammabar - 1.0); // store for computing saving
    }


    // Set prior value
    void set_gamma(double gammaBar_, double rBar_){
        gammabar = gammaBar_;
        rbar = rBar_;
        eta_gamma1 = init_xi_gamma1();
        eta_gamma2 = init_xi_gamma2();
        eta_gamma1_gradient = 0;
        eta_gamma2_gradient = 0;
        compute_viUS_gamma();
    }

    // Set posterior (initialization)
    void  set_posterior(double xi_eta1, double xi_eta2){
        eta_gamma1 = xi_eta1;
        eta_gamma2 = xi_eta2;
        eta_gamma1_gradient = 0;
        eta_gamma2_gradient = 0;
        compute_viUS_gamma();
    }

    // Input: posterior parameters and step size kappa
    void update_gamma(const double & xi_gamma1, const double & xi_gamma2, double kappa)
    {
        eta_gamma1_gradient = 0.8 * eta_gamma1_gradient +  (xi_gamma1 - eta_gamma1);
        eta_gamma2_gradient = 0.8 * eta_gamma2_gradient +  (xi_gamma2 - eta_gamma2);
        // eta_gamma1_gradient =  (xi_gamma1 - eta_gamma1);
        // eta_gamma2_gradient =  (xi_gamma2 - eta_gamma2);
        eta_gamma1 += kappa * eta_gamma1_gradient;
        eta_gamma2 += kappa * eta_gamma2_gradient;
        compute_viUS_gamma();
    }

    void compute_viUS_gamma()
    {
        viU_gamma = -(eta_gamma2 + 1) / eta_gamma1;
    }


};

typedef shared_ptr<plRelationClass> plr_ptr;

#endif
