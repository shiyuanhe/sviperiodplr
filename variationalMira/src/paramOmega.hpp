#ifndef OMEGA_CLASS
#define OMEGA_CLASS

#include "include.hpp"

class paramOmegaClass{
public:
    // Omega (Wishart distribution)
    mat OmegaBar, OmegaBarInv; // hyperparam
    double nBar;                  //hyperparam
    double eta_Omega2, eta_Omega2_gradient;
    mat eta_Omega1, eta_Omega1_gradient;
    mat viU_Omega;
    int p; // matrix dimesion

    paramOmegaClass(){
        // init p
    }

    mat init_xi_Omega1(){
        return (-0.5 * nBar * OmegaBarInv);  // store for computing saving
    }

    double init_xi_Omega2(){
        return 0.5 * (nBar - p - 1); // store for computing saving
    }



    // Set prior
    void set_omega(const mat & OmegaBar_, const double & nBar_){
        nBar = nBar_;
        OmegaBar = OmegaBar_;
        OmegaBarInv = inv(OmegaBar);
        p = OmegaBar.n_cols;


        eta_Omega1 = init_xi_Omega1();
        eta_Omega2 = init_xi_Omega2();

        eta_Omega1_gradient = zeros<mat>(p,p);
        eta_Omega2_gradient = 0.0;
        compute_viUS_Omega();
    }

    // Directly set posterior. (initialization)
    void set_posterior(mat xi_eta1,
                       double xi_eta2){
        eta_Omega1 = xi_eta1;
        eta_Omega2 = xi_eta2;
        eta_Omega1_gradient = zeros<mat>(p,p);
        eta_Omega2_gradient = 0.0;
        compute_viUS_Omega();
    }

    // Input: posterior parameters and step size kappa
    void updateOmega(const mat & xi_Omega1, const double & xi_Omega2, double kappa)
    {
        eta_Omega1_gradient = 0.8 * eta_Omega1_gradient +  (xi_Omega1 - eta_Omega1);
        eta_Omega2_gradient = 0.8 * eta_Omega2_gradient +  (xi_Omega2 - eta_Omega2);
        // eta_Omega1_gradient =  (xi_Omega1 - eta_Omega1);
        // eta_Omega2_gradient =  (xi_Omega2 - eta_Omega2);
        eta_Omega1 +=  kappa * eta_Omega1_gradient;
        eta_Omega2 +=  kappa * eta_Omega2_gradient;
        compute_viUS_Omega();
    }



    void compute_viUS_Omega()
    {
        viU_Omega = (2.0 * eta_Omega2 + p + 1.0) * (-0.5) *
            inv(eta_Omega1);
    }
};

typedef shared_ptr<plRelationClass> plr_ptr;

#endif
