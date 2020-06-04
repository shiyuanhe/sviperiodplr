#ifndef SQUARED_EXP_HEADER
#define SQUARED_EXP_HEADER

#include "RcppArmadillo.h"
#include <cmath>
using namespace std;
using namespace arma;


// All data sample in matrix columns

class squaredExp{
public:
    squaredExp(){
        beta = vec(2, fill::zeros);
        theta = exp(beta);
    }

    void set_beta(vec beta_){
        beta = beta_;
        theta = exp(beta);
    }


    //k = 0 covariance exp(beta[0]) * exp( - distSq / exp(beta[1]) )
    //k = 1 partial deriv wrt beta1, the value is the same as k = 0
    //k = 2 partial deriv wrt beta2,
    mat kernelVar(const vec & X, const int k = 0);

    mat kernelCov(const vec & X1, const vec & X2);

    // The first order partial direvative w.r.t X1(cj,i)
    //mat kernelCov_partial_cj(const mat & X1, const mat & X2, int cj);

    // The second order partial direvative w.r.t X1(cj,i)
    //mat kernelCov_partial2O_cj(const mat & X1, const mat & X2, int cj);

    // The second order partial direvative w.r.t X1(cj,i)
    //mat kernelCov_partial2OCross(const mat & X1, const mat & X2,
                                 //int c1, int c2);


    inline double compute_kij(const double & X1colI, const double & X2colJ){
        double k_ij;
        k_ij = X1colI - X2colJ;
        k_ij = theta(0) * exp( - k_ij * k_ij / theta(1));
        return k_ij;
    }

private:
    vec theta, beta;
};


#endif
