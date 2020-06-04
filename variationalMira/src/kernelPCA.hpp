#ifndef PCA_KERNEL_HEADER
#define PCA_KERNEL_HEADER

#include "RcppArmadillo.h"
#include <cmath>
using namespace std;
using namespace arma;


// Kernel is trained from FPCA.
// The cov between 0.5 and [0,1] is measured,
// which is recorded on a 1000 grid points `beta`.

class pcaKernel{
public:
    pcaKernel(){
    }

    void set_beta(vec beta_){
        beta = beta_;
    }


    // Input: sorted time points X
    // OUtput: the covariance matrix
    mat kernelVar(const vec & X);


private:
    vec beta;
};


#endif
