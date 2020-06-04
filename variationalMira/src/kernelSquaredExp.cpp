#include "kernelSquaredExp.hpp"


//k = 0 covariance exp(beta[0]) * exp( - distSq / exp(beta[1]) )
//k = 1 partial deriv wrt beta1, the value is the same as k = 0
//k = 2 partial deriv wrt beta2,
mat squaredExp::kernelVar(const vec & X, const int k){
    double distSq, k_ij;
    int i, j, nSample = X.n_elem;
    mat kMat(nSample, nSample, fill::zeros);

    for(i = 0; i < nSample; i++)
        for(j = i; j < nSample; j++){
            k_ij = compute_kij(X(i), X(j));
            if(k == 2){
                distSq = X(i) -  X(j);
                distSq = distSq * distSq;
                k_ij *= distSq / theta(1);
            }
            kMat(i,j) = k_ij;
            kMat(j,i) = k_ij;
        }

        return kMat;
}


// Compute the squared exponential kernel between X1 and X2
mat squaredExp::kernelCov(const vec & X1, const vec & X2){
    int i, j, nSample1, nSample2;
    nSample1 =  X1.n_elem;
    nSample2 = X2.n_elem;
    mat kMat(nSample1, nSample2, fill::zeros);

    for(i = 0; i < nSample1; i++)
        for(j = 0; j < nSample2; j++)
            kMat(i,j) = compute_kij(X1(i), X2(j));

    return kMat;
}

// // The first order partial direvative w.r.t X1(cj,i)
// mat squaredExp::kernelCov_partial_cj(const mat & X1, const mat & X2, int cj){
//     double  k_ij;
//     int i, j, nSample1, nSample2;
//     nSample1 =  X1.n_cols;
//     nSample2 = X2.n_cols;
//     mat kMat(nSample1, nSample2, fill::zeros);
//
//     for(i = 0; i < nSample1; i++)
//         for(j = 0; j < nSample2; j++){
//             k_ij = compute_kij(X1.col(i), X2.col(j));
//             k_ij *= -2.0 / theta(1) * (X1(cj,i) - X2(cj, j));
//             kMat(i,j) = k_ij;
//         }
//
//     return kMat;
// }
//
// // The second order partial direvative w.r.t X1(cj,i)
// mat squaredExp::kernelCov_partial2O_cj(const mat & X1, const mat & X2, int cj){
//     double k_ij, e_ij;
//     int i, j, nSample1, nSample2;
//     nSample1 =  X1.n_cols;
//     nSample2 = X2.n_cols;
//     mat kMat(nSample1, nSample2, fill::zeros);
//
//     for(i = 0; i < nSample1; i++)
//         for(j = 0; j < nSample2; j++){
//             k_ij = compute_kij(X1.col(i), X2.col(j));
//             e_ij = -2.0 / theta(1) * (X1(cj, i) - X2(cj, j));
//             k_ij = k_ij * e_ij * e_ij - 2.0 / theta(1) * k_ij;
//             kMat(i,j) = k_ij;
//         }
//     return kMat;
// }
//
//
// // The second order partial direvative w.r.t X1(cj,i)
// mat squaredExp::kernelCov_partial2OCross(const mat & X1, const mat & X2,
//                                          int c1, int c2){
//     double k_ij, e_ij, d_ij;
//     int i, j, nSample1, nSample2;
//     nSample1 = X1.n_cols;
//     nSample2 = X2.n_cols;
//     mat kMat(nSample1, nSample2, fill::zeros);
//
//     for(i = 0; i < nSample1; i++)
//         for(j = 0; j < nSample2; j++){
//             k_ij = compute_kij(X1.col(i), X2.col(j));
//             e_ij = -2.0 / theta(1) * (X1(c1, i) - X2(c1, j));
//             d_ij = -2.0 / theta(1) * (X1(c2, i) - X2(c2, j));
//             kMat(i,j) = k_ij * e_ij * d_ij;
//         }
//     return kMat;
// }
//








