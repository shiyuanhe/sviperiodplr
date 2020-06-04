#include "kernelPCA.hpp"

// Input: sorted time points X
// OUtput: the covariance matrix
mat pcaKernel::kernelVar(const vec & X){
    int bI, eI, mI, pcaI;
    mat covMat(X.n_elem, X.n_elem, fill::zeros);
    bI = 0;
    double diffT;
    for(mI = 0; mI < X.n_elem; mI++){
        while(X(bI) - X(mI) < -0.5) bI++;
        eI = bI;
        while((eI < X.n_elem) && (X(eI) - X(mI) < 0.5)){
            diffT = X(eI) - X(mI) +0.5;
            pcaI = floor(diffT * 1000);
            pcaI = (pcaI<0)?0:pcaI;
            pcaI = (pcaI>999)?999:pcaI;
            covMat(mI, eI) = beta(pcaI);
            eI++;
        }
    }
    covMat = (covMat + covMat.t())/2;
    vec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, covMat);
    eigval(find(eigval < 0)).zeros();
    eigval = sqrt(eigval);
    eigvec.each_row() %= eigval.t();
    covMat = eigvec * eigvec.t();
    return covMat;
}

