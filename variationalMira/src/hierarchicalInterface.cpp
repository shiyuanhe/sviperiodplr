#include "hierarchical.hpp"


// 1. Set Freq
// 2. Set nBand
// 3. Set Omega
// 4. Set PLR(alpha)
// 5. Set Gamma
// 6. Set Cov Mat(kernel)
// 7. Push Data

// Preparation

void hierarchicalPE::push_data(vector<vec> tt, vector<vec> yy,
                               vector<vec> ss, vec nSample)
{
    if(!haveFreq)
        throw(runtime_error("Frequency not set yet!"));
    if(!havePLR)
        throw(runtime_error("Please set PLR by set_alpha()!"));
    if(!haveNBand)
        throw(runtime_error("Please set band number by set_nBand()!"));
    if(!haveKernel)
        throw(runtime_error("Please set kernel parameters by set_nBand()!"));
    if(nSample.size() != nBand)
        throw(runtime_error("Band number does not match!"));

    miraClass tmp(tt, yy, ss, nSample);
    tmp.set_globalFreq(fMin, fMax, nfSeq);
    tmp.set_globalTheta(viU_Theta);
    tmp.set_globalPL(plrVec);
    tmp.set_globalCov(coreCov);
    allObs.push_back(tmp);
}



void hierarchicalPE::set_globalFreq(double fmin_, double fmax_, size_t nfseq_){
    fMin = fmin_;
    fMax = fmax_;
    nfSeq = nfseq_;
    haveFreq = true;
}

void hierarchicalPE::set_nBand(int nBand_)
{
    nBand = nBand_;
    //Iset = ones<uvec>(nBand);
    Jset = ones<uvec>( (BASIS_DOF - 1) * nBand);

    size_t bI, cI, posI;
    for (bI = 0; bI < nBand; bI++)
    {
        //Iset(bI) = BASIS_DOF * bI;
        for (cI = 1; cI < BASIS_DOF; cI++){
            posI = (BASIS_DOF - 1) * bI + cI - 1;
            Jset(posI) = BASIS_DOF * bI + cI;
        }

    }
    viU_Theta = make_shared<mat>();
    *viU_Theta = zeros<mat>(BASIS_DOF*nBand, BASIS_DOF*nBand);
    haveNBand = true;
}



void hierarchicalPE::set_Omega(mat OmegaBar_, int nBar_)
{
    paramOmega.set_omega(OmegaBar_, nBar_);
    haveOmega = true;
}



void hierarchicalPE::set_alpha(vector<vec> alphaBar_,
                               double deltaBar_)
{
    if(!haveNBand)
        throw(runtime_error("Please specify the number of bands first!") );
    if(alphaBar_.size() != nBand)
        throw(runtime_error("Band number does not match!"));

    plrVec.resize(nBand);
    size_t bI;
    for(bI = 0; bI < nBand; bI++){
        plrVec.at(bI) = make_shared<plRelationClass>();
        plrVec.at(bI)->set_alpha(alphaBar_.at(bI), deltaBar_);
    }

    havePLR = true;
}


// Each beta in the vector corresponds to one band kernel.
// typeI = 1: squared exponential kernel
//            pcaKernel, otherwise

void hierarchicalPE::set_kernel(vector<vec> beta, int typeI)
{
    if(!haveNBand)
        throw(runtime_error("Please specify the number of bands first!") );
    if(beta.size() != nBand)
        throw(runtime_error("Band number does not match!"));

    coreCov.resize(nBand);
    int i;
    if(typeI == 1){
        for(i = 0; i < nBand; i++)
            coreCov.at(i) = make_shared<covSquaredExp>(beta.at(i));
    }else{
        for(i = 0; i < nBand; i++)
            coreCov.at(i) = make_shared<covPCA>(beta.at(i));
    }

    kernelBetaDim = coreCov.at(0)->getBetaDim();
    kernelType = typeI;
    haveKernel = true;
}

void hierarchicalPE::set_kernelScaling(double s_){
    kernelScaling = s_;
}


void hierarchicalPE::set_gamma(vec gammaBar_, double rBar_)
{
    if(!haveNBand)
        throw(runtime_error("Please specify the number of bands first!") );
    if(gammaBar_.size() != nBand)
        throw(runtime_error("Band number does not match!"));

    paramGamma.resize(nBand);
    size_t bI;
    for(bI = 0; bI < nBand; bI++)
        paramGamma.at(bI).set_gamma(gammaBar_(bI), rBar_);

    haveGamma = true;
}

// Set Init value after setting prior.
// Their accumulated gradient in the past iteration
// will be set to zero.
void hierarchicalPE::set_eta_Omega(mat xi_eta1,
                                   double xi_eta2){
    if(!haveOmega)
        throw(runtime_error("Omega Prior not set yet!"));

    paramOmega.set_posterior(xi_eta1, xi_eta2);
}

void hierarchicalPE::set_eta_alpha(vector<mat> xi_eta1,
                                   vector<vec> xi_eta2){
    if(!havePLR)
        throw(runtime_error("Please set PLR by set_alpha()!"));
    size_t bI;
    for(bI = 0; bI < nBand; bI++)
        plrVec.at(bI)->set_posterior(xi_eta1.at(bI),
                  xi_eta2.at(bI));

}

void hierarchicalPE::set_eta_gamma(vec xi_eta1,
                                   vec xi_eta2){
    if(!haveGamma)
        throw(runtime_error("Gamma prior not set yet!"));
    size_t bI;
    for(bI = 0; bI < nBand; bI++)
        paramGamma.at(bI).set_posterior(xi_eta1(bI),
                      xi_eta2(bI));

}




mat hierarchicalPE::get_PLValues(){
    int mI;
    mat PLValues = zeros<mat>(2*nBand + 2, nObs);
    Progress p(nObs, true);
    update_viU_Theta();
    // pragma omp parallel for schedule(static)
    for (mI = 0; mI < nObs; mI++)
    {
        //allObs.at(mI).set_globalTheta(viU_Theta);
        allObs.at(mI).compute();
        PLValues.col(mI) = allObs.at(mI).getPLValues();
        p.increment();
    }
    return PLValues;
}

// Get the theta vector for all Miras.
mat hierarchicalPE::get_ThetaAll(){
    int mI;
    mat thetaAll = zeros<mat>(3*nBand, nObs);
    for (mI = 0; mI < nObs; mI++)
        thetaAll.col(mI) = allObs.at(mI).get_theta();
    return thetaAll;
}

// Sample and Compute
// Extract gradient
mat hierarchicalPE::get_qf(size_t miraI){
    // allObs.at(miraI).set_globalBeta(viU_Theta);
    // allObs.at(miraI).set_globalPL(viU_alpha, viS_alpha);
    // allObs.at(miraI).compute();
    return allObs.at(miraI).get_qfSeq();
}

mat hierarchicalPE::get_fit(size_t miraI, vec tSeq){
    // allObs.at(miraI).set_globalBeta(viU_Theta);
    // allObs.at(miraI).set_globalPL(viU_alpha, viS_alpha);
    // allObs.at(miraI).compute();
    return allObs.at(miraI).get_fit(tSeq);
}


// get fitted parameters
Rcpp::List hierarchicalPE::get_eta_alpha1(){
    Rcpp::List res(nBand);
    for(size_t i = 0; i < nBand; i++){
        res.at(i) = Rcpp::wrap(plrVec.at(i)->eta_alpha1);
    }
    return res;
}

Rcpp::List hierarchicalPE::get_eta_alpha2(){
    Rcpp::List res(nBand);
    for(size_t i = 0; i < nBand; i++){
        res.at(i) = Rcpp::wrap(plrVec.at(i)->eta_alpha2);
    }
    return res;
}

mat hierarchicalPE::get_eta_Omega1(){
    return paramOmega.eta_Omega1;
}

double hierarchicalPE::get_eta_Omega2(){
    return paramOmega.eta_Omega2;
}

vec hierarchicalPE::get_eta_gamma1(){
    vec eta_gamma1(nBand);
    for(size_t i = 0; i < nBand; i++)
        eta_gamma1(i) = paramGamma.at(i).eta_gamma1;
    return eta_gamma1;
}

vec hierarchicalPE::get_eta_gamma2(){
    vec eta_gamma2(nBand);
    for(size_t i = 0; i < nBand; i++)
        eta_gamma2(i) = paramGamma.at(i).eta_gamma2;
    return eta_gamma2;
}

