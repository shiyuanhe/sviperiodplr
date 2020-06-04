#include "hierarchical.hpp"
#include <RcppArmadilloExtensions/sample.h>

hierarchicalPE::hierarchicalPE():
    viU_Theta(make_shared<mat>(1,1))
{
    haveFreq = false;
    haveOmega = false;
    haveNBand = false;
    havePLR = false;
    haveGamma = false;
    miniBatchSize = 8;
    //kernelScaling = 1;
}

void hierarchicalPE::compute(size_t maxIter, double scale1, double scale2)
{
    if(scale2 > 0) throw(runtime_error("scale2 should be negative!"));
    if(scale1 < 10) throw(runtime_error("scale1 is too small!"));

    init();
    //Rcpp::Rcout << "a" << std::endl;
    size_t  iterTotal = 0, startObsI = 0;
    // iter = 0, epochI, maxEpoch = 0;
    // maxEpoch = std::ceil(nObs / static_cast<double>(miniBatchSize) );

    Progress p(maxIter, true); //maxEpoch * maxIter
    initTrace(maxIter); //maxEpoch * maxIter

    double kappa;
    while (iterTotal < maxIter)
    {
        // startObsI = 0; // mini-batch starting index
        //for(epochI = 0; epochI < maxEpoch; epochI++){
        kappa = pow(scale1 + iterTotal, scale2);
        // Rcpp::Rcout << startObsI << std::endl;
        computeMiniBatch(startObsI);
        //Rcpp::Rcout << "Mini computeMiniBatch" << std::endl;

        updateOmega(kappa);
        //Rcpp::Rcout << "Mini updateOmega" << std::endl;

        updateGamma(kappa);
        //Rcpp::Rcout << "Mini updateGamma" << std::endl;

        updateAlpha(kappa);
        //Rcpp::Rcout << "Mini updateAlpha" << std::endl;

        // updateKernel(kappa);
        //Rcpp::Rcout << "Mini updateKernel" << std::endl;

        //startObsI += miniBatchSize;
        recordTrace(iterTotal);
        //Rcpp::Rcout << "Mini recordTrace" << std::endl;

        iterTotal++;
        p.increment();
        //}
        //iter++;
    }
}


void hierarchicalPE::initTrace(size_t containerSize ){
    alphaTrace = zeros<mat>(PLR_BASIS_DOF * nBand, containerSize);
    gammaTrace = mat(nBand, containerSize);

    kernelTrace = zeros<mat>(kernelBetaDim * nBand, containerSize);
    negloglikeTrace = zeros<vec>(containerSize);
}

void hierarchicalPE::recordTrace(size_t iter){
    span sp;
    for(size_t bI = 0; bI < nBand; bI++){
        sp = span(bI * PLR_BASIS_DOF, (bI + 1) * PLR_BASIS_DOF - 1);
        alphaTrace(sp, iter) = plrVec.at(bI)->viU_alpha;
        gammaTrace(bI, iter) = paramGamma.at(bI).viU_gamma;
        if(kernelType == 1){
            sp = span(bI * kernelBetaDim, (bI + 1) * kernelBetaDim - 1);
            kernelTrace(sp, iter) = coreCov.at(bI)->getCovBeta();
        }
    }
    negloglikeTrace(iter) = sum(negLogLike);
}

void hierarchicalPE::init(){
    check_ready();
    nObs = allObs.size();
    samplingSequence = arma::regspace<uvec>(0, nObs-1);
    //omp_set_num_threads(NUM_PARA_THREAD);
    rx.resize(NUM_PARA_THREAD);
    for(size_t i = 0; i < NUM_PARA_THREAD; i++){
        rx.at(i).seed(100);
        rx.at(i).split(NUM_PARA_THREAD, i);
    }

    // Initialize container for mini-batch result
    viU_thetaVecCollect = mat(BASIS_DOF * nBand, miniBatchSize);
    viS_thetaVecCollect.resize(miniBatchSize);
    viS_thetaVecSum = mat(BASIS_DOF * nBand, BASIS_DOF * nBand, fill::zeros);
    freqSample = zeros<vec>(miniBatchSize);
    negLogLike = freqSample;
    plBasisCollect = mat(PLR_BASIS_DOF, miniBatchSize);
    if(kernelType == 1)
        kernelGradientCollect = mat(kernelBetaDim * nBand, miniBatchSize);
    else
        kernelGradientCollect = zeros<mat>(1,1);

}

void hierarchicalPE::check_ready(){

    if(!haveFreq)
        throw(runtime_error("Frequency not set yet!"));
    if(!havePLR)
        throw(runtime_error("Please set PLR by set_alpha()!"));
    if(!haveNBand)
        throw(runtime_error("Please set band number by sen_nBand()!"));
    if(!haveOmega)
        throw(runtime_error("Omega Prior not set yet!"));
    if(!haveGamma)
        throw(runtime_error("Gamma prior not set yet!"));
}

void hierarchicalPE::computeMiniBatch(size_t startIndex){
    int mI;
    double f;

    update_viU_Theta();

    size_t rank = 0; //omp_get_thread_num();

    // coreCov is not thread safe!!
    // pragma omp parallel for schedule(static)
    uvec currentI = Rcpp::RcppArmadillo::sample(samplingSequence, miniBatchSize, false);
    // Rcpp::Rcout << currentI << std::endl << std::endl;

    for (mI = 0; mI < miniBatchSize; mI++)
    {

        size_t sampleI = currentI(mI); //startIndex + mI; //uniformInt(samplingMechanism);
        //sampleI %= nObs;
        // Rcpp::Rcout << sampleI << std::endl;
        // Only need to update global Theta precision matrix.
        // The PLR is updated by pointer automatically.

        //Rcpp::Rcout << 2 << std::endl;
        allObs.at(sampleI).compute();
        //Rcpp::Rcout << 3 << std::endl;

        // The VI expectation is difficult to evaluate.
        // We use sampled f, theta, and PL basis.
        f = allObs.at(sampleI).sampleFreqComputeExpect(rx.at(rank));

        //Rcpp::Rcout << 4 << std::endl;
        viU_thetaVecCollect.col(mI) = allObs.at(sampleI).get_viUthetaVec();

        //Rcpp::Rcout << 5 << std::endl;
        viS_thetaVecCollect.at(mI) = allObs.at(sampleI).get_viSthetaVec();

        //Rcpp::Rcout << 6 << std::endl;
        freqSample(mI) = f;

        // Assume PL basis is the same for all bands.
        plBasisCollect.col(mI) = plrVec.at(0)->computePLBasis(f);
        //Rcpp::Rcout << 7 << std::endl;

        // Update for kernel gradient only
        // if(kernelType == 1){
            // kernelGradientCollect.col(mI) = allObs.at(sampleI).get_covDeriv();
            //negLogLike(mI) = allObs.at(sampleI).gpNegativeLogLik();
        //}
    }
    // ofstream myfile;
    // myfile.open("/Users/heshiyuan/Data/Mira/freq.txt", ios::app);

    viS_thetaVecSum.zeros();
    //vec tmp;
    for(mI = 0; mI < miniBatchSize; mI++)
        viS_thetaVecSum += viS_thetaVecCollect.at(mI);

    // myfile.close();

}

// Pull the global precicesion matrix of m and beta into one matrix
// Because we have employed pointer, the value in each star obs
// will be automatically updated.
void hierarchicalPE::update_viU_Theta(){
    (*viU_Theta).submat(Jset, Jset) = paramOmega.viU_Omega;
    for(size_t bI = 0; bI < nBand; bI++)
        (*viU_Theta)(BASIS_DOF*bI, BASIS_DOF*bI) = paramGamma.at(bI).viU_gamma;
}


// Compute the posterior parameters and update Omega
void hierarchicalPE::updateOmega(double kappa)
{
    mat xi_Omega1, xi_Omega1_Update;
    double xi_Omega2;
    xi_Omega1 = paramOmega.init_xi_Omega1();
    xi_Omega2 = paramOmega.init_xi_Omega2();

    // Extract the component for beta.
    xi_Omega1_Update = viU_thetaVecCollect.rows(Jset);
    xi_Omega1_Update = xi_Omega1_Update * xi_Omega1_Update.t();
    xi_Omega1_Update += viS_thetaVecSum.submat(Jset, Jset);

    xi_Omega1 += (-0.5 * nObs / static_cast<double>(miniBatchSize))
        * xi_Omega1_Update;
    xi_Omega2 += 0.5 * nObs;
    paramOmega.updateOmega(xi_Omega1, xi_Omega2, kappa);

}




// Compute expected value of the posteriror canonical parameter
void hierarchicalPE::updateGamma(double kappa)
{
    size_t bI;
    double  tgamma1;

    double xi_gamma1, xi_gamma2;
    vec viU_alpha;
    mat viS_alpha;
    rowvec tmp;
    for (bI = 0; bI < nBand; bI++)
    {
        //if(bI == 0) continue;
        xi_gamma1 = paramGamma.at(bI).init_xi_gamma1();
        xi_gamma2 = paramGamma.at(bI).init_xi_gamma2();
        viU_alpha = plrVec.at(bI)->viU_alpha;
        viS_alpha = plrVec.at(bI)->viS_alpha;

        tmp = viU_thetaVecCollect.row(bI * BASIS_DOF) - viU_alpha.t() * plBasisCollect;
        tgamma1 = sum(tmp % tmp);
        tgamma1 += dot(plBasisCollect * plBasisCollect.t(), viS_alpha);
        tgamma1 += viS_thetaVecSum(bI * BASIS_DOF, bI * BASIS_DOF);

        xi_gamma1 -= tgamma1 * 0.5 * nObs / static_cast<double>(miniBatchSize);
        xi_gamma2 += 0.5 * nObs;
        paramGamma.at(bI).update_gamma(xi_gamma1, xi_gamma2, kappa);
    }

}

void hierarchicalPE::updateAlpha(double kappa)
{
    mat xi_alpha1;
    vec xi_alpha2;
    double ctmp;
    // mat cholEta1;
    // vec viU_alpha;
    size_t bI;
    for (bI = 0; bI < nBand; bI++)
    {
        // if(bI == 0) continue;

        // init parameters with prior only
        xi_alpha1 = plrVec.at(bI)->init_eta_alpha1();
        ctmp =  paramGamma.at(bI).viU_gamma
            * nObs / static_cast<double>(miniBatchSize);
        xi_alpha1 += -0.5 * ctmp * (plBasisCollect * plBasisCollect.t());

        xi_alpha2 = plrVec.at(bI)->init_eta_alpha2(); //deltaBar * alphaBar.at(bI);
        ctmp = paramGamma.at(bI).viU_gamma * nObs
            / static_cast<double>(miniBatchSize);
        xi_alpha2 += ctmp * (plBasisCollect
                                 * viU_thetaVecCollect.row(bI * BASIS_DOF).t());

        plrVec.at(bI)->updateAlpha(xi_alpha1, xi_alpha2, kappa);

        // Additional to compute the expectation of xi_alpha
        // if(bI == 1){
        //     ofstream myfile;
        //     myfile.open("/Users/heshiyuan/Data/Mira/alpha.txt", ios::app);
        //     cholEta1 = chol((-2.0) * xi_alpha1, "lower"); //
        //     viU_alpha = solve(trimatl(cholEta1), xi_alpha2);
        //     viU_alpha = solve(trimatu(cholEta1.t()), viU_alpha);
        //     for(int k = 0; k < 3; k++){
        //         myfile << viU_alpha(k) << " ";
        //     }
        //     myfile << std::endl;
        //     myfile.close();
        //
        // }
    }

}


void hierarchicalPE::updateKernel(double kappa){
    int bI;
    vec gradBetaAll, bandCoef;
    gradBetaAll = mean(kernelGradientCollect, 1) * nObs;
    // Rcpp::Rcout << norm(gradBetaAll) << std::endl;
    gradBetaAll *= kernelScaling;
    for(bI = 0; bI < nBand; bI++){
        bandCoef = gradBetaAll.rows(bI * kernelBetaDim, (bI + 1) * kernelBetaDim -1);
        coreCov.at(bI)->update_beta(bandCoef, kappa);
    }
}


