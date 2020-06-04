#include "miraStar.hpp"




miraClass::miraClass(vector<vec> tt_, vector<vec> yy_,
                     vector<vec> ss_, vec nSample_)
    :component_basis(basisClass(nSample_))
{
    haveFreq = false;
    havePLR = false;
    haveData = false;
    haveThetaPrior = false;
    haveCov = false;

    tt = tt_;
    yy = yy_;
    ss = ss_;
    nSample = nSample_;
    nBand = nSample.size();
    haveData = true;

    thetaPriorMean = zeros<vec>(nBand * BASIS_DOF);
}



void miraClass::init()
{
    if(!haveData) throw(runtime_error("Please provide data first!"));
    if(!havePLR) throw(runtime_error("The PLR not found!"));
    if(!haveFreq) throw(runtime_error("Frequency range not specified!"));
    if(!haveThetaPrior) throw(runtime_error("Theta Prior not found!"));
    if(!haveCov) throw(runtime_error("Cov Class not found!"));

    thetaPriorMean = zeros<vec>(nBand * BASIS_DOF);
}


void miraClass::compute()
{
    init();
    //freq = component_freq.fetch_reset();
    vector<covSquaredExp> coreCovConcur;
    shared_ptr<covSquaredExp> tmpPtr;

    for(size_t bI = 0; bI < coreCov.size(); bI++){
        tmpPtr = std::dynamic_pointer_cast<covSquaredExp>(coreCov.at(bI));
        coreCovConcur.push_back( *tmpPtr );
    }

    mat gfA = zeros<mat>(PLR_BASIS_DOF, PLR_BASIS_DOF);
    for (size_t bI = 0; bI < nBand; bI++){
        // terms related to the PL relation
        gfA += (corePLR.at(bI)->viS_alpha ) *
            (*Omega_theta)(bI * BASIS_DOF, bI * BASIS_DOF);
    }

    int nfI, nfTotal;
    nfTotal = component_freq.get_nfSeq();
// #ifdef _OPENMP
// omp_set_num_threads(NUM_PARA_THREAD);
// #endif
//#pragma omp parallel for schedule(static) \
//   firstprivate(thetaPriorMean, component_basis, coreCovConcur) \
//   private(nfI) \
//    shared(gfA, tt, ss, yy, corePLR, Omega_theta)
    for(nfI = 0; nfI < nfTotal; nfI++)
    {
        double freq;
        freq = component_freq.get_freq_atI(nfI);
        // ************************************ //
        // computeUpdateCore(freq);
        // ************************************ //

        // ************************************ //
        // set_freq(double ff_)
        // ************************************ //
        size_t bI;
        vec tmpT;

        for (bI = 0; bI < nBand; bI++){
            thetaPriorMean(bI * BASIS_DOF ) =
                corePLR.at(bI)->individual_priorMean(freq);

            if (notEnoughBandSample(nSample.at(bI)) ) continue;
            tmpT = tt.at(bI) * freq;
            component_basis.generate_basis(tmpT, bI);
            coreCovConcur.at(bI).setupData(tmpT, ss.at(bI));
        }



        // ************************************ //
        // update_theta()
        // ************************************ //

        // natural parameters for the posterior
        mat xi1_theta;
        vec xi2_theta;

        mat tmpC;
        vec tmpY;
        // prior value
        xi1_theta = (*Omega_theta);
        xi2_theta = (*Omega_theta) * thetaPriorMean;

        // data update to the canonical parameters
        for (bI = 0; bI < nBand; bI++)
        {
            if (notEnoughBandSample(nSample.at(bI)) ) continue;
            tmpC = component_basis.basisMats.at(bI);

            // Sigma = LL^T, to get L^{-1}C and L^{-1}y
            tmpC = coreCovConcur.at(bI).halfActOnMat(tmpC);
            tmpY = coreCovConcur.at(bI).halfActOnVector(yy.at(bI));

            auto colSeq = component_basis.get_thetaSpan(bI);
            xi1_theta(colSeq, colSeq) += tmpC.t()* tmpC;
            xi2_theta(colSeq) += tmpC.t()* tmpY;
        }

        xi1_theta *= (-0.5);
        component_basis.update_theta(xi1_theta, xi2_theta);



        // ************************************ //
        // qfValue = update_Qf();
        // ************************************ //

        double cqf;
        mat tmpMat;
        vec tmpy;

        // two terms related the periodical basis.
        // -0.25 * (eta2)^T(eta1)^{-1}(eta2)
        cqf = component_basis.qf1;

        // The entropy - 0.5 *logdet(-eta1)
        // q(f) \propto exp(g(f) + Entropy)
        cqf += component_basis.qf2;

        // A = \sum_b E(gamma_b)* E(\alpha_b \alpha_b^T)
        // mat gfA = zeros<mat>(PLR_BASIS_DOF, PLR_BASIS_DOF);
        for (bI = 0; bI < nBand; bI++){
            // terms related to the PL relation
            // gfA += (corePLR.at(bI)->viS_alpha ) *
            //     (*Omega_theta)(bI * BASIS_DOF, bI * BASIS_DOF);

            // Remaining sample term
            // y^T\Sigma^{-1}y + logdet Sigma
            if (notEnoughBandSample(nSample.at(bI))) continue;
            tmpy = coreCovConcur.at(bI).halfActOnVector(yy.at(bI));
            cqf += -0.5 * dot(tmpy, tmpy);
            cqf += -0.5 * coreCovConcur.at(bI).logdet();
        }

        // The term <didi^T, A>
        vec di = corePLR.at(0)->computePLBasis(freq);
        cqf -=  0.5 * dot(di, gfA*di);
        cqf -= 0.5 * dot(thetaPriorMean,
                         (*Omega_theta) * thetaPriorMean);

        // ************************************ //
        // component_freq.set_current_qfValue(qfValue );
        // ************************************ //

        component_freq.set_qfValue_atI(cqf, nfI);

        // freq = component_freq.fetch_next();
    }
    // After a whole range of computation,
    // normalize to proper density.
    component_freq.normalize_freq();
}

// Set new frequency, and compute q(theta|f)
void miraClass::computeUpdateCore(double freq){
    set_freq(freq);
    update_theta();
}

// Use the frequency already selected
void miraClass::computeUpdateCore(){
    set_freq(freqValue);
    update_theta();
}

// double miraClass::getMiraModelEvidence(){
//     double res = component_freq.getMiraModelEvidence();
//     int dfT = sum(nSample) + (*Omega_theta).n_cols;
//     res += -(dfT) / 2.0 * log(2*datum::pi);
//     res -= 1.5 * (*Omega_theta).n_cols;
//
//     double logDetTheta, signT;
//     log_det(logDetTheta, signT, (*Omega_theta) );
//     res += logDetTheta * 0.5;
//     return  res;
// }

// Set the value of the frequency, update basis, and update covariance mat.
void miraClass::set_freq(double ff_)
{
    freqValue = ff_;
    size_t bI;
    vec tmpT;
    for (bI = 0; bI < nBand; bI++){
        thetaPriorMean(bI * BASIS_DOF ) =
            corePLR.at(bI)->individual_priorMean(freqValue);

        if (notEnoughBandSample(nSample.at(bI)) ) continue;
        tmpT = tt.at(bI) * freqValue;
        component_basis.generate_basis(tmpT, bI);
        coreCov.at(bI)->setupData(tmpT, ss.at(bI));
    }
    //Rcpp::Rcout << thetaPriorMean << endl;
}

// Compute the posterior canonical parameters for theta.(mean mag, beta)
void miraClass::update_theta()
{
    // natural parameters for the posterior
    mat xi1_theta;
    vec xi2_theta;

    size_t bI;
    mat tmpC;
    vec tmpY;
    // prior value
    xi1_theta = (*Omega_theta);
    xi2_theta = (*Omega_theta) * thetaPriorMean;
    //Rcpp::Rcout << "f" << std::endl;

    // data update to the canonical parameters
    for (bI = 0; bI < nBand; bI++)
    {
        if (notEnoughBandSample(nSample.at(bI)) ) continue;
        tmpC = component_basis.basisMats.at(bI);

        // Sigma = LL^T, to get L^{-1}C and L^{-1}y
        tmpC = coreCov.at(bI)->halfActOnMat(tmpC);
        tmpY = coreCov.at(bI)->halfActOnVector(yy.at(bI));
        //Rcpp::Rcout << "a" << std::endl;

        auto colSeq = component_basis.get_thetaSpan(bI);
        xi1_theta(colSeq, colSeq) += tmpC.t()* tmpC;
        xi2_theta(colSeq) += tmpC.t()* tmpY;
    }
    //Rcpp::Rcout << "b" << std::endl;

    xi1_theta *= (-0.5);
    component_basis.update_theta(xi1_theta, xi2_theta);
    //Rcpp::Rcout << "d" << std::endl;

}


// Compute q(f_i) for current f_i
double miraClass::update_Qf()
{
    size_t bI;
    double cqf;
    mat tmpMat;
    vec tmpy;

    // two terms related the periodical basis.
    // -0.25 * (eta2)^T(eta1)^{-1}(eta2)
    cqf = component_basis.qf1;

    // The entropy - 0.5 *logdet(-eta1)
    // q(f) \propto exp(g(f) + Entropy)
    cqf += component_basis.qf2;

    // A = \sum_b E(gamma_b)* E(\alpha_b \alpha_b^T)
    mat gfA = zeros<mat>(PLR_BASIS_DOF, PLR_BASIS_DOF);
    for (bI = 0; bI < nBand; bI++){
        // terms related to the PL relation
        gfA += (corePLR.at(bI)->viS_alpha ) *
            (*Omega_theta)(bI * BASIS_DOF, bI * BASIS_DOF);

        // Remaining sample term
        // y^T\Sigma^{-1}y + logdet Sigma
        if (notEnoughBandSample(nSample.at(bI))) continue;
        tmpy = coreCov.at(bI)->halfActOnVector(yy.at(bI));
        cqf += -0.5 * dot(tmpy, tmpy);
        cqf += -0.5 * coreCov.at(bI)->logdet();
    }

    // The term <didi^T, A>
    vec di = corePLR.at(0)->computePLBasis(freqValue);
    cqf -=  0.5 * dot(di, gfA*di);
    cqf -= 0.5 * dot(thetaPriorMean,
                     (*Omega_theta) * thetaPriorMean);

    return cqf;
}


//Sample one frequency value via VI q(f),
// and compute the expecation of q(theta | f)
double miraClass::sampleFreqComputeExpect(trng::yarn2 &r)
{
    double freq = component_freq.sample_freq(r);

    //Rcpp::Rcout << "Get Freq " << freq << std::endl;
    computeUpdateCore(freq);

    // ofstream myfile;
    // myfile.open("/Users/heshiyuan/Data/Mira/prior1.txt", ios::app);
    // myfile << freq << " " << thetaPriorMean(3) << " " << corePLR.at(1)->viU_alpha(0) << std::endl;
    // myfile.close();

    return freq;
}


// Sinusoid basis coef for all bands.
vec miraClass::get_viUBetaAll(){
    return component_basis.get_viUBetaAll();

}

// Get covariance derivative for the negative log-likelihood
// Only to squared exponetial class
vec miraClass::get_covDeriv(){
    vec result, theta1M;
    mat covCoreMat, theta2M;

    size_t bI, subI;
    mat tmpC;
    vec tmpY;
    // prior value
    result = zeros<vec>(nBand * 3);

    subI = 0;
    // data update to the canonical parameters
    for (bI = 0; bI < nBand; bI++)
    {
        if (notEnoughBandSample(nSample.at(bI)) ) continue;
        tmpC = component_basis.basisMats.at(bI);

        // Sigma = LL^T, to get L^{-1}C and L^{-1}y
        tmpC = coreCov.at(bI)->inverseOnMat(tmpC);
        tmpY = coreCov.at(bI)->inverseOnVector(yy.at(bI));

        auto colSeq = component_basis.get_thetaSpan(bI);
        theta1M = component_basis.viU_theta(colSeq);
        theta2M = component_basis.viS_theta(colSeq, colSeq);
        theta2M += theta1M * theta1M.t();

        covCoreMat = (tmpY - 2.0 * (tmpC * theta1M) ) * tmpY.t();
        covCoreMat +=  tmpC * theta2M * tmpC.t();
        covCoreMat -=  coreCov.at(bI)->getInverse();

        result(subI + 0) = dot(covCoreMat, coreCov.at(bI)->covPartial_j(0));
        result(subI + 1) = dot(covCoreMat, coreCov.at(bI)->covPartial_j(1));
        result(subI + 2) = dot(covCoreMat, coreCov.at(bI)->covPartial_j(2));

        subI += 3;
    }
    result *= -0.5;

    return result;
}

// This function serves to provide initial theta estimation.
// The expected negative loglikelihood under variational parameter
double miraClass::gpNegativeLogLik(){
    double result = 0;
    int bI;
    mat tmpC, theta2M;
    vec tmpy, theta1M;

    for (bI = 0; bI < nBand; bI++){
        if (notEnoughBandSample(nSample.at(bI)) ) continue;
        tmpC = component_basis.basisMats.at(bI);
        tmpC = coreCov.at(bI)->halfActOnMat(tmpC);

        tmpy = yy.at(bI);// - tmpC * theta1M;
        tmpy = coreCov.at(bI)->halfActOnVector(tmpy);

        auto colSeq = component_basis.get_thetaSpan(bI);
        theta1M = component_basis.viU_theta(colSeq);

        theta2M = component_basis.viS_theta(colSeq, colSeq);
        theta2M += theta1M * theta1M.t();

        result += 0.5 * dot(tmpy, tmpy);
        result -= dot(tmpC.t() * tmpy, theta1M);
        //Rcpp::Rcout << result << " ";
        result += 0.5 * dot(tmpC.t() * tmpC, theta2M);
        //Rcpp::Rcout << result << " ";
        result += 0.5 * coreCov.at(bI)->logdet();
        //Rcpp::Rcout << result << " " << std::endl;

    }

    return result;
}



