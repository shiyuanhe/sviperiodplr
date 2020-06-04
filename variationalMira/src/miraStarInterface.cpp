#include "miraStar.hpp"



void miraClass::set_globalPL(vector<plr_ptr> plrVec_)
{
    if(!haveData) throw(runtime_error("Please provide data first!"));
    if(plrVec_.size() != nBand)
        throw(runtime_error("PLR band number does not match obs band number!"));
    corePLR = plrVec_;
    havePLR = true;
}



void miraClass::Rset_alpha(vector<vec> alphaBar_, double deltaBar_)
{
    if(!haveData) throw(runtime_error("Please provide data first!"));
    corePLR.resize(nBand);
    size_t bI;
    for(bI = 0; bI < nBand; bI++){
        corePLR.at(bI) = make_shared<plRelationClass>();
        corePLR.at(bI)->set_alpha(alphaBar_.at(bI), deltaBar_);
    }

    havePLR = true;
}


void miraClass::set_globalTheta(mat_ptr Omega_theta_)
{
    if(!haveData) throw(runtime_error("Please provide data first!"));
    Omega_theta = Omega_theta_;
    haveThetaPrior = true;
}



void miraClass::Rset_globalCov(vector<vec> beta, int typeI){
    if(!haveData) throw(runtime_error("Please provide data first!"));
    // if(coreCov_.size() != nBand)
    //     throw(runtime_error("PLR band number does not match obs band number!"));
    coreCov.resize(nBand);
    int i;
    if(typeI == 1){
        for(i = 0; i < nBand; i++)
            coreCov.at(i) = make_shared<covSquaredExp>(beta.at(i));
    }else{
        for(i = 0; i < nBand; i++)
            coreCov.at(i) = make_shared<covPCA>(beta.at(i));
    }

    haveCov = true;
}


void miraClass::Rupdate_globalCov(vec betaNew){
    int i, betaDim;
    betaDim = coreCov.at(0)->getBetaDim();
    for(i = 0; i < nBand; i++){
        arma::span sp = span(i*betaDim, (i+1)*betaDim - 1);
        coreCov.at(i)->set_beta(betaNew.rows(sp));
    }
}

void miraClass::set_globalCov(vector<cov_ptr> covPtrVec_){
    if(!haveData) throw(runtime_error("Please provide data first!"));
    if(covPtrVec_.size() != nBand)
        throw(runtime_error("Cov vector does not match obs band number!"));

    coreCov = covPtrVec_;
    haveCov = true;
}




void miraClass::set_globalFreq(double freqMin_, double freqMax_,
                               size_t nSeqf_)
{
    if(!haveData) throw(runtime_error("Please provide data first!"));
    component_freq.set_freq_param(freqMin_, freqMax_, nSeqf_);
    haveFreq = true;
}




// Return f, f_sigma, and m, m_sigma for all bands.
// f is the value of posterior mode.
vec miraClass::getPLValues(){
    vec res1, res2;
    res1 = component_freq.get_estimatedPeriod();
    set_freq(1/res1(0));
    update_theta();


    res2 = component_basis.get_estimatedMeanMag();
    return join_cols(res1, res2);
}



mat miraClass::get_fit(vec tSeqNew){
    // double fHat;
    // uword imax = qfSeqLong.index_max();
    // fHat = 1/fSeqLong(imax);

    double fHat;

    fHat = 1/component_freq.get_estimatedPeriod()(0);
    set_freq(fHat);


    update_theta();
    vec phaseLong = tSeqNew * fHat;

    mat basisLong = ones<mat>(tSeqNew.size(), 3);
    basisLong.col(1) = cos(2*PI*phaseLong);
    basisLong.col(2) = sin(2*PI*phaseLong);

    mat result(tSeqNew.size(), nBand, fill::zeros);
    size_t bI;
    vec residOld, residNew;
    mat covNewOld;

    for(bI = 0; bI < nBand; bI++){
        if (notEnoughBandSample(nSample.at(bI)) ) continue;

        auto colSel = component_basis.get_thetaSpan(bI);
        result.col(bI) = basisLong *
            component_basis.viU_theta(colSel);

        residOld = yy.at(bI) - component_basis.get_fittedValue(bI);

        // compute residuals
        residOld = coreCov.at(bI)->inverseOnVector(residOld);
        covNewOld = coreCov.at(bI)->computeCrossCov(phaseLong);
        residNew = covNewOld * residOld;
        // predict residuals on the dense sequence

        result.col(bI) += residNew;
    }

    return result;
}


vec miraClass::get_theta(){
    double fHat;

    fHat = 1/component_freq.get_estimatedPeriod()(0);
    set_freq(fHat);
    update_theta();

    vec result(BASIS_DOF * nBand, fill::zeros);
    size_t bI;

    for(bI = 0; bI < nBand; bI++){
        if (notEnoughBandSample(nSample.at(bI)) ) continue;
        auto colSel = component_basis.get_thetaSpan(bI);
        result(colSel) =
            component_basis.viU_theta(colSel);
    }

    return result;
}
