#include "hierarchical.hpp"
#include "LombScargle.hpp"
#include "hierarchicalLS.hpp"

RCPP_MODULE(miraClass){

    Rcpp::class_<plRelationClass>("plRelationClass")
    .constructor()
    .method("set_alpha", &plRelationClass::set_alpha)
    .method("updateAlpha", &plRelationClass::updateAlpha)
    .method("compute_viUS_alpha", &plRelationClass::compute_viUS_alpha)
    .method("individual_priorMean", &plRelationClass::individual_priorMean)
    .method("computePLBasis", &plRelationClass::computePLBasis)
    .method("individual_priorVar", &plRelationClass::individual_priorVar)
    ;

    Rcpp::class_<basisClass>("basisClass")
        .constructor<vec>()
        .method("generate_basis", &basisClass::generate_basis)
        .method("get_fittedValue", &basisClass::get_fittedValue)
        .method("update_theta", &basisClass::update_theta)
        .method("get_estimatedMeanMag", &basisClass::get_estimatedMeanMag)
        .method("get_viUBetaAll", &basisClass::get_viUBetaAll)
    ;


    Rcpp::class_<miraClass>("miraClass")
        .constructor<vector<vec>, vector<vec>, vector<vec>, vec>()
        .method("compute", &miraClass::compute)
    //.method("computeUpdateCore", &miraClass::computeUpdateCore)
      .method("set_globalFreq", &miraClass::set_globalFreq)
      .method("Rset_globalCov", &miraClass::Rset_globalCov)
      .method("Rset_alpha", &miraClass::Rset_alpha)
      .method("Rset_globalTheta", &miraClass::Rset_globalTheta)
      .method("Rupdate_globalCov", &miraClass::Rupdate_globalCov)
      .method("get_qfSeq", &miraClass::get_qfSeq)
      .method("getPLValues", &miraClass::getPLValues)
      .method("get_viUBetaAll", &miraClass::get_viUBetaAll)
      .method("get_fit", &miraClass::get_fit)
      .method("get_covDeriv", &miraClass::get_covDeriv)
      .method("gpNegativeLogLik", &miraClass::gpNegativeLogLik)
    ;

    // Rcpp::class_<miraStarOutlier>("outlierClass")
    //     .derives<miraClass>("miraClass")
    //     .constructor<vector<vec>, vector<vec>, vector<vec>, vec>()
    //     .method("Rset_globalPrior", &miraStarOutlier::Rset_globalPrior)
    //     .method("compute2", &miraStarOutlier::compute)
    //     .method("objFun", &miraStarOutlier::objFun)
    //     .method("gradFun", &miraStarOutlier::gradFun)
    //     .method("hessFun", &miraStarOutlier::hessFun)
    // ;

    Rcpp::class_<hierarchicalPE>("hierarchicalPE")
        .constructor()
        .method("compute", &hierarchicalPE::compute)
        .method("set_globalFreq", &hierarchicalPE::set_globalFreq)
        .method("set_Omega", &hierarchicalPE::set_Omega)
        .method("set_nBand", &hierarchicalPE::set_nBand)
        .method("set_alpha", &hierarchicalPE::set_alpha)
        .method("set_gamma", &hierarchicalPE::set_gamma)
        .method("set_kernel", &hierarchicalPE::set_kernel)
        .method("set_kernelScaling", &hierarchicalPE::set_kernelScaling)

    .method("set_eta_Omega", &hierarchicalPE::set_eta_Omega)
    .method("set_eta_alpha", &hierarchicalPE::set_eta_alpha)
    .method("set_eta_gamma", &hierarchicalPE::set_eta_gamma)

    // .method("set_eta_alpha", &hierarchicalPE::set_eta_alpha)
    // .method("set_eta_Omega", &hierarchicalPE::set_eta_Omega)
    // .method("set_eta_gamma", &hierarchicalPE::set_eta_gamma)
       .method("push_data", &hierarchicalPE::push_data)
       .method("get_PLValues", &hierarchicalPE::get_PLValues)
       .method("get_gammaTrace", &hierarchicalPE::get_gammaTrace)
       .method("get_alphaTrace", &hierarchicalPE::get_alphaTrace)
       .method("get_kernelTrace", &hierarchicalPE::get_kernelTrace)
       .method("get_negloglikeTrace", &hierarchicalPE::get_negloglikeTrace)
       .method("get_eta_alpha1", &hierarchicalPE::get_eta_alpha1)
       .method("get_eta_alpha2", &hierarchicalPE::get_eta_alpha2)
       .method("get_eta_Omega1", &hierarchicalPE::get_eta_Omega1)
       .method("get_eta_Omega2", &hierarchicalPE::get_eta_Omega2)
       .method("get_eta_gamma1", &hierarchicalPE::get_eta_gamma1)
       .method("get_eta_gamma2", &hierarchicalPE::get_eta_gamma2)
       .method("get_qf", &hierarchicalPE::get_qf)
       .method("get_fit", &hierarchicalPE::get_fit)
       .method("get_ThetaAll", &hierarchicalPE::get_ThetaAll)

    .method("gpObjective", &hierarchicalPE::gpObjective)
    .method("gpGradient", &hierarchicalPE::gpGradient)
    ;


    Rcpp::class_<lombScargle>("lombScargle")
        .constructor<vector<vec>, vector<vec>, vector<vec>, vec>()
        .method("compute", &lombScargle::compute)
        .method("fitBandI", &lombScargle::fitBandI)
        .method("set_freq", &lombScargle::set_freq)
        .method("gpNegativeLogLik", &lombScargle::gpNegativeLogLik)
        .method("get_covDeriv", &lombScargle::get_covDeriv)
    ;

    Rcpp::class_<hierarchicalLS>("hierarchicalLS")
        .constructor()
        .method("compute", &hierarchicalLS::compute)
        .method("push_data", &hierarchicalLS::push_data)
        .method("setGPBand", &hierarchicalLS::setGPBand)
        .method("set_globalFreq", &hierarchicalLS::set_globalFreq)
        .method("gpNegativeLogLik", &hierarchicalLS::gpNegativeLogLik)
        .method("get_covDeriv", &hierarchicalLS::get_covDeriv)
        .method("set_scaling", &hierarchicalLS::set_scaling)
        .method("getFreqMag", &hierarchicalLS::getFreqMag)
        .method("set_gpUse", &hierarchicalLS::set_gpUse)
        .method("set_individualFreq", &hierarchicalLS::set_individualFreq)
    ;



}























