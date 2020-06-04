#ifndef COMP_FREQ
#define COMP_FREQ

#include "include.hpp"

class freqClass{
public:

    // Return the next frequency value for computation.
    // At the end of iteration, return -1;
    double fetch_next(){
        currentI++;
        if(currentI < n_fSeq)
            return fSeq(currentI);
        else
            return -1;
    }

    // Return the first frequency value
    double fetch_reset(){
        currentI = -1;
        return fetch_next();
    }

    // Assign q(f) at current position
    void set_current_qfValue(double qf){
        qfSeq(currentI) = qf;
        normalized = false;
    }

    // Assign E( p(f, theta| ...) ) at current position
    // void set_current_modelSel(double qf){
    //     modelSel(currentI) = qf;
    // }
    //double getMiraModelEvidence();

    mat get_qfSeq(){ return join_rows(fSeqLong, qfSeqLong); }

    // Sample according to VI q(f)
    double sample_freq(trng::yarn2 &r);
    void set_freq_param(double freqMin_, double freqMax_,
                        size_t nSeqf_);


    // Return the estimated period (posterior mode)
    // and its uncertainty
    vec get_estimatedPeriod();
    // Linear interpolation for more accurate result.
    void normalize_freq();

    double get_freq_atI(int i){
        return fSeq.at(i);
    }

    int get_nfSeq(){
        return n_fSeq;
    }

    void set_qfValue_atI(double qvf_, int i){
        qfSeq(i) = qvf_;
        normalized = false;
    }

private:
    // f
    size_t n_fSeq; // length of freq sequence
    size_t currentI;
    //deltaFreq for integration increment
    double freqMin, freqMax, deltaFreq;
    vec fSeq, qfSeq; // q(f) sequence evaluated at fSeq
    vec fSeqLong, qfSeqLong; // for linear interpolation
    bool normalized;

    // for model selection, not being used for now
    // vec modelSel, modelSelLong;



};



#endif
