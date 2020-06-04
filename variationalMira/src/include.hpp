#ifndef BASIC_INCLUDE
#define BASIC_INCLUDE

// #define USE_RCPP_ARMADILLO
// #include "optim.hpp"

#include <vector>
#include <RcppArmadillo.h>
#include <trng/uniform_int_dist.hpp>
#include <trng/discrete_dist.hpp>
#include <trng/mt19937_64.hpp>
#include <trng/yarn2.hpp>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <memory>
#include <stdexcept>
#include <fstream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace arma;
using std::vector;
using namespace std;




#define NUM_PARA_THREAD 4
typedef shared_ptr<mat> mat_ptr;



inline bool notEnoughBandSample(int sampleSize){
    return (sampleSize < 3);
}

#endif
