#define ARMA_DONT_PRINT_ERRORS
#include <string>
#include <armadillo>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;
using namespace arma;


// These three functions are dealing with file names or read lines
string rm_front_spaces(string in_string) {
    string out_string = in_string;
    while (out_string.substr(0,1) == " ") {
        out_string.erase(0,1);
    }
    return(out_string);
}

string trim_path_from_file_name(string in_string) {
    string out_string = in_string;
    size_t slash_pos = in_string.find_last_of("/\\");
    return(out_string.substr(slash_pos+1));
}

string rm_last_if_slash(string in_string) {
    string out_string = in_string;
    if (in_string.substr(in_string.size()-1) == "/")
        out_string = out_string.substr(0,in_string.size()-1);
    return(in_string);
}

// Define the class: mgpscale
class mgpscale {
    
protected:
    const double PI = 3.1415926536,
        epsilon = 1e-4,
        C1 = 0.001,
        C2 = 0.8;
    double largeV; // A very large value for initialization
    double mean_m1, mean_m2, mean_m3, mean_m4; // prior I-band, JHK bands mean magnitude
    double sig_m1_sq; // prior variance of mean IJHK
    double sig_b1_sq; // prior variance of amplitude & phase
    double period;
    bool flag_inv; // flag of matrix invertible
    mat spc, spc_1, spc_2; // freqs, log-likelihood, opt theta1, opt theta2, curvature
    vec Y, //measurements
    theta; // {theta1 and theta2}
    mat Ht, // Design matrix of the periodic component (transpose)
    sig_gamma, // prior covariance matrix
    HtBH; // Ht * sig_gamma * H
    mat K, // squared exponential kernel
    partial_K, // first order partial derivative of K
    Kc, // squared exponential kernel + observation uncertainties
    Ky, // combined variance matrix
    Ky_inv, // inverse of Ky
    Ky_chol; // Cholesky factorization
    double Ky_logdet; // log (|Ky|), natual log here
    vec alpha; // Ky_inv * Y, an intermediate notation
    // scaling relations
    double Aij_0 = -1.895, Aij_1 = 1.830, Pij_0 = 0.462, Pij_1 = -0.144,
        Aih_0 = -1.662, Aih_1 = 1.675, Pih_0 = 0.493, Pih_1 = -0.155,
        Aik_0 = -1.913, Aik_1 = 1.848, Pik_0 = 0.557, Pik_1 = -0.175;
    
    // methods
    mat compute_Ht(vec);  // get the design matrix
    void compute_kernal(int, int, int, vec &); // compute K or partial derivative of K
    void compute_Kc(vec &, int); // compute Kc
    
public:
    int n_obs, n_band;
    int n_obs_1, n_obs_2, n_obs_3, n_obs_4;
    int n_trials;
    double freq; // trial frequencies
    double tot1=0, tot2=0, tot3=0, tot4=0;
    vec MJD, // date, maybe HJH-2450000 or any reference time in the units of day
    bands, // numeric labels for bands
    mag, // all magnitude
    sig_mag, // uncertainties of measurements
    mag_1,
    sig_1;
    const double least_error = 0.01; // in mag; least measurement uncertainties to avoid unrealistic input
    mat band_info;
    
    // methods
    double compute_mloglik(vec); // compute the negative log-likehood
    vec compute_Dmloglik(vec); // compute the derivative of the above
    double BFGSopt(vec &, mat &); // optimization over theta_1 and theta_2 using the BFGS method
    void set_freq(double); // set the trial frequencies
    void freq_est(vec, int); // compute model for each trial frequency
    // vec check_bad(int);
    vec check_breaks(double);
    void load_input(string); // filename of the light curve table
    void spc_output(string); // write the results to a file
    void res_output(string); ///// added by Zhenfeng, 03-27-2019
};

// define the methods
mat mgpscale::compute_Ht(vec t) {
    int n, ntot=0, bcount=2;
    double A, dP;
    mat X = mat(n_obs, n_band+2, fill::zeros);
    vec phase = 2 * PI * freq * t;
    n = band_info(1,0);
    if (n > 0) {
        X(span(0,n-1),0) = cos(phase(span(0,n-1)));
        X(span(0,n-1),1) = sin(phase(span(0,n-1)));
        X(span(0,n-1),2) = vec(n, fill::ones);
        bcount = bcount + 1;
    }
    ntot = ntot + n;
    n = band_info(1,1);
    if (n > 0) {
        A = Aij_0 + Aij_1 * log10(period);
        A = 1. / A;
        dP = Pij_0 + Pij_1 * log10(period);
        dP = dP * 2 * PI;
        X(span(ntot, ntot+n-1),0) = A*cos(phase(span(ntot, ntot+n-1)) - dP);
        X(span(ntot, ntot+n-1),1) = A*sin(phase(span(ntot, ntot+n-1)) - dP);
        X(span(ntot, ntot+n-1),bcount) = vec(n, fill::ones);
        bcount = bcount + 1;
    }
    ntot = ntot + n;
    n = band_info(1,2);
    if (n > 0) {
        A = Aih_0 + Aih_1 * log10(period);
        A = 1. / A;
        dP = Pih_0 + Pih_1 * log10(period);
        dP = dP * 2 * PI;
        X(span(ntot, ntot+n-1),0) = A*cos(phase(span(ntot, ntot+n-1)) - dP);
        X(span(ntot, ntot+n-1),1) = A*sin(phase(span(ntot, ntot+n-1)) - dP);
        X(span(ntot, ntot+n-1),bcount) = vec(n, fill::ones);
        bcount = bcount + 1;
    }
    ntot = ntot + n;
    n = band_info(1,3);
    if (n > 0) {
        A = Aik_0 + Aik_1 * log10(period);
        A = 1. / A;
        dP = Pik_0 + Pik_1 * log10(period);
        dP = dP * 2 * PI;
        X(span(ntot, ntot+n-1),0) = A*cos(phase(span(ntot, ntot+n-1)) - dP);
        X(span(ntot, ntot+n-1),1) = A*sin(phase(span(ntot, ntot+n-1)) - dP);
        X(span(ntot, ntot+n-1),bcount) = vec(n, fill::ones);
        bcount = bcount + 1;
    }
    return X;
}

void mgpscale::compute_kernal(int p1, int p2, int k, vec & t) {
    int i, j;
    double tmpE, tmpE_deriv, tmpK, tmpDiff, eTheta1, eTheta2;
    K = mat(n_obs, n_obs, fill::zeros);
    partial_K = K;
    eTheta1 = exp(theta(p1));
    eTheta2 = exp(theta(p2));
    for (i = 0; i < n_obs; i++) {
        for (j = i; j < n_obs; j++) {
            tmpDiff = t(i) - t(j);
            tmpE = tmpDiff * tmpDiff / (2 * eTheta2 * eTheta2);
            if (k == 0) {
                tmpK = eTheta1 * eTheta1 * exp(-tmpE);
                if (bands(i) == 1 && bands(j) > 1) {
                    tmpK = 0;
                }
                if (bands(j) == 1 && bands(i) > 1) {
                    tmpK = 0;
                }
                K(i, j) = tmpK;
                K(j, i) = tmpK;
            } else if (k == 1) {
                tmpK = 2 * eTheta1 * exp(-tmpE) * eTheta1;
                if (bands(i) == 1 && bands(j) > 1) {
                    tmpK = 0;
                }
                if (bands(j) == 1 && bands(i) > 1) {
                    tmpK = 0;
                }
                partial_K(i, j) = tmpK;
                partial_K(j, i) = tmpK;
            } else if (k == 2) {
                tmpK = eTheta1 * eTheta1 * exp(-tmpE);
                tmpE_deriv = 2 * tmpE;
                tmpK *= tmpE_deriv;
                if (bands(i) == 1 && bands(j) > 1) {
                    tmpK = 0;
                }
                if (bands(j) == 1 && bands(i) > 1) {
                    tmpK = 0;
                }
                partial_K(i, j) = tmpK;
                partial_K(j, i) = tmpK;
            } else {
                throw " >> Invalid value of k for kernal derivative";
            }
        }
    }
}

void mgpscale::compute_Kc(vec & t, int k) {
    double tmp;
    compute_kernal(0, 1, 0, t);
    Kc = K;
    for (int i = 0; i < n_obs; i++) {
        tmp = max(sig_mag(i), least_error);
        Kc(i, i) += tmp * tmp;
    }
    Ky = Kc + HtBH;
    flag_inv = chol(Ky_chol, Ky);
    if (flag_inv) {
        vec Ky_diag = Ky_chol.diag();
        Ky_logdet = 2 * sum(log(Ky_diag));
        alpha = solve(trimatl(Ky_chol.t()), Y);
        alpha = solve(trimatu(Ky_chol), alpha);
    }
}

// Compute -1 * Q (refer to the paper https://ui.adsabs.harvard.edu/#abs/2016AJ....152..164H)
double mgpscale::compute_mloglik(vec theta_) {
    double mloglik;
    theta = theta_;
    compute_Kc(MJD, 0);
    if (flag_inv) {
        mloglik = 0.5 * (as_scalar(Y.t() * alpha) + Ky_logdet + n_obs * log(2 * PI));
    } else {
        mloglik = largeV;
    }
    return mloglik;
}

vec mgpscale::compute_Dmloglik(vec theta_) {
    vec partialDeri(2, fill::zeros);
    theta = theta_;
    mat tempI = eye(n_obs, n_obs);
    Ky_inv = solve(trimatl(Ky_chol.t()), tempI);
    Ky_inv = solve(trimatu(Ky_chol), Ky_inv);
    mat tempA = alpha * alpha.t() - Ky_inv;
    for (int i = 0; i < 2; i++) {
        compute_kernal(0, 1, i+1, MJD);
        partialDeri(i) = -0.5 * trace(tempA * partial_K);
    }
    return partialDeri;
}

double mgpscale::BFGSopt(vec & theta0, mat & H0) {
    vec deltaFk, deltaFkp1, Xk, Xkp1, pk, sk, yk;
    mat Hk,tmpOut;
    double rho, Fk,Fkp1;
    double armijoRight, curvL,curvR, tmpIn;
    double beta1,stepsize,alow,ahigh;
    int iter,iterMax;
    iter = 0;
    iterMax = 50;
    beta1 = 0.6;
    Hk = H0;
    Xk = theta0;
    stepsize = 1;
    // the covariance matrix, only once
    Fk = compute_mloglik(Xk);
    // if(!flag_inv) return largeV;
    Fkp1 = Fk;
    deltaFk = compute_Dmloglik(Xk);
    while(norm(deltaFk)>epsilon && iter<iterMax){
        iter++;
        // search direction
        pk = -Hk*deltaFk;
        // Wolfe Condition
        stepsize = 1;
        alow = 0;
        ahigh = 10.0;
        tmpIn = dot(deltaFk,pk);
        while(true){
            if(ahigh-alow<1e-3) stepsize = alow;
            Xkp1 = Xk + stepsize*pk;
            Fkp1 = compute_mloglik(Xkp1);
            //check the armijo condition
            armijoRight =  Fk + C1*stepsize*tmpIn;
            if(Fkp1 > armijoRight || !flag_inv){
                ahigh = stepsize;
                beta1 = 0.66;
            }
            else{
                // check the curvature condition
                deltaFkp1 = compute_Dmloglik(Xkp1);
                if(ahigh-alow<0.0001) break;
                curvR = tmpIn;
                curvL = dot(deltaFkp1,pk);
                if(fabs(curvL)<fabs(C2*curvR)){
                    break;
                }else if(curvL > abs(curvR)){
                    ahigh = stepsize;
                    beta1 = 0.66;
                }else{
                    alow = stepsize;
                    beta1 = 0.33;
                }
            }
            if(ahigh-alow<1e-3) break;
            //new try step
            stepsize = alow + (ahigh-alow) * beta1;
        }
        //after fixing alpha
        sk = Xkp1 - Xk;
        if(!flag_inv) return largeV;
        deltaFkp1 = compute_Dmloglik(Xkp1);
        yk = deltaFkp1 - deltaFk;
        // update H
        rho = 1/dot(yk,sk);
        tmpOut = eye<mat>(2,2) - rho*sk*yk.t();
        Hk = tmpOut *Hk * tmpOut.t() +rho*sk*sk.t();
        // k = k + 1
        Xk = Xkp1;
        Fk = Fkp1;
        deltaFk = deltaFkp1;
    }
    theta0 = Xk;
    H0 = Hk;
    // cout << iter << "-----";
    return Fkp1;
}

void mgpscale::set_freq(double f) {
    freq = f;
    period = 1. / freq;
    Ht = compute_Ht(MJD);
    HtBH = Ht * sig_gamma * Ht.t();
}

void mgpscale::freq_est(vec all_freqs, int restart_sep) {
    n_trials = all_freqs.size();
    double y_min, y_opt;
    bool restart;
    int i, j, k,
    n_trial_theta = 66; // Yes, 26 do not work for some lc
    spc = mat(n_trials, 8, fill::zeros);
    spc_1 = spc;
    spc_2 = spc;
    vec theta0 = vec(2, fill::zeros),
        theta0_opt,
        theta0_fix;
    mat H0, // inverse of approx. Hessian
    H0_opt, H0_pre, H0_fix,
    theta_guess = mat(n_trial_theta, 2, fill::randu);
    theta0_fix = theta0;
    H0_fix = eye<mat> (2, 2);
    
    theta_guess = (theta_guess - 0.5) * 20;
    theta_guess(0, 0) = 1;
    theta_guess(0, 1) = 1;
    theta_guess(1, 0) = 8;
    theta_guess(1, 1) = 1;
    theta_guess(2, 0) = 1;
    theta_guess(2, 1) = 8;
    theta_guess(3, 0) = 8;
    theta_guess(3, 1) = 8;
    theta_guess(4, 0) = -0.5;
    theta_guess(4, 1) = 4;
    theta_guess(5, 0) = -2.5;
    theta_guess(5, 1) = 5;
    
    for (k = 0; k < n_trials; k++) {
        // cout << k << " / " << n_trials << endl;
        if ((k % restart_sep) == 0) {
            restart = true;
        } else {
            restart = false;
        }
        set_freq(all_freqs(k));
        if (restart) {
            y_opt = largeV;
            for (i = 0; i < n_trial_theta; i ++) {
                theta0(0) = theta_guess(i, 0);
                theta0(1) = theta_guess(i, 1);
                H0 = eye<mat> (2, 2);
                y_min = BFGSopt(theta0, H0);
                if (y_min < y_opt) {
                    y_opt = y_min;
                    H0_opt = H0;
                    theta0_opt = theta0;
                }
            }
        }
        y_min = BFGSopt(theta0_opt, H0_opt);
        spc_1(k, 0) = freq;
        spc_1(k, 1) = -1 * y_min;
        spc_1(k, 2) = theta0_opt(0);
        spc_1(k, 3) = theta0_opt(1);
        spc_1(k, 4) = H0_opt(0, 0);
        spc_1(k, 5) = H0_opt(0, 1);
        spc_1(k, 6) = H0_opt(1, 0);
        spc_1(k, 7) = H0_opt(1, 1);
    }
    for (k = n_trials - 1; k >= 0; k--) {
        if ((k % restart_sep) == 0) {
            restart = true;
        } else {
            restart = false;
        }
        set_freq(all_freqs(k));
        if (restart) {
            theta0_opt(0) = spc_1(k, 2);
            theta0_opt(1) = spc_1(k, 3);
            H0_opt(0, 0) = spc_1(k, 4);
            H0_opt(0, 1) = spc_1(k, 5);
            H0_opt(1, 0) = spc_1(k, 6);
            H0_opt(1, 1) = spc_1(k, 7);
        } else if (k == n_trials - 1) {
            y_opt = largeV;
            for (i = 0; i < n_trial_theta; i ++) {
                theta0(0) = theta_guess(i, 0);
                theta0(1) = theta_guess(i, 1);
                H0 = eye<mat> (2, 2);
                y_min = BFGSopt(theta0, H0);
                if (y_min < y_opt) {
                    y_opt = y_min;
                    H0_opt = H0;
                    theta0_opt = theta0;
                }
            }
        }
        y_min = BFGSopt(theta0_opt, H0_opt);
        spc_2(k, 0) = freq;
        spc_2(k, 1) = -1 * y_min;
        spc_2(k, 2) = theta0_opt(0);
        spc_2(k, 3) = theta0_opt(1);
        spc_2(k, 4) = H0_opt(0, 0);
        spc_2(k, 5) = H0_opt(0, 1);
        spc_2(k, 6) = H0_opt(1, 0);
        spc_2(k, 7) = H0_opt(1, 1);
    }
    
    for (k = 0; k < n_trials; k++) {
        if (spc_1(k, 1) >= spc_2(k, 1)) {
            spc(k, 0) = spc_1(k, 0);
            spc(k, 1) = spc_1(k, 1);
            spc(k, 2) = spc_1(k, 2);
            spc(k, 3) = spc_1(k, 3);
            spc(k, 4) = spc_1(k, 4);
            spc(k, 5) = spc_1(k, 5);
            spc(k, 6) = spc_1(k, 6);
            spc(k, 7) = spc_1(k, 7);
        } else {
            spc(k, 0) = spc_2(k, 0);
            spc(k, 1) = spc_2(k, 1);
            spc(k, 2) = spc_2(k, 2);
            spc(k, 3) = spc_2(k, 3);
            spc(k, 4) = spc_2(k, 4);
            spc(k, 5) = spc_2(k, 5);
            spc(k, 6) = spc_2(k, 6);
            spc(k, 7) = spc_2(k, 7);
        }
    }
    // Check bad points and fix the initial guess of H0
    int last_idx_size = -1;
    for (j = 0; j < 30; j++) {
        double diff_val = 0.2;
        int inc_length = 1;
        vec idx = check_breaks(diff_val);
        vec idx_2 = check_breaks(diff_val*10);
        vec idx_refit = vec(idx.size() * (2*inc_length + 1));
        int jk = 0, ji, jj, max_idx = 0;
        
        if (last_idx_size <= int (idx.size())) {
            break;
        }
        last_idx_size = idx.size();
        
        // cout << idx_2.size() << endl;
        if (idx_2.size() > 10)
            break;
        
        if (idx.size() > 5 && max(spc(span(0, n_trials-1),1)) > 0)
            break;
        
        if (idx.size() > 0) {
            for (int ji = 0; ji < int (idx.size()); ji++) {
                for (jj = -inc_length; jj <= inc_length; jj++) {
                    idx_refit(jk) = idx(ji) + jj;
                    jk++;
                }
            }
            idx_refit = unique(idx_refit);
            int n_refit = idx_refit.size();
            cout << "  Refit for the " << j+1 << "th iteration:";
            cout << "   Fixing at " << n_refit << " positions..." << endl;
            
            H0_fix(0, 0) = median(spc(span(0, n_trials-1),4));
            H0_fix(0, 1) = median(spc(span(0, n_trials-1),5));
            H0_fix(1, 0) = median(spc(span(0, n_trials-1),6));
            H0_fix(1, 1) = median(spc(span(0, n_trials-1),7));
            theta_guess(0, 0) = median(spc(span(0, n_trials-1),2));
            theta_guess(0, 1) = median(spc(span(0, n_trials-1),3));
            theta_guess(1, 0) = (randu()-0.5) * j;
            theta_guess(1, 1) = (randu()-0.5) * j;
            
            vec tmp_vec = spc(span::all,1);
            double max_val = max(tmp_vec);
            for (int i_max = 0; i_max < n_trials; i_max++) {
                if (tmp_vec(i_max) == max_val)
                    max_idx = i_max;
            }
            theta_guess(2, 0) = spc(max_idx, 2);
            theta_guess(2, 1) = spc(max_idx, 3);
            
            
            for (jk = 0; jk < n_refit; jk++) {
                if (idx_refit(jk) >= 0 && idx_refit(jk) < n_trials) {
                    set_freq(all_freqs(idx_refit(jk)));
                    if (jk >= 0) {
                        y_opt = largeV;
                        y_min = largeV;
                        for (ji = 0; ji < 3; ji++) {
                            theta0(0) = theta_guess(ji, 0);
                            theta0(1) = theta_guess(ji, 1);
                            H0 = H0_fix;
                            y_min = BFGSopt(theta0, H0);
                            if (y_min < y_opt) {
                                y_opt = y_min;
                                H0_opt = H0; // Find a good Hessian Matrix
                                theta0_opt = theta0;
                            }
                        }
                        if (idx_refit(jk) > 1) {
                            theta0(0) = spc(idx_refit(jk)-1,2);
                            theta0(1) = spc(idx_refit(jk)-1,3);
                        }
                        y_min = BFGSopt(theta0, H0_opt);
                        if (y_min < y_opt) {
                            y_opt = y_min;
                            H0_opt = H0;
                            theta0_opt = theta0;
                        }
                    }
                    y_min = BFGSopt(theta0_opt, H0_opt);
                    spc(idx_refit(jk), 0) = freq;
                    spc(idx_refit(jk), 1) = -1 * y_min;
                    spc(idx_refit(jk), 2) = theta0_opt(0);
                    spc(idx_refit(jk), 3) = theta0_opt(1);
                    spc(idx_refit(jk), 4) = H0_opt(0, 0);
                    spc(idx_refit(jk), 5) = H0_opt(0, 1);
                    spc(idx_refit(jk), 6) = H0_opt(1, 0);
                    spc(idx_refit(jk), 7) = H0_opt(1, 1);
                }
            }
        } else {
            break;
        }
    }
}

vec mgpscale::check_breaks(double diff_val) {
    vec tmp_idx(n_trials), tmp_diff;
    int k, j = 0;
    tmp_diff = spc(span(2, n_trials-1), 1) - 2 * spc(span(1, n_trials-2), 1) + spc(span(0, n_trials-3), 1);
    tmp_diff = abs(tmp_diff);
    for (k = 0; k < n_trials - 2; k++) {
        if (tmp_diff(k) > diff_val) {
            tmp_idx(j) = k;
            j++;
        }
    }
    j--;
    if (j > 0) {
        tmp_idx = tmp_idx.head(j);
        return tmp_idx;
    } else {
        vec no_idx;
        return no_idx;
    }
}

void mgpscale::load_input(string f_lc) {
    string line;
    
    n_obs = 0;
    ifstream infile1(f_lc);
    while (getline(infile1, line))
        ++n_obs;
    infile1.close();
    
    MJD = vec(n_obs);
    bands = vec(n_obs);
    mag = vec(n_obs);
    sig_mag = vec(n_obs);
    
    std::vector<string> bands_tmp;
    string band_tmp;
    
    int iline;
    ifstream lcfile(f_lc);
    if (lcfile.is_open()) {
        for (iline=0; iline < n_obs; iline++) {
            getline(lcfile, line);
            istringstream in(line);
            in >> MJD(iline) >> band_tmp >> mag(iline) >> sig_mag(iline);
            bands_tmp.push_back(band_tmp);
        }
        lcfile.close();
        
        for (size_t i = 0; i < n_obs; i++) {
            if (bands_tmp.at(i).at(1) == 'I') {
                bands(i) = 1.;
            }
            if (bands_tmp.at(i).at(1) == 'J') {
                bands(i) = 2.;
            }
            if (bands_tmp.at(i).at(1) == 'H') {
                bands(i) = 3.;
            }
            if (bands_tmp.at(i).at(1) == 'K') {
                bands(i) = 4.;
            }
        }
        
        //for (size_t i = 0; i < n_obs; i++) {
        //  cout << bands.at(i) << endl;
        //}
        
        // already sorted by band, compute the number of measurements and mean mag in each band
        for (int i=0; i<n_obs; i++) {
            if (bands(i) == 1.) tot1 = tot1 + mag(i);
            if (bands(i) == 2.) tot2 = tot2 + mag(i);
            if (bands(i) == 3.) tot3 = tot3 + mag(i);
            if (bands(i) == 4.) tot4 = tot4 + mag(i);
        }
        n_obs_1 = sum(bands == 1.);
        n_obs_2 = sum(bands == 2.);
        n_obs_3 = sum(bands == 3.);
        n_obs_4 = sum(bands == 4.);
        mean_m1 = tot1 / n_obs_1;
        mean_m2 = tot2 / n_obs_2;
        mean_m3 = tot3 / n_obs_3;
        mean_m4 = tot4 / n_obs_4;
        vec meanmags = vec(n_obs, fill::zeros);
        for (int i=0; i < n_obs_1; i++) meanmags(i) = mean_m1;
        for (int i=0; i < n_obs_2; i++) meanmags(i+n_obs_1) = mean_m2;
        for (int i=0; i < n_obs_3; i++) meanmags(i+n_obs_1+n_obs_2) = mean_m3;
        for (int i=0; i < n_obs_4; i++) meanmags(i+n_obs_1+n_obs_2+n_obs_3) = mean_m4;
        Y = mag - meanmags;
        band_info = mat(2, 4, fill::zeros);
        band_info(0, 0) = 1;
        band_info(0, 1) = 2;
        band_info(0, 2) = 3;
        band_info(0, 3) = 4;
        band_info(1,0) = n_obs_1;
        band_info(1,1) = n_obs_2;
        band_info(1,2) = n_obs_3;
        band_info(1,3) = n_obs_4;
        n_band = sum(band_info(1, span(0,3)) > 0);
        sig_gamma = mat(2+n_band, 2+n_band, fill::zeros);
        for (int i=0; i<2+n_band; i++) {
            sig_gamma(i, i) = 5;
        }
        largeV = 1.0e30;
    } else {
        cout << "Cannot open file " << f_lc << "\n";
        abort();
    }
}

void mgpscale::spc_output(string f_output) {
    ofstream outfile (f_output);
    if (outfile.is_open()){
        for (int i = 0; i < n_trials; i++)
        {
            if (spc(i,0) != 0) {
                outfile << fixed << setprecision(6) << spc(i,0)
                        << "  " << spc(i,1)
                        << "  " << spc(i,2)
                        << "  " << spc(i,3)
                // << "  " << spc(i,4)
                // << "  " << spc(i,5)
                // << "  " << spc(i,6)
                // << "  " << spc(i,7)
                   <<endl;
            }
        }
        outfile.close();
    }
}

void mgpscale::res_output(string f_output) {
    ofstream outfile (f_output);
    
    uvec bad_ids = find(spc(span(),0) >= 0.01);
    /*
     for (uword row_id = 0; row_id < spc.n_rows; row_id++) {
     if (spc(row_id,0) >= 0.01) {
     bad_ids.push_back(row_id);
     }
     }
     spc.shed_row(bad_ids[0],bad_ids[bad_ids.n_slices-1]);
     */
    if (outfile.is_open()){
        uword i = spc(span(0, bad_ids[0]), 1).index_max();
        //cout << i << '\n';
        
        //for (int i = 0; i < bad_ids[0]; i++) {
        if (spc(i,0) != 0) {
            outfile << fixed << setprecision(6) << spc(i,0)
                    << "  " << spc(i,1)
                    << "  " << spc(i,2)
                    << "  " << spc(i,3)
            // << "  " << spc(i,4)
            // << "  " << spc(i,5)
            // << "  " << spc(i,6)
            // << "  " << spc(i,7)
               <<endl;
        }
        //}
        outfile.close();
    }
}

int main(int argc, char* argv[ ]) {
    //string output_dir = "./gp_spectra/";
    string output_extension = ".dat";
    int restart_sep = 50;
    vec all_freqs, freqs_hp, freqs_lp;
    double p1 = 2000.;
    double p2 = 100.;
    double p3 = 10.;
    double hp_resolution = 1.e-5; // in frequency
    double lp_resolution = 0.2; // in days
    int n_freqs_hp, n_freqs_lp;
    n_freqs_hp = int((1./p2 - 1./p1) / hp_resolution);
    n_freqs_lp = int((p2 - p3) / lp_resolution);
    all_freqs = vec(n_freqs_hp + n_freqs_lp);
    for (int i = 0; i < n_freqs_hp; i++) {
        all_freqs(i) = 1./p2 - double(i) * hp_resolution;
    }
    for (int i = 0; i < n_freqs_lp; i++) {
        all_freqs(i + n_freqs_hp) = 1./(p2 - double(i) * lp_resolution);
    }
    all_freqs = unique(all_freqs);
    
    string f_lc;
    string f_output, f_lc_trim, shell_command;
    //output_dir = rm_last_if_slash(output_dir);
    //shell_command = "mkdir -p " + output_dir;
    //const char *ts_command = shell_command.c_str();
    //system(ts_command);
    f_lc = argv[1];
    f_output = argv[2] + output_extension;
    mgpscale lc;
    lc.load_input(f_lc);
    lc.freq_est(all_freqs, restart_sep);
    //f_lc_trim = trim_path_from_file_name(f_lc);
    //f_output = output_dir + "/" + f_lc_trim + output_extension;
    if (argc == 4) {
        lc.res_output(f_output);
    } else {
        lc.spc_output(f_output);
    }
    return 0;
}