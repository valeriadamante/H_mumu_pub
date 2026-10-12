#include <boost/math/special_functions/erf.hpp>
#include <cstdint>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

struct CrystalBall{
    double pi=3.14159;
    double sqrtPiOver2=sqrt(pi/2.0);
    double sqrt2=sqrt(2.0);
    double m;
    double s;
    double a;
    double n;
    double B;
    double C;
    double D;
    double N;
    double NA;
    double Ns;
    double NC;
    double F;
    double G;
    double k;
    double cdfMa;
    double cdfPa;
CrystalBall():m(0),s(1),a(10),n(10){
    init();
}
CrystalBall(double mean, double sigma, double alpha, double n)
    :m(mean),s(sigma),a(alpha),n(n){
    init();
}
void init(){
    double fa = fabs(a);
    double ex = exp(-fa*fa/2);
    double A  = pow(n/fa, n) * ex;
    double C1 = n/fa/(n-1) * ex;
    double D1 = 2 * sqrtPiOver2 * erf(fa/sqrt2);
    B = n/fa-fa;
    C = (D1+2*C1)/C1;
    D = (D1+2*C1)/2;
    N = 1.0/s/(D1+2*C1);
    k = 1.0/(n-1);
    NA = N*A;
    Ns = N*s;
    NC = Ns*C1;
    F = 1-fa*fa/n;
    G = s*n/fa;
    cdfMa = cdf(m-a*s);
    cdfPa = cdf(m+a*s);
}
double pdf(double x) const{
    double d=(x-m)/s;
    if(d<-a) return NA*pow(B-d, -n);
    if(d>a) return NA*pow(B+d, -n);
    return N*exp(-d*d/2);
}
double pdf(double x, double ks, double dm) const{
    double d=(x-m-dm)/(s*ks);
    if(d<-a) return NA/ks*pow(B-d, -n);
    if(d>a) return NA/ks*pow(B+d, -n);
    return N/ks*exp(-d*d/2);

}
double cdf(double x) const{
    double d = (x-m)/s;
    if(d<-a) return NC / pow(F-s*d/G, n-1);
    if(d>a) return NC * (C - pow(F+s*d/G, 1-n) );
    return Ns * (D - sqrtPiOver2 * erf(-d/sqrt2));
}
double invcdf(double u) const{
    if(u<cdfMa) return m + G*(F - pow(NC/u, k));
    if(u>cdfPa) return m - G*(F - pow(C-u/NC, -k) );
    return m - sqrt2 * s * boost::math::erf_inv((D - u/Ns )/sqrtPiOver2);
}
};

double get_rndm(double eta, double phi, float nL, int evtNumber, int lumiNumber, bool use_VXBS) {
    // obtain parameters from correctionlib
    const auto &cs = use_VXBS ? *cset_vxbs : *cset;
    double mean = cs.at("cb_params")->evaluate({abs(eta), nL, 0});
    double sigma = cs.at("cb_params")->evaluate({abs(eta), nL, 1});
    double n = cs.at("cb_params")->evaluate({abs(eta), nL, 2});
    double alpha = cs.at("cb_params")->evaluate({abs(eta), nL, 3});

    // instantiate CB and get random number following the CB
    CrystalBall cb(mean, sigma, alpha, n);
    double rndm = cs.at("RandomSmearing")->evaluate({(int)evtNumber, (int)lumiNumber, phi});

    return cb.invcdf(rndm);
}


double get_std(double pt, double eta, float nL, bool use_VXBS) {

    // obtain paramters from correctionlib
    const auto &cs = use_VXBS ? *cset_vxbs : *cset;
    double param_0 = cs.at("poly_params")->evaluate({abs(eta), nL, 0});
    double param_1 = cs.at("poly_params")->evaluate({abs(eta), nL, 1});
    double param_2 = cs.at("poly_params")->evaluate({abs(eta), nL, 2});

    // calculate value and return max(0, val)
    double sigma = param_0 + param_1 * pt + param_2 * pt*pt;
    if (sigma < 0) sigma = 0;
    return sigma;
}


double get_k(double eta, string var, bool use_VXBS) {

    // obtain parameters from correctionlib
    const auto &cs = use_VXBS ? *cset_vxbs : *cset;
    double k_data = cs.at("k_data")->evaluate({abs(eta), var});
    double k_mc = cs.at("k_mc")->evaluate({abs(eta), var});

    // calculate residual smearing factor
    // return 0 if smearing in MC already larger than in data
    double k = 0;
    if (k_mc < k_data) k = sqrt(k_data*k_data - k_mc*k_mc);
    return k;
}


double pt_resol(double pt, double eta, double phi, float nL, int evtNumber, int lumiNumber, bool use_VXBS, double low_pt_threshold = 20) {
    // load correction values
    double rndm = (double) get_rndm(eta, phi, nL, evtNumber, lumiNumber, use_VXBS);
    double std = (double) get_std(pt, eta, nL, use_VXBS);
    double k = (double) get_k(eta, "nom",use_VXBS);

    // calculate corrected value and return original value if a parameter is nan
    double ptc = pt * ( 1 + k * std * rndm);
    if (isnan(ptc)) ptc = pt;
    if(ptc / pt > 2 || ptc / pt < 0.1 || ptc < 0 || pt < low_pt_threshold || pt > 200){
	ptc = pt;
    }
    // TODO: Understand why for evts with pT < threshold the pt_corr is set to one
    return ptc;
}

double pt_resol_var(double pt_woresol, double pt_wresol, double eta, string updn, bool use_VXBS){

    double k = (double) get_k(eta, "nom", use_VXBS);

    if (k==0) return pt_wresol;
    const auto &cs = use_VXBS ? *cset_vxbs : *cset;
    double k_unc = cs.at("k_mc")->evaluate({abs(eta), "stat"});

    double std_x_rndm = (pt_wresol / pt_woresol - 1) / k;

    double pt_var = pt_wresol;

    if (updn=="up"){
        pt_var = pt_woresol * (1 + (k+k_unc) * std_x_rndm);
    }
    else if (updn=="dn"){
        pt_var = pt_woresol * (1 + (k-k_unc) * std_x_rndm);
    }
    else {
        cout << "ERROR: updn must be 'up' or 'dn'" << endl;
    }
    if(pt_var / pt_woresol > 2 || pt_var / pt_woresol < 0.1 || pt_var < 0){
        pt_var = pt_woresol;
    }

    return pt_var;
}

double pt_scale(bool is_data, double pt, double eta, double phi, int charge,  bool use_VXBS, double low_pt_threshold = 20) {

    // use right correction
    string dtmc = "mc";
    if (is_data) dtmc = "data";
    const auto &cs = use_VXBS ? *cset_vxbs : *cset;
    double a = cs.at("a_"+dtmc)->evaluate({eta, phi, "nom"});
    double m = cs.at("m_"+dtmc)->evaluate({eta, phi, "nom"});
    if(pt < low_pt_threshold)
	    return pt;

    return 1. / (m/pt + charge * a);
}


double pt_scale_var(double pt, double eta, double phi, int charge, string updn, bool use_VXBS) {
    const auto &cs = use_VXBS ? *cset_vxbs : *cset;
    double stat_a = cs.at("a_mc")->evaluate({eta, phi, "stat"});
    double stat_m = cs.at("m_mc")->evaluate({eta, phi, "stat"});
    double stat_rho = cs.at("m_mc")->evaluate({eta, phi, "rho_stat"});

    double unc = pt*pt*sqrt(stat_m*stat_m / (pt*pt) + stat_a*stat_a + 2*charge*stat_rho*stat_m/pt*stat_a);

    double pt_var = pt;

    if (updn=="up"){
        pt_var = pt + unc;
    }
    else if (updn=="dn"){
        pt_var = pt - unc;
    }

    return pt_var;
}
