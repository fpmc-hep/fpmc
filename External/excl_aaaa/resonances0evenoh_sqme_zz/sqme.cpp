#include<iostream>
#include<fstream>
#include<math.h>
#include"helicity_amplitudes.h"
using namespace std;

namespace resonances0evenoh_zz {
  const double alpha_em = 1./137.036; // EM coupling at zero momentum (on shell scheme)

  void
  me_SM( void (*me)(double,double ,double *, double *, int), double s,double t, double *re, double*im, int exclude_loops )
  {
    // This routine computes the complex SM amplitude
    // The first argument can be any of the helicity amplitudes Mpppp,Mppmm,Mpmpm,Mpmmp,Mpppm

    // SM fermion content: (e,mu,tau,u,c,t,d,s,b)
    // SM_weight equals (number of colors) * (el. charge)^4
    // SM masses in GeV
  
    const double SM_weight [9] = {1, 1, 1, 16./27., 16./27., 16./27., 1./27.,  1./27., 1./27.};
    const double SM_masses [9] = {0.5e-3,0.105,1.77,0.0023,1.28,173.07,0.0048,0.095,4.18};

    double d_re;
    double d_im;

    *re = *im = 0.;

    for ( int i=0; i<=8; i++ ) {
      me(s/(4*SM_masses[i]*SM_masses[i]),t/(4*SM_masses[i]*SM_masses[i]), &d_re, &d_im, exclude_loops);
      *re += d_re * SM_weight[i];
      *im += d_im * SM_weight[i];
    }

    // Add also the W contribution

    const double mW=80.385;  // W mass in GeV

    if      ( me==Mpppp_fermion ) Mpppp_vector(s/(4*mW*mW),t/(4*mW*mW), &d_re, &d_im, exclude_loops);
    else if ( me==Mppmm_fermion ) Mppmm_vector(s/(4*mW*mW),t/(4*mW*mW), &d_re, &d_im, exclude_loops);
    else if ( me==Mpmpm_fermion ) Mpmpm_vector(s/(4*mW*mW),t/(4*mW*mW), &d_re, &d_im, exclude_loops);
    else if ( me==Mpmmp_fermion ) Mpmmp_vector(s/(4*mW*mW),t/(4*mW*mW), &d_re, &d_im, exclude_loops);
    else if ( me==Mpppm_fermion ) Mpppm_vector(s/(4*mW*mW),t/(4*mW*mW), &d_re, &d_im, exclude_loops);

    *re += d_re;
    *im += d_im;
    
    *re *= - 8 * alpha_em*alpha_em;
    *im *= - 8 * alpha_em*alpha_em;
    // the factor of - 8 is needed because of the conventions in
    // Costantini, DeTollis, Pistoni
  }

  /////////////////////////////// Z Z final state from OZ operator

  // Computes the exotic spin 0 even  squared matrix element and the SM interference
  double
  sqme( double s,double t, int exclude_loops_SM, int exclude_loops_EX, double m, double c, double ch, double w_c, double aa )
  {
    double mass = m;
    double f0 = c;
    double f0H = ch;
    double w_const = w_c;
    double a2 = aa;

    if ( s<0 || t >0 || t<-s ) {
      cout<<"Invalid domain. Valid range is s>=0 and -s<=t<=0"<<endl;
      return 0;
    }
    double re_ex;
    double im_ex;
    double value=0;

    // Mpppp:

    // the exotic matrix element:
    MZZOHpppp_spin0even(s, t, mass, f0, f0H, w_const,a2, &re_ex, &im_ex);

    value += re_ex*(re_ex) + im_ex*(im_ex);

    // repeat for the other helicities

    // Mppmm:
    MZZOHppmm_spin0even(s, t, mass, f0, f0H, w_const,a2, &re_ex, &im_ex);
    value += re_ex*re_ex + im_ex*im_ex;

    // Mpmmp:
    MZZOHpmmp_spin0even(s, t, mass, f0, f0H, w_const,a2, &re_ex, &im_ex);
    value += re_ex*re_ex + im_ex*im_ex;

    // Mpmpm:
    MZZOHpmpm_spin0even(s, t, mass, f0, f0H, w_const,a2, &re_ex, &im_ex);
    value += re_ex*re_ex + im_ex*im_ex;

    // Mpppm
    MZZOHpppm_spin0even(s, t, mass, f0, f0H, w_const,a2, &re_ex, &im_ex);
    value += re_ex*re_ex + im_ex*im_ex;
 
    // Mpp00
    MZZOHpp00_spin0even(s, t, mass, f0, f0H, w_const,a2, &re_ex, &im_ex);
    value += re_ex*re_ex + im_ex*im_ex;

    return 0.5*value;
  }
}
