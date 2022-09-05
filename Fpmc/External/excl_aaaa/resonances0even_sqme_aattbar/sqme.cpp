#include<iostream>
#include<fstream>
#include<math.h>
#include"helicity_amplitudes.h"
using namespace std;

namespace resonances0even_aattbar {
    
const double alpha_em = 1./137.036; // EM coupling at zero momentum (on shell scheme)



void me_SM(void (*me)(double,double ,double *, double *, int),
    double s,double t, double *re, double*im, int exclude_loops){
  // This routine computes the complex SM amplitude
  // The first argument can be any of the helicity amplitudes Mpppp,Mppmm,Mpmpm,Mpmmp,Mpppm


  // SM fermion content: (e,mu,tau,u,c,t,d,s,b)
  // SM_weight equals (number of colors) * (el. charge)^4  
  // SM masses in GeV
  
  const double SM_weight [9] = {1, 1, 1, 16./27., 16./27., 16./27., 1./27.,  1./27., 1./27.};
  const double SM_masses [9] = {0.5e-3,0.105,1.77,0.0023,1.28,173.07,0.0048,0.095,4.18};

  double d_re;
  double d_im;
  

  *re=0;
  *im=0;

  for (int i=0;i<=8;i++){
    me(s/(4*SM_masses[i]*SM_masses[i]),t/(4*SM_masses[i]*SM_masses[i]), &d_re, &d_im, exclude_loops);
    *re += d_re * SM_weight[i];
    *im += d_im * SM_weight[i];
  }

  //  cout<<*re<<"  "<<*im<<endl;


  // Add also the W contribution

  const double mW=80.385;  // W mass in GeV
 
 
  if (me==Mpppp_fermion){
    Mpppp_vector(s/(4*mW*mW),t/(4*mW*mW), &d_re, &d_im, exclude_loops);
  }
  else if (me==Mppmm_fermion){
    Mppmm_vector(s/(4*mW*mW),t/(4*mW*mW), &d_re, &d_im, exclude_loops);
  }
  else if (me==Mpmpm_fermion){
    Mpmpm_vector(s/(4*mW*mW),t/(4*mW*mW), &d_re, &d_im, exclude_loops);
  }
  else if (me==Mpmmp_fermion){
    Mpmmp_vector(s/(4*mW*mW),t/(4*mW*mW), &d_re, &d_im, exclude_loops);
  }
  else if (me==Mpppm_fermion){
    Mpppm_vector(s/(4*mW*mW),t/(4*mW*mW), &d_re, &d_im, exclude_loops);  
  }
   
  *re += d_re; 
  *im += d_im;
    
  *re *= - 8 * alpha_em*alpha_em;
  *im *= - 8 * alpha_em*alpha_em;
  // the factor of - 8 is needed because of the conventions in
  // Costantini, DeTollis, Pistoni
}

/////////////////////////////// yy->phi->ttbar

// Computes the exotic spin 0 even  squared matrix element and the SM interference
double sqme(double s,double t, int exclude_loops_SM, int exclude_loops_EX, double m, double c, double ch, double w_c, double aa){


  double mass=m;
  double f0=c; //diphoton-phi coupling
  double f0top=c; //phi-ttbar ocupling
//  double w_const=w_c; //width
  double a2=aa; //amplitude square, averaged over spin, summed over polarizations
  double mtop = 172.76;
  //f0top = mass/(4*M_PI);
  f0 = mass/sqrt(4*M_PI);
//  double width_photon = m;
  double width_photon = m*m*m / ( 4*M_PI*f0*f0 ); //a->diphoton width
  double w_const = 3/(8*M_PI)*mtop*mtop/(f0top*f0top)*mass*pow((1-4*mtop*mtop/(mass*mass)),3./2.); //width assumption for ttbar paper
  w_const += width_photon;

  if (s<2*mtop*mtop || t > 0 || t<-s ){
    cout<<"Invalid domain. Valid range is s>=2mtop^2 and -s<=t<=0"<<endl;
    return 0;
  }
  
  // read in input data

  double re_ex;
  double im_ex;

  double value=0;


  return 3.0*mtop*mtop/((f0*f0top)*(f0*f0top))*s*s*(s-4*mtop*mtop)/((s-mass*mass)*(s-mass*mass)+mass*mass*w_const*w_const);
}

}
