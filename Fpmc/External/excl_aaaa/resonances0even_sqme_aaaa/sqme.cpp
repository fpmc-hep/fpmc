#include<iostream>
#include<fstream>
#include<math.h>
#include"helicity_amplitudes.h"
using namespace std;

namespace resonances0even_aaaa {
    
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

    *re *= 8 * alpha_em*alpha_em;
    *im *= 8 * alpha_em*alpha_em;

    // the factor of 8 is needed because of the conventions in
    // Costantini, DeTollis, Pistoni
  }

  // Computes the exotic spin 0 even  squared matrix element and the SM interference
  double sqme(double s,double t, int exclude_loops_SM, int exclude_loops_EX, double m, double c, double w_c, double aa){
    double mass = m;
    double f0 = c;
    double w_const = w_c;
    double a2 = aa;
    double mtop = 172.76;
    double f0top = mass/(4*M_PI);
    if (exclude_loops_EX == 2)//even
    {
      w_const = 3/(8*M_PI)*mtop*mtop/(f0top*f0top)*mass*pow((1-4*mtop*mtop/(mass*mass)),3./2.); //width assumption for ttbar paper
    }

    if (exclude_loops_EX == 3)//odd
    {
      w_const = 3/(8*M_PI)*mtop*mtop/(f0top*f0top)*mass*sqrt(1-4*mtop*mtop/(mass*mass)); //width assumption for ttbar paper
    }

    if (s<0 || t >0 || t<-s ){
      cout<<"Invalid domain. Valid range is s>=0 and -s<=t<=0"<<endl;
      return 0;
    }

    //M.S now parameters 
    /*
    // read in input data
    ifstream data;

    data.open ("./data_resonance");

    data >> mass;
    data >> f0;
    data >> width;

    data.close();
    */

    // reduced Mandelstam variables
    //double sred = s/(4 * mass * mass);
    //double tred = t/(4 * mass * mass); 

    double re_ex;
    double im_ex;
    double re_SM;
    double im_SM;

    double value=0;
    
    // Mpppp:

    // the exotic matrix element:
    Mpppp_spin0even(s, t, mass, f0, w_const,a2, &re_ex, &im_ex);
    re_ex *=-1.;
    im_ex *=-1.;
    // the factor of -1 is needed to match the convention of the resonance with the one of 
    // Costantini, DeTollis, Pistoni

    // the SM matrix element:
    me_SM(Mpppp_fermion,s,t,&re_SM,&im_SM, exclude_loops_SM);

    value += re_ex*(re_ex+2*re_SM) + im_ex*(im_ex+2*im_SM) ;
  
    // repeat for the other helicities

    // Mppmm:

    Mppmm_spin0even(s, t, mass, f0, w_const,a2, &re_ex, &im_ex);
    re_ex *=-1.;
    im_ex *=-1.;

    me_SM(Mppmm_fermion,s,t,&re_SM,&im_SM, exclude_loops_SM);

    value += re_ex*(re_ex+2*re_SM) + im_ex*(im_ex+2*im_SM) ;

    // Mpmmp:

    Mpmmp_spin0even(s, t, mass, f0, w_const,a2, &re_ex, &im_ex);
    re_ex *=-1.;
    im_ex *=-1.;

    me_SM(Mpmmp_fermion,s,t,&re_SM,&im_SM, exclude_loops_SM);

    value += re_ex*(re_ex+2*re_SM) + im_ex*(im_ex+2*im_SM) ;

    // Mpmpm:

    Mpmpm_spin0even(s, t, mass, f0, w_const,a2, &re_ex, &im_ex);
    re_ex *=-1.;
    im_ex *=-1.;

    me_SM(Mpmpm_fermion,s,t,&re_SM,&im_SM, exclude_loops_SM);

    value += re_ex*(re_ex+2*re_SM) + im_ex*(im_ex+2*im_SM) ;

    // Mpppm

    Mpppm_spin0even(s, t, mass, f0, w_const,a2, &re_ex, &im_ex);
    re_ex *=-1.;
    im_ex *=-1.;

    me_SM(Mpppm_fermion,s,t,&re_SM,&im_SM, exclude_loops_SM);

    value += 4* (  re_ex*(re_ex+2*re_SM) + im_ex*(im_ex+2*im_SM) );

    return 0.5*value;
  }
}  //namespace resonances0even_aaaa
