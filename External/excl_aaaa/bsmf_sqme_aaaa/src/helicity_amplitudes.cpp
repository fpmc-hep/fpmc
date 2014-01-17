// Computes different helicity amplitudes as defined in 
// Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787 

#include<math.h>
#include"./functions.h"



void Mxxxx_fermion(double x, double y, double * re, double * im){
  // some auxilliary function used in Mpppp, Mpmpm, Mpmmp.


  *re=1;
  *im=0;

  double z = - x - y;
  
  double temp;
 
  temp = 2 * ( y*y + z*z ) / ( x * x ) - 2/x;
  *re += temp * ( ReT(y) + ReT(z) );
  *im += temp * ( ImT(y) + ImT(z) );

  temp =  1/(2 * x * y) - 1/y ; 
  *re += temp * ReI(x,y); 
  *im += temp * ImI(x,y);

  temp =  1/(2 * x * z) - 1/z ; 
  *re += temp * ReI(x,z); 
  *im += temp * ImI(x,z);

  temp =  4/x +1/y +1/z + 1/(2*z*y)  -  2 * ( y*y + z*z ) / ( x * x );
  *re += temp * ReI(y,z);
  *im += temp * ImI(y,z);

  temp = 2* (y-z)/x;
  *re += temp * ( ReB(y) - ReB(z) );
  *im += temp * ( ImB(y) - ImB(z) );

  return;

};

void Mpppp_fermion(double sred, double tred, double *re, double *im){
  // M++++ from Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787 

  Mxxxx_fermion(sred,tred,re,im);

  return;

};

void Mpmmp_fermion(double sred, double tred, double *re, double *im){
  // M+--+ from Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787 

  Mxxxx_fermion(tred,sred,re,im);

  return;

};

void Mpmpm_fermion(double sred, double tred, double *re, double *im){
// M+-+- from Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787 

  Mxxxx_fermion(-tred-sred,tred,re,im);

  return;

};

void Mpppm_fermion(double sred, double tred, double * re, double * im){
  // M+--- from Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787 

  double temp;

  double ured=-sred-tred;

  *re=-1;
  *im=0;

  temp=-1/sred-1/tred-1/ured;
  *re += temp*( ReT(sred) + ReT(tred) + ReT(ured) );
  *im += temp*( ImT(sred) + ImT(tred) + ImT(ured) );
  
  temp = 1/ured + 1/( 2 * sred * tred );
  *re += temp*ReI(sred,tred);
  *im += temp*ImI(sred,tred);

  temp = 1/tred + 1/( 2 * sred * ured );
  *re += temp*ReI(sred,ured);
  *im += temp*ImI(sred,ured);
  
  temp = 1/sred + 1/( 2 * tred * ured );
  *re += temp*ReI(tred,ured);
  *im += temp*ImI(tred,ured);

  return ;

};

void Mppmm_fermion(double sred, double tred, double * re, double * im){
  // M++-- from Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787 

  double temp;
  double ured = -sred-tred;


  *re=-1;
  *im=0;
    
  temp = 1/( 2 * sred * tred );
  *re += temp*ReI(sred,tred);
  *im += temp*ImI(sred,tred);

  temp = 1/( 2 * sred * ured );
  *re += temp*ReI(sred,ured);
  *im += temp*ImI(sred,ured);
  
  temp = 1/( 2 * tred * ured );
  *re += temp*ReI(tred,ured);
  *im += temp*ImI(tred,ured);

  return;

};



void Mxxxx_vector(double x, double y, double * re, double * im){

  // some auxilliary function used in Mpppp, Mpmpm, Mpmmp.


  *re=-1.5;
  *im=0;

  double z = - x - y;
  
  double temp;


  temp = -3* (y-z)/x;
  *re += temp * ( ReB(y) - ReB(z) );
  *im += temp * ( ImB(y) - ImB(z) );

 
  temp = -1/x*(8*x-3-6*y*z/x);
  *re += temp * ( ReT(y) + ReT(z) );
  *im += temp * ( ImT(y) + ImT(z) );

  temp =  1/x*(8*x-6-6*y*z/x)-4*(x-0.25)*(x-0.75)/(y*z);
  *re += temp * ReI(y,z); 
  *im += temp * ImI(y,z);


  temp = -4*(x-0.25)*(x-0.75)/(x*y); ; 
  *re += temp * ReI(x,y); 
  *im += temp * ImI(x,y);

  temp = -4*(x-0.25)*(x-0.75)/(x*z); ; 
  *re += temp * ReI(x,z); 
  *im += temp * ImI(x,z);


  return;


  
};

void Mpppp_vector(double sred, double tred, double *re, double *im){

  Mxxxx_vector(sred,tred,re,im);

  return;

};

void Mpmmp_vector(double sred, double tred, double *re, double *im){

  Mxxxx_vector(tred,sred,re,im);

  return;

};

void Mpmpm_vector(double sred, double tred, double *re, double *im){

  Mxxxx_vector(-tred-sred,tred,re,im);

  return;

};

void Mpppm_vector(double sred, double tred, double * re, double * im){
  Mpppm_fermion(sred,tred,re,im);
  *re *= -1.5;
  *im *= -1.5;

  return;
};

void Mppmm_vector(double sred, double tred, double * re, double * im){

  Mppmm_fermion(sred,tred,re,im);
  *re *= -1.5;
  *im *= -1.5;

  return;

};
