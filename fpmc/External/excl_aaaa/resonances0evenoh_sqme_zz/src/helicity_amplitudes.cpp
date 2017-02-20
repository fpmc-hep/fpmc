// Computes different helicity amplitudes as defined in 
// Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787 

#include<math.h>
#include"./functions.h"
#include<iostream>
#include<fstream>


const double PI = 4*atan(1);
const double mW = 80.4 ;
const double mZ = 91.2 ;
const double mh = 125. ;
const int low=1;
const int high=2;
const int forward=3;
const int backward=4;




int limits(double sred,double tred, double ured){

  double const shigh=pow(10.,9.);

  if ( sred <= 0.001 ) return low; // EFT limit
  else if ( ( sred <= 10. && -tred < 0.0001*sred ) ||
      ( sred >  10. && sred <= shigh && -tred < 0.001)  ||
      ( sred > shigh && -tred < 1.) ) return forward; // forward limit
  else if ( ( sred <= 10.  && -ured < 0.0001*sred ) ||
      ( sred >  10.  && -ured < 0.001 ) ||
      ( sred > shigh && -ured < 1.)) return backward;  // backward limit
  else if ( sred > shigh ) return high; // high energy limit
  else return 0; // no limit

  // explanation: 
  // for sred>shigh, optimal value to switch from HE limit to forward limit is |tred|=1
  // (at these value both limits are somewhat bad but quickly converge at either side)
  // for sred<shigh, only switch from exact result to forward limit at |t| < 0.001 for better accuracy  
};


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

void Mpppp_fermion(double sred, double tred, double *re, double *im, int exclude_loops){
  // M++++ from Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787 

  double ured=-sred-tred;

  if(exclude_loops==1||exclude_loops==3) {*re=0; *im=0;}
  else{ 
    int region=limits(sred,tred,ured);
    if( region==low ){ // EFT limit
      *re= -4.*(4.*(-1./36.)  +3.*(7./90.) )  *sred*sred ;  
      *im=0;  }
    else if( region == forward || region == backward )
      {            // Forward and backward limit 
  *re= 1./(2.* sred*sred)*( 2.* sred*sred+(-2.*sred+4.*sred*sred)*ReB(sred)
          +(2.*sred-8.*sred*sred)*ReB(-sred) +(-1.+2.*sred)*ReT(sred)+(-1.-2.*sred+4.*sred*sred)*ReT(-sred) )  ;
  *im= 1./(2.* sred*sred)*(               (-2.*sred+4.*sred*sred)*ImB(sred)
            +(2.*sred-8.*sred*sred)*ImB(-sred) +(-1.+2.*sred)*ImT(sred)+(-1.-2.*sred+4.*sred*sred)*ImT(-sred) )  ;
      }
    else if( region== high ) 
      {            // high energy limit
  *re = 1.+(tred-ured)/sred * log(tred/ured)+(tred*tred+ured*ured) / (2.*sred*sred)*(pow( log(tred/ured) , 2 )+PI*PI);
  *im = 0;
      }
    else{   
      Mxxxx_fermion(sred,tred,re,im);
    };
    
  } 
  
  return;

};

void Mpmmp_fermion(double sred, double tred, double *re, double *im, int exclude_loops){
  // M+--+ from Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787 

  double ured=-sred-tred;
  
  if(exclude_loops==1||exclude_loops==3) {*re=0; *im=0;}
  else
    {
      int region = limits(sred,tred,ured); 
      if( region == low ){ // EFT limit
  *re= -4.*(4.*(-1./36.)  +3.*(7./90.) )  *tred*tred ;
  *im=0;  
      }
      else if( region==forward )
  {                // Forward limit 
    *re=0.; *im=0.;
  }
      else if( region == backward )
  {                // Backward limit 
    *re= 1./(2.* sred*sred)*( 2.* sred*sred+(2.*sred+4.*sred*sred)*ReB(-sred)+(-2.*sred-8.*sred*sred)*ReB(sred) +(-1.-2.*sred)*ReT(-sred)+(-1.+2.*sred+4.*sred*sred)*ReT(sred) )  ;
    *im= 1./(2.* sred*sred)*(               (2.*sred+4.*sred*sred)*ImB(-sred)+(-2.*sred-8.*sred*sred)*ImB(sred) +(-1.-2.*sred)*ImT(-sred)+(-1.+2.*sred+4.*sred*sred)*ImT(sred) )  ;
  }
      else if( region == high ) 
  {            // high energy limit
    *re = 1. + (sred-ured)/tred * log(-sred/ured)+(sred*sred+ured*ured) / (2.*tred*tred)*pow( log(-sred/ured) , 2 );
    *im = -PI*( (sred-ured)/tred + (sred*sred+ured*ured) / (tred*tred)*log(-sred/ured));
  }
      else{  Mxxxx_fermion(tred,sred,re,im); }     
    }
  return;
 
};

void Mpmpm_fermion(double sred, double tred, double *re, double *im, int exclude_loops){
// M+-+- from Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787 

  double ured=-sred-tred;

  if(exclude_loops==1||exclude_loops==3) {*re=0; *im=0;}
  else{
    int region = limits(sred,tred,ured);     
    if( region == low ){ // EFT limit
      *re= -4.*(4.*(-1./36.)  +3.*(7./90.) )  *ured*ured ;
      *im=0;  }
    else if( region == forward )
      {                // Forward limit 
  *re= 1./(2.* sred*sred)*( 2.* sred*sred+(2.*sred+4.*sred*sred)*ReB(-sred)+(-2.*sred-8.*sred*sred)*ReB(sred) +(-1.-2.*sred)*ReT(-sred)+(-1.+2.*sred+4.*sred*sred)*ReT(sred) )  ;
  *im= 1./(2.* sred*sred)*(               (2.*sred+4.*sred*sred)*ImB(-sred)+(-2.*sred-8.*sred*sred)*ImB(sred) +(-1.-2.*sred)*ImT(-sred)+(-1.+2.*sred+4.*sred*sred)*ImT(sred) )  ;    
      }
    else if( region == backward )
      {                // Backward limit 
  *re=0.; *im=0.;
      }
    else if( region == high ) 
      {            // high energy limit
  *re = 1. + (tred-sred)/ured * log(-tred/sred)+(sred*sred+tred*tred) / (2.*ured*ured)*pow( log(-tred/sred) , 2 );
  *im = PI*( (tred-sred)/ured + (sred*sred+tred*tred) / (ured*ured)*log(-tred/sred));
      }   
    else{   
      Mxxxx_fermion(ured,tred,re,im);     
    }
    
  }
  return;
  
};

void Mpppm_fermion(double sred, double tred, double * re, double * im, int exclude_loops){
  // M+--- from Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787 

  double temp;
  double ured=-sred-tred;


  if(exclude_loops==1||exclude_loops==3) {*re=0; *im=0;}
  else{
    int region = limits (sred, tred, ured);    
    if( region == low ){ // EFT limit
      *re= 0.; *im=0.;}
    else if( region == forward || region == backward )
      {                // Forward and backward limit 
  *re= 0.; *im=0.;}
      else if (region == high ) 
      {            // high energy limit
  *re=-1;
  *im=0;
      }
    else{
      
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
    }
    
  }


  return ;

};

void Mppmm_fermion(double sred, double tred, double * re, double * im, int exclude_loops){
  // M++-- from Costantini, DeTollis, Pistoni; Nuovo Cim. A2 (1971) 733-787 

  double temp;
  double ured = -sred-tred;

  if(exclude_loops==1||exclude_loops==3) {*re=0; *im=0;}
  else{ 
    int region = limits (sred, tred, ured);
    if( region == low ){ // EFT limit
      *re= -4.*(4.*(-1./36.)  +(7./90.) )*(sred*sred+tred*tred+ured*ured); *im=0.;}
    else if( region == forward || region == backward )
      {                // Forward and backward limit 
  *re=1./(2.*sred*sred)*( -2.*sred*sred-2.*sred*ReB(sred)+2.*sred*ReB(-sred)-ReT(sred)-ReT(-sred)  ) ; 
  *im=1./(2.*sred*sred)*(              -2.*sred*ImB(sred)+2.*sred*ImB(-sred)-ImT(sred)-ImT(-sred)  );
      }
    else if ( region == high ) 
      {            // high energy limit
  *re=-1;
  *im=0;
      }
    else{ 
      
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
      
    }
    
  }


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

void Mpppp_vector(double sred, double tred, double *re, double *im, int exclude_loops){

  double ured=-sred-tred;

  if(exclude_loops==2||exclude_loops==3) {*re=0; *im=0;}
  else{ 
    int region = limits (sred, tred, ured);
    if( region == low ){ // EFT limit
      *re= -4.*(4.*(-5./32.)  +3.*(27./40.) )  *sred*sred ;
      *im=0;  }
    else if ( region == forward || region == backward )
      {           // Forward and backward limit 
  *re=-3./2.+8.*(sred-0.25)*(sred-0.75)/sred*ReB(sred)+
    (-8.*(sred-0.25)*(sred-0.75)/sred+3.)*ReB(-sred)+
    4.*(sred-0.25)*(sred-0.75)/(sred*sred)*ReT(sred)+
    (4.*(sred-0.25)*(sred-0.75)/(sred*sred)-(8.*sred-3.)/sred)*ReT(-sred);
  *im=8.*(sred-0.25)*(sred-0.75)/sred*ImB(sred)+
    (-8.*(sred-0.25)*(sred-0.75)/sred+3.)*ImB(-sred)+
    4.*(sred-0.25)*(sred-0.75)/(sred*sred)*ImT(sred)+
    (4.*(sred-0.25)*(sred-0.75)/(sred*sred)-(8.*sred-3.)/sred)*ImT(-sred);
      }
    else if (region == high)
      {          // high energy limit
  *re = - 1.*( 1.5 +  
      1.5* (ured-tred)/sred * log(ured/tred) +
      2. * ( 1. - 0.75 * tred*ured / (sred*sred) ) * ( pow( log(ured/tred),2) + PI*PI ) +
      2. * sred * sred * ( log(4.*sred)*log(-4.*tred)/(sred*tred)+
               log(4.*sred)*log(-4.*ured)/(sred*ured)+
               log(-4.*ured)*log(-4.*tred)/(ured*tred) )
      );
  *im =  ( 2. * PI * sred*sred *( log(-4.*ured)  / (sred * ured) +
          log(-4.*tred)  / (sred * tred ) ) 
      )  ;
      }
    else{   Mxxxx_vector(sred,tred,re,im);  
    }
  }
  return;

};

void Mpmmp_vector(double sred, double tred, double *re, double *im, int exclude_loops){

  double ured=-sred-tred;

  if(exclude_loops==2||exclude_loops==3) {*re=0; *im=0;}
  else{
    int region = limits (sred, tred, ured); 
    if( region == low ){ // EFT limit
      *re= -4.*(4.*(-5./32.)  +3.*(27./40.) )  *tred*tred ;
      *im=0;  }
    else if( region == forward )
      {                // Forward limit 
        *re=0.; *im=0.;
      }
    else if( region == backward )
      {                // Backward limit 
  *re=-3./2.-8.*(-sred-0.25)*(-sred-0.75)/sred*ReB(-sred)+
    (8.*(-sred-0.25)*(-sred-0.75)/sred+3.)*ReB(sred)+
    4.*(-sred-0.25)*(-sred-0.75)/(sred*sred)*ReT(-sred)+
    (4.*(-sred-0.25)*(-sred-0.75)/(sred*sred)+(-8.*sred-3.)/sred)*ReT(sred);
  *im=-8.*(-sred-0.25)*(-sred-0.75)/sred*ImB(-sred)+
    (8.*(-sred-0.25)*(-sred-0.75)/sred+3.)*ImB(sred)+
    4.*(-sred-0.25)*(-sred-0.75)/(sred*sred)*ImT(-sred)+
    (4.*(-sred-0.25)*(-sred-0.75)/(sred*sred)+(-8.*sred-3.)/sred)*ImT(sred);       
      }
    else if (region == high)
      {  // high energy limit
  *re = - ( 1.5 +  
      1.5* (ured-sred)/tred * log(-ured/sred) +
      2. * ( 1. - 0.75 * sred * ured / (tred*tred) ) *  pow( log( - ured/sred),2)   +
      2. * tred * tred * ( log(4.*sred)*log(-4.*tred)/(sred*tred)+
               log(4.*sred)*log(-4.*ured)/(sred*ured)+
               log(-4.*ured)*log(-4.*tred)/(ured*tred) )  
      );
  
  *im = - ( 1.5 * (sred-ured)/tred * (- PI)+ 
      2. * ( 1. - 0.75 * sred*ured / (tred*tred) ) * PI * 2. * log(-ured/sred) + 
      2. * (-PI) * tred*tred*  ( log(-4.*ured) / (ured*sred) +
               log(-4.*tred) / (tred*sred)  )
      )  ;
      }
    else{    Mxxxx_vector(tred,sred,re,im); }
  }
  return;

};

void Mpmpm_vector(double sred, double tred, double *re, double *im, int exclude_loops){

  double ured=-tred-sred;

  if(exclude_loops==2||exclude_loops==3) {*re=0; *im=0;}
  else{ 
    int region = limits (sred,tred,ured);
    if( region == low )
      {     // EFT limit 
  *re= -4.*(4.*(-5./32.)  +3.*(27./40.) )  *ured*ured ;
  *im=0;  
      }
    else if( region == forward)
      {                // Forward limit 
  *re=-3./2.-8.*(-sred-0.25)*(-sred-0.75)/sred*ReB(-sred)+
    (8.*(-sred-0.25)*(-sred-0.75)/sred+3.)*ReB(sred)+
    4.*(-sred-0.25)*(-sred-0.75)/(sred*sred)*ReT(-sred)+
    (4.*(-sred-0.25)*(-sred-0.75)/(sred*sred)+(-8.*sred-3.)/sred)*ReT(sred);
  *im=-8.*(-sred-0.25)*(-sred-0.75)/sred*ImB(-sred)+
    (8.*(-sred-0.25)*(-sred-0.75)/sred+3.)*ImB(sred)+
    4.*(-sred-0.25)*(-sred-0.75)/(sred*sred)*ImT(-sred)+
    (4.*(-sred-0.25)*(-sred-0.75)/(sred*sred)+(-8.*sred-3.)/sred)*ImT(sred);       
      }
    else if( region == backward )
      {                // Backward limit 
  *re=0.; *im=0.;
      }
    else if ( region == high )
      {         // high energy limit
  *re = - ( 1.5 +  
      1.5* (tred-sred)/ured * log(-tred/sred) +
      2. * ( 1. - 0.75 * sred * tred / (ured*ured) ) *  pow( log( - tred/sred),2)   +
      2. * ured * ured * ( log(4.*sred)*log(-4.*tred)/(sred*tred)+
                log(4.*sred)*log(-4.*ured)/(sred*ured)+
                log(-4.*ured)*log(-4.*tred)/(ured*tred) )
      );
  *im = - ( 1.5 * (sred-tred)/ured * (- PI)+ 
      2. * ( 1. - 0.75 * sred*tred / (ured*ured) ) * PI * 2. * log(-tred/sred) + 
      2. * (-PI) * ured*ured*  ( log(-4.*ured) / (ured*sred) +
               log(-4.*tred) / (tred*sred)  )
      )  ;  
      }
    else{   Mxxxx_vector(ured,tred,re,im); }
    
  }
  return;

};

void Mpppm_vector(double sred, double tred, double * re, double * im, int exclude_loops){

double ured=-tred-sred;


if(exclude_loops==2||exclude_loops==3) {*re=0; *im=0;}
else{ 

  //if(sred<0.001){ // EFT limit
  //    *re= 0.; *im=0.;}
  //else if(sred<10000. && sred>0.001 && (-tred<0.0001*sred||-ured<0.0001*sred ))
  //{                // Forward and backward limit 
  //            *re= 0.; *im=0.;}
  //else{

   Mpppm_fermion(sred,tred,re,im,exclude_loops);
  *re *= -1.5;
  *im *= -1.5;
  //}

    }

  return;
};

void Mppmm_vector(double sred, double tred, double * re, double * im, int exclude_loops){

  double ured = -sred-tred;

if(exclude_loops==2||exclude_loops==3) {*re=0; *im=0;}
else{   

  //if(sred<0.001){ // EFT limit 
  //    *re= -4.*(4.*(-5./32.)  +(27./40.) )*(sred*sred+tred*tred+ured*ured); *im=0.;}
  //else if(sred<10000. && sred>0.001 && (-tred<0.0001*sred||-ured<0.0001*sred ))
  //{                // Forward and backward limit 
  //*re=1./(2.*sred*sred)*( -2.*sred*sred-2.*sred*ReB(sred)+2.*sred*ReB(-sred)-ReT(sred)-ReT(-sred)  ) ; 
  //*im=1./(2.*sred*sred)*(              -2.*sred*ImB(sred)+2.*sred*ImB(-sred)-ImT(sred)-ImT(-sred)  );
  //*re *= -1.5;
  //*im *= -1.5;
  //}
  //else{ 

  Mppmm_fermion(sred,tred,re,im,exclude_loops);
  *re *= -1.5;
  *im *= -1.5;
  // }

   }

  return;

};


/// Neutral resonances

// Define the s-dependent width
double width_gen(double s, double m, double f0, double w_const, double a2){
 return 1./m*( a2*s*s/ (4*PI*f0*f0)) +w_const;
};

void Mxxxx_spin0even(double x, double y, double m, double f0, double w_const, double a2, double * re, double * im){
  // some auxiliary function used in Mpppp, Mpmpm, Mpmmp.

  *re=0;
  *im=0;

  double z = - x - y;

  *re += -4./(f0*f0) * x*x/((x-m*m)*(x-m*m) +m*m*width_gen(x,m,f0,w_const,a2)*width_gen(x,m,f0,w_const,a2) ) * (x-m*m);
  *im += -4./(f0*f0) * x*x/((x-m*m)*(x-m*m) +m*m*width_gen(x,m,f0,w_const,a2)*width_gen(x,m,f0,w_const,a2) ) *(- m*width_gen(x,m,f0,w_const,a2));

  return;
};

void Mpppp_spin0even(double s, double t, double m, double f0, double w_const, double a2, double *re, double *im){

double u=-s-t;
   Mxxxx_spin0even(s,t,m,f0,w_const,a2,re,im);
//std::cout <<"re pppp normal" <<"\t" << *re  << std::endl;
  return;
};

void Mpmmp_spin0even(double s, double t, double m, double f0, double w_const, double a2, double *re, double *im){
double u=-s-t;
   Mxxxx_spin0even(t,s,m,f0,w_const,a2,re,im);

  return;
};

void Mpmpm_spin0even(double s, double t, double m, double f0, double w_const, double a2, double *re, double *im){
double u=-s-t;
   Mxxxx_spin0even(u,t,m,f0,w_const,a2,re,im);
  return;
};

void Mppmm_spin0even(double s, double t, double m, double f0, double w_const, double a2, double *re, double *im){
  double u=-s-t;

  *re=0;
  *im=0;

  *re += -4./(f0*f0) * s*s/((s-m*m)*(s-m*m) +m*m*width_gen(s,m,f0,w_const,a2)*width_gen(s,m,f0,w_const,a2) ) * (s-m*m);
  *re += -4./(f0*f0) * t*t/((t-m*m)*(t-m*m) +m*m*width_gen(s,m,f0,w_const,a2)*width_gen(s,m,f0,w_const,a2) ) * (t-m*m);
  *re += -4./(f0*f0) * u*u/((u-m*m)*(u-m*m) +m*m*width_gen(s,m,f0,w_const,a2)*width_gen(s,m,f0,w_const,a2) ) * (u-m*m);

  *im += -4./(f0*f0) * s*s/((s-m*m)*(s-m*m) +m*m*width_gen(s,m,f0,w_const,a2)*width_gen(s,m,f0,w_const,a2) ) *(- m*width_gen(s,m,f0,w_const,a2));
  *im += -4./(f0*f0) * t*t/((t-m*m)*(t-m*m) +m*m*width_gen(s,m,f0,w_const,a2)*width_gen(s,m,f0,w_const,a2) ) *(- m*width_gen(s,m,f0,w_const,a2));
  *im += -4./(f0*f0) * u*u/((u-m*m)*(u-m*m) +m*m*width_gen(s,m,f0,w_const,a2)*width_gen(s,m,f0,w_const,a2) ) *(- m*width_gen(s,m,f0,w_const,a2));

//std::cout <<"Mppmm " <<"\t" << *re  << std::endl;
  return ;
};

void Mpppm_spin0even(double s, double t, double m, double f0, double w_const, double a2, double *re, double *im){
  double u=-s-t;
   *re=0;
   *im=0;

  return ;

};

void Mpppp_eft(double zeta1, double zeta2, double s, double t, double *re, double *im){
double u=-s-t;
*re= -1./4.*(4.*zeta1+3*zeta2)*s*s   ;
*im=0;
  return;
};

void Mpmmp_eft(double zeta1, double zeta2, double s, double t, double *re, double *im){
double u=-s-t;
*re= -1./4.*(4.*zeta1+3*zeta2)*t*t   ;
*im=0;
  return;
};

void Mpmpm_eft(double zeta1, double zeta2, double s, double t, double *re, double *im){
double u=-s-t;
*re= -1./4.*(4.*zeta1+3*zeta2)*u*u   ;
*im=0;
  return;
};

void Mpppm_eft(double zeta1, double zeta2, double s, double t, double *re, double *im){
double u=-s-t;
*re= 0 ;
*im=0;
  return;
};

void Mppmm_eft(double zeta1, double zeta2, double s, double t, double *re, double *im){
double u=-s-t;
*re= -1./4.*(4.*zeta1+zeta2)*(s*s+t*t+u*u)   ;
*im=0;
  return;
};

/// Z Z final state from OZ operator
/// Z Z final state from OZ operator

void MZZxxxx_spin0even(double x, double y, double m, double f0, double f0Z, double w_const, double a2, double * re, double * im){
  // some auxiliary function used in Mpppp, Mpmpm, Mpmmp.

  *re=0;
  *im=0;

  double z = - x - y;

  *re += -4./(f0*f0Z) * x*(x-2*mZ*mZ)/((x-m*m)*(x-m*m) +m*m*width_gen(x,m,f0,w_const,a2)*width_gen(x,m,f0,w_const,a2) ) * (x-m*m);
  *im += -4./(f0*f0Z) * x*(x-2*mZ*mZ)/((x-m*m)*(x-m*m) +m*m*width_gen(x,m,f0,w_const,a2)*width_gen(x,m,f0,w_const,a2) ) *(- m*width_gen(x,m,f0,w_const,a2));

  return;
};

void MZZpppp_spin0even(double s, double t, double m, double f0, double f0Z, double w_const, double a2, double *re, double *im){

double u=-s-t;
   MZZxxxx_spin0even(s,t,m,f0,f0Z,w_const,a2,re,im);
//std::cout <<"Mpppp " <<"\t" << *re  << std::endl;
  return;
};

void MZZpmmp_spin0even(double s, double t, double m, double f0, double f0Z, double w_const, double a2, double *re, double *im){
double u=-s-t;
   *re=0;
   *im=0;
//std::cout <<"Mpmmp " <<"\t" << *re  << std::endl;
  return;
};

void MZZpmpm_spin0even(double s, double t, double m, double f0, double f0Z, double w_const, double a2, double *re, double *im){
double u=-s-t;
   *re=0;
   *im=0;
//std::cout <<"Mpmpm " <<"\t" << *re  << std::endl;
  return;
};

void MZZppmm_spin0even(double s, double t, double m, double f0, double f0Z, double w_const, double a2, double *re, double *im){
  double u=-s-t;

  *re=0;
  *im=0;

  *re += -4./(f0*f0Z) * s*(s-2*mZ*mZ)/((s-m*m)*(s-m*m) +m*m*width_gen(s,m,f0,w_const,a2)*width_gen(s,m,f0,w_const,a2) ) * (s-m*m);

  *im += -4./(f0*f0Z) * s*(s-2*mZ*mZ)/((s-m*m)*(s-m*m) +m*m*width_gen(s,m,f0,w_const,a2)*width_gen(s,m,f0,w_const,a2) ) *(- m*width_gen(s,m,f0,w_const,a2));

//std::cout <<"Mppmm " <<"\t" << *re  << std::endl;
  return ;
};

void MZZpppm_spin0even(double s, double t, double m, double f0, double f0Z, double w_const, double a2, double *re, double *im){
  double u=-s-t;
   *re=0;
   *im=0;

  return ;
};

void MZZpp00_spin0even(double s, double t, double m, double f0, double f0Z, double w_const, double a2, double *re, double *im){
  double u=-s-t;

  *re=0;
  *im=0;

  *re += -4./(f0*f0Z) * s*(2*mZ*mZ)/((s-m*m)*(s-m*m) +m*m*width_gen(s,m,f0,w_const,a2)*width_gen(s,m,f0,w_const,a2) ) * (s-m*m);
 
  *im += -4./(f0*f0Z) * s*(2*mZ*mZ)/((s-m*m)*(s-m*m) +m*m*width_gen(s,m,f0,w_const,a2)*width_gen(s,m,f0,w_const,a2) ) *(- m*width_gen(s,m,f0,w_const,a2));
 
//std::cout <<"Mpp00 " <<"\t" << *re  << std::endl;
  return ;
};

/// Z gam final state from OZgam operator

void MZgxxxx_spin0even(double x, double y, double m, double f0, double f0Zg, double w_const, double a2, double * re, double * im){
  // some auxiliary function used in Mpppp, Mpmpm, Mpmmp.

  *re=0;
  *im=0;

  double z = - x - y;

  *re += -4./(f0*f0Zg) * x*(x-mZ*mZ)/((x-m*m)*(x-m*m) +m*m*width_gen(x,m,f0,w_const,a2)*width_gen(x,m,f0,w_const,a2) ) * (x-m*m);
  *im += -4./(f0*f0Zg) * x*(x-mZ*mZ)/((x-m*m)*(x-m*m) +m*m*width_gen(x,m,f0,w_const,a2)*width_gen(x,m,f0,w_const,a2) ) *(- m*width_gen(x,m,f0,w_const,a2));

  return;
};

void MZgpppp_spin0even(double s, double t, double m, double f0, double f0Zg, double w_const, double a2, double *re, double *im){

double u=-s-t;
   MZgxxxx_spin0even(s,t,m,f0,f0Zg,w_const,a2,re,im);
//std::cout <<"Mpppp " <<"\t" << *re  << std::endl;
  return;
};

void MZgpmmp_spin0even(double s, double t, double m, double f0, double f0Zg, double w_const, double a2, double *re, double *im){
double u=-s-t;
    *re=0;
    *im=0;
//std::cout <<"Mpmmp " <<"\t" << *re  << std::endl;
  return;
};

void MZgpmpm_spin0even(double s, double t, double m, double f0, double f0Zg, double w_const, double a2, double *re, double *im){
double u=-s-t;
   *re=0;
   *im=0;
//std::cout <<"Mpmpm " <<"\t" << *re  << std::endl;
  return;
};

void MZgppmm_spin0even(double s, double t, double m, double f0, double f0Zg, double w_const, double a2, double *re, double *im){
  double u=-s-t;

  *re=0;
  *im=0;

  *re += -4./(f0*f0Zg) * s*(s-mZ*mZ)/((s-m*m)*(s-m*m) +m*m*width_gen(s,m,f0,w_const,a2)*width_gen(s,m,f0,w_const,a2) ) * (s-m*m);
 

  *im += -4./(f0*f0Zg) * s*(s-mZ*mZ)/((s-m*m)*(s-m*m) +m*m*width_gen(s,m,f0,w_const,a2)*width_gen(s,m,f0,w_const,a2) ) *(- m*width_gen(s,m,f0,w_const,a2));
 

//std::cout <<"Mppmm " <<"\t" << *re  << std::endl;
  return ;
};

void MZgpppm_spin0even(double s, double t, double m, double f0, double f0Zg, double w_const, double a2, double *re, double *im){
  double u=-s-t;
   *re=0;
   *im=0;

  return ;
};

void MZgpp00_spin0even(double s, double t, double m, double f0, double f0Zg, double w_const, double a2, double *re, double *im){
  double u=-s-t;

  *re=0;
  *im=0;


//std::cout <<"Mpp00 " <<"\t" << *re  << std::endl;
  return ;
};

/// W+W- final state from OW operator
void MWWxxxx_spin0even(double x, double y, double m, double f0, double f0W, double w_const, double a2, double * re, double * im){
  // some auxiliary function used in Mpppp, Mpmpm, Mpmmp.

  *re=0;
  *im=0;

  double z = - x - y;

  *re += -4./(f0*f0W) * x*(x-2*mW*mW)/((x-m*m)*(x-m*m) +m*m*width_gen(x,m,f0,w_const,a2)*width_gen(x,m,f0,w_const,a2) ) * (x-m*m);
  *im += -4./(f0*f0W) * x*(x-2*mW*mW)/((x-m*m)*(x-m*m) +m*m*width_gen(x,m,f0,w_const,a2)*width_gen(x,m,f0,w_const,a2) ) *(- m*width_gen(x,m,f0,w_const,a2));

  return;
};

void MWWpppp_spin0even(double s, double t, double m, double f0, double f0W, double w_const, double a2, double *re, double *im){

double u=-s-t;
   MWWxxxx_spin0even(s,t,m,f0,f0W,w_const,a2,re,im);
//std::cout <<"Mpppp " <<"\t" << *re  << std::endl;
  return;
};

void MWWpmmp_spin0even(double s, double t, double m, double f0, double f0W, double w_const, double a2, double *re, double *im){
double u=-s-t;
   *re=0;
   *im=0;
//std::cout <<"Mpmmp " <<"\t" << *re  << std::endl;
  return;
};

void MWWpmpm_spin0even(double s, double t, double m, double f0, double f0W, double w_const, double a2, double *re, double *im){
double u=-s-t;
   *re=0;
   *im=0;
//std::cout <<"Mpmpm " <<"\t" << *re  << std::endl;
  return;
};

void MWWppmm_spin0even(double s, double t, double m, double f0, double f0W, double w_const, double a2, double *re, double *im){
  double u=-s-t;

  *re=0;
  *im=0;

  *re += -4./(f0*f0W) * s*(s-2*mW*mW)/((s-m*m)*(s-m*m) +m*m*width_gen(s,m,f0,w_const,a2)*width_gen(s,m,f0,w_const,a2) ) * (s-m*m);
  
  *im += -4./(f0*f0W) * s*(s-2*mW*mW)/((s-m*m)*(s-m*m) +m*m*width_gen(s,m,f0,w_const,a2)*width_gen(s,m,f0,w_const,a2) ) *(- m*width_gen(s,m,f0,w_const,a2));
  

//std::cout <<"Mppmm " <<"\t" << *re  << std::endl;
  return ;
};

void MWWpppm_spin0even(double s, double t, double m, double f0, double f0W, double w_const, double a2, double *re, double *im){
  double u=-s-t;
   *re=0;
   *im=0;

  return ;
};

void MWWpp00_spin0even(double s, double t, double m, double f0, double f0W, double w_const, double a2, double *re, double *im){
  double u=-s-t;

  *re=0;
  *im=0;

  *re += -4./(f0*f0W) * s*(2*mW*mW)/((s-m*m)*(s-m*m) +m*m*width_gen(s,m,f0,w_const,a2)*width_gen(s,m,f0,w_const,a2) ) * (s-m*m);
 
  *im += -4./(f0*f0W) * s*(2*mW*mW)/((s-m*m)*(s-m*m) +m*m*width_gen(s,m,f0,w_const,a2)*width_gen(s,m,f0,w_const,a2) ) *(- m*width_gen(s,m,f0,w_const,a2));
 

//std::cout <<"Mpp00 " <<"\t" << *re  << std::endl;
  return ;
};
// Spin 2

void Mxxxx_spin2(double x, double y, double m, double f0, double w_const, double a2, double * re, double * im){
  // some auxiliary function used in Mpppp, Mpmpm, Mpmmp.

  *re=0;
  *im=0;

  double z = - x - y;

  *re += -1./(f0*f0) * x*x/((y-m*m)*(y-m*m) +m*m*width_gen(y,m,f0,w_const,a2)*width_gen(y,m,f0,w_const,a2) ) * (y-m*m);
  *im += -1./(f0*f0) * x*x/((y-m*m)*(y-m*m) +m*m*width_gen(y,m,f0,w_const,a2)*width_gen(y,m,f0,w_const,a2) ) * (- m*width_gen(y,m,f0,w_const,a2));
  *re += -1./(f0*f0) * x*x/((z-m*m)*(z-m*m) +m*m*width_gen(z,m,f0,w_const,a2)*width_gen(z,m,f0,w_const,a2) ) * (z-m*m);
  *im += -1./(f0*f0) * x*x/((z-m*m)*(z-m*m) +m*m*width_gen(z,m,f0,w_const,a2)*width_gen(z,m,f0,w_const,a2) ) * (- m*width_gen(z,m,f0,w_const,a2));

  return;
};

void Mpppp_spin2(double s, double t, double m, double f0, double w_const, double a2, double *re, double *im){

double u=-s-t;
   Mxxxx_spin2(s,t,m,f0,w_const,a2,re,im);
//std::cout <<"re pppp normal" <<"\t" << *re  << std::endl;
  return;
};

void Mpmmp_spin2(double s, double t, double m, double f0, double w_const, double a2, double *re, double *im){
double u=-s-t;
   Mxxxx_spin2(t,s,m,f0,w_const,a2,re,im);

  return;
};

void Mpmpm_spin2(double s, double t, double m, double f0, double w_const, double a2, double *re, double *im){
double u=-s-t;
   Mxxxx_spin2(u,t,m,f0,w_const,a2,re,im);
  return;
};

void Mppmm_spin2(double s, double t, double m, double f0, double w_const, double a2, double *re, double *im){
  double u=-s-t;

  *re=0;
  *im=0;

  
  return ;
};

void Mpppm_spin2(double s, double t, double m, double f0, double w_const, double a2, double *re, double *im){
  double u=-s-t;
   *re=0;
   *im=0;

  return ;

};

/// gg final state from OG operator

void Mggxxxx_spin0even(double x, double y, double m, double f0, double f0g, double w_const, double a2, double * re, double * im){
  // some auxiliary function used in Mpppp, Mpmpm, Mpmmp.

  *re=0;
  *im=0;

  double z = - x - y;

  *re += -4./(f0*f0g) * x*x/((x-m*m)*(x-m*m) +m*m*width_gen(x,m,f0,w_const,a2)*width_gen(x,m,f0,w_const,a2) ) * (x-m*m);
  *im += -4./(f0*f0g) * x*x/((x-m*m)*(x-m*m) +m*m*width_gen(x,m,f0,w_const,a2)*width_gen(x,m,f0,w_const,a2) ) *(- m*width_gen(x,m,f0,w_const,a2));

  return;
};

void Mggpppp_spin0even(double s, double t, double m, double f0, double f0g, double w_const, double a2, double *re, double *im){

double u=-s-t;
   Mggxxxx_spin0even(s,t,m,f0,f0g,w_const,a2,re,im);
//std::cout <<"Mpppp " <<"\t" << *re  << std::endl;
  return;
};

void Mggpmmp_spin0even(double s, double t, double m, double f0, double f0g, double w_const, double a2, double *re, double *im){
double u=-s-t;
   *re=0;
   *im=0;
//std::cout <<"Mpmmp " <<"\t" << *re  << std::endl;
  return;
};

void Mggpmpm_spin0even(double s, double t, double m, double f0, double f0g, double w_const, double a2, double *re, double *im){
double u=-s-t;
   *re=0;
   *im=0;
//std::cout <<"Mpmpm " <<"\t" << *re  << std::endl;
  return;
};

void Mggppmm_spin0even(double s, double t, double m, double f0, double f0g, double w_const, double a2, double *re, double *im){
  double u=-s-t;

  *re=0;
  *im=0;

  *re += -4./(f0*f0g) * s*s/((s-m*m)*(s-m*m) +m*m*width_gen(s,m,f0,w_const,a2)*width_gen(s,m,f0,w_const,a2) ) * (s-m*m);

  *im += -4./(f0*f0g) * s*s/((s-m*m)*(s-m*m) +m*m*width_gen(s,m,f0,w_const,a2)*width_gen(s,m,f0,w_const,a2) ) *(- m*width_gen(s,m,f0,w_const,a2));

//std::cout <<"Mppmm " <<"\t" << *re  << std::endl;
  return ;
};

void Mggpppm_spin0even(double s, double t, double m, double f0, double f0g, double w_const, double a2, double *re, double *im){
  double u=-s-t;
   *re=0;
   *im=0;

  return ;
};

/// hh final state from OH operator

void MhhOHxxxx_spin0even(double x, double y, double m, double f0, double f0H, double w_const, double a2, double * re, double * im){
  // some auxiliary function used in Mpppp, Mpmpm, Mpmmp.

  *re=0;
  *im=0;

  double z = - x - y;

  *re += 4./(f0*f0H) * x*(mh*mh)/((x-m*m)*(x-m*m) +m*m*width_gen(x,m,f0,w_const,a2)*width_gen(x,m,f0,w_const,a2) ) * (x-m*m);
  *im += 4./(f0*f0H) * x*(mh*mh)/((x-m*m)*(x-m*m) +m*m*width_gen(x,m,f0,w_const,a2)*width_gen(x,m,f0,w_const,a2) ) *(- m*width_gen(x,m,f0,w_const,a2));

  return;
};

void MhhOHpppp_spin0even(double s, double t, double m, double f0, double f0H, double w_const, double a2, double *re, double *im){

double u=-s-t;
   MhhOHxxxx_spin0even(s,t,m,f0,f0H,w_const,a2,re,im);
//std::cout <<"Mpppp " <<"\t" << *re  << std::endl;
  return;
};

void MhhOHpmmp_spin0even(double s, double t, double m, double f0, double f0H, double w_const, double a2, double *re, double *im){
double u=-s-t;
    *re=0;
   *im=0;
//std::cout <<"Mpmmp " <<"\t" << *re  << std::endl;
  return;
};

void MhhOHpmpm_spin0even(double s, double t, double m, double f0, double f0H, double w_const, double a2, double *re, double *im){
double u=-s-t;
   *re=0;
   *im=0;
//std::cout <<"Mpmpm " <<"\t" << *re  << std::endl;
  return;
};

void MhhOHppmm_spin0even(double s, double t, double m, double f0, double f0H, double w_const, double a2, double *re, double *im){
  double u=-s-t;

  *re=0;
  *im=0;

  *re += 4./(f0*f0H) * s*(mh*mh)/((s-m*m)*(s-m*m) +m*m*width_gen(s,m,f0,w_const,a2)*width_gen(s,m,f0,w_const,a2) ) * (s-m*m);
 
  *im += 4./(f0*f0H) * s*(mh*mh)/((s-m*m)*(s-m*m) +m*m*width_gen(s,m,f0,w_const,a2)*width_gen(s,m,f0,w_const,a2) ) *(- m*width_gen(s,m,f0,w_const,a2));

//std::cout <<"Mppmm " <<"\t" << *re  << std::endl;
  return ;
};

void MhhOHpppm_spin0even(double s, double t, double m, double f0, double f0H, double w_const, double a2, double *re, double *im){
  double u=-s-t;
   *re=0;
   *im=0;

  return ;
};

void MhhOHpp00_spin0even(double s, double t, double m, double f0, double f0H, double w_const, double a2, double *re, double *im){
  double u=-s-t;

  *re=0;
  *im=0;

  *re += 2./(f0*f0H) * s*(s-2*mh*mh)/((s-m*m)*(s-m*m) +m*m*width_gen(s,m,f0,w_const,a2)*width_gen(s,m,f0,w_const,a2) ) * (s-m*m);

  *im += 2./(f0*f0H) * s*(s-2*mh*mh)/((s-m*m)*(s-m*m) +m*m*width_gen(s,m,f0,w_const,a2)*width_gen(s,m,f0,w_const,a2) ) *(- m*width_gen(s,m,f0,w_const,a2));

//std::cout <<"Mpp00 " <<"\t" << *re  << std::endl;
  return ;
};

//EFT for AA->AZ

void Mpppp_eft_AZ(double zeta1, double s, double t, double *re, double *im){
double u=-s-t;
*re= 2.*zeta1*s*(s-mZ*mZ);
*im=0;
  return;
};

void Mpmmp_eft_AZ(double zeta1, double s, double t, double *re, double *im){
double u=-s-t;
*re= 2*zeta1*s/(s-mZ*mZ)*t*t   ;
*im=0;
  return;
};

void Mpmpm_eft_AZ(double zeta1, double s, double t, double *re, double *im){
double u=-s-t;
*re= 2*zeta1*s/(s-mZ*mZ)*u*u   ;
*im=0;
  return;
};

void Mppmm_eft_AZ(double zeta1, double s, double t, double *re, double *im){
double u=-s-t;
*re= 2*zeta1*(s*(s-mZ*mZ)+s/(s-mZ*mZ)*(t*t+u*u));
*im=0;
  return;
};

void Mppmp_eft_AZ(double zeta1, double s, double t, double *re, double *im){
double u=-s-t;
*re= 2*zeta1*(s*(s-mZ*mZ)+s/(s-mZ*mZ)*(t*t+u*u));
*im=0;
  return;
};

void Mpppm_eft_AZ(double zeta1, double s, double t, double *re, double *im){
double u=-s-t;
*re= 0;
*im=0;
  return;
};


/// Z Z final state from OH operator

void MZZOHxxxx_spin0even(double x, double y, double m, double f0, double f0H, double w_const, double a2, double * re, double * im){
  // some auxiliary function used in Mpppp, Mpmpm, Mpmmp.

  *re=0;
  *im=0;

  double z = - x - y;

  *re += 4./(f0*f0H) * x*(mZ*mZ)/((x-m*m)*(x-m*m) +m*m*width_gen(x,m,f0,w_const,a2)*width_gen(x,m,f0,w_const,a2) ) * (x-m*m);
  *im += 4./(f0*f0H) * x*(mZ*mZ)/((x-m*m)*(x-m*m) +m*m*width_gen(x,m,f0,w_const,a2)*width_gen(x,m,f0,w_const,a2) ) *(- m*width_gen(x,m,f0,w_const,a2));

  return;
};

void MZZOHpppp_spin0even(double s, double t, double m, double f0, double f0H, double w_const, double a2, double *re, double *im){

double u=-s-t;
   MZZOHxxxx_spin0even(s,t,m,f0,f0H,w_const,a2,re,im);
//std::cout <<"Mpppp " <<"\t" << *re  << std::endl;
  return;
};

void MZZOHpmmp_spin0even(double s, double t, double m, double f0, double f0H, double w_const, double a2, double *re, double *im){
double u=-s-t;
    *re=0;
   *im=0;
//std::cout <<"Mpmmp " <<"\t" << *re  << std::endl;
  return;
};

void MZZOHpmpm_spin0even(double s, double t, double m, double f0, double f0H, double w_const, double a2, double *re, double *im){
double u=-s-t;
   *re=0;
   *im=0;
//std::cout <<"Mpmpm " <<"\t" << *re  << std::endl;
  return;
};

void MZZOHppmm_spin0even(double s, double t, double m, double f0, double f0H, double w_const, double a2, double *re, double *im){
  double u=-s-t;

  *re=0;
  *im=0;

  *re += 4./(f0*f0H) * s*(mZ*mZ)/((s-m*m)*(s-m*m) +m*m*width_gen(s,m,f0,w_const,a2)*width_gen(s,m,f0,w_const,a2) ) * (s-m*m);
  
  *im += 4./(f0*f0H) * s*(mZ*mZ)/((s-m*m)*(s-m*m) +m*m*width_gen(s,m,f0,w_const,a2)*width_gen(s,m,f0,w_const,a2) ) *(- m*width_gen(s,m,f0,w_const,a2));

//std::cout <<"Mppmm " <<"\t" << *re  << std::endl;
  return ;
};

void MZZOHpppm_spin0even(double s, double t, double m, double f0, double f0H, double w_const, double a2, double *re, double *im){
  double u=-s-t;
   *re=0;
   *im=0;

  return ;
};

void MZZOHpp00_spin0even(double s, double t, double m, double f0, double f0H, double w_const, double a2, double *re, double *im){
  double u=-s-t;

  *re=0;
  *im=0;

  *re += 2./(f0*f0H) * s*(s-2*mZ*mZ)/((s-m*m)*(s-m*m) +m*m*width_gen(s,m,f0,w_const,a2)*width_gen(s,m,f0,w_const,a2) ) * (s-m*m);
 
  *im += 2./(f0*f0H) * s*(s-2*mZ*mZ)/((s-m*m)*(s-m*m) +m*m*width_gen(s,m,f0,w_const,a2)*width_gen(s,m,f0,w_const,a2) ) *(- m*width_gen(s,m,f0,w_const,a2));
 

//std::cout <<"Mpp00 " <<"\t" << *re  << std::endl;
  return ;
};
