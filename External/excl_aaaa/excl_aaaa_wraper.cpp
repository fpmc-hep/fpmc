#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define DEBUG 0

 //---M.Saimpert 01/2014--matthias.saimpert@cern.ch--------//
 //--------Coding of the exclusive aa->aa process----------//
 //--formulas from G. von Gersdorff (gersdorff@gmail.com)--// 
 //--formulas from S. Fichet  sylvain.fichet@gmail.com--)--// 
 //modification of the comphep external module used for aQGC/

 #ifdef __cplusplus /* If this is a C++ compiler, use C linkage */
 extern "C" {
 #endif

// // routines called by fpmc 
 //SM
 void sm_sqme_aaaa_c__(double *_amp2, double *_s, double *_t, int *_exclude_loops);

 void sm_sqme_aaaa_c_(double *_amp2, double *_s, double *_t, int *_exclude_loops);

 //Exotic fermions
 void bsmf_sqme_aaaa_c__(double *_amp2, double *_s, double *_t, 
                         int *_exclude_loops_SM, int *_exclude_loops_EX,
                         double *_m, double *_q, double *_n);

 void bsmf_sqme_aaaa_c_(double *_amp2, double *_s, double *_t, 
                        int *_exclude_loops_SM, int *_exclude_loops_EX,
                        double *_m, double *_q, double *_n);
 
 //Exotic vectors
 void bsmv_sqme_aaaa_c__(double *_amp2, double *_s, double *_t,
                         int *_exclude_loops_SM, int *_exclude_loops_EX,
                         double *_m, double *_q, double *_n);

 void bsmv_sqme_aaaa_c_(double *_amp2, double *_s, double *_t,
                        int *_exclude_loops_SM, int *_exclude_loops_EX,
                        double *_m, double *_q, double *_n);
 #ifdef __cplusplus
 }
 #endif
 


 
// forward declaration
   namespace sm_aaaa {
      extern double sqme(double, double, int);
   }; // namespace sm_aaaa

   namespace bsmf_aaaa {
      extern double sqme(double, double, int, int, double, double, double);
   }; // namespace bsmf_aaaa

   namespace bsmv_aaaa {
      extern double sqme(double, double, int, int, double, double, double);
   }; // namespace bsmv_aaaa





//wrapper for g77
void sm_sqme_aaaa_c__(double *_amp2, double *_s, double *_t, int *_exclude_loops) {

    sm_sqme_aaaa_c_(_amp2, _s, _t, _exclude_loops);
 }

void sm_sqme_aaaa_c_(double *_amp2, double *_s, double *_t, int *_exclude_loops) {
//////////////////////////////////////////////////////////////////// 
// SM aaaa including:
//  fermion+W loop (exclude_loops=0), no fermions (1), no W (2)
//////////////////////////////////////////////////////////////////// 

   double amp2 = sm_aaaa::sqme(*_s, *_t, *_exclude_loops);
   *_amp2 = amp2;
}



//wrapper for g77
void bsmf_sqme_aaaa_c__(double *_amp2, double *_s, double *_t, 
                   int *_exclude_loops_SM, int *_exclude_loops_EX,
                   double *_m, double *_q, double *_n) {

    bsmf_sqme_aaaa_c_(_amp2, _s, _t, _exclude_loops_SM, _exclude_loops_EX, _m, _q, _n);
 }

void bsmf_sqme_aaaa_c_(double *_amp2, double *_s, double *_t,
                  int *_exclude_loops_SM, int *_exclude_loops_EX,
                  double *_m, double *_q, double *_n) {
//////////////////////////////////////////////////////////////////// 
// BSM aaaa with exotic fermions of mass m charge q multiplicity n
// including interference with SM (SM fermion loops excluded is default)
//////////////////////////////////////////////////////////////////// 

   double amp2 = bsmf_aaaa::sqme(*_s, *_t, *_exclude_loops_SM, *_exclude_loops_EX, *_m, *_q, *_n);
   *_amp2 = amp2;
}



//wrapper for g77
void bsmv_sqme_aaaa_c__(double *_amp2, double *_s, double *_t,
                   int *_exclude_loops_SM, int *_exclude_loops_EX,
                   double *_m, double *_q, double *_n) {

    bsmv_sqme_aaaa_c_(_amp2, _s, _t, _exclude_loops_SM, _exclude_loops_EX, _m, _q, _n);
 }

void bsmv_sqme_aaaa_c_(double *_amp2, double *_s, double *_t,
                  int *_exclude_loops_SM, int *_exclude_loops_EX,
                  double *_m, double *_q, double *_n) {
//////////////////////////////////////////////////////////////////// 
// BSM aaaa with exotic vector of mass m charge q multiplicity n
// including interferences with SM
//////////////////////////////////////////////////////////////////// 

   double amp2 = bsmv_aaaa::sqme(*_s, *_t, *_exclude_loops_SM, *_exclude_loops_EX, *_m, *_q, *_n);
   *_amp2 = amp2;
}  
