#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define DEBUG 0

//---M.Saimpert 01/2014--matthias.saimpert@cern.ch--------//
//--------Coding of the exclusive aa->aa process----------//
//--formulas from G. von Gersdorff (gersdorff@gmail.com)--// 
//--formulas from S. Fichet  sylvain.fichet@gmail.com--)--// 
//modification of the comphep external module used for aaww/aazz//

#ifdef __cplusplus /* If this is a C++ compiler, use C linkage */
extern "C" {
#endif
  // routines called by fpmc 

  void sm_sqme_aaaa_c_(double* _amp2, double* _s, double* _t, int* _exclude_loops); //SM
  void bsmf_sqme_aaaa_c_(double* _amp2, double* _s, double* _t, int* _exclude_loops_SM, int* _exclude_loops_EX, double* _m, double* _q, double* _n); //Exotic fermions
  void bsmv_sqme_aaaa_c_(double* _amp2, double* _s, double* _t, int* _exclude_loops_SM, int* _exclude_loops_EX, double* _m, double* _q, double* _n); //Exotic vectors
  void resonances0even_sqme_aaaa_c_(double* _amp2, double* _s, double* _t, int* _exclude_loops_SM, int* _exclude_loops_EX, double* _m, double* _c, double* _w_c, double* _aa); //Spin0even neutral resonances
  void resonances0even_sqme_zz_c_(double* _amp2, double* _s, double* _t, int* _exclude_loops_SM, int* _exclude_loops_EX, double* _m, double* _c, double* _cz, double* _w_c, double* _aa); //Spin0even neutral resonances
  void resonances0evenoh_sqme_zz_c_(double* _amp2, double* _s, double* _t, int* _exclude_loops_SM, int* _exclude_loops_EX, double* _m, double* _c, double* _cz, double* _w_c, double* _aa); //Spin0even neutral resonances
  void resonances0even_sqme_ww_c_(double* _amp2, double* _s, double* _t, int* _exclude_loops_SM, int* _exclude_loops_EX, double* _m, double* _c, double* _cw, double* _w_c, double* _aa); //Spin0even neutral resonances
  void resonances0even_sqme_hh_c_(double* _amp2, double* _s, double* _t, int* _exclude_loops_SM, int* _exclude_loops_EX, double* _m, double* _c, double* _cw, double* _w_c, double* _aa); //Spin0even neutral resonances
  void resonances0even_sqme_gluglu_c_(double* _amp2, double* _s, double* _t, int* _exclude_loops_SM, int* _exclude_loops_EX, double* _m, double* _c, double* _cw, double* _w_c, double* _aa); //Spin0even neutral resonances
  void resonances0even_sqme_az_c_(double* _amp2, double* _s, double* _t, int* _exclude_loops_SM, int* _exclude_loops_EX, double* _m, double* _c, double* _czg, double* _w_c, double* _aa); //Spin0even neutral resonances
  void resonances2_sqme_aaaa_c_(double* _amp2, double* _s, double* _t, int* _exclude_loops_SM, int* _exclude_loops_EX, double* _m, double* _c, double* _w_c, double* _aa); //Spin2 neutral resonances
  void eft_sqme_aaaa_c_(double* _amp2, double* _s, double* _t, int* _exclude_loops_SM, double* _z1, double* _z2, double* _cutoff); //EFT limit
  void eft_sqme_aaaz_c_( double*, double*, double*, int*, double*, double*, double* );

#ifdef __cplusplus
}
#endif

//////////////////////////////////////////////////////////////////// 
// SM aaaa including:
//  fermion+W loop (exclude_loops=0), no fermions (1), no W (2)
//////////////////////////////////////////////////////////////////// 
namespace sm_aaaa { extern double sqme(double, double, int); }
void sm_sqme_aaaa_c_(double* _amp2, double* _s, double* _t, int* _exclude_loops) {
  *_amp2 = sm_aaaa::sqme(*_s, *_t, *_exclude_loops);
}
#define sm_sqme_aaaa_c__ sm_sqme_aaaa_c_ //wrapper for g77

//////////////////////////////////////////////////////////////////// 
// BSM aaaa with exotic fermions of mass m charge q multiplicity n
// including interference with SM (SM fermion loops excluded is default)
//////////////////////////////////////////////////////////////////// 
namespace bsmf_aaaa { extern double sqme(double, double, int, int, double, double, double); }
void bsmf_sqme_aaaa_c_(double* _amp2, double* _s, double* _t, int* _exclude_loops_SM, int* _exclude_loops_EX, double* _m, double* _q, double* _n) {
  *_amp2 = bsmf_aaaa::sqme(*_s, *_t, *_exclude_loops_SM, *_exclude_loops_EX, *_m, *_q, *_n);
}
#define bsmf_sqme_aaaa_c__ bsmf_sqme_aaaa_c_ //wrapper for g77

//////////////////////////////////////////////////////////////////// 
// BSM aaaa with exotic vector of mass m charge q multiplicity n
// including interferences with SM
//////////////////////////////////////////////////////////////////// 
namespace bsmv_aaaa { extern double sqme(double, double, int, int, double, double, double); }
void bsmv_sqme_aaaa_c_(double* _amp2, double* _s, double* _t, int* _exclude_loops_SM, int* _exclude_loops_EX, double* _m, double* _q, double* _n) {
  *_amp2 = bsmv_aaaa::sqme(*_s, *_t, *_exclude_loops_SM, *_exclude_loops_EX, *_m, *_q, *_n);
}
#define bsmv_sqme_aaaa_c__ bsmv_sqme_aaaa_c_ //wrapper for g77

//////////////////////////////////////////////////////////////////// 
// BSM aaaa with neutral resonance of mass m coupling c width w
// including interferences with SM
//////////////////////////////////////////////////////////////////// 
namespace resonances0even_aaaa { extern double sqme(double, double, int, int, double, double, double, double); }
void resonances0even_sqme_aaaa_c_(double* _amp2, double* _s, double* _t, int* _exclude_loops_SM, int* _exclude_loops_EX, double* _m, double* _c, double* _w_c, double* _aa) {
  *_amp2 = resonances0even_aaaa::sqme(*_s, *_t, *_exclude_loops_SM, *_exclude_loops_EX, *_m, *_c, *_w_c, *_aa);
}
#define resonances0even_sqme_aaaa_c__ resonances0even_sqme_aaaa_c_ //wrapper for g77
 
//////////////////////////////////////////////////////////////////// 
// BSM aaaa with neutral resonance of mass m coupling c width w
// including interferences with SM
//////////////////////////////////////////////////////////////////// 
namespace resonances0even_zz { extern double sqme(double, double, int, int, double, double, double, double, double); }
void resonances0even_sqme_zz_c_(double* _amp2, double* _s, double* _t, int* _exclude_loops_SM, int* _exclude_loops_EX, double* _m, double* _c, double* _cz, double* _w_c, double* _aa) {
  *_amp2 = resonances0even_zz::sqme(*_s, *_t, *_exclude_loops_SM, *_exclude_loops_EX, *_m, *_c, *_cz, *_w_c, *_aa);
}
#define resonances0even_sqme_zz_c__ resonances0even_sqme_zz_c_ //wrapper for g77

namespace resonances0evenoh_zz { extern double sqme(double, double, int, int, double, double, double, double, double); }
void resonances0evenoh_sqme_zz_c_( double* _amp2, double* _s, double* _t, int* _exclude_loops_SM, int* _exclude_loops_EX, double* _m, double* _c, double* _cz, double* _w_c, double* _aa ) {
  *_amp2 = resonances0evenoh_zz::sqme( *_s, *_t, *_exclude_loops_SM, *_exclude_loops_EX, *_m, *_c, *_cz, *_w_c, *_aa );
}
#define resonances0evenoh_sqme_zz_c__ resonances0evenoh_sqme_zz_c_ //wrapper for g77

////////////////////////////////////////////////////////////////////
// BSM aaaa with neutral resonance of mass m coupling c width w
// including interferences with SM
////////////////////////////////////////////////////////////////////
namespace resonances0even_ww { extern double sqme(double, double, int, int, double, double, double, double, double); }
void resonances0even_sqme_ww_c_(double* _amp2, double* _s, double* _t, int* _exclude_loops_SM, int* _exclude_loops_EX, double* _m, double* _c, double* _cw, double* _w_c, double* _aa) {
  *_amp2 = resonances0even_ww::sqme(*_s, *_t, *_exclude_loops_SM, *_exclude_loops_EX, *_m, *_c, *_cw, *_w_c, *_aa);
}
#define resonances0even_sqme_ww_c__ resonances0even_sqme_ww_c_ //wrapper for g77

////////////////////////////////////////////////////////////////////
// BSM aahh with neutral resonance of mass m coupling c width w
// including interferences with SM
////////////////////////////////////////////////////////////////////
namespace resonances0even_hh { extern double sqme(double, double, int, int, double, double, double, double, double); }
void resonances0even_sqme_hh_c_( double* _amp2, double* _s, double* _t, int* _exclude_loops_SM, int* _exclude_loops_EX, double* _m, double* _c, double* _cw, double* _w_c, double* _aa ) {
  *_amp2 = resonances0even_hh::sqme( *_s, *_t, *_exclude_loops_SM, *_exclude_loops_EX, *_m, *_c, *_cw, *_w_c, *_aa );
}
#define resonances0even_sqme_hh_c__ resonances0even_sqme_hh_c_ //wrapper for g77

////////////////////////////////////////////////////////////////////
// BSM gluon-gluon with neutral resonance of mass m coupling c width w
// including interferences with SM
////////////////////////////////////////////////////////////////////
namespace resonances0even_gluglu { extern double sqme(double, double, int, int, double, double, double, double, double); }
void resonances0even_sqme_gluglu_c_( double* _amp2, double* _s, double* _t, int* _exclude_loops_SM, int* _exclude_loops_EX, double* _m, double* _c, double* _cw, double* _w_c, double* _aa ) {
  *_amp2 = resonances0even_gluglu::sqme( *_s, *_t, *_exclude_loops_SM, *_exclude_loops_EX, *_m, *_c, *_cw, *_w_c, *_aa );
}
#define resonances0even_sqme_gluglu_c__ resonances0even_sqme_gluglu_c_ //wrapper for g77

////////////////////////////////////////////////////////////////////
// BSM aaaa with neutral resonance of mass m coupling c width w
// including interferences with SM
////////////////////////////////////////////////////////////////////
namespace resonances0even_az { extern double sqme(double, double, int, int, double, double, double, double, double); }
void resonances0even_sqme_az_c_(double* _amp2, double* _s, double* _t, int* _exclude_loops_SM, int* _exclude_loops_EX, double* _m, double* _c, double* _czg, double* _w_c, double* _aa) {
  *_amp2 = resonances0even_az::sqme(*_s, *_t, *_exclude_loops_SM, *_exclude_loops_EX, *_m, *_c, *_czg, *_w_c, *_aa);
}
#define resonances0even_sqme_az_c__ resonances0even_sqme_az_c_ //wrapper for g77

//////////////////////////////////////////////////////////////////// 
// BSM aaaa with neutral resonance of mass m coupling c width w
// including interferences with SM
//////////////////////////////////////////////////////////////////// 
namespace resonances2_aaaa { extern double sqme(double, double, int, int, double, double, double, double); }
void resonances2_sqme_aaaa_c_(double* _amp2, double* _s, double* _t, int* _exclude_loops_SM, int* _exclude_loops_EX, double* _m, double* _c, double* _w_c, double* _aa) {
  *_amp2 = resonances2_aaaa::sqme(*_s, *_t, *_exclude_loops_SM, *_exclude_loops_EX, *_m, *_c, *_w_c, *_aa);
}
#define resonances2_sqme_aaaa_c__ resonances2_sqme_aaaa_c_ //wrapper for g77

/////////////////////////////////////////////////////////////////// 
// Anomalous aaaa in the EFT limit parametrized by z1 and z2
// including interferences with SM
//  _cutoff_scale means the scale cutoff - if < 0, no formfactor used
//////////////////////////////////////////////////////////////////// 
namespace eft_aaaa { extern double sqme(double, double, int, double, double); }
void eft_sqme_aaaa_c_(double* _amp2, double* _s, double* _t, int* _exclude_loops_SM, double* _z1, double* _z2, double* _cutoff) {
  // apply cutoff  for anomalous couplings 
  //R.S. Gupta arXiv:1111.3354 [hep-ph] 
  double fact = (*_cutoff > 0 ) ? 1/(1+pow(*_s/ *_cutoff / *_cutoff, 2)) : 1;
  *_amp2 = eft_aaaa::sqme( *_s, *_t, *_exclude_loops_SM, *_z1*fact, *_z2*fact );
}

/////////////////////////////////////////////////////////////////// 
// Anomalous aaaz in the EFT limit parametrized by z1 and z2
// including interferences with SM
//  _cutoff_scale means the scale cutoff - if < 0, no formfactor used
//////////////////////////////////////////////////////////////////// 
namespace eft_aaaz { extern double sqme( double, double, int, double, double ); }
void eft_sqme_aaaz_c_( double* _amp2, double* _s, double* _t, int* _exclude_loops_SM, double* _z1, double* _z2, double* _cutoff ) {
  // apply cutoff  for anomalous couplings 
  //R.S. Gupta arXiv:1111.3354 [hep-ph] 
  double fact = (*_cutoff > 0 ) ? 1/(1+pow(*_s/ *_cutoff / *_cutoff, 2)) : 1;
  *_amp2 = eft_aaaz::sqme( *_s, *_t, *_exclude_loops_SM, *_z1*fact, *_z2*fact );
}

#define eft_sqme_aaaa_c__ eft_sqme_aaaa_c_ //wrapper for g77

