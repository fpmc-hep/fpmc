//==========================================================================
// This file has been automatically generated for C++
// MadGraph5_aMC@NLO v. 2.7.3, 2020-06-21
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef Parameters_AAttbar_UFO_H
#define Parameters_AAttbar_UFO_H

#include <complex> 

#include "read_slha.h"
using namespace std; 

class Parameters_AAttbar_UFO
{
  public:

    static Parameters_AAttbar_UFO * getInstance(); 

    // Define "zero"
    double zero, ZERO; 
    // Model parameters independent of aS
    double mdl_WH, mdl_WW, mdl_WZ, mdl_WT, mdl_ymtau, mdl_ymm, mdl_yme,
        mdl_ymt, mdl_ymb, mdl_ymc, mdl_yms, mdl_ymup, mdl_ymdo, aS, mdl_Gf,
        aEWM1, mdl_MH, mdl_MZ, mdl_MTA, mdl_MMU, mdl_Me, mdl_MT, mdl_MB,
        mdl_MC, mdl_MS, mdl_MU, mdl_MD, mdl_xi6, mdl_xi5, mdl_xi4, mdl_xi3,
        mdl_xi2, mdl_xi1, mdl_cabi, mdl_cos__cabi, mdl_sin__cabi,
        mdl_MZ__exp__2, mdl_MZ__exp__4, mdl_sqrt__2, mdl_MH__exp__2, mdl_aEW,
        mdl_MW, mdl_sqrt__aEW, mdl_ee, mdl_MW__exp__2, mdl_sw2, mdl_cw,
        mdl_sqrt__sw2, mdl_sw, mdl_g1, mdl_gw, mdl_vev, mdl_vev__exp__2,
        mdl_lam, mdl_yb, mdl_yc, mdl_ydo, mdl_ye, mdl_ym, mdl_ys, mdl_yt,
        mdl_ytau, mdl_yup, mdl_muH, mdl_ee__exp__2, mdl_sw__exp__2,
        mdl_cw__exp__2;
    std::complex<double> mdl_CKM1x1, mdl_CKM1x2, mdl_CKM1x3, mdl_CKM2x1,
        mdl_CKM2x2, mdl_CKM2x3, mdl_CKM3x1, mdl_CKM3x2, mdl_CKM3x3,
        mdl_conjg__CKM1x1, mdl_conjg__CKM2x1, mdl_conjg__CKM3x1,
        mdl_conjg__CKM1x2, mdl_conjg__CKM2x2, mdl_conjg__CKM3x2,
        mdl_conjg__CKM1x3, mdl_conjg__CKM2x3, mdl_conjg__CKM3x3, mdl_complexi,
        mdl_I1b11, mdl_I1b12, mdl_I1b13, mdl_I1b21, mdl_I1b22, mdl_I1b23,
        mdl_I1b31, mdl_I1b32, mdl_I1b33, mdl_I2b11, mdl_I2b12, mdl_I2b13,
        mdl_I2b21, mdl_I2b22, mdl_I2b23, mdl_I2b31, mdl_I2b32, mdl_I2b33,
        mdl_I3b11, mdl_I3b12, mdl_I3b13, mdl_I3b21, mdl_I3b22, mdl_I3b23,
        mdl_I3b31, mdl_I3b32, mdl_I3b33, mdl_I4b11, mdl_I4b12, mdl_I4b13,
        mdl_I4b21, mdl_I4b22, mdl_I4b23, mdl_I4b31, mdl_I4b32, mdl_I4b33;
    // Model parameters dependent on aS
    double mdl_sqrt__aS, G, mdl_G__exp__2; 
    // Model couplings independent of aS
    std::complex<double> GC_100, GC_101, GC_102, GC_103, GC_104, GC_106; 
    // Model couplings dependent on aS


    // Set parameters that are unchanged during the run
    void setIndependentParameters(SLHAReader& slha, double* couplings, double m_top); 
    // Set couplings that are unchanged during the run
    void setIndependentCouplings(); 
    // Set parameters that are changed event by event
    void setDependentParameters(); 
    // Set couplings that are changed event by event
    void setDependentCouplings(); 

    // Print parameters that are unchanged during the run
    void printIndependentParameters(); 
    // Print couplings that are unchanged during the run
    void printIndependentCouplings(); 
    // Print parameters that are changed event by event
    void printDependentParameters(); 
    // Print couplings that are changed event by event
    void printDependentCouplings(); 


  private:
    static Parameters_AAttbar_UFO * instance; 
}; 

#endif  // Parameters_AAttbar_UFO_H

