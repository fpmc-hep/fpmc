//==========================================================================
// This file has been automatically generated for C++ by
// MadGraph5_aMC@NLO v. 2.7.3, 2020-06-21
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include <iostream> 
#include <iomanip> 
#include "Parameters_AAttbar_UFO.h"

// Initialize static instance
Parameters_AAttbar_UFO * Parameters_AAttbar_UFO::instance = 0; 

// Function to get static instance - only one instance per program
Parameters_AAttbar_UFO * Parameters_AAttbar_UFO::getInstance()
{
  if (instance == 0)
    instance = new Parameters_AAttbar_UFO(); 

  return instance; 
}

void Parameters_AAttbar_UFO::setIndependentParameters(SLHAReader& slha, double* couplings, double m_top)
{
  // Define "zero"
  zero = 0; 
  ZERO = 0; 
  // Prepare a vector for indices
  vector<int> indices(2, 0); 
  mdl_WH = slha.get_block_entry("decay", 25, 4.070000e-03); 
  mdl_WW = slha.get_block_entry("decay", 24, 2.085000e+00); 
  mdl_WZ = slha.get_block_entry("decay", 23, 2.495200e+00); 
  mdl_WT = slha.get_block_entry("decay", 6, 1.508336e+00); 
  mdl_ymtau = slha.get_block_entry("yukawa", 15, 1.777000e+00); 
  mdl_ymm = slha.get_block_entry("yukawa", 13, 1.056600e-01); 
  mdl_yme = slha.get_block_entry("yukawa", 11, 5.110000e-04); 
  mdl_ymt = slha.get_block_entry("yukawa", 6, 1.720000e+02); 
  mdl_ymb = slha.get_block_entry("yukawa", 5, 4.700000e+00); 
  mdl_ymc = slha.get_block_entry("yukawa", 4, 1.270000e+00); 
  mdl_yms = slha.get_block_entry("yukawa", 3, 1.010000e-01); 
  mdl_ymup = slha.get_block_entry("yukawa", 2, 2.550000e-03); 
  mdl_ymdo = slha.get_block_entry("yukawa", 1, 5.040000e-03); 
  aS = slha.get_block_entry("sminputs", 3, 1.184000e-01); 
  mdl_Gf = slha.get_block_entry("sminputs", 2, 1.166370e-05); 
  aEWM1 = slha.get_block_entry("sminputs", 1, 1.279000e+02); 
  mdl_MH = slha.get_block_entry("mass", 25, 1.250000e+02); 
  mdl_MZ = slha.get_block_entry("mass", 23, 9.118760e+01); 
  mdl_MTA = slha.get_block_entry("mass", 15, 1.777000e+00); 
  mdl_MMU = slha.get_block_entry("mass", 13, 1.056600e-01); 
  mdl_Me = slha.get_block_entry("mass", 11, 5.110000e-04); 
  mdl_MT = slha.get_block_entry("mass", 6, 1.720000e+02); 
  mdl_MB = slha.get_block_entry("mass", 5, 4.700000e+00); 
  mdl_MC = slha.get_block_entry("mass", 4, 1.270000e+00); 
  mdl_MS = slha.get_block_entry("mass", 3, 1.010000e-01); 
  mdl_MU = slha.get_block_entry("mass", 2, 2.550000e-03); 
  mdl_MD = slha.get_block_entry("mass", 1, 5.040000e-03); 

  mdl_xi6 = couplings[5]; 
  mdl_xi5 = couplings[4]; 
  mdl_xi4 = couplings[3]; 
  mdl_xi3 = couplings[2]; 
  mdl_xi2 = couplings[1]; 
  mdl_xi1 = couplings[0];

  mdl_cabi = slha.get_block_entry("ckmblock", 1, 2.277360e-01); 
  mdl_cos__cabi = cos(mdl_cabi); 
  mdl_CKM1x1 = mdl_cos__cabi; 
  mdl_sin__cabi = sin(mdl_cabi); 
  mdl_CKM1x2 = mdl_sin__cabi; 
  mdl_CKM1x3 = 0.; 
  mdl_CKM2x1 = -mdl_sin__cabi; 
  mdl_CKM2x2 = mdl_cos__cabi; 
  mdl_CKM2x3 = 0.; 
  mdl_CKM3x1 = 0.; 
  mdl_CKM3x2 = 0.; 
  mdl_CKM3x3 = 1.; 
  mdl_MZ__exp__2 = ((mdl_MZ) * (mdl_MZ)); 
  mdl_MZ__exp__4 = ((mdl_MZ) * (mdl_MZ) * (mdl_MZ) * (mdl_MZ)); 
  mdl_sqrt__2 = sqrt(2.); 
  mdl_MH__exp__2 = ((mdl_MH) * (mdl_MH)); 
  mdl_conjg__CKM1x1 = conj(mdl_CKM1x1); 
  mdl_conjg__CKM2x1 = conj(mdl_CKM2x1); 
  mdl_conjg__CKM3x1 = conj(mdl_CKM3x1); 
  mdl_conjg__CKM1x2 = conj(mdl_CKM1x2); 
  mdl_conjg__CKM2x2 = conj(mdl_CKM2x2); 
  mdl_conjg__CKM3x2 = conj(mdl_CKM3x2); 
  mdl_conjg__CKM1x3 = conj(mdl_CKM1x3); 
  mdl_conjg__CKM2x3 = conj(mdl_CKM2x3); 
  mdl_conjg__CKM3x3 = conj(mdl_CKM3x3); 
  mdl_complexi = std::complex<double> (0., 1.); 
  mdl_aEW = 1./aEWM1; 
  mdl_MW = sqrt(mdl_MZ__exp__2/2. + sqrt(mdl_MZ__exp__4/4. - (mdl_aEW * M_PI *
      mdl_MZ__exp__2)/(mdl_Gf * mdl_sqrt__2)));
  mdl_sqrt__aEW = sqrt(mdl_aEW); 
  mdl_ee = 2. * mdl_sqrt__aEW * sqrt(M_PI); 
  mdl_MW__exp__2 = ((mdl_MW) * (mdl_MW)); 
  mdl_sw2 = 1. - mdl_MW__exp__2/mdl_MZ__exp__2; 
  mdl_cw = sqrt(1. - mdl_sw2); 
  mdl_sqrt__sw2 = sqrt(mdl_sw2); 
  mdl_sw = mdl_sqrt__sw2; 
  mdl_g1 = mdl_ee/mdl_cw; 
  mdl_gw = mdl_ee/mdl_sw; 
  mdl_vev = (2. * mdl_MW * mdl_sw)/mdl_ee; 
  mdl_vev__exp__2 = ((mdl_vev) * (mdl_vev)); 
  mdl_lam = mdl_MH__exp__2/(2. * mdl_vev__exp__2); 
  mdl_yb = (mdl_ymb * mdl_sqrt__2)/mdl_vev; 
  mdl_yc = (mdl_ymc * mdl_sqrt__2)/mdl_vev; 
  mdl_ydo = (mdl_ymdo * mdl_sqrt__2)/mdl_vev; 
  mdl_ye = (mdl_yme * mdl_sqrt__2)/mdl_vev; 
  mdl_ym = (mdl_ymm * mdl_sqrt__2)/mdl_vev; 
  mdl_ys = (mdl_yms * mdl_sqrt__2)/mdl_vev; 
  mdl_yt = (mdl_ymt * mdl_sqrt__2)/mdl_vev; 
  mdl_ytau = (mdl_ymtau * mdl_sqrt__2)/mdl_vev; 
  mdl_yup = (mdl_ymup * mdl_sqrt__2)/mdl_vev; 
  mdl_muH = sqrt(mdl_lam * mdl_vev__exp__2); 
  mdl_I1b11 = mdl_ydo * mdl_conjg__CKM1x1; 
  mdl_I1b12 = mdl_ydo * mdl_conjg__CKM2x1; 
  mdl_I1b13 = mdl_ydo * mdl_conjg__CKM3x1; 
  mdl_I1b21 = mdl_ys * mdl_conjg__CKM1x2; 
  mdl_I1b22 = mdl_ys * mdl_conjg__CKM2x2; 
  mdl_I1b23 = mdl_ys * mdl_conjg__CKM3x2; 
  mdl_I1b31 = mdl_yb * mdl_conjg__CKM1x3; 
  mdl_I1b32 = mdl_yb * mdl_conjg__CKM2x3; 
  mdl_I1b33 = mdl_yb * mdl_conjg__CKM3x3; 
  mdl_I2b11 = mdl_yup * mdl_conjg__CKM1x1; 
  mdl_I2b12 = mdl_yc * mdl_conjg__CKM2x1; 
  mdl_I2b13 = mdl_yt * mdl_conjg__CKM3x1; 
  mdl_I2b21 = mdl_yup * mdl_conjg__CKM1x2; 
  mdl_I2b22 = mdl_yc * mdl_conjg__CKM2x2; 
  mdl_I2b23 = mdl_yt * mdl_conjg__CKM3x2; 
  mdl_I2b31 = mdl_yup * mdl_conjg__CKM1x3; 
  mdl_I2b32 = mdl_yc * mdl_conjg__CKM2x3; 
  mdl_I2b33 = mdl_yt * mdl_conjg__CKM3x3; 
  mdl_I3b11 = mdl_CKM1x1 * mdl_yup; 
  mdl_I3b12 = mdl_CKM1x2 * mdl_yup; 
  mdl_I3b13 = mdl_CKM1x3 * mdl_yup; 
  mdl_I3b21 = mdl_CKM2x1 * mdl_yc; 
  mdl_I3b22 = mdl_CKM2x2 * mdl_yc; 
  mdl_I3b23 = mdl_CKM2x3 * mdl_yc; 
  mdl_I3b31 = mdl_CKM3x1 * mdl_yt; 
  mdl_I3b32 = mdl_CKM3x2 * mdl_yt; 
  mdl_I3b33 = mdl_CKM3x3 * mdl_yt; 
  mdl_I4b11 = mdl_CKM1x1 * mdl_ydo; 
  mdl_I4b12 = mdl_CKM1x2 * mdl_ys; 
  mdl_I4b13 = mdl_CKM1x3 * mdl_yb; 
  mdl_I4b21 = mdl_CKM2x1 * mdl_ydo; 
  mdl_I4b22 = mdl_CKM2x2 * mdl_ys; 
  mdl_I4b23 = mdl_CKM2x3 * mdl_yb; 
  mdl_I4b31 = mdl_CKM3x1 * mdl_ydo; 
  mdl_I4b32 = mdl_CKM3x2 * mdl_ys; 
  mdl_I4b33 = mdl_CKM3x3 * mdl_yb; 
  mdl_ee__exp__2 = ((mdl_ee) * (mdl_ee)); 
  mdl_sw__exp__2 = ((mdl_sw) * (mdl_sw)); 
  mdl_cw__exp__2 = ((mdl_cw) * (mdl_cw)); 
}
void Parameters_AAttbar_UFO::setIndependentCouplings()
{
  GC_100 = 4. * mdl_complexi * mdl_MT * mdl_xi1; 
  GC_101 = -2. * mdl_MT * mdl_xi2; 
  GC_102 = -2. * mdl_complexi * mdl_MT * mdl_xi3; 
  GC_103 = 4. * mdl_MT * mdl_xi4; 
  GC_104 = mdl_complexi * mdl_xi5; 
  GC_106 = mdl_complexi * mdl_xi6; 
}
void Parameters_AAttbar_UFO::setDependentParameters()
{
  mdl_sqrt__aS = sqrt(aS); 
  G = 2. * mdl_sqrt__aS * sqrt(M_PI); 
  mdl_G__exp__2 = ((G) * (G)); 
}
void Parameters_AAttbar_UFO::setDependentCouplings()
{

}

// Routines for printing out parameters
void Parameters_AAttbar_UFO::printIndependentParameters()
{
  cout <<  "AAttbar_UFO model parameters independent of event kinematics:" <<
      endl;
  cout << setw(20) <<  "mdl_WH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WH << endl;
  cout << setw(20) <<  "mdl_WW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WW << endl;
  cout << setw(20) <<  "mdl_WZ " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WZ << endl;
  cout << setw(20) <<  "mdl_WT " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WT << endl;
  cout << setw(20) <<  "mdl_ymtau " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ymtau << endl;
  cout << setw(20) <<  "mdl_ymm " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ymm << endl;
  cout << setw(20) <<  "mdl_yme " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_yme << endl;
  cout << setw(20) <<  "mdl_ymt " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ymt << endl;
  cout << setw(20) <<  "mdl_ymb " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ymb << endl;
  cout << setw(20) <<  "mdl_ymc " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ymc << endl;
  cout << setw(20) <<  "mdl_yms " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_yms << endl;
  cout << setw(20) <<  "mdl_ymup " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ymup << endl;
  cout << setw(20) <<  "mdl_ymdo " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ymdo << endl;
  cout << setw(20) <<  "aS " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << aS << endl;
  cout << setw(20) <<  "mdl_Gf " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_Gf << endl;
  cout << setw(20) <<  "aEWM1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << aEWM1 << endl;
  cout << setw(20) <<  "mdl_MH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MH << endl;
  cout << setw(20) <<  "mdl_MZ " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MZ << endl;
  cout << setw(20) <<  "mdl_MTA " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MTA << endl;
  cout << setw(20) <<  "mdl_MMU " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MMU << endl;
  cout << setw(20) <<  "mdl_Me " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_Me << endl;
  cout << setw(20) <<  "mdl_MT " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MT << endl;
  cout << setw(20) <<  "mdl_MB " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MB << endl;
  cout << setw(20) <<  "mdl_MC " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MC << endl;
  cout << setw(20) <<  "mdl_MS " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MS << endl;
  cout << setw(20) <<  "mdl_MU " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MU << endl;
  cout << setw(20) <<  "mdl_MD " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MD << endl;
  cout << setw(20) <<  "mdl_xi6 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_xi6 << endl;
  cout << setw(20) <<  "mdl_xi5 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_xi5 << endl;
  cout << setw(20) <<  "mdl_xi4 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_xi4 << endl;
  cout << setw(20) <<  "mdl_xi3 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_xi3 << endl;
  cout << setw(20) <<  "mdl_xi2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_xi2 << endl;
  cout << setw(20) <<  "mdl_xi1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_xi1 << endl;
  cout << setw(20) <<  "mdl_cabi " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_cabi << endl;
  cout << setw(20) <<  "mdl_cos__cabi " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_cos__cabi << endl;
  cout << setw(20) <<  "mdl_CKM1x1 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_CKM1x1 << endl;
  cout << setw(20) <<  "mdl_sin__cabi " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_sin__cabi << endl;
  cout << setw(20) <<  "mdl_CKM1x2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_CKM1x2 << endl;
  cout << setw(20) <<  "mdl_CKM1x3 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_CKM1x3 << endl;
  cout << setw(20) <<  "mdl_CKM2x1 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_CKM2x1 << endl;
  cout << setw(20) <<  "mdl_CKM2x2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_CKM2x2 << endl;
  cout << setw(20) <<  "mdl_CKM2x3 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_CKM2x3 << endl;
  cout << setw(20) <<  "mdl_CKM3x1 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_CKM3x1 << endl;
  cout << setw(20) <<  "mdl_CKM3x2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_CKM3x2 << endl;
  cout << setw(20) <<  "mdl_CKM3x3 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_CKM3x3 << endl;
  cout << setw(20) <<  "mdl_MZ__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MZ__exp__2 << endl;
  cout << setw(20) <<  "mdl_MZ__exp__4 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MZ__exp__4 << endl;
  cout << setw(20) <<  "mdl_sqrt__2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_sqrt__2 << endl;
  cout << setw(20) <<  "mdl_MH__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MH__exp__2 << endl;
  cout << setw(20) <<  "mdl_conjg__CKM1x1 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_conjg__CKM1x1 << endl;
  cout << setw(20) <<  "mdl_conjg__CKM2x1 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_conjg__CKM2x1 << endl;
  cout << setw(20) <<  "mdl_conjg__CKM3x1 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_conjg__CKM3x1 << endl;
  cout << setw(20) <<  "mdl_conjg__CKM1x2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_conjg__CKM1x2 << endl;
  cout << setw(20) <<  "mdl_conjg__CKM2x2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_conjg__CKM2x2 << endl;
  cout << setw(20) <<  "mdl_conjg__CKM3x2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_conjg__CKM3x2 << endl;
  cout << setw(20) <<  "mdl_conjg__CKM1x3 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_conjg__CKM1x3 << endl;
  cout << setw(20) <<  "mdl_conjg__CKM2x3 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_conjg__CKM2x3 << endl;
  cout << setw(20) <<  "mdl_conjg__CKM3x3 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_conjg__CKM3x3 << endl;
  cout << setw(20) <<  "mdl_complexi " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_complexi << endl;
  cout << setw(20) <<  "mdl_aEW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_aEW << endl;
  cout << setw(20) <<  "mdl_MW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MW << endl;
  cout << setw(20) <<  "mdl_sqrt__aEW " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_sqrt__aEW << endl;
  cout << setw(20) <<  "mdl_ee " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ee << endl;
  cout << setw(20) <<  "mdl_MW__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MW__exp__2 << endl;
  cout << setw(20) <<  "mdl_sw2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_sw2 << endl;
  cout << setw(20) <<  "mdl_cw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_cw << endl;
  cout << setw(20) <<  "mdl_sqrt__sw2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_sqrt__sw2 << endl;
  cout << setw(20) <<  "mdl_sw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_sw << endl;
  cout << setw(20) <<  "mdl_g1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_g1 << endl;
  cout << setw(20) <<  "mdl_gw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_gw << endl;
  cout << setw(20) <<  "mdl_vev " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_vev << endl;
  cout << setw(20) <<  "mdl_vev__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_vev__exp__2 << endl;
  cout << setw(20) <<  "mdl_lam " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_lam << endl;
  cout << setw(20) <<  "mdl_yb " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_yb << endl;
  cout << setw(20) <<  "mdl_yc " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_yc << endl;
  cout << setw(20) <<  "mdl_ydo " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ydo << endl;
  cout << setw(20) <<  "mdl_ye " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ye << endl;
  cout << setw(20) <<  "mdl_ym " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ym << endl;
  cout << setw(20) <<  "mdl_ys " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ys << endl;
  cout << setw(20) <<  "mdl_yt " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_yt << endl;
  cout << setw(20) <<  "mdl_ytau " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ytau << endl;
  cout << setw(20) <<  "mdl_yup " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_yup << endl;
  cout << setw(20) <<  "mdl_muH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_muH << endl;
  cout << setw(20) <<  "mdl_I1b11 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I1b11 << endl;
  cout << setw(20) <<  "mdl_I1b12 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I1b12 << endl;
  cout << setw(20) <<  "mdl_I1b13 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I1b13 << endl;
  cout << setw(20) <<  "mdl_I1b21 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I1b21 << endl;
  cout << setw(20) <<  "mdl_I1b22 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I1b22 << endl;
  cout << setw(20) <<  "mdl_I1b23 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I1b23 << endl;
  cout << setw(20) <<  "mdl_I1b31 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I1b31 << endl;
  cout << setw(20) <<  "mdl_I1b32 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I1b32 << endl;
  cout << setw(20) <<  "mdl_I1b33 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I1b33 << endl;
  cout << setw(20) <<  "mdl_I2b11 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I2b11 << endl;
  cout << setw(20) <<  "mdl_I2b12 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I2b12 << endl;
  cout << setw(20) <<  "mdl_I2b13 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I2b13 << endl;
  cout << setw(20) <<  "mdl_I2b21 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I2b21 << endl;
  cout << setw(20) <<  "mdl_I2b22 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I2b22 << endl;
  cout << setw(20) <<  "mdl_I2b23 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I2b23 << endl;
  cout << setw(20) <<  "mdl_I2b31 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I2b31 << endl;
  cout << setw(20) <<  "mdl_I2b32 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I2b32 << endl;
  cout << setw(20) <<  "mdl_I2b33 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I2b33 << endl;
  cout << setw(20) <<  "mdl_I3b11 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I3b11 << endl;
  cout << setw(20) <<  "mdl_I3b12 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I3b12 << endl;
  cout << setw(20) <<  "mdl_I3b13 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I3b13 << endl;
  cout << setw(20) <<  "mdl_I3b21 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I3b21 << endl;
  cout << setw(20) <<  "mdl_I3b22 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I3b22 << endl;
  cout << setw(20) <<  "mdl_I3b23 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I3b23 << endl;
  cout << setw(20) <<  "mdl_I3b31 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I3b31 << endl;
  cout << setw(20) <<  "mdl_I3b32 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I3b32 << endl;
  cout << setw(20) <<  "mdl_I3b33 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I3b33 << endl;
  cout << setw(20) <<  "mdl_I4b11 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I4b11 << endl;
  cout << setw(20) <<  "mdl_I4b12 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I4b12 << endl;
  cout << setw(20) <<  "mdl_I4b13 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I4b13 << endl;
  cout << setw(20) <<  "mdl_I4b21 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I4b21 << endl;
  cout << setw(20) <<  "mdl_I4b22 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I4b22 << endl;
  cout << setw(20) <<  "mdl_I4b23 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I4b23 << endl;
  cout << setw(20) <<  "mdl_I4b31 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I4b31 << endl;
  cout << setw(20) <<  "mdl_I4b32 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I4b32 << endl;
  cout << setw(20) <<  "mdl_I4b33 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I4b33 << endl;
  cout << setw(20) <<  "mdl_ee__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_ee__exp__2 << endl;
  cout << setw(20) <<  "mdl_sw__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_sw__exp__2 << endl;
  cout << setw(20) <<  "mdl_cw__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_cw__exp__2 << endl;
}
void Parameters_AAttbar_UFO::printIndependentCouplings()
{
  cout <<  "AAttbar_UFO model couplings independent of event kinematics:" <<
      endl;
  cout << setw(20) <<  "GC_100 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_100 << endl;
  cout << setw(20) <<  "GC_101 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_101 << endl;
  cout << setw(20) <<  "GC_102 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_102 << endl;
  cout << setw(20) <<  "GC_103 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_103 << endl;
  cout << setw(20) <<  "GC_104 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_104 << endl;
  cout << setw(20) <<  "GC_106 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_106 << endl;
}
void Parameters_AAttbar_UFO::printDependentParameters()
{
  cout <<  "AAttbar_UFO model parameters dependent on event kinematics:" <<
      endl;
  cout << setw(20) <<  "mdl_sqrt__aS " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_sqrt__aS << endl;
  cout << setw(20) <<  "G " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << G << endl;
  cout << setw(20) <<  "mdl_G__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_G__exp__2 << endl;
}
void Parameters_AAttbar_UFO::printDependentCouplings()
{
  cout <<  "AAttbar_UFO model couplings dependent on event kinematics:" <<
      endl;

}


