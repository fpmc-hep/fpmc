#include "../MG_standalone_eft_ttbar/SubProcesses/P1_Sigma_AAttbar_UFO_aa_ttx/EFT_AAttbarProcess.h"
#include "herwig.h"
#include <iostream>
#include <math.h>
#include <unistd.h>
#define GetCurrentDir getcwd
using namespace std;

namespace eft_aattbar {
// Computes the squared matrix element
double sqme(double *_s, double *_t,
            double *_xi1, double *_xi2, double *_xi3, double *_xi4,
            double *_xi5, double *_xi6, double *_m_top) {

  // We compute the sqme (Lorentz-invariant) in the reference frame of the ttbar system
  // We choose phi = 0 for the outgoing particles
  // The z-axis in this frame is the same as in the lab, as photons are collinear
  // A simple boost along the z-axis is performed

  // Get Mandelstam variables
  double s = *_s;
  double t = *_t;
  double m_top = *_m_top;
  double u = -t -s + 2*m_top*m_top;

  double sqrts = sqrt(s);
  const double phi = 0;

  // EFT couplings
  double xi1 = *_xi1;
  double xi2 = *_xi2;
  double xi3 = *_xi3;
  double xi4 = *_xi4;
  double xi5 = *_xi5;
  double xi6 = *_xi6;
  double couplings[6] = {xi1, xi2, xi3, xi4, xi5, xi6};

  double beta=sqrt(1-4*m_top*m_top/s);
  double costheta = (t-u)/beta/s;

  // Energy and momentum for outgoing tops
  double e = sqrts/2;
  double p = sqrt(e*e - m_top*m_top);
  // Cartesian components of the momentum
  double pz = p*costheta;
  double pt = p*sqrt(1-costheta*costheta);
  double px = pt*cos(phi);
  double py = pt*sin(phi);

  // Convert momenta to MG desired format
  // MG: (E,px,py,pz)
  double p1[4] = {sqrts/2,0,0,sqrts/2};
  double p2[4] = {sqrts/2,0,0,-sqrts/2};
  double p3[4] = {e,px,py,pz};
  double p4[4] = {e,-px,-py,-pz};
  vector<double *> p_vec = {p1, p2, p3, p4};

  // Create a process object
  EFT_AAttbarProcess process;

  // Read param_card and set parameters
  process.initProc("External/excl_aaaa/MG_standalone_eft_ttbar/Cards/param_card.dat",
                   couplings,m_top);

  // Set momenta for this event
  process.setMomenta(p_vec);

  // Evaluate matrix element
  process.sigmaKin();

  const double *matrix_elements = process.getMatrixElements();

  // cout << "Momenta:" << endl;
  // for(int i=0;i < process.nexternal; i++)
  // 	cout << setw(4) << i+1
  // << setiosflags(ios::scientific) << setw(14) << p[i][0]
  // << setiosflags(ios::scientific) << setw(14) << p[i][1]
  // << setiosflags(ios::scientific) << setw(14) << p[i][2]
  // << setiosflags(ios::scientific) << setw(14) << p[i][3] << endl;
  // cout << "
  // -----------------------------------------------------------------------------"
  // << endl;

  // Display matrix elements
  // for (int i = 0; i < process.nprocesses; i++)
    // cout << " Matrix element = "
         //	 << setiosflags(ios::fixed) << setprecision(17)
         // << matrix_elements[i] << " GeV^" << -(2 * process.nexternal - 8)
         // << endl;

  // cout << " -------------------------------------------------------------------"
  //         "----------"
       // << endl;

  return matrix_elements[0];
}
} // namespace eft_aattbar