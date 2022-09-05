//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.7.3, 2020-06-21
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "HelAmps_AAttbar_UFO.h"
#include <complex> 
#include <cmath> 
#include <iostream> 
#include <cstdlib> 
using namespace std; 

namespace MG5_AAttbar_UFO 
{

void ixxxxx(double p[4], double fmass, int nhel, int nsf, complex<double> fi[6])
{
  complex<double> chi[2]; 
  double sf[2], sfomega[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int ip, im, nh; 
  fi[0] = complex<double> (-p[0] * nsf, -p[3] * nsf); 
  fi[1] = complex<double> (-p[1] * nsf, -p[2] * nsf); 
  nh = nhel * nsf; 
  if (fmass != 0.0)
  {
    pp = min(p[0], sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3])); 
    if (pp == 0.0)
    {
      sqm[0] = sqrt(std::abs(fmass)); 
      sqm[1] = Sgn(sqm[0], fmass); 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      fi[2] = ip * sqm[ip]; 
      fi[3] = im * nsf * sqm[ip]; 
      fi[4] = ip * nsf * sqm[im]; 
      fi[5] = im * sqm[im]; 
    }
    else
    {
      sf[0] = (1 + nsf + (1 - nsf) * nh) * 0.5; 
      sf[1] = (1 + nsf - (1 - nsf) * nh) * 0.5; 
      omega[0] = sqrt(p[0] + pp); 
      omega[1] = fmass/omega[0]; 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      sfomega[0] = sf[0] * omega[ip]; 
      sfomega[1] = sf[1] * omega[im]; 
      pp3 = max(pp + p[3], 0.0); 
      chi[0] = complex<double> (sqrt(pp3 * 0.5/pp), 0); 
      if (pp3 == 0.0)
      {
        chi[1] = complex<double> (-nh, 0); 
      }
      else
      {
        chi[1] = complex<double> (nh * p[1], p[2])/sqrt(2.0 * pp * pp3); 
      }
      fi[2] = sfomega[0] * chi[im]; 
      fi[3] = sfomega[0] * chi[ip]; 
      fi[4] = sfomega[1] * chi[im]; 
      fi[5] = sfomega[1] * chi[ip]; 
    }
  }
  else
  {
    if (p[1] == 0.0 and p[2] == 0.0 and p[3] < 0.0)
    {
      sqp0p3 = 0.0; 
    }
    else
    {
      sqp0p3 = sqrt(max(p[0] + p[3], 0.0)) * nsf; 
    }
    chi[0] = complex<double> (sqp0p3, 0.0); 
    if (sqp0p3 == 0.0)
    {
      chi[1] = complex<double> (-nhel * sqrt(2.0 * p[0]), 0.0); 
    }
    else
    {
      chi[1] = complex<double> (nh * p[1], p[2])/sqp0p3; 
    }
    if (nh == 1)
    {
      fi[2] = complex<double> (0.0, 0.0); 
      fi[3] = complex<double> (0.0, 0.0); 
      fi[4] = chi[0]; 
      fi[5] = chi[1]; 
    }
    else
    {
      fi[2] = chi[1]; 
      fi[3] = chi[0]; 
      fi[4] = complex<double> (0.0, 0.0); 
      fi[5] = complex<double> (0.0, 0.0); 
    }
  }
  return; 
}


double Sgn(double a, double b)
{
  return (b < 0)? - abs(a):abs(a); 
}


void txxxxx(double p[4], double tmass, int nhel, int nst, complex<double>
    tc[18])
{
  complex<double> ft[6][4], ep[4], em[4], e0[4]; 
  double pt, pt2, pp, pzpt, emp, sqh, sqs; 
  int i, j; 

  sqh = sqrt(0.5); 
  sqs = sqrt(0.5/3); 

  pt2 = p[1] * p[1] + p[2] * p[2]; 
  pp = min(p[0], sqrt(pt2 + p[3] * p[3])); 
  pt = min(pp, sqrt(pt2)); 

  ft[4][0] = complex<double> (p[0] * nst, p[3] * nst); 
  ft[5][0] = complex<double> (p[1] * nst, p[2] * nst); 

  // construct eps+
  if(nhel >= 0)
  {
    if(pp == 0)
    {
      ep[0] = complex<double> (0, 0); 
      ep[1] = complex<double> (-sqh, 0); 
      ep[2] = complex<double> (0, nst * sqh); 
      ep[3] = complex<double> (0, 0); 
    }
    else
    {
      ep[0] = complex<double> (0, 0); 
      ep[3] = complex<double> (pt/pp * sqh, 0); 

      if(pt != 0)
      {
        pzpt = p[3]/(pp * pt) * sqh; 
        ep[1] = complex<double> (-p[1] * pzpt, -nst * p[2]/pt * sqh); 
        ep[2] = complex<double> (-p[2] * pzpt, nst * p[1]/pt * sqh); 
      }
      else
      {
        ep[1] = complex<double> (-sqh, 0); 
        ep[2] = complex<double> (0, nst * Sgn(sqh, p[3])); 
      }
    }

  }

  // construct eps-
  if(nhel <= 0)
  {
    if(pp == 0)
    {
      em[0] = complex<double> (0, 0); 
      em[1] = complex<double> (sqh, 0); 
      em[2] = complex<double> (0, nst * sqh); 
      em[3] = complex<double> (0, 0); 
    }
    else
    {
      em[0] = complex<double> (0, 0); 
      em[3] = complex<double> (-pt/pp * sqh, 0); 

      if(pt != 0)
      {
        pzpt = -p[3]/(pp * pt) * sqh; 
        em[1] = complex<double> (-p[1] * pzpt, -nst * p[2]/pt * sqh); 
        em[2] = complex<double> (-p[2] * pzpt, nst * p[1]/pt * sqh); 
      }
      else
      {
        em[1] = complex<double> (sqh, 0); 
        em[2] = complex<double> (0, nst * Sgn(sqh, p[3])); 
      }
    }
  }

  // construct eps0
  if(std::labs(nhel) <= 1)
  {
    if(pp == 0)
    {
      e0[0] = complex<double> (0, 0); 
      e0[1] = complex<double> (0, 0); 
      e0[2] = complex<double> (0, 0); 
      e0[3] = complex<double> (1, 0); 
    }
    else
    {
      emp = p[0]/(tmass * pp); 
      e0[0] = complex<double> (pp/tmass, 0); 
      e0[3] = complex<double> (p[3] * emp, 0); 

      if(pt != 0)
      {
        e0[1] = complex<double> (p[1] * emp, 0); 
        e0[2] = complex<double> (p[2] * emp, 0); 
      }
      else
      {
        e0[1] = complex<double> (0, 0); 
        e0[2] = complex<double> (0, 0); 
      }
    }
  }

  if(nhel == 2)
  {
    for(j = 0; j < 4; j++ )
    {
      for(i = 0; i < 4; i++ )
        ft[i][j] = ep[i] * ep[j]; 
    }
  }
  else if(nhel == -2)
  {
    for(j = 0; j < 4; j++ )
    {
      for(i = 0; i < 4; i++ )
        ft[i][j] = em[i] * em[j]; 
    }
  }
  else if(tmass == 0)
  {
    for(j = 0; j < 4; j++ )
    {
      for(i = 0; i < 4; i++ )
        ft[i][j] = 0; 
    }
  }
  else if(tmass != 0)
  {
    if(nhel == 1)
    {
      for(j = 0; j < 4; j++ )
      {
        for(i = 0; i < 4; i++ )
          ft[i][j] = sqh * (ep[i] * e0[j] + e0[i] * ep[j]); 
      }
    }
    else if(nhel == 0)
    {
      for(j = 0; j < 4; j++ )
      {
        for(i = 0; i < 4; i++ )
          ft[i][j] = sqs * (ep[i] * em[j] + em[i] * ep[j]
         + 2.0 * e0[i] * e0[j]); 
      }
    }
    else if(nhel == -1)
    {
      for(j = 0; j < 4; j++ )
      {
        for(i = 0; i < 4; i++ )
          ft[i][j] = sqh * (em[i] * e0[j] + e0[i] * em[j]); 
      }
    }
    else
    {
      std::cerr <<  "Invalid helicity in txxxxx.\n"; 
      std::exit(1); 
    }
  }

  tc[0] = ft[4][0]; 
  tc[1] = ft[5][0]; 

  for(j = 0; j < 4; j++ )
  {
    for(i = 0; i < 4; i++ )
      tc[j * 4 + i + 2] = ft[j][i]; 
  }
}

void vxxxxx(double p[4], double vmass, int nhel, int nsv, complex<double> vc[6])
{
  double hel, hel0, pt, pt2, pp, pzpt, emp, sqh; 
  int nsvahl; 
  sqh = sqrt(0.5); 
  hel = double(nhel); 
  nsvahl = nsv * std::abs(hel); 
  pt2 = (p[1] * p[1]) + (p[2] * p[2]); 
  pp = min(p[0], sqrt(pt2 + (p[3] * p[3]))); 
  pt = min(pp, sqrt(pt2)); 
  vc[0] = complex<double> (p[0] * nsv, p[3] * nsv); 
  vc[1] = complex<double> (p[1] * nsv, p[2] * nsv); 
  if (vmass != 0.0)
  {
    hel0 = 1.0 - std::abs(hel); 
    if(pp == 0.0)
    {
      vc[2] = complex<double> (0.0, 0.0); 
      vc[3] = complex<double> (-hel * sqh, 0.0); 
      vc[4] = complex<double> (0.0, nsvahl * sqh); 
      vc[5] = complex<double> (hel0, 0.0); 
    }
    else
    {
      emp = p[0]/(vmass * pp); 
      vc[2] = complex<double> (hel0 * pp/vmass, 0.0); 
      vc[5] = complex<double> (hel0 * p[3] * emp + hel * pt/pp * sqh, 0.0); 
      if (pt != 0.0)
      {
        pzpt = p[3]/(pp * pt) * sqh * hel; 
        vc[3] = complex<double> (hel0 * p[1] * emp - p[1] * pzpt, -nsvahl *
            p[2]/pt * sqh);
        vc[4] = complex<double> (hel0 * p[2] * emp - p[2] * pzpt, nsvahl *
            p[1]/pt * sqh);
      }
      else
      {
        vc[3] = complex<double> (-hel * sqh, 0.0); 
        vc[4] = complex<double> (0.0, nsvahl * Sgn(sqh, p[3])); 
      }
    }
  }
  else
  {
    pp = p[0]; 
    pt = sqrt((p[1] * p[1]) + (p[2] * p[2])); 
    vc[2] = complex<double> (0.0, 0.0); 
    vc[5] = complex<double> (hel * pt/pp * sqh, 0.0); 
    if (pt != 0.0)
    {
      pzpt = p[3]/(pp * pt) * sqh * hel; 
      vc[3] = complex<double> (-p[1] * pzpt, -nsv * p[2]/pt * sqh); 
      vc[4] = complex<double> (-p[2] * pzpt, nsv * p[1]/pt * sqh); 
    }
    else
    {
      vc[3] = complex<double> (-hel * sqh, 0.0); 
      vc[4] = complex<double> (0.0, nsv * Sgn(sqh, p[3])); 
    }
  }
  return; 
}

void sxxxxx(double p[4], int nss, complex<double> sc[3])
{
  sc[2] = complex<double> (1.00, 0.00); 
  sc[0] = complex<double> (p[0] * nss, p[3] * nss); 
  sc[1] = complex<double> (p[1] * nss, p[2] * nss); 
  return; 
}

void oxxxxx(double p[4], double fmass, int nhel, int nsf, complex<double> fo[6])
{
  complex<double> chi[2]; 
  double sf[2], sfomeg[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int nh, ip, im; 
  fo[0] = complex<double> (p[0] * nsf, p[3] * nsf); 
  fo[1] = complex<double> (p[1] * nsf, p[2] * nsf); 
  nh = nhel * nsf; 
  if (fmass != 0.000)
  {
    pp = min(p[0], sqrt((p[1] * p[1]) + (p[2] * p[2]) + (p[3] * p[3]))); 
    if (pp == 0.000)
    {
      sqm[0] = sqrt(std::abs(fmass)); 
      sqm[1] = Sgn(sqm[0], fmass); 
      ip = -((1 - nh)/2) * nhel; 
      im = (1 + nh)/2 * nhel; 
      fo[2] = im * sqm[std::abs(ip)]; 
      fo[3] = ip * nsf * sqm[std::abs(ip)]; 
      fo[4] = im * nsf * sqm[std::abs(im)]; 
      fo[5] = ip * sqm[std::abs(im)]; 
    }
    else
    {
      pp = min(p[0], sqrt((p[1] * p[1]) + (p[2] * p[2]) + (p[3] * p[3]))); 
      sf[0] = double(1 + nsf + (1 - nsf) * nh) * 0.5; 
      sf[1] = double(1 + nsf - (1 - nsf) * nh) * 0.5; 
      omega[0] = sqrt(p[0] + pp); 
      omega[1] = fmass/omega[0]; 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      sfomeg[0] = sf[0] * omega[ip]; 
      sfomeg[1] = sf[1] * omega[im]; 
      pp3 = max(pp + p[3], 0.00); 
      chi[0] = complex<double> (sqrt(pp3 * 0.5/pp), 0.00); 
      if (pp3 == 0.00)
      {
        chi[1] = complex<double> (-nh, 0.00); 
      }
      else
      {
        chi[1] = complex<double> (nh * p[1], -p[2])/sqrt(2.0 * pp * pp3); 
      }
      fo[2] = sfomeg[1] * chi[im]; 
      fo[3] = sfomeg[1] * chi[ip]; 
      fo[4] = sfomeg[0] * chi[im]; 
      fo[5] = sfomeg[0] * chi[ip]; 
    }
  }
  else
  {
    if((p[1] == 0.00) and (p[2] == 0.00) and (p[3] < 0.00))
    {
      sqp0p3 = 0.00; 
    }
    else
    {
      sqp0p3 = sqrt(max(p[0] + p[3], 0.00)) * nsf; 
    }
    chi[0] = complex<double> (sqp0p3, 0.00); 
    if(sqp0p3 == 0.000)
    {
      chi[1] = complex<double> (-nhel, 0.00) * sqrt(2.0 * p[0]); 
    }
    else
    {
      chi[1] = complex<double> (nh * p[1], -p[2])/sqp0p3; 
    }
    if(nh == 1)
    {
      fo[2] = chi[0]; 
      fo[3] = chi[1]; 
      fo[4] = complex<double> (0.00, 0.00); 
      fo[5] = complex<double> (0.00, 0.00); 
    }
    else
    {
      fo[2] = complex<double> (0.00, 0.00); 
      fo[3] = complex<double> (0.00, 0.00); 
      fo[4] = chi[1]; 
      fo[5] = chi[0]; 
    }
  }
  return; 
}

void FFVV8_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> V4[], std::complex<double>
    COUP, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP5; 
  std::complex<double> TMP12; 
  std::complex<double> TMP1; 
  double P1[4]; 
  std::complex<double> TMP0; 
  std::complex<double> TMP7; 
  double P3[4]; 
  std::complex<double> TMP10; 
  std::complex<double> TMP6; 
  std::complex<double> TMP4; 
  double P4[4]; 
  std::complex<double> TMP11; 
  std::complex<double> TMP9; 
  std::complex<double> TMP3; 
  std::complex<double> TMP8; 
  P1[0] = F1[0].real(); 
  P1[1] = F1[1].real(); 
  P1[2] = F1[1].imag(); 
  P1[3] = F1[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  P4[0] = V4[0].real(); 
  P4[1] = V4[1].real(); 
  P4[2] = V4[1].imag(); 
  P4[3] = V4[0].imag(); 
  TMP1 = (V3[2] * P4[0] - V3[3] * P4[1] - V3[4] * P4[2] - V3[5] * P4[3]); 
  TMP0 = (V4[2] * P3[0] - V4[3] * P3[1] - V4[4] * P3[2] - V4[5] * P3[3]); 
  TMP9 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      (F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])) +
      (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
      F1[5] * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] + V3[5])))));
  TMP8 = (F1[2] * (F2[4] * (P4[0] + P4[3]) + F2[5] * (P4[1] + cI * (P4[2]))) +
      (F1[3] * (F2[4] * (P4[1] - cI * (P4[2])) + F2[5] * (P4[0] - P4[3])) +
      (F1[4] * (F2[2] * (P4[0] - P4[3]) - F2[3] * (P4[1] + cI * (P4[2]))) +
      F1[5] * (F2[2] * (+cI * (P4[2]) - P4[1]) + F2[3] * (P4[0] + P4[3])))));
  TMP5 = (V4[2] * P1[0] - V4[3] * P1[1] - V4[4] * P1[2] - V4[5] * P1[3]); 
  TMP4 = (P3[0] * P4[0] - P3[1] * P4[1] - P3[2] * P4[2] - P3[3] * P4[3]); 
  TMP7 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP6 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      (F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])) +
      (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
      F1[5] * (F2[2] * (+cI * (P3[2]) - P3[1]) + F2[3] * (P3[0] + P3[3])))));
  TMP11 = (F1[2] * (F2[4] * (V4[2] + V4[5]) + F2[5] * (V4[3] + cI * (V4[4]))) +
      (F1[3] * (F2[4] * (V4[3] - cI * (V4[4])) + F2[5] * (V4[2] - V4[5])) +
      (F1[4] * (F2[2] * (V4[2] - V4[5]) - F2[3] * (V4[3] + cI * (V4[4]))) +
      F1[5] * (F2[2] * (+cI * (V4[4]) - V4[3]) + F2[3] * (V4[2] + V4[5])))));
  TMP10 = (P4[0] * P1[0] - P4[1] * P1[1] - P4[2] * P1[2] - P4[3] * P1[3]); 
  TMP3 = (V4[2] * V3[2] - V4[3] * V3[3] - V4[4] * V3[4] - V4[5] * V3[5]); 
  TMP12 = (P3[0] * P1[0] - P3[1] * P1[1] - P3[2] * P1[2] - P3[3] * P1[3]); 
  vertex = COUP * (TMP0 * (-1.) * (+cI * (TMP7 * TMP8 + TMP9 * TMP10)) + (TMP1
      * (-1.) * (+cI * (TMP5 * TMP6 + TMP11 * TMP12)) + (TMP3 * (+cI * (TMP6 *
      TMP10 + TMP8 * TMP12)) + TMP4 * (+cI * (TMP5 * TMP9 + TMP7 * TMP11)))));
}


void FFVV11_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> V4[], std::complex<double>
    COUP, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP15; 
  std::complex<double> TMP1; 
  std::complex<double> TMP0; 
  double P3[4]; 
  std::complex<double> TMP16; 
  std::complex<double> TMP4; 
  double P4[4]; 
  std::complex<double> TMP3; 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  P4[0] = V4[0].real(); 
  P4[1] = V4[1].real(); 
  P4[2] = V4[1].imag(); 
  P4[3] = V4[0].imag(); 
  TMP15 = (F2[2] * F1[2] + F2[3] * F1[3]); 
  TMP4 = (P3[0] * P4[0] - P3[1] * P4[1] - P3[2] * P4[2] - P3[3] * P4[3]); 
  TMP16 = (F2[4] * F1[4] + F2[5] * F1[5]); 
  TMP1 = (V3[2] * P4[0] - V3[3] * P4[1] - V3[4] * P4[2] - V3[5] * P4[3]); 
  TMP0 = (V4[2] * P3[0] - V4[3] * P3[1] - V4[4] * P3[2] - V4[5] * P3[3]); 
  TMP3 = (V4[2] * V3[2] - V4[3] * V3[3] - V4[4] * V3[4] - V4[5] * V3[5]); 
  vertex = COUP * (TMP0 * TMP1 * (-cI * (TMP15) + cI * (TMP16)) + TMP3 * TMP4 *
      (-cI * (TMP16) + cI * (TMP15)));
}


void FFVV7_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> V4[], std::complex<double>
    COUP, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP2; 
  double P3[4]; 
  std::complex<double> TMP14; 
  double P4[4]; 
  std::complex<double> TMP13; 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  P4[0] = V4[0].real(); 
  P4[1] = V4[1].real(); 
  P4[2] = V4[1].imag(); 
  P4[3] = V4[0].imag(); 
  TMP14 = (-1.) * (P3[0] * (P4[1] * (V4[5] * V3[4] - V4[4] * V3[5]) + (P4[2] *
      (V4[3] * V3[5] - V4[5] * V3[3]) + P4[3] * (V4[4] * V3[3] - V4[3] *
      V3[4]))) + (P3[1] * (P4[0] * (V4[4] * V3[5] - V4[5] * V3[4]) + (P4[2] *
      (V4[5] * V3[2] - V4[2] * V3[5]) + P4[3] * (V4[2] * V3[4] - V4[4] *
      V3[2]))) + (P3[2] * (P4[0] * (V4[5] * V3[3] - V4[3] * V3[5]) + (P4[1] *
      (V4[2] * V3[5] - V4[5] * V3[2]) + P4[3] * (V4[3] * V3[2] - V4[2] *
      V3[3]))) + P3[3] * (P4[0] * (V4[3] * V3[4] - V4[4] * V3[3]) + (P4[1] *
      (V4[4] * V3[2] - V4[2] * V3[4]) + P4[2] * (V4[2] * V3[3] - V4[3] *
      V3[2]))))));
  TMP13 = (-1.) * (P3[0] * (P4[1] * (V4[4] * V3[5] - V4[5] * V3[4]) + (P4[2] *
      (V4[5] * V3[3] - V4[3] * V3[5]) + P4[3] * (V4[3] * V3[4] - V4[4] *
      V3[3]))) + (P3[1] * (P4[0] * (V4[5] * V3[4] - V4[4] * V3[5]) + (P4[2] *
      (V4[2] * V3[5] - V4[5] * V3[2]) + P4[3] * (V4[4] * V3[2] - V4[2] *
      V3[4]))) + (P3[2] * (P4[0] * (V4[3] * V3[5] - V4[5] * V3[3]) + (P4[1] *
      (V4[5] * V3[2] - V4[2] * V3[5]) + P4[3] * (V4[2] * V3[3] - V4[3] *
      V3[2]))) + P3[3] * (P4[0] * (V4[4] * V3[3] - V4[3] * V3[4]) + (P4[1] *
      (V4[2] * V3[4] - V4[4] * V3[2]) + P4[2] * (V4[3] * V3[2] - V4[2] *
      V3[3]))))));
  TMP2 = (F2[2] * F1[2] + F2[3] * F1[3] + F2[4] * F1[4] + F2[5] * F1[5]); 
  vertex = COUP * TMP2 * (-cI * (TMP14) + cI * (TMP13)); 
}


void FFVV12_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> V4[], std::complex<double>
    COUP, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP1; 
  std::complex<double> TMP10; 
  double P3[4]; 
  std::complex<double> TMP21; 
  std::complex<double> TMP5; 
  std::complex<double> TMP18; 
  std::complex<double> TMP20; 
  std::complex<double> TMP24; 
  std::complex<double> TMP19; 
  std::complex<double> TMP3; 
  std::complex<double> TMP12; 
  double P1[4]; 
  std::complex<double> TMP23; 
  std::complex<double> TMP7; 
  double P4[4]; 
  std::complex<double> TMP22; 
  std::complex<double> TMP0; 
  std::complex<double> TMP17; 
  std::complex<double> TMP4; 
  P1[0] = F1[0].real(); 
  P1[1] = F1[1].real(); 
  P1[2] = F1[1].imag(); 
  P1[3] = F1[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  P4[0] = V4[0].real(); 
  P4[1] = V4[1].real(); 
  P4[2] = V4[1].imag(); 
  P4[3] = V4[0].imag(); 
  TMP24 = (F1[4] * (F2[2] * (V4[2] - V4[5]) - F2[3] * (V4[3] + cI * (V4[4]))) +
      F1[5] * (F2[2] * (+cI * (V4[4]) - V4[3]) + F2[3] * (V4[2] + V4[5])));
  TMP20 = (F1[2] * (F2[4] * (V4[2] + V4[5]) + F2[5] * (V4[3] + cI * (V4[4]))) +
      F1[3] * (F2[4] * (V4[3] - cI * (V4[4])) + F2[5] * (V4[2] - V4[5])));
  TMP21 = (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
      F1[5] * (F2[2] * (+cI * (P3[2]) - P3[1]) + F2[3] * (P3[0] + P3[3])));
  TMP22 = (F1[4] * (F2[2] * (P4[0] - P4[3]) - F2[3] * (P4[1] + cI * (P4[2]))) +
      F1[5] * (F2[2] * (+cI * (P4[2]) - P4[1]) + F2[3] * (P4[0] + P4[3])));
  TMP23 = (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
      F1[5] * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] + V3[5])));
  TMP10 = (P4[0] * P1[0] - P4[1] * P1[1] - P4[2] * P1[2] - P4[3] * P1[3]); 
  TMP17 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
      F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
  TMP19 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
  TMP18 = (F1[2] * (F2[4] * (P4[0] + P4[3]) + F2[5] * (P4[1] + cI * (P4[2]))) +
      F1[3] * (F2[4] * (P4[1] - cI * (P4[2])) + F2[5] * (P4[0] - P4[3])));
  TMP5 = (V4[2] * P1[0] - V4[3] * P1[1] - V4[4] * P1[2] - V4[5] * P1[3]); 
  TMP4 = (P3[0] * P4[0] - P3[1] * P4[1] - P3[2] * P4[2] - P3[3] * P4[3]); 
  TMP7 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP1 = (V3[2] * P4[0] - V3[3] * P4[1] - V3[4] * P4[2] - V3[5] * P4[3]); 
  TMP0 = (V4[2] * P3[0] - V4[3] * P3[1] - V4[4] * P3[2] - V4[5] * P3[3]); 
  TMP3 = (V4[2] * V3[2] - V4[3] * V3[3] - V4[4] * V3[4] - V4[5] * V3[5]); 
  TMP12 = (P3[0] * P1[0] - P3[1] * P1[1] - P3[2] * P1[2] - P3[3] * P1[3]); 
  vertex = COUP * (TMP0 * (TMP10 * (-cI * (TMP19) + cI * (TMP23)) + TMP7 * (-cI
      * (TMP18) + cI * (TMP22))) + (TMP1 * (TMP12 * (-cI * (TMP20) + cI *
      (TMP24)) + TMP5 * (-cI * (TMP17) + cI * (TMP21))) + (TMP3 * (TMP10 * (-cI
      * (TMP21) + cI * (TMP17)) + TMP12 * (-cI * (TMP22) + cI * (TMP18))) +
      TMP4 * (TMP5 * (-cI * (TMP23) + cI * (TMP19)) + TMP7 * (-cI * (TMP24) +
      cI * (TMP20))))));
}


void FFVV10_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> V4[], std::complex<double>
    COUP, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP15; 
  double P3[4]; 
  std::complex<double> TMP16; 
  std::complex<double> TMP14; 
  double P4[4]; 
  std::complex<double> TMP13; 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  P4[0] = V4[0].real(); 
  P4[1] = V4[1].real(); 
  P4[2] = V4[1].imag(); 
  P4[3] = V4[0].imag(); 
  TMP15 = (F2[2] * F1[2] + F2[3] * F1[3]); 
  TMP14 = (-1.) * (P3[0] * (P4[1] * (V4[5] * V3[4] - V4[4] * V3[5]) + (P4[2] *
      (V4[3] * V3[5] - V4[5] * V3[3]) + P4[3] * (V4[4] * V3[3] - V4[3] *
      V3[4]))) + (P3[1] * (P4[0] * (V4[4] * V3[5] - V4[5] * V3[4]) + (P4[2] *
      (V4[5] * V3[2] - V4[2] * V3[5]) + P4[3] * (V4[2] * V3[4] - V4[4] *
      V3[2]))) + (P3[2] * (P4[0] * (V4[5] * V3[3] - V4[3] * V3[5]) + (P4[1] *
      (V4[2] * V3[5] - V4[5] * V3[2]) + P4[3] * (V4[3] * V3[2] - V4[2] *
      V3[3]))) + P3[3] * (P4[0] * (V4[3] * V3[4] - V4[4] * V3[3]) + (P4[1] *
      (V4[4] * V3[2] - V4[2] * V3[4]) + P4[2] * (V4[2] * V3[3] - V4[3] *
      V3[2]))))));
  TMP16 = (F2[4] * F1[4] + F2[5] * F1[5]); 
  TMP13 = (-1.) * (P3[0] * (P4[1] * (V4[4] * V3[5] - V4[5] * V3[4]) + (P4[2] *
      (V4[5] * V3[3] - V4[3] * V3[5]) + P4[3] * (V4[3] * V3[4] - V4[4] *
      V3[3]))) + (P3[1] * (P4[0] * (V4[5] * V3[4] - V4[4] * V3[5]) + (P4[2] *
      (V4[2] * V3[5] - V4[5] * V3[2]) + P4[3] * (V4[4] * V3[2] - V4[2] *
      V3[4]))) + (P3[2] * (P4[0] * (V4[3] * V3[5] - V4[5] * V3[3]) + (P4[1] *
      (V4[5] * V3[2] - V4[2] * V3[5]) + P4[3] * (V4[2] * V3[3] - V4[3] *
      V3[2]))) + P3[3] * (P4[0] * (V4[4] * V3[3] - V4[3] * V3[4]) + (P4[1] *
      (V4[2] * V3[4] - V4[4] * V3[2]) + P4[2] * (V4[3] * V3[2] - V4[2] *
      V3[3]))))));
  vertex = COUP * (TMP13 * (-cI * (TMP16) + cI * (TMP15)) + TMP14 * (-cI *
      (TMP15) + cI * (TMP16)));
}

void FFVV10_11_12_7_8_9_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> V4[], std::complex<double>
    COUP1, std::complex<double> COUP2, std::complex<double> COUP3,
    std::complex<double> COUP4, std::complex<double> COUP5,
    std::complex<double> COUP6, std::complex<double> & vertex)
{
  double P3[4]; 
  double P4[4]; 
  std::complex<double> tmp; 
  FFVV10_0(F1, F2, V3, V4, COUP1, vertex); 
  FFVV11_0(F1, F2, V3, V4, COUP2, tmp); 
  vertex = vertex + tmp; 
  FFVV12_0(F1, F2, V3, V4, COUP3, tmp); 
  vertex = vertex + tmp; 
  FFVV7_0(F1, F2, V3, V4, COUP4, tmp); 
  vertex = vertex + tmp; 
  FFVV8_0(F1, F2, V3, V4, COUP5, tmp); 
  vertex = vertex + tmp; 
  FFVV9_0(F1, F2, V3, V4, COUP6, tmp); 
  vertex = vertex + tmp; 
}

void FFVV9_0(std::complex<double> F1[], std::complex<double> F2[],
    std::complex<double> V3[], std::complex<double> V4[], std::complex<double>
    COUP, std::complex<double> & vertex)
{
  static std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP2; 
  std::complex<double> TMP1; 
  std::complex<double> TMP0; 
  double P3[4]; 
  std::complex<double> TMP4; 
  double P4[4]; 
  std::complex<double> TMP3; 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  P4[0] = V4[0].real(); 
  P4[1] = V4[1].real(); 
  P4[2] = V4[1].imag(); 
  P4[3] = V4[0].imag(); 
  TMP4 = (P3[0] * P4[0] - P3[1] * P4[1] - P3[2] * P4[2] - P3[3] * P4[3]); 
  TMP1 = (V3[2] * P4[0] - V3[3] * P4[1] - V3[4] * P4[2] - V3[5] * P4[3]); 
  TMP0 = (V4[2] * P3[0] - V4[3] * P3[1] - V4[4] * P3[2] - V4[5] * P3[3]); 
  TMP3 = (V4[2] * V3[2] - V4[3] * V3[3] - V4[4] * V3[4] - V4[5] * V3[5]); 
  TMP2 = (F2[2] * F1[2] + F2[3] * F1[3] + F2[4] * F1[4] + F2[5] * F1[5]); 
  vertex = COUP * TMP2 * (-cI * (TMP0 * TMP1) + cI * (TMP3 * TMP4)); 
}


}  // end namespace $(namespace)s_AAttbar_UFO

