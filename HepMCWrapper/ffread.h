#ifndef HepMCWrapper_ffread_h
#define HepMCWrapper_ffread_h

#include <stdio.h>

#ifdef __cplusplus
extern "C"
{
#endif
  void fpmc_welcome_();
  void ffread_();

  extern struct {
    float r_mass, w_mass, h_mass, t_mass;
    float mst1, msb1;
    float ecms;
    float yj_min, yj_max, pt_min, pt_max, em_min, em_max;
    float dkappa; // ATGC
    float acw, a0w, a0z, acz; // AQGC
    float a1a, a2a;
    float aam, aaq, aan, aaf0, aaf0z, aaf0w, aaf0zg, aaw, aaa2;
    float chide_x, chide_xp, chide_s2;
    float xi1_min, xi1_max, xi2_min, xi2_max;
    float chide_gap_min, chide_gap_max;
    float kmr2_q2_cut, kmr2_surv, kmr2_scale;
  } myffread1_;
  extern struct {
    float dlambda, anomcutoff;
    float yww_min, yww_max;
    float q2ww_min, q2ww_max;
  } myffread2_;
  extern struct {
    int output, outputlhe;
    float maxev;
    int iproc, nflux, nrn1, nrn2, ifit, isoftm, zion, aion;
    float bmin;
    int aaanom, aaexotic;
    int chide_iglu, kmr2_delta;
  } myffread3_;
  extern struct { char hadr; } cc0_;
  extern struct { char typepr[3]; } cc1_;
  extern struct { char typint[3]; } cc2_;
  extern struct { char part1[4]; } cc3_;
  extern struct { char part2[4]; } cc4_;
  extern struct { int modpdf1, modpdf2; } cc5_;
  extern struct { char ntname[32], chide_path[32]; } cyfflong1_;
  extern struct { char lhe_file[64]; } cyfflong2_;
#ifdef __cplusplus
}
#endif

void ffread_();
void ffinit();

#endif
