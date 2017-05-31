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

void ffread_()
{
  fpmc_welcome_();
  // FFC default initialiszation
  myffread1_.r_mass = 0.;
  myffread1_.w_mass = 80.425;
  myffread1_.h_mass = 125.0;
  myffread1_.t_mass = 174.3;
  myffread1_.mst1 = 250.;
  myffread1_.msb1 = 250.;
  myffread1_.ecms = 14.e3;
  myffread1_.yj_min = -6.;
  myffread1_.yj_max = 6.;
  myffread1_.pt_min = 0.;
  myffread1_.pt_max = 1.e8;
  myffread1_.em_min = 10.;
  myffread1_.em_max = 1.e8;
  myffread1_.dkappa = 0.;
  myffread1_.acw = 0.;
  myffread1_.a0w = 0.;
  myffread1_.a0z = 0.;
  myffread1_.acz = 0.;
  myffread1_.a1a = 0.;
  myffread1_.a2a = 0.;
  myffread1_.aam = 0.;
  myffread1_.aaq = 0.;
  myffread1_.aan = 0.;
  myffread1_.aaf0 = 0.;
  myffread1_.aaf0z = 0.;
  myffread1_.aaf0w = 0.;
  myffread1_.aaf0zg = 0.;
  myffread1_.aaw = 0.;
  myffread1_.aaa2 = 0.;
  myffread1_.chide_x = -1.;
  myffread1_.chide_xp = -1.;
  myffread1_.chide_s2 = -1.;
  myffread1_.xi1_min = -1.;
  myffread1_.xi1_max = -1.;
  myffread1_.xi2_min = -1.;
  myffread1_.xi2_max = -1.;
  myffread1_.chide_gap_min = 0.;
  myffread1_.chide_gap_max = 0.;
  myffread1_.kmr2_q2_cut = 2.;
  myffread1_.kmr2_surv = 0.3;
  myffread1_.kmr2_scale = 0.618;
  //
  myffread2_.dlambda = 0.;
  myffread2_.anomcutoff = -1.;
  myffread2_.yww_min = 0.;
  myffread2_.yww_max = 0.1;
  myffread2_.q2ww_min = 0.;
  myffread2_.q2ww_max = 4.;
  //
  myffread3_.output = 1;
  myffread3_.outputlhe = 0;
  myffread3_.maxev = 1000;
  myffread3_.iproc = 16010;
  myffread3_.nflux = 15;
  myffread3_.nrn1 = 33799;
  myffread3_.nrn2 = 11799;
  myffread3_.ifit = 10;
  myffread3_.isoftm = 1;
  myffread3_.zion = 1;
  myffread3_.aion = 1;
  myffread3_.bmin = 1.;
  myffread3_.aaanom = 0;
  myffread3_.aaexotic = 0;
  myffread3_.chide_iglu = -1;
  myffread3_.kmr2_delta = 1;
  //
  cc0_.hadr = 'Y';
  //
  sprintf(cc1_.typepr, "EXC");
  //
  sprintf(cc2_.typint, "QED");
  //
  sprintf(cc3_.part1, "E+");
  //
  sprintf(cc4_.part2, "E+");
  //
  cc5_.modpdf1 = -1;
  cc5_.modpdf2 = -1;
  //
  sprintf(cyfflong1_.ntname, "tmpntuple.ntp");
  sprintf(cyfflong1_.chide_path, "External/CHIDe/Data/");
  //
  sprintf(cyfflong2_.lhe_file, "FPMC.lhe");
}
