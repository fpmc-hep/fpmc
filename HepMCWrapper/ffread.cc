#include "ffread.h"

void ffread_()
{
  fpmc_welcome_();
  ffinit();
}

void ffread()
{
  
}

void ffinit()
{
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
