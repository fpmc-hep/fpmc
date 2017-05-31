#ifndef FPMC_INC
#define FPMC_INC

//------------------------- FPMC ------------------------------
//
#ifdef __cplusplus
extern "C" {
#endif
  void hwmodini_();
#define hwmodini hwmodini_
  void hwchek_();
#define hwchek hwchek_
  void hwfxer_( int* );
#define hwfxer hwfxer_

//      DOUBLE PRECISION XPQNRM(-6:6),FMCONV,XMSCUT
//      DOUBLE PRECISION AMASS,PMASS
//      CHARACTER*3 TYPEPR, TYPINT
      
//---Normalizations of partonic PDFs inside the pomeron:
//      COMMON /PDFNRM/ XPQNRM
      
//---Atomic and proton numbers of the ion
  extern struct {
    double AION,ZION,RBMIN;
  } ion_;
#define ion ion_
      
//---Flags carrying info whether the process is inclusive or exclusive:
//      PARAMETER(PI=3.141592654D0,FMCONV=0.197,XMSCUT=50.)
//      PARAMETER(AMASS=931.D-3,PMASS=938D-3)
  extern struct {
    char TYPEPR[3], TYPINT[3];
  } prtype_;
#define prtype prtype_
      
//---Model parameters
  extern struct {
    double GAPSPR,PROSPR,CDFFAC;
    int    ISOFTM,NFLUX;
  } xsect_;
#define xsect xsect_

//---PDF parameter
  extern struct {
    int IFITPDF;
  } pdfs_;
#define pdfs pdfs_

//c---FLUX - POMERON + REGEON PARAMETERS
//      DOUBLE PRECISION alphaP,Bpom,alphaPp, zh1
//      DOUBLE PRECISION alphaR,alphaRp,Breg,Cr
//      COMMON /PRPARAMS/ alphaP,Bpom,alphaPp,alphaR,alphaRp,Breg,Cr,zh1

//--- Anomalous coupling variables for anomalous AA->WW, AA-> ZZ, AA-> AA
  extern struct {
    double D_KAPPA, LAMBDA, ANOMCUTOFF, A0W, ACW, A0Z, ACZ,
           A1A, A2A;
    int    AAANOM;
  } aaanomal_;
#define aaanomal aaanomal_
 
//--- Variables for exotic  AA->AA
  extern struct {
    double AAM, AAQ, AAN, AAF0, AAW, AAA2;
    int    AAEXOTIC;
  } aaexotical_;
#define aaexotical aaexotical_
      
//--- CHIDe
  extern struct {
    double CHIDeX, CHIDeXP, CHIDeS2, CHIDeS;
    int    CHIDeIGLU;
  } chidefpmc_;
#define chidefpmc chidefpmc_

//--- KMR2
  extern struct {
    double KMR2Q2CUT, KMR2SURV, KMR2SCALE;
    int    KMR2DELTA;
  } kmr2fpmc_;
#define kmr2fpmc kmr2fpmc_

//--- temporal variables
//      DOUBLE PRECISION TMAXMAX, TMINMIN
//      DOUBLE PRECISION SMAXMAX, SMINMIN
//      INTEGER ISSET, MYDEBUG
//      COMMON /OKTEMP/ ICOUNT
//      INTEGER ICOUNT

//-------------------------------------------------------------
//
  void hwaend_() {}
  void fpmc_var_ini_( int* );
#define fpmc_var_ini fpmc_var_ini_
  void fpmc_welcome_();
#define fpmc_welcome fpmc_welcome_

  extern struct {
    float  URMASS,UWMASS,UHMASS,UTMASS,UMST1,UMSB1,UECMS,
           UYJMIN, UYJMAX,UPTMIN, UPTMAX, UEMMIN, UEMMAX, UDKAPPA,
           UACW, UA0W, UA0Z, UACZ, UA1A, UA2A,
           UAAM, UAAQ, UAAN, UAAF0, UAAW, UAAA2,
           UCHIDeX, UCHIDeXP, UCHIDeS2,
           UXI1Min, UXI1Max, UXI2Min, UXI2Max,
           UCHIDeGapMin, UCHIDeGapMax,
           UKMR2Q2CUT, UKMR2SURV, UKMR2SCALE;
  } myffread1_;
#define myffread1 myffread1_
 
  extern struct {
    float dlambda, anom_cutoff, yww_min, yww_max,
          q2ww_min, q2ww_max;
  } myffread2_;
#define myffread2 myffread2_

  extern struct {
    int UOUTPUT, UOUTPUTLHE, UMAXEV, UIPROC, UNFLUX, UNRN1,
        UNRN2, UIFIT, UISOFTM, UZION, UAION;
    float UBMIN;
    int UAAANOM, UAAEXOTIC,
        UCHIDeIGLU, UKMR2DELTA;
  } myffread3_;
#define myffread3 myffread3_

  extern struct {
    char UHADR[1];
  } cc0_;
#define cc0 cc0_

  extern struct {
    char UTYPEPR[3];
  } cc1_;
#define cc1 cc1_

  extern struct {
    char UTYPINT[3];
  } cc2_;
#define cc2 cc2_

  extern struct {
    char UPART1[4];
  } cc3_;
#define cc3 cc3_

  extern struct {
    char UPART2[4];
  } cc4_;
#define cc4 cc4_

  extern struct {
    int UMODPDF1, UMODPDF2;
  } cc5_;
#define cc5 cc5_

  extern struct {
    char UNTNAME[32], UCHIDePATH[32];
  } cyfflong1_;
#define cyfflong1 cyfflong1_

#ifdef __cplusplus
}
#endif

#endif
