#ifndef FPMC_INC
#define FPMC_INC

//------------------------- FPMC ------------------------------
//
#ifdef __cplusplus
extern "C"
{
#endif
  void hwmodini_();
#define hwmodini hwmodini_
  void hwchek_();
#define hwchek hwchek_
  void hwfxer_( int* );
#define hwfxer hwfxer_

  // PARAMETER(PI=3.141592654D0,FMCONV=0.197,XMSCUT=50.)
  // PARAMETER(AMASS=931.D-3,PMASS=938D-3)
  extern double pi_; // 3.141592654
  extern double fmconv_; // 0.197
  extern double xmscut_; // 50.
  extern double amass_; // 931e-3
  extern double pmass_; // 938e-3

  //--- normalizations of partonic PDFs inside the pomeron:
  typedef struct {
    double XPQNRM[12];
  } pdfnrm_t; // warning: XPQNRM(-6:6)
  //extern pdfnrm_t pdfnrm_;
      
  //--- atomic and proton numbers of the ion
  typedef struct {
    double AION, ZION, RBMIN;
  } ion_t;
  extern ion_t ion_;
      
  //--- flags carrying info whether the process is inclusive or exclusive:
  typedef struct {
    char TYPEPR[3], TYPINT[3];
  } prtype_t;
  extern prtype_t prtype_;

  //--- model parameters
  typedef struct {
    double GAPSPR,PROSPR,CDFFAC;
    int    ISOFTM,NFLUX;
  } xsect_t;
  extern xsect_t xsect_;

  //--- PDF parameter
  typedef struct {
    int IFITPDF;
  } pdfs_t;
  extern pdfs_t pdfs_;

  //--- FLUX - POMERON + REGEON PARAMETERS
  typedef struct {
    double alphaP, Bpom, alphaPp, alphaR, alphaRp, Breg, Cr, zh1;
  } prparams_t;
  //extern prparams_t prparams_;

  //--- anomalous coupling variables for anomalous AA->WW, AA-> ZZ, AA-> AA
  typedef struct {
    double D_KAPPA, LAMBDA, ANOMCUTOFF, A0W, ACW, A0Z, ACZ,
           A1A, A2A, XI1TTBAR, XI2TTBAR, XI3TTBAR, 
           XI4TTBAR, XI5TTBAR, XI6TTBAR;
    int    AAANOM;
  } aaanomal_t;
  extern aaanomal_t aaanomal_;
 
  //--- variables for exotic  AA->AA
  //--- and in the future AA->WW and AA->AZ
  //--- (C. Baldenegro, 05-2016)
  typedef struct {
    double AAM, AAQ, AAN,
           AAF0, AAF0Z, AAF0W, AAF0ZG, AAF0G, AAF0H,
           AAW, AAA2;
    int    AAEXOTIC;
  } aaexotical_t;
  extern aaexotical_t aaexotical_;
      
  //--- CHIDe
  typedef struct {
    double CHIDeX, CHIDeXP, CHIDeS2, CHIDeS;
    int    CHIDeIGLU;
  } chidefpmc_t;
  extern chidefpmc_t chidefpmc_;

  //--- KMR2
  typedef struct {
    double KMR2Q2CUT, KMR2SURV, KMR2SCALE;
    int    KMR2DELTA;
  } kmr2fpmc_t;
  extern kmr2fpmc_t kmr2fpmc_;

  //--- temporal variables
  typedef struct {
    int ICOUNT;
  } oktemp_t;
  //extern oktemp_t oktemp_;

  //--- custom definition of subroutine hwaend()
  void hwaend_();

  void fpmc_welcome_();
#define fpmc_welcome fpmc_welcome_
  void fpmc_var_ini_( int* );
#define fpmc_var_ini fpmc_var_ini_

  //-------------------------------------------------------------
  // outdated common blocks
  //-------------------------------------------------------------

  typedef struct {
    float  RMASS, WMASS, HMASS, TMASS, MST1, MSB1, ECMS,
           YJMIN, YJMAX, PTMIN, PTMAX, EMMIN, EMMAX,
           DKAPPA,
           ACW, A0W, A0Z, ACZ, A1A, A2A,
           AAM, AAQ, AAN, AAF0, AAW, AAA2,
           CHIDeX, CHIDeXP, CHIDeS2,
           XI1Min, XI1Max, XI2Min, XI2Max,
           CHIDeGapMin, CHIDeGapMax,
           KMR2Q2CT, KMR2SRV, KMR2SCALE,
           XI1TTBAR, XI2TTBAR, XI3TTBAR, 
           XI4TTBAR, XI5TTBAR, XI6TTBAR;
  } myffread1_t;
  //extern myffread1_t myffread1_;

  typedef struct {
    float dlambda, anom_cutoff, yww_min, yww_max,
          q2ww_min, q2ww_max;
  } myffread2_t;
  //extern myffread2_t myffread2_;

  typedef struct {
    int OUTPUT, OUTPUTLHE, MAXEV, IPROC, NFLUX, NRN1,
        NRN2, IFIT, ISOFTM, ZION, AION;
    float BMIN;
    int AAANOM, AAEXOTIC,
        CHIDeIGL, KMR2DELTA;
  } myffread3_t;
  //extern myffread3_t myffread3_;

  typedef struct {
    char HADR[1];
  } cc0_t;
  //extern cc0_t cc0_;

  typedef struct {
    char TYPEPR[3];
  } cc1_t;
  //extern cc1_t cc1_;

  typedef struct {
    char TYPINT[3];
  } cc2_t;
  //extern cc2_t cc2_;

  typedef struct {
    char PART1[4];
  } cc3_t;
  //extern cc3_t cc3_;

  typedef struct {
    char PART2[4];
  } cc4_t;
  //extern cc4_t cc4_;

  typedef struct {
    int MODPDF1, MODPDF2;
  } cc5_t;
  //extern cc5_t cc5_;

  typedef struct {
    char NTNAME[32], CHIDePATH[32];
  } cyfflong1_t;
  //extern cyfflong1_t cyfflong1_;

  typedef struct {
    char dgdtab1[50], dgdtab2[50], dgdtab3[50], dgdtab4[50], sudatab[50];
  } chidepath_t;
  extern chidepath_t chidepath_;

  typedef struct {
    char CHIDe_PATH[500];
  } chideenv_t;
  extern chideenv_t chideenv_;

#ifdef __cplusplus
}
#endif

#endif
