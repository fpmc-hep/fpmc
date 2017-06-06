#ifndef HERWIG_INC
#define HERWIG_INC

//-------------------------HERWIG common block -------------------------------------
// COMMON block stuff, that doesn't come with the HerwigWrapper6_4.h ....

/*C Arrays for particle properties (NMXRES = max no of particles defined)
      PARAMETER(NMXRES=500)
      COMMON/HWPROP/RLTIM(0:NMXRES),RMASS(0:NMXRES),RSPIN(0:NMXRES),
     & ICHRG(0:NMXRES),IDPDG(0:NMXRES),IFLAV(0:NMXRES),NRES,
     & VTOCDK(0:NMXRES),VTORDK(0:NMXRES),
     & QORQQB(0:NMXRES),QBORQQ(0:NMXRES) */

static const int nmxres = 500+1; // we need NMXRES+1 entries ...
static const int nmxsud = 1024; // max number of entries for Sudakov lookup table
static const int nmxcdk = 4000;
static const int nmxhep = 4000;
static const int hepevt_size = 4000; // check in HerwigWrapper
static const int modmax = 50;
static const int maxhrp = 100;
static const int imaxch = 20;

static const int NPROC = 117;
static const int MAXMS = 100;
static const int NPSIMP = 16;
static const double SMALL = 1.e-20;

#ifdef __cplusplus
extern "C"
{
#endif

  typedef struct {
    double AFCH[2][16], ALPHEM, B1LIM, BETAF, BTCLM, CAFAC, CFFAC, CLMAX, CLPOW, CLSMR[2], CSPEED;
    double ENSOF, ETAMIX, F0MIX, F1MIX, F2MIX, GAMH, GAMW, GAMZ, GAMZP, GEV2NB;
    double H1MIX, PDIQK, PGSMX, PGSPL[4], PHIMIX, PIFAC, PRSOF, PSPLT[2], PTRMS, PXRMS;
    double QCDL3, QCDL5, QCDLAM, QDIQK, QFCH[16], QG, QSPAC, QV;
    double SCABI, SWEIN, TMTOP, VFCH[2][16], VCKM[3][3], VGCUT, VQCUT, VPCUT;
    double ZBINM, EFFMIN, OMHMIX, ET2MIX, PH3MIX, GCUTME;
    int IOPREM, IPRINT, ISPAC, LRSUD, LWSUD;
    int MODPDF[2], NBTRY, NCOLO, NCTRY, NDTRY, NETRY, NFLAV, NGSPL, NSTRU, NSTRY, NZBIN;
    int IOP4JT[2], NPRFMT;
    int AZSOFT, AZSPIN; // logical
    int CLDIR[2];
    int HARDME, NOSPAC, PRNDEC, PRVTX, SOFTME, ZPRIME, PRNDEF, PRNTEX, PRNWEB; // logical
  } hwpram_t;
  extern hwpram_t hwpram_;

  typedef struct {
    double AVWGT, EVWGT, GAMWT, TLOUT, WBIGST, WGTMAX, WGTSUM, WSQSUM;
    int IDHW[nmxhep], IERROR, ISTAT, LWEVT, MAXER, MAXPR;
    int NOWGT; // logical
    int NRN[2], NUMER, NUMERU, NWGTS;
    int GENSOF; // logical
  } hwevnt_t;
  extern hwevnt_t hwevnt_;

  typedef struct { char PART1[8], PART2[8]; } hwbmch_t;
  extern hwbmch_t hwbmch_;

  typedef struct {
    double EBEAM1, EBEAM2, PBEAM1, PBEAM2;
    int IPROC, MAXEV;
  } hwproc_t;
  extern hwproc_t hwproc_;

  typedef struct {
    int NDECSY, NSEARCH, LRDEC, LWDEC;
    int SYSPIN, THREEB, FOURB; // logical
    char TAUDEC[6];
  } hwdspn_t;
  extern hwdspn_t hwdspn_;

  typedef struct {
    char AUTPDF[2][20];
    char BDECAY[4];
  } hwprch_t;
  extern hwprch_t hwprch_;

    
  typedef struct {
    double RLTIM[nmxres], RMASS[nmxres], RSPIN[nmxres];
    int ICHRG[nmxres], IDPDG[nmxres],IFLAV[nmxres], NRES;
    int VTOCDK[nmxres], VTORDK[nmxres], QORQQB[nmxres], QBORQQ[nmxres];    
  } hwprop_t;
  extern hwprop_t hwprop_;

  //--- parameters for Sudakov form factors
  typedef struct {
    double ACCUR, QEV[6][nmxsud],SUD[6][nmxsud];
    int INTER, NQEV, NSUD, SUDORD;
  } hwusud_t;
  extern hwusud_t hwusud_;

  //--- new parameters for version 6.203
  typedef struct {
    double ABWGT, ABWSUM, AVABW;
    int NNEGWT,NNEGEV,NEGWTS;
  } hw6203_t;
  extern hw6203_t hw6203_;

  //--- parameters for minimum bias/soft underlying event
  typedef struct {
    double PMBN1,PMBN2,PMBN3,PMBK1,PMBK2,PMBM1,PMBM2,PMBP1,PMBP2,PMBP3;
  } hwminb_t;
  extern hwminb_t hwminb_;

  //--- variables controling mixing and vertex information
  // VTXPIP should have been a 5-vector, problems with NAG compiler
  typedef struct {
    double EXAG,GEV2MM,HBAR,PLTCUT,VMIN2,VTXPIP[5],XMIX[2],XMRCT[2],YMIX[2],YMRCT[2];
    int IOPDKL,MAXDKL,MIXING,PIPSMR;
  } hwdist_t;
  extern hwdist_t hwdist_;

  typedef struct {
    double CLDKWT[nmxcdk],CTHRPW[12][12],PRECO,RESN[12][12], RMIN[12][12];
    int LOCN[12][12],NCLDK[nmxcdk], NRECO,CLRECO;
  } hwuclu_t;
  extern hwuclu_t hwuclu_; 

  //--- weights used in cluster decays
  typedef struct {
    double REPWT[5][5][4],SNGWT,DECWT,QWT[3],PWT[12],SWTEF[nmxres];
  } hwuwts_t;
  extern hwuwts_t hwuwts_;

  //--- other new parameters for version 6.2
  typedef struct {
    double VIPWID[3], DXRCYL,DXZMAX,DXRSPH;
    int WZRFR,FIX4JT,IMSSM,IHIGGS,PARITY,LRSUSY;
  } hw6202_t;
  extern hw6202_t hw6202_;

  typedef struct {
    double ALPFAC, BRHIG[12], ENHANC[12], GAMMAX, RHOHEP[hepevt_size][3];
    int IOPHIG, MODBOS[modmax];
  } hwbosc_t;
  extern hwbosc_t hwbosc_;

  typedef struct {
    double ASFIXD, CLQ[6][7], COSS, COSTH, CTMAX;
    double DISF[2][13], EMLST, EMMAX, EMMIN, EMPOW, EMSCA, EPOLN[3];
    double GCOEF[7], GPOLN;
    double OMEGA0, PHOMAS, PPOLN[3];
    double PTMAX, PTMIN, PTPOW;
    double Q2MAX, Q2MIN, Q2POW, Q2WWMN, Q2WWMX, QLIM;
    double SINS, THMAX, Y4JT, TMNISR, TQWT, XX[2], XLMIN, XXMIN;
    double YBMAX, YBMIN, YJMAX, YJMIN, YWWMAX, YWWMIN, WHMIN, ZJMAX, ZMXISR;
    int IAPHIG, IBRN[2], IBSH, ICO[10], IDCMF, IDN[10], IFLMAX, IFLMIN, IHPRO, IPRO;
    int MAPQ[6], MAXFL, BGSHAT, COLISR;
    int FSTEVT, FSTWGT, GENEV, HVFCEN, TPOL, DURHAM;
  } hwhard_t;
  extern hwhard_t hwhard_;

  //--- other HERWIG branching, event and hard subprocess common blocks
  typedef struct {
    double ANOMSC[2][2],HARDST,PTINT[2][3],XFACT;
    int INHAD,JNHAD,NSPAC[7],ISLENT,BREIT,FROST,USECMF;
  } hwbrch_t;
  extern hwbrch_t hwbrch_;

  typedef struct {
    int PRESPL;
  } hw6500_t;
  extern hw6500_t hw6500_;

  //--- R-Parity violating parameters and colours
  typedef struct {
    double LAMDA1[3][3][3],LAMDA2[3][3][3],LAMDA3[3][3][3];
    int HRDCOL[5][2],RPARTY,COLUPD;
  } hwrpar_t;
  extern hwrpar_t hwrpar_;

  //--- SUSY parameters
  typedef struct {
    double TANB,ALPHAH,COSBPA,SINBPA,COSBMA,SINBMA,COSA,SINA,COSB,SINB,COTB,
           ZMIXSS[4][4],ZMXNSS[4][4],ZSGNSS[4],
           LFCH[16],RFCH[16],SLFCH[4][16],SRFCH[4][16],
           WMXUSS[2][2],WMXVSS[2][2],WSGNSS[2],
           QMIXSS[2][2][6],LMIXSS[2][2][6],
           THETAT,THETAB,THETAL,ATSS,ABSS,ALSS,MUSS,FACTSS,
           GHWWSS[3],GHZZSS[3],GHDDSS[4],GHUUSS[4],GHWHSS[3],GHSQSS[2][2][6][4],
           XLMNSS,RMMNSS,DMSSM,
           SENHNC[24],SSPARITY;
    int SUSYIN;
  } hwsusy_t;
  extern hwsusy_t hwsusy_;

  //--- common block for Les Houches interface to store information we need
  typedef struct {
    double LHWGT[maxhrp],LHWGTS[maxhrp],LHMXSM,LHXSCT[maxhrp],LHXERR[maxhrp],LHXMAX[maxhrp];
    int LHIWGT,LHNEVT,ITYPLH,LHSOFT,LHGLSF;
  } hwgupr_t;
  extern hwgupr_t hwgupr_;

  //--- new parameters for version 6.3
  typedef struct {
    double MJJMIN,CHNPRB[imaxch];
    int IOPSTP,IOPSH,OPTM,CHON[imaxch];
  } hw6300_t;
  extern hw6300_t hw6300_;

  //-------------------------- JIMMY COMMON BLOCK -------------------------------
/*
      DOUBLE PRECISION YGAMMA, JMZMIN, JMRAD, PTJIM
      DOUBLE PRECISION PHAD, JMU2, JMV2, SMALL, JMARRY
c     JMARRY is the array storing gamma-p xsec at various z, & 
c	max weight for each z
      DOUBLE PRECISION TOTSCAT, NLOST

      INTEGER MAXMS, NPSIMP, MSFLAG, JMPTYP, JCMVAR, NPROC
      LOGICAL ANOMOFF

      PARAMETER( NPROC = 117 )
      PARAMETER( MAXMS  = 100  )  ! Maximum multiple scatters
      PARAMETER( NPSIMP = 16 )  ! No. of Simpson rule (YBJ)
C                                 intervals (must be even)
      PARAMETER( SMALL  = 1.0D-20  )
      INTEGER JMOUT, JMBUG, FN_TYPE, NSCAT, JMUEO, MAXMSTRY
      PARAMETER(JMOUT = 6)
      COMMON / JMPARM /  PTJIM, YGAMMA, JMZMIN, JMRAD(264)
     &     ,PHAD, JMU2, JMV2, JMARRY( 6+MAXMS,0:NPSIMP )
     &     ,NLOST, TOTSCAT, ANOMOFF, JCMVAR, JMUEO
     &     ,JMPTYP(NPROC), JMBUG, FN_TYPE, MSFLAG, MAXMSTRY
      DOUBLE PRECISION JMPROC, JMVETO
      COMMON / JMEVNT/ JMPROC(NPROC)
     &,        JMVETO(2,13), NSCAT
*/

  typedef struct {
    double PTJIM,YGAMMA,JMZMIN,JMRAD[264],PHAD,JMU2,JMV2,JMARRY[NPSIMP+1][6+MAXMS],
      NLOST,TOTSCAT;
    int ANAMOFF,JCMVAR,JMUEO,JMPTYP[NPROC],JMBUG,FN_TYPE,MSFLAG,MAXMSTRY;
  } jmparm_t;
  extern jmparm_t jmparm_;

  typedef struct {
    double JMPROC[NPROC],JMVETO[13][2];
    int NSCAT;
  } jmevnt_t;
  extern jmevnt_t jmevnt_;

  //------------------------------ JIMMY functions -------------------------------------------------
  void jimmin_();
#define jimmin jimmin_
  void jminit_();
#define jminit jminit_
  double hwmsct_( int* );
#define hwmsct hwmsct_
  void jmefin_();
#define jmefin jmefin_

  //------------------------------ LHAPDF functions -------------------------------------------------

  //--- MLM-style matching for aMC@NLO
  typedef struct {
    int n_match, max_multiplicity_flag;
    double matching_scale;
  } hwmatchpram_t;
  extern hwmatchpram_t hwmatchpram_;

  void hwmatch_( int* );
#define hwmatch hwmatch_
  void hwhdecay_();
#define hwhdecay hwhdecay_

  //---------------------------------------------------------------
  void hwuidt_( int* iopt, int* ipdg, int* iwig, char nwig[8] );
#define hwuidt hwuidt_
  double hwualf_( int *mode, double* scale );
#define hwualf hwualf_
  double hwuaem_( double* scale );
#define hwuaem hwuaem_

  //---------------------------------------------------------------
  void hweini_(); // initialise elementary process
#define hweini hweini_
  void hwuine_(); // initialise event
#define hwuine hwuine_
  void hwepro_(); // generate HERWIG hard subprocess
#define hwepro hwepro_
  void hwigin_(); // initialise other common blocks
#define hwigin hwigin_
  void hwusta_( const char*, int ); // make any particle stable
#define hwusta hwusta_
  void hwuinc_(); // compute parameter-dependent constants
#define hwuinc hwuinc_
  void hwbgen_(); // generate parton cascades
#define hwbgen hwbgen_
  void hwdhob_(); // do heavy object decays
#define hwdhob hwdhob_
  void hwcfor_(); // do cluster hadronization
#define hwcfor hwcfor_
  void hwcdec_(); // do cluster decay
#define hwcdec hwcdec_
  void hwdhad_(); // do unstable particle decays
#define hwdhad hwdhad_
  void hwdhvy_(); // do heavy flavour decays
#define hwdhvy hwdhvy_
  void hwmevt_(); // add soft underlying event if needed
#define hwmevt hwmevt_
  void hwufne_(); // event generation completed, wrap up event, ...
#define hwufne hwufne_

#ifdef __cplusplus
}
#endif
//---------------------------------------------------------------

#endif
