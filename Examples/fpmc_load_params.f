C-----------------------------------------------------------------------
C * 21/06/2018 L. Forthomme
C * Set parameters via FF cards
C-----------------------------------------------------------------------
      SUBROUTINE FPMC_LOAD_PARAMS
c---Common block
      INCLUDE 'HERWIG65.INC'
      INCLUDE 'fpmc.inc'
      INCLUDE 'CHIDe.inc'
c---User's declarations
      DOUBLE PRECISION CMSENR
      INTEGER N
c---External
      INCLUDE 'ffcard.inc'

c---Apply the setup parameters - either default/ffcard

c---Events in the run
      MAXEV = UMAXEV

c---Beam particles (In DPE as implemented by POMWIG E+ stands for initial state proton)
      PART1=UPART1
      PART2=UPART2

c---Set cms energy
      CMSENR = UECMS    !set by a card

c---Beam momenta
      PBEAM1=CMSENR/2d0
      PBEAM2=PBEAM1

c---Process type : EXClusive or INClusive, QCD or QED when relevant
c   It is important to set those right
c      TYPEPR='INC'
c      TYPINT='QCD'
      TYPEPR = UTYPEPR
      TYPINT = UTYPINT

c---Setting the hard subprocess (see manual)
      IPROC = UIPROC

c---Option to include hadronization and showering effects
      ANSWER=UHADR

c---Initialization of other Herwig common blocks
c...The default parameters has to be changed after this call
      CALL HWIGIN
cc      SUSYIN=.TRUE.

c---Random number genrator initializaiton
      NRN(1) = UNRN1 ! set again later in the code
      NRN(2) = UNRN2 ! set again later in the code

C ... increase inefficiency tolerance
      EFFMIN = 1d-6

C---Initialize LHAPDF
      AUTPDF(1) = "HWLHAPDF"
      MODPDF(1) = UMODPDF1
      AUTPDF(2) = "HWLHAPDF"
      MODPDF(2) = UMODPDF2

c---User's default kinematic parameters
c... beam momentum transfer range ( Q2 = |t| )
      RMASS(201) = UHMASS
      RMASS(6)   = UTMASS ! top mass
      RMASS(198) = UWMASS ! W mass
      RMASS(406) = UMST1 ! stop1 mass
      RMASS(405) = UMSB1 ! stop2 mass
c      RMASS(4)   = 1.5 ! charm mass
c      RMASS(5)   = 4.5 ! bottom mass

c stop mass
c      RMASS(449)=10000.D0
c      RMASS(401)=10000.D0
c      RMASS(402)=10000.D0
c      RMASS(403)=10000.D0
c      RMASS(404)=10000.D0
c      RMASS(405)=10000.D0
c      RMASS(407)=10000.D0
c      RMASS(408)=10000.D0
c      RMASS(409)=10000.D0
c      RMASS(410)=10000.D0
c      RMASS(411)=10000.D0

c      RMASS(413)=10000.D0
c      RMASS(414)=10000.D0
c      RMASS(415)=10000.D0
c      RMASS(416)=10000.D0
c      RMASS(417)=10000.D0
c      RMASS(418)=10000.D0
c      RMASS(419)=10000.D0
c      RMASS(420)=10000.D0
c      RMASS(421)=10000.D0
c      RMASS(422)=10000.D0
c      RMASS(423)=10000.D0
c      RMASS(424)=10000.D0
c      RMASS(406)=393.D0
c      RMASS(412)=393.D0

      Q2WWMN=UQ2WWMN
      Q2WWMX=UQ2WWMX

c      Q2WWMN=0
c      Q2WWMX=4

c... beam momentum loss range
c      YWWMIN=0.d0
c      YWWMAX=1.D0
      YWWMIN=UYWWMIN
      YWWMAX=UYWWMAX

c--- central products : rapidity, pT or mass cuts
c      YJMAX = 5d0
c---JC      YJMAX = 4d0
      YJMAX = UYJMAX
      YJMIN = UYJMIN

c      PTMIN = 25d0
c---JC
      PTMIN=UPTMIN
      PTMAX=UPTMAX

c masse des leptons / pas photons, pas jets
c---JC      EMMIN=0d0
      EMMIN=UEMMIN

c---Choosing the flux : 9,10 is Cox-Forshaw (DPE)
c                       11 is Bialas-Landshoff (exc. DPE) or Boonekamp et al (inc. DPE)
c                       12 is Cahn, Jackson (heavy ions, QED)
c                       13 is Drees, Ellis, Zeppenfled (heavy ions, QED)
c                       14 is Papageorgiu (proton, QED)
c                       15 is Budnev flux (RECOMMENDED for proton, QED)
c                       16 is KMR flux (tables from L.Lonnblad or ExHume)
      NFLUX = UNFLUX

c---Choosing pdf to use:
c                       according to hep/ph 0609291
c                       10:  H1
c                       20:  Zeus
c                       30:  H1Zeus combined
c                       older versions for compat., see h1qcd.f
c                        2:  NLO fit as in H1
c                        5:  LO  fit as in H1
c                        8:  extended version of H1
      IFITPDF = UIFIT

c ---Anomalous coupling settings
      AAANOM = UAAANOM
      D_KAPPA = UDKAPPA
      LAMBDA = UDLAMBDA
      A0W = UA0W
      ACW = UACW
      A0Z = UA0Z
      ACZ = UACZ
      A1A = UA1A
      A2A = UA2A
      ANOMCUTOFF = UANOMCUTOFF

c ---Ttbar EFT settings
      XI1TTBAR = UXI1TTBAR
      XI2TTBAR = UXI2TTBAR
      XI3TTBAR = UXI3TTBAR
      XI4TTBAR = UXI4TTBAR
      XI5TTBAR = UXI5TTBAR
      XI6TTBAR = UXI6TTBAR

c ---Exotic excl AA settings
      AAEXOTIC = UAAEXOTIC
      AAM = UAAM
      AAQ = UAAQ
      AAN = UAAN
      AAF0 = UAAF0

c ----aa->phi->zz exclusive production coupling
      AAF0Z = UAAF0Z
c ----aa->phi->zz exclusive production coupling
      AAF0W = UAAF0W
c ----aa->phi->az exclusive production coupling
      AAF0ZG = UAAF0ZG
      AAW = UAAW
      AAA2 = UAAA2

C ----CHIDe Model
      CHIDeIGLU = UCHIDeIGLU
      CHIDeX   =  UCHIDeX
      CHIDeXP  =  UCHIDeXP
      CHIDeS2  =  UCHIDeS2
      CHIDeS   =  UECMS*UECMS
      XI1MIN = UXI1MIN
      XI1MAX = UXI1MAX
      XI2MIN = UXI2MIN
      XI2MAX = UXI2MAX
      CHIDeGapMin = UCHIDeGapMin
      CHIDeGapMax = UCHIDeGapMax
      CHIDePATH = UCHIDePATH

C ----KMR2
      KMR2Q2CUT=UKMR2Q2CUT
      KMR2SURV=UKMR2SURV
      KMR2SCALE=UKMR2SCALE
      KMR2DELTA=UKMR2DELTA

c---Set normalisations for Exclusive processes
c... option for soft corrections : ISOFTM = 0 : no correction
c                                           1 : simple factor, see below
c                                           2 : KMR low mass diffractive
c                                           3 : effective opacity model
      ISOFTM = UISOFTM

c---Events printed
      MAXPR=10

c---Set the number of nucleons - for QED photon flux
      ZION = UZION
      AION = UAION
      RBMIN = UBMIN

      write(*,*) ''
      write(*,*) 'USER SETTINGS'
      write(*,*) '-------------'
      write(*,*) 'OUTPUT    = ',UOUTPUT
      write(*,*) 'OUTPUTLHE = ',UOUTPUTLHE
      write(*,*) 'NTNAME    = ',UNTNAME
      write(*,*) 'LHE FILE  = ',ULHEFILE
      write(*,*) 'MAXEV     = ',MAXEV
      write(*,*) 'TYPEPR    = ',TYPEPR
      write(*,*) 'TYPEINT   = ',TYPINT
      write(*,*) 'HADR      = ', ANSWER
      write(*,*) 'PART1     = ', PART1
      write(*,*) 'PART2     = ', PART2
      write(*,*) 'NRN(1)    = ',NRN(1)
      write(*,*) 'NRN(2)    = ',NRN(2)
      write(*,*) 'WMASS     = ',UWMASS
      write(*,*) 'TMASS     = ',UTMASS
      write(*,*) 'HMASS     = ',UHMASS
      write(*,*) 'MSTOP1    = ',UMST1
      write(*,*) 'MSBOY1    = ',UMSB1
      write(*,*) 'YJMIN     = ', YJMIN
      write(*,*) 'YJMAX     = ',YJMAX
      write(*,*) 'PTMIN     = ',PTMIN
      write(*,*) 'EMMIN     = ',EMMIN
      write(*,*) 'IFIT      = ',IFITPDF
      write(*,*) 'ISOFTM    = ',ISOFTM
      write(*,*) 'Q2WWMN    = ',Q2WWMN
      write(*,*) 'Q2WWMX    = ',Q2WWMX
      write(*,*) 'YWWMIN    = ',YWWMIN
      write(*,*) 'YWWMAX    = ',YWWMAX
      write(*,*) 'ZION      = ',ZION
      write(*,*) 'AION      = ',AION
      write(*,*) 'BMIN      = ',RBMIN
      write(*,*) 'AAANOM    = ',AAANOM
      write(*,*) 'DKAPPA    = ',D_KAPPA
      write(*,*) 'DLAMBDA   = ',LAMBDA
      write(*,*) 'ANOMCUTOFF= ',ANOMCUTOFF
      write(*,*) 'AAEXOTIC  = ',AAEXOTIC
      write(*,*) 'AUTPDF(1) = ',AUTPDF(1)
      write(*,*) 'MODPDF(1) = ',MODPDF(1)
      write(*,*) 'AUTPDF(2) = ',AUTPDF(2)
      write(*,*) 'MODPDF(2) = ',MODPDF(2)
      write(*,*)

      END

