C-----------------------------------------------------------------------
C * 15/11/2006 O. Kepka
C * Multifunctional module with cms fast simulation
C * Set parameters via FF cards
C-----------------------------------------------------------------------
      PROGRAM FPMC
c---Common block
      INCLUDE 'HERWIG65.INC'
      INCLUDE 'fpmc.inc'
      INCLUDE 'CHIDe.inc'
c---User's declarations
      DOUBLE PRECISION CMSENR
      INTEGER N
      CHARACTER ANSWER
      LOGICAL READCARD
c---External
      INCLUDE 'ffcard.inc'

c---Initialize FPMC setup parameters and read ffcards if READCARD = true
      READCARD = .TRUE.
      CALL FPMC_VAR_INI(READCARD)

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
c      ANSWER='Y' ! say 'N' to skip this part
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
      ANOMCUTOFF = UANOMCUTOFF

C ... begin R.S.
C     CHIDe Model
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

C ... end R.S.


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

      write(*,*) ''
      write(*,*) 'USER SETTINGS'
      write(*,*) '-------------'
      write(*,*) 'OUTPUT   = ',UOUTPUT
      write(*,*) 'NTNAME   = ',UNTNAME
      write(*,*) 'MAXEV    = ',MAXEV
      write(*,*) 'TYPEPR   = ',TYPEPR
      write(*,*) 'TYPEINT  = ',TYPINT
      write(*,*) 'HADR     = ', ANSWER
      write(*,*) 'PART1    = ', PART1
      write(*,*) 'PART2    = ', PART2
      write(*,*) 'NRN(1)   = ',NRN(1)
      write(*,*) 'NRN(2)   = ',NRN(2)
      write(*,*) 'WMASS    = ',UWMASS
      write(*,*) 'TMASS    = ',UTMASS
      write(*,*) 'HMASS    = ',UHMASS
      write(*,*) 'MSTOP1   = ',UMST1
      write(*,*) 'MSBOY1   = ',UMSB1
      write(*,*) 'YJMIN    = ', YJMIN
      write(*,*) 'YJMAX    = ',YJMAX
      write(*,*) 'PTMIN    = ',PTMIN
      write(*,*) 'EMMIN    = ',EMMIN
      write(*,*) 'IFIT     = ',IFITPDF
      write(*,*) 'ISOFTM   = ',ISOFTM
      write(*,*) 'Q2WWMN   = ',Q2WWMN
      write(*,*) 'Q2WWMX   = ',Q2WWMX
      write(*,*) 'YWWMIN   = ',YWWMIN
      write(*,*) 'YWWMAX   = ',YWWMAX
      write(*,*) 'ZION     = ',ZION
      write(*,*) 'AION     = ',AION
      write(*,*) 'AAANOM = ', AAANOM
      write(*,*) 'DKAPPA = ', D_KAPPA
      write(*,*) 'DLAMBDA = ', LAMBDA
      write(*,*) 'ANOMCUTOFF= ', ANOMCUTOFF
      write(*,*)

c---Initialize model/pdf dependant parameters
      CALL HWMODINI

c---User's initial calculations
      IF(UOUTPUT.NE.0) CALL HWABEG

c---Compute parameter dependent constants
      CALL HWUINC

c---Check POMWIG Settings + Initialisations for consistency
      CALL HWCHEK

c---Call HWUSTA to make any particle stable
      CALL HWUSTA('PI0     ')

c---Initialize elementary process
      CALL HWEINI

c---Initialize event record fixing : this will replace the beam
c   electrons by protons, radiated photons by pomerons/reggeons etc
      CALL HWFXER(.TRUE.,IPROC)
c---Write initialization of lhe file
      CALL LHEINI
c---Loop over events
      DO 100 N=1,MAXEV
c...Initialize event
         CALL HWUINE
c...Generate hard subprocesses
         CALL HWEPRO
c...Include showering and hadronization
         IF (ANSWER.EQ.'Y') THEN
            CALL HWBGEN
            CALL HWDHOB
            CALL HWCFOR
            CALL HWCDEC
            CALL HWDHAD
            CALL HWDHVY
            CALL HWMEVT
         END IF
c...Finish event
         CALL HWUFNE
c...Fix event record (i.e. restore correct intermediate states); print result
         CALL HWFXER(.FALSE.,IPROC)
         IF(N.LE.MAXPR) THEN
           PRINT*, ' '
           PRINT*, ' '
           PRINT*, ' '
           PRINT*, ' '
           PRINT*, 'AFTER EVENT RECORD FIXING:'
           CALL HWUEPR
         ENDIF
         CALL LHEEVT
c...User's event analysis
         IF(UOUTPUT.NE.0) CALL HWANAL
 100  CONTINUE
      CALL LHEEND
c---Terminate elementary process
      CALL HWEFIN

c---User's terminal calculations
      CALL HWAEND
      STOP
      END
C-----------------------------------------------------------------------
C * 07/03/2003, Tibor Kucs
C * User's routine for initialization
C-----------------------------------------------------------------------
      SUBROUTINE HWABEG
      INCLUDE 'HERWIG65.INC'
      INCLUDE 'fpmc.inc'
      INCLUDE 'ffcard.inc'
      INTEGER NWPAWC
      REAL*4 HMEMOR

      RETURN
      END


C----------------------------------------------------------------------
C * 07/11/2003 Maarten B.
C * User's routine to analyze data from event.
C * Distributions relevant for DPE scattering. Assumptions:
C    - initial protons are index 1 and 2
C    - scattered protons are index 5 and 7 (from 1 and 2 resp.)
C    - pomerons are 4 and 6 (radiated from index 1 and 2 resp.)
C    - gluons are 8 and 9 (from index 4 and 6 resp.)
C    - Higgs/Z/... is 10
C    - di-parton is 11 and 12
C-----------------------------------------------------------------------
      SUBROUTINE HWANAL
      INCLUDE 'HERWIG65.INC'
      INCLUDE 'fpmc.inc'
      INCLUDE 'ffcard.inc'

      integer IERR
      REAL V1SQ,V2SQ,XI1,XI2,OXI1,OXI2,WEIGHT,XG1,XG2
      REAL MT1SQ,MT2SQ,PTH,MTHSQ,YH
      REAL PHI1,PHI2,DPHI
      real xg1b, xg2b
      common /remnant/ xg1b,xg2b

      IF(IERROR.NE.0) RETURN

      WEIGHT=1d0

      V1SQ=PHEP(1,5)**2+PHEP(2,5)**2
      V2SQ=PHEP(1,7)**2+PHEP(2,7)**2
      XI1=1.-ABS(PHEP(3,5)/PHEP(3,1))
      XI2=1.-ABS(PHEP(3,7)/PHEP(3,2))
      OXI1=1.-XI1
      OXI2=1.-XI2

      PHI1 = DACOS(PHEP(1,5)/SQRT(V1SQ))
      IF(PHEP(2,5).LT.0d0) PHI1 = -PHI1
      PHI2 = DACOS(PHEP(1,7)/SQRT(V2SQ))
      IF(PHEP(2,7).LT.0d0) PHI2 = -PHI2
      DPHI=PHI2-PHI1
      IF(DPHI.LT.-PI) DPHI=DPHI+2*PI
      IF(DPHI.GT.+PI) DPHI=DPHI-2*PI
      DPHI=ABS(DPHI)*180/3.14159265359

      XG1=ABS(PHEP(4,8)/PHEP(4,4))
      XG2=ABS(PHEP(4,9)/PHEP(4,6))

      xg1b=xg1
      xg2b=xg2

      IF(IERR.NE.0) return

      RETURN
      END
C----------------------------------------------------------------------
C * 07/11/2003 Maarten B.
C * User's routine for terminal calculations, histogram output, etc.
C-----------------------------------------------------------------------
      SUBROUTINE HWAEND
      INCLUDE 'HERWIG65.INC'
      INCLUDE 'fpmc.inc'
      INCLUDE 'ffcard.inc'
      INTEGER ICYCLE

      ICYCLE = 2

      ! print the cross section

      PRINT *, ''
      PRINT *, '========Final summary========='
      PRINT *, 'Cross section[pb]=', 1000.*AVWGT , ' NEVENTS=', UMAXEV

      RETURN
      END

