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
      CALL FPMC_LOAD_PARAMS

c---Events printed
      MAXPR=10

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

