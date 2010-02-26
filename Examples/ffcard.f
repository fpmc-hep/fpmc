C-----------------------------------------------------------------------
C * 10/06/2008 O. Kepka, oldrich.kepka@cern.ch
C * Routine initialize FPMC parameters to default
C * if REACARD = .TRUE. it allows to read setup from external format
C * free card (FFC) with
C * ./module < datacard.ffcard
C----------------------------------------------------------------------- 

      SUBROUTINE FPMC_VAR_INI(READCARD)
      IMPLICIT NONE
      INCLUDE 'ffcard.inc'

      INTEGER SPACE(1000)
      COMMON/CFREAD/SPACE
      LOGICAL READCARD

      CALL FPMC_WELCOME

c---FFC default initialiszation
      UNTNAME    =    'tmpntuple.ntp'
      UMAXEV     =    1000
      UTYPEPR    =    'EXC'
      UTYPINT    =    'QED'
      UHADR      =    'Y'
      UPART1     =    'E+'
      UPART2     =    'E+'
      UIPROC     =    16010
      UHMASS     =    120.
      UNFLUX     =    15
      UNRN1      =    33799
      UNRN2      =    11799
      UWMASS     =    80.425
      UTMASS     =    174.3
      UMST1      =    250.
      UMSB1      =    250.
      UECMS      =    14.E3
      UPTMIN     =    0.
      UPTMAX     =    1.E8
      UYJMIN     =    -6.
      UYJMAX     =    6.
      UEMMIN     =    10.
      UEMMAX     =    1.E8
      UIFIT      =    10
      UISOFTM    =    1
      UZION      =    1
      UAION      =    1
      UAAANOM  =    0       ! original HERWIG AAWW formula
      UDKAPPA    =    0.
      UDLAMBDA   =    0.
      UA0W       =    0.
      UACW       =    0.
      UA0Z       =    0.
      UACZ       =    0.
      UANOMCUTOFF =  -1.   ! if < 0 turns formfactor off
      UQ2WWMN    =    0.
      UQ2WWMX    =    4.
      UYWWMIN    =    0.
      UYWWMAX    =    0.1

C ... begin R.S.
C     CHIDe Model
      UCHIDeIGLU = 4   ! Impact factor parameterisation
      UCHIDeX = 1      ! Scaling upper limit of Sudakov factor integration
      UCHIDeXP = 1     ! Scaling lower limit of Sudakov factor integration
      UCHIDeS2 = 0.075 ! Gap Survival Probability
C ... end R.S.

      IF(READCARD) THEN
         !---Key reading      
         PRINT *, "Reading datacard ..."
         call FFINIT(1000)
         !---JC read config file
         CALL FFSET('SIZE', 32)
         call FFKEY('MAXEV',UMAXEV,1,'integer')
         call FFKEY('TYPEPR',UTYPEPR,3,'mixed')
         call FFKEY('TYPINT',UTYPINT,3,'mixed')
         call FFKEY('HADR',UHADR,1,'mixed')
         call FFKEY('PART1',UPART1,1,'mixed')
         call FFKEY('PART2',UPART2,1,'mixed')
         call FFKEY('ECMS',UECMS,1,'real')
         call FFKEY('IPROC',UIPROC,1,'integer')
         call FFKEY('NFLUX',UNFLUX,1,'integer')
         call FFKEY('NRN1',UNRN1,1,'integer')
         call FFKEY('NRN2',UNRN2,1,'integer')
         call FFKEY('IFIT',UIFIT,1,'integer')
         call FFKEY('IZION',UZION,1,'integer')
         call FFKEY('IAION',UAION,1,'integer')
         call FFKEY('WMASS',UWMASS,1,'real')
         call FFKEY('TMASS',UTMASS,1,'real')
         call FFKEY('HMASS',UHMASS,1,'real')
         call FFKEY('MSTOP1',UMST1,1,'real')
         call FFKEY('MSBOT1',UMSB1,1,'real')
         call FFKEY('YJMAX',UYJMAX,1,'real')
         call FFKEY('YJMIN',UYJMIN,1,'real')
         call FFKEY('PTMIN',UPTMIN,1,'real')
         call FFKEY('PTMAX',UPTMAX,1,'real')
         call FFKEY('EMMIN',UEMMIN,1,'real')
         call FFKEY('EMMAX',UEMMAX,1,'real')
         call FFKEY('Q2WWMN',UQ2WWMN,1,'real')
         call FFKEY('Q2WWMX',UQ2WWMX,1,'real')
         call FFKEY('YWWMIN',UYWWMIN,1,'real')
         call FFKEY('YWWMAX',UYWWMAX,1,'real')
         call FFKEY('ISOFTM',UISOFTM,1,'integer')
         call FFKEY('AAANOM', UAAANOM,1,'integer')
         call FFKEY('DKAPPA', UDKAPPA,1,'real')
         call FFKEY('LAMBDA', UDLAMBDA,1,'real')
         call FFKEY('A0W', UA0W,1,'real')
         call FFKEY('ACW', UACW,1,'real')
         call FFKEY('A0Z', UA0Z,1,'real')
         call FFKEY('ACZ', UACZ,1,'real')
         call FFKEY('ANOMCUTOFF', UANOMCUTOFF,1,'real')

C ... begin R.S.
         call FFKEY('CHIDeIGLU', UCHIDeIGLU, 1, 'integer')
         call FFKEY('CHIDeX',  UCHIDeX, 1, 'real')
         call FFKEY('CHIDeXP', UCHIDeXp, 1, 'real')
         call FFKEY('CHIDeS2', UCHIDeS2, 1, 'real')
C ... end R.S.

         UNTNAME = ''
         call FFKEY('NTNAME',UNTNAME,32,'mixed')
         call FFGO
      ELSE
         PRINT *, "Reading datacard not requested, using default",
     & " parameters."
      ENDIF    

      END SUBROUTINE

      SUBROUTINE FPMC_WELCOME

      WRITE (*,*) ""
      WRITE (*,*) ""
      WRITE (*,*) "     FPMC - Forward Physics Monte Carlo v1.0 beta",
     &   " 12 Jun 2008"
      WRITE (*,*) ""
      WRITE (*,*) "     Authors: M. Boonekamp, O. Kepka, V. Juranek, M.
     &Rangel"
      WRITE (*,*) "              C. Royon"
      WRITE (*,*) "     www.cern.ch/fpmc"
      WRITE (*,*) ""
      WRITE (*,*) ""

      END SUBROUTINE
