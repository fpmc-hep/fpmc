C-----------------------------------------------------------------------
C * 13/03/2008 Oldrich Kepka 
C * Common blocks and declarations needed for the use of Format 
C * Free data cards. They are used in the module.f, module_cms.f, etc.
C----------------------------------------------------------------------- 
      integer UMAXEV, UIPROC, UNFLUX, UNRN1, UNRN2, UIFIT, UISOFTM,
     $ UAAANOM, UZION, UAION,
     $ UCHIDeIGLU
     
      real URMASS,UWMASS,UHMASS,UTMASS,UMST1,UMSB1,UECMS,
     $     UYJMIN,UYJMAX,UPTMIN,UPTMAX,UEMMIN, UEMMAX, UDKAPPA,
     $   UDLAMBDA,
     $   UANOMCUTOFF, UQ2WWMN, UQ2WWMX,UYWWMIN, UYWWMAX, UA0W, UACW,
     $   UA0Z, UACZ,
     $   UCHIDeX, UCHIDeXP, UCHIDeS2 ! R.S.
     
      COMMON /MYFFREAD1/ URMASS,UWMASS,UHMASS,UTMASS,UMST1,UMSB1,UECMS,
     $     UYJMIN, UYJMAX,UPTMIN, UPTMAX, UEMMIN, UEMMAX, UDKAPPA,
     $     UACW, UA0W, UA0Z, UACZ,
     $     UCHIDeX, UCHIDeXP, UCHIDeS2 ! R.S.
     
      COMMON /MYFFREAD2/ UDLAMBDA, UANOMCUTOFF, UYWWMIN, UYWWMAX,
     $ UQ2WWMN, UQ2WWMX

      COMMON /MYFFREAD3/ UMAXEV,UIPROC,UNFLUX,UNRN1,UNRN2,UIFIT, 
     $     UISOFTM, UZION, UAION, UAAANOM,
     $     UCHIDeIGLU ! R.S.

      character*1 UHADR
      COMMON /CC0/ UHADR

      character*3 UTYPEPR, UTYPINT
      COMMON /CC1/ UTYPEPR
      COMMON /CC2/ UTYPINT

      character*4 UPART1, UPART2
      COMMON /CC3/ UPART1
      COMMON /CC4/ UPART2

      character*32 UNTNAME
      COMMON /CYFFLONG1/ UNTNAME


