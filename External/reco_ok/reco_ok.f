* OK: 1/12/2006 interface to the pxcone jet algo. 
* 1) convert hep to structure of pxcone(pythia typ)
* 2) do the clustering
* 3) fill the appropriate structure in ntuples      

      SUBROUTINE MYRECO(IERR)
      !IERR = 0 OK
      !IERR!= 0 error

      include 'HERWIG65.INC'
      integer nfoundjet

      
C ... PXCONE variables, two pxcone algorithms 
      INTEGER  ITKDM,MXTRK
      ! ... ctyrvektor, max. pocet castic v jednom eventu
      PARAMETER  (ITKDM=4,MXTRK=4000)
      INTEGER  MXJET, MXTRAK, MXPROT
      PARAMETER  (MXJET=20,MXTRAK=4000,MXPROT=500)
      INTEGER  IPASS (MXTRAK),IJMUL (MXJET)
      INTEGER  NTRAK,MODE, NJET, IERR
      DOUBLE PRECISION  PTRAK(ITKDM,MXTRK),PJET(5,MXJET)
      DOUBLE PRECISION  CONER, EPSLON, OVLIM


      
      
c ... particles on generator level     
       integer ngenmax
       parameter(ngenmax=1000)
       integer ngen
       real px(ngenmax),py(ngenmax),pz(ngenmax)
       real e(ngenmax),rm(ngenmax)
       integer id(ngenmax)
       
       common /gener/ngen,
     &     px,py,pz,e,rm,id

c ... jet structure (as in simul)
      integer nrecmax
      parameter(nrecmax=100)
      integer nrec
      real typrec(nrecmax),pxrec(nrecmax),pyrec(nrecmax)
      real pzrec(nrecmax)
      real erec(nrecmax),qrec(nrecmax),eemrec(nrecmax),
     &ehadrec(nrecmax)
      real etrarec(nrecmax),widrec(nrecmax),var6(nrecmax)
      real ctagrec(nrecmax),btag1rec(nrecmax),btag2rec(nrecmax)
      integer ntrarec(nrecmax)
      real sum15ec(nrecmax),sum15hc(nrecmax), sum40(nrecmax)
      real sum40ec(nrecmax)
      
      common /recons1/ nrec,
     & typrec,pxrec,pyrec,pzrec,erec,qrec,eemrec,ehadrec,
     & etrarec,widrec,var6,
     & ctagrec,btag1rec,btag2rec,
     & ntrarec,
     & sum15ec,sum15hc, sum40, sum40ec

c ... missing et
       real circul,transen,transmass,ptmisstx,ptmissty
       real ptmissx,ptmissy,pznu1,pznu2
       
        common /global/ circul,transen,transmass,ptmisstx,ptmissty,
     & ptmissx,ptmissy,pznu1,pznu2


      data nevhep /0/

c ... variables for selecting leptons
      real*8 pt, pttresh
      integer absid
      integer nleptons
      INTEGER MXLEP
      PARAMETER (MXLEP=20)
      REAL*8 PLEP(5,MXLEP)
c ... local variables for misset
      real*8 misspx, misspy


c ... loop variables
      INTEGER I, J, N, IPART
      DOUBLE PRECISION RAP

      integer idiv


C...HEPEVT commonblock.
c      PARAMETER (NMXHEP=4000)
c      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
c     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
c      DOUBLE PRECISION PHEP,VHEP
c      SAVE /HEPEVT/
      
c------------------------------------------------------------------- 

C ... PXCONE parameters ...
      MODE    = 2        !Snow mass scheme
      CONER   = 0.7
      EPSLON  = 4
c     EPSLON  = 8
      OVLIM   = 0.5d0

c------------------------------------------------------------------- 
c particles on the generator level
      call vzero(px,ngenmax)
      call vzero(py,ngenmax)
      call vzero(pz,ngenmax)
      call vzero(e,ngenmax)
      call vzero(rm,ngenmax)
      call vzero(id,ngenmax)
        
        
      N=NHEP
      IPART=0
      do 1515 I=1,N
          IF(ISTHEP(I).EQ.1) THEN
             IPART=IPART+1
             px(IPART)=sngl(PHEP(1,I))
             py(IPART)=sngl(PHEP(2,I))
             pz(IPART)=sngl(PHEP(3,I))
             e(IPART) =sngl(PHEP(4,I))
             rm(IPART)=sngl(PHEP(5,I))
             id(IPART)=IDHEP(I)
c         rm(i)=sngl(PTRAK(5,i))
c          print '(A,4F8.2)','gen:',px(i),py(i),pz(i),e(i)
c         print '(A,F8.2,I6)','rm id :',p(i,5),k(i,2)
          ENDIF
1515  continue
      ngen=IPART
          
c-------------------------------------------------------------------          
      NTRAK=0
      N=NHEP
c     print '(A,I6)', 'NHEP', NHEP
      DO 180 I=1,N
          !if particle in the final state
          IF(ISTHEP(I).EQ.1) THEN
             RAP = DABS( 0.5*DLOG( (PHEP(4,I)+PHEP(3,I))/
     .                (PHEP(4,I)-PHEP(3,I)) ) )
c             print*,rap
             IF(RAP.LT.4) THEN
                NTRAK=NTRAK+1
                DO 160 J=1,5
                  PTRAK(J,NTRAK)=PHEP(J,I)
 160            CONTINUE
             ENDIF     
               
c          print '(A,4F8.2)','phep:',PHEP(1,i),PHEP(2,i),
c     .                              PHEP(3,i),PHEP(4,i)
               
          ENDIF     
  180 CONTINUE       


c-------------------------------------------------------------------           
      !clustering
      
      NJET=0                          !midpoint alg.
      CALL PXCONE(MODE,NTRAK,ITKDM,PTRAK,CONER,EPSLON,OVLIM,MXJET,
     +            NJET,PJET,IPASS,IJMUL,IERR)

      IF(IERR.NE.0) THEN
         print *, 'PXCONE did not converge'
         RETURN
      ENDIF   



      idiv = MOD(nevhep, 1000)
      if(nevhep.lt.100) then 
        print *,' *** cone_interface nevhep :', nevhep
      else  if ( idiv.EQ.0 ) then
            print *,' ***  cone_interface nevhep :', nevhep
      ENDIF


c...jets common blocks
      call vzero(typrec,nrecmax)
      call vzero(pxrec,nrecmax)
      call vzero(pyrec,nrecmax)
      call vzero(pzrec,nrecmax)
      call vzero(erec,nrecmax)
      call vzero(qrec,nrecmax)
      call vzero(eemrec,nrecmax)
      call vzero(ehadrec,nrecmax)
      call vzero(etrarec,nrecmax)
      call vzero(widrec,nrecmax)
      call vzero(ctagrec,nrecmax)
      call vzero(btag1rec,nrecmax)
      call vzero(btag2rec,nrecmax)
      call vzero(ntrarec,nrecmax)
      call vzero(sum15ec,nrecmax)
      call vzero(sum15hc,nrecmax)
      call vzero(sum40,nrecmax)
      call vzero(sum40ec,nrecmax)

      nrec=0


c --- reconstruction of the leptons
c will only copy all leptons above 
      pttresh=5d0
      nleptons=0

      misspx = 0
      misspy = 0
      
      do 1600 i=1, ngen
         absid=abs(id(i))
         if( (absid.eq.11).or.(absid.eq.13).or.(absid.eq.15) ) then
            ! it is lepton

            pt= sqrt( px(i)*px(i) + py(i)*py(i) )
            nleptons=nleptons+1

            if(pt.gt.pttresh) then
               nrec = nrec+1
               ! identity - charge not distingished
               ! save leptons
c               PLEP(1, nleptons) = px(i)
c               PLEP(2, nleptons) = py(i)
c               PLEP(3, nleptons) = pz(i)
c               PLEP(4, nleptons) = e(i)
c               PLEP(5, nleptons) = rm(i)

               ! fill reconstructed object block
               pxrec(nrec) = px(i)
               pyrec(nrec) = py(i)
               pzrec(nrec) = pz(i)
               erec(nrec)  = e(i)

               ! identify type of leptons
               if(abs( id(i) ).eq.11 ) then 
                  typrec(nrec) = 1
               else if (abs( id(i)).eq.13 ) then
                   typrec(nrec) = 2
               else if (abs( id(i)).eq.15 ) then 
                   typrec(nrec) = 3
               end if

             endif
          endif   

          ! missing Et
          if( (absid.eq.12).or.(absid.eq.14).or.(absid.eq.16) ) then
             misspx =  misspx + px(i)
             misspy =  misspy + py(i)
          endif   
          
1600  continue
c      print '(A,I6)', 'nleptons', nleptons

c --- look for misidentified leptons as jets

      !assumtion all jet hadronic
c      print '(A,I6)', 'njet', njet
c      do 1550 i=1, njet
c         typrec(i)= 4
c         eemrec(i)= 0
c         pxrec(i) = PJET(1,i)
c         pyrec(i) = PJET(2,i)
c         pzrec(i) = PJET(3,i)
c        erec(i)  = PJET(4,i)
cc         print '(A, 4F8.2)', 'pjet x y z e :', PJET(1,i),PJET(2,i),
c     .                                         PJET(3,i),PJET(4,i)
c         print *,''
c 1550 continue         
 
      !assumtion all jet hadronic
c      print '(A,I6)', 'njet', njet
      do 1550 i=1, njet
         nrec=nrec+1
         typrec(nrec)= 4
         eemrec(nrec)= 0
         pxrec(nrec) = PJET(1,i)
         pyrec(nrec) = PJET(2,i)
         pzrec(nrec) = PJET(3,i)
         erec(nrec)  = PJET(4,i)
c         print '(A, 4F8.2)', 'pjet x y z e :', PJET(1,i),PJET(2,i),
c     .                                         PJET(3,i),PJET(4,i)
c         print *,''
 1550 continue         


      ! assing the missing Et in the common block
      ptmissx = misspx
      ptmissy = misspy

c --- for debug use

c      do 1700 i=1, nrec
c            print '(A,4F8.2)','rec:',pxrec(i),pyrec(i),pzrec(i),erec(i)
c            print '(A,F8.2)','typrec :', typrec(i)
c1700  continue            
      
       
       RETURN
       END

      

