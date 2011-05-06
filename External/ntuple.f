c////////////////////////////////////////////////
c// @purpose: Simple ntuple
c//
c// @date:    20/02/2011 
c// @author:  O. Kepka (kepkao@fzu.cz)
c////////////////////////////////////////////////


c______________________________________________________________________________ 
c___ create ntuple

      SUBROUTINE NTINIT
c______________________________________________________________________________ 

      REAL hmemor
      COMMON/pawc/hmemor(1000000)
 
***** common block - generator level particles 
       integer ngenmax
       parameter(ngenmax=1000)
       integer ngen
       real px(ngenmax),py(ngenmax),pz(ngenmax)
       real e(ngenmax),m(ngenmax)
       integer id(ngenmax)
       
       common /gener/ngen,
     &     px,py,pz,e,m,id

***** common block - bjorken-x of the parton from pomeron         
       common /remnant/xg1b,xg2b


c***** jets 
      integer njetmax
      parameter(njetmax=30)
      real
     &     pxjet(njetmax),pyjet(njetmax),
     &     pzjet(njetmax),ejet(njetmax)
      common/jets/
     &     njetscel,
     &     pxjet,pyjet,
     &     pzjet,ejet


      INTEGER MODE, ICYCLE, ISTAT


*
      print '(A)',' hropen ...'


      call hbnt(777,'ntuple',' ')
      
      print *,'Creating the ntuple 777'


       call hbname(777,'remnant',xg1b,'xg1,'//
     & 'xg2')

       
     
       call hbname(777,'gener',ngen,'ngen[0,1000]:I,'//
     &     'px(ngen),py(ngen),pz(ngen),e(ngen),rm(ngen),id(ngen)')


      END 
      
c______________________________________________________________________________ 
      
      SUBROUTINE NTEND
c______________________________________________________________________________ 

*---Close the ntuple:
      INTEGER idbg,igener,irad,ifrad
      COMMON /CONST1/ idbg,igener,irad,ifrad

      call hrout(777,icycle,' ')

      END

