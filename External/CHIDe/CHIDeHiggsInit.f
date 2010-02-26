C   ============================================================================   C  
C   SUBROUTINE INIT:                                                               C
C   Initialization of parameters and constants                                     C
C   Content: value for gap survival, possibily to change the scales in the         C
C            the Sudakov form factor by changing the value of x, values for cuts,  C
C            parameters for VEGAS,...                                              C
C   =============================================================================  C
      subroutine CHIDeHiggsInit(mH_,mt_,s_,iglu_,x_,xp_,s2_)
      implicit none
      double precision mH_,mt_,s_,x_,xp_,s2_
      integer iglu_ 
      integer ncall1,ncall2,itmx1,itmx2
      integer ncall,itmx,nprn,ndev
      integer grid,ihist,formfac,event,sudaf,scale,evolv,interpol
      integer Vgrid,exp,Model,moode,Lcut,Lenergy,fullvertex,gaps
      integer i,j,iglu
      double precision mh
      double precision xl,xu,acc
      double precision kmax,pi,s,ncolor,nb,Gf,mt,mp
      double precision c,ael,gq
      double precision acut,bcut
      double precision hacutd,lacutd,hbcutd,lbcutd,ycut
      double precision hygap,lygap
      double precision y(50,100),xhis(50,100),count(50)
      double precision ave(50),nbhis(50)
      double precision facwgt
      double precision mmu2,s2
      double precision x,xp
      double precision dgdgrid(1:1001,1:1001),sudagrid(1:1001,1:1001)

      common/bveg1/ncall,itmx,nprn,ndev,xl(11),xu(11),acc
      common/paramint/ncall1,ncall2,itmx1,itmx2   
      common/const/pi,s,kmax,ncolor,gq,nb,Gf,mt,mp
      common/Ifactor/c,ael,iglu
      common/switch/grid,ihist,formfac,event,sudaf,scale,evolv,interpol,
     &Vgrid,exp,Model,moode,Lcut,Lenergy,fullvertex,gaps 
      common/histo/nbhis,y,xhis,count,facwgt,ave
      common/cuts/acut,bcut,hacutd,lacutd,hbcutd,lbcutd,ycut,
     &hygap,lygap
      common/uncert/x,xp
      common/sudaD/mmu2,s2
      common/cinterpol/dgdgrid,sudagrid
      common/higgs/mh


      mh=mH_

      exp=2
      Lenergy=1
      Lcut=3
      Model=2

      fullvertex=2
      formfac=2
      event=2
      iglu=iglu_       
      sudaf=1
      evolv=1
      interpol=1
      scale=2
      gaps=1

C     Constants
      pi=atan(1.d0)*4
      s=s_

      ncolor=3.
C      nb=0.3894*1.E9
      nb=389379.D0
      Gf=1.16639E-5
      mt=mt_
      mp=0.938d0

C    Gap survival
c      if(gaps.EQ.1)then
c       if(exp.EQ.1)then 
c        s2=0.15d0
c       elseif(exp.EQ.2)then
c        s2=0.15d0/2.
c       endif
c      elseif(gaps.EQ.2)then
c        s2=1.d0
c      endif  
      s2 = s2_ 

C   Vectors
      kmax=20.

C   Scale in Sudakov: uncertainties evaluation
      x=x_
      xp=xp_ 

C      x=1.d0
C      xp=1.d0 


      print*, " I=================================I"
      print*, " I  Initialisation of CHIDe model  I"
      print*, " I=================================I" 
      print*, " I  Using parameters:" 
      print*, " I    Higgs mass : ", mh
      print*, " I    Top mass : ", mt
      print*, " I    sqrt(s) = ", sqrt(s)
      print*, " I    iglu = ", iglu
      print*, " I    gap survival = ", s2
      print*, " I    upper sudakov limit = ", x
      print*, " I    lower sudakov limit = ", xp
      print*, " I=================================I" 
      print*, ""
      
C   Cuts for LCWF
C   Typical cuts alpha,beta < 0.1 lost by the p+: 3 quarks so 0.3
      bcut=0.3d0
      acut=0.3d0

C     Cuts TEVATRON and LHC
      if(exp.EQ.1)then
       ycut=2.5                           ! Cut on rapidity: ycut on jets = 2.5
       hygap=5.9                          ! Size of the gap: jets
       lygap=3.6                          !    3.6 < eta(gap) < 5.9
       hacutd=0.08                        ! Higher cut on transfered energy from the p+
       lacutd=0.03                        ! Lower cut on transfered energy from the p+
      elseif(exp.EQ.2)then
       if(Lcut.EQ.1)then
       hacutd=0.02d0
       lacutd=0.002d0
       hbcutd=0.02d0
       lbcutd=0.002d0
       ycut=1.d0
       elseif(Lcut.EQ.2)then
       hacutd=0.018d0
       lacutd=0.005d0
       hbcutd=0.014d0
       lbcutd=0.004d0
       ycut=1.75d0
       ycut=0.06d0
       elseif(Lcut.EQ.3)then              ! id No cuts 
       hacutd=1.d0
       lacutd=0.d0
       hbcutd=1.d0
       lbcutd=0.d0
       ycut=20.d0          
       endif
      endif

C     Parameters from CH paper
      c=-0.41d0
      ael=0.86d0
      
C     Vegas parameters
      ncall1=20000000
      itmx1=5
      ncall2=1000000
      itmx2=10
 

 
C   Initialization of histo
      call hinit('    dsigma')
      call histart(1,0.d0,2.d0,50)             !fort.8 k1
      call histart(4,0.d0,2.d0,50)             !fort.11 k3
      call histart(5,0.d0,2.*pi,20)             !fort.12 theta3
      call histart(6,0.d0,2.d0,50)              !fort.13 k
      call histart(7,0.d0,2.*pi,20)             !fort.14 theta
      call histart(8,0.d0,2.d0,50)             !fort.15 k'
      call histart(9,0.d0,2.*pi,20)             !fort.16 theta'
      call histart(10,0.2d0,1.d0,50)            !fort.17 b1
      call histart(11,0.2d0,1.d0,50)            !fort.18 a3
      call histart(12,0.d0,10.d0,50)       
C      call histart(13,0.d0,10d0,50)           
      call histart(14,-10.d0,10d0,50)           !fort.21 Higgs rapidity   
      call histart(15,-10.d0,10d0,50)           !fort.22 rapidity p1
      call histart(16,-10.d0,10.d0,50)          !fort.23 rapidity p2  


      if(formfac.EQ.2)then
C   Chose of dgdtab as a function of iglu
C   All tabs are made from [interpol.f]
        if(iglu.EQ.1)then
         open(unit = 99, file = "External/CHIDe/Higgsdgdtab1.d",
     &        status = "old")
        elseif(iglu.EQ.2)then
         open(unit = 99, file = "External/CHIDe/Higgsdgdtab2.d", 
     &        status = "old")
        elseif(iglu.EQ.3)then
         open(unit = 99, file = "External/CHIDe/Higgsdgdtab3.d",
     &        status = "old")
        elseif(iglu.EQ.4)then
         open(unit = 99, file = "External/CHIDe/Higgsdgdtab4.d",
     &        status = "old")
        endif
C     Read sudatab 
       open(unit = 98, file = "External/CHIDe/Higgssudatab.d",
     & status = "old")
C     Read the tab and put it in the common
       do i = 0, 1000
        do j = 0, 1000
         read(99, *, end=999) dgdgrid(i+1,j+1)
         read(98,*,end=999) sudagrid(i+1,j+1)
        end do
       end do
 999  continue  

      elseif(formfac.EQ.1)then
C     Read only sudatab 
       open(unit = 98, file = "External/CHIDe/Higgssudatab.d", 
     &      status = "old")
     
C     Read the tab and put it in the common
       do i = 0, 1000
        do j = 0, 1000
        read(98,*,end=111) sudagrid(i+1,j+1)   
       end do
       end do
 111  continue  
      endif
      end
