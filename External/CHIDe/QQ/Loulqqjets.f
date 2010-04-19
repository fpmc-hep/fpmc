C   ============================================================================   C
c   ----------->       Mofification to qq bar production    <-------------------   C
C   Jets production from Igor                                                      C
C   ifort -O2 -pc80 -mp1 -parallel -fpe:3 -xP -vec-report:0 -par-report:0          C
C   Louljet.f Histo.f vegas11.f dgd2007.f dgd2008.f dgdforward.f dcadredo.f        C
C   makeLHE.f                                                                      C
C   ============================================================================   C
C                                                                                  C
C   *Louljet.f    Evaluate the cross section of dijet central exclusive production C
C                 via 2-gluons exchange in nb.                                     C
C   *Histo.f      Histogram different variables and need ihist=1.                  C
C   *vegas11.f    Modified version of VEGAS used to evaluate multidimensionnal     C
C                 phase space integrals.                                           C
C   *dgd2007.f    2007 version of the unintegrated gluon densities. This factor    C
C                 doesn't go to zero when incoming momenta become soft.            C
C   *dgd2008.f    2008 version of the unintegrated gluon densities. This factor    C
C                 goes to zero when incoming momenta become soft.                  C 
C   *dgdforward.f Unintegrated gluon densities in the forward direction. Skewness  C
C                 is taken into account by x -> 0.41*x and we used a prefactor     C
C                 that introduced the non-zero transverse momentum transfer.       C
C                 Dgdforward need the file 'grv98lo.grid'                          C
C   *dcadredo.f   A 1d integral program used to evaluate the Sudakov form factor   C
C   *Durhamlike.f A external subroutine that mimic Durham calculation with 2       C
C                 possibilities: the reference curve with the approx ki << k or    C
C                 'Real Durham' with a different Sudakov form factor and the same  C
C                 approx ki << k                                                   C 
C   ----------------------------------------------------------------------------   C
C   !The integration is very weighty for the Cpu!                                  C
C   1. Speed the integration using interpolation for dgdforward and the Durham     C
C   Sudakov form factor with running as -> need dgdtab1.d, dgdtab2.d, dgdtab3.d,   C
C   dgdtab4.d for each value of iglu in dgdforward and sudatab.d.                  C
C   2. Help Vegas using preformated grid (cf: switch Vgrid) -> need grid.j with    C
C   j=[60,67] and LHCgrid.i with i=[70,78].                                        C
C   ----------------------------------------------------------------------------   C
C   STRUCTURE and POSSIBILITIES                                                    C
C    *Main part: Included the call to VEGAS and the do loop on minimum             C
C              transverse energy.                                                  C
C    *Subroutine init: Constant value, opening of needed files, initialisation of  C
C              histograms and switches between possibilities                       C 
C          - Experiment: exp=1,2 or 3 for TEVatron Run I, Run II or LHC            C
C              Also choose the cuts applied on physical quantities                 C 
C          - Use of a new or a preformated grid: Vgrid=1,2                         C
C          - Ingredients: Model=1,2 or 3 Choose all pieces, use the reference      C
C              pieces, use Durham like ingredients i.e. k<<ki                      C
C              If Model=1, you have to choose:                                     C
C              + A perturbative calculation or not (cuts on the gluons propagators)C
C              + A impact factor and its options                                   C
C              + A Sudakov form factor, its scales and if alpha_s is running       C
C    *Subroutine interinit: open and read dgdtab and sudatab for interpolation     C  
C    *Subroutine output: format the output                                         C
C              fort.2 -> cross section as a function of minimum transverse energy  C
C              fort.30 -> Mean values of the variables of integration from Histo.f C
C    *Function dsigma: compute the cross section from the function integ           C
C              Content: kinematic definition, cuts, rapidity definitions, phase    C
C              space factor, regge factor, gluon/gluon mandelstam variables        C
C              and make the histogram.                                             C
C    *Function integ: analytic part of the integrand                               C
C              Content: integrand, different defintions for the impact factor and  C
C              for the Sudakov form factor, definition of the coupling constants   C
C              One Integ function for each diagrams                                C
C    *Subroutine and function: several functions and definitions                   C
C          - Alpha strong                                                          C
C          - LCWF impact factor, elastic and diffractif                            C
C          - UgD forward and prefactor                                             C
C          - Sudakov form factor                                                   C
C          - Interpolation bilinear for UgD (dgdpol) and Sudakov (sudapol)         C
C          - Addition, substraction and multiplication of 4-vectors                C
C    ---------------------------------------------------------------------------   C
C    Made by Alice Dechambre                                           June 2008   C 
C    ===========================================================================   C  


C   ============================================================================   C
C   MAIN PART:                                                                     C
C   call initialisation, do loop on minimum transverse energy                      C   
C   Call VEGAS in 3 different ways depending on the experiment wanted:             C
C    TEVatron Run II, LHC, TEVatron Run I                                          C
C    Each call to VEGAS is divided in one part producing the grid from the         C
C    begining (Vgrid=1) and the other part used a preformated grid (Vgrid=2)       C
C   Call output                                                                    C
C   ============================================================================   C
      program louljets
	implicit none
      integer ncall,itmx,nprn,ndev,it,ndo
      integer ndmx,mds
      integer ncall1,ncall2,itmx1,itmx2,ihist
      integer event,exp,sudaf,grid,formfac,scale,splash
      integer Model,evolv,Vgrid,Lcut,interpol,moode,lhe,nscale
      integer iglu,j,i
      integer update
      double precision Etjetmin,Etjet0,Etjetmax,step
      double precision dsigma,chi2a
      double precision xl,xu,acc,si,swgt,schi,xi
      double precision alph
      double precision avgi,sd,weight,c,ael
      double precision Etcut,hygap,lygap,hbcutd,lbcutd,ycut
      double precision ylcut,hacutd,lacutd,yclucut,propcut,yyj
      double precision Etmin,as
      double precision Etjetveto,Rjj,Mjj,Mx,Rt,Mgg,Mj
      
      double precision mm,mq
       
      external dsigma

      common/bveg1/ncall,itmx,nprn,ndev,xl(11),xu(11),acc
      common/bveg2/it,ndo,si,swgt,schi,xi(50,11)
      common/bveg3/alph,ndmx,mds
      common/paramint/ncall1,ncall2,itmx1,itmx2,ihist
      common/point/Etjetmin,Etjet0,Etjetmax,step
      common/out/avgi,sd,weight,update
      common/switch/event,exp,sudaf,grid,formfac,
     &scale,Model,evolv,Vgrid,Lcut,interpol,moode,lhe,splash,nscale
      common/cut/Etcut,hygap,lygap,hbcutd,lbcutd,ycut,
     &ylcut,hacutd,lacutd,Etmin,yclucut,propcut,yyj
      common/Ifactor/c,ael,iglu
      common/algo/Etjetveto,Rjj,Mjj,Mx,Rt,Mgg,Mj
      common/mquark/mq

      call init
      call interinit

      update=0

C     grid.j: Stored grid file for Vegas and TEVatron RunII
C     LHCgrid.i: Stored grid file for Vegas and LHC
       j=60
       i=70

C     Variation of the minimum transverse energy
C       do mm=100.,200.,10.
       mm=120.
       mq=mm/2.   
C       do Etjetmin=Etjet0,Etjetmax,step
C       do Etjetmin=10.,30.,20.
        if(exp.EQ.2)Etjetmin=10.
        if(exp.GE.3)Etjetmin=50.

C     Splashout prescription
      if(splash.EQ.1)then
       Etcut=Etjetmin
       Etjetveto=10.
       print*,'Jets Etmin=',Etjetveto
      elseif(splash.EQ.2)then
       Etcut=Etjetmin/0.75d0
      elseif(splash.EQ.3)then
       Etcut=Etjetmin/0.8d0
      elseif(splash.EQ.4)then 
       Etcut=(as(Etjetmin**2)/2.+1.)*Etjetmin+1.
      endif
      
      print*,'Emin=',Etcut,' GeV'
      print*,'Mqq=',mm,' GeV'

C   ------------------------------ CDF RUN II ---------------------------------   C
      if(exp.EQ.2)then

C     New Grid from Vegas
      if(Vgrid.EQ.1)then
       nprn=0
       ncall=ncall1  
       itmx=itmx1 
       ihist=0
       if(Etjetmin.EQ.5.)then
         call VEGAS(11,dsigma,avgi,sd,chi2a)
        ncall=ncall2
        itmx=itmx2
        ihist=1
         call VEGAS1(11,dsigma,avgi,sd,chi2a)
       elseif(Etjetmin.GT.5.)then
        ncall=ncall2
        itmx=itmx2
        ihist=0
         call VEGAS1(11,dsigma,avgi,sd,chi2a)
       endif

C      Used the stored grid from Vegas in fort.(60+j)
       elseif(Vgrid.EQ.2)then
        read(j,*)ndo,xi
        nprn=0
        ncall=ncall2
        itmx=itmx2
        ihist=1
         call VEGAS1(11,dsigma,avgi,sd,chi2a)
        j=j+1
        if(it.EQ.1)write(82,*)ndo,xi
       endif
C   ---------------------------------- LHC -------------------------------------   C
      elseif(exp.GE.3)then

C     New Grid from Vegas
      if(Vgrid.EQ.1)then  
       nprn=0
       ncall=ncall1  
       itmx=itmx1 
       ihist=0
       if(Etjetmin.EQ.20.)then
         call VEGAS(11,dsigma,avgi,sd,chi2a)
        ncall=ncall2
        itmx=itmx2
        ihist=0
         call VEGAS1(11,dsigma,avgi,sd,chi2a)
       elseif(Etjetmin.GT.20)then
        ncall=ncall2
        itmx=itmx2
        ihist=0
         call VEGAS1(11,dsigma,avgi,sd,chi2a)
       endif

C      Used the stored grid from Vegas in fort.(70+i)
       elseif(Vgrid.EQ.2)then
        read(i,*)ndo,xi
        nprn=0
        ncall=ncall2
        itmx=itmx2
        ihist=1
         call VEGAS1(11,dsigma,avgi,sd,chi2a)
        i=i+1
       endif
C   -------------------------------- CDF RUN I -------------------------------   C
      elseif(exp.EQ.1)then
C     RunI: 3 initialisation of the grid
       grid=1
       nprn=0
       ncall=ncall1  
       itmx=itmx1 
       ihist=0
       if(Etmin.EQ.25.)then
         call VEGAS(11,dsigma,avgi,sd,chi2a)
        ncall=ncall2
        itmx=itmx2
        ihist=0
        grid=2
         call VEGAS1(11,dsigma,avgi,sd,chi2a)
       elseif(Etmin.GT.25.)then
        grid=2  
        ncall=ncall2
        itmx=itmx2
        ihist=0
         call VEGAS1(11,dsigma,avgi,sd,chi2a)
       endif  
C   ----------------------------------------------------------------------------   C 
      endif
      call output
C      end do
      end 


C   ============================================================================   C  
C   SUBROUTINE INIT:                                                               C
C   Initialization of parameters and constants                                     C
C   Content: value for gap survival, possibily to change the scales in the         C
C            the Sudakov form factor by changing the value of x, values for cuts,  C
C            parameters for VEGAS,...                                              C
C   =============================================================================  C
      subroutine init
	implicit none
      integer ncall1,ncall2,itmx1,itmx2,ihist
      integer ncall,itmx,nprn,ndev
      integer event,exp,sudaf,grid,formfac,scale,Vgrid,Lcut
      integer Model,evolv,interpol,moode,lhe,splash,nscale
      integer iglu
      integer pout1,pout2,pout3,pout4,sout1,sout2,sout3,sout4
      double precision xl,xu,acc 
      double precision pi,s,ncolor,gg,gq,kmax,nb,mp
      double precision y(50,100),xhis(50,100),count(50)
      double precision ave(50),nbhis(50)
      double precision facwgt,ael
      double precision etcut,hygap,lygap,hbcutd,lbcutd,ycut,yclucut
      double precision ylcut,hacutd,lacutd,Etmin,propcut,yyj
      double precision Etjetmin,Etjet0,Etjetmax,step
      double precision c
      double precision mmu2,s2
      double precision xmax,amax,amin
      double precision x,xp

      common/bveg1/ncall,itmx,nprn,ndev,xl(11),xu(11),acc
      common/paramint/ncall1,ncall2,itmx1,itmx2,ihist 
      common/point/Etjetmin,Etjet0,Etjetmax,step
      common/histo/nbhis,y,xhis,count,facwgt,ave
      common/const/pi,s,ncolor,gg,gq,kmax,nb,mp
      common/cut/Etcut,hygap,lygap,hbcutd,lbcutd,ycut,
     &ylcut,hacutd,lacutd,Etmin,yclucut,propcut,yyj
      common/switch/event,exp,sudaf,grid,formfac,
     &scale,Model,evolv,Vgrid,Lcut,interpol,moode,lhe,splash,nscale
      common/Ifactor/c,ael,iglu
      common/sudaD/mmu2,s2
      common/BotE/amax,xmax,amin
      common/uncert/x,xp
      common/logs/pout1,pout2,pout3,pout4,sout1,sout2,sout3,sout4

      print*,'Experiment and cuts?'
      print*,'      CDF runI=1, CDF runII=2,' 
      print*,'      LHC 14 TeV=3, LHC 10 TeV=4, LHC 7 TeV=5'
      read*,exp
      if(exp.GE.3)then
       print*,'LHC set of cut?'
       print*,'     FP420+FP220=1, LHC Higgs=2'
       read*,Lcut
      endif
      print*,'Model?'
      print*,'     Manual=1,'
      print*,'     Automatic: Reference=2, Durham like=3, Forshaw=4'
      read*,Model
      print*,'Splash-Out: run Pythia or use a parametrization?'
      print*,'      Pythia=1, Parametrization=2'
      read*,lhe
      if(lhe.EQ.1)then
       splash=1
      elseif(lhe.EQ.2)then
       print*,'Splash-out, Etpartons=?'
       print*,'        Etjet=1, Etjet/0.75=2, Etjet/0.8=3, Etjet~(as)=4'
       read*,splash 
      endif   
         
C   --------------------------- All Ingredients -------------------------------   C
      if(Model.EQ.1)then
       print*,'                                                        '  
       print*,'--------------------------------------------------------'
       print*,'-> Please, choose your ingredients'
       print*,'--------------------------------------------------------'
       print*,'                                                        '

       print*,'New grid for Vegas or saved grid?'
       print*,'      New=1, Old=2'
       read*,Vgrid
       print*,'Pertubative propagators, propagators cut?'
       print*,'     No cut=-10., 1GeV= 1.'
       read*,propcut
       print*,'Form Factor?'
       print*,'     CH Form Factor=1,' 
       print*,'     Igor Unintegrated Gluon Distribution=2'
       read*,formfac
       if(formfac.EQ.1)then
         print*,'     Elastic=1, Diffractif=2, SD with proton split=3'
       elseif(formfac.EQ.2)then
         print*,'     dgd2007=1, dgd2008=2, dgdforward=3'
       endif
       read*,event
       if(event.EQ.3)then
         print*,'Iglu?'
         print*,'     GRV(2000)=1, Matching(2002)=2, GRV(2005)=3,'
         print*,'     Smaller NPregion=4'
         read*,iglu
       endif
       print*,'Sudakov form factor?'
       print*,'      Durham=1, DDT=2, DLA=3, Exhume=4,' 
       print*,'      CDHI=5, without=7'
       read*,sudaf
       if(sudaf.LE.3.OR.sudaf.EQ.4)then
         print*,'Running alpha_s?'
         print*,'      Running=1, Fixed=2'
         read*,evolv
         if(evolv.EQ.1)then
          interpol=1
         elseif(evolv.EQ.2)then
          interpol=2
         endif 
         print*,'Scale in Sudakov?'
         print*,'      k^2=1, (k+ki)^2=2'
         read*,scale
         print*,'Number of scales in Sudakov?'
         print*,'      One scale=1, Two scales=2'
         read*,nscale
         if(nscale.EQ.2)interpol=2 
       endif
C   -------------------------- Automatic Model ---------------------------------   C
      elseif(Model.EQ.2)then
       Vgrid=2   
       propcut=-10.d0
       formfac=2
       event=3
       iglu=4
       sudaf=1
       evolv=1
       interpol=1
       scale=2
       nscale=1
      elseif(Model.EQ.3)then
       Vgrid=2   
       propcut=-10.
       print*,'k>>k_i or Durham like?'
       print*,'      k>>k_i=1, Durham like=2'
       read*,moode
       formfac=2
       event=3
       iglu=4
       sudaf=1
       evolv=1
       interpol=1
       nscale=1
       if(nscale.EQ.2)interpol=2
      elseif(Model.EQ.4)then
       Vgrid=2   
       propcut=-10.
       formfac=2
       event=3
       iglu=4
       sudaf=4
       evolv=1
       interpol=1
       scale=2
       nscale=1
      endif
C   ---------------------------------------------------------------------------   C
C     Constants
      pi=atan(1.d0)*4
      if(exp.EQ.2)then  
      s=1960.**2
      elseif(exp.EQ.1)then
      s=1800.**2
      elseif(exp.EQ.3)then
      s=14000.**2
      elseif(exp.EQ.4)then
      s=10000.**2 
      elseif(exp.EQ.5)then
      s=7000.**2 
      endif     

      nb=0.3894*1.E6
      mp=0.938d0
      ncolor=3.
      pi=atan(1.d0)*4

C    Gap survival
      if(exp.EQ.1.or.exp.EQ.2)then
       s2=0.15d0
      elseif(exp.GE.3)then
       s2=0.15/2.d0
      endif 

C     Alpha soft
      if(formfac.EQ.1)then
       c=0.61d0
       ael=0.88d0
      elseif(formfac.EQ.2)then
       ael=1.d0
      endif

C   Vectors
      kmax=50.

C   Scale in Sudakov: uncertainties evaluation
       x=1./2.d0
       xp=1./2.d0 

C   New Reference: 2010 results 
C       x=1./2.d0
C       xp=1./2.d0
C   Old Reference: pre-2009 results 
C       x=2.d0
C       xp=1.d0

C       if(sudaf.EQ.4)then
C        x=1.d0
C       endif
 
C   Cuts run I and run II
      if(exp.EQ.1)then
      ycut=2.4
      ylcut=-4.2
      hygap=5.9
      lygap=2.4
      hacutd=0.095
      lacutd=0.035
      hbcutd=0.03
      lbcutd=0.01
      elseif(exp.EQ.2)then
       ycut=2.5                           ! Cut on rapidity: ycut on jets = 2.5
       hygap=5.9                          ! Size of the gap: jets
       lygap=3.6                          !    3.6 < eta(gap) < 5.9
       hacutd=0.08                        ! Higher cut on transfered energy from the p+
       lacutd=0.03                        ! Lower cut on transfered energy from the p+
       yyj=-0.5d0
      elseif(exp.GE.3)then
       if(Lcut.EQ.1)then                  !FP420+FP220
       hacutd=0.2
       lacutd=0.002
       hbcutd=0.2
       lbcutd=0.002
       ycut=1.
       elseif(Lcut.EQ.2)then
       hacutd=0.018
       lacutd=0.005
       hbcutd=0.014
       lbcutd=0.004
       ycut=1.75
       yclucut=0.06
       endif
      endif

C     Parameters of integration for Vegas
      ncall1=500000
      itmx1=10
      ncall2=2000000
      itmx2=10

C   Counter
      pout1=0
      pout2=0
      pout3=0
      pout4=0
      sout1=0
      sout2=0
      sout3=0
      sout4=0

C     Points for the do loop
      if(exp.EQ.2)then
       Etjet0=5.
       Etjetmax=40.
       step=5.
      elseif(exp.GE.3)then
       Etjet0=20.
       Etjetmax=100.
       step=10.
      endif 

C     Open the grid for VEGAS
C     Unit 60 for TEVATRON and Unit 70 for LHC
      if(Vgrid.EQ.2)then
       if(exp.EQ.2)then
        open(unit = 60, file = "grid.60", status = "old")
        open(unit = 61, file = "grid.61", status = "old")
        open(unit = 62, file = "grid.62", status = "old")
        open(unit = 63, file = "grid.63", status = "old")
        open(unit = 64, file = "grid.64", status = "old")
        open(unit = 65, file = "grid.65", status = "old")
        open(unit = 66, file = "grid.66", status = "old")
        open(unit = 67, file = "grid.67", status = "old")
       elseif(exp.GE.3)then
        open(unit = 70, file = "LHCgrid.70", status = "old")
        open(unit = 71, file = "LHCgrid.71", status = "old")
        open(unit = 72, file = "LHCgrid.72", status = "old")
        open(unit = 73, file = "LHCgrid.73", status = "old")
        open(unit = 74, file = "LHCgrid.74", status = "old")
        open(unit = 75, file = "LHCgrid.75", status = "old")
        open(unit = 76, file = "LHCgrid.76", status = "old")
        open(unit = 77, file = "LHCgrid.77", status = "old")
        open(unit = 78, file = "LHCgrid.78", status = "old")
        open(unit = 79, file = "LHCgrid.79", status = "old")
       endif
      endif
 
C     Initialization of histo
      call hinit('    dsigma')
      call histart(1,0.d0,25.d0,50)                  !fort.8 k1²
      call histart(2,0.d0,25.d0,50)                  !fort.9 k2²
      call histart(3,0.d0,2.*pi,10)                       !fort.10 theta1
      call histart(4,0.d0,25.d0,50)                  !fort.11 k3²
      call histart(5,0.d0,2.*pi,10)                       !fort.12 theta3
      call histart(6,0.d0,5.d0,50)                  !fort.13 k²
      call histart(7,0.d0,2.*pi,20)                       !fort.14 theta
      call histart(8,0.d0,25.d0,50)                       !fort.15 k'²
      call histart(9,0.d0,2.*pi,20)                       !fort.16 theta'
      call histart(10,0.d0,0.2d0,50)                       !fort.17 b1
      call histart(11,0.d0,0.2d0,50)                       !fort.18 b1+b2
      call histart(12,0.d0,0.2d0,50)                       !fort.19 a1
      call histart(13,0.d0,0.2d0,50)                       !fort.20 a1+a2
      call histart(14,-4.d0,4.d0,50)                    !fort.21 y(jet1)
      call histart(15,-4.d0,4.d0,50)                    !fort.22 y(jet2)
      call histart(16,-16.d0,16.d0,50)                    !fort.23 y(p1)  
      call histart(17,-16.d0,16.d0,50)                    !fort.24 y(p2)
      call histart(18,-10000.d0,10000.d0,50)              !fort.25 tgg
      call histart(19,-10000.d0,10000.d0,50)              !fort.26 ugg
      call histart(20,-10000.d0,10000.d0,50)              !fort.27 sgg
      call histart(21,0.d0,2.*kmax,20)                 !fort.28 (k+k1)²
      call histart(22,0.d0,2.*kmax,20)                 !fort.29 (k+k3)²
      call histart(24,0.d0,2.*kmax,20)                 !fort.31 (k+k2)²
      call histart(25,0.d0,2.*kmax,20)                 !fort.32 (k+k2-k3)²
      call histart(26,-3.d0,3.d0,50)                      !fort.33 (yj1+yj2)/2
      call histart(27,0.d0,60.d0,50)                       !fort.34 (Etj1+Etj2)/2
      call histart(28,0.d0,5.d0,50)                       !fort.35 Rjj=Mjj/Mgg
      call histart(29,0.d0,1.d0,50)                       !fort.36 x=b1/b2
      call histart(30,20.d0,180.d0,50)                    !fort.37 Mjj
      call histart(31,20.d0,180.d0,50)                    !fort.38 Mx
      call histart(32,0.d0,1.d0,50)                       !fort.39 Rt=Mj/Mgg
      call histart(33,20.d0,180.d0,50)                    !fort.40 Mj
      call histart(34,20.d0,180.d0,50)                    !fort.41 Mgg
      call histart(35,0.d0,8.d0,50)                       !fort.42 nj = number of jets
      call histart(36,-3.d0,3.d0,50)                    !fort.43 y(jet1)+y(jet2)/2
      call histart(37,-4.d0,4.d0,50)                    !fort.44 y_cluster
      call histart(38,-3.d0,3.d0,50)                    !fort.45 y of leading jet
      call histart(39,-3.d0,3.d0,50)                    !fort.46 y of second jet
      end

C   ============================================================================   C
C   SUBROUTINE INTERINIT:                                                          C  
C   Open and read dgdtab and sudatab for interpolation                             C
C   ============================================================================   C
      subroutine interinit
        implicit none
      integer i,j,iglu
      integer event,exp,sudaf,grid,formfac,splash,nscale
      integer scale,Model,evolv,Vgrid,Lcut,interpol,moode,lhe
      double precision c,ael
      double precision dgdgrid(1:1001,1:1001),sudagrid(1:1001,1:1001)

      common/switch/event,exp,sudaf,grid,formfac,
     &scale,Model,evolv,Vgrid,Lcut,interpol,moode,lhe,splash,nscale
      common/Ifactor/c,ael,iglu
      common/cinterpol/dgdgrid,sudagrid

      if(formfac.EQ.2)then
C     Chose of dgdtab as a function of iglu
       if(iglu.EQ.1)then
        open(unit = 99, file = "dgdtab1.d", status = "old")
       elseif(iglu.EQ.2)then
        open(unit = 99, file = "dgdtab2.d", status = "old")
       elseif(iglu.EQ.3)then
        open(unit = 99, file = "dgdtab3.d", status = "old") 
       elseif(iglu.EQ.4)then
        open(unit = 99, file = "dgdtab4.d", status = "old")
       endif

C     Open sudatab
       open(unit = 98, file = "sudatab.d", status = "old")
C    Durham prescription for Delta=qt/(qt+mu)
C    Also used in the before 2010 results i.e. before the Forshaw paper
       if(Model.EQ.3.AND.moode.EQ.2)then
        open(unit = 98, file = "DurhamSudatab.d", status = "old")
       endif

  
C     Read the tab and put it in the common
       do i = 0, 1000
        do j = 0, 1000
         read(99, *) dgdgrid(i+1,j+1)
         read(98,*) sudagrid(i+1,j+1)
        end do
       end do
      continue

      elseif(formfac.EQ.1)then
C     Read only sudatab 
       open(unit = 98, file = "sudatab.d", status = "old")
       do i = 0, 1000
        do j = 0, 1000
        read(98,*,end=111) sudagrid(i+1,j+1)   
       end do
       end do
 111  continue  
      endif
      end
  
C   ============================================================================   C
C   SUBROUTINE OUTPUT:                                                             C
C   format the output                                                              C
C   fort.2 -> cross section as a function of minimum transverse energy             C
C   fort.30 -> Mean values of the variables of integration from Histo.f            C  
C   ============================================================================   C
      subroutine output
        implicit none
      integer ncall,itmx,nprn,ndev,it,ndo
      integer ncall1,ncall2,itmx1,itmx2,ihist
      integer event,exp,sudaf,grid,formfac,scale,splash,nscale
      integer iglu,Model,evolv,Vgrid,Lcut,interpol,moode,lhe
      integer update
      integer pout1,pout2,pout3,pout4,sout1,sout2,sout3,sout4
      double precision avgi,sd,weight
      double precision pi,s,ncolor,gg,gq,kmax,nb,mp
      double precision xl,xu,acc,si,swgt,schi,xi
      double precision y(50,100),xhis(50,100),count(50)
      double precision ave(50),nbhis(50)
      double precision facwgt
      double precision etcut,hygap,lygap,hbcutd,lbcutd,ycut,yclucut
      double precision ylcut,hacutd,lacutd,Etmin,propcut,yyj
      double precision c,ael
      double precision x,xp
      double precision bmax,xmax,bmin
      double precision Etjetmin,Etjet0,Etjetmax,step

      common/bveg1/ncall,itmx,nprn,ndev,xl(11),xu(11),acc
      common/bveg2/it,ndo,si,swgt,schi,xi(50,11)
      common/paramint/ncall1,ncall2,itmx1,itmx2,ihist 
      common/out/avgi,sd,weight,update
      common/const/pi,s,ncolor,gg,gq,kmax,nb,mp
      common/cut/Etcut,hygap,lygap,hbcutd,lbcutd,ycut,
     &ylcut,hacutd,lacutd,Etmin,yclucut,propcut,yyj
      common/histo/nbhis,y,xhis,count,facwgt,ave
      common/switch/event,exp,sudaf,grid,formfac,
     &scale,Model,evolv,Vgrid,Lcut,interpol,moode,lhe,splash,nscale
      common/Ifactor/c,ael,iglu
      common/uncert/x,xp
      common/BotE/bmax,xmax,bmin
      common/point/Etjetmin,Etjet0,Etjetmax,step
      common/logs/pout1,pout2,pout3,pout4,sout1,sout2,sout3,sout4

      open(unit = 30, file = "LoulJets.log", status = "unknown")
      
      write(30,*)'                '
      write(30,*)'RESULTS:                '
      write(30,*)'    Estimate of sigma=',avgi,' +-',sd
      write(30,*)'    with Etjetmin=',Etjetmin
      write(30,*)'    with Etmin=',Etcut   
      write(30,*)'                '

      if(exp.EQ.1)then
      write(30,*)'   CDF runI parameters:'
      elseif(exp.EQ.2)then
      write(30,*)'   CDF runII parameters:'
       elseif(exp.EQ.3)then
      write(30,*)'   LHC 14 TeV expected parameters:'
       elseif(exp.EQ.4)then
      write(30,*)'   LHC 10 TeV expected parameters:'
       elseif(exp.EQ.5)then
      write(30,*)'   LHC 7 TeV expected parameters:'
      endif
      write(30,*)'                     '
      write(30,*)'LOGS: ' 
      write(30,*)'    SUDA Number of point out of grid:'
      write(30,*)'       ',sout1,' points lower than p1min'
      write(30,*)'       ',sout2,' points lower than p2min'
      write(30,*)'       ',sout3,' points bigger than p1max'
      write(30,*)'       ',sout4,' points bigger than p2max'
      write(30,*)'    DGD Number of point out of grid:'
      write(30,*)'       ',pout1,' points lower than p1min'
      write(30,*)'       ',pout2,' points lower than p2min'
      write(30,*)'       ',pout3,' points bigger than p1max'
      write(30,*)'       ',pout4,' points bigger than p2max'      

      print*,'Estimate of sigma=',avgi,' +-',sd


C  -------------------------------- Output --------------------------------------  C
      write(2,*),Etjetmin,avgi,sd
C  ------------------------------------------------------------------------------  C  

      call hisdump(1,'        k1','    dsigma')
      call hisdump(2,'        k2','    dsigma')
      call hisdump(3,'    theta1','    dsigma')
      call hisdump(4,'        k3','    dsigma')
      call hisdump(5,'    theta3','    dsigma')
      call hisdump(6,'         k','    dsigma')
      call hisdump(7,'    thetak','    dsigma')
      call hisdump(8,'        kp','    dsigma')
      call hisdump(9,'   thetakp','    dsigma')
      call hisdump(10,'        b1','    dsigma')
      call hisdump(11,'     b1+b2','    dsigma')
      call hisdump(12,'        a1','    dsigma')
      call hisdump(13,'     a1+a2','    dsigma')
      call hisdump(14,'       yj1','    dsigma')
      call hisdump(15,'       yj2','    dsigma')
      call hisdump(16,'       yp1','    dsigma')
      call hisdump(17,'       yp2','    dsigma')
      call hisdump(18,'       tgg','    dsigma')
      call hisdump(19,'       ugg','    dsigma')
      call hisdump(20,'       sgg','    dsigma')
      call hisdump(21,' (k+k1)**2','    dsigma')
      call hisdump(22,' (k+k3)**2','    dsigma')
      call hisdump(24,' (k+k2)**2','    dsigma')
      call hisdump(25,'(k+k2-k3)2','    dsigma')
      call hisdump(26,' yj1+yj2/2','    dsigma')
      call hisdump(27,' Et1+Ej2/2','    dsigma')
      call hisdump(28,'       Rjj','    dsigma')
      call hisdump(29,'   x=b1/b2','    dsigma')
      call hisdump(30,'       Mjj','    dsigma')
      call hisdump(31,'        Mx','    dsigma')
      call hisdump(32,'        Rt','    dsigma')
      call hisdump(33,'        Mj','    dsigma')
      call hisdump(34,'       Mgg','    dsigma')
      call hisdump(35,'        nj','    dsigma')
      call hisdump(36,'   yj1-yj2','    dsigma')
      call hisdump(37,' y_cluster','    dsigma')
      call hisdump(38,' Leading j','    dsigma')
      call hisdump(39,' Second  j','    dsigma')

      end

 

C   ============================================================================   C 
C    FUNCTION DSIGMA:                                                              C
C    compute the cross section from the function integ                             C
C          Content: kinematic definition, cuts, rapidity definitions, phase        C
C          space factor, regge factor, gluon/gluon Mandelstam variables and call   C 
C          the integrant                                                           C 
C    dsigma -> Vegas                                                               C
C   ============================================================================   C
      function dsigma(x,wgt)
 	implicit none
      integer it,ndo  
      integer ncall1,ncall2,itmx1,itmx2,ihist
      integer event,exp,sudaf,grid,formfac,Lcut,nscale
      integer scale,Model,evolv,Vgrid,interpol,moode,lhe,splash
      integer jetveto,update
      integer iglu
      integer ik,ik1,ik2
      double precision si,swgt,schi,xi
      double precision x(11),dsigma,wgt,dot,dotdiff,dotsum
      double precision IntegM1,IntegM2,IntegM12
      double precision IntegdurhamM1,IntegdurhamM2,IntegdurhamM12
      double precision k1(2),ak1,theta1,k2(2),ak2,theta2
      double precision k3(2),ak3,theta3
      double precision k(2),ak,thetak,kp(2),akp,thetap
      double precision mk(2),mkp(2),mk1(2),mk2(2),mk3(2)
      double precision skk1(2),skpk1(2),skk3(2),skpk3(2)
      double precision mskk1(2),mskpk1(2),mskk3(2),mskpk3(2)
      double precision skk2(2),skpk2(2),mskk2(2),mskpk2(2)
      double precision skk2k3(2),mskk2k3(2),skpk2k3(2),mskpk2k3(2)
      double precision b1,b2,a1,a2,abcut,hbcutd,lbcutd
      double precision sgg,tgg,ugg,Sudaphi,phi
      double precision ael,c
      double precision kmax,pi,s,ncolor,gg,gq,nb,mp
      double precision weight,weightsuda,weightphi
      double precision fact,factsp,jac,factp,regge
      double precision as,Cf
      double precision gap1,gap2,gap3,gap4
      double precision yp1,yp2,hygap,lygap,yj1,yj2,ycut,yclu
      double precision ylcut,hacutd,lacutd,Etmin,yclucut,propcut,yyj
      double precision Ejet1,Ejet2,Etcut,Etran
      double precision T1,T2
      double precision a1t,a3t
      double precision mmu2,s2
      double precision xmax,amax,amin
      double precision normk1,normk3,normk,normkp,normk2,xb1b2
      double precision sskk1,sskk2,sskk3,sskk2k3
      double precision Etjetveto,Rjj,Mjj,Mx,Rt,Mgg,Mj
      double precision avgi,sd
      double precision nj
      double precision Lj,Sj

      double precision gqq,mq

      external IntegM1,IntegdurhamM1
      external IntegM2,IntegdurhamM2
      external IntegM12,IntegdurhamM12

      common/bveg2/it,ndo,si,swgt,schi,xi(50,11)
      common/paramint/ncall1,ncall2,itmx1,itmx2,ihist
      common/Energy/Ejet1,Ejet2    
      common/const/pi,s,ncolor,gg,gq,kmax,nb,mp
      common/ab/a1,a2,b1,b2
      common/integrand/tgg,sgg,ugg,Sudaphi,phi
      common/Ifactor/c,ael,iglu
      common/cut/Etcut,hygap,lygap,hbcutd,lbcutd,ycut,
     &           ylcut,hacutd,lacutd,Etmin,yclucut,propcut,yyj
      common/switch/event,exp,sudaf,grid,formfac,scale,Model,evolv,
     &              Vgrid,Lcut,interpol,moode,lhe,splash,nscale
      common/sudaD/mmu2,s2
      common/BotE/amax,xmax,amin
      common/lhe/k,k1,k2,k3
      common/out/avgi,sd,weight,update
      common/algo/Etjetveto,Rjj,Mjj,Mx,Rt,Mgg,Mj
      common/kinematic/skk1,skpk1,skk3,skpk3,mskk1,mskpk1,mskk3,mskpk3,
     &                 mk3,skk2,skpk2,mskk2,mskpk2,mk2,skk2k3,
     &                 mskk2k3,skpk2k3,mskpk2k3,mk,mkp,mk1
      common/mquark/mq

       dsigma=0.d0
       jac=1.d0  
  
C     Cinematic for momentum

      ak=kmax*x(1)
      thetak=2.*pi*x(2)
      jac=jac*2.*pi*ak*kmax
      k(1)=ak*cos(thetak)
      k(2)=ak*sin(thetak)

       akp=kmax*x(3)
       thetap=2.*pi*x(4)
       jac=jac*2.*pi*akp*kmax
       kp(1)=akp*cos(thetap)
       kp(2)=akp*sin(thetap)

      ak1=kmax*x(5)
      theta1=2.*pi*x(6)
      jac=jac*2.*pi*ak1*kmax
      k1(1)=ak1*cos(theta1)
      k1(2)=ak1*sin(theta1)

      ak2=5.+10.*kmax*x(7)
      theta2=0.d0
      jac=jac*2.*pi*ak2*kmax*10.
      k2(1)=ak2
      k2(2)=0.d0

      ak3=kmax*x(8)
      theta3=2.*pi*x(9)
      jac=jac*2.*pi*ak3*kmax
      k3(1)=ak3*cos(theta3)
      k3(2)=ak3*sin(theta3)

      call kinematics(k,kp,k1,k2,k3)
   
C     beta et alpha + symetric cuts on
      abcut=0.00000d0    
      b1=abcut+x(10)*(1.-2.*abcut)    
      b2=abcut+x(11)*(1.-2.*abcut)
      jac=jac*(1.-2.*abcut)**2

      a1=(dotdiff(k1,k2,k1,k2)-mq**2)/b1/s
      a2=(dotdiff(k2,k3,k2,k3)-mq**2)/b2/s

      if(a1.LT.abcut.or.a1.GT.(1.-2.*abcut))goto 666
      if(a2.LT.abcut.or.a2.GT.(1.-2.*abcut))goto 666
      if(a1+a2.GT.1.d0.or.b1+b2.GT.1.d0)goto 666

C     Gluon gluon MANDELSTAM variables
      sgg=dot(k2,k2)*(b1+b2)**2/b1/b2
      tgg=-dot(k2,k2)*(b1+b2)/b1
      ugg=-dot(k2,k2)*(b1+b2)/b2

C   Rapidities: antiproton denotes with p2
C   Leading jet is Lj and the second jet is Sj   
      yp1=0.5d0*log((1.-b1-b2)**2*s/(dot(k1,k1)+mp**2))
      yp2=-0.5d0*log((1.-a1-a2)**2*s/(dot(k3,k3)+mp**2))
      yclu=0.5d0*log((b1+b2)/(a1+a2))
      yj1=0.5d0*log(b1/a1)
      yj2=0.5d0*log(b2/a2)

      if(Ejet1.GT.Ejet2)then
       Lj=yj1
       Sj=yj2
      elseif(Ejet2.GT.Ejet1)then
       Lj=yj2
       Sj=yj1
      endif 

      if(yp1.LT.0.d0.or.yp2.GT.0.d0)goto 666

C     Transverse energy of the two jets
      Ejet1=sqrt(dotdiff(k1,k2,k1,k2)-mq**2)
      Ejet2=sqrt(dotdiff(k2,k3,k2,k3)-mq**2)

C     Gap
      gap1=yp1-yj1
      gap2=yp1-yj2
      gap3=-yp2+yj1
      gap4=-yp2+yj2

C   ----------------------------------------------------------------------------   C
C   EXPERIMENTAL CUTS: CDF RUNI
C   ----------------------------------------------------------------------------   C
      if(exp.EQ.1)then

C     Cuts on transverse energy of the two jets
      if(Ejet1.LT.Etcut.or.Ejet2.LT.Etcut)goto 666

      if(grid.EQ.1.or.grid.EQ.2)then
C     Cuts on alpha and beta ordered
      if(a1+a2.GT.hacutd.or.a1+a2.LT.lacutd)goto 666
      if(b1+b2.GT.hbcutd.or.b1+b2.LT.lbcutd)goto 666 

C     Cuts on rapidity
      if(yj1.GT.ycut.or.yj1.LT.ylcut)goto 666    
      if(yj2.GT.ycut.or.yj2.LT.ylcut)goto 666
      endif
    
      if(grid.EQ.2)then
C     Cuts on momentum transfer |t|
      Etran=dot(k3,k3)
      if(Etran.GT.1.)goto 666
      endif

C     Cuts on the gap
      if(min(gap1,gap2).LT.lygap)goto 666
      if(min(gap3,gap4).LT.lygap)goto 666
C      if(max(gap1,gap2).GT.hygap)goto 666
C      if(max(gap3,gap4).GT.hygap)goto 666

C   ----------------------------------------------------------------------------   C
C   EXPERIMENTAL CUTS: CDF RUNII
C   ----------------------------------------------------------------------------   C
      elseif(exp.EQ.2)then

CC      if(Lj.GT.yyj.AND.Sj.GT.yyj)goto 666  
 
C     Cuts on transverse energy of the two jets
      if(Ejet1.LT.Etcut.or.Ejet2.LT.Etcut)goto 666

C     Cuts on alpha and beta ordered
      if(a1+a2.GT.hacutd.or.a1+a2.LT.lacutd)goto 666

C     Cuts on rapidity
      if(abs(yj1).GT.ycut.or.abs(yj2).GT.ycut)goto 666  

C     Cuts on the gap
      if(min(gap1,gap2).LT.lygap)goto 666
      if(min(gap3,gap4).LT.lygap)goto 666

C      if(Lcut.EQ.2)then
C       if(max(gap1,gap2).GT.hygap)goto 666
C       if(max(gap3,gap4).GT.hygap)goto 666
C      endif 

C   ---------------------------------------------------------------------------   C
C   EXPERIMENTAL CUTS: Expected LHC on ATLAS
C   Lcut=1        M.G. Albrow et al. CERN-LHCC-2005-025, Jun 2005.
C   Lcut=2        JHEP 0710 (2007) 090,[arXiv:0709.3035 [hep-ph]].
C   ---------------------------------------------------------------------------   C
      elseif(exp.GE.3)then

C     Cuts on transverse energy of the two jets
        if(Ejet1.LT.Etcut.or.Ejet2.LT.Etcut)goto 666

C     Cuts on the mass of the central system
        if(Lcut.EQ.1)then
         if(sgg.LT.50.**2)goto 666
        elseif(Lcut.EQ.2)then
         if(sgg.LE.80.**2)goto 666 
        endif

C     Cuts on alpha and beta ordered
        if(a1+a2.GT.hacutd.or.a1+a2.LT.lacutd)goto 666
        if(b1+b2.GT.hbcutd.or.b1+b2.LT.lbcutd)goto 666 

C     Cuts on rapidity
        if(abs(yj1).GT.ycut.or.abs(yj2).GT.ycut)goto 666
        if(Lcut.EQ.2)then
         if(yclu.GE.yclucut)goto 666
        endif
      endif
C   ----------------------------------------------------------------------------   C 

C     Test Igor
      T1=dotsum(k,k1,kp,k1)*dotsum(k,k3,kp,k3)+dotsum(k,k1,k,k3)
     &*dotsum(kp,k1,kp,k3)-dotsum(k,k1,kp,k3)*dotsum(k,k3,kp,k1)
      T2=dotsum(k,k1,kp,k1)*dotsum(k,k3,kp,k3)-dotsum(k,k1,k,k3)
     &*dotsum(kp,k1,kp,k3)+dotsum(k,k1,kp,k3)*dotsum(k,k3,kp,k1)

C     Coefficient, function and phase space
C     First pick up the corresponding value of gg(as) from 
C     subroutine and ael,aDD from sigmael.f
      factsp=1./16./(2.*pi)**8/b1/b2
C      factp=81./2.              !3 quarks in p+ and 1/2 because of gluons
      factp=81.                  !3 quarks in p+ and NO 1/2 because qqbar
      gq=sqrt(4.*pi*ael)
C      gg=sqrt(4.*pi*as(sgg))
      gqq=sqrt(4.*pi*as(sgg))    !Sub-process coupling constant
      fact=gqq**4*gq**8*nb/4./pi**4

C     Regge factor
      if(formfac.EQ.1)then
      a1t=0.09d0-0.3d0*dot(k1,k1)/(1.-(b1+b2))
      a3t=0.09d0-0.3d0*dot(k3,k3)/(1.-(a2+a1))
      regge=((1./(a2+a1))**a1t*(1./(b1+b2))**a3t)**2
      else
       regge=1.d0
      endif
      
C   ============================================================================   C
C                        Cross section = m1-m2+2*m1m2                              C
C   ============================================================================   C
      if(Model.LT.3.OR.Model.GT.3)then
       dsigma=jac*fact*factsp*factp*s2*regge*
     &        IntegM1(a1,a2,b1,b2,k1,k2,k3,k,kp)
C     &       IntegM2(a1,a2,b1,b2,k1,k2,k3,k,kp)
C     &       (IntegM12(a1,a2,b1,b2,k1,k2,k3,k,kp)
C     &       +IntegM12(a1,a2,b1,b2,k1,k2,k3,kp,k))

      elseif(Model.EQ.3)then
       dsigma=jac*fact*factsp*factp*s2*regge*
     &        IntegdurhamM1(a1,a2,b1,b2,k1,k2,k3,k,kp)
C     &       IntegdurhamM2(a1,a2,b1,b2,k1,k2,k3,k,kp)
C     &       (IntegdurhamM12(a1,a2,b1,b2,k1,k2,k3,k,kp)
C     &       +IntegdurhamM12(a1,a2,b1,b2,k1,k2,k3,kp,k))
      endif 

       if(abs(dsigma).lt.1.d-25)dsigma=0.d0
C   ============================================================================   C

C      print*,"regge=",regge
C      print*,'Fact,factsp=',factp,factsp
C      print*,'jac=',jac
C      print*,'IntegM1=',IntegM1(a1,a2,b1,b2,k1,k2,k3,k,kp)
C      print*,'IntegM2=',IntegM2(a1,a2,b1,b2,k1,k2,k3,k,kp)
C      print*,'IntegM12=',IntegM12(a1,a2,b1,b2,k1,k2,k3,k,kp)
C      print*,'************************************************'



C   -----------------------------------------------------------------------------  C

C   Events weight
    
      weight=wgt*dsigma
      weightsuda=wgt*Sudaphi*dsigma
      weightphi=wgt*phi*dsigma

C   ----------------------------------------------------------------------------   C 
C   LHE -> input Pythia
C   Cuts and useful infos are in loulib/makeLHE.f and JetsAlgo/src/MidPointAlgo.f  
C   jetveto = 2 means jets don't pass the cuts and return sigma=0 

      if(lhe.EQ.1)then
        call makeLHE
        call pythia
        call midpointalgo(jetveto,nj)
        if(jetveto.EQ.2)goto 666 
      endif
C   ----------------------------------------------------------------------------   C
C     Makes the histograms
      if(ihist.EQ.1) then
      normk1=sqrt(dot(k1,k1))
      normk3=sqrt(dot(k3,k3))
      normk2=sqrt(dot(k2,k2))
      normk=sqrt(dot(k,k))
      normkp=sqrt(dot(kp,kp))
      call histdo(1,normk1,weight)
      call histdo(2,normk2,weight)
      call histdo(3,theta1,weight)
      call histdo(4,normk3,weight)
      call histdo(5,theta3,weight)
      call histdo(6,normk,weight)
      call histdo(7,thetak,weight)
      call histdo(8,normkp,weight)
      call histdo(9,thetap,weight)
      call histdo(10,b1,weight)
      call histdo(11,b1+b2,weight)
      call histdo(12,a1,weight)
      call histdo(13,a1+a2,weight)
      call histdo(14,yj1,weight)
      call histdo(15,-yj2,weight)
      call histdo(16,yp1,weight)
      call histdo(17,yp2,weight)
      call histdo(18,tgg,weight)
      call histdo(19,ugg,weight)
      call histdo(20,sgg,weight)
      sskk1=dotsum(k,k1,k,k1)
      sskk3=dotsum(k,k3,k,k3)
      sskk2=dotsum(k,k2,k,k2)
      sskk2k3=dot(skk2k3,skk2k3)
      xb1b2=b1/b2
      call histdo(21,sskk1,weight)
      call histdo(22,sskk3,weight)
      call histdo(24,sskk2,weight)
      call histdo(25,sskk2k3,weight)
      call histdo(26,(yj1+yj2)/2.d0,weight)
      call histdo(27,(Ejet1+Ejet2)/2.d0,weight)
      call histdo(28,Rjj,weight)
      call histdo(29,xb1b2,weight)
C      call histdo(30,Mjj,weight)
      call histdo(30,sqrt(sgg),weight)
      call histdo(31,Mx,weight)
      call histdo(32,Rt,weight)
      call histdo(33,Mj,weight)
      call histdo(34,Mgg,weight)
      call histdo(35,nj,weight)
      call histdo(36,(yj1-yj2),weight)
      call histdo(37,yclu,weight)
      call histdo(38,Lj,weight)
      call histdo(39,Sj,weight)
      endif

       update=update+1
       if(update.EQ.100)print*,'-> Working: 100 points done'
       if(update.EQ.1000)print*,'-> Working: 1000 points done'
       if(update.EQ.10000)print*,'-> Working: 10 000 points done'
       if(update.EQ.100000)print*,'-> Working: 100 000 points done'
 
C  666   enddo
 666   end

C   ============================================================================   C
C   FUNCTION INTEGM1:                                                              C
C   analytic part of the integrand                                                 C
C              Content: integrand, different defintions for the impact factor and  C
C              for the Sudakov form factor, definition of the coupling constants   C
C   ============================================================================   C
C   ==============================   Diagram 1   ===============================   C
C                --------                                   --------               C
C                 |   |                                       |   |                C
C                 |   |--                                   --|   |                C
C                 |   |--                                   --|   |                C
C                 |   |                                       |   |                C
C                --------                                   --------               C
C   ============================================================================   C
      function IntegM1(a1,a2,b1,b2,k1,k2,k3,k,kp)         
	implicit none
      integer nf,iglu
      integer event,exp,Sudaf,grid,formfac,splash,nscale
      integer scale,Model,evolv,Vgrid,Lcut,interpol,moode,lhe
      double precision dot,dotsum,dothalf,IntegM1,as
      double precision fc,fs
      double precision phi,dphi,ephi 
      double precision b1,b2,k2(2),a1,a2
      double precision k(2),kp(2),k1(2),k3(2)
      double precision mk(2),mkp(2),mk1(2),mk2(2),mk3(2)
      double precision skk1(2),skpk1(2),skk3(2),skpk3(2)
      double precision mskk1(2),mskpk1(2),mskk3(2),mskpk3(2)
      double precision skk2(2),skpk2(2),mskk2(2),mskpk2(2)
      double precision skk2k3(2),mskk2k3(2),skpk2k3(2),mskpk2k3(2)
      double precision I2,I3,I3p,I4,I5
      double precision Ij,C0,C2,Ij0,Ij2
      double precision sgg,tgg,ugg
      double precision I1,J1,J2
      double precision pi,s,ncolor,gg,gq,kmax,nb,mp,cf
      double precision Ejet1,Ejet2
      double precision Sphi,Sphib,Sudaphi
      double precision phip,phipb,Ifact
      double precision dgd2008,dgd07,dgdforward
      double precision bo,ao
      double precision c,ael
      double precision x,xp
      double precision Dsuda
      double precision Rg,prefac
      double precision amax,xmax,amin
      double precision Etcut,hygap,lygap,hbcutd,lbcutd,ycut
      double precision ylcut,hacutd,lacutd,Etmin,yclucut,propcut,yyj
      double precision dgdpol,Dsudapol

      double precision mq,bet

      common/const/pi,s,ncolor,gg,gq,kmax,nb,mp
      common/Energy/Ejet1,Ejet2
      common/integrand/tgg,sgg,ugg,Sudaphi,phi 
      common/flavor/nf
      common/switch/event,exp,Sudaf,grid,formfac,scale,Model,evolv,
     &              Vgrid,Lcut,interpol,moode,lhe,splash,nscale
      common/Ifactor/c,ael,iglu
      common/uncert/x,xp
      common/BotE/amax,xmax,amin
      common/cut/Etcut,hygap,lygap,hbcutd,lbcutd,ycut,
     &           ylcut,hacutd,lacutd,Etmin,yclucut,propcut,yyj
      common/kinematic/skk1,skpk1,skk3,skpk3,mskk1,mskpk1,mskk3,mskpk3,
     &                 mk3,skk2,skpk2,mskk2,mskpk2,mk2,skk2k3,
     &                 mskk2k3,skpk2k3,mskpk2k3,mk,mkp,mk1
      common/mquark/mq

C   Common factor
      Ifact=pi**2/gq**2

C   ----------------------  M1: Num and Den  ---------------------------------   C
CCCCCCCCCCC  GLUONS JETS
C     I1: numerator
c      I1=(1./2.)*(1.+(ugg**4/sgg**4)+(tgg**4/sgg**4))*(1./tgg+1./ugg)**2
C      J1=(1./2.)*(1.-(ugg**4/sgg**4)-(tgg**4/sgg**4))*(1./tgg+1./ugg)**2

C     Numerator
      I2=dotsum(k,k1,kp,k1)*dotsum(k,k3,kp,k3)
      J2=dotsum(k,k1,k,k3)*dotsum(kp,k1,kp,k3)
     &-dotsum(k,k1,kp,k3)*dotsum(k,k3,kp,k1)

C     Angular momentum dependence: Rules
C     J=0 then C2=0 or J=2 and then C0=0 
C      C0=I2+J2
C      C2=(I2-J2)

C      Ij=C0*(I1+J1)/2.+C2*(I1-J1)/2.
C       Ij0=C0*(I1+J1)/2.
C      Ij2=C2*(I1-J1)/2.
CCCCCCCCCC  GLUONS JETS

      bet=sqrt(1.-(4.*mq**2/sgg))
           
C     Numerator
      I1=(1./2.)*(2.+bet)**2/sgg**3
     &    *(1.+(tgg-mq**2)/(ugg-mq**2)+(ugg-mq**2)/(tgg-mq**2))**2
      J1=1./4./mq**2*(bet*sgg-2.*tgg+2.*ugg)/bet/sgg
     &    *((tgg-ugg-bet*sgg)/(tgg-mq**2)/(ugg-mq**2))**2

      I2=dotsum(k,k1,kp,k1)*dotsum(k,k3,kp,k3)
      J2=dotsum(k,k1,k,k3)*dotsum(kp,k1,kp,k3)
     &-dotsum(k,k1,kp,k3)*dotsum(k,k3,kp,k1)

      Ij=(I1+J1)*I2+(I1-J1)*J2

C     Denominator!
      I3=dot(k,k)*dotsum(k,k1,k,k1)*dotsum(k,k3,k,k3)
      I3p=dot(kp,kp)*dotsum(kp,k1,kp,k1)*dotsum(kp,k3,kp,k3)

      if(dot(k,k).LT.propcut.or.dot(kp,kp).LT.propcut)goto 666
      if(dotsum(k,k1,k,k1).LT.propcut)goto 666
      if(dotsum(k,k3,k,k3).LT.propcut)goto 666
      if(dotsum(kp,k1,kp,k1).LT.propcut)goto 666
      if(dotsum(kp,k3,kp,k3).LT.propcut)goto 666  


C   ----------------------  M1: Impact Factor  --------------------------------   C
      if(formfac.EQ.1)then

        if(event.EQ.1)then
         phip=ephi(k,mskk1)*ephi(skpk1,mkp)
         phipb=ephi(mk,skk3)*ephi(mskpk3,kp)
        elseif(event.EQ.2)then
         phip=dphi(k,mskk1,skpk1,mkp)
         phipb=dphi(mk,skk3,mskpk3,kp) 
        elseif(event.EQ.3)then
         phip=dphi(k,mskk1,skpk1,mkp)
         phipb=ephi(mk,skk3)*ephi(mskpk3,kp)
        endif
         
      elseif(formfac.EQ.2)then

        if(event.EQ.1)then
      phip=dgd07(0.41d0*(b1+b2),dothalf(k,k1,k,k1),dot(k1,k1),
     &dot(k,k1)+dot(k1,k1)/2.,iglu)
     &*dgd07(0.41d0*(b1+b2),dothalf(kp,k1,kp,k1)
     &,dot(k1,k1),dot(kp,k1)+dot(k1,k1)/2.,iglu)
        phipb=dgd07(0.41d0*(a1+a2),dothalf(k,k3,k,k3),dot(k3,k3),
     &dot(k,k3)+dot(k3,k3)/2.,iglu)
     &*dgd07(0.41d0*(a1+a2),dothalf(kp,k3,kp,k3)
     &,dot(k3,k3),dot(kp,k3)+dot(k3,k3)/2.,iglu)
        elseif(event.EQ.2)then        
      phip=dgd2008(0.41d0*(b1+b2),dot(k,k),dotsum(k,k1,k,k1),dot(k1,k1))
     &*dgd2008(0.41d0*(b1+b2),dot(kp,kp),dotsum(kp,k1,kp,k1),dot(k1,k1))
      phipb=dgd2008(0.41d0*(a1+a2),dot(k,k),dotsum(k,k3,k,k3)
     &,dot(k3,k3))*dgd2008(0.41d0*(a1+a2),dot(kp,kp)
     &,dotsum(kp,k3,kp,k3),dot(k3,k3))

        elseif(event.EQ.3)then
C     Exact definition of dgdforward -> Very slow
C      phip=dgdforward(0.41d0*(b1+b2),dothalf(k,k1,k,k1),iglu)
C     &*prefac(k,mskk1,dot(k1,k1),b1+b2)
C     &*dgdforward(0.41d0*(b1+b2),dothalf(kp,k1,kp,k1),iglu)
C     &*prefac(skpk1,mkp,dot(k1,k1),b1+b2)
C      phipb=dgdforward(0.41d0*(a1+a2),dothalf(k,k3,k,k3),iglu)
C     &*prefac(mk,skk3,dot(k3,k3),a1+a2) 
C     &*dgdforward(0.41d0*(a1+a2),dothalf(kp,k3,kp,k3),iglu)
C     &*prefac(mskpk3,kp,dot(k3,k3),a1+a2)

C    Interpolation of dgdforward
       phip=dgdpol(b1+b2,dothalf(k,k1,k,k1))
     &*prefac(k,mskk1,dot(k1,k1),b1+b2)
     &*dgdpol(b1+b2,dothalf(kp,k1,kp,k1))
     &*prefac(skpk1,mkp,dot(k1,k1),b1+b2)
      phipb=dgdpol(a1+a2,dothalf(k,k3,k,k3))
     &*prefac(mk,skk3,dot(k3,k3),a1+a2) 
     &*dgdpol(a1+a2,dothalf(kp,k3,kp,k3))
     &*prefac(mskpk3,kp,dot(k3,k3),a1+a2)
        endif
       endif

C   -------------------  M1:Sudakov Form Factor  -------------------------------   C 
C     Durham Sudakov factor
      if(sudaf.EQ.1)then
       if(interpol.EQ.2)then
         if(scale.EQ.1)then
       Sphi=Dsuda(dot(k2,k2)/x,dot(k,k)/xp)
       Sphib=Dsuda(dot(k2,k2)/x,dot(kp,kp)/xp)
         elseif(scale.EQ.2)then
           Sphi=sqrt(Dsuda(dot(k2,k2)/x,dot(skk1,skk1)/xp))
     &      *sqrt(Dsuda(dot(k2,k2)/x,dot(skk3,skk3)/xp))
           Sphib=sqrt(Dsuda(dot(k2,k2)/x,dot(skpk1,skpk1)/xp))
     &      *sqrt(Dsuda(dot(k2,k2)/x,dot(skpk3,skpk3)/xp))
         endif

       elseif(interpol.EQ.1)then
         if(scale.EQ.1)then
       Sphi=Dsudapol(dot(k2,k2)/x,dot(k,k)/xp)
       Sphib=Dsudapol(dot(k2,k2)/x,dot(kp,kp)/xp)
         elseif(scale.EQ.2)then
           Sphi=sqrt(Dsudapol(dot(k2,k2)/x,dot(skk1,skk1)/xp))
     &      *sqrt(Dsudapol(dot(k2,k2)/x,dot(skk3,skk3)/xp))
           Sphib=sqrt(Dsudapol(dot(k2,k2)/x,dot(skpk1,skpk1)/xp))
     &      *sqrt(Dsudapol(dot(k2,k2)/x,dot(skpk3,skpk3)/xp))

          endif
       endif  

C     DDT Sudakov form factor
      elseif(sudaf.EQ.2)then

         if(scale.EQ.1)then
       Sphi=Dsuda(dot(k2,k2)/x,dot(k,k)/xp)
       Sphib=Dsuda(dot(k2,k2)/x,dot(kp,kp)/xp)
         elseif(scale.EQ.2)then
           Sphi=sqrt(Dsuda(dot(k2,k2)/x,dot(skk1,skk1)/xp))
     &      *sqrt(Dsuda(dot(k2,k2)/x,dot(skk3,skk3)/xp))
           Sphib=sqrt(Dsuda(dot(k2,k2)/x,dot(skpk1,skpk1)/xp))
     &      *sqrt(Dsuda(dot(k2,k2)/x,dot(skpk3,skpk3)/xp))
         endif

C     DLA Sudakov form factor
      elseif(sudaf.EQ.3)then

         if(scale.EQ.1)then
       Sphi=Dsuda(dot(k2,k2)/x,dot(k,k)/xp)
       Sphib=Dsuda(dot(k2,k2)/x,dot(kp,kp)/xp)
         elseif(scale.EQ.2)then
           Sphi=sqrt(Dsuda(dot(k2,k2)/x,dot(skk1,skk1)/xp))
     &      *sqrt(Dsuda(dot(k2,k2)/x,dot(skk3,skk3)/xp))
           Sphib=sqrt(Dsuda(dot(k2,k2)/x,dot(skpk1,skpk1)/xp))
     &      *sqrt(Dsuda(dot(k2,k2)/x,dot(skpk3,skpk3)/xp))
         endif
                 

C   Exhume like: Dhuram group with sgg scale
      elseif(sudaf.EQ.4)then 
         
        if(scale.EQ.1)then
         Sphi=Dsuda(sgg/x,dot(k,k)/xp)
         Sphib=Dsuda(sgg/x,dot(kp,kp)/xp)
        elseif(scale.EQ.2)then 
         Sphi=sqrt(Dsudapol(sgg/x,dot(skk1,skk1)/xp))
     &      *sqrt(Dsudapol(sgg/x,dot(skk3,skk3)/xp))
         Sphib=sqrt(Dsudapol(sgg/x,dot(skpk1,skpk1)/xp))
     &      *sqrt(Dsudapol(sgg/x,dot(skpk3,skpk3)/xp)) 
        endif 
  
C   Naive Sudakov Factor      
      elseif(sudaf.EQ.5)then

       Sphi=Dexp(-3./pi/2.*as(dot(k2,k2)/x)
     &*((log(dot(skk3,skk3)/dot(k2,k2)/xp)
     &*log(dot(skk1,skk1)/dot(k2,k2)/xp))))
       Sphib=Dexp(-3./pi/2.*as(dot(k2,k2)/x)
     &*((log(dot(skpk3,skpk3)/dot(k2,k2)/xp)
     &*log(dot(skpk1,skpk1)/dot(k2,k2)/xp))))

C     No Sudakov Factor
      elseif(sudaf.EQ.7)then
       Sphi=1.d0
       Sphib=1.d0 
      endif

C   -------------------------------  Diag1=M1  -----------------------------------  C 
C     Sudakov en Impact
      if(formfac.EQ.1)then
       phi=phip*phipb
      elseif(formfac.EQ.2)then
       phi=Ifact**4*phip*phipb
      endif  
      Sudaphi=Sphi*Sphib

C    Color factor
C      Cf=1./64./ncolor**4
      Cf=1.


C    Function M1
      IntegM1=Cf*Ij*Sudaphi*phi/I3/I3p

C    Debug tools  
C      print*,'In Diag 1:'
C      print*,'phi=',phi,phip,phipb
C      print*,'in phip=',prefac(k,mskk1,dot(k1,k1),b1+b2)
C     &,k,mskk1,dot(k1,k1),b1+b2
C      print*,'Fact,Ifact=',Ifact,pi,gq
C      print*,'Suda=',Sudaphi,Sphi,Sphib
C      print*,'Deno, Numo',I3*I3p,Ij
C      print*,'********************************'    

      goto 66
666   IntegM1=0.
66    end



C   ============================================================================   C
C   FUNCTION INTEGM2:                                                              C
C   analytic part of the integrand                                                 C
C              Content: integrand, different defintions for the impact factor and  C
C              for the Sudakov form factor, definition of the coupling constants   C
C   ============================================================================   C
C   ==============================   Diagram 2   ===============================   C
C                --------                                   --------               C
C                 |   |                                       |   |                C
C                 |   |--                                   --|   |                C
C                 |------                                   ------|                C
C                 |   |                                       |   |                C
C                --------                                   --------               C
C   ============================================================================   C
      function IntegM2(a1,a2,b1,b2,k1,k2,k3,k,kp)         
	implicit none
      integer nf,iglu
      integer event,exp,Sudaf,grid,formfac,splash,nscale
      integer scale,Model,evolv,Vgrid,Lcut,interpol,moode,lhe
      double precision dot,dotsum,dothalf,IntegM2,as
      double precision fc,fs
      double precision phi,dphi,ephi 
      double precision b1,b2,a1,a2
      double precision k(2),kp(2),k1(2),k2(2),k3(2)
      double precision mk(2),mkp(2),mk1(2),mk2(2),mk3(2)
      double precision skk1(2),skpk1(2),skk3(2),skpk3(2)
      double precision mskk1(2),mskpk1(2),mskk3(2),mskpk3(2)
      double precision skk2(2),skpk2(2),mskk2(2),mskpk2(2)
      double precision skk2k3(2),mskk2k3(2),skpk2k3(2),mskpk2k3(2)
      double precision I2,I3,I3p,I4,I5
      double precision Ij,C0,C2,Ij0,Ij2
      double precision sgg,tgg,ugg
      double precision I1,J1,J2
      double precision pi,s,ncolor,gg,gq,kmax,nb,mp,cf
      double precision Ejet1,Ejet2
      double precision Sphi,Sphib,Sudaphi
      double precision phip,phipb,Ifact
      double precision dgd2008,dgd07,dgdforward
      double precision bo,ao
      double precision c,ael
      double precision x,xp
      double precision Dsuda
      double precision Rg,prefac
      double precision amax,xmax,amin
      double precision Etcut,hygap,lygap,hbcutd,lbcutd,ycut
      double precision ylcut,hacutd,lacutd,Etmin,yclucut,propcut,yyj
      double precision dgdpol,Dsudapol

      common/const/pi,s,ncolor,gg,gq,kmax,nb,mp
      common/Energy/Ejet1,Ejet2
      common/integrand/tgg,sgg,ugg,Sudaphi,phi 
      common/flavor/nf
      common/switch/event,exp,Sudaf,grid,formfac,scale,Model,evolv,
     &              Vgrid,Lcut,interpol,moode,lhe,splash,nscale
      common/Ifactor/c,ael,iglu
      common/uncert/x,xp
      common/BotE/amax,xmax,amin
      common/cut/Etcut,hygap,lygap,hbcutd,lbcutd,ycut,
     &           ylcut,hacutd,lacutd,Etmin,yclucut,propcut,yyj
      common/kinematic/skk1,skpk1,skk3,skpk3,mskk1,mskpk1,mskk3,mskpk3,
     &                 mk3,skk2,skpk2,mskk2,mskpk2,mk2,skk2k3,
     &                 mskk2k3,skpk2k3,mskpk2k3,mk,mkp,mk1


C   Common factor
      Ifact=pi**2/gq**2

C   ----------------------  M2: Den and Num  ----------------------------------   C
C     Ij: numerator
      I1=(fc(skk1,skk2)*fc(skpk1,skpk2))+fs(skk1,skk2,skpk1,skpk2)
      I2=(fc(mk,mskk2k3)*fc(mkp,mskpk2k3))+fs(mk,mskk2k3,mkp,mskpk2k3)
      Ij=16.*I1*I2/(sqrt(dot(k2,k2)))**4

C     Denominator!
      I3=dot(k,k)*dotsum(k,k1,k,k1)*dotsum(k,k2,k,k2)*dot(skk2k3,skk2k3)
      I3p=dot(kp,kp)*dotsum(kp,k1,kp,k1)
     &      *dotsum(kp,k2,kp,k2)*dot(skpk2k3,skpk2k3)

      if(dot(k,k).LT.propcut.or.dot(kp,kp).LT.propcut)goto 666
      if(dotsum(k,k1,k,k1).LT.propcut)goto 666
      if(dotsum(k,k2,k,k2).LT.propcut)goto 666
      if(dot(skk2k3,skk2k3).LT.propcut)goto 666
      if(dotsum(kp,k1,kp,k1).LT.propcut)goto 666
      if(dotsum(kp,k2,kp,k2).LT.propcut)goto 666
      if(dot(skpk2k3,skpk2k3).LT.propcut)goto 666

C   ----------------------  M2: Impact Factor  --------------------------------   C
      if(formfac.EQ.1)then

        if(event.EQ.1)then
         phip=ephi(k,mskk1)*ephi(skpk1,mkp)
         phipb=ephi(mskk2k3,skk2)*ephi(mskpk2,skpk2k3)
        elseif(event.EQ.2)then
         phip=dphi(k,mskk1,skpk1,mkp)
         phipb=dphi(mskk2k3,skk2,mskpk2,skpk2k3) 
        elseif(event.EQ.3)then
         phip=dphi(k,mskk1,skpk1,mkp)
         phipb=ephi(mskk2k3,skk2)*ephi(mskpk2,skpk2k3)
        endif
         
      elseif(formfac.EQ.2)then
C  TO BE CHECKED!
        if(event.EQ.1)then
      phip=dgd07(0.41d0*0.5d0*(b1+b2),dothalf(k,k1,k,k1),dot(k1,k1),
     &dot(k,k1)+dot(k1,k1)/2.,iglu)
     &*dgd07(0.41d0*0.5d0*(b1+b2),dothalf(kp,k1,kp,k1)
     &,dot(k1,k1),dot(kp,k1)+dot(k1,k1)/2.,iglu)
        phipb=dgd07(0.41d0*0.5d0*(a1+a2),dothalf(skk2,mk3,skk2,mk3)
     &,dot(k3,k3),dot(skk2,k3)+dot(k3,k3)/2.,iglu)
     &*dgd07(0.41d0*0.5d0*(a1+a2),dothalf(skpk2,mk3,skpk2,mk3)
     &,dot(k3,k3),dot(skpk2,k3)+dot(k3,k3)/2.,iglu)
        elseif(event.EQ.2)then        
      phip=dgd2008(0.41d0*0.5d0*(b1+b2),dot(k,k),
     & dotsum(k,k1,k,k1),dot(k1,k1))
     &*dgd2008(0.41d0*0.5d0*(b1+b2),dot(kp,kp),
     & dotsum(kp,k1,kp,k1),dot(k1,k1))
      phipb=dgd2008(0.41d0*0.5d0*(a1+a2),dot(skk2k3,skk2k3),
     & dotsum(k,k2,k,k2),dot(k3,k3))
     &*dgd2008(0.41d0*0.5d0*(a1+a2),dot(skpk2k3,skpk2k3)
     &,dotsum(kp,k2,kp,k2),dot(k3,k3))

        elseif(event.EQ.3)then
C    Interpolation of dgdforward
       phip=dgdpol(0.5d0*(b1+b2),dothalf(k,k1,k,k1))
     &*prefac(k,mskk1,dot(k1,k1),0.5d0*(b1+b2))
     &*dgdpol(0.5d0*(b1+b2),dothalf(kp,k1,kp,k1))
     &*prefac(skpk1,mkp,dot(k1,k1),0.5d0*(b1+b2))
      phipb=dgdpol(0.5d0*(a1+a2),dothalf(skk2,mk3,skk2,mk3))
     &*prefac(mskk2k3,skk2,dot(k3,k3),0.5*(a1+a2)) 
     &*dgdpol(0.5d0*(a1+a2),dothalf(skpk2,mk3,skpk2,mk3))
     &*prefac(mskpk2,skpk2k3,dot(k3,k3),0.5*(a1+a2))
        endif
       endif


C   ----------------------  Sudakov Form Factor  -------------------------------   C 
C     No Sudakov Factor
       Sphi=1.d0
       Sphib=1.d0 

C   -------------------------------  Dia2=M2  -----------------------------------  C
C     Sudakov en Impact
      if(formfac.EQ.1)then
       phi=phip*phipb
      elseif(formfac.EQ.2)then
       phi=Ifact**4*phip*phipb
      endif  
      Sudaphi=Sphi*Sphib

C    Color factor
      Cf=(ncolor**2-1.)/ncolor**2

C    Function INTEG
      IntegM2=Cf*Ij*Sudaphi*phi/I3/I3p

C    Debug tools  
C      print*,'In Diag 2:'
C      print*,'phi=',phi,skk2,mk3
C      print*,'Fact,Ifact=',Ifact,pi,gq
C      print*,'Suda=',Sudaphi,Sphi,Sphib
C      print*,'Deno, Numo=',Ij
C      print*,'********************************'
   
      goto 66
666   IntegM2=0.
66    end

C   ============================================================================   C
C   FUNCTION INTEGM12:                                                              C
C   analytic part of the integrand                                                 C
C              Content: integrand, different defintions for the impact factor and  C
C              for the Sudakov form factor, definition of the coupling constants   C
C   ============================================================================   C
C   ==============================   Diagram 12  ===============================   C
C                --------                                   --------               C
C                 |   |                                       |   |                C
C                 |   |--                                   --|   |                C
C                 |   |--                                   ------|                C
C                 |   |                                       |   |                C
C                --------                                   --------               C
C   ============================================================================   C
      function IntegM12(a1,a2,b1,b2,k1,k2,k3,k,kp)         
	implicit none
      integer nf,iglu
      integer event,exp,Sudaf,grid,formfac,splash,nscale
      integer scale,Model,evolv,Vgrid,Lcut,interpol,moode,lhe
      double precision dot,dotsum,dothalf,IntegM12,as
      double precision fc,fs
      double precision phi,dphi,ephi 
      double precision b1,b2,a1,a2
      double precision k(2),kp(2),k1(2),k2(2),k3(2)
      double precision mk(2),mkp(2),mk1(2),mk2(2),mk3(2)
      double precision skk1(2),skpk1(2),skk3(2),skpk3(2)
      double precision mskk1(2),mskpk1(2),mskk3(2),mskpk3(2)
      double precision skk2(2),skpk2(2),mskk2(2),mskpk2(2)
      double precision skk2k3(2),mskk2k3(2),skpk2k3(2),mskpk2k3(2)
      double precision I2,I3,I3p,I4,I5
      double precision Ij,C0,C2,Ij0,Ij2
      double precision sgg,tgg,ugg
      double precision I1,J1,J2
      double precision pi,s,ncolor,gg,gq,kmax,nb,mp,cf
      double precision Ejet1,Ejet2
      double precision Sphi,Sphib,Sudaphi
      double precision phip,phipb,Ifact
      double precision dgd2008,dgd07,dgdforward
      double precision Qp,Cff,m,mprime,mb,mbprime
      double precision bo,ao
      double precision c,ael
      double precision x,xp
      double precision Dsuda
      double precision Rg,prefac
      double precision amax,xmax,amin
      double precision Etcut,hygap,lygap,hbcutd,lbcutd,ycut
      double precision ylcut,hacutd,lacutd,Etmin,yclucut,propcut,yyj
      double precision dgdpol,Dsudapol

      common/const/pi,s,ncolor,gg,gq,kmax,nb,mp
      common/Energy/Ejet1,Ejet2
      common/integrand/tgg,sgg,ugg,Sudaphi,phi 
      common/flavor/nf
      common/switch/event,exp,Sudaf,grid,formfac,scale,Model,evolv,
     &              Vgrid,Lcut,interpol,moode,lhe,splash,nscale
      common/Ifactor/c,ael,iglu
      common/uncert/x,xp
      common/BotE/amax,xmax,amin
      common/cut/Etcut,hygap,lygap,hbcutd,lbcutd,ycut,
     &           ylcut,hacutd,lacutd,Etmin,yclucut,propcut,yyj
      common/kinematic/skk1,skpk1,skk3,skpk3,mskk1,mskpk1,mskk3,mskpk3,
     &                 mk3,skk2,skpk2,mskk2,mskpk2,mk2,skk2k3,
     &                 mskk2k3,skpk2k3,mskpk2k3,mk,mkp,mk1


      double precision tprefac1,tprefac2,tprefac3,tprefac4
      common/tmp/tprefac1,tprefac2,tprefac3,tprefac4

C   Common factor
      Ifact=pi**2/gq**2
 
C   Ij: numerator
      I1=(1.+ugg**2/sgg**2+tgg**2/sgg**2)
     &   *(fc(mk2,skk1)*fc(mk2,skk3)*fc(skpk1,skpk2)*fc(mkp,mskpk2k3)
     &   - fs(mk2,skk1,mk2,skk3)*fs(skpk1,skpk2,mkp,mskpk2k3))
      I2=(1.-ugg**2/sgg**2-tgg**2/sgg**2)
     &   *(fc(skpk1,skpk2)*fc(mkp,mskpk2k3)*fs(mk2,skk1,mk2,skk3)
     &   -fc(mk2,skk1)*fc(mk2,skk3)*fs(skpk1,skpk2,mkp,mskpk2k3))
      I4=(1.+ugg**2/sgg**2-tgg**2/sgg**2)
     &   *(fc(mk2,skk3)*fc(mkp,mskpk2k3)*fs(mk2,skk1,skpk1,skpk2)
     &   -fc(mk2,skk1)*fc(skpk1,skpk2)*fs(mk2,skk3,mkp,mskpk2k3))
      I5=(1.-ugg**2/sgg**2+tgg**2/sgg**2)
     &   *(fc(mk2,skk3)*fc(skpk1,skpk2)*fs(mk2,skk1,mkp,mskpk2k3)
     &   -fc(mk2,skk1)*fc(mkp,mskpk2k3)*fs(mk2,skk3,skpk1,skpk2))

      Ij=8.d0*(I1+I2+I4+I5)/(sqrt(dot(k2,k2)))**4

C     Denominator!
      I3=dot(k,k)*dotsum(k,k1,k,k1)*dotsum(k,k3,k,k3)
      I3p=dot(kp,kp)*dotsum(kp,k1,kp,k1)
     &      *dotsum(kp,k2,kp,k2)*dot(skpk2k3,skpk2k3)

      if(dot(k,k).LT.propcut.or.dot(kp,kp).LT.propcut)goto 666
      if(dotsum(k,k1,k,k1).LT.propcut)goto 666
      if(dotsum(k,k3,k,k3).LT.propcut)goto 666
      if(dotsum(kp,k1,kp,k1).LT.propcut)goto 666
      if(dotsum(kp,k2,kp,k2).LT.propcut)goto 666
      if(dot(skpk2k3,skpk2k3).LT.propcut)goto 666


C   ----------------------  M12:Impact Factor  --------------------------------   C
      if(formfac.EQ.1)then

        if(event.EQ.1)then
         phip=ephi(k,mskk1)*ephi(skpk1,mkp)
         phipb=ephi(mk,skk3)*ephi(mskpk2,skpk2k3)
        elseif(event.EQ.2)then
         phip=dphi(k,mskk1,skpk1,mkp)
         phipb=dphi(mk,skk3,mskpk2,skpk2k3) 
        elseif(event.EQ.3)then
         phip=dphi(k,mskk1,skpk1,mkp)
         phipb=ephi(mk,skk3)*ephi(mskpk2,skpk2k3)
        endif
         
      elseif(formfac.EQ.2)then
C  TO BE CHECKED!
        if(event.EQ.1)then
      phip=dgd07(0.41d0*(b1+b2),dothalf(k,k1,k,k1),dot(k1,k1),
     &dot(k,k1)+dot(k1,k1)/2.,iglu)
     &*dgd07(0.41d0*0.5d0*(b1+b2),dothalf(kp,k1,kp,k1)
     &,dot(k1,k1),dot(kp,k1)+dot(k1,k1)/2.,iglu)
        phipb=dgd07(0.41d0*(a1+a2),dothalf(k,k3,k,k3),dot(k3,k3),
     &dot(k,k3)+dot(k3,k3)/2.,iglu)
     &*dgd07(0.41d0*0.5d0*(a1+a2),dothalf(skpk2,mk3,skpk2,mk3)
     &,dot(k3,k3),dot(skpk2,k3)+dot(k3,k3)/2.,iglu)
        elseif(event.EQ.2)then        
      phip=dgd2008(0.41d0*(b1+b2),dot(k,k),
     & dotsum(k,k1,k,k1),dot(k1,k1))
     &*dgd2008(0.41d0*0.5d0*(b1+b2),dot(kp,kp),
     & dotsum(kp,k1,kp,k1),dot(k1,k1))
      phipb=dgd2008(0.41d0*(a1+a2),dot(k,k),dotsum(k,k3,k,k3)
     & ,dot(k3,k3))
     &*dgd2008(0.41d0*0.5d0*(a1+a2),dot(skpk2k3,skpk2k3)
     &,dotsum(kp,k2,kp,k2),dot(k3,k3))

        elseif(event.EQ.3)then
C    Interpolation of dgdforward
       phip=dgdpol((b1+b2),dothalf(k,k1,k,k1))
     &*prefac(k,mskk1,dot(k1,k1),(b1+b2))
     &*dgdpol(0.5d0*(b1+b2),dothalf(kp,k1,kp,k1))
     &*prefac(skpk1,mkp,dot(k1,k1),0.5d0*(b1+b2))
      phipb=dgdpol((a1+a2),dothalf(k,k3,k,k3))
     &*prefac(mk,skk3,dot(k3,k3),(a1+a2)) 
     &*dgdpol(0.5d0*(a1+a2),dothalf(skpk2,mk3,skpk2,mk3))
     &*prefac(mskpk2,skpk2k3,dot(k3,k3),0.5d0*(a1+a2))


      tprefac1=prefac(k,mskk1,dot(k1,k1),(b1+b2))
      tprefac2=prefac(skpk1,mkp,dot(k1,k1),0.5d0*(b1+b2))
      tprefac3=prefac(mk,skk3,dot(k3,k3),(a1+a2)) 
      tprefac4=prefac(mskpk2,skpk2k3,dot(k3,k3),0.5d0*(a1+a2))

      endif
      endif

C   -------------------  M12:Sudakov Form Factor  ------------------------------   C 
C     Durham Sudakov factor
      if(sudaf.EQ.1)then
       if(interpol.EQ.2)then
         if(scale.EQ.1)then
       Sphi=Dsuda(dot(k2,k2)/x,dot(k,k)/xp)
       Sphib=1.d0
         elseif(scale.EQ.2)then
           Sphi=sqrt(Dsuda(dot(k2,k2)/x,dotsum(k,k1,k,k1)/xp))
     &      *sqrt(Dsuda(dot(k2,k2)/x,dotsum(k,k3,k,k3)/xp))
           Sphib=1.d0
         endif

       elseif(interpol.EQ.1)then
         if(scale.EQ.1)then
       Sphi=Dsudapol(dot(k2,k2)/x,dot(k,k)/xp)
       Sphib=1.d0
         elseif(scale.EQ.2)then
           Sphi=sqrt(Dsudapol(dot(k2,k2)/x,dotsum(k,k1,k,k1)/xp))
     &      *sqrt(Dsudapol(dot(k2,k2)/x,dotsum(k,k3,k,k3)/xp))
           Sphib=1.d0
          endif
       endif  

C     DDT Sudakov form factor
      elseif(sudaf.EQ.2)then

         if(scale.EQ.1)then
       Sphi=Dsuda(dot(k2,k2)/x,dot(k,k)/xp)
       Sphib=1.d0
         elseif(scale.EQ.2)then
           Sphi=sqrt(Dsuda(dot(k2,k2)/x,dotsum(k,k1,k,k1)/xp))
     &      *sqrt(Dsuda(dot(k2,k2)/x,dotsum(k,k3,k,k3)/xp))
           Sphib=1.d0
         endif

C     DLA Sudakov form factor
      elseif(sudaf.EQ.3)then

         if(scale.EQ.1)then
       Sphi=Dsuda(dot(k2,k2)/x,dot(k,k)/xp)
       Sphib=1.d0
         elseif(scale.EQ.2)then
           Sphi=sqrt(Dsuda(dot(k2,k2)/x,dotsum(k,k1,k,k1)/xp))
     &      *sqrt(Dsuda(dot(k2,k2)/x,dotsum(k,k3,k,k3)/xp))
           Sphib=1.d0
         endif
                 
C     Naive Sudakov Factor      
      elseif(sudaf.EQ.5)then

       Sphi=Dexp(-3./pi/2.*as(dot(k2,k2)/x)
     &*((log(dotsum(k,k3,k,k3)/dot(k2,k2)/xp)
     &*log(dotsum(k,k1,k,k1)/dot(k2,k2)/xp))))
       Sphib=1.d0

C     No Sudakov Factor
      elseif(sudaf.EQ.7)then
       Sphi=1.d0
       Sphib=1.d0 
      endif

C   -------------------------------  Dia12=M12  ---------------------------------  C
C     Sudakov en Impact
      if(formfac.EQ.1)then
       phi=phip*phipb
      elseif(formfac.EQ.2)then
       phi=Ifact**4*phip*phipb
      endif 

      Sudaphi=Sphi*Sphib

C   Color factor
      Cf=-(ncolor**2-1.)/ncolor**2

C   Function INTEG
      IntegM12=Cf*Ij*Sudaphi*phi/I3/I3p

C   ------------------------------------------------------------------------------  C 

C    Debug tools  
C      print*,'In Diag 12:'
C      print*,'phi=',phi,phipb
C      print*,'Fact,Ifact=',Ifact,Cf
C      print*,'Suda=',Sudaphi,Sphi,Sphib
C      print*,'Deno, Numo=',Ij,I1,I2,I4,I5,1./I3,1./I3p
C      print*,'M12=',IntegM12
C      print*,'********************************'

      goto 66
666   IntegM12=0.
66    end


C   ============================================================================   C
C   SUBROUTINES and FUNCTIONS                                                      C
C   Several functions and definitions                                              C  
C          - Alpha strong                                                          C
C          - LCWF impact factor, elastic and diffractif                            C
C          - UgD forward and prefactor                                             C
C          - Sudakov form factor                                                   C
C          - Interpolation bilinear for UgD (dgdpol) and Sudakov (sudapol)         C
C          - Addition, substraction and multiplication of 4-vectors                C
C   ============================================================================   C
C   ----------------------------------------------------------------------------   C
C   SUBROUTINE KINEMATIC
C   ----------------------------------------------------------------------------   C
      subroutine kinematics(k,kp,k1,k2,k3)
       implicit none 
      double precision k(2),kp(2),k1(2),k2(2),k3(2)
      double precision mk(2),mkp(2),mk1(2),mk2(2),mk3(2)
      double precision skk1(2),skpk1(2),skk3(2),skpk3(2)
      double precision mskk1(2),mskpk1(2),mskk3(2),mskpk3(2)
      double precision skk2(2),skpk2(2),mskk2(2),mskpk2(2)
      double precision skk2k3(2),mskk2k3(2),skpk2k3(2),mskpk2k3(2)
      
      common/kinematic/skk1,skpk1,skk3,skpk3,mskk1,mskpk1,mskk3,mskpk3,
     &                 mk3,skk2,skpk2,mskk2,mskpk2,mk2,skk2k3,
     &                 mskk2k3,skpk2k3,mskpk2k3,mk,mkp,mk1

C   Usefulf sum of vectors
C   k+k1
      skk1(1)=k(1)+k1(1)
      skk1(2)=k(2)+k1(2)
      skpk1(1)=kp(1)+k1(1)
      skpk1(2)=kp(2)+k1(2)

      mskk1(1)=-k(1)-k1(1)
      mskk1(2)=-k(2)-k1(2)
      mskpk1(1)=-kp(1)-k1(1)
      mskpk1(2)=-kp(2)-k1(2)
      
      mk1(1)=-k1(1)
      mk1(2)=-k1(2)

C   k+k3
      skk3(1)=k(1)+k3(1)
      skk3(2)=k(2)+k3(2)
      skpk3(1)=kp(1)+k3(1)
      skpk3(2)=kp(2)+k3(2)

      mskk3(1)=-k(1)-k3(1)
      mskk3(2)=-k(2)-k3(2)
      mskpk3(1)=-kp(1)-k3(1)
      mskpk3(2)=-kp(2)-k3(2)

      mk3(1)=-k3(1)
      mk3(2)=-k3(2)

C   k+k2
      skk2(1)=k(1)+k2(1)
      skk2(2)=k(2)+k2(2)
      skpk2(1)=kp(1)+k2(1)
      skpk2(2)=kp(2)+k2(2)

      mskk2(1)=-k(1)-k2(1)
      mskk2(2)=-k(2)-k2(2)
      mskpk2(1)=-kp(1)-k2(1)
      mskpk2(2)=-kp(2)-k2(2)     

      mk2(1)=-k2(1)
      mk2(2)=-k2(2)

C   k+k2-k3
      skk2k3(1)=k(1)+k2(1)-k3(1)
      skk2k3(2)=k(2)+k2(2)-k3(2)
      skpk2k3(1)=kp(1)+k2(1)-k3(1)
      skpk2k3(2)=kp(2)+k2(2)-k3(2)      

      mskk2k3(1)=-k(1)-k2(1)+k3(1)
      mskk2k3(2)=-k(2)-k2(2)+k3(2)
      mskpk2k3(1)=-kp(1)-k2(1)+k3(1)
      mskpk2k3(2)=-kp(2)-k2(2)+k3(2)      
    
C   k and kp
      mk(1)=-k(1)
      mk(2)=-k(2)
      mkp(1)=-kp(1)
      mkp(2)=-kp(2)
    
      end

C   ----------------------------------------------------------------------------   C
C   SUBROUTINE ALPHA STRONG
C   Compute as in function of Q²:
C      -> choose a Q²
C      -> compute the corresponding lambda(nf)
C      -> compute as
C   ----------------------------------------------------------------------------   C
      function as(q) 
       implicit none 
      integer i,nf 
      double precision as
      double precision m(2:6)
      double precision smu,smc,sms,smb,smt 
      double precision lambda(2:6),la
      double precision q
      double precision pi,s,couleur,gg,gq,kmax,nb,mp

      common/const/pi,s,couleur,gg,gq,kmax,nb,mp
      common/flavor/nf
            
C     Constituent quarks masses and l(5) en GeV from particle data group 
      m(2)=0.312d0
      m(3)=0.45d0
      m(4)=1.5d0
      m(5)=4.5d0
      m(6)=173.4d0
      smu=m(2)**2
      sms=m(3)**2
      smc=m(4)**2
      smb=m(5)**2
      smt=m(6)**2
      lambda(5)=0.2d0
 
      do i=1,3 
      lambda(5-i)=lambda(6-i)*(m(6-i)
     &/lambda(6-i))**(2./(33.-2.*(5-i)))
      end do
      lambda(6)=lambda(5)**(23./21.)/(m(6)**(2./21.))
      continue

      if(q.GT.smt)then 
        nf=6
      elseif(q.GT.smb)then      
        nf=5
      elseif(q.GT.smc)then 
        nf=4
      elseif(q.GT.sms)then
        nf=3 
      elseif(q.GT.smu)then
        nf=2  
      else
      goto 666
      end if
      continue

      la=lambda(nf)
      as=12.*pi/(33.-2.*nf)/log(q/la**2)

      if(as.LT.0.or.as.GT.0.7d0)as=0.7d0

666   end    

C   ----------------------------------------------------------------------------   C
C   LCWF IMPACT FACTOR 
C   Elastic, diffractif, function f(x)~ effective slope B for the proton  
C   ----------------------------------------------------------------------------   C  
C   Elastic Impact factor
C   phi: Chosen from Jr and Hernandez paper
      function ephi(p1,p2)
	implicit none
      integer iglu
      double precision ephi,dot,dotsum,f
      double precision p1(2),p2(2),x1,x2,c,ael

      common/Ifactor/c,ael,iglu         
         
      x1=dotsum(p1,p2,p1,p2)
      x2=dot(p1,p1)+dot(p2,p2)-c*dot(p1,p2)

      ephi=f(x1)-f(x2)
C      if((dot(p1,p1)).LT.10E-8.or.(dot(p2,p2)).LT.10E-8)then
C      ephi=0.
C      endif 

      end
C   ----------------------------------------------------------------------------   C
C   Diffractif Impact factor
C   dphi from C&H p.479
      function dphi(p1,p2,p3,p4)
       implicit none
      double precision p1(2),p2(2),p3(2),p4(2),dphi
      double precision sp1p2(2),sp1p4(2),sp1p3(2)
      double precision f,E3,E2
      double precision x1,x2,x3,x4,x12,x14,x13
      double precision y12,y34,y23,y13,y14,y24
 
       x1=E2(p1)
       x2=E2(p2)
       x3=E2(p3)
       x4=E2(p4)

       sp1p2(1)=p1(1)+p2(1)
       sp1p2(2)=p1(2)+p2(2)
       sp1p4(1)=p1(1)+p4(1)
       sp1p4(2)=p1(2)+p4(2)
       sp1p3(1)=p1(1)+p3(1)
       sp1p3(2)=p1(2)+p3(2)

       x12=E2(sp1p2)
       x14=E2(sp1p4)
       x13=E2(sp1p3)

       y12=E3(p1,p2)
       y34=E3(p3,p4)
       y23=E3(p2,p3)
       y13=E3(p1,p3)
       y14=E3(p1,p4)
       y24=E3(p4,p2)

       dphi=1.d0-f(x1)-f(x2)-f(x3)-f(x4)
     &+2.d0*f(x12)+0.5d0*f(x14)+0.5d0*f(x13)
     &-f(y12)-f(y34)+0.5d0*f(y23)+0.5d0*f(y13)+0.5d0*f(y14)+0.5d0*f(y24)
       dphi=dphi/3.
      end
C  ------------------------------------------------------------------------------  C
C  Function in the Impact factor
      function f(x)
	implicit none
      double precision x,f
      f=(3.53d0+2.79d0*x)/(3.53d0+x)/(1.d0+(x/0.71d0))**2
      if(x.lt.0.)print*,'x<0'
      end

      function E2(p1)
       implicit none
      double precision E2,p1(2),pz(2)
      double precision E3
       
      pz(1)=0.
      pz(2)=0.
      E2=E3(p1,pz)
      end

      function E3(p1,p2)
        implicit none
      integer iglu
      double precision E3,p1(2),p2(2),p3(2)
      double precision c,dot,ael

      common/Ifactor/c,ael,iglu

      p3(1)=-p1(1)-p2(1)
      p3(2)=-p1(2)-p2(2)

      E3=dot(p1,p1)+dot(p2,p2)+dot(p3,p3)
     &-c*(dot(p1,p2)+dot(p2,p3)+dot(p1,p3))
      end

C   ----------------------------------------------------------------------------   C
C   UgD FORWARD PREFACTOR
C   ----------------------------------------------------------------------------   C
      function prefac(p1,p2,t,x)
         implicit none
      double precision p1(2),p2(2),t,x
      double precision aprim,xo,Bo
      double precision f,prefac,dot,tfac
      
      Bo=4.d0
      aprim=0.25d0
      xo=3.4d0*10.E-4

      prefac=2.*dot(p1,p1)*dot(p2,p2)/(dot(p1,p1)**2+dot(p2,p2)**2)
      prefac=prefac*Dexp(-0.5d0*(Bo+2.*aprim*log(xo/x))*t)

C      |t| dependance coming from hep-ph 0802.0177
C       tfac=Dexp(-4.*t)
C       prefac=2.*tfac
C     &*dot(p1,p1)*dot(p2,p2)/(dot(p1,p1)**2+dot(p2,p2)**2)

      return
      end 

C   ---------------------------------------------------------------------------   C
C   SUDAKOV FORM FACTOR
C   ---------------------------------------------------------------------------   C 
      function Dsuda(mu2,p2)
        implicit none
      integer IER
      integer event,exp,Sudaf,grid,formfac,splash,nscale
      integer scale,Model,evolv,Vgrid,Lcut,interpol,moode,lhe
      real*8 Dsuda,mu2,mmu2,p2
      real*8 arg,Iarg
      real*8 ERROR,dcadredo,prec
      double precision as,s2
      
      external dcadredo
      external Iarg

      common/sudaD/mmu2,s2
      common/switch/event,exp,Sudaf,grid,formfac,
     &scale,Model,evolv,Vgrid,Lcut,interpol,moode,lhe,splash,nscale
 
 
      prec=0.001d0
      mmu2=mu2

      if(p2.GE.mu2)then
        Dsuda=1.
      elseif(p2.LT.mu2.or.p2.EQ.mu2)then  
       arg=dcadredo(Iarg,p2,mu2,0.d0,prec,ERROR,IER)
      endif 

      if(evolv.EQ.1)then 
       Dsuda=Dexp(-arg)
      elseif(evolv.EQ.2)then
       Dsuda=Dexp(-as(mu2)*arg)
      endif

      end
C   ---------------------------------------------------------------------------   C
C   Integrand for dcadredo.f
      function Iarg(qt)
        implicit none
      integer nf,event,exp,Sudaf,grid,formfac,Lcut,splash
      integer scale,Model,evolv,Vgrid,interpol,moode,lhe,nscale
      double precision qt,mmu2,s2,nnu2
      double precision tgg,sgg,ugg,Sudaphi,phi
      real*8 pi,s,ncolor,gg,gq,kmax,nb,mp
      real*8 d,Iarg,Ia,Ib,Ilog,Iconst
      real*8 as
      double precision x,xp

      common/const/pi,s,ncolor,gg,gq,kmax,nb,mp
      common/flavor/nf
      common/integrand/tgg,sgg,ugg,Sudaphi,phi
      common/switch/event,exp,Sudaf,grid,formfac,
     &scale,Model,evolv,Vgrid,Lcut,interpol,moode,lhe,splash,nscale
      common/sudaD/mmu2,s2
      common/uncert/x,xp

C   Sudakov prescription: 1=Durham, 2=DDT, 3=DLA
C   Old prescription for delta -> PAPER and sudatab
C       d=sqrt(qt)/(sqrt(mmu2)+sqrt(qt)) 
C   NEW 2010 after Forshaw publication arXiv:0912.3280  
C       d=sqrt(qt)/sqrt(mmu2)


      if(nscale.EQ.1)then
       d=sqrt(qt)/sqrt(mmu2)
      elseif(nscale.EQ.2)then
       nnu2=sgg/x/4.  
       d=sqrt(qt)/sqrt(nnu2)
      endif


      if(sudaf.EQ.1.OR.sudaf.EQ.4)then
       Ia=(-36.0*log(d)-9.0*d**4-2.0*d**3*nf+24.0*d**
     &3+3.0*d**2*nf-54.0*d**2-3.0*d*nf+72.0*d+2.0*nf-33.0)/6.0

       Ilog=-6.*log(d)+(2.*nf-33.)/6.
       Iconst=(-9.0*d**4-2.0*d**3*nf+24.0*d**
     &3+3.0*d**2*nf-54.0*d**2-3.0*d*nf+72.0*d)/6.0

      elseif(sudaf.EQ.2)then
       d=qt/mmu2
       Ia=(2.*float(nf)-33.0-36.0*d-36.0*(d+1.0)**2*log(d)
     &+36.0*(d+1.0)**2*log(d+1.0))/6.0

       Ilog=(2.*float(nf)-33.0-36.0*log(d))/6.0
       Iconst=(-36.0*d-36.0*(d**2+2.*d)*log(d)
     &+36.0*(d+1.0)**2*log(d+1.0))/6.0

      elseif(sudaf.EQ.3)then
       Ia=-6.*log(d)

      elseif(Model.EQ.3.AND.moode.EQ.2.AND.sudaf.EQ.1)then
       d=sqrt(qt)/(sqrt(mmu2)+sqrt(qt))
       Ia=(-36.0*log(d)-9.0*d**4-2.0*d**3*nf+24.0*d**
     &3+3.0*d**2*nf-54.0*d**2-3.0*d*nf+72.0*d+2.0*nf-33.0)/6.0
      endif

      if(evolv.EQ.1)then 
       Ib=as(qt)/2./pi/qt
      elseif(evolv.EQ.2)then
       Ib=1./2./pi/qt
      endif
       
       Iarg=Ia*Ib
       
      end

C   ----------------------------------------------------------------------------   C
C   INTERPOLATION of UgD and SUDAKOV
C   ----------------------------------------------------------------------------   C
      function dgdpol(p1,p2)
        implicit none
      integer pout1,pout2,pout3,pout4,sout1,sout2,sout3,sout4    
      double precision p1,p2,p1log,p2log,dgdpol
      double precision dgdgrid(1:1001,1:1001),sudagrid(1:1001,1:1001)
      double precision lp1max,lp1min,lp2max,lp2min
      double precision step1,step2
      double precision ip1,ip2
      double precision lp1left,lp1right,lp2top,lp2bottom
      double precision f1,f2,f3,f4
      double precision t,u

      common/cinterpol/dgdgrid,sudagrid
      common/logs/pout1,pout2,pout3,pout4,sout1,sout2,sout3,sout4

      p1log=log(p1)
      p2log=log(p2)

C     Dimension and parameter of dgdtab
      lp1min=-14.d0
      lp1max=0.d0
      step1=(lp1max-lp1min)/1000.
      lp2min=-9.d0
      lp2max=8.d0
      step2=(lp2max-lp2min)/1000.

      if(p1log.LT.lp1min)pout1=pout1+1
      if(p2log.LT.lp2min)pout2=pout2+1
      if(p1log.GT.lp1max)pout3=pout3+1
      if(p2log.GT.lp2max)pout4=pout4+1    

      if(p1log.LT.lp1min)p1log=lp1min+step1
      if(p2log.LT.lp2min)p2log=lp2min+step2
      if(p1log.GT.lp1max)p1log=lp1max-step1
      if(p2log.GT.lp2max)p2log=lp2max-step2

C     Bililear interpolation as in Numerical Recipes
      ip1=(p1log-lp1min)/step1
      ip2=(p2log-lp2min)/step2

C     The point to the left is x(1)+ix*stepx
C     The point to the bottom is k(1)+ik*stepk
      lp1left=lp1min+int(ip1)*step1
      lp1right=lp1min+int(ip1+1)*step1
      lp2top=lp2min+int(ip2+1)*step2
      lp2bottom=lp2min+int(ip2)*step2
      
      f1=dgdgrid(int(ip1)+1,int(ip2)+1)
      f2=dgdgrid(int(ip1+1)+1,int(ip2)+1)
      f3=dgdgrid(int(ip1+1)+1,int(ip2+1)+1)
      f4=dgdgrid(int(ip1)+1,int(ip2+1)+1)

      t=(p1log-lp1left)/(lp1right-lp1left)
      u=(p2log-lp2bottom)/(lp2top-lp2bottom)

      dgdpol=(1.-t)*(1.-u)*f1+t*(1.-u)*f2+t*u*f3+(1.-t)*u*f4
   
      return
      end

C   ----------------------------------------------------------------------------   C 
      function Dsudapol(p1,p2)
        implicit none
      double precision p1,p2,p1log,p2log,Dsudapol
      double precision dgdgrid(1:1001,1:1001),sudagrid(1:1001,1:1001)
      double precision lp1max,lp1min,lp2max,lp2min
      double precision step1,step2
      double precision ip1,ip2
      double precision lp1left,lp1right,lp2top,lp2bottom
      double precision f1,f2,f3,f4
      double precision t,u

      common/cinterpol/dgdgrid,sudagrid

      p1log=log(p1)
      p2log=log(p2)

C   Dimension and parameter of dgdtab
C   Old grid notation
C      lp1min=-14.d0
C      lp1max=8.d0
C      step1=(lp1max-lp1min)/1000.
C      lp2min=-14.d0
C      lp2max=8.d0
C      step2=(lp2max-lp2min)/1000.

      lp1min=log(1.E-6)
      lp1max=log(2980.d0)
      step1=(lp1max-lp1min)/1000.
      lp2min=log(1.E-6)
      lp2max=log(2980.d0)
      step2=(lp2max-lp2min)/1000.


      if(p1log.LT.lp1min)p1log=lp1min+step1
      if(p2log.LT.lp2min)p2log=lp2min+step2
      if(p1log.GT.lp1max)p1log=lp1max-step1
      if(p2log.GT.lp2max)p2log=lp2max-step2

C     Bililear interpolation as in Numerical Recipes
      ip1=(p1log-lp1min)/step1
      ip2=(p2log-lp2min)/step2

C     The point to the left is x(1)+ix*stepx
C     The point to the bottom is k(1)+ik*stepk
      lp1left=lp1min+int(ip1)*step1
      lp1right=lp1min+int(ip1+1)*step1
      lp2top=lp2min+int(ip2+1)*step2
      lp2bottom=lp2min+int(ip2)*step2
      
      f1=sudagrid(int(ip1)+1,int(ip2)+1)
      f2=sudagrid(int(ip1+1)+1,int(ip2)+1)
      f3=sudagrid(int(ip1+1)+1,int(ip2+1)+1)
      f4=sudagrid(int(ip1)+1,int(ip2+1)+1)

      t=(p1log-lp1left)/(lp1right-lp1left)
      u=(p2log-lp2bottom)/(lp2top-lp2bottom)

      Dsudapol=(1.-t)*(1.-u)*f1+t*(1.-u)*f2+t*u*f3+(1.-t)*u*f4
   
      return
      end

C   ----------------------------------------------------------------------------   C  
C     OTHER FUNCTION for 4-VECTORS
C   ----------------------------------------------------------------------------   C
C   Definition of the scalar product
      function dot(p1,p2)
	implicit none
      double precision p1(2),p2(2),dot
      dot=p1(1)*p2(1)+p1(2)*p2(2)
      end
C   ----------------------------------------------------------------------------   C
C   Definition of the scalar product of two sum of two vectors
      function dotsum(p1,p2,p3,p4)
	implicit none
      double precision p1(2),p2(2),p3(2),p4(2),dotsum
      dotsum=(p1(1)+p2(1))*(p3(1)+p4(1))+(p1(2)+p2(2))*(p3(2)+p4(2))
      end
C   ----------------------------------------------------------------------------   C
C   Definition of the scalar product of two difference of two vectors
      function dotdiff(p1,p2,p3,p4)
	implicit none
      double precision p1(2),p2(2),p3(2),p4(2),dotdiff
      dotdiff=(p1(1)-p2(1))*(p3(1)-p4(1))+(p1(2)-p2(2))*(p3(2)-p4(2))
      end
C   ----------------------------------------------------------------------------   C
C     Definition of the scalar product of two difference in unintegrated 
C     distribution from Igor
      function dothalf(p1,p2,p3,p4)
        implicit none
      double precision p1(2),p2(2),p3(2),p4(2),dothalf
      dothalf=(p1(1)+0.5d0*p2(1))*(p3(1)+0.5d0*p4(1))
     &+(p1(2)+0.5d0*p2(2))*(p3(2)+0.5d0*p4(2))
      end
C   ----------------------------------------------------------------------------   C
C     Definition of the ~cos of the angle between two 4-vectors
C     [cf Dijet2,b3 and g3]
      function fc(p1,p2)
        implicit none
      double precision p1(2),p2(2),fc
      double precision dot
      fc=dot(p1,p2)
      end
C   ----------------------------------------------------------------------------   C
C     Definition of the ~sin of the angle between two 4-vectors
C     [cf Dijet2,b3 and g3]
      function fs(p1,p2,p3,p4)
        implicit none
      double precision p1(2),p2(2),p3(2),p4(2),fs
      double precision dot
      fs=dot(p1,p3)*dot(p2,p4)-(dot(p1,p4)*dot(p2,p3))
      end
C   ============================================================================   C
