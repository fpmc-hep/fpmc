C     ================================================================== 
C     Histogramming package
C     ==================================================================

C     Make the plot with 'n' the number of the plot, 'x' the variable
C     and 'wght' the weight 
C     'ave' is the average

      subroutine histdo(n,x,wght)
       implicit none
      integer n,ibin
      double precision facwgt,wght,x
      double precision y(50,100),xhis(50,100),count(50),ave(50)
      double precision nbhis(50)
      common/histo/nbhis,y,xhis,count,facwgt,ave

      count(n)=count(n)+wght
      ave(n)=ave(n)+x*wght
      ibin=0
      if(n.ge.1)then
c---
1     ibin=ibin+1
c--
      if(x.lt.xhis(n,ibin))then
      y(n,ibin)=y(n,ibin)+wght
      elseif(ibin.lt.nbhis(n))then
      goto 1
      else
      y(n,nbhis(n))=y(n,nbhis(n))+wght

      endif
c---
      else
c---
      if(x.gt.xhis(n,nbhis(n)))then
      y(n,nbhis(n))=y(n,nbhis(n))+wght
      return
      else
c--
2     ibin=ibin+1
c-
      if(x.lt.xhis(n,ibin))then
      y(n,ibin)=y(n,ibin)+wght
      endif
c-
      if(ibin.lt.nbhis(n))goto 2
      endif
c--
      endif
      return
      end

C     ====================================================================
C     Generating the histograms and the columns

      subroutine histart(n,xmin,xmax,nbin)
        implicit none
      integer i,n,nbin
      double precision facwgt
      double precision xmax,xmin
      double precision y(50,100),xhis(50,100),count(50)
      double precision ave(50),nbhis(50)
      common/histo/nbhis,y,xhis,count,facwgt,ave   

      count(n)=0.
      ave(n)=0.
      nbhis(n)=nbin
      do i=1,nbin
      y(n,i)=0.d0
      xhis(n,i)=xmin+(xmax-xmin)/float(nbin)*i
      end do 
      continue
      return
      end

C     ==================================================================
C     Initialization of the histograms 
 
      subroutine hinit(name)
       implicit none
      integer ncall,itmx,nprn,ndev
      integer ncall1,ncall2,itmx1,itmx2
      double precision pi,s,ncolor,gg,gq,kmax,nb,mp
      double precision xl,xu,acc
      double precision y(50,100),xhis(50,100),count(50)
      double precision ave(50),nbhis(50),facwgt

      common/bveg1/ncall,itmx,nprn,ndev,xl(11),xu(11),acc
      common/paramint/ncall1,ncall2,itmx1,itmx2 
      common/const/pi,s,ncolor,gg,gq,kmax,nb,mp
      common/histo/nbhis,y,xhis,count,facwgt,ave
      character*10 name

      write(30,100)'rs=',sqrt(s)
      write(30,*)'statistics Vegas ',ncall1,' X',itmx1
      write(30,*)'statistics Vegas1 ',ncall2,' X',itmx2 
     
      facwgt=1./itmx2
     
100	format(1x,1a,1f10.1,2x,1a,1f10.1,2x,1a,1f10.1,2x,1a,2f5.2) 
101	format(1x,1a,1i6,2x,1a,1i2) 

      return
      end
C     ==================================================================
C     Makes the columns of numbers
 
      subroutine hisdump(n,title,name)
       implicit none
      integer n,i,nbin
      double precision facwgt,xbin,res,x
      double precision y(50,100),xhis(50,100),count(50)
      double precision ave(50),nbhis(50)
      character*10 name,title
      common/histo/nbhis,y,xhis,count,facwgt,ave

      if(count(n).eq.0.)then
	return
	endif
      nbin=nbhis(n)
100	xbin=xhis(n,2)-xhis(n,1)
      if(xbin.eq.0)print*,'problem in hisdump'
      write(30,*)title
      write(30,*)'average= ',sngl(ave(n)/count(n))
      write(30,*)'total contrib. ',count(n)*facwgt
c	print*,title
c	print*,'average= ',sngl(ave(n)/count(n))
c	print*,'total contrib. ',count(n)*facwgt

      do i=1,nbin
      res=y(n,i)*facwgt/xbin
      x=xhis(n,i)-xbin/2.
      write(7+n,105)x,res
 105	format(1x,2g14.3)
      end do
      continue
      return
      end
C     =======================================================================
