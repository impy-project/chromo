c$Id: numrec.f,v 1.5 1999/01/18 09:57:10 ernst Exp $
C=======================================================================
C Routines taken from Numerical Recipes
C
C

c
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      implicit integer (i - n)
      implicit real*8 (a - h , o - z)
      INTEGER n,NMAXsp
      PARAMETER (NMAXsp=1000)
      REAL*8 yp1,ypn,x(NMAXsp),y(NMAXsp),y2(NMAXsp)
      INTEGER i,k
      REAL*8 p,qn,sig,un,u(NMAXsp)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END

      FUNCTION ran1(idum)
      implicit integer (i - n)
c     real*4 is on purpose, do not change!!!
      implicit real*4 (a - h , o - z)
      INTEGER*4 idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL*4 AM,EPS,RNMX
      real*8 ran1
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER*4 j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=dble(min(AM*iy,RNMX))
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .


