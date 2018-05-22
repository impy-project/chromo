c $Id: getmass.f,v 1.9 1999/01/18 09:57:02 ernst Exp $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function getmass(ssqrt,type)
c
c     Revision : 1.0 
c
cinput  ssqrt  :  maximum energy available
cinput  type   :  class of resonance defined in getrange
coutput getmass:  mass of the resonance
c
C  DETERMINES MASS OF A non-strange baryon RESONANCE.
c
c     function calls: massdist ranf
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      real*8 ssqrt,mdice,mm,mmax
      integer type
      real*8 ranf
      include 'comnorm.f'
      include 'comres.f'

c...  only n*,d,d* included in n_splint

      mmax=min(mresmax,ssqrt)

c get probability mdice
      call n_splint(x_norm,y_norm,n,mmax,mdice,type)
      if(mdice.eq.0)write(6,*)'getmass:mdice=',mdice
c for this probability choose mass mm
      call n_splint(y_norm,x_norm,n,mdice*ranf(0),mm,type)

      if(mm.lt.mresmin) then
c         write(6,*)'(W) getmass-error - m, mmin',mm,mresmin
           mm=mresmin
      else if(mm.gt.mresmax) then 
c         write(6,*)'(W) getmass-error - m, ,max',mm,mresmax
           mm=mresmax
      endif

      getmass=mm

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine norm_init
c
coutput : via common-block comnorm
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'comres.f'
      include 'comnorm.f'
      real*8 massdist
      real*8 x(n),y(n)
      integer i,j
c    linear weighting restored
c    in comnorm.f n is set to 400 again
      dx = (mresmax-mresmin)/dble(n-1)
      do 1 i=0,3
         y_norm(i,1) = 0.0
         y(1) = 0.0
         x_norm(i,1) = mresmin
         x(1) = mresmin
         do 2 j=2,n
            x_norm(i,j) = mresmin+dble(j-1)*dx ! (j-1)**2 for quad. weight
            x(j) = x_norm(i,j)
            y(j) = y(j-1) + (x(j)-x(j-1))
     &                      *(massdist(x(j-1),i)+massdist(x(j),i)
     &                  + 4.0*massdist(0.5*(x(j)+x(j-1)),i))/6.0 
            y_norm(i,j) = y(j)
2        continue
1     continue
      return
      end

c Numerical recipies:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  (modified)
      SUBROUTINE n_splint(xa,ya,n,x,y,m)
      implicit none
      INTEGER n,m
      real*8 x,y,xa(0:3,n),ya(0:3,n)
      INTEGER k,khi,klo
      real*8 a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(m,k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(m,khi)-xa(m,klo)
      if (h.eq.0.) then
       write(6,*) 'bad xa input in splint'
      end if
      a=(xa(m,khi)-x)/h
      y=a*ya(m,klo)+(1.0d0-a)*ya(m,khi)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software.

