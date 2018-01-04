c $Id: dwidth.f,v 1.11 2007/01/30 14:50:24 bleicher Exp $
C####C##1#########2#########3#########4#########5#########6#########7##
      real*8 function mmean (io,m0,g,mmin,mmax) 
c
cinput  io  : flag (see bleow)
cinput  m0  : pole mass
cinput  g   : nominal width 
cinput mmin : minimal mass
cinput mmax : maximal mass
c
c   io=0 : Yields average mass between {\rm mmin} and {\rm mmax}
c          according to a Breit-Wigner function with constant width {\rm g}
c          and pole {\rm m0}.\\
c   io=1 : Chooses a mass randomly between {\rm mmin} and {\rm mmax} 
c            according to a Breit-Wigner function with constant 
c            width {\rm g} and pole {\rm m0}.\\
c    else: Integral of a Breit-Wigner function from {\rm mmin} to {\rm mmax}.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none

      real*8 m0,g,mmin,mmax,x,i0,i1,inv,fmin,fmax,ranf,f,gcut
      parameter(gcut=1d-10)
      integer io
      logical errchk
      parameter (errchk=.true.)

      i0(x) =2.*g*atan( 2.* (x-m0)/ g  )
      i1(x) =.5*g**2*log( (x-m0)**2+g**2/4. ) + m0*i0(x)
      inv(x)=.5*g*tan( 0.5*x/g )+m0


c...check for some error conditions
      if(errchk)then
        if(mmin.gt.mmax)
     .    write(6,*)'mmean: mass range negative (mmin>mmax)'
     .          ,mmin,mmax
        if(g.le.gcut.and.(m0.gt.mmax.or.m0.lt.mmin))
     .    write(6,*)'mmean: narrow particle out of mass range'  
      end if

      if(io.eq.0)then
        if(g.le.gcut)then
          mmean=0d0
          if(mmin.le.m0.and.m0.le.mmax)mmean=1d0
        else
          mmean=(i1(mmax)-i1(mmin))/(i0(mmax)-i0(mmin))
        end if
      else if(io.eq.1)then
c... determin a mass in a given interval
        if(g.le.gcut)then
          mmean=max(mmin,min(mmax,m0))
        else
          fmin=i0(mmin)
          fmax=i0(mmax)
          f=fmin+(fmax-fmin)*ranf(0)
          mmean=inv(f)      
        end if
      else
            mmean=i0(mmax)-i0(mmin)
      end if      
      return
      end


C####C##1#########2#########3#########4#########5#########6#########7##
      subroutine getmas(m0,g0,i,iz,mmin,mmax,mrest,m)
c
cinput  m0   : pole mass of resonance
cinput  g0   : nominal width of resonance 
cinput  i    : resonance ID
cinput  iz   : iso3 of resonance
cinput  mmin : minimal mass
cinput  mmax : maximal mass
coutput m    : actual mass of the resonance
c
c       {\tt getmas} (not $\rightarrow$ {\tt getmass}) first chooses the 
c       mass {\tt m} of resonance {\tt i} between {\tt mmin} and {\tt mmax}
c       by a call of {\tt mmean}. Since {\tt mmean} only handles Breit-Wigners
c       with constant widths it follows 
c       a correction such that  {\tt m} is distributed according to mass 
c       dependent widths (corresponding to {\tt fbrwig(...,m,1)}).
c
C####C##1#########2#########3#########4#########5#########6#########7##
      implicit none
      include 'options.f'
      include 'comres.f'
      integer i,iz,nrej, nrejmax
      real*8 m,m0,g0,mmin,mmax,x,x0,gg,f,g,h,pi,al,alpha,ce,mmax2
      real*8 phi,k,k0,mrest
c...functions
      real*8 ranf,mmean,fbrwig,bwnorm,pcms
      parameter(pi=3.1415926535d0)
      parameter(alpha=3d0, ce=2d0, nrejmax=1000000)

c 'broadened' Breit-Wigner function with h(x0,al)=h(x0,1)
c normalised to alpha
      h(x,x0,gg,al)=al*0.5/pi*(al*gg)/((x-x0)**2+0.25*(al*gg)**2)

c cut-off for maximum resonance mass
      mmax2=min(mresmax,mmax)


      if(g0.lt.1d-4.or.CTOption(1).ne.0.or.CTOption(32).ne.0)then
        m=mmean(1,m0,g0,mmin,mmax2)
        return
      else
          nrej=0

c This is a Monte Carlo rejection method, where the invertable 
c BW-distribution with constant widths is used to limit the BW-distribution
c with mass-dep. widths whose inverse is not known analytically.

108     continue
        m=mmean(1,m0,alpha*g0,mmin,mmax2)
        if(m.gt.(mmax2+1d-8).or.m.lt.(mmin-1d-8))then
           write(*,*)'getmas (E): m outside (mmin,mmax2)',m,mmin,mmax2
           write(*,*)'called as getmas(',m0,g0,i,mmin,mmax,')'
           write(*,*)'Program stopped'
           stop 137
        endif
        if ((CTOption(25).eq.1).and.(mrest.gt.0.0)) then
           k=pcms(mmax2+mrest,mrest,m)
           k0=pcms(mmax2+mrest,mrest,mmin)
           phi = m*k / (mmin*k0)
        else
           phi = 1.0
        endif

c Breit-Wigner with mass dependent widths and phase space correction

            f=fbrwig(i,iz,m,1)*phi/bwnorm(i)
            g=ce*h(m,m0,g0,alpha)

            if(g.lt.f)then
c              write(*,*)'(W) getmas: C h(m) not limiting at m=',m
c              write(*,*)'->mass distribution of ',i,'might be corrupt'
            endif
          nrej=nrej+1
        if (nrej.le.nrejmax.and.(ranf(0)*g).gt.f) goto 108
        if (nrej.gt.nrejmax) then
c           write(*,*)'(W) getmas_space: too many rejections, m= ',m
c           write(*,*)'called with (',m0,g0,i,mmin,mmax,mrest,')'
c           write(*,*)'->mass distribution of ',i,' might be corrupt'
           m=mmean(1,m0,alpha*g0,mmin,mmax2)
        endif
        
      endif

      return
      end

C####C##1#########2#########3#########4#########5#########6#########7##
        real*8 function bwnorm(ires)
c
cinput  ires   : itype of resonance
c
c       This function calculates the integral of {\tt fbrwig}
c       between parameters {\tt mmin}(= 0~GeV) and {\tt mmax}(= 30~GeV) by
c       calling {\tt qsimp3} resp. by table lookup. It's value shall 
c       serve as the norm of the Breit-Wigner function of particle {\tt ires} 
c       with mass dependent width.
c
C####C##1#########2#########3#########4#########5#########6#########7##

        implicit none
        include 'comres.f'
        include 'comwid.f'
        include 'options.f'
        integer ires,iz,isoit,it
        real*8 mmin,mmax,pole,norm1,norm2
        real*8 widit,massit
        parameter(mmin=0d0,mmax=50d0)
        real*8 fbrwig
        external fbrwig

        if((CTOption(36).ne.0.or.CTOption(1).ne.0).and.wtabflg.gt.1)then
          bwnorm=1d0
          return
        endif
        
        it=iabs(ires)

      if (wtabflg.ge.2.and.CTOption(33).eq.0) then
c table lookup
         if (it.ge.minbar.and.it.le.maxbar) then
           bwnorm=bwbarnorm(it)
         else if (it.ge.minmes.and.it.le.maxmes) then
           bwnorm=bwmesnorm(it)
         else 
           write (6,*) '*** error(bwnorm) wrong id:',it
           bwnorm=1d0
         endif
        else
c calculate
           if (widit(it).gt.1d-3)then

           pole=massit(it)
c arbitrary value of iz 
             iz=isoit(it)

c     the integration is divided by the pole of the Breit-Wigner -
c     thus two integrations with pole as upper or lower boundary
c     respectively are necessary

             call qsimp3(fbrwig,mmin,pole,it,iz,norm1,-1)
             call qsimp3(fbrwig,pole,mmax,it,iz,norm2,+1)
             bwnorm=norm1+norm2
           else
             bwnorm=1d0
           endif
        endif

        return
        end


c no physics after these routines!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c(c) numerical receipies, adapted for f(idum1,idum2,x)
      SUBROUTINE qsimp3(func,a,b,idum1,idum2,s,flag)
c
c     Simpson integration via Numerical Receipies.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none

      include 'options.f'

      INTEGER JMAX,j,idum1,idum2,flag
      REAL*8 a,b,func,s,EPS
      REAL*8 os,ost,st
      PARAMETER (JMAX=100)
      external func
      if(b-a.le.1.d-4) then
         s=0.d0
         return
      endif

      EPS = 5d-3
      if (CTOption(35).eq.1) EPS=5d-4

      ost=-1.d30
      os= -1.d30
      do 11 j=1,JMAX
         if(flag.eq.-1) then
            call midsqu3(func,a,b,idum1,idum2,st,j)
         elseif(flag.eq.1) then
            call midsql3(func,a,b,idum1,idum2,st,j)
         endif
        s=(9.*st-ost)/8.

        if (abs(s-os).le.EPS*abs(os)) return
        os=s
        ost=st
11    continue

      write(6,*)  'too many steps in qsimp3, increase JMAX!'

      return      
      END


      SUBROUTINE midsqu3(funk,aa,bb,idum1,idum2,s,n)
c     modified midpoint rule; allows singuarity at upper limit
      implicit none
      integer idum1,idum2
      INTEGER n
      REAL*8 aa,bb,s,funk
      EXTERNAL funk
      INTEGER it,j
      REAL*8 ddel,del,sum,tnm,x,func,a,b
      func(x)=2.*x*funk(idum1,idum2,bb-x**2,1)
      b=sqrt(bb-aa)
      a=0.
      if (n.eq.1) then
        s=(b-a)*func(0.5*(a+b))
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.*tnm)
        ddel=del+del
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.
      endif
      return
      END

      SUBROUTINE midsql3(funk,aa,bb,idum1,idum2,s,n)
c     modified midpoint rule; allows singularity at lower limit
      implicit none
      integer idum1,idum2
      INTEGER n
      REAL*8 aa,bb,s,funk
      EXTERNAL funk
      INTEGER it,j
      REAL*8 ddel,del,sum,tnm,x,func,a,b
      func(x)=2.*x*funk(idum1,idum2,aa+x**2,1)
      b=sqrt(bb-aa)
      a=0.
      if (n.eq.1) then
        s=(b-a)*func(0.5*(a+b))
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.*tnm)
        ddel=del+del
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.
      endif
      return
      END


