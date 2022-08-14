c $Id: cascinit.f,v 1.13 2005/03/24 16:46:58 bleicher Exp $
      subroutine cascinit(ZZ,AA,nucleus)

      implicit none

      include 'coms.f'
      include 'options.f'
      include 'inputs.f'

      integer Z,A,i,getspin,fchg,ZZ,AA,antia, j,k, nrjc,izcnt,nucleus
      real*8 R,P,R2,ranf,sint,phi,nucrad,cost,drr
      real*8 xcm,ycm,zcm,pxcm,pycm,pzcm,densnorm,add,avd,ws
      real*8 r_long, r_short, rdef, ratio, alpha
c
c mass correction according to binding energy
      real*8 meff
ce bcorr to prevent recycling odd nuclei
        logical bcorr
        common /ini/ bcorr
        bcorr=.false.

C averaged density of one Gaussian in units of central value
      avd = 0.5**1.5 

      pxcm=0.d0
      pycm=0.d0
      pzcm=0.d0

      A=abs(AA)
      Z=abs(ZZ)
      antia=A/AA

      PT_AA(nucleus)=A

      if(A.gt.AAmax*MaxNTest.or.A.gt.AAmax*CTParam(67)) then
         write(6,*)'***(E): Mass ',A,' exceeds initialization limit!'
         write(6,*)'        -> increase parameter AAmax=',AAmax
         write(6,*)'           in include-file inputs.f '
         write(6,*)'        -> if test particles are used check also'
         write(6,*)'           if MaxNTest in inputs.f has to be '
         write(6,*)'           increased'
         write(6,*)'        -> uqmd aborting ... '
         stop 137
      endif

      if (A.eq.1*int(CTParam(67))) then  ! proton or neutron
        i = 1
        PT_iso3(i,nucleus) = antia*(-1+2*Z)
        PT_ityp(i,nucleus) = antia*1
        uid_cnt=uid_cnt+1
        PT_uid(i,nucleus)=uid_cnt
        PT_spin(i,nucleus) = 1
        PT_dectime(i,nucleus)=1.d31
        PT_charge(i,nucleus)=fchg(PT_iso3(i,nucleus),
     &                            PT_ityp(i,nucleus))
        PT_fmass(i,nucleus) = EMNUC
        PT_p0(i,nucleus)=sqrt(PT_fmass(i,nucleus)**2)
        return
      end if

C
C Initialisation of A nucleons in configuration space
C
      xcm=0.d0
      ycm=0.d0
      zcm=0.d0
      if (CTOption(24).eq.1) then
        R2 = nucrad(A) + 10.0
      else
        R2 = nucrad(A)    
      endif
      nrjc = 0
      if(CTOption(24).eq.2)then ! use fast method
        call nucfast(A,nucleus)
        do j=1,A
          PT_dectime(j,nucleus)=1.d31
          xcm=xcm+PT_rx(j,nucleus)
          ycm=ycm+PT_ry(j,nucleus)
          zcm=zcm+PT_rz(j,nucleus)
        enddo
      else
      do 1 j=1,A
        PT_dectime(j,nucleus)=1.d31
111     nrjc = nrjc+1
        R = R2*ranf(0)**(1./3.)
        cost = 1.-2.*ranf(0)
        sint = sqrt(1.-cost**2)
        phi = 2.*Pi*ranf(0)
        PT_r0(j,nucleus) = 0.d0
        PT_rx(j,nucleus) = R*sint*cos(phi)
        PT_ry(j,nucleus) = R*sint*sin(phi)
        PT_rz(j,nucleus) = R*cost
        if (CTParam(21).gt.0.0) then
           ratio = sqrt((1 + 4.0*CTParam(21)/3.0) / 
     $                  (1 - 2.0*CTParam(21)/3.0) )
           alpha = atan( sqrt(PT_rx(j,nucleus)*PT_rx(j,nucleus)
     $     + PT_ry(j,nucleus)*PT_ry(j,nucleus))/PT_rz(j,nucleus)) 
           r_short = nucrad(A)*(ratio**(-(1.0/3.0)))
           r_long  = nucrad(A)*(ratio**(2.0/3.0))
           rdef = r_short/sqrt(1-((r_long**2-r_short**2)/r_long**2)
     $            *cos(alpha)*cos(alpha))
        else
           rdef = nucrad(A)
        endif
      

        if (CTOption(24).eq.1) then
           WS = 1 / ( 1 + dexp( ( R - rdef ) / 0.545 ) )
           if (ranf(0) .gt. WS ) then
c               write (6,*) "rejected: ",R
              goto 111
           endif
        else
           do 11 k=1,j-1
              drr=(PT_rx(j,nucleus)-PT_rx(k,nucleus))**2
     &             +(PT_ry(j,nucleus)-PT_ry(k,nucleus))**2
     &             +(PT_rz(j,nucleus)-PT_rz(k,nucleus))**2
              if (drr.lt.2.6/CTParam(67)**(1./3.).and.
     &           nrjc.lt.CTParam(46)) goto 111
 11        continue
        endif
        xcm=xcm+PT_rx(j,nucleus)
        ycm=ycm+PT_ry(j,nucleus)
        zcm=zcm+PT_rz(j,nucleus)
1     continue
c end of fast init if statement
      endif
      if (nrjc.ge.CTParam(46)) then
c         write(6,*)'*** warning: initialisation corrupt '
         bcorr=.true.
      end if

      xcm = xcm/dble(A)
      ycm = ycm/dble(A)
      zcm = zcm/dble(A)
      do 13 j=1,A
        if (CTOption(24).ne.1) then
          PT_rx(j,nucleus) = PT_rx(j,nucleus)-xcm 
          PT_ry(j,nucleus) = PT_ry(j,nucleus)-ycm 
          PT_rz(j,nucleus) = PT_rz(j,nucleus)-zcm
        endif
        PT_rho(j,nucleus) = avd
13    continue

C local proton density in nucleus A,Z
      do 14 j=1,Z
        do 15 k=j+1,Z
          drr=(PT_rx(j,nucleus)-PT_rx(k,nucleus))**2
     &       +(PT_ry(j,nucleus)-PT_ry(k,nucleus))**2
     &       +(PT_rz(j,nucleus)-PT_rz(k,nucleus))**2
          add=exp(-(2.0*gw*drr))
          PT_rho(j,nucleus) = PT_rho(j,nucleus)+add
          PT_rho(k,nucleus) = PT_rho(k,nucleus)+add
15      continue
14    continue

C local neutron density in nucleus A,Z
      do 16 j=Z+1,A
        do 17 k=j+1,A
          drr=(PT_rx(j,nucleus)-PT_rx(k,nucleus))**2
     &       +(PT_ry(j,nucleus)-PT_ry(k,nucleus))**2
     &       +(PT_rz(j,nucleus)-PT_rz(k,nucleus))**2
          add=exp(-(2.0*gw*drr))
          PT_rho(j,nucleus) = PT_rho(j,nucleus)+add
          PT_rho(k,nucleus) = PT_rho(k,nucleus)+add
17      continue
16    continue

      densnorm = ((2.0*gw/pi)**1.5)/int(CTParam(67))
      do 18 j=1,A
        PT_rho(j,nucleus) = PT_rho(j,nucleus)*densnorm
        PT_pmax(j,nucleus) = hqc*(3.0*pi*pi*PT_rho(j,nucleus))**(1./3.)
18    continue

      izcnt=0
      do 12 j=1,A
         P = PT_pmax(j,nucleus)*ranf(0)**(1./3.)
         cost = 1.-2.*ranf(0)
         sint = sqrt(1.-cost**2)
         phi = 2.*Pi*ranf(0)
         PT_px(j,nucleus) = P*sint*cos(phi)
         PT_py(j,nucleus) = P*sint*sin(phi)
         PT_pz(j,nucleus) = P*cost
         pxcm=pxcm+PT_px(j,nucleus)
         pycm=pycm+PT_py(j,nucleus)
         pzcm=pzcm+PT_pz(j,nucleus)
         if (j.le.Z) then
            PT_iso3(j,nucleus)= 1*antia
            PT_charge(j,nucleus)=1*antia
         else
            PT_iso3(j,nucleus)= -(1*antia)
            PT_charge(j,nucleus)=0
         endif

         PT_spin(j,nucleus) = getspin(1,-1)
         PT_ityp(j,nucleus) = 1*antia
        uid_cnt=uid_cnt+1
        PT_uid(j,nucleus)=uid_cnt

12    continue

c perform CM-correction
      pxcm=pxcm/A
      pycm=pycm/A
      pzcm=pzcm/A
      do 2 i=1,A
         PT_px(i,nucleus)=PT_px(i,nucleus)-pxcm
         PT_py(i,nucleus)=PT_py(i,nucleus)-pycm
         PT_pz(i,nucleus)=PT_pz(i,nucleus)-pzcm
c effective masses for initial energy corr. (CTOption(11).eq.0)
         r=sqrt(PT_rx(i,nucleus)**2+PT_ry(i,nucleus)**2
     &         +PT_rz(i,nucleus)**2)
         p=sqrt(PT_px(i,nucleus)**2+PT_py(i,nucleus)**2
     &         +PT_pz(i,nucleus)**2)
         PT_fmass(i,nucleus) = meff(z,a,r,p)
         PT_p0(i,nucleus)=sqrt(PT_px(i,nucleus)**2+PT_py(i,nucleus)**2
     &                    +PT_pz(i,nucleus)**2+PT_fmass(i,nucleus)**2)
 2    continue
c end of CM-correction

      return
      end


      function nucrad(AA)
      implicit none
      real*8 nucrad, r_0
      integer A,AA
      include 'coms.f'
      include 'options.f'

      A=abs(AA)/CTParam(67)

c root mean square radius of nucleus of mass A 
c r_0 corresponding to rho0
      if (CTOption(24).eq.1) then
c root mean square radius of nucleus of mass A (Mayer-Kuckuck)
         nucrad = 1.128 * a**(1./3.) - 0.89 * a**(-(1./3.))
      else
         r_0 = (0.75/pi/rho0)**(1./3.) 
c subtract gaussian tails, for distributing centroids correctly
         nucrad = r_0*(0.5*(a + (a**(1./3.)-1.)**3.))**(1./3.)
      endif

      return
      end

      subroutine boostnuc(i1,i2,pin,b,dst)
      implicit none
      include 'coms.f'
      include 'options.f'
      integer i1,i2,i
      real*8 b,dst,ei,ti
      real*8 pin,beta,gamma

      do 1 i=i1,i2

      beta = pin/sqrt(pin**2+fmass(i)**2)
      gamma = 1.d0/sqrt(1.d0-beta**2)

c  Gallilei-Trafo in x-direction (impact parameter)
c  projectile hits at POSITIVE x
         rx(i) = rx(i) + b
c  distance between nuclei: projectile at NEGATIVE z for dst < 0
         if(CTOption(23).eq.0)then
           ti = r0(i)
           rz(i) = rz(i)/gamma+dst/gamma
         else
          rz(i) = (rz(i) + dst)
         end if


         Ei = p0(i)
         p0(i) = gamma*(p0(i) - beta*pz(i))
         pz(i) = gamma*(pz(i) - beta*Ei) 

 1    continue
      return
      end

      real*8 function meff(znuc,anuc,r,p)
c mean binding energy of a nucleon in a nucleus according to weizaecker
      implicit none
      include 'options.f'
      real*8 av,as,ac,aa,ap,mdef,r,p,e,EMNUC
      integer z,a,znuc,anuc
      parameter (av=0.01587,as=0.01834,ac=0.00071)
      parameter (aa=0.09286,ap=11.46,EMNUC=0.938)
      z=znuc/CTParam(67)
      a=anuc/CTParam(67)
      if(CTOption(11).ne.0.or.a.eq.1)then
        meff=EMNUC
        return
      end if
c...mass defect
      mdef=-(av*A)+as*A**0.66667+ac*z*z/a**0.33333+aa*(z-a/2.)**2/a
c...energy per nucleon = binding energy + nucleon mass 
      e=min(0d0,mdef/a)+EMNUC
      meff=sqrt(e**2-p**2)      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getnucleus(nucleus,offset)
c
c     Revision: 1.0 
c
cinput nucleus : 1=projectile, 2=target
cinput offset  : offset for location of nucleus in particle vectors
c
c output : via common blocks
c 
c This subroutine read in a nucleus which has been initialized
c by {\tt cascinit} and stored in the {\tt PT\_ *(i,nucleus)} arrays.
c The respective nucleus is then rotated randomly in configuration
c and momentum space to yield a new initial state.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'coms.f'
      include 'options.f'
      include 'inputs.f'

c local variables
      real*8 eul1,eul2,eul3,ceul1,seul1,ceul2,seul2,ceul3,seul3
      real*8 vecx0,vecy0,vecz0,vecx1,vecy1,vecz1,vecx2
      real*8 vecy2,vecz2
      integer k,nucleus,offset

c functions
      real*8 ranf

c ***
c *** rotation around Euler-angles
c ***
      if (CTParam(21).gt.0.0) then
         eul1 = 0.0
         eul2 = 0.0
         eul3 = 0.0
      else
         eul1 = ranf(0)*2.0*pi
         eul2 = ranf(0)*2.0*pi  
         eul3 = ranf(0)*2.0*pi
      endif
 
      ceul1 = cos(eul1)
      seul1 = sin(eul1)
      ceul2 = cos(eul2)
      seul2 = sin(eul2)
      ceul3 = cos(eul3)
      seul3 = sin(eul3)

      do 178 k=1,PT_AA(nucleus)

c rotate in configuration space

         vecx0 = PT_rx(k,nucleus)
         vecy0 = PT_ry(k,nucleus)
         vecz0 = PT_rz(k,nucleus)


         vecx1 =   ceul1*vecx0 + seul1*vecy0
         vecy1 = -(seul1*vecx0)+ ceul1*vecy0
         vecz1 =                                     vecz0

         vecx2 =   ceul2*vecx1               + seul2*vecz1
         vecy2 =                       vecy1
         vecz2 = -(seul2*vecx1)              + ceul2*vecz1

         vecx0 =   ceul3*vecx2 + seul3*vecy2
         vecy0 = -(seul3*vecx2)+ ceul3*vecy2
         vecz0 =                                     vecz2

         rx(k+offset) = vecx0 
         ry(k+offset) = vecy0
         rz(k+offset) = vecz0 

c rotate in momentum space

         vecx0 = PT_px(k,nucleus)
         vecy0 = PT_py(k,nucleus)
         vecz0 = PT_pz(k,nucleus)


         vecx1 =   ceul1*vecx0 + seul1*vecy0
         vecy1 = -(seul1*vecx0)+ ceul1*vecy0
         vecz1 =                                     vecz0

         vecx2 =   ceul2*vecx1               + seul2*vecz1
         vecy2 =                       vecy1
         vecz2 = -(seul2*vecx1)              + ceul2*vecz1

         vecx0 =   ceul3*vecx2 + seul3*vecy2
         vecy0 = -(seul3*vecx2) + ceul3*vecy2
         vecz0 =                                     vecz2

         px(k+offset) = vecx0 
         py(k+offset) = vecy0
         pz(k+offset) = vecz0 

c initialize the other quantum numbers

         iso3(k+offset)=PT_iso3(k,nucleus)
         ityp(k+offset)=PT_ityp(k,nucleus)
         uid(k+offset)=PT_uid(k,nucleus)
         spin(k+offset)=PT_spin(k,nucleus)
         dectime(k+offset)=PT_dectime(k,nucleus)
         charge(k+offset)=PT_charge(k,nucleus)
         fmass(k+offset)=PT_fmass(k,nucleus)
         r0(k+offset)=PT_r0(k,nucleus)
         p0(k+offset)=PT_p0(k,nucleus)

 178  continue
      
      return
      end

C####C##1#########2#########3#########4#########5#########6#########7##
      real*8 function rnfWSX(AA,zmin,zmax)
c  yields a $x^n$ distributet value for $x$ between mmin and mmax
cccccCcc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

      implicit none

      real*8 zmin,zmax,ranf,a,rr,yf,z,nucrad
      integer AA
      parameter (a=0.54d0)
      
      rr=nucrad(aa)

      if(zmax.lt.rr)then
        write(6,*)'rnfwsx: maximum radius seems too low'
        stop 137
      end if

 108  continue 
      z=zmin+ranf(0)*(zmax-zmin)
      yf=z*z/((zmax-zmin)**3)*0.5d0/(1d0+exp(z-rr)/a)
      rnfWSX=z 
      if(yf.gt.1d0)stop'rnfWSX: wrong normalisaton:' 
      if(yf.gt.ranf(0))    return
      goto 108                        
       
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine setonshell(i)
c
c     Revision : 1.0
c     This subroutine set particle i on-shell
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      include 'coms.f'
      integer i

      p0(i) = sqrt(px(i)**2+py(i)**2+pz(i)**2+fmass(i)**2)

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine nucfast(IA,JJ)
c
c     Performs very fast initialisation of nuclei
c     Adopted from QGSJET (S. Ostapchenko)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'inputs.f'
      real*8 nucrad

      am=0.545
      rad=nucrad(ia)/am
      cr1=1.+3./rad+6./rad**2+6./rad**3
      cr2=3./rad
      cr3=3./rad+6./rad**2

      DO I=1,IA
 1      zuk=ranf(0)*cr1-1.
        if(zuk.le.0.)then
          tt=rad*(ranf(0)**.3333-1.)
        elseif(zuk.le.cr2 )then
          tt=-log(ranf(0))
        elseif(zuk.lt.cr3 )then
          tt=-log(ranf(0))-log(ranf(0))
        else
          tt=-log(ranf(0))-log(ranf(0))-log(ranf(0))
        endif
        if(ranf(0).gt.1./(1.+exp(-abs(tt))))goto 1
        rim=tt+rad
        z=rim*(2D0*ranf(0)-1D0)
        rim=dsqrt(rim*rim-z*z)

        PT_r0(I,JJ)=0d0
        PT_rz(I,JJ)=Z * AM

 2      s1=2d0*ranf(0)-1d0
        s2=2d0*ranf(0)-1d0
        s3=s1*s1+s2*s2
        if(s3.gt.1d0)goto 2
        s3=dsqrt(s3)
        c=s1/s3
        s=s2/s3
        PT_rx(I,JJ)=rim * C * AM
        PT_ry(I,JJ)=rim * S * AM
      enddo

cc      if(debug.ge.3)then
cc        write (*,*) "nucleons"
cc        do i=1,ia
cc          write (*,'(i3,3g12.6)')i,PT_rx(i,JJ),PT_ry(i,JJ),PT_rz(i,JJ)
cc        enddo
cc        write (*,*)
cc      endif
      return
      END


