      subroutine sig_rpp2014(L,KT,SQS,SLOPE,SIGT,SIGEL,SIGINEL,RHO)
C-----------------------------------------------------------------------
C     implementation of the PDG RPP 2014 cross section fit
C     proton-, pion-, kaon-nucleon interactions
C      
c     projectile dependent parameters are stored in amp array
c     dimensions are: (beam,target,exchange mode)
c     cross section is used for interaction length in AIR
c     therefore proton and neutron cross sections are averaged.
c
C     Input:
c     L : beam id (1: proton, 2: pion, 3: kaon)
c     KT: target id (0: Nucleon, 1: proton, 2: neutron)
c     SQS: c.m. energy in GeV
c     SLOPE: fit does not include elastic slope, need input to calc
c            elastic and inelastic cross section
c     Output:      
c     SIGT,SIGEL,SIGINEL,RHO
c     cross sections and ratio of real and imaginary part of ela. amp.
C-----------------------------------------------------------------------
Cf2py double precision,intent(in) :: sqs,slope
Cf2py integer,intent(in) :: l
Cf2py double precision,intent(out) :: sigt,sigel,siginel,rho      
      implicit none
c     external types
      double precision SQS,SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO
      integer l,kt
c     commons
      include 'sib_debug_cmmn.inc'
      include 'sib_utl_cmmn.inc'
c     internal types
      double precision s,s0,sig,rho1,xi
      integer k,i
C     universal constants and parameters
      DOUBLE PRECISION M0,ETA1,ETA2,H
      DATA M0,ETA1,ETA2,H /2.076D0,0.412D0,0.5626D0,0.2838D0/
      DOUBLE PRECISION AMP(3,2,3)
c     hadron-proton
      DATA (AMP(1,1,i),i=1,3) /33.73, 13.67, 7.77 /
      DATA (AMP(2,1,i),i=1,3) /18.08, 10.44, 1.977 /
      DATA (AMP(3,1,i),i=1,3) /15.84, 5.12, 3.538 /
c     hadron-neutron
      DATA (AMP(1,2,i),i=1,3) /33.77, 14.05, 6.93 /
      DATA (AMP(2,2,i),i=1,3) /18.08, 10.44, 1.977 /
      DATA (AMP(3,2,i),i=1,3) /15.73, 4.81, 1.86 /
c     particle masses
      DOUBLE PRECISION XMA(3),XMB(2)
      DATA XMA /0.93827D0,0.13957D0,0.493667D0/
      DATA XMB /0.93827D0,0.939565D0/
      s = SQS**two
      sigt = zero
      rho = zero
      k = kt
 100  if(kt.eq.0.and.k.lt.2) k = k + 1
      s0=xma(l)+xmb(k)+M0
      s0=s0**two
      xi=s/s0
c     print *,'s,s0,xi',s,s0,xi
c     print *,'eta1,eta2,h,M0',eta1,eta2,h,M0
c     print *,'P,R1,R2',amp(l,k,1),amp(l,k,2),amp(l,k,3)
c     print *,H*log(xi)**two,amp(l,k,1),amp(l,k,2)*(ONE/xi)**eta1,
c     &        amp(l,k,3)*(one/xi)**eta2
      sig = H*log(xi)**two+amp(l,k,1)+amp(l,k,2)*(ONE/xi)**eta1
     &     +amp(l,k,3)*(one/xi)**eta2
c     print *,'sig',sig
c     print *,'pi,half,zero',pi,half,zero
c     print *,pi*h*log(xi),amp(l,k,2)*xi**(-eta1),tan(eta1*pi*half),
c     &        amp(l,k,3)*xi**(-eta2),(tan(pi*eta2*half)+EPS5)
      rho1 = pi*h*log(xi)-amp(l,k,2)*xi**(-eta1)*tan(eta1*pi*half)
     &     +amp(l,k,3)*xi**(-eta2)/(tan(pi*eta2*half)+EPS5)
c     print *,'rho:',rho1
      rho = rho + rho1/sig
      sigt = sigt + sig
c     write(6,*) 'l,k,sig,rho:',l,k,sig,rho
      if(kt.eq.0.and.k.lt.2) goto 100
      if(kt.eq.0) then
         sigt = sigt*half
         rho = rho*half
      endif
c     derive elastic and inelastic cross section
      sigel = sigt**2*(one+rho**2)/(16.D0*pi*slope*cmbarn)
      siginel = sigt-sigel
      IF(ndebug.gt.2)
     &  write(6,*) 'SIG_RPP2014: L,KT,SQS,SIGT,SIGEL,SIGINEL,SLOPE,RHO',
     &     L,KT,SQS,SIGT,SIGEL,SIGINEL,SLOPE,RHO
      end
