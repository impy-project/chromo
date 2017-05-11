      SUBROUTINE transfOnShell(ECM,XM1in,XM2in,XMAX,IMOD,P1,P2,LBAD)
C-----------------------------------------------------------------------
C     samples 2 --> 2 scattering that puts a particle on its mass shell
C
C     particle1 is along +z, always receives mass
C     particle2 is along -z, mass only sampled if both aquire mass
C
C     DEPENDS: slope-parameter in s_difmass
C
C     INPUT: ECM : center-of-mass energy of scattering particles
C            M1in  : mass of first particle
C            M2in  : mass of second particle
C            XMAX  : maximal mass that can be obtained
C            IMOD  : remnant or diffraction mode
C     
C     OUTPUT: P1,P2 : final state 4vectors in two-particle c.m.   \FR'14
C-----------------------------------------------------------------------
      IMPLICIT NONE
      
c     external types
      DOUBLE PRECISION ECM,XM1in,XM2in,XMAX,P1,P2
      DIMENSION P1(5),P2(5)
      INTEGER IMOD,LBAD
      
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      INCLUDE 'sib_cnt_cmmn.inc'
      INCLUDE 'sib_mass_cmmn.inc'

c     internal types
      DOUBLE PRECISION XMB2,XMT2,AXMX,S,X1,X2,ALX,SLOP0,SLOPE,DB,
     &     T,PTS,PZB2,PZT2,PT,PHI,XMB,XMT,S_RNDM,PTSWTCH
      DOUBLE PRECISION SLOP0_0,ASLOP,BSLOP
      INTEGER II
      DATA SLOP0_0 /6.5D0/        ! b (slope_ for Mx**2 > 5 GeV**2
      DATA ASLOP /31.10362D0/     ! fit to the slope parameter.
      DATA BSLOP /-15.29012D0/
      INCLUDE 'sib_utl_cmmn.inc'

      IF(NDEBUG.gt.3)
     &     WRITE(LUN,*) 'transfOnShell: called with (Ecm,M1,M2,XMAX):',
     &     ECM,XM1in,XM2in,XMAX
      
      XMB2 = XM1in**2
      XMT2 = XM2in**2     

      AXMX = LOG(XMAX)
      
      ITRY(6) = 0
      LBAD = 1

C     remnant pt parameters
c     distribution is: exp(-slope*t)
c     slope = aslop + bslop * log(Mx**2)
c     (by default same as in diff.
c      scale with paramterers 90 and 91)

c     diff. pt paramters
      ASLOP = PAR(133)
      BSLOP = PAR(134)
      SLOP0_0 = PAR(135)
      
      S = ECM*ECM
      X1 = ONE-(XMT2-XMB2)/S
      X2 = TWO-X1
      IF(X2.LT.EPS5) RETURN

 60   ITRY(6) = ITRY(6) + 1
      IF(ITRY(6).GT.NREJ(6)) RETURN
c     sample transverse momentum
      ALX = LOG(MAX(XMT2,XMB2))
c     set slope of pt distribution
      SELECT CASE(IMOD)
c     diffraction dissociation
      CASE(0)
         SLOP0 = SLOP0_0*PAR(93)
         SLOPE = MAX(SLOP0,ASLOP+BSLOP*ALX)
         PTSWTCH = ONE
c     remnant excitation
      CASE(1)
         IF(IPAR(57).eq.0)THEN
            ALX = ALX-LOG(AM2(13))
            SLOP0 = PAR(92)
            DB = (SLOP0-PAR(90))/AXMX
            SLOPE = MAX(SLOP0,PAR(90)+DB*PAR(91)*ALX)
         ELSE
            ALX = ALX-LOG(AM2(13))
            SLOP0 = PAR(92)
            SLOPE = MAX(SLOP0,PAR(90)+PAR(91)*ALX)
         ENDIF
         PTSWTCH = ONE
c     no pt
      CASE(3)
         PTSWTCH = ZERO
         SLOPE = ONE
      END SELECT
      IF(ndebug.gt.3)
     &     WRITE(LUN,*) 'transfOnShell: (SLOP0,SLOPE,log(M**2)):',
     &     SLOP0,SLOPE,ALX
      T = -dLOG(MAX(EPS10,S_RNDM(0)))/SLOPE
      PTS = T*X1*PTSWTCH
      PZB2 = S*HALF*HALF*X1*X1-XMB2-PTS
      PZT2 = S*HALF*HALF*X2*X2-XMT2-PTS
      IF(NDEBUG.gt.3) 
     &     WRITE(LUN,*) 'transfOnShell: (PTS,PZB2,PZT2):',PTS,PZB2,PZT2
c      IF (ABS(PZB2)-PZT2.GT.EPS10) GOTO 60
      IF (PZB2.lt.ZERO.or.PZT2.LT.ZERO) GOTO 60
      PT = dSQRT(PTS)
      PHI = TWOPI*S_RNDM(0)
      XMB = sqrt(XMB2)
      XMT = sqrt(XMT2)
      P2(4) = HALF*ECM*X2
      P2(3) = -dSQRT(PZT2)
      P2(1) = PT*dCOS(PHI)
      P2(2) = PT*dSIN(PHI)
      P2(5) = XMT

      P1(4) = HALF*ECM*X1
      P1(3) = dSQRT(PZB2)
      do ii = 1,2
         P1(ii) = -P2(ii)
      enddo
      P1(5) = XMB
      IF(NDEBUG.gt.3) THEN
          WRITE(LUN,*) 'transfOnShell: (P1):',(p1(ii),ii=1,5)
          WRITE(LUN,*) 'transfOnShell: (P2):',(p2(ii),ii=1,5)
       ENDIF
      LBAD = 0
      END


