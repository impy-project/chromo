
      FUNCTION XM2DIS(XM2MIN,XM2MAX,ALPHA)
C-----------------------------------------------------------------------
C     function that samples mass**2 from (1/M**2)**alpha
C     with alpha <= 1                                             
C     INPUT: Mmin**2 : minimal mass
C            Mmax**2 : maximal mass
C            alpha   : slope                                      \FR'14
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'sib_utl_cmmn.inc'            
      
c     reduced alpha
      ALPHArdc = 2*(ALPHA-ONE)
      AMIN = LOG(XM2MIN)
      AMAX = LOG(XM2MAX)
      ADLT = AMAX-AMIN
      IF(ABS(ALPHArdc).LT.EPS3)THEN
c     alpha = 1
         XRNDM = MAX(S_RNDM(0),EPS10)
         AX = AMIN+ADLT*XRNDM
         XM2DIS = EXP(AX)
      ELSEIF(ALPHArdc.LT.ZERO.and.ALPHA.gt.0)THEN
c     0 < alpha < 1
         XRNDM = MAX(S_RNDM(0),EPS10)
c     AX = AMAX-LOG(XRNDM)*ALPHArdc
         DX = XM2MAX**(ONE-ALPHA)*XRNDM + XM2MIN**(ONE-ALPHA)*(1-XRNDM)
         AX = LOG(DX)/(ONE-ALPHA)
         XM2DIS = EXP(AX)
      ELSEIF(ALPHArdc.GE.ONE)THEN
c     alpha >= 1
         ALPHAr = ONE-ALPHA
         XMINA = XM2MIN**ALPHAr
         XMAXA = XM2MAX**ALPHAr
         XDLT = XMAXA-XMINA
         XRNDM = MAX(S_RNDM(0),EPS10)
         Z = LOG(XMINA+XDLT*XRNDM)/ALPHAR
         XM2DIS = EXP(Z)
      ELSE
         WRITE(6,*) 'M2DIS: undefined exponent in mass distribution!'
         XA = -1.D0
         XA = LOG(XA)
         STOP
      ENDIF
      END
