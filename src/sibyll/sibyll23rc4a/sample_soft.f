
      SUBROUTINE SAMPLE_soft (STR_mass_min, X1,X2,PT)
C-----------------------------------------------------------------------
C...  Routine for the sampling the kinematical variables of sea quarks
C.     according to (1-x)**b / x**2
C.  INPUT:  STR_mass_min : minimal string mass ** 2 = x1 * x2 * s
C.          SLOPE : large x suppression exponent
C.  OUTPUT:  gluon 4momenta on parton stack (GeV)                /FR'14
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      include 'sib_nw_prm.inc'
      INCLUDE 'sib_run_cmmn.inc'
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_cqdis2_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
c$$$      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
c$$$      COMMON /S_CQDIS/ PPT0(35),ptflag
c$$$      COMMON /S_CQDIS2/ PPT02(44)
      INCLUDE 'sib_utl_cmmn.inc'      

      SLOPE = max(ONE,PAR(42))
      ZSOF = TWO*dLOG(STR_mass_min/SQS) ! minim. mass ~ x1 * x2
 50   XMIN = dEXP(ZSOF)
      axmin = ONE/xmin
 100  Z1 = -ONE*dLOG(axmin-(axmin-ONE)*s_rndm(0))
      x1 = dexp(z1)
      XR = dlog(ONE-X1) - dlog(ONE-xmin)
      if(SLOPE*XR.le.log(S_RNDM(0))) goto 100

 200  Z2 = -ONE*dLOG(axmin-(axmin-ONE)*s_rndm(0))
      X2 = dEXP(Z2)
      XR = dlog(ONE-X2) - dlog(ONE-dEXP(ZSOF))
      if(SLOPE*XR.le.dlog(S_RNDM(0))) goto 200     

      IF(Z1+Z2.LE.ZSOF) GOTO 50
      STR_mass2 = dsqrt(X1*X2*S)/TWO
      PPTT = PPT02(10)
 150  PT = PPTT*dSQRT(-dLOG(MAX(EPS10,S_RNDM(0))))
      IF(IPAR(3).eq.6)THEN
         XM = ZERO
         XM2 = XM**2
         RNDM = MAX(EPS10,S_RNDM(IFL))
         XMT = PPTT * dLOG(RNDM) - XM
         XMT2 = XMT**2
         PT = dSQRT(XMT2-XM2)
      ENDIF
      IF(PT.GT.PTmin) GOTO 150
      IF(PT.GE.STR_mass2) GOTO 150
      END


      SUBROUTINE SAMPLE_soft2 (STR_mass_min, X1,X2,PT)
C-----------------------------------------------------------------------
C...Routine for sampling the kinematical variables
C.  that characterize a soft cut pomeron (x1,x2, pT)
C.  from the differential cross section:
C.     d3sigma/(dx1 dx2 dpT)
C.      ~ 1/x_i**a .*. exp(-mT)
C.  INPUT: STR_mass_min : minimal string mass defined by kinematic limits
C.                        of the string fragmentation
C.  PAR:   PAR(42) : exponent a
C.  OUTPUT:  X1, X2, PT (GeV)
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      include 'sib_nw_prm.inc'
      INCLUDE 'sib_run_cmmn.inc'
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_cqdis2_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
c$$$      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
c$$$      COMMON /S_CQDIS/ PPT0(35),ptflag
c$$$      COMMON /S_CQDIS2/ PPT02(44)
      INCLUDE 'sib_utl_cmmn.inc'      

      SLOPE = PAR(42)
      ZSOF = TWO*dLOG(STR_mass_min/SQS) ! zmin
      zsof = zsof * slope
 100  Z1=ONE/SLOPE*(-zsof*s_rndm(0)+zsof)
      Z2=ONE/SLOPE*(-zsof*s_rndm(0)+zsof)
c      print *,'zsof,z1,z2',zsof,z1,z2
      IF(Z1+Z2.LE.ZSOF) GOTO 100
      X1=dEXP(Z1)
      X2=dEXP(Z2)
      STR_mass2 = sqrt(X1*X2*S)/TWO
      if(str_mass2.lt.0.9D0)goto 100
      PPTT = PPT02(10)
c      print *,'ptmin,str_mass:',ptmin,str_mass2
 150  PT = PPTT*dSQRT(-dLOG(MAX(EPS10,S_RNDM(0))))
      IF(IPAR(3).eq.6)THEN
         XM = ZERO
         XM2 = XM**2
         RNDM = MAX(EPS10,S_RNDM(IFL))
         XMT = PPTT * dLOG(RNDM) - XM
         XMT2 = XMT**2
         PT = dSQRT(XMT2-XM2)
      ENDIF
      IF(PT.GT.PTmin) GOTO 150
      IF(PT.GE.STR_mass2) GOTO 150
      PHI = TWOPI*S_RNDM(L)
      END

      SUBROUTINE SAMPLE_soft3 (STR_mass_min, X1,X2,PT)
C-----------------------------------------------------------------------
C...Routine for the sampling the kinematical variables
C.  that characterize a soft cut pomeron (x1,x2, pT)
C.  from the differential cross section:
C.     d3sigma/(dx1 dx2 dpT)
C.  INPUT:  L=1 incident proton, L=2  incident pi
C.          (soft strings identical for pi and p interactions)
C.  OUTPUT:  X1, X2, PT (GeV)
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      include 'sib_nw_prm.inc'
      INCLUDE 'sib_run_cmmn.inc'
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_cqdis2_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
c$$$      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
c$$$      COMMON /S_CQDIS/ PPT0(35),ptflag
c$$$      COMMON /S_CQDIS2/ PPT02(44)
      INCLUDE 'sib_utl_cmmn.inc'      

      SLOPE = max(ONE,PAR(42))
      ZSOF = TWO*dLOG(STR_mass_min/SQS) ! minim. mass ~ x1 * x2
 100  Z1=-ZSOF*S_RNDM(0)+ZSOF   ! sample envelope 1/x
      X1 = dEXP(Z1)
c      print *,'z1,x1:',z1,x1
      XR = dlog(ONE-X1) - dlog(ONE-dEXP(ZSOF))
c      print *,'ratio:',(1.-X1)/(1.-EXP(ZSOF)),(1.-X1),1.-EXP(ZSOF)
c      print *,'log ratio:',xr,log(1.-X1),log(1.-EXP(ZSOF))
      if(SLOPE*XR.le.dlog(S_RNDM(0))) goto 100

 200  Z2=-ZSOF*S_RNDM(0)+ZSOF   ! sample envelope 1/x
      X2 = dEXP(Z2)
      XR = dlog(ONE-X2) - dlog(ONE-dEXP(ZSOF))
      if(SLOPE*XR.le.dlog(S_RNDM(0))) goto 200     
c      print *,'zsof,z1,z2',zsof,z1,z2
      IF(Z1+Z2.LE.ZSOF) GOTO 100
      STR_mass2 = sqrt(X1*X2*S)/TWO
      PPTT = PPT02(10)
      IF(IPAR(3).eq.8) PPTT = PPT02(20)
 150  PT = PPTT*dSQRT(-dLOG(MAX(EPS10,S_RNDM(0))))
      IF(IPAR(3).ge.6)THEN
         XM = ZERO
         XM2 = XM**2
         RNDM = MAX(EPS10,S_RNDM(IFL))
         XMT = PPTT * dLOG(RNDM) - XM
         XMT2 = XMT**2
         PT = dSQRT(XMT2-XM2)
      ENDIF
      IF(PT.GT.PTmin) GOTO 150
      IF(PT.GE.STR_mass2) GOTO 150
      PHI = TWOPI*S_RNDM(L)
      END

      SUBROUTINE SAMPLE_soft5 (STR_mass_min, X1,X2,PT)
C-----------------------------------------------------------------------
C...Routine for the sampling the kinematical variables of sea quarks
C.     according to (1-x)**b / x**2
C.  INPUT:  STR_mass_min : minimal string mass ** 2 = x1 * x2 * s
C.          SLOPE : large x suppression exponent
C.  OUTPUT:  X1, X2, PT (GeV)                                   /FR'14
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      include 'sib_nw_prm.inc'
c      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, kb ,kt
      INCLUDE 'sib_run_cmmn.inc'
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_cqdis2_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
c$$$      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
c$$$      COMMON /S_CQDIS/ PPT0(35),ptflag
c$$$      COMMON /S_CQDIS2/ PPT02(44)
      INCLUDE 'sib_utl_cmmn.inc'      

      SLOPE = max(ONE,PAR(42))
      ZSOF = TWO*dLOG(STR_mass_min/SQS) ! minim. mass ~ x1 * x2
 50   XMIN = dEXP(ZSOF)
      axmin = ONE/xmin
 100  Z1 = -ONE*dLOG(axmin-(axmin-ONE)*s_rndm(0))
      x1 = dexp(z1)
      XR = dlog(ONE-X1) - dlog(ONE-xmin)
      if(SLOPE*XR.le.log(S_RNDM(0))) goto 100

 200  Z2 = -ONE*dLOG(axmin-(axmin-ONE)*s_rndm(0))
      X2 = dEXP(Z2)
      XR = dlog(ONE-X2) - dlog(ONE-dEXP(ZSOF))
      if(SLOPE*XR.le.dlog(S_RNDM(0))) goto 200     

      IF(Z1+Z2.LE.ZSOF) GOTO 50
      STR_mass2 = dsqrt(X1*X2*S)/TWO
      PPTT = PPT02(10)
      IF(IPAR(3).eq.8) PPTT = PPT02(20)
 150  PT = PPTT*dSQRT(-dLOG(MAX(EPS10,S_RNDM(0))))
      IF(IPAR(3).ge.6)THEN
         XM = ZERO
         XM2 = XM**2
         RNDM = MAX(EPS10,S_RNDM(IFL))
         XMT = PPTT * dLOG(RNDM) - XM
         XMT2 = XMT**2
         PT = dSQRT(XMT2-XM2)
      ENDIF
      IF(PT.GT.PTmin) GOTO 150
      IF(PT.GE.STR_mass2) GOTO 150
      END


      SUBROUTINE SAMPLE_soft6 (STR_mass_min, X1,X2,PT)
C-----------------------------------------------------------------------
C...Routine for the sampling the kinematical variables of sea quarks
C.     according to (1-x)**b / x
C.  INPUT:  STR_mass_min : minimal string mass ** 2 = x1 * x2 * s
C.          SLOPE : large x suppression exponent
C.  OUTPUT:  X1, X2, PT (GeV)                                   /FR'14
C-----------------------------------------------------------------------
Cf2py double precision, intent(in) :: STR_mass_min
Cf2py double precision, intent(out) :: X1
Cf2py double precision, intent(out) :: X2
Cf2py double precision, intent(out) :: PT

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      
      INCLUDE 'sib_debug_cmmn.inc'
      include 'sib_nw_prm.inc'
      INCLUDE 'sib_run_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      include 'sib_cqdis2_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'
      
      NOSLOPE = 0
      SLOPE = PAR(42)
      IF(SLOPE.lt.HALF) NOSLOPE = 1
      XMAX = 0.8D0
      ZSOF = TWO*LOG(STR_mass_min/SQS) ! minim. mass ~ x1 * x2
      XMINA = MAX(EXP(ZSOF),EPS10)
      AXMINA = ONE/XMINA
      IF(ndebug.gt.2)
     &write(lun,*) ' SAMPLE_soft: Mmin,ZSOF,XMINA,XMAX,SLOPE:',
     &     STR_mass_min,ZSOF,XMINA,XMAX,SLOPE
      
 100  X1 = XM2DIS(XMINA,XMAX,ONE) ! ~(1/x)**alpha
      IF(NOSLOPE.eq.1) goto 200
      XR = LOG(ONE-X1)-LOG(ONE-XMINA)
      IF(ndebug.gt.5)
     &write(lun,*) '  X1,XR,SLOPE*XR:',X1,XR,SLOPE*XR
      if(SLOPE*XR.le.LOG(S_RNDM(0))) goto 100

 200  X2 = XM2DIS(XMINA,XMAX,ONE) ! ~(1/x)**alpha
      IF(NOSLOPE.eq.1) goto 300
      XR = log(ONE-X2) - log(ONE-XMINA)
      IF(ndebug.gt.5)
     &write(lun,*) '  X2,XR,SLOPE*XR:',X2,XR,SLOPE*XR
      if(SLOPE*XR.le.log(S_RNDM(0))) goto 200

 300  Z1 = log(X1)
      Z2 = log(X2)
      IF(Z1+Z2.LE.ZSOF) GOTO 100     
      STR_mass2 = sqrt(X1*X2*S)/TWO
      PPTT = PPT02(10)
      IF(IPAR(3).eq.8) PPTT = PPT02(20)
      IF(ndebug.gt.2)
     &write(lun,*) ' SAMPLE_soft: PPTT,Mmin2,PTmin:',
     &PPTT,STR_mass2,PTmin
 150  PT = PPTT*SQRT(-LOG(MAX(EPS10,S_RNDM(0))))
      IF(IPAR(3).ge.6)THEN
         XM = ZERO
         XM2 = XM**2
         RNDM = MAX(EPS10,S_RNDM(0))
         XMT = PPTT * LOG(RNDM) - XM
         XMT2 = XMT**2
         PT = SQRT(XMT2-XM2)
      ENDIF
      IF(ndebug.gt.2)
     &write(lun,*) ' XM,XMT2,PT:',XM,XMT2,PT
      IF(PT.GT.PTmin) GOTO 150
      IF(PT.GE.STR_mass2) GOTO 150
      END
