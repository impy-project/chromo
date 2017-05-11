      SUBROUTINE PTDIS_4FLV (IFL,PX,PY)
C...Generate pT
Cf2py double precision,intent(out) :: px
Cf2py double precision,intent(out) :: py
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
c      COMMON /S_CQDIS2/ PPT02(44)
      include 'sib_cqdis2_cmmn.inc'
      include 'sib_cflafr_cmmn.inc'
      include 'sib_utl_cmmn.inc'
c      DATA ZERO,HALF,ONE,ONEHALF,TWO,THREE 
c     &     /0.D0,0.5D0,1.D0,1.5D0,2.D0,3.D0/
c      DATA TWOPI /6.2831853D0/
c      DATA EPS10 /1D-10/

      IF(IFL.eq.0)THEN
c     quark confinement pt
         PPTT = PAR(110)
         XM = 0.325
         XM2 = XM**2
         RNDM = MAX(EPS10,S_RNDM(IFL))
         XMT = PPTT * LOG(RNDM) - XM
         XMT2 = XMT**2
         PT = SQRT(XMT2-XM2)         
      ELSE
         IFLA = IABS(IFL)
         IFLA = MOD(IFLA,100)
         PPTT = PPT02(IFLA)
c     Gaussian distribution
         PT = PPTT*SQRT(-LOG(MAX(EPS10,S_RNDM(IFL))))
         IF (IPAR(3).GE.1) THEN
            IF(MOD(IFLA,10).NE.0) THEN
               XM = QMASS(IFL)
            ELSE
               XM = HALF        ! pomeron mass
               IF(IPAR(3).ge.6) XM = ZERO
            ENDIF
c     exponential transverse mass
            XM2 = XM**2
            RNDM = MAX(EPS10,S_RNDM(IFL))
            XMT = PPTT * LOG(RNDM) - XM
            XMT2 = XMT**2
            PT = SQRT(XMT2-XM2)
         ENDIF      
      ENDIF
      PHI= TWOPI*S_RNDM(IFL)
      PX=PT*COS(PHI)
      PY=PT*SIN(PHI)
      RETURN
      END


      SUBROUTINE PTSETUP_4FLV(ECM)
C     moved from sib_ndiff to seperate subroutine 
c     so that changes will affect diff. /FR'13
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      INCLUDE 'sib_cqdis2_cmmn.inc'

      SQS = ECM

c     NA22 piC retune
      PTU=.3+.08*dlog10(sqs/30.D0)
      PTS=.45+.08*dlog10(sqs/30.D0)
      PTQQ=.6+.08*dlog10(sqs/30.D0)
      PTPOM= .6+.08*dlog10(sqs/30.D0)
      if ( IPAR(3).eq.1 ) then 
c     pt0
         ptu=.15+.007*dlog10(sqs/20.D0)**2
         pts=.3+.007*dlog10(sqs/20.D0)**2
         ptqq=.3+.03*dlog10(sqs/20.D0)**2
         ptpom= .6+.08*dlog10(sqs/30.D0)
      elseif ( ipar(3).eq.2 ) then
C     pt1
         ptu=.15+.007*dlog10(sqs/20.D0)**2
         pts=.32+.007*dlog10(sqs/20.D0)**2
         ptqq=.4+.007*dlog10(sqs/20.D0)**2
         ptpom= .6+.08*dlog10(sqs/30.D0)
c     pt2
      elseif ( ipar(3).eq.3 ) then
         ptu=.17+.007*dlog10(sqs/20.D0)**2
         pts=.3+.007*dlog10(sqs/20.D0)**2
         ptqq=.3+.03*dlog10(sqs/20.D0)**2
         ptpom = .6+.08*dlog10(sqs/30.D0)
      elseif ( ipar(3).eq.5 ) then
         PTU=.16+.007*dlog10(sqs/20.D0)**2
         PTS=.28+.007*dlog10(sqs/20.D0)**2
         PTQQ= .3+.03*dlog10(sqs/20.D0)**2
         PTPOM = .23+.03*dlog10(sqs/20.D0)**2
      elseif ( IPAR(3).eq.6 ) then
         PTU=.16+.007*dlog10(sqs/20.D0)**2
         PTS=.28+.007*dlog10(sqs/20.D0)**2
         PTQQ= .3+.03*dlog10(sqs/20.D0)**2
         PTPOM = .23+.03*dlog10(sqs/20.D0)**2
      elseif ( IPAR(3).eq.7 ) then
         PTU= PAR(46) + .007*dlog10(sqs/20.D0)**2
         PTS= PAR(47) + .007*dlog10(sqs/20.D0)**2
         PTQQ= PAR(48) + .03*dlog10(sqs/20.D0)**2
         PTPOM = PAR(49) + .03*dlog10(sqs/20.D0)**2
      elseif ( IPAR(3).eq.8 ) then
         ASQS = MAX(log10(SQS/PAR(109)),ZERO)
         PTU= PAR(46) + PAR(68)*ASQS**2
         PTS= PAR(47) + PAR(70)*ASQS**2
         PTQQ= PAR(48) + PAR(69)*ASQS**2
         PTPOM = PAR(49) + PAR(51)*ASQS**2
         PTSEA = PAR(67) + PAR(52)*ASQS**2
      endif
      PPT02 (1) = PTU
      PPT02 (2) = PTU
      PPT02 (3) = PTS
c     valence pt
      PPT02 (10) = PTPOM
      DO J=11,33
         PPT02(J) = PTQQ
      ENDDO
c     soft minijet pt
      PPT02 (20) = PTSEA
c     sea quark pt
      PPT02 (30) = PAR(132)
c     charm pt
      IF(IPAR(16).eq.1)THEN
         PTCHM=0.5+.08*dlog10(sqs/30.D0)
         PTCHB=0.5+.08*dlog10(sqs/30.D0)
      ELSEIF(IPAR(16).eq.2)THEN
         PTCHM=0.3+.08*dlog10(sqs/30.D0)
         PTCHB=0.5+.08*dlog10(sqs/30.D0)
      ELSEIF(IPAR(16).eq.3)THEN
         PTCHM=0.3+.5*dlog10(sqs/30.D0)
         PTCHB=0.5+.5*dlog10(sqs/30.D0)
      ELSEIF(IPAR(16).eq.4)THEN
         PTCHM=0.365+.473*dlog10(sqs/30.D0)
         PTCHB=0.144+.473*dlog10(sqs/30.D0)
      ELSEIF(IPAR(16).eq.5)THEN
         PTCHM=0.3+.1*dlog10(sqs/30.D0)
         PTCHB=0.5+.1*dlog10(sqs/30.D0)
      ELSEIF(IPAR(16).eq.6)THEN
         PTCHM=0.303+.125*dlog10(sqs/30.D0)
         PTCHB=0.5+.125*dlog10(sqs/30.D0)
      ELSEIF(IPAR(16).eq.7)THEN
         PTCHM=0.308+.165*dlog10(sqs/30.D0)
         PTCHB=0.5+.165*dlog10(sqs/30.D0)
      ELSE
         PTCHM=1.D0+.08*dlog10(sqs/30.D0)
         PTCHB=1.5D0+.08*dlog10(sqs/30.D0)
      ENDIF
      PPT02(4) = PTCHM
      PPT02(14) = PTCHB
      PPT02(24) = PTCHB
      DO J=34,44
         PPT02(J) = PTCHB
      ENDDO
     
      IF(ndebug.gt.2)THEN
         WRITE(LUN,*)'PTSETUP_4FLV: (sqs,(u,d),s,diq,pom,cm,cb)',sqs
     +        ,ppt02(1),ppt02(3),ppt02(11), ppt02(10),ppt02(4),ppt02(34)
      ENDIF

      RETURN
      END
