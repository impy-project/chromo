
      DOUBLE PRECISION FUNCTION CHIDIS (KPARTin, IFL1, IFL2)
C...Generate CHI (fraction of energy of a hadron carried by 
C.                the valence quark, or diquark, as specified by IFL1)
C.  INPUT KPART = code of particle
C.        IFL1, IFL2 = codes of partons (3, 3bar of color)
C.........................................................
      IMPLICIT NONE
      SAVE
c     external types
      INTEGER KPARTIN,IFL1,IFL2
c     COMMONs
      include 'sib_nw_prm.inc'
      INCLUDE 'sib_run_cmmn.inc'
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_cpspl_cmmn.inc'
      INCLUDE 'sib_cutoff_cmmn.inc'
c     internal types
      DOUBLE PRECISION CUT,S_RNDM
      INTEGER KPART,IFQ
      kpart=IABS(kpartin)
      IFQ=IABS(IFL1)
      IF (IFQ.GT.10) IFQ=IABS(IFL2)
      CUT=2.D0*STR_mass_val/SQS
c     hyperon beam cut
      IF(kpart.gt.14) CUT=2.D0*STR_mass_val_hyp/SQS
100   CHIDIS=S_RNDM(0)**2
      if (chidis.lt.cut) goto 100
      if (chidis.gt.(1.D0-cut)) goto 100
      IF((CHIDIS**2/(CHIDIS**2+CUT**2))**0.5D0
     +   *(1.D0-CHIDIS)**CCHIK(IFQ,KPART).LT.S_RNDM(0)) GOTO 100
      CHIDIS = MAX(0.5D0*CUT,CHIDIS)
      CHIDIS = MIN(1.D0-CUT,CHIDIS)
      IF (IABS(IFL1).GT.10)  CHIDIS=1.D0-CHIDIS
      RETURN
      END

      FUNCTION QMASS(IFL)
C-----------------------------------------------------------------------
C...Return quark or diquark constituent masses
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      DIMENSION QMAS(4)
      DATA QMAS /0.325D0,0.325D0,0.5D0,1.5D0/
      IFLA = IABS(IFL)
      IFLA = MOD(IFLA,100)
      IF (IFLA .LE. 4)       THEN
         QMASS = QMAS(IFLA)
      ELSE
         QMA = QMAS(IFLA/10)
         QMB = QMAS(MOD(IFLA,10))
         QMASS = QMA+QMB
      ENDIF
      RETURN
      END
