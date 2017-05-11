
      INTEGER FUNCTION IMRG2HAD(IFLB1,IFLB2)
C     -----------------------------------------------------
C     function that merges two flavors into lightest hadron
C     -----------------------------------------------------
      IMPLICIT NONE
c     flavor merging array
      INCLUDE 'sib_kflv_cmmn.inc'
      INTEGER IFLB1,IFLB2,IFLA,IFLB,IFL1,IFL2
      IFLA = IFLB1
      IFLB = IFLB2
c     order by flavor, meson: antiquark-quark, baryon: quark-diquark
      IF(IFLB.lt.IFLA) call iswtch_lmnts(ifla,iflb)
c     if antibaryon switch again..
      IF(IFLB.lt.0) call iswtch_lmnts(ifla,iflb)
      IFL1 = IABS(IFLA)
      IFL2 = IABS(IFLB)
      IMRG2HAD = ISIGN(KFLV(IFL1,IFL2),IFLB)
      END
