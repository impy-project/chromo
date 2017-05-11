      SUBROUTINE REMOVE_PI0(XRATE,N1,N2)
C---------------------------------------------------------
C     routine to exchange pi0 on stack with charged pions
C     violating charge conservation.
C     final pions will be off-shell
C      
C     Input: exchange rate and stack positions inbetween
C     which pions shall be exchanged.
C     
C---------------------------------------------------------     
      IMPLICIT NONE
      SAVE
c     Commons
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_plist_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      include 'sib_utl_cmmn.inc'      
C     external types
      DOUBLE PRECISION XRATE
      INTEGER N1,N2
C     internals
      INTEGER I,LL,LA,IFPI0
      DOUBLE PRECISION S_RNDM
      IF(NDEBUG.gt.0)write(lun,*)'REMOVE_PI0: Rate,Mode:',xrate,IPAR(50)
C     select exchange model      
      SELECT CASE(IPAR(50))
C     flat (xf-independent) model
      CASE(1)         
C     stack loop     
         DO I=N1,N2
            LL = MOD(LLIST(I),10000)
            LA = IABS(LL)
c     IF(LA.eq.6)THEN
            IFPI0=(1-MIN(IABS(1-LA/6),1))*MAX(1-MOD(LA,6),0)
c     replace with pi+ or pi-
            LL=LL+IFPI0*(2-INT(MIN((TWO+XRATE)*S_RNDM(LA),THREE-EPS10)))
            LLIST(I) = LL
            IF(NDEBUG.gt.1)
     &           WRITE(LUN,*) 'REMOVE_PI0: LA,IFPI0,LNEW:',LA,IFPI0,LL
         ENDDO
      END SELECT      
      END
