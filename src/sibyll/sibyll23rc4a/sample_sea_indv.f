      subroutine sample_sea_indv(KRMNT,XMINA,XMINA_SEA,NSEA,
     &     XREM0,ALPHA,ASUP,XQMASS,XMAX,XX,IREJ)
      IMPLICIT NONE
      SAVE
Cf2py double precision, intent(in) :: XMINA,XMINA_SEA,XREM0,ALPHA,ASUP,XQMASS,XMAX
Cf2py integer, intent(in) :: Nsea,krmnt
Cf2py double precision, intent(out) :: xx

      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      INCLUDE 'sib_nw_prm.inc'
      INCLUDE 'sib_utl_cmmn.inc'
      INCLUDE 'sib_cnt_cmmn.inc'

      DOUBLE PRECISION XMINA,XMINA_SEA,XREM0,ALPHA,ASUP,XQMASS,XMAX
      INTEGER NSEA,KRMNT
      DOUBLE PRECISION XX
      DIMENSION XX(2*NW_max+2)
      INTEGER IREJ

      DOUBLE PRECISION XREM,XKIN,X1,X2,pt,S_RNDM,XQM
      INTEGER ICNT2,J,jj1,jj2
      DATA ICNT2 /0/
      
      IF(ndebug.gt.2)
     &write(lun,*)'SAMPLE_sea_tot2: called with ',
     &     '(KRMNT,XMINA,XMINA_SEA,NSEA,XREM0,ALPHA,ASUP,XQMASS,XMAX):', 
     &     KRMNT,XMINA,XMINA_SEA,NSEA,XREM0,ALPHA,ASUP,XQMASS,XMAX
      
      XREM = ZERO
      XKIN = 0.1D0
      XQM = XQMASS
      ITRY(4) = 0
      DO WHILE ( XREM .lt. XMINA )
         XREM = XREM0
         IF ( XREM .LT. 2*XMINA + Nsea*XMINA_SEA
     &        +XKIN*(ONEHALF-S_RNDM(Nsea)) ) THEN
            IREJ = 2            ! resample event
            RETURN
         ENDIF
         IF(ITRY(4).gt.Nsea/2*NREJ(4))THEN
            ICNT2 = ICNT2 + 1
            IF(ndebug.gt.2)THEN
               IF(ICNT2.le.5)THEN
                  write(lun,*)'SAMPLE_projectile: rejection!' 
                  write(lun,*)' reached max. no. of trials!', NREJ(4)
                  write(lun,*)' XREM0,N,XMIN:' ,XREM0,Nsea,XMINA_SEA
               ENDIF
               IF(ICNT2.eq.5) 
     &              write(lun,*)' last warning of this type..'
            ENDIF
            IREJ = IPAR(51)
            RETURN
         ENDIF
         DO j=1,Nsea/2
c     scale for interactions other thatn first if Nw>1
            IF(IPAR(75).eq.1.and.J.gt.1) XQM = XQM*PAR(118)
            call sample_sea(ALPHA,ASUP,XQM,XMAX,x1,x2,pt)
            jj1 = 2 + 2*(j-1) + 1
            IF(KRMNT.eq.0) jj1 = 4+2*(j-1) + 1
            jj2 = jj1 + 1
            XX(jj1) = x1
            XX(jj2) = x2
            XREM = XREM - XX(jj1) - XX(jj2)
            IF(NDEBUG.gt.2) 
     &           WRITE(LUN,*) 'x-frac: JW,X3,X4,XREM',
     &           J,XX(jj1),XX(jj2),XREM
         ENDDO
         ITRY(4) = ITRY(4) + 1
         IF(NDEBUG.gt.2) 
     &        WRITE(LUN,*) 'ISMPL,XREM0,XREM,XMINA,XMINSEA',
     &        ITRY(4),XREM0,XREM,XMINA,XMINA_SEA
      ENDDO
      XREM0 = XREM
      IREJ = 0
      END
