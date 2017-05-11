
      SUBROUTINE SAMPLE_SEA (ALPHA,ASUP,XMASS,XMAX,X1,X2,PT)
C-----------------------------------------------------------------------
C.    Routine that samples the kinematical variables of a sea quark pair.
C.  INPUT:  STR_mass_min : minimal string mass ** 2 = x1 * x2 * s
C.          ASUP : large x suppression exponent
C.  OUTPUT:  X1, X2, PT (GeV)                                   /FR'14
C-----------------------------------------------------------------------
Cf2py double precision, intent(in) :: ALPHA,ASUP,XMASS,XMAX
Cf2py double precision, intent(out) :: X1,X2,PT
      IMPLICIT NONE
      SAVE

c     include COMMONs
      INCLUDE 'sib_debug_cmmn.inc'
      include 'sib_nw_prm.inc'
      INCLUDE 'sib_run_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'

c     external type declarations
      DOUBLE PRECISION ALPHA,ASUP,XMASS,XMAX,X1,X2,PT

c     internal types
      DOUBLE PRECISION XMINA,XM2DIS,XR,SLOPE,S_RNDM
      
      IF(ndebug.gt.3)
     &write(lun,*) ' SAMPLE_sea: ALPHA,ASUP,QMASS,XMAX',
     &ALPHA,ASUP,XMASS,XMAX

c     min. momentum fraction for massive quarks
c     i.e. sample from 1/(x+x_min)
      XMINA = TWO*XMASS/SQS
      IF(ndebug.gt.3)
     &     write(lun,*) ' SAMPLE_sea: XMINA:',XMINA
c     exponent of large x suppression: (1-x)**b, b=0 or b>1
      IF(ABS(ASUP).lt.EPS3)THEN
c     b = 0 , no suppression, sample bare 1/(x+xmin)       
         X1 = XM2DIS(XMINA,XMAX,ALPHA) ! ~(1/x)**alpha
         X2 = XM2DIS(XMINA,XMAX,ALPHA) ! ~(1/x)**alpha
         
      ELSEIF(ASUP.ge.EPS3)THEN
c     b >= 1 , sample bare (1-x)**b/(x+xmin)
         SLOPE = MAX(ASUP,EPS3)
c     quark
 100     X1 = XM2DIS(XMINA,XMAX,ALPHA) ! ~(1/x)**alpha
         XR = LOG(ONE-X1)-LOG(ONE-XMINA)
         IF(ndebug.gt.4)
     &write(lun,*) '  X1,XR,SLOPE*XR:',X1,XR,SLOPE*XR
         if(SLOPE*XR.le.LOG(S_RNDM(0))) goto 100

c     anti-quark
 200     X2 = XM2DIS(XMINA,XMAX,ALPHA) ! ~(1/x)**alpha
         XR = log(ONE-X2)-log(ONE-XMINA)
         IF(ndebug.gt.4)
     &write(lun,*) '  X2,XR,SLOPE*XR:',X2,XR,SLOPE*XR
         if(SLOPE*XR.le.log(S_RNDM(0))) goto 200     
      ELSE
         WRITE(LUN,*) ' SAMPLE_sea: suppression exponent out of range.'
         WRITE(LUN,*) ' SAMPLE_sea: ASUP:',ASUP
         STOP
      ENDIF

c     sample pt
c     not yet implemented... to avoid problem with virtual partons
      pt = zero
      IF(ndebug.gt.3)
     &     write(lun,*) ' SAMPLE_sea: X1,X2,PT:',X1,X2,PT

      END
