      DOUBLE PRECISION FUNCTION FERMI(XARG,X0,XALPH)
C-----------------------------------------------------------------------
C     fermi function, used to smoothen samplings
C     f = 1/(1+exp((x-x0)/alpha))
C-----------------------------------------------------------------------
      IMPLICIT NONE
c     externals
      DOUBLE PRECISION XARG,X0,XALPH
c     COMMONs
      INCLUDE 'sib_utl_cmmn.inc'      
c     internals
      fermi=one+exp((xarg-x0)/xalph)
      fermi=ONE/fermi
      END
      
