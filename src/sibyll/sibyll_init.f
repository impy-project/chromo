    
      SUBROUTINE SIBINI(SEEDIN)

C-----------------------------------------------------------------------
C  SIB(YLL) INI(TIALIZATION)
C
C  FIRST INITIALIZATION OF SIBYLL PROGRAM PACKAGE.
C  THIS SUBROUTINE IS CALLED FROM START.
C  ARGUMENT:
C   SEED   : ANY INTEGER TO BE USED AS RANDOM GENERATOR SEED
C-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER          ASEED(3)
      INTEGER          SEEDIN
      SAVE
C-----------------------------------------------------------------------

C  init the random number generator
      Call INIT_RMMARD(SEEDIN)


      CALL SIBYLL_INI
      CALL SIGMA_INI
      CALL NUC_NUC_INI

      RETURN
      END