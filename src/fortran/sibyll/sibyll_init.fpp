
      SUBROUTINE SIBINI()

C-----------------------------------------------------------------------
C  SIB(YLL) INI(TIALIZATION)
C
C  FIRST INITIALIZATION OF SIBYLL PROGRAM PACKAGE.
C  THIS SUBROUTINE IS CALLED FROM START.
C-----------------------------------------------------------------------

      IMPLICIT NONE

C-----------------------------------------------------------------------

      CALL SIBYLL_INI
      CALL SIGMA_INI
      CALL NUC_NUC_INI

      RETURN
      END

      SUBROUTINE SIBHEP
C-----------------------------------------------------------------------
C  Convert SIBYLL output to HEPEVT common block
C-----------------------------------------------------------------------
      IMPLICIT NONE

#ifdef SIBYLL_21
      REAL P
#else
      DOUBLE PRECISION P
#endif
      INTEGER NP,LLIST,LLIST1,NP_max,NEVSIB
      DATA NEVSIB /0/
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP
      COMMON /S_PLIST1/ LLIST1(8000)
      INTEGER ICHP,ISTR,IBAR
#ifdef SIBYLL_21
      COMMON /S_CHP/ ICHP(49), ISTR(49), IBAR(49)
#else
      COMMON /S_CHP/ ICHP(99), ISTR(99), IBAR(99)
#endif
      INTEGER NEVHEP,NMXHEP,NHEP,ISTHEP,IDHEP,JMOHEP,JDAHEP
      DOUBLE PRECISION PHEP,VHEP
      PARAMETER (NMXHEP=NP_max)
      COMMON /HEPEVT/ NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &                JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),
     &                VHEP(4,NMXHEP)

      INTEGER ICHG
      COMMON /SCHG/  ICHG(NMXHEP)

      INTEGER I, ISIB_PID2PDG
      EXTERNAL ISIB_PID2PDG

      SAVE NEVSIB

      NHEP = NP
      NEVHEP = NEVSIB

      DO I=1,NP
         IF (ABS(LLIST(I)).LT.10000) THEN
            ISTHEP(I) = 1
         ELSE
            ISTHEP(I) = 2
         END IF
         ICHG(I) = ICHP(ABS(LLIST(I)))
         IDHEP(I) = ISIB_PID2PDG(MOD(LLIST(I),10000))
         JMOHEP(1,I) = LLIST1(I)
         JMOHEP(2,I) = LLIST1(I)
#ifdef SIBYLL_21
         PHEP(1,I) = DBLE(P(I,1))
         PHEP(2,I) = DBLE(P(I,2))
         PHEP(3,I) = DBLE(P(I,3))
         PHEP(4,I) = DBLE(P(I,4))
         PHEP(5,I) = DBLE(P(I,5))
#else
         PHEP(1,I) = P(I,1)
         PHEP(2,I) = P(I,2)
         PHEP(3,I) = P(I,3)
         PHEP(4,I) = P(I,4)
         PHEP(5,I) = P(I,5)
#endif
      END DO

      NEVSIB = NEVSIB + 1
      END
