    
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

      SUBROUTINE SIBHEP3
C-----------------------------------------------------------------------
C  Convert to HEPEVT common block for double precision sibyll
C
C-----------------------------------------------------------------------
      IMPLICIT NONE

      DOUBLE PRECISION P
      INTEGER NP,LLIST,LLIST1,NP_max,NEVSIB
      DATA NEVSIB /0/
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP
      COMMON /S_PLIST1/ LLIST1(8000)
      INTEGER ICHP,ISTR,IBAR
      COMMON /S_CHP/ ICHP(SIBYLL_TABLE_LENGTH), ISTR(SIBYLL_TABLE_LENGTH), IBAR(SIBYLL_TABLE_LENGTH)

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
         PHEP(1,I) = P(I,1)
         PHEP(2,I) = P(I,2)
         PHEP(3,I) = P(I,3)
         PHEP(4,I) = P(I,4)
         PHEP(5,I) = P(I,5)
      END DO

      NEVSIB = NEVSIB + 1
      END

#ifndef IMPY      
      DOUBLE PRECISION FUNCTION GASDEV(Idum)
C***********************************************************************
C     Gaussian deviation
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON /RNDMGAS/ ISET
      SAVE
      DATA ISET/0/      
      gasdev=idum
      IF (ISET.EQ.0) THEN
1       V1=2.D0*S_RNDM(0)-1.D0
        V2=2.D0*S_RNDM(1)-1.D0
        R=V1**2+V2**2
        IF(R.GE.1.D0)GO TO 1
        FAC=SQRT(-2.D0*LOG(R)/R)
        GSET=V1*FAC
        GASDEV=V2*FAC
        ISET=1
      ELSE
        GASDEV=GSET
        ISET=0
      ENDIF
      RETURN
      END
#endif      