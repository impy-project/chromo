
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

#ifdef SIBYLLSP
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
#ifdef SIBYLLSP
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
         ICHG(I) = sign(1, LLIST(I)) * ICHP(ABS(MOD(LLIST(I),10000)))
         IDHEP(I) = ISIB_PID2PDG(MOD(LLIST(I),10000))
         JMOHEP(1,I) = LLIST1(I)
         JMOHEP(2,I) = LLIST1(I)
#ifdef SIBYLLSP
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


C=======================================================================
      SUBROUTINE SIGMA_NUC_NUC(IA,IB,ECM,KINT)
C-----------------------------------------------------------------------
C.  wrapping for SIGMA_NUC in NUCLIB
C...Compute with a montecarlo method the "production"
C.  and "quasi-elastic" cross section for  
C.  a nucleus-nucleus interaction
C.  nucleon - nucleon cross section is taken from 
C.  the table calculated by SIBYLL_INI
C.
C.  INPUT : IA            = mass of projectile nucleus
C.          IB            = mass of target nucleus
C.          ECM          = c.m. energy
C.          KINT            = number  of interactions to generate
C.  OUTPUT : SIGMA (mbarn) = "production" cross section
C.           DSIGMA   "    = error
C.           SIGQE    "    = "quasi-elastic" cross section
C.           DSIGQE   "    = error
C.      in COMMON /NUCNUCSIG/ 
C.           additional output is in the common block  /CPROBAB/
C.           Prob(n_A), Prob(n_B), Prob(n_int)
C..........................................................................
#ifndef SIBYLLSP
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER(I-N)
#endif
      COMMON /NUCNUCSIG/ SIGPROD,DSIGPROD,SIGQE,DSIGQE,IBE,ITG
      DIMENSION SIGDIF(3)
      SAVE
      DATA NDB /0/
      
      DSIGPROD = 0.D0
      DSIGQE = 0.D0

      CALL SIB_SIGMA_HP(1,ECM,SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO)
      CALL SIGMA_MC(IA,IB,SIGINEL,SIGEL,KINT,SIGPROD,DSIGPROD,
     +     SIGQE,DSIGQE)
      IBE = IA
      ITG = IB
      IF(DSIGPROD/SIGPROD.gt.0.1D0)THEN
         IF( NDB.EQ.0 ) 
     +     PRINT*,'SIGMA_NUC_NUC: warning! large error in cross section'
         NDB = 1
      ENDIF
      RETURN
      END
