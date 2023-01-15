      SUBROUTINE NPYRNG( RVAL )
      IMPLICIT NONE

      EXTERNAL NPYNXT
      DOUBLE PRECISION RVAL
      integer*8 BITGEN
      COMMON /NPY/BITGEN
      SAVE

      CALL NPYNXT(RVAL, BITGEN)

      RETURN
      END

      SUBROUTINE RMMARD( RVEC,LENV,ISEQ )
      implicit none
      DOUBLE PRECISION RVEC(*), RVAL
      INTEGER ISEQ,LENV,IVEC

      DO IVEC = 1, LENV
        call NPYRNG(RVAL)
        RVEC(IVEC) = RVAL
      ENDDO

      RETURN
      END

C=======================================================================

      SUBROUTINE RM48( RVEC,LENV )

C-----------------------------------------------------------------------
C  R(ANDO)M (NUMBER GENERATOR FOR EVAPORATION MODULE/DPMJET)
C
C  THIS SUBROUTINE IS CALLED FROM ROUTINES OF EVAPORATION MODULE.
C  ARGUMENTS:
C   RVEC   = DOUBL.PREC. VECTOR FIELD TO BE FILLED WITH RANDOM NUMBERS
C   LENV   = LENGTH OF VECTOR (# OF RANDNUMBERS TO BE GENERATED)
C-----------------------------------------------------------------------

      implicit none
      DOUBLE PRECISION RVEC(*), RVAL
      INTEGER LENV,IVEC

      DO IVEC = 1, LENV
        call NPYRNG(RVAL)
        RVEC(IVEC) = RVAL
      ENDDO

      RETURN
      END

C=======================================================================

      DOUBLE PRECISION FUNCTION DRANF(dummy)

      implicit none
      double precision dummy

C-----------------------------------------------------------------------
C  (SIM)PLIFIED RANDOM NUMBER GENERATOR CALL TO NPYRNG
C-----------------------------------------------------------------------

      CALL NPYRNG(DRANF)

      RETURN
      END

C=======================================================================

      DOUBLE PRECISION FUNCTION SIMRND()

C-----------------------------------------------------------------------
C  (SIM)PLIFIED RANDOM NUMBER GENERATOR CALL TO NPYRNG
C-----------------------------------------------------------------------

      CALL NPYRNG(SIMRND)

      RETURN
      END

C=======================================================================

C====================       WRAPPERS     ===============================

C-----------------------------------------------------------------------
C  S(IBYLL) R(A)ND(O)M (GENERATOR)
C-----------------------------------------------------------------------
#ifdef SIBYLL_21
#define S_RNDM_RESULT REAL
#else
#define S_RNDM_RESULT DOUBLE PRECISION
#endif
      S_RNDM_RESULT FUNCTION S_RNDM()

      IMPLICIT NONE

      DOUBLE PRECISION SIMRND

#ifdef SIBYLL_21
555   S_RNDM = real(SIMRND())
      IF ((S_RNDM.LE.0E0).OR.(S_RNDM.GE.1E0)) GOTO 555     
#else
      S_RNDM = SIMRND()
#endif

      RETURN
      END

C-----------------------------------------------------------------------
C  PY(THIA) R(ANDOM GENERATOR)
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION PYR()
      IMPLICIT NONE

      CALL NPYRNG(PYR)

      RETURN
      END

      DOUBLE PRECISION FUNCTION RNDM()

C-----------------------------------------------------------------------
C  R(A)ND(O)M (GENERATOR FOR DPMJET)
C-----------------------------------------------------------------------

      IMPLICIT NONE

      CALL NPYRNG(RNDM)

      RETURN
      END

      DOUBLE PRECISION FUNCTION PSRAN()

C-----------------------------------------------------------------------
C  RAN(DOM GENERATOR FOR QGSJET)
C-----------------------------------------------------------------------

      IMPLICIT NONE

      CALL NPYRNG(PSRAN)

      RETURN
      END

      DOUBLE PRECISION FUNCTION RANF()

C-----------------------------------------------------------------------
C  RAN(DOM GENERATOR FOR URQMD
C-----------------------------------------------------------------------

      IMPLICIT NONE

      CALL NPYRNG(RANF)

      RETURN
      END


      DOUBLE PRECISION FUNCTION RLU()

C-----------------------------------------------------------------------
C  RLU  RANDOM GENERATOR FOR JETSET
C-----------------------------------------------------------------------

      IMPLICIT NONE

      CALL NPYRNG(RLU)

      RETURN
      END

      DOUBLE PRECISION FUNCTION DT_RNDM(VDUMMY)

C-----------------------------------------------------------------------
C  THIS FUNCTON IS CALLED FROM DPM_JET306 ROUTINES.
C-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION VDUMMY

      CALL NPYRNG(DT_RNDM)

      RETURN
      END
