      SUBROUTINE RMMAQD( ISEED,ISEQ,CHOPT )

C---------------------------------------------------------------------
C  SUBROUTINE FOR INITIALIZATION OF RMMARD
C  THESE ROUTINE RMMAQD IS A MODIFIED VERSION OF ROUTINE RMMAQ FROM
C  THE CERN LIBRARIES. DESCRIPTION OF ALGORITHM SEE:
C               http://consult.cern.ch/shortwrups/v113/top.html
C  FURTHER DETAILS SEE SUBR. RMMARD
C  ARGUMENTS:
C   ISEED  = SEED TO INITIALIZE A SEQUENCE (3 INTEGERS)
C   ISEQ   = # OF RANDOM SEQUENCE
C   CHOPT  = CHARACTER TO STEER INITIALIZE OPTIONS
C-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER          KSEQ
      PARAMETER        (KSEQ = 8)
      COMMON /CRRANMA3/CD,CINT,CM,TWOM24,TWOM48,MODCNS
      DOUBLE PRECISION CD,CINT,CM,TWOM24,TWOM48
      INTEGER          MODCNS

      COMMON /CRRANMA4/C,U,IJKL,I97,J97,NTOT,NTOT2,JSEQ
      DOUBLE PRECISION C(KSEQ),U(97,KSEQ),UNI
      INTEGER          IJKL(KSEQ),I97(KSEQ),J97(KSEQ),
     *                 NTOT(KSEQ),NTOT2(KSEQ),JSEQ

      DOUBLE PRECISION CC,S,T,UU(97)
      INTEGER          ISEED(3),I,IDUM,II,II97,IJ,IJ97,IORNDM,
     *                 ISEQ,J,JJ,K,KL,L,LOOP2,M,NITER
      CHARACTER        CHOPT*(*), CCHOPT*12
      LOGICAL          FIRST
      SAVE
      DATA             FIRST / .TRUE. /, IORNDM / 11 /, JSEQ / 1 /
C-----------------------------------------------------------------------

      IF ( FIRST ) THEN
        TWOM24 = 2.D0**(-24)
        TWOM48 = 2.D0**(-48)
        CD     = 7654321.D0*TWOM24
        CM     = 16777213.D0*TWOM24
        CINT   = 362436.D0*TWOM24
        MODCNS = 1000000000
        FIRST  = .FALSE.
      ENDIF

      CCHOPT = CHOPT
      IF ( CCHOPT .EQ. ' ' ) THEN
        ISEED(1) = 54217137
        ISEED(2) = 0
        ISEED(3) = 0
        CCHOPT   = 'S'
        JSEQ     = 1
      ENDIF

      IF     ( INDEX(CCHOPT,'S') .NE. 0 ) THEN
        IF ( ISEQ .GT. 0  .AND.  ISEQ .LE. KSEQ ) JSEQ = ISEQ
        IF ( INDEX(CCHOPT,'V') .NE. 0 ) THEN
          READ(IORNDM,'(3Z8)') IJKL(JSEQ),NTOT(JSEQ),NTOT2(JSEQ)
          READ(IORNDM,'(2Z8,Z16)') I97(JSEQ),J97(JSEQ),C(JSEQ)
          READ(IORNDM,'(24(4Z16,/),Z16)') U
          IJ = IJKL(JSEQ)/30082
          KL = IJKL(JSEQ) - 30082 * IJ
          I  = MOD(IJ/177, 177) + 2
          J  = MOD(IJ, 177)     + 2
          K  = MOD(KL/169, 178) + 1
          L  = MOD(KL, 169)
          CD =  7654321.D0 * TWOM24
          CM = 16777213.D0 * TWOM24
        ELSE
          IJKL(JSEQ)  = ISEED(1)
          NTOT(JSEQ)  = ISEED(2)
          NTOT2(JSEQ) = ISEED(3)
          IJ = IJKL(JSEQ) / 30082
          KL = IJKL(JSEQ) - 30082*IJ
          I  = MOD(IJ/177, 177) + 2
          J  = MOD(IJ, 177)     + 2
          K  = MOD(KL/169, 178) + 1
          L  = MOD(KL, 169)
          DO  II = 1, 97
            S = 0.D0
            T = 0.5D0
            DO  JJ = 1, 48
              M = MOD(MOD(I*J,179)*K, 179)
              I = J
              J = K
              K = M
              L = MOD(53*L+1, 169)
              IF ( MOD(L*M,64) .GE. 32 ) S = S + T
              T = 0.5D0 * T
            ENDDO
            UU(II) = S
          ENDDO
          CC    = CINT
          II97  = 97
          IJ97  = 33
C  COMPLETE INITIALIZATION BY SKIPPING (NTOT2*MODCNS+NTOT) RANDOMNUMBERS
          NITER = MODCNS
          DO  LOOP2 = 1, NTOT2(JSEQ)+1
            IF ( LOOP2 .GT. NTOT2(JSEQ) ) NITER = NTOT(JSEQ)
            DO  IDUM = 1, NITER
              UNI = UU(II97) - UU(IJ97)
              IF ( UNI .LT. 0.D0 ) UNI = UNI + 1.D0
              UU(II97) = UNI
              II97     = II97 - 1
              IF ( II97 .EQ. 0 ) II97 = 97
              IJ97     = IJ97 - 1
              IF ( IJ97 .EQ. 0 ) IJ97 = 97
              CC       = CC - CD
              IF ( CC .LT. 0.D0 ) CC  = CC + CM
            ENDDO
          ENDDO
          I97(JSEQ) = II97
          J97(JSEQ) = IJ97
          C(JSEQ)   = CC
          DO  JJ = 1, 97
            U(JJ,JSEQ) = UU(JJ)
          ENDDO
        ENDIF
      ELSEIF ( INDEX(CCHOPT,'R') .NE. 0 ) THEN
        IF ( ISEQ .GT. 0 ) THEN
          JSEQ = ISEQ
        ELSE
          ISEQ = JSEQ
        ENDIF
        IF ( INDEX(CCHOPT,'V') .NE. 0 ) THEN
          WRITE(IORNDM,'(3Z8)') IJKL(JSEQ),NTOT(JSEQ),NTOT2(JSEQ)
          WRITE(IORNDM,'(2Z8,Z16)') I97(JSEQ),J97(JSEQ),C(JSEQ)
          WRITE(IORNDM,'(24(4Z16,/),Z16)') U
        ELSE
          ISEED(1) = IJKL(JSEQ)
          ISEED(2) = NTOT(JSEQ)
          ISEED(3) = NTOT2(JSEQ)
        ENDIF
      ENDIF

      RETURN
      END

*-- Author :    D. HECK IK FZK KARLSRUHE   17/03/2003
C=======================================================================

      SUBROUTINE RMMARD( RVEC,LENV,ISEQ )

C-----------------------------------------------------------------------
C  R(ANDO)M (NUMBER GENERATOR OF) MAR(SAGLIA TYPE) D(OUBLE PRECISION)
C
C  THESE ROUTINES (RMMARD,RMMAQD) ARE MODIFIED VERSIONS OF ROUTINES
C  FROM THE CERN LIBRARIES. DESCRIPTION OF ALGORITHM SEE:
C   http://wwwasdoc.web.cern.ch/wwwasdoc/cernlib.html (v113)
C  IT HAS BEEN CHECKED THAT RESULTS ARE BIT-IDENTICAL WITH CERN
C  DOUBLE PRECISION RANDOM NUMBER GENERATOR RMM48, DESCRIBED IN
C   http://wwwasdoc.web.cern.ch/wwwasdoc/cernlib.html (v116)
C  ARGUMENTS:
C   RVEC   = DOUBLE PREC. VECTOR FIELD TO BE FILLED WITH RANDOM NUMBERS
C   LENV   = LENGTH OF VECTOR (# OF RANDNUMBERS TO BE GENERATED)
C   ISEQ   = # OF RANDOM SEQUENCE
C
C  VERSION OF D. HECK FOR DOUBLE PRECISION RANDOM NUMBERS.
C-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER          KSEQ
      PARAMETER        (KSEQ = 8)
      COMMON /CRRANMA3/CD,CINT,CM,TWOM24,TWOM48,MODCNS
      DOUBLE PRECISION CD,CINT,CM,TWOM24,TWOM48
      INTEGER          MODCNS

      COMMON /CRRANMA4/C,U,IJKL,I97,J97,NTOT,NTOT2,JSEQ
      DOUBLE PRECISION C(KSEQ),U(97,KSEQ),UNI
      INTEGER          IJKL(KSEQ),I97(KSEQ),J97(KSEQ),
     *                 NTOT(KSEQ),NTOT2(KSEQ),JSEQ

      DOUBLE PRECISION RVEC(*)
      INTEGER          ISEQ,IVEC,LENV
      SAVE
C-----------------------------------------------------------------------

      IF ( ISEQ .GT. 0  .AND.  ISEQ .LE. KSEQ ) JSEQ = ISEQ

      DO  IVEC = 1, LENV
        UNI = U(I97(JSEQ),JSEQ) - U(J97(JSEQ),JSEQ)
        IF ( UNI .LT. 0.D0 ) UNI = UNI + 1.D0
        U(I97(JSEQ),JSEQ) = UNI
        I97(JSEQ)  = I97(JSEQ) - 1
        IF ( I97(JSEQ) .EQ. 0 ) I97(JSEQ) = 97
        J97(JSEQ)  = J97(JSEQ) - 1
        IF ( J97(JSEQ) .EQ. 0 ) J97(JSEQ) = 97
        C(JSEQ)    = C(JSEQ) - CD
        IF ( C(JSEQ) .LT. 0.D0 ) C(JSEQ)  = C(JSEQ) + CM
        UNI        = UNI - C(JSEQ)
        IF ( UNI .LT. 0.D0 ) UNI = UNI + 1.D0
C  AN EXACT ZERO HERE IS VERY UNLIKELY, BUT LET''S BE SAFE.
        IF ( UNI .EQ. 0.D0 ) UNI = TWOM48
        RVEC(IVEC) = UNI
      ENDDO

      NTOT(JSEQ) = NTOT(JSEQ) + LENV
      IF ( NTOT(JSEQ) .GE. MODCNS )  THEN
        NTOT2(JSEQ) = NTOT2(JSEQ) + 1
        NTOT(JSEQ)  = NTOT(JSEQ) - MODCNS
      ENDIF

      RETURN
      END

*-- Author :    D. HECK IK FZK KARLSRUHE   18/03/2003
C=======================================================================

      SUBROUTINE RM48( RVEC,LENV )

C-----------------------------------------------------------------------
C  R(ANDO)M (NUMBER GENERATOR FOR EVAPORATION MODULE/DPMJET)
C
C  WE USE HERE A SIMPLIFIED FORM OF RMMARD WITH JSEQ=1.
C  THIS SUBROUTINE IS CALLED FROM ROUTINES OF EVAPORATION MODULE.
C  ARGUMENTS:
C   RVEC   = DOUBL.PREC. VECTOR FIELD TO BE FILLED WITH RANDOM NUMBERS
C   LENV   = LENGTH OF VECTOR (# OF RANDNUMBERS TO BE GENERATED)
C-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION RVEC(*)
      INTEGER          IVEC,LENV
      SAVE

      CALL RMMARD(RVEC,LENV,1)
      RETURN
      END


C-- Author :    D. HECK IK FZK KARLSRUHE   18/03/2003
C=======================================================================

      DOUBLE PRECISION FUNCTION SIMRND()

C-----------------------------------------------------------------------
C  (SIM)PLIFIED RANDOM NUMBER GENERATOR CALL TO RMMARD
C
C  SEE SUBROUT. RMMARD
C  WE USE HERE A SIMPLIFIED FORM OF RMMARD WITH JSEQ=1, LENV=1.
C  THIS FUNCTON IS CALLED FROM SIBYLL ROUTINES.
C-----------------------------------------------------------------------

      DOUBLE PRECISION RVEC(1)
      SAVE

      CALL RMMARD(RVEC,1,1)
      SIMRND = RVEC(1)

      RETURN
      END

C-- Author :    D. HECK IK FZK KARLSRUHE   18/03/2003
C=======================================================================

      DOUBLE PRECISION FUNCTION SIMRND2()

C-----------------------------------------------------------------------
C anfe  (SIM)PLIFIED RANDOM NUMBER GENERATOR CALL TO RMMARD's SEQUENCE 2
C
C  SEE SUBROUT. RMMARD
C-----------------------------------------------------------------------

      DOUBLE PRECISION RVEC(1)
      SAVE

      CALL RMMARD(RVEC,1,2)
      SIMRND2 = RVEC(1)

      RETURN
      END

C====================       WRAPPERS     ===============================
#ifdef SIBYLL_SP
      REAL FUNCTION S_RNDM()

C-----------------------------------------------------------------------
C  S(IBYLL) R(A)ND(O)M (GENERATOR)
C
C  SEE SUBROUT. RMMARD
C  WE USE HERE A SIMPLIFIED FORM OF RMMARD WITH JSEQ=1, LENV=1.
C  THIS FUNCTON IS CALLED FROM SIBYLL ROUTINES.
C-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION SIMRND

      S_RNDM = real(SIMRND())

      RETURN
      END
#else
      DOUBLE PRECISION FUNCTION S_RNDM()

C-----------------------------------------------------------------------
C  S(IBYLL) R(A)ND(O)M (GENERATOR)
C
C  SEE SUBROUT. RMMARD
C  WE USE HERE A SIMPLIFIED FORM OF RMMARD WITH JSEQ=1, LENV=1.
C  THIS FUNCTON IS CALLED FROM SIBYLL ROUTINES.
C-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION SIMRND

      S_RNDM = SIMRND()

      RETURN
      END
#endif


      DOUBLE PRECISION FUNCTION PYR()

C-----------------------------------------------------------------------
C  PY(THIA) R(ANDOM GENERATOR)
C
C  SEE SUBROUT. RMMARD
C  WE USE HERE A SIMPLIFIED FORM OF RMMARD WITH JSEQ=1, LENV=1.
C  THIS FUNCTON IS CALLED FROM SIBYLL ROUTINES.
C-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION SIMRND

      PYR = SIMRND()

      RETURN
      END

      DOUBLE PRECISION FUNCTION RNDM()

C-----------------------------------------------------------------------
C  R(A)ND(O)M (GENERATOR FOR DPMJET)
C
C  SEE SUBROUT. RMMARD
C  WE USE HERE A SIMPLIFIED FORM OF RMMARD WITH JSEQ=1, LENV=1.
C  THIS FUNCTON IS CALLED FROM DPMJET ROUTINES.
C-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION SIMRND

      RNDM = SIMRND()

      RETURN
      END

      DOUBLE PRECISION FUNCTION PSRAN()

C-----------------------------------------------------------------------
C  RAN(DOM GENERATOR FOR QGSJET)
C
C  SEE SUBROUT. RMMARD
C  WE USE HERE A SIMPLIFIED FORM OF RMMARD WITH JSEQ=1, LENV=1.
C  THIS FUNCTON IS CALLED FROM QGSJET01C AND QGSJET_II ROUTINES.
C* ARGUMENT:
C*  B10    =  DUMMY ARGUMENT
C-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION SIMRND

      PSRAN = SIMRND()

      RETURN
      END

      DOUBLE PRECISION FUNCTION RANF()

C-----------------------------------------------------------------------
C  RAN(DOM GENERATOR FOR URQMD
C-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION SIMRND

      RANF = SIMRND()

      RETURN
      END

      DOUBLE PRECISION FUNCTION RLU()

C-----------------------------------------------------------------------
C  RLU  RANDOM GENERATOR FOR JETSET
C-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION SIMRND

      RLU = SIMRND()

      RETURN
      END

      SUBROUTINE INIT_RMMARD(ISEEDIN)

C-----------------------------------------------------------------------
C  Initialization routine for rmmard
C
C* ARGUMENT:
C*  ISEEDIN    =  SEED
C-----------------------------------------------------------------------
      
      IMPLICIT NONE

      INTEGER ISEEDIN, ISEED(3), NSEQ, I
      PARAMETER( NSEQ = 2)

      ISEED(1) = ISEEDIN
      ISEED(2) = 0
      ISEED(3) = 0

      DO I = 1, NSEQ
        CALL RMMAQD(ISEED,I,'S')
        ISEED(1) = ISEED(1) + 1
      END DO

      RETURN
      END