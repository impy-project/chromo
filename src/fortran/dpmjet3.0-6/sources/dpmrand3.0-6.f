#ifndef CHROMO
************************************************************************
*                                                                      *
*                 7) Random number generator package                   *
*                                                                      *
*    THIS IS A PACKAGE CONTAINING A RANDOM NUMBER GENERATOR AND        *
*    SERVICE ROUTINES.                                                 *
*    THE ALGORITHM IS FROM                                             *
*      'TOWARD A UNVERSAL RANDOM NUMBER GENERATOR'                     *
*      G.MARSAGLIA, A.ZAMAN ;  FSU-SCRI-87-50                          *
*    IMPLEMENTATION BY K. HAHN  DEC. 88,                               *
*    THIS GENERATOR SHOULD NOT DEPEND ON THE HARD WARE ( IF A REAL HAS *
*    AT LEAST 24 SIGNIFICANT BITS IN INTERNAL REPRESENTATION ),        *
*    THE PERIOD IS ABOUT 2**144,                                       *
*    TIME FOR ONE CALL AT IBM-XT IS ABOUT 0.7 MILLISECONDS,            *
*    THE PACKAGE CONTAINS                                              *
*      FUNCTION DT_RNDM(I)                  : GENERATOR                *
*      SUBROUTINE DT_RNDMST(NA1,NA2,NA3,NB4): INITIALIZATION           *
*      SUBROUTINE DT_RNDMIN(U,C,CD,CM,I,J)  : PUT SEED TO GENERATOR    *
*      SUBROUTINE DT_RNDMOU(U,C,CD,CM,I,J)  : TAKE SEED FROM GENERATOR *
*      SUBROUTINE DT_RNDMTE(IO)             : TEST OF GENERATOR        *
*---                                                                   *
*    FUNCTION DT_RNDM(I)                                               *
*       GIVES UNIFORMLY DISTRIBUTED RANDOM NUMBERS  IN (0..1)          *
*       I  - DUMMY VARIABLE, NOT USED                                  *
*    SUBROUTINE DT_RNDMST(NA1,NA2,NA3,NB1)                             *
*       INITIALIZES THE GENERATOR, MUST BE CALLED BEFORE USING DT_RNDM *
*       NA1,NA2,NA3,NB1  - VALUES FOR INITIALIZING THE GENERATOR       *
*                          NA? MUST BE IN 1..178 AND NOT ALL 1         *
*                          12,34,56  ARE THE STANDARD VALUES           *
*                          NB1 MUST BE IN 1..168                       *
*                          78  IS THE STANDARD VALUE                   *
*    SUBROUTINE DT_RNDMIN(U,C,CD,CM,I,J)                               *
*       PUTS SEED TO GENERATOR ( BRINGS GENERATOR IN THE SAME STATUS   *
*       AS AFTER THE LAST DT_RNDMOU CALL )                             *
*       U(97),C,CD,CM,I,J  - SEED VALUES AS TAKEN FROM DT_RNDMOU       *
*    SUBROUTINE DT_RNDMOU(U,C,CD,CM,I,J)                               *
*       TAKES SEED FROM GENERATOR                                      *
*       U(97),C,CD,CM,I,J  - SEED VALUES                               *
*    SUBROUTINE DT_RNDMTE(IO)                                          *
*       TEST OF THE GENERATOR                                          *
*       IO     - DEFINES OUTPUT                                        *
*                  = 0  OUTPUT ONLY IF AN ERROR IS DETECTED            *
*                  = 1  OUTPUT INDEPENDEND ON AN ERROR                 *
*       DT_RNDMTE USES DT_RNDMIN AND DT_RNDMOU TO BRING GENERATOR TO   *
*       SAME STATUS                                                    *
*       AS BEFORE CALL OF DT_RNDMTE                                    *
************************************************************************
*$ CREATE DT_RNDM.FOR
*COPY DT_RNDM
*
*===rndm===============================================================*
*
      DOUBLE PRECISION FUNCTION DT_RNDM(VDUMMY)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

* random number generator
      COMMON /DTRAND/ U(97),C,CD,CM,I,J

* counter of calls to random number generator
* uncomment if needed
C     COMMON /DTRNCT/ IRNCT0,IRNCT1
C     LOGICAL LFIRST
C     DATA LFIRST /.TRUE./

* counter of calls to random number generator
* uncomment if needed
C     IF (LFIRST) THEN
C        IRNCT0 = 0
C        IRNCT1 = 0
C        LFIRST = .FALSE.
C     ENDIF
 100  CONTINUE
      DT_RNDM = U(I)-U(J)
      IF ( DT_RNDM.LT.0.0D0 ) DT_RNDM = DT_RNDM+1.0D0
      U(I) = DT_RNDM
      I    = I-1
      IF ( I.EQ.0 ) I = 97
      J    = J-1
      IF ( J.EQ.0 ) J = 97
      C    = C-CD
      IF ( C.LT.0.0D0 ) C = C+CM
      DT_RNDM = DT_RNDM-C
      IF ( DT_RNDM.LT.0.0D0 ) DT_RNDM = DT_RNDM+1.0D0

      IF ((DT_RNDM.EQ.0.D0).OR.(DT_RNDM.EQ.1.D0)) GOTO 100

* counter of calls to random number generator
* uncomment if needed
C     IRNCT0 = IRNCT0+1

      RETURN
      END

*$ CREATE DT_RNDMST.FOR
*COPY DT_RNDMST
*
*===rndmst=============================================================*
*
      SUBROUTINE DT_RNDMST(NA1,NA2,NA3,NB1)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

* random number generator
      COMMON /DTRAND/ U(97),C,CD,CM,I,J

      MA1 = NA1
      MA2 = NA2
      MA3 = NA3
      MB1 = NB1
      I   = 97
      J   = 33
      DO 20 II2 = 1,97
        S = 0
        T = 0.5D0
        DO 10 II1 = 1,24
          MAT  = MOD(MOD(MA1*MA2,179)*MA3,179)
          MA1  = MA2
          MA2  = MA3
          MA3  = MAT
          MB1  = MOD(53*MB1+1,169)
          IF ( MOD(MB1*MAT,64).GE.32 ) S = S+T
   10   T = 0.5D0*T
   20 U(II2) = S
      C  =   362436.0D0/16777216.0D0
      CD =  7654321.0D0/16777216.0D0
      CM = 16777213.0D0/16777216.0D0
      RETURN
      END

*$ CREATE DT_RNDMIN.FOR
*COPY DT_RNDMIN
*
*===rndmin=============================================================*
*
      SUBROUTINE DT_RNDMIN(UIN,CIN,CDIN,CMIN,IIN,JIN)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

* random number generator
      COMMON /DTRAND/ U(97),C,CD,CM,I,J

      DIMENSION UIN(97)

      DO 10 KKK = 1,97
   10 U(KKK) = UIN(KKK)
      C  = CIN
      CD = CDIN
      CM = CMIN
      I  = IIN
      J  = JIN

      RETURN
      END

*$ CREATE DT_RNDMOU.FOR
*COPY DT_RNDMOU
*
*===rndmou=============================================================*
*
      SUBROUTINE DT_RNDMOU(UOUT,COUT,CDOUT,CMOUT,IOUT,JOUT)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

* random number generator
      COMMON /DTRAND/ U(97),C,CD,CM,I,J

      DIMENSION UOUT(97)

      DO 10 KKK = 1,97
   10 UOUT(KKK) = U(KKK)
      COUT  = C
      CDOUT = CD
      CMOUT = CM
      IOUT  = I
      JOUT  = J

      RETURN
      END

*$ CREATE DT_RNDMTE.FOR
*COPY DT_RNDMTE
*
*===rndmte=============================================================*
*
      SUBROUTINE DT_RNDMTE(IO)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      DIMENSION UU(97),U(6),X(6),D(6)
      DATA U / 6533892.D0, 14220222.D0, 7275067.D0, 6172232.D0,
     +8354498.D0, 10633180.D0/

      CALL DT_RNDMOU(UU,CC,CCD,CCM,II,JJ)
      CALL DT_RNDMST(12,34,56,78)
      DO 10 II1 = 1,20000
   10 XX = DT_RNDM(XX)
      SD        = 0.0D0
      DO 20 II2 = 1,6
        X(II2)  = 4096.D0*(4096.D0*DT_RNDM(SD))
        D(II2)  = X(II2)-U(II2)
   20 SD = SD+D(II2)
      CALL DT_RNDMIN(UU,CC,CCD,CCM,II,JJ)
**sr 24.01.95
C     IF ( IO.EQ. 1.OR. SD.NE.0. 0) WRITE(6,500) (U(I),X(I),D(I),I=1,6)
      IF ((IO.EQ.1).OR.(SD.NE.0.0)) THEN
C        WRITE(6,1000)
 1000    FORMAT(/,/,1X,'DT_RNDMTE: Test of random-number generator...',
     &          ' passed')
      ENDIF
**
      RETURN
  500 FORMAT('  === TEST OF THE RANDOM-GENERATOR ===',/,
     &'    EXPECTED VALUE    CALCULATED VALUE     DIFFERENCE',/, 6(F17.
     &1,F20.1,F15.3,/), '  === END OF TEST ;',
     &'  GENERATOR HAS THE SAME STATUS AS BEFORE CALLING DT_RNDMTE')
      END

*$ CREATE PYR.FOR
*COPY PYR
*
*===pyr================================================================*
*
      DOUBLE PRECISION FUNCTION PYR(IDUMMY)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      DUMMY = DBLE(IDUMMY)
      PYR = DT_RNDM(DUMMY)

      RETURN
      END

*$ CREATE PHO_LHIST.FOR
*COPY PHO_LHIST
*
*===poluhi=============================================================*
*
      SUBROUTINE PHO_LHIST(I,X)
**

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      RETURN
      END

*$ CREATE PDFSET.FOR
*COPY PDFSET
*
C**********************************************************************
C
C   dummy subroutines, remove to link PDFLIB
C
C**********************************************************************
      SUBROUTINE PDFSET(PARAM,VALUE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PARAM(20),VALUE(20)
      CHARACTER*20 PARAM
      END

*$ CREATE STRUCTM.FOR
*COPY STRUCTM
*
      SUBROUTINE STRUCTM(XI,SCALE2,UV,DV,US,DS,SS,CS,BS,TS,GL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      END

*$ CREATE STRUCTP.FOR
*COPY STRUCTP
*
      SUBROUTINE STRUCTP(XI,SCALE2,P2,IP2,UV,DV,US,DS,SS,CS,BS,TS,GL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      END

*
#endif