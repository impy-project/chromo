C***********************************************************************
C
C    interface to PHOJET double precision random number generator 
C    for SIBYLL \FR'14
C
C***********************************************************************
      DOUBLE PRECISION FUNCTION S_RNDM(IDUMMY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DUMMY = dble(IDUMMY)
      S_RNDM= PHO_RNDM(DUMMY)
      END

C***********************************************************************
C
C    initialization routine for double precision random number generator
C    calls PHO_RNDIN \FR'14
C
C***********************************************************************
      SUBROUTINE RND_INI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /RNDMGAS/ ISET
      ISET = 0
      CALL PHO_RNDIN(12,34,56,78)
      END


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
C***********************************************************************
      

      DOUBLE PRECISION FUNCTION PHO_RNDM(DUMMY)
C***********************************************************************
C
C    random number generator
C
C    initialization by call to PHO_RNDIN needed!
C     
C    the algorithm is taken from
C      G.Marsaglia, A.Zaman: 'Toward a unversal random number generator'
C      Florida State Univ. preprint FSU-SCRI-87-70
C
C    implementation by K. Hahn (Dec. 88), changed to include possibility
C    of saving / reading generator registers to / from file (R.E. 10/98)
C
C    generator should not depend on the hardware (if a real has
C    at least 24 significant bits in internal representation),
C    the period is about 2**144,
C
C    internal registers:
C       U(97),C,CD,CM,I,J  - seed values as initialized in PHO_RNDIN
C
C
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      COMMON /PORAND/ U(97),C,CD,CM,I,J

 100  CONTINUE
      RNDMI = DUMMY
      RNDMI = U(I)-U(J)
      IF ( RNDMI.LT.0.D0 ) RNDMI = RNDMI+1.D0
      U(I) = RNDMI
      I    = I-1
      IF ( I.EQ.0 ) I = 97
      J    = J-1
      IF ( J.EQ.0 ) J = 97
      C    = C-CD
      IF ( C.LT.0.D0 ) C = C+CM
      RNDMI = RNDMI-C
      IF ( RNDMI.LT.0.D0 ) RNDMI = RNDMI+1.D0

      IF((ABS(RNDMI).LT.0.D0).OR.(ABS(RNDMI-1.D0).LT.1.D-10)) GOTO 100
      PHO_RNDM = RNDMI

      END


CDECK  ID>, PHO_RNDIN
      SUBROUTINE PHO_RNDIN(NA1,NA2,NA3,NB1)
C***********************************************************************
C
C     initialization of PHO_RNDM, has to be called before using PHO_RNDM
C
C     input:
C       NA1,NA2,NA3,NB1  - values for initializing the generator
C                          NA? must be in 1..178 and not all 1;
C                          12,34,56  are the standard values
C                          NB1 must be in 1..168;
C                          78  is the standard value
C
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      COMMON /PORAND/ U(97),C,CD,CM,I,J
      MA1 = NA1
      MA2 = NA2
      MA3 = NA3
      MB1 = NB1
      I   = 97
      J   = 33
      DO 20 II2 = 1,97
        S = 0.D0
        T = 0.5D0
        DO 10 II1 = 1,24
          MAT  = MOD(MOD(MA1*MA2,179)*MA3,179)
          MA1  = MA2
          MA2  = MA3
          MA3  = MAT
          MB1  = MOD(53*MB1+1,169)
          IF ( MOD(MB1*MAT,64).GE.32 ) S = S+T
          T    = 0.5D0*T
 10     CONTINUE
        U(II2) = S
 20   CONTINUE
      C  =   362436.D0/16777216.D0
      CD =  7654321.D0/16777216.D0
      CM = 16777213.D0/16777216.D0

      END


CDECK  ID>, PHO_RNDSI
      SUBROUTINE PHO_RNDSI(UIN,CIN,CDIN,CMIN,IIN,JIN)
C***********************************************************************
C
C     updates internal random number generator registers using
C     registers given as arguments
C
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      DIMENSION UIN(97)
      COMMON /PORAND/ U(97),C,CD,CM,I,J
      DO 10 KKK = 1,97
        U(KKK) = UIN(KKK)
 10   CONTINUE
      C  = CIN
      CD = CDIN
      CM = CMIN
      I  = IIN
      J  = JIN

      END


CDECK  ID>, PHO_RNDSO
      SUBROUTINE PHO_RNDSO(UOUT,COUT,CDOUT,CMOUT,IOUT,JOUT)
C***********************************************************************
C
C     copies internal registers from randon number generator
C     to arguments
C
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      DIMENSION UOUT(97)
      COMMON /PORAND/ U(97),C,CD,CM,I,J
      DO 10 KKK = 1,97
        UOUT(KKK) = U(KKK)
 10   CONTINUE
      COUT  = C
      CDOUT = CD
      CMOUT = CM
      IOUT  = I
      JOUT  = J

      END


CDECK  ID>, PHO_RNDTE
      SUBROUTINE PHO_RNDTE(IO)
C***********************************************************************
C
C     test of random number generator PHO_RNDM
C
C     input:
C       IO defines output
C           0  output only if an error is detected
C           1  output independend on an error
C
C     uses PHO_RNDSI and PHO_RNDSO to bring the random number generator
C     to same status as it had before the test run
C
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE


C  input/output channels
      INTEGER LI,LO
      COMMON /POINOU/ LI,LO

      DIMENSION UU(97)
      DIMENSION U(6),X(6),D(6)
      DATA U / 6533892.D0 , 14220222.D0 ,  7275067.D0 ,
     &         6172232.D0 ,  8354498.D0 , 10633180.D0 /

      CALL PHO_RNDSO(UU,CC,CCD,CCM,II,JJ)

      CALL PHO_RNDIN(12,34,56,78)
      DO 10 II1 = 1,20000
        XX      = PHO_RNDM(SD)
 10   CONTINUE

      SD        = 0.D0
      DO 20 II2 = 1,6
        X(II2)  = 4096.D0*(4096.D0*PHO_RNDM(XX))
        D(II2)  = X(II2)-U(II2)
        SD      = SD+ABS(D(II2))
 20   CONTINUE

      CALL PHO_RNDSI(UU,CC,CCD,CCM,II,JJ)

      IF ((IO.EQ.1).OR.(ABS(SD).GT.0.D-10)) THEN
        WRITE(LO,50) (U(I),X(I),D(I),I=1,6)
      ENDIF

 50   FORMAT(/,' PHO_RNDTE: test of the random number generator:',/,
     &  '    expected value    calculated value     difference',/,
     &  6(F17.1,F20.1,F15.3,/),
     &  ' generator has the same status as before calling PHO_RNDTE',/)

      END


CDECK  ID>, PHO_RNDST
      SUBROUTINE PHO_RNDST(MODE,FILENA)
C***********************************************************************
C
C     read / write random number generator status from / to file
C
C     input:    MODE        1   read registers from file
C                           2   dump registers to file
C
C               FILENA      file name
C
C***********************************************************************

      IMPLICIT NONE



      SAVE

      INTEGER       MODE
      CHARACTER*(*) FILENA


C  input/output channels
      INTEGER LI,LO
      COMMON /POINOU/ LI,LO


      DOUBLE PRECISION UU,CC,CCD,CCM
      DIMENSION UU(97)

      INTEGER I,II,JJ

      CHARACTER*80 CH_DUMMY

      IF(MODE.EQ.1) THEN

        WRITE(LO,'(/,1X,2A,A,/)') 'PHO_RNDST: ',
     &    'reading random number registers from file ',FILENA

        OPEN(12,FILE=FILENA,ERR=1010,STATUS='OLD')
        READ(12,*,ERR=1010) CH_DUMMY
        DO I=1,97
          READ(12,*,ERR=1010) UU(I)
        ENDDO
        READ(12,*,ERR=1010) CC
        READ(12,*,ERR=1010) CCD
        READ(12,*,ERR=1010) CCM
        READ(12,*,ERR=1010) II,JJ
        CLOSE(12)
        CALL PHO_RNDSI(UU,CC,CCD,CCM,II,JJ)

      ELSE IF(MODE.EQ.2) THEN

        WRITE(LO,'(/,1X,2A,A,/)') 'PHO_RNDST: ',
     &    'dumping random number registers to file ',FILENA

        OPEN(12,FILE=FILENA,ERR=1010,STATUS='UNKNOWN')
        CALL PHO_RNDSO(UU,CC,CCD,CCM,II,JJ)
        WRITE(12,'(1X,A)',ERR=1020) 'random number status registers:'
        DO I=1,97
          WRITE(12,'(1X,1P,E28.20)',ERR=1020) UU(I)
        ENDDO
        WRITE(12,'(1X,1P,E28.20)',ERR=1020) CC
        WRITE(12,'(1X,1P,E28.20)',ERR=1020) CCD
        WRITE(12,'(1X,1P,E28.20)',ERR=1020) CCM
        WRITE(12,'(1X,2I4)',ERR=1020) II,JJ
        CLOSE(12)

      ELSE

        WRITE(LO,'(/,1X,2A,I6,/)') 'PHO_RNDST: ',
     &    'called with invalid mode, nothing done (mode)',MODE

      ENDIF

      RETURN

 1010 CONTINUE
      WRITE(LO,'(1X,2A,A,/)') 'PHO_RNDST: ',
     &  'cannot open or read file ',FILENA
      RETURN

 1020 CONTINUE
      WRITE(LO,'(1X,A,A,/)') 'PHO_RNDST: ',
     &  'cannot open or write file ',FILENA
      RETURN
    
      END

C----------------------------------------
C standard generator
C----------------------------------------
      REAL FUNCTION S_RNDM_std(IDUMMY)
C...Generator  from the LUND montecarlo
C...Purpose: to generate random numbers uniformly distributed between
C...0 and 1, excluding the endpoints.
      COMMON/LUDATR/MRLU(6),RRLU(100)
      SAVE /LUDATR/
      EQUIVALENCE (MRLU1,MRLU(1)),(MRLU2,MRLU(2)),(MRLU3,MRLU(3)),
     &(MRLU4,MRLU(4)),(MRLU5,MRLU(5)),(MRLU6,MRLU(6)),
     &(RRLU98,RRLU(98)),(RRLU99,RRLU(99)),(RRLU00,RRLU(100))
 
C...  Initialize generation from given seed.
      S_RNDM_std = real(idummy)
      IF(MRLU2.EQ.0) THEN
        IF (MRLU1 .EQ. 0)  MRLU1 = 19780503    ! initial seed
        IJ=MOD(MRLU1/30082,31329)
        KL=MOD(MRLU1,30082)
        I=MOD(IJ/177,177)+2
        J=MOD(IJ,177)+2
        K=MOD(KL/169,178)+1
        L=MOD(KL,169)
        DO 110 II=1,97
        S=0.
        T=0.5
        DO 100 JJ=1,24
        M=MOD(MOD(I*J,179)*K,179)
        I=J
        J=K
        K=M
        L=MOD(53*L+1,169)
        IF(MOD(L*M,64).GE.32) S=S+T
        T=0.5*T
  100   CONTINUE
        RRLU(II)=S
  110   CONTINUE
        TWOM24=1.
        DO 120 I24=1,24
        TWOM24=0.5*TWOM24
  120   CONTINUE
        RRLU98=362436.*TWOM24
        RRLU99=7654321.*TWOM24
        RRLU00=16777213.*TWOM24
        MRLU2=1
        MRLU3=0
        MRLU4=97
        MRLU5=33
      ENDIF
 
C...Generate next random number.
  130 RUNI=RRLU(MRLU4)-RRLU(MRLU5)
      IF(RUNI.LT.0.) RUNI=RUNI+1.
      RRLU(MRLU4)=RUNI
      MRLU4=MRLU4-1
      IF(MRLU4.EQ.0) MRLU4=97
      MRLU5=MRLU5-1
      IF(MRLU5.EQ.0) MRLU5=97
      RRLU98=RRLU98-RRLU99
      IF(RRLU98.LT.0.) RRLU98=RRLU98+RRLU00
      RUNI=RUNI-RRLU98
      IF(RUNI.LT.0.) RUNI=RUNI+1.
      IF(RUNI.LE.0.OR.RUNI.GE.1.) GOTO 130
 
C...Update counters. Random number to output.
      MRLU3=MRLU3+1
      IF(MRLU3.EQ.1000000000) THEN
        MRLU2=MRLU2+1
        MRLU3=0
      ENDIF
      S_RNDM_std=RUNI

      END
