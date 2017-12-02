      SUBROUTINE NEWHIS(XLIM1,XLIM2,XLIM3,XLIMB,IBIN,IREFN)
C**********************************************************************
C
C     create new histogram
C
C     input:  XLIM1        lower limit
C             XLIM2        upper limit
C             XLIM3        bin size
C             XLIMB(*)     field to handle different bin sizes
C             IBIN         number of bins
C                          neg. number to delete all histograms
C
C     output: IREFN        reference number of histogram
C                          (failure: neg. number)
C
C**********************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
C
      DIMENSION XLIMB(*)
C
      PARAMETER (NHIS=200, NDIM=201)
      DOUBLE PRECISION HIST
      COMMON /HISTOS/ XLIMIT(3,NHIS),ISWI(NHIS),ENTRY(NHIS),
     &                HIST(5,NHIS,NDIM),HEVENT(NHIS),OVERF(NHIS),
     &                UNDERF(NHIS),IBINS(NHIS),
     &                IHMIN(NHIS),IHMAX(NHIS),IHISL
C
      PARAMETER (ZERO   =  0.,
     &           IZERO  =  0,
     &           ONE    =  1.,
     &           TWO    =  2.,
     &           OHALF  =  0.5,
     &           TINY   =  1.e-10)
C
C  delete all histograms
      IF(IBIN.LE.0) THEN
        IHISL = 0
        RETURN
      ENDIF
C
      IHIS = IHISL+1
      IF(IHIS.LE.NHIS) THEN
C  general initialization
        DO 75 I=1,NDIM
          HIST(2,IHIS,I) = ZERO
          HIST(3,IHIS,I) = ZERO
          HIST(4,IHIS,I) = ZERO
          HIST(5,IHIS,I) = ZERO
 75     CONTINUE
        ENTRY(IHIS) = ZERO
        HEVENT(IHIS) = ZERO
        OVERF(IHIS) = ZERO
        UNDERF(IHIS) = ZERO
C  field dimension check
        IF(IBIN.GE.NDIM) THEN
          WRITE(6,'(/1X,A,3I6)')
     &      'INIHIS:ERROR: too many bins (bins,limit,histgram ID)',
     &             IBIN,NDIM,IHIS
          STOP
        ENDIF
        IBINS(IHIS) = IBIN
        IHMIN(IHIS) = IBIN
        IHMAX(IHIS) = 0
        IREFN = IHIS
C  bin size and lower limit given by user
        IF(XLIM3.GT.ZERO) THEN
          XLIMIT(1,IHIS) = XLIM1
          XLIMIT(2,IHIS) = XLIM1+IBINS(IHIS)*XLIM3
          XLIMIT(3,IHIS) = XLIM3
          TMP = XLIM1-XLIM3
          DO 100 K=1,IBINS(IHIS)+1
            TMP = TMP + XLIM3
            HIST(1,IHIS,K) = TMP
 100      CONTINUE
          ISWI(IHIS) = 1
C  limits given by user
        ELSE IF(XLIM3.EQ.ZERO) THEN
          XLIMIT(1,IHIS) = XLIM1
          XLIMIT(2,IHIS) = XLIM2
          XLIMIT(3,IHIS) = (XLIM2-XLIM1)/DBLE(IBINS(IHIS))
          TMP = XLIM1-XLIMIT(3,IHIS)
          DO 200 K=1,IBINS(IHIS)+1
            TMP = TMP + XLIMIT(3,IHIS)
            HIST(1,IHIS,K) = TMP
 200      CONTINUE
          ISWI(IHIS) = 1
        ELSE
C  each bin given by user
          DO 300 K=1,IBIN+1
            HIST(1,IHIS,K) = XLIMB(K)
 300      CONTINUE
          ISWI(IHIS) = 2
        ENDIF
        IHISL = IHIS
      ELSE
        WRITE(6,'(/1X,A,I4)') 'NEWHIS:ERROR: TOO MANY HISTOGRAMS',IHIS
        IREFN = -1
        STOP
      ENDIF
      END
C
C
      SUBROUTINE FILHIS(X,Y,IHIS)
C**********************************************************************
C
C     fill value Y at position X into histogramm no IHIS
C
C     output: COMMON HISTOS
C
C**********************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
C
      PARAMETER (NHIS=200, NDIM=201)
      DOUBLE PRECISION HIST
      COMMON /HISTOS/ XLIMIT(3,NHIS),ISWI(NHIS),ENTRY(NHIS),
     &                HIST(5,NHIS,NDIM),HEVENT(NHIS),OVERF(NHIS),
     &                UNDERF(NHIS),IBINS(NHIS),
     &                IHMIN(NHIS),IHMAX(NHIS),IHISL
C
      PARAMETER (ZERO   =  0.,
     &           ONE    =  1.,
     &           TWO    =  2.,
     &           OHALF  =  0.5,
     &           TINY   =  1.e-10)
C
      IF((IHIS.LT.1).OR.(IHIS.GT.IHISL)) THEN
        IDUMMY = 0
        WRITE(6,'(/1X,A,I6)')
     &    'FILHIS:ERROR: histogram out of range ',IHIS
C  simulation of division by zero (to get the call of FILHIS)
        IHIS = IHIS/IDUMMY
        RETURN
      ENDIF

C  find relevant bin position
      IF(X.LT.HIST(1,IHIS,1)) THEN
        I1 = 0
        GOTO 100
      ELSE IF(X.GT.HIST(1,IHIS,IBINS(IHIS)+1)) THEN
        I1 = IBINS(IHIS)+1
        GOTO 100
      ENDIF

      IF(ISWI(IHIS).EQ.1) THEN
C  bin size not given explicitly
        I1 = (X-XLIMIT(1,IHIS))/XLIMIT(3,IHIS)+1
      ELSE IF(ISWI(IHIS).EQ.2) THEN
C  different bin sizes possible
C  binary search algorithm
        KMIN=0
        KMAX=IBINS(IHIS)+1
 300    CONTINUE
          IF((KMAX-KMIN).EQ.1) GOTO 400
          KK=(KMAX+KMIN)/2
          IF(X.LE.HIST(1,IHIS,KK)) THEN
            KMAX=KK
          ELSE
            KMIN=KK
          ENDIF
        GOTO 300
 400    CONTINUE
        I1 = KMIN
      ELSE
        WRITE(6,'(/1X,A,I6)')
     &    'FILHIS:ERROR: histogram not initialized',IHIS
        RETURN
      ENDIF
C
C  histogramming
 100  CONTINUE
      IF(I1.LE.0) THEN
        UNDERF(IHIS) = UNDERF(IHIS) + ONE
      ELSE IF(I1.LE.IBINS(IHIS)) THEN
        HIST(3,IHIS,I1) = HIST(3,IHIS,I1) + ONE
        HIST(4,IHIS,I1) = HIST(4,IHIS,I1) + Y
        ENTRY(IHIS) = ENTRY(IHIS) + ONE
        IF(I1.LT.IHMIN(IHIS)) IHMIN(IHIS) = I1
        IF(I1.GT.IHMAX(IHIS)) IHMAX(IHIS) = I1
      ELSE
        OVERF(IHIS) = OVERF(IHIS) + ONE
      ENDIF
      END
C
C
      SUBROUTINE ADDHIS(IHIS)
C**********************************************************************
C
C     add stored data to histogram no IHIS
C
C     output: COMMON HISTOS
C
C**********************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
C
      PARAMETER (ZERO   =  0.,
     &           ONE    =  1.,
     &           TWO    =  2.,
     &           OHALF  =  0.5,
     &           TINY   =  1.e-10)
C
      PARAMETER (NHIS=200, NDIM=201)
      DOUBLE PRECISION HIST
      COMMON /HISTOS/ XLIMIT(3,NHIS),ISWI(NHIS),ENTRY(NHIS),
     &                HIST(5,NHIS,NDIM),HEVENT(NHIS),OVERF(NHIS),
     &                UNDERF(NHIS),IBINS(NHIS),
     &                IHMIN(NHIS),IHMAX(NHIS),IHISL
C
      IF((IHIS.GT.0).AND.(IHIS.LE.IHISL)) THEN
        DO 300 K=IHMIN(IHIS),IHMAX(IHIS)
          HIST(2,IHIS,K) = HIST(2,IHIS,K) + HIST(4,IHIS,K)
          HIST(5,IHIS,K) = HIST(5,IHIS,K) + HIST(4,IHIS,K)**2
          HIST(4,IHIS,K) = ZERO
 300    CONTINUE
        HEVENT(IHIS) = HEVENT(IHIS)+ONE
        IHMIN(IHIS) = IBINS(IHIS)
        IHMAX(IHIS) = 0
      ELSE IF(IHIS.EQ.-1) THEN
        DO 100 I=1,IHISL
          DO 200 K=IHMIN(I),IHMAX(I)
            HIST(2,I,K) = HIST(2,I,K) + HIST(4,I,K)
            HIST(5,I,K) = HIST(5,I,K) + HIST(4,I,K)**2
            HIST(4,I,K) = ZERO
 200      CONTINUE
          HEVENT(I) = HEVENT(I)+ONE
          IHMIN(I) = IBINS(I)
          IHMAX(I) = 0
 100    CONTINUE
      ELSE
        WRITE(6,'(1X,A,I12)') 'ADDHIS:ERROR: invalid histogram',IHIS
      ENDIF
      END
C
C
      SUBROUTINE CLRHIS(IHIS)
C**********************************************************************
C
C     clear temporarily stored data of histogram no IHIS
C     (no effect after ADDHIS call)
C
C     output: COMMON HISTOS
C
C**********************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
C
      PARAMETER (ZERO   =  0.,
     &           ONE    =  1.,
     &           TWO    =  2.,
     &           OHALF  =  0.5,
     &           TINY   =  1.e-10)
C
      PARAMETER (NHIS=200, NDIM=201)
      DOUBLE PRECISION HIST
      COMMON /HISTOS/ XLIMIT(3,NHIS),ISWI(NHIS),ENTRY(NHIS),
     &                HIST(5,NHIS,NDIM),HEVENT(NHIS),OVERF(NHIS),
     &                UNDERF(NHIS),IBINS(NHIS),
     &                IHMIN(NHIS),IHMAX(NHIS),IHISL
C
      IF((IHIS.GT.0).AND.(IHIS.LE.IHISL)) THEN
        DO 300 K=IHMIN(IHIS),IHMAX(IHIS)
          HIST(4,IHIS,K) = ZERO
 300    CONTINUE
        IHMIN(IHIS) = IBINS(IHIS)
        IHMAX(IHIS) = 0
      ELSE IF(IHIS.EQ.-1) THEN
        DO 100 I=1,IHISL
          DO 200 K=IHMIN(I),IHMAX(I)
            HIST(4,I,K) = ZERO
 200      CONTINUE
          IHMIN(I) = IBINS(I)
          IHMAX(I) = 0
 100    CONTINUE
      ELSE
        WRITE(6,'(1X,A,I12)') 'CLRHIS:ERROR: invalid histogram',IHIS
      ENDIF
      END
C
C
      SUBROUTINE OUTHIS(I,ILOG,HEADER,FAC,INORM)
C**********************************************************************
C
C     write out histogram
C
C     input: COMMON HISTOS
C            I              histogram number
C            ILOG           0 linear axis
C                           1 log y axis
C                           2 log x axis
C                           3 log x and y axis
C            HEADER         header string
C            FAC            factor to change normalization
C            INORM          0 normalization per event and bin width
C                           1 normalization per entry and bin width
C                           2 normalization per bin entry
C                           3 normalization per event and bin width**2
C
C**********************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE

      PARAMETER (NHIS=200, NDIM=201)
      DOUBLE PRECISION HIST
      COMMON /HISTOS/ XLIMIT(3,NHIS),ISWI(NHIS),ENTRY(NHIS),
     &                HIST(5,NHIS,NDIM),HEVENT(NHIS),OVERF(NHIS),
     &                UNDERF(NHIS),IBINS(NHIS),
     &                IHMIN(NHIS),IHMAX(NHIS),IHISL

      PARAMETER (ZERO   =  0.,
     &           IZERO  =  0,
     &           ONE    =  1.,
     &           TWO    =  2.,
     &           OHALF  =  0.5,
     &           EPS    =  1.e-5,
     &           TINY   =  1.e-8)

      DIMENSION HI(NDIM),HE(NDIM)

      PARAMETER(NDIM2 = 2*NDIM)
      DIMENSION XX(NDIM2),YY(NDIM2)
      CHARACTER*72 HEADER

      AVERX  = ZERO
      AVERX2 = ZERO
      AVERY  = ZERO
      AVERY2 = ZERO
      SUM    = ZERO

      IF((I.LT.1).OR.(I.GT.NHIS)) THEN
        WRITE(6,'(/1X,A,I6)') 'OUTHIS:ERROR: histogram out of range ',I
        RETURN
      ELSE IF(ISWI(I).EQ.0) THEN
        WRITE(6,'(/1X,A,I6)')
     &    'OUTHIS:ERROR: histogram not initialized ',I
        RETURN
      ENDIF

      WRITE(6,'(//A,A,/,A,1X,76A)') '# ',HEADER,'#',('=',II=1,76)
      IF(ENTRY(I).LT.0.5D0) THEN
        WRITE(6,'(1X,A,2F10.0)')
     &    'OUTHIS:WARNING: no entries, underflows, overflows',
     &    UNDERF(I),OVERF(I)
        RETURN
      ENDIF

      IF(ILOG.GE.2) THEN
        DO 221 K=1,IBINS(I)+1
          HIST(1,I,K) = EXP(HIST(1,I,K))
 221    CONTINUE
        WRITE(6,'(A,5X,A)') '#','(logarithmic x axis assumed)'
      ENDIF

C  number of events
      FNORM  = MAX(HEVENT(I),ONE)

      DO 222 K=1,IBINS(I)
        HE(K) = HIST(5,I,K)/FNORM - (HIST(2,I,K)/FNORM)**2
        HE(K) = SQRT(MAX(ZERO,HE(K)))/SQRT(FNORM-ONE)
        IF(INORM.EQ.3) THEN
          HI(K) = HIST(2,I,K)/(HIST(1,I,K+1)**2-HIST(1,I,K)**2)*FAC
          HE(K) = HE(K)/(HIST(1,I,K+1)**2-HIST(1,I,K)**2)*FAC
        ELSE
          HI(K) = HIST(2,I,K)/(HIST(1,I,K+1)-HIST(1,I,K))*FAC
          HE(K) = HE(K)/(HIST(1,I,K+1)-HIST(1,I,K))*FAC
        ENDIF
        AVERY = AVERY + HI(K)
        AVERY2 = AVERY2 + HI(K)**2
        SUM = SUM + HI(K)
 222  CONTINUE

      IF(ABS(SUM).LT.EPS) SUM = ONE

      DO 223 K=1,IBINS(I)
        X = (HIST(1,I,K)+HIST(1,I,K+1))/TWO
        AVERX = AVERX   + X * HI(K)/SUM
        AVERX2 = AVERX2 + X**2 * HI(K)/SUM
 223  CONTINUE
      AVERX2 = AVERX2-AVERX**2

      AVERY  = AVERY/ENTRY(I)
      AVERY2 = AVERY2/ENTRY(I)-AVERY**2

      WRITE(6,'(A,5X,3(A,F10.0))') '#','entries ',ENTRY(I),
     &  '      underflows ',UNDERF(I),
     &  '      overflows ',OVERF(I)
      WRITE(6,'(A,5X,A,F10.0)')  
     &  '#','events                   ',HEVENT(I)
      WRITE(6,'(A,5X,A,1P2E12.4)') 
     &  '#','mean value / variance x  ',AVERX,AVERX2
      WRITE(6,'(A,5X,A,1P2E12.4)') 
     &  '#','mean value / variance y  ',AVERY,AVERY2
      WRITE(6,'(A,5X,A,1PE12.4)') 
     &  '#','external scaling factor  ',FAC
      IF(ENTRY(I).LT.0.5D0) RETURN

      IF(INORM.EQ.3)
     &  WRITE(6,'(A,5X,A)') '#','(values divided by squared bin width)'

      WRITE(6,'(A,1X,76A)') '#',('=',II=1,76)
      WRITE(6,'(A,2X,A)')
     &  '#',' X-center, X-min, X-max, Y-value, Y-error (stat.), entries'
      DO 333 K=1,IBINS(I)
        X1 = HIST(1,I,K)
        X2 = HIST(1,I,K+1)

C  suppress output of empty bins
        if((k.ne.1).and.(k.ne.IBINS(I))) then
          if(HIST(3,I,K-1)+HIST(3,I,K)+HIST(3,I,K+1).lt.0.1) goto 200
        endif

        if ((INORM.eq.0).or.(INORM.eq.3)) then
          WRITE(6,'(1X,1P6E12.4)') 
     &      (X1+X2)/2.,X1,X2,HI(K)/FNORM,HE(K),HIST(3,I,K)
        else
          WRITE(6,'(1X,1P6E12.4)') 
     &      (X1+X2)/2.,X1,X2,HI(K)/ENTRY(I),
     &      HI(K)/ENTRY(I)/SQRT(MAX(HIST(3,I,K),1.D0)),HIST(3,I,K)
        endif
 200    continue
        II = 2*K
        XX(II-1) = X1
        XX(II)   = X2
        IF((INORM.EQ.0).OR.(INORM.EQ.3)) THEN
          YY(II-1) = HI(K)/FNORM
          YY(II)   = YY(II-1)
        ELSE IF(INORM.EQ.1) THEN
          YY(II-1) = HI(K)/ENTRY(I)
          YY(II)   = YY(II-1)
        ELSE IF(INORM.EQ.2) THEN
          YY(II-1) = HI(K)*(HIST(1,I,K+1)-HIST(1,I,K))
     &      /FAC/MAX(HIST(3,I,K),1.D0)
          YY(II)   = YY(II-1)
        ELSE
          WRITE(6,'(/,1X,A,I5)') 'OUTHIS:ERROR: unsupported Inorm',
     &      INORM
        ENDIF
 333  CONTINUE
      WRITE(6,'(/1X,A,A69)') 'Preview: ',HEADER
      IBIN2 = 2.*IBINS(I)
      IF(ILOG.EQ.0) THEN
        CALL QGRAPH(0,0,IBIN2,1,XX,YY,YY)
      ELSE IF(ILOG.EQ.1) THEN
        CALL QGRAPH(0,1,IBIN2,1,XX,YY,YY)
      ELSE IF(ILOG.EQ.2) THEN
        CALL QGRAPH(1,0,IBIN2,1,XX,YY,YY)
      ELSE IF(ILOG.EQ.3) THEN
        CALL QGRAPH(1,1,IBIN2,1,XX,YY,YY)
      ELSE
        WRITE(6,'(/1X,A,I3)') 'OUTHIS:ERROR: invalid Ilog',ILOG
      ENDIF
C
      IF(ILOG.GE.2) THEN
        DO 224 K=1,IBINS(I)+1
          HIST(1,I,K) = LOG(HIST(1,I,K))
 224    CONTINUE
      ENDIF
C
      END
C
C
      SUBROUTINE QGRAPH(ILOGX,ILOGY,N,IARG,X,Y1,Y2)
C***********************************************************************
C
C     calculate quasi graphic picture with 25 lines and 79 columns
C     ranges will be chosen automatically
C
C     input     ILOGX      0   linear scale in x
C                          1   logarithmic scale in x
C               ILOGY      0   linear scale in y
C                          1   logarithmic scale in y
C               N          dimension of input fields
C               IARG       number of curves (fields) to plot
C               X          field of X
C               Y1         field of Y1
C               Y2         field of Y2
C
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
      PARAMETER (EPS = 1.D-10, DEPS=1.D-30)
C
      DIMENSION X(N),Y1(N),Y2(N)
      PARAMETER (IYRAST=5,IXRAST=10,IBREIT=79,IZEIL=20)
      CHARACTER SYMB(5)
      CHARACTER COL(0:149,0:49)
C
      DATA SYMB /'0','e','z','#','x'/
C
      ISPALT=IBREIT-10
C
C***  automatic range fitting
C
      XMAX=X(1)
      XMIN=X(1)
      DO 600 I=1,N
        XMAX=MAX(X(I),XMAX)
        XMIN=MIN(X(I),XMIN)
 600  CONTINUE
      IF(ILOGX.EQ.1) THEN
        IF(XMIN.LE.0.) THEN
          WRITE(6,'(/1X,A,2E12.5)') 
     &      'QGRAPH:ERROR: invalid range for log-x plot',XMIN,XMAX
          RETURN
        ENDIF
        XMAX = LOG(XMAX)
        XMIN = LOG(XMIN)
      ENDIF
      XZOOM=(XMAX-XMIN)/DBLE(ISPALT)
C
      ITEST=0
      DO 1100 K=0,IZEIL-1
        ITEST=ITEST+1
        IF(ITEST.EQ.IYRAST) THEN
          DO 1010 L=1,ISPALT-1
            COL(L,K)='-'
 1010     CONTINUE
          COL(ISPALT,K)='+'
          ITEST=0
          DO 1020 L=0,ISPALT-1,IXRAST
            COL(L,K)='+'
 1020     CONTINUE
        ELSE
          DO 1030 L=1,ISPALT-1
            COL(L,K)=' '
 1030     CONTINUE
          DO 1040 L=0,ISPALT-1,IXRAST
            COL(L,K)='|'
 1040     CONTINUE
          COL(ISPALT,K)='|'
        ENDIF
 1100 CONTINUE
C
C***  plot curve Y1
C
      YMAX=Y1(1)
      IF(ILOGY.EQ.0) THEN
        YMIN=Y1(1)
      ELSE
        YMIN=MAX(Y1(1),DEPS)
      ENDIF
      DO 500 I=1,N
        YMAX=MAX(Y1(I),YMAX)
        IF(ILOGY.EQ.0) THEN
          YMIN=MIN(Y1(I),YMIN)
        ELSE
          YMAX =MAX(Y1(I),YMAX)
          IF(Y1(I).GT.DEPS) THEN
            IF(YMIN.EQ.DEPS) THEN
              YMIN = Y1(I)/10.
            ELSE
              YMIN = MIN(Y1(I),YMIN)
            ENDIF
          ENDIF
        ENDIF
 500  CONTINUE
C
      IF(IARG.GT.1) THEN
        DO 550 I=1,N
          YMAX=MAX(Y2(I),YMAX)
          IF(ILOGY.EQ.0) THEN
            YMIN=MIN(Y2(I),YMIN)
          ELSE
            YMAX =MAX(Y2(I),YMAX)
            IF(Y2(I).GT.DEPS) THEN
              IF(YMIN.EQ.DEPS) THEN
                YMIN = Y2(I)/10.
              ELSE
                YMIN = MIN(Y2(I),YMIN)
              ENDIF
            ENDIF
          ENDIF
 550    CONTINUE
      ENDIF
C
      IF(ILOGY.EQ.1) THEN
        DO 560 I=1,N
          Y1(I) = MAX(Y1(I),YMIN)
 560    CONTINUE
        IF(IARG.GT.1) THEN
          DO 570 I=1,N
            Y2(I) = MAX(Y2(I),YMIN)
 570      CONTINUE
        ENDIF
      ENDIF
C
      IF(ILOGY.EQ.0) THEN
        YMA=(YMAX-YMIN)/20.0D0+YMAX
        YMI=YMIN-(YMAX-YMIN)/20.0D0
        YZOOM=(YMA-YMI)/DBLE(IZEIL)
      ELSE
        YMA=LOG(YMAX/YMIN)/20.0D0+LOG(YMAX)
        YMI=LOG(YMIN)-LOG(YMAX/YMIN)/20.0D0
        YZOOM=(YMA-YMI)/DBLE(IZEIL)
      ENDIF
C
      IF(YZOOM.LT.EPS) THEN
        WRITE(6,'(/1X,A,2E12.5)')
     &    'QGRAPH:WARNING: ymin = ymax, output suppressed',YMIN,YMAX
        RETURN
      ENDIF
C
C***  plot curve Y1
C
      ILAST=-1
      LLAST=-1
      DO 1200 K=1,N
        IF(ILOGX.EQ.0) THEN
          L=NINT((X(K)-XMIN)/XZOOM)
        ELSE
          L=NINT((LOG(X(K))-XMIN)/XZOOM)
        ENDIF
        IF(ILOGY.EQ.0) THEN
          I=NINT((YMA-Y1(K))/YZOOM)
        ELSE
          I=NINT((YMA-LOG(Y1(K)))/YZOOM)
        ENDIF
        IF(ILAST.GE.0) THEN
          LD = L-LLAST
          ID = I-ILAST
          DO 55 II=0,LD,SIGN(1,LD)
            DO 66 KK=0,ID,SIGN(1,ID)
              COL(II+LLAST,KK+ILAST)=SYMB(1)
 66         CONTINUE
 55       CONTINUE
        ELSE
          COL(L,I)=SYMB(1)
        ENDIF
        ILAST = I
        LLAST = L
 1200 CONTINUE
C
      IF(IARG.GT.1) THEN
C
C***  plot curve Y2
C
        DO 1250 K=1,N
          IF(ILOGX.EQ.0) THEN
            L=NINT((X(K)-XMIN)/XZOOM)
          ELSE
            L=NINT((LOG(X(K))-XMIN)/XZOOM)
          ENDIF
          IF(ILOGY.EQ.0) THEN
            I=NINT((YMA-Y2(K))/YZOOM)
          ELSE
            I=NINT((YMA-LOG(Y2(K)))/YZOOM)
          ENDIF
          COL(L,I)=SYMB(2)
 1250   CONTINUE
      ENDIF
C
C***  write it
C
      IF(ILOGX.EQ.1) WRITE(6,'(2X,A)') '(logarithmic x axis assumed)'
      IF(ILOGY.EQ.1) WRITE(6,'(2X,A)') '(logarithmic y axis)'
      WRITE(6,'(1X,79A)') ('-',I=1,IBREIT)
C
C***  write range of X
C
      XZOOM = (XMAX-XMIN)/DBLE(7)
      IF(ILOGX.EQ.0) THEN
        WRITE(6,120) (XZOOM*DBLE(I-1)+XMIN,I=1,7)
      ELSE
        WRITE(6,120) (EXP(XZOOM*DBLE(I-1)+XMIN),I=1,7)
      ENDIF
C
      DO 1300 K=0,IZEIL-1
        IF(ILOGY.EQ.0) THEN
          YPOS=YMA-((DBLE(K)+0.5D0)*YZOOM)
        ELSE
          YPOS=EXP(YMA-((DBLE(K)+0.5D0)*YZOOM))
        ENDIF
        WRITE(6,110) YPOS,(COL(I,K),I=0,ISPALT)
 110    FORMAT(1X,1PE9.2,70A1)
 1300 CONTINUE
C
C***  write range of X
C
      XZOOM = (XMAX-XMIN)/DBLE(7)
      IF(ILOGX.EQ.0) THEN
        WRITE(6,120) (XZOOM*DBLE(I-1)+XMIN,I=1,7)
      ELSE
        WRITE(6,120) (EXP(XZOOM*DBLE(I-1)+XMIN),I=1,7)
      ENDIF
      WRITE(6,'(1X,79A)') ('-',I=1,IBREIT)
 120  FORMAT(6X,7(1PE10.3))
      END
C
C
      SUBROUTINE MOMENT(IHIS)
C**********************************************************************
C
C     calculated normalized moments of distribution binned in
C                                      histogram IHIS
C
C**********************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
C
      PARAMETER (NHIS=200, NDIM=201)
      DOUBLE PRECISION HIST
      COMMON /HISTOS/ XLIMIT(3,NHIS),ISWI(NHIS),ENTRY(NHIS),
     &                HIST(5,NHIS,NDIM),HEVENT(NHIS),OVERF(NHIS),
     &                UNDERF(NHIS),IBINS(NHIS),
     &                IHMIN(NHIS),IHMAX(NHIS),IHISL
C
      PARAMETER (ZERO   =  0.,
     &           IZERO  =  0,
     &           ONE    =  1.,
     &           TWO    =  2.,
     &           OHALF  =  0.5,
     &           TINY   =  1.e-10)
      DIMENSION YM(5),C(5)
C
      DO 50 I=1,5
        YM(I) = ZERO
 50   CONTINUE
C
      SUM = ZERO
      DO 100 K=1,NDIM
        SUM = SUM + HIST(2,IHIS,K)
 100  CONTINUE
C
C  absolute moments
      DO 200 K=1,NDIM
        X = (DBLE(K-1)+OHALF)*XLIMIT(3,IHIS)+XLIMIT(1,IHIS)
        XX = X
        YM(1) = YM(1) + HIST(2,IHIS,K)*XX
        XX = XX*X
        YM(2) = YM(2) + HIST(2,IHIS,K)*XX
        XX = XX*X
        YM(3) = YM(3) + HIST(2,IHIS,K)*XX
        XX = XX*X
        YM(4) = YM(4) + HIST(2,IHIS,K)*XX
        XX = XX*X
        YM(5) = YM(5) + HIST(2,IHIS,K)*XX
 200  CONTINUE
C
C  normalized moments
      DO 300 K=1,5
        YM(K)=YM(K)/SUM
        C(K) = YM(K)/YM(1)**K
 300  CONTINUE
C
C  write results
      WRITE(6,'(/1X,A,/1X,A)') 'Absolute and scaled moments',
     &                         '==========================='
      DO 400 K=1,5
        WRITE(6,'(5X,A,I2,3X,A,E12.3,A,E12.3)')
     &    ' Rank ',K,' absolute ',YM(K),'          normalized ',C(K)
 400  CONTINUE
C
      END
