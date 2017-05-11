
C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION ZERO,HALF,ONE,ONEHALF,TWO,THREE,FOUR,EIGHT
      COMMON /SIB_NUM/ ZERO,HALF,ONE,ONEHALF,TWO,THREE,FOUR,EIGHT

      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
