C--------------------------------------------------------------------
C     SIBYLL common blocks containing event information       \FR'14
C--------------------------------------------------------------------

C     EVENT INFO COMMON
C     contains overall interaction properties, like
C     SQS : center-of-mass energy
C     S   :         "       "     squared
C     PTmin : low pt cut of QCD cross section, 
C             i.e. minimal pt of hard minijets
C     Xmin : low-x bound for PDFs, 
C            i.e. minimal momentum fraction of hard partons
C     Zmin : logarithm of that
C     KB : PID of beam hadron
C     KT() : PID of target
C     IAT : mass number of target
      PARAMETER (NWD_max = 20)
      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
      INTEGER KB,IAT
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NWD_max), IAT
