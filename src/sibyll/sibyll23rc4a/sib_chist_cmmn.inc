C     parameters that represent: NW: max. number of wounded nucleons,
C     NS,NH: max. number of soft and hard interactions
c      PARAMETER (NW_max = 20)
C     The COMMON block /S_CHIST/ contains information about the
C     the structure of the  generated event:
C     NWD   = number of wounded nucleons
C     NJET = total number of hard interactions
C     NSOF = total number of soft interactions
C     NNSOF (1:NW) = number of soft pomeron cuts in each interaction
C     NNJET (1:NW) = number of minijets produced in each interaction 
C     JDIF(1:NW) = diffraction code 
C                  0 : non-diff,
C                  1 : beam-diff
C                  2 : target-diff
C                  3 : double-diff
      INTEGER NNSOF,NNJET,JDIF,NWD,NJET,NSOF
      COMMON /S_CHIST/ NNSOF(NW_max),NNJET(NW_max),
     &     JDIF(NW_max),NWD,NJET,NSOF
