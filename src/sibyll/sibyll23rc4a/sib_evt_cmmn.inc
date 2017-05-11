C--------------------------------------------------------------------
C     SIBYLL common blocks containing event variables       \FR'14
C--------------------------------------------------------------------

C     parameters that represent: NW: max. number of wounded nucleons,
C     NS,NH: max. number of soft and hard interactions
      PARAMETER (NW_max = 20)
      PARAMETER (NS_max = 20, NH_max = 100)

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


      INTEGER IBMRDX,ITGRDX,IHMJDX,ISMJDX,ICSTDX
      COMMON /S_INDX/ IBMRDX(3),ITGRDX(NW_max,3),IHMJDX(NW_max*NH_max),
     &     ISMJDX(NW_max*NS_max),ICSTDX(2*NW_max,3)

      PARAMETER (NP_max=8000)
C     The final particle output is contained in COMMON /S_PLIST/    
C     NP           : number of final particles
C     P(1:NP, 1:5) : 4-momenta + masses of the final particles 
C     LLIST (1:NP) : codes of final particles
      DOUBLE PRECISION P
      INTEGER NP,LLIST
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP

C     array of pointer to mother particles of unstable particles
      INTEGER LLIST1
      COMMON /S_PLIST1/ LLIST1(NP_max)
      
      INTEGER NFORIG,NPORIG,IPFLAG,NINT
      COMMON /S_PARTO/ NFORIG(NP_max), NPORIG(NP_max), IPFLAG, NINT
