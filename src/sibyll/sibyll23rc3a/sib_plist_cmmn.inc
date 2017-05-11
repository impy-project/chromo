C     The final particle output is contained in COMMON /S_PLIST/    
C     NP           : number of final particles
C     P(1:NP, 1:5) : 4-momenta + masses of the final particles 
C     LLIST (1:NP) : codes of final particles
      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP
