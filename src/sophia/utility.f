      SUBROUTINE TOEVT()
C************************************************************************
C
C     Converts data in S_PLIST common block to
C     standard HEPEVT common block format, 
C     and S_PLIST1 and S_CHP to SCHG common block containing
C     charge and parent decayed particle information
C     
C************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

      SAVE
      INTEGER PDGCOD(49)
      DATA PDGCOD /
     &  22, -11, 11, -13, 13, 111, 211, -211, 321, -321, 130, 310, 2212,
     &  2112, 12, -12, 14, -14, -99999999, -99999999, 311, -311, 221, 
     &  331, 213, -213, 113, 323, -323, 313, -313, 223, 333, 3222, 3212,
     &  3112, 3322, 3312, 3122, 2224, 2214, 2114, 1114, 3224, 3214, 
     &  3114, 3324, 3314, 3334 /
       
      COMMON /S_PLIST/ P(2000,5), LLIST(2000), NP, Ideb
      COMMON /S_PLIST1/ LLIST1(2000)
      COMMON /S_CHP/ S_LIFE(49), ICHP(49), ISTR(49), IBAR(49)
        
C     Standard HEPEVT common block:
      PARAMETER (NMXHEP = 2000)
      PARAMETER (NMXHP2 = NMXHEP * 2, NMXHP4 = NMXHEP * 4, 
     & NMXHP5 = NMXHEP * 5)
      COMMON /HEPEVT/ NEVHEP, NHEP, ISTHEP(NMXHEP), IDHEP(NMXHEP),
     & JMOHEP(2, NMXHEP), JDAHEP(2, NMXHEP), 
     & PHEP(5, NMXHEP), VHEP(4, NMXHEP)

C     Initialization of HEVEVT variables
      DATA NEVHEP / 0 /
      DATA JMOHEP/ NMXHP2 * 0.0/, JDAHEP / NMXHP2 * 0.0 /
      DATA PHEP/ NMXHP5 * 0.0/, VHEP / NMXHP4 * 0.0 /

C     Common block SCHG with additional information:
C     ICHG - charge, 
C     DECPAR - position in the arrays (ISTHEP, PHEP, ...) of parent (decayed) particle     
      INTEGER ICHG, DECPAR
      COMMON /SCHG/ ICHG(NMXHEP), DECPAR(NMXHEP)

      INTEGER CID  ! id of the current (in the loop) particle
C     number added to decayed particle:
      INTEGER DNUM
      PARAMETER (DNUM = 10000)
       
      NEVHEP = NEVHEP + 1
      NHEP = NP
C     We do not transpose P, because we don't use PHEP         
C      PHEP = TRANSPOSE(P)
        
      DO I = 1, NHEP
        CID = LLIST(I)
        IF (ABS(CID) .GT. DNUM) THEN
            ISTHEP(I) = 2
            CID = CID - ISIGN(DNUM, CID)
        ELSE 
            ISTHEP(I) = 1
        END IF     
        IDHEP(I) = ISIGN(1, CID) * PDGCOD(ABS(CID))
        ICHG(I) = ICHP(ABS(CID))
        DECPAR(I) = LLIST1(I)
      ENDDO
  
      END      
      
          
      
C This subroutine should be removed:      
      SUBROUTINE PREPARE_EVENT_DATA()
C************************************************************************
C
C     Converts data in S_PLIST common block to convenient format.
C     The converted data are saved in EVENT_DATA common block.
C
C    
C     * Convert particles ids used in SOPHIA to PDG codes: LLIST -> PDG_IDS
C     * Transpose P(2000,5) to PARTICLE_MOMENTA(5, 2000)
C     * STATUS_CODES are set to 1 (final particles) to mimic
C     HEPEVT common block format
C     * At every call EVENT_NUMBER is incremented to count number of
C     events
C     
C     
C************************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

      SAVE
      INTEGER PDG_CODES(49)
      DATA PDG_CODES /
     &  22, -11, 11, -13, 13, 111, 211, -211, 321, -321, 130, 310, 2212,
     &  2112, 12, -12, 14, -14, -99999999, -99999999, 311, -311, 221, 
     &  331, 213, -213, 113, 323, -323, 313, -313, 223, 333, 3222, 3212,
     &  3112, 3322, 3312, 3122, 2224, 2214, 2114, 1114, 3224, 3214, 
     &  3114, 3324, 3314, 3334 /
     
      COMMON /S_PLIST/ P(2000,5), LLIST(2000), NP, Ideb
      COMMON /S_PLIST1/ LLIST1(2000)
      COMMON /S_CHP/ S_LIFE(49), ICHP(49), ISTR(49), IBAR(49)
      
      INTEGER CURRENT_ID
      INTEGER EVENT_NUMBER, NUMBER_OF_PARTICLES
      INTEGER STATUS_CODES(2000), PDG_IDS(2000)
      INTEGER PARENTS(2000), CHARGES(2000)
      DOUBLE PRECISION VERTICES(4, 2000)
      DOUBLE PRECISION JMOHEP(2, 2000), JDAHEP(2, 2000) 
      DATA STATUS_CODES / 2000 * 1 /, EVENT_NUMBER / 0 /
      DATA VERTICES / 8000 * 0 /
      COMMON /EVENT_DATA/ EVENT_NUMBER, NUMBER_OF_PARTICLES, 
     & STATUS_CODES, PDG_IDS, CHARGES, PARENTS,
     & VERTICES, JMOHEP, JDAHEP
     
      EVENT_NUMBER = EVENT_NUMBER + 1
      NUMBER_OF_PARTICLES = NP   


      
      DO I = 1, NP
        CURRENT_ID = LLIST(I)
        IF (ABS(CURRENT_ID) .GT. 10000) THEN
            STATUS_CODES(I) = 2
            CURRENT_ID = CURRENT_ID - ISIGN(10000, CURRENT_ID)
        ELSE 
            STATUS_CODES(I) = 1
        END IF
        CHARGES(I) = ICHP(ABS(CURRENT_ID))     
        PDG_IDS(I) = ISIGN(1, CURRENT_ID) * PDG_CODES(ABS(CURRENT_ID))
        PARENTS(I) = LLIST1(I)
      ENDDO

      END
