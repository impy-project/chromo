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
C     IPDG is the conversion table from sophia ID (index of the array) 
C     to pdg ID (value at the index)     
      INTEGER IPDG(49)
      DATA IPDG /
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
C     ICHG - charge
      INTEGER ICHG
      COMMON /SCHG/ ICHG(NMXHEP)

C     IDS is sophia ID of the current particle
      INTEGER IDS
C     number added to decayed particle:
      INTEGER IDEC
      PARAMETER (IDEC = 10000)

C     If KEEPDC = .TRUE. then LLIST1 has correct entries
C     otherwise LLIST1 contains not relevant data
      LOGICAL KEEPDC
      COMMON /EG_IO/ KEEPDC
       
      NEVHEP = NEVHEP + 1
      NHEP = NP

      DO I = 1, NHEP
C       Convert P to PHEP
        PHEP(1,I) = P(I,1)
        PHEP(2,I) = P(I,2)
        PHEP(3,I) = P(I,3)
        PHEP(4,I) = P(I,4)
        PHEP(5,I) = P(I,5)    

C       IDS - sophia ID        
        IDS = LLIST(I)
C       Normalize ID of decayed particle to sophia ID
C       by subtracting IDEC = 10000
C       and record the status ISTHEP:
C       ISTHEP = 1 for final particle 
C       ISTHEP = 2 for decayed particle        
        IF (ABS(IDS) .GT. IDEC) THEN
            ISTHEP(I) = 2
            IDS = IDS - ISIGN(IDEC, IDS)
        ELSE 
            ISTHEP(I) = 1
        END IF

C       Convert sophia ID to pdg ID of the particle           
        IDHEP(I) = ISIGN(1, IDS) * IPDG(ABS(IDS))
C       Record charge of the particle        
        ICHG(I) = ICHP(ABS(IDS))
        IF (KEEPDC) THEN
            JMOHEP(1, I) = LLIST1(I)
        END IF
      ENDDO

      END
      