      SUBROUTINE HSPLI (KF, KP1,KP2)
C...This subroutine splits one hadron of code KF
C.  into 2 partons of code KP1 and KP2
C.  KP1 refers to a color triplet [q or (qq)bar]         
C.  KP2 to a a color anti-triplet [qbar or (qq)]         
C.  allowed inputs:
C.  KF = 6:14 pi0,pi+-,k+-,k0L,k0s, p,n
C.     = -13,-14  pbar,nbar
C.     = 34 Sigma+
C.     = 35/39 Sigma0/Lambda0
C.     = 36    Sigma-
C.     = 37    Xi0
C.     = 38    Xi-
C-------------------------------------------------
      SAVE

      Select case(IABS(KF))

      Case (6)                  ! pi0
        R = S_RNDM(0)              
        IF (R.LE.0.)  THEN
           KP1 = 1                  
           KP2 = -1
        ELSE
          KP1 = 2
          KP2 = -2
        ENDIF
      Case (7)                  ! pi+
        KP1 = 1                  
        KP2 = -2
      Case (8)                  ! pi-
        KP1 = 2                  
        KP2 = -1
      Case (9)                  ! K+
        KP1 = 1                  
        KP2 = -3
      Case (10)                  ! K-
        KP1 = 3                  
        KP2 = -1
      Case (11,12)                  ! K0S/K0L
        KP1 = 2
        KP2 = -3
        IF (S_RNDM(0).GT. 0.5)  THEN
          KP1 = 3
          KP2 = -2
        ENDIF
      Case (13)                  ! p/pbar
        R = 6.*S_RNDM(0)            
        IF (R .LT.3.)       THEN
          KP1 = 1
          KP2 = 12
        ELSEIF (R .LT. 4.)  THEN
          KP1 = 1
          KP2 = 21
        ELSE
          KP1 = 2
          KP2 = 11
        ENDIF
      Case (14)                 ! n/nbar
        R = 6.*S_RNDM(0)                  
        IF (R .LT.3.)       THEN
           KP1 = 2
           KP2 = 12
        ELSEIF (R .LT. 4.)  THEN
          KP1 = 2
          KP2 = 21
        ELSE
          KP1 = 1
          KP2 = 22
        ENDIF
      Case (34)                 !Sigma+
        R = 6.*S_RNDM(0)                  
        IF (R .LT.3.)       THEN
           KP1 = 3
           KP2 = 11
        ELSEIF (R .LT. 4.)  THEN
          KP1 = 1
          KP2 = 31
        ELSE
          KP1 = 1
          KP2 = 13
        ENDIF
      Case (35,39)              !Sigma0/Lambda0
        R = 6.*S_RNDM(0)                  
        IF (R .LT.3.)       THEN
           KP1 = 3
           KP2 = 21
        ELSEIF (R .LT. 4.)  THEN
          KP1 = 1
          KP2 = 32
        ELSE
          KP1 = 2
          KP2 = 13
        ENDIF
      Case (36)                 !Sigma-
        R = 6.*S_RNDM(0)                  
        IF (R .LT.3.)       THEN
           KP1 = 3
           KP2 = 22
        ELSEIF (R .LT. 4.)  THEN
          KP1 = 2
          KP2 = 32
        ELSE
          KP1 = 2
          KP2 = 23
        ENDIF
      Case (37)                 !Xi0
        R = 6.*S_RNDM(0)                  
        IF (R .LT.3.)       THEN
           KP1 = 1
           KP2 = 33
        ELSEIF (R .LT. 4.)  THEN
          KP1 = 3
          KP2 = 13
        ELSE
          KP1 = 1
          KP2 = 33
        ENDIF
      Case (38)                 !Xi-
        R = 6.*S_RNDM(0)                  
        IF (R .LT.3.)       THEN
           KP1 = 2
           KP2 = 33
        ELSEIF (R .LT. 4.)  THEN
          KP1 = 3
          KP2 = 23
        ELSE
          KP1 = 2
          KP2 = 33
        ENDIF
      Case Default
C...Test for good input
        WRITE(6,*)
     &      'HSPLI : Routine entered with illegal particle code ',KF
        STOP ! This has to be replaced by the usual SIBYLL exception if exists.
      End Select

C if anti-baryon, invert valences
      IF (KF .LT. 0) THEN
        KPP = KP1
        KP1 = -KP2
        KP2 = -KPP
      ENDIF
      RETURN
      END