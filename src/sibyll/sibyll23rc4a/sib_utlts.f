      SUBROUTINE SIB_ALTRA(GA,BGX,BGY,BGZ,PCX,PCY,PCZ,EC,P,PX,PY,PZ,E)
C*********************************************************************
C
C    arbitrary Lorentz transformation
C
C     Input: GA : gamma factor
C            BG? : components of gamma * beta
C            PC?,EC : components of initial 4 vector
C
C     Output: P?,E : components of 4vector in final frame
C             P : 3-norm in final frame, a.k.a momentum
C
C     PHO_ALTRA taken from PHOJET /FR'14
C*********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      INCLUDE 'sib_debug_cmmn.inc'
c     consistency check: (gamma*beta)**2 = gamma**2 - 1
      BETGAM2 = BGX**2+BGY**2+BGZ**2
      xtst = 1.0-BETGAM2/GA**2-1.D0/GA**2
      IF(abs(xtst).gt.1.D-5) THEN
         WRITE(LUN,*) 'SIB_ALTRA: transf. inconsistent!'
         WRITE(LUN,*) 'SIB_ALTRA: input (GA,GABE):',GA,BGX,BGY,BGZ
      ENDIF
      IF(GA.LT.1.D0) THEN
         WRITE(LUN,*) 'SIB_ALTRA: you are joking right? GAMMA=',GA
         call sib_reject
      ENDIF
      EP=PCX*BGX+PCY*BGY+PCZ*BGZ
      PE=EP/(GA+1.D0)+EC
      PX=PCX+BGX*PE
      PY=PCY+BGY*PE
      PZ=PCZ+BGZ*PE
      P=SQRT(PX*PX+PY*PY+PZ*PZ)
      E=GA*EC+EP
      END


      SUBROUTINE SIB_TRANS(XO,YO,ZO,CDE,SDE,CFE,SFE,X,Y,Z)
C**********************************************************************
C
C  rotation of coordinate frame (1) de rotation around y axis
C                               (2) fe rotation around z axis
C  (inverse rotation to SIB_TRANI)
C
C     Input: ?0 : vector components in initial frame
C            C? : cosine of rotation angle
C            S? : sine of rotation angle
C            DE : angle of rotation around y axis 
C                 (polar angle in spherical coord.)
C            FE : angle of rotation around z axis 
C                 (azimuthal angle in spherical coord.)
C
C     Output: X,Y,Z: components of vector in rotated frame
C
C     PHO_TRANS taken from PHOJET \FR'14
C**********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      X= CDE*CFE*XO-SFE*YO+SDE*CFE*ZO
      Y= CDE*SFE*XO+CFE*YO+SDE*SFE*ZO
      Z=-SDE    *XO       +CDE    *ZO
      END


      SUBROUTINE SIB_TRANI(XO,YO,ZO,CDE,SDE,CFE,SFE,X,Y,Z)
C**********************************************************************
C
C  rotation of coordinate frame (1) -fe rotation around z axis
C                               (2) -de rotation around y axis
C  (inverse rotation to SIB_TRANS)
C
C     Input: ?0 : vector components in initial frame
C            C? : cosine of rotation angle
C            S? : sine of rotation angle
C            DE : angle of rotation around y axis 
C                 (polar angle in spherical coord.)
C            FE : angle of rotation around z axis 
C                 (azimuthal angle in spherical coord.)
C
C     Output: X,Y,Z: components of vector in rotated frame
C
C     PHO_TRANS taken from PHOJET \FR'14
C**********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      X= CDE*CFE*XO+CDE*SFE*YO-SDE*ZO
      Y=-SFE    *XO+CFE*    YO
      Z= SDE*CFE*XO+SDE*SFE*YO+CDE*ZO
      END


      SUBROUTINE SIROBO( NBEG, NEND, THE, PHI, DBEX, DBEY, DBEZ)
C **********************************************************************
C   THIS IS A SLIGHTLY ALTERED VERSION OF "LUROBO" [JETSET63.PYTHIA]   *
C SET TO WORK IN THE SIBYL ENVIROMENT. THE TRANSFORMATION IS PERFORMED *
C ON PARTICLES NUMBER FROM NBEG TO NEND. COMMON BLOCKS CHANGED.        *
C                                      TSS,   Oct '87                  *
C  modification  use directly BETA in double precision in input (PL)   *
C **********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      COMMON /S_PLIST/ P(8000,5), LLIST(8000), NP
      DIMENSION ROT(3,3),PV(3),DP(4)
      IF(THE**2+PHI**2 .LE. 1D-20) GO TO 131
C...ROTATE (TYPICALLY FROM Z AXIS TO DIRECTION THETA,PHI)
       ROT(1,1)=dCOS(THE)*dCOS(PHI)
       ROT(1,2)=-dSIN(PHI)
       ROT(1,3)=dSIN(THE)*dCOS(PHI)
       ROT(2,1)=dCOS(THE)*dSIN(PHI)
       ROT(2,2)=dCOS(PHI)
       ROT(2,3)=dSIN(THE)*dSIN(PHI)
       ROT(3,1)=-dSIN(THE)
       ROT(3,2)=0.D0
       ROT(3,3)=dCOS(THE)
       DO 120 I=NBEG,NEND
       DO 100 J=1,3
 100   PV(J)=P(I,J)
       DO 110 J=1,3
 110   P(I,J)=ROT(J,1)*PV(1)+ROT(J,2)*PV(2)+ROT(J,3)*PV(3)
 120   CONTINUE
 131    IF(DBEX**2+DBEY**2+DBEZ**2 .LE. 1D-20) GO TO 151
C...LORENTZ BOOST (TYPICALLY FROM REST TO MOMENTUM/ENERGY=BETA)
       DGA=1D0/DSQRT(1D0-DBEX**2-DBEY**2-DBEZ**2)
       DO 140 I=NBEG, NEND
       DO 130 J=1,4
 130   DP(J)=P(I,J)
       DBEP=DBEX*DP(1)+DBEY*DP(2)+DBEZ*DP(3)
       DGABEP=DGA*(DGA*DBEP/(1D0+DGA)+DP(4))
       P(I,1)=DP(1)+DGABEP*DBEX
       P(I,2)=DP(2)+DGABEP*DBEY
       P(I,3)=DP(3)+DGABEP*DBEZ
       P(I,4)=DGA*(DP(4)+DBEP)
 140   CONTINUE
 151   RETURN
      END



      subroutine iswtch_lmnts(ia,ib)
      itmp = ia
      ia = ib
      ib = itmp
      end

      subroutine swtch_lmnts(a,b)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      tmp = a
      a = b
      b = tmp
      end

      DOUBLE PRECISION FUNCTION PAWT(A,B,C)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA TWO,EPS10 /2.D0,1.D-10/
C...  c.m.s. Momentum in two particle decays
      PAWT = SQRT((A**2-(B+C)**2+EPS10)*(A**2-(B-C)**2))/(TWO*A)
      END


      SUBROUTINE HSPLI (KF,KP1,KP2)
C...This subroutine splits one hadron of code KF
C.  into 2 partons of code KP1 and KP2
C.  KP1 refers to a color triplet [q or (qq)bar]         
C.  KP2 to a a color anti-triplet [qbar or (qq)]         
C.  allowed inputs:
C.  KF = 6:14 pi0,pi+-,k+-,k0L,k0s, p,n
C.     = -13,-14  pbar,nbar
C.     = 34:39 Sig+, Sig0, Sig-, Xi0, Xi-, Lam0 
C.     = 49: Omega-
C.   \AF'13
C------------------------------------------------
Cf2py integer, intent(in) :: KF
Cf2py integer, intent(out) :: KP1
Cf2py integer, intent(out) :: KP2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'

      DIMENSION KPL(3)
      DATA KPL /1,2,3/

      Select case(IABS(KF))

      Case (6)                  ! pi0
         R = S_RNDM(0)              
         XBUG = 0.D0
         IF(IPAR(19).eq.1) XBUG = 0.5D0
         IF (R.LE.XBUG)  THEN
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
      Case (10)                 ! K-
         KP1 = 3                  
         KP2 = -1
      Case (11,12)              ! K0S/K0L
         KP1 = 2
         KP2 = -3
         IF (S_RNDM(0).GT. 0.5D0)  THEN
            KP1 = 3
            KP2 = -2
         ENDIF
      Case (21)                 ! K0
         KP1 = 2
         KP2 = -3
      Case (22)                 ! K0bar
         KP1 = 3
         KP2 = -2
      Case(33)                  ! phi
         KP1 = 3
         KP2 = -3
      Case (13,41)              ! p/pbar,delta+
         R = PAR(53)*S_RNDM(KF)
         IF (R .LT.3.D0)       THEN
            KP1 = 1
            KP2 = 12
         ELSEIF (R .LT. 4.D0)  THEN
            KP1 = 1
            KP2 = 21
         ELSE
            KP1 = 2
            KP2 = 11
         ENDIF
      Case (14,42)              ! n/nbar,delta0
         R = 6.D0*S_RNDM(0)                  
         IF (R .LT.3.D0)       THEN
            KP1 = 2
            KP2 = 12
         ELSEIF (R .LT. 4.D0)  THEN
            KP1 = 2
            KP2 = 21
         ELSE
            KP1 = 1
            KP2 = 22
         ENDIF
      Case (40)                 ! delta++
         KP1 = 1
         KP2 = 11
      Case (43)                 ! delta-
         KP1 = 2
         KP2 = 22
      Case (34)                 !Sigma+
         R = 6.D0*S_RNDM(0)                  
         IF (R .LT.3.D0)       THEN
            KP1 = 3
            KP2 = 11
         ELSEIF (R .LT. 4.D0)  THEN
            KP1 = 1
            KP2 = 31
         ELSE
            KP1 = 1
            KP2 = 13
         ENDIF
      Case (35,39)              !Sigma0/Lambda0     
c     all configurations equally likely --> Knuth shuffle
         do ii=3,2,-1
            jj=1+INT(ii*S_RNDM(0))
            IFL=KPL(jj)
            KPL(jj)=KPL(ii)
            KPL(ii)=IFL
         enddo
         KP1=KPL(1)
         KP2=KPL(2)*10+KPL(3)
      Case (36)                 !Sigma-
         R = 6.D0*S_RNDM(0)                  
         IF (R .LT.3.D0)       THEN
            KP1 = 3
            KP2 = 22
         ELSEIF (R .LT. 4.D0)  THEN
            KP1 = 2
            KP2 = 32
         ELSE
            KP1 = 2
            KP2 = 23
         ENDIF
      Case (37)                 !Xi0
         R = 6.D0*S_RNDM(0)                  
         IF (R .LT.3.D0)       THEN
            KP1 = 1
            KP2 = 33
         ELSEIF (R .LT. 4.D0)  THEN
            KP1 = 3
            KP2 = 13
         ELSE
            KP1 = 1
            KP2 = 33
         ENDIF
      Case (38)                 !Xi-
         R = 6.D0*S_RNDM(0)                  
         IF (R .LT.3.D0)       THEN
            KP1 = 2
            KP2 = 33
         ELSEIF (R .LT. 4.D0)  THEN
            KP1 = 3
            KP2 = 23
         ELSE
            KP1 = 2
            KP2 = 33
         ENDIF
      Case(49)                  ! Omega-
         KP1 = 3
         KP2 = 33

      Case(59)                  ! D+
         KP1 = 4
         KP2 = -2

      Case(60)                  ! D-
         KP1 = 2
         KP2 = -4

      Case(71)                  ! D0
         KP1 = 4
         KP2 = -1

      Case(72)                  ! D0bar
         KP1 = 1
         KP2 = -4

      Case(73)                  ! eta_c
         KP1 = 4
         KP2 = -4

      Case(74)                  ! Ds+
         KP1 = 4
         KP2 = -3

      Case(75)                  ! Ds-
         KP1 = 3
         KP2 = -4

      Case(76)                  ! Ds*+
         KP1 = 4
         KP2 = -3

      Case(77)                  ! Ds*-
         KP1 = 3
         KP2 = -4

      Case(78)                  ! D*+
         KP1 = 4
         KP2 = -2

      Case(79)                  ! D*-
         KP1 = 2
         KP2 = -4

      Case(80)                  ! D*0
         KP1 = 4
         KP2 = -1

      Case(81)                  ! D*0bar
         KP1 = 1
         KP2 = -4

      Case(83)                  ! J/psi
         KP1 = 4
         KP2 = -4

      Case Default
C...  Test for good input
         WRITE(LUN,*)
     &        'HSPLI : Routine entered with illegal particle code ',KF
         a = -1.D0
         a = dlog(a)
         STOP ! This has to be replaced by the usual SIBYLL exception if exists.
      End Select

C     if anti-baryon, invert valences
      IF (KF .LT. 0) THEN
         KPP = KP1
         KP1 = -KP2
         KP2 = -KPP
      ENDIF
      RETURN
      END

c$$$      subroutine invert_array (yy, xmin, dx, n, xnew, ymin, dy)
c$$$C..    This subroutine receives one   array
c$$$C      of n y values in input yy(1:n)
c$$$C      that correspond to  equispaced values of x_j = xmin + dx*(j-1)
c$$$C
c$$$C      and "reverse" the array returning an array of  x values
c$$$C      xnew (1:n) that  corresponds to equispaced values of y
c$$$C      The relation is assumed monotonous but can be 
c$$$C      increasing or decreasing
c$$$C..............................................................
c$$$      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c$$$      SAVE
c$$$      dimension  yy(n), xnew (n)
c$$$      ymin = yy(1)
c$$$      ymax = yy(n)
c$$$      dy = (ymax - ymin)/float(n-1)
c$$$      xnew (1) = xmin
c$$$      xnew (n) = xmin + dx*float(n-1)
c$$$      k0 = 1
c$$$      do j=2,n-1
c$$$         y = ymin + float(j-1)*dy 
c$$$         do k=k0,n
c$$$            if((yy(k) .gt. y) .eqv. (yy(n) .gt. yy(1))) goto 100
c$$$         enddo
c$$$100      y2 = yy(k)
c$$$         y1 = yy(k-1)
c$$$         k0 = k-1
c$$$         x1 = xmin + dx*float(k-2)
c$$$         x2 = x1+dx
c$$$         xnew (j)  = x1 + dx* (y-y1)/(y2-y1)
c$$$      enddo
c$$$      return
c$$$      end
c$$$
c$$$
c$$$      FUNCTION GASDEV(Idum)
c$$$C...Gaussian deviation
c$$$      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c$$$      SAVE
c$$$      SAVE GSET
c$$$      DATA ISET/0/
c$$$      IF (ISET.EQ.0) THEN
c$$$1       V1=2.*S_RNDM(0)-1.
c$$$        V2=2.*S_RNDM(0)-1.
c$$$        R=V1**2+V2**2
c$$$        IF(R.GE.1.)GO TO 1
c$$$        FAC=SQRT(-2.*LOG(R)/R)
c$$$        GSET=V1*FAC
c$$$        GASDEV=V2*FAC
c$$$        ISET=1
c$$$      ELSE
c$$$        GASDEV=GSET
c$$$        ISET=0
c$$$      ENDIF
c$$$      RETURN
c$$$      END
