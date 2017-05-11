      SUBROUTINE SIB_I4FLAV (IFL1, IFL2_A, IRNK, IFL2, KF)
C-----------------------------------------------------------------------
C.  This subroutine receives as input IFL1 the flavor code
C.  of a quark (antiquark) and  generates the antiquark (quark)
C.  of flavor code IFL2 that combine with the original parton
C.  to compose an hadron of code KF.
C.
C.  updated to 4 FLAVORS \FR'13
C.  Baryon sector is from jetset code
C.  assuming D*_s+- are J=1, only Charm=1 baryons
C.
C.  If (IFL2_A.NE.0) returns an hadron KF composed of IFL1 and IFL2_A
c
c     Input: IFL1 - flavor of first quark
c            IFL2_A - flavor of second quark ( if 0 randomly chosen ) 
c            IRNK - position in hadron chain
c     Output: IFL2 - flavor of second quark partner to be passed on
c             KF - final hadron
C-----------------------------------------------------------------------
Cf2py integer,intent(out) :: ifl2
Cf2py integer,intent(out) :: kf
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      SAVE

      INCLUDE 'sib_cflafr_cmmn.inc'
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'

      DIMENSION KFLA(4,4,2), CDIAG(16), KDIAG(8)
      DIMENSION KBAR(40), CFR(28), KFR(80)
      DATA KFLA /0,8,10,71,7,0,22,59,9,21,0,74,72,60,75,0, ! spin-zero mesons
     +     0,26,29,80,25,0,31,78,28,30,0,76,81,79,77,0/ ! spin-one mesons
      DATA CDIAG /0.5,0.25,0.5,0.25,1.,0.5,2.,1., ! spin-zero diagonal mesons
     +     0.5,0.,0.5,0.,1.,1.,2.,1./ ! spin-one diagonal mesons
      DATA KDIAG /6,23,24,73,27,32,33,83/
      DATA KBAR /13,14,34,35,36,37,38,84,85,86, ! jetset -> sibyll part. code map
     +     87,88,99,3*0,39,89,87,88, 
     +     40,41,42,43,44,45,46,47,48,49,      
     +     94,95,96,97,98,99,4*0/ ! spin-3/2 css baryon added to 1/2 css
      DATA CFR /0.75,0.,0.5,0.,0.,1.,0.1667,0.3333,0.0833,0.6667,0.1667,
     &     0.3333,-3.,1.,-2.,-2.,1.,0.,0.,-3.,1.,1.,1.,5*0./
      DATA KFR/0,16,17,19,100,104,109,115,0,26,27,29,122,126,131,137
     +  ,0,40,42,47,144,158,178,205,0,1,3,6,10,15,21,28,0,0,56,57,240,
     +  246,256,271,0,0,1,3,6,10,15,21,60,61,64,70,292,307,328,356,
     +  0,1,3,6,10,15,21,28,16*0/

      IF(NDEBUG.gt.6)
     &     WRITE(LUN,*)'  SIB_FLAV: input:',IFL1, IFL2_A, IRNK, IFL2, KF

c     set rho0 / ( omega, phi ) ratio, i.e. I=1 to I=0
c     default: 0.5, 0.0 ( phi only created from s-sbar)
      CDIAG(8+1) =  ONE-PAR(143) ! u-flavor, Prob. I=1 vs 0
      CDIAG(8+3) =  ONE-PAR(143) ! d-flavor, Prob. I=1 vs 0

      IARNK = IABS(IRNK)
      IFLA = IABS(IFL1)
      IFL2A = IFL2_A
      IF (IFL2A .NE. 0)  THEN
c     combine existing flavors to hadron
         IFL2A = MOD(IFL2A,100)
         IFL2 = IFL2A
         IFLB = IABS(IFL2A)
         MB = 0
         IF (IFLB .GT. 10)   MB=1
         IF (IFLA .GT. 10)   MB=2
      ELSE
c     sample new flavor
         MB = 2
         IF (IFLA .LT. 10)   THEN
             MB = 1
             IF ((ONE+PAR(1))*S_RNDM(0).LT. ONE)  MB=0
c     suppress baryons close to the string end
c     IPAR(55) defines largest forbidden rank
c     PAR(101) is the rejection probability
             IF (IPAR(54).eq.1)THEN
                IF(IARNK.le.IPAR(55).and.S_RNDM(0).lt.PAR(101)) MB=0
             ENDIF
         ENDIF
      ENDIF

 50   IF (MB .EQ. 0)  THEN
c     flavor open, sample from u,d,s,c
         IF (IFL2A.EQ.0)THEN
            IF(IPAR(69).eq.2)THEN
c     asymmetric between u,d
               IFL2 = MIN(2,1+INT((TWO+PAR(115))*S_RNDM(0)))
               IFLS = 3*(INT((2+PAR(2))*S_RNDM(0))/2)
               IFL2 = MAX(IFL2,IFLS)
               IFL2 = ISIGN(IFL2,-IFL1)
            ELSE
c     symmetric in u,d
               IFL2=ISIGN(1+INT((TWO+PAR(2))*S_RNDM(0)),-IFL1)
            ENDIF
            IF(IABS(IFL2).eq.3) THEN
               IF(S_RNDM(0).lt.PAR(24)*PAR(125))
     +              IFL2=ISIGN(4,-IFL1)
            ENDIF
         ENDIF
         IFLD = MAX(IFL1,IFL2)
         IFLE = MIN(IFL1,IFL2)
         GOTO 100
      ENDIF

C...Decide if the diquark must be split
      IF (MB .EQ. 2 .AND. IFLA .GT. 100)   THEN
         IFLA = MOD(IFLA,100)
           GOTO 200
      ENDIF
      IF (MB .EQ. 2 .AND. IFL2A .EQ. 0)   THEN
          IF (S_RNDM(0) .LT. PAR(8))  THEN
             MB = 0
             IFLG = MOD(IFL1,10)
             IFLH =(IFL1-IFLG)/10
             IF (S_RNDM(0) .GT. HALF)  THEN
                IFLDUM = IFLG
                IFLG = IFLH
                IFLH = IFLDUM
             ENDIF
             IFL11=IFLG
             IFL22=ISIGN(1+INT((TWO+PAR(2))*S_RNDM(0)),-IFL1)
             IFLD = MAX(IFL11,IFL22)
             IFLE = MIN(IFL11,IFL22)
             IFL2 = -IFLH*10+IFL22
             IF (S_RNDM(0) .GT. HALF)  IFL2 = IFL22*10-IFLH
             IFL2 = IFL2+ISIGN(100,IFL2)
          ENDIF
      ENDIF
       
C...Form a meson: consider spin and flavor mixing for the diagonal states
 100  IF (MB .EQ. 0)  THEN
         IF1 = IABS(IFLD)
         IF2 = IABS(IFLE)
         IFLC = MAX(IF1,IF2)
         KSP = INT(PAR(5)+S_RNDM(0))
         KSP = MIN(KSP,1)
         IF (IFLC.EQ.3)  KSP = INT(PAR(6)+S_RNDM(0))
         IF (IFLC.EQ.4)  KSP = INT(PAR(6)+S_RNDM(0))
         IF (IF1 .NE. IF2)   THEN
            KF = KFLA(IF1,IF2,KSP+1)
         ELSE
            R = S_RNDM(0)
            JF=1+INT(R+CDIAG(8*KSP+2*IF1-1))+
     +             INT(R+CDIAG(8*KSP+2*IF1))
            JF = MIN(JF,4)
            KF=KDIAG(JF+4*KSP)
c     suppress neutral pions
            IF(KF.eq.6)THEN
               IF(IPAR(82).eq.1.and.
     +              S_RNDM(kf).lt.PAR(137))then
                  IF(IFL2A.ne.0) goto 100
                  IF(IFLA.gt.10) mb = 2                  
                  GOTO 50
               endif
c     suppress neutral pions, depending on rank
               IF(IPAR(82).eq.2.and.S_RNDM(kf).lt.PAR(137).and.
     +              irnk.gt.0.and.irnk.lt.2) then
                  IF(IFL2A.ne.0) goto 100
                  IF(IFLA.gt.10) mb = 2
                  GOTO 50
               endif
            ENDIF
c     suppress rank1 (leading) omega
            IF(KF.eq.32)THEN
               IF(IPAR(83).ne.0.and.
     +              S_RNDM(kf).lt.PAR(138))then
                  IF(IFL2A.ne.0) goto 100
                  IF(IFLA.gt.10) mb = 2                  
                  GOTO 50
               endif
            ENDIF
         ENDIF
c         PRINT*,' I4FLAV returns :(IFL1,IFL2,LL)',IFL1,IFL2,KF
         RETURN
      ENDIF

C...Form a baryon
 200  IF (IFL2A .NE. 0)   THEN
         IF (MB .EQ. 1)  THEN
            IFLD = IFLA
            IFLE = IFLB/10
            IFLF = MOD(IFLB,10)
         ELSE
            IFLD = IFLB
            IFLE = IFLA/10
            IFLF = MOD(IFLA,10)
         ENDIF
         LFR = 3+2*((2*(IFLE-IFLF))/(1+IABS(IFLE-IFLF)))
         IF(IFLD.NE.IFLE.AND.IFLD.NE.IFLF)  LFR=LFR+1
      ELSE
 110     CONTINUE
         IF(MB.EQ.1)   THEN     ! generate diquark
            IFLD = IFLA
 120        IFLE = 1+INT((TWO+PAR(2)*PAR(3))*S_RNDM(0))
            IFLF = 1+INT((TWO+PAR(2)*PAR(3))*S_RNDM(0))          
            IF(IFLD.NE.4)THEN
               IF(IFLE.EQ.3)THEN 
                  IF(S_RNDM(0).lt.PAR(24)*PAR(125))
     +                 IFLE=4
               ENDIF
               IF(IFLF.EQ.3.and.IFLE.NE.4)THEN 
                  IF(S_RNDM(0).lt.PAR(24)*PAR(125))
     +                 IFLF=4
               ENDIF
            ENDIF
            IF(IFLE.GE.IFLF.AND.PAR(4).LT.S_RNDM(0))    GOTO 120
            IF(IFLE.LT.IFLF.AND.PAR(4)*S_RNDM(0).GT.ONE) GOTO 120     
            IFL2=ISIGN(10*IFLE+IFLF,IFL1)
         ELSE                   ! generate quark
            IF(IPAR(69).eq.2)THEN
c     asymmetric between u,d
               IFL2 = MIN(2,1+INT((TWO+PAR(115))*S_RNDM(0)))
               IFLS = 3*(INT((2+PAR(2))*S_RNDM(0))/2)
               IFL2 = MAX(IFL2,IFLS)
               IFL2 = ISIGN(IFL2,IFL1)
            ELSE
c     symmetric in u,d
               IFL2=ISIGN(1+INT((TWO+PAR(2))*S_RNDM(0)),IFL1)
            ENDIF
            IFLE=IFLA/10
            IFLF=MOD(IFLA,10)
            IF(IABS(IFL2).EQ.3.and.IFLF.ne.4.and.IFLE.ne.4) THEN
               IF(S_RNDM(0).lt.PAR(24)*PAR(125))
     +              IFL2=ISIGN(4,IFL1)
            ENDIF
            IFLD=IABS(IFL2)
         ENDIF
C...SU(6) factors for baryon formation
         LFR=3+2*((2*(IFLE-IFLF))/(1+IABS(IFLE-IFLF)))
         IF(IFLD.NE.IFLE.AND.IFLD.NE.IFLF)  LFR=LFR+1
         WT = CFR(2*LFR-1)+PAR(7)*CFR(2*LFR)
         IF(IFLE.LT.IFLF)  WT=WT/THREE
         IF (WT.LT.S_RNDM(0)) GOTO 110
      ENDIF

C...Form Baryon
      IFLG=MAX(IFLD,IFLE,IFLF)
      IFLI=MIN(IFLD,IFLE,IFLF)
      IFLH=IFLD+IFLE+IFLF-IFLG-IFLI
c      IF(IFLG+IFLH.gt.7) GOTO 200 ! forbid double charmed
      KSP=2+2*INT(ONE-CFR(2*LFR-1)+(CFR(2*LFR-1)+PAR(7)*
     1       CFR(2*LFR))*S_RNDM(0))

C...Distinguish Lambda- and Sigma- like particles
      IF (KSP.EQ.2.AND.IFLG.GT.IFLH.AND.IFLH.GT.IFLI)  THEN
         IF(IFLE.GT.IFLF.AND.IFLD.NE.IFLG) KSP=2+INT(0.75D0+S_RNDM(0))
         IF(IFLE.LT.IFLF.AND.IFLD.EQ.IFLG) KSP=3
         IF(IFLE.LT.IFLF.AND.IFLD.NE.IFLG) KSP=2+INT(0.25D0+S_RNDM(0))
      ENDIF
      KF=KFR(16*KSP-16+IFLG)+KFR(16*KSP-8+IFLH)+IFLI
      IF(KBAR(KF-40).eq.0)
     +     WRITE(LUN,*)'jetset code missing,flvs:',kf,IFLG,IFLH,IFLI
      KF=KBAR(KF-40)
      IF(KF.le.14)THEN
         IF(PAR(106).gt.S_RNDM(0).and.IARNK.le.IPAR(61)) KF=KF-13+51
     &        +2*INT(PAR(108)+S_RNDM(0))
      ENDIF
      KF=ISIGN(KF,IFL1)
      
      RETURN
      END
