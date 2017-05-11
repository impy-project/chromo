C*********************************************************************
 
C...PYDECY
C...Handles the decay of unstable particles.
 
      SUBROUTINE PYDECY(IP)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(4000,2),BRAT(4000),KFDP(4000,5)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/,/PYDAT3/
C...Local arrays.
      DIMENSION VDCY(4),KFLO(4),KFL1(4),PV(10,5),RORD(10),UE(3),BE(3),
     &WTCOR(10),PTAU(4),PCMTAU(4),DBETAU(3)
      CHARACTER CIDC*4
      DATA WTCOR/2D0,5D0,15D0,60D0,250D0,1500D0,1.2D4,1.2D5,150D0,16D0/
 
C...Functions: momentum in two-particle decays and four-product.
      PAWT(A,B,C)=SQRT((A**2-(B+C)**2)*(A**2-(B-C)**2))/(2D0*A)
      FOUR(I,J)=P(I,4)*P(J,4)-P(I,1)*P(J,1)-P(I,2)*P(J,2)-P(I,3)*P(J,3)
 
C...Initial values.
      NTRY=0
      NSAV=N
      KFA=IABS(K(IP,2))
      KFS=ISIGN(1,K(IP,2))
      KC=PYCOMP(KFA)
      MSTJ(92)=0

C...Choose lifetime and determine decay vertex.
      IF(K(IP,1).EQ.5) THEN
        V(IP,5)=0D0
      ELSEIF(K(IP,1).NE.4) THEN
        V(IP,5)=-PMAS(KC,4)*LOG(PYR(0))
      ENDIF
      DO 100 J=1,4
        VDCY(J)=V(IP,J)+V(IP,5)*P(IP,J)/P(IP,5)
  100 CONTINUE
 
C...Determine whether decay allowed or not.
      MOUT=0
      IF(MSTJ(22).EQ.2) THEN
        IF(PMAS(KC,4).GT.PARJ(71)) MOUT=1
      ELSEIF(MSTJ(22).EQ.3) THEN
        IF(VDCY(1)**2+VDCY(2)**2+VDCY(3)**2.GT.PARJ(72)**2) MOUT=1
      ELSEIF(MSTJ(22).EQ.4) THEN
        IF(VDCY(1)**2+VDCY(2)**2.GT.PARJ(73)**2) MOUT=1
        IF(ABS(VDCY(3)).GT.PARJ(74)) MOUT=1
      ENDIF
      IF(MOUT.EQ.1.AND.K(IP,1).NE.5) THEN
        K(IP,1)=4
        RETURN
      ENDIF
 
C...Interface to external tau decay library (for tau polarization).
      IF(KFA.EQ.15.AND.MSTJ(28).GE.1) THEN
 
C...Starting values for pointers and momenta.
        ITAU=IP
        DO 110 J=1,4
          PTAU(J)=P(ITAU,J)
          PCMTAU(J)=P(ITAU,J)
  110   CONTINUE
 
C...Iterate to find position and code of mother of tau.
        IMTAU=ITAU
  120   IMTAU=K(IMTAU,3)
 
        IF(IMTAU.EQ.0) THEN
C...If no known origin then impossible to do anything further.
          KFORIG=0
          IORIG=0
 
        ELSEIF(K(IMTAU,2).EQ.K(ITAU,2)) THEN
C...If tau -> tau + gamma then add gamma energy and loop.
          IF(K(K(IMTAU,4),2).EQ.22) THEN
            DO 130 J=1,4
              PCMTAU(J)=PCMTAU(J)+P(K(IMTAU,4),J)
  130       CONTINUE
          ELSEIF(K(K(IMTAU,5),2).EQ.22) THEN
            DO 140 J=1,4
              PCMTAU(J)=PCMTAU(J)+P(K(IMTAU,5),J)
  140       CONTINUE
          ENDIF
          GOTO 120
 
        ELSEIF(IABS(K(IMTAU,2)).GT.100) THEN
C...If coming from weak decay of hadron then W is not stored in record,
C...but can be reconstructed by adding neutrino momentum.
          KFORIG=-ISIGN(24,K(ITAU,2))
          IORIG=0
          DO 160 II=K(IMTAU,4),K(IMTAU,5)
            IF(K(II,2)*ISIGN(1,K(ITAU,2)).EQ.-16) THEN
              DO 150 J=1,4
                PCMTAU(J)=PCMTAU(J)+P(II,J)
  150         CONTINUE
            ENDIF
  160     CONTINUE
 
        ELSE
C...If coming from resonance decay then find latest copy of this
C...resonance (may not completely agree).
          KFORIG=K(IMTAU,2)
          IORIG=IMTAU
          DO 170 II=IMTAU+1,IP-1
            IF(K(II,2).EQ.KFORIG.AND.K(II,3).EQ.IORIG.AND.
     &      ABS(P(II,5)-P(IORIG,5)).LT.1D-5*P(IORIG,5)) IORIG=II
  170     CONTINUE
          DO 180 J=1,4
            PCMTAU(J)=P(IORIG,J)
  180     CONTINUE
        ENDIF
 
C...Boost tau to rest frame of production process (where known)
C...and rotate it to sit along +z axis.
        DO 190 J=1,3
          DBETAU(J)=PCMTAU(J)/PCMTAU(4)
  190   CONTINUE
        IF(KFORIG.NE.0) CALL PYROBO(ITAU,ITAU,0D0,0D0,-DBETAU(1),
     &  -DBETAU(2),-DBETAU(3))
        PHITAU=PYANGL(P(ITAU,1),P(ITAU,2))
        CALL PYROBO(ITAU,ITAU,0D0,-PHITAU,0D0,0D0,0D0)
        THETAU=PYANGL(P(ITAU,3),P(ITAU,1))
        CALL PYROBO(ITAU,ITAU,-THETAU,0D0,0D0,0D0,0D0)
 
C...Call tau decay routine (if meaningful) and fill extra info.
        IF(KFORIG.NE.0.OR.MSTJ(28).EQ.2) THEN
          CALL PYTAUD(ITAU,IORIG,KFORIG,NDECAY)
          DO 200 II=NSAV+1,NSAV+NDECAY
            K(II,1)=1
            K(II,3)=IP
            K(II,4)=0
            K(II,5)=0
  200     CONTINUE
          N=NSAV+NDECAY
        ENDIF
 
C...Boost back decay tau and decay products.
        DO 210 J=1,4
          P(ITAU,J)=PTAU(J)
  210   CONTINUE
        IF(KFORIG.NE.0.OR.MSTJ(28).EQ.2) THEN
          CALL PYROBO(NSAV+1,N,THETAU,PHITAU,0D0,0D0,0D0)
          IF(KFORIG.NE.0) CALL PYROBO(NSAV+1,N,0D0,0D0,DBETAU(1),
     &    DBETAU(2),DBETAU(3))
 
C...Skip past ordinary tau decay treatment.
          MMAT=0
          MBST=0
          ND=0
          GOTO 630
        ENDIF
      ENDIF
 
C...B-Bbar mixing: flip sign of meson appropriately.
      MMIX=0
      IF((KFA.EQ.511.OR.KFA.EQ.531).AND.MSTJ(26).GE.1) THEN
        XBBMIX=PARJ(76)
        IF(KFA.EQ.531) XBBMIX=PARJ(77)
        IF(SIN(0.5D0*XBBMIX*V(IP,5)/PMAS(KC,4))**2.GT.PYR(0)) MMIX=1
        IF(MMIX.EQ.1) KFS=-KFS
      ENDIF
 
C...Check existence of decay channels. Particle/antiparticle rules.
      KCA=KC
      IF(MDCY(KC,2).GT.0) THEN
        MDMDCY=MDME(MDCY(KC,2),2)
        IF(MDMDCY.GT.80.AND.MDMDCY.LE.90) KCA=MDMDCY
      ENDIF
      IF(MDCY(KCA,2).LE.0.OR.MDCY(KCA,3).LE.0) THEN
        CALL PYERRM(9,'(PYDECY:) no decay channel defined')
        RETURN
      ENDIF
      IF(MOD(KFA/1000,10).EQ.0.AND.KCA.EQ.85) KFS=-KFS
      IF(KCHG(KC,3).EQ.0) THEN
        KFSP=1
        KFSN=0
        IF(PYR(0).GT.0.5D0) KFS=-KFS
      ELSEIF(KFS.GT.0) THEN
        KFSP=1
        KFSN=0
      ELSE
        KFSP=0
        KFSN=1
      ENDIF
 
C...Sum branching ratios of allowed decay channels.
  220 NOPE=0
      BRSU=0D0
      DO 230 IDL=MDCY(KCA,2),MDCY(KCA,2)+MDCY(KCA,3)-1
        IF(MDME(IDL,1).NE.1.AND.KFSP*MDME(IDL,1).NE.2.AND.
     &  KFSN*MDME(IDL,1).NE.3) GOTO 230
        IF(MDME(IDL,2).GT.100) GOTO 230
        NOPE=NOPE+1
        BRSU=BRSU+BRAT(IDL)
  230 CONTINUE
      IF(NOPE.EQ.0) THEN
        CALL PYERRM(2,'(PYDECY:) all decay channels closed by user')
        RETURN
      ENDIF
 
C...Select decay channel among allowed ones.
  240 RBR=BRSU*PYR(0)
      IDL=MDCY(KCA,2)-1
  250 IDL=IDL+1
      IF(MDME(IDL,1).NE.1.AND.KFSP*MDME(IDL,1).NE.2.AND.
     &KFSN*MDME(IDL,1).NE.3) THEN
        IF(IDL.LT.MDCY(KCA,2)+MDCY(KCA,3)-1) GOTO 250
      ELSEIF(MDME(IDL,2).GT.100) THEN
        IF(IDL.LT.MDCY(KCA,2)+MDCY(KCA,3)-1) GOTO 250
      ELSE
        IDC=IDL
        RBR=RBR-BRAT(IDL)
        IF(IDL.LT.MDCY(KCA,2)+MDCY(KCA,3)-1.AND.RBR.GT.0D0) GOTO 250
      ENDIF
 
C...Start readout of decay channel: matrix element, reset counters.
      MMAT=MDME(IDC,2)
  260 NTRY=NTRY+1
      IF(MOD(NTRY,200).EQ.0) THEN
        WRITE(CIDC,'(I4)') IDC
        CALL PYERRM(4,'(PYDECY:) caught in loop for decay channel'//
     &  CIDC)
        GOTO 240
      ENDIF
      IF(NTRY.GT.1000) THEN
        CALL PYERRM(14,'(PYDECY:) caught in infinite loop')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
      I=N
      NP=0
      NQ=0
      MBST=0
      IF(MMAT.GE.11.AND.P(IP,4).GT.20D0*P(IP,5)) MBST=1
      DO 270 J=1,4
        PV(1,J)=0D0
        IF(MBST.EQ.0) PV(1,J)=P(IP,J)
  270 CONTINUE
      IF(MBST.EQ.1) PV(1,4)=P(IP,5)
      PV(1,5)=P(IP,5)
      PS=0D0
      PSQ=0D0
      MREM=0
      MHADDY=0
      IF(KFA.GT.80) MHADDY=1
C.. Random flavour and popcorn system memory.
      IRNDMO=0
      JTMO=0
      MSTU(121)=0
      MSTU(125)=10
 
C...Read out decay products. Convert to standard flavour code.
      JTMAX=5
      IF(MDME(IDC+1,2).EQ.101) JTMAX=10
      DO 280 JT=1,JTMAX
        IF(JT.LE.5) KP=KFDP(IDC,JT)
        IF(JT.GE.6) KP=KFDP(IDC+1,JT-5)
        IF(KP.EQ.0) GOTO 280
        KPA=IABS(KP)
        KCP=PYCOMP(KPA)
        IF(KPA.GT.80) MHADDY=1
        IF(KCHG(KCP,3).EQ.0.AND.KPA.NE.81.AND.KPA.NE.82) THEN
          KFP=KP
        ELSEIF(KPA.NE.81.AND.KPA.NE.82) THEN
          KFP=KFS*KP
        ELSEIF(KPA.EQ.81.AND.MOD(KFA/1000,10).EQ.0) THEN
          KFP=-KFS*MOD(KFA/10,10)
        ELSEIF(KPA.EQ.81.AND.MOD(KFA/100,10).GE.MOD(KFA/10,10)) THEN
          KFP=KFS*(100*MOD(KFA/10,100)+3)
        ELSEIF(KPA.EQ.81) THEN
          KFP=KFS*(1000*MOD(KFA/10,10)+100*MOD(KFA/100,10)+1)
        ELSEIF(KP.EQ.82) THEN
          CALL PYDCYK(-KFS*INT(1D0+(2D0+PARJ(2))*PYR(0)),0,KFP,KDUMP)
          IF(KFP.EQ.0) GOTO 260
          KFP=-KFP
          IRNDMO=1
          MSTJ(93)=1
          IF(PV(1,5).LT.PARJ(32)+2D0*PYMASS(KFP)) GOTO 260
        ELSEIF(KP.EQ.-82) THEN
          KFP=MSTU(124)
        ENDIF
        IF(KPA.EQ.81.OR.KPA.EQ.82) KCP=PYCOMP(KFP)
 
C...Add decay product to event record or to quark flavour list.
        KFPA=IABS(KFP)
        KQP=KCHG(KCP,2)
        IF(MMAT.GE.11.AND.MMAT.LE.30.AND.KQP.NE.0) THEN
          NQ=NQ+1
          KFLO(NQ)=KFP
C...set rndmflav popcorn system pointer
          IF(KP.EQ.82.AND.MSTU(121).GT.0) JTMO=NQ
          MSTJ(93)=2
          PSQ=PSQ+PYMASS(KFLO(NQ))
        ELSEIF((MMAT.EQ.42.OR.MMAT.EQ.43.OR.MMAT.EQ.48).AND.NP.EQ.3.AND.
     &    MOD(NQ,2).EQ.1) THEN
          NQ=NQ-1
          PS=PS-P(I,5)
          K(I,1)=1
          KFI=K(I,2)
          CALL PYKFDI(KFP,KFI,KFLDMP,K(I,2))
          IF(K(I,2).EQ.0) GOTO 260
          MSTJ(93)=1
          P(I,5)=PYMASS(K(I,2))
          PS=PS+P(I,5)
        ELSE
          I=I+1
          NP=NP+1
          IF(MMAT.NE.33.AND.KQP.NE.0) NQ=NQ+1
          IF(MMAT.EQ.33.AND.KQP.NE.0.AND.KQP.NE.2) NQ=NQ+1
          K(I,1)=1+MOD(NQ,2)
          IF(MMAT.EQ.4.AND.JT.LE.2.AND.KFP.EQ.21) K(I,1)=2
          IF(MMAT.EQ.4.AND.JT.EQ.3) K(I,1)=1
          K(I,2)=KFP
          K(I,3)=IP
          K(I,4)=0
          K(I,5)=0
          P(I,5)=PYMASS(KFP)
          PS=PS+P(I,5)
        ENDIF
  280 CONTINUE
 
C...Check masses for resonance decays.
      IF(MHADDY.EQ.0) THEN
        IF(PS+PARJ(64).GT.PV(1,5)) GOTO 240
      ENDIF
 
C...Choose decay multiplicity in phase space model.
  290 IF(MMAT.GE.11.AND.MMAT.LE.30) THEN
        PSP=PS
        CNDE=PARJ(61)*LOG(MAX((PV(1,5)-PS-PSQ)/PARJ(62),1.1D0))
        IF(MMAT.EQ.12) CNDE=CNDE+PARJ(63)
  300   NTRY=NTRY+1
C...Reset popcorn flags if new attempt. Re-select rndmflav if failed.
        IF(IRNDMO.EQ.0) THEN
           MSTU(121)=0
           JTMO=0
        ELSEIF(IRNDMO.EQ.1) THEN
           IRNDMO=2
        ELSE
           GOTO 260
        ENDIF
        IF(NTRY.GT.1000) THEN
          CALL PYERRM(14,'(PYDECY:) caught in infinite loop')
          IF(MSTU(21).GE.1) RETURN
        ENDIF
        IF(MMAT.LE.20) THEN
          GAUSS=SQRT(-2D0*CNDE*LOG(MAX(1D-10,PYR(0))))*
     &    SIN(PARU(2)*PYR(0))
          ND=0.5D0+0.5D0*NP+0.25D0*NQ+CNDE+GAUSS
          IF(ND.LT.NP+NQ/2.OR.ND.LT.2.OR.ND.GT.10) GOTO 300
          IF(MMAT.EQ.13.AND.ND.EQ.2) GOTO 300
          IF(MMAT.EQ.14.AND.ND.LE.3) GOTO 300
          IF(MMAT.EQ.15.AND.ND.LE.4) GOTO 300
        ELSE
          ND=MMAT-20
        ENDIF
C.. Set maximum popcorn meson number. Test rndmflav popcorn size.
        MSTU(125)=ND-NQ/2
        IF(MSTU(121).GT.MSTU(125)) GOTO 300
 
C...Form hadrons from flavour content.
        DO 310 JT=1,4
          KFL1(JT)=KFLO(JT)
  310   CONTINUE
        IF(ND.EQ.NP+NQ/2) GOTO 330
        DO 320 I=N+NP+1,N+ND-NQ/2
C.. Stick to started popcorn system, else pick side at random
          JT=JTMO
          IF(JT.EQ.0) JT=1+INT((NQ-1)*PYR(0))
          CALL PYDCYK(KFL1(JT),0,KFL2,K(I,2))
          IF(K(I,2).EQ.0) GOTO 300
          MSTU(125)=MSTU(125)-1
          JTMO=0
          IF(MSTU(121).GT.0) JTMO=JT
          KFL1(JT)=-KFL2
  320   CONTINUE
  330   JT=2
        JT2=3
        JT3=4
        IF(NQ.EQ.4.AND.PYR(0).LT.PARJ(66)) JT=4
        IF(JT.EQ.4.AND.ISIGN(1,KFL1(1)*(10-IABS(KFL1(1))))*
     &  ISIGN(1,KFL1(JT)*(10-IABS(KFL1(JT)))).GT.0) JT=3
        IF(JT.EQ.3) JT2=2
        IF(JT.EQ.4) JT3=2
        CALL PYDCYK(KFL1(1),KFL1(JT),KFLDMP,K(N+ND-NQ/2+1,2))
        IF(K(N+ND-NQ/2+1,2).EQ.0) GOTO 300
        IF(NQ.EQ.4) CALL PYDCYK(KFL1(JT2),KFL1(JT3),KFLDMP,K(N+ND,2))
        IF(NQ.EQ.4.AND.K(N+ND,2).EQ.0) GOTO 300
 
C...Check that sum of decay product masses not too large.
        PS=PSP
        DO 340 I=N+NP+1,N+ND
          K(I,1)=1
          K(I,3)=IP
          K(I,4)=0
          K(I,5)=0
          P(I,5)=PYMASS(K(I,2))
          PS=PS+P(I,5)
  340   CONTINUE
        IF(PS+PARJ(64).GT.PV(1,5)) GOTO 300
 
C...Rescale energy to subtract off spectator quark mass.
      ELSEIF((MMAT.EQ.31.OR.MMAT.EQ.33.OR.MMAT.EQ.44)
     &  .AND.NP.GE.3) THEN
        PS=PS-P(N+NP,5)
        PQT=(P(N+NP,5)+PARJ(65))/PV(1,5)
        DO 350 J=1,5
          P(N+NP,J)=PQT*PV(1,J)
          PV(1,J)=(1D0-PQT)*PV(1,J)
  350   CONTINUE
        IF(PS+PARJ(64).GT.PV(1,5)) GOTO 260
        ND=NP-1
        MREM=1
 
C...Fully specified final state: check mass broadening effects.
      ELSE
        IF(NP.GE.2.AND.PS+PARJ(64).GT.PV(1,5)) GOTO 260
        ND=NP
      ENDIF
 
C...Determine position of grandmother, number of sisters.
      NM=0
      KFAS=0
      MSGN=0
      IF(MMAT.EQ.3) THEN
        IM=K(IP,3)
        IF(IM.LT.0.OR.IM.GE.IP) IM=0
        IF(IM.NE.0) KFAM=IABS(K(IM,2))
        IF(IM.NE.0) THEN
          DO 360 IL=MAX(IP-2,IM+1),MIN(IP+2,N)
            IF(K(IL,3).EQ.IM) NM=NM+1
            IF(K(IL,3).EQ.IM.AND.IL.NE.IP) ISIS=IL
  360     CONTINUE
          IF(NM.NE.2.OR.KFAM.LE.100.OR.MOD(KFAM,10).NE.1.OR.
     &    MOD(KFAM/1000,10).NE.0) NM=0
          IF(NM.EQ.2) THEN
            KFAS=IABS(K(ISIS,2))
            IF((KFAS.LE.100.OR.MOD(KFAS,10).NE.1.OR.
     &      MOD(KFAS/1000,10).NE.0).AND.KFAS.NE.22) NM=0
          ENDIF
        ENDIF
      ENDIF
 
C...Kinematics of one-particle decays.
      IF(ND.EQ.1) THEN
        DO 370 J=1,4
          P(N+1,J)=P(IP,J)
  370   CONTINUE
        GOTO 630
      ENDIF
 
C...Calculate maximum weight ND-particle decay.
      PV(ND,5)=P(N+ND,5)
      IF(ND.GE.3) THEN
        WTMAX=1D0/WTCOR(ND-2)
        PMAX=PV(1,5)-PS+P(N+ND,5)
        PMIN=0D0
        DO 380 IL=ND-1,1,-1
          PMAX=PMAX+P(N+IL,5)
          PMIN=PMIN+P(N+IL+1,5)
          WTMAX=WTMAX*PAWT(PMAX,PMIN,P(N+IL,5))
  380   CONTINUE
      ENDIF
 
C...Find virtual gamma mass in Dalitz decay.
  390 IF(ND.EQ.2) THEN
      ELSEIF(MMAT.EQ.2) THEN
        PMES=4D0*PMAS(11,1)**2
        PMRHO2=PMAS(131,1)**2
        PGRHO2=PMAS(131,2)**2
  400   PMST=PMES*(P(IP,5)**2/PMES)**PYR(0)
        WT=(1+0.5D0*PMES/PMST)*SQRT(MAX(0D0,1D0-PMES/PMST))*
     &  (1D0-PMST/P(IP,5)**2)**3*(1D0+PGRHO2/PMRHO2)/
     &  ((1D0-PMST/PMRHO2)**2+PGRHO2/PMRHO2)
        IF(WT.LT.PYR(0)) GOTO 400
        PV(2,5)=MAX(2.00001D0*PMAS(11,1),SQRT(PMST))
 
C...M-generator gives weight. If rejected, try again.
      ELSE
  410   RORD(1)=1D0
        DO 440 IL1=2,ND-1
          RSAV=PYR(0)
          DO 420 IL2=IL1-1,1,-1
            IF(RSAV.LE.RORD(IL2)) GOTO 430
            RORD(IL2+1)=RORD(IL2)
  420     CONTINUE
  430     RORD(IL2+1)=RSAV
  440   CONTINUE
        RORD(ND)=0D0
        WT=1D0
        DO 450 IL=ND-1,1,-1
          PV(IL,5)=PV(IL+1,5)+P(N+IL,5)+(RORD(IL)-RORD(IL+1))*
     &    (PV(1,5)-PS)
          WT=WT*PAWT(PV(IL,5),PV(IL+1,5),P(N+IL,5))
  450   CONTINUE
        IF(WT.LT.PYR(0)*WTMAX) GOTO 410
      ENDIF
 
C...Perform two-particle decays in respective CM frame.
  460 DO 480 IL=1,ND-1
        PA=PAWT(PV(IL,5),PV(IL+1,5),P(N+IL,5))
        UE(3)=2D0*PYR(0)-1D0
        PHI=PARU(2)*PYR(0)
        UE(1)=SQRT(1D0-UE(3)**2)*COS(PHI)
        UE(2)=SQRT(1D0-UE(3)**2)*SIN(PHI)
        DO 470 J=1,3
          P(N+IL,J)=PA*UE(J)
          PV(IL+1,J)=-PA*UE(J)
  470   CONTINUE
        P(N+IL,4)=SQRT(PA**2+P(N+IL,5)**2)
        PV(IL+1,4)=SQRT(PA**2+PV(IL+1,5)**2)
  480 CONTINUE
 
C...Lorentz transform decay products to lab frame.
      DO 490 J=1,4
        P(N+ND,J)=PV(ND,J)
  490 CONTINUE
      DO 530 IL=ND-1,1,-1
        DO 500 J=1,3
          BE(J)=PV(IL,J)/PV(IL,4)
  500   CONTINUE
        GA=PV(IL,4)/PV(IL,5)
        DO 520 I=N+IL,N+ND
          BEP=BE(1)*P(I,1)+BE(2)*P(I,2)+BE(3)*P(I,3)
          DO 510 J=1,3
            P(I,J)=P(I,J)+GA*(GA*BEP/(1D0+GA)+P(I,4))*BE(J)
  510     CONTINUE
          P(I,4)=GA*(P(I,4)+BEP)
  520   CONTINUE
  530 CONTINUE
 
C...Check that no infinite loop in matrix element weight.
      NTRY=NTRY+1
      IF(NTRY.GT.800) GOTO 560
 
C...Matrix elements for omega and phi decays.
      IF(MMAT.EQ.1) THEN
        WT=(P(N+1,5)*P(N+2,5)*P(N+3,5))**2-(P(N+1,5)*FOUR(N+2,N+3))**2
     &  -(P(N+2,5)*FOUR(N+1,N+3))**2-(P(N+3,5)*FOUR(N+1,N+2))**2
     &  +2D0*FOUR(N+1,N+2)*FOUR(N+1,N+3)*FOUR(N+2,N+3)
        IF(MAX(WT*WTCOR(9)/P(IP,5)**6,0.001D0).LT.PYR(0)) GOTO 390
 
C...Matrix elements for pi0 or eta Dalitz decay to gamma e+ e-.
      ELSEIF(MMAT.EQ.2) THEN
        FOUR12=FOUR(N+1,N+2)
        FOUR13=FOUR(N+1,N+3)
        WT=(PMST-0.5D0*PMES)*(FOUR12**2+FOUR13**2)+
     &  PMES*(FOUR12*FOUR13+FOUR12**2+FOUR13**2)
        IF(WT.LT.PYR(0)*0.25D0*PMST*(P(IP,5)**2-PMST)**2) GOTO 460
 
C...Matrix element for S0 -> S1 + V1 -> S1 + S2 + S3 (S scalar,
C...V vector), of form cos**2(theta02) in V1 rest frame, and for
C...S0 -> gamma + V1 -> gamma + S2 + S3, of form sin**2(theta02).
      ELSEIF(MMAT.EQ.3.AND.NM.EQ.2) THEN
        FOUR10=FOUR(IP,IM)
        FOUR12=FOUR(IP,N+1)
        FOUR02=FOUR(IM,N+1)
        PMS1=P(IP,5)**2
        PMS0=P(IM,5)**2
        PMS2=P(N+1,5)**2
        IF(KFAS.NE.22) HNUM=(FOUR10*FOUR12-PMS1*FOUR02)**2
        IF(KFAS.EQ.22) HNUM=PMS1*(2D0*FOUR10*FOUR12*FOUR02-
     &  PMS1*FOUR02**2-PMS0*FOUR12**2-PMS2*FOUR10**2+PMS1*PMS0*PMS2)
        HNUM=MAX(1D-6*PMS1**2*PMS0*PMS2,HNUM)
        HDEN=(FOUR10**2-PMS1*PMS0)*(FOUR12**2-PMS1*PMS2)
        IF(HNUM.LT.PYR(0)*HDEN) GOTO 460
 
C...Matrix element for "onium" -> g + g + g or gamma + g + g.
      ELSEIF(MMAT.EQ.4) THEN
        HX1=2D0*FOUR(IP,N+1)/P(IP,5)**2
        HX2=2D0*FOUR(IP,N+2)/P(IP,5)**2
        HX3=2D0*FOUR(IP,N+3)/P(IP,5)**2
        WT=((1D0-HX1)/(HX2*HX3))**2+((1D0-HX2)/(HX1*HX3))**2+
     &  ((1D0-HX3)/(HX1*HX2))**2
        IF(WT.LT.2D0*PYR(0)) GOTO 390
        IF(K(IP+1,2).EQ.22.AND.(1D0-HX1)*P(IP,5)**2.LT.4D0*PARJ(32)**2)
     &  GOTO 390
 
C...Effective matrix element for nu spectrum in tau -> nu + hadrons.
      ELSEIF(MMAT.EQ.41) THEN
        HX1=2D0*FOUR(IP,N+1)/P(IP,5)**2
        HXM=MIN(0.75D0,2D0*(1D0-PS/P(IP,5)))
        IF(HX1*(3D0-2D0*HX1).LT.PYR(0)*HXM*(3D0-2D0*HXM)) GOTO 390
 
C...Matrix elements for weak decays (only semileptonic for c and b)
      ELSEIF((MMAT.EQ.42.OR.MMAT.EQ.43.OR.MMAT.EQ.44.OR.MMAT.EQ.48)
     &  .AND.ND.EQ.3) THEN
        IF(MBST.EQ.0) WT=FOUR(IP,N+1)*FOUR(N+2,N+3)
        IF(MBST.EQ.1) WT=P(IP,5)*P(N+1,4)*FOUR(N+2,N+3)
        IF(WT.LT.PYR(0)*P(IP,5)*PV(1,5)**3/WTCOR(10)) GOTO 390
      ELSEIF(MMAT.EQ.42.OR.MMAT.EQ.43.OR.MMAT.EQ.44.OR.MMAT.EQ.48) THEN
        DO 550 J=1,4
          P(N+NP+1,J)=0D0
          DO 540 IS=N+3,N+NP
            P(N+NP+1,J)=P(N+NP+1,J)+P(IS,J)
  540     CONTINUE
  550   CONTINUE
        IF(MBST.EQ.0) WT=FOUR(IP,N+1)*FOUR(N+2,N+NP+1)
        IF(MBST.EQ.1) WT=P(IP,5)*P(N+1,4)*FOUR(N+2,N+NP+1)
        IF(WT.LT.PYR(0)*P(IP,5)*PV(1,5)**3/WTCOR(10)) GOTO 390
      ENDIF
 
C...Scale back energy and reattach spectator.
  560 IF(MREM.EQ.1) THEN
        DO 570 J=1,5
          PV(1,J)=PV(1,J)/(1D0-PQT)
  570   CONTINUE
        ND=ND+1
        MREM=0
      ENDIF
 
C...Low invariant mass for system with spectator quark gives particle,
C...not two jets. Readjust momenta accordingly.
      IF(MMAT.EQ.31.AND.ND.EQ.3) THEN
        MSTJ(93)=1
        PM2=PYMASS(K(N+2,2))
        MSTJ(93)=1
        PM3=PYMASS(K(N+3,2))
        IF(P(N+2,5)**2+P(N+3,5)**2+2D0*FOUR(N+2,N+3).GE.
     &  (PARJ(32)+PM2+PM3)**2) GOTO 630
        K(N+2,1)=1
        KFTEMP=K(N+2,2)
        CALL PYKFDI(KFTEMP,K(N+3,2),KFLDMP,K(N+2,2))
        IF(K(N+2,2).EQ.0) GOTO 260
        P(N+2,5)=PYMASS(K(N+2,2))
        PS=P(N+1,5)+P(N+2,5)
        PV(2,5)=P(N+2,5)
        MMAT=0
        ND=2
        GOTO 460
      ELSEIF(MMAT.EQ.44) THEN
        MSTJ(93)=1
        PM3=PYMASS(K(N+3,2))
        MSTJ(93)=1
        PM4=PYMASS(K(N+4,2))
        IF(P(N+3,5)**2+P(N+4,5)**2+2D0*FOUR(N+3,N+4).GE.
     &  (PARJ(32)+PM3+PM4)**2) GOTO 600
        K(N+3,1)=1
        KFTEMP=K(N+3,2)
        CALL PYKFDI(KFTEMP,K(N+4,2),KFLDMP,K(N+3,2))
        IF(K(N+3,2).EQ.0) GOTO 260
        P(N+3,5)=PYMASS(K(N+3,2))
        DO 580 J=1,3
          P(N+3,J)=P(N+3,J)+P(N+4,J)
  580   CONTINUE
        P(N+3,4)=SQRT(P(N+3,1)**2+P(N+3,2)**2+P(N+3,3)**2+P(N+3,5)**2)
        HA=P(N+1,4)**2-P(N+2,4)**2
        HB=HA-(P(N+1,5)**2-P(N+2,5)**2)
        HC=(P(N+1,1)-P(N+2,1))**2+(P(N+1,2)-P(N+2,2))**2+
     &  (P(N+1,3)-P(N+2,3))**2
        HD=(PV(1,4)-P(N+3,4))**2
        HE=HA**2-2D0*HD*(P(N+1,4)**2+P(N+2,4)**2)+HD**2
        HF=HD*HC-HB**2
        HG=HD*HC-HA*HB
        HH=(SQRT(HG**2+HE*HF)-HG)/(2D0*HF)
        DO 590 J=1,3
          PCOR=HH*(P(N+1,J)-P(N+2,J))
          P(N+1,J)=P(N+1,J)+PCOR
          P(N+2,J)=P(N+2,J)-PCOR
  590   CONTINUE
        P(N+1,4)=SQRT(P(N+1,1)**2+P(N+1,2)**2+P(N+1,3)**2+P(N+1,5)**2)
        P(N+2,4)=SQRT(P(N+2,1)**2+P(N+2,2)**2+P(N+2,3)**2+P(N+2,5)**2)
        ND=ND-1
      ENDIF
 
C...Check invariant mass of W jets. May give one particle or start over.
  600 IF((MMAT.EQ.42.OR.MMAT.EQ.43.OR.MMAT.EQ.44.OR.MMAT.EQ.48)
     &.AND.IABS(K(N+1,2)).LT.10) THEN
        PMR=SQRT(MAX(0D0,P(N+1,5)**2+P(N+2,5)**2+2D0*FOUR(N+1,N+2)))
        MSTJ(93)=1
        PM1=PYMASS(K(N+1,2))
        MSTJ(93)=1
        PM2=PYMASS(K(N+2,2))
        IF(PMR.GT.PARJ(32)+PM1+PM2) GOTO 610
        KFLDUM=INT(1.5D0+PYR(0))
        CALL PYKFDI(K(N+1,2),-ISIGN(KFLDUM,K(N+1,2)),KFLDMP,KF1)
        CALL PYKFDI(K(N+2,2),-ISIGN(KFLDUM,K(N+2,2)),KFLDMP,KF2)
        IF(KF1.EQ.0.OR.KF2.EQ.0) GOTO 260
        PSM=PYMASS(KF1)+PYMASS(KF2)
        IF((MMAT.EQ.42.OR.MMAT.EQ.48).AND.PMR.GT.PARJ(64)+PSM) GOTO 610
        IF(MMAT.GE.43.AND.PMR.GT.0.2D0*PARJ(32)+PSM) GOTO 610
        IF(MMAT.EQ.48) GOTO 390
        IF(ND.EQ.4.OR.KFA.EQ.15) GOTO 260
        K(N+1,1)=1
        KFTEMP=K(N+1,2)
        CALL PYKFDI(KFTEMP,K(N+2,2),KFLDMP,K(N+1,2))
        IF(K(N+1,2).EQ.0) GOTO 260
        P(N+1,5)=PYMASS(K(N+1,2))
        K(N+2,2)=K(N+3,2)
        P(N+2,5)=P(N+3,5)
        PS=P(N+1,5)+P(N+2,5)
        IF(PS+PARJ(64).GT.PV(1,5)) GOTO 260
        PV(2,5)=P(N+3,5)
        MMAT=0
        ND=2
        GOTO 460
      ENDIF
 
C...Phase space decay of partons from W decay.
  610 IF((MMAT.EQ.42.OR.MMAT.EQ.48).AND.IABS(K(N+1,2)).LT.10) THEN
        KFLO(1)=K(N+1,2)
        KFLO(2)=K(N+2,2)
        K(N+1,1)=K(N+3,1)
        K(N+1,2)=K(N+3,2)
        DO 620 J=1,5
          PV(1,J)=P(N+1,J)+P(N+2,J)
          P(N+1,J)=P(N+3,J)
  620   CONTINUE
        PV(1,5)=PMR
        N=N+1
        NP=0
        NQ=2
        PS=0D0
        MSTJ(93)=2
        PSQ=PYMASS(KFLO(1))
        MSTJ(93)=2
        PSQ=PSQ+PYMASS(KFLO(2))
        MMAT=11
        GOTO 290
      ENDIF
 
C...Boost back for rapidly moving particle.
  630 N=N+ND
      IF(MBST.EQ.1) THEN
        DO 640 J=1,3
          BE(J)=P(IP,J)/P(IP,4)
  640   CONTINUE
        GA=P(IP,4)/P(IP,5)
        DO 660 I=NSAV+1,N
          BEP=BE(1)*P(I,1)+BE(2)*P(I,2)+BE(3)*P(I,3)
          DO 650 J=1,3
            P(I,J)=P(I,J)+GA*(GA*BEP/(1D0+GA)+P(I,4))*BE(J)
  650     CONTINUE
          P(I,4)=GA*(P(I,4)+BEP)
  660   CONTINUE
      ENDIF
 
C...Fill in position of decay vertex.
      DO 680 I=NSAV+1,N
        DO 670 J=1,4
          V(I,J)=VDCY(J)
  670   CONTINUE
        V(I,5)=0D0
  680 CONTINUE
 
C...Set up for parton shower evolution from jets.
      IF(MSTJ(23).GE.1.AND.MMAT.EQ.4.AND.K(NSAV+1,2).EQ.21) THEN
        K(NSAV+1,1)=3
        K(NSAV+2,1)=3
        K(NSAV+3,1)=3
        K(NSAV+1,4)=MSTU(5)*(NSAV+2)
        K(NSAV+1,5)=MSTU(5)*(NSAV+3)
        K(NSAV+2,4)=MSTU(5)*(NSAV+3)
        K(NSAV+2,5)=MSTU(5)*(NSAV+1)
        K(NSAV+3,4)=MSTU(5)*(NSAV+1)
        K(NSAV+3,5)=MSTU(5)*(NSAV+2)
        MSTJ(92)=-(NSAV+1)
      ELSEIF(MSTJ(23).GE.1.AND.MMAT.EQ.4) THEN
        K(NSAV+2,1)=3
        K(NSAV+3,1)=3
        K(NSAV+2,4)=MSTU(5)*(NSAV+3)
        K(NSAV+2,5)=MSTU(5)*(NSAV+3)
        K(NSAV+3,4)=MSTU(5)*(NSAV+2)
        K(NSAV+3,5)=MSTU(5)*(NSAV+2)
        MSTJ(92)=NSAV+2
      ELSEIF(MSTJ(23).GE.1.AND.(MMAT.EQ.32.OR.MMAT.EQ.44).AND.
     &  IABS(K(NSAV+1,2)).LE.10.AND.IABS(K(NSAV+2,2)).LE.10) THEN
        K(NSAV+1,1)=3
        K(NSAV+2,1)=3
        K(NSAV+1,4)=MSTU(5)*(NSAV+2)
        K(NSAV+1,5)=MSTU(5)*(NSAV+2)
        K(NSAV+2,4)=MSTU(5)*(NSAV+1)
        K(NSAV+2,5)=MSTU(5)*(NSAV+1)
        MSTJ(92)=NSAV+1
      ELSEIF(MSTJ(23).GE.1.AND.(MMAT.EQ.32.OR.MMAT.EQ.44).AND.
     &  IABS(K(NSAV+1,2)).LE.20.AND.IABS(K(NSAV+2,2)).LE.20) THEN
        MSTJ(92)=NSAV+1
      ELSEIF(MSTJ(23).GE.1.AND.MMAT.EQ.33.AND.IABS(K(NSAV+2,2)).EQ.21)
     &  THEN
        K(NSAV+1,1)=3
        K(NSAV+2,1)=3
        K(NSAV+3,1)=3
        KCP=PYCOMP(K(NSAV+1,2))
        KQP=KCHG(KCP,2)*ISIGN(1,K(NSAV+1,2))
        JCON=4
        IF(KQP.LT.0) JCON=5
        K(NSAV+1,JCON)=MSTU(5)*(NSAV+2)
        K(NSAV+2,9-JCON)=MSTU(5)*(NSAV+1)
        K(NSAV+2,JCON)=MSTU(5)*(NSAV+3)
        K(NSAV+3,9-JCON)=MSTU(5)*(NSAV+2)
        MSTJ(92)=NSAV+1
      ELSEIF(MSTJ(23).GE.1.AND.MMAT.EQ.33) THEN
        K(NSAV+1,1)=3
        K(NSAV+3,1)=3
        K(NSAV+1,4)=MSTU(5)*(NSAV+3)
        K(NSAV+1,5)=MSTU(5)*(NSAV+3)
        K(NSAV+3,4)=MSTU(5)*(NSAV+1)
        K(NSAV+3,5)=MSTU(5)*(NSAV+1)
        MSTJ(92)=NSAV+1
      ENDIF
 
C...Mark decayed particle; special option for B-Bbar mixing.
      IF(K(IP,1).EQ.5) K(IP,1)=15
      IF(K(IP,1).LE.10) K(IP,1)=11
      IF(MMIX.EQ.1.AND.MSTJ(26).EQ.2.AND.K(IP,1).EQ.11) K(IP,1)=12
      K(IP,4)=NSAV+1
      K(IP,5)=N
 
      RETURN
      END

C...PYDATA
C...Default values for switches and parameters,
C...and particle, decay and process data.
 
      BLOCK DATA PYDATA
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(4000,2),BRAT(4000),KFDP(4000,5)
      COMMON/PYDAT4/CHAF(500,2)
      CHARACTER CHAF*16
      COMMON/PYDATR/MRPY(6),RRPY(100)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      COMMON/PYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)
      COMMON/PYINT4/MWID(500),WIDS(500,5)
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
      COMMON/PYINT6/PROC(0:500)
      CHARACTER PROC*28
      COMMON/PYINT7/SIGT(0:6,0:6,0:5)
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      COMMON/PYSSMT/ZMIX(4,4),UMIX(2,2),VMIX(2,2),SMZ(4),SMW(2),
     &SFMIX(16,4)
      COMMON/PYBINS/IHIST(4),INDX(1000),BIN(20000)
      SAVE /PYDAT1/,/PYDAT2/,/PYDAT3/,/PYDAT4/,/PYDATR/,/PYSUBS/,
     &/PYPARS/,/PYINT1/,/PYINT2/,/PYINT3/,/PYINT4/,/PYINT5/,
     &/PYINT6/,/PYINT7/,/PYMSSM/,/PYSSMT/,/PYBINS/
 
C...PYDAT1, containing status codes and most parameters.
      DATA MSTU/
     &   0,    0,    0, 4000,10000,  500, 4000,    0,    0,    2,
     1   6,    1,    1,    0,    1,    1,    0,    0,    0,    0,
     2   2,   10,    0,    0,    1,   10,    0,    0,    0,    0,
     3   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     4   2,    2,    1,    4,    2,    1,    1,    0,    0,    0,
     5  25,   24,    0,    1,    0,    0,    0,    0,    0,    0,
     6   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     7  30*0,
     1   1,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     2   1,    5,    3,    5,    0,    0,    0,    0,    0,    0,
     &  80*0/
      DATA PARU/
     &  3.141592653589793D0, 6.283185307179586D0,
     &  0.197327D0, 5.06773D0, 0.389380D0, 2.56819D0,  4*0D0,
     1  0.001D0, 0.09D0, 0.01D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0,
     2  0D0,   0D0,   0D0,   0D0,  0D0,  0D0,  0D0,  0D0,  0D0,  0D0,
     3  0D0,   0D0,   0D0,   0D0,  0D0,  0D0,  0D0,  0D0,  0D0,  0D0,
     4  2.0D0,  1.0D0, 0.25D0,  2.5D0, 0.05D0,
     4  0D0,   0D0, 0.0001D0, 0D0,   0D0,
     5  2.5D0,1.5D0,7.0D0,1.0D0,0.5D0,2.0D0,3.2D0, 0D0, 0D0, 0D0,
     6  40*0D0,
     &  0.00729735D0, 0.232D0, 0.007764D0, 1.0D0, 1.16639D-5,
     &  0D0, 0D0, 0D0, 0D0,  0D0,
     1  0.20D0, 0.25D0, 1.0D0, 4.0D0, 10D0, 0D0, 0D0,  0D0, 0D0, 0D0,
     2 -0.693D0, -1.0D0, 0.387D0, 1.0D0, -0.08D0,
     2 -1.0D0,  1.0D0,  1.0D0,  1.0D0,  0D0,
     3  1.0D0,-1.0D0, 1.0D0,-1.0D0, 1.0D0,  0D0,  0D0, 0D0, 0D0, 0D0,
     4  5.0D0, 1.0D0, 1.0D0,  0D0, 1.0D0, 1.0D0,  0D0, 0D0, 0D0, 0D0,
     5  1.0D0, 0D0, 0D0, 0D0, 1000D0, 1.0D0, 1.0D0, 1.0D0, 1.0D0,0D0,
     6  1.0D0, 1.0D0, 1.0D0, 1.0D0, 1.0D0,  0D0,  0D0, 0D0, 0D0, 0D0,
     7  1.0D0, 1.0D0, 1.0D0, 1.0D0, 1.0D0, 1.0D0, 1.0D0, 0D0,0D0,0D0,
     8  1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0, 1.0D0, 1.0D0, 0D0,0D0,0D0,
     9  0D0,  0D0,  0D0,  0D0, 1.0D0,  0D0,  0D0, 0D0, 0D0, 0D0/
      DATA MSTJ/
     &  1,    3,    0,    0,    0,    0,    0,    0,    0,    0,
     1  4,    2,    0,    1,    0,    0,    0,    0,    0,    0,
     2  2,    1,    1,    2,    1,    2,    2,    0,    0,    0,
     3  0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     4  2,    2,    4,    2,    5,    3,    3,    0,    0,    3,
     5  0,    3,    0,    0,    0,    0,    0,    0,    0,    0,
     6  40*0,
     &  5,    2,    7,    5,    1,    1,    0,    2,    0,    2,
     1  0,    0,    0,    0,    1,    1,    0,    0,    0,    0,
     2  80*0/
      DATA PARJ/
     &  0.10D0, 0.30D0, 0.40D0, 0.05D0, 0.50D0,
     &  0.50D0, 0.50D0,   0.6D0,   1.2D0,   0.6D0,
     1  0.50D0,0.60D0,0.75D0, 0D0, 0D0, 0D0, 0D0, 1.0D0, 1.0D0, 0D0,
     2  0.36D0, 1.0D0,0.01D0, 2.0D0,1.0D0,0.4D0, 0D0, 0D0, 0D0, 0D0,
     3  0.10D0, 1.0D0, 0.8D0, 1.5D0,0D0,2.0D0,0.2D0,2.5D0,0.6D0,0D0,
     4  0.3D0, 0.58D0, 0.5D0, 0.9D0,0.5D0,1.0D0,1.0D0,1.0D0,0D0,0D0,
     5  0.77D0, 0.77D0, 0.77D0, -0.05D0, -0.005D0,
     5 -0.00001D0, -0.00001D0, -0.00001D0, 1.0D0, 0D0,
     6  4.5D0, 0.7D0, 0D0,0.003D0, 0.5D0, 0.5D0, 0D0, 0D0, 0D0, 0D0,
     7  10D0, 1000D0, 100D0, 1000D0, 0D0, 0.7D0,10D0, 0D0, 0D0, 0D0,
     8  0.29D0, 1.0D0, 1.0D0,  0D0,  10D0, 10D0, 0D0, 0D0, 0D0, 0D0,
     9  0.02D0, 1.0D0, 0.2D0,  0D0,  0D0,  0D0,  0D0, 0D0, 0D0, 0D0,
     &  0D0,  0D0,  0D0,  0D0,   0D0,   0D0,  0D0,  0D0,  0D0,  0D0,
     1  0D0,  0D0,  0D0,  0D0,   0D0,   0D0,  0D0,  0D0,  0D0,  0D0,
     2  1.0D0, 0.25D0,91.187D0,2.489D0, 0.01D0,
     2  2.0D0,  1.0D0, 0.25D0,0.002D0,   0D0,
     3  0D0, 0D0, 0D0, 0D0, 0.01D0, 0.99D0, 0D0, 0D0,  0.2D0,   0D0,
     4  60*0D0/
 
C...PYDAT2, with particle data and flavour treatment parameters.
      DATA (KCHG(I,1),I=   1, 500)/-1,2,-1,2,-1,2,-1,2,2*0,-3,0,-3,0,
     &-3,0,-3,6*0,3,9*0,3,2*0,3,0,-1,12*0,3,2*0,3,28*0,2,-1,20*0,4*3,
     &8*0,3*3,4*0,3*3,3*0,3*3,7*0,3*3,3*0,3*3,3*0,-2,-3,2*1,3*0,4,3*3,
     &6,2*-2,2*-3,0,2*1,2*0,2*3,-2,2*-3,2*0,-3,2*1,2*0,3,0,2*4,2*3,2*6,
     &3,2*1,2*0,2*3,2*0,4,2*3,2*6,2*3,6,2*-2,2*-3,0,-3,0,2*1,2*0,2*3,0,
     &3,2*-2,2*-3,2*0,2*-3,0,2*1,2*0,2*3,2*0,2*3,-2,2*-3,2*0,2*-3,2*0,
     &-3,2*0,2*3,4*0,2*3,2*0,2*3,2*0,2*3,4*0,2*3,2*0,2*3,3*0,3,2*0,3,0,
     &3,0,3,2*0,3,0,3,3*0,-1,2,-1,2,-1,2,-3,0,-3,0,-3,4*0,3,2*0,3,0,-1,
     &2,-1,2,-1,2,-3,0,-3,0,-3,0,-1,2,-3,164*0/
      DATA (KCHG(I,2),I=   1, 500)/8*1,12*0,2,16*0,2,1,113*0,-1,0,2*-1,
     &3*0,-1,4*0,2*-1,3*0,2*-1,4*0,-1,5*0,2*-1,4*0,2*-1,5*0,2*-1,6*0,
     &-1,7*0,2*-1,5*0,2*-1,6*0,2*-1,7*0,2*-1,8*0,-1,56*0,6*1,6*0,2,7*0,
     &6*1,6*0,2*1,165*0/
      DATA (KCHG(I,3),I=   1, 500)/8*1,2*0,8*1,5*0,1,9*0,1,2*0,1,0,2*1,
     &11*0,1,2*0,1,26*0,1,0,2*1,20*0,4*1,5*0,6*1,4*0,9*1,4*0,12*1,3*0,
     &102*1,2*0,2*1,2*0,4*1,2*0,6*1,2*0,8*1,3*0,1,0,2*1,0,3*1,0,4*1,
     &3*0,12*1,3*0,1,2*0,1,0,16*1,163*0/
      DATA (KCHG(I,4),I=   1, 293)/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
     &16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
     &37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,
     &58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,
     &79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,
     &100,110,111,113,115,130,210,211,213,215,220,221,223,225,310,311,
     &313,315,321,323,325,330,331,333,335,411,413,415,421,423,425,431,
     &433,435,440,441,443,445,511,513,515,521,523,525,531,533,535,541,
     &543,545,551,553,555,1103,1114,2101,2103,2110,2112,2114,2203,2210,
     &2212,2214,2224,3101,3103,3112,3114,3122,3201,3203,3212,3214,3222,
     &3224,3303,3312,3314,3322,3324,3334,4101,4103,4112,4114,4122,4132,
     &4201,4203,4212,4214,4222,4224,4232,4301,4303,4312,4314,4322,4324,
     &4332,4334,4403,4412,4414,4422,4424,4432,4434,4444,5101,5103,5112,
     &5114,5122,5132,5142,5201,5203,5212,5214,5222,5224,5232,5242,5301,
     &5303,5312,5314,5322,5324,5332,5334,5342,5401,5403,5412,5414,5422,
     &5424,5432,5434,5442,5444,5503,5512,5514,5522,5524,5532,5534,5542,
     &5544,5554,10111,10113,10211,10213,10221,10223,10311,10313,10321,
     &10323,10331,10333,10411,10413,10421,10423,10431,10433,10441,
     &10443,10511,10513,10521,10523,10531,10533,10541,10543,10551,
     &10553,20113,20213,20223,20313,20323,20333,20413,20423,20433/
      DATA (KCHG(I,4),I= 294, 500)/20443,20513,20523,20533,20543,20553,
     &100443,100553,1000001,1000002,1000003,1000004,1000005,1000006,
     &1000011,1000012,1000013,1000014,1000015,1000016,1000021,1000022,
     &1000023,1000024,1000025,1000035,1000037,1000039,2000001,2000002,
     &2000003,2000004,2000005,2000006,2000011,2000012,2000013,2000014,
     &2000015,2000016,4000001,4000002,4000011,4000012,163*0/
      DATA (PMAS(I,1),I=   1, 214)/0.0099D0,0.0056D0,0.199D0,1.35D0,
     &5D0,175D0,2*400D0,2*0D0,0.00051D0,0D0,0.10566D0,0D0,1.777D0,0D0,
     &400D0,5*0D0,91.187D0,80.33D0,80D0,6*0D0,500D0,900D0,500D0,
     &3*300D0,350D0,200D0,5000D0,10*0D0,3*100D0,3*200D0,26*0D0,1D0,2D0,
     &5D0,16*0D0,0.13498D0,0.7685D0,1.318D0,0.49767D0,0D0,0.13957D0,
     &0.7669D0,1.318D0,0D0,0.54745D0,0.78194D0,1.275D0,2*0.49767D0,
     &0.8961D0,1.432D0,0.4936D0,0.8916D0,1.425D0,0D0,0.95777D0,
     &1.0194D0,1.525D0,1.8693D0,2.01D0,2.46D0,1.8645D0,2.0067D0,2.46D0,
     &1.9685D0,2.1124D0,2.5735D0,0D0,2.9798D0,3.09688D0,3.5562D0,
     &5.2792D0,5.3248D0,5.83D0,5.2789D0,5.3248D0,5.83D0,5.3693D0,
     &5.4163D0,6.07D0,6.594D0,6.602D0,7.35D0,9.4D0,9.4603D0,9.9132D0,
     &0.77133D0,1.234D0,0.57933D0,0.77133D0,0D0,0.93957D0,1.233D0,
     &0.77133D0,0D0,0.93827D0,1.232D0,1.231D0,0.80473D0,0.92953D0,
     &1.19744D0,1.3872D0,1.11568D0,0.80473D0,0.92953D0,1.19255D0,
     &1.3837D0,1.18937D0,1.3828D0,1.09361D0,1.3213D0,1.535D0,1.3149D0,
     &1.5318D0,1.67245D0,1.96908D0,2.00808D0,2.4521D0,2.5D0,2.2849D0,
     &2.4703D0,1.96908D0,2.00808D0,2.4535D0,2.5D0,2.4529D0,2.5D0,
     &2.4656D0,2.15432D0,2.17967D0,2.55D0,2.63D0,2.55D0,2.63D0,2.704D0,
     &2.8D0,3.27531D0,3.59798D0,3.65648D0,3.59798D0,3.65648D0,
     &3.78663D0,3.82466D0,4.91594D0,5.38897D0,5.40145D0,5.8D0,5.81D0/
      DATA (PMAS(I,1),I= 215, 500)/5.641D0,5.84D0,7.00575D0,5.38897D0,
     &5.40145D0,5.8D0,5.81D0,5.8D0,5.81D0,5.84D0,7.00575D0,5.56725D0,
     &5.57536D0,5.96D0,5.97D0,5.96D0,5.97D0,6.12D0,6.13D0,7.19099D0,
     &6.67143D0,6.67397D0,7.03724D0,7.0485D0,7.03724D0,7.0485D0,
     &7.21101D0,7.219D0,8.30945D0,8.31325D0,10.07354D0,10.42272D0,
     &10.44144D0,10.42272D0,10.44144D0,10.60209D0,10.61426D0,
     &11.70767D0,11.71147D0,15.11061D0,0.9835D0,1.231D0,0.9835D0,
     &1.231D0,1D0,1.17D0,1.429D0,1.29D0,1.429D0,1.29D0,2*1.4D0,2.272D0,
     &2.424D0,2.272D0,2.424D0,2.5D0,2.536D0,3.4151D0,3.46D0,5.68D0,
     &5.73D0,5.68D0,5.73D0,5.92D0,5.97D0,7.25D0,7.3D0,9.8598D0,9.875D0,
     &2*1.23D0,1.282D0,2*1.402D0,1.427D0,2*2.372D0,2.56D0,3.5106D0,
     &2*5.78D0,6.02D0,7.3D0,9.8919D0,3.686D0,10.0233D0,32*500D0,
     &4*400D0,163*0D0/
      DATA (PMAS(I,2),I=   1, 500)/5*0D0,1.4D0,16*0D0,2.47833D0,
     &2.069D0,0.00295D0,6*0D0,14.67788D0,0D0,16.79392D0,8.45231D0,
     &4.93534D0,5.80468D0,19.1898D0,0.39162D0,417.35283D0,62*0D0,
     &0.151D0,0.107D0,3*0D0,0.149D0,0.107D0,2*0D0,0.00843D0,0.185D0,
     &2*0D0,0.0505D0,0.109D0,0D0,0.0498D0,0.098D0,0D0,0.0002D0,
     &0.00443D0,0.076D0,2*0D0,0.023D0,2*0D0,0.023D0,2*0D0,0.015D0,0D0,
     &0.0013D0,0D0,0.002D0,2*0D0,0.02D0,2*0D0,0.02D0,2*0D0,0.02D0,
     &2*0D0,0.02D0,4*0D0,0.12D0,4*0D0,0.12D0,3*0D0,2*0.12D0,3*0D0,
     &0.0394D0,4*0D0,0.036D0,0D0,0.0358D0,2*0D0,0.0099D0,0D0,0.0091D0,
     &74*0D0,0.06D0,0.142D0,0.06D0,0.142D0,0D0,0.36D0,0.287D0,0.09D0,
     &0.287D0,0.09D0,0.25D0,0.08D0,0.05D0,0.02D0,0.05D0,0.02D0,0.05D0,
     &0D0,0.014D0,0.01D0,8*0.05D0,0D0,0.01D0,2*0.4D0,0.025D0,2*0.174D0,
     &0.053D0,3*0.05D0,0.0009D0,4*0.05D0,3*0D0,19*1D0,0D0,7*1D0,0D0,
     &1D0,0D0,1D0,0D0,2.60511D0,2.60839D0,0.42904D0,0.41921D0,163*0D0/
      DATA (PMAS(I,3),I=   1, 500)/5*0D0,14D0,16*0D0,24.78326D0,
     &20.69D0,0.02954D0,6*0D0,146.77876D0,0D0,167.93924D0,84.52308D0,
     &49.35344D0,58.04675D0,191.89803D0,3.91624D0,4173.5283D0,62*0D0,
     &0.4D0,0.25D0,3*0D0,0.4D0,0.25D0,2*0D0,0.1D0,0.17D0,2*0D0,0.2D0,
     &0.12D0,0D0,0.2D0,0.12D0,0D0,0.002D0,0.015D0,0.2D0,2*0D0,0.12D0,
     &2*0D0,0.12D0,2*0D0,0.05D0,0D0,0.005D0,0D0,0.01D0,2*0D0,0.05D0,
     &2*0D0,0.05D0,2*0D0,0.05D0,2*0D0,0.05D0,4*0D0,0.14D0,4*0D0,0.14D0,
     &3*0D0,2*0.14D0,3*0D0,0.04D0,4*0D0,0.035D0,0D0,0.035D0,2*0D0,
     &0.05D0,0D0,0.05D0,74*0D0,0.05D0,0.25D0,0.05D0,0.25D0,0D0,0.2D0,
     &0.4D0,0.005D0,0.4D0,0.01D0,0.35D0,0.001D0,0.1D0,0.08D0,0.1D0,
     &0.08D0,0.1D0,0D0,0.05D0,0.02D0,6*0.1D0,0.05D0,0.1D0,0D0,0.02D0,
     &2*0.3D0,0.05D0,2*0.3D0,0.02D0,2*0.1D0,0.03D0,0.001D0,4*0.1D0,
     &3*0D0,19*10D0,0.00001D0,7*10D0,0.00001D0,10D0,0.00001D0,10D0,
     &0.00001D0,26.05109D0,26.08388D0,4.29043D0,4.19206D0,163*0D0/
      DATA (PMAS(I,4),I=   1, 500)/12*0D0,658654D0,0D0,0.0872D0,68*0D0,
     &0.1D0,0.387D0,16*0D0,0.00003D0,2*0D0,15500D0,0D0,7804.5D0,6*0D0,
     &26.762D0,3*0D0,3709D0,6*0D0,0.317D0,2*0D0,0.1244D0,2*0D0,0.14D0,
     &6*0D0,0.468D0,2*0D0,0.462D0,2*0D0,0.483D0,2*0D0,0.15D0,19*0D0,
     &44.34D0,0D0,78.88D0,4*0D0,23.96D0,2*0D0,49.1D0,0D0,87.1D0,0D0,
     &24.6D0,4*0D0,0.0618D0,0.029D0,6*0D0,0.106D0,6*0D0,0.019D0,2*0D0,
     &7*0.1D0,4*0D0,0.342D0,2*0.387D0,6*0D0,2*0.387D0,6*0D0,0.387D0,
     &0D0,0.387D0,2*0D0,8*0.387D0,0D0,9*0.387D0,83*0D0,163*0D0/
      DATA PARF/
     &  0.5D0,0.25D0, 0.5D0,0.25D0, 1D0, 0.5D0,  0D0,  0D0,  0D0, 0D0,
     1  0.5D0,  0D0, 0.5D0,  0D0,  1D0,  1D0,  0D0,  0D0,  0D0, 0D0,
     2  0.5D0,  0D0, 0.5D0,  0D0,  1D0,  1D0,  0D0,  0D0,  0D0, 0D0,
     3  0.5D0,  0D0, 0.5D0,  0D0,  1D0,  1D0,  0D0,  0D0,  0D0, 0D0,
     4  0.5D0,  0D0, 0.5D0,  0D0,  1D0,  1D0,  0D0,  0D0,  0D0, 0D0,
     5  0.5D0,  0D0, 0.5D0,  0D0,  1D0,  1D0,  0D0,  0D0,  0D0, 0D0,
     6  0.75D0, 0.5D0, 0D0,0.1667D0,0.0833D0,0.1667D0,0D0,0D0,0D0, 0D0,
     7  0D0,  0D0,  1D0,0.3333D0,0.6667D0,0.3333D0,0D0,0D0,0D0, 0D0,
     8  0D0,  0D0,  0D0,  0D0,  0D0,  0D0,  0D0,  0D0,  0D0, 0D0,
     9  0D0,  0D0,  0D0,  0D0,  0D0,  0D0,  0D0,  0D0,  0D0, 0D0,
     & 0.325D0,0.325D0,0.5D0,1.6D0, 5.0D0,  0D0,  0D0,  0D0,  0D0, 0D0,
     1 0D0,0.11D0,0.16D0,0.048D0,0.50D0,0.45D0,0.55D0,0.60D0,0D0,0D0,
     2 0.2D0, 0.1D0,  0D0,  0D0,  0D0,  0D0,  0D0,  0D0,  0D0, 0D0,
     3 60*0D0,
     4 0.2D0,  0.5D0,  8*0D0, 
     5 1800*0D0/
      DATA ((VCKM(I,J),J=1,4),I=1,4)/
     &  0.95113D0,  0.04884D0,  0.00003D0,  0.00000D0,
     &  0.04884D0,  0.94940D0,  0.00176D0,  0.00000D0,
     &  0.00003D0,  0.00176D0,  0.99821D0,  0.00000D0,
     &  0.00000D0,  0.00000D0,  0.00000D0,  1.00000D0/
 
C...PYDAT3, with particle decay parameters and data.
      DATA (MDCY(I,1),I=   1, 500)/5*0,3*1,6*0,1,0,1,5*0,3*1,6*0,1,0,
     &7*1,10*0,2*1,0,3*1,26*0,3*1,16*0,3*1,3*0,2*1,0,7*1,0,2*1,0,12*1,
     &0,18*1,0,1,4*0,1,3*0,2*1,2*0,3*1,2*0,4*1,0,5*1,2*0,4*1,2*0,5*1,
     &2*0,6*1,0,7*1,2*0,5*1,2*0,6*1,2*0,7*1,2*0,8*1,0,75*1,0,7*1,0,1,0,
     &1,0,4*1,163*0/
      DATA (MDCY(I,2),I=   1, 500)/1,9,17,25,33,41,54,64,2*0,74,78,80,
     &85,87,141,143,148,2*0,151,160,172,188,208,6*0,287,0,309,332,414,
     &494,521,524,525,10*0,534,539,0,544,564,588,26*0,606,607,611,16*0,
     &620,622,627,636,0,645,647,649,0,656,664,670,679,681,683,686,696,
     &702,705,0,716,722,733,739,802,805,813,874,876,884,917,919,0,923,
     &924,927,929,965,966,974,1010,1011,1019,1058,1059,1063,1094,1095,
     &1099,1100,1109,0,1111,4*0,1112,3*0,1115,1118,2*0,1119,1121,1124,
     &2*0,1128,1129,1132,1135,0,1138,1143,1145,1148,1150,2*0,1154,1155,
     &1156,1232,2*0,1236,1237,1238,1239,1240,2*0,1244,1245,1247,1248,
     &1250,1254,0,1255,1259,1263,1267,1271,1275,1279,2*0,1283,1284,
     &1285,1302,1311,2*0,1320,1321,1322,1323,1324,1333,2*0,1342,1343,
     &1344,1345,1346,1355,1356,2*0,1365,1374,1383,1392,1401,1410,1419,
     &1428,0,1437,1446,1455,1464,1473,1482,1491,1500,1509,1518,1519,
     &1520,1521,1522,1527,1530,1532,1537,1539,1544,1551,1555,1557,1559,
     &1561,1563,1565,1567,1569,1570,1572,1574,1576,1578,1580,1582,1584,
     &1586,1588,1589,1591,1593,1607,1609,1611,1615,1617,1619,1621,1623,
     &1625,1627,1629,1631,1633,1644,1658,1670,1682,1694,1706,1718,1731,
     &1742,1753,1764,1775,1786,1797,1858,1863,1965,2021,2139,2273,0,
     &2344,2360,2376,2392,2408,2424,2440,0,2455,0,2470,0,2485,2489,
     &2493,2496,163*0/
      DATA (MDCY(I,3),I=   1, 500)/5*8,13,2*10,2*0,4,2,5,2,54,2,5,3,
     &2*0,9,12,16,20,79,6*0,22,0,23,82,80,27,3,1,9,10*0,2*5,0,20,24,18,
     &26*0,1,4,9,16*0,2,5,2*9,0,2*2,7,0,8,6,9,2*2,3,10,6,3,11,0,6,11,6,
     &63,3,8,61,2,8,33,2,4,0,1,3,2,36,1,8,36,1,8,39,1,4,31,1,4,1,9,2,0,
     &1,4*0,3,3*0,3,1,2*0,2,3,4,2*0,1,3*3,0,5,2,3,2,4,2*0,2*1,76,4,2*0,
     &4*1,4,2*0,1,2,1,2,4,1,0,7*4,2*0,2*1,17,2*9,2*0,4*1,2*9,2*0,4*1,9,
     &1,9,2*0,8*9,0,9*9,4*1,5,3,2,5,2,5,7,4,7*2,1,9*2,1,2*2,14,2*2,4,
     &9*2,11,14,5*12,13,6*11,61,5,102,56,118,134,71,0,6*16,15,0,15,0,
     &15,0,2*4,3,2,163*0/
      DATA (MDME(I,1),I=   1,4000)/6*1,-1,7*1,-1,7*1,-1,7*1,-1,7*1,-1,
     &7*1,-1,1,-1,12*1,2*-1,8*1,2*-1,73*1,-1,2*1,-1,6*1,2*-1,7*1,2*-1,
     &3*1,-1,6*1,2*-1,6*1,2*-1,3*1,-1,3*1,-1,3*1,5*-1,3*1,-1,85*1,2*-1,
     &6*1,8*-1,3*1,-1,3*1,-1,3*1,5*-1,3*1,4*-1,197*1,2*-1,2*1,-1,20*1,
     &2*-1,6*1,2*-1,7*1,-1,3*1,-1,3*1,5*-1,3*1,-1,1,-1,6*1,2*-1,6*1,
     &2*-1,1892*1,1503*0/
      DATA (MDME(I,2),I=   1,4000)/43*102,4*0,102,0,4*53,3*102,4*0,102,
     &2*0,3*102,4*0,102,2*0,6*102,42,6*102,2*42,2*0,8*41,2*0,36*41,
     &8*102,0,102,0,102,2*0,21*102,8*32,8*0,16*32,21*0,62*53,8*32,14*0,
     &16*32,27*0,62*53,18*0,62*53,9*0,18*53,3*32,0,6*32,3*0,2*32,3*0,
     &2*32,7*0,8*32,12*0,16*32,6*0,8*32,8*0,12,2*42,2*11,9*42,0,2,3,
     &15*0,4*42,5*0,3,12*0,2,3*0,1,0,3,16*0,2*3,15*0,2*42,2*3,18*0,2*3,
     &3*0,1,11*0,22*42,41*0,2*3,9*0,16*42,45*0,3,10*0,10*42,20*0,2*13,
     &6*0,12,2*0,12,0,12,14*42,16*0,48,3*13,2*42,9*0,14*42,16*0,48,
     &3*13,2*42,9*0,14*42,19*0,48,3*13,2*42,6*0,2*11,28*42,5*0,32,3*0,
     &4*32,2*4,0,32,45*0,14*42,52*0,10*13,2*42,2*11,4*0,2*42,2*11,6*0,
     &2*42,2*11,0,2*42,2*11,2*42,2*11,2*42,2*11,2*42,2*11,2*42,2*11,
     &2*42,2*11,2*42,2*11,2*0,3*42,8*0,48,3*13,20*42,4*0,18*42,4*0,
     &9*42,0,162*42,50*0,2*12,17*0,2*32,33*0,12,9*0,32,2*0,12,11*0,
     &4*32,2*4,5*0,828*53,1515*0/
      DATA (BRAT(I)  ,I=   1, 418)/43*0D0,0.00003D0,0.00177D0,0.9982D0,
     &33*0D0,1D0,6*0D0,0.1783D0,0.1735D0,0.1131D0,0.2494D0,0.003D0,
     &0.09D0,0.0027D0,0.01D0,0.0014D0,0.0012D0,2*0.00025D0,0.0071D0,
     &0.012D0,0.0004D0,0.00075D0,0.00006D0,2*0.00078D0,0.0034D0,0.08D0,
     &0.011D0,0.0191D0,0.00006D0,0.005D0,0.0133D0,0.0067D0,0.0005D0,
     &0.0035D0,0.0006D0,0.0015D0,0.00021D0,0.0002D0,0.00075D0,0.0001D0,
     &0.0002D0,0.0011D0,3*0.0002D0,0.00022D0,0.0004D0,0.0001D0,
     &2*0.00205D0,2*0.00069D0,0.00025D0,0.00051D0,0.00025D0,35*0D0,
     &0.15403D0,0.11945D0,0.15402D0,0.11931D0,0.15215D0,3*0D0,
     &0.03357D0,0.0668D0,0.03357D0,0.0668D0,0.0335D0,0.0668D0,2*0D0,
     &0.32139D0,0.0165D0,2*0D0,0.0165D0,0.32067D0,2*0D0,0.00001D0,
     &0.00059D0,6*0D0,2*0.10814D0,0.10806D0,3*0D0,0.00031D0,0.04438D0,
     &0.88031D0,4*0D0,0.0002D0,0.05531D0,0D0,0.01838D0,0.00071D0,0D0,
     &0.00009D0,0.00032D0,62*0D0,0.14449D0,0.11223D0,0.14449D0,
     &0.11223D0,0.14443D0,0.05782D0,2*0D0,0.03172D0,0.06305D0,
     &0.03172D0,0.06305D0,0.03172D0,0.06305D0,8*0D0,0.24928D0,0.0128D0,
     &0.00001D0,0D0,0.0128D0,0.24882D0,0.00039D0,0D0,0.00001D0,
     &0.00046D0,0.22153D0,5*0D0,2*0.08464D0,0.08463D0,7*0D0,0.00005D0,
     &0.00097D0,5*0D0,0.00007D0,0D0,0.00049D0,0.00001D0,0.00006D0,
     &0.30591D0,0.68863D0,0D0,0.0038D0,66*0D0,0.00008D0,0.00167D0/
      DATA (BRAT(I)  ,I= 419, 722)/5*0D0,0.00013D0,0D0,0.00294D0,
     &0.00001D0,3*0D0,0.99517D0,63*0D0,0.00002D0,0.07231D0,2*0D0,
     &0.00001D0,0.00269D0,0D0,0.92497D0,18*0D0,0.0024D0,0.99483D0,
     &0.00278D0,1D0,3*0.21511D0,0.21478D0,2*0D0,2*0.06995D0,2*0D0,1D0,
     &3*0D0,0.95D0,0.05D0,3*0D0,4*0.25D0,16*0D0,4*0.25D0,20*0D0,1D0,
     &17*0D0,1D0,2*0.08D0,0.76D0,0.08D0,2*0.105D0,0.04D0,0.5D0,0.08D0,
     &0.14D0,0.01D0,0.015D0,0.005D0,0.988D0,0.012D0,0.998739D0,
     &0.00079D0,0.00038D0,0.000046D0,0.000045D0,2*0.34725D0,0.144D0,
     &0.104D0,0.0245D0,2*0.01225D0,0.0028D0,0.0057D0,0.2112D0,0.1256D0,
     &2*0.1939D0,2*0.1359D0,0.002D0,0.001D0,0.0006D0,0.999877D0,
     &0.000123D0,0.99955D0,0.00045D0,2*0.34725D0,0.144D0,0.104D0,
     &0.049D0,0.0028D0,0.0057D0,0.3923D0,0.321D0,0.2317D0,0.0478D0,
     &0.0049D0,0.0013D0,0.0003D0,0.0007D0,0.89D0,0.08693D0,0.0221D0,
     &0.00083D0,2*0.00007D0,0.564D0,0.282D0,0.072D0,0.028D0,0.023D0,
     &2*0.0115D0,0.005D0,0.003D0,0.6861D0,0.3139D0,2*0.5D0,0.665D0,
     &0.333D0,0.002D0,0.333D0,0.166D0,0.168D0,0.084D0,0.087D0,0.043D0,
     &0.059D0,2*0.029D0,0.002D0,0.6352D0,0.2116D0,0.0559D0,0.0173D0,
     &0.0482D0,0.0318D0,0.666D0,0.333D0,0.001D0,0.332D0,0.166D0,
     &0.168D0,0.084D0,0.086D0,0.043D0,0.059D0,2*0.029D0,2*0.002D0,
     &0.437D0,0.208D0,0.302D0,0.0302D0,0.0212D0,0.0016D0,0.48947D0/
      DATA (BRAT(I)  ,I= 723, 897)/0.34D0,3*0.043D0,0.027D0,0.0126D0,
     &0.0013D0,0.0003D0,0.00025D0,0.00008D0,0.444D0,2*0.222D0,0.104D0,
     &2*0.004D0,0.07D0,0.065D0,2*0.005D0,2*0.011D0,5*0.001D0,0.07D0,
     &0.065D0,2*0.005D0,2*0.011D0,5*0.001D0,0.026D0,0.019D0,0.066D0,
     &0.041D0,0.045D0,0.076D0,0.0073D0,2*0.0047D0,0.026D0,0.001D0,
     &0.0006D0,0.0066D0,0.005D0,2*0.003D0,2*0.0006D0,2*0.001D0,0.006D0,
     &0.005D0,0.012D0,0.0057D0,0.067D0,0.008D0,0.0022D0,0.027D0,
     &0.004D0,0.019D0,0.012D0,0.002D0,0.009D0,0.0218D0,0.001D0,0.022D0,
     &0.087D0,0.001D0,0.0019D0,0.0015D0,0.0028D0,0.683D0,0.306D0,
     &0.011D0,0.3D0,0.15D0,0.16D0,0.08D0,0.13D0,0.06D0,0.08D0,0.04D0,
     &0.034D0,0.027D0,2*0.002D0,2*0.004D0,2*0.002D0,0.034D0,0.027D0,
     &2*0.002D0,2*0.004D0,2*0.002D0,0.0365D0,0.045D0,0.073D0,0.062D0,
     &3*0.021D0,0.0061D0,0.015D0,0.025D0,0.0088D0,0.074D0,0.0109D0,
     &0.0041D0,0.002D0,0.0035D0,0.0011D0,0.001D0,0.0027D0,2*0.0016D0,
     &0.0018D0,0.011D0,0.0063D0,0.0052D0,0.018D0,0.016D0,0.0034D0,
     &0.0036D0,0.0009D0,0.0006D0,0.015D0,0.0923D0,0.018D0,0.022D0,
     &0.0077D0,0.009D0,0.0075D0,0.024D0,0.0085D0,0.067D0,0.0511D0,
     &0.017D0,0.0004D0,0.0028D0,0.619D0,0.381D0,0.3D0,0.15D0,0.16D0,
     &0.08D0,0.13D0,0.06D0,0.08D0,0.04D0,0.01D0,2*0.02D0,0.03D0,
     &2*0.005D0,2*0.02D0,0.03D0,2*0.005D0,0.015D0,0.037D0,0.028D0/
      DATA (BRAT(I)  ,I= 898,1063)/0.079D0,0.095D0,0.052D0,0.0078D0,
     &4*0.001D0,0.028D0,0.033D0,0.026D0,0.05D0,0.01D0,4*0.005D0,0.25D0,
     &0.0952D0,0.94D0,0.06D0,2*0.4D0,2*0.1D0,1D0,0.0602D0,0.0601D0,
     &0.8797D0,0.135D0,0.865D0,0.02D0,0.055D0,2*0.005D0,0.008D0,
     &0.012D0,0.02D0,0.055D0,2*0.005D0,0.008D0,0.012D0,0.01D0,0.03D0,
     &0.0035D0,0.011D0,0.0055D0,0.0042D0,0.009D0,0.018D0,0.015D0,
     &0.0185D0,0.0135D0,0.025D0,0.0004D0,0.0007D0,0.0008D0,0.0014D0,
     &0.0019D0,0.0025D0,0.4291D0,0.08D0,0.07D0,0.02D0,0.015D0,0.005D0,
     &1D0,0.3D0,0.15D0,0.16D0,0.08D0,0.13D0,0.06D0,0.08D0,0.04D0,
     &0.02D0,0.055D0,2*0.005D0,0.008D0,0.012D0,0.02D0,0.055D0,
     &2*0.005D0,0.008D0,0.012D0,0.01D0,0.03D0,0.0035D0,0.011D0,
     &0.0055D0,0.0042D0,0.009D0,0.018D0,0.015D0,0.0185D0,0.0135D0,
     &0.025D0,0.0004D0,0.0007D0,0.0008D0,0.0014D0,0.0019D0,0.0025D0,
     &0.4291D0,0.08D0,0.07D0,0.02D0,0.015D0,0.005D0,1D0,0.3D0,0.15D0,
     &0.16D0,0.08D0,0.13D0,0.06D0,0.08D0,0.04D0,0.02D0,0.055D0,
     &2*0.005D0,0.008D0,0.012D0,0.02D0,0.055D0,2*0.005D0,0.008D0,
     &0.012D0,0.01D0,0.03D0,0.0035D0,0.011D0,0.0055D0,0.0042D0,0.009D0,
     &0.018D0,0.015D0,0.0185D0,0.0135D0,0.025D0,2*0.0002D0,0.0007D0,
     &2*0.0004D0,0.0014D0,0.001D0,0.0009D0,0.0025D0,0.4291D0,0.08D0,
     &0.07D0,0.02D0,0.015D0,0.005D0,1D0,2*0.3D0,2*0.2D0,0.047D0/
      DATA (BRAT(I)  ,I=1064,1254)/0.122D0,0.006D0,0.012D0,0.035D0,
     &0.012D0,0.035D0,0.003D0,0.007D0,0.15D0,0.037D0,0.008D0,0.002D0,
     &0.05D0,0.015D0,0.003D0,0.001D0,0.014D0,0.042D0,0.014D0,0.042D0,
     &0.24D0,0.065D0,0.012D0,0.003D0,0.001D0,0.002D0,0.001D0,0.002D0,
     &0.014D0,0.003D0,1D0,2*0.3D0,2*0.2D0,1D0,0.0252D0,0.0248D0,
     &0.0267D0,0.015D0,0.045D0,0.015D0,0.045D0,0.7743D0,0.029D0,0.22D0,
     &0.78D0,1D0,0.331D0,0.663D0,0.006D0,0.663D0,0.331D0,0.006D0,1D0,
     &0.999D0,0.001D0,0.88D0,2*0.06D0,0.639D0,0.358D0,0.002D0,0.001D0,
     &1D0,0.88D0,2*0.06D0,0.516D0,0.483D0,0.001D0,0.88D0,2*0.06D0,
     &0.9988D0,0.0001D0,0.0006D0,0.0004D0,0.0001D0,0.667D0,0.333D0,
     &0.9954D0,0.0011D0,0.0035D0,0.333D0,0.667D0,0.676D0,0.234D0,
     &0.085D0,0.005D0,2*1D0,0.018D0,2*0.005D0,0.003D0,0.002D0,
     &2*0.006D0,0.018D0,2*0.005D0,0.003D0,0.002D0,2*0.006D0,0.0066D0,
     &0.025D0,0.016D0,0.0088D0,2*0.005D0,0.0058D0,0.005D0,0.0055D0,
     &4*0.004D0,2*0.002D0,2*0.004D0,0.003D0,0.002D0,2*0.003D0,
     &3*0.002D0,2*0.001D0,0.002D0,2*0.001D0,2*0.002D0,0.0013D0,
     &0.0018D0,5*0.001D0,4*0.003D0,2*0.005D0,2*0.002D0,2*0.001D0,
     &2*0.002D0,2*0.001D0,0.2432D0,0.057D0,2*0.035D0,0.15D0,2*0.075D0,
     &0.03D0,2*0.015D0,2*0.08D0,0.76D0,0.08D0,4*1D0,2*0.08D0,0.76D0,
     &0.08D0,1D0,2*0.5D0,1D0,2*0.5D0,2*0.08D0,0.76D0,0.08D0,1D0/
      DATA (BRAT(I)  ,I=1255,1447)/2*0.08D0,0.76D0,3*0.08D0,0.76D0,
     &3*0.08D0,0.76D0,3*0.08D0,0.76D0,3*0.08D0,0.76D0,3*0.08D0,0.76D0,
     &3*0.08D0,0.76D0,0.08D0,2*1D0,2*0.105D0,0.04D0,0.0077D0,0.02D0,
     &0.0235D0,0.0285D0,0.0435D0,0.0011D0,0.0022D0,0.0044D0,0.4291D0,
     &0.08D0,0.07D0,0.02D0,0.015D0,0.005D0,2*0.105D0,0.04D0,0.5D0,
     &0.08D0,0.14D0,0.01D0,0.015D0,0.005D0,2*0.105D0,0.04D0,0.5D0,
     &0.08D0,0.14D0,0.01D0,0.015D0,0.005D0,4*1D0,2*0.105D0,0.04D0,
     &0.5D0,0.08D0,0.14D0,0.01D0,0.015D0,0.005D0,2*0.105D0,0.04D0,
     &0.5D0,0.08D0,0.14D0,0.01D0,0.015D0,0.005D0,4*1D0,2*0.105D0,
     &0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,0.015D0,0.005D0,1D0,2*0.105D0,
     &0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,0.015D0,0.005D0,2*0.105D0,
     &0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,0.015D0,0.005D0,2*0.105D0,
     &0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,0.015D0,0.005D0,2*0.105D0,
     &0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,0.015D0,0.005D0,2*0.105D0,
     &0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,0.015D0,0.005D0,2*0.105D0,
     &0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,0.015D0,0.005D0,2*0.105D0,
     &0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,0.015D0,0.005D0,2*0.105D0,
     &0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,0.015D0,0.005D0,2*0.105D0,
     &0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,0.015D0,0.005D0,2*0.105D0,
     &0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,0.015D0,0.005D0,2*0.105D0/
      DATA (BRAT(I)  ,I=1448,1648)/0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,
     &0.015D0,0.005D0,2*0.105D0,0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,
     &0.015D0,0.005D0,2*0.105D0,0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,
     &0.015D0,0.005D0,2*0.105D0,0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,
     &0.015D0,0.005D0,2*0.105D0,0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,
     &0.015D0,0.005D0,2*0.105D0,0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,
     &0.015D0,0.005D0,2*0.105D0,0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,
     &0.015D0,0.005D0,2*0.105D0,0.04D0,0.5D0,0.08D0,0.14D0,0.01D0,
     &0.015D0,0.005D0,4*1D0,0.52D0,0.26D0,0.11D0,2*0.055D0,0.333D0,
     &0.334D0,0.333D0,0.667D0,0.333D0,0.28D0,0.14D0,0.313D0,0.157D0,
     &0.11D0,0.667D0,0.333D0,0.28D0,0.14D0,0.313D0,0.157D0,0.11D0,
     &0.36D0,0.18D0,0.03D0,2*0.015D0,2*0.2D0,4*0.25D0,0.667D0,0.333D0,
     &0.667D0,0.333D0,0.667D0,0.333D0,0.667D0,0.333D0,4*0.5D0,0.007D0,
     &0.993D0,1D0,0.667D0,0.333D0,0.667D0,0.333D0,0.667D0,0.333D0,
     &0.667D0,0.333D0,8*0.5D0,0.02D0,0.98D0,1D0,4*0.5D0,3*0.146D0,
     &3*0.05D0,0.15D0,2*0.05D0,4*0.024D0,0.066D0,0.667D0,0.333D0,
     &0.667D0,0.333D0,4*0.25D0,0.667D0,0.333D0,0.667D0,0.333D0,2*0.5D0,
     &0.273D0,0.727D0,0.667D0,0.333D0,0.667D0,0.333D0,4*0.5D0,0.35D0,
     &0.65D0,2*0.0083D0,0.1866D0,0.324D0,0.184D0,0.027D0,0.001D0,
     &0.093D0,0.087D0,0.078D0,0.0028D0,3*0.014D0,0.008D0,0.024D0/
      DATA (BRAT(I)  ,I=1649,4000)/0.008D0,0.024D0,0.425D0,0.02D0,
     &0.185D0,0.088D0,0.043D0,0.067D0,0.066D0,827*0D0,0.8516D0,
     &0.00539D0,0.04483D0,0.09819D0,0.85053D0,0.02152D0,0.02989D0,
     &0.09806D0,0.29439D0,0.10943D0,0.59618D0,0.38983D0,0.61017D0,
     &1503*0D0/
      DATA (KFDP(I,1),I=   1, 375)/21,22,23,4*-24,25,21,22,23,4*24,25,
     &21,22,23,4*-24,25,21,22,23,4*24,25,21,22,23,4*-24,25,21,22,23,
     &4*24,25,37,1000022,1000023,1000025,1000035,21,22,23,4*-24,25,
     &2*-37,21,22,23,4*24,25,2*37,22,23,-24,25,23,24,-12,22,23,-24,25,
     &23,24,-12,-14,48*16,22,23,-24,25,23,24,22,23,-24,25,-37,23,24,37,
     &1,2,3,4,5,6,7,8,21,1,2,3,4,5,6,7,8,11,13,15,17,1,2,3,4,5,6,7,8,
     &11,12,13,14,15,16,17,18,4*-1,4*-3,4*-5,4*-7,-11,-13,-15,-17,1,2,
     &3,4,5,6,7,8,11,13,15,17,21,2*22,23,24,1000022,2*1000023,
     &3*1000025,4*1000035,2*1000024,2*1000037,1000001,2000001,1000001,
     &-1000001,1000002,2000002,1000002,-1000002,1000003,2000003,
     &1000003,-1000003,1000004,2000004,1000004,-1000004,1000005,
     &2000005,1000005,-1000005,1000006,2000006,1000006,-1000006,
     &1000011,2000011,1000011,-1000011,1000012,2000012,1000012,
     &-1000012,1000013,2000013,1000013,-1000013,1000014,2000014,
     &1000014,-1000014,1000015,2000015,1000015,-1000015,1000016,
     &2000016,1000016,-1000016,1,2,3,4,5,6,7,8,11,12,13,14,15,16,17,18,
     &24,37,2*23,25,35,4*-1,4*-3,4*-5,4*-7,-11,-13,-15,-17,3*24,1,2,3,
     &4,5,6,7,8,11,13,15,17,21,2*22,23,24,23,25,36,1000022,2*1000023,
     &3*1000025,4*1000035,2*1000024,2*1000037,1000001,2000001,1000001,
     &-1000001,1000002,2000002,1000002,-1000002,1000003,2000003/
      DATA (KFDP(I,1),I= 376, 606)/1000003,-1000003,1000004,2000004,
     &1000004,-1000004,1000005,2000005,1000005,-1000005,1000006,
     &2000006,1000006,-1000006,1000011,2000011,1000011,-1000011,
     &1000012,2000012,1000012,-1000012,1000013,2000013,1000013,
     &-1000013,1000014,2000014,1000014,-1000014,1000015,2000015,
     &1000015,-1000015,1000016,2000016,1000016,-1000016,1,2,3,4,5,6,7,
     &8,11,13,15,17,21,2*22,23,24,23,1000022,2*1000023,3*1000025,
     &4*1000035,2*1000024,2*1000037,1000001,2000001,1000001,-1000001,
     &1000002,2000002,1000002,-1000002,1000003,2000003,1000003,
     &-1000003,1000004,2000004,1000004,-1000004,1000005,2000005,
     &1000005,-1000005,1000006,2000006,1000006,-1000006,1000011,
     &2000011,1000011,-1000011,1000012,2000012,1000012,-1000012,
     &1000013,2000013,1000013,-1000013,1000014,2000014,1000014,
     &-1000014,1000015,2000015,1000015,-1000015,1000016,2000016,
     &1000016,-1000016,-1,-3,-5,-7,-11,-13,-15,-17,24,2*1000022,
     &2*1000023,2*1000025,2*1000035,1000006,2000006,1000006,2000006,
     &-1000001,-1000003,-1000011,-1000013,-1000015,-2000015,5,6,21,2,1,
     &2,3,4,5,6,11,13,15,4,5,11,13,15,2*4,-11,-13,-15,2*24,2*52,1,2,3,
     &4,5,6,7,8,11,12,13,14,15,16,17,18,2*24,2*52,4*-1,4*-3,4*-5,4*-7,
     &-11,-13,-15,-17,22,23,1,2,3,4,5,6,7,8,11,12,13,14,15,16,17,18,82/
      DATA (KFDP(I,1),I= 607,1001)/-11,-13,2*2,-12,-14,-16,2*-2,2*-4,
     &-2,-4,2*22,211,111,221,13,11,213,-213,221,223,321,130,310,111,
     &331,111,211,-12,12,-14,14,211,111,22,-13,-11,2*211,213,113,221,
     &223,321,211,331,22,111,211,2*22,211,22,111,211,22,211,221,111,11,
     &211,111,2*211,321,130,310,221,111,211,111,130,310,321,2*311,321,
     &311,323,313,323,313,321,3*311,-13,3*211,12,14,311,2*321,311,321,
     &313,323,313,323,311,4*321,211,111,3*22,111,321,130,-213,113,213,
     &211,22,111,11,13,211,321,130,310,221,211,111,11*-11,11*-13,-311,
     &-313,-311,-313,-20313,2*-311,-313,-311,-313,2*111,2*221,2*331,
     &2*113,2*223,2*333,-311,-313,2*-321,211,-311,-321,333,-311,-313,
     &-321,211,2*-321,2*-311,-321,211,113,421,2*411,421,411,423,413,
     &423,413,421,411,8*-11,8*-13,-321,-323,-321,-323,-311,2*-313,-311,
     &-313,2*-311,-321,-10323,-321,-323,-321,-311,2*-313,211,111,333,
     &3*-321,-311,-313,-321,-313,310,333,211,2*-321,-311,-313,-311,211,
     &-321,3*-311,211,113,321,2*421,411,421,413,423,413,423,411,421,
     &-15,5*-11,5*-13,221,331,333,221,331,333,10221,211,213,211,213,
     &321,323,321,323,2212,221,331,333,221,2*2,2*431,421,411,423,413,
     &82,11,13,82,443,82,6*12,6*14,2*16,3*-411,3*-413,2*-411,2*-413,
     &2*441,2*443,2*20443,2*2,2*4,2,4,511,521,511,523,513,523,513,521,
     &511,6*12,6*14,2*16,3*-421,3*-423,2*-421,2*-423,2*441,2*443/
      DATA (KFDP(I,1),I=1002,1428)/2*20443,2*2,2*4,2,4,521,511,521,513,
     &523,513,523,511,521,6*12,6*14,2*16,3*-431,3*-433,2*-431,2*-433,
     &3*441,3*443,3*20443,2*2,2*4,2,4,531,521,511,523,513,16,2*4,2*12,
     &2*14,2*16,4*2,4*4,2*-11,2*-13,2*-1,2*-3,2*-11,2*-13,2*-1,541,511,
     &521,513,523,21,11,13,15,1,2,3,4,21,22,553,21,2112,2212,2*2112,
     &2212,2112,2*2212,2112,-12,3122,3212,3112,2212,2*2112,-12,2*3122,
     &3222,3112,2212,2112,2212,3122,3222,3212,3122,3112,-12,-14,-12,
     &3322,3312,2*3122,3212,3322,3312,3122,3322,3312,-12,2*4122,7*-11,
     &7*-13,2*2224,2*2212,2*2214,2*3122,2*3212,2*3214,5*3222,4*3224,
     &2*3322,3324,2*2224,7*2212,5*2214,2*2112,2*2114,2*3122,2*3212,
     &2*3214,2*3222,2*3224,4*2,3,2*2,1,2*2,-11,-13,2*2,4*4122,-11,-13,
     &2*2,3*4132,3*4232,-11,-13,2*2,4332,-11,-13,2*2,-11,-13,2*2,-11,
     &-13,2*2,-11,-13,2*2,-11,-13,2*2,-11,-13,2*2,-11,-13,2*2,2*5122,
     &-12,-14,-16,5*4122,441,443,20443,2*-2,2*-4,-2,-4,-12,-14,-16,
     &2*-2,2*-4,-2,-4,-12,-14,-16,2*-2,2*-4,-2,-4,4*5122,-12,-14,-16,
     &2*-2,2*-4,-2,-4,-12,-14,-16,2*-2,2*-4,-2,-4,2*5132,2*5232,-12,
     &-14,-16,2*-2,2*-4,-2,-4,5332,-12,-14,-16,2*-2,2*-4,-2,-4,-12,-14,
     &-16,2*-2,2*-4,-2,-4,-12,-14,-16,2*-2,2*-4,-2,-4,-12,-14,-16,2*-2,
     &2*-4,-2,-4,-12,-14,-16,2*-2,2*-4,-2,-4,-12,-14,-16,2*-2,2*-4,-2,
     &-4,-12,-14,-16,2*-2,2*-4,-2,-4,-12,-14,-16,2*-2,2*-4,-2,-4,-12/
      DATA (KFDP(I,1),I=1429,1710)/-14,-16,2*-2,2*-4,-2,-4,-12,-14,-16,
     &2*-2,2*-4,-2,-4,-12,-14,-16,2*-2,2*-4,-2,-4,-12,-14,-16,2*-2,
     &2*-4,-2,-4,-12,-14,-16,2*-2,2*-4,-2,-4,-12,-14,-16,2*-2,2*-4,-2,
     &-4,-12,-14,-16,2*-2,2*-4,-2,-4,-12,-14,-16,2*-2,2*-4,-2,-4,-12,
     &-14,-16,2*-2,2*-4,-2,-4,-12,-14,-16,2*-2,2*-4,-2,-4,221,223,221,
     &223,211,111,321,130,310,213,113,-213,321,311,321,311,323,313,
     &2*311,321,311,321,313,323,321,211,111,321,130,310,2*211,313,-313,
     &323,-323,421,411,423,413,411,421,413,423,411,421,423,413,443,
     &2*82,521,511,523,513,511,521,513,523,521,511,523,513,511,521,513,
     &523,553,2*21,213,-213,113,213,10211,10111,-10211,2*221,213,2*113,
     &-213,2*321,2*311,113,323,2*313,323,313,-313,323,-323,423,2*413,
     &2*423,413,443,82,523,2*513,2*523,2*513,523,553,21,11,13,82,4*443,
     &10441,20443,445,441,11,13,15,1,2,3,4,21,22,2*553,10551,20553,555,
     &1000039,-1000024,-1000037,1000022,1000023,1000025,1000035,
     &1000002,2000002,1000002,2000002,1000021,1000039,1000024,1000037,
     &1000022,1000023,1000025,1000035,1000001,2000001,1000001,2000001,
     &1000021,1000039,-1000024,-1000037,1000022,1000023,1000025,
     &1000035,1000004,2000004,1000004,2000004,1000021,1000039,1000024,
     &1000037,1000022,1000023,1000025,1000035,1000003,2000003,1000003,
     &2000003,1000021,1000039,-1000024,-1000037,1000022,1000023/
      DATA (KFDP(I,1),I=1711,1900)/1000025,1000035,1000006,2000006,
     &1000006,2000006,1000021,1000039,1000024,1000037,1000022,1000023,
     &1000025,1000035,1000005,2000005,1000005,2000005,1000021,1000022,
     &1000039,-1000024,-1000037,1000022,1000023,1000025,1000035,
     &1000012,2000012,1000012,2000012,1000039,1000024,1000037,1000022,
     &1000023,1000025,1000035,1000011,2000011,1000011,2000011,1000039,
     &-1000024,-1000037,1000022,1000023,1000025,1000035,1000014,
     &2000014,1000014,2000014,1000039,1000024,1000037,1000022,1000023,
     &1000025,1000035,1000013,2000013,1000013,2000013,1000039,-1000024,
     &-1000037,1000022,1000023,1000025,1000035,1000016,2000016,1000016,
     &2000016,1000039,1000024,1000037,1000022,1000023,1000025,1000035,
     &1000015,2000015,1000015,2000015,1000039,1000001,-1000001,2000001,
     &-2000001,1000002,-1000002,2000002,-2000002,1000003,-1000003,
     &2000003,-2000003,1000004,-1000004,2000004,-2000004,1000005,
     &-1000005,2000005,-2000005,1000006,-1000006,2000006,-2000006,
     &6*1000022,6*1000023,6*1000025,6*1000035,1000024,-1000024,1000024,
     &-1000024,1000024,-1000024,1000037,-1000037,1000037,-1000037,
     &1000037,-1000037,10*1000039,16*1000022,1000024,-1000024,1000024,
     &-1000024,1000024,-1000024,1000024,-1000024,1000024,-1000024,
     &1000024,-1000024,1000037,-1000037,1000037,-1000037,1000037/
      DATA (KFDP(I,1),I=1901,2095)/-1000037,1000037,-1000037,1000037,
     &-1000037,1000037,-1000037,1000024,-1000024,1000037,-1000037,
     &1000001,-1000001,2000001,-2000001,1000002,-1000002,2000002,
     &-2000002,1000003,-1000003,2000003,-2000003,1000004,-1000004,
     &2000004,-2000004,1000005,-1000005,2000005,-2000005,1000006,
     &-1000006,2000006,-2000006,1000011,-1000011,2000011,-2000011,
     &1000012,-1000012,2000012,-2000012,1000013,-1000013,2000013,
     &-2000013,1000014,-1000014,2000014,-2000014,1000015,-1000015,
     &2000015,-2000015,1000016,-1000016,2000016,-2000016,5*1000021,
     &2*1000039,6*1000022,6*1000023,6*1000025,6*1000035,1000022,
     &1000023,1000025,1000035,1000002,2000002,-1000001,-2000001,
     &1000004,2000004,-1000003,-2000003,1000006,2000006,-1000005,
     &-2000005,1000012,2000012,-1000011,-2000011,1000014,2000014,
     &-1000013,-2000013,1000016,2000016,-1000015,-2000015,2*1000021,
     &5*1000039,16*1000022,16*1000023,1000024,-1000024,1000024,
     &-1000024,1000024,-1000024,1000024,-1000024,1000024,-1000024,
     &1000024,-1000024,1000037,-1000037,1000037,-1000037,1000037,
     &-1000037,1000037,-1000037,1000037,-1000037,1000037,-1000037,
     &1000024,-1000024,1000037,-1000037,1000001,-1000001,2000001,
     &-2000001,1000002,-1000002,2000002,-2000002,1000003,-1000003/
      DATA (KFDP(I,1),I=2096,2323)/2000003,-2000003,1000004,-1000004,
     &2000004,-2000004,1000005,-1000005,2000005,-2000005,1000006,
     &-1000006,2000006,-2000006,1000011,-1000011,2000011,-2000011,
     &1000012,-1000012,2000012,-2000012,1000013,-1000013,2000013,
     &-2000013,1000014,-1000014,2000014,-2000014,1000015,-1000015,
     &2000015,-2000015,1000016,-1000016,2000016,-2000016,5*1000021,
     &5*1000039,16*1000022,16*1000023,16*1000025,1000024,-1000024,
     &1000024,-1000024,1000024,-1000024,1000024,-1000024,1000024,
     &-1000024,1000024,-1000024,1000037,-1000037,1000037,-1000037,
     &1000037,-1000037,1000037,-1000037,1000037,-1000037,1000037,
     &-1000037,1000024,-1000024,1000037,-1000037,1000001,-1000001,
     &2000001,-2000001,1000002,-1000002,2000002,-2000002,1000003,
     &-1000003,2000003,-2000003,1000004,-1000004,2000004,-2000004,
     &1000005,-1000005,2000005,-2000005,1000006,-1000006,2000006,
     &-2000006,1000011,-1000011,2000011,-2000011,1000012,-1000012,
     &2000012,-2000012,1000013,-1000013,2000013,-2000013,1000014,
     &-1000014,2000014,-2000014,1000015,-1000015,2000015,-2000015,
     &1000016,-1000016,2000016,-2000016,5*1000021,2*1000039,15*1000024,
     &6*1000022,6*1000023,6*1000025,6*1000035,1000022,1000023,1000025,
     &1000035,1000002,2000002,-1000001,-2000001,1000004,2000004/
      DATA (KFDP(I,1),I=2324,4000)/-1000003,-2000003,1000006,2000006,
     &-1000005,-2000005,1000012,2000012,-1000011,-2000011,1000014,
     &2000014,-1000013,-2000013,1000016,2000016,-1000015,-2000015,
     &2*1000021,1000039,-1000024,-1000037,1000022,1000023,1000025,
     &1000035,4*1000001,1000002,2000002,1000002,2000002,1000021,
     &1000039,1000024,1000037,1000022,1000023,1000025,1000035,
     &4*1000002,1000001,2000001,1000001,2000001,1000021,1000039,
     &-1000024,-1000037,1000022,1000023,1000025,1000035,4*1000003,
     &1000004,2000004,1000004,2000004,1000021,1000039,1000024,1000037,
     &1000022,1000023,1000025,1000035,4*1000004,1000003,2000003,
     &1000003,2000003,1000021,1000039,-1000024,-1000037,1000022,
     &1000023,1000025,1000035,4*1000005,1000006,2000006,1000006,
     &2000006,1000021,1000039,1000024,1000037,1000022,1000023,1000025,
     &1000035,4*1000006,1000005,2000005,1000005,2000005,1000021,
     &1000039,-1000024,-1000037,1000022,1000023,1000025,1000035,
     &4*1000011,1000012,2000012,1000012,2000012,1000039,-1000024,
     &-1000037,1000022,1000023,1000025,1000035,4*1000013,1000014,
     &2000014,1000014,2000014,1000039,-1000024,-1000037,1000022,
     &1000023,1000025,1000035,4*1000015,1000016,2000016,1000016,
     &2000016,21,22,23,-24,21,22,23,24,22,23,-24,23,24,1503*0/
      DATA (KFDP(I,2),I=   1, 337)/3*1,2,4,6,8,1,3*2,1,3,5,7,2,3*3,2,4,
     &6,8,3,3*4,1,3,5,7,4,3*5,2,4,6,8,5,3*6,1,3,5,7,6,5,4*1000006,3*7,
     &2,4,6,8,7,4,6,3*8,1,3,5,7,8,5,7,2*11,12,11,12,2*11,2*13,14,13,14,
     &13,11,13,-211,-213,-211,-213,-211,-213,-211,-213,2*-211,-321,
     &-323,-321,2*-323,3*-321,4*-211,-213,-211,-213,-211,-213,-211,
     &-213,-211,-213,3*-211,-213,4*-211,-323,-321,2*-211,2*-321,3*-211,
     &2*15,16,15,16,15,2*17,18,17,2*18,2*17,-1,-2,-3,-4,-5,-6,-7,-8,21,
     &-1,-2,-3,-4,-5,-6,-7,-8,-11,-13,-15,-17,-1,-2,-3,-4,-5,-6,-7,-8,
     &-11,-12,-13,-14,-15,-16,-17,-18,2,4,6,8,2,4,6,8,2,4,6,8,2,4,6,8,
     &12,14,16,18,-1,-2,-3,-4,-5,-6,-7,-8,-11,-13,-15,-17,21,22,2*23,
     &-24,2*1000022,1000023,1000022,1000023,1000025,1000022,1000023,
     &1000025,1000035,-1000024,-1000037,-1000024,-1000037,-1000001,
     &2*-2000001,2000001,-1000002,2*-2000002,2000002,-1000003,
     &2*-2000003,2000003,-1000004,2*-2000004,2000004,-1000005,
     &2*-2000005,2000005,-1000006,2*-2000006,2000006,-1000011,
     &2*-2000011,2000011,-1000012,2*-2000012,2000012,-1000013,
     &2*-2000013,2000013,-1000014,2*-2000014,2000014,-1000015,
     &2*-2000015,2000015,-1000016,2*-2000016,2000016,-1,-2,-3,-4,-5,-6,
     &-7,-8,-11,-12,-13,-14,-15,-16,-17,-18,-24,-37,22,25,2*36,2,4,6,8,
     &2,4,6,8,2,4,6,8,2,4,6,8,12,14,16,18,23,22,25,-1,-2,-3,-4,-5,-6/
      DATA (KFDP(I,2),I= 338, 524)/-7,-8,-11,-13,-15,-17,21,22,2*23,
     &-24,2*25,36,2*1000022,1000023,1000022,1000023,1000025,1000022,
     &1000023,1000025,1000035,-1000024,-1000037,-1000024,-1000037,
     &-1000001,2*-2000001,2000001,-1000002,2*-2000002,2000002,-1000003,
     &2*-2000003,2000003,-1000004,2*-2000004,2000004,-1000005,
     &2*-2000005,2000005,-1000006,2*-2000006,2000006,-1000011,
     &2*-2000011,2000011,-1000012,2*-2000012,2000012,-1000013,
     &2*-2000013,2000013,-1000014,2*-2000014,2000014,-1000015,
     &2*-2000015,2000015,-1000016,2*-2000016,2000016,-1,-2,-3,-4,-5,-6,
     &-7,-8,-11,-13,-15,-17,21,22,2*23,-24,25,2*1000022,1000023,
     &1000022,1000023,1000025,1000022,1000023,1000025,1000035,-1000024,
     &-1000037,-1000024,-1000037,-1000001,2*-2000001,2000001,-1000002,
     &2*-2000002,2000002,-1000003,2*-2000003,2000003,-1000004,
     &2*-2000004,2000004,-1000005,2*-2000005,2000005,-1000006,
     &2*-2000006,2000006,-1000011,2*-2000011,2000011,-1000012,
     &2*-2000012,2000012,-1000013,2*-2000013,2000013,-1000014,
     &2*-2000014,2000014,-1000015,2*-2000015,2000015,-1000016,
     &2*-2000016,2000016,2,4,6,8,12,14,16,18,25,1000024,1000037,
     &1000024,1000037,1000024,1000037,1000024,1000037,2*-1000005,
     &2*-2000005,1000002,1000004,1000012,1000014,2*1000016,-5,-6,21,11/
      DATA (KFDP(I,2),I= 525, 940)/-3,-4,-5,-6,-7,-8,-13,-15,-17,-4,-5,
     &-11,-13,-15,-5,-3,12,14,16,-24,-52,-24,-52,-1,-2,-3,-4,-5,-6,-7,
     &-8,-11,-12,-13,-14,-15,-16,-17,-18,23,51,23,51,2,4,6,8,2,4,6,8,2,
     &4,6,8,2,4,6,8,12,14,16,18,2*51,-1,-2,-3,-4,-5,-6,-7,-8,-11,-12,
     &-13,-14,-15,-16,-17,-18,-82,12,14,-1,-3,11,13,15,1,4,3,4,1,3,22,
     &11,-211,2*22,-13,-11,-211,211,111,211,-321,130,310,22,2*111,-211,
     &11,-11,13,-13,-211,111,22,14,12,111,22,111,3*211,-311,22,211,22,
     &111,-211,211,11,-211,13,22,-211,111,-211,22,111,-11,-211,111,
     &2*-211,-321,130,310,221,111,-211,111,2*0,-211,111,22,-211,111,
     &-211,111,-211,211,-213,113,223,221,14,111,211,111,-11,-13,211,
     &111,22,211,111,211,111,2*211,213,113,223,221,22,-211,111,113,223,
     &22,111,-321,310,211,111,2*-211,221,22,-11,-13,-211,-321,130,310,
     &221,-211,111,11*12,11*14,2*211,2*213,211,20213,2*321,2*323,211,
     &213,211,213,211,213,211,213,211,213,211,213,3*211,213,211,2*321,
     &8*211,2*113,3*211,111,22,211,111,211,111,4*211,8*12,8*14,2*211,
     &2*213,2*111,221,2*113,223,333,20213,211,2*321,323,2*311,313,-211,
     &111,113,2*211,321,2*211,311,321,310,211,-211,4*211,321,4*211,113,
     &2*211,-321,111,22,-211,111,-211,111,-211,211,-211,211,16,5*12,
     &5*14,3*211,3*213,211,2*111,2*113,2*-311,2*-313,-2112,3*321,323,
     &2*-1,22,111,321,311,321,311,-82,-11,-13,-82,22,-82,6*-11,6*-13/
      DATA (KFDP(I,2),I= 941,1318)/2*-15,211,213,20213,211,213,20213,
     &431,433,431,433,311,313,311,313,311,313,-1,-4,-3,-4,-1,-3,22,
     &-211,111,-211,111,-211,211,-211,211,6*-11,6*-13,2*-15,211,213,
     &20213,211,213,20213,431,433,431,433,321,323,321,323,321,323,-1,
     &-4,-3,-4,-1,-3,22,211,111,211,111,4*211,6*-11,6*-13,2*-15,211,
     &213,20213,211,213,20213,431,433,431,433,221,331,333,221,331,333,
     &221,331,333,-1,-4,-3,-4,-1,-3,22,-321,-311,-321,-311,-15,-3,-1,
     &2*-11,2*-13,2*-15,-1,-4,-3,-4,-3,-4,-1,-4,2*12,2*14,2,3,2,3,2*12,
     &2*14,2,1,22,411,421,411,421,21,-11,-13,-15,-1,-2,-3,-4,2*21,22,
     &21,2*-211,111,22,111,211,22,211,-211,11,2*-211,111,-211,111,22,
     &11,22,111,-211,211,111,211,22,211,111,211,-211,22,11,13,11,-211,
     &2*111,2*22,111,211,-321,-211,111,11,2*-211,7*12,7*14,-321,-323,
     &-311,-313,-311,-313,211,213,211,213,211,213,111,221,331,113,223,
     &111,221,113,223,321,323,321,-211,-213,111,221,331,113,223,333,
     &10221,111,221,331,113,223,211,213,211,213,321,323,321,323,321,
     &323,311,313,311,313,2*-1,-3,-1,2203,3201,3203,2203,2101,2103,12,
     &14,-1,-3,2*111,2*211,12,14,-1,-3,22,111,2*22,111,22,12,14,-1,-3,
     &22,12,14,-1,-3,12,14,-1,-3,12,14,-1,-3,12,14,-1,-3,12,14,-1,-3,
     &12,14,-1,-3,12,14,-1,-3,2*-211,11,13,15,-211,-213,-20213,-431,
     &-433,3*3122,1,4,3,4,1,3,11,13,15,1,4,3,4,1,3,11,13,15,1,4,3,4,1/
      DATA (KFDP(I,2),I=1319,1774)/3,2*111,2*211,11,13,15,1,4,3,4,1,3,
     &11,13,15,1,4,3,4,1,3,4*22,11,13,15,1,4,3,4,1,3,22,11,13,15,1,4,3,
     &4,1,3,11,13,15,1,4,3,4,1,3,11,13,15,1,4,3,4,1,3,11,13,15,1,4,3,4,
     &1,3,11,13,15,1,4,3,4,1,3,11,13,15,1,4,3,4,1,3,11,13,15,1,4,3,4,1,
     &3,11,13,15,1,4,3,4,1,3,11,13,15,1,4,3,4,1,3,11,13,15,1,4,3,4,1,3,
     &11,13,15,1,4,3,4,1,3,11,13,15,1,4,3,4,1,3,11,13,15,1,4,3,4,1,3,
     &11,13,15,1,4,3,4,1,3,11,13,15,1,4,3,4,1,3,11,13,15,1,4,3,4,1,3,
     &11,13,15,1,4,3,4,1,3,11,13,15,1,4,3,4,1,3,2*111,2*211,-211,111,
     &-321,130,310,-211,111,211,-211,111,-213,113,-211,111,223,211,111,
     &213,113,211,111,223,-211,111,-321,130,310,2*-211,-311,311,-321,
     &321,211,111,211,111,-211,111,-211,111,311,2*321,311,22,2*-82,
     &-211,111,-211,111,211,111,211,111,-321,-311,-321,-311,411,421,
     &411,421,22,2*21,-211,2*211,111,-211,111,2*211,111,-211,211,111,
     &211,-321,2*-311,-321,22,-211,111,211,111,-311,311,-321,321,211,
     &111,-211,111,321,311,22,-82,-211,111,211,111,-321,-311,411,421,
     &22,21,-11,-13,-82,211,111,221,111,4*22,-11,-13,-15,-1,-2,-3,-4,
     &2*21,211,111,3*22,1,2*2,4*1,2*-24,2*-37,1,2,2*1,4*2,2*24,2*37,2,
     &3,2*4,4*3,2*-24,2*-37,3,4,2*3,4*4,2*24,2*37,4,5,2*6,4*5,2*-24,
     &2*-37,5,6,2*5,4*6,2*24,2*37,6,4,11,2*12,4*11,2*-24,2*-37,12,2*11,
     &4*12,2*24,2*37,13,2*14,4*13,2*-24,2*-37,14,2*13,4*14,2*24,2*37/
      DATA (KFDP(I,2),I=1775,2218)/15,2*16,4*15,2*-24,2*-37,16,2*15,
     &4*16,2*24,2*37,21,-1,1,-1,1,-2,2,-2,2,-3,3,-3,3,-4,4,-4,4,-5,5,
     &-5,5,-6,6,-6,6,1,3,5,2,4,6,1,3,5,2,4,6,1,3,5,2,4,6,1,3,5,2,4,6,1,
     &-1,3,-3,5,-5,1,-1,3,-3,5,-5,22,23,25,35,36,22,23,25,35,36,22,23,
     &11,13,15,12,14,16,1,3,5,2,4,25,35,36,-24,24,11,-11,13,-13,15,-15,
     &1,-1,3,-3,-24,24,11,-11,13,-13,15,-15,1,-1,3,-3,-37,37,-37,37,-1,
     &1,-1,1,-2,2,-2,2,-3,3,-3,3,-4,4,-4,4,-5,5,-5,5,-6,6,-6,6,-11,11,
     &-11,11,-12,12,-12,12,-13,13,-13,13,-14,14,-14,14,-15,15,-15,15,
     &-16,16,-16,16,1,3,5,2,4,24,37,24,-11,-13,-15,-1,-3,24,-11,-13,
     &-15,-1,-3,24,-11,-13,-15,-1,-3,24,-11,-13,-15,-1,-3,4*37,2*-1,
     &2*2,2*-3,2*4,2*-5,2*6,2*-11,2*12,2*-13,2*14,2*-15,2*16,-1,-3,22,
     &23,25,35,36,22,23,11,13,15,12,14,16,1,3,5,2,4,25,35,36,22,23,11,
     &13,15,12,14,16,1,3,5,2,4,25,35,36,-24,24,11,-11,13,-13,15,-15,1,
     &-1,3,-3,-24,24,11,-11,13,-13,15,-15,1,-1,3,-3,-37,37,-37,37,-1,1,
     &-1,1,-2,2,-2,2,-3,3,-3,3,-4,4,-4,4,-5,5,-5,5,-6,6,-6,6,-11,11,
     &-11,11,-12,12,-12,12,-13,13,-13,13,-14,14,-14,14,-15,15,-15,15,
     &-16,16,-16,16,1,3,5,2,4,22,23,25,35,36,22,23,11,13,15,12,14,16,1,
     &3,5,2,4,25,35,36,22,23,11,13,15,12,14,16,1,3,5,2,4,25,35,36,22,
     &23,11,13,15,12,14,16,1,3,5,2,4,25,35,36,-24,24,11,-11,13,-13,15,
     &-15,1,-1,3,-3,-24,24,11,-11,13,-13,15,-15,1,-1,3,-3,-37,37,-37/
      DATA (KFDP(I,2),I=2219,4000)/37,-1,1,-1,1,-2,2,-2,2,-3,3,-3,3,-4,
     &4,-4,4,-5,5,-5,5,-6,6,-6,6,-11,11,-11,11,-12,12,-12,12,-13,13,
     &-13,13,-14,14,-14,14,-15,15,-15,15,-16,16,-16,16,1,3,5,2,4,24,37,
     &23,11,13,15,12,14,16,1,3,5,2,4,25,35,36,24,-11,-13,-15,-1,-3,24,
     &-11,-13,-15,-1,-3,24,-11,-13,-15,-1,-3,24,-11,-13,-15,-1,-3,4*37,
     &2*-1,2*2,2*-3,2*4,2*-5,2*6,2*-11,2*12,2*-13,2*14,2*-15,2*16,-1,
     &-3,1,2*2,4*1,23,25,35,36,2*-24,2*-37,1,2,2*1,4*2,23,25,35,36,
     &2*24,2*37,2,3,2*4,4*3,23,25,35,36,2*-24,2*-37,3,4,2*3,4*4,23,25,
     &35,36,2*24,2*37,4,5,2*6,4*5,23,25,35,36,2*-24,2*-37,5,6,2*5,4*6,
     &23,25,35,36,2*24,2*37,6,11,2*12,4*11,23,25,35,36,2*-24,2*-37,13,
     &2*14,4*13,23,25,35,36,2*-24,2*-37,15,2*16,4*15,23,25,35,36,2*-24,
     &2*-37,3*1,4*2,1,2*11,2*12,11,1503*0/
      DATA (KFDP(I,3),I=   1,1087)/79*0,14,6*0,2*16,2*0,6*111,310,130,
     &2*0,3*111,310,130,321,113,211,223,221,2*113,2*211,2*223,2*221,
     &2*113,221,2*113,2*213,-213,113,2*111,310,130,310,130,2*310,130,
     &470*0,4*3,4*4,1,4,3,2*2,0,-11,8*0,-211,5*0,2*111,211,-211,211,
     &-211,10*0,111,4*0,2*111,-211,-11,11,-13,22,111,3*0,22,3*0,111,
     &211,4*0,111,11*0,111,-211,6*0,-211,3*111,7*0,111,-211,5*0,2*221,
     &3*0,111,5*0,111,11*0,-311,-313,-311,-321,-313,-323,111,221,331,
     &113,223,-311,-313,-311,-321,-313,-323,111,221,331,113,223,22*0,
     &111,113,2*211,-211,-311,211,111,3*211,-211,7*211,7*0,111,-211,
     &111,-211,-321,-323,-311,-321,-313,-323,-211,-213,-321,-323,-311,
     &-321,-313,-323,-211,-213,22*0,111,113,-311,2*-211,211,-211,310,
     &-211,2*111,211,2*-211,-321,-211,2*211,-211,111,-211,2*211,6*0,
     &111,-211,111,-211,0,221,331,333,321,311,221,331,333,321,311,20*0,
     &3,13*0,-411,-413,-10413,-10411,-20413,-415,-411,-413,-10413,
     &-10411,-20413,-415,-411,-413,16*0,-4,-1,-4,-3,2*-2,5*0,111,-211,
     &111,-211,-421,-423,-10423,-10421,-20423,-425,-421,-423,-10423,
     &-10421,-20423,-425,-421,-423,16*0,-4,-1,-4,-3,2*-2,5*0,111,-211,
     &111,-211,-431,-433,-10433,-10431,-20433,-435,-431,-433,-10433,
     &-10431,-20433,-435,-431,-433,19*0,-4,-1,-4,-3,2*-2,8*0,441,443,
     &441,443,441,443,-4,-1,-4,-3,-4,-3,-4,-1,531,533,531,533,3,2,3,2/
      DATA (KFDP(I,3),I=1088,2186)/511,513,511,513,1,2,13*0,2*21,11*0,
     &2112,6*0,2212,12*0,2*3122,3212,10*0,3322,2*0,3122,3212,3214,2112,
     &2114,2212,2112,3122,3212,3214,2112,2114,2212,2112,52*0,3*3,1,6*0,
     &4*3,4*0,4*3,6*0,4*3,0,28*3,2*0,3*4122,8*0,4,1,4,3,2*2,4*4,1,4,3,
     &2*2,4*4,1,4,3,2*2,4*0,4*4,1,4,3,2*2,4*4,1,4,3,2*2,4*0,4*4,1,4,3,
     &2*2,0,4*4,1,4,3,2*2,4*4,1,4,3,2*2,4*4,1,4,3,2*2,4*4,1,4,3,2*2,
     &4*4,1,4,3,2*2,4*4,1,4,3,2*2,4*4,1,4,3,2*2,4*4,1,4,3,2*2,4*4,1,4,
     &3,2*2,4*4,1,4,3,2*2,4*4,1,4,3,2*2,4*4,1,4,3,2*2,4*4,1,4,3,2*2,
     &4*4,1,4,3,2*2,4*4,1,4,3,2*2,4*4,1,4,3,2*2,4*4,1,4,3,2*2,4*4,1,4,
     &3,2*2,31*0,211,111,45*0,-211,2*111,-211,3*111,-211,111,211,30*0,
     &-211,111,13*0,2*21,-211,111,167*0,-1,-3,-5,-2,-4,-6,-1,-3,-5,-2,
     &-4,-6,-1,-3,-5,-2,-4,-6,-1,-3,-5,-2,-4,-6,-2,2,-4,4,-6,6,-2,2,-4,
     &4,-6,6,12*0,-11,-13,-15,-12,-14,-16,-1,-3,-5,-2,-4,5*0,-12,12,
     &-14,14,-16,16,-2,2,-4,4,2*0,-12,12,-14,14,-16,16,-2,2,-4,4,52*0,
     &-1,-3,-5,-2,-4,3*0,12,14,16,2,4,0,12,14,16,2,4,0,12,14,16,2,4,0,
     &12,14,16,2,4,28*0,2,4,7*0,-11,-13,-15,-12,-14,-16,-1,-3,-5,-2,-4,
     &5*0,-11,-13,-15,-12,-14,-16,-1,-3,-5,-2,-4,5*0,-12,12,-14,14,-16,
     &16,-2,2,-4,4,2*0,-12,12,-14,14,-16,16,-2,2,-4,4,52*0,-1,-3,-5,-2,
     &-4,7*0,-11,-13,-15,-12,-14,-16,-1,-3,-5,-2,-4,5*0,-11,-13,-15,
     &-12,-14,-16,-1,-3,-5,-2,-4,5*0,-11,-13,-15,-12,-14,-16,-1,-3,-5/
      DATA (KFDP(I,3),I=2187,4000)/-2,-4,5*0,-12,12,-14,14,-16,16,-2,2,
     &-4,4,2*0,-12,12,-14,14,-16,16,-2,2,-4,4,52*0,-1,-3,-5,-2,-4,3*0,
     &-11,-13,-15,-12,-14,-16,-1,-3,-5,-2,-4,4*0,12,14,16,2,4,0,12,14,
     &16,2,4,0,12,14,16,2,4,0,12,14,16,2,4,28*0,2,4,1657*0/
      DATA (KFDP(I,4),I=   1,4000)/92*0,4*111,6*0,111,2*0,-211,0,-211,
     &3*0,111,2*-211,0,111,0,2*111,113,221,2*111,-213,-211,211,113,
     &6*111,310,2*130,470*0,13*81,41*0,-11,10*0,111,-211,4*0,111,62*0,
     &111,211,111,211,7*0,111,211,111,211,35*0,2*-211,2*111,211,111,
     &-211,2*211,2*-211,13*0,-211,111,-211,111,4*0,-211,111,-211,111,
     &34*0,111,-211,3*111,3*-211,2*111,3*-211,14*0,-321,-311,3*0,-321,
     &-311,20*0,-3,43*0,6*1,39*0,6*2,42*0,6*3,14*0,8*4,4*0,4*-5,4*0,
     &2*-5,67*0,-211,111,5*0,-211,111,52*0,2101,2103,2*2101,6*0,4*81,
     &4*0,4*81,6*0,4*81,0,28*81,13*0,6*2101,18*81,4*0,18*81,4*0,9*81,0,
     &162*81,31*0,-211,111,2450*0/
      DATA (KFDP(I,5),I=   1,4000)/94*0,2*111,17*0,111,7*0,2*111,0,
     &3*111,0,111,665*0,-211,2*111,-211,111,-211,111,65*0,111,-211,
     &3*111,-211,111,3127*0/
 
C...PYDAT4, with particle names (character strings).
      DATA (CHAF(I,1),I=   1, 190)/'d','u','s','c','b','t','b''','t''',
     &2*' ','e-','nu_e','mu-','nu_mu','tau-','nu_tau','tau''-',
     &'nu''_tau',2*' ','g','gamma','Z0','W+','h0',2*' ','reggeon',
     &'pomeron',2*' ','Z''0','Z"0','W''+','H0','A0','H+','eta_tech0',
     &'LQ_ue','R0',10*' ','pi_tech0','pi_tech+','pi''_tech0',
     &'rho_tech0','rho_tech+','omega_tech',24*' ','specflav',
     &'rndmflav','phasespa','c-hadron','b-hadron',5*' ','cluster',
     &'string','indep.','CMshower','SPHEaxis','THRUaxis','CLUSjet',
     &'CELLjet','table',' ','rho_diff0','pi0','rho0','a_20','K_L0',
     &'pi_diffr+','pi+','rho+','a_2+','omega_di','eta','omega','f_2',
     &'K_S0','K0','K*0','K*_20','K+','K*+','K*_2+','phi_diff','eta''',
     &'phi','f''_2','D+','D*+','D*_2+','D0','D*0','D*_20','D_s+',
     &'D*_s+','D*_2s+','J/psi_di','eta_c','J/psi','chi_2c','B0','B*0',
     &'B*_20','B+','B*+','B*_2+','B_s0','B*_s0','B*_2s0','B_c+',
     &'B*_c+','B*_2c+','eta_b','Upsilon','chi_2b','dd_1','Delta-',
     &'ud_0','ud_1','n_diffr0','n0','Delta0','uu_1','p_diffr+','p+',
     &'Delta+','Delta++','sd_0','sd_1','Sigma-','Sigma*-','Lambda0',
     &'su_0','su_1','Sigma0','Sigma*0','Sigma+','Sigma*+','ss_1','Xi-',
     &'Xi*-','Xi0','Xi*0','Omega-','cd_0','cd_1','Sigma_c0',
     &'Sigma*_c0','Lambda_c+','Xi_c0','cu_0','cu_1','Sigma_c+'/
      DATA (CHAF(I,1),I= 191, 317)/'Sigma*_c+','Sigma_c++',
     &'Sigma*_c++','Xi_c+','cs_0','cs_1','Xi''_c0','Xi*_c0','Xi''_c+',
     &'Xi*_c+','Omega_c0','Omega*_c0','cc_1','Xi_cc+','Xi*_cc+',
     &'Xi_cc++','Xi*_cc++','Omega_cc+','Omega*_cc+','Omega*_ccc++',
     &'bd_0','bd_1','Sigma_b-','Sigma*_b-','Lambda_b0','Xi_b-',
     &'Xi_bc0','bu_0','bu_1','Sigma_b0','Sigma*_b0','Sigma_b+',
     &'Sigma*_b+','Xi_b0','Xi_bc+','bs_0','bs_1','Xi''_b-','Xi*_b-',
     &'Xi''_b0','Xi*_b0','Omega_b-','Omega*_b-','Omega_bc0','bc_0',
     &'bc_1','Xi''_bc0','Xi*_bc0','Xi''_bc+','Xi*_bc+','Omega''_bc0',
     &'Omega*_bc0','Omega_bcc+','Omega*_bcc+','bb_1','Xi_bb-',
     &'Xi*_bb-','Xi_bb0','Xi*_bb0','Omega_bb-','Omega*_bb-',
     &'Omega_bbc0','Omega*_bbc0','Omega*_bbb-','a_00','b_10','a_0+',
     &'b_1+','f_0','h_1','K*_00','K_10','K*_0+','K_1+','f''_0','h''_1',
     &'D*_0+','D_1+','D*_00','D_10','D*_0s+','D_1s+','chi_0c','h_1c',
     &'B*_00','B_10','B*_0+','B_1+','B*_0s0','B_1s0','B*_0c+','B_1c+',
     &'chi_0b','h_1b','a_10','a_1+','f_1','K*_10','K*_1+','f''_1',
     &'D*_1+','D*_10','D*_1s+','chi_1c','B*_10','B*_1+','B*_1s0',
     &'B*_1c+','chi_1b','psi''','Upsilon''','~d_L','~u_L','~s_L',
     &'~c_L','~b_1','~t_1','~e_L-','~nu_eL','~mu_L-','~nu_muL',
     &'~tau_1-','~nu_tauL','~g','~chi_10','~chi_20','~chi_1+'/
      DATA (CHAF(I,1),I= 318, 500)/'~chi_30','~chi_40','~chi_2+',
     &'~gravitino','~d_R','~u_R','~s_R','~c_R','~b_2','~t_2','~e_R-',
     &'~nu_eR','~mu_R-','~nu_muR','~tau_2-','~nu_tauR','d*','u*','e*-',
     &'nu*_e0',163*' '/
      DATA (CHAF(I,2),I=   1, 206)/'dbar','ubar','sbar','cbar','bbar',
     &'tbar','b''bar','t''bar',2*' ','e+','nu_ebar','mu+','nu_mubar',
     &'tau+','nu_taubar','tau''+','nu''_taubar',5*' ','W-',9*' ',
     &'W''-',2*' ','H-',' ','LQ_uebar','Rbar0',11*' ','pi_tech-',2*' ',
     &'rho_tech-',26*' ','rndmflavbar',' ','c-hadronbar','b-hadronbar',
     &20*' ','pi_diffr-','pi-','rho-','a_2-',5*' ','Kbar0','K*bar0',
     &'K*_2bar0','K-','K*-','K*_2-',4*' ','D-','D*-','D*_2-','Dbar0',
     &'D*bar0','D*_2bar0','D_s-','D*_s-','D*_2s-',4*' ','Bbar0',
     &'B*bar0','B*_2bar0','B-','B*-','B*_2-','B_sbar0','B*_sbar0',
     &'B*_2sbar0','B_c-','B*_c-','B*_2c-',3*' ','dd_1bar','Deltabar+',
     &'ud_0bar','ud_1bar','n_diffrbar0','nbar0','Deltabar0','uu_1bar',
     &'p_diffrbar-','pbar-','Deltabar-','Deltabar--','sd_0bar',
     &'sd_1bar','Sigmabar+','Sigma*bar+','Lambdabar0','su_0bar',
     &'su_1bar','Sigmabar0','Sigma*bar0','Sigmabar-','Sigma*bar-',
     &'ss_1bar','Xibar+','Xi*bar+','Xibar0','Xi*bar0','Omegabar+',
     &'cd_0bar','cd_1bar','Sigma_cbar0','Sigma*_cbar0','Lambda_cbar-',
     &'Xi_cbar0','cu_0bar','cu_1bar','Sigma_cbar-','Sigma*_cbar-',
     &'Sigma_cbar--','Sigma*_cbar--','Xi_cbar-','cs_0bar','cs_1bar',
     &'Xi''_cbar0','Xi*_cbar0','Xi''_cbar-','Xi*_cbar-','Omega_cbar0',
     &'Omega*_cbar0','cc_1bar','Xi_ccbar-','Xi*_ccbar-','Xi_ccbar--'/
      DATA (CHAF(I,2),I= 207, 324)/'Xi*_ccbar--','Omega_ccbar-',
     &'Omega*_ccbar-','Omega*_cccbar-','bd_0bar','bd_1bar',
     &'Sigma_bbar+','Sigma*_bbar+','Lambda_bbar0','Xi_bbar+',
     &'Xi_bcbar0','bu_0bar','bu_1bar','Sigma_bbar0','Sigma*_bbar0',
     &'Sigma_bbar-','Sigma*_bbar-','Xi_bbar0','Xi_bcbar-','bs_0bar',
     &'bs_1bar','Xi''_bbar+','Xi*_bbar+','Xi''_bbar0','Xi*_bbar0',
     &'Omega_bbar+','Omega*_bbar+','Omega_bcbar0','bc_0bar','bc_1bar',
     &'Xi''_bcbar0','Xi*_bcbar0','Xi''_bcbar-','Xi*_bcbar-',
     &'Omega''_bcba','Omega*_bcbar0','Omega_bccbar-','Omega*_bccbar-',
     &'bb_1bar','Xi_bbbar+','Xi*_bbbar+','Xi_bbbar0','Xi*_bbbar0',
     &'Omega_bbbar+','Omega*_bbbar+','Omega_bbcbar0','Omega*_bbcbar0',
     &'Omega*_bbbbar+',2*' ','a_0-','b_1-',2*' ','K*_0bar0','K_1bar0',
     &'K*_0-','K_1-',2*' ','D*_0-','D_1-','D*_0bar0','D_1bar0',
     &'D*_0s-','D_1s-',2*' ','B*_0bar0','B_1bar0','B*_0-','B_1-',
     &'B*_0sbar0','B_1sbar0','B*_0c-','B_1c-',3*' ','a_1-',' ',
     &'K*_1bar0','K*_1-',' ','D*_1-','D*_1bar0','D*_1s-',' ',
     &'B*_1bar0','B*_1-','B*_1sbar0','B*_1c-',3*' ','~d_Lbar',
     &'~u_Lbar','~s_Lbar','~c_Lbar','~b_1bar','~t_1bar','~e_L+',
     &'~nu_eLbar','~mu_L+','~nu_muLbar','~tau_1+','~nu_tauLbar',3*' ',
     &'~chi_1-',2*' ','~chi_2-',' ','~d_Rbar','~u_Rbar','~s_Rbar'/
      DATA (CHAF(I,2),I= 325, 500)/'~c_Rbar','~b_2bar','~t_2bar',
     &'~e_R+','~nu_eRbar','~mu_R+','~nu_muRbar','~tau_2+',
     &'~nu_tauRbar','d*bar','u*bar','e*bar+','nu*_ebar0',163*' '/
 
C...PYDATR, with initial values for the random number generator.
      DATA MRPY/19780503,0,0,97,33,0/
 
C...Default values for allowed processes and kinematics constraints.
      DATA MSEL/1/
      DATA MSUB/500*0/
      DATA ((KFIN(I,J),J=-40,40),I=1,2)/16*0,4*1,4*0,6*1,5*0,5*1,0,
     &5*1,5*0,6*1,4*0,4*1,16*0,16*0,4*1,4*0,6*1,5*0,5*1,0,5*1,5*0,
     &6*1,4*0,4*1,16*0/
      DATA CKIN/
     &  2.0D0, -1.0D0,  0.0D0, -1.0D0,  1.0D0,
     &  1.0D0,  -10D0,   10D0,  -10D0,   10D0,
     1  -10D0,   10D0,  -10D0,   10D0,  -10D0,
     1   10D0, -1.0D0,  1.0D0, -1.0D0,  1.0D0,
     2  0.0D0,  1.0D0,  0.0D0,  1.0D0, -1.0D0,
     2  1.0D0, -1.0D0,  1.0D0,    0D0,    0D0,
     3  2.0D0, -1.0D0,    0D0,    0D0,  0.0D0,
     3 -1.0D0,  0.0D0, -1.0D0,  4.0D0, -1.0D0,
     4 12.0D0, -1.0D0, 12.0D0, -1.0D0, 12.0D0,
     4 -1.0D0, 12.0D0, -1.0D0,    0D0,    0D0,
     5  0.0D0, -1.0D0,  0.0D0, -1.0D0,  0.0D0,
     5 -1.0D0,    0D0,    0D0,    0D0,    0D0,
     6  140*0D0/
 
C...Default values for main switches and parameters. Reset information.
      DATA (MSTP(I),I=1,100)/
     &  3,    1,    2,    0,    0,    0,    0,    0,    0,    0,
     1  1,    0,    1,    0,    5,    0,    0,    0,    0,    0,
     2  1,    0,    1,    0,    0,    0,    0,    0,    0,    1,
     3  1,    2,    0,    1,    0,    2,    1,    5,    2,    0,
     4  1,    1,    3,    7,    3,    1,    1,    0,    1,    0,
     5  4,    1,    3,    1,    5,    1,    1,    6,    1,    7,
     6  1,    3,    2,    2,    1,    1,    2,    0,    0,    0,
     7  1,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     8  1,    1,  100,    0,    0,    0,    0,    0,    0,    0,
     9  1,    4,    1,    2,    0,    0,    0,    0,    0,    0/
      DATA (MSTP(I),I=101,200)/
     &  3,    1,    0,    0,    0,    0,    0,    0,    0,    0,
     1  1,    1,    1,    0,    0,    0,    0,    0,    0,    0,
     2  0,    1,    2,    1,    1,   50,    0,    0,   10,    0,
     3  0,    4,    0,    1,    0,    0,    0,    0,    0,    0,
     4  0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     5  0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     6  0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     7  0,    2,    0,    0,    0,    0,    0,    0,    0,    0,
     8  6,  115, 1998,   01,   27,    0,    0,    0,    0,    0,
     9  0,    0,    0,    0,    0,    0,    0,    0,    0,    0/
      DATA (PARP(I),I=1,100)/
     &  0.25D0,  10D0, 8*0D0,
     1  0D0,   0D0,  1.0D0, 0.01D0,  0.6D0,  1.0D0,  1.0D0, 3*0D0,
     2  10*0D0,
     3  1.5D0,2.0D0,0.075D0,1.0D0,0.2D0,0D0,2.0D0,0.70D0,0.006D0,0D0,
     4  0.02D0,2.0D0,0.10D0,1000D0,2054D0, 123D0, 246D0, 50D0, 2*0D0,
     5  1.0D0, 9*0D0,
     6  0.25D0, 1.0D0,0.25D0, 1.0D0, 2.0D0,1D-3, 4.0D0,1D-3,2*0D0,
     7  4.0D0, 0.25D0, 8*0D0,
     8  1.40D0,1.55D0,0.5D0, 0.2D0,0.33D0,0.66D0, 0.7D0, 0.5D0,2*0D0,
     9  0.44D0,0.20D0,2.0D0,1.0D0,0D0,3.0D0,1.0D0,0.75D0,0.44D0,2.0D0/
      DATA (PARP(I),I=101,200)/
     &  0.5D0, 0.28D0,  1.0D0, 0.8D0, 6*0D0,
     1  2.0D0, 3*0D0, 1.5D0, 0.5D0, 0.6D0, 2.5D0, 2.0D0, 1.0D0,
     2  1.0D0,  0.4D0, 8*0D0,
     3  0.01D0, 9*0D0,
     4  0.33333D0, 82D0, 1D0, 4D0, 200D0, 5*0D0,
     5  0D0,   0D0,   0D0,   0D0, 6*0D0,
     6  2.20D0, 23.6D0, 18.4D0, 11.5D0, 6*0D0,
     7  0D0,   0D0,   0D0,  1.0D0, 6*0D0,
     8  20*0D0/
      DATA MSTI/200*0/
      DATA PARI/200*0D0/
      DATA MINT/400*0/
      DATA VINT/400*0D0/
 
C...Constants for the generation of the various processes.
      DATA (ISET(I),I=1,100)/
     &  1,    1,    1,   -1,    3,   -1,   -1,    3,   -2,    2,
     1  2,    2,    2,    2,    2,    2,   -1,    2,    2,    2,
     2 -1,    2,    2,    2,    2,    2,   -1,    2,    2,    2,
     3  2,   -1,    2,    2,    2,    2,   -1,   -1,   -1,   -1,
     4 -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
     5 -1,   -1,    2,    2,   -1,   -1,   -1,    2,   -1,   -1,
     6 -1,   -1,   -1,   -1,   -1,   -1,   -1,    2,    2,    2,
     7  4,    4,    4,   -1,   -1,    4,    4,   -1,   -1,    2,
     8  2,    2,    2,    2,    2,    2,    2,    2,    2,   -2,
     9  0,    0,    0,    0,    0,    9,   -2,   -2,   -2,   -2/
      DATA (ISET(I),I=101,200)/
     & -1,    1,    1,   -2,   -2,    2,    2,    2,   -2,    2,
     1  2,    2,    2,    2,    2,   -1,   -1,   -1,   -2,   -2,
     2  5,    5,    5,    5,   -2,   -2,   -2,   -2,   -2,   -2,
     3 -1,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,
     4  1,    1,    1,    1,    1,   -2,    1,    1,    1,   -2,
     5  1,    1,    1,   -2,   -2,    1,    1,    1,   -2,   -2,
     6  2,    2,    2,    2,    2,    2,    2,    2,   -2,   -2,
     7  2,    2,    5,    5,   -2,    2,    2,    5,    5,   -2,
     8  5,    5,   -2,   -2,   -2,    5,    5,   -2,   -2,   -2,
     9  1,    1,    1,    2,   -2,   -2,   -2,   -2,   -2,   -2/
      DATA (ISET(I),I=201,300)/
     &  2,    2,    2,    2,    2,    2,    2,    2,    2,    2,
     1  2,    2,    2,    2,   -2,    2,    2,    2,    2,    2,
     2  2,    2,    2,    2,    2,    2,    2,    2,    2,    2,
     3  2,    2,    2,    2,    2,    2,    2,    2,    2,    2,
     4  2,    2,    2,    2,   -1,    2,    2,    2,    2,    2,
     5  2,    2,    2,    2,   -1,    2,   -1,    2,    2,   -2,
     6  2,    2,    2,    2,    2,   -1,   -1,   -1,   -1,   -1,
     7  2,    2,    2,    2,    2,    2,    2,    2,    2,    2,
     8 -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,
     9 -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2,   -2/
      DATA (ISET(I),I=301,500)/200*-2/
      DATA ((KFPR(I,J),J=1,2),I=1,50)/
     &  23,    0,   24,    0,   25,    0,   24,    0,   25,    0,
     &  24,    0,   23,    0,   25,    0,    0,    0,    0,    0,
     1   0,    0,    0,    0,   21,   21,   21,   22,   21,   23,
     1  21,   24,   21,   25,   22,   22,   22,   23,   22,   24,
     2  22,   25,   23,   23,   23,   24,   23,   25,   24,   24,
     2  24,   25,   25,   25,    0,   21,    0,   22,    0,   23,
     3   0,   24,    0,   25,    0,   21,    0,   22,    0,   23,
     3   0,   24,    0,   25,    0,   21,    0,   22,    0,   23,
     4   0,   24,    0,   25,    0,   21,    0,   22,    0,   23,
     4   0,   24,    0,   25,    0,   21,    0,   22,    0,   23/
      DATA ((KFPR(I,J),J=1,2),I=51,100)/
     5   0,   24,    0,   25,    0,    0,    0,    0,    0,    0,
     5   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     6   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     6   0,    0,    0,    0,   21,   21,   24,   24,   23,   24,
     7  23,   23,   24,   24,   23,   24,   23,   25,   22,   22,
     7  23,   23,   24,   24,   24,   25,   25,   25,    0,  211,
     8   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     8 443,   21,10441,   21,20443,   21,  445,   21,    0,    0,
     9   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     9   0,    0,    0,    0,    0,    0,    0,    0,    0,    0/
      DATA ((KFPR(I,J),J=1,2),I=101,150)/
     &  23,    0,   25,    0,   25,    0,    0,    0,    0,    0,
     & 443,   22,  443,   21,  443,   22,    0,    0,   22,   25,
     1  21,   25,    0,   25,   21,   25,   22,   22,   21,   22,
     1  22,   23,   23,   23,   24,   24,    0,    0,    0,    0,
     2  25,    6,   25,    6,   25,    0,   25,    0,    0,    0,
     2   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     3  23,    5,    0,    0,    0,    0,    0,    0,    0,    0,
     3   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     4  32,    0,   34,    0,   37,    0,   40,    0,   39,    0,
     4   0,    0, 4000001, 0, 4000002, 0,   38,    0,    0,    0/
      DATA ((KFPR(I,J),J=1,2),I=151,200)/
     5  35,    0,   35,    0,   35,    0,    0,    0,    0,    0,
     5  36,    0,   36,    0,   36,    0,    0,    0,    0,    0,
     6   6,   37,   39,    0,   39,   39,   39,   39,   11,    0,
     6  11,    0, 0, 4000001, 0, 4000002,    0,    0,    0,    0,
     7  23,   35,   24,   35,   35,    0,   35,    0,    0,    0,
     7  23,   36,   24,   36,   36,    0,   36,    0,    0,    0,
     8  35,    6,   35,    6,    0,    0,    0,    0,    0,    0,
     8  36,    6,   36,    6,    0,    0,    0,    0,    0,    0,
     9  54,    0,   55,    0,   56,    0,   11,    0,    0,    0,
     9   0,    0,    0,    0,    0,    0,    0,    0,    0,    0/
      DATA ((KFPR(I,J),J=1,2),I=201,240)/
     &  1000011,   1000011,   2000011,   2000011,   1000011,
     &  2000011,   1000013,   1000013,   2000013,   2000013,
     &  1000013,   2000013,   1000015,   1000015,   2000015,
     &  2000015,   1000015,   2000015,   1000011,   1000012,
     1  1000015,   1000016,   2000015,   1000016,   1000012,
     1  1000012,   1000016,   1000016,         0,         0,
     1  1000022,   1000022,   1000023,   1000023,   1000025,
     1  1000025,   1000035,   1000035,   1000022,   1000023,
     2  1000022,   1000025,   1000022,   1000035,   1000023,
     2  1000025,   1000023,   1000035,   1000025,   1000035,
     2  1000024,   1000024,   1000037,   1000037,   1000024,
     2  1000037,   1000022,   1000024,   1000023,   1000024,
     3  1000025,   1000024,   1000035,   1000024,   1000022,
     3  1000037,   1000023,   1000037,   1000025,   1000037,
     3  1000035,   1000037,   1000021,   1000022,   1000021,
     3  1000023,   1000021,   1000025,   1000021,   1000035/
      DATA ((KFPR(I,J),J=1,2),I=241,280)/
     4  1000021,   1000024,   1000021,   1000037,   1000021,
     4  1000021,   1000021,   1000021,         0,         0,
     4  1000002,   1000022,   2000002,   1000022,   1000002,
     4  1000023,   2000002,   1000023,   1000002,   1000025,
     5  2000002,   1000025,   1000002,   1000035,   2000002,
     5  1000035,   1000001,   1000024,   2000005,   1000024,
     5  1000001,   1000037,   2000005,   1000037,   1000002,
     5  1000021,   2000002,   1000021,         0,         0,
     6  1000006,   1000006,   2000006,   2000006,   1000006,
     6  2000006,   1000006,   1000006,   2000006,   2000006,
     6        0,         0,         0,         0,         0,
     6        0,         0,         0,         0,         0,
     7  1000002,   1000002,   2000002,   2000002,   1000002,
     7  2000002,   1000002,   1000002,   2000002,   2000002,
     7  1000002,   2000002,   1000002,   1000002,   2000002,
     7  2000002,   1000002,   1000002,   2000002,   2000002/
      DATA ((KFPR(I,J),J=1,2),I=281,500)/440*0/
      DATA COEF/10000*0D0/
      DATA (((ICOL(I,J,K),K=1,2),J=1,4),I=1,40)/
     &4,0,3,0,2,0,1,0,3,0,4,0,1,0,2,0,2,0,0,1,4,0,0,3,3,0,0,4,1,0,0,2,
     &3,0,0,4,1,4,3,2,4,0,0,3,4,2,1,3,2,0,4,1,4,0,2,3,4,0,3,4,2,0,1,2,
     &3,2,1,0,1,4,3,0,4,3,3,0,2,1,1,0,3,2,1,4,1,0,0,2,2,4,3,1,2,0,0,1,
     &3,2,1,4,1,4,3,2,4,2,1,3,4,2,1,3,3,4,4,3,1,2,2,1,2,0,3,1,2,0,0,0,
     &4,2,1,0,0,0,1,0,3,0,0,3,1,2,0,0,4,0,0,4,0,0,1,2,2,0,0,1,4,4,3,3,
     &2,2,1,1,4,4,3,3,3,3,4,4,1,1,2,2,3,2,1,3,1,2,0,0,4,2,1,4,0,0,1,2,
     &4,0,0,0,4,0,1,3,0,0,3,0,2,4,3,0,3,4,0,0,1,0,0,1,0,0,3,4,2,0,0,2,
     &3,0,0,0,1,0,0,0,0,0,3,0,2,0,0,0,2,0,3,1,2,0,0,0,3,2,1,0,1,0,0,0,
     &4,4,3,3,2,2,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
 
C...Treatment of resonances.
      DATA (MWID(I)  ,I=   1, 500)/5*0,3*1,8*0,1,5*0,3*1,6*0,1,0,7*1,
     &10*0,2*1,0,3*1,245*0,19*2,0,7*2,0,2,0,2,0,4*1,163*0/
 
C...Character constants: name of processes.
      DATA PROC(0)/                    'All included subprocesses   '/
      DATA (PROC(I),I=1,20)/
     &'f + fbar -> gamma*/Z0       ',  'f + fbar'' -> W+/-           ',
     &'f + fbar -> h0              ',  'gamma + W+/- -> W+/-        ',
     &'Z0 + Z0 -> h0               ',  'Z0 + W+/- -> W+/-           ',
     &'                            ',  'W+ + W- -> h0               ',
     &'                            ',  'f + f'' -> f + f'' (QFD)      ',
     1'f + f'' -> f + f'' (QCD)      ','f + fbar -> f'' + fbar''      ',
     1'f + fbar -> g + g           ',  'f + fbar -> g + gamma       ',
     1'f + fbar -> g + Z0          ',  'f + fbar'' -> g + W+/-       ',
     1'f + fbar -> g + h0          ',  'f + fbar -> gamma + gamma   ',
     1'f + fbar -> gamma + Z0      ',  'f + fbar'' -> gamma + W+/-   '/
      DATA (PROC(I),I=21,40)/
     2'f + fbar -> gamma + h0      ',  'f + fbar -> Z0 + Z0         ',
     2'f + fbar'' -> Z0 + W+/-      ', 'f + fbar -> Z0 + h0         ',
     2'f + fbar -> W+ + W-         ',  'f + fbar'' -> W+/- + h0      ',
     2'f + fbar -> h0 + h0         ',  'f + g -> f + g              ',
     2'f + g -> f + gamma          ',  'f + g -> f + Z0             ',
     3'f + g -> f'' + W+/-          ', 'f + g -> f + h0             ',
     3'f + gamma -> f + g          ',  'f + gamma -> f + gamma      ',
     3'f + gamma -> f + Z0         ',  'f + gamma -> f'' + W+/-      ',
     3'f + gamma -> f + h0         ',  'f + Z0 -> f + g             ',
     3'f + Z0 -> f + gamma         ',  'f + Z0 -> f + Z0            '/
      DATA (PROC(I),I=41,60)/
     4'f + Z0 -> f'' + W+/-         ', 'f + Z0 -> f + h0            ',
     4'f + W+/- -> f'' + g          ', 'f + W+/- -> f'' + gamma      ',
     4'f + W+/- -> f'' + Z0         ', 'f + W+/- -> f'' + W+/-       ',
     4'f + W+/- -> f'' + h0         ', 'f + h0 -> f + g             ',
     4'f + h0 -> f + gamma         ',  'f + h0 -> f + Z0            ',
     5'f + h0 -> f'' + W+/-         ', 'f + h0 -> f + h0            ',
     5'g + g -> f + fbar           ',  'g + gamma -> f + fbar       ',
     5'g + Z0 -> f + fbar          ',  'g + W+/- -> f + fbar''       ',
     5'g + h0 -> f + fbar          ',  'gamma + gamma -> f + fbar   ',
     5'gamma + Z0 -> f + fbar      ',  'gamma + W+/- -> f + fbar''   '/
      DATA (PROC(I),I=61,80)/
     6'gamma + h0 -> f + fbar      ',  'Z0 + Z0 -> f + fbar         ',
     6'Z0 + W+/- -> f + fbar''      ', 'Z0 + h0 -> f + fbar         ',
     6'W+ + W- -> f + fbar         ',  'W+/- + h0 -> f + fbar''      ',
     6'h0 + h0 -> f + fbar         ',  'g + g -> g + g              ',
     6'gamma + gamma -> W+ + W-    ',  'gamma + W+/- -> Z0 + W+/-   ',
     7'Z0 + Z0 -> Z0 + Z0          ',  'Z0 + Z0 -> W+ + W-          ',
     7'Z0 + W+/- -> Z0 + W+/-      ',  'Z0 + Z0 -> Z0 + h0          ',
     7'W+ + W- -> gamma + gamma    ',  'W+ + W- -> Z0 + Z0          ',
     7'W+/- + W+/- -> W+/- + W+/-  ',  'W+/- + h0 -> W+/- + h0      ',
     7'h0 + h0 -> h0 + h0          ',  'q + gamma -> q'' + pi+/-     '/
      DATA (PROC(I),I=81,100)/
     8'q + qbar -> Q + Qbar, mass  ',  'g + g -> Q + Qbar, massive  ',
     8'f + q -> f'' + Q, massive    ', 'g + gamma -> Q + Qbar, mass ',
     8'gamma + gamma -> F + Fbar, m',  'g + g -> J/Psi + g          ',
     8'g + g -> chi_0c + g         ',  'g + g -> chi_1c + g         ',
     8'g + g -> chi_2c + g         ',  '                            ',
     9'Elastic scattering          ',  'Single diffractive (XB)     ',
     9'Single diffractive (AX)     ',  'Double  diffractive         ',
     9'Low-pT scattering           ',  'Semihard QCD 2 -> 2         ',
     9'                            ',  '                            ',
     9'                            ',  '                            '/
      DATA (PROC(I),I=101,120)/
     &'g + g -> gamma*/Z0          ',  'g + g -> h0                 ',
     &'gamma + gamma -> h0         ',  '                            ',
     &'                            ',  'g + g -> J/Psi + gamma      ',
     &'gamma + g -> J/Psi + g      ',  'gamma+gamma -> J/Psi + gamma',
     &'                            ',  'f + fbar -> gamma + h0      ',
     1'f + fbar -> g + h0          ',  'q + g -> q + h0             ',
     1'g + g -> g + h0             ',  'g + g -> gamma + gamma      ',
     1'g + g -> g + gamma          ',  'g + g -> gamma + Z0         ',
     1'g + g -> Z0 + Z0            ',  'g + g -> W+ + W-            ',
     1'                            ',  '                            '/
      DATA (PROC(I),I=121,140)/
     2'g + g -> Q + Qbar + h0      ',  'q + qbar -> Q + Qbar + h0   ',
     2'f + f'' -> f + f'' + h0       ',
     2'f + f'' -> f" + f"'' + h0     ',
     2'                            ',  '                            ',
     2'                            ',  '                            ',
     2'                            ',  '                            ',
     3'g + g -> Z0 + q + qbar      ',  '                            ',
     3'                            ',  '                            ',
     3'                            ',  '                            ',
     3'                            ',  '                            ',
     3'                            ',  '                            '/
      DATA (PROC(I),I=141,160)/
     4'f + fbar -> gamma*/Z0/Z''0   ', 'f + fbar'' -> W''+/-          ',
     4'f + fbar'' -> H+/-           ', 'f + fbar'' -> R              ',
     4'q + l -> LQ                 ',  '                            ',
     4'd + g -> d*                 ',  'u + g -> u*                 ',
     4'g + g -> eta_techni         ',  '                            ',
     5'f + fbar -> H0              ',  'g + g -> H0                 ',
     5'gamma + gamma -> H0         ',  '                            ',
     5'                            ',  'f + fbar -> A0              ',
     5'g + g -> A0                 ',  'gamma + gamma -> A0         ',
     5'                            ',  '                            '/
      DATA (PROC(I),I=161,180)/
     6'f + g -> f'' + H+/-          ', 'q + g -> LQ + lbar          ',
     6'g + g -> LQ + LQbar         ',  'q + qbar -> LQ + LQbar      ',
     6'f + fbar -> f'' + fbar'' (g/Z)',
     6'f +fbar'' -> f" + fbar"'' (W) ',
     6'q + q'' -> q" + d*           ',  'q + q'' -> q" + u*           ',
     6'                            ',  '                            ',
     7'f + fbar -> Z0 + H0         ', 'f + fbar'' -> W+/- + H0      ',
     7'f + f'' -> f + f'' + H0       ',
     7'f + f'' -> f" + f"'' + H0     ',
     7'                            ',  'f + fbar -> Z0 + A0         ',
     7'f + fbar'' -> W+/- + A0      ',
     7'f + f'' -> f + f'' + A0       ',
     7'f + f'' -> f" + f"'' + A0     ',
     7'                            '/
      DATA (PROC(I),I=181,200)/
     8'g + g -> Q + Qbar + H0      ',  'q + qbar -> Q + Qbar + H0   ',
     8'                            ',  '                            ',
     8'                            ',  'g + g -> Q + Qbar + A0      ',
     8'q + qbar -> Q + Qbar + A0   ',  '                            ',
     8'                            ',  '                            ',
     9'f + fbar -> rho_tech0       ',  'f + f'' -> rho_tech+/-       ',
     9'f + fbar -> omega_tech0     ',  'f+fbar -> f''+fbar'' (technic)',
     9'                            ',  '                            ',
     9'                            ',  '                            ',
     9'                            ',  '                            '/
      DATA (PROC(I),I=201,220)/
     &'f + fbar -> ~e_L + ~e_Lbar  ',  'f + fbar -> ~e_R + ~e_Rbar  ',
     &'f + fbar -> ~e_R + ~e_Lbar  ',  'f + fbar -> ~mu_L + ~mu_Lbar',
     &'f + fbar -> ~mu_R + ~mu_Rbar',  'f + fbar -> ~mu_L + ~mu_Rbar',
     &'f+fbar -> ~tau_1 + ~tau_1bar',  'f+fbar -> ~tau_2 + ~tau_2bar',
     &'f+fbar -> ~tau_1 + ~tau_2bar',  'q + qbar'' -> ~l_L + ~nulbar ',
     1'q+qbar''-> ~tau_1 + ~nutaubar', 'q+qbar''-> ~tau_2 + ~nutaubar',
     1'f + fbar -> ~nul + ~nulbar  ',  'f+fbar -> ~nutau + ~nutaubar',
     1'                            ',  'f + fbar -> ~chi1 + ~chi1   ',
     1'f + fbar -> ~chi2 + ~chi2   ',  'f + fbar -> ~chi3 + ~chi3   ',
     1'f + fbar -> ~chi4 + ~chi4   ',  'f + fbar -> ~chi1 + ~chi2   '/
      DATA (PROC(I),I=221,240)/
     2'f + fbar -> ~chi1 + ~chi3   ',  'f + fbar -> ~chi1 + ~chi4   ',
     2'f + fbar -> ~chi2 + ~chi3   ',  'f + fbar -> ~chi2 + ~chi4   ',
     2'f + fbar -> ~chi3 + ~chi4   ',  'f+fbar -> ~chi+-1 + ~chi-+1 ',
     2'f+fbar -> ~chi+-2 + ~chi-+2 ',  'f+fbar -> ~chi+-1 + ~chi-+2 ',
     2'q + qbar'' -> ~chi1 + ~chi+-1', 'q + qbar'' -> ~chi2 + ~chi+-1',
     3'q + qbar'' -> ~chi3 + ~chi+-1', 'q + qbar'' -> ~chi4 + ~chi+-1',
     3'q + qbar'' -> ~chi1 + ~chi+-2', 'q + qbar'' -> ~chi2 + ~chi+-2',
     3'q + qbar'' -> ~chi3 + ~chi+-2', 'q + qbar'' -> ~chi4 + ~chi+-2',
     3'q + qbar -> ~chi1 + ~g      ',  'q + qbar -> ~chi2 + ~g      ',
     3'q + qbar -> ~chi3 + ~g      ',  'q + qbar -> ~chi4 + ~g      '/
      DATA (PROC(I),I=241,260)/
     4'q + qbar'' -> ~chi+-1 + ~g   ', 'q + qbar'' -> ~chi+-2 + ~g  ',
     4'q + qbar -> ~g + ~g         ',  'g + g -> ~g + ~g            ',
     4'                            ',  'qj + g -> ~qj_L + ~chi1     ',
     4'qj + g -> ~qj_R + ~chi1     ',  'qj + g -> ~qj_L + ~chi2     ',
     4'qj + g -> ~qj_R + ~chi2     ',  'qj + g -> ~qj_L + ~chi3     ',
     5'qj + g -> ~qj_R + ~chi3     ',  'qj + g -> ~qj_L + ~chi4     ',
     5'qj + g -> ~qj_R + ~chi4     ',  'qj + g -> ~qk_L + ~chi+-1   ',
     5'qj + g -> ~qk_R + ~chi+-1   ',  'qj + g -> ~qk_L + ~chi+-2   ',
     5'qj + g -> ~qk_R + ~chi+-2   ',  'qj + g -> ~qj_L + ~g        ',
     5'qj + g -> ~qj_R + ~g        ',  '                            '/
      DATA (PROC(I),I=261,280)/
     6'f + fbar -> ~t_1 + ~t_1bar  ',  'f + fbar -> ~t_2 + ~t_2bar  ',
     6'f + fbar -> ~t_1 + ~t_2bar  ',  'g + g -> ~t_1 + ~t_1bar     ',
     6'g + g -> ~t_2 + ~t_2bar     ',  '                            ',
     6'                            ',  '                            ',
     6'                            ',  '                            ',
     7'qi + qj -> ~qi_L + ~qj_L    ',  'qi + qj -> ~qi_R + ~qj_R    ',
     7'qi + qj -> ~qi_L + ~qj_R    ',  'qi+qjbar -> ~qi_L + ~qj_Lbar',
     7'qi+qjbar -> ~qi_R + ~qj_Rbar',  'qi+qjbar -> ~qi_L + ~qj_Rbar',
     7'f + fbar -> ~qi_L + ~qi_Lbar',  'f + fbar -> ~qi_R + ~qi_Rbar',
     7'g + g -> ~qi_L + ~qi_Lbar   ',  'g + g -> ~qi_R + ~qi_Rbar   '/
      DATA (PROC(I),I=281,500)/220*'                            '/
 
C...Cross sections and slope offsets.
      DATA SIGT/294*0D0/
 
C...Supersymmetry switches and parameters.
      DATA IMSS/0,
     &  0,  0,  0,  1,  0,  0,  0,  1,  0,  0,
     1  89*0/
      DATA RMSS/0D0,
     &  80D0,160D0,500D0,800D0,2D0,250D0,200D0,800D0,700D0,800D0,
     1  700D0,500D0,250D0,200D0,800D0,400D0,0D0,0.1D0,850D0,0.041D0,
     2   1D0,800D0,1D4,1D4,1D4,0D0,0D0,24D17,2*0D0,
     3  69*0D0/
 
C...Data for histogramming routines.
      DATA IHIST/1000,20000,55,1/
      DATA INDX/1000*0/
 
      END

C*********************************************************************
 
C...PYCOMP
C...Compress the standard KF codes for use in mass and decay arrays;
C...also checks whether a given code actually is defined.
 
      FUNCTION PYCOMP(KF)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYDAT1/,/PYDAT2/
C...Local arrays and saved data.
      DIMENSION KFORD(100:500),KCORD(101:500)
      SAVE KFORD,KCORD,NFORD,KFLAST,KCLAST
 
C...Whenever necessary reorder codes for faster search.
      IF(MSTU(20).EQ.0) THEN
        NFORD=100
        KFORD(100)=0
        DO 120 I=101,500
          KFA=KCHG(I,4)
          IF(KFA.LE.100) GOTO 120
          NFORD=NFORD+1
          DO 100 I1=NFORD-1,0,-1
            IF(KFA.GE.KFORD(I1)) GOTO 110
            KFORD(I1+1)=KFORD(I1)
            KCORD(I1+1)=KCORD(I1)
  100     CONTINUE
  110     KFORD(I1+1)=KFA
          KCORD(I1+1)=I
  120   CONTINUE
        MSTU(20)=1
        KFLAST=0
        KCLAST=0
      ENDIF
 
C...Fast action if same code as in latest call.
      IF(KF.EQ.KFLAST) THEN
        PYCOMP=KCLAST
        RETURN
      ENDIF
 
C...Starting values. Remove internal diquark flags.
      PYCOMP=0
      KFA=IABS(KF)
      IF(MOD(KFA/10,10).EQ.0.AND.KFA.LT.100000
     &     .AND.MOD(KFA/1000,10).GT.0) KFA=MOD(KFA,10000)
 
C...Simple cases: direct translation.
      IF(KFA.GT.KFORD(NFORD)) THEN
      ELSEIF(KFA.LE.100) THEN
        PYCOMP=KFA
 
C...Else binary search.
      ELSE
        IMIN=100
        IMAX=NFORD+1
  130   IAVG=(IMIN+IMAX)/2
        IF(KFORD(IAVG).GT.KFA) THEN
          IMAX=IAVG
          IF(IMAX.GT.IMIN+1) GOTO 130
        ELSEIF(KFORD(IAVG).LT.KFA) THEN
          IMIN=IAVG
          IF(IMAX.GT.IMIN+1) GOTO 130
        ELSE
          PYCOMP=KCORD(IAVG)
        ENDIF
      ENDIF
 
C...Check if antiparticle allowed.
      IF(PYCOMP.NE.0.AND.KF.LT.0) THEN
        IF(KCHG(PYCOMP,3).EQ.0) PYCOMP=0
      ENDIF
 
C...Save codes for possible future fast action.
      KFLAST=KF
      KCLAST=PYCOMP
 
      RETURN
      END

C*********************************************************************
 
C...PYDCYK
C...Handles flavour production in the decay of unstable particles
C...and small string clusters.
 
      SUBROUTINE PYDCYK(KFL1,KFL2,KFL3,KF)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYDAT1/,/PYDAT2/
 

C.. Call PYKFDI directly if no popcorn option is on
      IF(MSTJ(12).LT.2) THEN
         CALL PYKFDI(KFL1,KFL2,KFL3,KF)
         MSTU(124)=KFL3
         RETURN
      ENDIF
 
      KFL3=0
      KF=0
      IF(KFL1.EQ.0) RETURN
      KF1A=IABS(KFL1)
      KF2A=IABS(KFL2)
 
      NSTO=130
      NMAX=MIN(MSTU(125),10)
 
C.. Identify rank 0 cluster qq
      IRANK=1
      IF(KF1A.GT.10.AND.KF1A.LT.10000) IRANK=0
 
      IF(KF2A.GT.0)THEN
C.. Join jets: Fails if store not empty
         IF(MSTU(121).GT.0) THEN
            MSTU(121)=0
            RETURN
         ENDIF
         CALL PYKFDI(KFL1,KFL2,KFL3,KF)
      ELSEIF(KF1A.GT.10.AND.MSTU(121).GT.0)THEN
C.. Pick popcorn meson from store, return same qq, decrease store
         KF=MSTU(NSTO+MSTU(121))
         KFL3=-KFL1
         MSTU(121)=MSTU(121)-1
      ELSE
C.. Generate new flavour. Then done if no diquark is generated
  100    CALL PYKFDI(KFL1,0,KFL3,KF)
         IF(MSTU(121).EQ.-1) GOTO 100
         MSTU(124)=KFL3
         IF(KF.EQ.0.OR.IABS(KFL3).LE.10) RETURN
 
C.. Simple case if no dynamical popcorn suppressions are considered
         IF(MSTJ(12).LT.4) THEN
            IF(MSTU(121).EQ.0) RETURN
            NMES=1
            KFPREV=-KFL3
            CALL PYKFDI(KFPREV,0,KFL3,KFM)
C.. Due to eta+eta' suppr., a qq->M+qq attempt might end as qq->B+q
            IF(IABS(KFL3).LE.10)THEN
               KFL3=-KFPREV
               RETURN
            ENDIF
            GOTO 120
         ENDIF
 
C test output qq against fake Gamma, then return if no popcorn.
         GB=2D0
         IF(IRANK.NE.0)THEN
            CALL PYZDIS(1,2103,5D0,Z)
            GB=3D0*(1D0-Z)/Z
            IF(1D0-PARF(192)**GB.LT.PYR(0)) THEN
               MSTU(121)=0
               GOTO 100
            ENDIF
         ENDIF      
         IF(MSTU(121).EQ.0) RETURN
 
C..Set store size memory. Pick fake dynamical variables of qq.
         NMES=MSTU(121)
         CALL PYPTDI(1,PX3,PY3)
         X=1D0
         POPM=0D0
         G=GB
         POPG=GB
 
C.. Pick next popcorn meson, test with fake dynamical variables
  110    KFPREV=-KFL3
         PX1=-PX3
         PY1=-PY3
         CALL PYKFDI(KFPREV,0,KFL3,KFM)
         IF(MSTU(121).EQ.-1) GOTO 100
         CALL PYPTDI(KFL3,PX3,PY3)
         PM=PYMASS(KFM)**2+(PX1+PX3)**2+(PY1+PY3)**2
         CALL PYZDIS(KFPREV,KFL3,PM,Z)
         G=(1D0-Z)*(G+PM/Z)
         X=(1D0-Z)*X
 
         PTST=1D0
         GTST=1D0
         RTST=PYR(0)
         IF(MSTJ(12).GT.4)THEN
            POPMN=SQRT((1D0-X)*(G/X-GB))
            POPM=POPM+PMAS(PYCOMP(KFM),1)-PMAS(PYCOMP(KFM),3)
            PTST=EXP((POPM-POPMN)*PARF(193))
            POPM=POPMN
         ENDIF
         IF(IRANK.NE.0)THEN
            POPGN=X*GB
            GTST=(1D0-PARF(192)**POPGN)/(1D0-PARF(192)**POPG)
            POPG=POPGN
         ENDIF
         IF(RTST.GT.PTST*GTST)THEN
            MSTU(121)=0
            IF(RTST.GT.PTST) MSTU(121)=-1
            GOTO 100
         ENDIF
 
C.. Store meson
  120    IF(NMES.LE.NMAX) MSTU(NSTO+MSTU(121)+1)=KFM
         IF(MSTU(121).GT.0) GOTO 110
 
C.. Test accepted system size. If OK set global popcorn size variable.
         IF(NMES.GT.NMAX)THEN
            KF=0
            KFL3=0
            RETURN
         ENDIF
         MSTU(121)=NMES
      ENDIF
 
      RETURN
      END

C*********************************************************************
 
C...PYROBO
C...Performs rotations and boosts.
 
      SUBROUTINE PYROBO(IMI,IMA,THE,PHI,BEX,BEY,BEZ)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /PYJETS/,/PYDAT1/
C...Local arrays.
      DIMENSION ROT(3,3),PR(3),VR(3),DP(4),DV(4)
 
C...Find and check range of rotation/boost.
      IMIN=IMI
      IF(IMIN.LE.0) IMIN=1
      IF(MSTU(1).GT.0) IMIN=MSTU(1)
      IMAX=IMA
      IF(IMAX.LE.0) IMAX=N
      IF(MSTU(2).GT.0) IMAX=MSTU(2)
      IF(IMIN.GT.MSTU(4).OR.IMAX.GT.MSTU(4)) THEN
        CALL PYERRM(11,'(PYROBO:) range outside PYJETS memory')
        RETURN
      ENDIF
 
C...Optional resetting of V (when not set before.)
      IF(MSTU(33).NE.0) THEN
        DO 110 I=MIN(IMIN,MSTU(4)),MIN(IMAX,MSTU(4))
          DO 100 J=1,5
            V(I,J)=0D0
  100     CONTINUE
  110   CONTINUE
        MSTU(33)=0
      ENDIF
 
C...Rotate, typically from z axis to direction (theta,phi).
      IF(THE**2+PHI**2.GT.1D-20) THEN
        ROT(1,1)=COS(THE)*COS(PHI)
        ROT(1,2)=-SIN(PHI)
        ROT(1,3)=SIN(THE)*COS(PHI)
        ROT(2,1)=COS(THE)*SIN(PHI)
        ROT(2,2)=COS(PHI)
        ROT(2,3)=SIN(THE)*SIN(PHI)
        ROT(3,1)=-SIN(THE)
        ROT(3,2)=0D0
        ROT(3,3)=COS(THE)
        DO 140 I=IMIN,IMAX
          IF(K(I,1).LE.0) GOTO 140
          DO 120 J=1,3
            PR(J)=P(I,J)
            VR(J)=V(I,J)
  120     CONTINUE
          DO 130 J=1,3
            P(I,J)=ROT(J,1)*PR(1)+ROT(J,2)*PR(2)+ROT(J,3)*PR(3)
            V(I,J)=ROT(J,1)*VR(1)+ROT(J,2)*VR(2)+ROT(J,3)*VR(3)
  130     CONTINUE
  140   CONTINUE
      ENDIF
 
C...Boost, typically from rest to momentum/energy=beta.
      IF(BEX**2+BEY**2+BEZ**2.GT.1D-20) THEN
        DBX=BEX
        DBY=BEY
        DBZ=BEZ
        DB=SQRT(DBX**2+DBY**2+DBZ**2)
        EPS1=1D0-1D-12
        IF(DB.GT.EPS1) THEN
C...Rescale boost vector if too close to unity.
          CALL PYERRM(3,'(PYROBO:) boost vector too large')
          DBX=DBX*(EPS1/DB)
          DBY=DBY*(EPS1/DB)
          DBZ=DBZ*(EPS1/DB)
          DB=EPS1
        ENDIF
        DGA=1D0/SQRT(1D0-DB**2)
        DO 160 I=IMIN,IMAX
          IF(K(I,1).LE.0) GOTO 160
          DO 150 J=1,4
            DP(J)=P(I,J)
            DV(J)=V(I,J)
  150     CONTINUE
          DBP=DBX*DP(1)+DBY*DP(2)+DBZ*DP(3)
          DGABP=DGA*(DGA*DBP/(1D0+DGA)+DP(4))
          P(I,1)=DP(1)+DGABP*DBX
          P(I,2)=DP(2)+DGABP*DBY
          P(I,3)=DP(3)+DGABP*DBZ
          P(I,4)=DGA*(DP(4)+DBP)
          DBV=DBX*DV(1)+DBY*DV(2)+DBZ*DV(3)
          DGABV=DGA*(DGA*DBV/(1D0+DGA)+DV(4))
          V(I,1)=DV(1)+DGABV*DBX
          V(I,2)=DV(2)+DGABV*DBY
          V(I,3)=DV(3)+DGABV*DBZ
          V(I,4)=DGA*(DV(4)+DBV)
  160   CONTINUE
      ENDIF
 
      RETURN
      END

C*********************************************************************
 
C...PYMASS
C...Gives the mass of a particle/parton.
 
      FUNCTION PYMASS(KF)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYDAT1/,/PYDAT2/
 
C...Reset variables. Compressed code. Special case for popcorn diquarks.
      PYMASS=0D0
      KFA=IABS(KF)
      KC=PYCOMP(KF)
      IF(KC.EQ.0) THEN
        MSTJ(93)=0
        RETURN
      ENDIF
 
C...Guarantee use of constituent masses for internal checks.
      IF((MSTJ(93).EQ.1.OR.MSTJ(93).EQ.2).AND.
     &(KFA.LE.10.OR.MOD(KFA/10,10).EQ.0)) THEN
        PARF(106)=PMAS(6,1)
        PARF(107)=PMAS(7,1)
        PARF(108)=PMAS(8,1)
        IF(KFA.LE.10) THEN
          PYMASS=PARF(100+KFA)
          IF(MSTJ(93).EQ.2) PYMASS=MAX(0D0,PYMASS-PARF(121))
        ELSEIF(MSTJ(93).EQ.1) THEN
          PYMASS=PARF(100+MOD(KFA/1000,10))+PARF(100+MOD(KFA/100,10))
        ELSE
          PYMASS=MAX(0D0,PMAS(KC,1)-PARF(122)-2D0*PARF(112)/3D0)
        ENDIF
 
C...Other masses can be read directly off table.
      ELSE
        PYMASS=PMAS(KC,1)
      ENDIF
 
C...Optional mass broadening according to truncated Breit-Wigner
C...(either in m or in m^2).
      IF(MSTJ(24).GE.1.AND.PMAS(KC,2).GT.1D-4) THEN
        IF(MSTJ(24).EQ.1.OR.(MSTJ(24).EQ.2.AND.KFA.GT.100)) THEN
          PYMASS=PYMASS+0.5D0*PMAS(KC,2)*TAN((2D0*PYR(0)-1D0)*
     &    ATAN(2D0*PMAS(KC,3)/PMAS(KC,2)))
        ELSE
          PM0=PYMASS
          PMLOW=ATAN((MAX(0D0,PM0-PMAS(KC,3))**2-PM0**2)/
     &    (PM0*PMAS(KC,2)))
          PMUPP=ATAN(((PM0+PMAS(KC,3))**2-PM0**2)/(PM0*PMAS(KC,2)))
          PYMASS=SQRT(MAX(0D0,PM0**2+PM0*PMAS(KC,2)*TAN(PMLOW+
     &    (PMUPP-PMLOW)*PYR(0))))
        ENDIF
      ENDIF
      MSTJ(93)=0
 
      RETURN
      END

C********************************************************************
 
C...PYKFDI
C...Generates a new flavour pair and combines off a hadron
 
      SUBROUTINE PYKFDI(KFL1,KFL2,KFL3,KF)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYDAT1/,/PYDAT2/
C...Local arrays.
      DIMENSION PD(7)
 
      IF(MSTU(123).EQ.0.AND.MSTJ(12).GT.0)  CALL PYKFIN
 
C...Default flavour values. Input consistency checks.
      KF1A=IABS(KFL1)
      KF2A=IABS(KFL2)
      KFL3=0
      KF=0
      IF(KF1A.EQ.0) RETURN
      IF(KF2A.NE.0)THEN
        IF(KF1A.LE.10.AND.KF2A.LE.10.AND.KFL1*KFL2.GT.0) RETURN
        IF(KF1A.GT.10.AND.KF2A.GT.10) RETURN
        IF((KF1A.GT.10.OR.KF2A.GT.10).AND.KFL1*KFL2.LT.0) RETURN
      ENDIF
 
C...Check if tabulated flavour probabilities are to be used.
      IF(MSTJ(15).EQ.1) THEN
        IF(MSTJ(12).GE.5)  CALL PYERRM(29,
     &        '(PYKFDI:) Sorry, option MSTJ(15)=1 not available' //
     &        ' together with MSTJ(12)>=5 modification')
        KTAB1=-1
        IF(KF1A.GE.1.AND.KF1A.LE.6) KTAB1=KF1A
        KFL1A=MOD(KF1A/1000,10)
        KFL1B=MOD(KF1A/100,10)
        KFL1S=MOD(KF1A,10)
        IF(KFL1A.GE.1.AND.KFL1A.LE.4.AND.KFL1B.GE.1.AND.KFL1B.LE.4)
     &  KTAB1=6+KFL1A*(KFL1A-2)+2*KFL1B+(KFL1S-1)/2
        IF(KFL1A.GE.1.AND.KFL1A.LE.4.AND.KFL1A.EQ.KFL1B) KTAB1=KTAB1-1
        IF(KF1A.GE.1.AND.KF1A.LE.6) KFL1A=KF1A
        KTAB2=0
        IF(KF2A.NE.0) THEN
          KTAB2=-1
          IF(KF2A.GE.1.AND.KF2A.LE.6) KTAB2=KF2A
          KFL2A=MOD(KF2A/1000,10)
          KFL2B=MOD(KF2A/100,10)
          KFL2S=MOD(KF2A,10)
          IF(KFL2A.GE.1.AND.KFL2A.LE.4.AND.KFL2B.GE.1.AND.KFL2B.LE.4)
     &    KTAB2=6+KFL2A*(KFL2A-2)+2*KFL2B+(KFL2S-1)/2
          IF(KFL2A.GE.1.AND.KFL2A.LE.4.AND.KFL2A.EQ.KFL2B) KTAB2=KTAB2-1
        ENDIF
        IF(KTAB1.GE.0.AND.KTAB2.GE.0) GOTO 140
      ENDIF
 
C.. Recognize rank 0 diquark case
  100 IRANK=1
      KFDIQ=MAX(KF1A,KF2A)
      IF(KFDIQ.GT.10.AND.KFDIQ.LT.10000) IRANK=0
 
C.. Join two flavours to meson or baryon. Test for popcorn.
      IF(KF2A.GT.0)THEN
        MBARY=0
        IF(KFDIQ.GT.10) THEN
          IF(IRANK.EQ.0.AND.MSTJ(12).LT.5)
     &         CALL PYNMES(KFDIQ)
          IF(MSTU(121).NE.0) RETURN
          MBARY=2
        ENDIF
        KFQOLD=KF1A
        KFQVER=KF2A
        GOTO 130
      ENDIF
 
C.. Separate incoming flavours, curtain flavour consistency check
      KFIN=KFL1
      KFQOLD=KF1A
      KFQPOP=KF1A/10000
      IF(KF1A.GT.10)THEN
         KFIN=-KFL1
         KFL1A=MOD(KF1A/1000,10)
         KFL1B=MOD(KF1A/100,10)
         IF(IRANK.EQ.0)THEN
            QAWT=1D0
            IF(KFL1A.GE.3) QAWT=PARF(136+KFL1A/4)
            IF(KFL1B.GE.3) QAWT=QAWT/PARF(136+KFL1B/4)
            KFQPOP=KFL1A+(KFL1B-KFL1A)*INT(1D0/(QAWT+1D0)+PYR(0))
         ENDIF
         IF(KFQPOP.NE.KFL1B.AND.KFQPOP.NE.KFL1A) RETURN
         KFQOLD=KFL1A+KFL1B-KFQPOP
      ENDIF
 
C...Meson/baryon choice. Set number of mesons if starting a popcorn
C...system.
  110 MBARY=0
      IF(KF1A.LE.10.AND.MSTJ(12).GT.0)THEN
         IF(MSTU(121).EQ.-1.OR.(1D0+PARJ(1))*PYR(0).GT.1D0)THEN
            MBARY=1
            CALL PYNMES(0)
         ENDIF
      ELSEIF(KF1A.GT.10)THEN
         MBARY=2
         IF(IRANK.EQ.0) CALL PYNMES(KF1A)
         IF(MSTU(121).GT.0) MBARY=-1
      ENDIF
 
C..x->H+q: Choose single vertex quark. Jump to form hadron.
      IF(MBARY.EQ.0.OR.MBARY.EQ.2)THEN
         KFQVER=1+INT((2D0+PARJ(2))*PYR(0))
         KFL3=ISIGN(KFQVER,-KFIN)
         GOTO 130
      ENDIF
 
C..x->H+qq: (IDW=proper PARF position for diquark weights)
      IDW=160
C..   q->B+qq: Get curtain quark, different weights for q->B+B and
C..   q->B+M+...
      IF(MBARY.EQ.1)THEN
         IF(MSTU(121).EQ.0) IDW=150
         SQWT=PARF(IDW+1)
         IF(MSTU(121).GT.0) SQWT=SQWT*PARF(135)*PARF(138)**MSTU(121)
         KFQPOP=1+INT((2D0+SQWT)*PYR(0))
C..   Shift to s-curtain parameters if needed
         IF(KFQPOP.GE.3.AND.MSTJ(12).GE.5)THEN
            PARF(194)=PARF(138)*PARF(139)
            PARF(193)=PARJ(8)+PARJ(9)
         ENDIF
      ENDIF
 
C.. x->H+qq: Get vertex quark
      IF(MBARY.EQ.-1.AND.MSTJ(12).GE.5)THEN
         IDW=MSTU(122)
         MSTU(121)=MSTU(121)-1
         IF(IDW.EQ.170) THEN
            IF(MSTU(121).EQ.0)THEN
               IPOS=3*MIN(KFQPOP-1,2)+MIN(KFQOLD-1,2)
            ELSE
               IPOS=3*3+3*MAX(0,MIN(KFQPOP-2,1))+MIN(KFQOLD-1,2)
            ENDIF
         ELSE
            IF(MSTU(121).EQ.0)THEN
               IPOS=3*5+5*MIN(KFQPOP-1,3)+MIN(KFQOLD-1,4)
            ELSE
               IPOS=3*5+5*4+MIN(KFQOLD-1,4)
            ENDIF
         ENDIF
         IPOS=200+30*IPOS+1
 
         IMES=-1
         RMES=PYR(0)*PARF(194)
  120    IMES=IMES+1
         RMES=RMES-PARF(IPOS+IMES)
         IF(IMES.EQ.30) THEN
            MSTU(121)=-1
            KF=-111
            RETURN
         ENDIF
         IF(RMES.GT.0D0) GOTO 120
         KMUL=IMES/5
         KFJ=2*KMUL+1
         IF(KMUL.EQ.2) KFJ=10003
         IF(KMUL.EQ.3) KFJ=10001
         IF(KMUL.EQ.4) KFJ=20003
         IF(KMUL.EQ.5) KFJ=5
         IDIAG=0
         KFQVER=MOD(IMES,5)+1
         IF(KFQVER.GE.KFQOLD) KFQVER=KFQVER+1
         IF(KFQVER.GT.3)THEN
            IDIAG=KFQVER-3
            KFQVER=KFQOLD
         ENDIF
      ELSE
         IF(MBARY.EQ.-1) IDW=170
         SQWT=PARF(IDW+2)
         IF(KFQPOP.EQ.3) SQWT=PARF(IDW+3)
         IF(KFQPOP.GT.3) SQWT=PARF(IDW+3)*(1D0/PARF(IDW+5)+1D0)/2D0
         KFQVER=MIN(3,1+INT((2D0+SQWT)*PYR(0)))
         IF(KFQPOP.LT.3.AND.KFQVER.LT.3)THEN
            KFQVER=KFQPOP
            IF(PYR(0).GT.PARF(IDW+4)) KFQVER=3-KFQPOP
         ENDIF
      ENDIF
 
C..x->H+qq: form outgoing diquark with KFQPOP flag at 10000-pos
      KFLDS=3
      IF(KFQPOP.NE.KFQVER)THEN
         SWT=PARF(IDW+7)
         IF(KFQVER.EQ.3) SWT=PARF(IDW+6)
         IF(KFQPOP.GE.3) SWT=PARF(IDW+5)
         IF((1D0+SWT)*PYR(0).LT.1D0) KFLDS=1
      ENDIF
      KFDIQ=900*MAX(KFQVER,KFQPOP)+100*(KFQVER+KFQPOP)+KFLDS
     &      +10000*KFQPOP
      KFL3=ISIGN(KFDIQ,KFIN)
 
C..x->M+y: flavour for meson.
  130 IF(MBARY.LE.0)THEN
        KFLA=MAX(KFQOLD,KFQVER)
        KFLB=MIN(KFQOLD,KFQVER)
        KFS=ISIGN(1,KFL1)
        IF(KFLA.NE.KFQOLD) KFS=-KFS
C... Form meson, with spin and flavour mixing for diagonal states.
        IF(MBARY.EQ.-1.AND.MSTJ(12).GE.5)THEN
           IF(IDIAG.GT.0) KF=110*IDIAG+KFJ
           IF(IDIAG.EQ.0) KF=(100*KFLA+10*KFLB+KFJ)*KFS*(-1)**KFLA
           RETURN
        ENDIF
        IF(KFLA.LE.2) KMUL=INT(PARJ(11)+PYR(0))
        IF(KFLA.EQ.3) KMUL=INT(PARJ(12)+PYR(0))
        IF(KFLA.GE.4) KMUL=INT(PARJ(13)+PYR(0))
        IF(KMUL.EQ.0.AND.PARJ(14).GT.0D0)THEN
          IF(PYR(0).LT.PARJ(14)) KMUL=2
        ELSEIF(KMUL.EQ.1.AND.PARJ(15)+PARJ(16)+PARJ(17).GT.0D0)THEN
          RMUL=PYR(0)
          IF(RMUL.LT.PARJ(15)) KMUL=3
          IF(KMUL.EQ.1.AND.RMUL.LT.PARJ(15)+PARJ(16)) KMUL=4
          IF(KMUL.EQ.1.AND.RMUL.LT.PARJ(15)+PARJ(16)+PARJ(17)) KMUL=5
        ENDIF
        KFLS=3
        IF(KMUL.EQ.0.OR.KMUL.EQ.3) KFLS=1
        IF(KMUL.EQ.5) KFLS=5
        IF(KFLA.NE.KFLB)THEN
          KF=(100*KFLA+10*KFLB+KFLS)*KFS*(-1)**KFLA
        ELSE
          RMIX=PYR(0)
          IMIX=2*KFLA+10*KMUL
          IF(KFLA.LE.3) KF=110*(1+INT(RMIX+PARF(IMIX-1))+
     &    INT(RMIX+PARF(IMIX)))+KFLS
          IF(KFLA.GE.4) KF=110*KFLA+KFLS
        ENDIF
        IF(KMUL.EQ.2.OR.KMUL.EQ.3) KF=KF+ISIGN(10000,KF)
        IF(KMUL.EQ.4) KF=KF+ISIGN(20000,KF)
 
C..Optional extra suppression of eta and eta'.
C..Allow shift to qq->B+q in old version (set IRANK to 0)
        IF(KF.EQ.221.OR.KF.EQ.331)THEN
           IF(PYR(0).GT.PARJ(25+KF/300))THEN
              IF(KF2A.GT.0) GOTO 130
              IF(MSTJ(12).LT.4) IRANK=0
              GOTO 110
           ENDIF
        ENDIF
        MSTU(121)=0
 
C.. x->B+y: Flavour for baryon
      ELSE
        KFLA=KFQVER
        IF(KF1A.LE.10) KFLA=KFQOLD
        KFLB=MOD(KFDIQ/1000,10)
        KFLC=MOD(KFDIQ/100,10)
        KFLDS=MOD(KFDIQ,10)
        KFLD=MAX(KFLA,KFLB,KFLC)
        KFLF=MIN(KFLA,KFLB,KFLC)
        KFLE=KFLA+KFLB+KFLC-KFLD-KFLF
 
C...  SU(6) factors for formation of baryon.
        KBARY=3
        KDMAX=5
        KFLG=KFLB
        IF(KFLB.NE.KFLC)THEN
           KBARY=2*KFLDS-1
           KDMAX=1+KFLDS/2
           IF(KFLB.GT.2) KDMAX=KDMAX+2
        ENDIF
        IF(KFLA.NE.KFLB.AND.KFLA.NE.KFLC)THEN
           KBARY=KBARY+1
           KFLG=KFLA
        ENDIF
 
        SU6MAX=PARF(140+KDMAX)
        SU6DEC=PARJ(18)
        SU6S  =PARF(146)
        IF(MSTJ(12).GE.5.AND.IRANK.EQ.0) THEN
           SU6MAX=1D0
           SU6DEC=1D0
           SU6S  =1D0
        ENDIF
        SU6OCT=PARF(60+KBARY)
        IF(KFLG.GT.MAX(KFLA+KFLB-KFLG,2))THEN
           SU6OCT=SU6OCT*4*SU6S/(3*SU6S+1)
           IF(KBARY.EQ.2) SU6OCT=PARF(60+KBARY)*4/(3*SU6S+1)
        ELSE
           IF(KBARY.EQ.6) SU6OCT=SU6OCT*(3+SU6S)/(3*SU6S+1)
        ENDIF
        SU6WT=SU6OCT+SU6DEC*PARF(70+KBARY)
 
C..   SU(6) test. Old options enforce new baryon if q->B+qq is rejected.
        IF(SU6WT.LT.PYR(0)*SU6MAX.AND.KF2A.EQ.0)THEN
           MSTU(121)=0
           IF(MSTJ(12).LE.2.AND.MBARY.EQ.1) MSTU(121)=-1
           GOTO 110
        ENDIF
 
C.. Form baryon. Distinguish Lambda- and Sigmalike baryons.
        KSIG=1
        KFLS=2
        IF(SU6WT*PYR(0).GT.SU6OCT) KFLS=4
        IF(KFLS.EQ.2.AND.KFLD.GT.KFLE.AND.KFLE.GT.KFLF)THEN
          KSIG=KFLDS/3
          IF(KFLA.NE.KFLD) KSIG=INT(3*SU6S/(3*SU6S+KFLDS**2)+PYR(0))
        ENDIF
        KF=ISIGN(1000*KFLD+100*KFLE+10*KFLF+KFLS,KFL1)
        IF(KSIG.EQ.0) KF=ISIGN(1000*KFLD+100*KFLF+10*KFLE+KFLS,KFL1)
      ENDIF
      RETURN
 
C...Use tabulated probabilities to select new flavour and hadron.
  140 IF(KTAB2.EQ.0.AND.MSTJ(12).LE.0) THEN
        KT3L=1
        KT3U=6
      ELSEIF(KTAB2.EQ.0.AND.KTAB1.GE.7.AND.MSTJ(12).LE.1) THEN
        KT3L=1
        KT3U=6
      ELSEIF(KTAB2.EQ.0) THEN
        KT3L=1
        KT3U=22
      ELSE
        KT3L=KTAB2
        KT3U=KTAB2
      ENDIF
      RFL=0D0
      DO 160 KTS=0,2
        DO 150 KT3=KT3L,KT3U
          RFL=RFL+PARF(120+80*KTAB1+25*KTS+KT3)
  150   CONTINUE
  160 CONTINUE
      RFL=PYR(0)*RFL
      DO 180 KTS=0,2
        KTABS=KTS
        DO 170 KT3=KT3L,KT3U
          KTAB3=KT3
          RFL=RFL-PARF(120+80*KTAB1+25*KTS+KT3)
          IF(RFL.LE.0D0) GOTO 190
  170   CONTINUE
  180 CONTINUE
  190 CONTINUE
 
C...Reconstruct flavour of produced quark/diquark.
      IF(KTAB3.LE.6) THEN
        KFL3A=KTAB3
        KFL3B=0
        KFL3=ISIGN(KFL3A,KFL1*(2*KTAB1-13))
      ELSE
        KFL3A=1
        IF(KTAB3.GE.8) KFL3A=2
        IF(KTAB3.GE.11) KFL3A=3
        IF(KTAB3.GE.16) KFL3A=4
        KFL3B=(KTAB3-6-KFL3A*(KFL3A-2))/2
        KFL3=1000*KFL3A+100*KFL3B+1
        IF(KFL3A.EQ.KFL3B.OR.KTAB3.NE.6+KFL3A*(KFL3A-2)+2*KFL3B) KFL3=
     &  KFL3+2
        KFL3=ISIGN(KFL3,KFL1*(13-2*KTAB1))
      ENDIF
 
C...Reconstruct meson code.
      IF(KFL3A.EQ.KFL1A.AND.KFL3B.EQ.KFL1B.AND.(KFL3A.LE.3.OR.
     &KFL3B.NE.0)) THEN
        RFL=PYR(0)*(PARF(143+80*KTAB1+25*KTABS)+PARF(144+80*KTAB1+
     &  25*KTABS)+PARF(145+80*KTAB1+25*KTABS))
        KF=110+2*KTABS+1
        IF(RFL.GT.PARF(143+80*KTAB1+25*KTABS)) KF=220+2*KTABS+1
        IF(RFL.GT.PARF(143+80*KTAB1+25*KTABS)+PARF(144+80*KTAB1+
     &  25*KTABS)) KF=330+2*KTABS+1
      ELSEIF(KTAB1.LE.6.AND.KTAB3.LE.6) THEN
        KFLA=MAX(KTAB1,KTAB3)
        KFLB=MIN(KTAB1,KTAB3)
        KFS=ISIGN(1,KFL1)
        IF(KFLA.NE.KF1A) KFS=-KFS
        KF=(100*KFLA+10*KFLB+2*KTABS+1)*KFS*(-1)**KFLA
      ELSEIF(KTAB1.GE.7.AND.KTAB3.GE.7) THEN
        KFS=ISIGN(1,KFL1)
        IF(KFL1A.EQ.KFL3A) THEN
          KFLA=MAX(KFL1B,KFL3B)
          KFLB=MIN(KFL1B,KFL3B)
          IF(KFLA.NE.KFL1B) KFS=-KFS
        ELSEIF(KFL1A.EQ.KFL3B) THEN
          KFLA=KFL3A
          KFLB=KFL1B
          KFS=-KFS
        ELSEIF(KFL1B.EQ.KFL3A) THEN
          KFLA=KFL1A
          KFLB=KFL3B
        ELSEIF(KFL1B.EQ.KFL3B) THEN
          KFLA=MAX(KFL1A,KFL3A)
          KFLB=MIN(KFL1A,KFL3A)
          IF(KFLA.NE.KFL1A) KFS=-KFS
        ELSE
          CALL PYERRM(2,'(PYKFDI:) no matching flavours for qq -> qq')
          GOTO 100
        ENDIF
        KF=(100*KFLA+10*KFLB+2*KTABS+1)*KFS*(-1)**KFLA
 
C...Reconstruct baryon code.
      ELSE
        IF(KTAB1.GE.7) THEN
          KFLA=KFL3A
          KFLB=KFL1A
          KFLC=KFL1B
        ELSE
          KFLA=KFL1A
          KFLB=KFL3A
          KFLC=KFL3B
        ENDIF
        KFLD=MAX(KFLA,KFLB,KFLC)
        KFLF=MIN(KFLA,KFLB,KFLC)
        KFLE=KFLA+KFLB+KFLC-KFLD-KFLF
        IF(KTABS.EQ.0) KF=ISIGN(1000*KFLD+100*KFLF+10*KFLE+2,KFL1)
        IF(KTABS.GE.1) KF=ISIGN(1000*KFLD+100*KFLE+10*KFLF+2*KTABS,KFL1)
      ENDIF
 
C...Check that constructed flavour code is an allowed one.
      IF(KFL2.NE.0) KFL3=0
      KC=PYCOMP(KF)
      IF(KC.EQ.0) THEN
        CALL PYERRM(2,'(PYKFDI:) user-defined flavour probabilities '//
     &  'failed')
        GOTO 100
      ENDIF
 
      RETURN
      END


C*********************************************************************
 
C...PYERRM
C...Informs user of errors in program execution.
 
      SUBROUTINE PYERRM(MERR,CHMESS)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /PYJETS/,/PYDAT1/
C...Local character variable.
      CHARACTER CHMESS*(*)
 
C...Write first few warnings, then be silent.
      IF(MERR.LE.10) THEN
        MSTU(27)=MSTU(27)+1
        MSTU(28)=MERR
        IF(MSTU(25).EQ.1.AND.MSTU(27).LE.MSTU(26)) WRITE(MSTU(11),5000)
     &  MERR,MSTU(31),CHMESS
 
C...Write first few errors, then be silent or stop program.
      ELSEIF(MERR.LE.20) THEN
        MSTU(23)=MSTU(23)+1
        MSTU(24)=MERR-10
        IF(MSTU(21).GE.1.AND.MSTU(23).LE.MSTU(22)) WRITE(MSTU(11),5100)
     &  MERR-10,MSTU(31),CHMESS
        IF(MSTU(21).GE.2.AND.MSTU(23).GT.MSTU(22)) THEN
          WRITE(MSTU(11),5100) MERR-10,MSTU(31),CHMESS
          WRITE(MSTU(11),5200)
          IF(MERR.NE.17) CALL PYLIST(2)
          STOP
        ENDIF
 
C...Stop program in case of irreparable error.
      ELSE
        WRITE(MSTU(11),5300) MERR-20,MSTU(31),CHMESS
        STOP
      ENDIF
 
C...Formats for output.
 5000 FORMAT(/5X,'Advisory warning type',I2,' given after',I9,
     &' PYEXEC calls:'/5X,A)
 5100 FORMAT(/5X,'Error type',I2,' has occured after',I9,
     &' PYEXEC calls:'/5X,A)
 5200 FORMAT(5X,'Execution will be stopped after listing of last ',
     &'event!')
 5300 FORMAT(/5X,'Fatal error type',I2,' has occured after',I9,
     &' PYEXEC calls:'/5X,A/5X,'Execution will now be stopped!')
 
      RETURN
      END

C*********************************************************************
 
C...PYPTDI
C...Generates transverse momentum according to a Gaussian.
 
      SUBROUTINE PYPTDI(KFL,PX,PY)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /PYDAT1/
 
C...Generate p_T and azimuthal angle, gives p_x and p_y.
      KFLA=IABS(KFL)
      PT=PARJ(21)*SQRT(-LOG(MAX(1D-10,PYR(0))))
      IF(PARJ(23).GT.PYR(0)) PT=PARJ(24)*PT
      IF(MSTJ(91).EQ.1) PT=PARJ(22)*PT
      IF(KFLA.EQ.0.AND.MSTJ(13).LE.0) PT=0D0
      PHI=PARU(2)*PYR(0)
      PX=PT*COS(PHI)
      PY=PT*SIN(PHI)
 
      RETURN
      END

C*********************************************************************
 
C...PYZDIS
C...Generates the longitudinal splitting variable z.
 
      SUBROUTINE PYZDIS(KFL1,KFL2,PR,Z)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYDAT1/,/PYDAT2/
 
C...Check if heavy flavour fragmentation.
      KFLA=IABS(KFL1)
      KFLB=IABS(KFL2)
      KFLH=KFLA
      IF(KFLA.GE.10) KFLH=MOD(KFLA/1000,10)
 
C...Lund symmetric scaling function: determine parameters of shape.
      IF(MSTJ(11).EQ.1.OR.(MSTJ(11).EQ.3.AND.KFLH.LE.3).OR.
     &MSTJ(11).GE.4) THEN
        FA=PARJ(41)
        IF(MSTJ(91).EQ.1) FA=PARJ(43)
        IF(KFLB.GE.10) FA=FA+PARJ(45)
        FBB=PARJ(42)
        IF(MSTJ(91).EQ.1) FBB=PARJ(44)
        FB=FBB*PR
        FC=1D0
        IF(KFLA.GE.10) FC=FC-PARJ(45)
        IF(KFLB.GE.10) FC=FC+PARJ(45)
        IF(MSTJ(11).GE.4.AND.KFLH.GE.4.AND.KFLH.LE.5) THEN
          FRED=PARJ(46)
          IF(MSTJ(11).EQ.5.AND.KFLH.EQ.5) FRED=PARJ(47)
          FC=FC+FRED*FBB*PARF(100+KFLH)**2
        ELSEIF(MSTJ(11).GE.4.AND.KFLH.GE.6.AND.KFLH.LE.8) THEN
          FRED=PARJ(46)
          IF(MSTJ(11).EQ.5) FRED=PARJ(48)
          FC=FC+FRED*FBB*PMAS(KFLH,1)**2
        ENDIF
        MC=1
        IF(ABS(FC-1D0).GT.0.01D0) MC=2
 
C...Determine position of maximum. Special cases for a = 0 or a = c.
        IF(FA.LT.0.02D0) THEN
          MA=1
          ZMAX=1D0
          IF(FC.GT.FB) ZMAX=FB/FC
        ELSEIF(ABS(FC-FA).LT.0.01D0) THEN
          MA=2
          ZMAX=FB/(FB+FC)
        ELSE
          MA=3
          ZMAX=0.5D0*(FB+FC-SQRT((FB-FC)**2+4D0*FA*FB))/(FC-FA)
          IF(ZMAX.GT.0.9999D0.AND.FB.GT.100D0) ZMAX=MIN(ZMAX,1D0-FA/FB)
        ENDIF
 
C...Subdivide z range if distribution very peaked near endpoint.
        MMAX=2
        IF(ZMAX.LT.0.1D0) THEN
          MMAX=1
          ZDIV=2.75D0*ZMAX
          IF(MC.EQ.1) THEN
            FINT=1D0-LOG(ZDIV)
          ELSE
            ZDIVC=ZDIV**(1D0-FC)
            FINT=1D0+(1D0-1D0/ZDIVC)/(FC-1D0)
          ENDIF
        ELSEIF(ZMAX.GT.0.85D0.AND.FB.GT.1D0) THEN
          MMAX=3
          FSCB=SQRT(4D0+(FC/FB)**2)
          ZDIV=FSCB-1D0/ZMAX-(FC/FB)*LOG(ZMAX*0.5D0*(FSCB+FC/FB))
          IF(MA.GE.2) ZDIV=ZDIV+(FA/FB)*LOG(1D0-ZMAX)
          ZDIV=MIN(ZMAX,MAX(0D0,ZDIV))
          FINT=1D0+FB*(1D0-ZDIV)
        ENDIF
 
C...Choice of z, preweighted for peaks at low or high z.
  100   Z=PYR(0)
        FPRE=1D0
        IF(MMAX.EQ.1) THEN
          IF(FINT*PYR(0).LE.1D0) THEN
            Z=ZDIV*Z
          ELSEIF(MC.EQ.1) THEN
            Z=ZDIV**Z
            FPRE=ZDIV/Z
          ELSE
            Z=(ZDIVC+Z*(1D0-ZDIVC))**(1D0/(1D0-FC))
            FPRE=(ZDIV/Z)**FC
          ENDIF
        ELSEIF(MMAX.EQ.3) THEN
          IF(FINT*PYR(0).LE.1D0) THEN
            Z=ZDIV+LOG(Z)/FB
            FPRE=EXP(FB*(Z-ZDIV))
          ELSE
            Z=ZDIV+Z*(1D0-ZDIV)
          ENDIF
        ENDIF
 
C...Weighting according to correct formula.
        IF(Z.LE.0D0.OR.Z.GE.1D0) GOTO 100
        FEXP=FC*LOG(ZMAX/Z)+FB*(1D0/ZMAX-1D0/Z)
        IF(MA.GE.2) FEXP=FEXP+FA*LOG((1D0-Z)/(1D0-ZMAX))
        FVAL=EXP(MAX(-50D0,MIN(50D0,FEXP)))
        IF(FVAL.LT.PYR(0)*FPRE) GOTO 100
 
C...Generate z according to Field-Feynman, SLAC, (1-z)**c OR z**c.
      ELSE
        FC=PARJ(50+MAX(1,KFLH))
        IF(MSTJ(91).EQ.1) FC=PARJ(59)
  110   Z=PYR(0)
        IF(FC.GE.0D0.AND.FC.LE.1D0) THEN
          IF(FC.GT.PYR(0)) Z=1D0-Z**(1D0/3D0)
        ELSEIF(FC.GT.-1.AND.FC.LT.0D0) THEN
          IF(-4D0*FC*Z*(1D0-Z)**2.LT.PYR(0)*((1D0-Z)**2-FC*Z)**2)
     &    GOTO 110
        ELSE
          IF(FC.GT.0D0) Z=1D0-Z**(1D0/FC)
          IF(FC.LT.0D0) Z=Z**(-1D0/FC)
        ENDIF
      ENDIF
 
      RETURN
      END

C*********************************************************************
 
C...PYNMES
C...Generates number of popcorn mesons and stores some relevant
C...parameters.
 
      SUBROUTINE PYNMES(KFDIQ)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYDAT1/,/PYDAT2/
 
      MSTU(121)=0
      IF(MSTJ(12).LT.2) RETURN
 
C..Old version: Get 1 or 0 popcorn mesons
      IF(MSTJ(12).LT.5)THEN
         POPWT=PARF(131)
         IF(KFDIQ.NE.0) THEN
            KFDIQA=IABS(KFDIQ)
            KFA=MOD(KFDIQA/1000,10)
            KFB=MOD(KFDIQA/100,10)
            KFS=MOD(KFDIQA,10)
            POPWT=PARF(132)
            IF(KFA.EQ.3) POPWT=PARF(133)
            IF(KFB.EQ.3) POPWT=PARF(134)
            IF(KFS.EQ.1) POPWT=POPWT*SQRT(PARJ(4))
         ENDIF
         MSTU(121)=INT(POPWT/(1D0+POPWT)+PYR(0))
         RETURN
      ENDIF
 
C..New version: Store popcorn- or rank 0 diquark parameters
      MSTU(122)=170
      PARF(193)=PARJ(8)
      PARF(194)=PARF(139)
      IF(KFDIQ.NE.0) THEN
         MSTU(122)=180
         PARF(193)=PARJ(10)
         PARF(194)=PARF(140)
      ENDIF
      IF(PARF(194).LT.1D-5.OR.PARF(194).GT.1D0-1D-5) THEN
         IF(PARF(194).GT.1D0-1D-5) CALL PYERRM(9,
     &        '(PYNMES:) Neglecting too large popcorn possibility')
         RETURN
      ENDIF
 
C..New version: Get number of popcorn mesons
  100 RTST=PYR(0)
      MSTU(121)=-1
  110 MSTU(121)=MSTU(121)+1
      RTST=RTST/PARF(194)
      IF(RTST.LT.1D0) GOTO 110
      IF(KFDIQ.EQ.0.AND.PYR(0)*(2D0+PARF(135)).GT.
     &     (2D0+PARF(135)*PARF(138)**MSTU(121))) GOTO 100
      RETURN
      END

C*********************************************************************
 
C...PYKFIN
C...Precalculates a set of diquark and popcorn weights.
C.. (Results stored in order SU0,US0,SS1,UU1,SU1,US1,UD1)
 
      SUBROUTINE PYKFIN
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYDAT1/,/PYDAT2/
 
      DIMENSION SU6(12),SU6M(7)
 
      MSTU(123)=1
C..Curtain tunneling factor T(D,q)/T(ud0,u).
      IF(MSTJ(12).GE.5) THEN
         PMUD0=PYMASS(2101)
         PMUD1=PYMASS(2103)-PMUD0
         PMUS0=PYMASS(3201)-PMUD0
         PMUS1=PYMASS(3203)-PMUS0-PMUD0
         PMSS1=PYMASS(3303)-PMUS0-PMUD0
         PARF(151)=EXP(-(PARJ(9)+PARJ(8))*PMUS0-PARJ(9)*PARF(191))
         PARF(152)=EXP(-PARJ(8)*PMUS0)
         PARF(153)=EXP(-(PARJ(9)+PARJ(8))*PMSS1)*PARF(151)
         PARF(154)=EXP(-PARJ(8)*PMUD1)
         PARF(155)=EXP(-(PARJ(9)+PARJ(8))*PMUS1)*PARF(151)
         PARF(156)=EXP(-PARJ(8)*PMUS1)*PARF(152)
         PARF(157)=PARF(154)
      ELSE
         PAR2M=SQRT(PARJ(2))
         PAR3M=SQRT(PARJ(3))
         PAR4M=SQRT(PARJ(4))
         PARF(151)=PAR2M*PAR3M
         PARF(152)=PAR3M
         PARF(153)=PAR2M*PARJ(3)*PAR4M
         PARF(154)=PAR4M
         PARF(155)=PAR4M*PARF(151)
         PARF(156)=PAR4M*PARF(152)
         PARF(157)=PAR4M
      ENDIF
 
C.. Total tunneling factor tau(D,q)=T*vertex*spin.
      PARF(161)=PARF(151)
      PARF(162)=PARJ(2)*PARF(152)
      PARF(163)=PARJ(2)*6D0*PARF(153)
      PARF(164)=6D0*PARF(154)
      PARF(165)=3D0*PARF(155)
      PARF(166)=PARJ(2)*3D0*PARF(156)
      PARF(167)=3D0*PARF(157)
 
      DO 100 I=1,7
         PARF(150+I)=PARF(150+I)*PARF(160+I)
  100 CONTINUE
 
C..Modified SU(6) factors.
      PARF(146)=1D0
      IF(MSTJ(12).GE.5) PARF(146)=3D0*PARJ(18)/(2D0*PARJ(18)+1D0)
      IF(PARJ(18).LT.1D0-1D-5.AND.MSTJ(12).LT.5) CALL PYERRM(9,
     &     '(PYKFIN:) PARJ(18)<1 combined with 0<MSTJ(12)<5 option')
      DO 110 I=1,6
         SU6(I)=PARF(60+I)
         SU6(6+I)=SU6(I)*4*PARF(146)/(3*PARF(146)+1)
  110 CONTINUE
      SU6(8)=SU6(2)*4/(3*PARF(146)+1)
      SU6(6)=SU6(6)*(3+PARF(146))/(3*PARF(146)+1)
      DO 120 I=1,6
         SU6(I)=SU6(I)+PARJ(18)*PARF(70+I)
         SU6(6+I)=SU6(6+I)+PARJ(18)*PARF(70+I)
  120 CONTINUE
 
C..Total diquark quark*SU(6).
      PUD0=(2D0*SU6(1)+PARJ(2)*SU6(8))
      PARF(171)=(SU6(7)+SU6(2)+PARJ(2)*SU6(1))/PUD0
      PARF(172)=PARF(171)
      PARF(173)=(2D0*SU6(4)+PARJ(2)*SU6(3))/PUD0
      PARF(174)=(SU6(3)+SU6(4)+PARJ(2)*SU6(10))/PUD0
      PARF(175)=(SU6(11)+SU6(6)+PARJ(2)*SU6(5))/PUD0
      PARF(176)=PARF(175)
      PARF(177)=(2D0*SU6(5)+PARJ(2)*SU6(12))/PUD0
 
C..SU(6)max         q       q'     s,c,b
      SU6MUD =MAX(SU6(1) ,       SU6(8) )
      SU6M(7)=MAX(SU6(5) ,       SU6(12))
      SU6M(1)=MAX(SU6(7) ,SU6(2),SU6MUD )
      SU6M(4)=MAX(SU6(3) ,SU6(4),SU6(10))
      SU6M(5)=MAX(SU6(11),SU6(6),SU6M(7))
      SU6M(2)=SU6M(1)
      SU6M(3)=SU6M(4)
      SU6M(6)=SU6M(5)
 
      IF(MSTJ(12).GE.5)THEN
C..New version: tau for rank 0 diquark.
         PARF(181)=EXP(-PARJ(10)*PMUS0)
         PARF(182)=PARJ(2)*PARF(181)
         PARF(183)=6D0*PARJ(2)*EXP(-PARJ(10)*PMSS1)*PARF(181)
         PARF(184)=3D0*EXP(-PARJ(10)*PMUD1)
         PARF(185)=3D0*EXP(-PARJ(10)*PMUS1)*PARF(181)
         PARF(186)=PARJ(2)*PARF(185)
         PARF(187)=2D0*PARF(184)
 
C..New version: s/u curtain ratios.
         WU=1D0+PARF(167)+PARF(162)+PARF(166)+PARF(164)
         PARF(135)=(2D0*(PARF(161)+PARF(165))+PARF(163))/WU
         WU=1D0+PARF(187)+PARF(182)+PARF(186)+PARF(184)
         PARF(136)=(2D0*(PARF(181)+PARF(185))+PARF(183))/WU
         PARF(137)=(PARF(181)+PARF(185))*
     &        (2D0+PARF(183)/(2D0*PARF(185)))/WU
      ELSE
C..Old version: Shuffle PARJ(7) into tau
         PARF(162)=PARF(162)*PARJ(7)
         PARF(163)=PARF(163)*PARJ(7)
         PARF(166)=PARF(166)*PARJ(7)
 
C..Old version: s/u curtain ratios.
         WU=1D0+PARF(167)+PARF(162)+PARF(166)+PARF(164)
         PARF(135)=(2D0*(PARF(161)+PARF(165))+PARF(163))/WU
         PARF(136)=PARF(135)*PARJ(6)*PARF(161)/PARF(162)
         PARF(137)=(1D0+PARF(167))*(2D0+PARF(162))/WU
      ENDIF
 
C..Combine SU(6), SU(6)max, tau and T into proper products
      DO 140 I=1,7
         PARF(180+I)=PARF(180+I)*PARF(170+I)
         PARF(170+I)=PARF(170+I)*PARF(160+I)
         PARF(160+I)=PARF(160+I)*SU6M(I)/SU6MUD
         PARF(150+I)=PARF(150+I)*SU6M(I)/SU6MUD
  140 CONTINUE
 
C..Store SU(6)max, in order UD0,UD1,US0,US1,QQ1
      PARF(141)=SU6MUD
      PARF(142)=SU6M(7)
      PARF(143)=SU6M(1)
      PARF(144)=SU6M(5)
      PARF(145)=SU6M(3)
 
      IF(MSTJ(12).LT.5)THEN
C.. Old version: Resulting popcorn weights.
         PARF(138)=PARJ(6)
         WS=PARF(135)*PARF(138)
         WQ=WU*PARJ(5)/3D0
         PARF(132)=WQ*PARF(167)/PARF(157)
         PARF(133)=WQ*(PARF(166)/PARF(156)+WS*PARF(165)/PARF(155))/2D0
         PARF(134)=WQ*WS*PARF(163)/PARF(153)
         PARF(131)=WQ*((1D0+PARF(167))*(1D0+PARF(162)+WS*PARF(161))+
     &     PARF(164)+WS*PARF(163)/2D0)/
     &    ((1D0+PARF(157))*(1D0+2D0*PARF(152))+PARF(154)+PARF(153)/2D0)
      ELSE
C..New version: Store weights for popcorn mesons,
C..get prel. popcorn weights.
         DO 150 IPOS=201,1400
            PARF(IPOS)=0D0
  150    CONTINUE
         DO 160 I=138,140
            PARF(I)=0D0
  160    CONTINUE
         IPOS=200
         PARF(193)=PARJ(8)
         DO 240 MR=170,180,10
           IF(MR.EQ.180) PARF(193)=PARJ(10)
           SQWT=2D0*(PARF(MR+2)+PARF(MR+6))/(1D0+PARF(MR+7)+PARF(MR+4))
           QQWT=PARF(MR+4)/(1D0+PARF(MR+7)+PARF(MR+4))
           DO 230 NMES=0,1
             IF(NMES.EQ.1) SQWT=PARJ(2)
             DO 220 KFQPOP=1,4
               IF(MR.EQ.170.AND.KFQPOP.GT.3) GOTO 220
               IF(NMES.EQ.0.AND.KFQPOP.GE.3)THEN
                  SQWT=PARF(MR+3)/(PARF(MR+1)+PARF(MR+5))
                  QQWT=0.5D0
                  IF(MR.EQ.170) PARF(193)=PARJ(8)+PARJ(9)
                  IF(KFQPOP.EQ.4) SQWT=SQWT*(1D0/PARF(185)+1D0)/2D0
               ENDIF
               DO 210 KFQOLD =1,5
                  IF(MR.EQ.170.AND.KFQOLD.GT.3) GOTO 210
                  IF(MR*NMES.EQ.170.AND.KFQPOP.EQ.1) GOTO 210
                  IF(MR*NMES.EQ.180.AND.KFQPOP.NE.1) GOTO 210
                  WTTOT=0D0
                  WTFAIL=0D0
      DO 190 KMUL=0,5
         PJWT=PARJ(12+KMUL)
         IF(KMUL.EQ.0) PJWT=1D0-PARJ(14)
         IF(KMUL.EQ.1) PJWT=1D0-PARJ(15)-PARJ(16)-PARJ(17)
         IF(PJWT.LE.0D0) GOTO 190
         IF(PJWT.GT.1D0) PJWT=1D0
         IMES=5*KMUL
         IMIX=2*KFQOLD+10*KMUL
         KFJ=2*KMUL+1
         IF(KMUL.EQ.2) KFJ=10003
         IF(KMUL.EQ.3) KFJ=10001
         IF(KMUL.EQ.4) KFJ=20003
         IF(KMUL.EQ.5) KFJ=5
         DO 180 KFQVER =1,3
            KFLA=MAX(KFQOLD,KFQVER)
            KFLB=MIN(KFQOLD,KFQVER)
            SWT=PARJ(11+KFLA/3+KFLA/4)
            IF(KMUL.EQ.0.OR.KMUL.EQ.2) SWT=1D0-SWT
            SWT=SWT*PJWT
            QWT=SQWT/(2D0+SQWT)
            IF(KFQVER.LT.3)THEN
               IF(KFQVER.EQ.KFQPOP) QWT=(1D0-QWT)*QQWT
               IF(KFQVER.NE.KFQPOP) QWT=(1D0-QWT)*(1D0-QQWT)
            ENDIF
            IF(KFQVER.NE.KFQOLD)THEN
               IMES=IMES+1
               KFM=100*KFLA+10*KFLB+KFJ
               PMM=PMAS(PYCOMP(KFM),1)-PMAS(PYCOMP(KFM),3)
               PARF(IPOS+IMES)=QWT*SWT*EXP(-PARF(193)*PMM)
               WTTOT=WTTOT+PARF(IPOS+IMES)
            ELSE
               DO 170 ID=3,5
                  IF(ID.EQ.3) DWT=1D0-PARF(IMIX-1)
                  IF(ID.EQ.4) DWT=PARF(IMIX-1)-PARF(IMIX)
                  IF(ID.EQ.5) DWT=PARF(IMIX)
                  KFM=110*(ID-2)+KFJ
                  PMM=PMAS(PYCOMP(KFM),1)-PMAS(PYCOMP(KFM),3)
                  PARF(IPOS+5*KMUL+ID)=QWT*SWT*DWT*EXP(-PARF(193)*PMM)
                  IF(KMUL.EQ.0.AND.ID.GT.3) THEN
                     WTFAIL=WTFAIL+QWT*SWT*DWT*(1D0-PARJ(21+ID))
                     PARF(IPOS+5*KMUL+ID)=
     &                    PARF(IPOS+5*KMUL+ID)*PARJ(21+ID)
                  ENDIF
                  WTTOT=WTTOT+PARF(IPOS+5*KMUL+ID)
  170          CONTINUE
            ENDIF
  180    CONTINUE
  190 CONTINUE
                  DO 200 IMES=1,30
                     PARF(IPOS+IMES)=PARF(IPOS+IMES)/(1D0-WTFAIL)
  200             CONTINUE
                  IF(MR.EQ.180) PARF(140)=
     &                 MAX(PARF(140),WTTOT/(1D0-WTFAIL))
                  IF(MR.EQ.170) PARF(139-KFQPOP/3)=
     &                 MAX(PARF(139-KFQPOP/3),WTTOT/(1D0-WTFAIL))
                  IPOS=IPOS+30
  210           CONTINUE
  220         CONTINUE
  230       CONTINUE
  240    CONTINUE
         IF(PARF(139).GT.1D-10) PARF(138)=PARF(138)/PARF(139)
         MSTU(121)=0
 
         PARF(186)=PARF(186)/PARF(182)
         PARF(185)=PARF(185)/PARF(181)
      ENDIF
 
C..Recombine diquark weights to flavour and spin ratios
      DO 250 I=150,170,10
         WSWQ=(2D0*(PARF(I+1)+PARF(I+5))+PARF(I+3))/
     &        (1D0+PARF(I+7)+PARF(I+4)+PARF(I+2)+PARF(I+6))
         WSSWSQ=PARF(I+3)/(PARF(I+1)+PARF(I+5))
         WQSWQQ=2D0*(PARF(I+2)+PARF(I+6))/(1D0+PARF(I+7)+PARF(I+4))
         WUUWQQ=PARF(I+4)/(1D0+PARF(I+7)+PARF(I+4))
         PARF(I+5)=PARF(I+5)/PARF(I+1)
         PARF(I+6)=PARF(I+6)/PARF(I+2)
         PARF(I+1)=WSWQ
         PARF(I+2)=WQSWQQ
         PARF(I+3)=WSSWSQ
         PARF(I+4)=WUUWQQ
  250 CONTINUE
      RETURN
      END

C*********************************************************************
 
C...PYLIST
C...Gives program heading, or lists an event, or particle
C...data, or current parameter values.
 
      SUBROUTINE PYLIST(MLIST)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter statement to help give large particle numbers.
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KEXCIT=4000000)
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(4000,2),BRAT(4000),KFDP(4000,5)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/,/PYDAT3/
C...Local arrays, character variables and data.
      CHARACTER CHAP*16,CHAC*16,CHAN*16,CHAD(5)*16,CHDL(7)*4
      DIMENSION PS(6)
      DATA CHDL/'(())',' ','()','!!','<>','==','(==)'/
 
C...Initialization printout: version number and date of last change.
      IF(MLIST.EQ.0.OR.MSTU(12).EQ.1) THEN
        CALL PYLOGO
        MSTU(12)=0
        IF(MLIST.EQ.0) RETURN
      ENDIF
 
C...List event data, including additional lines after N.
      IF(MLIST.GE.1.AND.MLIST.LE.3) THEN
        IF(MLIST.EQ.1) WRITE(MSTU(11),5100)
        IF(MLIST.EQ.2) WRITE(MSTU(11),5200)
        IF(MLIST.EQ.3) WRITE(MSTU(11),5300)
        LMX=12
        IF(MLIST.GE.2) LMX=16
        ISTR=0
        IMAX=N
        IF(MSTU(2).GT.0) IMAX=MSTU(2)
        DO 120 I=MAX(1,MSTU(1)),MAX(IMAX,N+MAX(0,MSTU(3)))
          IF((I.GT.IMAX.AND.I.LE.N).OR.K(I,1).LT.0) GOTO 120
 
C...Get particle name, pad it and check it is not too long.
          CALL PYNAME(K(I,2),CHAP)
          LEN=0
          DO 100 LEM=1,16
            IF(CHAP(LEM:LEM).NE.' ') LEN=LEM
  100     CONTINUE
          MDL=(K(I,1)+19)/10
          LDL=0
          IF(MDL.EQ.2.OR.MDL.GE.8) THEN
            CHAC=CHAP
            IF(LEN.GT.LMX) CHAC(LMX:LMX)='?'
          ELSE
            LDL=1
            IF(MDL.EQ.1.OR.MDL.EQ.7) LDL=2
            IF(LEN.EQ.0) THEN
              CHAC=CHDL(MDL)(1:2*LDL)//' '
            ELSE
              CHAC=CHDL(MDL)(1:LDL)//CHAP(1:MIN(LEN,LMX-2*LDL))//
     &        CHDL(MDL)(LDL+1:2*LDL)//' '
              IF(LEN+2*LDL.GT.LMX) CHAC(LMX:LMX)='?'
            ENDIF
          ENDIF
 
C...Add information on string connection.
          IF(K(I,1).EQ.1.OR.K(I,1).EQ.2.OR.K(I,1).EQ.11.OR.K(I,1).EQ.12)
     &    THEN
            KC=PYCOMP(K(I,2))
            KCC=0
            IF(KC.NE.0) KCC=KCHG(KC,2)
            IF(IABS(K(I,2)).EQ.39) THEN
              IF(LEN+2*LDL+3.LE.LMX) CHAC(LMX-1:LMX-1)='X'
            ELSEIF(KCC.NE.0.AND.ISTR.EQ.0) THEN
              ISTR=1
              IF(LEN+2*LDL+3.LE.LMX) CHAC(LMX-1:LMX-1)='A'
            ELSEIF(KCC.NE.0.AND.(K(I,1).EQ.2.OR.K(I,1).EQ.12)) THEN
              IF(LEN+2*LDL+3.LE.LMX) CHAC(LMX-1:LMX-1)='I'
            ELSEIF(KCC.NE.0) THEN
              ISTR=0
              IF(LEN+2*LDL+3.LE.LMX) CHAC(LMX-1:LMX-1)='V'
            ENDIF
          ENDIF
 
C...Write data for particle/jet.
          IF(MLIST.EQ.1.AND.ABS(P(I,4)).LT.9999D0) THEN
            WRITE(MSTU(11),5400) I,CHAC(1:12),(K(I,J1),J1=1,3),
     &      (P(I,J2),J2=1,5)
          ELSEIF(MLIST.EQ.1.AND.ABS(P(I,4)).LT.99999D0) THEN
            WRITE(MSTU(11),5500) I,CHAC(1:12),(K(I,J1),J1=1,3),
     &      (P(I,J2),J2=1,5)
          ELSEIF(MLIST.EQ.1) THEN
            WRITE(MSTU(11),5600) I,CHAC(1:12),(K(I,J1),J1=1,3),
     &      (P(I,J2),J2=1,5)
          ELSEIF(MSTU(5).EQ.10000.AND.(K(I,1).EQ.3.OR.K(I,1).EQ.13.OR.
     &      K(I,1).EQ.14)) THEN
            WRITE(MSTU(11),5700) I,CHAC,(K(I,J1),J1=1,3),
     &      K(I,4)/100000000,MOD(K(I,4)/10000,10000),MOD(K(I,4),10000),
     &      K(I,5)/100000000,MOD(K(I,5)/10000,10000),MOD(K(I,5),10000),
     &      (P(I,J2),J2=1,5)
          ELSE
            WRITE(MSTU(11),5800) I,CHAC,(K(I,J1),J1=1,5),
     &      (P(I,J2),J2=1,5)
          ENDIF
          IF(MLIST.EQ.3) WRITE(MSTU(11),5900) (V(I,J),J=1,5)
 
C...Insert extra separator lines specified by user.
          IF(MSTU(70).GE.1) THEN
            ISEP=0
            DO 110 J=1,MIN(10,MSTU(70))
              IF(I.EQ.MSTU(70+J)) ISEP=1
  110       CONTINUE
            IF(ISEP.EQ.1.AND.MLIST.EQ.1) WRITE(MSTU(11),6000)
            IF(ISEP.EQ.1.AND.MLIST.GE.2) WRITE(MSTU(11),6100)
          ENDIF
  120   CONTINUE
 
C...Sum of charges and momenta.
        DO 130 J=1,6
          PS(J)=PYP(0,J)
  130   CONTINUE
        IF(MLIST.EQ.1.AND.ABS(PS(4)).LT.9999D0) THEN
          WRITE(MSTU(11),6200) PS(6),(PS(J),J=1,5)
        ELSEIF(MLIST.EQ.1.AND.ABS(PS(4)).LT.99999D0) THEN
          WRITE(MSTU(11),6300) PS(6),(PS(J),J=1,5)
        ELSEIF(MLIST.EQ.1) THEN
          WRITE(MSTU(11),6400) PS(6),(PS(J),J=1,5)
        ELSE
          WRITE(MSTU(11),6500) PS(6),(PS(J),J=1,5)
        ENDIF
 
C...Give simple list of KF codes defined in program.
      ELSEIF(MLIST.EQ.11) THEN
        WRITE(MSTU(11),6600)
        DO 140 KF=1,80
          CALL PYNAME(KF,CHAP)
          CALL PYNAME(-KF,CHAN)
          IF(CHAP.NE.' '.AND.CHAN.EQ.' ') WRITE(MSTU(11),6700) KF,CHAP
          IF(CHAN.NE.' ') WRITE(MSTU(11),6700) KF,CHAP,-KF,CHAN
  140   CONTINUE
        DO 170 KFLS=1,3,2
          DO 160 KFLA=1,5
            DO 150 KFLB=1,KFLA-(3-KFLS)/2
              KF=1000*KFLA+100*KFLB+KFLS
              CALL PYNAME(KF,CHAP)
              CALL PYNAME(-KF,CHAN)
              WRITE(MSTU(11),6700) KF,CHAP,-KF,CHAN
  150       CONTINUE
  160     CONTINUE
  170   CONTINUE
        KF=130
        CALL PYNAME(KF,CHAP)
        WRITE(MSTU(11),6700) KF,CHAP
        KF=310
        CALL PYNAME(KF,CHAP)
        WRITE(MSTU(11),6700) KF,CHAP
        DO 200 KMUL=0,5
          KFLS=3
          IF(KMUL.EQ.0.OR.KMUL.EQ.3) KFLS=1
          IF(KMUL.EQ.5) KFLS=5
          KFLR=0
          IF(KMUL.EQ.2.OR.KMUL.EQ.3) KFLR=1
          IF(KMUL.EQ.4) KFLR=2
          DO 190 KFLB=1,5
            DO 180 KFLC=1,KFLB-1
              KF=10000*KFLR+100*KFLB+10*KFLC+KFLS
              CALL PYNAME(KF,CHAP)
              CALL PYNAME(-KF,CHAN)
              WRITE(MSTU(11),6700) KF,CHAP,-KF,CHAN
  180       CONTINUE
            KF=10000*KFLR+110*KFLB+KFLS
            CALL PYNAME(KF,CHAP)
            WRITE(MSTU(11),6700) KF,CHAP
  190     CONTINUE
  200   CONTINUE
        KF=100443
        CALL PYNAME(KF,CHAP)
        WRITE(MSTU(11),6700) KF,CHAP
        KF=100553
        CALL PYNAME(KF,CHAP)
        WRITE(MSTU(11),6700) KF,CHAP
        DO 240 KFLSP=1,3
          KFLS=2+2*(KFLSP/3)
          DO 230 KFLA=1,5
            DO 220 KFLB=1,KFLA
              DO 210 KFLC=1,KFLB
                IF(KFLSP.EQ.1.AND.(KFLA.EQ.KFLB.OR.KFLB.EQ.KFLC))
     &          GOTO 210
                IF(KFLSP.EQ.2.AND.KFLA.EQ.KFLC) GOTO 210
                IF(KFLSP.EQ.1) KF=1000*KFLA+100*KFLC+10*KFLB+KFLS
                IF(KFLSP.GE.2) KF=1000*KFLA+100*KFLB+10*KFLC+KFLS
                CALL PYNAME(KF,CHAP)
                CALL PYNAME(-KF,CHAN)
                WRITE(MSTU(11),6700) KF,CHAP,-KF,CHAN
  210         CONTINUE
  220       CONTINUE
  230     CONTINUE
  240   CONTINUE
        DO 250 KF=KSUSY1+1,KSUSY1+40
          CALL PYNAME(KF,CHAP)
          CALL PYNAME(-KF,CHAN)
          IF(CHAP.NE.' '.AND.CHAN.EQ.' ') WRITE(MSTU(11),6700) KF,CHAP
          IF(CHAN.NE.' ') WRITE(MSTU(11),6700) KF,CHAP,-KF,CHAN
  250   CONTINUE
        DO 260 KF=KSUSY2+1,KSUSY2+40
          CALL PYNAME(KF,CHAP)
          CALL PYNAME(-KF,CHAN)
          IF(CHAP.NE.' '.AND.CHAN.EQ.' ') WRITE(MSTU(11),6700) KF,CHAP
          IF(CHAN.NE.' ') WRITE(MSTU(11),6700) KF,CHAP,-KF,CHAN
  260   CONTINUE
        DO 270 KF=KEXCIT+1,KEXCIT+40
          CALL PYNAME(KF,CHAP)
          CALL PYNAME(-KF,CHAN)
          IF(CHAP.NE.' '.AND.CHAN.EQ.' ') WRITE(MSTU(11),6700) KF,CHAP
          IF(CHAN.NE.' ') WRITE(MSTU(11),6700) KF,CHAP,-KF,CHAN
  270   CONTINUE
 
C...List parton/particle data table. Check whether to be listed.
      ELSEIF(MLIST.EQ.12) THEN
        WRITE(MSTU(11),6800)
        DO 300 KC=1,MSTU(6)
          KF=KCHG(KC,4)
          IF(KF.EQ.0) GOTO 300
          IF(KF.LT.MSTU(1).OR.(MSTU(2).GT.0.AND.KF.GT.MSTU(2)))
     &    GOTO 300
 
C...Find particle name and mass. Print information.
          CALL PYNAME(KF,CHAP)
          IF(KF.LE.100.AND.CHAP.EQ.' '.AND.MDCY(KC,2).EQ.0) GOTO 300
          CALL PYNAME(-KF,CHAN)
          WRITE(MSTU(11),6900) KF,KC,CHAP,CHAN,(KCHG(KC,J1),J1=1,3),
     &    (PMAS(KC,J2),J2=1,4),MDCY(KC,1)
 
C...Particle decay: channel number, branching ratios, matrix element,
C...decay products.
          DO 290 IDC=MDCY(KC,2),MDCY(KC,2)+MDCY(KC,3)-1
            DO 280 J=1,5
              CALL PYNAME(KFDP(IDC,J),CHAD(J))
  280       CONTINUE
            WRITE(MSTU(11),7000) IDC,MDME(IDC,1),MDME(IDC,2),BRAT(IDC),
     &      (CHAD(J),J=1,5)
  290     CONTINUE
  300   CONTINUE
 
C...List parameter value table.
      ELSEIF(MLIST.EQ.13) THEN
        WRITE(MSTU(11),7100)
        DO 310 I=1,200
          WRITE(MSTU(11),7200) I,MSTU(I),PARU(I),MSTJ(I),PARJ(I),PARF(I)
  310   CONTINUE
      ENDIF
 
C...Format statements for output on unit MSTU(11) (by default 6).
 5100 FORMAT(///28X,'Event listing (summary)'//4X,'I particle/jet KS',
     &5X,'KF  orig    p_x      p_y      p_z       E        m'/)
 5200 FORMAT(///28X,'Event listing (standard)'//4X,'I  particle/jet',
     &'  K(I,1)   K(I,2) K(I,3)     K(I,4)      K(I,5)       P(I,1)',
     &'       P(I,2)       P(I,3)       P(I,4)       P(I,5)'/)
 5300 FORMAT(///28X,'Event listing (with vertices)'//4X,'I  particle/j',
     &'et  K(I,1)   K(I,2) K(I,3)     K(I,4)      K(I,5)       P(I,1)',
     &'       P(I,2)       P(I,3)       P(I,4)       P(I,5)'/73X,
     &'V(I,1)       V(I,2)       V(I,3)       V(I,4)       V(I,5)'/)
 5400 FORMAT(1X,I4,1X,A12,1X,I2,I8,1X,I4,5F9.3)
 5500 FORMAT(1X,I4,1X,A12,1X,I2,I8,1X,I4,5F9.2)
 5600 FORMAT(1X,I4,1X,A12,1X,I2,I8,1X,I4,5F9.1)
 5700 FORMAT(1X,I4,2X,A16,1X,I3,1X,I9,1X,I4,2(3X,I1,2I4),5F13.5)
 5800 FORMAT(1X,I4,2X,A16,1X,I3,1X,I9,1X,I4,2(3X,I9),5F13.5)
 5900 FORMAT(66X,5(1X,F12.3))
 6000 FORMAT(1X,78('='))
 6100 FORMAT(1X,130('='))
 6200 FORMAT(19X,'sum:',F6.2,5X,5F9.3)
 6300 FORMAT(19X,'sum:',F6.2,5X,5F9.2)
 6400 FORMAT(19X,'sum:',F6.2,5X,5F9.1)
 6500 FORMAT(19X,'sum charge:',F6.2,3X,'sum momentum and inv. mass:',
     &5F13.5)
 6600 FORMAT(///20X,'List of KF codes in program'/)
 6700 FORMAT(4X,I9,4X,A16,6X,I9,4X,A16)
 6800 FORMAT(///30X,'Particle/parton data table'//8X,'KF',5X,'KC',4X,
     &'particle',8X,'antiparticle',6X,'chg  col  anti',8X,'mass',7X,
     &'width',7X,'w-cut',5X,'lifetime',1X,'decay'/11X,'IDC',1X,'on/off',
     &1X,'ME',3X,'Br.rat.',4X,'decay products')
 6900 FORMAT(/1X,I9,3X,I4,4X,A16,A16,3I5,1X,F12.5,2(1X,F11.5),
     &1X,1P,E13.5,3X,I2)
 7000 FORMAT(10X,I4,2X,I3,2X,I3,2X,F10.6,4X,5A16)
 7100 FORMAT(///20X,'Parameter value table'//4X,'I',3X,'MSTU(I)',
     &8X,'PARU(I)',3X,'MSTJ(I)',8X,'PARJ(I)',8X,'PARF(I)')
 7200 FORMAT(1X,I4,1X,I9,1X,F14.5,1X,I9,1X,F14.5,1X,F14.5)
 
      RETURN
      END

C*********************************************************************
 
C...PYANGL
C...Reconstructs an angle from given x and y coordinates.
 
      FUNCTION PYANGL(X,Y)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /PYDAT1/
 
      PYANGL=0D0
      R=SQRT(X**2+Y**2)
      IF(R.LT.1D-20) RETURN
      IF(ABS(X)/R.LT.0.8D0) THEN
        PYANGL=SIGN(ACOS(X/R),Y)
      ELSE
        PYANGL=ASIN(Y/R)
        IF(X.LT.0D0.AND.PYANGL.GE.0D0) THEN
          PYANGL=PARU(1)-PYANGL
        ELSEIF(X.LT.0D0) THEN
          PYANGL=-PARU(1)-PYANGL
        ENDIF
      ENDIF
 
      RETURN
      END

C*********************************************************************
 
C...PYLOGO
C...Writes a logo for the program.
 
      SUBROUTINE PYLOGO
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter for length of information block.
      PARAMETER (IREFER=17)
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      SAVE /PYDAT1/,/PYPARS/
C...Local arrays and character variables.
      INTEGER IDATI(6)
      CHARACTER MONTH(12)*3, LOGO(48)*32, REFER(2*IREFER)*36, LINE*79,
     &VERS*1, SUBV*3, DATE*2, YEAR*4, HOUR*2, MINU*2, SECO*2
 
C...Data on months, logo, titles, and references.
      DATA MONTH/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',
     &'Oct','Nov','Dec'/
      DATA (LOGO(J),J=1,19)/
     &'            *......*            ',
     &'       *:::!!:::::::::::*       ',
     &'    *::::::!!::::::::::::::*    ',
     &'  *::::::::!!::::::::::::::::*  ',
     &' *:::::::::!!:::::::::::::::::* ',
     &' *:::::::::!!:::::::::::::::::* ',
     &'  *::::::::!!::::::::::::::::*! ',
     &'    *::::::!!::::::::::::::* !! ',
     &'    !! *:::!!:::::::::::*    !! ',
     &'    !!     !* -><- *         !! ',
     &'    !!     !!                !! ',
     &'    !!     !!                !! ',
     &'    !!                       !! ',
     &'    !!        ep             !! ',
     &'    !!                       !! ',
     &'    !!                 pp    !! ',
     &'    !!   e+e-                !! ',
     &'    !!                       !! ',
     &'    !!                          '/
      DATA (LOGO(J),J=20,38)/
     &'Welcome to the Lund Monte Carlo!',
     &'                                ',
     &'PPP  Y   Y TTTTT H   H III   A  ',
     &'P  P  Y Y    T   H   H  I   A A ',
     &'PPP    Y     T   HHHHH  I  AAAAA',
     &'P      Y     T   H   H  I  A   A',
     &'P      Y     T   H   H III A   A',
     &'                                ',
     &'This is PYTHIA version x.xxx    ',
     &'Last date of change: xx xxx 199x',
     &'                                ',
     &'Now is xx xxx 199x at xx:xx:xx  ',
     &'                                ',
     &'Disclaimer: this program comes  ',
     &'without any guarantees. Beware  ',
     &'of errors and use common sense  ',
     &'when interpreting results.      ',
     &'                                ',
     &'Copyright T. Sjostrand (1997)   '/
      DATA (REFER(J),J=1,18)/
     &'An archive of program versions and d',
     &'ocumentation is found on the web:   ',
     &'http://www.thep.lu.se/tf2/staff/torb',
     &'jorn/Pythia.html                    ',
     &'                                    ',
     &'                                    ',
     &'When you cite this program, currentl',
     &'y the official reference is         ',
     &'T. Sjostrand, Computer Physics Commu',
     &'n. 82 (1994) 74.                    ',
     &'The supersymmetry extensions are des',
     &'cribed in                           ',
     &'S. Mrenna, Computer Physics Commun. ',
     &'101 (1997) 232                      ',
     &'Also remember that the program, to a',
     &' large extent, represents original  ',
     &'physics research. Other publications',
     &' of special relevance to your       '/
      DATA (REFER(J),J=19,2*IREFER)/
     &'studies may therefore deserve separa',
     &'te mention.                         ',
     &'                                    ',
     &'                                    ',
     &'Main author: Torbjorn Sjostrand; Dep',
     &'artment of Theoretical Physics 2,   ',
     &'  Lund University, Solvegatan 14A, S',
     &'-223 62 Lund, Sweden;               ',
     &'  phone: + 46 - 46 - 222 48 16; e-ma',
     &'il: torbjorn@thep.lu.se             ',
     &'SUSY author: Stephen Mrenna, Argonne',
     &' National Laboratory,               ',
     &'  9700 South Cass Avenue, Argonne, I',
     &'L 60439, USA;                       ',
     &'  phone: + 1 - 630 - 252 - 7615; e-m',
     &'ail: mrenna@hep.anl.gov             '/
 
C...Check that PYDATA linked.
      IF(MSTP(183)/10.NE.199) THEN
        WRITE(MSTU(11),'(1X,A)')
     &  'Error: PYDATA has not been linked.'
        WRITE(MSTU(11),'(1X,A)') 'Execution stopped!'
        STOP
 
C...Write current version number and current date+time.
      ELSE
        WRITE(VERS,'(I1)') MSTP(181)
        LOGO(28)(24:24)=VERS
        WRITE(SUBV,'(I3)') MSTP(182)
        LOGO(28)(26:28)=SUBV
        IF(MSTP(182).LT.100) LOGO(28)(26:26)='0'
        WRITE(DATE,'(I2)') MSTP(185)
        LOGO(29)(22:23)=DATE
        LOGO(29)(25:27)=MONTH(MSTP(184))
        WRITE(YEAR,'(I4)') MSTP(183)
        LOGO(29)(29:32)=YEAR
        CALL PYTIME(IDATI)
        IF(IDATI(1).LE.0) THEN
          LOGO(31)='                                '
        ELSE
          WRITE(DATE,'(I2)') IDATI(3)
          LOGO(31)(8:9)=DATE
          LOGO(31)(11:13)=MONTH(MAX(1,MIN(12,IDATI(2))))
          WRITE(YEAR,'(I4)') IDATI(1)
          LOGO(31)(15:18)=YEAR
          WRITE(HOUR,'(I2)') IDATI(4)
          LOGO(31)(23:24)=HOUR
          WRITE(MINU,'(I2)') IDATI(5)
          LOGO(31)(26:27)=MINU
          IF(IDATI(5).LT.10) LOGO(31)(26:26)='0'
          WRITE(SECO,'(I2)') IDATI(6)
          LOGO(31)(29:30)=SECO
          IF(IDATI(6).LT.10) LOGO(31)(29:29)='0'
        ENDIF
      ENDIF
 
C...Loop over lines in header. Define page feed and side borders.
      DO 100 ILIN=1,29+IREFER
        LINE=' '
        IF(ILIN.EQ.1) THEN
          LINE(1:1)='1'
        ELSE
          LINE(2:3)='**'
          LINE(78:79)='**'
        ENDIF
 
C...Separator lines and logos.
        IF(ILIN.EQ.2.OR.ILIN.EQ.3.OR.ILIN.GE.28+IREFER) THEN
          LINE(4:77)='***********************************************'//
     &    '***************************'
        ELSEIF(ILIN.GE.6.AND.ILIN.LE.24) THEN
          LINE(6:37)=LOGO(ILIN-5)
          LINE(44:75)=LOGO(ILIN+14)
        ELSEIF(ILIN.GE.26.AND.ILIN.LE.25+IREFER) THEN
          LINE(5:40)=REFER(2*ILIN-51)
          LINE(41:76)=REFER(2*ILIN-50)
        ENDIF
 
C...Write lines to appropriate unit.
        WRITE(MSTU(11),'(A79)') LINE
  100 CONTINUE
 
      RETURN
      END

C*********************************************************************
 
C...PYNAME
C...Gives the particle/parton name as a character string.
 
      SUBROUTINE PYNAME(KF,CHAU)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT4/CHAF(500,2)
      CHARACTER CHAF*16
      SAVE /PYDAT1/,/PYDAT2/,/PYDAT4/
C...Local character variable.
      CHARACTER CHAU*16
 
C...Read out code with distinction particle/antiparticle.
      CHAU=' '
      KC=PYCOMP(KF)
      IF(KC.NE.0) CHAU=CHAF(KC,(3-ISIGN(1,KF))/2)
 
 
      RETURN
      END

C*********************************************************************
 
C...PYP
C...Provides various real-valued event related data.
 
      FUNCTION PYP(I,J)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/
C...Local array.
      DIMENSION PSUM(4)
 
C...Set default value. For I = 0 sum of momenta or charges,
C...or invariant mass of system.
      PYP=0D0
      IF(I.LT.0.OR.I.GT.MSTU(4).OR.J.LE.0) THEN
      ELSEIF(I.EQ.0.AND.J.LE.4) THEN
        DO 100 I1=1,N
          IF(K(I1,1).GT.0.AND.K(I1,1).LE.10) PYP=PYP+P(I1,J)
  100   CONTINUE
      ELSEIF(I.EQ.0.AND.J.EQ.5) THEN
        DO 120 J1=1,4
          PSUM(J1)=0D0
          DO 110 I1=1,N
            IF(K(I1,1).GT.0.AND.K(I1,1).LE.10) PSUM(J1)=PSUM(J1)+
     &      P(I1,J1)
  110     CONTINUE
  120   CONTINUE
        PYP=SQRT(MAX(0D0,PSUM(4)**2-PSUM(1)**2-PSUM(2)**2-PSUM(3)**2))
      ELSEIF(I.EQ.0.AND.J.EQ.6) THEN
        DO 130 I1=1,N
          IF(K(I1,1).GT.0.AND.K(I1,1).LE.10) PYP=PYP+PYCHGE(K(I1,2))/3D0
  130   CONTINUE
      ELSEIF(I.EQ.0) THEN
 
C...Direct readout of P matrix.
      ELSEIF(J.LE.5) THEN
        PYP=P(I,J)
 
C...Charge, total momentum, transverse momentum, transverse mass.
      ELSEIF(J.LE.12) THEN
        IF(J.EQ.6) PYP=PYCHGE(K(I,2))/3D0
        IF(J.EQ.7.OR.J.EQ.8) PYP=P(I,1)**2+P(I,2)**2+P(I,3)**2
        IF(J.EQ.9.OR.J.EQ.10) PYP=P(I,1)**2+P(I,2)**2
        IF(J.EQ.11.OR.J.EQ.12) PYP=P(I,5)**2+P(I,1)**2+P(I,2)**2
        IF(J.EQ.8.OR.J.EQ.10.OR.J.EQ.12) PYP=SQRT(PYP)
 
C...Theta and phi angle in radians or degrees.
      ELSEIF(J.LE.16) THEN
        IF(J.LE.14) PYP=PYANGL(P(I,3),SQRT(P(I,1)**2+P(I,2)**2))
        IF(J.GE.15) PYP=PYANGL(P(I,1),P(I,2))
        IF(J.EQ.14.OR.J.EQ.16) PYP=PYP*180D0/PARU(1)
 
C...True rapidity, rapidity with pion mass, pseudorapidity.
      ELSEIF(J.LE.19) THEN
        PMR=0D0
        IF(J.EQ.17) PMR=P(I,5)
        IF(J.EQ.18) PMR=PYMASS(211)
        PR=MAX(1D-20,PMR**2+P(I,1)**2+P(I,2)**2)
        PYP=SIGN(LOG(MIN((SQRT(PR+P(I,3)**2)+ABS(P(I,3)))/SQRT(PR),
     &  1D20)),P(I,3))
 
C...Energy and momentum fractions (only to be used in CM frame).
      ELSEIF(J.LE.25) THEN
        IF(J.EQ.20) PYP=2D0*SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2)/PARU(21)
        IF(J.EQ.21) PYP=2D0*P(I,3)/PARU(21)
        IF(J.EQ.22) PYP=2D0*SQRT(P(I,1)**2+P(I,2)**2)/PARU(21)
        IF(J.EQ.23) PYP=2D0*P(I,4)/PARU(21)
        IF(J.EQ.24) PYP=(P(I,4)+P(I,3))/PARU(21)
        IF(J.EQ.25) PYP=(P(I,4)-P(I,3))/PARU(21)
      ENDIF
 
      RETURN
      END

C*********************************************************************
 
C...PYTAUD
C...Dummy routine, to be replaced by user, to handle the decay of a
C...polarized tau lepton.
C...Input:
C...ITAU is the position where the decaying tau is stored in /PYJETS/.
C...IORIG is the position where the mother of the tau is stored;
C...     is 0 when the mother is not stored.
C...KFORIG is the flavour of the mother of the tau;
C...     is 0 when the mother is not known.
C...Note that IORIG=0 does not necessarily imply KFORIG=0;
C...     e.g. in B hadron semileptonic decays the W  propagator
C...     is not explicitly stored but the W code is still unambiguous.
C...Output:
C...NDECAY is the number of decay products in the current tau decay.
C...These decay products should be added to the /PYJETS/ common block,
C...in positions N+1 through N+NDECAY. For each product I you must
C...give the flavour codes K(I,2) and the five-momenta P(I,1), P(I,2),
C...P(I,3), P(I,4) and P(I,5). The rest will be stored automatically.
 
      SUBROUTINE PYTAUD(ITAU,IORIG,KFORIG,NDECAY)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /PYJETS/,/PYDAT1/
 
C...Stop program if this routine is ever called.
C...You should not copy these lines to your own routine.
      NDECAY=ITAU+IORIG+KFORIG
      WRITE(MSTU(11),5000)
      IF(PYR(0).LT.10D0) STOP
 
C...Format for error printout.
 5000 FORMAT(1X,'Error: you did not link your PYTAUD routine ',
     &'correctly.'/1X,'Dummy routine in PYTHIA file called instead.'/
     &1X,'Execution stopped!')
 
      RETURN
      END

C*********************************************************************
 
C...PYTIME
C...Finds current date and time.
C...Since this task is not standardized in Fortran 77, the routine
C...is dummy, to be replaced by the user. Examples are given for
C...the Fortran 90 routine and DEC Fortran 77, and what to do if
C...you do not have access to suitable routines.
 
      SUBROUTINE PYTIME(IDATI)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP
      CHARACTER*8 ATIME
C...Local array.
      INTEGER IDATI(6),IDTEMP(3)
 
C...Example 0: if you do not have suitable routines.
      DO 100 J=1,6
      IDATI(J)=0
  100 CONTINUE
 
C...Example 1: Fortran 90 routine.
C      INTEGER IVAL(8)
C      CALL DATE_AND_TIME(VALUES=IVAL)
C      IDATI(1)=IVAL(1)
C      IDATI(2)=IVAL(2)
C      IDATI(3)=IVAL(3)
C      IDATI(4)=IVAL(5)
C      IDATI(5)=IVAL(6)
C      IDATI(6)=IVAL(7)
 
C...Example 2: DEC Fortran 77.
C      CALL IDATE(IMON,IDAY,IYEAR)
C      IDATI(1)=1900+IYEAR
C      IDATI(2)=IMON
C      IDATI(3)=IDAY
C      CALL ITIME(IHOUR,IMIN,ISEC)
C      IDATI(4)=IHOUR
C      IDATI(5)=IMIN
C      IDATI(6)=ISEC
 
C...Example 3: DEC Fortran
C      CALL IDATE(IMON,IDAY,IYEAR)
C      IDATI(1)=1900+IYEAR
C      IDATI(2)=IMON
C      IDATI(3)=IDAY
C      CALL TIME(ATIME)
C      IHOUR=0
C      IMIN=0
C      ISEC=0
C      READ(ATIME(1:2),'(I2)') IHOUR
C      READ(ATIME(4:5),'(I2)') IMIN
C      READ(ATIME(7:8),'(I2)') ISEC
C      IDATI(4)=IHOUR
C      IDATI(5)=IMIN
C      IDATI(6)=ISEC
 
C...Example 4: GNU LINUX libU77.
C      CALL IDATE(IDTEMP)
C      IDATI(1)=IDTEMP(3)
C      IDATI(2)=IDTEMP(2)
C      IDATI(3)=IDTEMP(1)
C      CALL ITIME(IDTEMP)
C      IDATI(4)=IDTEMP(1)
C      IDATI(5)=IDTEMP(2)
C      IDATI(6)=IDTEMP(3)
 
      RETURN
      END

C*********************************************************************
 
C...PYCHGE
C...Gives three times the charge for a particle/parton.
 
      FUNCTION PYCHGE(KF)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYDAT2/
 
C...Read out charge and change sign for antiparticle.
      PYCHGE=0
      KC=PYCOMP(KF)
      IF(KC.NE.0) PYCHGE=KCHG(KC,1)*ISIGN(1,KF)
 
      RETURN
      END


C*********************************************************************
 
C...PYR
C...Generates random numbers uniformly distributed between
C...0 and 1, excluding the endpoints.
C     linked to SIBYLL random number generator \FR'14
 
C      DOUBLE PRECISION FUNCTION PYR(IDUMMY)
C      IMPLICIT NONE
C      double precision s_rndm
C      integer idummy
C      pyr = s_rndm(idummy)
C      end
