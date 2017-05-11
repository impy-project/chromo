C=======================================================================
C          SSSSSS   IIIIIII  BBBBB   YY      YY   L        L
C         S            I     B    B    YY  YY     L        L
C          SSSSS       I     BBBBB       YY       L        L
C               S      I     B    B      YY       L        L
C         SSSSSS    IIIIIII  BBBBB       YY       LLLLLLL  LLLLLLL
C=======================================================================
C  Code for SIBYLL:  hadronic interaction Monte Carlo event generator
C=======================================================================
C   Version 2.3rc1     (August-8-2014)
C	with CHARM production
C
C       By   Eun-Joo Ahn
C            Ralph Engel
C            R.S. Fletcher
C            T.K. Gaisser
C            Paolo Lipari
C            Felix Riehn
C            Todor Stanev
C
C-----------------------------------------------------------------------
C***  Please  have people who want this code contact one of the authors.
C***  Please report any problems.       ****
C
C      For a correct copy contact:
C                sein@fnal.gov
C                ralph.engel@kit.edu
C                gaisser@bartol.udel.edu
C                stanev@bartol.udel.edu
C                paolo.lipari@roma1.infn.it
C                felix.riehn@kit.edu
C
C-----------------------------------------------------------------------

      SUBROUTINE SIBYLL (K_beam, IATARG, Ecm)
C-----------------------------------------------------------------------
C...Main routine for the production of hadronic events,
C.  generates an inelastic hadronic interaction of 
C.  a `projectile particle' of code K_beam with a 
C.  target nucleus of mass number A = IATARG (integer)
C.  IATARG = 0 is an "air" nucleus  (superposition of oxygen and nitrogen)
C.  with c.m. energy for the hadron-nucleon system Ecm (GeV)
C.  
C.  Allowed values of K_beam: 7,8,9,10,11,12,13,14,-13,-14
C.                            pi+-,K+-,KL,KS,p,n,pbar,nbar
C. 
C.  The output is contained in COMMON /S_PLIST/ that contains:
C.
C.     NP           number of final particles
C.     P(1:NP, 1:5) 4-momenta + masses of the final particles 
C.     LLIST (1:NP) codes of final particles.
C.  the reaction is studied in the c.m. of  hadron-nucleon system
C.
C.  The COMMON block /S_CHIST/ contains information about the
C.  the structure of the  generated event:
C.    NW   = number of wounded nucleons
C.    NJET = total number of hard interactions
C.    NSOF = total number of soft interactions
C.    NNSOF (1:NW) = number of soft pomeron cuts in each interaction
C.    NNJET (1:NW) = number of minijets produced in each interaction 
C.    XJ1 (1:Index) = x1  for each string
C.    XJ2 (1:Index) = x2   "   "     "
C.    PTJET (1:Index) = pT   "   "     "
C.    NNPJET (1:Index) = total number of particles in each string
C.    NNPSTR (1:2*NW) = number of particles in each `beam string'
C.    JDIF(1:NW) = diffraction code
C----------------------------------------------------------------------
      SAVE

      COMMON /S_PLIST/ P(8000,5), LLIST(8000), NP
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, kb ,kt
      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      PARAMETER (NW_max = 20)
      PARAMETER (NS_max = 20, NH_max = 80)
      PARAMETER (NJ_max = (NS_max+NH_max)*NW_max)
      COMMON /S_CHIST/ X1J(NJ_max),X2J(NJ_max),
     &    X1JSUM(NW_max),X2JSUM(NW_max),PTJET(NJ_max),PHIJET(NJ_max),
     &    NNPJET(NJ_max),NNPSTR(2*NW_max),NNSOF(NW_max),NNJET(NW_max),
     &    JDIF(NW_max),NW,NJET,NSOF
      COMMON /S_CCSTR/ X1(2*NW_max),X2(2*NW_max),
     &    PXB(2*NW_max),PYB(2*NW_max),PXT(2*NW_max),PYT(2*NW_max),
     &    IFLB(2*NW_max),IFLT(2*NW_max)
      COMMON /S_CLDIF/ LDIFF
      COMMON /S_CQDIS/ PPT0 (35),ptflag
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)

      DIMENSION LL(39)
      DATA LL /5*0,7*2,2*1,19*0,6*1/
      DATA FOX /0.257/

      if(Ndebug.gt.1) 
     &  WRITE(LUN,*)' SIBYLL: called with (K_beam,IATARG,Ecm):',
     &  K_beam,IATARG,Ecm

      kb = K_beam
      SQS = Ecm
      S = SQS*SQS

      Ncall = Ncall+1

 100  CONTINUE


      NP = 0
      NJET = 0
      NSOF = 0
      IATARGET = IATARG
      
C...Generate an 'air' interaction by choosing Nitrogen or Oxygen

      IF (IATARGET .EQ. 0) THEN
          R = S_RNDM(0)
          IATARGET = 14
          IF (R .LT. FOX)  IATARGET = 16
      ENDIF
      L = LL(IABS(KB))

C...Generate number NW wounded nucleons, and diffraction code.

1000  CALL SIB_START_EV (Ecm, L, IATARGET, NW, JDIF)

C...limits on simulation of pure diffraction dissociation
      IF((LDIFF.NE.0).and.(NW.EQ.1)) THEN 
         IF((LDIFF.EQ.-1) .AND. (JDIF(1).NE.0) ) GOTO 1000
         IF((LDIFF.EQ. 1) .AND. ((JDIF(1).NE.0).AND.(JDIF(1).NE.3)))
     +     GOTO 1000
         IF((LDIFF.EQ. 5) .AND. (JDIF(1).EQ.2)) GOTO 1000
         IF((LDIFF.GE. 2) .AND. (LDIFF.LE.4)) THEN
           JDIF(1) = LDIFF-1
         ENDIF
      ENDIF

C...`soft increase of pT'
      IF(IPAR(17).gt.0)THEN
         CALL PTSETUP_4FLV(ECM)
      ELSE
         CALL PTSETUP(ECM)
      ENDIF

C...Diffractive/non-diffractive interactions

      IF((NW.EQ.1).and.(JDIF(1).NE.0)) THEN
        CALL SIB_DIFF (KB, JDIF(1), Ecm, 1, IREJ)
      ELSE
        CALL SIB_NDIFF (KB, IATARGET, Ecm, 1, IREJ)
      ENDIF

      IF (IREJ.NE.0) THEN
        if(Ndebug.gt.2) WRITE(LUN,*)
     &    'SIBYLL: rejection (Ecm,Ncall,Nw,JDIF):',Ecm,Ncall,NW,JDIF(1)
        GOTO 100
      ENDIF

      do J=1,NP
         if (P(J,4).lt.0.0 ) then
            IF(NDEBUG.ne.0)
     +           WRITE(LUN,*)'negative energy particle!' , P(J,4)
            goto 100
         endif
      enddo

C...Check energy-momentum conservation

      CALL PFsum(1,NP,Esum,PXsum,PYsum,PZsum,NF)
      IF (ABS(Esum/(0.5*Ecm*FLOAT(NW+1)) - 1.) .GT. 1.E-03)  THEN
         IF(NDEBUG.ne.0)THEN
            WRITE(6,*)' SIBYLL: energy not conserved (L,call): '
     &           ,L,Ncall
            WRITE(LUN,*)' SIBYLL: energy not conserved (L,call): '
     &           ,L,Ncall
            WRITE(LUN,*) ' sqs_inp = ', Ecm, ' sqs_out = ', Esum
            CALL SIB_LIST(LUN)
            WRITE(LUN,*) ' SIBYLL: event rejected'
         ENDIF
         goto 100
      ENDIF

C...list final state particles
      if(Ndebug.gt.10) call sib_list(6)

      END

  
      SUBROUTINE SIB_NDIFF (K_beam, IATARGET, Ecm, Irec, IREJ)
C----------------------------------------------------------------------
C...Non-diffractive or multiple non-diff./diff. interactions
C.    Irec  flag to avoid recursive calls of SIB_DIFF and SIB_NDIFF
C----------------------------------------------------------------------
      SAVE

      COMMON /S_PLIST/ P(8000,5), LLIST(8000), NP
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, kb ,kt
      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)

      PARAMETER (NW_max = 20)
      PARAMETER (NS_max = 20, NH_max = 80)
      PARAMETER (NJ_max = (NS_max+NH_max)*NW_max)
      COMMON /S_CHIST/ X1J(NJ_max),X2J(NJ_max),
     &    X1JSUM(NW_max),X2JSUM(NW_max),PTJET(NJ_max),PHIJET(NJ_max),
     &    NNPJET(NJ_max),NNPSTR(2*NW_max),NNSOF(NW_max),NNJET(NW_max),
     &    JDIF(NW_max),NW,NJET,NSOF
      COMMON /S_CCSTR/ X1(2*NW_max),X2(2*NW_max),
     &    PXB(2*NW_max),PYB(2*NW_max),PXT(2*NW_max),PYT(2*NW_max),
     &    IFLB(2*NW_max),IFLT(2*NW_max)

      COMMON /S_CLDIF/ LDIFF
      COMMON /S_CUTOFF/ STR_mass_val, STR_mass_val_hyp, STR_mass_sea

      DIMENSION X2JET(NW_max),BET(2*NW_max),GAM(2*NW_max),EE(2*NW_max)

      DIMENSION QMAS(43),LL(39)
      DATA QMAS
     &   /2*0.35,0.6,1.5,6*0.,2*1.1,1.25,1.8,6*0.,1.25,1.1,1.25,1.8,6*0,
     &    2*1.25,1.5,2.0,6*0,1.8,1.8,2.0/
      DATA LL /5*0,7*2,2*1,19*0,6*1/
      DATA PI /3.14159/
      DOUBLE PRECISION DTHETA,DPHI,DPSTRTOT,DPT2STR,cod,sid,cof,sif
      DOUBLE PRECISION DBET(3),DGABE(4),DPSTR(5),PC(4),dP1(4),dptot1,
     &     PCt(4),dP2(4),dp1tot_t,dp0(4),dpb(5),dpt(5),dxmt
      COMMON /S_PARTO/ nforig(8000), nporig(8000), Ipflag, Nint
      CHARACTER*25 CWRN

      if(Ndebug.gt.1) 
     &  print *,' SIB_NDIFF: called with (K_beam,IATARGET,Ecm,Irec):',
     &  K_beam,IATARGET,Ecm,Irec

      if(Irec.eq.1) IPFLAG = 1

      IREJ = 1

      NP_0    = NP
      SQS_0   = SQS

      SQS   = Ecm
      S     = SQS*SQS

*        print *,' current NP,QSQ (-3) ',NP,SQS,NP_0

C...energy-dependent transverse momentum cutoff
c...EJA correction 2007.03.27
      PTmin = PAR(10)+PAR(11)*EXP(PAR(12)*SQRT(LOG(S)))
      XMIN = 4.*PTmin**2/S
      ZMIN = LOG(XMIN)

c     store default diq / q ratio
      PAR1def= PAR(1)

2000  CONTINUE
      PAR(1)= PAR1def 
C...sample multiple interaction configuration
*        print *,' current NP,QSQ (-2a) ',NP,SQS,NP_0

      L = LL(IABS(K_beam))
      DO I=1,NW
        if(JDIF(I).eq.0) then
          CALL CUT_PRO(L, SQS, PTmin, NNSOF(I), NNJET(I))
        else
          NNSOF(I) = 1
          NNJET(I) = 0
        endif
      ENDDO

*        print *,' current NP,QSQ (-2b) ',NP,SQS,NP_0

C...sample x values

      ITRY = 0
3000  CONTINUE
      ITRY = ITRY+1
      IF(ITRY.GT.5) GOTO 2000
      NP = NP_0
      NJET = 0
      NSOF = 0
      Nall = 0
      X1JET = 0.
      DO JW=1,NW
C...hard sea-sea interactions
         X2JET(JW) = 0.
         X1JSUM(JW) = 0.
         X2JSUM(JW) = 0.
         DO JJ=1,NNJET(JW)
           Nall = Nall+1
           NJET = NJET+1
*          print *,' Ncall,JW,NW,Njet,NNJET(JW),Nall',
*    &       Ncall,JW,NW,Njet,NNJET(JW),Nall
           CALL SAMPLE_hard (L,X1J(Nall),X2J(Nall),PTJET(Nall))
*        print *,' current NP,QSQ (-2c) ',NP,SQS,NP_0
           X1JET = X1JET + X1J(Nall)
           X2JET(JW) = X2JET(JW)+X2J(Nall)
           if(Ndebug.gt.2)
     &       print *,' SIB_NDIFF: hard JJ,JW,X1JET,X2JET(JW):',
     &       JJ,JW,X1JET,X2JET(JW)
           IF ((X2JET(JW).GT.0.9).OR.(X1JET.GT.0.9)) then
             if(Ndebug.gt.2) print *,
     &         ' SIB_NDIFF: not enough phase space (Ncall,Njet):',
     &         Ncall,Njet
             GOTO 3000
           ENDIF
           X1JSUM(JW) = X1JSUM(JW)+X1J(Nall)
           X2JSUM(JW) = X2JSUM(JW)+X2J(Nall)
         ENDDO

C...soft sea-sea interactions 
         NSOF_JW = 0
         DO JJ=1,NNSOF(JW)-1
c	print*,'inside jj=nnsof(jw)'
            IF(IPAR(28).eq.1)then
c               CALL SAMPLE_soft2 (STR_mass_sea,X1S,X2S,PTSOF)
               
            elseif(IPAR(28).eq.2)then
c     (1-x)**b / x distribution
               CALL SAMPLE_soft3 (STR_mass_sea,X1S,X2S,PTSOF)
            elseif(IPAR(28).eq.3)then
c               CALL SAMPLE_soft5 (STR_mass_sea,X1S,X2S,PTSOF)
               CWRN = 'SAMPLE_soft5 not incld. !'
               call sib_abrt(cwrn)
            else
               CALL SAMPLE_soft (STR_mass_sea,X1S,X2S,PTSOF)
            endif
*        print *,' current NP,QSQ (-2d) ',NP,SQS,NP_0
           IF ((X2JET(JW)+X2S.LT.0.9).AND.(X1JET+X1S.LT.0.9)) THEN
             NSOF = NSOF+1
             Nall = Nall+1
*            print *,' Ncall,JW,NW,Nsof,NNSOF(JW),Nall',
*    &         Ncall,JW,NW,Nsof,NNSOF(JW),Nall
             NSOF_JW = NSOF_JW+1
             X1J(Nall) = X1S
             X2J(Nall) = X2S
             PTjet(Nall) = PTsof
             X1JSUM(JW) = X1JSUM(JW)+X1S
             X2JSUM(JW) = X2JSUM(JW)+X2S
             X1JET = X1JET + X1S
             X2JET(JW) = X2JET(JW)+X2S
           ENDIF
           if(Ndebug.gt.2)
     &       print *,' SIB_NDIFF: soft JJ,JW,X1JET,X2JET(JW):',
     &       JJ,JW,X1JET,X2JET(JW)
         ENDDO

         NNSOF(JW) = NSOF_JW+1
 3500    CONTINUE
      ENDDO

*        print *,' current NP,QSQ (-1) ',NP,SQS,NP_0

C...Prepare 2*NW valence/sea color strings.

      CALL BEAM_SPLIT (K_beam, NW, X1, IFLB, X1JET, LXBAD )
      IF (LXBAD .EQ. 1) then
        if(Ndebug.gt.2) 
     &        WRITE(LUN,*)' BEAM_SPLIT: rejection (Ncall):',Ncall
        NP    = NP_0
        SQS   = SQS_0
        S     = SQS*SQS
        return
      ENDIF
*        print *,' current NP,QSQ (-1a) ',NP,SQS,NP_0
      DO J=1,NW

         J1=2*(J-1)+1
         J2=J1+1
*        print *,' J,J1,J2,NW ',J,J1,J2,NW
         KT=13
         IF (IATARGET .GT. 1)  KT = 13+INT(1.9999*S_RNDM(0))
         if (Irec.ne.1.and.IATARGET.eq.1 ) KT=6    !fix baryon number problem in diffraction
         CALL HSPLI (KT,IFLT(J2),IFLT(J1))
         IF(IFLT(J2).eq.0.and.IFLT(J1).eq.0)THEN
            WRITE(LUN,*) 'SIB_NDIFF: HSPLI rejection!'
            WRITE(LUN,*) ' RUN (SQS,KB,KT):',SQS,KB,KT
            STOP
         ENDIF

*        XMINA = 2.*STR_mass_val/(SQS*(1.-X2JET(J)))
         XMINA = 1./(SQS*(1.-X2JET(J)))**2
C        XMINA = 2.*0.20/(SQS*(1.-X2JET(J)))  ! change RSF. 5-92
         CHI=CHIDIS (KT,IFLT(J2),IFLT(J1))
         XVAL=1.-X2JET(J)
         IF (XVAL.LT.XMINA) GOTO 3000
         X2(J2) = MAX(CHI*XVAL,XMINA)
         X2(J2) = MIN(X2(J2),XVAL-XMINA)
         X2(J1) = XVAL-X2(J2)
      ENDDO

C...Generates primordial pT for the partons 
*        print *,' current NP,QSQ (-1b) ',NP,SQS,NP_0
      IF(IPAR(17).eq.0) THEN
         DO J=1,NW
            J1 = 2*(J-1)+1
            J2 = J1+1
            CALL PTDIS (10,PXT(J1),PYT(J1),0)
            if (j.eq.1) then
               CALL PTDIS (10,PXB(J2),PYB(J2),0)
            else
               CALL PTDIS (IFLB(J2),PXB(J2),PYB(J2),0)
            endif
            PXB(J1) = -PXB(J2)
            PYB(J1) = -PYB(J2)
            PXT(J2) = -PXT(J1)
            PYT(J2) = -PYT(J1)
         ENDDO
      ELSE
         DO J=1,NW
            J1 = 2*(J-1)+1
            J2 = J1+1
            CALL PTDIS_4FLV (10,PXT(J1),PYT(J1))
            if (j.eq.1) then
               CALL PTDIS_4FLV (10,PXB(J2),PYB(J2))
            else
               CALL PTDIS_4FLV (IFLB(J2),PXB(J2),PYB(J2))
            endif
            PXB(J1) = -PXB(J2)
            PYB(J1) = -PYB(J2)
            PXT(J2) = -PXT(J1)
            PYT(J2) = -PYT(J1)
         ENDDO
      ENDIF

*        print *,' current NP,QSQ (-1c) ',NP,SQS,NP_0
C...Check consistency of kinematics

      DO J=1,2*NW
         EE(J) = SQS*SQRT(X1(J)*X2(J))
         XM1 = SQRT(PXB(J)**2+PYB(J)**2+QMAS(IABS(IFLB(J)))**2)
         XM2 = SQRT(PXT(J)**2+PYT(J)**2+QMAS(IABS(IFLT(J)))**2)
*        print *,' current NP,QSQ (-1d) ',NP,SQS,NP_0
*        print *,' J,IFLB(J),IFLT(J),NW ',J,IFLB(J),IFLT(J),NW
         IF (EE(J) .LT. XM1+XM2+0.3)  GOTO 2000
      ENDDO

C...Fragmentation of soft/hard sea color strings
      if( NSOF+NJET.gt.0 .and. Irec.eq.1) then
         PAR(1)= PAR(15)
      else
         PAR(1)= PAR(14)
      endif

*     print *,' current NP,SQS (0)',NP,SQS,NP_0
      Ipflag_old= Ipflag
      DO I=1,Nall
        NOLD=NP
        if( Irec .eq. 1 ) then 
           if( I.le.NJET ) then
              Ipflag= Ipflag*100
              Nint= I
           else
              Ipflag= Ipflag*10
              Nint= I-NJET
           endif
        endif
        CALL JET_FRAG (I)
        NNPJET (I) = NP-NOLD
*       print *,' current NP,SQS (1)',NP,SQS,NP_0
        Ipflag= Ipflag_old
        Nint= 0
      ENDDO
c     valence and hard scattering are coupled, reset PAR1 after valence strings!!
c      PAR(1) = PAR1def

C...Fragment the 2*NW valence/sea color strings
      IF(IPAR(4).eq.1)THEN
         PAR(1)= PAR(17)
      ENDIF
      DO JW=1,NW
        if((Irec.eq.1).and.(JDIF(JW).ne.0)) then
          J1 = 2*JW-1
          J2 = J1+1
          X1D = X1(J1)+X1(J2)
          X2D = X2(J1)+X2(J2)
          EE (J1) = SQS*SQRT(X1D*X2D)
          BET(J1) = (X1D-X2D)/(X1D+X2D)
          GAM(J1) = (X1D+X2D)/(2.*SQRT(X1D*X2D))
          if(JW.eq.1) then
            KD = K_beam
          else
            IF(IPAR(18).eq.1) THEN 
               KD = 6
            ELSEIF(IPAR(18).eq.2) THEN  !  diagonal sea states only
               IS = -1 + 2.*INT(1.9999*S_RNDM(0))
               IFL1 = IS*(INT((2.+0.3)*S_RNDM(0))+1)
               IF(IPAR(17).gt.0)THEN
                  CALL SIB_I4FLAV(IFL1,-IFL1,IDUM,KD)
               ELSE
                  CALL SIB_IFLAV(IFL1,-IFL1,IDUM,KD,0)
               ENDIF
            ELSEIF(IPAR(18).eq.3) THEN  
c           allow all sea states, introduces strangeness violation
               IS = -1 + 2.*INT(1.9999*S_RNDM(0))
               IFL1 = IS*(INT((2.+0.3)*S_RNDM(0))+1)
               IF(IPAR(17).gt.0)THEN
                  CALL SIB_I4FLAV(IFL1,0,IDUM,KD)
               ELSE
                  CALL SIB_IFLAV(IFL1,0,IDUM,KD,0)
               ENDIF
            ELSEIF(IPAR(18).eq.4) THEN  
c           use sea from beam_split
               IF(IPAR(17).gt.0)THEN
                  CALL SIB_I4FLAV(IFLB(J1),IFLB(J2),IDUM,KD)
               ELSE
                  CALL SIB_IFLAV(IFLB(J1),IFLB(J2),IDUM,KD,0)
               ENDIF
            ELSE
               KD = 9
            ENDIF
          endif
          Nold = NP
          call SIB_DIFF(KD, JDIF(JW), EE(J1), 0, IREJ)
          if(IREJ.ne.0) 
     &         WRITE(LUN,*)' SIB_NDIFF: SIB_DIFF rejection:',Ncall
          DO K=NOLD+1,NP
            PZ = P(K,3)
            P(K,3) = GAM(J1)*(PZ+BET(J1)*P(K,4))
            P(K,4) = GAM(J1)*(P(K,4)+BET(J1)*PZ)
          ENDDO
          NNPSTR(J1) = NP-Nold
          NNPSTR(J2) = 0
        else
          DO J=2*JW-1,2*JW
            Nint= J
            EE (J) = SQS*SQRT(X1(J)*X2(J))
            BET(J) = (X1(J)-X2(J))/(X1(J)+X2(J))
            GAM(J) = (X1(J)+X2(J))/(2.*SQRT(X1(J)*X2(J)))
            PAR24_def = PAR(24)
            IF(IPAR(15).eq.2) PAR(24) = PAR(25)*EXP(-PAR(26)/EE(J))
            IF(IPAR(15).eq.3) PAR(24) = PAR(25)*EXP(-PAR(26)/EE(J))
            IF(IPAR(15).eq.4) PAR(24) = PAR(25)*EXP(-PAR(26)/EE(J))
            IF(IPAR(15).eq.5) PAR(24) = PAR(25)*EXP(-PAR(26)/EE(J))
            IF(IPAR(15).eq.6) PAR(24) = PAR(25)*EXP(-PAR(26)/EE(J))
            IF(IPAR(15).eq.8) PAR(24) = PAR(25)*EXP(-PAR(26)/EE(J))
            IF(IPAR(15).eq.9) PAR(24) = PAR(25)*EXP(-PAR(26)/EE(J))
            IF(IPAR(15).eq.10) PAR(24) = PAR(25)*EXP(-PAR(26)/EE(J))
            NOLD=NP
            IF(IPAR(38).ne.0.and.Irec.eq.1)THEN
C...  rotate strings instead of attaching all pt to string end hadrons
               PX1 = 0.
               PY1 = 0.
               PX2 = 0.
               PY2 = 0.
c     parton momenta
c     beam side, use off-shell partons
               dpb(1) = PXB(J)
               dpb(2) = PYB(J)
               dpb(3) = 0.5*SQS*X1(J)
c               dpb(5) = QMAS(IABS(IFLB(J)))
               dpb(5) = 0.D0
               dpb(4) = 0.5*SQS*X1(J)
c     target side
               dpt(1) = PXT(J)
               dpt(2) = PYT(J)
               dpt(3) = -0.5*SQS*X2(J)
c               dpt(5) = QMAS(IABS(IFLT(J)))
               dpt(5) = 0.D0
               dpt(4) = 0.5*SQS*X2(J)
c     string momentum in hadron-hadron rest frame
               do jj=1,4
                  dPSTR(jj) = dpb(jj)+dpt(jj)
               enddo
               dPT2STR = dPSTR(1)**2 + dPSTR(2)**2
c               dxmt = dPstr(4)**2-dPstr(3)**2
c               dPstr(5)=dsqrt(dxmt-dpt2str)
               dxmt = (dPstr(4)+dPstr(3))*(dPstr(4)-dPstr(3))
               dPstr(5)=(dsqrt(dxmt)+dsqrt(dPT2str))
     &              *(dsqrt(dxmt)-dsqrt(dPT2str))
               dPstr(5)=dsqrt(dPstr(5))
               EE(j) = dPstr(5)
            ELSE
c     assign pt to hadrons at string end (old model)
               PX1 = PXB(J)
               PY1 = PYB(J)
               PX2 = PXT(J)
               PY2 = PYT(J)
            ENDIF
            IF(IPAR(17).eq.1)THEN
               CALL STRING_FRAG_4FLV
     &      (EE(J),IFLB(J),IFLT(J),PXB(J),PYB(J),PXT(J),PYT(J),IFBAD,1)
            ELSE
c               CALL STRING_FRAG
c     &      (EE(J),IFLB(J),IFLT(J),PXB(J),PYB(J),PXT(J),PYT(J),IFBAD,1)
            CWRN = '3flv routines removed    '
            CALL SIB_ABRT(LUN,CWRN)
            ENDIF
            PAR(24) = PAR24_def
            Nint= 0
            IF (IFBAD .EQ. 1) then
              if(Ndebug.gt.2) 
     &          WRITE(LUN,*)' STRING_FRAG: rejection (Ncall):',Ncall
              GOTO 2000
            ENDIF
            IF(IPAR(38).ge.1.and.Irec.eq.1)THEN
               do jj=1,3
                  dbet(jj)=dpstr(jj)/dpstr(4)
                  dgabe(jj) = dpstr(jj)/dpstr(5)
               enddo
               dgabe(4) = dpstr(4)/dpstr(5) ! gamma
c     energy, momentum check
               IF(ndebug.gt.4)THEN
                  write(LUN,*) '*************************'
                  write(LUN,*)
     &         'total momentum of string final state before rot/boost'
                  CALL PFsum(NOLD+1,NP,Esum,PXsum,PYsum,PZsum,NF)
                  write(LUN,*)'string mass:',dpstr(5)
                  write(LUN,*)'pz:',PZsum
                  write(LUN,*)'px, delta px:',Pxsum,Pxsum-PXT(J)-PXB(J)
                  write(LUN,*)'py, delta py:',Pysum,Pysum-PYT(J)-PYB(J)
                  write(LUN,*)'en:',esum
                  write(LUN,*)'beta',dbet
                  write(LUN,*)'gamma*beta:',dgabe
                  write(LUN,*)'string momentum:',dpstr
               ENDIF
               IF(IPAR(38).eq.2)THEN
c...  rotate string by scattering angle in hadron-hadron rest frame
c     (incorrect since the lorentz transformation does not preserve angles)
                  theta = dacos(dPSTR(3)/dsqrt(dPSTRtot))
                  if(dPSTR(3).lt.0.) theta = PI-theta
                  phi = dacos(dPSTR(1)/dsqrt(dPT2STR))
                  IF(ndebug.gt.4) 
     &            write(LUN,*)'rotating string (theta,phi):',theta,phi
                  CALL SIROBO (NOLD+1,NP,theta,phi,0.D0,0.D0,0.D0)
               ELSEIF(IPAR(38).ge.3)then
c...  rotate string by angles in string rest frame
c     parton momentum in string frame
c     always choose parton in +z direction for angle definition
c                  dgam = gam(j)
                  IF(ndebug.gt.4) 
     &                 write(6,*)'initial 4 vector of beam parton:',dpb
                  call SIB_ALTRA(dgabe(4),-dgabe(1),-dgabe(2),-dgabe(3),
     &                 dpb(1),dpb(2),dpb(3),dpb(4),
     &                 dPTOT1,PC(1),PC(2),PC(3),PC(4))
                  IF(ndebug.gt.4) then
                     write(6,*)'boosted 4 vector (p,px,py,pz,e):',
     &                    dPTOT1,PC
                     print*,'target parton 4 momentum for cross check:',
     &                    dpt
                     print*,'boosting 4 vector..'
                     call SIB_ALTRA(dgabe(4),-dgabe(1),-dgabe(2),
     &                    -dgabe(3),dpt(1),dpt(2),dpt(3),dpt(4),
     &                    dP1TOT_t,pct(1),pct(2),pct(3),pct(4))
                     print*,'boosted 4 vector (p,px,py,pz,e):',
     &                    dP1TOT_t,pct
                  endif
c     rotation factors
                  COD= PC(3)/dPTOT1
                  SID= DSQRT(PC(1)**2+PC(2)**2)/dPTOT1
                  COF=1.D0
                  SIF=0.D0
                  IF(dPTOT1*SID.GT.1.D-5) THEN
                     COF=PC(1)/(SID*dPTOT1)
                     SIF=PC(2)/(SID*dPTOT1)
                     ANORF=DSQRT(COF*COF+SIF*SIF)
                     COF=COF/ANORF
                     SIF=SIF/ANORF
                  ENDIF
                  if(ndebug.gt.4)
     &                 print*,'rotation factors: cth,sth,cph,sph',
     &                 cod,sid,cof,sif
c     rotate string final state
                  DO K=NOLD+1,NP
                     do ii=1,3
                        dP0(ii)=P(K,ii)
                     enddo
                     call SIB_TRANI(dP0(1),dP0(2),dP0(3),cod,sid,cof,sif
     &                    ,dP1(1),dP1(2),dP1(3))
                     do ii=1,3
                        P(K,ii)=dP1(ii)
                     enddo
                  ENDDO
               ENDIF
c...  boost string along x,y,z into hadron-hadron restframe
c     sign depends on z component of parton-parton 4 momentum
               DO K=NOLD+1,NP
                  IF(IPAR(38).eq.4)then
c...  ralph's double prec. boost routine
                     do ii=1,4
                        dP0(ii)=P(K,ii)
                     enddo
                     call SIB_ALTRA(dgabe(4),dgabe(1),dgabe(2),
     &                    dgabe(3),dp0(1),dp0(2),dp0(3),dp0(4),
     &                    dP1TOT_t,dp1(1),dp1(2),dp1(3),dp1(4))
                     do ii=1,4
                        P(K,ii)=dP1(ii)
                     enddo
                  else
c...  not so smart boost
                     BETP = dbetz*P(k,3) + dbety*P(K,2) + dbetx*P(K,1)
                     do ii=1,3
                        P(K,ii)=gam(j)*P(K,ii)+dgabe(ii)*P(K,4)
                     enddo
                     P(K,4) = GAM(J)*(P(K,4)+BETP)
                  ENDIF
               ENDDO
            else
c...  boost along z into hadron-hadron restframe 
c     (parton pt is put into string end hadrons)
               DO K=NOLD+1,NP
                  PZ = P(K,3)
                  P(K,3) = GAM(J)*(PZ+BET(J)*P(K,4))
                  P(K,4) = GAM(J)*(P(K,4)+BET(J)*PZ)
               ENDDO
            endif
            NNPSTR(J) = NP-NOLD
          ENDDO
        endif
      ENDDO
      PAR(1)= PAR1def 
      IREJ = 0
      SQS   = SQS_0
      S     = SQS*SQS

      if(Ndebug.gt.2)
     &  WRITE(LUN,*)'SIB_NDIFF: generated interactions (Ns,Nh):',
     &  NSOF+NW,NJET

      END


      SUBROUTINE SIBNUC (IAB, IAT, SQS)
C-----------------------------------------------------------------------
C.  Routine that generates the interaction of a nucleus of
C.  mass number IAB with a  target nucleus  of mass IAT
C.  (IAT=0 : air).
C.  SQS (GeV) is the  center of mass energy of each
C.  nucleon - nucleon cross section
C-----------------------------------------------------------------------
      SAVE
      COMMON /S_PLIST/ P(8000,5), LLIST(8000), NP
      COMMON /S_PLNUC/ PA(5,40000), LLA(40000), NPA
      COMMON /S_MASS1/ AM(99), AM2(99)
      COMMON /CKFRAG/ KODFRAG
      PARAMETER (IAMAX=56)
      COMMON /CNUCMS/ B, BMAX, NTRY, NA, NB, NI, NAEL, NBEL
     +         ,JJA(IAMAX), JJB(IAMAX), JJINT(IAMAX,IAMAX)
     +         ,JJAEL(IAMAX), JJBEL(IAMAX)            
      COMMON /FRAGMENTS/ PPP(3,60)
      DIMENSION SIGDIF(3)
      DIMENSION IAF(60)
      DATA RPOX /0.3624/

C...Target mass
      IF (IAT .EQ. 0) THEN
         IATARGET = 14 + 2*INT((1.+RPOX)*S_RNDM(0))
      ELSE
         IATARGET = IAT
      ENDIF
       
C...Single nucleon (proton) case

      IF (IAB .EQ. 1)  THEN
         NPA = 0
         CALL SIBYLL (13,IATARGET, SQS)
         CALL DECSIB
         DO J=1,NP
            LA = IABS(LLIST(J))
            IF (LA .LT. 10000)  THEN
               NPA = NPA + 1
               LLA(NPA) = LLIST(J)
               DO K=1,5
                  PA(K,NPA) = P(J,K)
               ENDDO
            ENDIF
         ENDDO
         RETURN
      ENDIF


C...Nuclei

      CALL SIB_SIGMA_HP(1,SQS,SIGT,SIGEL,SIG0,SIGDIF,SLOPE,RHO)
      CALL INT_NUC (IATARGET, IAB, SIG0, SIGEL) 

C...fragment spectator nucleons
      NBT = NB + NBEL
      IF (KODFRAG .EQ. 1)  THEN
          CALL FRAGM1(IAB,NBT, NF, IAF)
      ELSE IF(KODFRAG .EQ. 2)  THEN
          CALL FRAGM2(IAB,NBT, NF, IAF)
      ELSE 
          CALL FRAGM (IATARGET, IAB, NBT,B, NF, IAF)
      ENDIF
     
C...Spectator fragments
      NPA = 0
      DO J=1,NF
         NPA = NPA+1
         if(NPA.gt.40000) then
           write(6,'(1x,a,2i8)') 
     &       'SIBNUC: no space left in S_PLNUC (NPA,NF)',NPA,NF
           NPA = NPA-1
           return
         endif
         LLA(NPA) = 1000+IAF(J)
         PA(1,NPA) = 0.
         PA(2,NPA) = 0.
         PA(3,NPA) = SQS/2.
         PA(4,NPA) = SQS/2.
         PA(5,NPA) = FLOAT(IAF(J))*0.5*(AM(13)+AM(14))
      ENDDO

C...Elastically scattered fragments
      DO J=1,NBEL
         NPA = NPA+1
         if(NPA.gt.40000) then
           write(6,'(1x,a,2i8)') 
     &       'SIBNUC: no space left in S_PLNUC (NPA,NBEL)',NPA,NBEL
           NPA = NPA-1
           return
         endif
         LLA(NPA) = 1001
         PA(1,NPA) = 0.
         PA(2,NPA) = 0.
         PA(3,NPA) = SQS/2.
         PA(4,NPA) = SQS/2.
         PA(5,NPA) = 0.5*(AM(13)+AM(14))
      ENDDO

C...Superimpose NB  nucleon interactions
      DO JJ=1,NB
          CALL SIBYLL (13,IATARGET, SQS)
          CALL DECSIB
          DO J=1,NP
             LA = IABS(LLIST(J))
             IF (LA .LT. 10000)   THEN
                NPA = NPA + 1
                if(NPA.gt.40000) then
                  write(6,'(1x,a,2i8)') 
     &              'SIBNUC: no space left in S_PLNUC (NPA,NP)',NPA,NP
                  NPA = NPA-1
                  return
                endif
                LLA(NPA) = LLIST(J)
                DO K=1,5
                    PA(K,NPA) = P(J,K)
                ENDDO
            ENDIF
         ENDDO
      ENDDO

      END



      FUNCTION CHIDIS (KPARTin, IFL1, IFL2)
C...Generate CHI (fraction of energy of a hadron carried by 
C.                the valence quark, or diquark, as specified by IFL1)
C.  INPUT KPART = code of particle
C.        IFL1, IFL2 = codes of partons (3, 3bar of color)
C.........................................................
      SAVE
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, kb ,kt
      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      COMMON /S_CPSPL/ CCHIK(3,39)
      COMMON /S_CUTOFF/ STR_mass_val, STR_mass_val_hyp, STR_mass_sea

      kpart=IABS(kpartin)
      IFQ=IABS(IFL1)
      IF (IFQ.GT.10) IFQ=IABS(IFL2)
      CUT=2.*STR_mass_val/SQS
c     hyperon beam cut
      IF(kpart.gt.14) CUT=2.*STR_mass_val_hyp/SQS
100   CHIDIS=S_RNDM(0)**2
      if (chidis.lt.cut) goto 100
      if (chidis.gt.(1.-cut)) goto 100
      IF((CHIDIS**2/(CHIDIS**2+CUT**2))**0.5
     +   *(1.-CHIDIS)**CCHIK(IFQ,KPART).LT.S_RNDM(0)) GOTO 100
      CHIDIS = MAX(0.5*CUT,CHIDIS)
      CHIDIS = MIN(1.-CUT,CHIDIS)
      IF (IABS(IFL1).GT.10)  CHIDIS=1.-CHIDIS
      RETURN
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
C.   \AF'13
C------------------------------------------------
      SAVE
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)

      DIMENSION KPL(3)
      DATA KPL /1,2,3/

      Select case(IABS(KF))

      Case (6)                  ! pi0
         R = S_RNDM(0)              
         XBUG = 0.
         IF(IPAR(19).eq.1) XBUG = 0.5
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
         IF (S_RNDM(0).GT. 0.5)  THEN
            KP1 = 3
            KP2 = -2
         ENDIF
      Case (13)                 ! p/pbar
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
C...  Test for good input
         WRITE(6,*)
     &        'HSPLI : Routine entered with illegal particle code ', KF
         KP1 = 0
         KP2 = 0

      End Select

C     if anti-baryon, invert valences
      IF (KF .LT. 0) THEN
         KPP = KP1
         KP1 = -KP2
         KP2 = -KPP
      ENDIF
      RETURN
      END



      SUBROUTINE JET_FRAG (Index)
C-----------------------------------------------------------------------
C.   Fragmentation of a jet-jet system
C.   Input : kinematical variables of a jet-jet system, 
C.           taken from /S_CHIST/
C-----------------------------------------------------------------------
      SAVE

      REAL*8 DX1J, DX2J, DBETJ
      COMMON /S_PLIST/ P(8000,5), LLIST(8000), NP
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, kb ,kt
      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      PARAMETER (NW_max = 20)
      PARAMETER (NS_max = 20, NH_max = 80)
      PARAMETER (NJ_max = (NS_max+NH_max)*NW_max)
      COMMON /S_CHIST/ X1J(NJ_max),X2J(NJ_max),
     &    X1JSUM(NW_max),X2JSUM(NW_max),PTJET(NJ_max),PHIJET(NJ_max),
     &    NNPJET(NJ_max),NNPSTR(2*NW_max),NNSOF(NW_max),NNJET(NW_max),
     &    JDIF(NW_max),NW,NJET,NSOF
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)
      DATA PGG /1./
      CHARACTER*25 CWRN

      if(Ndebug.gt.2) then
        WRITE(LUN,*)' JET_FRAG: called for entry (I,NP):',Index,NP
        WRITE(LUN,*)' JET_FRAG: (X1J,X2J,PTjet):',X1J(Index),X2J(Index),
     &    PTjet(Index)
      endif

      E0 = SQRT(S*X1J(Index)*X2J(Index))
      TH = ASIN(MIN(0.999999,2.*PTJET(Index)/E0))
      FI = 6.283185*S_RNDM(0)

C...  charm setup
      PAR24_def = PAR(24)
      SELECT CASE(IPAR(15))
      CASE(2,3)
         PAR(24) = PAR(25)*EXP(-PAR(26)/E0)
      CASE(4)
         PAR(24) = PAR(27)*EXP(-PAR(26)/E0)
      CASE(5)
         PAR(24) = PAR(27)*EXP(-PAR(26)/E0)
         PAR(29) = PAR(27)*EXP(-PAR(28)/E0)
      CASE(6,8,9,11) 
         PAR(24) = PAR(27)*EXP(-PAR(28)/E0)
      CASE(7)
         PAR(24) = PAR(27)
      CASE(10)
         IF(INDEX.LT.NJET+1) THEN
c     only increase perturbative charm production
            PAR(24) = PAR(27)*EXP(-PAR(28)/E0) 
         ELSE
c     non pert. parameters
            PAR(24) = PAR(25)*EXP(-PAR(26)/E0) 
         ENDIF
      END SELECT

C...  strange setup
      PAR2_1_def = PAR(2)
      PAR3_def = PAR(3)
      SELECT CASE(IPAR(42))
c     change to constant value 
      CASE(1)
         PAR(2) = PAR(72)
c     change according to string mass, saturating
      CASE(2)
         PAR(2) = PAR(72)*EXP(-PAR(73)/E0)
c     change strange diq fraction as well
      CASE(3)
         PAR(2) = PAR(72)       ! P_s / P_ud
         PAR(3) = PAR(73)       ! P_us / P_ud
      END SELECT

C...  vector setup
      PAR5_def = PAR(5)
      PAR6_def = PAR(6)
      SELECT CASE (IPAR(43))
c     change vector rate and kaon vector rate
      CASE(1)    
         PAR(5) = PAR(74)       ! P_vec
         PAR(6) = PAR(74)       ! P_K* from K
      END SELECT

      NOLD = NP
      IF ( (E0.LT.8.) .OR. (S_RNDM(0).GT.PGG)) THEN
C...  'leading' strange fraction
         PAR2_2_def = PAR(2)
         IF(IPAR(39).eq.2) PAR(2) = PAR(66)
         IS = -1 + 2.*INT(1.9999*S_RNDM(0))
 100     IFL1 = IS*(INT((2.+PAR(2))*S_RNDM(0))+1)
         XM = 2.*QMASS(IFL1)+0.3
         if(E0.LE.XM) GOTO 100
         PAR(2) = PAR2_2_def
         IF(IPAR(17).eq.1)THEN
            IF(IABS(IFL1).eq.3)THEN
               IF(S_RNDM(IFL1).lt.PAR(24))IFL1 = IS*4
               XM = 2.*QMASS(IFL1)+0.3
               if(E0.LE.XM) GOTO 100
            ENDIF
            CALL STRING_FRAG_4FLV (E0,IFL1,-IFL1,0.,0.,0.,0.,IFBAD,0)
         ELSE
c            CALL STRING_FRAG (E0,IFL1,-IFL1,0.,0.,0.,0.,IFBAD,0)
            CWRN = '3flv routines removed    '
            CALL SIB_ABRT(CWRN)
         ENDIF
         if(IFBAD.ne.0) WRITE(LUN,*)
     &     ' JET_FRAG: rejection in STRING_FRAG (IFL,E0):',IFL1,E0,
     &     ' strg.ratio,charm ratio:', PAR(2),PAR(24)
      ELSE
         IF(IPAR(17).eq.1)THEN
            CALL GG_FRAG_4FLV(E0)
         ELSE
            CWRN = '3flv routines removed    '
            CALL SIB_ABRT(LUN,CWRN)
c            CALL GG_FRAG(E0) 
         ENDIF
      ENDIF
      DX1J = X1J(Index)
      DX2J = X2J(Index)
      DBETJ = (DX1J-DX2J)/(DX1J+DX2J)
      CALL SIROBO (NOLD+1,NP,TH,FI,0.D0,0.D0,DBETJ)

      if(Ndebug.gt.2) 
     &     WRITE(LUN,*)' JET_FRAG: particles produced:',NP-NOLD
      PAR(24) = PAR24_def
      PAR(2) = PAR2_1_def
      PAR(5) = PAR5_def
      PAR(6) = PAR6_def
      PAR(3) = PAR3_def 
      END


      SUBROUTINE STRING_FRAG_4FLV
     +     (E0,IFL1,IFL2,PX1,PY1,PX2,PY2,IFBAD,IFQRK)
C-----------------------------------------------------------------------
C.  This routine fragments a string of energy E0
C.  the ends of the strings  have flavors IFL1 and IFL2
C.  the particles produced are in the  jet-jet frame
C.  with IFL1 going in the +z direction
C.     E0 = total energy in jet-jet system
C.  This version consider also a primordial pT attached
C.  to the ends of the string PX1,PY1,  PX2,PY2
C.  OUTPUT:  IFBAD =1  kinematically impossible decay
c	2010.03.11 ifqrk - leading quark flag
c	1 in valence quark, 0 in others
c
c      Modified Nov. 91.  RSF and TSS to fragment symmetrically
c      ie forward and backward are fragmented as leading.
c      Change- Dec. 92  RSF.  call to ptdis moved- to use flavor
c      of NEW quark in fragmentation.
c
c     includes 4 FLAVORS \FR'13
C-----------------------------------------------------------------------
      SAVE

      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      COMMON /S_PLIST/ P(8000,5), LLIST(8000), NP
      COMMON /S_MASS1/ AM(99), AM2(99)
      DIMENSION WW(2,2), PTOT(4), PX(3),PY(3),IFL(3),ILEAD(2)
      DIMENSION LPOINT(8000), PMQ(3)
      LOGICAL LRANK
      DATA LRANK/.true./
      COMMON /S_PARTO/ NFORIG(8000), NPORIG(8000), IPFLAG, NINT
      COMMON /S_ZLIST/ ZLIST(8000)
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)
      COMMON /S_CZDIS/ FAin, FB0in
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, kb ,kt
      
      IF(Ndebug.gt.2) THEN
         WRITE(LUN,*)
     &        ' STRING_FRAG_4FLV: called with ',
     &        '(E0,IFL1,IFL2,PX1,PY1,PX2,PY2,IVAL)',
     &        E0,IFL1,IFL2,PX1,PY1,PX2,PY2,IFQRK
         WRITE(LUN,*)' STRING_FRAG_4FLV: NP before fragmentation:',NP
      ENDIF

c...  remember initial values
c     strange fraction
      par2_def = PAR(2)
c     vector model
      IPAR11_def = IPAR(11)
c     vector fraction
      PAR5_def = PAR(5)      
c     charm fraction
      PAR24_def = PAR(24)
c     popcorn fraction
      PAR8_def = PAR(8)
      
C...initialise
      NTRY = 0
      IFBAD = 0
 200  NTRY = NTRY + 1

c     reset parameters after rejection
      PAR(2) = PAR2_def
      PAR(5) = PAR5_def
      PAR(24) = PAR24_def
      IPAR(11) = IPAR11_def
      PAR(8) = PAR8_def

      IF (NTRY .GT. 50)  THEN
         IFBAD = 1
         RETURN
      ENDIF
      I = NP
      DO K=1,2
         WW(K,1) = 1.
         WW(K,2) = 0.
      ENDDO
      PX(1) = PX1
      PY(1) = PY1
      PX(2) = PX2
      PY(2) = PY2
      PX(3) = 0.
      PY(3) = 0.
      PTOT (1) = PX1+PX2
      PTOT (2) = PY1+PY2
      PTOT (3) = 0.
      PTOT (4) = E0
      IFL(1) = IFL1
      IFL(2) = IFL2
      PMQ(1) = QMASS(IFL(1))
      PMQ(2) = QMASS(IFL(2))

      ILEAD(1) = 0
      ILEAD(2) = 0
      IBLEAD = 0
      IF(IABS(IFQRK).eq.1) THEN
         ILEAD(1) = 1
         ILEAD(2) = 1
      ENDIF
      IF(IPAR(20).eq.0.and.IFQRK.gt.0) GOTO 300
C
C      SET FLAG FOR GENERATION OF LEADING PARTICLES. 
C      "AND" IS FOR PPBAR ( DIQUARK AT BOTH ENDS)
C      "OR" IS FOR PP, PPI, ( DIQUARK AT ONE END.)
C
      IF (IABS(IFL1) .GT. 10 .AND. IABS(IFL2) .GT. 10)  THEN
         IBLEAD = 2
         I = I+1
         JT = 1.5+S_RNDM(0)
         GOTO 350
      ENDIF         
      IF (IABS(IFL1) .GT. 10 .OR. IABS(IFL2) .GT. 10)  THEN
         IBLEAD = 1
         I = I+1
         JT = 2
         IF (IABS(IFL2) .GT. 10) JT = 1
         GOTO 350
      ENDIF         

C...produce new particle: side, pT
 300  continue
      I=I+1
      if(i.gt.8000) then
        write(LUN,'(1x,a,i8)') 
     &        'STRING_FRAG_4FLV: no space left in S_PLIST:',I
        write(LUN,'(1x,a,i8)') 
     &       ' STRING_FRAG_4FLV: called with ',
     &       '(E0,IFL1,IFL2,PX1,PY1,PX2,PY2,IVAL)',
     &       E0,IFL1,IFL2,PX1,PY1,PX2,PY2,IFQRK
        stop   
      endif
      IF (IBLEAD .GT. 0)  THEN
         JT = 3 - JT   
         GO TO 350              
      ENDIF
c     
 349  continue
      JT=1.5+S_RNDM(0)                            
 350  JR=3-JT
      LPOINT(I) = JT

      nporig(I)= Ipflag*2 + Nint
      IF(ILEAD(JT).eq.1) nporig(I)= -1 * nporig(I)
      nforig(I) = 0

 555  CONTINUE
c
c.... CHARM config
c
      charmPARdef=PAR(24)
      IF(IPAR(15).lt.9)THEN
c     no s->c
         PAR(24) = 0.
         IF (IFQRK.EQ.1) THEN
c     ifqrk = 1 (valence quark attatched) 
            IF(IPAR(15).ge.1) THEN
c     enforce s->c at string end
               IF(ILEAD(JT).eq.1) PAR(24)=charmPARdef
c     produce charm in all strings
               IF(IPAR(15).eq.8) PAR(24)=charmPARdef
            ELSE
c     compatibility to broken version
               PAR(24)=charmPARdef
            ENDIF
         ELSE
c     no val. quark at string end or diff
            PAR(24)=charmPARdef
         ENDIF
      ENDIF
c
C.... Vector meson config
c
c...  switch off for proton beam
      IF(IPAR(31).eq.1)then
c         print*,'ipar11,ipar11def,1-kb/13,kb',ipar(11),ipar11_def,
c     +        max((1-iabs(kb)/13),0),kb
         IPAR(11) = IPAR(11)*max((1-iabs(kb)/13),0) ! meson beam only
      endif
c     increase vec.meson ratio for leading quarks
      IF(IABS(IFQRK).eq.1)THEN
         IF(IPAR(11).le.-5.and.IPAR(11).ge.-7
     &        .and.ilead(jt).eq.1)
     &        PAR(5) = 9.
         
c     increase vec.meson ratio for diff.
         IF(IFQRK.eq.-1.and.IPAR(11).le.-4.and.IPAR(11).ge.-7)
     &        PAR(5) = 9.

c     increase vec.meson ratio for leading particle in str. diff. (lvec16)
         IF(IFQRK.eq.-1.and.IPAR(11).le.-11.and.ILEAD(JT).EQ.1)
     &        PAR(5) = 99.
      ENDIF

c...  suppress leading charm for pion and kaon beams
      IF(IPAR(15).eq.11)then
         IF((1-IABS(KB)/13)*ILEAD(JT).gt.0) PAR(24)=0.
      ENDIF

C...  suppress rank-1 baryon through popcorn
      IF(IBLEAD .GT. 0.and.abs(ifl(jt)).gt.10
     &     .and.abs(ifl(3)).lt.10) PAR(8)=PAR(63)*PAR(8)

C...  leading strange/charm
      IF(ILEAD(JT).eq.1.and.IPAR(39).gt.0) PAR(2) = PAR(65)

C..   scale valence string end charm for assoc. prod.      
      IF(IPAR(41).eq.1)THEN
         IF(ILEAD(JT).eq.1.and.IFQRK.eq.1) PAR(24) = PAR(71)*PAR(24)
      ENDIF


C...particle ID and pt.

      CALL SIB_I4FLAV (IFL(JT), 0, IFL(3), LLIST(I))

c     reset strange fraction
      PAR(2) = PAR2_def
c     reset vec.meson production
      PAR(5) = PAR5_def
c     reset charm fraction
      PAR(24) = PAR24_def
c     reset popcorn
      PAR(8) = par8_def

c     replace leading pi0 by rho0's
      IF(IABS(IFQRK).eq.1)THEN
         IF(ABS(IPAR(11)).ge.2.and.IPAR(11).ge.-3)THEN
            IF(ilead(jt).EQ.1) then 
               IF(ABS(LLIST(I)).EQ.6) THEN
                  LLIST(I) = 27*isign(1,LLIST(I))
               endif
            endif
        
c     replace leading pi0 in string diff by rho0's (lvec15)
         ELSEIF(IFQRK.eq.-1.and.IPAR(11).eq.-10)THEN
            IF(ILEAD(JT).EQ.1) THEN 
               IF(ABS(LLIST(I)).EQ.6) THEN
                  LLIST(I) = 27*isign(1,LLIST(I))
               ENDIF
            ENDIF
c     replace leading pi0 in string diff by rho0's 
c     in addition to increased leading vec.meson ratio (lvec20)
         ELSEIF(IFQRK.eq.-1.and.IPAR(11).eq.-15)THEN
            IF(ILEAD(JT).EQ.1) THEN 
               IF(ABS(LLIST(I)).EQ.6) THEN
                  LLIST(I) = 27*isign(1,LLIST(I))
               ENDIF
            ENDIF     
c     replace leading omega in string diff by rho0's 
c     in addition to increased leading vec.meson ratio (lvec21)
         ELSEIF(IFQRK.eq.-1.and.IPAR(11).eq.-16)THEN
            IF(ILEAD(JT).EQ.1) THEN 
               IF(ABS(LLIST(I)).EQ.32) 
     &              LLIST(I) = 27*isign(1,LLIST(I))
            ENDIF     
c     replace leading omega in string diff by rho0's 
c     suppress pi0 in diff. strings
c     in addition to increased leading vec.meson ratio (lvec22)
         ELSEIF(IFQRK.eq.-1.and.IPAR(11).eq.-17)THEN
            IF(ILEAD(JT).EQ.1) THEN 
c     print*,'replacing leading omega with rho0'
               IF(ABS(LLIST(I)).EQ.32)
     &              LLIST(I) = 27*isign(1,LLIST(I))
            ENDIF
            IF(LLIST(I).EQ.6) then
c     print*,'pi0 found! start again.. '
               GOTO 555
            endif

c     replace all for diff.
         ELSEIF(IFQRK.eq.-1.and.ipar(11).lt.0.and.
     &           ipar(11).ge.-3) then
            IF(ABS(LLIST(I)).EQ.6)  LLIST(I) = 27*isign(1,LLIST(I))

c     increased vec.meson ratio and replace pi0 with rho0 in str.diff
         ELSEIF(IFQRK.eq.-1.and.ipar(11).eq.-7) then
            IF(ABS(LLIST(I)).EQ.6)  LLIST(I) = 27*isign(1,LLIST(I))  

c     replace leading pi's by vec.mesons, iso-spin conserving
         ELSEIF(IPAR(11).eq.-8.and.IPAR(11).lt.0)THEN
            PAR(5) = 9.
            IF(ilead(jt).EQ.1.and.INT((PAR(5)+1)*S_RNDM(0)).gt.1) then 
               IF(ABS(LLIST(I)).EQ.6) LLIST(I) = 27*isign(1,LLIST(I))
               IF(ABS(LLIST(I)).EQ.7) LLIST(I) = 25*isign(1,LLIST(I))
c     IF(ABS(LLIST(I)).EQ.8) LLIST(I) = 26*isign(1,LLIST(I))
            endif

c     replace almost all for diff.
         ELSEIF(IFQRK.eq.-1.and.ipar(11).eq.-8.and.ipar(11).lt.0) then
            PAR(5) = 9.
            if( INT((PAR(5)+1)*S_RNDM(0)).gt.1 ) then
               IF(ABS(LLIST(I)).EQ.6)  LLIST(I) = 27*isign(1,LLIST(I))
               IF(ABS(LLIST(I)).EQ.7) LLIST(I) = 25*isign(1,LLIST(I))
            endif
      
c     replace leading pi0's by vec.mesons
         ELSEIF(IPAR(11).eq.-9.and.IPAR(11).lt.0)THEN
            PCHF = 0.1
            IF(ilead(jt).EQ.1.and.ABS(LLIST(I)).EQ.6) 
     &           LLIST(I) = 27*isign(1,LLIST(I))
            if(ilead(jt).EQ.1.and.ABS(LLIST(I)).EQ.7)then
               if(S_RNDM(0).lt.PCHF) LLIST(I) = 25*isign(1,LLIST(I))
            endif        

c     replace for string diff.
         ELSEIF(IFQRK.eq.-1.and.ipar(11).eq.-9) then
            IF(ABS(LLIST(I)).EQ.6) 
     &           LLIST(I) = 27*isign(1,LLIST(I))
            if(ABS(LLIST(I)).EQ.7)then
               if(S_RNDM(0).lt.PCHF) 
     &              LLIST(I) = 25*isign(1,LLIST(I))
            endif
         ELSE
            CONTINUE
         ENDIF
      ENDIF

c     reset vec.meson ratio
      PAR(5) = 0.3
      IF(IABS(IFQRK).eq.1) ILEAD(JT) = 0
      
      PMQ(3) = QMASS(IFL(3))
      P(I,5) = AM(IABS(LLIST(I)))
      CALL PTDIS_4FLV (IFL(3), PX(3),PY(3))

C...fill transverse momentum
      P(I,1) = PX(JT) + PX(3)
      P(I,2) = PY(JT) + PY(3)
      XMT2 = P(I,5)**2+P(I,1)**2+P(I,2)**2

C...test end of fragmentation

      WREM2 = PTOT(4)**2-PTOT(1)**2-PTOT(2)**2-PTOT(3)**2
c      IF (WREM2 .LT. 0.1)  GOTO 200
      IF (WREM2 .LT. 0.1)  GOTO 200
      WMIN = PMQ(1)+PMQ(2)+2.*PMQ(3)+ PAR(59) + (2.*S_RNDM(0)-1.)*0.2
      IF (WREM2 .LT. WMIN**2)    Then		!   goto 400
         if (abs(ifl(3)).ne.3.and.ABS(IFL(3)).ne.4) GOTO 400
         goto 200
      endif

C...Choose z
      IF(IABS(IFQRK).eq.1) THEN
c     valence strings: ( str.diff and non diff. )
         IF(IPAR(11).EQ.1) THEN
c     use hard distribution for leading quarks ( no exchange )
            IF(ILEAD(JT).eq.1) THEN
               Z = ZBLEAD (IABS(LLIST(I)))
            ELSE
               IF(IBLEAD .GT. 0.and.abs(ifl(jt)).gt.10)  THEN
                  Z = ZBLEAD (IABS(LLIST(I)))   
                  IBLEAD = IBLEAD - 1
               ELSE
                  Z = ZDIS_4FLV (IFL(3),ifl(jt),XMT2)
               ENDIF
            ENDIF
c     use hard frag. for leading particles
         ELSEIF(IPAR(11).ge.3.or.IPAR(11).eq.-3.or.IPAR(11).eq.-6
     &           .or.IPAR(11).eq.-7) THEN
            IF(ILEAD(jt).eq.1) THEN
               Z = ZBLEAD (IABS(LLIST(I)))
            ELSE
               IF(IBLEAD .GT. 0.and.abs(ifl(jt)).gt.10)  THEN
                  Z = ZBLEAD (IABS(LLIST(I)))   
                  IBLEAD = IBLEAD - 1
               ELSE
                  Z = ZDIS_4FLV (IFL(3),ifl(jt),XMT2)
               ENDIF
            ENDIF
         ELSEIF(IPAR(11).EQ.-11) THEN
c     very hard leading frag. for diff and non. diff val. strings (lvec16)
            IF(IBLEAD .GT. 0.and.abs(ifl(jt)).gt.10) THEN
               Z = ZBLEAD (IABS(LLIST(I)))
            ELSEIF(ILEAD(jt).eq.1)THEN
               Z = 1. - ZDISN(1)
            ELSE
               Z = ZDIS_4FLV (IFL(3),ifl(jt),XMT2)
            ENDIF

         ELSEIF(IPAR(11).EQ.-12.OR.IPAR(11).LE.-15)THEN
c     very hard leading frag. for diff. val. strings only (lvec17)
            IF(IBLEAD .GT. 0.and.abs(ifl(jt)).gt.10) THEN
               Z = ZBLEAD (IABS(LLIST(I)))
            ELSEIF(ILEAD(jt).eq.1.and.IFQRK.eq.-1)THEN
               Z = 1. - ZDISN(1)
            ELSE
               Z = ZDIS_4FLV (IFL(3),ifl(jt),XMT2)
            ENDIF

         ELSEIF(IPAR(11).EQ.-13.AND.IFQRK.eq.-1) THEN
c     hard leading frag. for diff. val. strings only (lvec18)
            IF(IBLEAD .GT. 0.and.abs(ifl(jt)).gt.10) THEN
               Z = ZBLEAD (IABS(LLIST(I)))
            ELSEIF(ILEAD(jt).eq.1)THEN
               Z = S_RNDM(JT)
            ELSE
               Z = ZDIS_4FLV (IFL(3),ifl(jt),XMT2)
            ENDIF
         ELSEIF(IPAR(11).EQ.-14.AND.IFQRK.eq.-1) THEN
c     hard leading frag. for diff. AND ndiff. val. strings (lvec19)
            IF(IBLEAD .GT. 0.and.abs(ifl(jt)).gt.10) THEN
               Z = ZBLEAD (IABS(LLIST(I)))
            ELSEIF(ILEAD(jt).eq.1)THEN
               Z = S_RNDM(JT)
            ELSE
               Z = ZDIS_4FLV (IFL(3),ifl(jt),XMT2)
            ENDIF
            
         ELSE

c     hard leading baryons only ( standard )
            IF(IBLEAD .GT. 0.and.abs(ifl(jt)).gt.10)  THEN
               IF(IPAR(20).eq.3)THEN
c     use lund function with different parameters for leading baryon
                  fa_def = FAin
                  fb_def = FB0in
                  FAin = PAR(57)
                  FB0in = PAR(58)
                  Z = ZDIS_4FLV(IFL(3),IFL(JT),XMT2)
c     set parameters to initial values again
                  FAin = fa_def
                  FB0in = fb_def
               ELSE
                  Z = ZBLEAD (IABS(LLIST(I)))
               ENDIF
               IBLEAD = IBLEAD - 1
            ELSE
               Z = ZDIS_4FLV (IFL(3),ifl(jt),XMT2)
            ENDIF
         ENDIF
      ELSE
c     non valence string
         IF (IBLEAD .GT. 0.and.abs(ifl(jt)).gt.10)  THEN
C     Special frag. for leading Baryon only
            Z = ZBLEAD (IABS(LLIST(I)))   
            IBLEAD = IBLEAD - 1
         ELSE
            Z = ZDIS_4FLV (IFL(3),ifl(jt),XMT2)
         ENDIF
      ENDIF
      IF(IFQRK.eq.1) ILEAD(JT) = 0

      WW(JT,2) = Z*WW(JT,1)
      WW(JR,2) = XMT2/(WW(JT,2)*E0**2)

      P(I,3) = WW(1,2)*0.5*E0 - WW(2,2)*0.5*E0
      P(I,4) = WW(1,2)*0.5*E0 + WW(2,2)*0.5*E0

      DO J=1,4
         PTOT (J) = PTOT(J) - P(I,J)
      ENDDO
      DO K=1,2
         WW(K,1) = WW(K,1) - WW(K,2)
      ENDDO

C...Reset pT and flavor at ends of the string
      PX(JT) = -PX(3)
      PY(JT) = -PY(3)
      IFL(JT) =-IFL(3)
      PMQ(JT) = PMQ(3)

      GOTO 300

C...Final two hadrons
 400  IAFL1 = IABS(mod(IFL(JR),100))
      IAFL2 = IABS(mod(IFL(3),100))
      IF(IPAR(40).eq.0)THEN
c     reject two diquarks, two anti-diquarks AND diquark anti-diquark pairs
         IF (IAFL1*IAFL2 .GT. 100)  GOTO 200 
      ELSE
c     ONLY reject two diquarks or two anti-diquarks (unphysical) 
c     AND KEEP diquark anti-diquark pairs 
         IF (mod(IFL(JR),100)*mod(IFL(3),100).GT.100) GOTO 200 
      ENDIF

      IF ((IAFL1/10.eq.4.or.mod(IAFL1,10).eq.4)
     +     .and.(IAFL2/10.eq.4.or.mod(IAFL2,10).eq.4))
     +     GOTO 200             ! reject two charm quarks

C.... Vector meson config
c     increase vec.meson ratio for diff.
      IF(IFQRK.eq.-1.and.IPAR(11).le.-4.and.IPAR(11).gt.-8) PAR(5) = 9.
c     increase vec.meson ratio for leading quarks in valence interactions
      IF(IABS(IFQRK).eq.1.and.IPAR(11).le.-5.and.ilead(jr).eq.1
     &     .and.IPAR(11).gt.-8) PAR(5) = 9.
      
 666  CALL SIB_I4FLAV (IFL(JR), -IFL(3), IFLA, LLIST(I+1))

      nporig(I+1)= Ipflag*2 + Nint
      IF(ILEAD(1).eq.1.or.ILEAD(2).eq.1) nporig(I+1)= -1 * nporig(I+1)
      
c     replace all for diff.
      IF(IABS(IFQRK).EQ.1)THEN
         IF(IFQRK.eq.-1.and.ipar(11).lt.0
     &        .and.ipar(11).ge.-3) then
            IF(ABS(LLIST(I+1)).EQ.6) LLIST(I+1) = 27*isign(1,LLIST(I+1))
         endif
c     replace all for leading val.
         IF(ipar(11).le.-2.and.ipar(11).ge.-3) then
            if( ilead(jr).eq.1 ) then
               IF(ABS(LLIST(I+1)).EQ.6)
     &              LLIST(I+1) = 27*isign(1,LLIST(I+1))
            endif
         endif

c     increased vec.meson ratio and replace pi0 with rho0
         IF(IFQRK.eq.-1.and.ipar(11).eq.-7) then
            IF(ABS(LLIST(I+1)).EQ.6) LLIST(I+1) = 27*isign(1,LLIST(I+1))
c     IF(ABS(LLIST(I+1)).EQ.7)  LLIST(I+1) = 25*isign(1,LLIST(I+1))
         endif
         
c     replace all for diff. ( same as lvec6 but for rhop as well )
c     reset vec.meson ratio
         IF(IFQRK.eq.-1.and.ipar(11).eq.-8) then
            PAR(5) = 9.
            if( INT((PAR(5)+1)*S_RNDM(0)).gt.1 ) then
               IF(ABS(LLIST(I+1)).EQ.6)
     &              LLIST(I+1) = 27*isign(1,LLIST(I+1))
               IF(ABS(LLIST(I+1)).EQ.7)
     &              LLIST(I+1) = 25*isign(1,LLIST(I+1))
            endif
         endif
c     replace leading pseudoscalar by vector
         IF(ipar(11).eq.-8.and.ilead(jr).eq.1) then
            PAR(5) = 9.
            if( INT((PAR(5)+1)*S_RNDM(0)).gt.1 ) then
               IF(ABS(LLIST(I+1)).EQ.6) 
     &              LLIST(I+1) = 27*isign(1,LLIST(I+1))
               IF(ABS(LLIST(I+1)).EQ.7)
     &              LLIST(I+1) = 25*isign(1,LLIST(I+1))
            endif
         endif
         
c     replace all pi0 for string diff.( same as lvec7 but for rhop as well )
         IF(IFQRK.eq.-1.and.ipar(11).eq.-9) then
            if(ABS(LLIST(I+1)).EQ.6) LLIST(I+1) = 27*isign(1,LLIST(I+1))
         endif
c     replace leading pi0 by vector
         IF(ipar(11).eq.-9.and.ILEAD(JR).eq.1) then
            if(ABS(LLIST(I+1)).EQ.6) LLIST(I+1) = 27*isign(1,LLIST(I+1))
         endif

c     replace leading omega in string diff by rho0's 
c     suppress pi0 in diff. strings
c     in addition to increased leading vec.meson ratio (lvec22)
         IF(IFQRK.eq.-1.and.IPAR(11).eq.-17)THEN
            IF(ABS(LLIST(I+1)).EQ.6)THEN
c     print*,'found pi0, restarting..'
               GOTO 666
            ENDIF
         ENDIF
         ILEAD(JR)= 0
      ENDIF
c     reset vec.mes. ratio
      PAR(5) = PAR5_def
      PAR(24) = charmPARdef
      IPAR(11) = IPAR11_def

      P(I,1)   = PX(JT)+PX(3)      
      P(I,2)   = PY(JT)+PY(3)
      LPOINT(I) = JT
      I1 = I+1
      nforig(I1) = 0      
      P(I1,5) = AM(IABS(LLIST(I1)))
      P(I1,1) = PX(JR)-PX(3)      
      P(I1,2) = PY(JR)-PY(3)   
      LPOINT(I1) = JR 
      XM1 = P(I,5)**2+P(I,1)**2+P(I,2)**2
      XM2 = P(I1,5)**2+P(I1,1)**2+P(I1,2)**2
      IF (SQRT(XM1)+SQRT(XM2) .GT. SQRT(WREM2)) GOTO 200
      WREM = SQRT(WREM2)

c...RE & EJA fix
      PT2 = (P(I,1)+P(I1,1))**2+(P(I,2)+P(I1,2))**2
      WREMPT = sqrt(WREM2+PT2)
      EA1 = (WREM2+XM1-XM2+PT2)/(2.*WREMPT)
    
      PA2 = (EA1**2-XM1)
      if (pa2.gt.0)  then
            PA = SQRT(PA2)
      else
            goto 200
      endif
      BA = PTOT(3)/PTOT(4)
      GA = PTOT(4)/WREMPT
      SGN = FLOAT(3-2*JT)
      P(I,3) = GA*(BA*EA1+SGN*PA)
      P(I,4) = GA*(EA1+BA*SGN*PA)
      P(I+1,3) = PTOT(3)-P(I,3)
      P(I+1,4) = PTOT(4)-P(I,4)


      NA= NP+1
      NP=I+1
         
C...reorder  particles along chain (in rank)
      IF (LRANK)  THEN
      N1 = NA-1
      N2 = 0
      DO J=NA,NP
         IF(P(J,4).lt.0) THEN
            NP=NA-1
            GOTO 200            ! negative energy bug 'fix'
         ENDIF
         IF(LPOINT(J) .EQ. 2)  THEN
            N2=N2+1
            LLIST (NP+N2) = LLIST(J)
            nporig(NP+N2) = nporig(J)
            nforig(NP+N2) = 0
            DO K=1,5
               P(NP+N2,K)=P(J,K)
            ENDDO
         ELSE
            N1= N1+1
            IF (N1.LT.J)   THEN
               LLIST(N1) = LLIST(J)
               nporig(N1) = nporig(J)
               nforig(N1) = nforig(J)
               DO K=1,5
                  P(N1,K) = P(J,K)
               ENDDO
            ENDIF
         ENDIF
      ENDDO
      JJ=N1
      DO J=NP+N2,NP+1,-1
         JJ= JJ+1
         LLIST(JJ) = LLIST(J)
         nporig(JJ) = nporig(J)
         nforig(JJ) = nforig(J)
         DO K=1,5
            P(JJ,K) = P(J,K)
         ENDDO
      ENDDO
      ENDIF

      if(Ndebug.gt.2)
     &     WRITE(LUN,*)' STRING_FRAG_4FLV: NP after fragmentation:',NP

      END


      FUNCTION ZDIS (IFL1,ifl2, XMT2,lc)
c	includes charmed fragmentation
c	if lc = 0 - no charm, use default fragmentation function
c	   lc = 1 - charm, use Peterson/SLAC (zmefn)
C...z distribution
      SAVE
      COMMON /S_CZDIS/ FAin, FB0in
      COMMON /S_CZDISs/ FAs1, fAs2
      COMMON /S_CZDISc/ ZDMAX, EPSI
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, kb ,kt
      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)

      IF (lc .eq. 1) then
 90      z = s_rndm(0)
         tcp = zmefn(z,epsi)/zdmax
         if (tcp .lt. s_rndm(0)) goto 90
         zdis = z
      ELSE                      ! original
      fa=fain
      fb0=fb0in
c-----hard scattering fragmentation
      IF(IPAR(6).eq.2)THEN
         fa= PAR(18)
         fb0= PAR(19)
      ENDIF   
CDH   correction  may 10-1996
      if (iabs(kb).ge.13) then   ! baryons only
          if (abs(ifl2).eq.3)  fa=fain+fas2
          if (abs(ifl1).eq.3)  fa=fain+fas1
      endif
      FB = FB0*XMT2
      IF(FA.GT.0.01.AND.ABS(FA-1.)/FB.LE.0.01) ZMAX=FB/(1.+FB)+
     +  (1.-FA)*FB**2/(1.+FB)**3
      IF(FA.GT.0.01.AND.ABS(FA-1.)/FB.GT.0.01) ZMAX=0.5*(1.+FB-
     +  SQRT((1.-FB)**2+4.*FA*FB))/(1.-FA)
      IF(ZMAX.LT.0.1)  ZDIV=2.75*ZMAX
      IF(ZMAX.GT.0.85) 
     +     ZDIV=ZMAX-0.6/FB**2+(FA/FB)*ALOG((0.01+FA)/FB)
C...Choice if z, preweighted for peaks at low or high z
100   Z=S_RNDM(0)
      IDIV=1
      FPRE=1.
      IF (ZMAX.LT.0.1)  THEN
         IF(1..LT.S_RNDM(0)*(1.-ALOG(ZDIV)))  IDIV=2
         IF (IDIV.EQ.1)  Z=ZDIV*Z
         IF (IDIV.EQ.2)  Z=ZDIV**Z
         IF (IDIV.EQ.2)  FPRE=ZDIV/Z
      ELSEIF (ZMAX.GT.0.85)  THEN
         IF(1..LT.S_RNDM(0)*(FB*(1.-ZDIV)+1.)) IDIV=2
         IF (IDIV.EQ.1)  Z=ZDIV+ALOG(Z)/FB
         IF (IDIV.EQ.1)  FPRE=EXP(FB*(Z-ZDIV))
         IF (IDIV.EQ.2)  Z=ZDIV+Z*(1.-ZDIV)
      ENDIF
C...weighting according to the correct formula
      IF (Z.LE.FB/(50.+FB).OR.Z.GE.1.)  GOTO 100
      FVAL=(ZMAX/Z)*EXP(FB*(1./ZMAX-1./Z))
      IF(FA.GT.0.01)  FVAL=((1.-Z)/(1.-ZMAX))**FA*FVAL
      IF(FVAL.LT.S_RNDM(0)*FPRE)  GOTO 100
      ZDIS=Z


	endif

	return
	end	

      FUNCTION ZDIS_4FLV (IFL1,IFL2, XMT2)
C...z distribution
c     includes charmed fragmentation (Peterson/SLAC)
      SAVE
      COMMON /S_CZDIS/ FAin, FB0in
      COMMON /S_CZDISs/ FAs1, fAs2
      COMMON /S_CZDISc/ ZDMAX, EPSI
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, kb ,kt
      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)

      IAFL1 = IABS(mod(IFL1,100))
      IAFL2 = IABS(mod(IFL2,100))
c     SLAC-Peterson fragmentation function for charm
      IF ((IAFL1/10.eq.4.or.mod(IAFL1,10).eq.4)
     +     .or.(IAFL2/10.eq.4.or.mod(IAFL2,10).eq.4))THEN
 90      z = s_rndm(0)
         tcp = zmefn(z,epsi)/zdmax
         if (tcp .lt. s_rndm(0)) goto 90
         zdis_4flv = z
      else
c     original lund function, non charm
         fa=fain                ! lund parameter a
         fb0=fb0in              ! lund parameter b
c-----hard scattering fragmentation
         IF(IPAR(6).eq.2)THEN
            fa= PAR(18)
            fb0= PAR(19)
         ENDIF   
C     DH   correction  may 10-1996
         if (iabs(kb).ge.13) then ! baryons only
            if (iafl2.eq.3)  fa=fain+fas2
            if (iafl1.eq.3)  fa=fain+fas1
         endif
c     special parameters for baryon fragmentation
c     similar to pythia
         IF((IAFL1+IAFL2).gt.10.and.
     &        (IPAR(36).eq.1.or.IPAR(20).eq.3))then
            fa = fain + PAR(45)
            fb0 = PAR(60)
         ENDIF        
         FB = FB0*XMT2
         IF(FA.GT.0.01.AND.ABS(FA-1.)/FB.LE.0.01) ZMAX=FB/(1.+FB)+
     +        (1.-FA)*FB**2/(1.+FB)**3
         IF(FA.GT.0.01.AND.ABS(FA-1.)/FB.GT.0.01) ZMAX=0.5*(1.+FB-
     +        SQRT((1.-FB)**2+4.*FA*FB))/(1.-FA)
         IF(ZMAX.LT.0.1)  ZDIV=2.75*ZMAX
         IF(ZMAX.GT.0.85) 
     +        ZDIV=ZMAX-0.6/FB**2+(FA/FB)*ALOG((0.01+FA)/FB)
C...  Choice if z, preweighted for peaks at low or high z
 100     Z=S_RNDM(0)
         IDIV=1
         FPRE=1.
         IF (ZMAX.LT.0.1)  THEN
            IF(1..LT.S_RNDM(0)*(1.-ALOG(ZDIV)))  IDIV=2
            IF (IDIV.EQ.1)  Z=ZDIV*Z
            IF (IDIV.EQ.2)  Z=ZDIV**Z
            IF (IDIV.EQ.2)  FPRE=ZDIV/Z
         ELSEIF (ZMAX.GT.0.85)  THEN
            IF(1..LT.S_RNDM(0)*(FB*(1.-ZDIV)+1.)) IDIV=2
            IF (IDIV.EQ.1)  Z=ZDIV+ALOG(Z)/FB
            IF (IDIV.EQ.1)  FPRE=EXP(FB*(Z-ZDIV))
            IF (IDIV.EQ.2)  Z=ZDIV+Z*(1.-ZDIV)
         ENDIF
C...weighting according to the correct formula
         IF (Z.LE.FB/(50.+FB).OR.Z.GE.1.)  GOTO 100
         FVAL=(ZMAX/Z)*EXP(FB*(1./ZMAX-1./Z))
         IF(FA.GT.0.01)  FVAL=((1.-Z)/(1.-ZMAX))**FA*FVAL
         IF(FVAL.LT.S_RNDM(0)*FPRE)  GOTO 100
         ZDIS_4FLV=Z
         
      ENDIF
      
      RETURN
      END	
      
      SUBROUTINE ZNORMAL
C...  normalisation for Peterson/SLAC frag. func
      COMMON /S_CZDISc/ ZDMAX, EPSI
      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
c     get the maximum zmefn value first for normalisation
      jmax = 1000
      zdmax = 1.e-10
      DO j = 1, jmax
         z = float(j)/float(jmax+1)
         zdmax = max(zdmax, zmefn(z,epsi))
      enddo
      WRITE(LUN,*)'ZDMAX,EPS:',zdmax, epsi
      RETURN
      END	
	
      FUNCTION ZMEFN(z,eps)
C...  Peterson/SLAC frag. func
      zmefn = (z*(1.-z**(-1)-eps/(1.-z))**2)**(-1)
      RETURN
      END
	

      FUNCTION ZBLEAD (LB)
C...fragmentation function for leading baryon
C.  simple form:  f(z) = a + x**b
C   INPUT : LB = particle code.
C..................................................
      SAVE
      COMMON /S_CZLEAD/ CLEAD, FLEAD
c      COMMON /S_SZLEAD/ CLEADs, FLEADs
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)
      COMMON /S_CHP/ ICHP(99), ISTR(99), IBAR(99)

      IF(IPAR(30).eq.1)THEN
C     Sibyll 2.1 hard fragmentation function

      IC = ICHP(Lb)*ISIGN(1,Lb)

      if (lb.ge.34.and.lb.le.39)  then  ! Lambda's and Sigma's
  665               ZBLEAD = S_RNDM(0)
                if (zblead.le..01) goto 665
c          zblead=zdisn(1) ! blead**2   ! soft
      else if (ic.eq.0)     then
          zblead=zdisn(1)   ! blead**2   !soft
      else if (ic.eq.1)  then  ! fast protons only
            if (abs(lb).eq.13) then
              IF (S_RNDM(0) .LT. CLEAD)  THEN
  666               ZBLEAD = S_RNDM(0)
                if (zblead.le..01) goto 666
              ELSE
                  zblead=1.-zdisn(1)  ! zblead**2   !hard
              ENDIF
            continue
           else
               zblead=zdisn(1)  ! zblead**2   !hard
           endif   
      else if (ic.eq.2)  then  ! fast delta++
          zblead=1.- zdisn(1)  ! (zblead)**.3333
      else
               zblead=S_RNDM(0) ! zdisn(1)     !hard
      endif

      RETURN
C     flat baryon fragmentation function
      ELSE
 999     zblead = s_rndm(0)
         if (zblead .le. 0.01) goto 999     
         RETURN
      ENDIF
      END


      FUNCTION ZDISN (n)
C...Generate (1-x)**n
      SAVE
666   rmin=1.1
      do i=1,n+1
         R1=S_RNDM(0)
         IF (R1.LE.RMIN) RMIN=R1
      ENDDO
      ZDISn=RMIN
      if (zdisn.le..01) goto 666
      if (zdisn.ge..99) goto 666
      END


      SUBROUTINE GG_FRAG_4FLV (E0)
C...This routine fragments a  gluon-gluon system
C.  of mass E0 (GeV)
C.  the particles produced are in the  jet-jet frame
C.  oriented along the z axis
C...........................................................
      SAVE
      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      COMMON /S_PLIST/ P(8000,5), LLIST(8000), NP
      COMMON /S_MASS1/ AM(99), AM2(99)
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)
      DIMENSION WW(2,2),PTOT(4),PX(3),PY(3),IFL(3),PMQ(3)
      COMMON /S_PARTO/ NFORIG(8000), NPORIG(8000), IPFLAG, NINT
      if(Ndebug.gt.2) then
         WRITE(LUN,*)' GG_FRAG_4FLV: called with (E0)',
     &        E0
         WRITE(LUN,*)' GG_FRAG_4FLV: NP before fragmentation:',NP
      endif

C...  'leading' strange fraction
      PAR2_def = PAR(2)
      IF(IPAR(39).eq.2) PAR(2) = PAR(66)

C...Generate the 'forward' leading particle.
100   I = NP+1
      I0 = -1 + 2.*INT(1.9999*S_RNDM(0))
      CALL SIB_I4FLAV(I0,0,IFL1, LDUM)
      CALL SIB_I4FLAV(IFL1,0,IFL2, LLIST(I))
      CALL PTDIS_4FLV(IFL1,PX1,PY1)
      CALL PTDIS_4FLV(IFL2,PX2,PY2)
      P(I,1) = PX1+PX2
      P(I,2) = PY1+PY2
      P(I,5) = AM(IABS(LLIST(I)))
      XM1 = P(I,5)**2+P(I,1)**2+P(I,2)**2
      Z1 = ZDIS_4FLV (IFL1,1,0.25*XM1)
      Z2 = ZDIS_4FLV (IFL2,1,0.25*XM1)
      T1  = 4.*XM1/(E0*E0*(Z1+Z2))
      P(I,4) = 0.25*E0*(Z1+Z2 + T1)
      P(I,3) = 0.25*E0*(Z1+Z2 - T1)

      nforig(I)= 0
      nporig(I)= Ipflag*3 + Nint


C...Generate the 'backward' leading particle.
      I = I+1
      CALL SIB_I4FLAV(-I0,0,IFL3, LDUM)
      CALL SIB_I4FLAV(IFL3,0,IFL4, LLIST(I))
      CALL PTDIS_4FLV(IFL3,PX3,PY3)
      CALL PTDIS_4FLV(IFL4,PX4,PY4)
      P(I,1) = PX3+PX4
      P(I,2) = PY3+PY4
      P(I,5) = AM(IABS(LLIST(I)))
      XM2 = P(I,5)**2+P(I,1)**2+P(I,2)**2
      Z3 = ZDIS_4FLV (IFL3,1,0.25*XM2)
      Z4 = ZDIS_4FLV (IFL4,1,0.25*XM2)
      T2  = 4.*XM2/(E0*E0*(Z3+Z4))
      P(I,4) = 0.25*E0*( Z3+Z4 + T2)
      P(I,3) = 0.25*E0*(-Z3-Z4 + T2)

      nforig(I)= 0
      nporig(I)= Ipflag*3 + Nint

c     reset strange fraction
      PAR(2) = PAR2_def
      
C...Fragment the two remaning strings
      N0 = 0
      DO KS=1,2
      
      NTRY = 0
200      NTRY = NTRY+1
      I = NP+2+N0
      IF (NTRY .GT. 30)  GOTO 100

      IF (KS .EQ. 1)  THEN
         WW(1,1) = 0.5 * (1 - Z1 - 0.5*T2) 
         WW(2,1) = 0.5 * (1 - Z3 - 0.5*T1)
         PX(1) = -PX1
         PY(1) = -PY1
         PX(2) = -PX3
         PY(2) = -PY3
         IFL(1) = -IFL1
         IFL(2) = -IFL3
      ELSE
         WW(1,1) = 0.5 * (1 - Z2 - 0.5*T2) 
         WW(2,1) = 0.5 * (1 - Z4 - 0.5*T1)
         PX(1) = -PX2
         PY(1) = -PY2
         PX(2) = -PX4
         PY(2) = -PY4
         IFL(1) = -IFL2
         IFL(2) = -IFL4
      ENDIF
      PX(3) = 0.
      PY(3) = 0.
      PTOT (1) = PX(1)+PX(2)
      PTOT (2) = PY(1)+PY(2)
      PTOT (3) = 0.5*E0*(WW(1,1)-WW(2,1))
      PTOT (4) = 0.5*E0*(WW(1,1)+WW(2,1))

      PMQ(1) = QMASS(IFL(1))
      PMQ(2) = QMASS(IFL(2))

C...produce new particle: side, pT
300      I=I+1
      if(i.gt.8000) then
        write(6,'(1x,a,i8)') 
     &    'GG_FRAG_4FLV: no space left in S_PLIST:',I
        stop   
      endif
      nforig(I)= 0
      nporig(I)= Ipflag*2 + Nint

      JT=1.5+S_RNDM(0)
      JR=3-JT
c      CALL PTDIS (IFL(JT), PX(3),PY(3))

C...particle ID
      CALL SIB_I4FLAV (IFL(JT), 0, IFL(3), LLIST(I))
      PMQ(3) = QMASS(IFL(3))
      P(I,5) = AM(IABS(LLIST(I)))

      CALL PTDIS_4FLV (IFL(3), PX(3),PY(3))
      
C...test end of fragmentation
      WREM2 = PTOT(4)**2-PTOT(1)**2-PTOT(2)**2-PTOT(3)**2
      IF (WREM2 .LT. 0.1)  GOTO 200
      WMIN = PMQ(1)+PMQ(2)+2.*PMQ(3)+1.1 + (2.*S_RNDM(0)-1.)*0.2
      IF (WREM2 .LT. WMIN**2)THEN
         GOTO 400
      ENDIF

C...fill transverse momentum
      P(I,1) = PX(JT) + PX(3)
      P(I,2) = PY(JT) + PY(3)

C...Choose z
      XMT2 = P(I,5)**2+P(I,1)**2+P(I,2)**2
      Z = ZDIS_4FLV (ifl(3),IFL(JT), XMT2)
      
      WW(JT,2) = Z*WW(JT,1)
      WW(JR,2) = XMT2/(WW(JT,2)*E0**2)

      P(I,3) = WW(1,2)*0.5*E0 - WW(2,2)*0.5*E0
      P(I,4) = WW(1,2)*0.5*E0 + WW(2,2)*0.5*E0

      DO J=1,4
         PTOT (J) = PTOT(J) - P(I,J)
      ENDDO
      DO K=1,2
         WW(K,1) = WW(K,1) - WW(K,2)
      ENDDO

C...Reset pT and flavor at ends of the string
      PX(JT) = -PX(3)
      PY(JT) = -PY(3)
      IFL(JT) =-IFL(3)
      PMQ(JT) = PMQ(3)
      GOTO 300

C...Final two hadrons
 400  IAFL1 = mod(IABS(IFL(JR)),100)
      IAFL2 = mod(IABS(IFL(3)),100)
      IF (IAFL1*IAFL2 .GT. 100)  GOTO 200 ! reject two diquarks
      IF ((IAFL1/10.eq.4.or.mod(IAFL1,10).eq.4)
     +     .and.(IAFL2/10.eq.4.or.mod(IAFL2,10).eq.4))
     +     GOTO 200             ! reject two charm quarks
      
      CALL SIB_I4FLAV (IFL(JR), -IFL(3), IFLA, LLIST(I+1))
      P(I+1,5) = AM(IABS(LLIST(I+1)))
      P(I,1)   = PX(JT)+PX(3)      
      P(I,2)   = PY(JT)+PY(3)      
      nporig(I)= Ipflag*2 + Nint
      I1 = I+1
      nporig(I1)= Ipflag*2 + Nint
      P(I1,1) = PX(JR)-PX(3)      
      P(I1,2) = PY(JR)-PY(3)      
      XM1 = P(I,5)**2+P(I,1)**2+P(I,2)**2
      XM2 = P(I1,5)**2+P(I1,1)**2+P(I1,2)**2
      IF (SQRT(XM1)+SQRT(XM2) .GT. SQRT(WREM2)) GOTO 200
      if (ptot(4).le.0) goto 200
      PT2 = (P(I,1)+P(I1,1))**2+(P(I,2)+P(I1,2))**2
      WREMPT = sqrt(WREM2+PT2)
      EA1 = (WREM2+XM1-XM2+PT2)/(2.*WREMPT)
      PA2 = (EA1**2-XM1)
      if (PA2.ge.0.0) then
        PA = SQRT(PA2)
      else
         goto 200
      endif
      BA = PTOT(3)/PTOT(4)
      GA = PTOT(4)/WREMPT
      SGN = FLOAT(3-2*JT)
      P(I,3) = GA*(BA*EA1+SGN*PA)
      P(I,4) = GA*(EA1+BA*SGN*PA)
      P(I+1,3) = PTOT(3)-P(I,3)
      P(I+1,4) = PTOT(4)-P(I,4)
      N0 = I-NP-1
      ENDDO                  ! loop on two `remaining strings'
      NP = I+1
      IF(Ndebug.gt.2) then
         WRITE(LUN,*)' GG_FRAG_4FLV: NP after fragmentation:',NP
      ENDIF
      RETURN
      END

      FUNCTION QMASS(IFL)
C-----------------------------------------------------------------------
C...Return quark or diquark constituent masses
C     extended to 4 flavors /FR'13
C-----------------------------------------------------------------------
      SAVE
      DIMENSION QMAS(4)
      DATA QMAS /0.325,0.325,0.5,1.5/
      IFLA = IABS(IFL)
      IFLA = MOD(IFLA,100)
      IF (IFLA .LE. 4)       THEN
         QMASS = QMAS(IFLA)
      ELSE
         QMA = QMAS(IFLA/10)
         QMB = QMAS(MOD(IFLA,10))
         QMASS = QMA+QMB
      ENDIF
      RETURN
      END


      SUBROUTINE SIB_IFLAV (IFL1,IFL2_A, IFL2, KF,iqch)
C-----------------------------------------------------------------------
C.  This subroutine receives as input IFL1 the flavor code
C.  of a quark (antiquark) and  generates the antiquark (quark)
C.  of flavor code IFL2 that combine with the original parton
C.  to compose an hadron of code KF. ONLY 3 FLAVORS
C.  If (IFL2_A.NE.0) returns an hadron KF composed of IFL1 and IFL2_A
c	2010.05.18 reinstating suppression of s->c for valance quark
c	iqch = 0 is valence quark: no s->c
c	iqch = 1 s->c indiscriminately
c
c     Input: IFL1 - flavor of first quark
c            IFL2_A - flavor of second quark ( if 0 randomly chosen ) 
c            IQCH - charm flag
c     Output: IFL2 - flavor of second quark partner to be passed on
c             KF - final hadron
C-----------------------------------------------------------------------
      SAVE

      COMMON /S_CFLAFR/ PAR(80), IPAR(80)

      DIMENSION KFLA(3,3,2), CDIAG(12), KDIAG(6)
      DIMENSION KBAR(30), CFR(12), KFR(80), KCHMES(33)
      DATA KFLA /0,8,10,7,0,22,9,21,0,0,26,29,25,0,31,28,30,0/
      DATA CDIAG /0.5,0.25,0.5,0.25,1.,0.5,0.5,0.,0.5,0.,1.,1./
      DATA KDIAG /6,23,24,27,32,33/
      DATA KBAR /13,14,34,35,36,37,38,9*0,39,3*0,40,41,42,43,44,
     +             45,46,47,48,49/
      DATA CFR /0.75,0.,0.5,0.,0.,1.,0.1667,0.3333,0.0833,0.6667,
     +            0.1667,0.3333/
      DATA KFR/0,16,17,19,100,104,109,115,0,26,27,29,122,126,131,137
     +  ,0,40,42,47,144,158,178,205,0,1,3,6,10,15,21,28,0,0,56,57,240,
     +  246,256,271,0,0,1,3,6,10,15,21,60,61,64,70,292,307,328,356,
     +  0,1,3,6,10,15,21,28,16*0/
      DATA KCHMES /8*0,72,71,10*0,60,59,73,4*0,81,80,79,78,0,83/

      IFLA = IABS(IFL1)
      IFL2A = IFL2_A
      IF (IFL2A .NE. 0)  THEN
         IFL2A = MOD(IFL2A,100)
         IFL2 = IFL2A
         IFLB = IABS(IFL2A)
         MB = 0
         IF (IFLB .GT. 10)   MB=1
         IF (IFLA .GT. 10)   MB=2
      ELSE
          MB = 2
         IF (IFLA .LT. 10)   THEN
             MB = 1
             IF ((1.+PAR(1))*S_RNDM(0).LT. 1.)  MB=0
         ENDIF
      ENDIF
      
      IF (MB .EQ. 0)  THEN
         IF (IFL2A.EQ.0)  
     +        IFL2=ISIGN(1+INT((2.+PAR(2))*S_RNDM(0)),-IFL1)
         IFLD = MAX(IFL1,IFL2)
         IFLE = MIN(IFL1,IFL2)
         GOTO 100
      ENDIF

C...Decide if the diquark must be split
      IF (MB .EQ. 2 .AND. IFLA .GT. 100)   THEN
         IFLA = MOD(IFLA,100)
           GOTO 200
      ENDIF
      IF (MB .EQ. 2 .AND. IFLA .EQ. 0)   THEN
          IF (S_RNDM(0) .LT. PAR(8))  THEN
             MB = 0
             IFLG = MOD(IFL1,10)
             IFLH =(IFL1-IFLG)/10
             IF (S_RNDM(0) .GT. 0.5)  THEN
                IFLDUM = IFLG
                IFLG = IFLH
                IFLH = IFLDUM
             ENDIF
             IFL11=IFLG
             IFL22=ISIGN(1+INT((2.+PAR(2))*S_RNDM(0)),-IFL1)
             IFLD = MAX(IFL11,IFL22)
             IFLE = MIN(IFL11,IFL22)
             IFL2 = -IFLH*10+IFL22
             IF (S_RNDM(0) .GT. 0.5)  IFL2 = IFL22*10-IFLH
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
         IF (IF1 .NE. IF2)   THEN
            KF = KFLA(IF1,IF2,KSP+1)
         ELSE
            R = S_RNDM(0)
            JF=1+INT(R+CDIAG(6*KSP+2*IF1-1))+
     +             INT(R+CDIAG(6*KSP+2*IF1))
            JF = MIN(JF,3)
            KF=KDIAG(JF+3*KSP)
         ENDIF

         IF(IPAR(15).ne.9)THEN
c. no c correlation, simple branching ratio
            if (iqch .gt. 0.and.KF.ne.24) then
               if (if1 .eq. 3 .or. if2 .eq. 3) then
                  if (s_rndm(0) .lt. PAR(24)) kf = kf + 50
               endif
            endif

         ELSE
c. c conservation for mesons
            IF(IQCH.eq.0.or.KF.eq.24) RETURN
            IF(mod(if1*if2,3).ne.0) RETURN
            if (iqch .gt. 0) then
               IF(ifla.eq.3) RETURN
               if (s_rndm(0) .lt. PAR(24)) then
                  kf = KCHMES( KF )
                  IF(IF1.ne.IF2) iqch = -1 * iqch
               ENDIF
            else
               kf = KCHMES( KF )
               IF(IF1.ne.IF2) iqch = -1 * iqch
            endif
         ENDIF

         RETURN
      ENDIF

C...Form a baryon
200      IF (IFL2A .NE. 0)   THEN
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
110          CONTINUE
          IF(MB.EQ.1)   THEN            ! generate diquark
             IFLD = IFLA
120             IFLE = 1+INT((2.+PAR(2)*PAR(3))*S_RNDM(0))          
             IFLF = 1+INT((2.+PAR(2)*PAR(3))*S_RNDM(0))          
             IF(IFLE.GE.IFLF.AND.PAR(4).LT.S_RNDM(0))    GOTO 120
             IF(IFLE.LT.IFLF.AND.PAR(4)*S_RNDM(0).GT.1.) GOTO 120     
             IFL2=ISIGN(10*IFLE+IFLF,IFL1)
          ELSE                  ! generate quark
             IFL2=ISIGN(1+INT((2.+PAR(2))*S_RNDM(0)),IFL1)
             IFLD=IABS(IFL2)
             IFLE=IFLA/10
             IFLF=MOD(IFLA,10)
          ENDIF
C...SU(6) factors for baryon formation
             LFR=3+2*((2*(IFLE-IFLF))/(1+IABS(IFLE-IFLF)))
          IF(IFLD.NE.IFLE.AND.IFLD.NE.IFLF)  LFR=LFR+1
          WT = CFR(2*LFR-1)+PAR(7)*CFR(2*LFR)
          IF(IFLE.LT.IFLF)   WT=WT/3.
          IF (WT.LT.S_RNDM(0)) GOTO 110
      ENDIF

C...Form Baryon
      IFLG=MAX(IFLD,IFLE,IFLF)
      IFLI=MIN(IFLD,IFLE,IFLF)
      IFLH=IFLD+IFLE+IFLF-IFLG-IFLI
      KSP=2+2*INT(1.-CFR(2*LFR-1)+(CFR(2*LFR-1)+PAR(7)*
     1       CFR(2*LFR))*S_RNDM(0))

C...Distinguish Lambda- and Sigma- like particles
      IF (KSP.EQ.2.AND.IFLG.GT.IFLH.AND.IFLH.GT.IFLI)  THEN
      IF(IFLE.GT.IFLF.AND.IFLD.NE.IFLG) KSP=2+INT(0.75+S_RNDM(0))
       IF(IFLE.LT.IFLF.AND.IFLD.EQ.IFLG) KSP=3
       IF(IFLE.LT.IFLF.AND.IFLD.NE.IFLG) KSP=2+INT(0.25+S_RNDM(0))
      ENDIF
      KF=KFR(16*KSP-16+IFLG)+KFR(16*KSP-8+IFLH)+IFLI
      KF=ISIGN(KBAR(KF-40),IFL1)

      IF(IPAR(15).ne.9)THEN
c... c/s ratio for baryon 
         if (iqch .gt. 0) then
            if (iabs(ifld) .eq. 3 .or. iabs(ifle) .eq. 3 .or. 
     +           iabs(iflf) .eq. 3) then   
               IF(IPAR(15).eq.5) THEN
                  if (s_rndm(0) .lt. PAR(29)) kf = abs(kf) + 50
               ELSE
                  if (s_rndm(0) .lt. PAR(24)) kf = abs(kf) + 50
               ENDIF
               kf = isign(kf, ifl1)
            endif 	 
         endif
         kf = isign(kf, ifl1)
      ELSE
c... c conservation for baryons
         IF(IQCH.eq.0) RETURN
         IF(mod(ifld*ifle*iflf,3).ne.0) RETURN
         if (iqch .gt. 0) then
            IF((ifla/10).eq.3.or.mod(ifla,10).eq.3) RETURN
            if (s_rndm(0) .lt. PAR(24))then 
               kf = abs(kf) + 50
               iqch = -1 * iqch
            ENDIF
         else
            kf = abs(kf) + 50
            iqch = -1 * iqch
         endif
         kf = isign(kf, ifl1)
      ENDIF

      RETURN
      END

      SUBROUTINE SIB_I4FLAV (IFL1,IFL2_A, IFL2, KF)
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
c            IQCH - charm flag
c     Output: IFL2 - flavor of second quark partner to be passed on
c             KF - final hadron
C-----------------------------------------------------------------------
      SAVE

      COMMON /S_CFLAFR/ PAR(80), IPAR(80)
      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      DIMENSION KFLA(4,4,2), CDIAG(16), KDIAG(8)
      DIMENSION KBAR(40), CFR(28), KFR(80)
      DATA KFLA /0,8,10,71,7,0,22,59,9,21,0,74,72,60,75,0, ! spons
     +     0,26,29,80,25,0,31,78,28,30,0,76,81,79,77,0/ ! spin-
      DATA CDIAG /0.5,0.25,0.5,0.25,1.,0.5,2.,1., ! spin-zero dons
     +     0.5,0.,0.5,0.,1.,1.,2.,1./ ! spin-one diagonal meson
      DATA KDIAG /6,23,24,73,27,32,33,83/
      DATA KBAR /13,14,34,35,36,37,38,84,85,86, ! jetset->sibyl
     +     87,88,99,3*0,39,89,87,88, 
     +     40,41,42,43,44,45,46,47,48,49,      
     +     94,95,96,97,98,99,4*0/ ! spin-3/2 css baryon 
      DATA CFR /0.75,0.,0.5,0.,0.,1.,0.1667,0.3333,0.0833,0.6667,0.1667,
     &     0.3333,-3.,1.,-2.,-2.,1.,0.,0.,-3.,1.,1.,1.,5*0./
      DATA KFR/0,16,17,19,100,104,109,115,0,26,27,29,122,126,131,137
     +  ,0,40,42,47,144,158,178,205,0,1,3,6,10,15,21,28,0,0,56,57,240,
     +  246,256,271,0,0,1,3,6,10,15,21,60,61,64,70,292,307,328,356,
     +  0,1,3,6,10,15,21,28,16*0/
     
c      PRINT*,' I4FLAV called with :(IFL1,IFL2_A)',IFL1,IFL2_A

      IFLA = IABS(IFL1)
      IFL2A = IFL2_A
      IF (IFL2A .NE. 0)  THEN
         IFL2A = MOD(IFL2A,100)
         IFL2 = IFL2A
         IFLB = IABS(IFL2A)
         MB = 0
         IF (IFLB .GT. 10)   MB=1
         IF (IFLA .GT. 10)   MB=2
      ELSE
          MB = 2
         IF (IFLA .LT. 10)   THEN
             MB = 1
             IF ((1.+PAR(1))*S_RNDM(0).LT. 1.)  MB=0
         ENDIF
      ENDIF
      
      IF (MB .EQ. 0)  THEN
         IF (IFL2A.EQ.0)  
     +        IFL2=ISIGN(1+INT((2.+PAR(2))*S_RNDM(0)),-IFL1)
         IF(IFL2A.EQ.0.and.IABS(IFL2).eq.3) THEN
            IF(S_RNDM(0).lt.PAR(24))
     +           IFL2=ISIGN(4,-IFL1)
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
             IF (S_RNDM(0) .GT. 0.5)  THEN
                IFLDUM = IFLG
                IFLG = IFLH
                IFLH = IFLDUM
             ENDIF
             IFL11=IFLG
             IFL22=ISIGN(1+INT((2.+PAR(2))*S_RNDM(0)),-IFL1)
             IFLD = MAX(IFL11,IFL22)
             IFLE = MIN(IFL11,IFL22)
             IFL2 = -IFLH*10+IFL22
             IF (S_RNDM(0) .GT. 0.5)  IFL2 = IFL22*10-IFLH
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
 120        IFLE = 1+INT((2.+PAR(2)*PAR(3))*S_RNDM(0))
            IFLF = 1+INT((2.+PAR(2)*PAR(3))*S_RNDM(0))          
            IF(IFLD.NE.4)THEN
               IF(IFLE.EQ.3)THEN 
                  IF(S_RNDM(0).lt.PAR(24))
     +                 IFLE=4
               ENDIF
               IF(IFLF.EQ.3.and.IFLE.NE.4)THEN 
                  IF(S_RNDM(0).lt.PAR(24))
     +                 IFLF=4
               ENDIF
            ENDIF
            IF(IFLE.GE.IFLF.AND.PAR(4).LT.S_RNDM(0))    GOTO 120
            IF(IFLE.LT.IFLF.AND.PAR(4)*S_RNDM(0).GT.1.) GOTO 120     
            IFL2=ISIGN(10*IFLE+IFLF,IFL1)
         ELSE                   ! generate quark
            IFL2=ISIGN(1+INT((2.+PAR(2))*S_RNDM(0)),IFL1)
            IFLE=IFLA/10
            IFLF=MOD(IFLA,10)
            IF(IABS(IFL2).EQ.3.and.IFLF.ne.4.and.IFLE.ne.4) THEN
               IF(S_RNDM(0).lt.PAR(24))
     +              IFL2=ISIGN(4,IFL1)
            ENDIF
            IFLD=IABS(IFL2)
         ENDIF
C...SU(6) factors for baryon formation
         LFR=3+2*((2*(IFLE-IFLF))/(1+IABS(IFLE-IFLF)))
         IF(IFLD.NE.IFLE.AND.IFLD.NE.IFLF)  LFR=LFR+1
         WT = CFR(2*LFR-1)+PAR(7)*CFR(2*LFR)
         IF(IFLE.LT.IFLF)   WT=WT/3.
         IF (WT.LT.S_RNDM(0)) GOTO 110
      ENDIF

C...Form Baryon
      IFLG=MAX(IFLD,IFLE,IFLF)
      IFLI=MIN(IFLD,IFLE,IFLF)
      IFLH=IFLD+IFLE+IFLF-IFLG-IFLI
      KSP=2+2*INT(1.-CFR(2*LFR-1)+(CFR(2*LFR-1)+PAR(7)*
     1       CFR(2*LFR))*S_RNDM(0))

C...Distinguish Lambda- and Sigma- like particles
      IF (KSP.EQ.2.AND.IFLG.GT.IFLH.AND.IFLH.GT.IFLI)  THEN
         IF(IFLE.GT.IFLF.AND.IFLD.NE.IFLG) KSP=2+INT(0.75+S_RNDM(0))
         IF(IFLE.LT.IFLF.AND.IFLD.EQ.IFLG) KSP=3
         IF(IFLE.LT.IFLF.AND.IFLD.NE.IFLG) KSP=2+INT(0.25+S_RNDM(0))
      ENDIF
      KF=KFR(16*KSP-16+IFLG)+KFR(16*KSP-8+IFLH)+IFLI
      IF(KF-40.le.0) 
     &     WRITE(LUN,*)'missing particle code?',
     &     KF-40,KF,IFLG,IFLH,IFLI,mod(23,100)
      IF(KBAR(KF-40).eq.0)
     +     WRITE(LUN,*)'jetset code missing,flvs:',kf,IFLG,IFLH,IFLI
      KF=ISIGN(KBAR(KF-40),IFL1)
      
      RETURN
      END


      SUBROUTINE PTDIS (IFL,PX,PY,LPTL)
C...Generate pT
      SAVE
      COMMON /S_CQDIS/ PPT0(35),ptflag
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)
      PPTT = PPT0(IABS(IFL))
      LL = abs( LPTL )
c.  force larger value for charm ptls
c   charm mesons
      IF (LL .GT. 49 .AND. LL .LE. 83) PPTT = PPT0(34)
c   charm baryons
      IF (LL .GT. 83) PPTT = PPT0(35)
c     Gaussian distribution
      PT = PPTT*SQRT(-ALOG(MAX(1E-10,S_RNDM(IFL))))
      IF (IPAR(3).GE.1) THEN
         IF(IABS(IFL).NE.10) THEN
            XM = QMASS(IFL)
            IF(IPAR(3).ge.6.and.LL.gt.49)THEN
               XM = 1.5         ! charmed meson mass
               IF(LL.gt.83) XM = 2. ! charmed baryon mass
            ENDIF
         ELSE
            XM = .5             ! pomeron mass
            IF(IPAR(3).ge.6) XM = 0.
         ENDIF
c     exponential transverse mass
         XM2 = XM**2
         RNDM = MAX(1E-10,S_RNDM(IFL))
         XMT = PPTT * ALOG(RNDM) - XM
         XMT2 = XMT**2
         PT = SQRT(XMT2-XM2)
      ENDIF      
      PHI= 6.2831853*S_RNDM(IFL)
      PX=PT*COS(PHI)
      PY=PT*SIN(PHI)
      RETURN
      END


      SUBROUTINE PTDIS_4FLV (IFL,PX,PY)
C...Generate pT
      SAVE
      COMMON /S_CQDIS2/ PPT02(44)
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)
      IFLA = IABS(IFL)
      IFLA = MOD(IFLA,100)
      PPTT = PPT02(IFLA)
c     Gaussian distribution
      PT = PPTT*SQRT(-ALOG(MAX(1E-10,S_RNDM(IFL))))
      IF (IPAR(3).GE.1) THEN
         IF(IFLA.NE.10) THEN
            XM = QMASS(IFLA)
         ELSE
            XM = .5             ! pomeron mass
            IF(IPAR(3).ge.6) XM = 0.
         ENDIF
c     exponential transverse mass
         XM2 = XM**2
         RNDM = MAX(1E-10,S_RNDM(IFL))
         XMT = PPTT * ALOG(RNDM) - XM
         XMT2 = XMT**2
         PT = SQRT(XMT2-XM2)
      ENDIF      
      PHI= 6.2831853*S_RNDM(IFL)
      PX=PT*COS(PHI)
      PY=PT*SIN(PHI)
      RETURN
      END


      SUBROUTINE PTSETUP(ECM)
C Setting ptflag = 0 will result in
C     underestimating the P_t at high energies.
C     moved from sib_ndiff to seperate subroutine 
c     so that changes will affect diff. /FR'13

      COMMON /S_DEBUG/ Ncall, Ndebug, Lun      
      COMMON /S_CQDIS/ PPT0(35),ptflag
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)

      SQS = ECM
c     original mean of gaussian distributions
c      ptu=.3+.08*log10(sqs/30.)
c      pts=.45+.08*log10(sqs/30.)
c      ptqq=.6+.08*log10(sqs/30.)
     
c     alternate set of parameters used for 
c     exp. transverse mass distribution /FR
c      PTU=.15+.007*log10(sqs/20.)**2
c      PTS=.3+.007*log10(sqs/20.)**2
c      PTQQ=.3+.03*log10(sqs/20.)**2
c     NA22 piC retune
      PTU=.3+.08*log10(sqs/30.)
      PTS=.45+.08*log10(sqs/30.)
      PTQQ=.6+.08*log10(sqs/30.)
      PTPOM= .6+.08*log10(sqs/30.)
      if ( IPAR(3).eq.1 ) then 
c     pt0
         ptu=.15+.007*log10(sqs/20.)**2
         pts=.3+.007*log10(sqs/20.)**2
         ptqq=.3+.03*log10(sqs/20.)**2
         ptpom= .6+.08*log10(sqs/30.)
      elseif ( ipar(3).eq.2 ) then
C     pt1
         ptu=.15+.007*log10(sqs/20.)**2
         pts=.32+.007*log10(sqs/20.)**2
         ptqq=.4+.007*log10(sqs/20.)**2
         ptpom= .6+.08*log10(sqs/30.)
c     pt2
      elseif ( ipar(3).eq.3 ) then
         ptu=.17+.007*log10(sqs/20.)**2
         pts=.3+.007*log10(sqs/20.)**2
         ptqq=.3+.03*log10(sqs/20.)**2
         ptpom = .6+.08*log10(sqs/30.)
      elseif ( ipar(3).eq.5 ) then
         PTU=.16+.007*log10(sqs/20.)**2
         PTS=.28+.007*log10(sqs/20.)**2
         PTQQ= .3+.03*log10(sqs/20.)**2
         PTPOM = .23+.03*log10(sqs/20.)**2
      elseif ( IPAR(3).eq.6 ) then
         PTU=.16+.007*log10(sqs/20.)**2
         PTS=.28+.007*log10(sqs/20.)**2
         PTQQ= .3+.03*log10(sqs/20.)**2
         PTPOM = .23+.03*log10(sqs/20.)**2
      elseif ( IPAR(3).eq.7 ) then
         PTU= PAR(46) + .007*log10(sqs/20.)**2
         PTS= PAR(47) + .007*log10(sqs/20.)**2
         PTQQ= PAR(48) + .03*log10(sqs/20.)**2
         PTPOM = PAR(49) + .03*log10(sqs/20.)**2
      endif
      PPT0 (1) = PTU
      PPT0 (2) = PTU
      PPT0 (3) = PTS
c      PPT0 (10) = .6+.08*log10(sqs/30.)
      PPT0 (10) = PTPOM
      DO J=11,33
         PPT0(J) = PTQQ
      ENDDO
c     charm pt
      IF(IPAR(16).eq.1)THEN
         PPT0(34)=0.5+.08*log10(sqs/30.)
         PPT0(35)=0.5+.08*log10(sqs/30.)
      ELSEIF(IPAR(16).eq.2)THEN
         PPT0(34)=0.3+.08*log10(sqs/30.)
         PPT0(35)=0.5+.08*log10(sqs/30.)
      ELSEIF(IPAR(16).eq.3)THEN
         PPT0(34)=0.3+.5*log10(sqs/30.)
         PPT0(35)=0.5+.5*log10(sqs/30.)
      ELSEIF(IPAR(16).eq.4)THEN
         PPT0(34)=0.365+.473*log10(sqs/30.)
         PPT0(35)=0.144+.473*log10(sqs/30.)
      ELSEIF(IPAR(16).eq.5)THEN
         PPT0(34)=0.3+.1*log10(sqs/30.)
         PPT0(35)=0.5+.1*log10(sqs/30.)
      ELSEIF(IPAR(16).eq.6)THEN
         PPT0(34)=0.303+.125*log10(sqs/30.)
         PPT0(35)=0.5+.125*log10(sqs/30.)
      ELSEIF(IPAR(16).eq.7)THEN
         PPT0(34)=0.308+.165*log10(sqs/30.)
         PPT0(35)=0.5+.165*log10(sqs/30.)
      ELSE
         PPT0(34)=1.0+.08*log10(sqs/30.)
         PPT0(35)=1.5+.08*log10(sqs/30.)
      ENDIF
      
      IF(ndebug.gt.2)THEN
         print*,'PTSETUP: ((u,d),s,diq,pom,cm,cb)',ppt0(1),ppt0(3),
     +        ppt0(11), ppt0(10),ppt0(34),ppt0(35)
      ENDIF

      RETURN
      END

      SUBROUTINE PTSETUP_4FLV(ECM)
C     moved from sib_ndiff to seperate subroutine 
c     so that changes will affect diff. /FR'13

      COMMON /S_DEBUG/ Ncall, Ndebug, Lun    
      COMMON /S_CQDIS2/ PPT0(44)  
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)

      SQS = ECM

c     NA22 piC retune
      PTU=.3+.08*log10(sqs/30.)
      PTS=.45+.08*log10(sqs/30.)
      PTQQ=.6+.08*log10(sqs/30.)
      PTPOM= .6+.08*log10(sqs/30.)
      if ( IPAR(3).eq.1 ) then 
c     pt0
         ptu=.15+.007*log10(sqs/20.)**2
         pts=.3+.007*log10(sqs/20.)**2
         ptqq=.3+.03*log10(sqs/20.)**2
         ptpom= .6+.08*log10(sqs/30.)
      elseif ( ipar(3).eq.2 ) then
C     pt1
         ptu=.15+.007*log10(sqs/20.)**2
         pts=.32+.007*log10(sqs/20.)**2
         ptqq=.4+.007*log10(sqs/20.)**2
         ptpom= .6+.08*log10(sqs/30.)
c     pt2
      elseif ( ipar(3).eq.3 ) then
         ptu=.17+.007*log10(sqs/20.)**2
         pts=.3+.007*log10(sqs/20.)**2
         ptqq=.3+.03*log10(sqs/20.)**2
         ptpom = .6+.08*log10(sqs/30.)
      elseif ( ipar(3).eq.5 ) then
         PTU=.16+.007*log10(sqs/20.)**2
         PTS=.28+.007*log10(sqs/20.)**2
         PTQQ= .3+.03*log10(sqs/20.)**2
         PTPOM = .23+.03*log10(sqs/20.)**2
      elseif ( IPAR(3).eq.6 ) then
         PTU=.16+.007*log10(sqs/20.)**2
         PTS=.28+.007*log10(sqs/20.)**2
         PTQQ= .3+.03*log10(sqs/20.)**2
         PTPOM = .23+.03*log10(sqs/20.)**2
      elseif ( IPAR(3).eq.7 ) then
         PTU= PAR(46) + .007*log10(sqs/20.)**2
         PTS= PAR(47) + .007*log10(sqs/20.)**2
         PTQQ= PAR(48) + .03*log10(sqs/20.)**2
         PTPOM = PAR(49) + .03*log10(sqs/20.)**2
      elseif ( IPAR(3).eq.8 ) then
         PTU= PAR(46) + PAR(68)*log10(sqs/20.)**2
         PTS= PAR(47) + PAR(70)*log10(sqs/20.)**2
         PTQQ= PAR(48) + PAR(69)*log10(sqs/20.)**2
         PTPOM = PAR(49) + PAR(69)*log10(sqs/20.)**2
         PTSEA = PAR(67) + PAR(69)*log10(sqs/20.)**2
      endif
      PPT0 (1) = PTU
      PPT0 (2) = PTU
      PPT0 (3) = PTS
c      PPT0 (10) = .6+.08*log10(sqs/30.)
      PPT0 (10) = PTPOM
      DO J=11,33
         PPT0(J) = PTQQ
      ENDDO
      PPT0 (20) = PTSEA
c     charm pt
      IF(IPAR(16).eq.1)THEN
         PTCHM=0.5+.08*log10(sqs/30.)
         PTCHB=0.5+.08*log10(sqs/30.)
      ELSEIF(IPAR(16).eq.2)THEN
         PTCHM=0.3+.08*log10(sqs/30.)
         PTCHB=0.5+.08*log10(sqs/30.)
      ELSEIF(IPAR(16).eq.3)THEN
         PTCHM=0.3+.5*log10(sqs/30.)
         PTCHB=0.5+.5*log10(sqs/30.)
      ELSEIF(IPAR(16).eq.4)THEN
         PTCHM=0.365+.473*log10(sqs/30.)
         PTCHB=0.144+.473*log10(sqs/30.)
      ELSEIF(IPAR(16).eq.5)THEN
         PTCHM=0.3+.1*log10(sqs/30.)
         PTCHB=0.5+.1*log10(sqs/30.)
      ELSEIF(IPAR(16).eq.6)THEN
         PTCHM=0.303+.125*log10(sqs/30.)
         PTCHB=0.5+.125*log10(sqs/30.)
      ELSEIF(IPAR(16).eq.7)THEN
         PTCHM=0.308+.165*log10(sqs/30.)
         PTCHB=0.5+.165*log10(sqs/30.)
      ELSE
         PTCHM=1.0+.08*log10(sqs/30.)
         PTCHB=1.5+.08*log10(sqs/30.)
      ENDIF
      PPT0(4) = PTCHM
      PPT0(14) = PTCHB
      PPT0(24) = PTCHB
      DO J=34,44
         PPT0(J) = PTCHB
      ENDDO
     
      IF(ndebug.gt.2)THEN
         print*,'PTSETUP2: ((u,d),s,diq,pom,cm,cb)',ppt0(1),ppt0(3),
     +        ppt0(11), ppt0(10),ppt0(4),ppt0(34)
      ENDIF

      RETURN
      END

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
      SAVE
      COMMON /S_PLIST/ PLIST(8000,5), LLIST(8000), NP
      DIMENSION ROT(3,3),PV(3)
      DOUBLE PRECISION DP(4),DBEX,DBEY,DBEZ,DGA,DBEP,DGABEP
      IF(THE**2+PHI**2 .LE. 1E-20) GO TO 131
C...ROTATE (TYPICALLY FROM Z AXIS TO DIRECTION THETA,PHI)
       ROT(1,1)=COS(THE)*COS(PHI)
       ROT(1,2)=-SIN(PHI)
       ROT(1,3)=SIN(THE)*COS(PHI)
       ROT(2,1)=COS(THE)*SIN(PHI)
       ROT(2,2)=COS(PHI)
       ROT(2,3)=SIN(THE)*SIN(PHI)
       ROT(3,1)=-SIN(THE)
       ROT(3,2)=0.
       ROT(3,3)=COS(THE)
       DO 120 I=NBEG,NEND
       DO 100 J=1,3
 100   PV(J)=PLIST(I,J)
       DO 110 J=1,3
 110   PLIST(I,J)=ROT(J,1)*PV(1)+ROT(J,2)*PV(2)+ROT(J,3)*PV(3)
 120   CONTINUE
 131    IF(DBEX**2+DBEY**2+DBEZ**2 .LE. 1D-20) GO TO 151
C...LORENTZ BOOST (TYPICALLY FROM REST TO MOMENTUM/ENERGY=BETA)
       DGA=1D0/DSQRT(1D0-DBEX**2-DBEY**2-DBEZ**2)
       DO 140 I=NBEG, NEND
       DO 130 J=1,4
 130   DP(J)=PLIST(I,J)
       DBEP=DBEX*DP(1)+DBEY*DP(2)+DBEZ*DP(3)
       DGABEP=DGA*(DGA*DBEP/(1D0+DGA)+DP(4))
       PLIST(I,1)=DP(1)+DGABEP*DBEX
       PLIST(I,2)=DP(2)+DGABEP*DBEY
       PLIST(I,3)=DP(3)+DGABEP*DBEZ
       PLIST(I,4)=DGA*(DP(4)+DBEP)
 140   CONTINUE
 151   RETURN
      END


      SUBROUTINE BEAM_SPLIT (L, NW, XX, IFL, XJET, LXBAD)
C...This subroutine split a hadron of code L
C.  into 2*NW partons, each of energy XX(j) and
C.  flavor IFL.  The minimum fractional energy of 
C.  each parton is X_min = 2*STR_mass/sqrt(s)
C.  
C.  Variable qmas changed to STR_mass to agree with name in SIBYLL
C.      and added to calling sequenceto insure symmetry.
C.     Also a factor of (1-xjet) is added to the def. of xmin for nw=1
C.                               RSF  Apr-2-92
C---------------------------------------------------------------------
      SAVE

      PARAMETER (NW_max = 20)
      DIMENSION XX(2*NW_max), IFL(2*NW_max)

      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, kb ,kt
      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      COMMON /S_CUTOFF/ STR_mass_val, STR_mass_val_hyp, STR_mass_sea

      DATA AC /-0.2761856692/             ! log(2) - gamma(Eulero)
      DATA GAMMA /2./
      DATA NBAD / 0 /

      STR_mass = STR_mass_val
      IF(ABS(L).gt.14) STR_mass = STR_mass_val_hyp
c-------
c  New code to handle low energy p nuc problem.
c------
      LXBAD = 0
      XMIN = 2.*STR_mass/SQS
      IF (1.-XJET .LT. FLOAT(2*NW)*XMIN)  THEN
         NBAD = NBAD + 1
         LXBAD = 1
         IF(ndebug.gt.0)THEN
            IF (NBAD .LE. 20) THEN
               WRITE (LUN, *) 
     &              'BEAM_SPLIT: kinematically forbidden situation'
               WRITE (LUN, 5)  NBAD, SQS, XJET, NW
            ENDIF
 5          FORMAT(1X,'NBAD = ',I3,3X,'sqs = ',E10.3,
     &           3X, 'x_jet = ', F9.3, 3X, ' NW = ',I2)
            IF (NBAD .eq. 20) THEN
               WRITE (LUN, *) 
     &              ' BEAM_SPLIT : Last warning about bad splittings '
               WRITE(LUN, *)' The energy threshold is probably too low.'
            ENDIF
         ENDIF
         RETURN
      ENDIF

      IF (NW .EQ. 1)  THEN
         XVAL = 1.-XJET
         GOTO 200
      ENDIF

C...Choose total energy of sea partons
      N = 2*(NW-1)
      Z1 = LOG(FLOAT(N))
      Z2 = LOG(0.5*SQS*(1.-XJET)/STR_mass-2.)
100   R=S_RNDM(0)
      Z=(Z1+AC)*(1.+R*(((Z2+AC)/(Z1+AC))**N-1.))**(1./FLOAT(N))-AC
      XSEA = XMIN*EXP(Z)
      IF ( (1.-XSEA)**GAMMA .LT. S_RNDM(0)) GOTO 100
C...Split the energy  of sea partons among the different partons
      XREM = XSEA - FLOAT(N)*XMIN
      DO J=3,N+1
         XA = XREM*S_RNDM(0)
         XREM = XREM - XA
*        print *,' BEAM_SPLIT: XX index ',J
         XX(J) = XMIN + XA
      ENDDO
*     print *,' BEAM_SPLIT: XX index ',N+2
      XX(N+2) = XMIN + XREM
      XVAL = 1.-XSEA-XJET
C...Flavor of sea partons
      DO J=1,N/2
         J1 =  3 + (J-1)*2 
*        print *,' BEAM_SPLIT: flavour indices ',J1,J1+1
         IFL(J1) = INT(1.+1.99*S_RNDM(0))
         IFL(J1+1) = -IFL(J1)
      ENDDO
C...Prepare the valence partons
200   CALL HSPLI (L,IFL(1),IFL(2))
      IF(IFL(1).eq.0.and.IFL(2).eq.0)THEN
         WRITE(LUN,*) 'BEAM_SPLIT: HSPLI rejection!'
         WRITE(LUN,*) ' RUN (SQS,KB,KT):',SQS,KB,KT
         STOP
      ENDIF
      CHI = CHIDIS(L,IFL(1),IFL(2))
      XX(1) = MAX(CHI*XVAL,XMIN)
      XX(1) = MIN(XX(1),XVAL-XMIN)
C      FOR MESONS, SPLIT ENERGY SYMETRICALLY.
C????? SPLIT K'S WITH ENERGY TO S QUARK? 
C
      if (abs(l).le.12.and.S_RNDM(0).le.0.5) xx(1)=XVAL-XX(1)
      XX(2) = XVAL-XX(1)

      END

      SUBROUTINE pho_cpcini(Nrows,Number,List)
C***********************************************************************
C
C     initialization of particle hash table
C
C     input:   Number     vector with Nrows entries according to PDG
C                         convention
C
C     output:  List       vector with hash table
C
C     (this code is based on the function initpns written by
C      Gerry Lynch, LBL, January 1990)
C
C***********************************************************************
      IMPLICIT NONE
      SAVE
      integer Number(*),List(*),Nrows
      Integer Nin,Nout,Ip,I

      do I = 1,577
        List(I) = 0
      enddo
C    Loop over all of the elements in the Number vector

        Do 500 Ip = 1,Nrows
            Nin = Number(Ip)

C    Calculate a list number for this particle id number
            If(Nin.Gt.99999.or.Nin.Le.0) Then
                 Nout = -1
            Else If(Nin.Le.577) Then
                 Nout = Nin
            Else
                 Nout = Mod(Nin,577)
            End If

 200        continue

            If(Nout.Lt.0) Then
C    Count the bad entries
c                Write(6,'(1x,a,i10)')
c     &            'pho_cpcini: invalid particle ID',Nin
                Go to 500
            End If
            If(List(Nout).eq.0) Then
                List(Nout) = Ip
            Else
                If(Nin.eq.Number(List(Nout))) Then
c                  Write(6,'(1x,a,i10)')
c     &              'pho_cpcini: double particle ID',Nin
                End If
                Nout = Nout + 5
                If(Nout.Gt.577) Nout = Mod(Nout, 577)

                Go to 200
            End If
 500      Continue

      END

C--------------------------------------------------------------------------
C    CODE OF ANALYSIS (not needed to generate events)
C--------------------------------------------------------------------------

      SUBROUTINE PFsum(N1,N2,ETOT,PXT,PYT,PZT,NF)
C...Return the energy,px,py,pz and the number of stable
C.  particles in the list between N1 and N2
      SAVE
      COMMON /S_PLIST/ P(8000,5), LLIST(8000), NP
      NF=0
      ETOT=0.
      PXT=0.
      PYT=0.
      PZT=0.
      DO J=N1,N2
         L = LLIST(J)     
         IF (IABS(L) .LT. 10000)  THEN
           NF = NF+1
           ETOT = ETOT + ABS( P(J,4) )
           PXT = PXT + P(J,1)
           PYT = PYT + P(J,2)
           PZT = PZT + P(J,3)
         ENDIF
      ENDDO
      RETURN
      END


      SUBROUTINE QNUM (JQ,JS,JC,JB,JBA, NC, NF)
C...Return the quantum numbers of one event
C.  JQ = charge, JB = baryon number, JS = strangeness, JC = charmedness
C.  JBA = (number of baryons+antibaryons)
C.  NC  = number of charged particles
C.  NF  = number of final particles
C..................................................
      SAVE
      COMMON /S_PLIST/ P(8000,5), LLIST(8000), NP
      COMMON /S_CHP/ ICHP(99), ISTR(99), IBAR(99)
      COMMON /S_CHM/ ICHM(99)
      JQ = 0
      JB = 0
      JS = 0
      JC = 0
      JBA= 0
      NC = 0
      NF = 0
      DO J=1,NP
          L = LLIST(J)
          LL = IABS(L)
          IF (LL .LT. 10000)  THEN
              IF(ICHP(LL) .NE. 0) NC = NC + 1
              NF = NF + 1
              JQ = JQ + ICHP(LL)*ISIGN(1,L)
              JB = JB + IBAR(LL)*ISIGN(1,L)
              JBA= JBA+ IBAR(LL)
              JS = JS + ISTR(LL)*ISIGN(1,L)
              JC = JC + ICHM(LL)*ISIGN(1,L)
          ENDIF
      ENDDO
      RETURN
      END


      SUBROUTINE SIB_LIST(LUN)
C-----------------------------------------------------------------------
C...This routine prints the event record for the
C.  current event on unit LUN
C-----------------------------------------------------------------------
      SAVE

      COMMON /S_PLIST/ P(8000,5), LLIST(8000), NP
      COMMON /S_PLIST1/ LLIST1(8000)
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, kb ,kt
      COMMON /S_DEBUG/ Ncall, Ndebug, Lunn
      PARAMETER (NW_max = 20)
      PARAMETER (NS_max = 20, NH_max = 80)
      PARAMETER (NJ_max = (NS_max+NH_max)*NW_max)
      COMMON /S_CHIST/ X1J(NJ_max),X2J(NJ_max),
     &    X1JSUM(NW_max),X2JSUM(NW_max),PTJET(NJ_max),PHIJET(NJ_max),
     &    NNPJET(NJ_max),NNPSTR(2*NW_max),NNSOF(NW_max),NNJET(NW_max),
     &    JDIF(NW_max),NW,NJET,NSOF
      COMMON /S_CCSTR/ X1(2*NW_max),X2(2*NW_max),
     &    PXB(2*NW_max),PYB(2*NW_max),PXT(2*NW_max),PYT(2*NW_max),
     &    IFLB(2*NW_max),IFLT(2*NW_max)
      CHARACTER*6 NAMP
      COMMON /S_CNAM/ NAMP (0:99)
      COMMON /S_CHP/ ICHP(99), ISTR(99), IBAR(99)
      COMMON /S_CHM/ ICHM(99)
      COMMON /S_PARTO/ NFORIG(8000), NPORIG(8000), IPFLAG, NINT

      CHARACTER CODE*18
      CHARACTER*18 NAMDIF(0:3)
      DATA NAMDIF /'Non-diff. event   ',
     &  'Beam diffraction  ','Target diffraction','Double diffraction'/

      WRITE (LUN,*)
      WRITE (LUN, *) ' Event record '
      if(NW.eq.1) WRITE (LUN,*) '  ',NAMDIF(JDIF(1))
      WRITE (LUN,*) '  N_w/N_s/N_j = ', NW, NSOF, NJET
      WRITE (LUN,100)

C...Print particle list
      ichar = 0
      ibary = 0
      ichmd = 0
      istrg = 0
      DO J=1,NP
        L = MOD(LLIST(J),10000)
        CODE = '                  '
        CODE(1:6) = NAMP(IABS(L))
        IF (L .LT. 0) CODE(7:9) = 'bar'
        IF(IABS(LLIST(J)) .GT. 10000)   CODE(10:10) = '*'
        WRITE (LUN,120) J, CODE, LLIST1(J), NPORIG(J), (P(J,K),K=1,4)
        if(abs(LLIST(J)).LT.10000) then
          ichar = ichar+sign(1,l)*ICHP(iabs(l))
          ibary = ibary+sign(1,l)*IBAR(iabs(l))
          ichmd = ichmd+sign(1,l)*ICHM(iabs(l))
          istrg = istrg+sign(1,l)*ISTR(iabs(l))
        endif
      ENDDO
      CALL PFsum(1,NP,Esum,PXsum,PYsum,PZsum,NF)
      WRITE(LUN,140) PXsum,PYsum,PZsum,Esum
100      FORMAT(3X,'N  Particle',12X,'Ori',2x,'Proc',6x,'PX',9x,'PY',9x,
     +         'PZ',9x,'E', /, 3X,75('-'))
120      FORMAT(1X,I4,1X,A18,1X,I4,1X,I4,2X,2(F9.3,2X),2(E9.3,2X))
140      FORMAT(3X,75('-'),/,1X,'Tot = ',29X,2(F9.3,2X),G9.3,2X,E9.3)
      write(LUN,'(1x,a,i3,3x,a,i3))') 'Total charge:',ichar,
     &  'baryon number:',ibary
      write(LUN,'(1x,a,i3,3x,a,i3))') 'Total strangeness:',istrg,
     &  'charm number:',ichmd

      END



      SUBROUTINE KCODE (J,CODE,NC)
C...Produce the code for parton J
C.  Input K, Output CODE, NC=number of characters
C..................................................
      SAVE
      CHARACTER*5 CODE
      CHARACTER*1 NAMQ(3)
      DATA NAMQ /'U','D','S'/
      CODE = '     '
      IF(J.EQ.0)  THEN
         CODE(1:3) = 'GLU'
         NC = 3
         RETURN
      ENDIF
      JA = IABS(J)
      J1 = MOD(JA,10)
      J2 = (JA-J1)/10
      IF(JA .GT. 10) THEN
         CODE(1:1) = NAMQ(J2)
         CODE(2:2) = NAMQ(J1)
         NC = 2
      ELSE
         CODE(1:1) = NAMQ(J1)
         NC = 1      
      ENDIF
      IF (J .LT. 0)  THEN
         CODE(NC+1:NC+3) = 'bar'
         NC = NC+3
      ENDIF
      RETURN
      END

      SUBROUTINE SIB_ABRT(CWARN)
      CHARACTER*25 CWARN
      
      WRITE (6,*) CWARN , ' aborting..'
      a = -1.
      a = log(a)
      STOP

      END

C----------------------------------------------------------------------------
C  Code for sampling
C-----------------------------------------------------------------------------

      SUBROUTINE SAMPLE_soft (STR_mass_min, X1,X2,PT)
C-----------------------------------------------------------------------
C...Routine for the sampling the kinematical variables
C.  that characterize a soft cut pomeron (x1,x2, pT)
C.  from the differential cross section:
C.     d3sigma/(dx1 dx2 dpT)
C.  INPUT:  L=1 incident proton, L=2  incident pi
C.          (soft strings identical for pi and p interactions)
C.  OUTPUT:  X1, X2, PT (GeV)
C-----------------------------------------------------------------------
      SAVE

      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, kb ,kt
      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      COMMON /S_CQDIS/ PPT0(35),ptflag
      COMMON /S_CQDIS2/ PPT02(44)
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)

      ZSOF = 2.*LOG(STR_mass_min/SQS)
 100  Z1=-ZSOF*S_RNDM(0)+ZSOF
      Z2=-ZSOF*S_RNDM(0)+ZSOF
      IF(Z1+Z2.LE.ZSOF) GOTO 100
      X1=EXP(Z1)
      X2=EXP(Z2)
      STR_mass2 = sqrt(X1*X2*S)/2.
      PPTT = PPT0(10)
      IF(IPAR(17).GT.0)then
         PPTT = PPT02(10)
         IF(IPAR(3).eq.8) PPTT = PPT02(20)
      endif
 150  PT = PPTT*SQRT(-ALOG(MAX(1E-10,S_RNDM(0))))
      IF(IPAR(3).ge.6)THEN
         XM = .0
         XM2 = XM**2
         RNDM = MAX(1E-10,S_RNDM(IFL))
         XMT = PPTT * ALOG(RNDM) - XM
         XMT2 = XMT**2
         PT = SQRT(XMT2-XM2)
      ENDIF
      IF(PT.GT.PTmin) GOTO 150
      IF(PT.GE.STR_mass2) GOTO 150

      END

      SUBROUTINE SAMPLE_soft3 (STR_mass_min, X1,X2,PT)
C-----------------------------------------------------------------------
C...Routine for the sampling the kinematical variables
C.  that characterize a soft cut pomeron (x1,x2, pT)
C.  from the differential cross section:
C.     d3sigma/(dx1 dx2 dpT)
C.  INPUT:  L=1 incident proton, L=2  incident pi
C.          (soft strings identical for pi and p interactions)
C.  OUTPUT:  X1, X2, PT (GeV)
C-----------------------------------------------------------------------
      SAVE

      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, kb ,kt
      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      COMMON /S_CQDIS/ PPT0(35),ptflag
      COMMON /S_CQDIS2/ PPT02(44)
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)

      SLOPE = max(1.,PAR(42))
      ZSOF = 2.*LOG(STR_mass_min/SQS) ! minim. mass ~ x1 * x2
 100  Z1=-ZSOF*S_RNDM(0)+ZSOF   ! sample envelope 1/x
      X1 = EXP(Z1)
c      print *,'z1,x1:',z1,x1
      XR = log(1.-X1) - log(1.-EXP(ZSOF))
c      print *,'ratio:',(1.-X1)/(1.-EXP(ZSOF)),(1.-X1),1.-EXP(ZSOF)
c      print *,'log ratio:',xr,log(1.-X1),log(1.-EXP(ZSOF))
      if(SLOPE*XR.le.log(S_RNDM(0))) goto 100

 200  Z2=-ZSOF*S_RNDM(0)+ZSOF   ! sample envelope 1/x
      X2 = EXP(Z2)
      XR = log(1.-X2) - log(1.-EXP(ZSOF))
      if(SLOPE*XR.le.log(S_RNDM(0))) goto 200     
c      print *,'zsof,z1,z2',zsof,z1,z2
      IF(Z1+Z2.LE.ZSOF) GOTO 100
      STR_mass2 = sqrt(X1*X2*S)/2.
      PPTT = PPT0(10)
      IF(IPAR(17).GT.0)then
         PPTT = PPT02(10)
         IF(IPAR(3).eq.8) PPTT = PPT02(20)
      endif
 150  PT = PPTT*SQRT(-ALOG(MAX(1E-10,S_RNDM(0))))
      IF(IPAR(3).ge.6)THEN
         XM = .0
         XM2 = XM**2
         RNDM = MAX(1E-10,S_RNDM(IFL))
         XMT = PPTT * ALOG(RNDM) - XM
         XMT2 = XMT**2
         PT = SQRT(XMT2-XM2)
      ENDIF
      IF(PT.GT.PTmin) GOTO 150
      IF(PT.GE.STR_mass2) GOTO 150

      END


      SUBROUTINE SAMPLE_hard (L, X1,X2,PT)
C-----------------------------------------------------------------------
C...Routine for the sampling the kinematical variables 
C.  that determine a  jet-jet  system (x1,x2, pT) 
C.  from the differential cross section:
C.     d3sigma/(dx1 dx2 dpT)
C.  This version assumes the `single parton approximation'
C.  INPUT:  L=1 incident proton, L=2  incident pi
C.  OUTPUT:  X1, X2, PT (GeV)
C-----------------------------------------------------------------------
      SAVE

      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, kb ,kt
      COMMON /S_DEBUG/ Ncall, Ndebug, Lun

100   Z1=ZSAMPLE (ZMIN,L)
      Z2=ZSAMPLE (ZMIN,1)
      SIG=1.-XMIN*EXP(-Z1-Z2)
      IF (SIG .LT. S_RNDM(0))  GOTO 100      
      X1=EXP(Z1)
      X2=EXP(Z2)
      IF (X1.gt.0.9.or.X2.gt.0.9) GOTO 100
      Q2=PTmin**2/(1.-S_RNDM(0)*SIG)
      PT=SQRT(Q2*(1.-Q2/(S*X1*X2)))
      RETURN
      END


      FUNCTION ZSAMPLE (ZMIN,L)
C...This function returns as output a value z=log(x)
C.  distributed as f(x) = g(x) + 4/9 *(q(x) + qbar(x))
C.  from a minimum value ZMIN to 0,
C.  for a proton (L=1) or a pi (L=2)
C.  needs to be initialised with: CALL ZSAMPLE_INI
C.....................................................
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)
      COMMON /S_CZGEN/ XA(2),XB(2),XMAX,ZA(2),ZB(2),ZMAX,
     +     DX(2),DZ(2),APART(2),FFA(2),FFB(2),
     +     DFX(2),DFZ(2),XX(200,2),ZZ(200,2),FFX(200,2),FFZ(200,2),
     +     NX,NZ
      parameter (b=0.268)
      parameter (bpi=3.7)
      parameter (cpi=0.698)
      SAVE

      F = PART_INT(ZMIN,L)*S_RNDM(0)
      IF (F .GE. FFA(L))  THEN
         IF(IPAR(8).EQ.0)THEN
            ZSAMPLE = ZA(L) - (F-FFA(L))/APART(L)
         ELSE
            if(L.eq.1) then
               ZSAMPLE = -1./b * LOG( 1. - F / APART(L) ) 
            else
               ZSAMPLE = -1. * ( (F - cpi)/APART(L) ) ** (1/bpi)
            endif
         ENDIF
      ELSE IF (F .GE. FFB(L))  THEN
         JF = (F-FFB(L))/DFZ(L) + 1
         JF = min(JF,199)
         F0 = FFB(L) + DFZ(L)*FLOAT(JF-1)
         T = (F-F0)/DFZ(L)
         ZSAMPLE = ZZ(JF,L)*(1.-T)+ZZ(JF+1,L)*T
      ELSE
         JF = F/DFX(L)+1
         JF = min(JF,199)
         F0 = DFX(L)*FLOAT(JF-1)
         T = (F-F0)/DFX(L)
         X = XX(JF,L)*(1.-T)+XX(JF+1,L)*T
         ZSAMPLE = LOG(X)
      ENDIF

      RETURN
      END


      FUNCTION PART_INT (ZMIN,L)
C...This function returns as output the integral of
C.  the parton structure function:
C.     f(x) = g(x) + 4/9 *(q(x) + qbar(x))
C.  from xmin = exp(zmin) to 1 
C.  for a proton (L=1) or a pi (L=2)
C.  needs to be initialised with: CALL ZSAMPLE_INI
C.....................................................
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)
      COMMON /S_CZGEN/ XA(2),XB(2),XMAX,ZA(2),ZB(2),ZMAX,
     +     DX(2),DZ(2),APART(2),FFA(2),FFB(2),
     +     DFX(2),DFZ(2),XX(200,2),ZZ(200,2),FFX(200,2),FFZ(200,2),
     +     NX,NZ
      parameter (b=0.268)
      parameter (bpi=3.7)
      parameter (cpi=0.698)
      SAVE

      IF (ZMIN .LT. ZA(L))  THEN
         IF(IPAR(8).EQ.0)THEN
            PART_INT = FFA(L) + APART(L) * (ZA(L) - ZMIN)
         ELSE
            if(L.eq.1) then
               PART_INT = APART(L) * ( 1. - EXP(-b*ZMIN) ) 
            else
               PART_INT = APART(L) * ( -ZMIN ) ** bpi + cpi
            endif
         ENDIF
      ELSE IF (ZMIN .LT. ZB(L)) THEN
         JZ = (ZB(L)-ZMIN)/DZ(L)+1
         JZ = min(JZ,199)
         Z0 = ZB(L)-DZ(L)*FLOAT(JZ-1)
         T = (Z0-ZMIN)/DZ(L)
         PART_INT = FFZ(JZ,L)*(1.-T) + FFZ(JZ+1,L)*T

      ELSE
         X = EXP(ZMIN)
         JX = (XMAX-X)/DX(L)+1
         JX = min(JX,199)
         X0 = XMAX-DX(L)*FLOAT(JX-1)
         T = (X0-X)/DX(L)
         PART_INT = FFX(JX,L)*(1.-T) + FFX(JX+1,L)*T
      
      ENDIF
      RETURN
      END


      BLOCK DATA PDFINI
C..   tabled parton distribution function
c     Proton: GRV98LO , Eur.Phys.J. C5(1998) 461-470
c     Pion:   GRV91 , Z. Phys. C53, 651-655 (1992)
      
      SAVE
      COMMON /S_CZGEN/ XA(2),XB(2),XMAX,ZA(2),ZB(2),ZMAX,
     +     DX(2),DZ(2),APART(2),FFA(2),FFB(2),
     +     DFX(2),DFZ(2),XX(200,2),ZZ(200,2),FFX(200,2),FFZ(200,2),
     +     NX,NZ
      DATA XA /1e-06,0.0001/
      DATA XB /0.1,0.1/
      DATA XMAX /0.800000011921/
      DATA ZMAX /-0.223143532872/
      DATA NX /200/
      DATA NZ /200/
      DATA ZA /-13.8155,-9.21034/
      DATA ZB /-2.30259,-2.30259/
      DATA DX /0.00351759,0.00351759/
      DATA DZ /0.0578539,0.0347123/
      DATA DFX /0.00952501,0.00847474/
      DATA DFZ /1.93863,0.326082/
      DATA APART /-9.80215,0.0178207/
      DATA FFA /387.684,66.5767/
      DATA FFB /1.89548,1.68647/
      
      DATA (FFX(K,1),K=1,200 ) /
     &0.000E+00,6.380E-05,1.315E-04,2.034E-04,2.795E-04,
     &3.601E-04,4.454E-04,5.356E-04,6.309E-04,7.315E-04,
     &8.377E-04,9.497E-04,1.068E-03,1.192E-03,1.323E-03,
     &1.460E-03,1.605E-03,1.756E-03,1.916E-03,2.083E-03,
     &2.258E-03,2.441E-03,2.633E-03,2.835E-03,3.045E-03,
     &3.265E-03,3.496E-03,3.736E-03,3.988E-03,4.250E-03,
     &4.524E-03,4.810E-03,5.108E-03,5.418E-03,5.742E-03,
     &6.078E-03,6.429E-03,6.794E-03,7.174E-03,7.570E-03,
     &7.981E-03,8.408E-03,8.852E-03,9.313E-03,9.793E-03,
     &1.029E-02,1.081E-02,1.134E-02,1.190E-02,1.247E-02,
     &1.307E-02,1.369E-02,1.433E-02,1.500E-02,1.568E-02,
     &1.640E-02,1.714E-02,1.790E-02,1.869E-02,1.951E-02,
     &2.035E-02,2.123E-02,2.213E-02,2.307E-02,2.403E-02,
     &2.503E-02,2.607E-02,2.713E-02,2.823E-02,2.937E-02,
     &3.054E-02,3.176E-02,3.301E-02,3.430E-02,3.563E-02,
     &3.701E-02,3.842E-02,3.989E-02,4.139E-02,4.295E-02,
     &4.455E-02,4.620E-02,4.791E-02,4.966E-02,5.147E-02,
     &5.334E-02,5.526E-02,5.724E-02,5.927E-02,6.137E-02,
     &6.353E-02,6.576E-02,6.805E-02,7.041E-02,7.284E-02,
     &7.534E-02,7.791E-02,8.056E-02,8.329E-02,8.609E-02,
     &8.898E-02,9.195E-02,9.500E-02,9.814E-02,1.014E-01,
     &1.047E-01,1.081E-01,1.116E-01,1.153E-01,1.190E-01,
     &1.228E-01,1.267E-01,1.308E-01,1.350E-01,1.392E-01,
     &1.436E-01,1.481E-01,1.528E-01,1.575E-01,1.624E-01,
     &1.674E-01,1.725E-01,1.778E-01,1.832E-01,1.888E-01,
     &1.946E-01,2.005E-01,2.066E-01,2.128E-01,2.193E-01,
     &2.259E-01,2.327E-01,2.397E-01,2.469E-01,2.543E-01,
     &2.619E-01,2.698E-01,2.778E-01,2.862E-01,2.947E-01,
     &3.035E-01,3.125E-01,3.218E-01,3.314E-01,3.413E-01,
     &3.514E-01,3.618E-01,3.726E-01,3.836E-01,3.950E-01,
     &4.067E-01,4.188E-01,4.312E-01,4.440E-01,4.572E-01,
     &4.708E-01,4.848E-01,4.992E-01,5.141E-01,5.294E-01,
     &5.452E-01,5.615E-01,5.783E-01,5.956E-01,6.134E-01,
     &6.319E-01,6.509E-01,6.706E-01,6.909E-01,7.118E-01,
     &7.334E-01,7.558E-01,7.789E-01,8.029E-01,8.276E-01,
     &8.532E-01,8.797E-01,9.072E-01,9.356E-01,9.650E-01,
     &9.956E-01,1.027E+00,1.060E+00,1.094E+00,1.130E+00,
     &1.167E+00,1.205E+00,1.245E+00,1.287E+00,1.331E+00,
     &1.376E+00,1.423E+00,1.473E+00,1.525E+00,1.579E+00,
     &1.636E+00,1.696E+00,1.759E+00,1.826E+00,1.895E+00/
      
      DATA (FFX(K,2),K=1,200 ) /
     &0.000E+00,7.266E-04,1.470E-03,2.231E-03,3.009E-03,
     &3.805E-03,4.619E-03,5.450E-03,6.300E-03,7.168E-03,
     &8.055E-03,8.961E-03,9.886E-03,1.083E-02,1.179E-02,
     &1.278E-02,1.378E-02,1.481E-02,1.585E-02,1.692E-02,
     &1.800E-02,1.911E-02,2.024E-02,2.139E-02,2.256E-02,
     &2.376E-02,2.498E-02,2.622E-02,2.748E-02,2.877E-02,
     &3.008E-02,3.142E-02,3.278E-02,3.416E-02,3.557E-02,
     &3.701E-02,3.847E-02,3.996E-02,4.147E-02,4.301E-02,
     &4.458E-02,4.617E-02,4.779E-02,4.945E-02,5.112E-02,
     &5.283E-02,5.457E-02,5.634E-02,5.813E-02,5.996E-02,
     &6.182E-02,6.371E-02,6.563E-02,6.759E-02,6.957E-02,
     &7.159E-02,7.365E-02,7.573E-02,7.786E-02,8.001E-02,
     &8.221E-02,8.443E-02,8.670E-02,8.900E-02,9.134E-02,
     &9.372E-02,9.614E-02,9.860E-02,1.011E-01,1.036E-01,
     &1.062E-01,1.088E-01,1.115E-01,1.142E-01,1.170E-01,
     &1.197E-01,1.226E-01,1.255E-01,1.284E-01,1.314E-01,
     &1.344E-01,1.375E-01,1.406E-01,1.438E-01,1.470E-01,
     &1.503E-01,1.536E-01,1.570E-01,1.605E-01,1.640E-01,
     &1.675E-01,1.712E-01,1.748E-01,1.786E-01,1.824E-01,
     &1.862E-01,1.901E-01,1.941E-01,1.982E-01,2.023E-01,
     &2.065E-01,2.107E-01,2.151E-01,2.195E-01,2.239E-01,
     &2.285E-01,2.331E-01,2.378E-01,2.426E-01,2.474E-01,
     &2.524E-01,2.574E-01,2.625E-01,2.677E-01,2.730E-01,
     &2.784E-01,2.839E-01,2.895E-01,2.951E-01,3.009E-01,
     &3.068E-01,3.128E-01,3.189E-01,3.251E-01,3.314E-01,
     &3.378E-01,3.443E-01,3.510E-01,3.578E-01,3.647E-01,
     &3.717E-01,3.789E-01,3.862E-01,3.937E-01,4.012E-01,
     &4.090E-01,4.169E-01,4.249E-01,4.331E-01,4.415E-01,
     &4.500E-01,4.587E-01,4.676E-01,4.767E-01,4.859E-01,
     &4.954E-01,5.050E-01,5.148E-01,5.249E-01,5.352E-01,
     &5.457E-01,5.564E-01,5.674E-01,5.786E-01,5.901E-01,
     &6.019E-01,6.139E-01,6.262E-01,6.388E-01,6.517E-01,
     &6.649E-01,6.785E-01,6.923E-01,7.066E-01,7.212E-01,
     &7.362E-01,7.516E-01,7.673E-01,7.836E-01,8.002E-01,
     &8.174E-01,8.350E-01,8.532E-01,8.718E-01,8.911E-01,
     &9.109E-01,9.313E-01,9.524E-01,9.742E-01,9.966E-01,
     &1.020E+00,1.044E+00,1.069E+00,1.094E+00,1.121E+00,
     &1.149E+00,1.177E+00,1.207E+00,1.238E+00,1.271E+00,
     &1.304E+00,1.339E+00,1.376E+00,1.414E+00,1.454E+00,
     &1.496E+00,1.540E+00,1.586E+00,1.635E+00,1.686E+00/
      
      DATA (FFZ(K,1),K=1,200 ) /
     &1.895E+00,2.014E+00,2.137E+00,2.263E+00,2.393E+00,
     &2.527E+00,2.665E+00,2.807E+00,2.953E+00,3.103E+00,
     &3.257E+00,3.417E+00,3.580E+00,3.748E+00,3.921E+00,
     &4.098E+00,4.281E+00,4.469E+00,4.663E+00,4.861E+00,
     &5.065E+00,5.274E+00,5.489E+00,5.710E+00,5.937E+00,
     &6.170E+00,6.409E+00,6.654E+00,6.906E+00,7.164E+00,
     &7.430E+00,7.702E+00,7.981E+00,8.267E+00,8.561E+00,
     &8.862E+00,9.171E+00,9.487E+00,9.811E+00,1.014E+01,
     &1.048E+01,1.083E+01,1.119E+01,1.156E+01,1.193E+01,
     &1.232E+01,1.271E+01,1.311E+01,1.352E+01,1.395E+01,
     &1.438E+01,1.482E+01,1.527E+01,1.573E+01,1.621E+01,
     &1.669E+01,1.718E+01,1.769E+01,1.821E+01,1.874E+01,
     &1.928E+01,1.983E+01,2.040E+01,2.097E+01,2.156E+01,
     &2.217E+01,2.278E+01,2.341E+01,2.406E+01,2.471E+01,
     &2.539E+01,2.607E+01,2.677E+01,2.749E+01,2.822E+01,
     &2.896E+01,2.973E+01,3.050E+01,3.130E+01,3.211E+01,
     &3.293E+01,3.378E+01,3.464E+01,3.552E+01,3.642E+01,
     &3.733E+01,3.827E+01,3.922E+01,4.020E+01,4.119E+01,
     &4.220E+01,4.323E+01,4.429E+01,4.536E+01,4.646E+01,
     &4.758E+01,4.872E+01,4.988E+01,5.106E+01,5.227E+01,
     &5.350E+01,5.476E+01,5.604E+01,5.735E+01,5.868E+01,
     &6.003E+01,6.142E+01,6.282E+01,6.426E+01,6.572E+01,
     &6.721E+01,6.873E+01,7.028E+01,7.186E+01,7.346E+01,
     &7.510E+01,7.677E+01,7.847E+01,8.020E+01,8.196E+01,
     &8.375E+01,8.558E+01,8.744E+01,8.934E+01,9.127E+01,
     &9.324E+01,9.524E+01,9.728E+01,9.936E+01,1.015E+02,
     &1.036E+02,1.058E+02,1.080E+02,1.103E+02,1.126E+02,
     &1.150E+02,1.174E+02,1.198E+02,1.223E+02,1.248E+02,
     &1.274E+02,1.300E+02,1.327E+02,1.354E+02,1.381E+02,
     &1.409E+02,1.438E+02,1.467E+02,1.496E+02,1.526E+02,
     &1.557E+02,1.588E+02,1.619E+02,1.652E+02,1.684E+02,
     &1.718E+02,1.751E+02,1.786E+02,1.821E+02,1.856E+02,
     &1.892E+02,1.929E+02,1.967E+02,2.005E+02,2.043E+02,
     &2.083E+02,2.122E+02,2.163E+02,2.204E+02,2.246E+02,
     &2.289E+02,2.332E+02,2.376E+02,2.421E+02,2.467E+02,
     &2.513E+02,2.560E+02,2.608E+02,2.656E+02,2.706E+02,
     &2.756E+02,2.807E+02,2.859E+02,2.911E+02,2.965E+02,
     &3.019E+02,3.074E+02,3.130E+02,3.187E+02,3.245E+02,
     &3.304E+02,3.364E+02,3.425E+02,3.486E+02,3.549E+02,
     &3.612E+02,3.677E+02,3.743E+02,3.809E+02,3.877E+02/
      
      DATA (FFZ(K,2),K=1,200 ) /
     &1.686E+00,1.738E+00,1.791E+00,1.844E+00,1.899E+00,
     &1.955E+00,2.011E+00,2.069E+00,2.128E+00,2.188E+00,
     &2.249E+00,2.311E+00,2.374E+00,2.438E+00,2.504E+00,
     &2.570E+00,2.638E+00,2.708E+00,2.778E+00,2.850E+00,
     &2.923E+00,2.997E+00,3.072E+00,3.149E+00,3.228E+00,
     &3.307E+00,3.388E+00,3.471E+00,3.555E+00,3.640E+00,
     &3.727E+00,3.815E+00,3.905E+00,3.997E+00,4.090E+00,
     &4.184E+00,4.281E+00,4.378E+00,4.478E+00,4.579E+00,
     &4.682E+00,4.787E+00,4.893E+00,5.002E+00,5.112E+00,
     &5.224E+00,5.337E+00,5.453E+00,5.571E+00,5.690E+00,
     &5.811E+00,5.935E+00,6.060E+00,6.188E+00,6.317E+00,
     &6.449E+00,6.583E+00,6.719E+00,6.857E+00,6.997E+00,
     &7.139E+00,7.284E+00,7.431E+00,7.580E+00,7.732E+00,
     &7.886E+00,8.042E+00,8.201E+00,8.363E+00,8.526E+00,
     &8.693E+00,8.862E+00,9.033E+00,9.207E+00,9.384E+00,
     &9.563E+00,9.746E+00,9.930E+00,1.012E+01,1.031E+01,
     &1.050E+01,1.070E+01,1.090E+01,1.110E+01,1.130E+01,
     &1.151E+01,1.172E+01,1.194E+01,1.215E+01,1.237E+01,
     &1.260E+01,1.283E+01,1.306E+01,1.329E+01,1.353E+01,
     &1.377E+01,1.401E+01,1.426E+01,1.451E+01,1.476E+01,
     &1.502E+01,1.528E+01,1.554E+01,1.581E+01,1.608E+01,
     &1.636E+01,1.664E+01,1.692E+01,1.721E+01,1.750E+01,
     &1.780E+01,1.810E+01,1.840E+01,1.871E+01,1.902E+01,
     &1.934E+01,1.966E+01,1.998E+01,2.031E+01,2.065E+01,
     &2.098E+01,2.133E+01,2.167E+01,2.203E+01,2.238E+01,
     &2.274E+01,2.311E+01,2.348E+01,2.385E+01,2.423E+01,
     &2.462E+01,2.501E+01,2.541E+01,2.581E+01,2.621E+01,
     &2.662E+01,2.704E+01,2.746E+01,2.789E+01,2.832E+01,
     &2.875E+01,2.920E+01,2.965E+01,3.010E+01,3.056E+01,
     &3.103E+01,3.150E+01,3.198E+01,3.246E+01,3.295E+01,
     &3.344E+01,3.395E+01,3.445E+01,3.497E+01,3.549E+01,
     &3.601E+01,3.655E+01,3.709E+01,3.763E+01,3.819E+01,
     &3.875E+01,3.931E+01,3.989E+01,4.047E+01,4.105E+01,
     &4.165E+01,4.225E+01,4.286E+01,4.347E+01,4.410E+01,
     &4.473E+01,4.537E+01,4.601E+01,4.666E+01,4.732E+01,
     &4.799E+01,4.867E+01,4.935E+01,5.005E+01,5.075E+01,
     &5.146E+01,5.217E+01,5.290E+01,5.363E+01,5.437E+01,
     &5.512E+01,5.588E+01,5.665E+01,5.743E+01,5.821E+01,
     &5.901E+01,5.981E+01,6.062E+01,6.145E+01,6.228E+01,
     &6.312E+01,6.397E+01,6.483E+01,6.570E+01,6.658E+01/
      
      DATA (XX(K,1),K=1,200 ) /
     &8.000E-01,6.472E-01,5.944E-01,5.597E-01,5.335E-01,
     &5.121E-01,4.941E-01,4.785E-01,4.647E-01,4.522E-01,
     &4.409E-01,4.306E-01,4.210E-01,4.122E-01,4.039E-01,
     &3.961E-01,3.887E-01,3.817E-01,3.751E-01,3.688E-01,
     &3.628E-01,3.571E-01,3.516E-01,3.463E-01,3.413E-01,
     &3.365E-01,3.318E-01,3.273E-01,3.230E-01,3.188E-01,
     &3.147E-01,3.108E-01,3.070E-01,3.033E-01,2.998E-01,
     &2.963E-01,2.929E-01,2.896E-01,2.864E-01,2.833E-01,
     &2.802E-01,2.773E-01,2.744E-01,2.715E-01,2.688E-01,
     &2.661E-01,2.634E-01,2.608E-01,2.583E-01,2.558E-01,
     &2.534E-01,2.510E-01,2.487E-01,2.464E-01,2.442E-01,
     &2.420E-01,2.398E-01,2.377E-01,2.356E-01,2.336E-01,
     &2.316E-01,2.296E-01,2.277E-01,2.257E-01,2.239E-01,
     &2.220E-01,2.202E-01,2.184E-01,2.167E-01,2.150E-01,
     &2.132E-01,2.116E-01,2.099E-01,2.083E-01,2.067E-01,
     &2.051E-01,2.036E-01,2.020E-01,2.005E-01,1.990E-01,
     &1.976E-01,1.961E-01,1.947E-01,1.933E-01,1.919E-01,
     &1.905E-01,1.891E-01,1.878E-01,1.865E-01,1.852E-01,
     &1.839E-01,1.826E-01,1.814E-01,1.801E-01,1.789E-01,
     &1.777E-01,1.765E-01,1.753E-01,1.741E-01,1.730E-01,
     &1.718E-01,1.707E-01,1.696E-01,1.685E-01,1.674E-01,
     &1.663E-01,1.653E-01,1.642E-01,1.632E-01,1.622E-01,
     &1.611E-01,1.601E-01,1.591E-01,1.581E-01,1.572E-01,
     &1.562E-01,1.552E-01,1.543E-01,1.534E-01,1.524E-01,
     &1.515E-01,1.506E-01,1.497E-01,1.488E-01,1.479E-01,
     &1.471E-01,1.462E-01,1.453E-01,1.445E-01,1.437E-01,
     &1.428E-01,1.420E-01,1.412E-01,1.404E-01,1.396E-01,
     &1.388E-01,1.380E-01,1.372E-01,1.365E-01,1.357E-01,
     &1.349E-01,1.342E-01,1.335E-01,1.327E-01,1.320E-01,
     &1.313E-01,1.306E-01,1.299E-01,1.292E-01,1.284E-01,
     &1.278E-01,1.271E-01,1.264E-01,1.257E-01,1.251E-01,
     &1.244E-01,1.237E-01,1.231E-01,1.224E-01,1.218E-01,
     &1.212E-01,1.205E-01,1.199E-01,1.193E-01,1.187E-01,
     &1.181E-01,1.175E-01,1.169E-01,1.163E-01,1.157E-01,
     &1.151E-01,1.145E-01,1.139E-01,1.134E-01,1.128E-01,
     &1.123E-01,1.117E-01,1.112E-01,1.106E-01,1.101E-01,
     &1.095E-01,1.090E-01,1.085E-01,1.079E-01,1.074E-01,
     &1.069E-01,1.064E-01,1.059E-01,1.054E-01,1.049E-01,
     &1.044E-01,1.039E-01,1.034E-01,1.029E-01,1.024E-01,
     &1.019E-01,1.014E-01,1.010E-01,1.005E-01,1.000E-01/
      
      DATA (XX(K,2),K=1,200 ) /
     &8.000E-01,7.632E-01,7.331E-01,7.073E-01,6.846E-01,
     &6.643E-01,6.458E-01,6.289E-01,6.132E-01,5.986E-01,
     &5.849E-01,5.721E-01,5.600E-01,5.485E-01,5.376E-01,
     &5.272E-01,5.172E-01,5.077E-01,4.986E-01,4.899E-01,
     &4.815E-01,4.734E-01,4.656E-01,4.581E-01,4.508E-01,
     &4.438E-01,4.370E-01,4.304E-01,4.240E-01,4.178E-01,
     &4.118E-01,4.059E-01,4.002E-01,3.947E-01,3.893E-01,
     &3.840E-01,3.789E-01,3.739E-01,3.690E-01,3.643E-01,
     &3.597E-01,3.551E-01,3.507E-01,3.464E-01,3.421E-01,
     &3.380E-01,3.340E-01,3.300E-01,3.261E-01,3.223E-01,
     &3.186E-01,3.150E-01,3.114E-01,3.079E-01,3.045E-01,
     &3.011E-01,2.978E-01,2.945E-01,2.914E-01,2.883E-01,
     &2.852E-01,2.822E-01,2.792E-01,2.763E-01,2.735E-01,
     &2.707E-01,2.679E-01,2.652E-01,2.625E-01,2.599E-01,
     &2.574E-01,2.548E-01,2.523E-01,2.499E-01,2.475E-01,
     &2.451E-01,2.428E-01,2.405E-01,2.382E-01,2.360E-01,
     &2.338E-01,2.316E-01,2.295E-01,2.274E-01,2.254E-01,
     &2.233E-01,2.213E-01,2.193E-01,2.174E-01,2.155E-01,
     &2.136E-01,2.117E-01,2.099E-01,2.081E-01,2.063E-01,
     &2.045E-01,2.028E-01,2.011E-01,1.994E-01,1.977E-01,
     &1.961E-01,1.944E-01,1.929E-01,1.913E-01,1.897E-01,
     &1.882E-01,1.867E-01,1.851E-01,1.837E-01,1.822E-01,
     &1.808E-01,1.793E-01,1.779E-01,1.765E-01,1.752E-01,
     &1.738E-01,1.725E-01,1.711E-01,1.698E-01,1.686E-01,
     &1.673E-01,1.660E-01,1.648E-01,1.635E-01,1.623E-01,
     &1.611E-01,1.599E-01,1.588E-01,1.576E-01,1.564E-01,
     &1.553E-01,1.542E-01,1.531E-01,1.520E-01,1.509E-01,
     &1.498E-01,1.488E-01,1.477E-01,1.467E-01,1.457E-01,
     &1.447E-01,1.437E-01,1.427E-01,1.417E-01,1.407E-01,
     &1.398E-01,1.388E-01,1.379E-01,1.369E-01,1.360E-01,
     &1.351E-01,1.342E-01,1.333E-01,1.324E-01,1.316E-01,
     &1.307E-01,1.299E-01,1.290E-01,1.282E-01,1.273E-01,
     &1.265E-01,1.257E-01,1.249E-01,1.241E-01,1.233E-01,
     &1.225E-01,1.218E-01,1.210E-01,1.203E-01,1.195E-01,
     &1.188E-01,1.180E-01,1.173E-01,1.166E-01,1.159E-01,
     &1.152E-01,1.144E-01,1.138E-01,1.131E-01,1.124E-01,
     &1.117E-01,1.110E-01,1.104E-01,1.097E-01,1.091E-01,
     &1.084E-01,1.078E-01,1.072E-01,1.065E-01,1.059E-01,
     &1.053E-01,1.047E-01,1.041E-01,1.035E-01,1.029E-01,
     &1.023E-01,1.017E-01,1.012E-01,1.006E-01,1.000E-01/
      
      DATA (ZZ(K,1),K=1,200 ) /
     &-2.303E+00,-3.084E+00,-3.649E+00,-4.098E+00,
     &-4.472E+00,-4.795E+00,-5.080E+00,-5.335E+00,
     &-5.568E+00,-5.781E+00,-5.978E+00,-6.161E+00,
     &-6.333E+00,-6.494E+00,-6.647E+00,-6.792E+00,
     &-6.929E+00,-7.060E+00,-7.186E+00,-7.306E+00,
     &-7.421E+00,-7.532E+00,-7.639E+00,-7.742E+00,
     &-7.842E+00,-7.938E+00,-8.031E+00,-8.122E+00,
     &-8.210E+00,-8.295E+00,-8.378E+00,-8.459E+00,
     &-8.538E+00,-8.614E+00,-8.689E+00,-8.762E+00,
     &-8.834E+00,-8.904E+00,-8.972E+00,-9.039E+00,
     &-9.104E+00,-9.168E+00,-9.231E+00,-9.293E+00,
     &-9.353E+00,-9.412E+00,-9.470E+00,-9.528E+00,
     &-9.584E+00,-9.639E+00,-9.693E+00,-9.746E+00,
     &-9.799E+00,-9.851E+00,-9.901E+00,-9.951E+00,
     &-1.000E+01,-1.005E+01,-1.010E+01,-1.014E+01,
     &-1.019E+01,-1.024E+01,-1.028E+01,-1.033E+01,
     &-1.037E+01,-1.041E+01,-1.046E+01,-1.050E+01,
     &-1.054E+01,-1.058E+01,-1.062E+01,-1.066E+01,
     &-1.070E+01,-1.074E+01,-1.078E+01,-1.082E+01,
     &-1.086E+01,-1.089E+01,-1.093E+01,-1.097E+01,
     &-1.101E+01,-1.104E+01,-1.108E+01,-1.111E+01,
     &-1.115E+01,-1.118E+01,-1.122E+01,-1.125E+01,
     &-1.128E+01,-1.132E+01,-1.135E+01,-1.138E+01,
     &-1.141E+01,-1.145E+01,-1.148E+01,-1.151E+01,
     &-1.154E+01,-1.157E+01,-1.160E+01,-1.163E+01,
     &-1.166E+01,-1.169E+01,-1.172E+01,-1.175E+01,
     &-1.178E+01,-1.181E+01,-1.184E+01,-1.186E+01,
     &-1.189E+01,-1.192E+01,-1.195E+01,-1.198E+01,
     &-1.200E+01,-1.203E+01,-1.206E+01,-1.208E+01,
     &-1.211E+01,-1.214E+01,-1.216E+01,-1.219E+01,
     &-1.221E+01,-1.224E+01,-1.226E+01,-1.229E+01,
     &-1.231E+01,-1.234E+01,-1.236E+01,-1.239E+01,
     &-1.241E+01,-1.244E+01,-1.246E+01,-1.248E+01,
     &-1.251E+01,-1.253E+01,-1.255E+01,-1.258E+01,
     &-1.260E+01,-1.262E+01,-1.264E+01,-1.267E+01,
     &-1.269E+01,-1.271E+01,-1.273E+01,-1.276E+01,
     &-1.278E+01,-1.280E+01,-1.282E+01,-1.284E+01,
     &-1.286E+01,-1.289E+01,-1.291E+01,-1.293E+01,
     &-1.295E+01,-1.297E+01,-1.299E+01,-1.301E+01,
     &-1.303E+01,-1.305E+01,-1.307E+01,-1.309E+01,
     &-1.311E+01,-1.313E+01,-1.315E+01,-1.317E+01,
     &-1.319E+01,-1.321E+01,-1.323E+01,-1.325E+01,
     &-1.327E+01,-1.329E+01,-1.330E+01,-1.332E+01,
     &-1.334E+01,-1.336E+01,-1.338E+01,-1.340E+01,
     &-1.342E+01,-1.343E+01,-1.345E+01,-1.347E+01,
     &-1.349E+01,-1.351E+01,-1.352E+01,-1.354E+01,
     &-1.356E+01,-1.358E+01,-1.360E+01,-1.361E+01,
     &-1.363E+01,-1.365E+01,-1.366E+01,-1.368E+01,
     &-1.370E+01,-1.372E+01,-1.373E+01,-1.375E+01,
     &-1.377E+01,-1.378E+01,-1.380E+01,-1.382E+01/
      
      DATA (ZZ(K,2),K=1,200 ) /
     &-2.303E+00,-2.512E+00,-2.700E+00,-2.871E+00,
     &-3.029E+00,-3.175E+00,-3.310E+00,-3.438E+00,
     &-3.557E+00,-3.670E+00,-3.778E+00,-3.880E+00,
     &-3.977E+00,-4.070E+00,-4.159E+00,-4.245E+00,
     &-4.328E+00,-4.407E+00,-4.484E+00,-4.558E+00,
     &-4.630E+00,-4.699E+00,-4.767E+00,-4.832E+00,
     &-4.896E+00,-4.958E+00,-5.019E+00,-5.078E+00,
     &-5.135E+00,-5.191E+00,-5.246E+00,-5.300E+00,
     &-5.352E+00,-5.403E+00,-5.453E+00,-5.503E+00,
     &-5.551E+00,-5.598E+00,-5.645E+00,-5.690E+00,
     &-5.735E+00,-5.779E+00,-5.822E+00,-5.864E+00,
     &-5.906E+00,-5.947E+00,-5.988E+00,-6.027E+00,
     &-6.067E+00,-6.105E+00,-6.143E+00,-6.181E+00,
     &-6.217E+00,-6.254E+00,-6.290E+00,-6.325E+00,
     &-6.360E+00,-6.394E+00,-6.428E+00,-6.462E+00,
     &-6.495E+00,-6.528E+00,-6.560E+00,-6.592E+00,
     &-6.624E+00,-6.655E+00,-6.686E+00,-6.716E+00,
     &-6.746E+00,-6.776E+00,-6.805E+00,-6.835E+00,
     &-6.863E+00,-6.892E+00,-6.920E+00,-6.948E+00,
     &-6.976E+00,-7.003E+00,-7.030E+00,-7.057E+00,
     &-7.084E+00,-7.110E+00,-7.136E+00,-7.162E+00,
     &-7.188E+00,-7.213E+00,-7.238E+00,-7.263E+00,
     &-7.288E+00,-7.312E+00,-7.336E+00,-7.360E+00,
     &-7.384E+00,-7.408E+00,-7.431E+00,-7.455E+00,
     &-7.478E+00,-7.501E+00,-7.523E+00,-7.546E+00,
     &-7.568E+00,-7.590E+00,-7.612E+00,-7.634E+00,
     &-7.656E+00,-7.677E+00,-7.698E+00,-7.720E+00,
     &-7.741E+00,-7.761E+00,-7.782E+00,-7.803E+00,
     &-7.823E+00,-7.843E+00,-7.863E+00,-7.883E+00,
     &-7.903E+00,-7.923E+00,-7.943E+00,-7.962E+00,
     &-7.981E+00,-8.001E+00,-8.020E+00,-8.039E+00,
     &-8.057E+00,-8.076E+00,-8.095E+00,-8.113E+00,
     &-8.132E+00,-8.150E+00,-8.168E+00,-8.186E+00,
     &-8.204E+00,-8.222E+00,-8.239E+00,-8.257E+00,
     &-8.274E+00,-8.292E+00,-8.309E+00,-8.326E+00,
     &-8.343E+00,-8.360E+00,-8.377E+00,-8.394E+00,
     &-8.411E+00,-8.427E+00,-8.444E+00,-8.460E+00,
     &-8.476E+00,-8.493E+00,-8.509E+00,-8.525E+00,
     &-8.541E+00,-8.557E+00,-8.572E+00,-8.588E+00,
     &-8.604E+00,-8.619E+00,-8.635E+00,-8.650E+00,
     &-8.666E+00,-8.681E+00,-8.696E+00,-8.711E+00,
     &-8.726E+00,-8.741E+00,-8.756E+00,-8.771E+00,
     &-8.786E+00,-8.800E+00,-8.815E+00,-8.829E+00,
     &-8.844E+00,-8.858E+00,-8.872E+00,-8.887E+00,
     &-8.901E+00,-8.915E+00,-8.929E+00,-8.943E+00,
     &-8.957E+00,-8.971E+00,-8.985E+00,-8.998E+00,
     &-9.012E+00,-9.026E+00,-9.039E+00,-9.053E+00,
     &-9.066E+00,-9.080E+00,-9.093E+00,-9.106E+00,
     &-9.119E+00,-9.133E+00,-9.146E+00,-9.159E+00,
     &-9.172E+00,-9.185E+00,-9.197E+00,-9.210E+00/
      END

      BLOCK DATA PARAM_INI
C-----------------------------------------------------------------------
C....This block data contains default values
C.   of the parameters used in fragmentation
C-----------------------------------------------------------------------
      SAVE

      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      COMMON /S_CZDIS/ FA, FB0
      COMMON /S_CZDISs/ FAs1, fAs2
      COMMON /S_CZDISc/ ZDMAX, EPSI
      COMMON /S_CZLEAD/ CLEAD, FLEAD
      COMMON /S_CPSPL/ CCHIK(3,39)
      COMMON /S_CQDIS/ PPT0 (35),ptflag
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)
      COMMON /S_CUTOFF/ STR_mass_val, STR_mass_val_hyp, STR_mass_sea
      COMMON /CKFRAG/ KODFRAG
      COMMON /S_DIFMAss/ XM2MIN(3),ALXMIN(3),SLOP0,ASLOP,BSLOP,XMASS(2)

C...  default output unit
      DATA Lun /6/

c...new fragmentation for charmed particles
      DATA EPSI /2.0/

C...mass cutoff for soft strings
      data STR_mass_val /.35/ 
      data STR_mass_val_hyp /.4/ 
      data STR_mass_sea /1./ 
C...Longitudinal Fragmentation function
      DATA FA /0.5/, FB0 /0.8/
C...Longitudinal Fragmentation function for leading baryons
      DATA CLEAD  /0.6/, FLEAD  /0.6/
c     strange fragmentation
      data FAs1 /3./, fAs2 /3./	
C...pT of sea partons
      DATA PTFLAG /1./
      DATA PPT0 /0.30,0.30,0.450,30*0.60,0.3,0.9/
C...Splitting parameters
      DATA CCHIK /15*0,21*2.,6*3.,57*0,18*3./
C...Parameters of flavor formation 
c     last in use: 74
      DATA PAR /0.04,0.3,0.3,0.14,0.3,0.3,0.15,0.,7.0,0., ! 10
     &     3*0.,0.04,0.04,0.04,0.04, 0.5,0.8, 0.5,        ! 20
     &     0.8,6.,0.5, 0.004,5*0.,0.7,                    ! 30
     &     2*0.,0.1,0.,3.,0.35,0.,0.5,2*0,                ! 40
     &     1.,2.,0.,.99,0.,0.3,.45,.6,.6,0.6,             ! 50
     &     1.,0.,6.,0.2,4*0.,1.1,0.8,                     ! 60
     &     0.33,3.,1.,0.25,0.3,0.3,0.6,0.007,0.03,0.007,  ! 70
     &     1.,0.3,0.,0.3,6*0./                            ! 80
c     last in use: 43
      DATA IPAR /9*0,1,0,2,8*0,20*0,40*0/

C...Fragmentation of nuclei
      DATA KODFRAG /0/
C...Debug label and event counter
      DATA Ndebug /0/
      DATA Ncall /0/
C...Diffractive mass parameters
      DATA XM2MIN /1.5, 0.2, 0.6/                  ! M_x**2(min) GeV**2
      DATA ALXMIN /0.405465,-1.6094379,-0.5108256/ ! log[M_x**2(min)]
      DATA SLOP0 /6.5/                 ! b (slope_ for Mx**2 > 5 GeV**2
      DATA ASLOP /31.10362/            ! fit to the slope parameter.
      DATA BSLOP /-15.29012/


      END



      SUBROUTINE PARAM_PRINT(LUN)
      SAVE

      COMMON /S_CZDIS/ FA, FB0
      COMMON /S_CZLEAD/ CLEAD, FLEAD
      COMMON /S_CPSPL/ CCHIK(3,39)
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, kb ,kt
      COMMON /S_DEBUG/ Ncall, Ndebug, Lunn
      COMMON /S_CQDIS/ PPT0 (35),ptflag
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)

      WRITE (LUN, 25)
25      FORMAT( //,1x,40('-'), /
     +   ' SIBYLL MONTE CARLO PROGRAM. Version 2.2f',/,1x,40('-'),/
     +   ' List of parameters: ' )

      WRITE (LUN, 31) FA, FB0
31      FORMAT (' Parameters of longitudinal fragmentation: ', /,
     +          '  f(z) = (1-z)**a * exp(-b * mt**2/z) ', /,
     +          '  a = ', f9.3, 3x, ' b = ', f9.3, ' GeV**-2' )
      WRITE (LUN, 32) CLEAD, 1./FLEAD-1.
32      FORMAT (' Parameters of leading fragmentation: ', /,
     +   '  f(z) = c + (1-z)**a ', /,
     +   '  c = ',f9.3,3x,' a = ',f9.3) 

      WRITE (LUN, 35) PPT0(1), PPT0(3), PPT0(11),ppt0(10)
35      FORMAT (' <pT> of sea partons ', /,
     +   2x,'<pT>(u/d) ',F8.3,2x,'<pT>(s) ',f8.3,2x,'<pT>(qq) ',f8.3,
     +     2x,'<pT>(val) ',f8.3)

      WRITE (LUN, 120) (PAR(K),K=1,16)
120      FORMAT (1x, 'Parameters of flavor formation: ',/,
     +   3x,'PAR(1) = Prob(qq)/Prob(q) =              ',F10.2,/,
     +   3x,'PAR(2) = Prob(s)/Prob(u)  =              ',F10.2,/,
     +   3x,'PAR(3) = Prob(us)/Prob(ud) =             ',F10.2,/,
     +   3x,'PAR(4) = Prob(ud_0)/Prob(ud_1) =         ',F10.2,/,
     +   3x,'PAR(5) = Prob(Vector)/Prob(Scalar) =     ',F10.2,/,
     +   3x,'PAR(6) = Prob(K*)/Prob(K) =              ',F10.2,/,
     +   3x,'PAR(7) = Prob(spin 3/2)/Prob(spin=1/2) = ',F10.2,/,
     +   3x,'PAR(8) = Prob(B-M-Bbar)/Prob(B-Bbar) =   ',F10.2,/,
     +   3x,'PAR(9) = Phase space suppression of MI = ',F10.2,/,
     +   3x,'PAR(10)= Low-energy limit for pt cutoff= ',F10.2,/,
     +   3x,'PAR(11)= Pt cutoff factor for exp      = ',F10.2,/,
     +   3x,'PAR(12)= Pt cutoff factor for exp      = ',F10.2,/,
     +   3x,'PAR(13)= max. mass in diffraction      = ',F10.2,/,
     +   3x,'PAR(14)= Prob(qq)/Prob(q) std. value   = ',F10.2,/,
     +   3x,'PAR(15)= Prob(qq)/Prob(q) in hard jets = ',F10.2,/,
     +   3x,'PAR(16)= Prob(qq)/Prob(q) in diff.     = ',F10.2)

      WRITE (LUN, 130) (IPAR(K),K=1,16)
130      FORMAT (1x, 'Model switches: ',/,
     +   3x,'IPAR(1) = not used =                     ',I4,/,
     +   3x,'IPAR(2) = not used =                     ',I4,/,
     +   3x,'IPAR(3) = exponential pt =               ',I4,/,
     +   3x,'IPAR(4) = decouple qq/q in val. strings= ',I4,/,
     +   3x,'IPAR(5) = decouple qq/q in hm. diff. =   ',I4,/,
     +   3x,'IPAR(6) = decouple qq/q in hard strings= ',I4,/,
     +   3x,'IPAR(7) = remnant (not implemented yet)= ',I4,/,
     +   3x,'IPAR(8) = jet kinematic pdf set (DO/GRV)=',I4,/,
     +   3x,'IPAR(9) = smear lowest diff. mass =      ',I4,/,
     +   3x,'IPAR(10)= high mass diff. mode (d:ON)=   ',I4,/,
     +   3x,'IPAR(11)= leading vec. meson prod. model=',I4,/,
     +   3x,'IPAR(12)= inel. screening in pAir=       ',I4,/,
     +   3x,'IPAR(13)= decouple qq/q in val. strings= ',I4,/,
     +   3x,'IPAR(14)= fireball model =               ',I4,/,
     +   3x,'IPAR(15)= intrinsic charm fix =          ',I4,/,
     +   3x,'IPAR(16) = not used =                    ',I4)

      WRITE (LUN, 40)
      WRITE (LUN, 41) CCHIK (1,13), CCHIK(2,13)
40      FORMAT(' Parameters of hadron splitting ' )
41      FORMAT('   p -> [(ud) u] splitting: alpha = ', F10.3, /,
     +         '   p -> [(uu) d] splitting: alpha = ', F10.3 )

      END



C-----------------------------------------------------------------------
C  Code for diffraction
C-----------------------------------------------------------------------

c EJA: change in more detail now

      SUBROUTINE SIB_DIFF (L0, JDIF1, Ecm, Irec, IREJ)
C-----------------------------------------------------------------------
C...diffraction dissociation
C.  INPUT L0 = index of "beam particle"
C.             the target is assumed to be a proton.
C.    JDIF1 = 1  "beam diffraction"
C.          = 2  "target diffraction"
C.          = 3  "double diffraction"
C     Irec  flag to avoid recursive calls of SIB_DIFF and SIB_NDIFF
C-----------------------------------------------------------------------
      SAVE

      COMMON /S_PLIST/ P(8000,5), LLIST(8000), NP
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, kb ,kt
      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      COMMON /S_MASS1/ AM(99), AM2(99)
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)
      DIMENSION P0(5)
      DIMENSION KK(39)
      DATA PI /3.1415926/
      DATA KK /5*0,7*2,2*1,19*0,6*1/
      COMMON /S_PARTO/ NFORIG(8000), NPORIG(8000), IPFLAG, NINT
      COMMON /S_DIFMAss/ XM2MIN(3),ALXMIN(3),SLOP0,ASLOP,BSLOP,XMASS(2)
c      DIMENSION XM2MIN(3), ALXMIN(3)
c      DATA XM2MIN /1.5, 0.2, 0.6/                  ! M_x**2(min) GeV**2
c      DATA ALXMIN /0.405465,-1.6094379,-0.5108256/ ! log[M_x**2(min)]
c      DATA SLOP0 /6.5/                 ! b (slope_ for Mx**2 > 5 GeV**2
c      DATA ASLOP /31.10362/            ! fit to the slope parameter.
c      DATA BSLOP /-15.29012/

      if(Ndebug.gt.1) 
     &  print *,' SIB_DIFF: called with (L0,JDIF1,Ecm):',
     &  L0,JDIF1,Ecm

      if(Irec.eq.1) Ipflag= -1

 20   IREJ = 1
      LA = IABS(L0)
      XM2MAX = PAR(13)*Ecm*Ecm

C...Double diffraction
      IF (JDIF1 .EQ. 3)   THEN
         K = KK(LA)
         AL = LOG(XM2MAX/XM2MIN(K))
         ALX = ALXMIN(K) + AL*S_RNDM(0)
         XMB2 = EXP(ALX)
         XMB = SQRT (XMB2)
         AL = LOG(XM2MAX/XM2MIN(1))
         ALX = ALXMIN(1) + AL*S_RNDM(0)
         XMT2 = EXP(ALX)
         XMT = SQRT (XMT2)
         X1 = 1.+(XMB2-XMT2)/(Ecm*Ecm)
         X2 = 2.-X1
         SLOPE = MAX(SLOP0, ASLOP+BSLOP*ALX)
50       T = -LOG(S_RNDM(0))/SLOPE
         PT = SQRT(T)
         PZ1 = 0.25*Ecm*Ecm*X1*X1-XMB2-PT*PT
         PZ2 = 0.25*Ecm*Ecm*X2*X2-XMT2-PT*PT
         IF (PZ1.LT.0. .OR. PZ2.LT.0.)   GOTO 50
         PHI = PI*S_RNDM(0)
         P0(5) = XMB
         P0(4) = 0.5*Ecm*X1
         P0(1) = PT*COS(PHI)
         P0(2) = PT*SIN(PHI)
         P0(3) = SQRT(PZ1)
         XMASS(1) = XMB
         CALL DIFDEC (L0, Irec, IDBAD, P0)
         IF(IDBAD.eq.1)goto 20
         P0(5) = XMT
         P0(4) = 0.5*Ecm*X2
         P0(1) = -P0(1)
         P0(2) = -P0(2)
         P0(3) = -SQRT(PZ2)
         Ipflag= -2
         XMASS(2) = XMT
         CALL DIFDEC (13, Irec, IDBAD, P0)
         IF(IDBAD.eq.1)goto 20
         IREJ = 0
         RETURN
      ENDIF

C...Single diffraction
      IF (JDIF1.EQ. 1)  THEN
         K = KK(LA)
         EM  = AM(13)
         EM2 = AM2(13)
         L = 13
         ZD = -1.
      ELSE
         K = 1
         EM  = AM(LA)
         EM2 = AM2(LA)
         L = L0
         L0 = 13
         ZD = +1.
      ENDIF
C...Generate the mass of the diffracted system Mx (1/Mx**2 distribution)
      AL = LOG(XM2MAX/XM2MIN(K))
      ALX = ALXMIN(K) + AL*S_RNDM(0)
      XM2 = EXP(ALX)

c... added part
      X = XM2/XM2MAX*PAR(13)
      IF (X.GT.PAR(13)-0.05) THEN
        PRO = 0.5*(1.+(X-PAR(13))/0.05)
        IF (S_RNDM(0).LT.PRO) X = 2.*PAR(13)-X
        XM2 = XM2MAX*X/PAR(13)
      ENDIF
c...

      XM = SQRT (XM2)
      XMB = XM
      XMT = XM
      XMASS(1) = XMB
      XMASS(2) = XMT
C...Generate the Kinematics of the pseudoelastic hadron
      X = 1.-(XM2-EM2)/(Ecm*Ecm)
      NP = NP+1
      P(NP,4) = 0.5*Ecm*X
      SLOPE = MAX(SLOP0, ASLOP+BSLOP*ALX)
60    T = -LOG(MAX(1.E-10,S_RNDM(0)))/SLOPE
      PT = SQRT(T*X)
      PZ2 = P(NP,4)**2-EM2 - PT*PT
      IF (PZ2 .LT.0.)   GOTO 60
      PHI = PI*S_RNDM(0)
      P(NP,3) = SQRT(PZ2)*ZD
      P(NP,1) = PT*COS(PHI)
      P(NP,2) = PT*SIN(PHI)
      P(NP,5) = EM
      LLIST(NP) = L
      NPORIG(NP) = IPFLAG
C...Generating the hadronic system recoling against the produced particle
      P0(5) = SQRT(XM2)
      P0(4) = 0.5*Ecm*(2.-X)
      DO J=1,3
         P0(J) = -P(NP,J)
      ENDDO
      CALL DIFDEC (L0, Irec, IDBAD, P0)
      IF(IDBAD.eq.1)goto 20

      IREJ = 0

      END


      SUBROUTINE DIFDEC (L0, Irec, IBAD, P0)
C-----------------------------------------------------------------------
C..."decay" of an excited state with the quantum numbers
C.   of particle L0 and the 5-momentum P0
C.   - low energy: phase space decay (fire ball model)
C.   - intermediate energy: one-string decay (longitudinal phase space)
C.   - high energy: pomeron-hadron scattering (multi-string model) 
C-----------------------------------------------------------------------
      SAVE

      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      COMMON /S_PLIST/ P(8000,5), LLIST(8000), NP
      COMMON /S_MASS1/ AM(99), AM2(99)
      COMMON /S_CHP/ ICHP(99), ISTR(99), IBAR(99)
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)
      COMMON /S_PARTO/ NFORIG(8000), NPORIG(8000), IPFLAG, NINT
      DIMENSION P0(5), LL(10), PD(10,5), BE(3), LCON(6:39), LRES(6:39),
     +     LRES1(6:39)
      DATA EMIN /0.7/
      DATA EMIN2 /10./
      DATA LCON /7,6,6,11,11,9,9,14,13,19*0,35,34,35,38,37,39/
      DATA LRES /26,27,27,11,11,9,9,14,13,19*0,35,34,35,38,37,39/ 
      DATA LRES1 /27,25,26,11,11,9,9,14,13,19*0,35,34,35,38,37,39/
      DATA PCHEX /0.33/            ! probability of charge exchange
      DATA PRES /0.7/            ! probability of forming a resonance
      CHARACTER CWRN*(25)

      LA = IABS(L0)
      DELTAE = P0(5) - AM(LA)
      EMIN = PAR(30)
      PAR1def= PAR(1)
      if(Irec.gt.0) PAR(1)= PAR(16)
C...pomeron-hadron scattering (pi0 is used instead of pomeron)
      IF ((IPAR(10).gt.0).and.(Irec.gt.0).and.(DELTAE.gt.EMIN2))  THEN
         if(Ndebug.gt.2) 
     &   print *,'  DIFDEC: central (DELTAE):', DELTAE
         N1 = NP+1
         if(irec.gt.0.and.ipar(5).eq.1) par(1)= par(15)
 50      CONTINUE
         IPFLAG= IPFLAG*100
         CALL SIB_NDIFF(L0, 1, P0(5), 0, IREJ) ! ori
         IF(IREJ.NE.0) THEN
            NP = N1-1
            GOTO 50
         ENDIF
         PAR(1) = PAR1def
         DO J=1,3
            BE(J)=P0(J)/P0(4)
         ENDDO
         GA=P0(4)/P0(5)
         if(P0(3).lt.0.) then
           do i=N1,NP
             P(I,3) = -P(I,3)
           enddo
         endif
         DO I=N1,NP
            BEP=BE(1)*P(I,1)+BE(2)*P(I,2)+BE(3)*P(I,3)
            DO J=1,3
               P(I,J)=P(I,J)+GA*(GA*BEP/(1.+GA)+P(I,4))*BE(J)
            ENDDO
            P(I,4)=GA*(P(I,4)+BEP)
         ENDDO

C..."string-like" decay
      ELSE IF (DELTAE .GT. EMIN)  THEN            
 10      IF(NDEBUG.gt.2) 
     &   WRITE(LUN,*)'  DIFDEC: string-like, (DELTAE):',DELTAE
           N1 = NP+1
         CALL HSPLI(L0,IFL1,IFL2)
         IF(IFL1.eq.0.and.IFL2.eq.0)THEN
            WRITE(LUN,*) 'DIFDEC: HSPLI rejection in string-like!'
            WRITE(LUN,*) ' RUN (SQS,KB,KT):',SQS,KB,KT
            STOP
         ENDIF
         IF (P0(3) .GT. 0..and.L0.gt.0)  THEN
            IFLA = IFL2
            IFL2 = IFL1
            IFL1 = IFLA
         ENDIF
         PAR24_def = PAR(24)
         SELECT CASE(IPAR(15))
         CASE(2,4,5,6,8,9,10,11)
            PAR(24) = PAR(25)*EXP(-PAR(26)/P0(5))
         CASE(7)
            PAR(24) = PAR(25)
         END SELECT
         IPFLAG = IPFLAG*10
         IF(IPAR(17).eq.1)THEN
            CALL STRING_FRAG_4FLV 
     +           (P0(5), IFL1, IFL2, 0.,0.,0.,0.,IBAD,-1)
         ELSE
c            CALL STRING_FRAG (P0(5), IFL1, IFL2, 0.,0.,0.,0.,IFBAD,-1)
            CWRN = '3flv routines removed    '
            CALL SIB_ABRT(LUN,CWRN)
         ENDIF
         PAR(24) = PAR24_def
         IF (IBAD .EQ. 1)then
            if(ndebug.gt.1) 
     &           WRITE(LUN,*)' SIB_DIFF: string-frag rejection!'
            NTRYs = NTRYs + 1
            NP = N1-1
            IBAD = 0
            IF(NTRYs.gt.5)then ! resample diff. mass              
               PAR(1) = PAR1def
               NP = 0
               IBAD = 1
               RETURN
            endif
            GOTO 10
         ENDIF
         DO J=1,3
            BE(J)=P0(J)/P0(4)
         ENDDO
         GA=P0(4)/P0(5)
         DO I=N1,NP
            BEP=BE(1)*P(I,1)+BE(2)*P(I,2)+BE(3)*P(I,3)
            DO J=1,3
               P(I,J)=P(I,J)+GA*(GA*BEP/(1.+GA)+P(I,4))*BE(J)
            ENDDO
            P(I,4)=GA*(P(I,4)+BEP)
         ENDDO

C...Phase space decay of the excited state
      ELSEIF(DELTAE.GT.AM(7)+0.02)THEN
         if(Ndebug.gt.2) 
     &     WRITE(LUN,*)'  DIFDEC: fireball, (DELTAE):',DELTAE
         IF(IPAR(14).GT.0.and.IPAR(14).NE.7)THEN
            IF(IPAR(14).eq.5) PCHEX = 0.
            NPI=0
            IRES = 0
            IF (S_RNDM(0).LT.PRES) THEN
               IF (LA.LT.9) THEN
c     if kinematically possible produce rho0 in charge exchange
                  LL(1) = LRES(LA)
                  DELTAE = P0(5) -  AM(LRES(LA))
                  IF (DELTAE.GT.AM(7)+0.02) GOTO 100
               ENDIF
            ENDIF
c     switch charge exchange on/off
            IF( S_RNDM(0).LT.PCHEX)THEN
               LL(1) = LCON(LA)*ISIGN(1,L0)
               IF( (L0 .EQ. 6) .OR. (L0 .EQ. 11) )
     .              LL(1) = LL(1)+INT(1.99999*S_RNDM(0))
            ELSE
               LL(1) = L0
            ENDIF
            
            DELTAE = P0(5) -  AM(LA)
 100        AV = 2.*SQRT(DELTAE)
            LA1 = IABS(LL(1))
            NPI = AV*(2.+0.5*GASDEV(LA))
            IF (IPAR(14).EQ.6)THEN
               IF(NPI.LT.1.OR.NPI.GT.9.OR.AM(LA1)+NPI*AM(7)+0.02
     .              .GT.P0(5))  GOTO 100
            ELSE
               IF(NPI.LT.0.OR.NPI.GT.9.OR.AM(LA1)+NPI*AM(7)+0.02
     .              .GT.P0(5))  GOTO 100
            ENDIF
c     create resonances inside fireball..
            IF(IPAR(14).ge.2
     +           .and.DELTAE.GE.AM(LA1)+AM(27)+(NPI-1)*AM(7)+0.02)
     +           IRES = 1
            IF(IPAR(14).ge.3.and.DELTAE.GE.AM(LA1)+NPI*AM(27)+0.02) 
     +           IRES=3
            JQQ = ICHP(LA)*ISIGN(1,L0)-
     .           ICHP(IABS(LL(1)))*ISIGN(1,LL(1))  
 120        JQTOT = 0.
            DO K=2,NPI
               LL(K) = 6+INT(S_RNDM(0)*2.99999)
c     suppress pi0 in fireball
               IF(IPAR(14).ge.4)
     +              LL(K) = 7+INT(S_RNDM(0)*1.99999)
c     IF(IRES.EQ.1.and.S_RNDM(LA).LT.0.5)
               IF(IRES.EQ.1) THEN
                  LL(K) = 27-INT(S_RNDM(0)*2.99999)
                  IRES = 2
               ENDIF
               IF(IRES.EQ.3)
     +              LL(K) = 27-INT(S_RNDM(0)*2.99999)
               JQTOT = JQTOT + ICHP(LL(K))
            ENDDO
            JQR = JQQ-JQTOT
            IF (JQR.LT.-1.OR.JQR.GT.1)  GOTO 120
            LL(NPI+1) = 6+JQR
            IF (LL(NPI+1) .EQ. 5)  LL(NPI+1)=8
            CALL DECPAR (0,P0,NPI+1,LL, PD)
            DO J=1,NPI+1
               NP = NP+1
               LLIST(NP) = LL(J)
               nporig(NP)= Ipflag*2
               DO K=1,5
                  P(NP,K) = PD(J,K)
               ENDDO
            ENDDO

         ELSEIF (IPAR(14).EQ.7.AND.LA.LT.9) THEN
c     all diff states go to resonances for pi beam ..
            NPI=0
            IRES = 0
            LL(1) = LRES1(LA)
            DELTAE = P0(5) -  AM(LL(1))
            IF( DELTAE.LT.AM(7)+0.02) GOTO 222
            IF( S_RNDM(0).LT.PAR(31))THEN
               LL(1) = LRES1(LCON(LA))
               IF( (L0 .EQ. 6) .OR. (L0 .EQ. 11) )
     .              LL(1) = LRES1(IABS(L0)+INT(1.99999*S_RNDM(0)))
            ENDIF
 300        AV = 2.*SQRT(DELTAE)
            LA1 = IABS(LL(1))
            NPI = AV*(2.+0.5*GASDEV(LA))
            IF(PAR(32).ne.0.0) NPI = AV*(PAR(32)+0.5*GASDEV(LA))
            IF(NPI.LT.0.OR.NPI.GT.9.OR.AM(LA1)+NPI*AM(7)+0.02
     .           .GT.P0(5))  GOTO 300
c     create resonances inside fireball..
c            IRES=3
            JQQ = ICHP(LA)*ISIGN(1,L0)-
     .           ICHP(IABS(LL(1)))*ISIGN(1,LL(1))  
 320        JQTOT = 0.
            DO K=2,NPI
               LL(K) = 6+INT(S_RNDM(0)*2.99999)
c     suppress pi0 in fireball
c               IF(IPAR(14).ge.4)
c     +              LL(K) = 7+INT(S_RNDM(0)*1.99999)
c     IF(IRES.EQ.1.and.S_RNDM(LA).LT.0.5)
c               LL(K) = 27-INT(S_RNDM(0)*2.99999)
               JQTOT = JQTOT + ICHP(LL(K))
            ENDDO
            JQR = JQQ-JQTOT
            IF (JQR.LT.-1.OR.JQR.GT.1)  GOTO 320
            LL(NPI+1) = 6+JQR
            IF (LL(NPI+1) .EQ. 5)  LL(NPI+1)=8
            CALL DECPAR (0,P0,NPI+1,LL, PD)
            DO J=1,NPI+1
               NP = NP+1
               LLIST(NP) = LL(J)
               nporig(NP)= Ipflag*2
               DO K=1,5
                  P(NP,K) = PD(J,K)
               ENDDO
            ENDDO

         ELSEIF (IPAR(14).LE.-1) THEN
C...  generalized fireball model
            IF(Ndebug.gt.2) 
     &           PRINT *,'  DIFDEC: using generalized fireball!'
            CALL FIREBALL_4FLV(L0,P0,IFBAD)
            IF(IFBAD.eq.1)THEN
               IF(ndebug.ne.0)THEN
                  IF(NRJECT.le.10)THEN
                  WRITE(LUN,*)' DIFDEC: warning: fireball rejection! ',
     &                    'diff. mass to low to dissociate beam!'
                  WRITE(LUN,*)' DIFDEC: m_Beam, DELTAE ,AM(7)+0.02 : ', 
     &                    AM(LA),DELTAE,'>',AM(7)+0.02
                  ENDIF
                  IF(NRJECT.eq.10) 
     &            WRITE(LUN,*)' this was the last warning.. good luck!'
               ENDIF
               NRJECT = NRJECT + 1
               NP = 0
               IBAD = 1
               RETURN
            ENDIF

         ELSE
 222        IF(IPAR(14).EQ.7)  DELTAE = P0(5) - AM(LA)
            AV = 2.*SQRT(DELTAE)
 200        NPI = AV*(1.+0.5*GASDEV(0))
c            print *,'npi:',npi,'av',av,'p05',p0(5),am(la),deltae
            IF(NPI.LE.0.OR.NPI.GT.9.OR.AM(LA)+NPI*AM(7)+0.02
     .           .GT.P0(5))  GOTO 200
            IF (S_RNDM(0).LT.PCHEX)  THEN
               LL(NPI+1) = LCON(LA)*ISIGN(1,L0)
               IF( (L0 .EQ. 6) .OR. (L0 .EQ. 11) )
     .              LL(NPI+1) = LL(NPI+1)+INT(1.99999*S_RNDM(0))
            ELSE
               LL(NPI+1) = L0
            ENDIF
            JQQ = ICHP(LA)*ISIGN(1,L0)-
     .           ICHP(IABS(LL(NPI+1)))*ISIGN(1,LL(NPI+1))  
 220        JQTOT = 0.
            DO K=1,NPI-1
               LL(K) = 6+INT(S_RNDM(0)*2.99999)
               JQTOT = JQTOT + ICHP(LL(K))
            ENDDO
            JQR = JQQ-JQTOT
            IF (JQR.LT.-1.OR.JQR.GT.1)  GOTO 220
            LL(NPI) = 6+JQR
            IF (LL(NPI) .EQ. 5)  LL(NPI)=8
            CALL DECPAR (0,P0,NPI+1,LL, PD)
            DO J=1,NPI+1
               NP = NP+1
               LLIST(NP) = LL(J)
               NPORIG(NP) = IPFLAG*2
               DO K=1,5
                  P(NP,K) = PD(J,K)
               ENDDO
            ENDDO
         ENDIF
      ELSE
         NTRY=NTRY+1
         IF(ndebug.gt.0)THEN
         IF(NTRY.le.20)THEN
            WRITE(LUN,*) '  DIFDEC rejection! ',
     &           'diff. mass to low to dissociate beam!'
            WRITE(LUN,*) '  DIFDEC: m_Beam, DELTAE : ', AM(LA),DELTAE
            IF(NTRY.eq.20) WRITE(LUN,*) '  DIFDEC: last warning...'
         ENDIF
         ENDIF
         NP = 0
         IBAD = 1
         PAR(1) = PAR1def
         RETURN
      ENDIF
      PAR(1) = PAR1def
      END

      SUBROUTINE FIREBALL_4FLV(L0,P0,IREJ)
C-----------------------------------------------------------------------
C... "decay" of an excited state with the quantum numbers
C.   of particle L0 and the 5-momentum P0
C.   4 flavor generalization /FR'13
C-----------------------------------------------------------------------
      SAVE
      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      COMMON /S_PLIST/ P(8000,5), LLIST(8000), NP
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, kb ,kt
      COMMON /S_MASS1/ AM(99), AM2(99)
      COMMON /S_CHP/ ICHP(99), ISTR(99), IBAR(99)
      COMMON /S_CNAM/ NAMP (0:99)
      CHARACTER NAMP*6
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)
      COMMON /S_PARTO/ NFORIG(8000), NPORIG(8000), IPFLAG, NINT
      DIMENSION P0(5), LL(10), PD(10,5), IFL(3), INONLEAD(2)
      DIMENSION LRESCHEX(6:39), LRES(6:39), LCON(6:43), LPIC(-1:1)

c     charge exchange map
      DATA LCON /7,6,6,22,21,9,9,14,13,4*0,20,19,10,9,23,24,27,27,25,
     &     31,30,29,28,32,33,35,34,35,38,37,39,41,42,41,42/
c     pion charge conversion map
      DATA LPIC /8,6,7/
c     charge exchange to resonances map
      DATA LRESCHEX /26,27,27,30,31,9,9,42,41,19*0,45,44,45,48,47,39/ 
c     resonance excitation map
      DATA LRES /27,25,26,28,29,9,9,41,42,19*0,44,45,46,47,48,39/

c...  charge exchange reaction rate
c      DATA PCHEX /0.33/
      PCHEX = PAR(61)

c...  suppression of high mass particles in fireball
c     xmpsuppr = prob. accepting additional proton
      XMPSUPPR=PAR(33)
      IF(ABS(XMPSUPPR).lt.1.e-3) THEN
         WRITE(LUN,*)
     &        'Error: too low mass suppression in 4 flv fireball!'
         WRITE(LUN,*)
     &        'Probably PAR(33)/IPAR(14) not properly set, aborting..'
         STOP
      ENDIF
      XTEMPH=(AM(6)-AM(13))/LOG(XMPSUPPR)

      IF(Ndebug.gt.1) THEN
         WRITE(LUN,*)' FIRBALL_4FLV: called with (L0,P0):',
     &        L0,P0
         WRITE(LUN,*)' 2nd Proton rejection prob.:',XMPSUPPR
         WRITE(LUN,*)' fireball temperature:',XTEMPH
         WRITE(LUN,*)' charge exchange prob.:',PCHEX
      ENDIF

c...  special vector resonance treatment for meson projectiles
c     i.e. spin exchange probability
      PAR5def = PAR(5)
      IF(IPAR(14).eq.-2.and.abs(kb).lt.13)THEN
         PAR(5)=PAR(34)
      ENDIF

      NTRY=0
 100  NTRY=NTRY+1
      IF(NTRY.GT.20)THEN
         WRITE(LUN,*)' FIRBALL_4FLV: unable to sample 4flv fireball!'
         WRITE(LUN,*)' lacking rejection mechanism, abort..'
         xa=-1.
         xa = log(xa)
         STOP
c         RETURN
      ENDIF

      LA = ABS(L0)
      ISGN = ISIGN(1,L0)
      DELTAE = P0(5) - AM(LA)
      IF(DELTAE.lt.AM(6)+0.02)THEN
         IREJ = 1
         IF(ndebug.ne.0)
     &   WRITE(LUN,*)'FIRBALL_4FLV:  too low mass!! aborting...',IREJ
c         xa=-1.
c         xa=log(xa)
c         stop        
         RETURN
      ENDIF
      AV = 2.*SQRT(DELTAE)

c...  select number of particles in fireball
c     at least two
 200  NPI = AV*(1.+PAR(38)*GASDEV(0))
      IF(Ndebug.gt.3) 
     &     WRITE(LUN,*)'AV,NPI,P0(5),limit',
     &     AV,NPI,P0(5),AM(LA)+(NPI-1)*AM(7)+0.02
      IF(NPI.LE.1.OR.NPI.GT.9.OR.AM(LA)+(NPI-1)*AM(6)+0.02
     .     .GT.P0(5)) GOTO 200
      IF(Ndebug.gt.2) 
     &     WRITE(LUN,*)' FIRBALL_4FLV: No. of particles sampled. ',
     &     '(NPI,DELTAE,NTRY):',NPI,DELTAE,NTRY

c...  sample particle list      
      NTRYL=0
 210  CONTINUE
c...  special vector resonance treatment with meson projectile
      IF(IPAR(14).eq.-3.and.LA.lt.13)THEN
c     form resonance from meson beam
         IF(NTRY.GT.5) GOTO 211
         I=1
         IF(PCHEX.gt.S_RNDM(LA))THEN
            LL(I)=LRESCHEX(LA)
            CALL HSPLI(LCON(LA),IFL1,IFL2)
            IF(IFL1.eq.0.and.IFL2.eq.0)THEN
               WRITE(LUN,*) 'FIREBALL_4FLV: HSPLI rejection in ',
     &              'resonance-charge exchange'
               WRITE(LUN,*) ' RUN (SQS,KB,KT):',SQS,KB,KT
               WRITE(LUN,*) ' CONFIG: (LA,LCON,LRESCHEX)',
     &              LA,LCON(LA),LRESCHEX(LA)
               STOP
            ENDIF
            IFL(1)=IFL1
            IFL(2)=IFL2
         ELSE
            LL(I)=LRES(LA)
            CALL HSPLI(L0,IFL1,IFL2)
            IF(IFL1.eq.0.and.IFL2.eq.0)THEN
               WRITE(LUN,*) 'FIREBALL_4FLV: HSPLI rejection!'
               WRITE(LUN,*) ' RUN (SQS,KB,KT):',SQS,KB,KT
               STOP
            ENDIF
            IFL(1)=-IFL1
            IFL(2)=-IFL2
         ENDIF
         WREM = P0(5)-AM(ABS(LL(1)))
         WREM2 = AM2(ABS(LL(1)))
         INONLEAD(1)=1
         INONLEAD(2)=1
      ELSE
c...  baryon projectile
c     first two particles defined by charge exchange
         I=1
         IF(PCHEX.gt.S_RNDM(LA))THEN
            L1=LCON(LA)
            if(la.eq.42) l1 = l1 + 2 * int(2.*S_RNDM(LA))
            LL(I)=L1*ISGN
c            print*,'charge exchange!',ISGN*LA,'->',L1
         ELSE
            L1=LA
            LL(I)=LA*ISGN
         ENDIF
c     determine remaining charge
         IDQ=ICHP(LA)*ISGN-ICHP(L1)*ISIGN(1,LL(I))
         LL(I+1)=LPIC(IDQ)      ! compensate with meson
         IF(NPI.eq.2) GOTO 300
c     split last hadron again to start hadron chain
 211     CALL HSPLI (LL(I+1),IFL(1),IFL(2))
         IF(IFL(1).eq.0.and.IFL(2).eq.0)THEN
            WRITE(LUN,*) 'FIREBALL_4FLV: HSPLI rejection in chain head'
            WRITE(LUN,*) ' RUN (SQS,KB,KT):',SQS,KB,KT
            STOP
         ENDIF

         IF(Ndebug.gt.3) 
     &        WRITE(LUN,*)' FIRBALL_4FLV: Input hadron split. ',
     &        '(L0,IFL1,IFL2):',LL(I+1),IFL(1),IFL(2)
         WREM = P0(5)
         WREM2 = AM2(ABS(LL(1)))
         INONLEAD(1)=0
         INONLEAD(2)=0
      ENDIF

      IF(NTRYL.gt.20) GOTO 100
      NTRYL=NTRYL+1

 230  I=I+1    
      JT=1.5+S_RNDM(0)                            
      JR=3-JT
      NTRYS=0
      IFLB=IFL(JT)
 240  CALL SIB_I4FLAV (IFL(JT), 0, IFL(3), LL(I))
      IF(NTRYS.gt.50) GOTO 210    
      NTRYS=NTRYs+1
      W=EXP(-AM(ABS(LL(I)))/XTEMPH)
      IF(Ndebug.gt.5) 
     &     WRITE(LUN,*)' FIRBALL_4FLV: flavor added: ',
     &     '(I,NTRYS,LL(I),IFL3,W):',I,NTRYS,LL(I),IFL(3),W
      IF(W.LT.S_RNDM(I).and.INONLEAD(JT).eq.1) GOTO 240

c...  kinematic limits...     
      WREM = WREM-AM(ABS(LL(I)))
      WREM2_2=WREM2+2*SQRT(WREM2)*AM(ABS(LL(I)))+AM2(ABS(LL(I)))
      IF(Ndebug.gt.5) 
     &  WRITE(LUN,*)' FIRBALL_4FLV: kinematic limits: ',
     &  '(I,NTRYS,P05**2,WREM2):',I,NTRYS,P0(5)**2,WREM2_2
      IF(WREM2_2+0.2*S_RNDM(I).ge.P0(5)**2) GOTO 240
      WREM2=WREM2_2
      IF(Ndebug.gt.3) 
     &     WRITE(LUN,*)
     &     ' FIRBALL_4FLV: Hadron added: (KF,NAMP,I,NONlead,WRME2)',
     &     LL(I),NAMP(ABS(LL(I))),I,INONLEAD(JT),WREM2

      IFL(JT)=-IFL(3)
      INONLEAD(JT)=1
      IF(I.lt.NPI-1) GOTO 230
      IF(ABS(IFL(JT)).gt.3.and.ABS(IFL(JR)).gt.3) THEN
         IFL(JT)=IFLB
         GOTO 240
      ENDIF

c...  close list
      I=I+1
      NTRYC=0
c$$$      IAFL1 = IABS(mod(IFL(JR),100))
c$$$      IAFL2 = IABS(mod(IFL(jt),100))
c$$$      IF ((IAFL1/10.eq.4.or.mod(IAFL1,10).eq.4)
c$$$     +     .and.(IAFL2/10.eq.4.or.mod(IAFL2,10).eq.4))
c$$$     +     GOTO 100             ! reject two charm quarks
c$$$      IF(IAFL1*IAFL2.GT.100)  GOTO 100
 250  CALL SIB_I4FLAV (IFL(JT), IFL(JR), IFL(3), LL(I))
      IF(NTRYC.gt.10) GOTO 210
      NTRYC=NTRYC+1
      WREM2_2=WREM2+2*SQRT(WREM2)*AM(ABS(LL(I)))+AM2(ABS(LL(I)))
      IF(Ndebug.gt.5) 
     &     WRITE(LUN,*)' FIRBALL_4FLV: closing List: (IFL1,IFL2,KF,',
     &     'NAMP,I,NTRYC,WREM2)',
     &     IFL(JT),IFL(JR),LL(I),NAMP(ABS(LL(I))),I,NTRYC,WREM2_2

      IF(WREM2_2+0.2*S_RNDM(I).ge.P0(5)**2) GOTO 250

 300  IF(Ndebug.gt.2) 
     &     WRITE(LUN,*)
     &     ' FIRBALL_4FLV: flavors sampled. (LL,NPI,WREM,NTRYL):',
     &     LL,NPI,WREM,NTRYL

c...  fill phasespace
      CALL DECPAR (0,P0,NPI,LL,PD)
      DO J=1,NPI
         NP = NP+1
         LLIST(NP) = LL(J)
         NPORIG(NP) = IPFLAG*2
         DO K=1,5
            P(NP,K) = PD(J,K)
         ENDDO
      ENDDO
      PAR(5)=PAR5def
      IREJ = 0
      RETURN
      END


      SUBROUTINE CUT_PRO (L, SQS, PTmin, NSOFR, NJETR)
C-----------------------------------------------------------------------
C...Generate a number of soft/hard (jet-)pairs for a 'projectile'
C.  (K=1:p),(K=2:pi) interacting with a nucleon at sqrt(s)=SQS(GeV)
C-----------------------------------------------------------------------
      SAVE

      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)
      PARAMETER (NS_max = 20, NH_max = 80)
      COMMON /S_CCSIG/ SSIG(61,3), PJETC(0:NS_max,0:NH_max,61,2),
     &     SSIGN(61,3), SSIGNSD(61,3), ALINT(61,3),
     &     ASQSMIN, ASQSMAX, DASQS, NSQS
      COMMON /S_CUTOFF/ STR_mass_val, STR_mass_val_hyp, STR_mass_sea

      K = L
      if(K.eq.3) K = 2

      AL = LOG10 (SQS)
      IF (AL .LT. ASQSMIN)  THEN
          WRITE(LUN,*)  ' CUT_PRO:  low sqrt(s) ', SQS
          NSOFR = 1
          NJETR = 0
          RETURN
      ENDIF
      IF (AL .GT. ASQSMAX)  THEN
          WRITE(LUN,*)  ' CUT_PRO:  sqrt(s) out of bounds ', SQS
          NJETR = 0
          RETURN
      ENDIF

      J1 = (AL - ASQSMIN)/DASQS + 1
      J1 = MIN(J1,60)
      J1 = MAX(J1,1)
      J2 = J1+1
      T = (AL-ASQSMIN)/DASQS - FLOAT(J1-1)

      R = 0.9999*S_RNDM(0)
      DO I=0,NS_max
        DO J=0,NH_max
          IF (R.LT.(1.-T)*PJETC(I,J,J1,K)+T*PJETC(I,J,J2,K)) GOTO 100
        ENDDO
      ENDDO
100   CONTINUE

C...phase space limitation

 120  CONTINUE
      XM = FLOAT(2*I)*STR_mass_sea + FLOAT(2*J)*PTmin
      PACC = EXP(PAR(9)*(2.-XM)/SQS)
      IF(S_RNDM(0).GT.PACC) THEN
        IF(I+J.GT.1) THEN
          IF(I.GT.0) THEN
            I = I-1
            GOTO 120
          ELSE IF(J.GT.0) THEN
            J = J-1
            GOTO 120
          ENDIF
        ENDIF
      ENDIF

      NSOFR = I
      NJETR = J

      if(Ndebug.gt.2) 
     &     WRITE(LUN,*)' CUT_PRO: (L,SQS,PTmin,Ns,Nh)',K,SQS,PTmin,I,J

      END


C===========================================================================
C  Code for initialization
C===========================================================================


      SUBROUTINE SIBYLL_INI
C-----------------------------------------------------------------------
C...Initialization routine for SYBILL 
C.  
C.  the routine fills the COMMON block /CCSIG/ that contains
C.  important information for the generation of events
C.
C     PARAMETER (NS_max = 20, NH_max = 80)
C     COMMON /S_CCSIG/ SSIG(61,3), PJETC(0:NS_max,0:NH_max,61,2),
C    &    SSIGN(61,3),SSIGNSD(61,3) ALINT(61,3), ASQSMIN, ASQSMAX, DASQS, NSQS
C.
C.  NSQS = number of energy points  (61 is current version)
C.  ASQSMIN = log_10 [sqrt(s) GeV]   minimum value
C.  ASQSMIN = log_10 [sqrt(s) GeV]   maximum value
C.  DASQS   = step  in log_10[sqrt(s)]
C.            DASQS = (ASQSMAX - ASQSMIN)/(NSQS-1)
C.
C.  SSIG(J,1) inelastic cross section for pp interaction
C.            at energy: sqrt(s)(GeV) = 10**[ASQSMIN+DASQS*(J-1)]
C.  SSIG(J,2)  inelastic cross section for pi-p interaction
C.  SSIGN(J,1) inelastic cross section for p-Air interaction
C.  SSIGN(J,2) inelastic cross section for pi-Air interaction
C.
C.  PJETC(n_s,n_j,J,1) Cumulative  probability distribution
C.                 for the production of n_s soft interactions and
C.                 n_j (n_j=0:30) jet pairs at sqrt(s) labeled 
C.                 by J, for p-p interaction
C.  PJETC(n_s,n_j,J,2) Same as above for pi-p interaction
C.  ALINT(J,1)   proton-air  interaction length (g cm-2)
C.  ALINT(J,2)   pi-air  interaction length (g cm-2)
C-----------------------------------------------------------------------
      SAVE

      WRITE(*,100)
 100  FORMAT(' ','====================================================',
     *     /,' ','|                                                  |',
     *     /,' ','|                 S I B Y L L  2.3rc1              |',
     *     /,' ','|                                                  |',
     *     /,' ','|         HADRONIC INTERACTION MONTE CARLO         |',
     *     /,' ','|                        BY                        |',
     *     /,' ','|            Eun-Joo AHN, Felix RIEHN              |',
     *     /,' ','|     R. ENGEL, R.S. FLETCHER, T.K. GAISSER        |',
     *     /,' ','|               P. LIPARI, T. STANEV               |',
     *     /,' ','|                                                  |',
     *     /,' ','| Publication to be cited when using this program: |',
     *     /,' ','| R. Engel et al., Proc. 26th ICRC, 1 (1999) 415   |',
     *     /,' ','|                                                  |',
     *     /,' ','| last tampered with by: F. Riehn (Aug. 8th 2014)  |',
     *     /,' ','|     --> isvhecri 2014 release candidate <--      |',
     *     /,' ','====================================================',
     *     /)

      CALL PAR_INI
      CALL JET_INI
      CALL BLOCK_INI
      CALL NUC_GEOM_INI
      CALL SIG_AIR_INI
c...  charm frag. normalisation
      call znormal

      END

      subroutine dec_ini
C-----------------------------------------------------------------------
C     decay initialization routine
C     sets which particles should decay and wich should be stable
C-----------------------------------------------------------------------
      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)
      COMMON /S_CSYDEC/ CBR(223), KDEC(1338), LBARP(99), IDB(99)
      
      WRITE(LUN,*)'-----------------------------------------'
      WRITE(LUN,*)'SIBYLL DECAYINI: setting particle decays!'
      WRITE(LUN,*)'-----------------------------------------'
C...  Definition of stable particles
      DO J=4,12
         IDB(J) = -abs(IDB(J))
      ENDDO
c----------------------------------------------------------
c     if the folowing is commented out then all particles
c     except leptons,protons and neutrons are UNSTABLE
c----------------------------------------------------------
c     all particles with t<0.3e-10s are considered unstable
c     i.e. all the mesons from K0s onwards(K0l is stable)
c----------------------------------------------------------
C     K0s stable
      WRITE(LUN,*)'making K0s stable..'
      IDB(12) = -abs(IDB(12))

C     Lambda/Anti-lambda stable
      WRITE(LUN,*)'making LAMBDA stable..'
      IDB(39) = -abs(IDB(39))

c     Sigmas stable
      WRITE(LUN,*)'making SIGMAs stable..'
      do i=34,36
         IDB(i) = -abs(IDB(i))
      enddo
C     Eta stable
c      WRITE(LUN,*)'making eta stable..'
c      IDB(23) = -abs(IDB(23))
      WRITE(LUN,*)'------------------------------------------'
      end

      SUBROUTINE PAR_INI
C------------------------------------------------------------------
C.    parameter config: 
C     aug14.2strgJet3.0vecJet1diqEnhc2chmRe6.3
C------------------------------------------------------------------
      SAVE
      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)
      COMMON /S_CZDIS/ FA, FB0
      COMMON /S_CZDISs/ FAs1, fAs2
      COMMON /S_CUTOFF/ STR_mass_val, STR_mass_val_hyp, STR_mass_sea
      COMMON /S_CPSPL/ CCHIK(3,39)
      PARAMETER ( NPARFIT = 22 )
      COMMON /XSCTN_FIT/ PARS( 50 , 2 )
      DOUBLE PRECISION PARS
      DATA (PARS(K,1),K=    1,NPARFIT) /
     &3.9223E+01,4.2055E+01,5.0913E-02,-4.0000E-01,2.000E-01,5.0000E-01,
     &0.0000E+00,6.0000E-01,9.0000E-02,1.0000E+00,2.0000E+00,3.2327E+00,
     &2.5000E-01,5.4000E-01,1.0000E+00,-8.8000E-01,5.400E-01,5.0000E-01,
     &9.0000E-01,5.4000E-01,6.5000E-02,9.0000E-01/

      DATA (PARS(K,2),K=    1,NPARFIT) /
     &2.0590E+01,9.6579E+01,5.6069E-02,-7.6393E-01,2.000E-01,5.0000E-01,
     &0.0000E+00,6.0000E-01,9.0000E-02,1.0000E+00,2.0000E+00,2.9191E+00,
     &2.5000E-01,5.4000E-01,1.0000E+00,-8.8000E-01,5.400E-01,5.4895E-01,
     &9.0000E-01,5.4000E-01,6.5000E-02,9.0000E-01/

C...  Model parameters
      PAR( 1 ) =  0.04
      PAR( 2 ) =  0.3
      PAR( 3 ) =  0.3
      PAR( 4 ) =  0.14
      PAR( 5 ) =  0.3
      PAR( 6 ) =  0.3
      PAR( 7 ) =  0.15
      PAR( 8 ) =  0.15
      PAR( 9 ) =  7.0
      PAR( 10 ) =  1.0
      PAR( 11 ) =  0.065
      PAR( 12 ) =  0.9
      PAR( 13 ) =  0.2
      PAR( 14 ) =  0.06
      PAR( 15 ) =  0.1
      PAR( 16 ) =  0.04
      PAR( 17 ) =  0.04
      PAR( 18 ) =  0.5
      PAR( 19 ) =  0.8
      PAR( 20 ) =  0.8
      PAR( 21 ) =  0.8
      PAR( 22 ) =  4.0
      PAR( 23 ) =  0.5
      PAR( 24 ) =  0.004
      PAR( 25 ) =  0.007
      PAR( 26 ) =  25.0
      PAR( 27 ) =  0.075
      PAR( 28 ) =  10.0
      PAR( 29 ) =  0.0
      PAR( 30 ) =  2.0
      PAR( 31 ) =  0.33
      PAR( 32 ) =  0.0
      PAR( 33 ) =  0.1
      PAR( 34 ) =  0.7
      PAR( 35 ) =  0.0
      PAR( 36 ) =  0.35
      PAR( 37 ) =  0.0
      PAR( 38 ) =  0.5
      PAR( 39 ) =  0.0
      PAR( 40 ) =  0.0
      PAR( 41 ) =  1.0
      PAR( 42 ) =  3.0
      PAR( 43 ) =  0.0
      PAR( 44 ) =  0.99
      PAR( 45 ) =  1.0
      PAR( 46 ) =  0.16
      PAR( 47 ) =  0.35
      PAR( 48 ) =  0.4
      PAR( 49 ) =  0.15
      PAR( 50 ) =  0.6
      PAR( 51 ) =  1.0
      PAR( 52 ) =  0.0
      PAR( 53 ) =  6.0
      PAR( 54 ) =  0.2
      PAR( 55 ) =  0.0
      PAR( 56 ) =  0.0
      PAR( 57 ) =  -0.9
      PAR( 58 ) =  0.8
      PAR( 59 ) =  0.4
      PAR( 60 ) =  0.8
      PAR( 61 ) =  0.33
      PAR( 62 ) =  3.0
      PAR( 63 ) =  1.0
      PAR( 64 ) =  0.25
      PAR( 65 ) =  0.3
      PAR( 66 ) =  0.3
      PAR( 67 ) =  0.6
      PAR( 68 ) =  0.006
      PAR( 69 ) =  0.01
      PAR( 70 ) =  0.007
      PAR( 71 ) =  3.0
      PAR( 72 ) =  0.3
      PAR( 73 ) =  0.5
      PAR( 74 ) =  0.6
      PAR( 75 ) =  0.0
      PAR( 76 ) =  0.0
      PAR( 77 ) =  0.0
      PAR( 78 ) =  0.0
      PAR( 79 ) =  0.0
      PAR( 80 ) =  0.0
      IPAR( 1 ) =  1
      IPAR( 2 ) =  0
      IPAR( 3 ) =  8
      IPAR( 4 ) =  0
      IPAR( 5 ) =  1
      IPAR( 6 ) =  0
      IPAR( 7 ) =  0
      IPAR( 8 ) =  1
      IPAR( 9 ) =  1
      IPAR( 10 ) =  1
      IPAR( 11 ) =  -16
      IPAR( 12 ) =  1
      IPAR( 13 ) =  0
      IPAR( 14 ) =  -2
      IPAR( 15 ) =  9
      IPAR( 16 ) =  7
      IPAR( 17 ) =  1
      IPAR( 18 ) =  4
      IPAR( 19 ) =  1
      IPAR( 20 ) =  3
      IPAR( 21 ) =  0
      IPAR( 22 ) =  0
      IPAR( 23 ) =  0
      IPAR( 24 ) =  0
      IPAR( 25 ) =  0
      IPAR( 26 ) =  0
      IPAR( 27 ) =  0
      IPAR( 28 ) =  2
      IPAR( 29 ) =  0
      IPAR( 30 ) =  0
      IPAR( 31 ) =  1
      IPAR( 32 ) =  0
      IPAR( 33 ) =  0
      IPAR( 34 ) =  0
      IPAR( 35 ) =  0
      IPAR( 36 ) =  1
      IPAR( 37 ) =  0
      IPAR( 38 ) =  0
      IPAR( 39 ) =  0
      IPAR( 40 ) =  0
      IPAR( 41 ) =  1
      IPAR( 42 ) =  3
      IPAR( 43 ) =  1
      IPAR( 44 ) =  0
      IPAR( 45 ) =  0
      IPAR( 46 ) =  0
      IPAR( 47 ) =  0
      IPAR( 48 ) =  0
      IPAR( 49 ) =  0
      IPAR( 50 ) =  0
      IPAR( 51 ) =  0
      IPAR( 52 ) =  0
      IPAR( 53 ) =  0
      IPAR( 54 ) =  0
      IPAR( 55 ) =  0
      IPAR( 56 ) =  0
      IPAR( 57 ) =  0
      IPAR( 58 ) =  0
      IPAR( 59 ) =  0
      IPAR( 60 ) =  0
      IPAR( 61 ) =  0
      IPAR( 62 ) =  0
      IPAR( 63 ) =  0
      IPAR( 64 ) =  0
      IPAR( 65 ) =  0
      IPAR( 66 ) =  0
      IPAR( 67 ) =  0
      IPAR( 68 ) =  0
      IPAR( 69 ) =  0
      IPAR( 70 ) =  0
      IPAR( 71 ) =  0
      IPAR( 72 ) =  0
      IPAR( 73 ) =  0
      IPAR( 74 ) =  0
      IPAR( 75 ) =  0
      IPAR( 76 ) =  0
      IPAR( 77 ) =  0
      IPAR( 78 ) =  0
      IPAR( 79 ) =  0
	  IPAR( 80 ) =  0

C...energy dependence of PTmin
c     pt_cut offset
      PAR(10) = PARS(10 , 1)
c     lambda
      PAR(11) = PARS(21 , 1)
c     c parameter
      PAR(12) = PARS(22 , 1)

C...fragmentation function
      FA = PAR(20)
      FB0 = PAR(21)

C...Strange fragmentation function
      FAs1 = PAR(35)
      FAs2 = PAR(35)

C...mass cut
      STR_mass_val = PAR(36) 
      STR_mass_sea = PAR(41)

C...  valence quark distribution function
c     large x suppression
      do i=1,3                  ! quark flavors
         CCHIK(i,13)=PAR(62)
         CCHIK(i,14)=PAR(62)
      enddo

C...leading baryon fragmentation function
c     hard proton mixing
      CLEAD = PAR(50)

      END
      

      SUBROUTINE JET_INI
C-----------------------------------------------------------------------
C...Compute table of cross sections, and table of probability
C.  for the production of multiple soft and hard interactions
C.
C.  The output of this routine  is the COMMON block /S_CCSIG/
C.  that contains  the cross sections h-p, h-Air, and the 
C.  cumulative probability of NS soft and NH hard interactions
C-----------------------------------------------------------------------
      SAVE

      PARAMETER (NS_max = 20, NH_max = 80)
      COMMON /S_CCSIG/ SSIG(61,3), PJETC(0:NS_max,0:NH_max,61,2),
     &     SSIGN(61,3), SSIGNSD(61,3), ALINT(61,3),
     &     ASQSMIN, ASQSMAX, DASQS, NSQS
      COMMON /S_CCSIG2/ SSIG_TOT(61,3),SSIG_SD1(61,3),SSIG_SD2(61,3),
     &    SSIG_DD(61,3),SSIG_B(61,3),SSIG_RHO(61,3)

      COMMON /S_CCSIG3/ SSIG_SD1LM(61,3),SSIG_SD1HM(61,3),
     &     SSIG_SD2LM(61,3),SSIG_SD2HM(61,3),
     &     SSIG_DDLM(61,3),SSIG_DDHM(61,3)

      DIMENSION Pjet(0:NS_max,0:NH_max)
      DIMENSION SIG_df(3),SIG_df2(3,2),SIGDIF(3),SIGDIF_pi(3),
     &          PS_tab(61),PH_tab(61),PT_tab(61)

      COMMON /S_DEBUG/ Ncall, Ndebug, Lun

C...spacing in energy for table of cross sections.

      NSQS = 61
      ASQSMIN = 1.
      ASQSMAX = 7.
      DASQS = (ASQSMAX-ASQSMIN)/FLOAT(NSQS-1)

C...initialization of proton and pion tables

      DO KK=1,2

         WRITE(LUN,'(2(/,1X,A,A))') 
     &     'Table: J, sqs,  PT_cut,  SIG_tot,  SIG_inel,  B_el,  ',
     &     'rho,  <n_s>,  <n_h>,  SIG_SD,  SD1_lm,  SD1_hm',
     &     '-----------------------------------------------------',
     &     '----------------------------------------------'

         JINT = KK
         DO J=1, NSQS
           ASQS = ASQSMIN + DASQS*FLOAT(J-1)
           SQS = 10.**ASQS

           CALL SIB_SIG (JINT, SQS, PTmin,
     &                   SIG_tot, SIG_inel, SIG_df, SIG_df2, B_el, Pjet)

C...low-energy interpolation with data-parametrizations
           call SIB_HADCSL(JINT,SQS,
     &                     SIGTOT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO)
           if(SQS.le.100.) then
             SIG_TOT  = SIGTOT
             SIG_inel = SIGINEL
             B_EL     = SLOPE
           else if(SQS.le.1000.) then
             Xi = log(SQS/100.)/2.30258509299405
             SIG_TOT  = Xi*SIG_TOT+(1.-Xi)*SIGTOT
             SIG_inel = Xi*SIG_inel+(1.-Xi)*SIGINEL
             B_EL     = Xi*B_EL+(1.-Xi)*SLOPE
           endif

           SSIG_TOT(J,KK) = SIG_TOT
           SSIG(J,KK)     = SIG_inel
           SSIG_SD1(J,KK) = SIGDIF(1)
           SSIG_SD2(J,KK) = SIGDIF(2)
           SSIG_DD(J,KK)  = SIG_df(3)
           SSIG_B(J,KK)   = B_EL
           SSIG_RHO(J,KK) = RHO

           SSIG_SD1LM(J,KK) = SIG_df2(1,1)
           SSIG_SD1HM(J,KK) = SIG_df2(1,2)
           SSIG_SD2LM(J,KK) = SIG_df2(2,1)
           SSIG_SD2HM(J,KK) = SIG_df2(2,2)
           SSIG_DDLM(J,KK) = SIG_df2(3,1)
           SSIG_DDHM(J,KK) = SIG_df2(3,2)

           PSUM = 0.
           PH = 0.
           PS = 0.
           DO NS=0,NS_max
             DO NJ=0,NH_max

               PS = PS+FLOAT(NS)*Pjet(NS,NJ)
               PH = PH+FLOAT(NJ)*Pjet(NS,NJ)

               PSUM = PSUM+Pjet(NS,NJ)
               PJETC(NS,NJ,J,KK) = PSUM

             ENDDO
           ENDDO
           PS_tab(J) = PS
           PH_tab(J) = PH
           PT_tab(J) = PTmin

           WRITE(LUN,'(3X,I2,1P,E12.3,0P,4F8.2,6F8.3)') 
     &       JINT,SQS,PTmin,SIG_tot,SIG_inel,B_el,RHO,PS,PH
     &          ,SIGDIF(1)+SIGDIF(2),SIG_df2(1,1),SIG_df2(1,2)

         ENDDO
      ENDDO

C...initialization of kaon tables

      JINT = 3

      WRITE(LUN,'(2(/,1X,A,A))') 
     &  'Table: J, sqs,  PT_cut,  SIG_tot,  SIG_inel,  B_el,  ',
     &  'rho,  <n_s>,  <n_h>',
     &  '-----------------------------------------------------',
     &  '-------------------'
      DO J=1, NSQS
        ASQS = ASQSMIN + DASQS*FLOAT(J-1)
        SQS = 10.**ASQS
C...use pion cross section rescaled for high-energy extrapolation
        SIG_tot   = SSIG_TOT(J,2)
        SIG_inel  = SSIG(J,2)
        SIG_df(1) = SSIG_SD1(J,2)
        SIG_df(2) = SSIG_SD2(J,2)
        SIG_df(3) = SSIG_DD(J,2)
        B_el = SSIG_B(J,2)
        PTmin = PT_tab(J)
        PS = PS_tab(J)
        PH = PH_tab(J)

C...low-energy interpolation with data-parametrizations
        call SIB_HADCSL(2,SQS,
     &                  SIGTOT_pi,SIGEL_pi,SIGINEL,SIGDIF_pi,SLOPE,RHO)
        call SIB_HADCSL(3,SQS,
     &                  SIGTOT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO)
        SIG_el    = (SIGEL/SIGEL_pi)*(SIG_TOT-SIG_inel)
        SIG_TOT   = (SIGTOT/SIGTOT_pi)*SIG_TOT
        SIG_inel  = SIG_TOT-SIG_el
        SIG_df(3) = (SIGDIF(3)/SIGDIF_pi(3))*SIG_df(3)
        if(SQS.le.100.) then
          SIG_TOT  = SIGTOT
          SIG_inel = SIGINEL
          B_EL     = SLOPE
        else if(SQS.le.1000.) then
          Xi = log(SQS/100.)/2.30258509299405
          SIG_TOT  = Xi*SIG_TOT+(1.-Xi)*SIGTOT
          SIG_inel = Xi*SIG_inel+(1.-Xi)*SIGINEL
          B_EL     = Xi*B_EL+(1.-Xi)*SLOPE
        endif

        SSIG_TOT(J,3) = SIG_TOT
        SSIG(J,3)     = SIG_inel
        SSIG_SD1(J,3) = SIGDIF(1)
        SSIG_SD2(J,3) = SIGDIF(2)
        SSIG_DD(J,3)  = SIG_df(3)
        SSIG_B(J,3)   = B_EL
        SSIG_RHO(J,3) = RHO

        WRITE(LUN,'(3X,I2,1P,E12.3,0P,4F8.2,3F8.3)') 
     &    JINT,SQS,PTmin,SIG_tot,SIG_inel,B_el,RHO,PS,PH

      ENDDO


      END


      SUBROUTINE INI_WRITE (LUN)
C-----------------------------------------------------------------------
C   This subroutine prints on unit LUN
C   a table of the cross sections  used in the program
C   and of the average number of hard interactions, and the average
C   number of wounded nucleons in a hadron-air interaction
C-----------------------------------------------------------------------
      SAVE

      PARAMETER (NS_max = 20, NH_max = 80)
      COMMON /S_CCSIG/ SSIG(61,3), PJETC(0:NS_max,0:NH_max,61,2),
     &     SSIGN(61,3), SSIGNSD(61,3), ALINT(61,3),
     &     ASQSMIN, ASQSMAX, DASQS, NSQS
      DIMENSION PJ(2),PS(2),PW(2)

      DATA ATARG /14.514/

*      CALL PARAM_PRINT(LUN)
      WRITE (LUN, 10)
      WRITE (LUN, 15)
      WRITE (LUN, 16)
      WRITE (LUN, 18)
10    FORMAT(//,' Table of cross sections, and average number',
     &         ' of minijets and wounded nucleons ')
15    FORMAT('        [sqrt(s) in GeV, cross sections in mbarn]. ')
16    FORMAT(' sqrt(s) sig(pp) sig(pA) <n_s> <n_j> <n_w>',
     &    ' sig(pip) sig(piA) <n_s> <n_j> <n_w>')
18    FORMAT(1X,77('-') )
      DO J=1,61,1
         SQS = 10.**(ASQSMIN + DASQS*FLOAT(J-1))

         DO K=1,2

           PW(K) = ATARG*SSIG(J,K)/SSIGN(J,K)

           PJ(K) = 0.
           PS(K) = 0.
           DO NS=0,NS_max
             DO NJ=0,NH_max
               IF(NJ.GT.0) THEN
                 PROB = PJETC(NS,NJ,J,K) - PJETC(NS,NJ-1,J,K)
               ELSE IF(NS.GT.0) THEN
                 PROB = PJETC(NS,NJ,J,K) - PJETC(NS-1,NH_max,J,K)
               ELSE
                 PROB = 0.
               ENDIF
               PJ(K) = PJ(K)+FLOAT(NJ)*PROB
               PS(K) = PS(K)+FLOAT(NS)*PROB
             ENDDO
           ENDDO

         ENDDO

         WRITE(LUN,20) SQS,SSIG(J,1),SSIGN(J,1),PS(1),PJ(1),PW(1)
     &                      ,SSIG(J,2),SSIGN(J,2),PS(2),PJ(2),PW(2)

      ENDDO
      WRITE (LUN, 18)
20    FORMAT (1X,E8.2, 2(2F7.1,1X,3F6.2,1X))

      END

C*************************************************************************
C=========================================================================
C. UTILITIES ROUTINES
C=========================================================================
C***********************************************************************

C=======================================================================
C. Code for the wounded nucleon distribution
C=======================================================================


      SUBROUTINE SIB_START_EV (SQS, L, IA, NW, JDIF)
C-----------------------------------------------------------------------
C...Beginning of a SIBYLL interaction 
C.
C.  add l.m. Glauber SD cross section for pAir  13/FR
C.
C.  INPUT : SQS = c.m.s. energy (GeV)
C.          L = 1:proton, 2:charged pion
C.          IA = mass of target nucleon
C. 
C.  OUTPUT: NW    = number of wounded nucleons
C.          JDIF(JW)  = diffraction code    !!!! changed to field !!!!
C.                  (0 : non-diffractive interaction)
C.                  (1 : forward diffraction)
C.                  (2 : backward diffraction)
C.                  (3 : double diffraction)
C.
C-----------------------------------------------------------------------
      SAVE

      PARAMETER (NW_max = 20)
      DIMENSION JDIF(NW_max)
      COMMON /S_CNCM0/ B, BMAX, NTRY, NA
      COMMON /S_DIFMAss/ XM2MIN(3),ALXMIN(3),SLOP0,ASLOP,BSLOP,XMASS(2)
      COMMON /GLAUB_SCR/ XI_MAX, ALAM(61)
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)
      DIMENSION SIGDIF(3)

C...sample number of wounded nucleons
      CALL SIB_SIGMA_HP(L,SQS,SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO) 
      IF (IA .GT. 1)  THEN
         IF(IPAR(12).NE.0)THEN
            CALL SIB_SIGMA_HAIR(L,SQS,SIGprod,SIGbdif)
C     2channel low-mass (coherent) diffraction?
            IF(S_RNDM(L).LT.SIGbdif/SIGprod)THEN
               NW = 1
               JDIF(1) = 1
               RETURN
            ENDIF
         ENDIF
c     sample number of wounded nucleons
         CALL INT_H_NUC (IA, SIGT, SLOPE, RHO) 
      ELSE
         NA = 1
      ENDIF      
      NW = NA

C...new treatment of diffraction 
      IF(IA.GT.1) THEN
c     nuclear target
         IF(NW.eq.1)THEN
            IF(IPAR(12).NE.0)THEN
c     high mass (incoherent) diffraction?
               S = SQS ** 2
               PF =(1.-LOG(S*XI_MAX/XM2MIN(L))/LOG(S*PAR(13)/XM2MIN(L)))
     &              * SIGDIF(1)/SIGINEL
               PB = SIGDIF(2)/SIGINEL
               PD = SIGDIF(3)/SIGINEL
            ELSE
               PF = SIGDIF(1)/SIGINEL
               PB = SIGDIF(2)/SIGINEL
               PD = SIGDIF(3)/SIGINEL
            ENDIF
         ELSE
            IF(IPAR(12).EQ.1)THEN
c     all interactions with Nw>1 are non-diff.
               DO K=1, NW
                  JDIF(K) = 0
               ENDDO
               RETURN
            ELSE
c     some Nw>1 are attached by diff. 
               PF = SIGDIF(1)/SIGINEL
               PB = SIGDIF(2)/SIGINEL
               PD = SIGDIF(3)/SIGINEL
            ENDIF
         ENDIF
      ELSE
         PF = SIGDIF(1)/SIGINEL
         PB = SIGDIF(2)/SIGINEL
         PD = SIGDIF(3)/SIGINEL
      ENDIF
      P0 = 1.-PF-PB-PD
      P1 = P0 + PF
      P2 = P1 + PB
      DO K=1, NW
         R = S_RNDM(0)
         IF (R .LT. P0)  THEN
            JDIF(K) = 0
         ELSE IF (R .LT. P1)  THEN
            JDIF(K) = 1
         ELSE IF (R .LT. P2)  THEN
            JDIF(K) = 2
         ELSE 
            JDIF(K) = 3
         ENDIF
      ENDDO
      
      END



      SUBROUTINE INT_H_NUC (IA, SIGT, SLOPE, RHO) 
C...Compute with a montecarlo method the "multiple interaction structure"
C.  of an hadron-nucleus collision.
C.  
C.
C.  INPUT : IA               = mass of target nucleus
C.          SIGT (mbarn)     = total hp cross section
C.          SLOPE (GeV**-2)  = slope of hp elastic scattering
C.          RHO              = real/imaginary part of forward elastic
C.                             scattering amplitude
C.
C.  OUTPUT : in COMMON block /CNCMS0/
C.           B = impact parameter (fm)
C.           BMAX = maximum impact parameter for generation
C.           NTRY = number of "trials" before one interaction
C.           NA = number of wounded nucleons in A
C. Author : P.Lipari  (may 1993)
C---------------------------------------------------------------------------
      SAVE
      PARAMETER (IAMAX=56)
      COMMON /S_CNCM0/ B, BMAX, NTRY, NA
      DIMENSION XA(IAMAX), YA(IAMAX)
      DATA PI /3.1415926/
      DATA CMBARN /0.389385/
      CC = SIGT/(4.*PI*SLOPE*CMBARN)         
      DEN = 2.*SLOPE*CMBARN*0.1
      BMAX = 10.                             ! fm
      NTRY = 0
      CALL NUC_CONF (IA, XA, YA)
1000  B = BMAX*SQRT(S_RNDM(0))
      PHI = 2.*PI*S_RNDM(0)
      BX = B*COS(PHI)
      BY = B*SIN(PHI)
      NTRY = NTRY+1
      NA = 0
      DO JA=1,IA
         S = (XA(JA)-BX)**2 + (YA(JA)-BY)**2
         F = EXP(-S/DEN)
         PEL = CC*CC*(1.+RHO*RHO)*F*F
         PINEL  = 2.*CC*F-PEL
         R = S_RNDM(0)
         IF (R .LT. PINEL)  THEN
            NA = NA + 1
         ENDIF
      ENDDO
      IF (NA .EQ. 0)  GOTO 1000
      RETURN
      END


C==========================================================================
C. Cross sections
C==========================================================================


      SUBROUTINE SIG_AIR_INI 
C-----------------------------------------------------------------------
C...Initialize the cross section and interaction lengths on air
C.  (this version initializes p-air, pi-air, and K-air cross sections)
C.
C.  also calculates the low mass beam diffraction cross section in hAir \FR
C.  using the same lambda for all hadrons
C-----------------------------------------------------------------------
      SAVE

      PARAMETER (NS_max = 20, NH_max = 80)
      COMMON /S_CCSIG/ SSIG(61,3), PJETC(0:NS_max,0:NH_max,61,2),
     &     SSIGN(61,3), SSIGNSD(61,3), ALINT(61,3),
     &     ASQSMIN, ASQSMAX, DASQS, NSQS
      COMMON /GLAUB_SCR/ XI_MAX , ALAM(61)
      COMMON /S_CFLAFR/ PAR(80), IPAR(80)
      COMMON /S_DEBUG/ Ncall, Ndebug, Lun

      DIMENSION SIGDIF(3)

      DATA AVOG /6.0221367E-04/

      ATARGET = 14.514

      IF ( IPAR(12).GT.0 ) THEN
      WRITE(LUN,*) '==================================================='
      WRITE(LUN,*) 'SIG_AIR_INI:'
      WRITE(LUN,*) 'using Goulianos param. for res.coupling..'
      XI_MAX = 0.02
      WRITE(LUN,*) 'low mass Xi_max: ' , XI_MAX
      WRITE(LUN,*) '==================================================='
      ENDIF
C...particle loop (p, pi, K)
      DO K=1,3
         
         WRITE(LUN,'(2(/,1X,A,A))') 
     &        'Table: J, sqs,  SIGtot,   SIGprod,  SIG_SD,  Lambda'
         WRITE(LUN,*) 
     &      '---------------------------------------------------'

        DO J=1,NSQS

           ASQS = ASQSMIN + DASQS*FLOAT(J-1)
           SQS = 10.**ASQS

           IF (K.EQ.1) THEN
c     Goulianos param. from GAP-2012-056, Mx**2s = 0.02
c     against PDG elastic cross section
              CALL SIB_HADCS1(K,SQS,SIGT1,SIGEL1,SIGINEL1,SLOPE1,RHO1)
              SIGEFF = 0.68*(1.+36./SQS**2)*log(0.6+XI_MAX/1.5*SQS**2)
              ALAM(J) = SQRT(SIGEFF/SIGEL1)
           ENDIF

           CALL SIB_SIGMA_HP(K,SQS,SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO)
           CALL SIG_H_AIR2
     &          (SIGT, SLOPE, RHO, ALAM(J),
     &          SSIGT, SSIGEL, SSIGQE, SSIGSD, SSIGQSD)
c           IF(IPAR(12).EQ.0)
c     &          CALL SIG_H_AIR(SIGT, SLOPE, RHO, SSIGT, SSIGEL, SSIGQE)

           WRITE(LUN,'(1X,1P,5E12.3)') 
     &          SQS,SSIGT,SSIGT-SSIGQE,SSIGQSD,ALAM(J)
C  particle production cross section
           SSIGN(J,K) = SSIGT-SSIGQE
           SSIGNSD(J,K) = SSIGQSD
           ALINT(J,K) = 1./(AVOG*SSIGn(j,K)/ATARGET)
        ENDDO
      ENDDO

      WRITE(LUN,*) '==================================================='
      WRITE(LUN,'(1X,A)') 
     &  'SIG_AIR_INI: NUCLIB interaction lengths (p-air, pi-air, K-air)'
      DO J=1,NSQS
         ASQS = ASQSMIN + DASQS*FLOAT(J-1)
         SQS = 10.**ASQS
         WRITE(LUN,'(1X,1P,4E12.3)')SQS,ALINT(J,1),ALINT(J,2),ALINT(J,3)
      ENDDO
      
      END


      SUBROUTINE SIB_SIGMA_HP(L,SQS,SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO)
C-----------------------------------------------------------------------
C     Hadron-proton cross sections, taken from interpolation table
C     calculated by SIBYLL_INI
C
C     input:       L     1      proton-proton
C                        2      pi-proton
C                        3      K-proton
C                  SQS   sqrt(s)
C
C     output:      SIGT       total cross section (mb)
C                  SIGEL      elastic cross section (mb)
C                  SIGINEL    inelastic cross section (mb)
C                  SIGDIF     diffraction dissociation CS (mb)
C                  SLOPE      elastic slope parameter (GeV^-2)
C                  RHO        real/imaginary part of forward amplitude
C-----------------------------------------------------------------------
      SAVE

      DIMENSION SIGDIF(3)

      PARAMETER (NS_max = 20, NH_max = 80)
      COMMON /S_CCSIG/ SSIG(61,3), PJETC(0:NS_max,0:NH_max,61,2),
     &     SSIGN(61,3), SSIGNSD(61,3), ALINT(61,3),
     &     ASQSMIN, ASQSMAX, DASQS, NSQS
      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      COMMON /S_CCSIG2/ SSIG_TOT(61,3),SSIG_SD1(61,3),SSIG_SD2(61,3),
     &    SSIG_DD(61,3),SSIG_B(61,3),SSIG_RHO(61,3)
      DIMENSION LL(39)
      DATA LL /5*0,3*2,4*3,2*1,19*0,6*1/
Cf2py intent(out) SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO

      IF(NSQS.LE.0) THEN
        WRITE(LUN,'(//,1X,A)') 
     &    'SIB_SIGMA_HP: interpolation table not initialized.'
        STOP
      ENDIF
      IF(IABS(L).gt.39)THEN
         WRITE(LUN,'(//,1X,A)')     
     &        'SIB_SIGMA_HAIR: unknown beam particle!',L
         STOP
      ENDIF
      IF(L.GT.3) L=LL(IABS(L))
      IF(L.EQ.0)THEN
         WRITE(LUN,'(//,1X,A)')     
     &        'SIB_SIGMA_HAIR: unknown beam particle!', L
         STOP
      ENDIF
        
      AL = LOG10(SQS)
      J1 = (AL - 1.)*10. + 1
      if((j1.lt.1).or.(j1.ge.NSQS)) then
        write (LUN,'(1x,a,i3,1p,e12.3)') 
     &    'SIB_SIGMA_HP: energy out of range ',L,sqs
        J1 = min(J1,NSQS-1)
        J1 = max(J1,1)
      endif
      T = (AL-1.)*10. - FLOAT(J1-1)
      SIGT    = SSIG_TOT(J1,L)*(1.-T) + SSIG_TOT(J1+1,L)*T
      SIGINEL = SSIG(J1,L)*(1.-T) + SSIG(J1+1,L)*T
      SIGEL   = SIGT-SIGINEL
      SIGDIF(1) = SSIG_SD1(J1,L)*(1.-T) + SSIG_SD1(J1+1,L)*T
      SIGDIF(2) = SSIG_SD2(J1,L)*(1.-T) + SSIG_SD2(J1+1,L)*T
      SIGDIF(3) = SSIG_DD(J1,L)*(1.-T) + SSIG_DD(J1+1,L)*T
      SLOPE   = SSIG_B(J1,L) *(1.-T) + SSIG_B(J1+1,L)*T
      RHO     = SSIG_RHO(J1,L) *(1.-T) + SSIG_RHO(J1+1,L)*T

      END


      SUBROUTINE SIB_SIGMA_HP2
     +     (L,SQS,SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO)
C-----------------------------------------------------------------------
C     Hadron-proton cross sections, taken from interpolation table
C     calculated by SIBYLL_INI
C
C     input:       L     1      proton-proton
C                        2      pi-proton
C                        3      K-proton
C                  SQS   sqrt(s)
C
C     output:      SIGT       total cross section (mb)
C                  SIGEL      elastic cross section (mb)
C                  SIGINEL    inelastic cross section (mb)
C                  SIGDIF     diffraction dissociation CS (mb)
C                             split in high and low mass !!
C                             ( taken from S_CCSIG3 )
C                  SLOPE      elastic slope parameter (GeV^-2)
C                  RHO        real/imaginary part of forward amplitude
C-----------------------------------------------------------------------
      SAVE

      DIMENSION SIGDIF(3,2)

      PARAMETER (NS_max = 20, NH_max = 80)
      COMMON /S_CCSIG/ SSIG(61,3), PJETC(0:NS_max,0:NH_max,61,2),
     &     SSIGN(61,3), SSIGNSD(61,3), ALINT(61,3),
     &     ASQSMIN, ASQSMAX, DASQS, NSQS
      COMMON /S_CCSIG2/ SSIG_TOT(61,3),SSIG_SD1(61,3),SSIG_SD2(61,3),
     &    SSIG_DD(61,3),SSIG_B(61,3),SSIG_RHO(61,3)
      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      COMMON /S_CCSIG3/ SSIG_SD1LM(61,3),SSIG_SD1HM(61,3),
     &     SSIG_SD2LM(61,3),SSIG_SD2HM(61,3),
     &     SSIG_DDLM(61,3),SSIG_DDHM(61,3)
Cf2py intent(out) SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO

      IF(NSQS.LE.0) THEN
        WRITE(LUN,'(//,1X,A)') 
     &    'SIB_SIGMA_HP2: interpolation table not initialized.'
        STOP
      ENDIF
        
      AL = LOG10(SQS)
      J1 = (AL - 1.)*10. + 1
      if((j1.lt.1).or.(j1.ge.NSQS)) then
        write (LUN,'(1x,a,i3,1p,e12.3)') 
     &    'SIB_SIGMA_HP2: energy out of range ',L,sqs
        J1 = min(J1,NSQS-1)
        J1 = max(J1,1)
      endif
      T = (AL-1.)*10. - FLOAT(J1-1)
      SIGT    = SSIG_TOT(J1,L)*(1.-T) + SSIG_TOT(J1+1,L)*T
      SIGINEL = SSIG(J1,L)*(1.-T) + SSIG(J1+1,L)*T
      SIGEL   = SIGT-SIGINEL
      SIGDIF(1,1) = SSIG_SD1LM(J1,L)*(1.-T) + SSIG_SD1LM(J1+1,L)*T
      SIGDIF(1,2) = SSIG_SD1HM(J1,L)*(1.-T) + SSIG_SD1HM(J1+1,L)*T
      SIGDIF(2,1) = SSIG_SD2LM(J1,L)*(1.-T) + SSIG_SD2LM(J1+1,L)*T
      SIGDIF(2,2) = SSIG_SD2HM(J1,L)*(1.-T) + SSIG_SD2HM(J1+1,L)*T
      SIGDIF(3,1) = SSIG_DDLM(J1,L)*(1.-T) + SSIG_DDLM(J1+1,L)*T
      SIGDIF(3,2) = SSIG_DDHM(J1,L)*(1.-T) + SSIG_DDHM(J1+1,L)*T
      SLOPE   = SSIG_B(J1,L) *(1.-T) + SSIG_B(J1+1,L)*T
      RHO     = SSIG_RHO(J1,L) *(1.-T) + SSIG_RHO(J1+1,L)*T

      END


      SUBROUTINE SIB_SIGMA_HAIR (L,SQS,SIGprod,SIGbdif) 
C-----------------------------------------------------------------------
C     Hadron-air cross sections, taken from interpolation table
C     calculated by SIBYLL_INI
C
C     input:       L     1      proton-air
C                        2      pi-air
C                        3      K-air
C                  SQS   sqrt(s)
C
C     output:      SIGprod    particle production cross section (mb)
C                  SIGbdif    q.ela and ela beam diff. cross section
C-----------------------------------------------------------------------
      SAVE
Cf2py intent(out) SIGprod, SIGbdif
      PARAMETER (NS_max = 20, NH_max = 80)
      COMMON /S_CCSIG/ SSIG(61,3), PJETC(0:NS_max,0:NH_max,61,2),
     &     SSIGN(61,3), SSIGNSD(61,3), ALINT(61,3),
     &     ASQSMIN, ASQSMAX, DASQS, NSQS
      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      DIMENSION LL(39)
      DATA LL /5*0,3*2,4*3,2*1,19*0,6*1/

      IF(NSQS.LE.0) THEN
        WRITE(LUN,'(//,1X,A)') 
     &    'SIB_SIGMA_HAIR: interpolation table not initialized.'
        STOP
      ENDIF
      IF(IABS(L).gt.39)THEN
         WRITE(LUN,'(//,1X,A)')     
     &        'SIB_SIGMA_HAIR: unknown beam particle!',L
         STOP
      ENDIF
      IF(L.GT.3) L=LL(IABS(L))
      IF(L.EQ.0)THEN
         WRITE(LUN,'(//,1X,A)')     
     &        'SIB_SIGMA_HAIR: unknown beam particle!',L
         STOP
      ENDIF

      AL = LOG10(SQS)
      J1 = (AL - 1.)*10. + 1
      if((j1.lt.1).or.(j1.ge.NSQS)) then
        write(LUN,'(1x,a,i3,1p,e12.3)') 
     &    'SIB_SIGMA_HAIR: energy out of range ',L,sqs
        J1 = min(J1,NSQS-1)
        J1 = max(J1,1)
      endif
      T = (AL-1.)*10. - FLOAT(J1-1)
      SIGprod = SSIGN(J1,L)*(1.-T) + SSIGN(J1+1,L)*T
      SIGbdif = SSIGNSD(J1,L)*(1.-T) + SSIGNSD(J1+1,L)*T

      END


      SUBROUTINE SIB_HADCSL(L,ECM,SIGTOT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO)
C-----------------------------------------------------------------------
C     low-energy cross section parametrizations (target always proton)
C
C     input:   L           beam particle: (1 - proton,
C                                          2 - pion,
C                                          3 - kaon)
C                          target is always proton
C              ECM         c.m. energy (GeV)
C
C     output:  SIGTOT      total cross section (mb)
C              SIGEL       elastic cross section (mb)
C              SIGDIF      diffractive cross section (sd-1,sd-2,dd, mb)
C              SLOPE       forward elastic slope (GeV**-2)
C              RHO         real/imaginary part of elastic amplitude
C-----------------------------------------------------------------------
      SAVE
      DIMENSION SIGDIF(3)

      COMMON /S_CFLAFR/ PAR(80), IPAR(80)

C  proton-proton cross section as reference
      CALL SIB_HADCS1(1,ECM,SIGTOT,SIGEL,SIGINEL,SLOPE,RHO)

C  parametrization for diffraction
      Xi_min = 1.5/(ECM*ECM)
      Xi_max = PAR(13)
      SIGeff = SIGEL
      call SIB_HADCS2(ECM,Xi_min,Xi_max,SIGeff,SIGDIF)

      if(L.eq.1) return

C  regge motivated rescaling of diffraction dissociation
      sigtot_pp = SIGTOT
      sigel_pp  = SIGEL
      slope_pp  = SLOPE
      CALL SIB_HADCS1(L,ECM,SIGTOT,SIGEL,SIGINEL,SLOPE,RHO)
      SIGDIF(1) = slope_pp/SLOPE*SIGTOT/sigtot_pp*SIGDIF(1)
      SIGDIF(2) = slope_pp/SLOPE*SIGEL/sigel_pp*SIGDIF(2)
      SIGDIF(3) = SIGTOT/sigtot_pp*SIGDIF(3)

      END


      SUBROUTINE SIB_HADCS1(L,ECM,SIGTOT,SIGEL,SIGINEL,SLOPE,RHO)
C-----------------------------------------------------------------------
C     low-energy cross section parametrizations
C
C     input:   L           beam particle: (1 - proton,
C                                          2 - pion,
C                                          3 - kaon)
C                          target is always proton
C              ECM         c.m. energy (GeV)
C
C     output:  SIGTOT      total cross section (mb)
C              SIGEL       elastic cross section (mb)
C              SIGDIF      diffractive cross section (sd-1,sd-2,dd, mb)
C              SLOPE       forward elastic slope (GeV**-2)
C              RHO         real/imaginary part of elastic amplitude
C
C     comments:
C     - low-energy data interpolation uses PDG fits from 1992
C     - slopes from ???, new fit to pp data
C     - high-energy extrapolation by Donnachie-Landshoff like fit made
C       by PDG 1996
C     - analytic extension of amplitude to calculate rho
C-----------------------------------------------------------------------
      SAVE

      DIMENSION TPDG92(7,2,6),TPDG96(9,6),BURQ83(3,6),XMA(6)

      DATA TPDG92  /
     &  3.D0, 2100.D0, 48.D0, 0.D0, 1.D0, 0.522D0, -4.51D0,
     &  3.D0, 2100.D0, 11.9D0, 26.9D0, -1.21D0, 0.169D0, -1.85D0,
     &  5.D0, 2100.D0, 38.4D0, 77.6D0, -0.64D0, 0.26D0, -1.2D0,
     &  5.D0, 2100.D0, 10.2D0, 52.7D0, -1.16D0, 0.125D0, -1.28D0,
     &  4.D0, 340.D0,  16.4D0, 19.3D0, -0.42D0, 0.19D0, 0.D0,
     &  4.D0, 340.D0,  0.D0, 11.4D0, -0.4D0, 0.079D0, 0.D0,
     &  2.5D0, 370.D0, 33.D0, 14.D0, -1.36D0, 0.456D0, -4.03D0,
     &  2.5D0, 370.D0, 1.76D0, 11.2D0, -0.64D0, 0.043D0, 0.D0,
     &  2.D0, 310.D0,  18.1D0, 0.D0, 1.D0, 0.26D0, -1.D0,
     &  2.D0, 310.D0,  5.D0, 8.1D0, -1.8D0, 0.16D0, -1.3D0,
     &  3.D0, 310.D0,  32.1D0, 0.D0, 1.D0, 0.66D0, -5.6D0,
     &  3.D0, 310.D0,  7.3D0, 0.D0, 1.D0, 0.29D0, -2.4D0  /

      DATA TPDG96  /
     &  50.D0, 22.D0,0.079D0,0.25D0,0.D0,
     &         77.15D0,-21.05D0,0.46D0,0.9D0,
     &  50.D0, 22.D0,0.079D0,0.25D0,0.D0,
     &         77.15D0,21.05D0,0.46D0,0.9D0,
     &  10.D0, 13.70,0.079D0,0.25D0,0.D0,
     &         31.85D0,-4.05D0,0.45D0,0.9D0,
     &  10.D0, 13.70,0.079D0,0.25D0,0.D0,
     &         31.85D0,4.05D0,0.45D0,0.9D0,
     &  10.D0, 12.20,0.079D0,0.25D0,0.D0,
     &         17.35D0,-9.05D0,0.50D0,0.9D0,
     &  10.D0, 12.20,0.079D0,0.25D0,0.D0,
     &         17.35D0,9.05D0,0.50D0,0.9D0  /

      DATA BURQ83 /
     &  8.557D0,  0.00D0, 0.574D0,
     &  11.13D0,  7.23D0, 0.30D0,
     &  9.11D0,  -0.73D0, 0.28D0,
     &  9.11D0,   0.65D0, 0.28D0,
     &  8.55D0,  -5.98D0, 0.28D0,
     &  8.55D0,   1.60D0, 0.28D0  /

      DATA XMA / 2*0.93956563, 2*0.13956995, 2*0.493677 /
      DATA GEV2MB /0.389365/
      DATA PI /3.14159265358979/

C  find index
      IF(L.eq.1) THEN
        K = 1                            ! p p
      ELSE IF(L.eq.2) THEN
        K = 3                            ! pi+ p
*       K = 4                            ! pi- p
      ELSE IF(L.eq.3) THEN
        K = 5                            ! K+ p
*       K = 6                            ! K- p
      ELSE
        GOTO 100
      ENDIF

C  calculate lab momentum
      SS = ECM**2
      E1 = (SS-XMA(1)**2-XMA(K)**2)/(2.*XMA(1))
      PL = SQRT((E1-XMA(K))*(E1+XMA(K)))
      PLL = LOG(PL)

C  check against lower limit
      IF(ECM.LE.XMA(1)+XMA(K)) GOTO 200

      XP  = TPDG96(2,K)*SS**TPDG96(3,K)
      YP  = TPDG96(6,K)/SS**TPDG96(8,K)
      YM  = TPDG96(7,K)/SS**TPDG96(8,K)

      PHR = TAN(PI/2.*(1.-TPDG96(8,K)))
      PHP = TAN(PI/2.*(1.+TPDG96(3,K)))
      RHO = (-YP/PHR + YM*PHR - XP/PHP)/(YP+YM+XP)

      SLOPE = BURQ83(1,K)+BURQ83(2,K)/SQRT(PL)+BURQ83(3,K)*PLL

C  select energy range and interpolation method
      IF(PL.LT.TPDG96(1,K)) THEN
        SIGTOT = TPDG92(3,1,K)+TPDG92(4,1,K)*PL**TPDG92(5,1,K)
     &          + TPDG92(6,1,K)*PLL**2+TPDG92(7,1,K)*PLL
        SIGEL  = TPDG92(3,2,K)+TPDG92(4,2,K)*PL**TPDG92(5,2,K)
     &          + TPDG92(6,2,K)*PLL**2+TPDG92(7,2,K)*PLL
      ELSE IF(PL.LT.TPDG92(2,1,K)) THEN
        SIGTO1 = TPDG92(3,1,K)+TPDG92(4,1,K)*PL**TPDG92(5,1,K)
     &          + TPDG92(6,1,K)*PLL**2+TPDG92(7,1,K)*PLL
        SIGEL1 = TPDG92(3,2,K)+TPDG92(4,2,K)*PL**TPDG92(5,2,K)
     &          + TPDG92(6,2,K)*PLL**2+TPDG92(7,2,K)*PLL
        SIGTO2 = YP+YM+XP
        SIGEL2 = SIGTO2**2/(16.*PI*SLOPE*GEV2MB)*(1.+RHO**2)
        X2 = LOG(PL/TPDG96(1,K))/LOG(TPDG92(2,1,K)/TPDG96(1,K))
        X1 = 1. - X2
        SIGTOT = SIGTO2*X2 + SIGTO1*X1
        SIGEL  = SIGEL2*X2 + SIGEL1*X1
      ELSE
        SIGTOT = YP+YM+XP
        SIGEL  = SIGTOT**2/(16.*PI*SLOPE*GEV2MB)*(1.+RHO**2)
      ENDIF
      SIGINEL = SIGTOT-SIGEL

      RETURN

 100  CONTINUE
        WRITE(6,'(1X,2A,2I7)') 'SIB_HADCSL: ',
     &    'invalid beam particle: ',L
        RETURN

 200  CONTINUE
        WRITE(6,'(1X,2A,1P,E12.4)') 'SIB_HADCSL: ',
     &    'energy too small (Ecm): ',ECM

      END


      SUBROUTINE SIB_HADCS2(SQS,Xi_min,Xi_max,SIGeff,SIGDIF)
C-----------------------------------------------------------------------
C   cross section for diffraction dissociation 
C
C   - single diffraction dissociation:
C     Goulianos' parametrization (Ref: PL B358 (1995) 379)
C   - double diffration dissociation: simple scaling model using 
C     single diff. cross section
C
C     in addition rescaling for different particles is applied using
C     internal rescaling tables (not implemented yet)
C
C     input:     SQS         c.m. energy (GeV)
C                Xi_min      min. diff mass (squared) = Xi_min*SQS**2
C                Xi_max      max. diff mass (squared) = Xi_max*SQS**2
C                SIGeff      effective cross section for DD scaling
C
C     output:    sig_sd1     cross section for diss. of particle 1 (mb)
C                sig_sd2     cross section for diss. of particle 2 (mb)
C                sig_dd      cross section for diss. of both particles
C-----------------------------------------------------------------------
      SAVE

      DIMENSION SIGDIF(3)
      DOUBLE PRECISION Xpos1(96),Xwgh1(96),Xpos2(96),Xwgh2(96)
      DOUBLE PRECISION xil,xiu,tl,tu

C  model parameters
      DATA delta    / 0.104 /
      DATA alphap   / 0.25 /
      DATA beta0    / 6.56 /
      DATA gpom0    / 1.21 /
      DATA xm_p     / 0.938 /
      DATA x_rad2   / 0.71 /

C  integration precision
      DATA Ngau1    / 32 /
      DATA Ngau2    / 32 /

      DATA PI /3.14159265358979/
      DATA GEV2MB /0.389365/


      SIGDIF(1) = 0.
      SIGDIF(2) = 0.
      SIGDIF(3) = 0.

      XIL = LOG(Xi_min)
      XIU = LOG(Xi_max)

      if(XIL.ge.XIU) return

      SS = SQS*SQS
      xm4_p2 = 4.*xm_p**2
      fac = beta0**2/(16.*PI)

      t1 = -5.
      t2 = 0.
      tl = x_rad2/3./(1.-t1/x_rad2)**3
      tu = x_rad2/3./(1.-t2/x_rad2)**3

C  flux renormalization and cross section for pp/ppbar case

      Xnorm  = 0.

      xil = log(1.5/SS)
      xiu = log(0.1)

      IF(xiu.LE.xil) goto 1000

      CALL SIB_GAUSET(xil,xiu,Ngau1,xpos1,xwgh1)
      CALL SIB_GAUSET(tl,tu,Ngau2,xpos2,xwgh2)

      do i1=1,Ngau1

        xi = exp(xpos1(i1))
        w_xi = Xwgh1(i1)

        do i2=1,Ngau2

          tt = x_rad2-x_rad2*(x_rad2/(3.*xpos2(i2)))**(1./3.)

          alpha_t =  1.+delta+alphap*tt
          f2_t = ((xm4_p2-2.8*tt)/(xm4_p2-tt))**2
            
          Xnorm = Xnorm
     &      + f2_t*xi**(2.-2.*alpha_t)*Xwgh2(i2)*w_xi

        enddo
      enddo   

      Xnorm = Xnorm*fac

 1000 continue

      XIL = LOG(Xi_min)
      XIU = LOG(Xi_max)

      T1 = -5.
      T2 = 0.

      TL = x_rad2/3./(1.-t1/x_rad2)**3
      TU = x_rad2/3./(1.-t2/x_rad2)**3

C  single diffraction diss. cross section 

      CSdiff = 0.

      CALL SIB_GAUSET(XIL,XIU,NGAU1,XPOS1,XWGH1)
      CALL SIB_GAUSET(TL,TU,NGAU2,XPOS2,XWGH2)

      do i1=1,Ngau1

        xi = exp(xpos1(i1))
        w_xi = Xwgh1(i1)*beta0*gpom0*(xi*ss)**delta

        do i2=1,Ngau2

          tt = x_rad2-x_rad2*(x_rad2/(3.*xpos2(i2)))**(1./3.)

          alpha_t =  1.+delta+alphap*tt
          f2_t = ((xm4_p2-2.8*tt)/(xm4_p2-tt))**2

          CSdiff = CSdiff 
     &      + f2_t*xi**(2.-2.*alpha_t)*Xwgh2(i2)*w_xi

        enddo
      enddo

      CSdiff = CSdiff*fac*GEV2MB/MAX(1.,Xnorm)

*     write(6,'(1x,1p,4e14.3)') 
*    &  sqrt(SS),Xnorm,2.*CSdiff*MAX(1.,Xnorm),2.*CSdiff

      SIGDIF(1) = CSdiff
      SIGDIF(2) = CSdiff

C  double diff. dissociation from simple probability consideration
*     Pdiff = 0.5-sqrt(0.25-CSdiff/SIGeff)
      Pdiff = CSdiff/SIGeff
      SIGDIF(3) = Pdiff*Pdiff*SIGeff

      END


      SUBROUTINE DECSIB 
C-----------------------------------------------------------------------
C...Decay all unstable particle in Sibyll
C.  decayed particle have the code increased by 10000
C
C   changed to allow for multiple calls to DECSIB in one event
C-----------------------------------------------------------------------
      SAVE
      COMMON /S_CSYDEC/ CBR(223), KDEC(1338), LBARP(99), IDB(99)
      COMMON /S_PLIST/ P(8000,5), LLIST(8000), NP
      COMMON /S_PLIST1/ LLIST1(8000)
      COMMON /S_MASS1/ AM(99), AM2(99)
      DIMENSION P0(5), LL(10), PD(10,5)
      COMMON /S_PARTO/ NFORIG(8000), NPORIG(8000), IPFLAG, NINT
      NN = 1
      DO J=1,NP
         LLIST1(J) = 0
      ENDDO
      DO WHILE (NN .LE. NP)
         L= LLIST(NN)
         LA = IABS(L)
         if(LA.lt.100) then
           IF (IDB(LA) .GT. 0)  THEN
              DO K=1,5
                P0(K) = P(NN,K)
              ENDDO
              CALL DECPAR (L,P0,ND,LL,PD)
              LLIST(NN) = LLIST(NN)+ISIGN(10000,LLIST(NN))
              DO J=1,ND
                NP = NP+1
                if(NP.gt.8000) then
                  write(6,'(1x,a,2i8)') 
     &              'DECSIB: no space left in S_PLIST (NP,ND):',NP,ND
                  NP = NP-1
                  return
                endif
                DO K=1,5
                  P(NP,K) = PD(J,K)
                ENDDO
                LLIST(NP)=LL(J)
                LLIST1(NP)=NN
                NPORIG(NP)= NPORIG(NN)
                NFORIG(NP) = L
              ENDDO
           ENDIF
         endif
         NN = NN+1
      ENDDO

c	call sib_list(20)

      END



      SUBROUTINE DECPAR (LA,P0,ND,LL,P)
C-----------------------------------------------------------------------
C...This subroutine generates the decay of a particle
C.  with ID = LA, and 5-momentum P0(1:5)
C.  into ND particles of 5-momenta P(j,1:5) (j=1:ND)
C.
C.  If the initial particle code is LA=0
C.  then ND and LL(1:ND) are considered as  input and
C.  the routine generates a phase space decay into ND
C.  particles of codes LL(1:nd)
C.
C.  june 1992
C.  This version  contains the decay of polarized muons
C.  The muon codes are  L =  4 : mu+ R
C.                          -4 : mu+ L
C.                           5 : mu- L
C.                          -5 : mu- R
C-----------------------------------------------------------------------
      SAVE
      COMMON /S_CSYDEC/ CBR(223), KDEC(1338), LBARP(99), IDB(99)
      COMMON /S_MASS1/ AM(99), AM2(99)
      DIMENSION P0(5), LL(10), P(10,5)
      DIMENSION PV(10,5), RORD(10), UE(3),BE(3), FACN(3:10)
      DATA FACN /2.,5.,15.,60.,250.,1500.,12000.,120000./
      DATA PI /3.1415926/
      
C...c.m.s. Momentum in two particle decays
      PAWT(A,B,C) = SQRT((A**2-(B+C)**2+1.e-5)*(A**2-(B-C)**2))/(2.*A)

C...Phase space decay into the particles in the list
      IF (LA .EQ. 0)  THEN
          MAT = 0
          MBST = 0
          PS = 0.
          DO J=1,ND
CDH          following statements corrected by D.H. dec 20.,1995
             P (J,5) = AM(IABS(LL(J)))
             PV(J,5) = AM(IABS(LL(J)))
             PS = PS+P(J,5)
          ENDDO
          DO J=1,4
             PV(1,J) = P0(J)
          ENDDO
          PV(1,5) = P0(5)
          GOTO 140
      ENDIF
         
C...Choose decay channel
      L = IABS(LA)
      ND=0
      IDC = IDB(L)-1
      IF (IDC+1 .LE.0)  RETURN
      RBR = S_RNDM(0)
110   IDC=IDC+1
      IF(RBR.GT.CBR(IDC))  GOTO 110

      KD =6*(IDC-1)+1
      ND = KDEC(KD)
      MAT= KDEC(KD+1)
      MBST=0
      IF (MAT .GT.0 .AND. P0(4) .GT. 20*P0(5)) MBST=1
      IF (MAT .GT.0 .AND. MBST .EQ. 0) 
     +        BETA = SQRT(P0(1)**2+P0(2)**2+P0(3)**2)/P0(4)
      PS = 0.
      DO J=1,ND
         LL(J) = KDEC(KD+1+J)
         P(J,5)  = AM(LL(J))
         PV(J,5) = AM(LL(J))
         PS = PS + P(J,5)
      ENDDO
      DO J=1,4
         PV(1,J) = 0.
         IF (MBST .EQ. 0)  PV(1,J) = P0(J)
      ENDDO
      IF (MBST .EQ. 1)  PV(1,4) = P0(5)
      PV(1,5) = P0(5)

140   IF (ND .EQ. 2) GOTO 280

      IF (ND .EQ. 1)  THEN
         DO J=1,4
            P(1,J) = P0(J)
         ENDDO
         RETURN
      ENDIF

C...Calculate maximum weight for ND-particle decay
      WWTMAX = 1./FACN(ND)      
      PMAX=PV(1,5)-PS+P(ND,5)
      PMIN=0.
      DO IL=ND-1,1,-1
         PMAX = PMAX+P(IL,5)
         PMIN = PMIN+P(IL+1,5)
         WWTMAX = WWTMAX*PAWT(PMAX,PMIN,P(IL,5))
      ENDDO

C...generation of the masses, compute weight, if rejected try again
240   RORD(1) = 1.
      DO 260 IL1=2,ND-1
      RSAV = S_RNDM(0)
      DO 250 IL2=IL1-1,1,-1
      IF(RSAV.LE.RORD(IL2))   GOTO 260
250     RORD(IL2+1)=RORD(IL2)
260     RORD(IL2+1)=RSAV
      RORD(ND) = 0.
      WT = 1.      
      DO 270 IL=ND-1,1,-1
      PV(IL,5)=PV(IL+1,5)+P(IL,5)+(RORD(IL)-RORD(IL+1))*(PV(1,5)-PS)
270   WT=WT*PAWT(PV(IL,5),PV(IL+1,5),P(IL,5))
      IF (WT.LT.S_RNDM(0)*WWTMAX)   GOTO 240

C...Perform two particle decays in respective cm frame
280   DO 300 IL=1,ND-1
      PA=PAWT(PV(IL,5),PV(IL+1,5),P(IL,5))
      UE(3)=2.*S_RNDM(0)-1.
      PHI=2.*PI*S_RNDM(0)
      UT = SQRT(1.-UE(3)**2)
      UE(1) = UT*COS(PHI)
      UE(2) = UT*SIN(PHI)
      DO 290 J=1,3
      P(IL,J)=PA*UE(J)
290   PV(IL+1,J)=-PA*UE(J)
      P(IL,4)=SQRT(PA**2+P(IL,5)**2)
300   PV(IL+1,4)=SQRT(PA**2+PV(IL+1,5)**2)

C...Lorentz transform decay products to lab frame
      DO 310 J=1,4
310   P(ND,J)=PV(ND,J)
      DO 340 IL=ND-1,1,-1
      DO 320 J=1,3
320   BE(J)=PV(IL,J)/PV(IL,4)
      GA=PV(IL,4)/PV(IL,5)
      DO 340 I=IL,ND
      BEP = BE(1)*P(I,1)+BE(2)*P(I,2)+BE(3)*P(I,3)
      DO 330 J=1,3
330   P(I,J)=P(I,J)+GA*(GA*BEP/(1.+GA)+P(I,4))*BE(J)
340   P(I,4)=GA*(P(I,4)+BEP)
      
C...Weak decays
      IF (MAT .EQ. 1)  THEN
         F1=P(2,4)*P(3,4)-P(2,1)*P(3,1)-P(2,2)*P(3,2)-P(2,3)*P(3,3)      
         IF (MBST.EQ.1)  THEN
C        WT = P0(5)*P(1,4)*F1
            WT = P0(5)*(P(1,4)+FLOAT(LA/L)*P(1,3))*F1
         ENDIF
         IF (MBST.EQ.0)  THEN  
            WT=F1*(P(1,4)*P0(4)-P(1,1)*P0(1)-P(1,2)*P0(2)-P(1,3)*P0(3))
            IF(L.lt.50)
     +           WT= WT-FLOAT(LA/L)*(P0(4)*BETA*P(1,4)-P0(4)*P(1,3))*F1
         ENDIF
         WTMAX = P0(5)**4/8.
         IF(WT.LT.S_RNDM(0)*WTMAX)   GOTO 240
      ENDIF

C...Boost back for rapidly moving particle
      IF (MBST .EQ. 1)   THEN
         DO 440 J=1,3
440      BE(J)=P0(J)/P0(4)
         GA= P0(4)/P0(5)
         DO 460 I=1,ND
         BEP=BE(1)*P(I,1)+BE(2)*P(I,2)+BE(3)*P(I,3)
         DO 450 J=1,3
450         P(I,J)=P(I,J)+GA*(GA*BEP/(1.+GA)+P(I,4))*BE(J)
460         P(I,4)=GA*(P(I,4)+BEP)
      ENDIF

C...labels for antiparticle decay
      IF (LA .LT. 0 .AND. L .GT. 18)  THEN
           DO J=1,ND
            LL(J) = LBARP(LL(J))
         ENDDO
      ENDIF

      RETURN
      END


      BLOCK DATA DATDEC
C-----------------------------------------------------------------------
C...initialization of SIBYLL particle data
C-----------------------------------------------------------------------
      SAVE

      COMMON /S_CSYDEC/ CBR(223), KDEC(1338), LBARP(99), IDB(99)
      COMMON /S_MASS1/ AM(99), AM2(99)
      COMMON /S_CHP/ ICHP(99), ISTR(99), IBAR(99)
      COMMON /S_CHM/ ICHM(99)
      COMMON /S_CNAM/ NAMP (0:99)
      CHARACTER NAMP*6
      DATA CBR /3*1.,0.,1.,1.,0.6354,0.8422,0.8981,0.9157,0.9492,1.,
     +   0.6354,0.8422,0.8981,0.9157,0.9492,1.,0.1965,0.3224,0.4579,
     +   0.5934,0.7967,1.,0.6925,1.,3*0.,0.5,1.,0.5,1.,
     +   0.3941,0.7197,0.9470,0.9930,1.,0.,0.4460,0.6530,0.9470,0.9770,
     +   0.9980,4*1.,0.6670,1.,9*0.,0.6670,1.,0.6670,1.,0.6670,1.,
     +   0.8940,0.9830,1.,0.4930,0.8340,0.9870,1.,0.5160,5*1.,0.6410,1.,
     +   1.,0.67,1.,0.33,1.,1.,0.88,0.94,1.,0.88,0.94,1.,0.88,0.94,1.,
     +   0.33,1.,0.67,1.,0.678,0.914,1., 0.217,0.398,0.506,0.595,0.684,
     +	 0.768,0.852,0.923,0.976,1., 0.217,0.398,0.506,0.595,0.684,
     +	 0.768,0.852,0.923,0.976,1., 0.2490,0.4604,0.5338,0.5703,0.7440,
     +	 0.7840,0.8460,0.8880,0.9230,0.9650,1., 0.2490,0.4604,0.5338,
     +	 0.5703,0.7440,0.7840,0.8460,0.8880,0.9230,0.9650,1.,
     +	 0.1666,0.3332,0.4998,0.6664,0.8330,1., 0.6770,0.9840,1.,
     +	 0.6770,0.9840,1., 0.6190,1., 0.6190,1., 0.0602, 0.1203, 1.,
     +	 3*1., 0.06,0.08,0.14,0.16,0.73,0.855,0.98,1., 0.08,0.16,0.92,1,
     +	 0.2335,0.4283,0.6446,0.7099,0.8080,0.9080,0.9380,0.9540,
     +	 0.9840,1., 3*1., 0.5,1., 0.5,1., 0.08,0.16,0.92,1. ,0.942,1.,
     +   0.942,1., 0.2493,0.4061,0.5602,0.6860,0.7608,0.8305,0.8818,
     +   0.9277,0.9691,1., 0.2493,0.4061,0.5602,0.6860,0.7608,0.8305,
     +   0.8818,0.9277,0.9691,1./
      DATA AM / 0.,2*0.511E-3, 2*0.10566, 0.13497, 2*0.13957,
     +   2*0.49368, 2*0.49761, 0.93827, 0.93957, 4*0.,0.93827,
     +   0.93957, 2*0.49761, 0.54785,0.95766,2*0.76690,0.76850,
     +   2*0.89166,2*0.89600,0.78265,1.01946,1.18937,1.19264,
     +   1.19745,1.31486,1.32171,1.11568,1.23100,1.23500,
     +   1.23400,1.23300,1.38280,1.38370,1.38720,
     +   1.53180,1.53500,1.67245, 9*0., 2*1.86926, 10*0., 2*1.86484,
     +	 2.9803, 2*1.9685, 2*2.1123, 2*2.01027, 2*2.00697, 0, 3.09692,
     +	 2.45402, 2.4529, 2.45376, 2.4679, 2.4710, 2.28646, 4*0.,2.5184,
     +	 2.5175, 2.5180, 2.6466, 2.6461, 2.6975 /
      DATA AM2 /0.,2*2.61121E-07,2*0.011164,0.018217,0.019480,
     + 0.019480,0.243720,0.243720,0.247616,0.247616,0.880351,0.882792,
     + 0.000000,0.000000,0.000000,0.000000,0.880351,0.882792,0.247616,
     + 0.247616,0.300140,0.917113,0.588136,0.588136,0.590592,0.795058,
     + 0.795058,0.802816,0.802816,0.612541,1.039299,1.414601,1.422390,
     + 1.433887,1.728857,1.746917,1.244742,1.515361,1.525225,1.522765,
     + 1.520289,1.912136,1.914626,1.924324,2.346411,2.356225,2.797022,
     + 9*0., 2*3.49414, 10*0., 2*3.477628, 8.882188, 2*3.8750,2*4.4618,
     + 2*4.041186,2*4.027928, 0, 9.590914, 6.022214, 6.016718, 6.020938,
     + 6.09053, 6.105841, 5.227899, 4*0., 6.342339, 6.337806, 6.340323,
     + 7.004492, 7.001845, 7.276506 /
      DATA IDB /
     +    0,0,0,1,2,3,5,6,7,13,19,25,8*0,30,32,34,40,46,47,48,49,60,62,
     +    64,66,69,73,75,76,77,78,79,81,82,84,86,87,90,93,96,98,100,
     +	  9*0,103,113,10*0,123,134,145,204,214,200,202,151,154,157,159,0,
     +	  161,164,165,166,167,175,179,4*0,189,190,191,192,194,196 /
      DATA KDEC /
     + 3,1,15,2,18,0,3,1,16,3,17,0,2,0,1,1,8*0,2,0,4,17,0,0,2,0,5,18,0,
     + 0,2,0,4,17,0,0,2,0,7,6,0,0,3,0,7,7,8,0,3,0,7,6,6,0,3,1,17,4,6,0,
     + 3,1,15,2,6,0,2,0,5,18,0,0,2,0,8,6,0,0,3,0,8,8,7,0,3,0,8,6,6,0,3,
     + 1,18,5,6,0,3,1,16,3,6,0,3,0,6,6,6,0,3,0,7,8,6,0,3,1,18,5,7,0,3,
     + 1,17,4,8,0,3,1,16,3,7,0,3,1,15,2,8,0,2,0,7,8,0,0,2,0,6,6,20*0,1,
     + 0,11,3*0,1,0,12,0,0,0,1,0,11,0,0,0,1,0,12,0,0,0,2,0,1,1,0,0,3,0,
     + 6,6,6,0,3,0,7,8,6,0,3,0,1,7,8,0,3,0,1,3,2,7*0,3,0,7,8,23,0,3,0,6
     + ,6,23,0,2,0,1,27,0,0,2,0,1,32,0,0,2,0,1,1,0,0,3,0,6,6,6,0,2,0,7,
     + 6,0,0,2,0,8,6,0,0,2,0,7,8,0,0,2,0,21,7,0,0,2,0,9,6,0,0,54*0,2,0,
     + 22,8,0,0,2,0,10,6,0,0,2,0,9,8,0,0,2,0,21,6,0,0,2,0,10,7,0,0,
     + 2,0,22,6,0,0,3,0,7,8,6,0,2,0,1,6,0,0,2,0,7,8,0,0,2,0,9,10,0,
     + 0,2,0,11,12,0,0,3,0,7,
     + 8,6,0,2,0,1,23,0,0,2,0,13,6,0,0,2,0,14,7,0,0,2,0,39,1,0,0,2,
     + 0,14,8,0,0,2,0,39,6,0,0,2,0,39,8,0,0,2,0,13,8,0,0,2,0,
     + 14,6,0,0,2,0,13,7,0,0,2,0,13,6,
     + 0,0,2,0,14,7,0,0,2,0,13,8,0,0,2,0,14,6,0,0,2,0,14,8,0,0,2,0,
     + 39,7,0,0,2,0,34,6,0,0,2,0,35,7,0,0,2,0,39,6,0,0,2,0,34,8,0,0,
     + 2,0,36,7,0,0,2,0,39,8,0,0,2,
     + 0,35,8,0,0,2,0,36,6,0,0,2,0,37,6,0,0,2,0,38,7,0,0,2,0,
     + 37,8,0,0,2,0,38,6,0,0,2,0,39,10,0,0,2,0,37,8,0,0,2,0,38,6,0,0,
     + 3,0,22,7,6,0,3,0,22,9,22,0,2,0,22,7,0,0,3,1,2,15,22,0,3,1,4,17,
     + 22,0,3,1,2,15,31,0,3,1,4,17,31,0,2,0,31,25,0,0,3,0,33,7,6,0,
     + 3,0,10,7,7,0,
     + 3,0,21,8,6,0,3,0,21,10,21,0,2,0,21,8,0,0,3,1,3,16,21,0,3,1,5,18,
     + 21,0,3,1,3,16,30,0,3,1,5,18,30,0,2,0,30,26,0,0,3,0,33,8,6,0,
     + 3,0,9,8,8,0,
     + 2,0,29,7,0,0,2,0,31,6,0,0,2,0,22,6,0,0,2,0,10,7,0,0,2,0,31,27,0,
     + 0,2,0,30,27,0,0,2,0,29,25,0,0,3,1,2,15,10,0,3,1,2,15,29,0,
     + 3,1,4,17,10,0,3,1,4,17,29,0,
     + 2,0,28,8,0,0,2,0,30,6,0,0,2,0,21,6,0,0,2,0,9,8,0,0,2,0,30,27,0,
     + 0,2,0,31,27,0,0,2,0,28,26,0,0,3,1,3,16,9,0,3,1,3,16,28,0,
     + 3,1,5,18,9,0,3,1,5,18,28,0,
     + 3,0,6,21,22,0,3,0,6,9,10,0,3,0,23,6,6,0,3,0,23,7,8,0,3,0,24,6,6,
     + 0,3,0,24,7,8,0,
     + 2,0,71,7,0,0,2,0,59,6,0,0,2,0,59,1,0,0,
     + 2,0,72,8,0,0,2,0,60,6,0,0,2,0,60,1,0,0,
     + 2,0,71,6,0,0,2,0,71,1,0,0,2,0,72,6,0,0,2,0,72,1,0,0,
     + 2,0,2,3,0,0,2,0,4,5,0,0,3,0,6,7,8,0,
     + 2,0,89,7,0,0,2,0,89,6,0,0,2,0,89,8,0,0,
     + 3,1,2,15,22,0,3,1,2,15,33,0,3,1,4,17,22,0,3,1,4,17,33,0,2,0,7,22,
     + 0,0,2,0,9,22,0,0,2,0,7,33,0,0,2,0,9,33,0,0,
     + 3,1,2,15,10,0,3,1,4,17,10,0,2,0,7,10,0,0,2,0,9,10,0,0,
     + 3,0,7,10,13,0,3,0,7,22,14,0,3,0,7,8,13,0,3,0,9,10,13,0,3,0,9,22,
     + 14,0,3,0,22,8,40,0,3,1,2,15,39,0,3,1,2,15,14,0,3,1,4,17,39,0,3,
     + 1,4,17,14,0,
     + 2,0,89,7,0,0,2,0,89,6,0,0,2,0,89,8,0,0,
     + 2,0,87,6,0,0,2,0,87,1,0,0,2,0,88,6,0,0,2,0,88,1,0,0,
     + 3,1,2,15,10,0,3,1,4,17,10,0,2,0,7,10,0,0,2,0,9,10,0,0 ,
     + 2,0,74,1,0,0 ,2,0,74,6,0,0 , 2,0,75,1,0,0 ,2,0,75,6,0,0, 
     + 2,0,23,25,0,0, 4,0,9,10,7,6, 3,0,9,10,7,0, 2,0,33,7,0,0, 
     + 3,1,23,2,15,0, 3,1,33,2,15,0, 2,0,23,7,0,0, 4,0,12,10,7,7,
     + 2,0,9,12,0,0, 4,0,7,8,7,8, 2,0,23,26,0,0, 4,0,10,9,8,6, 
     + 3,0,10,9,8,0, 2,0,33,8,0,0, 3,1,23,3,16,0, 3,1,33,3,16,0,
     + 2,0,23,8,0,0, 4,0,12,9,8,8, 2,0,10,12,0,0, 4,0,7,8,7,8/
      DATA LBARP/1,3,2,5,4,6,8,7,10,9,11,12,-13,-14,16,15,18,17,13,14,
     +  22,21,23,24,26,25,27,29,28,31,30,32,33,-34,-35,-36,-37,-38,-39,
     +  -40,-41,-42,-43,-44,-45,-46,-47,-48,-49,9*0,60,59,10*0,72,71,
     +	73,75,74,77,76,79,78,81,80,0,83,-84,-85,-86,-87,-88,-89,4*0,-94,
     +  -95,-96,-97,-98,-99 /
      DATA ICHP /0,1,-1,1,-1,0,1,-1,1,-1,0,0,1,0,4*0,-1,0,4*0,
     +    1,-1,0,1,-1,4*0,1,0,-1,0,-1,0,2,1,0,-1,1,0,-1,0,-1,-1,
     +	  9*0,1,-1,10*0,0,0,0,1,-1,1,-1,1,-1,0,0,0, ! charmed mesons
     + 0,2,1,0,1,0,1,4*0,2,1,0,1,0,0 / ! charmed baryons

      DATA ISTR /8*0,-1,+1,-1,-1,8*0,-1,+1,5*0,-1,+1,-1,+1,2*0, ! mesons
     +           3*1,2*2,1,4*0,3*1,2*2,3, ! baryons
     +           9*0,2*0,10*0,2*0,0,-1,1,-1,1,2*0,2*0,0,0, ! charmed mesons
     +		 3*0,2*1,0,4*0,3*0,2*1,2 / ! charmed baryons
      DATA IBAR /12*0,2*1,4*0,2*-1,13*0,16*1,9*0,2*0,10*0,2*0,0,4*0,
     + 		 2*0,2*0,0,0,6*1,4*0,6*1 /
      DATA ICHM /58*0,1,-1,10*0,1,-1,0,1,-1,+1,-1,1,-1,1,-1,0,0,
     +     6*1,4*0,6*1/
      DATA NAMP /
     +     '     ','gam   ','e+','e-','mu+','mu-','pi0',
     +     'pi+','pi-','k+', 'k-', 'k0l','k0s',
     +     'p', 'n', 'nue', 'nueb', 'num', 'numb', 'pbar', 'nbar',
     +     'k0', 'k0b', 'eta', 'etap', 'rho+', 'rho-','rho0',
     +     'k*+','k*-','k*0','k*0b','omeg', 'phi', 'SIG+', 'SIG0',
     +     'SIG-','XI0','XI-','LAM','DELT++','DELT+','DELT0','DELT-',
     +     'SIG*+','SIG*0','SIG*-', 'XI*0', 'XI*-', 'OME-',
     +	    9*'     ', 'D+', 'D-',10*'     ', 'D0', 'D0b', 'eta_c', 
     +	    'D_s+','D_s-','D*_s+','D*_s-','D*+', 'D*-', 'D*0', 'D*0b',
     +     '     ', 'J/psi',
     +	    'SIGc++', 'SIGc+', 'SIGc0','XI_c+','XI_c0','LAM_c+',
     +	    4*'     ', 'SIc*++','SIGc*0','SIGc*-', 'XI_c*+', 'XI_c*0',
     +	    'OME_c0'  /
      END
C->
      SUBROUTINE DECPR (LUN)
C...Print on unit LUN the list of particles and decay channels
      SAVE
      COMMON /S_CSYDEC/ CBR(223), KDEC(1338), LBARP(99), IDB(99)
      COMMON /S_MASS1/ AM(99), AM2(99)
      COMMON /S_CNAM/ NAMP (0:99)
      CHARACTER NAMP*6
      DIMENSION LL(4)
      
      DO L=1,99
         IDC = IDB(L)-1
         NC = 0
         WRITE (LUN,10) L,NAMP(L), AM(L)
         IF(IDC+1 .GT. 0)  THEN
            CB = 0.
110         IDC=IDC+1
            NC = NC+1
            CBOLD = CB
            CB = CBR(IDC)
            BR = CB-CBOLD
            KD = 6*(IDC-1)+1
            ND = KDEC(KD)
            MAT= KDEC(KD+1)
            DO J=1,ND
              LL(J) = KDEC(KD+1+J)
            ENDDO
            WRITE (LUN,15) NC,BR,ND,MAT, (NAMP(LL(J)),J=1,ND)
            IF (CB .LT. 1.)  GOTO 110
         ENDIF
      ENDDO
      RETURN
10    FORMAT(1X,I3,2X,A6,3X,F10.4)
15    FORMAT(5X,I2,2X,F9.4,I4,I4,2X,3(A6,2X))
      END



      SUBROUTINE DEC_DEBUG (L,P0, ND, LL, PD)
      SAVE
      COMMON /S_CNAM/ NAMP (0:99)
      CHARACTER*6 NAMP
      DIMENSION P0(5), LL(10), PD(10,5)
      ETOT = 0.
      DO J=1,ND
         ETOT = ETOT + PD(J,4)
      ENDDO
      WRITE(*,*)  NAMP(IABS(L)),' -> ', (NAMP(IABS(LL(J))),J=1,ND)
      WRITE(*,*)  ' Ei, Ef = ', P0(4), ETOT, ' L = ', L
      RETURN
      END



      SUBROUTINE SIB_SIG(Jint,SIB_SQS,SIB_PTmin,SIB_SIG_tot,
     &                 SIB_SIG_ine,SIB_diff,SIB_diff2,SIB_B_el,SIB_PJET)
C-----------------------------------------------------------------------
C
C...SIBYLL 2.1 cross sections 
C
C   input parameter: SIB_SQS   c.m.s. energy (GeV)
C                    Jint      1 p-p cross sections
C                              2 pi-p cross sections
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE

      PARAMETER (NS_max = 20, NH_max = 80)
      REAL SIB_PJET(0:NS_max,0:NH_max)
      REAL SIB_SQS,SIB_PTmin,
     &     SIB_SIG_ine,SIB_SIG_tot,SIB_diff(3),SIB_diff2(3,2),SIB_B_el


      COMMON /SIGMAS/SQS,SIGTOT,SIGEL,SIGINE,
     &               SIGSD1(2),SIGSD2(2),SIGDD(2),
     &               SLOPE,SLOPEc,RHO,PROB(0:NS_max,0:NH_max),SIGSUM


      COMMON /PROFILE/XNUS2,XMUS2,XNUSPI2,
     &                XNUH2,XMUH2,XNUHPI2,
     &                ENHPP,ENHPIP,al1,be1,al2,be2

      COMMON /S_CHDCNV/ABR(2,400),ABP(2,400),ABH(2,400),DB,NB

      PARAMETER ( NPARFIT = 22 )
      DIMENSION XI(50)
      COMMON /XSCTN_FIT/PARS(50,2)
      DOUBLE PRECISION PARS

      DIMENSION SIG_BRN(3)
      DIMENSION SIG_dif_1(2),SIG_dif_2(2),SIG_dd(2)

      DIMENSION IHAR(2)
      COMMON /QCD_XSCTN/SIGQCD(61,2),INIT
      DOUBLE PRECISION SIGQCD
      DATA (SIGQCD(K,1),K=    1,   61) /
     &8.4663E-02,1.8246E-01,3.3880E-01,5.6845E-01,8.8686E-01,1.3116E+00,
     &1.8626E+00,2.5645E+00,3.4445E+00,4.5343E+00,5.8715E+00,7.4962E+00,
     &9.4579E+00,1.1811E+01,1.4620E+01,1.7955E+01,2.1890E+01,2.6522E+01,
     &3.1952E+01,3.8303E+01,4.5704E+01,5.4307E+01,6.4284E+01,7.5818E+01,
     &8.9121E+01,1.0447E+02,1.2213E+02,1.4240E+02,1.6562E+02,1.9221E+02,
     &2.2260E+02,2.5733E+02,2.9694E+02,3.4207E+02,3.9348E+02,4.5194E+02,
     &5.1838E+02,5.9376E+02,6.7921E+02,7.7609E+02,8.8578E+02,1.0099E+03,
     &1.1504E+03,1.3090E+03,1.4882E+03,1.6903E+03,1.9183E+03,2.1754E+03,
     &2.4650E+03,2.7912E+03,3.1582E+03,3.5707E+03,4.0341E+03,4.5538E+03,
     &5.1360E+03,5.7883E+03,6.5193E+03,7.3358E+03,8.2428E+03,9.2498E+03,
     &1.0369E+04/
      DATA (SIGQCD(K,2),K=    1,   61) /
     &1.5665E-01,2.8800E-01,4.7863E-01,7.4235E-01,1.0949E+00,1.5547E+00,
     &2.1433E+00,2.8859E+00,3.8118E+00,4.9547E+00,6.3534E+00,8.0525E+00,
     &1.0103E+01,1.2563E+01,1.5498E+01,1.8986E+01,2.3111E+01,2.7971E+01,
     &3.3678E+01,4.0358E+01,4.8154E+01,5.7228E+01,6.7762E+01,7.9965E+01,
     &9.4071E+01,1.1034E+02,1.2909E+02,1.5063E+02,1.7536E+02,2.0370E+02,
     &2.3613E+02,2.7321E+02,3.1553E+02,3.6379E+02,4.1875E+02,4.8129E+02,
     &5.5238E+02,6.3311E+02,7.2470E+02,8.2854E+02,9.4614E+02,1.0792E+03,
     &1.2298E+03,1.3999E+03,1.5920E+03,1.8089E+03,2.0534E+03,2.3291E+03,
     &2.6396E+03,2.9892E+03,3.3825E+03,3.8248E+03,4.3219E+03,4.8803E+03,
     &5.5072E+03,6.2109E+03,7.0001E+03,7.8849E+03,8.8764E+03,9.9871E+03,
     &1.1231E+04/

      DATA CMBARN /0.389385/
      DATA PI /3.1415926/
      DATA INIT /0/

      IF(INIT.EQ.0) THEN
*        CALL HAR_INI
        CALL FACT_INI
        IHAR(1) = 0
        IHAR(2) = 0
        INIT = 1
      ENDIF

      ECM = SIB_SQS

      IF(JINT.EQ.1) THEN
c     K = 1 , proton
         DO K=1,NPARFIT
            XI(K) = PARS(K,1)
         ENDDO

      ELSE IF(JINT.EQ.2) THEN
c     K = 2 , pion
         DO K=1,NPARFIT
            XI(K) = PARS(K,2)
         ENDDO

      ENDIF

      XNUS2   = XI(12)
      XMUS2   = XI(13)
      XNUSPI2 = XI(14)

      XNUH2   = XI(15)
      XMUH2   = XI(16)
      XNUHPI2 = XI(17)

      CALL HAD_CONV(IABS(JINT))

      PTCUT = XI(10)+XI(21)*EXP(XI(22)*SQRT(2.D0*LOG(ECM)))
      INDX = abs(JINT)
      IHAR(INDX) = IHAR(INDX)+1
      SIGHAR = SIGQCD(IHAR(INDX),INDX)

      S = ECM**2

      BREG =  ABS(XI(18)) + XI(19)*LOG(S)
      BPOM =  ABS(XI(12)) + XI(13)*LOG(S)
      IK = ABS(JINT)
      DO JB=1,NB
        B = DB*FLOAT(JB-1)
        ABR(IK,JB) = 2./(8.*PI*BREG)*EXP(-B**2/(4.*BREG))
        ABP(IK,JB) = 2./(8.*PI*BPOM)*EXP(-B**2/(4.*BPOM))
      ENDDO

C  reggeon
      SIGSR = ABS(XI(2))*S**(-ABS(XI(4)))
      SIG_BRN(1) = SIGSR/CMBARN
C  pomeron (soft part)
      SIGSP = ABS(XI(1))*S**ABS(XI(3))
      SIG_BRN(2) = SIGSP/CMBARN
C  pomeron (hard part)
      SIG_BRN(3) = SIGHAR/CMBARN

C  2x2 channel low-mass model and separate high-mass diffraction
        
      al1 = XI(5)
      be1 = XI(6)
      al2 = al1
      be2 = be1
      EnhPP  = XI(9)
      EnhPiP = EnhPP

      CALL SIG_JET_3 (SIG_brn,JINT,SIG_tot,SIG_ela,SIG_ine,SIG_sum,
     &                SIG_dif_1,SIG_dif_2,SIG_dd,B_el,PROB)

      SIGTOT = SIG_tot*CMBARN
      SIGINE = SIG_ine*CMBARN
      SIGSUM = SIG_sum*CMBARN
      SIGELc = SIGTOT-SIGINE
      SIGEL  = SIG_ela*CMBARN
      SIGSD1(1) = SIG_dif_1(1)*CMBARN
      SIGSD1(2) = SIG_dif_1(2)*CMBARN
      SIGSD2(1) = SIG_dif_2(1)*CMBARN
      SIGSD2(2) = SIG_dif_2(2)*CMBARN
      SIGDD(1)  = SIG_dd(1)*CMBARN
      SIGDD(2)  = SIG_dd(2)*CMBARN
      SLOPE  = B_EL
      SLOPEc = SIG_tot**2/(16.*Pi*SIG_ela)

      DE = ABS(SIGEL+SIGINE-SIGTOT)/SIGTOT
      IF(DE.GT.0.01) THEN
        print *,'SIBSIG:      Ecm: ',ECM
        print *,'          SIGTOT: ',SIGTOT
        print *,'        SIGEL1/2: ',SIGEL,SIGELc
        print *,'        SLOPE1/2: ',SLOPE,SLOPEc
        print *,'        SIGDIF 1: ',SIGSD1
        print *,'        SIGDIF 2: ',SIGSD2
        print *,'         SIGDDIF: ',SIGDD
        print *,'      SUM-SIGTOT: ',SIGEL+SIGINE-SIGTOT
      ENDIF

C  SIBYLL interface to single precision

      SIB_PTmin   = PTCUT
      SIB_SIG_tot = SIGTOT
      SIB_SIG_ine = SIGINE
      SIB_diff(1) = SIGSD1(1)+SIGSD1(2)
      SIB_diff(2) = SIGSD2(1)+SIGSD2(2)
      SIB_diff(3) = SIGDD(1)+SIGDD(2)
      SIB_B_el    = SLOPE
      DO I=0,NS_max
        DO K=0,NH_max
          SIB_PJET(I,K) = PROB(I,K)
        ENDDO
      ENDDO
c     full diff. cross section 
c     ( ( b.single , t.single , double ) , ( low mass , high mass ) ) 
      SIB_diff2(1,1) = SIGSD1(1)
      SIB_diff2(1,2) = SIGSD1(2)
      SIB_diff2(2,1) = SIGSD2(1)
      SIB_diff2(2,2) = SIGSD2(2)
      SIB_diff2(3,1) = SIGDD(1)
      SIB_diff2(3,2) = SIGDD(2)
      END


      SUBROUTINE SIG_JET_3 (SIG_brn, JINT, SIG_TOT, SIG_ELA, 
     &        SIG_INE, SIG_sum, SIG_DIF1, SIG_DIF2, SIG_DD, B_EL, P_int)
C-----------------------------------------------------------------------
C
C...This subroutine  receives in INPUT:
C.       SIG_brn (GeV-2)  Born graph cross sections
C.       JINT (1 = pp interaction)    (2 pi-p interaction)
C.       neg. value: without calculation of interaction probabilities
C.
C.  and returns as output:
C.       SIG_???  , B_el
C.       and P_int(0:NS_max,0:NH_max)   interaction probabilities
C
C   two x two -channel approximation for diffraction
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE

      DIMENSION SIG_brn(3)
      PARAMETER (NS_max = 20, NH_max = 80)

      COMMON /S_CFACT/ FACT(0:NH_max), CO_BIN(0:NH_max,0:NH_max)
      COMMON /S_CHDCNV/ABR(2,400),ABP(2,400),ABH(2,400),DB,NB

      COMMON /PROFILE/XNUS2,XMUS2,XNUSPI2,
     &                XNUH2,XMUH2,XNUHPI2,
     &                EnhPP,EnhPiP,al1,be1,al2,be2

      DIMENSION SIG_DIF1(2),SIG_DIF2(2),SIG_DD(2),
     &          P_int(0:NS_max,0:NH_max)
      DATA PI /3.1415926/

      DO J=0,NH_max
        DO I=0,NS_max
          P_int(I,J) = 0.
        ENDDO
      ENDDO

      ga1 = sqrt(al1*al1+be1*be1)
      ga2 = sqrt(al2*al2+be2*be2)

      fe_a_1  = (1.+al1/ga1)/2.
      fe_a_2  = (1.-al1/ga1)/2.
      fd_a_1  = sqrt(1.-(al1/ga1)**2)/2.
      fd_a_2  = -fd_a_1

      fe_b_1  = (1.+al2/ga2)/2.
      fe_b_2  = (1.-al2/ga2)/2.
      fd_b_1  = sqrt(1.-(al2/ga2)**2)/2.
      fd_b_2  = -fd_b_1

      fe_11 = fe_a_1*fe_b_1
      fe_22 = fe_a_2*fe_b_2
      fe_12 = fe_a_1*fe_b_2
      fe_21 = fe_a_2*fe_b_1

      fd_a_11 = fd_a_1*fe_b_1
      fd_a_22 = fd_a_2*fe_b_2
      fd_a_12 = fd_a_1*fe_b_2
      fd_a_21 = fd_a_2*fe_b_1

      fd_b_11 = fe_a_1*fd_b_1
      fd_b_22 = fe_a_2*fd_b_2
      fd_b_12 = fe_a_1*fd_b_2
      fd_b_21 = fe_a_2*fd_b_1

      fdd_11 = fd_a_1*fd_b_1
      fdd_22 = fd_a_2*fd_b_2
      fdd_12 = fd_a_1*fd_b_2
      fdd_21 = fd_a_2*fd_b_1


      sum_abs = 0.
      sum_tot = 0.
      sum_ela = 0.
      sum_sd_a = 0.
      sum_sd_b = 0.
      sum_dd  = 0.
      sum_B   = 0.

      IK = ABS(JINT)
      if(JINT.GT.0) then
        I0MAX = NS_max
        J0MAX = NH_max
      ELSE
        I0MAX = 1
        J0MAX = 1
      ENDIF
      SIG_REG = SIG_BRN(1)
      SIG_POM = SIG_BRN(2)
      SIG_HAR = SIG_BRN(3)

      DO JB=1,NB

         B = DB*FLOAT(JB-1)

         ABREG = ABR(IK,JB)
         ABPOM = ABP(IK,JB)
         ABHAR = ABH(IK,JB)

         chi2_soft = ABREG*SIG_REG+ABPOM*SIG_POM
         chi2_soft_11 = (1.-al1+ga1)*(1.-al2+ga2)*chi2_soft
         chi2_soft_22 = (1.-al1-ga1)*(1.-al2-ga2)*chi2_soft
         chi2_soft_12 = (1.-al1+ga1)*(1.-al2-ga2)*chi2_soft
         chi2_soft_21 = (1.-al1-ga1)*(1.-al2+ga2)*chi2_soft

         chi2_hard = ABHAR*SIG_HAR
         chi2_hard_11 = (1.-al1+ga1)*(1.-al2+ga2)*chi2_hard
         chi2_hard_22 = (1.-al1-ga1)*(1.-al2-ga2)*chi2_hard
         chi2_hard_12 = (1.-al1+ga1)*(1.-al2-ga2)*chi2_hard
         chi2_hard_21 = (1.-al1-ga1)*(1.-al2+ga2)*chi2_hard
          

         ef_11  = exp(-0.5*(chi2_soft_11+chi2_hard_11))
         ef_22  = exp(-0.5*(chi2_soft_22+chi2_hard_22))
         ef_12 = exp(-0.5*(chi2_soft_12+chi2_hard_12))
         ef_21 = exp(-0.5*(chi2_soft_21+chi2_hard_21))

         esf_11  = ef_11**2
         esf_22  = ef_22**2
         esf_12  = ef_12**2
         esf_21  = ef_21**2

         F_ine = B*(1. - fe_11*esf_11 - fe_12*esf_12 
     &                 - fe_21*esf_21 - fe_22*esf_22)
         F_tot = 1. - fe_11*ef_11 - fe_12*ef_12
     &              - fe_21*ef_21 - fe_22*ef_22
         F_ela = B*F_tot**2
         F_tot = B*F_tot

         F_sd_a = B*(fd_a_11*ef_11 + fd_a_12*ef_12
     &              + fd_a_21*ef_21 + fd_a_22*ef_22)**2
         F_sd_b = B*(fd_b_11*ef_11 + fd_b_12*ef_12
     &              + fd_b_21*ef_21 + fd_b_22*ef_22)**2
         F_dd  = B*(fdd_11*ef_11 + fdd_12*ef_12
     &              + fdd_21*ef_21 + fdd_22*ef_22)**2

         sum_abs = sum_abs+F_ine
         sum_tot = sum_tot+F_tot
         sum_ela = sum_ela+F_ela

         sum_sd_a = sum_sd_a+F_sd_a
         sum_sd_b = sum_sd_b+F_sd_b
         sum_dd  = sum_dd +F_dd 

         sum_B   = sum_b+B**2*F_tot

         fac_11 = B*esf_11
         fac_22 = B*esf_22
         fac_12 = B*esf_12
         fac_21 = B*esf_21
         soft_rec_11 = 1./chi2_soft_11
         soft_rec_22 = 1./chi2_soft_22
         soft_rec_12 = 1./chi2_soft_12
         soft_rec_21 = 1./chi2_soft_21
         chi2_hard_11 = max(chi2_hard_11,1.d-10)
         chi2_hard_22 = max(chi2_hard_22,1.d-10)
         chi2_hard_12 = max(chi2_hard_12,1.d-10)
         chi2_hard_21 = max(chi2_hard_21,1.d-10)
         DO I=0,I0MAX
           soft_rec_11 = soft_rec_11*chi2_soft_11
           soft_rec_22 = soft_rec_22*chi2_soft_22
           soft_rec_12 = soft_rec_12*chi2_soft_12
           soft_rec_21 = soft_rec_21*chi2_soft_21
           hard_rec_11 = 1./chi2_hard_11
           hard_rec_22 = 1./chi2_hard_22
           hard_rec_12 = 1./chi2_hard_12
           hard_rec_21 = 1./chi2_hard_21
           DO J=0,J0MAX
             hard_rec_11 = hard_rec_11*chi2_hard_11
             hard_rec_22 = hard_rec_22*chi2_hard_22
             hard_rec_12 = hard_rec_12*chi2_hard_12
             hard_rec_21 = hard_rec_21*chi2_hard_21
             P_int(I,J) = P_int(I,J) 
     &                + fe_11*soft_rec_11*hard_rec_11*fac_11
     &                + fe_22*soft_rec_22*hard_rec_22*fac_22
     &                + fe_12*soft_rec_12*hard_rec_12*fac_12
     &                + fe_21*soft_rec_21*hard_rec_21*fac_21
           ENDDO
         ENDDO

      ENDDO

      SIG_abs  = SUM_abs*2.*PI*DB
      SIG_tot  = SUM_tot*4.*PI*DB
      SIG_ela  = SUM_ela*2.*PI*DB
      SIG_dif1(1) = SUM_sd_a*2.*PI*DB
      SIG_dif2(1) = SUM_sd_b*2.*PI*DB
      SIG_dd(1)   = SUM_dd*2.*PI*DB
      SIG_ine  = SIG_abs + SIG_dif1(1) + SIG_dif2(1) + SIG_dd(1)
      B_EL     = sum_B/SUM_tot/2.

      SA = 0.
      P_int(0,0) = 0.
      DO I=0,I0MAX
        DO J=0,J0MAX
          fac = FACT(I)*FACT(J)
          P_int(I,J) = P_int(I,J)/fac
          SA = SA + P_int(I,J)
        ENDDO
      ENDDO

      SIG_hmsd = EnhPP*(P_int(1,0)+P_int(0,1))*2.*PI*DB
      SIG_hmdd = be1**2*SIG_hmsd + be2**2*SIG_hmsd
     &          + EnhPP**2*P_int(1,1)*2.*PI*DB

      SIG_dif1(2) = SIG_hmsd
      SIG_dif2(2) = SIG_hmsd
      SIG_dd(2)   = SIG_hmdd

      SIG_sum = SA*2.*PI*DB

      DO I=0,I0MAX
        DO J=0,J0MAX
          P_int(I,J) = P_int(I,J)/SA
        ENDDO
      ENDDO

      END
C
C
      SUBROUTINE HAD_CONV(JINT)
C-----------------------------------------------------------------------
C
C...Convolution of hadrons profile
C.  [function A(b) of Durand and Pi]
C.  precalculate and put  in COMMON block
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
C
      COMMON /S_CHDCNV/ABR(2,400),ABP(2,400),ABH(2,400),DB,NB

      DOUBLE PRECISION NU2, MU2, NUPI2, NU, MU, NUPI
      COMMON /S_CH0CNV/ NU2, MU2, NUPI2, NU, MU, NUPI

C
      COMMON /PROFILE/XNUS2,XMUS2,XNUSPI2,
     &                XNUH2,XMUH2,XNUHPI2,
     &                ENHPP,ENHPIP,al1,be1,al2,be2

C...integration constants
      BMAX = 50.
      NB  = 400
      DB = BMAX/FLOAT(NB)

C  soft reggeon interactions

      NU2   = XNUS2
      MU2   = XMUS2
      NUPI2 = XNUSPI2

      NU = SQRT(NU2)
      MU = SQRT(ABS(MU2))
      NUPI = SQRT(NUPI2) 

      DO JB=1,NB
         B = DB*FLOAT(JB-1)
         IF(JINT.EQ.1) THEN
           ABR(JINT,JB) = A_PP(B)
         ELSE
           ABR(JINT,JB) = A_PIP(B)
         ENDIF
      ENDDO 

C  soft pomeron interactions

      NU2   = XNUS2
      MU2   = XMUS2
      NUPI2 = XNUSPI2

      NU = SQRT(NU2)
      MU = SQRT(ABS(MU2))
      NUPI = SQRT(NUPI2)

      DO JB=1,NB
         B = DB*FLOAT(JB-1)
         IF(JINT.EQ.1) THEN
           ABP(JINT,JB) = A_PP(B)
         ELSE
           ABP(JINT,JB) = A_PIP(B)
         ENDIF
      ENDDO

C  hard pomeron interactions

      NU2   = XNUH2
      MU2   = XMUH2
      NUPI2 = XNUHPI2

      NU = SQRT(NU2)
      MU = SQRT(ABS(MU2))
      NUPI = SQRT(NUPI2)

      DB = BMAX/FLOAT(NB)
      DO JB=1,NB
         B = DB*FLOAT(JB-1)
         IF(JINT.EQ.1) THEN
           ABH(JINT,JB) = A_PP(B)
         ELSE
           ABH(JINT,JB) = A_PIP(B)
         ENDIF
      ENDDO

      END
C
C
      DOUBLE PRECISION FUNCTION A_pp (b)
C...Convolution of parton distribution for pp interaction
      IMPLICIT DOUBLE PRECISION (A-Z)
      SAVE
C
      COMMON /S_CH0CNV/ NU2, MU2, NUPI2, NU, MU, NUPI
      data pi / 3.1415926/

      ETA = NU2/MU2
 
      IF(ETA.LT.0.D0) THEN
   
        c = nu**5/(96.*pi)
        if (b .gt. 0.0001D0)  then
           A_pp = c*b**3 * bessk (3, b*nu)
        else
           A_pp = nu**2/(12.*pi)
        endif

      ELSE

        X = B*NU
        Y = B*MU
        C = NU2/(12.*PI)/(1.-ETA)**2
        IF(X.GT.0.0001D0) THEN
          A_PP = C*(1./8.*X**3*BESSK(3,X)
     &          -3./2.*ETA/(1.-ETA)*X**2*BESSK(2,X)
     &          +9*ETA**2/(1.-ETA)**2*X*BESSK1(X)
     &          -24*ETA**3/(1.-ETA)**3*(BESSK0(X)-BESSK0(Y))
     &          +3.*ETA**3/(1.-ETA)**2*Y*BESSK1(Y))
        ELSE
          A_PP = C*(1./8.*8.
     &          -3./2.*ETA/(1.-ETA)*2.
     &          +9*ETA**2/(1.-ETA)**2*1.
     &          -24*ETA**3/(1.-ETA)**3*LOG(MU/NU)
     &          +3.*ETA**3/(1.-ETA)**2*1.)
        ENDIF

      ENDIF

      END
C
C
      DOUBLE PRECISION FUNCTION A_pip (b)
C...Convolution of parton distribution for pip interaction
      IMPLICIT DOUBLE PRECISION (A-Z)
      SAVE
C
      COMMON /S_CH0CNV/ NU2, MU2, NUPI2, NU, MU, NUPI
      data pi / 3.1415926/

      eta = nu2/nupi2
      c = nu2/(2.*pi) * 1./(1.-eta)

      if (b .gt. 0.0001D0)  then
         b1 = b*nu
         b2 = b*nupi
         f1 = 0.5*b1 * bessk1(b1)
         f2 = eta/(1.-eta)*(bessk0(b2)- bessk0(b1))
         A_pip = c*(f1+f2)
      else
         A_pip = c*(0.5 + eta/(1.-eta)*log(nu/nupi))
      endif
      return
      end
C
C
C----------------------------------------------------------------------------
C  Bessel functions
C----------------------------------------------------------------------------
C
      FUNCTION BESSK0(X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
C
*     REAL*8 Y,P1,P2,P3,P4,P5,P6,P7,
*    *    Q1,Q2,Q3,Q4,Q5,Q6,Q7
      DATA P1,P2,P3,P4,P5,P6,P7/-0.57721566D0,0.42278420D0,
     *    0.23069756D0,0.3488590D-1,0.262698D-2,0.10750D-3,0.74D-5/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,-0.7832358D-1,
     * 0.2189568D-1,-0.1062446D-1,0.587872D-2,-0.251540D-2,0.53208D-3/
      IF (X.LE.2.0) THEN
        Y=X*X/4.0
        BESSK0=(-LOG(X/2.0)*BESSI0(X))+(P1+Y*(P2+Y*(P3+
     *        Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        Y=(2.0/X)
        BESSK0=(EXP(-X)/SQRT(X))*(Q1+Y*(Q2+Y*(Q3+
     *        Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END
C
C
      FUNCTION BESSK1(X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
C
*     REAL*8 Y,P1,P2,P3,P4,P5,P6,P7,
*    *    Q1,Q2,Q3,Q4,Q5,Q6,Q7
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,0.15443144D0,-0.67278579D0,
     *    -0.18156897D0,-0.1919402D-1,-0.110404D-2,-0.4686D-4/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,0.23498619D0,
     *    -0.3655620D-1,0.1504268D-1,-0.780353D-2,0.325614D-2,
     *    -0.68245D-3/
      IF (X.LE.2.0) THEN
        Y=X*X/4.0
        BESSK1=(LOG(X/2.0)*BESSI1(X))+(1.0/X)*(P1+Y*(P2+
     *      Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        Y=2.0/X
        BESSK1=(EXP(-X)/SQRT(X))*(Q1+Y*(Q2+Y*(Q3+
     *      Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END
C
C
      FUNCTION BESSK(N,X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
C
      IF (N.LT.2) PAUSE 'bad argument N in BESSK'
      TOX=2.0/X
      BKM=BESSK0(X)
      BK=BESSK1(X)
      DO 11 J=1,N-1
        BKP=BKM+J*TOX*BK
        BKM=BK
        BK=BKP
11    CONTINUE
      BESSK=BK
      RETURN
      END
C
C
      FUNCTION BESSI0(X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
C
*     REAL*8 Y,P1,P2,P3,P4,P5,P6,P7,
*    *    Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,3.5156229D0,3.0899424D0,
     *    1.2067492D0,
     *    0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,
     *    0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,
     *    0.2635537D-1,-0.1647633D-1,0.392377D-2/
      IF (ABS(X).LT.3.75) THEN
        Y=(X/3.75)**2
        BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
        AX=ABS(X)
        Y=3.75/AX
        BESSI0=(EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4
     *      +Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END
C
C
      FUNCTION BESSI1(X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
C
*     REAL*8 Y,P1,P2,P3,P4,P5,P6,P7,
*    *    Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
      DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,
     *    0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1,
     *    -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1,
     *    -0.2895312D-1,0.1787654D-1,-0.420059D-2/
      IF (ABS(X).LT.3.75) THEN
        Y=(X/3.75)**2
        BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        AX=ABS(X)
        Y=3.75/AX
        BESSI1=(EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+
     *      Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END


      SUBROUTINE FACT_INI
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE

      PARAMETER (NS_max = 20, NH_max = 80)
      COMMON /S_CFACT/ FACT(0:NH_max), CO_BIN(0:NH_max,0:NH_max)

      FACT(0) = 1.
      DO J=1,NH_max
         FACT(J) = FACT(J-1)*FLOAT(J)
      ENDDO
      DO J=0,NH_max
         DO K=0,J
            CO_BIN(J,K) = FACT(J)/(FACT(K)*FACT(J-K))
         ENDDO
      ENDDO

      END


      SUBROUTINE SIB_GAUSET(AX,BX,NX,Z,W)
C-----------------------------------------------------------------------
C
C     N-point gauss zeros and weights for the interval (AX,BX) are
C           stored in  arrays Z and W respectively.
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      SAVE
      COMMON /GQCOM/A(273),X(273),KTAB(96)
      DIMENSION Z(NX),W(NX)
      DATA INIT/0/
C
      ALPHA=0.5*(BX+AX)
      BETA=0.5*(BX-AX)
      N=NX
*
*  the N=1 case:
      IF(N.NE.1) GO TO 1
      Z(1)=ALPHA
      W(1)=BX-AX
      RETURN
*
*  the Gauss cases:
    1 IF((N.LE.16).AND.(N.GT.1)) GO TO 2
      IF(N.EQ.20) GO TO 2
      IF(N.EQ.24) GO TO 2
      IF(N.EQ.32) GO TO 2
      IF(N.EQ.40) GO TO 2
      IF(N.EQ.48) GO TO 2
      IF(N.EQ.64) GO TO 2
      IF(N.EQ.80) GO TO 2
      IF(N.EQ.96) GO TO 2
*
*  the extended Gauss cases:
      IF((N/96)*96.EQ.N) GO TO 3
*
C  jump to center of intervall intrgration:
      GO TO 100
*
C  get Gauss point array
*
    2 CALL PO106BD
C     -print out message
*     IF(INIT.LE.20)THEN
*       INIT=init+1
*       WRITE (6,*) ' initialization of Gauss int. N=',N
*     ENDIF
C  extract real points
      K=KTAB(N)
      M=N/2
      DO 21 J=1,M
C       extract values from big array
        JTAB=K-1+J
        WTEMP=BETA*A(JTAB)
        DELTA=BETA*X(JTAB)
C       store them backward
        Z(J)=ALPHA-DELTA
        W(J)=WTEMP
C       store them forward
        JP=N+1-J
        Z(JP)=ALPHA+DELTA
        W(JP)=WTEMP
   21 CONTINUE
C     store central point (odd N)
      IF((N-M-M).EQ.0) RETURN
      Z(M+1)=ALPHA
      JMID=K+M
      W(M+1)=BETA*A(JMID)
      RETURN
C
C  get ND96 times chained 96 Gauss point array
C
    3 CALL PO106BD
C  print out message
      IF(INIT.LE.20)THEN
        INIT=init+1
        WRITE (6,*) ' initialization of extended Gauss int. N=',N
      ENDIF
C     -extract real points
      K=KTAB(96)
      ND96=N/96
      DO 31 J=1,48
C       extract values from big array
        JTAB=K-1+J
        WTEMP=BETA*A(JTAB)
        DELTA=BETA*X(JTAB)
        WTeMP=WTEMP/ND96
        DeLTA=DELTA/ND96
        DO 32 JD96=0,ND96-1
          ZCNTR= (ALPHA-BETA)+ BETA*FLOAT(2*JD96+1)/FLOAT(ND96)
C         store them backward
          Z(J+JD96*96)=ZCNTR-DELTA
          W(J+JD96*96)=WTEMP
C         store them forward
          JP=96+1-J
          Z(JP+JD96*96)=ZCNTR+DELTA
          W(JP+JD96*96)=WTEMP
   32   CONTINUE
   31 CONTINUE
      RETURN
*
C  the center of intervall cases:
  100 CONTINUE
C  print out message
      IF(INIT.LE.20)THEN
        INIT=init+1
        WRITE (6,*) ' init. of center of intervall int. N=',N
      ENDIF
C  put in constant weight and equally spaced central points
      N=IABS(N)
      DO 111 IN=1,N
        WIN=(BX-AX)/FLOAT(N)
        Z(IN)=AX  + (FLOAT(IN)-.5)*WIN
  111 W(IN)=WIN
      RETURN
      END
C
C
      SUBROUTINE PO106BD
C-----------------------------------------------------------------------
C
C     store big arrays needed for Gauss integral, CERNLIB D106BD
C     (arrays A,X,ITAB copied on B,Y,LTAB)
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      SAVE
      COMMON /GQCOM/ B(273),Y(273),LTAB(96)
      DIMENSION      A(273),X(273),KTAB(96)
C
C-----TABLE OF INITIAL SUBSCRIPTS FOR N=2(1)16(4)96
      DATA KTAB(2)/1/
      DATA KTAB(3)/2/
      DATA KTAB(4)/4/
      DATA KTAB(5)/6/
      DATA KTAB(6)/9/
      DATA KTAB(7)/12/
      DATA KTAB(8)/16/
      DATA KTAB(9)/20/
      DATA KTAB(10)/25/
      DATA KTAB(11)/30/
      DATA KTAB(12)/36/
      DATA KTAB(13)/42/
      DATA KTAB(14)/49/
      DATA KTAB(15)/56/
      DATA KTAB(16)/64/
      DATA KTAB(20)/72/
      DATA KTAB(24)/82/
      DATA KTAB(28)/82/
      DATA KTAB(32)/94/
      DATA KTAB(36)/94/
      DATA KTAB(40)/110/
      DATA KTAB(44)/110/
      DATA KTAB(48)/130/
      DATA KTAB(52)/130/
      DATA KTAB(56)/130/
      DATA KTAB(60)/130/
      DATA KTAB(64)/154/
      DATA KTAB(68)/154/
      DATA KTAB(72)/154/
      DATA KTAB(76)/154/
      DATA KTAB(80)/186/
      DATA KTAB(84)/186/
      DATA KTAB(88)/186/
      DATA KTAB(92)/186/
      DATA KTAB(96)/226/
C
C-----TABLE OF ABSCISSAE (X) AND WEIGHTS (A) FOR INTERVAL (-1,+1).
C
C-----N=2
      DATA X(1)/0.577350269189626D0  /, A(1)/1.000000000000000D0  /
C-----N=3
      DATA X(2)/0.774596669241483D0  /, A(2)/0.555555555555556D0  /
      DATA X(3)/0.000000000000000D0  /, A(3)/0.888888888888889D0  /
C-----N=4
      DATA X(4)/0.861136311594053D0  /, A(4)/0.347854845137454D0  /
      DATA X(5)/0.339981043584856D0  /, A(5)/0.652145154862546D0  /
C-----N=5
      DATA X(6)/0.906179845938664D0  /, A(6)/0.236926885056189D0  /
      DATA X(7)/0.538469310105683D0  /, A(7)/0.478628670499366D0  /
      DATA X(8)/0.000000000000000D0  /, A(8)/0.568888888888889D0  /
C-----N=6
      DATA X(9)/0.932469514203152D0  /, A(9)/0.171324492379170D0  /
      DATA X(10)/0.661209386466265D0 /, A(10)/0.360761573048139D0 /
      DATA X(11)/0.238619186083197D0 /, A(11)/0.467913934572691D0 /
C-----N=7
      DATA X(12)/0.949107912342759D0 /, A(12)/0.129484966168870D0 /
      DATA X(13)/0.741531185599394D0 /, A(13)/0.279705391489277D0 /
      DATA X(14)/0.405845151377397D0 /, A(14)/0.381830050505119D0 /
      DATA X(15)/0.000000000000000D0 /, A(15)/0.417959183673469D0 /
C-----N=8
      DATA X(16)/0.960289856497536D0 /, A(16)/0.101228536290376D0 /
      DATA X(17)/0.796666477413627D0 /, A(17)/0.222381034453374D0 /
      DATA X(18)/0.525532409916329D0 /, A(18)/0.313706645877887D0 /
      DATA X(19)/0.183434642495650D0 /, A(19)/0.362683783378362D0 /
C-----N=9
      DATA X(20)/0.968160239507626D0 /, A(20)/0.081274388361574D0 /
      DATA X(21)/0.836031107326636D0 /, A(21)/0.180648160694857D0 /
      DATA X(22)/0.613371432700590D0 /, A(22)/0.260610696402935D0 /
      DATA X(23)/0.324253423403809D0 /, A(23)/0.312347077040003D0 /
      DATA X(24)/0.000000000000000D0 /, A(24)/0.330239355001260D0 /
C-----N=10
      DATA X(25)/0.973906528517172D0 /, A(25)/0.066671344308688D0 /
      DATA X(26)/0.865063366688985D0 /, A(26)/0.149451349150581D0 /
      DATA X(27)/0.679409568299024D0 /, A(27)/0.219086362515982D0 /
      DATA X(28)/0.433395394129247D0 /, A(28)/0.269266719309996D0 /
      DATA X(29)/0.148874338981631D0 /, A(29)/0.295524224714753D0 /
C-----N=11
      DATA X(30)/0.978228658146057D0 /, A(30)/0.055668567116174D0 /
      DATA X(31)/0.887062599768095D0 /, A(31)/0.125580369464905D0 /
      DATA X(32)/0.730152005574049D0 /, A(32)/0.186290210927734D0 /
      DATA X(33)/0.519096129206812D0 /, A(33)/0.233193764591990D0 /
      DATA X(34)/0.269543155952345D0 /, A(34)/0.262804544510247D0 /
      DATA X(35)/0.000000000000000D0 /, A(35)/0.272925086777901D0 /
C-----N=12
      DATA X(36)/0.981560634246719D0 /, A(36)/0.047175336386512D0 /
      DATA X(37)/0.904117256370475D0 /, A(37)/0.106939325995318D0 /
      DATA X(38)/0.769902674194305D0 /, A(38)/0.160078328543346D0 /
      DATA X(39)/0.587317954286617D0 /, A(39)/0.203167426723066D0 /
      DATA X(40)/0.367831498998180D0 /, A(40)/0.233492536538355D0 /
      DATA X(41)/0.125233408511469D0 /, A(41)/0.249147045813403D0 /
C-----N=13
      DATA X(42)/0.984183054718588D0 /, A(42)/0.040484004765316D0 /
      DATA X(43)/0.917598399222978D0 /, A(43)/0.092121499837728D0 /
      DATA X(44)/0.801578090733310D0 /, A(44)/0.138873510219787D0 /
      DATA X(45)/0.642349339440340D0 /, A(45)/0.178145980761946D0 /
      DATA X(46)/0.448492751036447D0 /, A(46)/0.207816047536889D0 /
      DATA X(47)/0.230458315955135D0 /, A(47)/0.226283180262897D0 /
      DATA X(48)/0.000000000000000D0 /, A(48)/0.232551553230874D0 /
C-----N=14
      DATA X(49)/0.986283808696812D0 /, A(49)/0.035119460331752D0 /
      DATA X(50)/0.928434883663574D0 /, A(50)/0.080158087159760D0 /
      DATA X(51)/0.827201315069765D0 /, A(51)/0.121518570687903D0 /
      DATA X(52)/0.687292904811685D0 /, A(52)/0.157203167158194D0 /
      DATA X(53)/0.515248636358154D0 /, A(53)/0.185538397477938D0 /
      DATA X(54)/0.319112368927890D0 /, A(54)/0.205198463721296D0 /
      DATA X(55)/0.108054948707344D0 /, A(55)/0.215263853463158D0 /
C-----N=15
      DATA X(56)/0.987992518020485D0 /, A(56)/0.030753241996117D0 /
      DATA X(57)/0.937273392400706D0 /, A(57)/0.070366047488108D0 /
      DATA X(58)/0.848206583410427D0 /, A(58)/0.107159220467172D0 /
      DATA X(59)/0.724417731360170D0 /, A(59)/0.139570677926154D0 /
      DATA X(60)/0.570972172608539D0 /, A(60)/0.166269205816994D0 /
      DATA X(61)/0.394151347077563D0 /, A(61)/0.186161000015562D0 /
      DATA X(62)/0.201194093997435D0 /, A(62)/0.198431485327111D0 /
      DATA X(63)/0.000000000000000D0 /, A(63)/0.202578241925561D0 /
C-----N=16
      DATA X(64)/0.989400934991650D0 /, A(64)/0.027152459411754D0 /
      DATA X(65)/0.944575023073233D0 /, A(65)/0.062253523938648D0 /
      DATA X(66)/0.865631202387832D0 /, A(66)/0.095158511682493D0 /
      DATA X(67)/0.755404408355003D0 /, A(67)/0.124628971255534D0 /
      DATA X(68)/0.617876244402644D0 /, A(68)/0.149595988816577D0 /
      DATA X(69)/0.458016777657227D0 /, A(69)/0.169156519395003D0 /
      DATA X(70)/0.281603550779259D0 /, A(70)/0.182603415044924D0 /
      DATA X(71)/0.095012509837637D0 /, A(71)/0.189450610455069D0 /
C-----N=20
      DATA X(72)/0.993128599185094D0 /, A(72)/0.017614007139152D0 /
      DATA X(73)/0.963971927277913D0 /, A(73)/0.040601429800386D0 /
      DATA X(74)/0.912234428251325D0 /, A(74)/0.062672048334109D0 /
      DATA X(75)/0.839116971822218D0 /, A(75)/0.083276741576704D0 /
      DATA X(76)/0.746331906460150D0 /, A(76)/0.101930119817240D0 /
      DATA X(77)/0.636053680726515D0 /, A(77)/0.118194531961518D0 /
      DATA X(78)/0.510867001950827D0 /, A(78)/0.131688638449176D0 /
      DATA X(79)/0.373706088715419D0 /, A(79)/0.142096109318382D0 /
      DATA X(80)/0.227785851141645D0 /, A(80)/0.149172986472603D0 /
      DATA X(81)/0.076526521133497D0 /, A(81)/0.152753387130725D0 /
C-----N=24
      DATA X(82)/0.995187219997021D0 /, A(82)/0.012341229799987D0 /
      DATA X(83)/0.974728555971309D0 /, A(83)/0.028531388628933D0 /
      DATA X(84)/0.938274552002732D0 /, A(84)/0.044277438817419D0 /
      DATA X(85)/0.886415527004401D0 /, A(85)/0.059298584915436D0 /
      DATA X(86)/0.820001985973902D0 /, A(86)/0.073346481411080D0 /
      DATA X(87)/0.740124191578554D0 /, A(87)/0.086190161531953D0 /
      DATA X(88)/0.648093651936975D0 /, A(88)/0.097618652104113D0 /
      DATA X(89)/0.545421471388839D0 /, A(89)/0.107444270115965D0 /
      DATA X(90)/0.433793507626045D0 /, A(90)/0.115505668053725D0 /
      DATA X(91)/0.315042679696163D0 /, A(91)/0.121670472927803D0 /
      DATA X(92)/0.191118867473616D0 /, A(92)/0.125837456346828D0 /
      DATA X(93)/0.064056892862605D0 /, A(93)/0.127938195346752D0 /
C-----N=32
      DATA X(94)/0.997263861849481D0 /, A(94)/0.007018610009470D0 /
      DATA X(95)/0.985611511545268D0 /, A(95)/0.016274394730905D0 /
      DATA X(96)/0.964762255587506D0 /, A(96)/0.025392065309262D0 /
      DATA X(97)/0.934906075937739D0 /, A(97)/0.034273862913021D0 /
      DATA X(98)/0.896321155766052D0 /, A(98)/0.042835898022226D0 /
      DATA X(99)/0.849367613732569D0 /, A(99)/0.050998059262376D0 /
      DATA X(100)/0.794483795967942D0/, A(100)/0.058684093478535D0/
      DATA X(101)/0.732182118740289D0/, A(101)/0.065822222776361D0/
      DATA X(102)/0.663044266930215D0/, A(102)/0.072345794108848D0/
      DATA X(103)/0.587715757240762D0/, A(103)/0.078193895787070D0/
      DATA X(104)/0.506899908932229D0/, A(104)/0.083311924226946D0/
      DATA X(105)/0.421351276130635D0/, A(105)/0.087652093004403D0/
      DATA X(106)/0.331868602282127D0/, A(106)/0.091173878695763D0/
      DATA X(107)/0.239287362252137D0/, A(107)/0.093844399080804D0/
      DATA X(108)/0.144471961582796D0/, A(108)/0.095638720079274D0/
      DATA X(109)/0.048307665687738D0/, A(109)/0.096540088514727D0/
C-----N=40
      DATA X(110)/0.998237709710559D0/, A(110)/0.004521277098533D0/
      DATA X(111)/0.990726238699457D0/, A(111)/0.010498284531152D0/
      DATA X(112)/0.977259949983774D0/, A(112)/0.016421058381907D0/
      DATA X(113)/0.957916819213791D0/, A(113)/0.022245849194166D0/
      DATA X(114)/0.932812808278676D0/, A(114)/0.027937006980023D0/
      DATA X(115)/0.902098806968874D0/, A(115)/0.033460195282547D0/
      DATA X(116)/0.865959503212259D0/, A(116)/0.038782167974472D0/
      DATA X(117)/0.824612230833311D0/, A(117)/0.043870908185673D0/
      DATA X(118)/0.778305651426519D0/, A(118)/0.048695807635072D0/
      DATA X(119)/0.727318255189927D0/, A(119)/0.053227846983936D0/
      DATA X(120)/0.671956684614179D0/, A(120)/0.057439769099391D0/
      DATA X(121)/0.612553889667980D0/, A(121)/0.061306242492928D0/
      DATA X(122)/0.549467125095128D0/, A(122)/0.064804013456601D0/
      DATA X(123)/0.483075801686178D0/, A(123)/0.067912045815233D0/
      DATA X(124)/0.413779204371605D0/, A(124)/0.070611647391286D0/
      DATA X(125)/0.341994090825758D0/, A(125)/0.072886582395804D0/
      DATA X(126)/0.268152185007253D0/, A(126)/0.074723169057968D0/
      DATA X(127)/0.192697580701371D0/, A(127)/0.076110361900626D0/
      DATA X(128)/0.116084070675255D0/, A(128)/0.077039818164247D0/
      DATA X(129)/0.038772417506050D0/, A(129)/0.077505947978424D0/
C-----N=48
      DATA X(130)/0.998771007252426D0/, A(130)/0.003153346052305D0/
      DATA X(131)/0.993530172266350D0/, A(131)/0.007327553901276D0/
      DATA X(132)/0.984124583722826D0/, A(132)/0.011477234579234D0/
      DATA X(133)/0.970591592546247D0/, A(133)/0.015579315722943D0/
      DATA X(134)/0.952987703160430D0/, A(134)/0.019616160457355D0/
      DATA X(135)/0.931386690706554D0/, A(135)/0.023570760839324D0/
      DATA X(136)/0.905879136715569D0/, A(136)/0.027426509708356D0/
      DATA X(137)/0.876572020274247D0/, A(137)/0.031167227832798D0/
      DATA X(138)/0.843588261624393D0/, A(138)/0.034777222564770D0/
      DATA X(139)/0.807066204029442D0/, A(139)/0.038241351065830D0/
      DATA X(140)/0.767159032515740D0/, A(140)/0.041545082943464D0/
      DATA X(141)/0.724034130923814D0/, A(141)/0.044674560856694D0/
      DATA X(142)/0.677872379632663D0/, A(142)/0.047616658492490D0/
      DATA X(143)/0.628867396776513D0/, A(143)/0.050359035553854D0/
      DATA X(144)/0.577224726083972D0/, A(144)/0.052890189485193D0/
      DATA X(145)/0.523160974722233D0/, A(145)/0.055199503699984D0/
      DATA X(146)/0.466902904750958D0/, A(146)/0.057277292100403D0/
      DATA X(147)/0.408686481990716D0/, A(147)/0.059114839698395D0/
      DATA X(148)/0.348755886292160D0/, A(148)/0.060704439165893D0/
      DATA X(149)/0.287362487355455D0/, A(149)/0.062039423159892D0/
      DATA X(150)/0.224763790394689D0/, A(150)/0.063114192286254D0/
      DATA X(151)/0.161222356068891D0/, A(151)/0.063924238584648D0/
      DATA X(152)/0.097004699209462D0/, A(152)/0.064466164435950D0/
      DATA X(153)/0.032380170962869D0/, A(153)/0.064737696812683D0/
C-----N=64
      DATA X(154)/0.999305041735772D0/, A(154)/0.001783280721696D0/
      DATA X(155)/0.996340116771955D0/, A(155)/0.004147033260562D0/
      DATA X(156)/0.991013371476744D0/, A(156)/0.006504457968978D0/
      DATA X(157)/0.983336253884625D0/, A(157)/0.008846759826363D0/
      DATA X(158)/0.973326827789910D0/, A(158)/0.011168139460131D0/
      DATA X(159)/0.961008799652053D0/, A(159)/0.013463047896718D0/
      DATA X(160)/0.946411374858402D0/, A(160)/0.015726030476024D0/
      DATA X(161)/0.929569172131939D0/, A(161)/0.017951715775697D0/
      DATA X(162)/0.910522137078502D0/, A(162)/0.020134823153530D0/
      DATA X(163)/0.889315445995114D0/, A(163)/0.022270173808383D0/
      DATA X(164)/0.865999398154092D0/, A(164)/0.024352702568710D0/
      DATA X(165)/0.840629296252580D0/, A(165)/0.026377469715054D0/
      DATA X(166)/0.813265315122797D0/, A(166)/0.028339672614259D0/
      DATA X(167)/0.783972358943341D0/, A(167)/0.030234657072402D0/
      DATA X(168)/0.752819907260531D0/, A(168)/0.032057928354851D0/
      DATA X(169)/0.719881850171610D0/, A(169)/0.033805161837141D0/
      DATA X(170)/0.685236313054233D0/, A(170)/0.035472213256882D0/
      DATA X(171)/0.648965471254657D0/, A(171)/0.037055128540240D0/
      DATA X(172)/0.611155355172393D0/, A(172)/0.038550153178615D0/
      DATA X(173)/0.571895646202634D0/, A(173)/0.039953741132720D0/
      DATA X(174)/0.531279464019894D0/, A(174)/0.041262563242623D0/
      DATA X(175)/0.489403145707052D0/, A(175)/0.042473515123653D0/
      DATA X(176)/0.446366017253464D0/, A(176)/0.043583724529323D0/
      DATA X(177)/0.402270157963991D0/, A(177)/0.044590558163756D0/
      DATA X(178)/0.357220158337668D0/, A(178)/0.045491627927418D0/
      DATA X(179)/0.311322871990210D0/, A(179)/0.046284796581314D0/
      DATA X(180)/0.264687162208767D0/, A(180)/0.046968182816210D0/
      DATA X(181)/0.217423643740007D0/, A(181)/0.047540165714830D0/
      DATA X(182)/0.169644420423992D0/, A(182)/0.047999388596458D0/
      DATA X(183)/0.121462819296120D0/, A(183)/0.048344762234802D0/
      DATA X(184)/0.072993121787799D0/, A(184)/0.048575467441503D0/
      DATA X(185)/0.024350292663424D0/, A(185)/0.048690957009139D0/
C-----N=80
      DATA X(186)/0.999553822651630D0/, A(186)/0.001144950003186D0/
      DATA X(187)/0.997649864398237D0/, A(187)/0.002663533589512D0/
      DATA X(188)/0.994227540965688D0/, A(188)/0.004180313124694D0/
      DATA X(189)/0.989291302499755D0/, A(189)/0.005690922451403D0/
      DATA X(190)/0.982848572738629D0/, A(190)/0.007192904768117D0/
      DATA X(191)/0.974909140585727D0/, A(191)/0.008683945269260D0/
      DATA X(192)/0.965485089043799D0/, A(192)/0.010161766041103D0/
      DATA X(193)/0.954590766343634D0/, A(193)/0.011624114120797D0/
      DATA X(194)/0.942242761309872D0/, A(194)/0.013068761592401D0/
      DATA X(195)/0.928459877172445D0/, A(195)/0.014493508040509D0/
      DATA X(196)/0.913263102571757D0/, A(196)/0.015896183583725D0/
      DATA X(197)/0.896675579438770D0/, A(197)/0.017274652056269D0/
      DATA X(198)/0.878722567678213D0/, A(198)/0.018626814208299D0/
      DATA X(199)/0.859431406663111D0/, A(199)/0.019950610878141D0/
      DATA X(200)/0.838831473580255D0/, A(200)/0.021244026115782D0/
      DATA X(201)/0.816954138681463D0/, A(201)/0.022505090246332D0/
      DATA X(202)/0.793832717504605D0/, A(202)/0.023731882865930D0/
      DATA X(203)/0.769502420135041D0/, A(203)/0.024922535764115D0/
      DATA X(204)/0.744000297583597D0/, A(204)/0.026075235767565D0/
      DATA X(205)/0.717365185362099D0/, A(205)/0.027188227500486D0/
      DATA X(206)/0.689637644342027D0/, A(206)/0.028259816057276D0/
      DATA X(207)/0.660859898986119D0/, A(207)/0.029288369583267D0/
      DATA X(208)/0.631075773046871D0/, A(208)/0.030272321759557D0/
      DATA X(209)/0.600330622829751D0/, A(209)/0.031210174188114D0/
      DATA X(210)/0.568671268122709D0/, A(210)/0.032100498673487D0/
      DATA X(211)/0.536145920897131D0/, A(211)/0.032941939397645D0/
      DATA X(212)/0.502804111888784D0/, A(212)/0.033733214984611D0/
      DATA X(213)/0.468696615170544D0/, A(213)/0.034473120451753D0/
      DATA X(214)/0.433875370831756D0/, A(214)/0.035160529044747D0/
      DATA X(215)/0.398393405881969D0/, A(215)/0.035794393953416D0/
      DATA X(216)/0.362304753499487D0/, A(216)/0.036373749905835D0/
      DATA X(217)/0.325664370747701D0/, A(217)/0.036897714638276D0/
      DATA X(218)/0.288528054884511D0/, A(218)/0.037365490238730D0/
      DATA X(219)/0.250952358392272D0/, A(219)/0.037776364362001D0/
      DATA X(220)/0.212994502857666D0/, A(220)/0.038129711314477D0/
      DATA X(221)/0.174712291832646D0/, A(221)/0.038424993006959D0/
      DATA X(222)/0.136164022809143D0/, A(222)/0.038661759774076D0/
      DATA X(223)/0.097408398441584D0/, A(223)/0.038839651059051D0/
      DATA X(224)/0.058504437152420D0/, A(224)/0.038958395962769D0/
      DATA X(225)/0.019511383256793D0/, A(225)/0.039017813656306D0/
C-----N=96
      DATA X(226)/0.999689503883230D0/, A(226)/0.000796792065552D0/
      DATA X(227)/0.998364375863181D0/, A(227)/0.001853960788946D0/
      DATA X(228)/0.995981842987209D0/, A(228)/0.002910731817934D0/
      DATA X(229)/0.992543900323762D0/, A(229)/0.003964554338444D0/
      DATA X(230)/0.988054126329623D0/, A(230)/0.005014202742927D0/
      DATA X(231)/0.982517263563014D0/, A(231)/0.006058545504235D0/
      DATA X(232)/0.975939174585136D0/, A(232)/0.007096470791153D0/
      DATA X(233)/0.968326828463264D0/, A(233)/0.008126876925698D0/
      DATA X(234)/0.959688291448742D0/, A(234)/0.009148671230783D0/
      DATA X(235)/0.950032717784437D0/, A(235)/0.010160770535008D0/
      DATA X(236)/0.939370339752755D0/, A(236)/0.011162102099838D0/
      DATA X(237)/0.927712456722308D0/, A(237)/0.012151604671088D0/
      DATA X(238)/0.915071423120898D0/, A(238)/0.013128229566961D0/
      DATA X(239)/0.901460635315852D0/, A(239)/0.014090941772314D0/
      DATA X(240)/0.886894517402420D0/, A(240)/0.015038721026994D0/
      DATA X(241)/0.871388505909296D0/, A(241)/0.015970562902562D0/
      DATA X(242)/0.854959033434601D0/, A(242)/0.016885479864245D0/
      DATA X(243)/0.837623511228187D0/, A(243)/0.017782502316045D0/
      DATA X(244)/0.819400310737931D0/, A(244)/0.018660679627411D0/
      DATA X(245)/0.800308744139140D0/, A(245)/0.019519081140145D0/
      DATA X(246)/0.780369043867433D0/, A(246)/0.020356797154333D0/
      DATA X(247)/0.759602341176647D0/, A(247)/0.021172939892191D0/
      DATA X(248)/0.738030643744400D0/, A(248)/0.021966644438744D0/
      DATA X(249)/0.715676812348967D0/, A(249)/0.022737069658329D0/
      DATA X(250)/0.692564536642171D0/, A(250)/0.023483399085926D0/
      DATA X(251)/0.668718310043916D0/, A(251)/0.024204841792364D0/
      DATA X(252)/0.644163403784967D0/, A(252)/0.024900633222483D0/
      DATA X(253)/0.618925840125468D0/, A(253)/0.025570036005349D0/
      DATA X(254)/0.593032364777572D0/, A(254)/0.026212340735672D0/
      DATA X(255)/0.566510418561397D0/, A(255)/0.026826866725591D0/
      DATA X(256)/0.539388108324357D0/, A(256)/0.027412962726029D0/
      DATA X(257)/0.511694177154667D0/, A(257)/0.027970007616848D0/
      DATA X(258)/0.483457973920596D0/, A(258)/0.028497411065085D0/
      DATA X(259)/0.454709422167743D0/, A(259)/0.028994614150555D0/
      DATA X(260)/0.425478988407300D0/, A(260)/0.029461089958167D0/
      DATA X(261)/0.395797649828908D0/, A(261)/0.029896344136328D0/
      DATA X(262)/0.365696861472313D0/, A(262)/0.030299915420827D0/
      DATA X(263)/0.335208522892625D0/, A(263)/0.030671376123669D0/
      DATA X(264)/0.304364944354496D0/, A(264)/0.031010332586313D0/
      DATA X(265)/0.273198812591049D0/, A(265)/0.031316425596861D0/
      DATA X(266)/0.241743156163840D0/, A(266)/0.031589330770727D0/
      DATA X(267)/0.210031310460567D0/, A(267)/0.031828758894411D0/
      DATA X(268)/0.178096882367618D0/, A(268)/0.032034456231992D0/
      DATA X(269)/0.145973714654896D0/, A(269)/0.032206204794030D0/
      DATA X(270)/0.113695850110665D0/, A(270)/0.032343822568575D0/
      DATA X(271)/0.081297495464425D0/, A(271)/0.032447163714064D0/
      DATA X(272)/0.048812985136049D0/, A(272)/0.032516118713868D0/
      DATA X(273)/0.016276744849602D0/, A(273)/0.032550614492363D0/
      DATA IBD/0/
      IF(IBD.NE.0) RETURN
      IBD=1
      DO 10 I=1,273
        B(I) = A(I)
10      Y(I) = X(I)
      DO 20 I=1,96
20      LTAB(I) = KTAB(I)
      RETURN
      END

C==========================================================================
C.  Library of programs for the generation of nucleus-nucleus interactions
C.  and the study of nucleus-induced cosmic ray showers
C.
C.  September 2001  changes in FPNI, and SIGMA_INI,
C.                  new SIGMA_PP, SIGMA_PPI, SIGMA_KP  (R. Engel)
C.
C.  may  1996       small bug  corrected by Dieter Heck in NUC_CONF 
C.
C.  march 1996      small modification to the superposition code
C.
C.  February 1996   change to FPNI to give an interaction length
C.                   also  at very low energy  
C.
C.  Version 1.01  september 1995 
C.       (small corrections P.L.)
C.       the random number generator is called as S_RNDM(0)
C.  ------------------------------------------------------
C.  Version 1.00  April 1992
C.
C.  Authors:
C.
C.     J. Engel
C.     T.K Gaisser
C.     P.Lipari
C.     T. Stanev
C. 
C.  This set of routines  when used in  the simulation of cosmic ray
C.  showers have only three  "contact points" with the "external world"
C.
C.    (i) SUBROUTINE NUC_NUC_INI
C.        (no  calling arguments)         
C.         to be called once during general initialization
C.    (ii) SUBROUTINE HEAVY (IA, E0)
C.         where IA (integer) is the mass number of the projectile
C.         nucleus  and E0 (TeV) is the energy per nucleon
C.         The output (positions of first interaction for the IA
C.         nucleons of the projectile) is  contained in  the common block:
C.           COMMON /C1STNC/ XX0(60),XX(60),YY(60),AX(60),AY(60)
C.         In detail:
C.             XX0(j)   (g cm-2) =  position of interaction
C.             XX(j) (mm)    x-distance from shower axis
C.             YY(j) (mm)    y-distance from shower axis
C.             AX(j) (radiants)  Theta_x with respect to original direction
C.             AY(j) (radiants)  Theta_y with respect to original direction
C.      
C.    (iii)  FUNCTION FPNI (E,L)
C.           Interaction length in air.
C.           E (TeV) is the energy of the particle, L is the particle
C.           code (NOTE: "Sibyll" codes are used : L =1-18) 
C.           WANRNING : The nucleus-nucleus cross section
C.           tabulated in the program are "matched" to the p-Air
C.           cross section calculated  with this FUNCTION, in other words 
C.           they are both calculated with the same input pp cross section
C==========================================================================

      SUBROUTINE NUC_NUC_INI
C...Initialization for the generation of nucleus-nucleus interactions
C.  INPUT : E0 (TeV) Energy per nucleon of the beam nucleus
C........................................................................
      SAVE

      CALL NUC_GEOM_INI                       ! nucleus profiles
      CALL SIGMA_INI                          ! initialize pp cross sections

      RETURN
      END


      SUBROUTINE FRAGM1 (IA,NW, NF, IAF)
C...Nuclear Fragmentation 
C.  total dissolution of nucleus
C..........................................
      SAVE

      DIMENSION IAF(60)
      NF = IA-NW
      DO J=1,NF
         IAF(J) = 1
      ENDDO
      RETURN
      END
C->
      SUBROUTINE FRAGM2 (IA,NW, NF, IAF)
C...Nuclear Fragmentation 
C.  Spectator in one single fragment 
C..........................................
      SAVE

      DIMENSION IAF(60)
      IF (IA-NW .GT. 0)  THEN
         NF = 1
         IAF(1) = IA-NW
      ELSE
         NF = 0
      ENDIF
      RETURN
      END

C====================================================================
C...Code of fragmentation  of spectator nucleons
C.  based on Jon Engel  abrasion-ablation algorithms
C====================================================================

      BLOCK DATA FRAG_DATA
      SAVE

C...Data for the fragmentation of  nucleus  projectiles
      COMMON /FRAGMOD/A(10,10,20),AE(10,10,20),ERES(10,10),NFLAGG(10,10)
      DATA (NFLAGG(I, 1),I=1,10)  / 
     +    0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
      DATA (NFLAGG(I, 2),I=1,10)  / 
     +    0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
      DATA (NFLAGG(I, 3),I=1,10)  / 
     +    0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
      DATA (NFLAGG(I, 4),I=1,10)  / 
     +    0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
      DATA (NFLAGG(I, 5),I=1,10)  / 
     +    0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
      DATA (NFLAGG(I, 6),I=1,10)  / 
     +    0,  0,  0,  0,  0,  0,  0,  1,  1,  1 /
      DATA (NFLAGG(I, 7),I=1,10)  / 
     +    1,  1,  1,  1,  1,  1,  1,  1,  1,  1 /
      DATA (NFLAGG(I, 8),I=1,10)  / 
     +    1,  1,  1,  1,  1,  1,  1,  1,  1,  1 /
      DATA (NFLAGG(I, 9),I=1,10)  / 
     +    1,  1,  1,  1,  1,  1,  1,  1,  1,  1 /
      DATA (NFLAGG(I,10),I=1,10)  / 
     +    1,  1,  1,  1,  1,  1,  1,  1,  1,  1 /
      DATA (A(I, 1, 1),I=1,10)  / 
     +  .438E-01,.172    ,.283    ,.511    ,.715    ,.920    ,1.19    ,
     +  1.37    ,1.65    ,2.14     /
      DATA (A(I, 1, 2),I=1,10)  / 
     +  .147E-01,.249E-01,.439E-01,.592E-01,.776E-01,.886E-01,.108    ,
     +  .117    ,.126    ,.128     /
      DATA (A(I, 1, 3),I=1,10)  / 
     +  .216E-02,.627E-02,.834E-02,.108E-01,.144E-01,.152E-01,.196E-01,
     +  .200E-01,.210E-01,.224E-01 /
      DATA (A(I, 1, 4),I=1,10)  / 
     +  .593E-01,.653E-01,.116    ,.145    ,.184    ,.204    ,.234    ,
     +  .257    ,.271    ,.248     /
      DATA (A(I, 1, 5),I=1,10)  / 
     +  .000E+00,.918E-02,.362E-02,.805E-02,.436E-02,.728E-02,.466E-02,
     +  .707E-02,.932E-02,.130E-01 /
      DATA (A(I, 1, 6),I=1,10)  / 
     +  .000E+00,.180E-02,.247E-02,.208E-02,.224E-02,.214E-02,.226E-02,
     +  .233E-02,.230E-02,.194E-02 /
      DATA (A(I, 1, 7),I=1,10)  / 
     +  .000E+00,.106E-02,.703E-03,.687E-03,.739E-03,.674E-03,.819E-03,
     +  .768E-03,.756E-03,.720E-03 /
      DATA (A(I, 1, 8),I=1,10)  / 
     +  .000E+00,.000E+00,.188E-02,.130E-02,.138E-02,.117E-02,.124E-02,
     +  .119E-02,.111E-02,.829E-03 /
      DATA (A(I, 1, 9),I=1,10)  / 
     +  .000E+00,.000E+00,.302E-03,.258E-03,.249E-03,.208E-03,.248E-03,
     +  .222E-03,.210E-03,.187E-03 /
      DATA (A(I, 1,10),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.235E-03,.222E-03,.172E-03,.181E-03,
     +  .166E-03,.152E-03,.124E-03 /
      DATA (A(I, 1,11),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.238E-03,.179E-03,.145E-03,.156E-03,
     +  .138E-03,.129E-03,.111E-03 /
      DATA (A(I, 1,12),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.368E-03,.400E-03,.255E-03,.262E-03,
     +  .221E-03,.182E-03,.112E-03 /
      DATA (A(I, 1,13),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.753E-04,.712E-04,.527E-04,
     +  .537E-04,.538E-04,.487E-04 /
      DATA (A(I, 1,14),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.103E-03,.589E-04,.578E-04,
     +  .468E-04,.385E-04,.269E-04 /
      DATA (A(I, 1,15),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.444E-04,.372E-04,
     +  .318E-04,.284E-04,.218E-04 /
      DATA (A(I, 1,16),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.487E-04,.473E-04,
     +  .338E-04,.243E-04,.122E-04 /
      DATA (A(I, 1,17),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.121E-04,.117E-04,
     +  .932E-05,.792E-05,.583E-05 /
      DATA (A(I, 1,18),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.147E-04,
     +  .101E-04,.756E-05,.496E-05 /
      DATA (A(I, 1,19),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.755E-05,
     +  .612E-05,.505E-05,.341E-05 /
      DATA (A(I, 1,20),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
     +  .630E-05,.444E-05,.282E-05 /
      DATA (A(I, 2, 1),I=1,10)  / 
     +  .269    ,.510    ,.738    ,1.12    ,1.46    ,1.83    ,2.22    ,
     +  2.57    ,3.00    ,3.67     /
      DATA (A(I, 2, 2),I=1,10)  / 
     +  .121    ,.133    ,.190    ,.234    ,.293    ,.332    ,.395    ,
     +  .431    ,.468    ,.502     /
      DATA (A(I, 2, 3),I=1,10)  / 
     +  .227E-01,.374E-01,.474E-01,.578E-01,.722E-01,.794E-01,.960E-01,
     +  .102    ,.110    ,.120     /
      DATA (A(I, 2, 4),I=1,10)  / 
     +  .287    ,.196    ,.270    ,.314    ,.373    ,.408    ,.462    ,
     +  .498    ,.529    ,.523     /
      DATA (A(I, 2, 5),I=1,10)  / 
     +  .000E+00,.433E-01,.218E-01,.384E-01,.263E-01,.385E-01,.298E-01,
     +  .405E-01,.504E-01,.671E-01 /
      DATA (A(I, 2, 6),I=1,10)  / 
     +  .000E+00,.151E-01,.177E-01,.159E-01,.173E-01,.173E-01,.187E-01,
     +  .196E-01,.201E-01,.191E-01 /
      DATA (A(I, 2, 7),I=1,10)  / 
     +  .000E+00,.457E-02,.607E-02,.610E-02,.677E-02,.670E-02,.784E-02,
     +  .787E-02,.806E-02,.803E-02 /
      DATA (A(I, 2, 8),I=1,10)  / 
     +  .000E+00,.000E+00,.702E-02,.536E-02,.558E-02,.510E-02,.554E-02,
     +  .546E-02,.538E-02,.489E-02 /
      DATA (A(I, 2, 9),I=1,10)  / 
     +  .000E+00,.000E+00,.190E-02,.199E-02,.205E-02,.191E-02,.221E-02,
     +  .214E-02,.213E-02,.204E-02 /
      DATA (A(I, 2,10),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.226E-02,.219E-02,.195E-02,.208E-02,
     +  .204E-02,.203E-02,.194E-02 /
      DATA (A(I, 2,11),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.213E-02,.195E-02,.175E-02,.191E-02,
     +  .183E-02,.179E-02,.166E-02 /
      DATA (A(I, 2,12),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.588E-03,.186E-02,.137E-02,.141E-02,
     +  .128E-02,.117E-02,.947E-03 /
      DATA (A(I, 2,13),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.554E-03,.562E-03,.454E-03,
     +  .485E-03,.505E-03,.509E-03 /
      DATA (A(I, 2,14),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.490E-03,.533E-03,.531E-03,
     +  .476E-03,.437E-03,.369E-03 /
      DATA (A(I, 2,15),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.427E-03,.382E-03,
     +  .358E-03,.340E-03,.294E-03 /
      DATA (A(I, 2,16),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.239E-03,.298E-03,
     +  .238E-03,.196E-03,.134E-03 /
      DATA (A(I, 2,17),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.299E-04,.893E-04,
     +  .796E-04,.744E-04,.683E-04 /
      DATA (A(I, 2,18),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.127E-03,
     +  .107E-03,.916E-04,.720E-04 /
      DATA (A(I, 2,19),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.397E-04,
     +  .630E-04,.565E-04,.461E-04 /
      DATA (A(I, 2,20),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
     +  .511E-04,.459E-04,.402E-04 /
      DATA (A(I, 3, 1),I=1,10)  / 
     +  .708    ,1.02    ,1.41    ,1.91    ,2.42    ,3.00    ,3.53    ,
     +  4.09    ,4.71    ,5.57     /
      DATA (A(I, 3, 2),I=1,10)  / 
     +  .397    ,.410    ,.539    ,.648    ,.795    ,.910    ,1.06    ,
     +  1.17    ,1.29    ,1.42     /
      DATA (A(I, 3, 3),I=1,10)  / 
     +  .845E-01,.122    ,.157    ,.190    ,.232    ,.262    ,.307    ,
     +  .335    ,.366    ,.402     /
      DATA (A(I, 3, 4),I=1,10)  / 
     +  .210    ,.379    ,.450    ,.490    ,.574    ,.636    ,.709    ,
     +  .769    ,.820    ,.849     /
      DATA (A(I, 3, 5),I=1,10)  / 
     +  .000E+00,.102    ,.675E-01,.104    ,.858E-01,.115    ,.102    ,
     +  .129    ,.154    ,.194     /
      DATA (A(I, 3, 6),I=1,10)  / 
     +  .000E+00,.392E-01,.615E-01,.593E-01,.649E-01,.674E-01,.735E-01,
     +  .779E-01,.817E-01,.828E-01 /
      DATA (A(I, 3, 7),I=1,10)  / 
     +  .000E+00,.539E-02,.222E-01,.238E-01,.269E-01,.280E-01,.320E-01,
     +  .334E-01,.350E-01,.361E-01 /
      DATA (A(I, 3, 8),I=1,10)  / 
     +  .000E+00,.000E+00,.838E-02,.130E-01,.133E-01,.131E-01,.141E-01,
     +  .144E-01,.149E-01,.152E-01 /
      DATA (A(I, 3, 9),I=1,10)  / 
     +  .000E+00,.000E+00,.228E-02,.647E-02,.688E-02,.687E-02,.772E-02,
     +  .786E-02,.811E-02,.824E-02 /
      DATA (A(I, 3,10),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.664E-02,.828E-02,.802E-02,.845E-02,
     +  .869E-02,.902E-02,.930E-02 /
      DATA (A(I, 3,11),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.338E-02,.735E-02,.710E-02,.767E-02,
     +  .767E-02,.776E-02,.756E-02 /
      DATA (A(I, 3,12),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.280E-03,.262E-02,.349E-02,.342E-02,
     +  .322E-02,.312E-02,.291E-02 /
      DATA (A(I, 3,13),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.618E-03,.161E-02,.138E-02,
     +  .148E-02,.155E-02,.166E-02 /
      DATA (A(I, 3,14),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.313E-03,.128E-02,.161E-02,
     +  .150E-02,.144E-02,.134E-02 /
      DATA (A(I, 3,15),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.645E-03,.118E-02,
     +  .115E-02,.111E-02,.103E-02 /
      DATA (A(I, 3,16),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.117E-03,.497E-03,
     +  .581E-03,.501E-03,.401E-03 /
      DATA (A(I, 3,17),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.115E-04,.997E-04,
     +  .202E-03,.203E-03,.206E-03 /
      DATA (A(I, 3,18),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.877E-04,
     +  .242E-03,.263E-03,.226E-03 /
      DATA (A(I, 3,19),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.158E-04,
     +  .881E-04,.152E-03,.136E-03 /
      DATA (A(I, 3,20),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
     +  .358E-04,.997E-04,.117E-03 /
      DATA (A(I, 4, 1),I=1,10)  / 
     +  .945    ,1.29    ,1.40    ,1.98    ,2.73    ,3.17    ,3.77    ,
     +  4.29    ,4.78    ,5.54     /
      DATA (A(I, 4, 2),I=1,10)  / 
     +  .581    ,.599    ,.645    ,.839    ,1.10    ,1.25    ,1.47    ,
     +  1.64    ,1.78    ,1.99     /
      DATA (A(I, 4, 3),I=1,10)  / 
     +  .127    ,.182    ,.202    ,.264    ,.344    ,.387    ,.455    ,
     +  .504    ,.549    ,.611     /
      DATA (A(I, 4, 4),I=1,10)  / 
     +  .183    ,.464    ,.351    ,.444    ,.642    ,.659    ,.772    ,
     +  .830    ,.882    ,.930     /
      DATA (A(I, 4, 5),I=1,10)  / 
     +  .000E+00,.122    ,.803E-01,.136    ,.134    ,.173    ,.164    ,
     +  .203    ,.239    ,.300     /
      DATA (A(I, 4, 6),I=1,10)  / 
     +  .000E+00,.393E-01,.766E-01,.872E-01,.108    ,.111    ,.123    ,
     +  .132    ,.139    ,.145     /
      DATA (A(I, 4, 7),I=1,10)  / 
     +  .000E+00,.416E-02,.289E-01,.360E-01,.454E-01,.477E-01,.549E-01,
     +  .583E-01,.618E-01,.654E-01 /
      DATA (A(I, 4, 8),I=1,10)  / 
     +  .000E+00,.000E+00,.761E-02,.157E-01,.214E-01,.205E-01,.233E-01,
     +  .241E-01,.255E-01,.271E-01 /
      DATA (A(I, 4, 9),I=1,10)  / 
     +  .000E+00,.000E+00,.238E-02,.803E-02,.123E-01,.123E-01,.140E-01,
     +  .145E-01,.153E-01,.160E-01 /
      DATA (A(I, 4,10),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.695E-02,.150E-01,.154E-01,.166E-01,
     +  .172E-01,.181E-01,.192E-01 /
      DATA (A(I, 4,11),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.355E-02,.104E-01,.143E-01,.156E-01,
     +  .158E-01,.164E-01,.165E-01 /
      DATA (A(I, 4,12),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.112E-03,.276E-02,.568E-02,.736E-02,
     +  .684E-02,.691E-02,.661E-02 /
      DATA (A(I, 4,13),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.740E-03,.222E-02,.339E-02,
     +  .352E-02,.382E-02,.409E-02 /
      DATA (A(I, 4,14),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.369E-03,.160E-02,.322E-02,
     +  .375E-02,.375E-02,.355E-02 /
      DATA (A(I, 4,15),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.750E-03,.190E-02,
     +  .298E-02,.319E-02,.299E-02 /
      DATA (A(I, 4,16),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.260E-03,.673E-03,
     +  .117E-02,.156E-02,.126E-02 /
      DATA (A(I, 4,17),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.283E-05,.131E-03,
     +  .363E-03,.618E-03,.690E-03 /
      DATA (A(I, 4,18),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.205E-03,
     +  .378E-03,.709E-03,.844E-03 /
      DATA (A(I, 4,19),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.654E-05,
     +  .150E-03,.341E-03,.527E-03 /
      DATA (A(I, 4,20),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
     +  .957E-04,.197E-03,.406E-03 /
      DATA (A(I, 5, 1),I=1,10)  / 
     +  1.16    ,1.70    ,2.19    ,2.79    ,3.33    ,3.90    ,4.49    ,
     +  5.07    ,5.66    ,6.38     /
      DATA (A(I, 5, 2),I=1,10)  / 
     +  .779    ,.899    ,1.09    ,1.28    ,1.51    ,1.71    ,1.96    ,
     +  2.18    ,2.39    ,2.62     /
      DATA (A(I, 5, 3),I=1,10)  / 
     +  .167    ,.263    ,.334    ,.408    ,.482    ,.548    ,.632    ,
     +  .700    ,.767    ,.840     /
      DATA (A(I, 5, 4),I=1,10)  / 
     +  .203    ,.565    ,.845    ,.867    ,.906    ,.961    ,1.08    ,
     +  1.13    ,1.21    ,1.25     /
      DATA (A(I, 5, 5),I=1,10)  / 
     +  .000E+00,.129    ,.152    ,.237    ,.208    ,.268    ,.258    ,
     +  .312    ,.368    ,.450     /
      DATA (A(I, 5, 6),I=1,10)  / 
     +  .000E+00,.460E-01,.126    ,.174    ,.182    ,.188    ,.208    ,
     +  .219    ,.233    ,.239     /
      DATA (A(I, 5, 7),I=1,10)  / 
     +  .000E+00,.289E-02,.380E-01,.611E-01,.788E-01,.845E-01,.974E-01,
     +  .103    ,.111    ,.117     /
      DATA (A(I, 5, 8),I=1,10)  / 
     +  .000E+00,.000E+00,.137E-01,.223E-01,.374E-01,.436E-01,.488E-01,
     +  .488E-01,.524E-01,.547E-01 /
      DATA (A(I, 5, 9),I=1,10)  / 
     +  .000E+00,.000E+00,.162E-02,.114E-01,.198E-01,.263E-01,.315E-01,
     +  .323E-01,.348E-01,.364E-01 /
      DATA (A(I, 5,10),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.149E-01,.240E-01,.320E-01,.428E-01,
     +  .436E-01,.469E-01,.493E-01 /
      DATA (A(I, 5,11),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.562E-02,.194E-01,.290E-01,.408E-01,
     +  .460E-01,.492E-01,.500E-01 /
      DATA (A(I, 5,12),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.476E-04,.106E-01,.134E-01,.191E-01,
     +  .227E-01,.264E-01,.253E-01 /
      DATA (A(I, 5,13),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.281E-02,.679E-02,.879E-02,
     +  .123E-01,.165E-01,.190E-01 /
      DATA (A(I, 5,14),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.542E-04,.847E-02,.125E-01,
     +  .144E-01,.173E-01,.192E-01 /
      DATA (A(I, 5,15),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.652E-02,.982E-02,
     +  .129E-01,.159E-01,.192E-01 /
      DATA (A(I, 5,16),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.109E-03,.688E-02,
     +  .751E-02,.845E-02,.905E-02 /
      DATA (A(I, 5,17),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.823E-06,.237E-02,
     +  .318E-02,.446E-02,.569E-02 /
      DATA (A(I, 5,18),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.604E-03,
     +  .610E-02,.673E-02,.827E-02 /
      DATA (A(I, 5,19),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.716E-06,
     +  .412E-02,.519E-02,.617E-02 /
      DATA (A(I, 5,20),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
     +  .710E-03,.543E-02,.674E-02 /
      DATA (A(I, 6, 1),I=1,10)  / 
     +  1.36    ,2.08    ,2.67    ,3.30    ,3.94    ,4.62    ,5.18    ,
     +  3.60    ,3.64    ,3.95     /
      DATA (A(I, 6, 2),I=1,10)  / 
     +  1.07    ,1.33    ,1.58    ,1.82    ,2.10    ,2.44    ,2.74    ,
     +  1.78    ,1.73    ,1.80     /
      DATA (A(I, 6, 3),I=1,10)  / 
     +  .158    ,.276    ,.402    ,.506    ,.609    ,.700    ,.802    ,
     +  .638    ,.629    ,.658     /
      DATA (A(I, 6, 4),I=1,10)  / 
     +  .308    ,.739    ,1.02    ,1.12    ,1.26    ,1.35    ,1.57    ,
     +  1.94    ,1.71    ,1.55     /
      DATA (A(I, 6, 5),I=1,10)  / 
     +  .000E+00,.217    ,.183    ,.324    ,.276    ,.395    ,.393    ,
     +  .558    ,.602    ,.681     /
      DATA (A(I, 6, 6),I=1,10)  / 
     +  .000E+00,.658E-01,.251    ,.267    ,.299    ,.326    ,.386    ,
     +  .452    ,.475    ,.409     /
      DATA (A(I, 6, 7),I=1,10)  / 
     +  .000E+00,.198E-02,.774E-01,.136    ,.149    ,.164    ,.187    ,
     +  .210    ,.238    ,.256     /
      DATA (A(I, 6, 8),I=1,10)  / 
     +  .000E+00,.000E+00,.290E-01,.122    ,.139    ,.128    ,.129    ,
     +  .137    ,.147    ,.167     /
      DATA (A(I, 6, 9),I=1,10)  / 
     +  .000E+00,.000E+00,.699E-03,.617E-01,.750E-01,.801E-01,.905E-01,
     +  .974E-01,.105    ,.122     /
      DATA (A(I, 6,10),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.310E-01,.112    ,.127    ,.140    ,
     +  .143    ,.155    ,.176     /
      DATA (A(I, 6,11),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.277E-02,.889E-01,.143    ,.150    ,
     +  .175    ,.184    ,.208     /
      DATA (A(I, 6,12),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.202E-04,.343E-01,.959E-01,.109    ,
     +  .115    ,.112    ,.116     /
      DATA (A(I, 6,13),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.186E-02,.435E-01,.512E-01,
     +  .744E-01,.856E-01,.103     /
      DATA (A(I, 6,14),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.144E-04,.427E-01,.786E-01,
     +  .911E-01,.993E-01,.108     /
      DATA (A(I, 6,15),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.466E-02,.518E-01,
     +  .848E-01,.109    ,.119     /
      DATA (A(I, 6,16),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.655E-05,.330E-01,
     +  .586E-01,.617E-01,.594E-01 /
      DATA (A(I, 6,17),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.228E-06,.328E-02,
     +  .190E-01,.301E-01,.454E-01 /
      DATA (A(I, 6,18),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.218E-04,
     +  .272E-01,.501E-01,.707E-01 /
      DATA (A(I, 6,19),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.146E-06,
     +  .441E-02,.378E-01,.556E-01 /
      DATA (A(I, 6,20),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
     +  .160E-03,.204E-01,.679E-01 /
      DATA (A(I, 7, 1),I=1,10)  / 
     +  .522    ,.862    ,1.14    ,1.40    ,1.70    ,1.94    ,2.26    ,
     +  2.48    ,2.72    ,3.95     /
      DATA (A(I, 7, 2),I=1,10)  / 
     +  .314    ,.450    ,.588    ,.692    ,.834    ,.936    ,1.09    ,
     +  1.18    ,1.28    ,1.80     /
      DATA (A(I, 7, 3),I=1,10)  / 
     +  .814E-01,.147    ,.189    ,.226    ,.272    ,.302    ,.351    ,
     +  .378    ,.406    ,.658     /
      DATA (A(I, 7, 4),I=1,10)  / 
     +  .252    ,.864    ,1.01    ,.851    ,.837    ,.774    ,.763    ,
     +  .757    ,.748    ,1.55     /
      DATA (A(I, 7, 5),I=1,10)  / 
     +  .000E+00,.225    ,.180    ,.276    ,.193    ,.240    ,.190    ,
     +  .228    ,.259    ,.681     /
      DATA (A(I, 7, 6),I=1,10)  / 
     +  .000E+00,.485E-01,.272    ,.273    ,.253    ,.216    ,.206    ,
     +  .197    ,.191    ,.409     /
      DATA (A(I, 7, 7),I=1,10)  / 
     +  .000E+00,.137E-02,.752E-01,.137    ,.152    ,.134    ,.125    ,
     +  .119    ,.116    ,.256     /
      DATA (A(I, 7, 8),I=1,10)  / 
     +  .000E+00,.000E+00,.220E-01,.155    ,.175    ,.155    ,.116    ,
     +  .977E-01,.858E-01,.167     /
      DATA (A(I, 7, 9),I=1,10)  / 
     +  .000E+00,.000E+00,.326E-03,.695E-01,.881E-01,.106    ,.897E-01,
     +  .782E-01,.706E-01,.122     /
      DATA (A(I, 7,10),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.261E-01,.124    ,.131    ,.156    ,
     +  .141    ,.121    ,.176     /
      DATA (A(I, 7,11),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.785E-03,.864E-01,.130    ,.170    ,
     +  .182    ,.172    ,.208     /
      DATA (A(I, 7,12),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.896E-05,.225E-01,.105    ,.126    ,
     +  .126    ,.135    ,.116     /
      DATA (A(I, 7,13),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.542E-03,.427E-01,.553E-01,
     +  .744E-01,.980E-01,.103     /
      DATA (A(I, 7,14),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.515E-05,.377E-01,.831E-01,
     +  .985E-01,.104    ,.108     /
      DATA (A(I, 7,15),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.285E-02,.495E-01,
     +  .871E-01,.106    ,.119     /
      DATA (A(I, 7,16),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.110E-05,.284E-01,
     +  .588E-01,.657E-01,.594E-01 /
      DATA (A(I, 7,17),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.722E-07,.176E-02,
     +  .170E-01,.305E-01,.454E-01 /
      DATA (A(I, 7,18),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.148E-05,
     +  .213E-01,.492E-01,.707E-01 /
      DATA (A(I, 7,19),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.323E-07,
     +  .722E-02,.359E-01,.556E-01 /
      DATA (A(I, 7,20),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
     +  .461E-05,.155E-01,.679E-01 /
      DATA (A(I, 8, 1),I=1,10)  / 
     +  .630    ,.974    ,1.29    ,1.58    ,1.89    ,2.16    ,2.49    ,
     +  2.75    ,3.02    ,3.95     /
      DATA (A(I, 8, 2),I=1,10)  / 
     +  .328    ,.459    ,.613    ,.735    ,.879    ,.994    ,1.15    ,
     +  1.27    ,1.38    ,1.80     /
      DATA (A(I, 8, 3),I=1,10)  / 
     +  .748E-01,.121    ,.164    ,.197    ,.235    ,.265    ,.310    ,
     +  .339    ,.370    ,.658     /
      DATA (A(I, 8, 4),I=1,10)  / 
     +  .194    ,.211    ,.337    ,.344    ,.339    ,.351    ,.390    ,
     +  .419    ,.442    ,1.55     /
      DATA (A(I, 8, 5),I=1,10)  / 
     +  .000E+00,.869E-01,.725E-01,.113    ,.810E-01,.106    ,.951E-01,
     +  .120    ,.143    ,.681     /
      DATA (A(I, 8, 6),I=1,10)  / 
     +  .000E+00,.288E-01,.102    ,.922E-01,.857E-01,.845E-01,.932E-01,
     +  .983E-01,.102    ,.409     /
      DATA (A(I, 8, 7),I=1,10)  / 
     +  .000E+00,.668E-03,.533E-01,.575E-01,.493E-01,.482E-01,.539E-01,
     +  .558E-01,.582E-01,.256     /
      DATA (A(I, 8, 8),I=1,10)  / 
     +  .000E+00,.000E+00,.205E-01,.808E-01,.510E-01,.409E-01,.406E-01,
     +  .394E-01,.389E-01,.167     /
      DATA (A(I, 8, 9),I=1,10)  / 
     +  .000E+00,.000E+00,.999E-04,.647E-01,.385E-01,.325E-01,.325E-01,
     +  .316E-01,.314E-01,.122     /
      DATA (A(I, 8,10),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.169E-01,.834E-01,.611E-01,.565E-01,
     +  .533E-01,.519E-01,.176     /
      DATA (A(I, 8,11),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.107E-03,.769E-01,.922E-01,.805E-01,
     +  .745E-01,.711E-01,.208     /
      DATA (A(I, 8,12),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.180E-05,.143E-01,.983E-01,.775E-01,
     +  .627E-01,.541E-01,.116     /
      DATA (A(I, 8,13),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.157E-04,.346E-01,.507E-01,
     +  .479E-01,.455E-01,.103     /
      DATA (A(I, 8,14),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.752E-06,.248E-01,.721E-01,
     +  .728E-01,.611E-01,.108     /
      DATA (A(I, 8,15),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.686E-04,.356E-01,
     +  .731E-01,.791E-01,.119     /
      DATA (A(I, 8,16),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.838E-07,.151E-01,
     +  .470E-01,.567E-01,.594E-01 /
      DATA (A(I, 8,17),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.759E-08,.400E-04,
     +  .193E-01,.313E-01,.454E-01 /
      DATA (A(I, 8,18),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.385E-07,
     +  .921E-02,.353E-01,.707E-01 /
      DATA (A(I, 8,19),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.219E-08,
     +  .348E-03,.226E-01,.556E-01 /
      DATA (A(I, 8,20),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
     +  .212E-07,.149E-01,.679E-01 /
      DATA (A(I, 9, 1),I=1,10)  / 
     +  .736    ,1.13    ,1.49    ,1.82    ,2.20    ,2.49    ,2.86    ,
     +  3.17    ,3.49    ,3.95     /
      DATA (A(I, 9, 2),I=1,10)  / 
     +  .339    ,.492    ,.658    ,.789    ,.958    ,1.08    ,1.25    ,
     +  1.37    ,1.50    ,1.80     /
      DATA (A(I, 9, 3),I=1,10)  / 
     +  .680E-01,.110    ,.150    ,.180    ,.222    ,.247    ,.289    ,
     +  .318    ,.349    ,.658     /
      DATA (A(I, 9, 4),I=1,10)  / 
     +  .110    ,.104    ,.157    ,.156    ,.210    ,.205    ,.246    ,
     +  .274    ,.300    ,1.55     /
      DATA (A(I, 9, 5),I=1,10)  / 
     +  .000E+00,.379E-01,.347E-01,.477E-01,.486E-01,.576E-01,.569E-01,
     +  .732E-01,.893E-01,.681     /
      DATA (A(I, 9, 6),I=1,10)  / 
     +  .000E+00,.223E-01,.354E-01,.312E-01,.436E-01,.400E-01,.489E-01,
     +  .548E-01,.600E-01,.409     /
      DATA (A(I, 9, 7),I=1,10)  / 
     +  .000E+00,.338E-03,.149E-01,.142E-01,.215E-01,.188E-01,.248E-01,
     +  .278E-01,.307E-01,.256     /
      DATA (A(I, 9, 8),I=1,10)  / 
     +  .000E+00,.000E+00,.553E-02,.862E-02,.150E-01,.106E-01,.145E-01,
     +  .165E-01,.181E-01,.167     /
      DATA (A(I, 9, 9),I=1,10)  / 
     +  .000E+00,.000E+00,.375E-04,.641E-02,.111E-01,.792E-02,.112E-01,
     +  .127E-01,.140E-01,.122     /
      DATA (A(I, 9,10),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.112E-01,.200E-01,.127E-01,.176E-01,
     +  .200E-01,.220E-01,.176     /
      DATA (A(I, 9,11),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.244E-04,.261E-01,.162E-01,.232E-01,
     +  .263E-01,.287E-01,.208     /
      DATA (A(I, 9,12),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.455E-06,.635E-02,.121E-01,.186E-01,
     +  .201E-01,.207E-01,.116     /
      DATA (A(I, 9,13),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.146E-05,.922E-02,.116E-01,
     +  .145E-01,.165E-01,.103     /
      DATA (A(I, 9,14),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.135E-06,.128E-01,.202E-01,
     +  .215E-01,.220E-01,.108     /
      DATA (A(I, 9,15),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.237E-05,.229E-01,
     +  .259E-01,.271E-01,.119     /
      DATA (A(I, 9,16),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.100E-07,.534E-02,
     +  .210E-01,.193E-01,.594E-01 /
      DATA (A(I, 9,17),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.915E-09,.847E-06,
     +  .119E-01,.125E-01,.454E-01 /
      DATA (A(I, 9,18),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.298E-08,
     +  .101E-01,.242E-01,.707E-01 /
      DATA (A(I, 9,19),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.196E-09,
     +  .243E-05,.234E-01,.556E-01 /
      DATA (A(I, 9,20),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
     +  .575E-09,.364E-02,.679E-01 /
      DATA (A(I,10, 1),I=1,10)  / 
     +  .959    ,1.46    ,1.92    ,2.34    ,2.80    ,3.24    ,3.64    ,
     +  4.05    ,4.48    ,3.95     /
      DATA (A(I,10, 2),I=1,10)  / 
     +  .343    ,.516    ,.692    ,.836    ,1.01    ,1.16    ,1.31    ,
     +  1.46    ,1.61    ,1.80     /
      DATA (A(I,10, 3),I=1,10)  / 
     +  .512E-01,.837E-01,.115    ,.138    ,.169    ,.195    ,.220    ,
     +  .245    ,.270    ,.658     /
      DATA (A(I,10, 4),I=1,10)  / 
     +  .274E-01,.361E-01,.510E-01,.562E-01,.703E-01,.828E-01,.877E-01,
     +  .996E-01,.111    ,1.55     /
      DATA (A(I,10, 5),I=1,10)  / 
     +  .000E+00,.850E-02,.875E-02,.118E-01,.124E-01,.170E-01,.154E-01,
     +  .194E-01,.237E-01,.681     /
      DATA (A(I,10, 6),I=1,10)  / 
     +  .000E+00,.345E-02,.519E-02,.533E-02,.691E-02,.842E-02,.844E-02,
     +  .987E-02,.113E-01,.409     /
      DATA (A(I,10, 7),I=1,10)  / 
     +  .000E+00,.722E-04,.130E-02,.135E-02,.189E-02,.240E-02,.235E-02,
     +  .281E-02,.331E-02,.256     /
      DATA (A(I,10, 8),I=1,10)  / 
     +  .000E+00,.000E+00,.283E-03,.272E-03,.394E-03,.557E-03,.480E-03,
     +  .616E-03,.775E-03,.167     /
      DATA (A(I,10, 9),I=1,10)  / 
     +  .000E+00,.000E+00,.457E-05,.122E-03,.192E-03,.275E-03,.225E-03,
     +  .292E-03,.373E-03,.122     /
      DATA (A(I,10,10),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.119E-03,.185E-03,.278E-03,.201E-03,
     +  .274E-03,.364E-03,.176     /
      DATA (A(I,10,11),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.140E-05,.129E-03,.200E-03,.137E-03,
     +  .188E-03,.252E-03,.208     /
      DATA (A(I,10,12),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.207E-07,.307E-04,.518E-04,.278E-04,
     +  .421E-04,.608E-04,.116     /
      DATA (A(I,10,13),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.306E-07,.252E-04,.111E-04,
     +  .188E-04,.295E-04,.103     /
      DATA (A(I,10,14),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.321E-08,.220E-04,.104E-04,
     +  .162E-04,.243E-04,.108     /
      DATA (A(I,10,15),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.770E-08,.632E-05,
     +  .105E-04,.162E-04,.119     /
      DATA (A(I,10,16),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.117E-09,.199E-05,
     +  .321E-05,.492E-05,.594E-01 /
      DATA (A(I,10,17),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.888E-11,.323E-09,
     +  .106E-05,.192E-05,.454E-01 /
      DATA (A(I,10,18),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.174E-10,
     +  .131E-05,.218E-05,.707E-01 /
      DATA (A(I,10,19),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.994E-12,
     +  .233E-09,.104E-05,.556E-01 /
      DATA (A(I,10,20),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
     +  .144E-11,.724E-06,.679E-01 /
      DATA (AE(I, 1, 1),I=1,10)  / 
     +  7.27    ,6.29    ,7.76    ,6.70    ,8.17    ,7.34    ,8.70    ,
     +  8.02    ,7.37    ,6.18     /
      DATA (AE(I, 1, 2),I=1,10)  / 
     +  7.41    ,7.52    ,8.14    ,8.20    ,8.96    ,9.05    ,9.96    ,
     +  10.0    ,10.1    ,9.86     /
      DATA (AE(I, 1, 3),I=1,10)  / 
     +  7.72    ,7.69    ,9.17    ,8.99    ,10.6    ,10.5    ,12.1    ,
     +  12.1    ,12.0    ,11.5     /
      DATA (AE(I, 1, 4),I=1,10)  / 
     +  7.90    ,8.48    ,9.50    ,9.94    ,10.8    ,11.4    ,12.2    ,
     +  12.8    ,13.3    ,13.8     /
      DATA (AE(I, 1, 5),I=1,10)  / 
     +  .000E+00,8.52    ,9.59    ,10.1    ,11.1    ,11.8    ,12.7    ,
     +  13.3    ,13.8    ,14.4     /
      DATA (AE(I, 1, 6),I=1,10)  / 
     +  .000E+00,9.00    ,10.7    ,11.7    ,13.2    ,14.2    ,15.6    ,
     +  16.5    ,17.3    ,18.0     /
      DATA (AE(I, 1, 7),I=1,10)  / 
     +  .000E+00,9.01    ,11.1    ,11.9    ,14.3    ,15.0    ,17.4    ,
     +  18.0    ,18.6    ,18.8     /
      DATA (AE(I, 1, 8),I=1,10)  / 
     +  .000E+00,.000E+00,11.2    ,12.4    ,14.5    ,15.7    ,17.6    ,
     +  18.8    ,19.9    ,20.9     /
      DATA (AE(I, 1, 9),I=1,10)  / 
     +  .000E+00,.000E+00,11.4    ,12.7    ,15.5    ,16.6    ,19.3    ,
     +  20.2    ,21.1    ,21.7     /
      DATA (AE(I, 1,10),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,13.2    ,15.8    ,17.3    ,19.9    ,
     +  21.2    ,22.4    ,23.2     /
      DATA (AE(I, 1,11),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,13.2    ,16.3    ,17.8    ,20.8    ,
     +  22.1    ,23.3    ,24.2     /
      DATA (AE(I, 1,12),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,13.4    ,16.2    ,18.2    ,21.0    ,
     +  22.8    ,24.4    ,25.9     /
      DATA (AE(I, 1,13),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,16.5    ,18.4    ,21.6    ,
     +  23.2    ,24.8    ,26.2     /
      DATA (AE(I, 1,14),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,16.7    ,19.0    ,22.3    ,
     +  24.3    ,26.1    ,27.4     /
      DATA (AE(I, 1,15),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,19.1    ,22.8    ,
     +  24.7    ,26.6    ,28.2     /
      DATA (AE(I, 1,16),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,19.2    ,23.0    ,
     +  25.3    ,27.5    ,29.5     /
      DATA (AE(I, 1,17),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,19.6    ,23.3    ,
     +  25.6    ,27.8    ,29.6     /
      DATA (AE(I, 1,18),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,23.6    ,
     +  26.2    ,28.5    ,30.4     /
      DATA (AE(I, 1,19),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,23.7    ,
     +  26.3    ,28.8    ,31.0     /
      DATA (AE(I, 1,20),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
     +  26.5    ,29.2    ,31.5     /
      DATA (AE(I, 2, 1),I=1,10)  / 
     +  8.74    ,8.16    ,9.25    ,8.45    ,9.46    ,8.90    ,9.83    ,
     +  9.38    ,8.96    ,8.15     /
      DATA (AE(I, 2, 2),I=1,10)  / 
     +  8.96    ,9.30    ,9.95    ,10.0    ,10.8    ,10.9    ,11.7    ,
     +  11.8    ,11.9    ,11.8     /
      DATA (AE(I, 2, 3),I=1,10)  / 
     +  9.44    ,9.66    ,11.0    ,11.0    ,12.3    ,12.5    ,13.7    ,
     +  13.9    ,14.0    ,13.8     /
      DATA (AE(I, 2, 4),I=1,10)  / 
     +  8.86    ,9.81    ,10.8    ,11.2    ,12.0    ,12.6    ,13.4    ,
     +  14.0    ,14.5    ,15.1     /
      DATA (AE(I, 2, 5),I=1,10)  / 
     +  .000E+00,10.2    ,11.4    ,12.0    ,12.9    ,13.6    ,14.5    ,
     +  15.1    ,15.7    ,16.3     /
      DATA (AE(I, 2, 6),I=1,10)  / 
     +  .000E+00,10.7    ,12.5    ,13.5    ,15.1    ,16.0    ,17.5    ,
     +  18.3    ,19.2    ,19.9     /
      DATA (AE(I, 2, 7),I=1,10)  / 
     +  .000E+00,11.5    ,12.9    ,13.9    ,16.1    ,17.0    ,19.1    ,
     +  19.8    ,20.6    ,21.0     /
      DATA (AE(I, 2, 8),I=1,10)  / 
     +  .000E+00,.000E+00,12.4    ,13.8    ,15.9    ,17.2    ,19.1    ,
     +  20.3    ,21.4    ,22.3     /
      DATA (AE(I, 2, 9),I=1,10)  / 
     +  .000E+00,.000E+00,13.4    ,14.5    ,17.1    ,18.3    ,20.9    ,
     +  21.9    ,23.0    ,23.7     /
      DATA (AE(I, 2,10),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,14.9    ,17.5    ,19.1    ,21.6    ,
     +  22.9    ,24.1    ,25.0     /
      DATA (AE(I, 2,11),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,15.0    ,18.0    ,19.6    ,22.4    ,
     +  23.8    ,25.2    ,26.2     /
      DATA (AE(I, 2,12),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,16.2    ,17.3    ,19.4    ,22.2    ,
     +  24.0    ,25.7    ,27.2     /
      DATA (AE(I, 2,13),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,17.8    ,19.8    ,22.9    ,
     +  24.6    ,26.2    ,27.7     /
      DATA (AE(I, 2,14),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,19.1    ,20.4    ,23.7    ,
     +  25.7    ,27.6    ,29.1     /
      DATA (AE(I, 2,15),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,20.5    ,24.1    ,
     +  26.1    ,28.1    ,29.9     /
      DATA (AE(I, 2,16),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,20.9    ,23.9    ,
     +  26.4    ,28.7    ,30.7     /
      DATA (AE(I, 2,17),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,22.4    ,24.2    ,
     +  26.7    ,29.0    ,30.9     /
      DATA (AE(I, 2,18),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,24.8    ,
     +  27.3    ,29.7    ,31.8     /
      DATA (AE(I, 2,19),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,26.1    ,
     +  27.3    ,29.9    ,32.3     /
      DATA (AE(I, 2,20),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
     +  27.4    ,30.1    ,32.6     /
      DATA (AE(I, 3, 1),I=1,10)  / 
     +  11.0    ,11.0    ,11.7    ,11.3    ,11.9    ,11.4    ,12.1    ,
     +  11.7    ,11.5    ,11.0     /
      DATA (AE(I, 3, 2),I=1,10)  / 
     +  11.2    ,12.0    ,12.7    ,12.9    ,13.6    ,13.7    ,14.4    ,
     +  14.6    ,14.7    ,14.6     /
      DATA (AE(I, 3, 3),I=1,10)  / 
     +  12.1    ,12.6    ,13.7    ,13.9    ,15.0    ,15.2    ,16.3    ,
     +  16.5    ,16.7    ,16.7     /
      DATA (AE(I, 3, 4),I=1,10)  / 
     +  12.6    ,11.3    ,12.4    ,13.0    ,13.8    ,14.2    ,15.0    ,
     +  15.6    ,16.1    ,16.6     /
      DATA (AE(I, 3, 5),I=1,10)  / 
     +  .000E+00,12.6    ,13.7    ,14.4    ,15.3    ,16.0    ,16.8    ,
     +  17.5    ,18.1    ,18.6     /
      DATA (AE(I, 3, 6),I=1,10)  / 
     +  .000E+00,14.0    ,14.6    ,15.8    ,17.4    ,18.4    ,19.8    ,
     +  20.6    ,21.5    ,22.2     /
      DATA (AE(I, 3, 7),I=1,10)  / 
     +  .000E+00,16.0    ,15.2    ,16.3    ,18.3    ,19.3    ,21.1    ,
     +  22.0    ,22.8    ,23.5     /
      DATA (AE(I, 3, 8),I=1,10)  / 
     +  .000E+00,.000E+00,15.6    ,15.1    ,17.2    ,18.6    ,20.6    ,
     +  21.8    ,22.9    ,23.8     /
      DATA (AE(I, 3, 9),I=1,10)  / 
     +  .000E+00,.000E+00,17.8    ,16.3    ,18.8    ,20.1    ,22.5    ,
     +  23.6    ,24.7    ,25.6     /
      DATA (AE(I, 3,10),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,17.5    ,19.0    ,20.7    ,23.1    ,
     +  24.5    ,25.8    ,26.8     /
      DATA (AE(I, 3,11),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,19.2    ,19.4    ,21.1    ,23.8    ,
     +  25.4    ,26.8    ,28.0     /
      DATA (AE(I, 3,12),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,20.7    ,19.6    ,19.7    ,22.4    ,
     +  24.4    ,26.2    ,27.9     /
      DATA (AE(I, 3,13),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,21.6    ,20.4    ,23.2    ,
     +  25.1    ,26.9    ,28.5     /
      DATA (AE(I, 3,14),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,23.5    ,22.0    ,23.8    ,
     +  26.1    ,28.1    ,29.9     /
      DATA (AE(I, 3,15),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,23.7    ,24.2    ,
     +  26.3    ,28.5    ,30.4     /
      DATA (AE(I, 3,16),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,25.4    ,24.8    ,
     +  25.6    ,28.1    ,30.5     /
      DATA (AE(I, 3,17),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,26.9    ,26.8    ,
     +  26.1    ,28.4    ,30.8     /
      DATA (AE(I, 3,18),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,28.8    ,
     +  27.6    ,29.0    ,31.5     /
      DATA (AE(I, 3,19),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,30.5    ,
     +  29.2    ,28.9    ,31.5     /
      DATA (AE(I, 3,20),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
     +  31.0    ,30.0    ,31.7     /
      DATA (AE(I, 4, 1),I=1,10)  / 
     +  13.0    ,13.2    ,14.8    ,14.2    ,14.2    ,14.1    ,14.5    ,
     +  14.4    ,14.3    ,14.0     /
      DATA (AE(I, 4, 2),I=1,10)  / 
     +  13.5    ,14.5    ,16.1    ,15.9    ,16.0    ,16.3    ,16.8    ,
     +  17.0    ,17.1    ,17.2     /
      DATA (AE(I, 4, 3),I=1,10)  / 
     +  14.9    ,15.3    ,17.2    ,17.1    ,17.5    ,17.8    ,18.6    ,
     +  18.9    ,19.1    ,19.3     /
      DATA (AE(I, 4, 4),I=1,10)  / 
     +  15.1    ,13.5    ,16.4    ,16.7    ,16.4    ,17.3    ,17.8    ,
     +  18.5    ,19.0    ,19.6     /
      DATA (AE(I, 4, 5),I=1,10)  / 
     +  .000E+00,15.6    ,17.5    ,17.7    ,17.8    ,18.6    ,19.2    ,
     +  19.9    ,20.3    ,21.1     /
      DATA (AE(I, 4, 6),I=1,10)  / 
     +  .000E+00,18.0    ,18.4    ,19.2    ,19.8    ,20.9    ,22.0    ,
     +  23.1    ,23.6    ,24.7     /
      DATA (AE(I, 4, 7),I=1,10)  / 
     +  .000E+00,27.4    ,19.1    ,19.8    ,20.7    ,21.8    ,23.2    ,
     +  24.4    ,24.9    ,25.9     /
      DATA (AE(I, 4, 8),I=1,10)  / 
     +  .000E+00,.000E+00,18.9    ,18.9    ,19.3    ,21.1    ,22.5    ,
     +  24.0    ,24.7    ,26.0     /
      DATA (AE(I, 4, 9),I=1,10)  / 
     +  .000E+00,.000E+00,21.1    ,19.7    ,20.7    ,22.3    ,24.0    ,
     +  25.6    ,26.3    ,27.7     /
      DATA (AE(I, 4,10),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,21.0    ,21.1    ,22.9    ,24.6    ,
     +  26.5    ,27.3    ,29.0     /
      DATA (AE(I, 4,11),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,21.3    ,22.4    ,23.1    ,25.0    ,
     +  27.1    ,27.9    ,29.8     /
      DATA (AE(I, 4,12),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,36.6    ,21.5    ,22.2    ,23.1    ,
     +  25.6    ,26.8    ,29.1     /
      DATA (AE(I, 4,13),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,22.9    ,23.1    ,23.7    ,
     +  26.2    ,27.3    ,29.6     /
      DATA (AE(I, 4,14),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,30.5    ,23.6    ,25.0    ,
     +  26.9    ,28.2    ,30.7     /
      DATA (AE(I, 4,15),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,25.4    ,26.2    ,
     +  27.2    ,28.3    ,31.0     /
      DATA (AE(I, 4,16),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,24.5    ,25.9    ,
     +  27.4    ,27.6    ,30.7     /
      DATA (AE(I, 4,17),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,43.3    ,28.4    ,
     +  27.5    ,27.9    ,30.9     /
      DATA (AE(I, 4,18),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,27.2    ,
     +  29.1    ,29.0    ,31.4     /
      DATA (AE(I, 4,19),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,51.3    ,
     +  30.6    ,29.5    ,31.4     /
      DATA (AE(I, 4,20),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
     +  28.8    ,30.6    ,32.4     /
      DATA (AE(I, 5, 1),I=1,10)  / 
     +  15.0    ,14.9    ,15.5    ,15.4    ,15.9    ,15.8    ,16.2    ,
     +  16.2    ,16.1    ,15.9     /
      DATA (AE(I, 5, 2),I=1,10)  / 
     +  15.4    ,16.1    ,17.0    ,17.4    ,18.0    ,18.2    ,18.7    ,
     +  18.9    ,19.0    ,19.1     /
      DATA (AE(I, 5, 3),I=1,10)  / 
     +  17.1    ,17.2    ,18.3    ,18.7    ,19.3    ,19.6    ,20.3    ,
     +  20.6    ,20.8    ,20.9     /
      DATA (AE(I, 5, 4),I=1,10)  / 
     +  14.7    ,14.8    ,15.0    ,16.0    ,17.0    ,17.7    ,18.1    ,
     +  19.0    ,19.4    ,20.0     /
      DATA (AE(I, 5, 5),I=1,10)  / 
     +  .000E+00,16.7    ,17.6    ,18.1    ,18.6    ,19.2    ,19.7    ,
     +  20.4    ,20.8    ,21.2     /
      DATA (AE(I, 5, 6),I=1,10)  / 
     +  .000E+00,17.8    ,18.2    ,19.2    ,20.0    ,21.0    ,21.9    ,
     +  23.0    ,23.6    ,24.3     /
      DATA (AE(I, 5, 7),I=1,10)  / 
     +  .000E+00,35.2    ,18.9    ,20.3    ,20.6    ,21.5    ,22.6    ,
     +  23.7    ,24.2    ,24.7     /
      DATA (AE(I, 5, 8),I=1,10)  / 
     +  .000E+00,.000E+00,16.4    ,18.9    ,18.8    ,19.6    ,20.7    ,
     +  22.3    ,23.1    ,23.9     /
      DATA (AE(I, 5, 9),I=1,10)  / 
     +  .000E+00,.000E+00,33.9    ,19.8    ,20.3    ,20.7    ,21.9    ,
     +  23.4    ,24.1    ,24.8     /
      DATA (AE(I, 5,10),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,18.0    ,20.0    ,21.4    ,22.0    ,
     +  23.8    ,24.6    ,25.4     /
      DATA (AE(I, 5,11),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,26.4    ,20.4    ,21.2    ,22.3    ,
     +  23.8    ,24.7    ,25.5     /
      DATA (AE(I, 5,12),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,41.7    ,18.2    ,19.8    ,21.1    ,
     +  22.6    ,23.4    ,24.6     /
      DATA (AE(I, 5,13),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,22.5    ,20.0    ,21.7    ,
     +  22.8    ,23.7    ,24.7     /
      DATA (AE(I, 5,14),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,54.1    ,19.9    ,21.9    ,
     +  23.2    ,24.3    ,25.3     /
      DATA (AE(I, 5,15),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,21.2    ,22.2    ,
     +  23.6    ,24.9    ,25.5     /
      DATA (AE(I, 5,16),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,44.9    ,21.9    ,
     +  23.8    ,25.2    ,25.6     /
      DATA (AE(I, 5,17),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,47.8    ,22.7    ,
     +  23.8    ,24.9    ,26.3     /
      DATA (AE(I, 5,18),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,35.5    ,
     +  23.9    ,25.9    ,26.6     /
      DATA (AE(I, 5,19),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,64.3    ,
     +  24.1    ,25.7    ,27.1     /
      DATA (AE(I, 5,20),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
     +  34.0    ,25.7    ,27.7     /
      DATA (AE(I, 6, 1),I=1,10)  / 
     +  16.6    ,16.5    ,16.8    ,16.7    ,17.0    ,16.5    ,16.7    ,
     +  18.3    ,18.9    ,19.0     /
      DATA (AE(I, 6, 2),I=1,10)  / 
     +  16.2    ,16.6    ,17.2    ,17.4    ,17.9    ,17.4    ,17.7    ,
     +  20.7    ,22.0    ,22.6     /
      DATA (AE(I, 6, 3),I=1,10)  / 
     +  18.9    ,18.7    ,18.8    ,18.6    ,18.9    ,18.6    ,18.9    ,
     +  21.0    ,22.3    ,22.9     /
      DATA (AE(I, 6, 4),I=1,10)  / 
     +  18.3    ,12.7    ,14.2    ,15.0    ,15.7    ,16.1    ,16.3    ,
     +  16.5    ,17.9    ,19.0     /
      DATA (AE(I, 6, 5),I=1,10)  / 
     +  .000E+00,15.7    ,15.1    ,15.3    ,16.5    ,16.4    ,16.4    ,
     +  17.0    ,18.3    ,19.4     /
      DATA (AE(I, 6, 6),I=1,10)  / 
     +  .000E+00,22.9    ,14.9    ,15.2    ,16.2    ,16.9    ,17.4    ,
     +  18.2    ,19.5    ,21.1     /
      DATA (AE(I, 6, 7),I=1,10)  / 
     +  .000E+00,40.7    ,18.4    ,15.9    ,17.1    ,17.7    ,18.9    ,
     +  19.5    ,20.3    ,21.1     /
      DATA (AE(I, 6, 8),I=1,10)  / 
     +  .000E+00,.000E+00,23.3    ,16.2    ,16.3    ,17.3    ,18.7    ,
     +  19.5    ,20.3    ,21.1     /
      DATA (AE(I, 6, 9),I=1,10)  / 
     +  .000E+00,.000E+00,49.2    ,19.0    ,19.1    ,19.4    ,20.2    ,
     +  20.8    ,21.6    ,22.0     /
      DATA (AE(I, 6,10),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,27.2    ,21.2    ,20.8    ,21.4    ,
     +  22.3    ,22.8    ,23.3     /
      DATA (AE(I, 6,11),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,45.6    ,25.0    ,22.8    ,23.9    ,
     +  23.6    ,24.3    ,24.4     /
      DATA (AE(I, 6,12),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,45.8    ,29.7    ,25.1    ,25.3    ,
     +  25.3    ,26.0    ,26.3     /
      DATA (AE(I, 6,13),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,42.7    ,29.0    ,28.0    ,
     +  27.0    ,27.2    ,27.6     /
      DATA (AE(I, 6,14),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,62.0    ,32.0    ,30.0    ,
     +  29.8    ,29.5    ,29.6     /
      DATA (AE(I, 6,15),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,44.5    ,34.4    ,
     +  32.7    ,31.5    ,31.8     /
      DATA (AE(I, 6,16),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,75.6    ,37.1    ,
     +  34.6    ,34.4    ,34.4     /
      DATA (AE(I, 6,17),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,51.2    ,45.2    ,
     +  39.0    ,37.5    ,36.4     /
      DATA (AE(I, 6,18),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,74.9    ,
     +  42.3    ,39.9    ,38.3     /
      DATA (AE(I, 6,19),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,69.5    ,
     +  50.7    ,42.3    ,41.4     /
      DATA (AE(I, 6,20),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
     +  66.3    ,48.0    ,43.4     /
      DATA (AE(I, 7, 1),I=1,10)  / 
     +  27.0    ,25.8    ,26.3    ,26.2    ,26.7    ,26.7    ,27.1    ,
     +  27.1    ,27.2    ,19.0     /
      DATA (AE(I, 7, 2),I=1,10)  / 
     +  29.1    ,28.9    ,29.7    ,30.3    ,31.0    ,31.4    ,32.0    ,
     +  32.3    ,32.7    ,22.6     /
      DATA (AE(I, 7, 3),I=1,10)  / 
     +  31.6    ,29.7    ,30.9    ,31.4    ,32.5    ,33.1    ,34.0    ,
     +  34.6    ,35.1    ,22.9     /
      DATA (AE(I, 7, 4),I=1,10)  / 
     +  27.4    ,19.9    ,20.8    ,22.8    ,24.6    ,26.4    ,28.2    ,
     +  29.6    ,30.8    ,19.0     /
      DATA (AE(I, 7, 5),I=1,10)  / 
     +  .000E+00,24.6    ,24.1    ,25.0    ,27.2    ,28.7    ,30.7    ,
     +  31.8    ,32.9    ,19.4     /
      DATA (AE(I, 7, 6),I=1,10)  / 
     +  .000E+00,35.6    ,25.2    ,25.6    ,27.9    ,30.4    ,32.7    ,
     +  34.6    ,36.3    ,21.1     /
      DATA (AE(I, 7, 7),I=1,10)  / 
     +  .000E+00,45.4    ,30.9    ,28.2    ,29.0    ,31.2    ,34.0    ,
     +  35.8    ,37.4    ,21.1     /
      DATA (AE(I, 7, 8),I=1,10)  / 
     +  .000E+00,.000E+00,38.2    ,29.6    ,29.4    ,30.3    ,33.2    ,
     +  35.5    ,37.6    ,21.1     /
      DATA (AE(I, 7, 9),I=1,10)  / 
     +  .000E+00,.000E+00,59.3    ,34.5    ,33.7    ,32.9    ,35.4    ,
     +  37.6    ,39.6    ,22.0     /
      DATA (AE(I, 7,10),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,44.5    ,37.8    ,37.5    ,37.2    ,
     +  39.0    ,41.4    ,23.3     /
      DATA (AE(I, 7,11),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,67.0    ,43.6    ,42.0    ,40.8    ,
     +  41.4    ,43.0    ,24.4     /
      DATA (AE(I, 7,12),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,49.9    ,50.9    ,44.6    ,43.9    ,
     +  44.2    ,44.2    ,26.3     /
      DATA (AE(I, 7,13),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,67.2    ,50.5    ,48.7    ,
     +  48.1    ,47.2    ,27.6     /
      DATA (AE(I, 7,14),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,68.1    ,55.2    ,52.3    ,
     +  51.5    ,51.6    ,29.6     /
      DATA (AE(I, 7,15),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,68.7    ,58.6    ,
     +  56.5    ,55.7    ,31.8     /
      DATA (AE(I, 7,16),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,89.3    ,62.9    ,
     +  60.0    ,59.1    ,34.4     /
      DATA (AE(I, 7,17),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,56.0    ,72.9    ,
     +  66.3    ,64.2    ,36.4     /
      DATA (AE(I, 7,18),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,105.    ,
     +  71.3    ,68.3    ,38.3     /
      DATA (AE(I, 7,19),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,73.4    ,
     +  76.8    ,72.4    ,41.4     /
      DATA (AE(I, 7,20),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
     +  107.    ,79.9    ,43.4     /
      DATA (AE(I, 8, 1),I=1,10)  / 
     +  35.5    ,35.3    ,35.7    ,35.7    ,36.3    ,36.3    ,36.7    ,
     +  36.7    ,36.7    ,19.0     /
      DATA (AE(I, 8, 2),I=1,10)  / 
     +  40.6    ,41.4    ,41.9    ,42.3    ,43.2    ,43.5    ,44.0    ,
     +  44.3    ,44.5    ,22.6     /
      DATA (AE(I, 8, 3),I=1,10)  / 
     +  45.4    ,45.7    ,46.4    ,47.0    ,48.1    ,48.7    ,49.4    ,
     +  49.8    ,50.2    ,22.9     /
      DATA (AE(I, 8, 4),I=1,10)  / 
     +  43.9    ,44.3    ,43.4    ,45.1    ,47.3    ,48.7    ,49.6    ,
     +  50.5    ,51.3    ,19.0     /
      DATA (AE(I, 8, 5),I=1,10)  / 
     +  .000E+00,49.3    ,49.6    ,50.5    ,53.2    ,54.2    ,55.4    ,
     +  56.1    ,56.8    ,19.4     /
      DATA (AE(I, 8, 6),I=1,10)  / 
     +  .000E+00,59.1    ,53.0    ,55.4    ,58.0    ,60.0    ,61.2    ,
     +  62.5    ,63.6    ,21.1     /
      DATA (AE(I, 8, 7),I=1,10)  / 
     +  .000E+00,54.5    ,57.1    ,59.2    ,62.3    ,64.4    ,66.0    ,
     +  67.3    ,68.5    ,21.1     /
      DATA (AE(I, 8, 8),I=1,10)  / 
     +  .000E+00,.000E+00,65.9    ,62.1    ,65.1    ,67.6    ,69.4    ,
     +  71.1    ,72.6    ,21.1     /
      DATA (AE(I, 8, 9),I=1,10)  / 
     +  .000E+00,.000E+00,72.2    ,67.1    ,70.5    ,73.1    ,75.1    ,
     +  76.8    ,78.4    ,22.0     /
      DATA (AE(I, 8,10),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,80.1    ,75.0    ,78.0    ,80.0    ,
     +  82.1    ,83.9    ,23.3     /
      DATA (AE(I, 8,11),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,94.5    ,82.2    ,82.8    ,85.1    ,
     +  87.3    ,89.2    ,24.4     /
      DATA (AE(I, 8,12),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,56.8    ,92.5    ,87.2    ,89.4    ,
     +  91.9    ,94.1    ,26.3     /
      DATA (AE(I, 8,13),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,116.    ,96.2    ,94.4    ,
     +  97.0    ,99.2    ,27.6     /
      DATA (AE(I, 8,14),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,78.1    ,104.    ,102.    ,
     +  102.    ,105.    ,29.6     /
      DATA (AE(I, 8,15),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,128.    ,111.    ,
     +  109.    ,110.    ,31.8     /
      DATA (AE(I, 8,16),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,104.    ,118.    ,
     +  117.    ,115.    ,34.4     /
      DATA (AE(I, 8,17),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,64.4    ,138.    ,
     +  124.    ,122.    ,36.4     /
      DATA (AE(I, 8,18),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,133.    ,
     +  133.    ,132.    ,38.3     /
      DATA (AE(I, 8,19),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,83.6    ,
     +  146.    ,139.    ,41.4     /
      DATA (AE(I, 8,20),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
     +  166.    ,147.    ,43.4     /
      DATA (AE(I, 9, 1),I=1,10)  / 
     +  43.3    ,43.2    ,43.6    ,43.8    ,44.1    ,44.3    ,44.7    ,
     +  44.8    ,44.8    ,19.0     /
      DATA (AE(I, 9, 2),I=1,10)  / 
     +  50.9    ,51.4    ,52.0    ,52.6    ,53.1    ,53.6    ,54.2    ,
     +  54.5    ,54.7    ,22.6     /
      DATA (AE(I, 9, 3),I=1,10)  / 
     +  58.0    ,58.4    ,59.3    ,60.1    ,60.7    ,61.5    ,62.3    ,
     +  62.7    ,63.1    ,22.9     /
      DATA (AE(I, 9, 4),I=1,10)  / 
     +  62.0    ,63.9    ,63.7    ,65.7    ,65.5    ,67.5    ,68.2    ,
     +  68.9    ,69.7    ,19.0     /
      DATA (AE(I, 9, 5),I=1,10)  / 
     +  .000E+00,72.2    ,72.5    ,74.2    ,74.2    ,76.1    ,77.0    ,
     +  77.8    ,78.6    ,19.4     /
      DATA (AE(I, 9, 6),I=1,10)  / 
     +  .000E+00,80.4    ,80.5    ,83.1    ,83.0    ,85.5    ,86.8    ,
     +  88.1    ,89.2    ,21.1     /
      DATA (AE(I, 9, 7),I=1,10)  / 
     +  .000E+00,63.4    ,88.5    ,91.3    ,91.1    ,94.0    ,95.8    ,
     +  97.3    ,98.6    ,21.1     /
      DATA (AE(I, 9, 8),I=1,10)  / 
     +  .000E+00,.000E+00,98.8    ,98.6    ,97.8    ,102.    ,104.    ,
     +  106.    ,108.    ,21.1     /
      DATA (AE(I, 9, 9),I=1,10)  / 
     +  .000E+00,.000E+00,84.1    ,107.    ,107.    ,111.    ,113.    ,
     +  116.    ,117.    ,22.0     /
      DATA (AE(I, 9,10),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,116.    ,115.    ,119.    ,122.    ,
     +  125.    ,127.    ,23.3     /
      DATA (AE(I, 9,11),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,111.    ,123.    ,127.    ,131.    ,
     +  134.    ,137.    ,24.4     /
      DATA (AE(I, 9,12),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,65.6    ,136.    ,135.    ,140.    ,
     +  143.    ,146.    ,26.3     /
      DATA (AE(I, 9,13),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,146.    ,144.    ,149.    ,
     +  152.    ,155.    ,27.6     /
      DATA (AE(I, 9,14),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,88.7    ,152.    ,158.    ,
     +  162.    ,165.    ,29.6     /
      DATA (AE(I, 9,15),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,181.    ,167.    ,
     +  171.    ,174.    ,31.8     /
      DATA (AE(I, 9,16),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,117.    ,174.    ,
     +  180.    ,183.    ,34.4     /
      DATA (AE(I, 9,17),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,72.0    ,201.    ,
     +  189.    ,192.    ,36.4     /
      DATA (AE(I, 9,18),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,151.    ,
     +  198.    ,201.    ,38.3     /
      DATA (AE(I, 9,19),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,95.2    ,
     +  220.    ,210.    ,41.4     /
      DATA (AE(I, 9,20),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
     +  192.    ,217.    ,43.4     /
      DATA (AE(I,10, 1),I=1,10)  / 
     +  62.1    ,62.1    ,62.6    ,62.9    ,63.3    ,63.3    ,64.0    ,
     +  64.0    ,64.0    ,19.0     /
      DATA (AE(I,10, 2),I=1,10)  / 
     +  75.1    ,75.4    ,76.3    ,76.8    ,77.6    ,77.9    ,78.8    ,
     +  79.0    ,79.3    ,22.6     /
      DATA (AE(I,10, 3),I=1,10)  / 
     +  87.5    ,88.3    ,89.4    ,90.2    ,91.3    ,91.9    ,93.0    ,
     +  93.5    ,93.9    ,22.9     /
      DATA (AE(I,10, 4),I=1,10)  / 
     +  104.    ,104.    ,105.    ,106.    ,107.    ,108.    ,109.    ,
     +  110.    ,110.    ,19.0     /
      DATA (AE(I,10, 5),I=1,10)  / 
     +  .000E+00,122.    ,122.    ,123.    ,124.    ,125.    ,126.    ,
     +  127.    ,128.    ,19.4     /
      DATA (AE(I,10, 6),I=1,10)  / 
     +  .000E+00,138.    ,139.    ,140.    ,142.    ,143.    ,144.    ,
     +  146.    ,147.    ,21.1     /
      DATA (AE(I,10, 7),I=1,10)  / 
     +  .000E+00,85.3    ,158.    ,159.    ,161.    ,162.    ,164.    ,
     +  166.    ,167.    ,21.1     /
      DATA (AE(I,10, 8),I=1,10)  / 
     +  .000E+00,.000E+00,176.    ,177.    ,179.    ,181.    ,183.    ,
     +  184.    ,186.    ,21.1     /
      DATA (AE(I,10, 9),I=1,10)  / 
     +  .000E+00,.000E+00,114.    ,199.    ,201.    ,202.    ,205.    ,
     +  206.    ,207.    ,22.0     /
      DATA (AE(I,10,10),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,218.    ,219.    ,220.    ,224.    ,
     +  225.    ,226.    ,23.3     /
      DATA (AE(I,10,11),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,150.    ,238.    ,238.    ,243.    ,
     +  244.    ,245.    ,24.4     /
      DATA (AE(I,10,12),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,85.8    ,255.    ,255.    ,261.    ,
     +  262.    ,263.    ,26.3     /
      DATA (AE(I,10,13),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,195.    ,272.    ,279.    ,
     +  279.    ,280.    ,27.6     /
      DATA (AE(I,10,14),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,115.    ,290.    ,296.    ,
     +  297.    ,298.    ,29.6     /
      DATA (AE(I,10,15),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,263.    ,313.    ,
     +  314.    ,315.    ,31.8     /
      DATA (AE(I,10,16),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,150.    ,330.    ,
     +  331.    ,332.    ,34.4     /
      DATA (AE(I,10,17),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,90.0    ,319.    ,
     +  349.    ,349.    ,36.4     /
      DATA (AE(I,10,18),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,196.    ,
     +  366.    ,367.    ,38.3     /
      DATA (AE(I,10,19),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,122.    ,
     +  387.    ,384.    ,41.4     /
      DATA (AE(I,10,20),I=1,10)  / 
     +  .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
     +  247.    ,401.    ,43.4     /
      DATA (ERES(I, 1),I=1,10)  / 10*0./
      DATA (ERES(I, 2),I=1,10)  / 10*0./
      DATA (ERES(I, 3),I=1,10)  / 10*0./
      DATA (ERES(I, 4),I=1,10)  / 10*0./
      DATA (ERES(I, 5),I=1,10)  / 10*0./
      DATA (ERES(I, 6),I=1,10)  / 
     +     0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,
     +     2.780,   2.880,   2.890 /
      DATA (ERES(I, 7),I=1,10)  / 
     +     1.500,   2.460,   2.510,   2.610,   2.700,   2.920,   3.070,
     +     3.200,   3.330,   2.890 /
      DATA (ERES(I, 8),I=1,10)  / 
     +     4.470,   4.350,   4.390,   4.550,   4.660,   4.890,   4.980,
     +     5.100,   5.220,   2.890 /
      DATA (ERES(I, 9),I=1,10)  / 
     +     7.480,   7.380,   7.370,   7.480,   7.510,   7.630,   7.660,
     +     7.750,   7.820,   2.890 /
      DATA (ERES(I,10),I=1,10)  / 
     +    15.270,  15.190,  15.200,  15.370,  15.380,  15.430,  15.540,
     +    15.590,  15.630,   2.890 /
      END
C->
      SUBROUTINE FRAGM (IAT,IAP, NW,B, NF, IAF)
C...Nuclear Fragmentation, Abrasion-ablation model, 
C...Based on Jon Engel's routines ABRABL 
C...This most recent version adds for all prefragment
C...masses > 10 the model calculation for the fragment
C...mass distribution and the energy carried by the fragment
C...of W. Friedmann
C...The average values are used to implement the model
C...in the montecarlo fashion / TSS, Dec '91
C.
C.  INPUT: IAP = mass of incident nucleus
C.         IAT = mass of target   nucleus
C.         NW = number of wounded nucleons in the beam nucleus
C.         B  = impact parameter in the interaction
C.     
C.  OUTPUT : NF = number of fragments  of the spectator nucleus
C.           IAF(1:NF) = mass number of each fragment
C.           PF(3,60) in common block /FRAGMENTS/ contains
C.           the three momentum components (MeV/c) of each
C.           fragment in the projectile frame
C..............................................................
      SAVE

      COMMON /FRAGMENTS/ PPP(3,60)
      COMMON /FRAGMOD/A(10,10,20),AE(10,10,20),ERES(10,10),NFLAGG(10,10)
      DIMENSION IAF(60)
      DIMENSION AA(10), EAA(10) 
      DATA AA/10.,15.,20.,25.,30.,35.,40.,45.,50.,56./
      DATA EAA/1.,2.,4.,6.,8.,10.,12.,16.,20.,30/
      AP=IAP
      AT=IAT
      NPF = IAP - NW
      IF (NPF .EQ. 0) THEN
         NF = 0
         RETURN
      ENDIF

      EB = ESTAR(AP,AT, B)
      EBP = ESTARP (NPF, NW)
C CONTRIBUTION TO E* FROM ENERGY DEPOSITED BY SECONDARIES
      EB = EB + EBP
C TOTAL E* IS THE SUM OF THE TWO COMPONENTS

C.....Prefragment transverse momentum (MeV/nucleon)...
            FK = FERMK(AP)
C FERMI MOMENTUM OF THE PROJECTILE NUCLEUS
            IF (NW .LT. IAP) THEN
            SIG = FK*SQRT(NW*NPF/(AP-1.))/3.162
C GAUSSIAN SIGMA IN ALL THREE DIRECTION
            ELSE
            SIG = FK/3.162
C THIS IS NOT CORRECT, TOO LARGE !!!!!!!!!!!!!!
            ENDIF
             PPFX = SIG*GASDEV(0)/NPF
             PPFY = SIG*GASDEV(0)/NPF
C THREE MOMENTUM COMPONENTS PER NUCLEON FOR THE PREFRAGMENT

C.............Crude model for small prefragment mass .......
            IF (NPF .LT. 10) THEN
                 CALL EVAP(NPF, EB, EPS, NNUC, NALP)
C                  EPS IS THE KINETIC ENERGY CARRIED BY THE EVAPORATED NUCLEONS
               ETOT = 938. + EPS
                 PP = SQRT((ETOT*ETOT - 8.79844E5)/3.)
C                  AVERAGE MOMENTUM OF EVAPORATED NUCLEONS IN EACH DIRECTION
                 NUC = NPF - NNUC - 4*NALP
                 NF = 0
                 IF (NUC .GT. 0) THEN
                    NF = NF + 1
                    IAF(NF) = NUC
                    PPP(1,NF) = NUC*PPFX
                    PPP(2,NF) = NUC*PPFY
                 ENDIF
                 IF (NALP .NE. 0) THEN
                 DO I=1,NALP
                   NF = NF + 1
                    IAF(NF) = 4
                   CALL SINCO(S1,C1)
                   CALL SINCO(S2,C2)
                   PXE = 4.*PP*S1*S2
                   PYE = 4.*PP*S1*C2
                   PPP(1,NF) = 4.*PPFX + PXE
                   PPP(2,NF) = 4.*PPFY + PYE
                   PPP(1,1) = PPP(1,1) - PXE
                   PPP(2,1) = PPP(2,1) - PYE
                 ENDDO
                 ENDIF
                 IF (NNUC .NE. 0) THEN
                 DO I=1,NNUC
                    NF = NF + 1
                    IAF(NF) = 1
                    CALL SINCO(S1,C1)
                    CALL SINCO(S2,C2)
                    PXE = PP*S1*S2
                    PYE = PP*S1*C2
                    PPP(1,NF) = 4.*PPFX + PXE
                    PPP(2,NF) = 4.*PPFY + PYE
                    PPP(1,1) = PPP(1,1) - PXE
                    PPP(2,1) = PPP(2,1) - PYE
                 ENDDO
                 ENDIF
                 RETURN
            ENDIF

C.........More refined model calculation .............
      JA = NPF/5 -1
      IF (JA .LT. 10) THEN
      IF ((NPF - AA(JA)) .GT. (AA(JA+1)-NPF)) JA = JA + 1
      ENDIF
      ARAT = FLOAT(NPF)/AA(JA)
      DO J=1,10
      IF (EB .LT. EAA(J)) GO TO 29
      ENDDO
      JE = 10
      GO TO 39
   29      JE = J
   39      IF (JE .GT. 1 .AND. JE .NE. 10) THEN
      IF ((EB - EAA(J-1)) .LT. (EAA(J)-EB)) JE = J - 1
      ENDIF
      ERAT = EB/EAA(JE)
        IF (EB .LT. 1.) THEN
        ERAT = EB
        ENDIF
C INTERPOLATE BETWEEN EB=0. (NOTHING HAPPENS) AND EB = 1. MeV

         IF (JA .EQ. 10 .AND. JE .GT. 6) THEN
            WRITE(*,*)' JA=',JA,',   JE=',JE
         ENDIF
   43    ESUM = 0.
      NSUM = 0
      JF = 0
      DO J=20,1,-1
      FR =  A(JA, JE, J)*ARAT*ERAT
      N1 = 1 + FR
      FR1 = FR/FLOAT(N1)
      DO K=1, N1
      IF (S_RNDM(0) .LT. FR1) THEN
      JF = JF + 1
      IAF(JF) = J
      NSUM = NSUM + J
      EKIN = ERAT*AE(JA,JE, J)
         IF (EKIN .GT. 0.) THEN
         ESUM = ESUM + EKIN
         ETOT = 938.*IAF(JF) + EKIN
         PP = SQRT(2.*(ETOT*ETOT - IAF(JF)**2*8.79844E5)/3.)
         CALL SINCO(S1,C1)
         CALL SINCO(S2,C2)
         PPP(1,JF) = PP*S1*S2 + IAF(JF)*PPFX
         PPP(2,JF) = PP*S1*C2 + IAF(JF)*PPFY
         ENDIF
         IF (NSUM .GT. NPF) THEN
C         WRITE(*,*)' WARNING, NSUM=', NSUM,',  NPF=',NPF
C         WRITE(*,*)'  ARAT =', ARAT
         GO TO 43
        ELSE
        IF (NSUM .EQ. NPF) THEN
        GO TO 44
        ENDIF
        ENDIF
      ENDIF
      ENDDO
      ENDDO
      IF (NFLAGG(JA,JE) .EQ. 0) THEN
C 'THE RESIDUE' IS A NUCLEAR FRAGMENT
      JF = JF + 1
      IAF(JF) = NPF - NSUM
      F1 = NPF*EB - ESUM
      IF (F1 .LT. 0.) F1 = 0.
C GIVE THE REST OF EB TO THE FRAGMENT
      EKIN = F1
         IF (EKIN .GT. 0.) THEN
         ETOT = 938.*IAF(JF) + EKIN
         PP = SQRT(2.*(ETOT*ETOT - IAF(JF)**2*8.79844E5)/3.)
         CALL SINCO(S1,C1)
         CALL SINCO(S2,C2)
         PPP(1,JF) = PP*S1*S2 + IAF(JF)*PPFX
         PPP(2,JF) = PP*S1*C2 + IAF(JF)*PPFY
         ENDIF
      ELSE
C 'THE RESIDUE' CONSISTS OF SPECTATOR NUCLEONS
      N1 = NPF - NSUM
      DO K=1,N1
      JF = JF + 1
      IAF(JF) = 1
      EKIN = ERAT*ERES(JA,JE)
         IF (EKIN .GT. 0.) THEN
         ETOT = 938.*IAF(JF) + EKIN
         PP = SQRT(2.*(ETOT*ETOT - IAF(JF)**2*8.79844E5)/3.)
         CALL SINCO(S1,C1)
         CALL SINCO(S2,C2)
         PPP(1,JF) = PP*S1*S2 + PPFX
         PPP(2,JF) = PP*S1*C2 + PPFY
         ENDIF
      ENDDO
      ENDIF
  44  NF = JF
      RETURN
      END
C->
      FUNCTION ESTARP (NPF, NW)
C CONTRIBUTION TO E* FROM ENERGY DEPOSITED BY SECONDARIES
C VERY NAIVE VERSION INCORPORATING HUEFFNER'S IDEAS
      SAVE

      APF = NPF
      F1 = 15.3/APF**0.666666666
C AVERAGE KINETIC ENERGY/NUCLEON IN PREFRAGMENT (MeV)
C PER PATHLENGTH EQUAL TO THE PREFRAGMENT RADIUS
      ESTARP = 0.
      DO I=1,NW
      IF (S_RNDM(0) .GT. 0.5) THEN
      F2 = F1*RDIS(0)
      ESTARP = ESTARP + F2
      ENDIF
      ENDDO
C SAMPLE RANDOMLY PER WOUNDED NUCLEON, x NW
      RETURN
      END
      
      function rdis(Idum)
      SAVE

      dimension probr(20)
      data probr/
     *      0.10000, 0.15748, 0.21778, 0.28605, 0.36060,
     *      0.43815, 0.51892, 0.60631, 0.70002, 0.79325,
     *      0.88863, 0.98686, 1.10129, 1.21202, 1.32932,
     *      1.44890, 1.57048, 1.70139, 1.83417, 2.00000/
      nr = 20.*S_RNDM(0) + 1
      if (nr .eq. 1) then
      f1 = 0.
      else
      f1 = probr(nr-1)
      endif
      dr = probr(nr) - f1
      rdis = f1 + dr*S_RNDM(0)
      return
      end


      function estar(ap,at,b)
      implicit real*8(a-h,o-z)
      SAVE

      real*4 ap,at,b,estar
      sigma=4.5  !total n-n cross section in fm**2
      rt=.82*at**.3333 !target radius
      rp=.82*ap**.3333 !projectile radius
      alpha=rt**2/rp**2
      beta=b**2/rt**2
      f=at*sigma/(3.14159*rt**2)
      alf = log(f)
      alalf = log(alpha)
      gfac=0
      gfac1=0
      s1=0.
      s2=0.
      s3=0.      
      ii=1
      do n=0,10 ! This limit may not need to be so high.
         if(n.ge.2) then
            gfac1=gfac
            gfac=gfac+log(float(n)) 
         endif
         g0=n*alf -n*beta*alpha/(n+alpha)+alalf
         g1=g0-log(alpha+n)-gfac
         g2=(n+2)*log(f)-(n+2)*beta*alpha/(n+2+alpha) 
     >      +log(n+2+alpha+beta*alpha**2)-3*log(n+2+alpha)-gfac
         g3=g0-2*log(n+alpha)-gfac1
         ii=-ii
         s1=s1+ii*exp(g1)
         s2=s2+ii*exp(g2)
         if(n.ge.1) s3=s3+ii*exp(g3)
      enddo

      pb=s1
      e1b=197.**2/(2*938.*rp**2*pb) *s2
c      a=b*(s3/pb-1)
c      a=-b*s3/pb
c      e2b=-.5* 938. * (41./(ap**.333))**2 * a**2 /(197.**2)
c      estar=e1b+e2b
      estar = e1b
      return
      end

      subroutine evap(npf,eb,eps,nnuc,nalp)
      SAVE

      eps=7.5+sqrt(8*eb)
      n=min(npf*int(eb/eps),npf)
      nalp=n/5
      nnuc=n-4*nalp
      return
      end
C->
      FUNCTION FERMK(A)
      SAVE

      DIMENSION AA(6), FK(6)
      DATA AA/4., 6., 12., 24., 40., 57./
      DATA FK/130.,169.,221.,235.,251.,260./
      DO I=2,4
      IF (A .LT. AA(I)) GO TO 25
      ENDDO
      I = 5
   25      F11 = AA(I-1)
      F12 = AA(I)
      F13 = AA(I+1)
      F21 = FK(I-1)
      F22 = FK(I)
      F23 = FK(I+1)
      FERMK = QUAD_INT(A,F11,F12,F13, F21,F22,F23)
      RETURN
      END

C========================================================================
C. Multiple interaction structure
C========================================================================

      SUBROUTINE INT_NUC (IA, IB, SIG0, SIGEL) 
C...Compute with a montecarlo code  the  "multiple interaction structure"
C.  of a nucleus-nucleus interaction
C.
C.  INPUT : IA            = mass of target nucleus
C.          IB            = mass of projectile nucleus
C.          SIG0 (mbarn)  = inelastic pp cross section
C.          SIGEL(mbarn)  = elastic pp cross section
C.
C.  OUTPUT : in common block /CNUCMS/
C.           B = impact parameter (fm)
C.           BMAX = maximum impact parameter for generation
C.           NTRY = number of "trials" before one interaction
C.           NA = number of wounded nucleons in A
C.           NB =    "        "        "     in B
C.           NI = number of nucleon-nucleon inelastic interactions 
C.           NAEL = number of elastically scattered nucleons in  A 
C.           NBEL =    "         "           "          "    in  B
C.           JJA(J)  [J=1:IA]   = number of inelastic interactions 
C.                                of J-th nucleon of nucleus A
C.           JJB(J)  [J=1:IB]   = number of inelastic interactions 
C.                                of J-th nucleon of nucleus B
C.           JJAEL(J)  [J=1:IA]   = number of elastic interactions 
C.                                of J-th nucleon of nucleus A
C.           JJBEL(J)  [J=1:IB]   = number of elastic interactions 
C.                                of J-th nucleon of nucleus B
C.           JJINT(J,K)  [J=1:NB, K=1:NA]  (0 = no interaction) 
C.                                         (1 = interaction )
C.                                         between nucleon J of A and K of B
C-----------------------------------------------------------------------------
      SAVE

      PARAMETER (IAMAX=56)
      COMMON /CNUCMS/ B, BMAX, NTRY, NA, NB, NI, NAEL, NBEL
     +         ,JJA(IAMAX), JJB(IAMAX), JJINT(IAMAX,IAMAX)
     +         ,JJAEL(IAMAX), JJBEL(IAMAX)
      DIMENSION XA(IAMAX), YA(IAMAX), XB(IAMAX), YB(IAMAX)
      DATA PI /3.1415926/
      SIGT = SIG0 + SIGEL
      R2  = 0.1 * SIG0/PI
      R2T = 0.1 * SIGT/PI
      BMAX = 15.                             ! fm
      NTRY = 0
      CALL NUC_CONF (IA, XA, YA)
      CALL NUC_CONF (IB, XB, YB)
      NI = 0
      NIEL = 0
      DO JA=1,IA
         JJA(JA) = 0
         JJAEL(JA) = 0
      ENDDO
      DO JB=1,IB
         JJB(JB) = 0
         JJBEL(JB) = 0
         DO JA=1,IA
            JJINT(JB,JA) = 0
         ENDDO
      ENDDO
1000  B = BMAX*SQRT(S_RNDM(0))
      PHI = 2.*PI*S_RNDM(0)
      BX = B*COS(PHI)
      BY = B*SIN(PHI)
      NTRY = NTRY+1
      DO JA=1,IA
         DO JB=1,IB
            S = (XA(JA)-XB(JB)-BX)**2 + (YA(JA)-YB(JB)-BY)**2
            IF (S .LT. R2)  THEN
               NI = NI + 1
               JJA(JA) = JJA(JA)+1
               JJB(JB) = JJB(JB)+1
               JJINT(JB,JA) = 1
            ELSE IF (S .LT. R2T)  THEN
               NIEL = NIEL + 1
               JJAEL(JA) = JJAEL(JA)+1
               JJBEL(JB) = JJBEL(JB)+1
            ENDIF
         ENDDO
      ENDDO
      IF (NI + NIEL .EQ. 0)  GOTO 1000
      NA = 0
      NB = 0
      NAEL = 0
      NBEL = 0
      DO JA=1,IA
         IF (JJA(JA) .GT. 0)  THEN
            NA = NA + 1
         ELSE
            IF (JJAEL(JA) .GT. 0)  NAEL = NAEL+1
         ENDIF
      ENDDO
      DO JB=1,IB
         IF (JJB(JB) .GT. 0)  THEN
            NB = NB + 1
         ELSE
            IF (JJBEL(JB) .GT. 0)  NBEL = NBEL+1
         ENDIF
      ENDDO
      RETURN
      END

       SUBROUTINE NUC_CONF (IA, XX, YY)
C...This routine generates the configuration  of a nucleus 
C.  need an initialization call to NUC_GEOM_INI
C.
C.  INPUT  : IA = mass number of the nucleus
C.  OUTPUT : XX(1:IA), YY(1:IA) (fm) = position in impact parameter
C.                                     space of the IA nucleons
C...................................................................
      SAVE

      PARAMETER (IAMAX=56)
      DIMENSION XX(IAMAX), YY(IAMAX)
      PARAMETER (NB=401)
      COMMON /CPROFA/ ZMIN, DZ, BBZ(NB,IAMAX)
      DATA PI /3.1415926/
      DO J=1,IA
         Z = S_RNDM(0)
         JZ = INT((Z-ZMIN)/DZ)+1
CDH
         JZ = MIN(JZ,400)
         T = (Z-ZMIN)/DZ - FLOAT(JZ-1)
         B = BBZ(JZ,IA)*(1.-T) + BBZ(JZ+1,IA)*T
         PHI = 2.*PI*S_RNDM(0)
         XX(J) = B*COS(PHI)
         YY(J) = B*SIN(PHI)
      ENDDO
      RETURN
      END

      SUBROUTINE NUC_GEOM_INI
C...Initialize all nucleus profiles
      SAVE

      PARAMETER (NB=401)
      PARAMETER (IAMAX=56)
      COMMON /CPROF/ DB, BMAX, BB(NB), TB(NB), A
      COMMON /CPROFA/ ZMIN, DZ, BBZ(NB,IAMAX)
      DIMENSION FFB(NB), GGB(NB)
      DATA PI /3.1415926/
      CALL SHELL_INI
      CALL WOOD_SAXON_INI
      DO IA= 2,IAMAX
           JA = IA
         CALL NUC_PROFIL(JA)
         DO K=1,NB
           FFB(K) = BB(K)*TB(K) * (2.*PI)
         ENDDO            
         GGB(1) = 0.
         GGB(NB) = 1.
         DO K=2,NB-1
           GGB(K) = GGB(K-1) + FFB(K-1)*DB
         ENDDO            
         CALL INVERT_ARRAY(GGB,0.,DB,NB, BBZ(1,IA), ZMIN, DZ)
      ENDDO
      RETURN
      END

      SUBROUTINE NUC_PROFIL (JA)
C...Compute the profile function T(b)
C.  normalised as INT[d2b T(b) = 1]
C.  INPUT : JA = integer mass number of nucleus
C...............................................
      SAVE

      PARAMETER (NB=401)
      EXTERNAL DENSA
      REAL DENSA
      COMMON /CC01/  B
      COMMON /CCDA/ JJA
      COMMON /CPROF/ DB, BMAX, BB(NB), TB(NB), A
      BMAX = 7.5
      DB = BMAX/FLOAT(NB-1)
      JJA = JA
      A = JA
      DO JB=1,NB
        B = DB*FLOAT(JB-1)
        BB(JB) = B
        IF (JA .LE. 18)  THEN
            TB(JB) = PROFNUC (B, JA)
         ELSE
            TB(JB) = 2.*GAUSS (DENSA,0.,BMAX)
         ENDIF
      ENDDO
      RETURN
      END

      SUBROUTINE NUC1_PROFIL (AA)
C...Compute the profile function T(b)
C.  normalised as INT[d2b T(b) = 1]
C.  INPUT : AA = mass number of nucleus
C...............................................
      SAVE

      PARAMETER (NB=401)
      EXTERNAL DENSA
      REAL DENSA
      COMMON /CC01/  B
      COMMON /CPROF/ DB, BMAX, BB(NB), TB(NB), A
      A = AA
      IA1 = INT(AA)
      IA2 = IA1 + 1
      U = AA - FLOAT(IA1)
      BMAX = 7.5
      DB = BMAX/FLOAT(NB-1)
      DO JB=1,NB
         B = DB*FLOAT(JB-1)
         BB(JB) = B
         IF (A .LE. 18.)  THEN
             T1 = PROFNUC (B, IA1)
             T2 = PROFNUC (B, IA2)
          ELSE
             JJA = IA1
             T1 = 2.*GAUSS (DENSA,0.,BMAX)
             JJA = IA2
             T2 = 2.*GAUSS (DENSA,0.,BMAX)
          ENDIF
          TB(JB) = (1.-U)*T1  + U*T2
      ENDDO
      RETURN
      END

C===========================================================================
C.   Code about nuclear densities
C===========================================================================

      FUNCTION DENS_NUC (R, JA)
C....Nuclear density (normalised to 1)
C.   for a nucleus of mass number JA
C.   INPUT R = radial coordinate  (fm)
C.         JA = integer mass number
C.  OUTPUT (fm**-3)
C--------------------------------------------------------
      SAVE

      COMMON /CWOOD/ RR0(19:56), AA0(19:56), CC0(19:56)
      IF (JA .GT. 18)  THEN
         DENS_NUC = WOOD_SAXON(R,JA)
      ELSE IF (JA .NE. 4)  THEN
         DENS_NUC = HELIUM(R)
      ELSE
         DENS_NUC = SHELL(R,JA)
      ENDIF
      RETURN
      END

      FUNCTION WOOD_SAXON (R, JA) 
C....Wood-Saxon nuclear density (normalised to 1)
C.   for a nucleus of mass number A.
C.   INPUT R =  (fm)
C.         JA = mass number
C.   OUTPUT (fm**-3)
C------------------------------------------------------
      SAVE

      COMMON /CWOOD/ RR0(19:56), AA0(19:56), CC0(19:56)
      WOOD_SAXON = CC0(JA)/(1.+EXP((R-RR0(JA))/AA0(JA)))
      RETURN
      END      

      FUNCTION HELIUM (R)
C... Helium density from Barrett and Jackson
C.   INPUT R = r coordinate (fm)
C.   OUTPUT (fm**-3)
C........................................................
      SAVE

      DATA R0 /0.964/, CA /0.322/   ! fm
      DATA W /0.517/, CC /5.993224E-02/
      HELIUM = CC*(1.+W*(R/R0)**2)/(1. + EXP((R-R0)/CA))
      RETURN
      END

      FUNCTION SHELL (R,JA)
C...Density in the shell model
      COMMON /CSHELL/ RR0(18), RR02(18)
      SAVE

      DATA PI /3.1415926/
      R0 = RR0(JA)
      C1 = MIN(1.,4./FLOAT(JA))
      CS = 1./(R0**3*PI**(1.5))
      CP = 2.*CS/3.
      FS = EXP(-(R/R0)**2)
      FP = (R/R0)**2 * FS
      SHELL = C1*CS*FS + (1.-C1)*CP*FP
      RETURN
      END

      FUNCTION PROFNUC (B, JA)
C...This function return
C.  the profile T(b) for a nucleus of mass number A
C.  INPUT B = impact parameter (GeV**-1)
C.        JA = integer mass number
C.  OUTPUT  (fm**-2)
C.
C.  The  density of the nucleus is the `shell model density'
C.  the parameter r0 must beinitialized in the common block
C.............................................................
      SAVE

      COMMON /CSHELL/ RR0(18), RR02(18)
      DATA PI /3.1415926/
      B2 = B*B
      ARG = B2/RR02(JA)
      TS = EXP(-ARG)
      TP = TS*(2.*B2+RR02(JA))/(3.*RR02(JA))
      CS = MIN(1.,4./FLOAT(JA))
      PROFNUC = (CS*TS + (1.-CS)*TP)/(PI*RR02(JA))
      RETURN
      END

      SUBROUTINE SHELL_INI
C...Initialize the parameter  of the shell model
C.  for the nuclei with    6 < A < 18
C..............................................
      SAVE

      COMMON /CSHELL/ RR0(18), RR02(18)
      DIMENSION RR(18)
C...Data on Sqrt[<r**2>]  in fermi
      DATA RR /0.81,2.095,1.88,1.674, -1.,2.56,2.41,-1.,2.519,2.45
     +          ,2.37, 2.460, 2.440, 2.54, 2.58, 2.718, 2.662,2.789 /
      DO JA=1,18
         A = FLOAT(JA)
         RMED = RR(JA)
         IF (RMED .LE. 0.)   RMED = 0.5*(RR(JA-1) + RR(JA+1))
         C = MAX(1.5,(5./2. - 4./A) )
         R0 = RMED/SQRT(C)
         RR0 (JA) = R0
         RR02(JA) = R0*R0
      ENDDO
      RETURN
      END
C->
      SUBROUTINE WOOD_SAXON_INI
      COMMON /CWOOD/ RR0(19:56), AA0(19:56), CC0(19:56)
      SAVE

      DATA PI /3.1415926/
C...Wood-Saxon parameters from  table 6.2   of Barrett and Jackson
      RR0 (19) = 2.59
      AA0 (19) = 0.564
      RR0 (20) = 2.74
      AA0 (20) = 0.569
      RR0 (22) = 2.782
      AA0 (22) = 0.549
      RR0 (24) = 2.99
      AA0 (24) = 0.548
      RR0 (27) = 2.84
      AA0 (27) = 0.569
      RR0 (28) = 3.14
      AA0 (28) = 0.537
      RR0 (29) = 3.77
      AA0 (29) = 0.52
      RR0 (48) = 3.912
      AA0 (48) = 0.5234
      RR0 (56) = 3.98
      AA0 (56) = 0.569
      DO J=19, 56
         IF (RR0(J) .LE. 0.)  THEN
            RR0(J) = 1.05*FLOAT(J)**0.333333
            AA0(J) = 0.545
         ENDIF
         CC0(J)=3./(4.*PI*RR0(J)**3)/(1.+((AA0(J)*PI)/RR0(J))**2)
      ENDDO
      RETURN
      END

      FUNCTION DENSA (Z)
C....Woods Saxon nuclear density (normalised to 1)
C.   for a nucleus of mass number A.
C.   INPUT z = z coordinate (fm)
C.         JA = integer mass number
C.         B (in common /CC01/)  impact parameter  (fm)
C.  OUTPUT (fm**-3)
C--------------------------------------------------------
      SAVE

      COMMON /CC01/  B
      COMMON /CCDA/ JA
      COMMON /CWOOD/ RR0(19:56), AA0(19:56), CC0(19:56)
      R = SQRT (Z*Z + B*B)
      DENSA = CC0(JA)/(1.+EXP((R-RR0(JA))/AA0(JA)))
      RETURN
      END

C==========================================================================
C. Cross sections
C==========================================================================

      SUBROUTINE SIGMA_AIR (IB,SIG0,SIGEL,NINT,
     +                            SIGMA,DSIGMA,SIGQE,DSIGQE)
C...Compute with a montecarlo method the "production"
C.  and "quasi-elastic" cross section for  
C.  a nucleus-air  interaction 
C.
C.  INPUT : IB            = mass of projectile nucleus
C.          SIG0 (mbarn)  = inelastic pp cross section
C.          NINT            = number  of interactions to generate
C.  OUTPUT : SIGMA (mbarn) = "production" cross section
C.           DSIGMA   "    = error
C.           SIGQE    "    = "quasi-elastic" cross section
C.           DSIGQE   "    = error
C.           additional output is in the common block  /CPROBAB/
C..........................................................................
      SAVE

      PARAMETER (IAMAX=56)
      PARAMETER (IAMAX2=3136)          ! IAMAX*IAMAX
      COMMON  /CPROBAB/ PROBA(IAMAX), DPROBA(IAMAX), 
     +   PROBB(IAMAX), DPROBB(IAMAX), PROBI(IAMAX2), DPROBI(IAMAX2),
     +   P1AEL(0:IAMAX),DP1AEL(0:IAMAX),P1BEL(0:IAMAX), DP1BEL(0:IAMAX),
     +   P2AEL(0:IAMAX),DP2AEL(0:IAMAX),P2BEL(0:IAMAX), DP2BEL(0:IAMAX)
      COMMON /CNUCMS/ B, BMAX, NTRY, NA, NB, NI, NAEL, NBEL
     +         ,JJA(IAMAX), JJB(IAMAX), JJINT(IAMAX,IAMAX)
     +         ,JJAEL(IAMAX), JJBEL(IAMAX)
      DIMENSION  MMA(0:IAMAX), MMB(0:IAMAX), MMI(0:IAMAX2)
      DIMENSION  M1AEL(0:IAMAX), M1BEL(0:IAMAX)
      DIMENSION  M2AEL(0:IAMAX), M2BEL(0:IAMAX)
      DATA WOX /0.346/
      DATA PI /3.1415926/
      R2 = 0.1 * SIG0/PI
      BMAX = 15.                             ! fm
      SIGMA0 = PI*BMAX*BMAX*10.              ! mbarn
      IA = 16
      DO J=1,IA
         MMA(J) = 0
         M1AEL(J) = 0
         M2AEL(J) = 0
      ENDDO
      DO J=1,IB
         MMB(J) = 0
         M1BEL(J) = 0
         M2BEL(J) = 0
      ENDDO
      DO J=1,IA*IB
         MMI(J) = 0
      ENDDO
      NN = 0
      M = 0
      DO KK=1,NINT
         IA = 14 + 2*INT((1.+WOX)*S_RNDM(0))
         CALL INT_NUC (IA, IB, SIG0, SIGEL) 
         NN = NN + NTRY
         MMI(NI) = MMI(NI) + 1
         MMA(NA) = MMA(NA)+1
         MMB(NB) = MMB(NB)+1
         IF (NI .GT. 0)  THEN
            M = M+1
            M1AEL(NAEL) = M1AEL(NAEL)+1
            M1BEL(NBEL) = M1BEL(NBEL)+1
         ELSE
            M2AEL(NAEL) = M2AEL(NAEL)+1
            M2BEL(NBEL) = M2BEL(NBEL)+1
         ENDIF
      ENDDO
      MQE = NINT - M
      SIGMA  = SIGMA0 * FLOAT(M)/FLOAT(NN)
      DSIGMA = SIGMA0 * SQRT(FLOAT(M))/FLOAT(NN)
      SIGQE  = SIGMA0 * FLOAT(MQE)/FLOAT(NN)
      DSIGQE = SIGMA0 * SQRT(FLOAT(MQE))/FLOAT(NN)
      DO J=1,IA
         PROBA(J) = FLOAT(MMA(J))/FLOAT(M)
         DPROBA(J) = SQRT(FLOAT(MMA(J)))/FLOAT(M)
      ENDDO
      DO J=1,IB
         PROBB(J) = FLOAT(MMB(J))/FLOAT(M)
         DPROBB(J) = SQRT(FLOAT(MMB(J)))/FLOAT(M)
      ENDDO
      DO J=1,IA*IB
         PROBI(J) = FLOAT(MMI(J))/FLOAT(M)
         DPROBI(J) = SQRT(FLOAT(MMI(J)))/FLOAT(M)
      ENDDO
      DO J=0,IA
         P1AEL(J) = FLOAT(M1AEL(J))/FLOAT(M)
         DP1AEL(J) = SQRT(FLOAT(M1AEL(J)))/FLOAT(M)
         P2AEL(J) = FLOAT(M2AEL(J))/FLOAT(MQE)
         DP2AEL(J) = SQRT(FLOAT(M2AEL(J)))/FLOAT(MQE)
      ENDDO
      DO J=0,IB
         P1BEL(J) = FLOAT(M1BEL(J))/FLOAT(M)
         DP1BEL(J) = SQRT(FLOAT(M1BEL(J)))/FLOAT(M)
         P2BEL(J) = FLOAT(M2BEL(J))/FLOAT(MQE)
         DP2BEL(J) = SQRT(FLOAT(M2BEL(J)))/FLOAT(MQE)
      ENDDO
      RETURN
      END
C->
      SUBROUTINE SIGMA_MC (IA,IB,SIG0,SIGEL,NINT,
     +                            SIGMA,DSIGMA,SIGQE,DSIGQE)
C...Compute with a montecarlo method the "production"
C.  and "quasi-elastic" cross section for  
C.  a nucleus-nucleus interaction
C.
C.  INPUT : IA            = mass of target nucleus
C.          IB            = mass of projectile nucleus
C.          SIG0 (mbarn)  = inelastic pp cross section
C.          NINT            = number  of interactions to generate
C.  OUTPUT : SIGMA (mbarn) = "production" cross section
C.           DSIGMA   "    = error
C.           SIGQE    "    = "quasi-elastic" cross section
C.           DSIGQE   "    = error
C.           additional output is in the common block  /CPROBAB/
C.           Prob(n_A), Prob(n_B), Prob(n_int)
C..........................................................................
      SAVE

      PARAMETER (IAMAX=56)
      PARAMETER (IAMAX2=3136)          ! IAMAX*IAMAX
      COMMON  /CPROBAB/ PROBA(IAMAX), DPROBA(IAMAX), 
     +   PROBB(IAMAX), DPROBB(IAMAX), PROBI(IAMAX2), DPROBI(IAMAX2),
     +   P1AEL(0:IAMAX),DP1AEL(0:IAMAX),P1BEL(0:IAMAX), DP1BEL(0:IAMAX),
     +   P2AEL(0:IAMAX),DP2AEL(0:IAMAX),P2BEL(0:IAMAX), DP2BEL(0:IAMAX)
      COMMON /CNUCMS/ B, BMAX, NTRY, NA, NB, NI, NAEL, NBEL
     +         ,JJA(IAMAX), JJB(IAMAX), JJINT(IAMAX,IAMAX)
     +         ,JJAEL(IAMAX), JJBEL(IAMAX)
      DIMENSION  MMA(0:IAMAX), MMB(0:IAMAX), MMI(0:IAMAX2)
      DIMENSION  M1AEL(0:IAMAX), M1BEL(0:IAMAX)
      DIMENSION  M2AEL(0:IAMAX), M2BEL(0:IAMAX)
      DATA PI /3.1415926/
      R2 = 0.1 * SIG0/PI
      BMAX = 15.                             ! fm
      SIGMA0 = PI*BMAX*BMAX*10.              ! mbarn
      DO J=1,IA
         MMA(J) = 0
         M1AEL(J) = 0
         M2AEL(J) = 0
      ENDDO
      DO J=1,IB
         MMB(J) = 0
         M1BEL(J) = 0
         M2BEL(J) = 0
      ENDDO
      DO J=1,IA*IB
         MMI(J) = 0
      ENDDO
      NN = 0
      M = 0
      DO KK=1,NINT
         CALL INT_NUC (IA, IB, SIG0, SIGEL) 
         NN = NN + NTRY
         MMI(NI) = MMI(NI) + 1
         MMA(NA) = MMA(NA)+1
         MMB(NB) = MMB(NB)+1
         IF (NI .GT. 0)  THEN
            M = M+1
            M1AEL(NAEL) = M1AEL(NAEL)+1
            M1BEL(NBEL) = M1BEL(NBEL)+1
         ELSE
            M2AEL(NAEL) = M2AEL(NAEL)+1
            M2BEL(NBEL) = M2BEL(NBEL)+1
         ENDIF
      ENDDO
      MQE = NINT - M
      SIGMA  = SIGMA0 * FLOAT(M)/FLOAT(NN)
      DSIGMA = SIGMA0 * SQRT(FLOAT(M))/FLOAT(NN)
      SIGQE  = SIGMA0 * FLOAT(MQE)/FLOAT(NN)
      DSIGQE = SIGMA0 * SQRT(FLOAT(MQE))/FLOAT(NN)
      DO J=1,IA
         PROBA(J) = FLOAT(MMA(J))/FLOAT(M)
         DPROBA(J) = SQRT(FLOAT(MMA(J)))/FLOAT(M)
      ENDDO
      DO J=1,IB
         PROBB(J) = FLOAT(MMB(J))/FLOAT(M)
         DPROBB(J) = SQRT(FLOAT(MMB(J)))/FLOAT(M)
      ENDDO
      DO J=1,IA*IB
         PROBI(J) = FLOAT(MMI(J))/FLOAT(M)
         DPROBI(J) = SQRT(FLOAT(MMI(J)))/FLOAT(M)
      ENDDO
      DO J=0,IA
         P1AEL(J) = FLOAT(M1AEL(J))/FLOAT(M)
         DP1AEL(J) = SQRT(FLOAT(M1AEL(J)))/FLOAT(M)
         P2AEL(J) = FLOAT(M2AEL(J))/FLOAT(MQE)
         DP2AEL(J) = SQRT(FLOAT(M2AEL(J)))/FLOAT(MQE)
      ENDDO
      DO J=0,IB
         P1BEL(J) = FLOAT(M1BEL(J))/FLOAT(M)
         DP1BEL(J) = SQRT(FLOAT(M1BEL(J)))/FLOAT(M)
         P2BEL(J) = FLOAT(M2BEL(J))/FLOAT(MQE)
         DP2BEL(J) = SQRT(FLOAT(M2BEL(J)))/FLOAT(MQE)
      ENDDO
      RETURN
      END

C=============================================================
C.  Cross sections
C=============================================================

      SUBROUTINE SIG_H_AIR2
     +     (SSIG,SLOPE,ALPHA,ALAM,SIGT,SIGEL,SIGQE,SIGSD,SIGQSD)
C**********************************************************************
C...Subroutine to compute hadron-air cross sections
C.  according to:
C.  R.J. Glauber and G.Matthiae  Nucl.Phys. B21, 135, (1970)
C.
C.  Air is a linear combination of Nitrogen and oxygen
C.
C.  INPUT :  SSIG  (mbarn) total pp cross section
C.           SLOPE (GeV**-2)  elastic scattering slope for pp
C.           ALPHA    real/imaginary part of the forward pp elastic
C.                                               scattering amplitude
C.  OUTPUT : SIGT  = Total cross section
C.           SIGEL = Elastic cross section
C.           SIGQEL  = Elastic + Quasi elastic cross section
C.           SIGSD   = single diff. cross section (beam) 
C.           SIGQSD  = Elastic + Quasi elastic SD cross section (beam)
C.
C.  ALSO including interface from single precision in SIBYLL to
C.       double precision in GLAUBER2
C......................................................................
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      DATA  FOX /0.257/
      DOUBLE PRECISION DSSIG,DSLOPE,DALPHA,DALAM
      DOUBLE PRECISION SIG1,SIGEL1,SIGQE1,SIGSD1,SIGQSD1
      DOUBLE PRECISION SIG2,SIGEL2,SIGQE2,SIGSD2,SIGQSD2
      DSSIG = SSIG
      DSLOPE = SLOPE
      DALPHA = ALPHA
      DALAM = ALAM
      CALL GLAUBER2
     +  (14,DSSIG,DSLOPE,DALPHA,DALAM,SIG1,SIGEL1,SIGQE1,SIGSD1,SIGQSD1)
      CALL GLAUBER2
     +  (16,DSSIG,DSLOPE,DALPHA,DALAM,SIG2,SIGEL2,SIGQE2,SIGSD2,SIGQSD2)
      SIGT  = (1.-FOX)*SIG1   + FOX*SIG2
      SIGEL = (1.-FOX)*SIGEL1 + FOX*SIGEL2
      SIGQE = (1.-FOX)*SIGQE1 + FOX*SIGQE2
      SIGSD = (1.-FOX)*SIGSD1 + FOX*SIGSD2
      SIGQSD = (1.-FOX)*SIGQSD1 + FOX*SIGQSD2
      RETURN
      END


      SUBROUTINE GLAUBER2
     +     (JA,SSIG,SLOPE,ALPHA,ALAM,SIGT,SIGEL,SIGQEL,SIGSD,SIGQSD)
C-----------------------------------------------------------------------
C...Subroutine to compute hadron-Nucleus cross sections
C.  according to:
C.  R.J. Glauber and G.Matthiae  Nucl.Phys. B21, 135, (1970)
C.
C.  This formulas assume that the target nucleus  density is
C.  modeled by a shell-model form.  A reasonable range of models
C.  is  4 < JA < 18
C.
C.  This is a modified version with a two-channel model for inelastic
C.  intermediate states of low mass (R. Engel 2012/03/26)
C.
C.  INPUT :  A = mass number of the nucleus
C.           SSIG  (mbarn) total pp cross section
C.           SLOPE (GeV**-2)  elastic scattering slope for pp
C.           ALAM  enhancement factor (sqrt of sigma_sd1/sigma_ela)
C.           ALPHA    real/imaginary part of the forward pp elastic
C.                                               scattering amplitude
C.  OUTPUT : SIGT  = Total cross section
C.           SIGEL = Elastic cross section
C.           SIGQEL  = Elastic + Quasi elastic cross section
C.           SIGSD = single diff. cross section
C.           SIGQSD = Quasi single diff. cross section
C.
C. Internally  everything is computed in GeV (length = GeV**-1)
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CA0SH/ R0, R02
      SAVE
      COMPLEX  ZS1, ZS2, ZP1, ZP2, Z1, Z2, OM12
      DIMENSION RR(18)

      DATA CMBARN /0.389385D0/
      DATA PI /3.1415926D0/
      DATA BMAX /100.D0/            ! GeV**-1
      DATA NB /500/
C...data on Sqrt[<r**2>] (fm). (A=5,8 are not correct).
C   From Barett and Jackson
      DATA RR /0.81,2.095,1.88,1.674, 2.56,2.56,2.41,2.5,2.519,2.45
     +          ,2.37, 2.460, 2.440, 2.54, 2.58, 2.718, 2.662,2.789 /

      A = FLOAT(JA)
C...Parameter of shell model density
      R0 = RR(JA)/0.197D0/SQRT(5.D0/2.D0 - 4.D0/A)    ! GeV**-1
      R02 = R0*R0

      SIG1 = (1.D0+ALAM) * SSIG/CMBARN            ! GeV**-2
      SIG2 = (1.D0-ALAM) * SSIG/CMBARN
      SIG12 = SQRT((1.D0+ALAM)*(1.D0-ALAM)) * SSIG/CMBARN
      DB = BMAX/FLOAT(NB)
      SUM  = 0.D0
      SUM1 = 0.D0
      SUM2 = 0.D0
      SUM3 = 0.D0
      SUM4 = 0.D0
      DO JB=1,NB

        B = DB*(FLOAT(JB)-0.5D0)

        GS1 = GLAUBGS_D (B,SLOPE, SIG1)
        XS1 = (1.D0- GS1)
        YS1 = GS1*ALPHA
        ZS1 = CMPLX(XS1,YS1)

        GP1 = GLAUBGP_D (B,SLOPE, SIG1)
        XP1 = (1.D0- GP1)
        YP1 = GP1*ALPHA
        ZP1 = CMPLX(XP1,YP1)

        Z1 = ZS1**4.D0 * ZP1**(A-4.D0)

        GS2 = GLAUBGS_D (B,SLOPE, SIG2)
        XS2 = (1.D0- GS2)
        YS2 = GS2*ALPHA
        ZS2 = CMPLX(XS2,YS2)

        GP2 = GLAUBGP_D (B,SLOPE, SIG2)
        XP2 = (1.D0- GP2)
        YP2 = GP2*ALPHA
        ZP2 = CMPLX(XP2,YP2)

        Z2 = ZS2**4.D0 * ZP2**(A-4.D0)

        XZ = 0.5D0 * REAL(Z1+Z2)
        YZ = 0.5D0 * AIMAG(Z1+Z2)

        XZ2 = 0.5D0 * REAL(Z2-Z1)
        YZ2 = 0.5D0 * AIMAG(Z2-Z1)

        SUM = SUM + (1.D0-XZ)*B

        SUM1 = SUM1 + ((1.D0-XZ)**2 + YZ**2)*B

        SUM3 = SUM3 + (XZ2**2 + YZ2**2)*B

        OMS1 = OMEGAS_D(B,SIG1,SLOPE,ALPHA)
        OMS2 = OMEGAS_D(B,SIG2,SLOPE,ALPHA)
        OMS12 = OMEGAS_D(B,SIG12,SLOPE,ALPHA)

        OMP1 = OMEGAP_D(B,SIG1,SLOPE,ALPHA)
        OMP2 = OMEGAP_D(B,SIG2,SLOPE,ALPHA)
        OMP12 = OMEGAP_D(B,SIG12,SLOPE,ALPHA)

        OM1 = (1.D0 - 2.D0*GS1 + OMS1)**4.D0
     &      * (1.D0 - 2.D0*GP1 + OMP1)**(A-4.D0)
        OM2 = (1.D0 - 2.D0*GS2 + OMS2)**4.D0
     &      * (1.D0 - 2.D0*GP2 + OMP2)**(A-4.D0)
        OM12 = (1.D0 - GS1*CMPLX(1.D0,ALPHA) - GS2*CMPLX(1.D0,-ALPHA)
     &               + OMS12)**4.D0
     &       * (1.D0 - GP1*CMPLX(1.D0,ALPHA) - GP2*CMPLX(1.D0,-ALPHA)
     &               + OMP12)**(A-4.D0)
        SUM2 = SUM2 + (1.D0-2.D0*XZ + (OM1+OM2)/4.D0
     &                 + REAL(OM12)/2.D0)*B
        SUM4 = SUM4 + ((OM1+OM2)/4.D0
     &                 - REAL(OM12)/2.D0)*B

      ENDDO

      SIGT =   SUM  * DB * 4.D0*PI * CMBARN
      SIGEL =  SUM1 * DB * 2.D0*PI * CMBARN
      SIGQEL = SUM2 * DB * 2.D0*PI * CMBARN
      SIGSD =  SUM3 * DB * 2.D0*PI * CMBARN
      SIGQSD = SUM4 * DB * 2.D0*PI * CMBARN
      END


      FUNCTION GLAUBGS_D (B,SLOPE, SIG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CA0SH/ A0, A02
      SAVE
      DATA PI /3.1415926/
      GAMMA2 = A02/4. + 0.5*SLOPE
      ARG = B**2/(4.*GAMMA2)
      GLAUBGS_D = SIG/(8.*PI*GAMMA2) * EXP(-ARG)
      RETURN
      END


      FUNCTION GLAUBGP_D (B,SLOPE, SIG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CA0SH/ A0, A02
      SAVE
      DATA PI /3.1415926/
      GAMMA2 = A02/4. + 0.5*SLOPE
      ARG = B**2/(4.*GAMMA2)
      C1 = 1.- A02/(6.*GAMMA2)*(1.-ARG)
      GLAUBGP_D = SIG/(8.*PI*GAMMA2) *  C1 * EXP(-ARG)
      RETURN
      END


      FUNCTION OMEGAS_D (B, SIG, SLOPE, RHO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CA0SH/ A0, A02
      SAVE
      DATA PI /3.1415926/
      ETA2 = 0.25*(A02 + SLOPE)
      F02 = SIG*SIG*(1.+RHO*RHO)/(16.*PI**2)
      ARG = -B*B/(4.*ETA2)
      OMEGAS_D = F02/(4.*ETA2*SLOPE) *EXP(ARG)
      RETURN
      END


      FUNCTION OMEGAP_D (B, SIG, SLOPE, RHO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CA0SH/ A0, A02
      SAVE
      DATA PI /3.1415926/
      ETA2 = 0.25*(A02 + SLOPE)
      F02 = SIG*SIG*(1.+RHO*RHO)/(16.*PI**2)
      ARG = -B*B/(4.*ETA2)
      OMEGAP_D=F02/(4.*ETA2*SLOPE)*(1.-A02/(6.*ETA2)*(1.+ARG))*EXP(ARG)
      RETURN
      END

C------------------------------------------------------------------------
C.  Fit of Block and Cahn to pp and pbar-p cross sections
C------------------------------------------------------------------------

      SUBROUTINE BLOCK(SQS,SIG1,SIG2,SLOP1,SLOP2,
     +                 RHO1,RHO2,SIGEL1,SIGEL2)
C...p-p and pbar-p cross sections
C.  Parametrization of  Block and Cahn
C
C.  INPUT  : SQS   (GeV)  = c.m. energy
C.  
C.  OUPUT : SIG1 (mbarn)    = pp  total  cross section 
C.          SLOP1 (GeV**2)  = slope of elastic scattering
C.          RHO1            = Real/Imaginary part of the amplitude
C.                            for forward elastic  scattering (pp)
C.          SIGEL1 (mbarn)  = pp  elastic scattering  cross section
C.          [1 -> 2   : pp -> pbar p]
C-----------------------------------------------------------------------
      SAVE

      DATA PI /3.1415926/
      DATA CMBARN /0.389385/
      S = SQS*SQS
      CALL FPLUS  (S, FR, FI)
      CALL FMINUS (S, GR, GI)
      SIG1 = FI-GI
      SIG2 = FI+GI
      RHO1 = (FR-GR)/(FI-GI)
      RHO2 = (FR+GR)/(FI+GI)
      CALL SSLOPE (S, BP, BM)
      SLOP1 = BP - GI/FI*(BM-BP)
      SLOP2 = BP + GI/FI*(BM-BP)
      SIGEL1 = SIG1**2*(1.+RHO1**2)/(16.*PI*SLOP1)/CMBARN
      SIGEL2 = SIG2**2*(1.+RHO2**2)/(16.*PI*SLOP2)/CMBARN
      RETURN
      END

      SUBROUTINE FPLUS (S, FR, FI)
      SAVE

      COMMON /BLOCKC/ AA, BETA, S0, CC, AMU, DD, ALPHA, A0
      COMPLEX Z1, Z2, Z3
      DATA PI /3.1415926/
      F1 = LOG(S/S0)
      Z1 = CMPLX(F1,-PI/2.)
      Z1 = Z1*Z1
      Z2 = 1. + A0*Z1
      Z3 = Z1/Z2
      F2 = CC*S**(AMU-1.)
      F3 = 0.5*PI*(1.-AMU)
      FI = AA + F2*COS(F3) + BETA*REAL(Z3)
      FR = -BETA*AIMAG(Z3)+F2*SIN(F3)
      RETURN
      END

      SUBROUTINE FMINUS (S, FR, FI)
      SAVE

      COMMON /BLOCKC/ AA, BETA, S0, CC, AMU, DD, ALPHA, A0
      DATA PI /3.1415926/
      F1 = S**(ALPHA-1.)
      F2 = 0.5*PI*(1.-ALPHA)
      FR = -DD*F1*COS(F2)
      FI = -DD*F1*SIN(F2)
      RETURN
      END

      SUBROUTINE SSLOPE (S, BP, BM)
      SAVE

      COMMON /BLOCKD/ CP, DP, EP, CM, DM
      AL = LOG(S)
      BP = CP + DP*AL + EP*AL*AL
      BM = CM + DM*AL
      RETURN
      END

      SUBROUTINE BLOCK_INI
C...Parameters of fit IFIT=1 of Block and Cahn
      SAVE

      COMMON /BLOCKC/ AA, BETA, S0, CC, AMU, DD, ALPHA, A0
      COMMON /BLOCKD/ CP, DP, EP, CM, DM
      AA = 41.74
      BETA = 0.66
      S0 = 338.5
      CC = 0.
      AMU = 0.
      DD = -39.37
      ALPHA = 0.48
      A0 = 0.
      CP = 10.90
      DP = -0.08
      EP = 0.043
      CM = 23.27
      DM = 0.93
      RETURN
      END

C=============================================================
C.  Nucleus-nucleus cross sections
C=============================================================

      SUBROUTINE SIGNUC_INI (IA,E0)
C...This subroutine receives in INPUT E0 (TeV)
C.  energy per nucleon and computes the cross sections
C.  and interactions lengths for  all nuclei
C.  with A  between 2 and IA
C.  The output is contained in common block /CLENNN/
C.
C.  Attention: the tabulated cross sections are obtained with
C.  new p-p cross sections as used in SIBYLL 2x,
C.  in addition field dimensions changed (RE 04/2000)
C.
C........................................................
      COMMON /CLENNN/ SSIGNUC(60), ALNUC(60)
      DIMENSION SIGMA(6,56), SIGQE(6,56)
      DIMENSION AA(6)
      DATA NE /6/, AMIN /1./, DA /1./
      DATA AA /1.,2.,3.,4.,5.,6./
      DATA AVOG /6.0221367E-04/
      DATA ATARGET /14.514/               ! effective masss of air
C...Data on `inelastic-production' nucleus-air cross section
      DATA (SIGMA(J, 2),J=1,6) /
     &3.936E+02,4.276E+02,5.021E+02,5.985E+02,6.767E+02,7.639E+02/
      DATA (SIGMA(J, 3),J=1,6) /
     &4.566E+02,5.029E+02,5.746E+02,6.725E+02,7.578E+02,8.563E+02/
      DATA (SIGMA(J, 4),J=1,6) /
     &4.800E+02,5.300E+02,6.020E+02,6.965E+02,7.993E+02,8.805E+02/
      DATA (SIGMA(J, 5),J=1,6) /
     &5.936E+02,6.330E+02,7.166E+02,8.050E+02,9.138E+02,1.039E+03/
      DATA (SIGMA(J, 6),J=1,6) /
     &6.960E+02,7.343E+02,8.397E+02,9.276E+02,1.053E+03,1.158E+03/
      DATA (SIGMA(J, 7),J=1,6) /
     &7.046E+02,7.528E+02,8.420E+02,9.438E+02,1.076E+03,1.174E+03/
      DATA (SIGMA(J, 8),J=1,6) /
     &7.492E+02,7.967E+02,9.013E+02,1.016E+03,1.111E+03,1.263E+03/
      DATA (SIGMA(J, 9),J=1,6) /
     &7.992E+02,8.358E+02,9.335E+02,1.058E+03,1.188E+03,1.287E+03/
      DATA (SIGMA(J, 10),J=1,6) /
     &7.957E+02,8.597E+02,9.500E+02,1.083E+03,1.206E+03,1.298E+03/
      DATA (SIGMA(J, 11),J=1,6) /
     &8.130E+02,8.696E+02,9.477E+02,1.051E+03,1.153E+03,1.287E+03/
      DATA (SIGMA(J, 12),J=1,6) /
     &8.590E+02,9.039E+02,1.000E+03,1.112E+03,1.256E+03,1.363E+03/
      DATA (SIGMA(J, 13),J=1,6) /
     &8.634E+02,9.138E+02,1.007E+03,1.122E+03,1.246E+03,1.367E+03/
      DATA (SIGMA(J, 14),J=1,6) /
     &9.321E+02,9.639E+02,1.052E+03,1.193E+03,1.281E+03,1.430E+03/
      DATA (SIGMA(J, 15),J=1,6) /
     &9.264E+02,1.006E+03,1.083E+03,1.197E+03,1.325E+03,1.456E+03/
      DATA (SIGMA(J, 16),J=1,6) /
     &9.858E+02,1.048E+03,1.124E+03,1.259E+03,1.414E+03,1.541E+03/
      DATA (SIGMA(J, 17),J=1,6) /
     &9.909E+02,1.042E+03,1.150E+03,1.254E+03,1.406E+03,1.560E+03/
      DATA (SIGMA(J, 18),J=1,6) /
     &1.051E+03,1.109E+03,1.205E+03,1.344E+03,1.426E+03,1.587E+03/
      DATA (SIGMA(J, 19),J=1,6) /
     &1.090E+03,1.147E+03,1.249E+03,1.392E+03,1.521E+03,1.656E+03/
      DATA (SIGMA(J, 20),J=1,6) /
     &1.125E+03,1.196E+03,1.310E+03,1.442E+03,1.600E+03,1.721E+03/
      DATA (SIGMA(J, 21),J=1,6) /
     &1.195E+03,1.201E+03,1.330E+03,1.465E+03,1.599E+03,1.757E+03/
      DATA (SIGMA(J, 22),J=1,6) /
     &1.156E+03,1.221E+03,1.301E+03,1.467E+03,1.602E+03,1.741E+03/
      DATA (SIGMA(J, 23),J=1,6) /
     &1.212E+03,1.266E+03,1.369E+03,1.501E+03,1.635E+03,1.811E+03/
      DATA (SIGMA(J, 24),J=1,6) /
     &1.219E+03,1.294E+03,1.411E+03,1.530E+03,1.688E+03,1.834E+03/
      DATA (SIGMA(J, 25),J=1,6) /
     &1.253E+03,1.290E+03,1.452E+03,1.554E+03,1.717E+03,1.843E+03/
      DATA (SIGMA(J, 26),J=1,6) /
     &1.271E+03,1.308E+03,1.425E+03,1.578E+03,1.748E+03,1.866E+03/
      DATA (SIGMA(J, 27),J=1,6) /
     &1.251E+03,1.295E+03,1.418E+03,1.563E+03,1.676E+03,1.850E+03/
      DATA (SIGMA(J, 28),J=1,6) /
     &1.287E+03,1.369E+03,1.466E+03,1.621E+03,1.751E+03,1.905E+03/
      DATA (SIGMA(J, 29),J=1,6) /
     &1.404E+03,1.510E+03,1.615E+03,1.747E+03,1.925E+03,2.056E+03/
      DATA (SIGMA(J, 30),J=1,6) /
     &1.348E+03,1.398E+03,1.524E+03,1.651E+03,1.813E+03,1.963E+03/
      DATA (SIGMA(J, 31),J=1,6) /
     &1.379E+03,1.438E+03,1.537E+03,1.683E+03,1.823E+03,1.962E+03/
      DATA (SIGMA(J, 32),J=1,6) /
     &1.371E+03,1.439E+03,1.552E+03,1.698E+03,1.846E+03,2.028E+03/
      DATA (SIGMA(J, 33),J=1,6) /
     &1.399E+03,1.477E+03,1.607E+03,1.743E+03,1.873E+03,2.026E+03/
      DATA (SIGMA(J, 34),J=1,6) /
     &1.407E+03,1.466E+03,1.581E+03,1.746E+03,1.881E+03,2.041E+03/
      DATA (SIGMA(J, 35),J=1,6) /
     &1.406E+03,1.501E+03,1.612E+03,1.749E+03,1.887E+03,2.069E+03/
      DATA (SIGMA(J, 36),J=1,6) /
     &1.465E+03,1.518E+03,1.624E+03,1.763E+03,1.924E+03,2.087E+03/
      DATA (SIGMA(J, 37),J=1,6) /
     &1.466E+03,1.526E+03,1.650E+03,1.830E+03,1.964E+03,2.109E+03/
      DATA (SIGMA(J, 38),J=1,6) /
     &1.469E+03,1.548E+03,1.668E+03,1.827E+03,1.974E+03,2.133E+03/
      DATA (SIGMA(J, 39),J=1,6) /
     &1.511E+03,1.567E+03,1.670E+03,1.808E+03,2.012E+03,2.143E+03/
      DATA (SIGMA(J, 40),J=1,6) /
     &1.513E+03,1.591E+03,1.679E+03,1.843E+03,2.008E+03,2.164E+03/
      DATA (SIGMA(J, 41),J=1,6) /
     &1.510E+03,1.589E+03,1.729E+03,1.878E+03,2.040E+03,2.144E+03/
      DATA (SIGMA(J, 42),J=1,6) /
     &1.546E+03,1.591E+03,1.744E+03,1.889E+03,2.025E+03,2.161E+03/
      DATA (SIGMA(J, 43),J=1,6) /
     &1.554E+03,1.625E+03,1.747E+03,1.893E+03,2.053E+03,2.200E+03/
      DATA (SIGMA(J, 44),J=1,6) /
     &1.552E+03,1.653E+03,1.738E+03,1.920E+03,2.046E+03,2.257E+03/
      DATA (SIGMA(J, 45),J=1,6) /
     &1.574E+03,1.625E+03,1.797E+03,1.937E+03,2.133E+03,2.275E+03/
      DATA (SIGMA(J, 46),J=1,6) /
     &1.592E+03,1.661E+03,1.795E+03,1.938E+03,2.087E+03,2.265E+03/
      DATA (SIGMA(J, 47),J=1,6) /
     &1.610E+03,1.677E+03,1.802E+03,1.961E+03,2.114E+03,2.262E+03/
      DATA (SIGMA(J, 48),J=1,6) /
     &1.636E+03,1.720E+03,1.839E+03,1.960E+03,2.129E+03,2.308E+03/
      DATA (SIGMA(J, 49),J=1,6) /
     &1.644E+03,1.725E+03,1.816E+03,1.981E+03,2.133E+03,2.339E+03/
      DATA (SIGMA(J, 50),J=1,6) /
     &1.645E+03,1.710E+03,1.865E+03,1.990E+03,2.214E+03,2.317E+03/
      DATA (SIGMA(J, 51),J=1,6) /
     &1.671E+03,1.719E+03,1.871E+03,2.057E+03,2.217E+03,2.388E+03/
      DATA (SIGMA(J, 52),J=1,6) /
     &1.668E+03,1.762E+03,1.877E+03,2.012E+03,2.206E+03,2.394E+03/
      DATA (SIGMA(J, 53),J=1,6) /
     &1.694E+03,1.786E+03,1.886E+03,2.038E+03,2.217E+03,2.408E+03/
      DATA (SIGMA(J, 54),J=1,6) /
     &1.716E+03,1.791E+03,1.941E+03,2.051E+03,2.201E+03,2.419E+03/
      DATA (SIGMA(J, 55),J=1,6) /
     &1.729E+03,1.819E+03,1.947E+03,2.064E+03,2.245E+03,2.431E+03/
      DATA (SIGMA(J, 56),J=1,6) /
     &1.760E+03,1.850E+03,1.958E+03,2.138E+03,2.290E+03,2.440E+03/
C...Data on `quasi-elastic' nucleus-air cross section
      DATA (SIGQE(J, 2),J=1,6) /
     &4.002E+01,4.131E+01,5.277E+01,8.935E+01,1.366E+02,1.941E+02/
      DATA (SIGQE(J, 3),J=1,6) /
     &4.340E+01,4.119E+01,5.503E+01,8.936E+01,1.404E+02,2.010E+02/
      DATA (SIGQE(J, 4),J=1,6) /
     &4.089E+01,3.946E+01,5.845E+01,8.997E+01,1.443E+02,1.942E+02/
      DATA (SIGQE(J, 5),J=1,6) /
     &4.674E+01,4.371E+01,6.512E+01,1.070E+02,1.602E+02,2.173E+02/
      DATA (SIGQE(J, 6),J=1,6) /
     &5.166E+01,4.870E+01,7.074E+01,1.110E+02,1.755E+02,2.371E+02/
      DATA (SIGQE(J, 7),J=1,6) /
     &4.874E+01,4.712E+01,7.104E+01,1.091E+02,1.727E+02,2.284E+02/
      DATA (SIGQE(J, 8),J=1,6) /
     &5.260E+01,5.158E+01,6.752E+01,1.122E+02,1.829E+02,2.378E+02/
      DATA (SIGQE(J, 9),J=1,6) /
     &5.492E+01,5.099E+01,7.059E+01,1.202E+02,1.737E+02,2.441E+02/
      DATA (SIGQE(J, 10),J=1,6) /
     &5.523E+01,4.927E+01,7.239E+01,1.225E+02,1.827E+02,2.491E+02/
      DATA (SIGQE(J, 11),J=1,6) /
     &5.061E+01,4.750E+01,7.002E+01,1.198E+02,1.772E+02,2.414E+02/
      DATA (SIGQE(J, 12),J=1,6) /
     &5.396E+01,5.170E+01,7.055E+01,1.153E+02,1.937E+02,2.544E+02/
      DATA (SIGQE(J, 13),J=1,6) /
     &5.044E+01,5.124E+01,7.360E+01,1.228E+02,1.799E+02,2.504E+02/
      DATA (SIGQE(J, 14),J=1,6) /
     &5.592E+01,4.903E+01,7.205E+01,1.162E+02,1.790E+02,2.533E+02/
      DATA (SIGQE(J, 15),J=1,6) /
     &5.330E+01,5.318E+01,7.377E+01,1.232E+02,1.940E+02,2.609E+02/
      DATA (SIGQE(J, 16),J=1,6) /
     &5.550E+01,5.341E+01,7.676E+01,1.276E+02,2.021E+02,2.508E+02/
      DATA (SIGQE(J, 17),J=1,6) /
     &5.867E+01,5.543E+01,7.889E+01,1.269E+02,1.894E+02,2.667E+02/
      DATA (SIGQE(J, 18),J=1,6) /
     &5.755E+01,5.386E+01,7.800E+01,1.322E+02,2.006E+02,2.785E+02/
      DATA (SIGQE(J, 19),J=1,6) /
     &6.176E+01,6.446E+01,8.212E+01,1.311E+02,2.023E+02,2.750E+02/
      DATA (SIGQE(J, 20),J=1,6) /
     &6.196E+01,6.324E+01,8.657E+01,1.458E+02,2.091E+02,2.876E+02/
      DATA (SIGQE(J, 21),J=1,6) /
     &6.649E+01,6.427E+01,8.414E+01,1.279E+02,2.071E+02,2.905E+02/
      DATA (SIGQE(J, 22),J=1,6) /
     &6.494E+01,5.181E+01,8.334E+01,1.438E+02,2.120E+02,2.881E+02/
      DATA (SIGQE(J, 23),J=1,6) /
     &5.991E+01,5.688E+01,8.417E+01,1.398E+02,2.192E+02,2.728E+02/
      DATA (SIGQE(J, 24),J=1,6) /
     &6.742E+01,6.296E+01,8.719E+01,1.417E+02,2.150E+02,2.815E+02/
      DATA (SIGQE(J, 25),J=1,6) /
     &6.357E+01,6.148E+01,8.514E+01,1.419E+02,2.218E+02,2.874E+02/
      DATA (SIGQE(J, 26),J=1,6) /
     &6.775E+01,6.150E+01,9.240E+01,1.467E+02,2.149E+02,2.972E+02/
      DATA (SIGQE(J, 27),J=1,6) /
     &6.210E+01,6.245E+01,7.905E+01,1.467E+02,2.069E+02,3.029E+02/
      DATA (SIGQE(J, 28),J=1,6) /
     &6.886E+01,6.138E+01,8.515E+01,1.396E+02,2.182E+02,2.989E+02/
      DATA (SIGQE(J, 29),J=1,6) /
     &7.329E+01,6.651E+01,9.455E+01,1.492E+02,2.292E+02,3.234E+02/
      DATA (SIGQE(J, 30),J=1,6) /
     &6.665E+01,6.710E+01,9.592E+01,1.412E+02,2.238E+02,2.877E+02/
      DATA (SIGQE(J, 31),J=1,6) /
     &6.981E+01,5.773E+01,8.913E+01,1.519E+02,2.315E+02,2.943E+02/
      DATA (SIGQE(J, 32),J=1,6) /
     &6.626E+01,6.984E+01,8.322E+01,1.398E+02,2.305E+02,2.988E+02/
      DATA (SIGQE(J, 33),J=1,6) /
     &6.871E+01,6.684E+01,9.337E+01,1.497E+02,2.186E+02,2.996E+02/
      DATA (SIGQE(J, 34),J=1,6) /
     &6.876E+01,6.505E+01,9.522E+01,1.465E+02,2.258E+02,3.064E+02/
      DATA (SIGQE(J, 35),J=1,6) /
     &6.674E+01,6.646E+01,8.897E+01,1.559E+02,2.313E+02,3.045E+02/
      DATA (SIGQE(J, 36),J=1,6) /
     &6.841E+01,6.359E+01,9.195E+01,1.577E+02,2.173E+02,2.959E+02/
      DATA (SIGQE(J, 37),J=1,6) /
     &6.762E+01,6.873E+01,9.512E+01,1.437E+02,2.231E+02,3.016E+02/
      DATA (SIGQE(J, 38),J=1,6) /
     &7.021E+01,6.484E+01,9.780E+01,1.552E+02,2.304E+02,3.080E+02/
      DATA (SIGQE(J, 39),J=1,6) /
     &6.955E+01,7.334E+01,9.661E+01,1.602E+02,2.273E+02,3.193E+02/
      DATA (SIGQE(J, 40),J=1,6) /
     &6.533E+01,6.682E+01,9.137E+01,1.546E+02,2.393E+02,3.044E+02/
      DATA (SIGQE(J, 41),J=1,6) /
     &6.952E+01,6.810E+01,9.831E+01,1.571E+02,2.450E+02,3.087E+02/
      DATA (SIGQE(J, 42),J=1,6) /
     &6.577E+01,6.096E+01,9.797E+01,1.532E+02,2.391E+02,3.297E+02/
      DATA (SIGQE(J, 43),J=1,6) /
     &7.203E+01,7.122E+01,9.700E+01,1.615E+02,2.405E+02,3.189E+02/
      DATA (SIGQE(J, 44),J=1,6) /
     &6.652E+01,7.105E+01,9.380E+01,1.593E+02,2.398E+02,3.142E+02/
      DATA (SIGQE(J, 45),J=1,6) /
     &7.277E+01,7.158E+01,9.596E+01,1.582E+02,2.343E+02,3.247E+02/
      DATA (SIGQE(J, 46),J=1,6) /
     &7.256E+01,6.919E+01,9.289E+01,1.617E+02,2.383E+02,3.357E+02/
      DATA (SIGQE(J, 47),J=1,6) /
     &7.995E+01,7.389E+01,1.033E+02,1.666E+02,2.284E+02,3.290E+02/
      DATA (SIGQE(J, 48),J=1,6) /
     &7.566E+01,7.109E+01,9.538E+01,1.576E+02,2.352E+02,3.070E+02/
      DATA (SIGQE(J, 49),J=1,6) /
     &7.672E+01,6.647E+01,1.053E+02,1.604E+02,2.384E+02,3.132E+02/
      DATA (SIGQE(J, 50),J=1,6) /
     &7.121E+01,7.144E+01,9.875E+01,1.660E+02,2.408E+02,3.337E+02/
      DATA (SIGQE(J, 51),J=1,6) /
     &7.745E+01,7.760E+01,9.869E+01,1.591E+02,2.302E+02,3.216E+02/
      DATA (SIGQE(J, 52),J=1,6) /
     &7.076E+01,7.458E+01,1.031E+02,1.610E+02,2.616E+02,3.236E+02/
      DATA (SIGQE(J, 53),J=1,6) /
     &7.429E+01,7.034E+01,9.657E+01,1.621E+02,2.491E+02,3.098E+02/
      DATA (SIGQE(J, 54),J=1,6) /
     &6.706E+01,7.424E+01,9.635E+01,1.658E+02,2.383E+02,3.301E+02/
      DATA (SIGQE(J, 55),J=1,6) /
     &7.373E+01,6.657E+01,9.151E+01,1.582E+02,2.497E+02,3.543E+02/
      DATA (SIGQE(J, 56),J=1,6) /
     &7.335E+01,6.709E+01,1.011E+02,1.671E+02,2.426E+02,3.327E+02/

      ASQS = 0.5*LOG10(1.876E+03*E0)
      JE = MIN(INT((ASQS-AMIN)/DA)+1,NE-2)
      DO JA=2,IA
         ABEAM = FLOAT(JA)
         S1 = QUAD_INT(ASQS, AA(JE),AA(JE+1),AA(JE+2),
     +                   SIGMA(JE,JA),SIGMA(JE+1,JA),SIGMA(JE+2,JA))
         S2 = QUAD_INT(ASQS, AA(JE),AA(JE+1),AA(JE+2),
     +                   SIGQE(JE,JA),SIGQE(JE+1,JA),SIGQE(JE+2,JA))
         SSIGNUC(JA) = S1 + S2
         ALNUC(JA) = ATARGET/(AVOG*SSIGNUC(JA))
      ENDDO
      ALNUC(1) = FPNI(E0, 13)
      SSIGNUC(1) = ATARGET/(AVOG*ALNUC(1))

      RETURN
      END


C=======================================================================
C.  General utilities
C.======================================================================

      FUNCTION QUAD_INT (R,X0,X1,X2,V0,V1,V2)
C...Quadratic interpolation
      SAVE

      R0=R-X0
      R1=R-X1
      R2=R-X2
      S0=X0-X1
      S1=X0-X2
      S2=X1-X2
      QUAD_INT = V0*R1*R2/(S0*S1)-V1*R0*R2/(S0*S2)+V2*R0*R1/(S1*S2)
      RETURN
      END

      FUNCTION GAUSS (FUN, A,B)
C...Returns the  8 points Gauss-Legendre integral
C.  of function FUN from A to B
C...........................................................
      SAVE

      DIMENSION X(8), W(8)
      DATA X / .0950125098, .2816035507, .4580167776, .6178762444
     1          ,.7554044083, .8656312023, .9445750230, .9894009349/
      DATA W / .1894506104, .1826034150, .1691565193, .1495959888
     1          ,.1246289712, .0951585116, .0622535239, .0271524594/
      XM = 0.5*(B+A)
      XR = 0.5*(B-A)
      SS = 0.
      DO J=1,8
        DX = XR*X(J)
        SS = SS + W(J) * (FUN(XM+DX) + FUN(XM-DX))
      ENDDO
      GAUSS = XR*SS
      RETURN
      END

      subroutine invert_array (yy, xmin, dx, n, xnew, ymin, dy)
C..    This subroutine receives one   array
C      of n y values in input yy(1:n)
C      that correspond to  equispaced values of x_j = xmin + dx*(j-1)
C
C      and "reverse" the array returning an array of  x values
C      xnew (1:n) that  corresponds to equispaced values of y
C      The relation is assumed monotonous but can be 
C      increasing or decreasing
C..............................................................
      SAVE

      dimension  yy(n), xnew (n)
      ymin = yy(1)
      ymax = yy(n)
      dy = (ymax - ymin)/float(n-1)
      xnew (1) = xmin
      xnew (n) = xmin + dx*float(n-1)
      k0 = 1
      do j=2,n-1
         y = ymin + float(j-1)*dy 
         do k=k0,n
            if((yy(k) .gt. y) .eqv. (yy(n) .gt. yy(1))) goto 100
         enddo
100      y2 = yy(k)
         y1 = yy(k-1)
         k0 = k-1
         x1 = xmin + dx*float(k-2)
         x2 = x1+dx
         xnew (j)  = x1 + dx* (y-y1)/(y2-y1)
      enddo
      return
      end
C->
      SUBROUTINE SINCO(S,C)
      SAVE

      DATA PI /3.1415926/
      F = 2.*PI*S_RNDM(0)
      C = COS (F)
      S = SIN (F)
      RETURN
      END



C=============================================================
C.  Cross sections for cascade calculations (FPNI)
C=============================================================


      
      SUBROUTINE SIGMA_PP (E0, SIGT, SIGEL, SIGINEL, SLOPE, RHO) 
C-----------------------------------------------------------------------
C...p-p cross sections
C.
C.  this routine serves the purpose to calculate cascades with different 
C.  cross sections
C.
C. INPUT: E0 = Laboratory Energy  (TeV)
C. 
C. OUTPUT: SIGT = total cross section
C.         SIGEL = elastic cross section
C.         SIGINEL = inelastic cross section
C.         SLOPE = slope of elastic scattering (GeV**-2)
C.         RHO = Imaginary/Real part of forward elastic amplitude
C.   
C.  (old cross section tables end at 10^6 GeV)
C-----------------------------------------------------------------------
      SAVE

      DIMENSION SSIG0(51)
      DIMENSION SIGDIF(3)
      DATA ICSPA /0/
      DATA PI /3.1415926/
      DATA CMBARN /0.389385/

C...p-p inelastic cross sections (mbarn)
      DATA (SSIG0(J),J=1,51) /
     +      32.05,    32.06,    32.08,    32.13,    32.22,    32.36,
     +      32.56,    32.85,    33.24,    33.75,    34.37,    35.14,
     +      36.05,    37.12,    38.37,    39.78,    41.36,    43.13,
     +      45.07,    47.18,    49.47,    51.91,    54.54,    57.28,
     +      60.15,    63.15,    66.28,    69.48,    72.80,    76.22,
     +      79.71,    83.27,    86.87,    90.55,    94.26,    98.05,
     +     101.89,   105.75,   109.71,   113.65,   117.60,   121.55,
     +     125.53,   129.56,   133.60,   137.70,   141.77,   145.84,
     +     149.92,   154.02,   158.15/


      SQS = SQRT(2000.*0.938*E0)      

*  old standard NUCLIB/SIBYLL model

      IF(ICSPA.EQ.-1) THEN

        AL = LOG10(SQS)
        if(AL.le.1.) then
          SIGINEL = SSIG0(1)
        else
          J1 = (AL - 1.)*10. + 1
          J1 = min(J1,50)
          T = (AL-1.)*10. - FLOAT(J1-1)
          SIGINEL = SSIG0(J1)*(1.-T) + SSIG0(J1+1)*T
        endif
        CALL BLOCK(SQS,SIGT1,SIGT2,SLOP1,SLOP2,RHO1,RHO2,SIGEL1,SIGEL2)
        R = SIGEL1/SIGT1
        RHO = RHO1
        SIGT  = SIGINEL/(1.-R)
        SIGEL = SIGINEL*R/(1.-R)
        SLOPE = SIGT**2/(SIGEL * 16.*PI) * (1.+RHO1**2) /CMBARN

*  cross section as calculated in SIBYLL

      ELSE IF(ICSPA.EQ.0) THEN

        CALL SIB_SIGMA_HP(1,SQS,SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO)

*  Donnachie-Landshoff  (sig-tot)

      ELSE IF(ICSPA.EQ.1) THEN

        CALL BLOCK(SQS,SIGT1,SIGT2,SLOP1,SLOP2,RHO1,RHO2,
     +             SIGEL1,SIGEL2)
        R = SIGEL1/SIGT1
        RHO = RHO1

        DELDL = 0.0808
        EPSDL = -0.4525
        S = SQS*SQS
        SIGT = 21.7*S**DELDL+56.08*S**EPSDL
        SIGEL = R*SIGT
        SIGINEL = SIGT-SIGEL
        SLOPE = SIGT**2/(SIGEL * 16.*PI) * (1.+RHO**2) /CMBARN

*  Donnachie-Landshoff (sig-tot and sig-el)

      ELSE IF(ICSPA.EQ.2) THEN

        DELDL = 0.0808
        EPSDL = -0.4525
        S = SQS*SQS
        SIGT = 21.7*S**DELDL+56.08*S**EPSDL
        IMODEL = 1
        IF(IMODEL.EQ.1) THEN
          ALPHAP = 0.25D0
          SLOPE = 8.5D0+2.D0*ALPHAP*LOG(S)
        ELSE IF(IMODEL.EQ.2) THEN
          ALPHAP = 0.3D0
          SLOPE = 8.D0+2.D0*ALPHAP*LOG(S)
        ENDIF
        SIGEL = SIGT**2/(16.D0*PI*SLOPE*CMBARN)
        SIGINEL = SIGT-SIGEL
        RHO = 0.

*  geometrical scaling with Donnachie-Landshoff sig-tot

      ELSE IF(ICSPA.EQ.3) THEN

        R = 0.17D0

        DELDL = 0.0808
        EPSDL = -0.4525
        S = SQS*SQS
        SIGT = 21.7*S**DELDL+56.08*S**EPSDL

        SIGEL = R*SIGT
        SIGINEL = SIGT-SIGEL
        SLOPE = SIGT**2/(16*PI*SIGEL)/CMBARN
        RHO = 0.

      ENDIF

      RETURN
      END


      SUBROUTINE SIGMA_PIP (E0, SIGT, SIGEL, SIGINEL, SLOPE, RHO) 
C-----------------------------------------------------------------------
C...pi-p cross sections
C.
C.  this routine serves the purpose to calculate cascades with different 
C.  cross sections
C.
C. INPUT: E0 = Laboratory Energy  (TeV)
C. 
C. OUTPUT: SIGT = total cross section
C.         SIGEL = elastic cross section
C.         SIGINEL = inelastic cross section
C.         SLOPE = slope of elastic scattering (GeV**-2)
C.         RHO = Imaginary/Real part of forward elastic amplitude
C.
C.  (old cross section tables end at 10^6 GeV)
C-----------------------------------------------------------------------
      SAVE

      DIMENSION SSIG0(51)
      DIMENSION SIGDIF(3)
      DATA ICSPA /0/
      DATA PI /3.1415926/
      DATA CMBARN /0.389385/

C...pi-p inelastic cross sections (mbarn)
      DATA (SSIG0(J),J=1,51) /
     +      20.76,    20.78,    20.81,    20.88,    20.98,    21.13,
     +      21.33,    21.61,    21.96,    22.39,    22.92,    23.56,
     +      24.31,    25.18,    26.18,    27.32,    28.60,    30.04,
     +      31.64,    33.40,    35.34,    37.43,    39.72,    42.16,
     +      44.77,    47.56,    50.53,    53.66,    56.99,    60.50,
     +      64.17,    68.03,    72.05,    76.27,    80.67,    85.27,
     +      90.08,    95.04,   100.27,   105.65,   111.21,   116.94,
     +     122.87,   129.03,   135.37,   141.93,   148.62,   155.49,
     +     162.48,   169.60,   176.94/

      SQS = SQRT(2000.*0.938*E0)      

*  old standard NUCLIB/SIBYLL model

      IF(ICSPA.EQ.-1) THEN

        AL = LOG10(SQS)
        if(AL.le.1.) then
          SIGINEL = SSIG0(1)
        else
          J1 = (AL - 1.)*10. + 1
          J1 = min(J1,50)
          T = (AL-1.)*10. - FLOAT(J1-1)
          SIGINEL = SSIG0(J1)*(1.-T) + SSIG0(J1+1)*T
        endif
        CALL BLOCK(SQS,SIGT1,SIGT2,SLOP1,SLOP2,RHO1,RHO2,SIGEL1,SIGEL2)
        R = SIGEL1/SIGT1
        RHO = RHO1
        SIGT  = SIGINEL/(1.-R)
        SIGEL = SIGINEL*R/(1.-R)
        SLOPE = SIGT**2/(SIGEL * 16.*PI) * (1.+RHO1**2) /CMBARN

*  cross section as calculated in SIBYLL

      ELSE IF(ICSPA.EQ.0) THEN

        CALL SIB_SIGMA_HP(2,SQS,SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO)

*  Donnachie-Landshoff  (sig-tot)

      ELSE IF(ICSPA.EQ.1) THEN

        CALL BLOCK(SQS,SIGT1,SIGT2,SLOP1,SLOP2,RHO1,RHO2,
     +             SIGEL1,SIGEL2)
        R = SIGEL1/SIGT1
        RHO = RHO1

        DELDL = 0.0808
        EPSDL = -0.4525
        S = SQS*SQS
        SIGT = 13.63*S**DELDL+(36.02+27.56)/2.*S**EPSDL
        SIGEL = R*SIGT
        SIGINEL = SIGT-SIGEL
        SLOPE = SIGT**2/(SIGEL * 16.*PI) * (1.+RHO**2) /CMBARN

*  Donnachie-Landshoff (sig-tot and sig-el)

      ELSE IF(ICSPA.EQ.2) THEN

        DELDL = 0.0808
        EPSDL = -0.4525
        S = SQS*SQS
        SIGT = 13.63*S**DELDL+(36.02+27.56)/2.*S**EPSDL
        IMODEL = 1
        IF(IMODEL.EQ.1) THEN
          ALPHAP = 0.25D0
          SLOPE = 8.5D0+2.D0*ALPHAP*LOG(S)
        ELSE IF(IMODEL.EQ.2) THEN
          ALPHAP = 0.3D0
          SLOPE = 8.D0+2.D0*ALPHAP*LOG(S)
        ENDIF
        SIGEL = SIGT**2/(16.D0*PI*SLOPE*CMBARN)
        SIGINEL = SIGT-SIGEL
        RHO = 0.

*  geometrical scaling with Donnachie-Landshoff sig-tot

      ELSE IF(ICSPA.EQ.3) THEN

        R = 0.17D0

        DELDL = 0.0808
        EPSDL = -0.4525
        S = SQS*SQS
        SIGT = 13.63*S**DELDL+(36.02+27.56)/2.*S**EPSDL

        SIGEL = R*SIGT
        SIGINEL = SIGT-SIGEL
        SLOPE = SIGT**2/(16*PI*SIGEL)/CMBARN
        RHO = 0.

      ENDIF

      RETURN
      END


      SUBROUTINE SIGMA_KP (E0, SIGT, SIGEL, SIGINEL, SLOPE, RHO) 
C-----------------------------------------------------------------------
C...K-p cross sections
C.
C.  this routine serves the purpose to calculate cascades with different 
C.  cross sections
C.
C.  if old cross sections are selected then sigma_pi = sigma_K
C.
C. INPUT: E0 = Laboratory Energy  (TeV)
C. 
C. OUTPUT: SIGT = total cross section
C.         SIGEL = elastic cross section
C.         SIGINEL = inelastic cross section
C.         SLOPE = slope of elastic scattering (GeV**-2)
C.         RHO = Imaginary/Real part of forward elastic amplitude
C.
C.  (old cross section tables end at 10^6 GeV)
C-----------------------------------------------------------------------
      SAVE

      DIMENSION SSIG0(51)
      DIMENSION SIGDIF(3)
      DATA ICSPA /0/
      DATA PI /3.1415926/
      DATA CMBARN /0.389385/

C...pi-p inelastic cross sections (mbarn)
      DATA (SSIG0(J),J=1,51) /
     +      20.76,    20.78,    20.81,    20.88,    20.98,    21.13,
     +      21.33,    21.61,    21.96,    22.39,    22.92,    23.56,
     +      24.31,    25.18,    26.18,    27.32,    28.60,    30.04,
     +      31.64,    33.40,    35.34,    37.43,    39.72,    42.16,
     +      44.77,    47.56,    50.53,    53.66,    56.99,    60.50,
     +      64.17,    68.03,    72.05,    76.27,    80.67,    85.27,
     +      90.08,    95.04,   100.27,   105.65,   111.21,   116.94,
     +     122.87,   129.03,   135.37,   141.93,   148.62,   155.49,
     +     162.48,   169.60,   176.94/

      SQS = SQRT(2000.*0.938*E0)      

*  old standard NUCLIB/SIBYLL model

      IF(ICSPA.EQ.-1) THEN

        AL = LOG10(SQS)
        if(AL.le.1.) then
          SIGINEL = SSIG0(1)
        else
          J1 = (AL - 1.)*10. + 1
          J1 = min(J1,50)
          T = (AL-1.)*10. - FLOAT(J1-1)
          SIGINEL = SSIG0(J1)*(1.-T) + SSIG0(J1+1)*T
        endif
        CALL BLOCK(SQS,SIGT1,SIGT2,SLOP1,SLOP2,RHO1,RHO2,SIGEL1,SIGEL2)
        R = SIGEL1/SIGT1
        RHO = RHO1
        SIGT  = SIGINEL/(1.-R)
        SIGEL = SIGINEL*R/(1.-R)
        SLOPE = SIGT**2/(SIGEL * 16.*PI) * (1.+RHO1**2) /CMBARN

*  cross section as calculated in SIBYLL

      ELSE IF(ICSPA.EQ.0) THEN

        CALL SIB_SIGMA_HP(3,SQS,SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO)

*  Donnachie-Landshoff  (sig-tot)

      ELSE IF(ICSPA.EQ.1) THEN

        CALL BLOCK(SQS,SIGT1,SIGT2,SLOP1,SLOP2,RHO1,RHO2,
     +             SIGEL1,SIGEL2)
        R = SIGEL1/SIGT1
        RHO = RHO1

        DELDL = 0.0808
        EPSDL = -0.4525
        S = SQS*SQS
        SIGT = 11.82*S**DELDL+(26.36+ 8.15)/2.*S**EPSDL
        SIGEL = R*SIGT
        SIGINEL = SIGT-SIGEL
        SLOPE = SIGT**2/(SIGEL * 16.*PI) * (1.+RHO**2) /CMBARN

*  Donnachie-Landshoff (sig-tot and sig-el)

      ELSE IF(ICSPA.EQ.2) THEN

        DELDL = 0.0808
        EPSDL = -0.4525
        S = SQS*SQS
        SIGT = 11.82*S**DELDL+(26.36+ 8.15)/2.*S**EPSDL
        IMODEL = 1
        IF(IMODEL.EQ.1) THEN
          ALPHAP = 0.25D0
          SLOPE = 8.5D0+2.D0*ALPHAP*LOG(S)
        ELSE IF(IMODEL.EQ.2) THEN
          ALPHAP = 0.3D0
          SLOPE = 8.D0+2.D0*ALPHAP*LOG(S)
        ENDIF
        SIGEL = SIGT**2/(16.D0*PI*SLOPE*CMBARN)
        SIGINEL = SIGT-SIGEL
        RHO = 0.

*  geometrical scaling with Donnachie-Landshoff sig-tot

      ELSE IF(ICSPA.EQ.3) THEN

        R = 0.17D0

        DELDL = 0.0808
        EPSDL = -0.4525
        S = SQS*SQS
        SIGT = 11.82*S**DELDL+(26.36+ 8.15)/2.*S**EPSDL

        SIGEL = R*SIGT
        SIGINEL = SIGT-SIGEL
        SLOPE = SIGT**2/(16*PI*SIGEL)/CMBARN
        RHO = 0.

      ENDIF

      RETURN
      END



      SUBROUTINE SIGMA_INI 
C-----------------------------------------------------------------------
C.  Initialize the cross section and interaction lengths on air
C-----------------------------------------------------------------------
      SAVE

      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      COMMON /CSAIR/ ASQSMIN, ASQSMAX, DASQS,
     &               SSIG0(61,3),SSIGA(61,3),ALINT(61,3),NSQS
      COMMON /GLAUB_SCR/ XI_MAX , ALAM(61)
      DATA AVOG /6.0221367E-04/
      ATARGET = 14.514

      CALL BLOCK_INI

C...Loop on c.m. energy 
      NSQS = 61
      SQSMIN = 10.
      SQSMAX = 1.E+07
      ASQSMIN = LOG10(SQSMIN)
      ASQSMAX = LOG10(SQSMAX)
      DASQS = (ASQSMAX-ASQSMIN)/FLOAT(NSQS-1)
      DO J=1,NSQS
         ASQS = ASQSMIN + DASQS*FLOAT(J-1)
         SQS = 10.**ASQS
         E0 = SQS*SQS/(2.*0.938) * 1.E-03       ! TeV
C...p-air
         CALL SIGMA_PP (E0, SIGT, SIGEL, SIGINEL, SLOPE, RHO) 
         CALL SIG_H_AIR2 (SIGT, SLOPE, RHO, ALAM(J), SSIGT, SSIGEL, 
     +        SSIGQE, SSIGSD, SIGQSD)
         SSIGA(J,1) = SSIGT-SSIGQE
         SSIG0(J,1) = SIGINEL
         ALINT(J,1) = 1./(AVOG*SSIGA(J,1)/ATARGET)
C...pi-air
         CALL SIGMA_PIP (E0, SIGT, SIGEL, SIGINEL, SLOPE, RHO) 
         CALL  SIG_H_AIR2 (SIGT, SLOPE, RHO, ALAM(J), SSIGT, SSIGEL,
     +         SSIGQE, SSIGSD, SIGQSD)
         SSIGA(J,2) = SSIGT-SSIGQE
         SSIG0(J,2) = SIGINEL
         ALINT(J,2) = 1./(AVOG*SSIGA(J,2)/ATARGET)
C...K-air
         CALL SIGMA_KP (E0, SIGT, SIGEL, SIGINEL, SLOPE, RHO) 
         CALL  SIG_H_AIR2 (SIGT, SLOPE, RHO, ALAM(J), SSIGT, SSIGEL, 
     +        SSIGQE, SSIGSD, SIGQSD)
         SSIGA(J,3) = SSIGT-SSIGQE
         SSIG0(J,3) = SIGINEL
         ALINT(J,3) = 1./(AVOG*SSIGA(J,3)/ATARGET)
      ENDDO

      WRITE(LUN,'(1X,A)') 
     &  'SIGMA_INI: NUCLIB interaction lengths (p-air, pi-air, K-air)'
      DO J=1,NSQS
         SQS = 10.**(ASQSMIN + DASQS*FLOAT(J-1))
         WRITE(LUN,'(1X,1P,4E12.3)')SQS,ALINT(J,1),ALINT(J,2),ALINT(J,3)
      ENDDO

      RETURN
      END


      FUNCTION FPNI (E,Linp)
C-----------------------------------------------------------------------
C...This function  returns the interaction length 
C.  of an hadronic particle travelling in air
C.
C.  INPUT:   E (TeV)   particle energy
C.           Linp      particle code
C.  OUTPUT:  FPNI      (g cm-2)
C-----------------------------------------------------------------------
      SAVE

      COMMON /CSAIR/ ASQSMIN, ASQSMAX, DASQS,
     &               SSIG0(61,3),SSIGA(61,3),ALINT(61,3),NSQS

      DIMENSION KK(6:14)
      DATA KK /3*2, 4*3, 2*1/

      SQS = SQRT(2000.*E*0.937)                        ! GeV
      AL = LOG10 (SQS)
      L = abs(Linp)
      IF (AL .LE. ASQSMIN)  THEN
         FPNI = ALINT(1,KK(L))
      ELSE
         T = (AL-ASQSMIN)/DASQS
         J = INT(T)
         J = MIN(J,NSQS-2)
         T = T-FLOAT(J)
         FPNI = ((1.-T)*ALINT(J+1,KK(L)) + T*ALINT(J+2,KK(L)))
      ENDIF

      RETURN
      END



      SUBROUTINE INT_LEN_INI
C-----------------------------------------------------------------------
C...Initialize the interaction lengths from NUCLIB
C-----------------------------------------------------------------------
      SAVE

      CALL NUC_GEOM_INI                 ! nucleus profiles
      CALL SIGMA_INI                    ! initialize cross sections

      RETURN
      END

      SUBROUTINE PDG_INI
C----------------------------------------------------------------
C     PDG conversion blocks \FR'13
C----------------------------------------------------------------
      SAVE
      COMMON /S_DEBUG/ Ncall, Ndebug, lun
      PARAMETER ( ID_PDG_MAX = 99 )
      COMMON /S_PDG2PID/ ID_PDG_LIST(ID_PDG_MAX),ID_LIST(577)
      DATA ID_PDG_LIST /22,-11,11,-13,13,111,211,-211,321,-321, !10
     &     130,310,2212,2112,12,-12,14,-14,-2212,-2112,         !20
     &     311,-311,221,331,213,-213,113,10321,-10321,10311,    !30
     &     -10311,223,333,3222,3212,3112,3322,3312,3122,2224,   !40
     &     2214,2114,1114,3224,3214,3114,3324,3314,3334,0,      !50
     &     8*0,411,-411,10*0,                                   !70
     &     421,-421,441,431,-431,433,-433,413,-413,423,         !80
     &     -423,0,443,4222,4212,4112,4232,4132,4122,0,          !90
     &     0,0,0,4224,4214,4114,4324,4314,4332/

      IF(Ndebug.gt.2)
     & WRITE(lun,*) 'INITIALIZING PDG TABLES..'
      CALL pho_cpcini(ID_pdg_max,ID_pdg_list,ID_list)
      
      END

      INTEGER FUNCTION ISIB_PDG2PID(Npdg)
C----------------------------------------------------------------
C     conversion of PDG standard particle code to SIBYLL internal
C
C     input:     Npdg        PDG particle number
C     output:    sib_pdg2pid internal particle id
C
C     based on similar phojet function \FR'13
C----------------------------------------------------------------
      SAVE
      COMMON /S_PDG2PID/ IPID_PDG_LIST(99),ID_LIST(577)
      COMMON /S_DEBUG/ Ncall, Ndebug, lun
      DIMENSION LBAR(99)
      DATA LBAR /1,3,2,5,4,6,8,7,10,9,11,12,-13,-14,16,15,18,17,13,14,
     +  22,21,23,24,26,25,27,29,28,31,30,32,33,-34,-35,-36,-37,-38,-39,
     +  -40,-41,-42,-43,-44,-45,-46,-47,-48,-49,9*0,60,59,10*0,72,71,
     +  73,75,74,77,76,79,78,81,80,0,83,-84,-85,-86,-87,-88,-89,4*0,-94,
     +  -95,-96,-97,-98,-99 /

      Nin = abs(Npdg)
      if((Nin.gt.99999).or.(Nin.eq.0)) then
C  invalid particle number
        if(ndebug.gt.5) write(lun,'(1x,A,I10)')
     &    'isib_pdg2pid: invalid PDG ID number ',Npdg
        isib_pdg2pid = 0
        return
      else If(Nin.le.577) then
C  simple case
        Nout = Nin
      else
C  use hash algorithm
        Nout = mod(Nin,577)
      endif
 100  continue

C  particle not in table
      if(ID_list(Nout).Eq.0) then
        if(ndebug.ge.0) write(lun,'(1x,A,I10)')
     &    'isib_pdg2pid: particle not in table ',Npdg
        isib_pdg2pid = 0
        return
      endif

      if(IPID_pdg_list(ID_list(Nout)).eq.Nin) then
C  particle ID found
        isib_pdg2pid = ID_list(Nout)
        if (NPDG.lt.0) isib_pdg2pid = lbar( isib_pdg2pid )
        return
      else
C  increment and try again
        Nout = Nout + 5
        If(Nout.gt.577) Nout = Mod(Nout,577)
        goto 100
      endif

      END
            
      FUNCTION GASDEV(Idum)
C...Gaussian deviation
      SAVE
      SAVE GSET
      DATA ISET/0/
      IF (ISET.EQ.0) THEN
1       V1=2.*S_RNDM(0)-1.
        V2=2.*S_RNDM(0)-1.
        R=V1**2+V2**2
        IF(R.GE.1.) GO TO 1
        FAC=SQRT(-2.*LOG(R)/R)
        GSET=V1*FAC
        GASDEV=V2*FAC
        ISET=1
      ELSE
        GASDEV=GSET
        ISET=0
      ENDIF
      RETURN
      END