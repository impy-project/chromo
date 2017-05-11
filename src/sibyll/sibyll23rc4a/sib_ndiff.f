
      SUBROUTINE SIB_NDIFF(K_beam, NW, Ecm, Irec, IREJ)
C-----------------------------------------------------------------------
C     routine that samples and fragments a non-diffractive interaction
C
C     3 stages: 0: setup
C               1: sampling of event structure (number of parton interactions)
C                  (labeled as 2000)
C               2: sampling of kinematics
C                  (labeled as 3000)
C               3: fragmentation
C-----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE

c     external types
      DOUBLE PRECISION ECM
      INTEGER K_beam, NW, Irec, IREJ

c     COMMONs
      INCLUDE 'sib_debug_cmmn.inc'
      include 'sib_nw_prm.inc'
      INCLUDE 'sib_run_cmmn.inc'
      include 'sib_plist_cmmn.inc'
      include 'sib_parto_cmmn.inc'
      include 'sib_chist_cmmn.inc'
      include 'sib_int_prm.inc'
      INCLUDE 'sib_indx_cmmn.inc'
      INCLUDE 'sib_mass_cmmn.inc'
      INCLUDE 'sib_cnt_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'

c     internal type declarations
      DOUBLE PRECISION X2JET,SQS_0,PZ,E2,PAWT,xnsof,xnjet,xjdif,x1jet,
     &     Esum,PXsum,PYsum,PZsum
      DIMENSION X2JET(NW_max)
      INTEGER LL,LXBAD,NP_0,NPP_0,NPP0_0,J,JJ,I,KBA,L,NPP_1,NPP0_1,
     &     IREFout,IREF,nj,ns,nv,II,Idm,LPID,NF,NPP,NPP0
      DIMENSION LL(39)      
      DATA LL /5*0,7*2,2*1,19*0,6*1/

C..   setup stage
      IREJ = 1
c     default return point is kinematic sampling stage
      LXBAD = 3

c     remember initial setup
      NP_0    = NP
      SQS_0   = SQS
c     remember position on parton stack
      call GET_NPP(NPP_0,NPP0_0)

c     set interaction properties
c      IF(Irec.ne.1) CALL INI_EVENT(ECM,K_beam,Idm,Irec)

      IF(ndebug.gt.0)then
         IF(Irec.eq.0)THEN
            WRITE(LUN,*) 
     &           'SIB_NDIFF: recursive call with (ecm,kb,kt,np,jdif):',
     &           ecm,k_beam,kt(1),(jdif(j),j=1,NW),NP
         ELSE
            WRITE(LUN,*)' SIB_NDIFF: regular call with (ECM,KB,NW,KT,',
     &           'JDIF,NP):',ecm,k_beam,NW,(kt(ii),ii=1,NW),
     &           (jdif(j),j=1,NW),NP
         ENDIF
      ENDIF
      
 2000 CONTINUE

c     reset parton stack
      call ini_prtn_stck(NPP_0,NPP0_0)

C...  sample multiple interaction configuration
      KBA = IABS(K_beam)
      L = LL(KBA)
      DO I=1,NW
        if(JDIF(I).eq.0) then
           CALL CUT_PRO(L, SQS, PTmin, NNSOF(I), NNJET(I))
        else
          NNSOF(I) = 1
          NNJET(I) = 0
        endif
c     add incoming target particles
        PZ = PAWT(SQS,AM(KBA),AM(KT(I)))
        E2 = SQRT(PZ**2+AM2(KT(I)))
        call add_prtn(ZERO,ZERO,-PZ,E2,AM(KT(I)),KT(I),-2,0,IREFout)

c     add interactions
        xjdif = dble(jdif(I))
        xnjet = dble(nnjet(I))
        xnsof = dble(nnsof(I))
        call add_prtn(xnsof,xnjet,xjdif,sqs,zero,I,-1,IREFout,IREF)
c     write parton stack index to interaction index
        IINTDX(I) = IREF
      ENDDO
c     remember state of parton stack
      call GET_NPP(NPP_1,NPP0_1)

C...  kinematic sampling stage

C...  sample x values
      ITRY(1) = 0
 3000 CONTINUE
      ITRY(1) = ITRY(1)+1
      IF(ITRY(1).GT.NREJ(1)) THEN 
c         NCALL = NCALL + 1
         GOTO 2000
      ENDIF
      NP = NP_0
      call ini_prtn_stck(NPP_1,NPP0_1)

      call sample_minijet(L,NW,NNJET,NNSOF,NJET,NSOF,x1jet,x2jet,lxbad)
      SELECT CASE(LXBAD)
      CASE (3)
c     reject kinematics
         GOTO 3000
      CASE (2) 
c     reject kinematics and event structure
c         NCALL = NCALL + 1
         GOTO 2000
      CASE (1)
c     reject entire event
         if(Ndebug.gt.0) 
     &        WRITE(LUN,*)' SIB_NDIFF: minijet rejection (Ncall):',Ncall
c     restore initial state
         NP    = NP_0
         call ini_prtn_stck(NPP_0,NPP0_0)
         SQS   = SQS_0
         S     = SQS*SQS
         RETURN
      END SELECT


C...  Prepare 2*NW valence/sea color strings and/or remnant.

c     default return point, jump back to sampling interaction structure
c      LXBAD = 2
      call SAMPLE_rmnt(K_beam,NW,X1Jet,X2JET,Irec,LXBAD)
      SELECT CASE(LXBAD)
      CASE (3)
c     reject kinematics
         GOTO 3000
      CASE (2) 
c     reject kinematics and event structure
c         NCALL = NCALL + 1
         GOTO 2000
      CASE (1)
c     reject entire event
         if(Ndebug.gt.0) 
     &   WRITE(LUN,*)' SIB_NDIFF: rmnt rejection (Ncall,NW):',Ncall,NW
c     restore initial state
         NP    = NP_0
         call ini_prtn_stck(NPP_0,NPP0_0)
         SQS   = SQS_0
         S     = SQS*SQS
         RETURN
      END SELECT

C     Check parton final state..
      call GET_NPP(NPP,NPP0)
      CALL PPsum(1,NPP,Esum,PXsum,PYsum,PZsum,NF)
      IF(ABS(Esum/(HALF*Ecm*DBLE(NW+1))-ONE).GT.EPS3)THEN
         WRITE(LUN,*) ' SIB_NDIFF: energy not conserved! : ',Ncall
         WRITE(LUN,*) ' sqs_inp = ', Ecm, ' sqs_out = ', Esum
         call prnt_prtn_stck
         WRITE(LUN,*) ' SIB_NDIFF: event rejected! ',
     &        'partons do not conserve energy'
         WRITE(LUN,*)' (Ncall,NW,NPP,NJET,NSOF):',Ncall,NW,NPP,NJET,NSOF
         call sib_reject
      ENDIF
      IF(NDEBUG.gt.0) THEN
         IF(NDEBUG.gt.1) call prnt_prtn_stck
         WRITE(LUN,*) ' SIB_NDIFF: entering fragmentation stage...'
      ENDIF

C...  Fragmentation stage
      nj = 0
      ns = 0
      nv = 0
      II = NPP0_0+1
      DO WHILE (II.gt.0)
c     default return point: reject event if fragmentation fails
         LXBAD = 1         
c     loop over level0 partons
         call ITR_LVL0_PRTN(II,JJ,LPID)
c     read interaction
         call rd_int(jj,Idm,iiflag)
         SELECT CASE(LPID)

C...  Fragmentation of soft/hard sea color strings
         CASE(100)
            nj = nj + 1
            ipflag = 100
            Nint = nj
            call frag_minijet(jj,LXBAD)
            IF(LXBAD.ne.0) RETURN

         CASE(10)
            ns = ns + 1
            ipflag = 10
            Nint = ns
            call frag_minijet(jj,LXBAD)
            IF(LXBAD.ne.0) RETURN

C...  fragment 'valence' strings
         CASE(1)
            nv = nv + 1
            Nint = nv
            ipflag = 1
            call frag_vlnce(jj,LXBAD)
            IF(LXBAD.ne.0) RETURN

C...  fragment remnants
         CASE(2,-2)
            CALL EXCTDEC(JJ,LXBAD)
            IF(LXBAD.ne.0) RETURN

C...  fragment incoherent diffraction
         CASE(-10,-20,-30)
            CALL frag_inchrnt_diff(jj,lxbad)
            IF(LXBAD.ne.0) RETURN

         END SELECT
      ENDDO
      IREJ = 0
      
      END
