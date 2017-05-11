      SUBROUTINE EXCT_RMNT(JW,KRMNT,IREJ)
C-----------------------------------------------------------------------
C     routine to produce massive excitations of beam and/or target \FR'14
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_mass_cmmn.inc'
      include 'sib_nw_prm.inc'
      include 'sib_int_prm.inc'
      INCLUDE 'sib_run_cmmn.inc'
      INCLUDE 'sib_indx_cmmn.inc'
      INCLUDE 'sib_rmnt_cmmn.inc'
      INCLUDE 'sib_chp_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      INCLUDE 'sib_chist_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'
      INCLUDE 'sib_cnt_cmmn.inc'
      INCLUDE 'sib_difmass_cmmn.inc'

      DIMENSION P1(5),P2(5),P1n(5),P2n(5),PBM1(5),PBM2(5),PBM(5),
     &     PTG1(5),PTG2(5),PTG(5),PTT(5),GABE(4)

      DIMENSION LL(39)    
      DATA LL /5*0,7*2,2*1,19*0,6*1/

      
c     default return point, beam and target sampling
c      IREJ = 1

      IF(NDEBUG.gt.2)
     &     WRITE(LUN,*) 'EXCT_RMNT: input (JW,KRMNT,IREJ)', 
     &     JW,KRMNT,IREJ

      IF(NDEBUG.gt.3)then
         write(LUN,*) 'beam remnant index: (lvl0,flv1,flv2)' , IBMRDX
         write(LUN,*) '1st central string index: (lvl0,bm,tg)' , 
     &        (ICSTDX(2*(JW-1)+1,ii),ii=1,3)
         write(LUN,*) '2nd central string index: (lvl0,bm,tg)' ,
     &        (ICSTDX(2*(JW-1)+2,ii),ii=1,3)

         write(LUN,*) 'target remnant index: (lvl0,flv1,flv2)' ,
     &        (ITGRDX(JW,ii),ii=1,3)
      ENDIF

      ITRY(5) = 0

C...  select indices depending on configuration
C     krmnt = 0 : no excitation on either side
c           = 1 : beam side excited remnant
c           = 2 : target side
c           = 3 : both sides

c     write remnant configuration to remnant common
      IRMNT(JW) = KRMNT
      SELECT CASE(KRMNT)
      CASE(1)
c     beam side remnant only
c     proto-remnant position: IBMRDX(1)
c     partons in : IBMRDX(2:3)
         IBM1 = IBMRDX(2)
         IBM2 = IBMRDX(3)
c     target side to transfer energy from: 
c     (sofar always choose valence pair)
         ITG1 = ICSTDX(2*(JW-1)+1,3)
         ITG2 = ICSTDX(2*(JW-1)+2,3)
c     beam-side partons to go into central strings
         IBMST1 = ICSTDX(2*(JW-1)+1,2)
         IBMST2 = ICSTDX(2*(JW-1)+2,2)
c     target-side partons to go into central strings
         ITGST1 = ITG1
         ITGST2 = ITG2

      CASE(2)
c     target side remnant only
c     proto-remnant in ITGRDX(JW,1)
         ITG1 = ITGRDX(JW,2)
         ITG2 = ITGRDX(JW,3)
c     transfer energy from beam remnant or 
c     central strings with valence quarks
c     in ICSTDX(JW+0:1,2)
c     means no beam remnant --> get from valence strings
         IBM1 = ICSTDX(2*(JW-1)+1,2)
         IBM2 = ICSTDX(2*(JW-1)+2,2)
c     beam-side partons to go into central strings
         IBMST1 = IBM1
         IBMST2 = IBM2
c     target-side partons to go into central strings
         ITGST1 = ICSTDX(2*(JW-1)+1,3)
         ITGST2 = ICSTDX(2*(JW-1)+2,3)

      CASE(3)
c     beam and target side remnant
c     transfer energy from pairs in rmnt or central strings
c     listed in I?RDX and ICSTDX()
         IBM1 = IBMRDX(2)
         IBM2 = IBMRDX(3)
         ITG1 = ITGRDX(JW,2)
         ITG2 = ITGRDX(JW,3)

      CASE(0)
c     no excited remnant case, jump straight to central strings..
         GOTO 100

      END SELECT

      IF(NDEBUG.gt.3)then
         write(lun,*) 'beam parton1:',IBM1
         write(lun,*) 'beam parton2:',IBM2
         write(lun,*) 'target parton1:',ITG1
         write(lun,*) 'target parton2:',ITG2
      endif

c     save status of parton stack
      call get_npp(npld,np0ld)

 10   ITRY(5) = ITRY(5) + 1
      IF(ITRY(5).GT.NREJ(5))THEN
         IF(NDEBUG.gt.2) 
     &        WRITE(LUN,*) '  EXCT_RMNT: no. of trials exceeded, ',
     &        NREJ(5), 'resample minijets ...' , IREJ
         RETURN
      ENDIF
c     reset parton stack after rmnt mass rejection
      call ini_prtn_stck(npld,np0ld)

C..   construct 4momenta of proto-remnants
c     index of beam remnant on stack: IBMRDX(1)

C..   center-of-mass energy of parton system (s hat)
c     calculated in hadron-hadron frame
c     for first interaction (jw=1) partons are massless and collinear (sum pt=0)
c     in this case ecm = SQS*SQRT(XB*XT), xb,t=x1+x2 
c      for jw>1 beam partons may have already acquired mass and additional pt
c     therefore ecm = sqs*sqrt(xb*xt) + corr.
c     IRDX: index of remnant on parton stack
c      SHAT = S*XB*XT+XM2+(XT/XB)*XMT2

c     with 4momenta of partons on stack, momentum fractions are obsolete
c     center-of-mass energy is simply: shat = (pbm+ptg)**2

c     construct total 4momentum
c     add beam-side parton momenta, in had.-had. frame
      call rd_prtn_4vec(IBM1,PBM1,IFL,Idm)
      call rd_prtn_4vec(IBM2,PBM2,IFL,IBMH)
      call add_4vecs(PBM1,PBM2,PBM)
      
c     target-side parton momenta, in had.-had. frame
      call rd_prtn_4vec(ITG1,PTG1,IFL,Idm)
      call rd_prtn_4vec(ITG2,PTG2,IFL,Idm)
      call add_4vecs(PTG1,PTG2,PTG)
      
c     add beam and target side to get total 4momentum
      call add_4vecs(PBM,PTG,PTT)
      shat = PTT(5)**2
      ecm = PTT(5)
c     catch virtual remnants
      IF(PTT(5).LT.ZERO) THEN
         IF(NDEBUG.GT.2)THEN
            WRITE(LUN,*) ' EXCT_RMNT: too little mass left (Shat):',
     &           shat
            WRITE(LUN,*) '        resample minijets...'
         ENDIF
         LREJ = 2
         RETURN                 ! resample minijets
      ENDIF


      IF(NDEBUG.GT.2) WRITE(LUN,*) ' EXCT_RMNT: try no.',ITRY(5)
      IF(NDEBUG.GT.3)THEN
         write(LUN,*) '4momenta before scattering:'
         write(LUN,*) ' PBM1:' , (PBM1(jj),jj=1,5)
         write(LUN,*) ' PBM2:' , (PBM2(jj),jj=1,5)
         write(LUN,*) ' PBM:' , (PBM(jj),jj=1,5)

         write(LUN,*) ' PTG1:' , (PTG1(jj),jj=1,5)
         write(LUN,*) ' PTG2:' , (PTG2(jj),jj=1,5)
         write(LUN,*) ' PTG:' , (PTG(jj),jj=1,5)

         write(LUN,*) ' PTT:' , (PTT(jj),jj=1,5)
      ENDIF

      IF(NDEBUG.gt.2)
     &     WRITE(LUN,*)' EXCT_RMNT: Shat:',shat

      XMFRAC = PAR(81)
      XSFRAC = PAR(82)

c     exponent of remnant mass distribution (1/Mx**2)**alpha
c     by default: alpha = 1
c     different for baryons and mesons
c      ALPHA = PAR(98)

C..   Sample masses
      SELECT CASE(KRMNT)
      CASE(1)
         XM2MAX = MIN(XSFRAC*S,XMFRAC*AM2(IABS(KB)))
         XM2MAX = MAX(XM2MAX,ONE)

c     mass of target-side: 0
         XMT = ZERO
         XMT2 = ZERO
c     get remnant mass
c     (might have received mass from prior interaction)
         call get_mass2(IBMRDX(1),XM2)
c     allowing excitation to fallback to beam means min.
c     mass is beam mass, or more exact smallest mass of hadrons 
c     with flavors in remnant
         IF(IPAR(64).eq.1)THEN
c     remnant mass can also decrease through interactions
            Xmsqmin = AM2(IABS(KB))
         ELSE
c     remnant mass only increased by multiple interactions..
            Xmsqmin = MAX(AM2(IABS(KB)),XM2)
         ENDIF
C     select exponent from COMMON
         ALPHA = XRMEX(LL(IABS(KB)))
c     sample beam mass
         XMB2 = XM2DIS(XMSQMIN,XM2MAX,ALPHA)
         IF(NDEBUG.gt.2)
     &        WRITE(LUN,*)' EXCT_RMNT: XM2min,XM2max,ALPHA,XM2:',
     &        Xmsqmin,XM2max,ALPHA,XMB2
c     check if resonance or massive hadron has to be formed
         call sel_res(XMB2,KRB,IBMRDX(1),IBMH)
         XMB = sqrt(XMB2)


      CASE(2)
c     target side mass
         XM2MAX = MIN(XSFRAC*S,XMFRAC*AM2(IABS(KT(JW))))
         XM2MAX = MAX(XM2MAX,ONE)

         XMB = ZERO
         XMB2 = ZERO
         Xmsqmin = AM2(KT(JW))
C     select exponent from COMMON
         ALPHA = XRMEX(LL(IABS(KT(JW))))
         XMT2 = XM2DIS(XMSQMIN,XM2MAX,ALPHA)
         IF(NDEBUG.gt.2)
     &        WRITE(LUN,*)' EXCT_RMNT: XM2min,XM2max,ALPHA,XM2:',
     &        Xmsqmin,XM2max,ALPHA,XMT2

c     check if resonance or massive hadron has to be formed
         call sel_res(XMT2,KRT(JW),ITGRDX(JW,1),ITGH)
         XMT = sqrt(XMT2)


      CASE(3)
         XM2MAX = MIN(XSFRAC*S,XMFRAC*AM2(IABS(KB)))
         XM2MAX = MAX(XM2MAX,ONE)

         call get_mass2(IBMRDX(1),xm2)
         IF(IPAR(64).eq.1)THEN
c     remnant mass can also decrease through interactions
            Xmsqmin = AM2(IABS(KB))
         ELSE
c     remnant mass only increased by multiple interactions..
            Xmsqmin = MAX(AM2(IABS(KB)),XM2)
         ENDIF
C     select exponent from COMMON
         ALPHA = XRMEX(LL(IABS(KB)))        
         XMB2 = XM2DIS(XMSQMIN,XM2MAX,ALPHA)
         IF(NDEBUG.gt.2)
     &        WRITE(LUN,*)' EXCT_RMNT: XM2min,XM2max,ALPHA,XM2:',
     &        Xmsqmin,XM2max,ALPHA,XMB2

c     check if resonance or massive hadron has to be formed
         call sel_res(XMB2,KRB,IBMRDX(1),IBMH)
         XMB = SQRT(XMB2)
         
c     target always nucleon
         XM2MAX = MIN(XSFRAC*S,XMFRAC*AM2(IABS(KT(JW))))
         XM2MAX = MAX(XM2MAX,ONE)

         Xmsqmin = AM2(IABS(KT(JW)))
C     select exponent from COMMON
         ALPHA = XRMEX(LL(IABS(KT(JW))))        
         XMT2 = XM2DIS(XMSQMIN,XM2MAX,ALPHA)
         IF(NDEBUG.gt.2)
     &        WRITE(LUN,*)' EXCT_RMNT: XM2min,XM2max,ALPHA,XM2:',
     &        Xmsqmin,XM2max,ALPHA,XMT2

c     check if resonance or massive hadron has to be formed
         call sel_res(XMT2,KRT(JW),ITGRDX(JW,1),ITGH)
         XMT = SQRT(XMT2)

      END SELECT
c     write excitation mass to output common
      XRMASS(1) = XMB
      XRMASS(2) = XMT

c     minimal mass requirement
c      IF(SHAT.lt.XMB2+XMT2+0.3) GOTO 10
      IF(SHAT.lt.XMB2+XMT2+TWO*XMB*XMT+0.3D0) GOTO 10

C     transfer cm energy to mass of particle in parton-parton cm
      call transfOnShell(ECM,XMB,XMT,XM2MAX,1,P1,P2,IBD)
      IF(IBD.eq.1) THEN
         IF(NDEBUG.gt.2) WRITE(LUN,*) ' EXCT_RMNT: excitation rejected!'
         RETURN
      ENDIF

C...  Boost 4momenta to hadron-hadron center-of-mass
c     along z only if initial partons do not carry transverse momentum
c     (cancels between val1 and val2)
c     with multiple nucleons interacting beam val partons can aquire 
c     transverse momentum from the target. in this case need arbitrary boost
      do k = 1,4
         gabe(k) = PTT(k)/PTT(5)
      enddo
      call SIB_ALTRA(gabe(4), gabe(1), gabe(2), gabe(3),
     &     P1(1),p1(2),p1(3),p1(4),
     &     P1TOT,p1n(1),p1n(2),p1n(3),p1n(4))
      p1n(5)=p1(5)
      call SIB_ALTRA(gabe(4), gabe(1), gabe(2), gabe(3),
     &     P2(1),p2(2),p2(3),p2(4),
     &     P2TOT,p2n(1),p2n(2),p2n(3),p2n(4))
      p2n(5)=p2(5)

C...  Calculate new 4momentum of partons in had.-had. frame
c     P1,P2: momenta after scattering in parton-parton cm.
c     P1n,P2n: momenta after scattering in had.-had. cm
c     PBM1,2: momenta of beam partons in had.-had. before scattering
c     PTG1,2: momenta of target partons in had.-had. before scattering
c     PBM: combined momentum of all beam partons before scattering
c     PTG: combined momentum of all target partons before scattering

c     energy and z component
      do ii=3,4
         PBM1(ii) = PBM1(ii)*P1n(ii)/PBM(ii)
         PBM2(ii) = PBM2(ii)*P1n(ii)/PBM(ii)

         PTG1(ii) = PTG1(ii)*abs(P2n(ii)/PTG(ii))
         PTG2(ii) = PTG2(ii)*abs(P2n(ii)/PTG(ii))
      enddo

c     if transverse momentum prior to interaction zero then
c     assign transverse momentum of partons according to random fraction
      IF(ABS(PBM(1)).LT.EPS10.or.ABS(PBM(2)).LT.EPS10)THEN
         do ii = 1,2
            XI = S_RNDM(JW)
            PBM1(ii) = XI*P1N(ii)
            PBM2(ii) = (ONE-XI)*P1N(ii)
         enddo
      ELSE
         do ii=1,2
            PBM1(ii) = PBM1(ii)*P1n(ii)/PBM(ii)
            PBM2(ii) = PBM2(ii)*P1n(ii)/PBM(ii)
         enddo         
      ENDIF

      IF(ABS(PTG(1)).LT.EPS10.or.ABS(PTG(2)).LT.EPS10)THEN
         do ii=1,2
            XI = S_RNDM(JW)
            PTG1(ii) = XI*P2N(ii)
            PTG2(ii) = (ONE-XI)*P2N(ii)
         enddo
      ELSE
         do ii=1,2
            PTG1(ii) = PTG1(ii)*P2n(ii)/PTG(ii)
            PTG2(ii) = PTG2(ii)*P2n(ii)/PTG(ii)
         enddo                  
      endif

      IF(NDEBUG.GT.3)THEN
         write(LUN,*) 'parton 4momenta after scattering:'
         write(LUN,*) ' PBM1:' , (PBM1(jj),jj=1,5)
         write(LUN,*) ' PBM2:' , (PBM2(jj),jj=1,5)
         write(LUN,*) ' sum:' , (PBM2(jj)+PBM1(jj),jj=1,5)
         write(LUN,*) ' PTG1:' , (PTG1(jj),jj=1,5)
         write(LUN,*) ' PTG2:' , (PTG2(jj),jj=1,5)
         write(LUN,*) ' sum:' , (PTG2(jj)+PTG1(jj),jj=1,5)
      ENDIF
      
C...  change parton 4momenta on stack
      call edt_prtn(IBM1,PBM1(1),PBM1(2),PBM1(3),PBM1(4),PBM1(5),Idm)
      call edt_prtn(IBM2,PBM2(1),PBM2(2),PBM2(3),PBM2(4),PBM2(5),Idm)

      call edt_prtn(ITG1,PTG1(1),PTG1(2),PTG1(3),PTG1(4),PTG1(5),Idm)
      call edt_prtn(ITG2,PTG2(1),PTG2(2),PTG2(3),PTG2(4),PTG2(5),Idm)
         
C...  add remnants
c     references are circular: 
c     rmnt --> parton1 --> parton2 --> lvl2 rmnt (hadron) --> rmnt
      SELECT CASE(KRMNT)
      CASE(1)
c     beam side remnant, add only if does not exist yet otherwise edit
         IF(IBMRDX(1).eq.0)THEN
            call add_prtn
     &           (P1n(1),P1n(2),P1n(3),P1n(4),P1n(5),2,0,IBM1,IBMRDX(1))
         ELSE
            call edt_prtn
     &           (IBMRDX(1),P1n(1),P1n(2),P1n(3),P1n(4),P1n(5),Iref)
         ENDIF
c     add beam hadron as hypothetical final state
         IF(IBMH.eq.0)THEN
            call add_prtn
     &         (P1n(1),P1n(2),P1n(3),P1n(4),P1n(5),KRB,2,IBMRDX(1),IBMH)
         else
            call edt_prtn
     &           (IBMH,P1n(1),P1n(2),P1n(3),P1n(4),P1n(5),Iref)
         ENDIF
c     add references rmnt --> parton1 etc
         call add_ref(IBMRDX(1),IBM1)
         call add_ref(IBM1,IBM2)
         call add_ref(IBM2,IBMH)

      CASE(2)
c     add target side remnant
         IF(ITGRDX(JW,1).eq.0)THEN
            call add_prtn
     &           (P2n(1),P2n(2),P2n(3),P2n(4),P2n(5),
     &           -2,0,0,ITGRDX(JW,1))
         else
            call edt_prtn
     &           (ITGRDX(JW,1),P2n(1),P2n(2),P2n(3),P2n(4),P2n(5),Iref)
         endif
         IF(ITGH.eq.0)then
c     add target hadron as hypothetical final state, always nucleon
            call add_prtn
     &           (P2n(1),P2n(2),P2n(3),P2n(4),P2n(5),
     &           KRT(JW),2,ITGRDX(JW,1),ITGH)
         else
            call edt_prtn
     &           (ITGH,P2n(1),P2n(2),P2n(3),P2n(4),P2n(5),Iref)
         endif

c     add references rmnt --> parton1 etc
         call add_ref(ITGRDX(JW,1),ITG1)
         call add_ref(ITG1,ITG2)
         call add_ref(ITG2,ITGH)

      CASE(3)
c     beam side remnant, add only if does not exist yet, otherwise edit
         IF(IBMRDX(1).eq.0)THEN
            call add_prtn
     &           (P1n(1),P1n(2),P1n(3),P1n(4),P1n(5),2,0,0,IBMRDX(1))
         ELSE
            call edt_prtn
     &           (IBMRDX(1),P1n(1),P1n(2),P1n(3),P1n(4),P1n(5),Iref)
         ENDIF
c     add beam hadron as hypothetical final state
         IF(IBMH.eq.0)then
            call add_prtn
     &         (P1n(1),P1n(2),P1n(3),P1n(4),P1n(5),KRB,2,IBMRDX(1),IBMH)
         else
            call edt_prtn
     &           (IBMH,P1n(1),P1n(2),P1n(3),P1n(4),P1n(5),Iref)
         endif
         call add_ref(IBMRDX(1),IBM1)
         call add_ref(IBM1,IBM2)
         call add_ref(IBM2,IBMH)

c     add target side remnant
         IF(ITGRDX(JW,1).eq.0)THEN
            call add_prtn
     &           (P2n(1),P2n(2),P2n(3),P2n(4),P2n(5),-2,0,0,Iref)
            ITGRDX(JW,1) = Iref
         else
            call edt_prtn
     &           (ITGRDX(JW,1),P2n(1),P2n(2),P2n(3),P2n(4),P2n(5),Iref)
         endif
         IF(ITGH.eq.0)then
c     add target hadron as hypothetical final state
            call add_prtn
     &           (P2n(1),P2n(2),P2n(3),P2n(4),P2n(5),
     &           KRT(JW),2,ITGRDX(JW,1),ITGH)
         else
            call edt_prtn
     &           (ITGH,P2n(1),P2n(2),P2n(3),P2n(4),P2n(5),Iref)
         endif
c     add references rmnt --> parton1 etc
         call add_ref(ITGRDX(JW,1),ITG1)
         call add_ref(ITG1,ITG2)
         call add_ref(ITG2,ITGH)

      END SELECT

 100  IF(JDIF(JW).ne.0.and.NWD.ne.1)THEN
c     incoherent diffraction case
c     add parton 4momenta to obtain c.m energy
         
c     beam side
         IBMST1 = ICSTDX(2*(JW-1)+1,2)
         IBMST2 = ICSTDX(2*(JW-1)+2,2)

c     target side
         ITGST1 = ICSTDX(2*(JW-1)+1,3)
         ITGST2 = ICSTDX(2*(JW-1)+2,3)
         
         call rd_prtn_4vec(IBMST1,PBM1,IFLB1,Idm)
         call rd_prtn_4vec(IBMST2,PBM2,IFLB2,Idm)
         call add_4vecs(PBM1,PBM2,PBM)
         call rd_prtn_4vec(ITGST1,PTG1,IFLT1,Idm)
         call rd_prtn_4vec(ITGST2,PTG2,IFLT2,Idm)
         call add_4vecs(PTG1,PTG2,PTG)
c     total 4momentum
         call add_4vecs(PBM,PTG,PTT)
c     add diffractive system to parton stack
c     references are: diff --> diff. hadron 
c     --> beam parton1 --> beam parton2 --> target parton1 etc
         call add_prtn_4vec(PTT,-10*JDIF(JW),0,IBMST1,Iref)
         call add_int_ref(Iref,IINTDX(JW))
c     both string indices point to diff. system
         ICSTDX(2*(JW-1)+1,1) = Iref
         ICSTDX(2*(JW-1)+2,1) = Iref
c     add diff. beam hadron to stack
c     model assumes remnant always excited in first interaction
         L0 = KB
c     if not first interaction or remnant excited, merge sea pair to hadron
         if(krmnt.ne.0.or.jw.ne.1) then       
            L0 = IMRG2HAD(IFLB1,IFLB2)
c     call sib_i4flav(IFLB1,IFLB2,Idm,Idm1,L0)
         endif
c     check kinematic limits
c     m2_max should be smaller than m2_min
         IREJ = 1
         EE = PTT(5)
         EE2 = PTT(5)**2
         K = 2-IBAR(IABS(L0))
         IF(JDIF(jw).gt.1)THEN
            deltae = ee-am(13)
            XMMIN=max(XM2MIN(1),(am(iabs(l0))+am(7)+0.02D0)**2)
         else
            deltae = EE-AM(IABS(L0))
            xmmin=max(XM2MIN(K),(am(iabs(l0))+am(7)+0.02D0)**2)
         endif
c         print *,'jw,jdif,nwd,l0,ifl1,ifl2,deltae,xmin,ee,xmax',
c     &        jw,jdif(jw),nwd,l0,ifl1,ifl2,deltae,xmmin,ee,par(13)*ee2
         IF(DELTAE.lt.AM(7)+0.02D0) THEN
            IF(ndebug.gt.2) 
     &           WRITE(lun,*) ' EXCT_RMNT: inchoherent diff. :',
     &           ' not enough mass left for excitation! (DELTAE,PION,',
     &           'IREJ,NCALL)',DELTAE,AM(7)+0.02D0,IREJ,NCALL
            RETURN
         ENDIF
         IF(PAR(13)*EE2.lt.XMMIN)THEN
            IF(ndebug.gt.2)
     &           WRITE(lun,*) ' EXCT_RMNT: inchoherent diff. :',
     &           ' not enough mass left for excitation! (min,max,',
     &           'IREJ,NCALL)',PAR(13)*EE2,XMMIN,IREJ,NCALL
            RETURN
         ENDIF
         call add_prtn_4vec(PTT,L0,2,IBMST1,Idhad)
         call add_ref(Iref,Idhad)
c     reset references of partons
         call add_ref(IBMST1,IBMST2)
         call add_ref(IBMST2,ITGST1)
         call add_ref(ITGST1,ITGST2)
         call add_ref(ITGST2,Iref)
         IF(ndebug.gt.2) THEN
            WRITE(LUN,*) ' EXCT_RMNT: incoherent diff. ',
     &           '(IDX,IDX2,JDIF,ECM,L0)',Iref,IDhad,JDIF(jw),PTT(5),L0
            WRITE(LUN,*) ' EXCT_RMNT: DELTAE,XM2MAX:',DELTAE,PAR(13)*EE2
         ENDIF
         IREJ = 0
         RETURN
      ENDIF

C...  add central strings to stack
c     partons designated for central strings 
c     are indexed in ICSTDX(JW,2:3)
c     pstr_j = p_j_bm + p_j_tg
c     string mass ** 2 = pstr_j ** 2
c     --> read momenta from stack, add beam and target side, 
c     references are set in a loop:
c     string --> beam-parton --> target-parton --> string
c     then write string 4momentum on stack
      IMRG = 0
      DO JJ=1,2
         ISTH = 0
         IBMST = ICSTDX(2*(JW-1)+JJ,2)
         ITGST = ICSTDX(2*(JW-1)+JJ,3)
         call rd_prtn_4vec(IBMST,PBM1,IFL1,Idm)
         call rd_prtn_4vec(ITGST,PTG1,IFL2,Idm)
         call add_4vecs(PBM1,PTG1,PTT)
c     transverse mass of string end partons (pt**2)
         call get_xmt2(IBMST,XMT12)
         call get_xmt2(ITGST,XMT22)
c     available mass for string
         EE = sqrt(PTT(4)**2-PTT(3)**2)
c     catch virtual strings
         if(PTT(5).lt.ZERO) then
            IREJ = 1
            IF(ndebug.gt.2)
     &           write(LUN,*)' EXCT_RMNT: virt. string (M):',EE
            IF(ndebug.gt.3)then
               call get_imass2(IBMST,xm2)
               write(LUN,*) 'PBM1:', (PBM1(j),j=1,5),xm2
               call get_imass2(ITGST,xm2)
               write(LUN,*) 'PTG1:', (Ptg1(j),j=1,5),xm2
               write(LUN,*) 'Ptot:', (Ptt(j),j=1,5)
            ENDIF
c               stop
            RETURN
         ENDIF
c     minimal string mass requierement
         if(EE.lt.sqrt(xmt12)+sqrt(xmt22)+PAR(123))then
            IF(IPAR(74).eq.1)THEN
c     try to form single meson, set merge flag
               IF(IABS(IFL1).gt.10.and.IABS(IFL2).gt.10) THEN
c     skip if two diquarks need merging..                  
                  IREJ = 1
                  RETURN
               ENDIF
               L0 = IMRG2HAD(IFL1,IFL2)
               IF(EE.gt.AM(IABS(L0))) then
                  IMRG = IMRG + JJ
                  call add_prtn_4vec(PTT,L0,2,IBMST,ISTH)
                  IF(ndebug.gt.2)then
                     write(lun,*)
     &                    ' EXCT_RMNT: c.string mass too low! ',
     &                    'merge into hadron..',l0
                  ENDIF
               ENDIF
            ELSE
               IF(ndebug.gt.2)then
                  write(lun,*)
     &                 ' EXCT_RMNT: c.string kinematic rejection!'
                  write(lun,*) ' EE,limit,XMT1,XMT2:',
     &                 EE,sqrt(xmt12)+sqrt(xmt22)+0.3D0,sqrt(xmt12),
     &                 sqrt(xmt22)
                  write(lun,*) ' return to momentum sampling..'
               endif
               IREJ = 1
               RETURN
            ENDIF
         ENDIF
c     add central string to stack, refering to beam-end parton
         call add_prtn_4vec(PTT,1,0,IBMST,Iref)
         ICSTDX(2*(JW-1)+JJ,1) = Iref
         call add_int_ref(Iref,IINTDX(JW))
c     add reference to target parton to beam parton
         call add_ref(IBMST,ITGST)
         IF(ISTH.ne.0) THEN
c     if string merged to hadron add reference corresponding reference            
            call add_ref(ITGST,ISTH)
            call add_ref(ISTH,IRef)
         ELSE
c     add reference to corresponding central string to target parton
            call add_ref(ITGST,Iref)
         ENDIF
      ENDDO
      
c     form single hadron from string if mass was too low ..
c     need to put hadron on shell by exchanging energy with other string            
      SELECT CASE(IMRG)
      CASE(1,2)
         IF(ndebug.gt.2)
     &        WRITE(LUN,*)' EXCT_RMNT: merging one string..',IMRG
c     one string merged
c     index of merged string and its last parton
         IMST = ICSTDX(2*(JW-1)+IMRG,1)
         IMST1 = ICSTDX(2*(JW-1)+IMRG,3)
c     index of ordinary string
         IMRGBAR = 3-IMRG
         ICST2 = ICSTDX(2*(JW-1)+IMRGBAR,1)
c     read 4momenta
         call rd_ref(IMST1,ISTH)
         call rd_prtn_4vec(ISTH,P1,L0,Iref)
c     string two
         call rd_prtn_4vec(ICST2,P2,IFL2,Idm)
c     cm energy
         call add_4vecs(P1,P2,PTT)
         IF(ndebug.gt.2)THEN
            write(lun,*)' EXCT_RMNT: string A :',(P1(i),i=1,5)
            write(lun,*)' EXCT_RMNT: string B :',(P2(i),i=1,5)
            write(lun,*)' EXCT_RMNT: total :',(PTT(i),i=1,5)
         ENDIF
         ecm = PTT(5)
         xm1 = AM(IABS(L0))
         xm2 = P2(5)
         call transfOnShell(ecm,xm1,xm2,one,3,P1n,P2n,LBD)
         IF(LBD.eq.1) THEN
            IF(NDEBUG.gt.2)
     &           WRITE(LUN,*)' EXCT_RMNT: mass transfer failed!'
            RETURN
         ENDIF
c     by definition p1n is along +z in string cm, need to invert if pzA < pzB
c         IF(P2(3).gt.P1(3)) call swtch_lmnts(P1N(3),P2N(3))

C..   rotate parton-parton axis onto string-string axis
c     therefore boost to parton-parton cm
c     to calc. rotation angles BEFORE interaction !
         do k = 1,4
            GABE(k) = PTT(k)/PTT(5)
         enddo         
         call SIB_ALTRA(gabe(4),-gabe(1),-gabe(2),-gabe(3),
     &        P1(1),p1(2),p1(3),p1(4),
     &        P1TOT,pbm1(1),pbm1(2),pbm1(3),pbm1(4))
c     rotation factors
         COD= PBM1(3)/P1TOT
         SID= DSQRT(PBM1(1)**2+PBM1(2)**2)/P1TOT
         COF=ONE
         SIF=ZERO
         IF(P1TOT*SID.GT.EPS5) THEN
            COF=PBM1(1)/(SID*P1TOT)
            SIF=PBM1(2)/(SID*P1TOT)
            ANORF=DSQRT(COF*COF+SIF*SIF)
            COF=COF/ANORF
            SIF=SIF/ANORF
         ENDIF
         IF(ndebug.gt.2)THEN
            write(lun,*)' EXCT_RMNT: momentum in cm:',(pbm1(i),i=1,5)
            write(lun,*)' EXCT_RMNT: rotation factors:',cod,sid,cof,sif
            write(lun,*)' EXCT_RMNT: rotation angles (theta,phi):',
     &           ACOS(cod),ACOS(cof),Asin(sid),Asin(sif)
            write(lun,*)' EXCT_RMNT: momentum:',
     &           sqrt(p1n(1)**2+p1n(2)**2+p1n(3)**2)
         ENDIF
c     rotate parton momenta after interaction, still in parton-parton frame
         call SIB_TRANI(P1N(1),P1N(2),P1N(3),cod,sid,cof,sif
     &        ,Px,Py,Pz)
         P1N(1)=Px
         P1N(2)=Py
         P1N(3)=Pz
         call SIB_TRANI(P2N(1),P2N(2),P2N(3),cod,sid,cof,sif
     &        ,Px,Py,Pz)
         P2N(1)=Px
         P2N(2)=Py
         P2N(3)=Pz
         IF(ndebug.gt.2) write(lun,*)' EXCT_RMNT: momentum*:',
     &        sqrt(p1n(1)**2+p1n(2)**2+p1n(3)**2)

c     boost back to hadron-hadron
         do k = 1,4
            gabe(k) = PTT(k)/PTT(5)
         enddo
         call SIB_ALTRA(gabe(4), gabe(1), gabe(2), gabe(3),
     &        P1n(1),p1n(2),p1n(3),p1n(4),
     &        P1TOT,p1(1),p1(2),p1(3),p1(4))
         p1(5)=p1n(5)
         call SIB_ALTRA(gabe(4), gabe(1), gabe(2), gabe(3),
     &        P2n(1),p2n(2),p2n(3),p2n(4),
     &        P2TOT,p2(1),p2(2),p2(3),p2(4))
         p2(5)=p2n(5)
         IF(ndebug.gt.2)THEN
            write(lun,*)' EXCT_RMNT: momenta after scattering:'
            write(lun,*)' EXCT_RMNT: hadron A :',(P1(i),i=1,5)
            write(lun,*)' EXCT_RMNT: string B :',(P2(i),i=1,5)
         ENDIF

c     edit partons on stack
         call edt_prtn
     &        (ISTH,P1(1),P1(2),P1(3),P1(4),P1(5),Iref)
         ICST11 = ICSTDX(2*(JW-1)+IMRG,2)
         call edt_prtn
     &        (IMST,P1(1),P1(2),P1(3),P1(4),P1(5),ICST11)
         ICST21 = ICSTDX(2*(JW-1)+IMRGBAR,2)
         call edt_prtn
     &        (ICST2,P2(1),P2(2),P2(3),P2(4),P2(5),ICST21)
         
      CASE(3)
         IF(ndebug.gt.2)
     &        WRITE(LUN,*)' EXCT_RMNT: merge both strings..'

c     both strings merged
c     index of merged string and its last parton
         IMST1 = ICSTDX(2*(JW-1)+1,1)
         IMST11 = ICSTDX(2*(JW-1)+1,3)
c     index of ordinary string
         IMST2 = ICSTDX(2*(JW-1)+2,1)
         IMST21 = ICSTDX(2*(JW-1)+2,3)
c     read 4momenta
         call rd_ref(IMST11,ISTH1)
         call rd_prtn_4vec(ISTH1,P1,L01,Iref)
c     string two
         call rd_ref(IMST21,ISTH2)
         call rd_prtn_4vec(ISTH2,P2,L02,Iref)
         xm1 = AM(IABS(L01))
         xm2 = AM(IABS(L02))
c     cm energy
         call add_4vecs(P1,P2,PTT)
         ecm = PTT(5)
         etot = ptt(4)
         IF(ndebug.gt.2)THEN
            write(lun,*)' EXCT_RMNT: string A :',(P1(i),i=1,5)
            write(lun,*)' EXCT_RMNT: string B :',(P2(i),i=1,5)
            write(lun,*)' EXCT_RMNT: total :',(PTT(i),i=1,5)
         ENDIF

         call transfOnShell(ecm,xm1,xm2,one,3,P1n,P2n,LBD)
         IF(LBD.eq.1) THEN
            IF(NDEBUG.gt.2)
     &           WRITE(LUN,*)' EXCT_RMNT: mass transfer failed!'
            RETURN
         ENDIF
c     by definition p1n is along +z in string cm, need to invert if pzA < pzB
c         IF(P2(3).gt.P1(3)) call swtch_lmnts(P1N(3),P2N(3))
c     rotate parton-parton axis onto string-string axis
c     boost to parton-parton cm to calc. rotation angles BEFORE interaction!
         do k = 1,4
            GABE(k) = PTT(k)/PTT(5)
         enddo
         call SIB_ALTRA(gabe(4),-gabe(1),-gabe(2),-gabe(3),
     &        P1(1),p1(2),p1(3),p1(4),
     &        P1TOT,pbm1(1),pbm1(2),pbm1(3),pbm1(4))
c     rotation factors
         COD= PBM1(3)/P1TOT
         SID= DSQRT(PBM1(1)**2+PBM1(2)**2)/P1TOT
         COF=ONE
         SIF=ZERO
         IF(P1TOT*SID.GT.EPS5) THEN
            COF=PBM1(1)/(SID*P1TOT)
            SIF=PBM1(2)/(SID*P1TOT)
            ANORF=DSQRT(COF*COF+SIF*SIF)
            COF=COF/ANORF
            SIF=SIF/ANORF
         ENDIF
c     rotate parton momenta after interaction
         call SIB_TRANI(P1N(1),P1N(2),P1N(3),cod,sid,cof,sif
     &        ,Px,Py,Pz)
         P1N(1)=Px
         P1N(2)=Py
         P1N(3)=Pz
         call SIB_TRANI(P2N(1),P2N(2),P2N(3),cod,sid,cof,sif
     &        ,Px,Py,Pz)
         P2N(1)=Px
         P2N(2)=Py
         P2N(3)=Pz

c     boost massive hadrons back to hadron-hadron
         call SIB_ALTRA(gabe(4), gabe(1), gabe(2), gabe(3),
     &        P1n(1),p1n(2),p1n(3),p1n(4),
     &        P1TOT,p1(1),p1(2),p1(3),p1(4))
         p1(5)=p1n(5)
         call SIB_ALTRA(gabe(4), gabe(1), gabe(2), gabe(3),
     &        P2n(1),p2n(2),p2n(3),p2n(4),
     &        P2TOT,p2(1),p2(2),p2(3),p2(4))
         p2(5)=p2n(5)
         IF(ndebug.gt.2)THEN
            write(lun,*)' EXCT_RMNT: hadron A :',(P1(i),i=1,5)
            write(lun,*)' EXCT_RMNT: hadron B :',(P2(i),i=1,5)
         ENDIF

c     edit partons on stack
         call edt_prtn
     &        (ISTH1,P1(1),P1(2),P1(3),P1(4),P1(5),Iref)
         ICST11 = ICSTDX(2*(JW-1)+1,2)
         call edt_prtn
     &        (IMST1,P1(1),P1(2),P1(3),P1(4),P1(5),ICST11)

         call edt_prtn
     &        (ISTH2,P2(1),P2(2),P2(3),P2(4),P2(5),Iref)
         ICST21 = ICSTDX(2*(JW-1)+2,2)
         call edt_prtn
     &        (IMST2,P2(1),P2(2),P2(3),P2(4),P2(5),ICST21)
         
      CASE default
         CONTINUE
      END SELECT

      IREJ = 0
      END

