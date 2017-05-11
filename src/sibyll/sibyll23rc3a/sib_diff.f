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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_plist_cmmn.inc'
      INCLUDE 'sib_difmass_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      INCLUDE 'sib_parto_cmmn.inc'
      INCLUDE 'sib_mass_cmmn.inc'
      INCLUDE 'sib_chp_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'
      
      DIMENSION P0(5),P1(5),P2(5)
      DIMENSION KK(39)
      DATA KK /5*0,7*2,2*1,6*0,13*2,6*1/
c      DIMENSION XM2MIN(3), ALXMIN(3)
c      DATA XM2MIN /1.5, 0.2, 0.6/                  ! M_x**2(min) GeV**2
c      DATA ALXMIN /0.405465,-1.6094379,-0.5108256/ ! log[M_x**2(min)]
c      DATA SLOP0 /6.5/                 ! b (slope_ for Mx**2 > 5 GeV**2
c      DATA ASLOP /31.10362/            ! fit to the slope parameter.
c      DATA BSLOP /-15.29012/

      if(Ndebug.gt.1) 
     &  WRITE(LUN,*)' SIB_DIFF: called with (L0,JDIF1,Ecm):',
     &  L0,JDIF1,Ecm

      if(Irec.eq.1) THEN
         Ipflag= -1
         IIFLAG = 1
c     add incoming target particles
         PZ = PAWT(ECM,AM(IABS(L0)),AM(13))
         E2 = SQRT(PZ**2+AM2(13))
         call add_prtn(ZERO,ZERO,-PZ,E2,AM(13),13,-2,0,IREFout)

c     add interactions
        xjdif = dble(jdif1)
        call add_prtn(zero,zero,xjdif,ecm,zero,1,-1,IREFout,IREF)
      ENDIF
      call GET_NPP(NPP_0,NPP0_0)

      IDBAD = 0
      NTRY = 0
 20   IREJ = 1
      call ini_prtn_stck(NPP_0,NPP0_0)

      IF(NTRY.gt.20*Irec) RETURN ! zero tolerance for recursive calls 
      NTRY = NTRY + 1
     
      LL = L0
      LA = IABS(L0)
      XM2MAX = PAR(13)*Ecm*Ecm
      if(Ndebug.gt.1) 
     &   WRITE(LUN,*)' SIB_DIFF: max diff. mass (M,Xi):',XM2MAX,PAR(13)
      
C...Double diffraction
      IF (JDIF1 .EQ. 3)   THEN
         K = MAX(1,2-IBAR(LA)-ISTR(LA)-ICHM(LA))
         IF(Irec.eq.1) K = KK(LA)
c     minimal mass if larger than particle mass plus one pion
         XMMIN = XM2MIN(K)
         IF(Irec.eq.0) XMMIN = MAX(XMMIN,(AM(LA)+AM(7)+0.02D0)**2)
         XMB2 = XM2DIS(XMMIN,XM2MAX,ONE)
         XMB = SQRT (XMB2)
         XMT2 = XM2DIS(XM2MIN(1),XM2MAX,ONE)
         XMT = SQRT (XMT2)
         call transfOnShell(ECM,XMB,XMT,XM2MAX,0,P1,P2,IBAD)
         IF(IBAD.ne.0) goto 20
         XMASS(1) = XMB
         IF(Irec.eq.1)THEN
c     add diffractive system to parton stack
            call add_prtn_4vec(P1,3,0,0,Iref)
            call add_int_ref(Iref,1)
            call add_prtn_4vec(P2,-3,0,0,Iref)
            call add_int_ref(Iref,1)
         ENDIF
         if(Ndebug.gt.1) 
     &        write(lun,*)'double-diff.: (kb,xmb,kt,xmt)',LL,xmb,13,xmt
         CALL DIFDEC (LL, Irec, IDBAD, P1)
         IF(IDBAD.eq.1)goto 20
         Ipflag= -2
         XMASS(2) = XMT
         CALL DIFDEC (13, Irec, IDBAD, P2)
         IF(IDBAD.eq.1)goto 20
         IREJ = 0
         RETURN
      ENDIF

C...Single diffraction
      IF (JDIF1.EQ. 1)  THEN
         K = MAX(1,2-IBAR(LA))
         IF(Irec.eq.1) K = KK(LA)
         EM  = AM(13)
         EM2 = AM2(13)
         L = 13
         ZD = -ONE
         if(Ndebug.gt.1) 
     &        write(lun,*)'single-diff. (beam): (kb)',LL
      ELSE
         K = 1
         EM  = AM(LA)
         EM2 = AM2(LA)
         L = LL
         LL = 13
         ZD = +ONE
         if(Ndebug.gt.1) 
     &        write(lun,*)'single-diff. (target): (kt)', LL

      ENDIF
C...Generate the mass of the diffracted system Mx (1/Mx**2 distribution)
      XMMIN = XM2MIN(K)
      IF(Irec.eq.0) XMMIN = MAX(XMMIN,(AM(LA)+AM(7)+0.02D0)**2)
      XM2 = XM2DIS(XMMIN,XM2MAX,ONE)
      ALX = log(XM2)
c... added part
      X = XM2/XM2MAX*PAR(13)
      IF (X.GT.PAR(13)-0.05D0) THEN
        PRO = HALF*(ONE+(X-PAR(13))/0.05D0)
        IF (S_RNDM(0).LT.PRO) X = TWO*PAR(13)-X
        XM2 = XM2MAX*X/PAR(13)
      ENDIF
c...

      XM = SQRT (XM2)
      XMB = XM
      XMT = XM
      XMASS(1) = XMB
      XMASS(2) = XMT

C..   kinematics
      call transfOnShell(ECM,XMB,EM,XM2MAX,0,P1,P2,IBAD)
      IF(IBAD.ne.0) goto 20

C...Generate the Kinematics of the pseudoelastic hadron
      NP = NP+1
      P(NP,4) = P2(4)
      P(NP,3) = abs(P2(3))*ZD
      P(NP,1) = p2(1)
      P(NP,2) = p2(2)
      P(NP,5) = EM
      LLIST(NP) = L
      NPORIG(NP) = IPFLAG
      niorig(NP) = iiflag

C...Generating the hadronic system recoiling against the produced particle
      P0(5) = SQRT(XM2)
      P0(4) = P1(4)
      DO J=1,3
         P0(J) = -P(NP,J)
      ENDDO
      IF(Irec.eq.1)THEN
c     add diffractive system to parton stack
         call add_prtn_4vec(P1,JDIF1,0,0,Iref)
         call add_int_ref(Iref,1)
         call add_prtn_4vec(P2,int(zd),0,0,Iref)
         call add_int_ref(Iref,1)
      ENDIF
      CALL DIFDEC (LL, Irec, IDBAD, P0)
      IF(IDBAD.eq.1)goto 20
      IREJ = 0

      END
