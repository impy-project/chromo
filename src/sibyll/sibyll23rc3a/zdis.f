C-----------------------------------------------------------------------
C     fragmentation functions in SIBYLL                        \FR'14
C-----------------------------------------------------------------------

      FUNCTION ZDIS_4FLV (IFL1,IFL2, XMT2)
C...z distribution
c     includes charmed fragmentation (Peterson/SLAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      INCLUDE 'sib_czdis_cmmn.inc'
      INCLUDE 'sib_czdiss_cmmn.inc'
      INCLUDE 'sib_czdisc_cmmn.inc'
      include 'sib_nw_prm.inc'
      INCLUDE 'sib_run_cmmn.inc'
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'
c$$$      COMMON /S_CZDIS/ FAin, FB0in
c$$$      COMMON /S_CZDISs/ FAs1, fAs2
c$$$      COMMON /S_CZDISc/ ZDMAX, EPSI
c$$$      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, kb ,kt
c$$$      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
c$$$      DATA ZERO,HALF,ONE,TWO,THREE /0.D0,0.5D0,1.D0,2.D0,3.D0/      

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
c     parameters for hard scattering (gluon) fragmentation
         IF(IPAR(6).eq.2)THEN
            fa= PAR(18)
            fb0= PAR(19)
         ENDIF   
c     special parameters for strange fragmentation
c     only active for baryon beams (or K0,K0bar)
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
         IF(FA.GT.0.01D0.AND.ABS(FA-ONE)/FB.LE.0.01D0)
     +         ZMAX=FB/(ONE+FB)+(ONE-FA)*FB**2/(ONE+FB)**3
         IF(FA.GT.0.01D0.AND.ABS(FA-ONE)/FB.GT.0.01D0)
     +     ZMAX=HALF*(ONE+FB-dSQRT((ONE-FB)**2+4.D0*FA*FB))/(ONE-FA)
         IF(ZMAX.LT.0.1D0)  ZDIV=2.75D0*ZMAX
         IF(ZMAX.GT.0.85D0) 
     +        ZDIV=ZMAX-0.6D0/FB**2+(FA/FB)*dLOG((0.01D0+FA)/FB)
C...  Choice if z, preweighted for peaks at low or high z
 100     Z=S_RNDM(0)
         IDIV=1
         FPRE=ONE
         IF (ZMAX.LT.0.1D0)  THEN
            IF(ONE.LT.S_RNDM(0)*(ONE-dLOG(ZDIV)))  IDIV=2
            IF (IDIV.EQ.1)  Z=ZDIV*Z
            IF (IDIV.EQ.2)  Z=ZDIV**Z
            IF (IDIV.EQ.2)  FPRE=ZDIV/Z
         ELSEIF (ZMAX.GT.0.85D0)  THEN
            IF(ONE.LT.S_RNDM(0)*(FB*(ONE-ZDIV)+ONE)) IDIV=2
            IF (IDIV.EQ.1)  Z=ZDIV+dLOG(Z)/FB
            IF (IDIV.EQ.1)  FPRE=dEXP(FB*(Z-ZDIV))
            IF (IDIV.EQ.2)  Z=ZDIV+Z*(ONE-ZDIV)
         ENDIF
C...weighting according to the correct formula
         IF (Z.LE.FB/(50.D0+FB).OR.Z.GE.ONE)  GOTO 100
         FVAL=(ZMAX/Z)*dEXP(FB*(ONE/ZMAX-ONE/Z))
         IF(FA.GT.0.01D0)  FVAL=((ONE-Z)/(ONE-ZMAX))**FA*FVAL
         IF(FVAL.LT.S_RNDM(0)*FPRE)  GOTO 100
         ZDIS_4FLV=Z
         
      ENDIF
      
      RETURN
      END
      
      SUBROUTINE ZNORMAL
C...  normalisation for Peterson/SLAC frag. func
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'sib_czdisc_cmmn.inc'
      INCLUDE 'sib_debug_cmmn.inc'
c$$$      COMMON /S_CZDISc/ ZDMAX, EPSI
c$$$      COMMON /S_DEBUG/ Ncall, Ndebug, Lun

c     get the maximum zmefn value first for normalisation
      jmax = 1000
      zdmax = 1.D-10
	
      DO j = 1, jmax
         z = dble(j)/dble(jmax+1)
         zdmax = max(zdmax, zmefn(z,epsi))
      enddo
      WRITE(LUN,*)'ZDMAX,EPS:',zdmax, epsi
      RETURN
      END
	
      FUNCTION ZMEFN(z,eps)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C...  Peterson/SLAC frag. func
      zmefn = (z*(1.D0-z**(-1)-eps/(1.D0-z))**2)**(-1)
      RETURN
      END
	

      FUNCTION ZBLEAD (LB)
C...fragmentation function for leading baryon
C.  simple form:  f(z) = a + x**b
C   INPUT : LB = particle code.
C..................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      INCLUDE 'sib_czlead_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      INCLUDE 'sib_chp_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'
c$$$      COMMON /S_CZLEAD/ CLEAD, FLEAD
c$$$      COMMON /S_CHP/ ICHP(99), ISTR(99), IBAR(99)
c$$$      DATA ZERO,HALF,ONE,TWO,THREE /0.D0,0.5D0,1.D0,2.D0,3.D0/

c      ncall = ncall + 1
c      print*,'leading baryon frag. called:',lb,ncall

C...  leading z lower bound
c     used for protons only in Sib21 (if ..)
c     used for all baryons alike in Sib22 (else..)    
      ZLMIN = PAR(55)
      ZSMR = PAR(56)

      IF(IPAR(30).ne.0)THEN
C     Sibyll 2.1 hard fragmentation function

      IC = ICHP(Lb)*ISIGN(1,Lb)
      
      if (lb.ge.34.and.lb.le.39)  then ! Lambda's and Sigma's
         IF(IPAR(35).eq.1)then
            zblead=zdisn(1)     ! zblead**2   !soft
         ELSE
 665        ZBLEAD = S_RNDM(0)
            if (zblead.le.0.01D0) goto 665
         ENDIF
c     zblead=zdisn(1) ! blead**2   ! soft
      else if (ic.eq.0)     then
         if(ipar(30).eq.2)then
 555        zblead = s_rndm(0)
            if (zblead .le. 0.01D0) goto 555     
         else
            zblead=zdisn(1)     ! blead**2   !soft
         endif
      else if (ic.eq.1)  then   ! fast protons only
         if (abs(lb).eq.13) then
 661        IF (S_RNDM(0) .LT. CLEAD)  THEN
 666           ZBLEAD = S_RNDM(0)
               if (zblead.le.0.01D0) goto 666
            ELSE
               zblead=ONE-zdisn(1) ! zblead**2   !hard
            ENDIF
c     truncated zblead to fix antiprotons
            if (zblead.le.ZLMIN+ZSMR*(ONE-TWO*S_RNDM(LB))) goto 661
         else
            zblead=zdisn(1)     ! zblead**2   !hard
         endif   
      else if (ic.eq.2)  then   ! fast delta++
         zblead=ONE- zdisn(1)    ! (zblead)**.3333
      else
         zblead=S_RNDM(0)       ! zdisn(1)     !hard
      endif
      RETURN
      ELSE
C...  Sein's flat baryon fragmentation function a.k.a. Sibyll 2.2
 999     zblead = s_rndm(0)
         if (zblead .le. 0.01D0) goto 999
c     truncated zblead to fix instring pair production (antiprotons)
         if (zblead.le.ZLMIN+ZSMR*(ONE-TWO*S_RNDM(LB))) goto 999
         RETURN
      ENDIF
      END


      FUNCTION ZDISN (n)
C...Generate (1-x)**n
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
666   rmin=1.1D0
      do i=1,n+1
         R1=S_RNDM(0)
         IF (R1.LE.RMIN) RMIN=R1
      ENDDO
      ZDISn=RMIN
      if (zdisn.le.0.01D0) goto 666
      if (zdisn.ge.0.99D0) goto 666
      END
