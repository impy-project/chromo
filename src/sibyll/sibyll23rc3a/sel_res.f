      subroutine sel_res(XM2in,KDin,IRDX,IKDH)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C--------------------------------------------------------------------
C     routine that checks if excitation should go into resonant state
C     or rather should fallback to on-shell beam hadron
C     Input: XM2in : squared excitation mass
C            KDin : projectile hadron code
C            IRDX : reference to remnant on stack
C     Output: adds hadron to stack
C             IKDH : parton stack index of final hadron
C--------------------------------------------------------------------
      include 'sib_debug_cmmn.inc'
      include 'sib_cflafr_cmmn.inc'
      include 'sib_mass_cmmn.inc'
      include 'sib_width_cmmn.inc'
      DIMENSION LRES(6:99,2)
      DATA (LRES(k,1),k=6,22)  /27,25,26,28,29,0,0,51,52,6*0,30,31/
      DATA (LRES(k,1),k=23,33) /23,24,25,26,27,28,29,30,31,27,27/
      DATA (LRES(k,1),k=34,49) /34,35,36,37,38,39,40,41,42,43,34,35,36,
     &     37,38,49/
      DATA (LRES(k,1),k=50,83) /0,51,52,53,54,4*0,78,79,10*0,80,81,73,
     &     74,75,76,77,78,79,80,81,0,83/
      DATA (LRES(k,1),k=84,99) /94,95,96,97,98,89,4*0,94,95,96,97,98,99/

      DATA (LRES(k,2),k=6,22)  /61,62,63,64,65,0,0,53,54,6*0,66,67/
      DATA (LRES(k,2),k=23,33) /61,61,62,63,61,64,65,66,67,61,61/
      DATA (LRES(k,2),k=34,49) /34,35,36,37,38,39,40,41,42,43,44,45,46,
     &     47,48,49/
      DATA (LRES(k,2),k=50,83) /0,51,52,53,54,4*0,78,79,10*0,80,81,73,74
     &     ,75,76,77,78,79,80,81,0,83/
      DATA (LRES(k,2),k=84,99) /94,95,96,97,98,89,4*0,94,95,96,97,98,99/

      XM2 = XM2in
      XM1 = sqrt(XM2)
      KD = KDin

C     thresholds
c     fallback threshold
      EMIN1 = PAR(76)
      
c     resonance threshold
      EMIN2 = PAR(77)
      
c     parton stack index of incoming hadron
      IKDH = 0
      
c     if too low, fallback on beam
      IF(ndebug.gt.3)
     &     write(lun,*)' SEL_RES: input (XM2in,KDin):',XM2,KD
      DELTAE = XM1-AM(ABS(KD))
      IF(ndebug.gt.3)then
         write(lun,*)' SEL_RES: DELTAE,EMIN1,EMIN2',deltae,emin1,emin2
         write(lun,*)' SEL_RES: XM,XM1,XM2',
     &        XM1,emin1+AM(ABS(KD)),emin2+AM(ABS(KD))
      endif
      IF(DELTAE.LT.EMIN1)THEN
c     fallback to beam region
         KDH = kd
         XM1 = AM(abs(kd))
         XM2 = AM2(abs(kd))

      ELSEIF(DELTAE.LT.EMIN2)THEN
c     form resonance
         II = 1
         KDH = KD
         DO WHILE (II.le.2.and.KDH.eq.KD)
            KDD = IABS(KD)
c     K0s and K0l projection on K0 and K0bar
            IF(KDD.eq.11.or.KDD.eq.12)KDD=21+INT((TWO-EPS10)*S_RNDM(KD))
            IL = LRES(KDD,II)
            IF(ndebug.gt.3)
     &           write(lun,*) ' SEL_RES: res. select (KD,II,ISPN,IL):',
     &           KD,II,ISPN,IL
            IF(IL.eq.0) write(lun,*) 'KD:' , KD
c     sample probability for resonance to occur at this mass
c     from the relativistic breit-wigner dist.
c     scale widths to artificially increase or decrease resonance occurence
            XWDTH = PAR(94)*AW2(IL)
            PRES = BREIT_WIGNER(XM2,AM2(IL),XWDTH)
            IF(ndebug.gt.3)
     &           write(lun,*)
     &           ' SEL_RES: res. proposal (AM2,AW2,Prob.):',
     &           AM2(IL),XWDTH,PRES
            IF(S_RNDM(ii).lt.PRES) KDH = ISIGN(IL,KD)            
            II = II + 1
         ENDDO
c     no resonance selected, fallback to beam or phasespace decay?
         IF(IPAR(59).eq.1.and.KDH.eq.KD)THEN
c     distinguish regions in deltaE
            IF(DELTAE.LT.EMIN1)THEN
c     fallback to beam
               XM1 = AM(abs(kdh))
               XM2 = AM2(abs(kdh))           
            ELSE
               KDH = 0
            ENDIF
         ELSE
c     case where resonance has been selected 
c     or no overlap between resonance and phasespace region exists
c     set mass to pole masses of selected particles
            XM1 = AM(abs(kdh))
            XM2 = AM2(abs(kdh))
         ENDIF
      ELSE
c     neither resonance nor fallback
         KDH = 0
      ENDIF
      IF(KDH.ne.0)THEN
c     add new beam hadron to stack
         XM2in = XM2
         call add_prtn
     &        (ZERO,ZERO,ZERO,ZERO,XM1,KDH,2,IRDX,IKDH)
      endif
      IF(ndebug.gt.3)
     &     write(lun,*)' SEL_RES: output (XM2in,KDin,KDH):',XM2,KD,KDH
      end


      DOUBLE PRECISION FUNCTION BREIT_WIGNER(S,XM2,XWDTH2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     peak set to one
      x1 = (s-xm2)**2+xm2*xwdth2
      breit_wigner = xm2*xwdth2/x1
      end

      DOUBLE PRECISION FUNCTION TBREIT_WIGNER(S,XM2,XWDTH2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     breit-wigner truncated at 2*gamma from peak
C     peak set to one
      DATA N /10/
      DATA ZERO,ONE /0.D0,1.D0/
      XMLOW = MAX(XM2-N*XWDTH2,ZERO)
      XCUT = sign(one,S-XMLOW)
      XCUT = MAX(XCUT,ZERO)
      x1 = (S-xm2)**2+xm2*xwdth2
      tbreit_wigner = xcut * xm2*xwdth2/x1
      
      end
