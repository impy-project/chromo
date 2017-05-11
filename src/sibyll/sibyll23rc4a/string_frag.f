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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_plist_cmmn.inc'
      INCLUDE 'sib_zlist_cmmn.inc'
      INCLUDE 'sib_mass_cmmn.inc'
      include 'sib_nw_prm.inc'
      INCLUDE 'sib_run_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      INCLUDE 'sib_parto_cmmn.inc'
      INCLUDE 'sib_czdis_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'
      INCLUDE 'sib_rnk_cmmn.inc'
      INCLUDE 'sib_chp_cmmn.inc'

      DIMENSION WW(2,2), PTOT(4), PX(3),PY(3),IFL(3),ILEAD(2)
      DIMENSION LPOINT(8000), PMQ(3), IRNK(2), LRES(6:99)
      LOGICAL LRANK
      DATA LRANK/.true./

      DATA LRES(6:39) /27,25,26,28,29,9,9,41,42,19*0,44,45,46,47,48,39/
      DATA LRES(40:49) /40,41,42,43,44,45,46,47,48,49/
      DATA LRES(50:83) /0,51,52,53,54,4*0,78,79,10*0,71,72,73,76,77,76,
     &     77,78,79,80,81,0,83/
      DATA LRES(84:99) /94,95,96,97,98,89,4*0,94,95,96,97,98,99/
      
      IF(Ndebug.gt.3) THEN
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
         WW(K,1) = ONE
         WW(K,2) = ZERO
         IRNK(K) = 0
      ENDDO
      PX(1) = PX1
      PY(1) = PY1
      PX(2) = PX2
      PY(2) = PY2
      PX(3) = ZERO
      PY(3) = ZERO
      PTOT (1) = PX1+PX2
      PTOT (2) = PY1+PY2
      PTOT (3) = ZERO
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
c     switch leading baryon fragmentation function on/off
      IF(IPAR(20).eq.0) GOTO 300
c     set flags for leading baryon
C
C      SET FLAG FOR GENERATION OF LEADING PARTICLES. 
C      "AND" IS FOR PPBAR ( DIQUARK AT BOTH ENDS)
C      "OR" IS FOR PP, PPI, ( DIQUARK AT ONE END.)
C
      IF (IABS(IFL1) .GT. 10 .AND. IABS(IFL2) .GT. 10)  THEN
         IBLEAD = 2
         I = I+1
         JT = 1.5D0+S_RNDM(0)
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
        write(6,'(1x,a,i8)') 
     &        'STRING_FRAG_4FLV: no space left in S_PLIST:',I
        call sib_reject
      endif
      IF (IBLEAD .GT. 0)  THEN
         JT = 3 - JT   
         GO TO 350              
      ENDIF
c     
 349  continue
c     choose side (1 or 2)
      JT=1.5D0+S_RNDM(0)                 
c     set 'other' side
 350  JR=3-JT
c     remember side particle was produced
      LPOINT(I) = JT
c     increase rank counter
      IRNK(JT) = ISIGN(ABS(IRNK(JT))+1,1-JT)
c     set particle rank
      LRNK(I) = IRNK(JT)

      nporig(I)= Ipflag*2 + Nint
      niorig(I)= iiflag
      IF(ILEAD(JT).eq.1) nporig(I)= -1 * nporig(I)
      nforig(I) = 0

 555  CONTINUE
c
c.... CHARM config
c
      charmPARdef=PAR(24)
      IF(IPAR(15).lt.9)THEN
c     no s->c
         PAR(24) = ZERO
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
c     increase vec.meson ratio for leading particle in str. diff.
      IF(IFQRK.eq.-1)THEN
         select case(ipar(66))
         case(1)
            IF(ILEAD(JT).EQ.1)THEN
              IF(IBAR(IABS(kb)).eq.0.or.IPAR(70).eq.1) PAR(5) = PAR(113)
           ENDIF
           
         case(2)
            IF(IBAR(IABS(kb)).eq.0.or.IPAR(70).eq.1) PAR(5) = PAR(113)
           
         case(3)
c     increase vector meson rate for meson beam
c     on beam side (rank+) only!            
            IF(ILEAD(JT).EQ.1)THEN               
               IF(IBAR(IABS(kb)).eq.0.and.IRNK(JT).gt.0)
     &              PAR(5) = PAR(113)
c     always incr. vector rate for diff. strings independent of beam type
               IF(IPAR(70).eq.1) PAR(5) = PAR(113)               
            ENDIF         
            
         END select
      endif
      
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
     &        PAR(5) = 9.D0
         
c     increase vec.meson ratio for diff.
         IF(IFQRK.eq.-1.and.IPAR(11).le.-4.and.IPAR(11).ge.-7)
     &        PAR(5) = 9.D0

c     increase vec.meson ratio for leading particle in str. diff. (lvec16)
         IF(IFQRK.eq.-1.and.IPAR(11).le.-11.and.ILEAD(JT).EQ.1)
     &        PAR(5) = 99.D0
      ENDIF

c...  suppress leading charm for pion and kaon beams
      IF(IPAR(15).eq.11)then
         IF((1-IABS(KB)/13)*ILEAD(JT).gt.0) PAR(24)=ZERO
      ENDIF

C...  suppress rank-1 baryon through popcorn
      IF(IBLEAD .GT. 0.and.abs(ifl(jt)).gt.10
     &     .and.abs(ifl(3)).lt.10) PAR(8)=PAR(63)*PAR(8)

C...  leading strange/charm
      IF(ILEAD(JT).eq.1.and.IPAR(39).gt.0) PAR(2) = PAR(65)

c     scale valence string end charm for assoc. prod.      
      IF(IPAR(41).eq.1)THEN
         IF(ILEAD(JT).eq.1.and.IFQRK.eq.1) PAR(24) = PAR(71)*PAR(24)
      ENDIF

c     suppress direct pi0 for meson projectiles
c     rate set by par( 137 )
      ipar82_def = ipar(82)
c     skip if baryon projectile or minijet (i.e no flavor attached)
      if(ibar(iabs(kb)).ne.0.or.ifqrk.eq.0) ipar(82) = 0

c     suppress direct omega for meson projectiles
c     rate set by par( 138 )
      ipar83_def = IPAR(83)
c     skip if baryon projectile or central string
      if(ibar(iabs(kb)).ne.0.or.(ifqrk.gt.0.and.IPAR(83).eq.2))
     &     IPAR(83) = 0

c     change rho0 / omega ratio
      PAR143_def = PAR(143)
      SELECT CASE(IPAR(81))
      CASE(1)
c     change if beam is meson
         if(ibar(iabs(kb)).eq.0) PAR(143) = PAR(144)
      CASE(2)
c     change if beam is meson, on meson side only         
         if(ibar(iabs(kb)).eq.0.and.IRNK(JT).gt.0) PAR(143) = PAR(144)         
      CASE(3)
c     change if beam is meson, on meson side only, for leading only
         if(ibar(iabs(kb)).eq.0.and.ISIGN(ILEAD(JT),IRNK(JT)).eq.1)
     &        PAR(143) = PAR(144)
      CASE(4)
c     change if beam is meson, on meson side only, for diff. strings only
         if(ibar(iabs(kb)).eq.0.and.IFQRK.eq.-1)
     &        PAR(143) = PAR(144)
      CASE(5)
c     change if beam is meson, for leading on meson side only and
c     for diff. strings only
         if(ibar(iabs(kb)).eq.0.and.IFQRK.eq.-1.and.
     &        ISIGN(ILEAD(JT),IRNK(JT)).eq.1) PAR(143) = PAR(144)

      END SELECT        
      
C...particle ID and pt.

      CALL SIB_I4FLAV (IFL(JT), 0, IRNK(JT), IFL(3), LLIST(I))

c     reset strange fraction
      PAR(2) = PAR2_def
c     reset vec.meson production
      PAR(5) = PAR5_def
c     reset charm fraction
      PAR(24) = PAR24_def
c     reset popcorn
      PAR(8) = par8_def

c     reset pi0 suppr.
      IPAR(82) = ipar82_def

c     reset omega suppr.
      IPAR(83) = ipar83_def

c     reset rho0 / omega ratio
      PAR(143) = PAR143_def
      
c     reject iso 0 spin 1 for meson projectiles
      IF(IBAR(IABS(KB)).eq.0)THEN
c     reject leading spin1,isospin singlett
         IF(ILEAD(JT).EQ.1.and.LLIST(I).eq.32.and.
     +        PAR(136).gt.S_RNDM(I)) LLIST(I) = 27
      endif
      
c     replace leading or all pi0 with rho0
      IF(IFQRK.eq.-1) THEN 
         IF(IPAR(67).eq.1)THEN
            IF(ILEAD(JT).EQ.1) THEN 
c     replace leading pi0 with rho0
               IF(ABS(LLIST(I)).EQ.6) LLIST(I) = 27*isign(1,LLIST(I))           
            ENDIF
         ELSEIF(IPAR(67).eq.2)THEN
c     replace all pi0 with rho0 for all beams
            IF(ABS(LLIST(I)).EQ.6) LLIST(I) = 27*isign(1,LLIST(I))
         ELSEIF(IPAR(67).eq.3)THEN
c     replace all pi0 with rho0 for meson beam only
            IF(IBAR(IABS(KB)).eq.0)THEN
               IF(ABS(LLIST(I)).EQ.6) LLIST(I) = 27*isign(1,LLIST(I))
            ENDIF
         ELSEIF(IPAR(67).eq.4)THEN
c     replace all pi0 with rho0 for meson beam only
c     replace some beam mesons with their vector partner
            IF(IBAR(IABS(KB)).eq.0)THEN
               IF(ABS(LLIST(I)).EQ.6) LLIST(I) = 27*isign(1,LLIST(I))
c     reject leading spin1,isospin singlett
               IF(ILEAD(JT).EQ.1.and.LLIST(I).eq.32.and.
     +              PAR(136).gt.S_RNDM(I)) LLIST(I) = 27
               IF(s_rndm(0).lt.PAR(120).and.LLIST(I).eq.KB) 
     &              LLIST(I) = LRES(LLIST(I))
            ENDIF
         ENDIF
      ENDIF

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
            PAR(5) = 9.D0
            IF(ilead(jt).EQ.1.and.INT((PAR(5)+1)*S_RNDM(0)).gt.1) then 
               IF(ABS(LLIST(I)).EQ.6) LLIST(I) = 27*isign(1,LLIST(I))
               IF(ABS(LLIST(I)).EQ.7) LLIST(I) = 25*isign(1,LLIST(I))
c     IF(ABS(LLIST(I)).EQ.8) LLIST(I) = 26*isign(1,LLIST(I))
            endif

c     replace almost all for diff.
         ELSEIF(IFQRK.eq.-1.and.ipar(11).eq.-8.and.ipar(11).lt.0) then
            PAR(5) = 9.D0
            if( INT((PAR(5)+1)*S_RNDM(0)).gt.1 ) then
               IF(ABS(LLIST(I)).EQ.6)  LLIST(I) = 27*isign(1,LLIST(I))
               IF(ABS(LLIST(I)).EQ.7) LLIST(I) = 25*isign(1,LLIST(I))
            endif
      
c     replace leading pi0's by vec.mesons
         ELSEIF(IPAR(11).eq.-9.and.IPAR(11).lt.0)THEN
            PCHF = 0.1D0
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
      PAR(5) = 0.3D0
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
      IF (WREM2 .LT. 0.1D0)  GOTO 200
c      WMIN = PMQ(1)+PMQ(2)+2.*PMQ(3)+ 1.1 + (2.*S_RNDM(0)-1.)*0.2
      WMIN=PMQ(1)+PMQ(2)+TWO*PMQ(3)+PAR(59)+(TWO*S_RNDM(0)-ONE)*0.2D0
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
               IF(IBLEAD .GT. 0.and.abs(ifl(jt)).gt.10
     &              .and.abs(ifl(3)).lt.10)  THEN
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
               IF(IBLEAD .GT. 0.and.abs(ifl(jt)).gt.10
     &              .and.abs(ifl(3)).lt.10)  THEN
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
               Z = ONE - ZDISN(1)
            ELSE
               Z = ZDIS_4FLV (IFL(3),ifl(jt),XMT2)
            ENDIF

         ELSEIF(IPAR(11).EQ.-12.OR.IPAR(11).LE.-15.or.IPAR(68).eq.1)THEN
c     very hard leading frag. for diff. val. strings only (lvec17)
            IF(IBLEAD .GT. 0.and.abs(ifl(jt)).gt.10) THEN
               Z = ZBLEAD (IABS(LLIST(I)))
            ELSEIF(ILEAD(jt).eq.1.and.IFQRK.eq.-1)THEN
               Z = ONE - ZDISN(1)
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
            IF(IBLEAD .GT. 0.and.abs(ifl(jt)).gt.10
     &           .and.abs(ifl(3)).lt.10)  THEN
c           print*,'calling zblead: i,id,jt,ncall', i,llist(i),jt,ncall
               IF(IPAR(20).eq.3)THEN
c     use lund function with different parameters for leading baryon
                  fa_def = FAin
                  fb_def = FB0in
                  FAin = PAR(57)
                  FB0in = PAR(58)
                  z = zdis_4flv(IFL(3),ifl(jt),xmt2)
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
         IF (IBLEAD .GT. 0.and.abs(ifl(jt)).gt.10
     &        .and.abs(ifl(3)).lt.10)  THEN
C     Special frag. for leading Baryon only
c            print*,'calling zblead: i,id,jt,ncall', i,llist(i),jt,ncall
            Z = ZBLEAD (IABS(LLIST(I)))   
            IBLEAD = IBLEAD - 1
         ELSE
            Z = ZDIS_4FLV (IFL(3),ifl(jt),XMT2)
         ENDIF
      ENDIF
      IF(IPAR(20).eq.2)IBLEAD = 2
      IF(IFQRK.eq.1) ILEAD(JT) = 0

      ZLIST(I) = Z
      WW(JT,2) = Z*WW(JT,1)
      WW(JR,2) = XMT2/(WW(JT,2)*E0**2)

      P(I,3) = WW(1,2)*0.5D0*E0 - WW(2,2)*0.5D0*E0
      P(I,4) = WW(1,2)*0.5D0*E0 + WW(2,2)*0.5D0*E0

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
c     increase vec.meson ration for diff.
      IF(IFQRK.eq.-1.and.IPAR(11).le.-4.and.IPAR(11).gt.-8) PAR(5) =9.D0
c     increase vec.meson ration for leading quarks in valence interactions
      IF(IABS(IFQRK).eq.1.and.IPAR(11).le.-5.and.ilead(jr).eq.1
     &     .and.IPAR(11).gt.-8) PAR(5) = 9.D0

c     suppress direct pi0 for meson projectiles
c     rate set by par( 137 )
 666  ipar82_def = ipar(82)
c     skip if baryon projectile
      if(ibar(iabs(kb)).ne.0.or.ifqrk.eq.0) ipar(82) = 0

c     suppress direct omega for meson projectiles
c     rate set by par( 138 )
      ipar83_def = IPAR(83)
c     skip if baryon projectile or central string     
      if(ibar(iabs(kb)).ne.0.or.(ifqrk.gt.0.and.IPAR(83).eq.2))
     &     IPAR(83) = 0

c     set current rank
      IRNK(JR)=ISIGN(ABS(IRNK(JR))+1,1-JR)
      
c     change rho0 / omega ratio
      SELECT CASE(IPAR(81))
      CASE(1)
c     change if beam is meson
         if(ibar(iabs(kb)).eq.0) PAR(143) = PAR(144)
      CASE(2)
c     change if beam is meson, on meson side only         
         if(ibar(iabs(kb)).eq.0.and.IRNK(JR).gt.0) PAR(143) = PAR(144)         
      CASE(3)
c     change if beam is meson, on meson side only, for leading only
         if(ibar(iabs(kb)).eq.0.and.ISIGN(ILEAD(JR),IRNK(JR)).eq.1)
     &        PAR(143) = PAR(144)
      CASE(4)
c     change if beam is meson, on meson side only, for diff. strings only
         if(ibar(iabs(kb)).eq.0.and.IFQRK.eq.-1)
     &        PAR(143) = PAR(144)
      CASE(5)
c     change if beam is meson, for leading on meson side only and
c     for diff. strings only
         if(ibar(iabs(kb)).eq.0.and.IFQRK.eq.-1.and.
     &        ISIGN(ILEAD(JR),IRNK(JR)).eq.1) PAR(143) = PAR(144)
         
      END SELECT

c     increase vec.meson ratio for leading particle in str. diff.
      select case(ipar(66))
      case(1)        
         IF(ILEAD(JT).EQ.1.and.IFQRK.eq.-1)THEN
            IF(IBAR(IABS(kb)).eq.0.or.IPAR(70).eq.1) PAR(5) = PAR(113)
         ENDIF
           
      case(2)
         IF(IFQRK.eq.-1)THEN
            IF(IBAR(IABS(kb)).eq.0.or.IPAR(70).eq.1) PAR(5) = PAR(113)
         ENDIF

      case(3)
c     increase vector meson rate for meson beam
c     on beam side (rank+) only!
         IF(IFQRK.eq.-1)THEN
            IF(ILEAD(JR).EQ.1)THEN               
               IF(IBAR(IABS(kb)).eq.0.and.IRNK(JR).gt.0)
     &              PAR(5) = PAR(113)
c     always incr. vector rate for diff. strings independent of beam type
               IF(IPAR(70).eq.1) PAR(5) = PAR(113)               
            ENDIF
         ENDIF
      END select

      
      CALL SIB_I4FLAV (IFL(JR), -IFL(3), IRNK(JR), IFLA, LLIST(I+1))

      ipar(82) = ipar82_def
      IPAR(83) = ipar83_def
      PAR(143) = PAR143_def
      
      nporig(I+1)= Ipflag*2 + Nint
      niorig(I+1)= iiflag
      IF(ILEAD(1).eq.1.or.ILEAD(2).eq.1) nporig(I+1)= -1 * nporig(I+1)

c     replace leading or all pi0 with rho0
      IF(IFQRK.eq.-1) THEN 
         IF(IPAR(67).eq.1)THEN
            IF(ILEAD(JR).EQ.1) THEN 
               IF(ABS(LLIST(I+1)).EQ.6) 
     &              LLIST(I+1) = 27*isign(1,LLIST(I+1))           
            ENDIF
         ELSEIF(IPAR(67).eq.2)THEN
            IF(ABS(LLIST(I+1)).EQ.6) LLIST(I+1) = 27*isign(1,LLIST(I+1))           
         ELSEIF(IPAR(67).eq.3)THEN
            IF(IBAR(IABS(KB)).eq.0)THEN
               IF(ABS(LLIST(I+1)).EQ.6)LLIST(I+1)=27*isign(1,LLIST(I+1))
            ENDIF
         ENDIF
      ENDIF
      
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
            PAR(5) = 9.D0
            if( INT((PAR(5)+1)*S_RNDM(0)).gt.1 ) then
               IF(ABS(LLIST(I+1)).EQ.6)
     &              LLIST(I+1) = 27*isign(1,LLIST(I+1))
               IF(ABS(LLIST(I+1)).EQ.7)
     &              LLIST(I+1) = 25*isign(1,LLIST(I+1))
            endif
         endif
c     replace leading pseudoscalar by vector
         IF(ipar(11).eq.-8.and.ilead(jr).eq.1) then
            PAR(5) = 9.D0
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

c     reject iso 0 spin 1 (omega) for meson projectiles
      IF(IBAR(IABS(KB)).eq.0)THEN
c     reject leading spin1,isospin singlett
         IF(ILEAD(JR).EQ.1.and.LLIST(I+1).eq.32.and.
     +        PAR(136).gt.S_RNDM(I)) LLIST(I+1) = 27
      endif
      
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
      LRNK(I1) = IRNK(JR)
      XM1 = P(I,5)**2+P(I,1)**2+P(I,2)**2
      XM2 = P(I1,5)**2+P(I1,1)**2+P(I1,2)**2
      IF (dSQRT(XM1)+dSQRT(XM2) .GT. dSQRT(WREM2)) GOTO 200

c...RE & EJA fix
      PT2 = (P(I,1)+P(I1,1))**2+(P(I,2)+P(I1,2))**2
      WREMPT = dsqrt(WREM2+PT2)
      EA1 = (WREM2+XM1-XM2+PT2)/(TWO*WREMPT)

      PA2 = (EA1**2-XM1)
      if (pa2.gt.ZERO)  then
            PA = dSQRT(PA2)
      else
            goto 200
      endif
      BA = PTOT(3)/PTOT(4)
      GA = PTOT(4)/WREMPT
      SGN = DBLE(3-2*JT)
      P(I,3) = GA*(BA*EA1+SGN*PA)
      P(I,4) = GA*(EA1+BA*SGN*PA)
      P(I+1,3) = PTOT(3)-P(I,3)
      P(I+1,4) = PTOT(4)-P(I,4)

c     mark as final hadrons
      ZLIST(I) = ZERO
      ZLIST(I+1) = ZERO

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
            LRNK(NP+N2) = LRNK(J)
            ZLIST (NP+N2) = ZLIST(J)
            nporig(NP+N2) = nporig(J)
            niorig(NP+N2) = niorig(J)
            nforig(NP+N2) = 0
            DO K=1,5
               P(NP+N2,K)=P(J,K)
            ENDDO
         ELSE
            N1= N1+1
            IF (N1.LT.J)   THEN
               LLIST(N1) = LLIST(J)
               LRNK(N1) = LRNK(J)
               ZLIST(N1) = ZLIST(J)
               nporig(N1) = nporig(J)
               niorig(N1) = niorig(J)
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
         LRNK(JJ) = LRNK(J)
         ZLIST(JJ) = ZLIST(J)
         nporig(JJ) = nporig(J)
         niorig(JJ) = niorig(J)
         nforig(JJ) = nforig(J)
         DO K=1,5
            P(JJ,K) = P(J,K)
         ENDDO
      ENDDO
      ENDIF

      if(Ndebug.gt.3)
     &  WRITE(LUN,*)' STRING_FRAG_4FLV: NP after fragmentation:',NP

      END



      SUBROUTINE GG_FRAG_4FLV (E0)
C...This routine fragments a  gluon-gluon system
C.  of mass E0 (GeV)
C.  the particles produced are in the  jet-jet frame
C.  oriented along the z axis
C...........................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_plist_cmmn.inc'
      INCLUDE 'sib_zlist_cmmn.inc'
      INCLUDE 'sib_mass_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      INCLUDE 'sib_parto_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'

c$$$      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
c$$$      COMMON /S_PLIST/ P(8000,5), LLIST(8000), NP
c$$$      COMMON /S_MASS1/ AM(99), AM2(99)
c$$$      COMMON /S_PARTO/ NFORIG(8000), NPORIG(8000), IPFLAG, NINT
c$$$      COMMON /S_ZLIST/ ZLIST(8000)
c$$$      DATA ZERO,HALF,ONE,TWO,THREE,FOUR /0.D0,0.5D0,1.D0,2.D0,3.D0,4.D0/
c$$$      DATA EPS3,EPS5,EPS8,EPS10 /1.D-3,1.D-5,1.D-8,1.D-10/

      DIMENSION WW(2,2),PTOT(4),PX(3),PY(3),IFL(3),PMQ(3)

      if(Ndebug.gt.3) then
        WRITE(LUN,*)
     &    ' GG_FRAG_4FLV: called with (E0)',
     &    E0
        WRITE(LUN,*)' GG_FRAG_4FLV: NP before fragmentation:',NP
      endif

C...  'leading' strange fraction
      PAR2_def = PAR(2)
      IF(IPAR(39).eq.2) PAR(2) = PAR(66)

C...Generate the 'forward' leading particle.
100   I = NP+1
      I0 = -1 + TWO*INT((TWO-EPS8)*S_RNDM(0))
c     dummy rank argument
      IDM = 5
      CALL SIB_I4FLAV(I0,0,IDM,IFL1, LDUM)
      CALL SIB_I4FLAV(IFL1,0,IDM,IFL2, LLIST(I))
      CALL PTDIS_4FLV(IFL1,PX1,PY1)
      CALL PTDIS_4FLV(IFL2,PX2,PY2)
      P(I,1) = PX1+PX2
      P(I,2) = PY1+PY2
      P(I,5) = AM(IABS(LLIST(I)))
      XM1 = P(I,5)**2+P(I,1)**2+P(I,2)**2
      Z1 = ZDIS_4FLV (IFL1,1,0.25D0*XM1)
      Z2 = ZDIS_4FLV (IFL2,1,0.25D0*XM1)
      T1  = FOUR*XM1/(E0*E0*(Z1+Z2))
      P(I,4) = 0.25D0*E0*(Z1+Z2 + T1)
      P(I,3) = 0.25D0*E0*(Z1+Z2 - T1)

      nforig(I)= 0
      nporig(I)= Ipflag*3 + Nint
      niorig(I)= iiflag
      ZLIST(I) = Z1 + Z2

C...Generate the 'backward' leading particle.
      I = I+1
      CALL SIB_I4FLAV(-I0,0,IDM,IFL3, LDUM)
      CALL SIB_I4FLAV(IFL3,0,IDM,IFL4, LLIST(I))
      CALL PTDIS_4FLV(IFL3,PX3,PY3)
      CALL PTDIS_4FLV(IFL4,PX4,PY4)
      P(I,1) = PX3+PX4
      P(I,2) = PY3+PY4
      P(I,5) = AM(IABS(LLIST(I)))
      XM2 = P(I,5)**2+P(I,1)**2+P(I,2)**2
      Z3 = ZDIS_4FLV (IFL3,1,0.25D0*XM2)
      Z4 = ZDIS_4FLV (IFL4,1,0.25D0*XM2)
      T2  = FOUR*XM2/(E0*E0*(Z3+Z4))
      P(I,4) = 0.25D0*E0*( Z3+Z4 + T2)
      P(I,3) = 0.25D0*E0*(-Z3-Z4 + T2)

      nforig(I)= 0
      nporig(I)= Ipflag*3 + Nint
      niorig(I)= iiflag
      ZLIST(I) = Z3 + Z4
c      PAR24def = PAR(24)
c     charm QCD fusion
c      IF(IPAR(17).eq.2) PAR(24) = 0.

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
         WW(1,1) = HALF * (1 - Z1 - HALF*T2) 
         WW(2,1) = HALF * (1 - Z3 - HALF*T1)
         PX(1) = -PX1
         PY(1) = -PY1
         PX(2) = -PX3
         PY(2) = -PY3
         IFL(1) = -IFL1
         IFL(2) = -IFL3
      ELSE
         WW(1,1) = HALF * (1 - Z2 - HALF*T2) 
         WW(2,1) = HALF * (1 - Z4 - HALF*T1)
         PX(1) = -PX2
         PY(1) = -PY2
         PX(2) = -PX4
         PY(2) = -PY4
         IFL(1) = -IFL2
         IFL(2) = -IFL4
      ENDIF
      PX(3) = ZERO
      PY(3) = ZERO
      PTOT (1) = PX(1)+PX(2)
      PTOT (2) = PY(1)+PY(2)
      PTOT (3) = HALF*E0*(WW(1,1)-WW(2,1))
      PTOT (4) = HALF*E0*(WW(1,1)+WW(2,1))

      PMQ(1) = QMASS(IFL(1))
      PMQ(2) = QMASS(IFL(2))

C...produce new particle: side, pT
300      I=I+1
      if(i.gt.8000) then
        write(6,'(1x,a,i8)') 
     &    'GG_FRAG: no space left in S_PLIST:',I
        call sib_reject   
      endif
      nforig(I)= 0
      nporig(I)= Ipflag*2 + Nint
      niorig(I)= iiflag

      JT=1.5D0+S_RNDM(0)
      JR=3-JT
c      CALL PTDIS (IFL(JT), PX(3),PY(3))

C...particle ID
      CALL SIB_I4FLAV (IFL(JT), 0, IDM, IFL(3), LLIST(I))
      PMQ(3) = QMASS(IFL(3))
      P(I,5) = AM(IABS(LLIST(I)))

      CALL PTDIS_4FLV (IFL(3), PX(3),PY(3))
      
C...test end of fragmentation
      WREM2 = PTOT(4)**2-PTOT(1)**2-PTOT(2)**2-PTOT(3)**2
      IF (WREM2 .LT. 0.1D0)  GOTO 200
      WMIN = PMQ(1)+PMQ(2)+TWO*PMQ(3)+1.1D0+(TWO*S_RNDM(0)-ONE)*0.2D0
      IF (WREM2 .LT. WMIN**2)THEN
         GOTO 400
      ENDIF

C...fill transverse momentum
      P(I,1) = PX(JT) + PX(3)
      P(I,2) = PY(JT) + PY(3)

C...Choose z
      XMT2 = P(I,5)**2+P(I,1)**2+P(I,2)**2
      Z = ZDIS_4FLV (ifl(3),IFL(JT), XMT2)

      ZLIST(I) = Z      
      WW(JT,2) = Z*WW(JT,1)
      WW(JR,2) = XMT2/(WW(JT,2)*E0**2)

      P(I,3) = WW(1,2)*HALF*E0 - WW(2,2)*HALF*E0
      P(I,4) = WW(1,2)*HALF*E0 + WW(2,2)*HALF*E0

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

      CALL SIB_I4FLAV (IFL(JR), -IFL(3), IDM, IFLA, LLIST(I+1))
      P(I+1,5) = AM(IABS(LLIST(I+1)))
      P(I,1)   = PX(JT)+PX(3)      
      P(I,2)   = PY(JT)+PY(3)      
      nporig(I)= Ipflag*2 + Nint
      niorig(I)= iiflag
      I1 = I+1
      nporig(I1)= Ipflag*2 + Nint
      niorig(I1)= iiflag
      P(I1,1) = PX(JR)-PX(3)      
      P(I1,2) = PY(JR)-PY(3)      
      XM1 = P(I,5)**2+P(I,1)**2+P(I,2)**2
      XM2 = P(I1,5)**2+P(I1,1)**2+P(I1,2)**2
      IF (dSQRT(XM1)+dSQRT(XM2) .GT. dSQRT(WREM2)) GOTO 200
      if (ptot(4).le.ZERO) goto 200
      PT2 = (P(I,1)+P(I1,1))**2+(P(I,2)+P(I1,2))**2
      WREMPT = dsqrt(WREM2+PT2)
      EA1 = (WREM2+XM1-XM2+PT2)/(TWO*WREMPT)
      PA2 = (EA1**2-XM1)
      if (PA2.ge.ZERO) then
        PA = dSQRT(PA2)
      else
         goto 200
      endif
      BA = PTOT(3)/PTOT(4)
      GA = PTOT(4)/WREMPT
      SGN = DBLE(3-2*JT)
      P(I,3) = GA*(BA*EA1+SGN*PA)
      P(I,4) = GA*(EA1+BA*SGN*PA)
      P(I+1,3) = PTOT(3)-P(I,3)
      P(I+1,4) = PTOT(4)-P(I,4)
      ZLIST(I) = ZERO
      ZLIST(I+1) = ZERO
      N0 = I-NP-1
      ENDDO                  ! loop on two `remaining strings'
      NP = I+1
c      PAR(24) = PAR24def 
      IF(Ndebug.gt.3) then
        WRITE(LUN,*)' GG_FRAG_4FLV: NP after fragmentation:',NP
      ENDIF
      RETURN
      END

