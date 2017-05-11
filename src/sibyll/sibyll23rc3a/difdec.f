
      SUBROUTINE DIFDEC (L0, Irec, IBAD, P0)
C-----------------------------------------------------------------------
C..."decay" of an excited state with the quantum numbers
C.   of particle L0 and the 5-momentum P0
C.   - low energy: phase space decay (fire ball model)
C.   - intermediate energy: one-string decay (longitudinal phase space)
C.   - high energy: pomeron-hadron scattering (multi-string model) 
C-----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE

c     external types
      INTEGER L0, Irec, IBAD
      DOUBLE PRECISION P0
      DIMENSION P0(5)
      
      include 'sib_nw_prm.inc'
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_plist_cmmn.inc'
      INCLUDE 'sib_parto_cmmn.inc'
      INCLUDE 'sib_chp_cmmn.inc'
      INCLUDE 'sib_run_cmmn.inc'
      INCLUDE 'sib_mass_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      
      INCLUDE 'sib_utl_cmmn.inc'

c     internal types
      INTEGER LL,LCON,LRES,LRES1,NTRYS,NRJECT,LA,N1,IREJ,I,J,IFLA,
     &     IFL1,IFL2,IFBAD,NPI,IRES,LA1,JQQ,JQTOT,K,JQR
      DOUBLE PRECISION PD,BE,EMIN,EMIN2,PCHEX,PRES,DELTAE,SQS_0,FERMI,
     &     PAR1_def,PAR24_def,PAR53_def,GA,BEP,S_RNDM,AV,GASDEV,PCXG,
     &     XI1,XI2,XSMR
      DIMENSION LL(10), PD(10,5), BE(3), LCON(6:39), LRES(6:39),
     +     LRES1(6:39)
      DATA EMIN /0.7D0/
      DATA EMIN2 /10.D0/

      DATA LCON /7,6,6,11,11,9,9,14,13,19*0,35,34,35,38,37,39/
      DATA LRES /26,27,27,11,11,9,9,14,13,19*0,35,34,35,38,37,39/ 
      DATA LRES1 /27,25,26,11,11,9,9,14,13,19*0,35,34,35,38,37,39/
      DATA PCHEX /0.33D0/            ! probability of charge exchange
      DATA PRES /0.7D0/         ! probability of forming a resonance
      DATA NRJECT /0/

      IF(NDEBUG.gt.2)
     &     WRITE(LUN,'(2X,A,1x,I2,1x,I2,/,5(2x,F8.3))')
     &     'DIFDEC: (L0,Irec,P0):',L0,Irec,(P0(i),i=1,5)
      
      
      NTRYS = 0

      LA = IABS(L0)
      DELTAE = P0(5) - AM(LA)
      IF(IBAR(LA).ne.0.or.IPAR(65).eq.0)THEN
c     baryons
         EMIN = PAR(30)
      ELSE
c     mesons
         EMIN = PAR(112)
      ENDIF
c      IBAD = 0
      PAR1_def= PAR(1)
      if(Irec.gt.0) PAR(1)= PAR(16)
c      XSMR = HALF
c     XI2=FERMI(DELTAE,EMIN2,XSMR)
c     XI1=FERMI(DELTAE,EMIN,XSMR)
      XSMR=PAR(131)*EMIN
      XI1=MAX((EMIN-DELTAE)/XSMR,ZERO)      
      XSMR=PAR(131)*EMIN2
      XI2=MAX((EMIN2-DELTAE)/XSMR,ZERO)
      if(Ndebug.gt.2) 
     &     WRITE(LUN,'(A29,2(2x,F5.2),2(2x,F8.3))')
     &     '  DIFDEC: EMIN1,EMIN2,XI1,XI2',
     &     EMIN,EMIN2,Xi1,Xi2
      
C...  pomeron-hadron scattering (pi0 is used instead of pomeron)      
      IF ((IPAR(10).gt.0).and.(Irec.gt.0).and.
     &     (DELTAE.gt.EMIN2.or.S_RNDM(LA).gt.XI2))  THEN
         if(Ndebug.gt.2) 
     &        WRITE(LUN,*)'  DIFDEC: central (L0,DELTAE,NP,XI):',
     &        L0,DELTAE,NP,XI2
         N1 = NP+1
         if(irec.gt.0.and.ipar(5).eq.1) par(1)= par(15)
 50      CONTINUE
         IPFLAG= IPFLAG*100
c     create subevent
         SQS_0 = SQS
         CALL INI_EVENT(P0(5),L0,6,0)
c     create L0 - pi0 interaction, pi0(pid=6) target
         CALL SIB_NDIFF(L0, 1, P0(5), 0, IREJ) ! ori
c     restore main event
         SQS = SQS_0         
         IF(IREJ.NE.0) THEN
            NP = N1-1
            GOTO 50
         ENDIF
         PAR(1) = PAR1_def
         DO J=1,3
            BE(J)=P0(J)/P0(4)
         ENDDO
         GA=P0(4)/P0(5)
         if(P0(3).lt.ZERO) then
           do i=N1,NP
             P(I,3) = -P(I,3)
           enddo
         endif
         DO I=N1,NP
            BEP=BE(1)*P(I,1)+BE(2)*P(I,2)+BE(3)*P(I,3)
            DO J=1,3
               P(I,J)=P(I,J)+GA*(GA*BEP/(ONE+GA)+P(I,4))*BE(J)
            ENDDO
            P(I,4)=GA*(P(I,4)+BEP)
         ENDDO

C..."string-like" decay
      ELSE IF (DELTAE .GT. EMIN .or. S_RNDM(LA).gt.XI1)  THEN            
         IF(NDEBUG.gt.2) 
     &        WRITE(LUN,'(2X,A,3(2x,F8.3))')
     &        'DIFDEC: string-like, (DELTAE,E0,central prob.):',
     &        DELTAE,P0(5),ONE-XI2
c     set charge exchange probability, i.e. prob for p* -> n + pip
         PAR53_def = PAR(53)
         PAR(53) = PAR(130)
         N1 = NP+1
         CALL HSPLI(L0,IFL1,IFL2)
         PAR(53) = PAR53_def
         IF (P0(3) .GT. ZERO.and.L0.gt.0)  THEN
            IFLA = IFL2
            IFL2 = IFL1
            IFL1 = IFLA
         ENDIF
c     randomize flavor orientation in string
         IF(IPAR(25).eq.1.and.S_RNDM(L0).gt.PAR(39))THEN
            IFLA = IFL2
            IFL2 = IFL1
            IFL1 = IFLA
         ENDIF
         PAR24_def = PAR(24)
         SELECT CASE(IPAR(15))
         CASE(2,4,5,6,8,9,10,11)
            PAR(24) = PAR(25)*dEXP(-PAR(26)/P0(5))
         CASE(7)
            PAR(24) = PAR(25)
         END SELECT
 10      CONTINUE
         IPFLAG = IPFLAG*10
         CALL STRING_FRAG_4FLV 
     +        (P0(5), IFL1, IFL2, ZERO,ZERO,ZERO,ZERO,IFBAD,-1)
         IF (IFBAD .EQ. 1)then
            if(ndebug.gt.1)
     &           WRITE(lun,*)' SIB_DIFF: string-frag rejection! ',
     &           '(M,NCALL)',P0(5),NCALL
            NTRYs = NTRYs + 1
            NP = N1-1
            IFBAD = 0
            IF(NTRYs.gt.5)then ! resample diff. mass              
               NP = 0
               IBAD = 1
               PAR(24) = PAR24_def
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
               P(I,J)=P(I,J)+GA*(GA*BEP/(ONE+GA)+P(I,4))*BE(J)
            ENDDO
            P(I,4)=GA*(P(I,4)+BEP)
         ENDDO
         PAR(24) = PAR24_def

C...Phase space decay of the excited state
      ELSEIF(DELTAE.GT.AM(7)+0.02D0)THEN
         if(Ndebug.gt.2) 
     &        WRITE(LUN,*)'  DIFDEC: fireball, (DELTAE,string prob.):',
     &        DELTAE,ONE-XI1
         IF(IPAR(14).GT.0.and.IPAR(14).NE.7)THEN
            IF(IPAR(14).eq.5) PCHEX = ZERO
            NPI=0
            IRES = 0
            IF (S_RNDM(0).LT.PRES) THEN
               IF (LA.LT.9) THEN
c     if kinematically possible produce rho0 in charge exchange
                  LL(1) = LRES(LA)
                  DELTAE = P0(5) -  AM(LRES(LA))
                  IF (DELTAE.GT.AM(7)+0.02D0) GOTO 100
               ENDIF
            ENDIF
c     switch charge exchange on/off
            IF( S_RNDM(0).LT.PCHEX)THEN
               LL(1) = LCON(LA)*ISIGN(1,L0)
               IF( (L0 .EQ. 6) .OR. (L0 .EQ. 11) )
     .              LL(1) = LL(1)+INT((TWO-EPS8)*S_RNDM(0))
            ELSE
               LL(1) = L0
            ENDIF
            
            DELTAE = P0(5) -  AM(LA)
 100        AV = TWO*dSQRT(DELTAE)
            LA1 = IABS(LL(1))
            NPI = INT(AV*(TWO+HALF*GASDEV(LA)))
            IF (IPAR(14).EQ.6)THEN
               IF(NPI.LT.1.OR.NPI.GT.9.OR.AM(LA1)+NPI*AM(7)+0.02D0
     .              .GT.P0(5))  GOTO 100
            ELSE
               IF(NPI.LT.0.OR.NPI.GT.9.OR.AM(LA1)+NPI*AM(7)+0.02D0
     .              .GT.P0(5))  GOTO 100
            ENDIF
c     create resonances inside fireball..
            IF(IPAR(14).ge.2
     +           .and.DELTAE.GE.AM(LA1)+AM(27)+(NPI-1)*AM(7)+0.02D0)
     +           IRES = 1
            IF(IPAR(14).ge.3.and.DELTAE.GE.AM(LA1)+NPI*AM(27)+0.02D0) 
     +           IRES=3
            JQQ = ICHP(LA)*ISIGN(1,L0)-
     .           ICHP(IABS(LL(1)))*ISIGN(1,LL(1))  
 120        JQTOT = 0
            DO K=2,NPI
               LL(K) = 6+INT(S_RNDM(0)*(THREE-EPS8))
c     suppress pi0 in fireball
               IF(IPAR(14).ge.4)
     +              LL(K) = 7+INT(S_RNDM(0)*(TWO-EPS8))
c     IF(IRES.EQ.1.and.S_RNDM(LA).LT.0.5)
               IF(IRES.EQ.1) THEN
                  LL(K) = 27-INT(S_RNDM(0)*(THREE-EPS8))
                  IRES = 2
               ENDIF
               IF(IRES.EQ.3)
     +              LL(K) = 27-INT(S_RNDM(0)*(THREE-EPS8))
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
               niorig(NP)= iiflag
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
            IF( DELTAE.LT.AM(7)+0.02D0) GOTO 222
            IF( S_RNDM(0).LT.PAR(31))THEN
               LL(1) = LRES1(LCON(LA))
               IF( (L0 .EQ. 6) .OR. (L0 .EQ. 11) )
     .              LL(1) = LRES1(IABS(L0)+INT((TWO-EPS8)*S_RNDM(0)))
            ENDIF
 300        AV = TWO*dSQRT(DELTAE)
            LA1 = IABS(LL(1))
            NPI = INT(AV*(TWO+HALF*GASDEV(LA)))
            IF(ABS(PAR(32)).gt.ZERO) 
     &           NPI = INT(AV*(PAR(32)+HALF*GASDEV(LA)))
            IF(NPI.LT.0.OR.NPI.GT.9.OR.AM(LA1)+NPI*AM(7)+0.02D0
     .           .GT.P0(5))  GOTO 300
c     create resonances inside fireball..
c            IRES=3
            JQQ = ICHP(LA)*ISIGN(1,L0)-
     .           ICHP(IABS(LL(1)))*ISIGN(1,LL(1))  
 320        JQTOT = 0
            DO K=2,NPI
               LL(K) = 6+INT(S_RNDM(0)*(THREE-EPS8))
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
               niorig(NP)= iiflag
               DO K=1,5
                  P(NP,K) = PD(J,K)
               ENDDO
            ENDDO

         ELSEIF (IPAR(14).LE.-1) THEN
C...  generalized fireball model
            IF(Ndebug.gt.2) 
     &           WRITE(LUN,*)'  DIFDEC: using generalized fireball!'
c     set charge exchange probability, 
c     i.e. prob for p* -> n + pip
            PCXG = PAR(61)
            CALL FIREBALL_4FLV(L0,P0,PCXG,IFBAD)
            IF(IFBAD.eq.1)THEN
               IF(ndebug.ne.0)THEN
                  IF(NRJECT.le.10)THEN
                     WRITE(LUN,*)
     &                    ' DIFDEC: warning: fireball rejection! ',
     &                    'diff. mass to low to dissociate beam!'
                     WRITE(LUN,*)
     &                ' DIFDEC: m_Beam, DELTAE ,AM(7)+0.02, NCALL: ', 
     &                AM(LA),DELTAE,'>',AM(7)+0.02D0,NCALL
                  ENDIF
                  IF(NRJECT.eq.10) 
     &             write(lun,*)' this was the last warning.. good luck!'
               ENDIF
               NRJECT = NRJECT + 1
               NP = 0
               IBAD = 1
               RETURN
            ENDIF

         ELSE
 222        IF(IPAR(14).EQ.7)  DELTAE = P0(5) - AM(LA)
            AV = TWO*dSQRT(DELTAE)
 200        NPI = INT(AV*(ONE+HALF*GASDEV(0)))
c            print *,'npi:',npi,'av',av,'p05',p0(5),am(la),deltae
            IF(NPI.LE.0.OR.NPI.GT.9.OR.AM(LA)+NPI*AM(7)+0.02D0
     .           .GT.P0(5))  GOTO 200
            IF (S_RNDM(0).LT.PCHEX)  THEN
               LL(NPI+1) = LCON(LA)*ISIGN(1,L0)
               IF( (L0 .EQ. 6) .OR. (L0 .EQ. 11) )
     .              LL(NPI+1) = LL(NPI+1)+INT((TWO-EPS8)*S_RNDM(0))
            ELSE
               LL(NPI+1) = L0
            ENDIF
            JQQ = ICHP(LA)*ISIGN(1,L0)-
     .           ICHP(IABS(LL(NPI+1)))*ISIGN(1,LL(NPI+1))  
 220        JQTOT = 0
            DO K=1,NPI-1
               LL(K) = 6+INT(S_RNDM(0)*(THREE-EPS8))
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
               niorig(NP)= iiflag
               DO K=1,5
                  P(NP,K) = PD(J,K)
               ENDDO
            ENDDO
         ENDIF
      ELSE
         IF(NRJECT.le.10)THEN
            WRITE(LUN,*) '  DIFDEC rejection! ',
     &           'diff. mass to low to dissociate beam!'
            WRITE(LUN,*) '  DIFDEC: LA, m_Beam, DELTAE, NCALL : ', 
     &           LA, AM(LA),DELTAE,'>',AM(7)+0.02D0,NCALL
            IF(Irec.ne.1) 
     &           WRITE(LUN,*) '   was recursive call! (ECM):',P0(5)
         ENDIF
         IF(NRJECT.eq.10) 
     &        write(lun,*)' this was the last warning.. good luck!'
         NRJECT = NRJECT + 1            
         NP = 0
         IBAD = 1
         RETURN
      ENDIF
      PAR(1) = PAR1_def
      END
