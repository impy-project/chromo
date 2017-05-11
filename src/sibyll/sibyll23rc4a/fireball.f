      SUBROUTINE FIREBALL_4FLV(L0,P0,PCHEXin,IREJ)
C-----------------------------------------------------------------------
C... "decay" of an excited state with the quantum numbers
C.   of particle L0 and the 5-momentum P0
C.   4 flavor generalization /FR'13
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_plist_cmmn.inc'
      INCLUDE 'sib_nw_prm.inc'
      INCLUDE 'sib_run_cmmn.inc'
      INCLUDE 'sib_mass_cmmn.inc'
      INCLUDE 'sib_chp_cmmn.inc'
      INCLUDE 'sib_cnam_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      INCLUDE 'sib_parto_cmmn.inc'

      DIMENSION P0(5), LL(10), PD(10,5), IFL(3), INONLEAD(2)
      DIMENSION LRESCHEX(6:99), LRES(6:99), LCON(6:99), LPIC(-1:1)

c     charge exchange map
      DATA LCON(6:33) /7,6,6,22,21,9,9,14,13,4*0,20,19,10,9,23,24,27,27,
     &     25,31,30,29,28,32,33/
      DATA LCON(34:49) /35,34,35,38,37,39,41,42,41,42,45,44,45,48,47,49/
      DATA LCON(50:83) /0,52,51,54,53,4*0,71,72,10*0,
     &     59,60,73,74,75,76,77,80,81,78,79,0,83/
      DATA LCON(84:99) /84,85,86,87,88,89,4*0,94,95,96,97,98,99/
c     pion charge conversion map
      DATA LPIC /8,6,7/
c     charge exchange to resonances map
      DATA LRESCHEX(6:33) /26,27,27,30,31,9,9,42,41,19*0/
      DATA LRESCHEX(34:39) /45,44,45,48,47,39/ 
      DATA LRESCHEX(40:49) /41,42,43,42,45,46,45,48,47,49/
      DATA LRESCHEX(50:83) /0,52,51,54,53,4*0,60,59,10*0,71,72,73,75,74,
     &     77,76,79,78,80,81,0,83/
      DATA LRESCHEX(84:99) /84,85,86,87,88,89,4*0,94,95,96,97,98,99/
c     resonance excitation map
      DATA LRES(6:39) /27,25,26,28,29,9,9,41,42,19*0,44,45,46,47,48,39/
      DATA LRES(40:49) /40,41,42,43,44,45,46,47,48,49/
      DATA LRES(50:83) /0,51,52,53,54,4*0,78,79,10*0,71,72,73,76,77,76,
     &     77,78,79,80,81,0,83/
      DATA LRES(84:99) /94,95,96,97,98,89,4*0,94,95,96,97,98,99/
      INCLUDE 'sib_utl_cmmn.inc'

c...  charge exchange reaction rate
c      DATA PCHEX /0.33/
c     default parameter: PAR(61)
      PCHEX = PCHEXin

c...  suppression of high mass particles in fireball
c     xmpsuppr = prob. accepting additional proton
      XMPSUPPR=PAR(33)
      IF(ABS(XMPSUPPR).lt.EPS3) THEN
         WRITE(LUN,*)
     &        'Error: too low mass suppression in 4 flv fireball!'
         WRITE(LUN,*)
     &        'Probably PAR(33)/IPAR(14) not properly set, aborting..'
         STOP
      ENDIF
      XTEMPH=(AM(6)-AM(13))/dLOG(XMPSUPPR)

      IF(Ndebug.gt.3) THEN
         WRITE(LUN,*)' FIRBALL_4FLV: called with (L0,P0):',
     &        L0,P0
         WRITE(LUN,*)' 2nd Proton rejection prob.:',XMPSUPPR
         WRITE(LUN,*)' fireball temperature:',XTEMPH
         WRITE(LUN,*)' charge exchange prob.:',PCHEX
         WRITE(LUN,*)' multiplicity width:',PAR(38)
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
         xa=-1.D0
         xa = dlog(xa)
         STOP
c         RETURN
      ENDIF

      LA = ABS(L0)
      ISGN = ISIGN(1,L0)
      DELTAE = P0(5) - AM(LA)
      IF(DELTAE.lt.AM(6)+0.02D0)THEN
         IREJ = 1
         IF(ndebug.gt.3)
     &    WRITE(LUN,*)'FIRBALL_4FLV:  too low mass!! aborting...',IREJ
c         xa=-1.
c         xa=log(xa)
c         stop        
         RETURN
      ENDIF
      AV = TWO*SQRT(DELTAE)

c...  select number of particles in fireball
c     at least two
 200  XRNDM = GASDEV(LA)
      NPI = INT(AV*(ONE+PAR(38)*XRNDM))
      XMMIN = AM(LA)+(NPI-1)*AM(6)+0.02D0
      IF(Ndebug.gt.3)
     &     WRITE(LUN,*)'NPI,av,rndm,xmin,delta',
     &     NPI,av,XRNDM,xmmin,P0(5)-XMMIN

      IF((NPI.LE.1).OR.(NPI.GT.9).OR.(P0(5).LT.XMMIN))THEN
         GOTO 200
      ENDIF
      IF(Ndebug.gt.3) 
     &  WRITE(LUN,*)' FIRBALL_4FLV: No. of particles sampled. ',
     &  '(NPI,DELTAE,NTRY):',NPI,DELTAE,NTRY

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
            IFL(1)=IFL1
            IFL(2)=IFL2
         ELSE
            LL(I)=LRES(LA)
            CALL HSPLI(L0,IFL1,IFL2)
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
            if(la.eq.42) l1 = l1 + 2 * int(TWO*S_RNDM(LA))
            LL(I)=L1*ISGN
c            WRITE(LUN,*)'charge exchange!',ISGN*LA,'->',L1
         ELSE
            L1=LA
            LL(I)=LA*ISGN
         ENDIF
c     determine remaining charge
         IDQ=ICHP(LA)*ISGN-ICHP(L1)*ISIGN(1,LL(I))
         IF(ABS(IDQ).gt.1) write(lun,*) 'LA,L1',LA,L1
         LL(I+1)=LPIC(IDQ)      ! compensate with meson
         IF(NPI.eq.2) GOTO 300
c     split last hadron again to start hadron chain
 211     CALL HSPLI (LL(I+1),IFL(1),IFL(2))

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
      JT=INT(ONEHALF+S_RNDM(I))
      JR=3-JT
      NTRYS=0
      IFLB=IFL(JT)
      IDM = 5
 240  CALL SIB_I4FLAV (IFL(JT), 0, IDM, IFL(3), LL(I))
      IF(NTRYS.gt.50) GOTO 210    
      NTRYS=NTRYs+1
      W=dEXP(-AM(ABS(LL(I)))/XTEMPH)
      IF(Ndebug.gt.4) 
     &  WRITE(LUN,*)' FIRBALL_4FLV: flavor added: ',
     &  '(I,NTRYS,LL(I),IFL3,W):',I,NTRYS,LL(I),IFL(3),W
      IF(W.LT.S_RNDM(I).and.INONLEAD(JT).eq.1) GOTO 240

c...  kinematic limits...     
      WREM = WREM-AM(ABS(LL(I)))
      WREM2_2=WREM2+2*dSQRT(WREM2)*AM(ABS(LL(I)))+AM2(ABS(LL(I)))
      IF(Ndebug.gt.4) 
     &  WRITE(LUN,*)' FIRBALL_4FLV: kinematic limits: ',
     &  '(I,NTRYS,P05**2,WREM2):',I,NTRYS,P0(5)**2,WREM2_2
      IF(WREM2_2+0.2D0*S_RNDM(I).ge.P0(5)**2) GOTO 240
      WREM2=WREM2_2
      IF(Ndebug.gt.3) 
     & WRITE(LUN,*)
     & ' FIRBALL_4FLV: Hadron added: (KF,NAMP,I,NONlead,WRME2)',
     & LL(I),NAMP(ABS(LL(I))),I,INONLEAD(JT),WREM2

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
 250  CALL SIB_I4FLAV (IFL(JT), IFL(JR), IDM, IFL(3), LL(I))
      IF(NTRYC.gt.10) GOTO 210
      NTRYC=NTRYC+1
      WREM2_2=WREM2+2*dSQRT(WREM2)*AM(ABS(LL(I)))+AM2(ABS(LL(I)))
      IF(Ndebug.gt.5) 
     & WRITE(LUN,*)' FIRBALL_4FLV: closing List: (IFL1,IFL2,KF,',
     &'NAMP,I,NTRYC,WREM2)',
     & IFL(JT),IFL(JR),LL(I),NAMP(ABS(LL(I))),I,NTRYC,WREM2_2

      IF(WREM2_2+0.2D0*S_RNDM(I).ge.P0(5)**2) GOTO 250

 300  IF(Ndebug.gt.3) 
     &     WRITE(LUN,*)
     &     ' FIRBALL_4FLV: flavors sampled. (NPI,LL,WREM,NTRYL):',
     &     NPI,(LL(ii),ii=1,NPI),WREM,NTRYL

c...  fill phasespace
      CALL DECPAR (0,P0,NPI,LL,PD)
      DO J=1,NPI
         NP = NP+1
         LLIST(NP) = LL(J)
         NPORIG(NP) = IPFLAG*2
         niorig(NP)= iiflag
         DO K=1,5
            P(NP,K) = PD(J,K)
         ENDDO
      ENDDO
      PAR(5)=PAR5def
      IREJ = 0
      RETURN
      END
