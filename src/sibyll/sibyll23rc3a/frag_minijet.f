
      SUBROUTINE FRAG_MINIJET(IDX,IBAD)
C-----------------------------------------------------------------------
C     routine that fragments a gluon - gluon system \FR'14
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER IDX,IBAD
           
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_plist_cmmn.inc'
      include 'sib_nw_prm.inc'
      INCLUDE 'sib_run_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      INCLUDE 'sib_chist_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'
      
      DOUBLE PRECISION PGG,PST,PBM,PTG,E0,PT2JET,PTJET,TH,FI,S_RNDM,
     &     PAR1_def,PAR24_def,PAR3_def,PAR2_1_def,PAR2_2_def,PAR5_def,
     &     PAR6_def,XM,QMASS,DBETJ
      DIMENSION PST(5),PBM(5),PTG(5)
      INTEGER IST,ITGST,IBMST,IPID,IFLB,IFLT,NOLD,IS,IFL1,IFBAD
      DATA PGG /1.D0/

C     read partons from stack
c     references are string --> bm-parton --> tg-parton
c     read string 4momentum from stack
      call rd_prtn_4vec(IDX,PST,IPID,IBMST)
      call rd_prtn_4vec(IBMST,PBM,IFLB,ITGST)
      call rd_prtn_4vec(ITGST,PTG,IFLT,IST)
      IF(IDX.ne.IST) then
         write(lun,*) 'FRAG_minijet: reference loop broken!' , IDX
         call sib_reject
      endif

C..   kinematic variables
      E0 = PST(5)            ! string mass
      PT2JET = PBM(1)**2 + PBM(2)**2
      PTJET = sqrt(PT2JET)
      TH = ASIN(MIN((ONE-EPS8),TWO*PTJET/E0))
c      FI = ASIN(MIN((ONE-EPS8),PBM(2)/PTJET))
      FI = TWOPI*S_RNDM(IDX)
c      TH = PST(1)
c      FI = PST(2)

      IF(NDEBUG.gt.1) WRITE(LUN,*)' FRAG_MINIJET: IDX,EE,IFLB,IFLT,PT',
     &     IDX,E0,IFLB,IFLT,PTJET,IBAD
      IF(NDEBUG.gt.1) WRITE(LUN,*)' FRAG_MINIJET: PTJET,TH,FI:',
     &     PTJET,TH,FI

C...  parameter setup

c     baryon production setup
      PAR1_def = PAR(1)
      if( NSOF+NJET.gt.0) then
         PAR(1)= PAR(15)
      else
         PAR(1)= PAR(14)
      endif
            
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
         WRITE(LUN,*)' FRAG_minijet: charm model not implemented!'
         call sib_reject
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
      IF ( (E0.LT.EIGHT) .OR. (S_RNDM(0).GT.PGG)) THEN
C...  one string case, q - qbar
C     set 'leading' strange fraction
         PAR2_2_def = PAR(2)
         IF(IPAR(39).eq.2) PAR(2) = PAR(66)
         IS = -1 + 2*INT((TWO-EPS8)*S_RNDM(0))
 100     IFL1 = IS*(INT((TWO+PAR(2))*S_RNDM(0))+1)
         XM = TWO*QMASS(IFL1)+0.3D0
         if(E0.LE.XM) GOTO 100
         PAR(2) = PAR2_2_def
         IF(IABS(IFL1).eq.3)THEN
            IF(S_RNDM(IFL1).lt.PAR(24)*PAR(125))IFL1 = IS*4
            XM = TWO*QMASS(IFL1)+0.3D0
            if(E0.LE.XM) GOTO 100
         ENDIF
c     charm QCD fusion, soft+hard
c            IF(IPAR(17).eq.2) PAR(24) = 0.
         CALL STRING_FRAG_4FLV 
     &        (E0,IFL1,-IFL1,ZERO,ZERO,ZERO,ZERO,IFBAD,0)
         if(IFBAD.gt.0) then
            IF(ndebug.gt.1)
     &       WRITE(LUN,*)
     &           ' JET_FRAG: rejection in STRING_FRAG (IFL,E0,NCALL):',
     &           IFL1,E0,NCALL
            PAR(24) = PAR24_def
            PAR(1) = PAR1_def
            PAR(2) = PAR2_1_def
            PAR(5) = PAR5_def
            PAR(6) = PAR6_def
            PAR(3) = PAR3_def 
            RETURN
         ENDIF
      ELSE
C...  two string case, gluon - gluon
         CALL GG_FRAG_4FLV(E0)
      ENDIF

c      DBETJ = (DX1J-DX2J)/(DX1J+DX2J)
      DBETJ = PST(3)/PST(4)
      CALL SIROBO (NOLD+1,NP,TH,FI,ZERO,ZERO,DBETJ)

      if(Ndebug.gt.1) WRITE(LUN,*)
     &     ' JET_FRAG: particles produced:',NP-NOLD
      PAR(24) = PAR24_def
      PAR(1) = PAR1_def
      PAR(2) = PAR2_1_def
      PAR(5) = PAR5_def
      PAR(6) = PAR6_def
      PAR(3) = PAR3_def 
      IBAD = 0
      END

