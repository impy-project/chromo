C-----------------------------------------------------------------------
C     parton level administration tools for SIBYLL                \FR'14
C-----------------------------------------------------------------------

C...  COMMON /S_PRTNS/ : parton stack
c     PP: 4momentum of parton, px,py,pz,energy,mass
c     LPID(1): parton id, i.e. flavor (u:1,d:2,s:3,c:4) for quarks
c     LPID(2): level of parton
c              fragmenting systems (strings,remnants) are marked as level0
c              partons that make up these systems are marked as level1
c     LPID(3): 'downward' reference 
c               pointer from level1 partons to their level0 parent
c     LPID(4): 'upward' reference 
c               pointer from level0 partons to their level-1 parent
c     LVL0IDX: index cache for level0 partons
c     NPP: total number of partons on stack
c     NPP0: number of level0 partons on stack


      SUBROUTINE ADD_PRTN(PX,PY,PZ,E,XMS,IPID,LVL,IREFin,IREFout)
C-------------------------------------------------------
C     routine to add a parton to the stack \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_prtns_cmmn.inc'
      PP(NPP+1,1) = PX
      PP(NPP+1,2) = PY
      PP(NPP+1,3) = PZ
      PP(NPP+1,4) = E
      PP(NPP+1,5) = XMS
      LPID(NPP+1,1) = IPID
      LPID(NPP+1,2) = LVL     
      LPID(NPP+1,3) = IREFin
      NPP = NPP + 1
c     level0 index
      IF(LVL.eq.0)THEN
         LVL0IDX(NPP0+1) = NPP
         NPP0 = NPP0 + 1
      ENDIF
      IREFout = NPP
      IF(NDEBUG.gt.6)THEN
         WRITE(LUN,*) ' ADD_PRTN: (#,PID,LEVEL,REF)',
     &        NPP,LPID(NPP,1),LPID(NPP,2),LPID(NPP,3)
         WRITE(LUN,*) '   4momentum:',(PP(NPP,JJ),JJ=1,5)
      ENDIF
      END


      SUBROUTINE ADD_PRTN_4VEC(PIN,IPID,LVL,IREFin,IREFout)
C----------------------------------------------------------
C     wrapper for ADD_PRTN to add 4momentum directly \FR'14
C----------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'sib_debug_cmmn.inc'
      DIMENSION PIN(5)
      call add_prtn
     &     (PIN(1),PIN(2),PIN(3),PIN(4),PIN(5),IPID,LVL,IREFin,IRF)
      IREFout = IRF
      END


      SUBROUTINE ADD_REF(IDX,Irefin)
C-------------------------------------------------------
C     routine to add a reference label to a particle
C     after it has been added to the stack        \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_prtns_cmmn.inc'
      
c      IF(LPID(IDX,3).ne.0)  WRITE(LUN,*)
c     &     'ADD_ref: warning particle already has defined reference,',
c     &     IDX,' overwritting..'
      IF(NDEBUG.gt.6)
     &WRITE(LUN,*) 'ADD_ref: (IDX,REFin)',IDX,Irefin
      LPID(IDX,3) = Irefin
      END


      SUBROUTINE RD_REF(IDX,Irefout)
C-------------------------------------------------------
C     routine to add a reference label to a particle
C     after it has been added to the stack        \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_prtns_cmmn.inc'

      Irefout = LPID(IDX,3)
      IF(NDEBUG.gt.6)
     &WRITE(LUN,*) 'RD_ref: (IDX,REFout)',IDX,Irefout
      END


      SUBROUTINE ADD_INT_REF(IDX,Irefin)
C-------------------------------------------------------
C     routine to add a reference label to an interaction
C     after it has been added to the stack        \FR'15
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_prtns_cmmn.inc'

      IF(NDEBUG.gt.6)
     &WRITE(LUN,*) 'ADD_int_ref: (IDX,REFin)',IDX,Irefin
      LPID(IDX,4) = Irefin
      END


      SUBROUTINE RD_INT(IDX,Irefout,Iout)
C-------------------------------------------------------
C     routine to add a reference label to an interaction
C     after it has been added to the stack        \FR'15
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_prtns_cmmn.inc'

      Irefout = LPID(IDX,4)
      IF(Irefout.ne.0) Iout = LPID(Irefout,1)
      IF(NDEBUG.gt.6)
     &WRITE(LUN,*) 'RD_int: (IDX,REFout,Iint)',IDX,Irefout,Iout
      END


      SUBROUTINE EDT_PRTN(IDX,PX,PY,PZ,EN,XMS,IREFout)
C-------------------------------------------------------
C     routine to edit a parton already on stack   \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_prtns_cmmn.inc'
      IF(NDEBUG.gt.6)THEN
         WRITE(LUN,*) ' EDT_PRTN: (#,PID,LEVEL,REF)',
     &        IDX,LPID(IDX,1),LPID(IDX,2),LPID(IDX,3)
         WRITE(LUN,*) '   initial 4momentum:',(PP(IDX,JJ),JJ=1,5)
      ENDIF
      PP(IDX,1) = PX
      PP(IDX,2) = PY
      PP(IDX,3) = PZ
      PP(IDX,4) = EN
      PP(IDX,5) = XMS
c     return reference to other partons
      IREFout = LPID(IDX,3)
      IF(NDEBUG.gt.6)
     &     WRITE(LUN,*) '   final 4momentum:',(PP(IDX,JJ),JJ=1,5)
      END


      SUBROUTINE RD_PRTN(IDX,PX,PY,PZ,EN,XMS,IFL,IREFout)
C-------------------------------------------------------
C     routine to read a parton from the stack     \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_prtns_cmmn.inc'
      IF(NDEBUG.gt.6)THEN
         WRITE(LUN,*) ' RD_PRTN: (#,PID,LEVEL,REF)',
     &        IDX,LPID(IDX,1),LPID(IDX,2),LPID(IDX,3)
         WRITE(LUN,*) '   4momentum:',(PP(IDX,JJ),JJ=1,5)
      ENDIF
      PX = PP(IDX,1)
      PY = PP(IDX,2)
      PZ = PP(IDX,3)
      EN = PP(IDX,4)
      XMS = PP(IDX,5)
      IFL = LPID(IDX,1)
c     return reference to other partons
      IREFout = LPID(IDX,3)
      END


      SUBROUTINE RD_PRTN_4VEC(IDX,Pin,IFL,IREFout)
C-------------------------------------------------------
C     routine to read a parton from the stack     \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_prtns_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'
      DIMENSION Pin(5)

      IF(IDX.EQ.0) THEN
         WRITE(LUN,*) 'RD_PRTN_4VEC: invalid index!',IDX
         xa = -ONE
         xa = log(xa)
         RETURN
      ELSE         
         do ii = 1,5
            PIN(ii) = PP(IDX,ii)
         enddo
         IFL = LPID(IDX,1)
c     return reference to other partons
         IREFout = LPID(IDX,3)
         IF(NDEBUG.gt.6)THEN
            WRITE(LUN,*) ' RD_PRTN: (#,PID,LEVEL,REF)',
     &           IDX,IFL,LPID(IDX,2),IREFout
            WRITE(LUN,*) '   4momentum:',(PIN(JJ),JJ=1,5)
         ENDIF

      ENDIF
      END


      SUBROUTINE ITR_LVL0_PRTN(JJ,IDX,LID)
C-------------------------------------------------------
C     routine that serves as iterator over the level0
C     partons on the stack                        \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_prtns_cmmn.inc'      
      IDX = LVL0IDX(JJ)
      IF(ndebug.gt.6)
     & WRITE(LUN,*) 'ITR_LVL0_PRTN: JJ,IDX',JJ,IDX
      LID = LPID(IDX,1)
      IF(JJ+1.gt.NPP0) THEN
         JJ = -1
         RETURN
      ELSE
         JJ = JJ + 1
      ENDIF      
      END


      SUBROUTINE INI_PRTN_STCK(NOLD,N0OLD)
C-------------------------------------------------------
C     reset parton stack                          \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_prtns_cmmn.inc'
      IF(NDEBUG.gt.6) WRITE(LUN,*) 'PRTN_STCK: reset .. '
      IF(NDEBUG.gt.6) WRITE(LUN,*) ' OLD STATE: NPP,NPP0',NPP,NPP0
      
      NPP = NOLD
      NPP0 = N0OLD

      IF(NDEBUG.gt.6) WRITE(LUN,*) ' NEW STATE: NPP,NPP0',NPP,NPP0

      END


      SUBROUTINE GET_NPP(NPPLD,NPP0LD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_prtns_cmmn.inc'
      NPPLD = NPP
      NPP0LD = NPP0
      END


      SUBROUTINE GET_LVL0(NPP0LD,IDXLIST)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_prtns_cmmn.inc'
      DIMENSION IDXLIST(NPP0_max)
      NPP0LD = NPP0
      IDXLIST = LVL0IDX
      END


      SUBROUTINE PRNT_PRTN_STCK
C-------------------------------------------------------
C     as the name suggests, prints the current state
C     of the parton stack                         
C     print unit is defined in S_DEBUG:LUN        \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_prtns_cmmn.inc'
      CHARACTER*5 CDE
      CHARACTER*9 CODE
      INCLUDE 'sib_cnam_cmmn.inc'

      WRITE (LUN,50) 
 50   FORMAT(3X,75('-'),/,21X,'SIBYLL PARTON LEVEL EVENT SUMMARY',21X,
     &     /,3X,75('-'),13('-'))

c     beam particles
      WRITE(LUN,*) '  BEAM PARTICLES'
 52   FORMAT(4X,'#',3X,'PID',2x,'LVL',2x,'REF',20x,'PX',9x,'PY',7x,
     +     'PZ',9x,'E',11X,'Mass', /, 3X,75('-'),13('-'))
      WRITE (LUN,52)
      DO J=1,NPP
         IF(LPID(J,2).eq.-2)then
            WRITE (LUN,60) J, (LPID(J,KK),KK=1,3), (PP(J,K),K=1,5)
         ENDIF
      ENDDO
c     level -2 format
 60   FORMAT(4(1X,I4),16X,2(F9.3,2X),2(E9.3,2X),F9.3)
 61   FORMAT(3X,75('-'),13('-'))
      WRITE(LUN,61)

c     interactions
      WRITE(LUN,*) '  INTERACTIONS'
 62   FORMAT(4X,'#',3X,'PID',2x,'LVL',2x,'REF',20x,'NSOF',9x,'NJET',7x,
     +     'JDIF',9x,'E',11X,'Mass', /, 3X,75('-'),13('-'))
      WRITE (LUN,62)
      DO J=1,NPP
         IF(LPID(J,2).eq.-1)then
            WRITE (LUN,63) J, (LPID(J,KK),KK=1,3), (PP(J,K),K=1,5)
         ENDIF
      ENDDO
c     level -1 format
 63   FORMAT(4(1X,I4),16X,4(F8.0,4X),F9.3)
 64   FORMAT(3X,75('-'),13('-'))
      WRITE(LUN,64)
      
c     partons
      WRITE (LUN,100)
      DO J=1,NPP
         IF(LPID(J,2).eq.0)then
            WRITE (LUN,120) J, (LPID(J,KK),KK=1,3), (PP(J,K),K=1,5)
         elseif(LPID(J,2).eq.1)then
            call kcode(LPID(J,1),cde,nc)
            WRITE (LUN,121) J, CDE(1:nc),(LPID(J,KK),KK=2,3), 
     &           (PP(J,K),K=1,5)
         elseif(LPID(J,2).eq.2)then
            CODE = '        '
            L = LPID(J,1)
            CODE(1:6) = NAMP(IABS(L))
            IF (L .LT. 0) CODE(7:9) = 'bar'
            WRITE (LUN,122) J, CODE,(LPID(J,KK),KK=2,3), (PP(J,K),K=1,5)
         endif
      ENDDO
      CALL PPsum(1,NPP,Esum,PXsum,PYsum,PZsum,NF)
      WRITE(LUN,140) PXsum,PYsum,PZsum,Esum
 
 100  FORMAT(4X,'#',3X,'PID',2x,'LVL',2x,'REF',20x,'PX',9x,'PY',7x,
     +     'PZ',9x,'E',11X,'Mass', /, 3X,75('-'),13('-'))
c     level 0 format
 120  FORMAT(4(1X,I4),16X,2(F9.3,2X),2(E9.3,2X),F9.3)
c     level 1 format
 121  FORMAT(3X,I4,1X,A5,2(1X,I4),16X,2(F9.3,2X),2(E9.3,2X),F9.3)
c     level 2 format
 122  FORMAT(6X,I4,1X,A9,2(1X,I4),12X,2(F9.3,2X),2(E9.3,2X),F9.3)
 140  FORMAT(3X,75('-'),13('-'),/,
     &     1X,'Tot = ',29X,2(F9.3,2X),G9.3,2X,E9.3)

      END


      SUBROUTINE PPsum(N1,N2,ETOT,PXT,PYT,PZT,NF)
C-------------------------------------------------------
C     Return the energy,px,py,pz of level0 partons 
C     in the list between N1 and N2
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'sib_prtns_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'

      NF=0
      ETOT=ZERO
      PXT=ZERO
      PYT=ZERO
      PZT=ZERO
      DO J=N1,N2
         IF (LPID(J,2) .EQ. 0)  THEN
           NF = NF+1
           ETOT = ETOT + ABS( PP(J,4) )
           PXT = PXT + PP(J,1)
           PYT = PYT + PP(J,2)
           PZT = PZT + PP(J,3)
         ENDIF
      ENDDO
      RETURN
      END

      SUBROUTINE FOUR_LENGTH(XP,XM2)
C-------------------------------------------------------
C     Calculate the length of a 4vector (+---)    \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XP(5)
      XM2 = XP(4)**2 - XP(1)**2 - XP(2)**2 - XP(3)**2
      END

      DOUBLE PRECISION FUNCTION FOUR(XP)
C-------------------------------------------------------
C     Calculate the length of a 4vector (+---)    \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XP(5)
      FOUR = XP(4)**2 - XP(1)**2 - XP(2)**2 - XP(3)**2
      END

      DOUBLE PRECISION FUNCTION CALC_INVM(XP1,XP2)
C-------------------------------------------------------
C     Calculate the invariant mass of two 4vectors FR'15
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XP1(5),XP2(5)
      CALC_INVM = (XP1(4)+ XP2(4))**2
      DO I=1,3
         CALC_INVM = CALC_INVM-(XP1(I)+XP2(I))**2         
      ENDDO
      CALC_INVM = SQRT(CALC_INVM)
      END

      
      SUBROUTINE GET_XMT2(IDX,XM2)
C-------------------------------------------------------
C     Calculate the transverse mass of a parton 
C     on the stack                                \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'sib_prtns_cmmn.inc'
      XM2 = PP(IDX,1)**2 + PP(IDX,2)**2 + PP(IDX,5)**2
      END

      SUBROUTINE GET_IMASS2(IDX,XM2)
C-------------------------------------------------------
C     Calculate the invariant mass squared of a parton 
C     on the stack  (+---)                        \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'sib_prtns_cmmn.inc'
      XM2 = PP(IDX,1)**2 + PP(IDX,2)**2 + PP(IDX,3)**2
      XM2 = PP(IDX,4)**2 - XM2
      END


      SUBROUTINE GET_MASS(IDX,XM)
C-------------------------------------------------------
C     read mass of parton on stack                \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'sib_prtns_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'

      IF(IDX.EQ.0) THEN 
         XM2 = ZERO
      else
         XM = PP(IDX,5)
      ENDIF
      END

      SUBROUTINE GET_MASS2(IDX,XM2)
C-------------------------------------------------------
C     read mass of parton on stack                \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'sib_prtns_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'

      IF(IDX.EQ.0) THEN 
         XM2 = ZERO
      else
         XM2 = PP(IDX,5)**2
      ENDIF
      END


      SUBROUTINE GET_VRTLTY(IDX,XX)
C-------------------------------------------------------
C     calculate virtuality of parton on stack     \FR'14
C     = on-shell mass - inv. mass
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'sib_prtns_cmmn.inc'
      COMMON /SIB_NUM/ ZERO,HALF,ONE,ONEHALF,TWO,THREE,FOUR,EIGHT
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      IF(IDX.EQ.0) XM2 = ZERO
      call get_imass2(IDX,xm2)
      XX = PP(IDX,5)**2-xm2
      END


      SUBROUTINE ADD_4VECS(P1,P2,POUT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'sib_utl_cmmn.inc'
      INCLUDE 'sib_debug_cmmn.inc'
      DIMENSION P1(5),P2(5),POUT(5)
      DO II=1,4
         POUT(II) = P1(II) + P2(II)
      ENDDO
      CALL FOUR_LENGTH(POUT,XM2)
      IF(XM2.LT.0)THEN
c     virtual particle
         POUT(5) = -ONE
         IF(NDEBUG.gt.6)then
            WRITE(LUN,*)
     &           'ADD_4VECS: resulting particle virtual!! (m**2):',XM2
            WRITE(LUN,*) 'p**2' , POUT(1)**2+POUT(2)**2+POUT(3)**2
            WRITE(LUN,*) 'E**2: ', POUT(4)**2
         ENDIF
      ELSE
         POUT(5) = sqrt(xm2)
      ENDIF
      END
