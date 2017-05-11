      SUBROUTINE FRAG_VLNCE(IDX,LBAD)
C-----------------------------------------------------------------------
C     routine that fragments a quark - quark system               \FR'14
C
C     INPUT: IDX : parton stack index of central string
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER IDX,LBAD
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_plist_cmmn.inc'
      include 'sib_nw_prm.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      INCLUDE 'sib_chist_cmmn.inc'
      INCLUDE 'sib_parto_cmmn.inc'
      INCLUDE 'sib_mass_cmmn.inc'    
      INCLUDE 'sib_utl_cmmn.inc'

      DOUBLE PRECISION PST,PBM,PTG,PSTH,P1,P2,GABE,EE,
     &     PAR1_def,PAR24_def,PX1,PY1,PX2,PY2,GAM,BET,P1TOT,P2TOT,
     &     SIF,COF,COD,SID,ANORF,PZ
      DIMENSION PST(5),PBM(5),PTG(5),PSTH(5),P1(4),P2(4),GABE(4)
      INTEGER LSTH,IPID,IBMST,ITGST,ISTH,IFLB,IFLT,IST,I,IFBAD,JJ,
     &     NOLD,II,K,J
      
      LBAD = 2
      LSTH = 0
      
c     references are:
c     string --> bm-parton --> tg-parton (--> merged string/hadron)
c     read string 4momentum from stack
      call rd_prtn_4vec(IDX,PST,IPID,IBMST)
      call rd_prtn_4vec(IBMST,PBM,IFLB,ITGST)
      call rd_prtn_4vec(ITGST,PTG,IFLT,ISTH)

C     kinematic variables
      EE = PST(5)            ! string mass
      
      IF(NDEBUG.gt.1) WRITE(LUN,*)' FRAG_VLNCE: IDX,EE,IFLB,IFLT',
     &     IDX,EE,IFLB,IFLT

      IF(IDX.ne.ISTH) then
c     read merged string and add hadron to final particle stack..
         call rd_prtn_4vec(ISTH,PstH,LSTH,IST)
         IF(NDEBUG.gt.1) WRITE(LUN,*)' FRAG_VLNCE: found merged string',
     &        LSTH,(PSTH(I),I=1,5)         
         IF(IDX.ne.IST) then
            write(lun,*) 'FRAG_VLNCE: reference loop broken!' , IDX
            call sib_reject
         endif
         NP = NP + 1
         DO I=1,4
            P(NP,I) = PST(I)
         ENDDO
         P(NP,5) = AM(IABS(LSTH))
         LLIST(NP) = LSTH
         NPORIG(NP) = IPFLAG*2+NINT
         niorig(NP) = iiflag
         LBAD = 0
         RETURN
      ENDIF

c     baryon production setup
      PAR1_def = PAR(1)
      if( NSOF+NJET.gt.0) then
         PAR(1)= PAR(15)
      else
         PAR(1)= PAR(14)
      endif

c     charm fractions in different parameterizations
      PAR24_def = PAR(24)
      SELECT CASE(IPAR(15))
      CASE(2,3,4,5,6,8,9,10,11)
         PAR(24) = PAR(25)*EXP(-PAR(26)/EE)
      END SELECT

      IF(NDEBUG.gt.1)
     &     WRITE(LUN,*)' FRAG_VLNCE: parameters (CHM,DIQ,STR,VEC,POP)',
     &     PAR(24),PAR(1),PAR(2),PAR(5),PAR(8)

      NOLD=NP
      SELECT CASE(IPAR(38))
      CASE(1,2)
C...  rotate strings instead of attaching all pt to string end hadrons
         PX1 = ZERO
         PY1 = ZERO
         PX2 = ZERO
         PY2 = ZERO
      CASE(0,3)
c     assign pt to hadrons at string end (old model)
         PX1 = PBM(1)
         PY1 = PBM(2)
         PX2 = PTG(1)
         PY2 = PTG(2)
         GAM = PST(4)/EE
         BET = PST(3)/PST(4)
      END SELECT

C...  fragment strings in string restframe
      CALL STRING_FRAG_4FLV
     &     (EE,IFLB,IFLT,PX1,PY1,PX2,PY2,IFBAD,1)

      PAR(24) = PAR24_def
      PAR(1) = PAR1_def
      Nint= 0
      IF (IFBAD .EQ. 1) then
         if(Ndebug.gt.1) 
     &        WRITE(LUN,*)' STRING_FRAG: rejection (Ncall):',Ncall
         RETURN
      ENDIF

C...  rotate and boost string
      SELECT CASE(IPAR(38))
      CASE(1,2)
C     boost quark momentum to string center-of-mass 
c     to calculate rotation angles in string center-of-mass
         do jj=1,3
            gabe(jj) = PST(jj)/PST(5)
         enddo
         GABE(4) = PST(4)/PST(5)
         call SIB_ALTRA(gabe(4),-gabe(1),-gabe(2),-gabe(3),
     &        PBM(1),pbm(2),pbm(3),pbm(4),
     &        P1TOT,p1(1),p1(2),p1(3),p1(4))
         call SIB_ALTRA(gabe(4),-gabe(1),-gabe(2),-gabe(3),
     &        PTG(1),pTG(2),ptg(3),ptg(4),
     &        P2TOT,p2(1),p2(2),p2(3),p2(4))

c     should be back-to-back...
         IF(ndebug.gt.1)THEN
            write(lun,*)
     &      ' FRAG_VLNCE: string c.m. momentum, parton 1 (Pabs,P(i)):' ,
     &           P1TOT, (P1(j),j=1,4)
            write(lun,*)
     &      ' FRAG_VLNCE: string c.m. momentum, parton 2 (Pabs,P(i)):' ,
     &           P2TOT, (P2(j),j=1,4)
            write(lun,*) 'partons should be back to back...'
         ENDIF
c     rotation factors
         COD= P1(3)/P1TOT
         SID= DSQRT(P1(1)**2+P1(2)**2)/P1TOT
         COF=ONE
         SIF=ZERO
         IF(P1TOT*SID.GT.EPS5) THEN
            COF=P1(1)/(SID*P1TOT)
            SIF=P1(2)/(SID*P1TOT)
            ANORF=DSQRT(COF*COF+SIF*SIF)
            COF=COF/ANORF
            SIF=SIF/ANORF
         ENDIF
c     rotate string final state
         DO K=NOLD+1,NP
            call SIB_TRANI(P(K,1),P(k,2),P(k,3),cod,sid,cof,sif
     &           ,P2(1),P2(2),P2(3))
            do ii=1,3
               P(K,ii)=P2(ii)
            enddo
         ENDDO
c     boost to hadron - hadron center-of-mass
         DO K=NOLD+1,NP
            call SIB_ALTRA(gabe(4),gabe(1),gabe(2),
     &           gabe(3),P(k,1),p(k,2),p(k,3),p(k,4),
     &           P1TOT,p2(1),p2(2),p2(3),p2(4))
            do ii=1,4
               P(K,ii)=P2(ii)
            enddo
         ENDDO
      CASE(0,3)
C...  boost string
         DO K=NOLD+1,NP
            PZ = P(K,3)
            P(K,3) = GAM*(PZ+BET*P(K,4))
            P(K,4) = GAM*(P(K,4)+BET*PZ)
         ENDDO
      END SELECT
      LBAD = 0
      END
