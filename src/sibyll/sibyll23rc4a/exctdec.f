
      SUBROUTINE EXCTDEC( IDX, LBAD)
C-----------------------------------------------------------------------
C     routine to fragment an excited system with known flavor via
C     resonance decay
C
C     
C-----------------------------------------------------------------------
      IMPLICIT NONE
c     external variables
      INTEGER IDX,LBAD
      
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_plist_cmmn.inc'
      INCLUDE 'sib_parto_cmmn.inc'
      INCLUDE 'sib_rnk_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      INCLUDE 'sib_mass_cmmn.inc'
      INCLUDE 'sib_chp_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'

c     local variables
      DOUBLE PRECISION P0,BE,PR1,PR2,PRH,GABE,P2,
     &     PAR2_def,PAR8_def,PAR24_def,DELTAE,PCXG,
     &     EMIN1,EMIN2,EMIN3,EMIN4,S_RNDM,GA,PTR,PTOT,P1TOT,PX,PY,
     &     COD,SID,COF,SIF,ANORF,BEP
      DIMENSION P0(5),BE(3),PR1(5),PR2(5),PRH(5),GABE(4),
     &     P2(5)
      INTEGER IPID,IR1DX,IFLR1,IR2DX,IFLR2,IRH,IRHPID,IR,
     &     KK,KD,IFAIL,N1,IFBAD,J,K,I
      
c      LBAD = 1

c     initial parameters
      PAR2_def = PAR(2)         ! ud/s rate
      PAR8_def = PAR(8)         ! popcorn rate
      PAR24_def = PAR(24)       ! c/s rate
      if(ndebug.gt.1)
     &     WRITE(LUN,*) ' EXCTDEC: IDX,IREJ',IDX,LBAD
      
c     read remnant 4momentum from stack
      call rd_prtn_4vec(IDX,P0,IPID,IR1DX)
      call rd_prtn_4vec(IR1DX,PR1,IFLR1,IR2DX)
      call rd_prtn_4vec(IR2DX,PR2,IFLR2,IRH)
      call rd_prtn_4vec(IRH,PRH,IRHPID,IR)
      IPFLAG = IPID
      IF(IDX.ne.IR)then
         write(lun,*) ' EXCTDEC: reference loop broken!' , IDX
         call sib_reject
      endif
      IF(NDEBUG.GT.2)THEN
         WRITE(LUN,*) ' EXCTDEC: P0:' , (P0(kk),kk=1,5)
         WRITE(LUN,*) ' EXCTDEC: PR1:' , (PR1(kk),kk=1,5)
         WRITE(LUN,*) ' EXCTDEC: PR2:' , (PR2(kk),kk=1,5)
         WRITE(LUN,*) ' EXCTDEC: PH:' , (PRH(kk),kk=1,5)
      ENDIF
      
C     identity of remnant
c     form hadron from flavors in remnant
c     (not preserving spin or isospin!)
c      call sib_i4flav(iflr1,iflr2,Idm, KD )
      KD = IRHPID

c     available kinetic energy
      DELTAE = P0(5)-AM(ABS(KD))
c     fallback region: 0 < DELTAE < EMIN1
      EMIN1 = PAR(76)
c     resonance region: EMIN1 < DELTAE < EMIN2
      EMIN2 = PAR(77)
c     phasespace decay region: EMIN2 < DELTAE < EMIN3
      EMIN3 = PAR(78)
c     string decay region: EMIN3 < DELTAE < EMIN4
      EMIN4 = PAR(79)

      IF(NDEBUG.gt.2)THEN
         WRITE(LUN,*) 
     &        ' EXCTDEC: MASS,IFL1,IFL2,PID',P0(5),IFLR1,IFLR2,KD
         WRITE(LUN,*) ' EXCTDEC: DELTAE,EMIN1,EMIN2,EMIN3',
     &        DELTAE,EMIN1,EMIN2,EMIN3
      ENDIF
      
c     strange quark rate
      IF(IPAR(48).eq.1)THEN
         PAR(2) = PAR(89)
      ENDIF

c     charm quark rate
      IF(IPAR(62).eq.1)THEN
         PAR(24) = PAR(107)
      ENDIF
     
c     popcorn rate in remnant
      IF(IPAR(56).eq.1)THEN
         PAR(8) = PAR(102)
      ENDIF     

      IF(DELTAE.lt.EMIN2)THEN
c     beam or resonance region
         IF(NDEBUG.gt.1) then 
            if(DELTAE.lt.EMIN1)then
               WRITE(LUN,*)' EXCTDEC: fallback to beam..'
            else
               WRITE(LUN,*)' EXCTDEC: forming resonance..'
            endif
         endif
         NP = NP + 1
         LLIST(NP) = KD
         NPORIG(NP) = IPFLAG
         LRNK(NP) = 0
         niorig(NP) = iiflag
         DO kk=1,5
            P(NP,KK) = P0(KK)
         ENDDO
         LBAD = 0
         PAR(2) = PAR2_def
         PAR(8) = PAR8_def
         PAR(24) = PAR24_def
         RETURN         

      ELSEIF(DELTAE.lt.EMIN3)THEN
c     phasespace decay region
         IF(NDEBUG.gt.1) WRITE(LUN,*)' EXCTDEC: phasespace decay ..'
         IPFLAG = IPID/iabs(IPID) + ISIGN(1000,IPID)
c     set charge exchange probability, 
c     i.e. prob for p* -> n + pip
         PCXG = PAR(99)
         CALL FIREBALL_4FLV(KD,P0,PCXG,IFAIL)
         PAR(2) = PAR2_def
         PAR(8) = PAR8_def
         PAR(24) = PAR24_def
         IF(IFAIL.eq.1) THEN
            IF(ndebug.gt.0)
     &           WRITE(6,*) ' RMNT_FRAG: remnant frag. rejection!'
            LBAD = 1
            RETURN
         ENDIF
         LBAD = 0
         RETURN

c      ELSEIF(DELTAE.lt.EMIN4)THEN
      ELSE
C     string fragmentation region
         IF(NDEBUG.gt.1) WRITE(LUN,*)' EXCTDEC: string decay ..'
         N1 = NP+1
         IPFLAG = IPFLAG + ISIGN(3000,IPID)
c     for meson remnant quark and anti-quark should be treated equally
c     therefor switch randomly
         IF(IBAR(ABS(KD)).eq.0.and.S_RNDM(KD).lt.HALF)
     &        call iswtch_lmnts(IFLR1,IFLR2)

c     turn remnant string around
         IF(IPAR(23).eq.1)THEN
            IF(S_RNDM(KD).gt.PAR(39))
     &           call iswtch_lmnts(IFLR1,IFLR2)
         ENDIF

         CALL STRING_FRAG_4FLV 
     +        (P0(5), IFLR2, IFLR1, ZERO,ZERO,ZERO,ZERO,IFBAD,1)
         IF (IFBAD .EQ. 1)THEN
            IF(ndebug.gt.0)
     &           WRITE(6,*) ' EXCTDEC: remnant frag. rejection!'
            LBAD = 1
            PAR(2) = PAR2_def
            PAR(8) = PAR8_def
            PAR(24) = PAR24_def
            RETURN
         ENDIF
         DO J=1,3
            BE(J)=P0(J)/P0(4)
            GABE(J)=P0(J)/P0(5)
         ENDDO
         GA=P0(4)/P0(5)
         GABE(4)=P0(4)/P0(5)
C...  rotate and boost string
         SELECT CASE(IPAR(38))
         CASE(1,3)
c     sample additional soft pt for remnant partons
            call ptdis_4flv(0,PX,PY)
            PTR = SQRT(PX**2+PY**2)
            PTOT = SQRT(FOUR*PTR**2+P0(5)**2)*HALF
c     rotation factors
            COD = HALF*P0(5)/PTOT
            SID = PTR/PTOT
c            COD= ONE/SQRT(ONE+FOUR*PTR**2/P0(5))
c            SID= TWO*PTR/P0(5)*COD
            COF=ONE
            SIF=ZERO
            IF(PTOT*SID.GT.EPS5) THEN
               COF=PX/(SID*PTOT)
               SIF=PY/(SID*PTOT)
               ANORF=DSQRT(COF*COF+SIF*SIF)
               COF=COF/ANORF
               SIF=SIF/ANORF
            ENDIF
            IF(ndebug.gt.3)THEN
            write(lun,*)' EXCTDEC: rotation factors (cod,sid,cof,sif):',
     &           cod,sid,cof,sif
            write(lun,*)' EXCTDEC: rotation angles (theta,phi):',
     &           ACOS(cod),ACOS(cof)
            ENDIF
c     rotate string final state
            DO K=N1,NP
               call SIB_TRANI(P(K,1),P(k,2),P(k,3),cod,sid,cof,sif
     &              ,P2(1),P2(2),P2(3))
               do j=1,3
                  P(K,j)=P2(j)
               enddo
            ENDDO
c     boost to hadron-hadron center-of-mass
            IF(ndebug.gt.3)
     &      write(lun,*) ' EXCTDEC: boost to had-had (gabe,gam):',
     &           (gabe(j),j=1,4)
            DO K=N1,NP
               NPORIG(K) = IPFLAG
               niorig(K) = iiflag
               call SIB_ALTRA(gabe(4),gabe(1),gabe(2),
     &              gabe(3),P(k,1),p(k,2),p(k,3),p(k,4),
     &              P1TOT,p2(1),p2(2),p2(3),p2(4))
               do j=1,4
                  P(K,j)=P2(j)
               enddo
            ENDDO
         CASE(0,2)
C...  boost string
            DO I=N1,NP
               NPORIG(I) = IPFLAG
               niorig(I) = iiflag
               BEP=BE(1)*P(I,1)+BE(2)*P(I,2)+BE(3)*P(I,3)
               DO J=1,3
                  P(I,J)=P(I,J)+GA*(GA*BEP/(ONE+GA)+P(I,4))*BE(J)
               ENDDO
               P(I,4)=GA*(P(I,4)+BEP)
            ENDDO
         END SELECT
      ENDIF
      LBAD = 0
      PAR(2) = PAR2_def
      PAR(8) = PAR8_def
      PAR(24) = PAR24_def
      END

