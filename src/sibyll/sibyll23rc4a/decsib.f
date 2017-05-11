      SUBROUTINE DECSIB 
C-----------------------------------------------------------------------
C...Decay all unstable particle in Sibyll
C.  decayed particle have the code increased by 10000
C
C   changed to allow for multiple calls to DECSIB in one event
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      INCLUDE 'sib_csydec_cmmn.inc'
      INCLUDE 'sib_plist_cmmn.inc'
      INCLUDE 'sib_plist1_cmmn.inc'
      INCLUDE 'sib_parto_cmmn.inc'
      INCLUDE 'sib_mass_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      INCLUDE 'sib_rnk_cmmn.inc'
      DIMENSION P0(5), LL(10), PD(10,5)

c     call pythia decay routine      
      IF(IPAR(44).eq.1) CALL PYDEC

c     decay with sibyll
      NN = 1
      IF(IPAR(44).ne.1)THEN
         DO J=1,NP
            LLIST1(J) = 0
         ENDDO
      ENDIF
      DO WHILE (NN .LE. NP)
         L= LLIST(NN)
         LA = IABS(L)
         if(LA.lt.100) then
           IF (IDB(LA) .GT. 0)  THEN
              DO K=1,5
                P0(K) = P(NN,K)
              ENDDO
              CALL DECPAR (L,P0,ND,LL,PD)
              LLIST(NN) = LLIST(NN)+ISIGN(10000,LLIST(NN))
              DO J=1,ND
                NP = NP+1
                if(NP.gt.8000) then
                  write(6,'(1x,a,2i8)') 
     &              'DECSIB: no space left in S_PLIST (NP,ND):',NP,ND
                  NP = NP-1
                  return
                endif
                DO K=1,5
                  P(NP,K) = PD(J,K)
                ENDDO
                LLIST(NP)=LL(J)
                LLIST1(NP)=NN
                LRNK(NP)=LRNK(NN)
                NPORIG(NP)= NPORIG(NN)
                niorig(NP)= NIORIG(NN)
                NFORIG(NP) = L
              ENDDO
           ENDIF
         endif
         NN = NN+1
      ENDDO

c	call sib_list(20)

      END
