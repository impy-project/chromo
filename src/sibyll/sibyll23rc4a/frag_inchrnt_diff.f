      SUBROUTINE FRAG_INCHRNT_DIFF(IDX,LBAD)
C-----------------------------------------------------------------------
C     routine that fragments a diffractive system               \FR'15
C
C     INPUT: IDX : parton stack index of 4momentum
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER IDX,LBAD
      
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_plist_cmmn.inc'
      INCLUDE 'sib_parto_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'

      DOUBLE PRECISION PST,PDIFF,GABE,P2,EE,P1TOT
      DIMENSION PST(5),PDIFF(5),GABE(4),P2(4)
      INTEGER IDIFF1,IDIFF,IPID,L0,JDIFF,NOLD,LXBAD,K,II
      LBAD = 2

c     references are diff --> diff.hadron --> bm-partons --> tg-partons
c     only diff and diff. hadron are read out
c     read diff 4momentum from stack
      call rd_prtn_4vec(IDX,PST,IPID,IDIFF1)
      call rd_prtn_4vec(IDIFF1,PDIFF,L0,IDIFF)
      
C     kinematic variables
      EE = PDIFF(5)             ! center of mass energy in diff. system
      
c     set diffraction code of system (1:beam,2:target,3:double)
      JDIFF = ABS(IPID)/10

      IF(NDEBUG.gt.1) WRITE(LUN,*)' FRAG_INCHRNT_DIFF: IDX,EE,L0',
     &     IDX,EE,L0

      IPFLAG = -1

      NOLD = NP

c     diffractive interaction in center-of-mass system of (sea,rmnt)-nuc
      CALL SIB_DIFF(L0,JDIFF,EE,0,LXBAD)
      IF(LXBAD.ne.0) THEN
         IF(NDEBUG.gt.1) 
     &        WRITE(LUN,*)' FRAG_INCHRNT_DIFF: fragmentation rejection'         
         RETURN
      ENDIF
      IF(NDEBUG.gt.1) 
     &     WRITE(LUN,*)' FRAG_INCHRNT_DIFF: particles before/after :',
     &     NOLD,NP

c     boost to hadron - hadron center-of-mass
      do ii=1,4
         gabe(ii) = PDIFF(ii)/PDIFF(5)
      enddo
      DO K=NOLD+1,NP
         call SIB_ALTRA(gabe(4),gabe(1),gabe(2),
     &        gabe(3),P(k,1),p(k,2),p(k,3),p(k,4),
     &        P1TOT,p2(1),p2(2),p2(3),p2(4))
         do ii=1,4
            P(K,ii)=P2(ii)
         enddo
      ENDDO

      LBAD = 0
      END
