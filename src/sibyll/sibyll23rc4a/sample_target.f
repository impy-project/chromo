      SUBROUTINE sample_target(NW,XCHG,KRMNT,XJET,Irec,IREJ)
C...Subroutine to sample valence and sea quark kinematic variables
C     on the target side
C.    fills IFLT,X2 and PXT,PYT
C.    1,2 are valence quarks, 3,4 are additional sea quarks
C.    transverse momentum is shared between the val. and sea pairs
C.    X and flv are exchanged occasionally, not pt so far
C-------------------------------------------------------------------      
      IMPLICIT NONE
      SAVE

      include 'sib_nw_prm.inc'
c     external types
      DOUBLE PRECISION XJET,XCHG
      DIMENSION XJET(NW_max)
      INTEGER KRMNT,NW,IREC,IREJ
      DIMENSION KRMNT(NW_max)

      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      INCLUDE 'sib_cutoff_cmmn.inc'
      INCLUDE 'sib_run_cmmn.inc'
      include 'sib_int_prm.inc'
      INCLUDE 'sib_indx_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'
      INCLUDE 'sib_rmnt_cmmn.inc'

c     internal types
      DOUBLE PRECISION XX,X2,PX,PXT,PY,PYT,PZ,PZ1,PZ2
      DIMENSION XX(2*NW_max+2),PX(2*NW_max+2),PY(2*NW_max+2)
      DIMENSION X2(4*NW_max),PXT(4*NW_max),PYT(4*NW_max)
      INTEGER IFL,IFLT,IREJ1,J,J1,J2,J3,J4,JJ,JJJ,JI,I,KID,Iref1,
     &     Iref,KID1
      DIMENSION IFL(2*NW_max+2),IFLT(4*NW_max)

      IREJ1 = 1

      IF(ndebug.gt.2) 
     +     WRITE(LUN,*)
     +     ' SAMPLE_target: NW,XCHG,LRMNT,XJET,IREC,IREJ',
     +     NW,XCHG,(KRMNT(j),j=1,NW),(XJET(j),j=1,NW),IREC,IREJ

      DO J=1,NW ! zero arrays
         j1 = 1+4*(j-1)
         j2 = j1 + 1
         j3 = j2 + 1
         j4 = j3 + 1
         X2(j1) = ZERO
         X2(j2) = ZERO
         X2(j3) = ZERO
         X2(j4) = ZERO
         PXT(j1) = ZERO
         PXT(j2) = ZERO
         PXT(j3) = ZERO
         PXT(j4) = ZERO
         PyT(j1) = ZERO
         PyT(j2) = ZERO
         PyT(j3) = ZERO
         PyT(j4) = ZERO
      ENDDO

      DO j=1,NW
c     read target id from event info 
         KID = KT(J)
c     reset rejection
         IREJ = IREJ1
c     always fills remnant partons into 1,2 and c.strings into 3,4
c     so far only one interaction possible (beam is always a single hadron!)        
         call sample_projectile
     +        (KID,1,KRMNT(j),XCHG,XJET(j),XX,PX,PY,IFL,KID1,IREJ)
         IF(IREJ.ne.0) RETURN

c     write to target variables
         do jj=3-2*KRMNT(j),4
            ji = jj+4*(j-1)
            IFLT(ji) = IFL(jj)
            X2(ji) = XX(jj)
            PXT(ji) = PX(jj)
            PYT(ji) = PY(jj)
         enddo

         IF(KRMNT(j).ne.0)THEN
c     by convention hadron is split such that diq is 2nd flv
c     for string frag routine argument flv1 is along +z, flv2 -z
c     by convention again flv2 in the remnant is passed to +z and flv1 to -z
c     therefor on the target side the flavors need to be switched such that
c     the diq is along -z
            j1 = 1+4*(j-1)
            j2 = j1 + 1
            call iswtch_lmnts(IFLT(j1),IFLT(j2))
         ENDIF

c     central strings
c     flavors need to be switched as well (strictly speaking color)
c     in dual-parton model: q : color , diq : anticolor
c     need to combine q with diq for color neutral system..
         j3 = 3+4*(j-1)
         j4 = j3 + 1
         call iswtch_lmnts(IFLT(j3),IFLT(j4))
         call swtch_lmnts(X2(j3),X2(j4))
         
c     reset remnant id 
c     might have changed in flavor exchange (actually color)...
         KRT(J) = KID1
      ENDDO

C..   write target partons to stack
      DO I=1,NW
         IF(KRMNT(I).ne.0)THEN
c     add proto-remnant
            j1 = 1+4*(i-1)
            j2 = j1 + 1
            call add_prtn(PXT(J1)+PXT(J2),PYT(J1)+PYT(J2),
     &           -HALF*SQS*(X2(J1)+X2(j2)),HALF*SQS*(X2(J1)+X2(j2)),
     &           ZERO,-2,0,0,Iref1)
            ITGRDX(I,1) = Iref1
            call add_int_ref(Iref1,IINTDX(I))
c     add quarks to stack
            do j = 1,2
               jj = 4*(i-1)+j
               jjj = 4*(i-1)+j + 2
               pz1 = (HALF*SQS*X2(JJ))**2
c               PZ1 = (HALF*SQS*X2(JJ))**2-PXT(JJ)**2-PYT(JJ)**2
               call add_prtn(PXT(JJ),PYT(JJ),-sqrt(pz1),
     &              HALF*SQS*X2(JJ),ZERO,IFLT(JJ),1,Iref1,Iref)
               ITGRDX(I,j+1) = Iref
               pz2 = (HALF*SQS*X2(JJj))**2
c               pz2 = (HALF*SQS*X2(JJj))**2-PXT(JJj)**2-PYT(JJj)**2
               call add_prtn(PXT(JJj),PYT(JJj),-sqrt(pz2),
     &              HALF*SQS*X2(JJj),ZERO,IFLT(JJj),1,0,Iref)
               ICSTDX(2*(I-1)+j,3) = Iref
            enddo
         else
            do j = 3,4
               jj = 4*(i-1)+j
               pz = (HALF*SQS*X2(JJ))**2
c               pz = (HALF*SQS*X2(JJ))**2-PXT(JJ)**2-PYT(JJ)**2
               call add_prtn(PXT(JJ),PYT(JJ),-sqrt(pz),
     &              HALF*SQS*X2(JJ),ZERO,IFLT(JJ),1,0,Iref)
               ICSTDX(2*(I-1)+(J-2),3) = Iref
            enddo
         ENDIF
      ENDDO
      IF(NDEBUG.GT.3) call prnt_prtn_stck

      IREJ = 0
      END
