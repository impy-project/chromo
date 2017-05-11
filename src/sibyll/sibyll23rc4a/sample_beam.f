      SUBROUTINE sample_beam(KID,NW,XCHG,KRMNT,XJET,IREJ)
C...Subroutine to sample valence and sea quark kinematics
C.    fills IFL?,X? and PX?,PY?
C.    1,2 are valence quarks, 3,4 are additional sea quarks
C.    transverse momentum is shared between the val. and sea pairs
C.    X and flv are exchanged occasionally
C-------------------------------------------------------------------      
      IMPLICIT NONE
      SAVE

      DOUBLE PRECISION XCHG,XJET
      INTEGER KID,NW,KRMNT,IREJ

      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      include 'sib_nw_prm.inc'
      include 'sib_int_prm.inc'
      INCLUDE 'sib_run_cmmn.inc'
      INCLUDE 'sib_indx_cmmn.inc'
      INCLUDE 'sib_rmnt_cmmn.inc'

      DOUBLE PRECISION X1,PXB,PYB
      DIMENSION X1(2*NW_max+2),PXB(2*NW_max+2),PYB(2*NW_max+2)
      INTEGER IFLB,KID1,J,J1,J2,J3,J4,Iref1,Iref,Idm
      DIMENSION IFLB(2*NW_max+2)

      INCLUDE 'sib_utl_cmmn.inc'      

c     default rejection
c     options are: 1: resample minijets (Xjet)..
c                  2: resample non-diff event (Ns,Nh)..
c                  3: resample event (Nw,diff,ndiff)..
      IREJ = 1

      IF(ndebug.gt.2) 
     +     WRITE(LUN,*)
     +     ' SAMPLE_beam: KID,NW,XCHG,KRMNT,XJET,IREJ',
     +     KID,NW,XCHG,KRMNT,XJET,IREJ

      call sample_projectile
     +     (KID,NW,KRMNT,XCHG,XJET,X1,PXB,PYB,IFLB,KID1,IREJ)    
      IF(IREJ.ne.0) RETURN

c     set remnant id to beam
c     will be changed if flavor is exchanged between central strings and remnant
      KRB = KID1        

C..   write beam partons to stack
c     order is: val1, val2, q, qbar etc
      IF(KRMNT.ne.0)THEN
         j1 = 1
         j2 = 2
c     add proto-remnant (still massless)
         call add_prtn(PXB(J1)+PXB(J2),PYB(J1)+PYB(J2),
     &        HALF*SQS*(X1(J1)+X1(J2)),
     &        HALF*SQS*(X1(J1)+X1(J2)),ZERO,2,0,0,Iref1)
         IBMRDX(1) = Iref1
c     beam remnant always associated with first interaction
         call add_int_ref(Iref1,IINTDX(1))
c     add quarks designated for remnant
         IF(KID.lt.0)THEN
c     if beam is antibaryon then hspli puts diq into 1st flv
c     need to switch to fit call to string frag routine 
c     such that diq is along +z
            call iswtch_lmnts(IFLB(j1),IFLB(j2))
         ENDIF
         call add_prtn(PXB(J1),PYB(J1),HALF*SQS*X1(J1),
     &        HALF*SQS*X1(J1),ZERO,IFLB(J1),1,Iref1,Iref)
         IBMRDX(2) = Iref
         call add_prtn(PXB(J2),PYB(J2),HALF*SQS*X1(J2),
     &        HALF*SQS*X1(J2),ZERO,IFLB(J2),1,Idm,Iref)
         IBMRDX(3) = Iref
      ENDIF
      DO j=1,NW
         j3 = 3+(j-1)*2
         j4 = j3+1
c     add sea quarks
         call add_prtn(PXB(J3),PYB(J3),HALF*SQS*X1(J3),
     &        HALF*SQS*X1(J3),ZERO,IFLB(J3),1,0,Iref)
         ICSTDX(2*(J-1)+1,2) = Iref
         call add_prtn(PXB(J4),PYB(J4),HALF*SQS*X1(J4),
     &        HALF*SQS*X1(J4),ZERO,IFLB(J4),1,0,Iref)
         ICSTDX(2*(J-1)+2,2) = Iref
c     add parton index to cache
      ENDDO
      IF(NDEBUG.GT.3) call prnt_prtn_stck

      IREJ = 0

      END
