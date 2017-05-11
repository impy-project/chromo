
      subroutine sample_minijet
     &     (L,NW,NNJET,NNSOF,NJET,NSOF,X1JET,X2JET,LBAD)
C-----------------------------------------------------------------------
C     routine to sample minijets
C     INPUT: L - hadron type (1:nucleon,2:pion or 3:kaon)
C            NW - number of hadron-nucleon interactions
C            NNJET(1:NW) - number of hard interactions per nucleon
C            NNSOF(1:NW) - number of soft interactions per nucleon
C     OUTPUT: X1JET - momentum fraction of beam in minijets
C             X2JET(1:NW) - momentum fraction of target in minijets
C     
C     in addition minijets are added to parton stack
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      include 'sib_int_prm.inc'
      include 'sib_nw_prm.inc'

      INCLUDE 'sib_run_cmmn.inc'
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      INCLUDE 'sib_cutoff_cmmn.inc'
      INCLUDE 'sib_mass_cmmn.inc'
c$$$      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
c$$$      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, kb ,kt
c$$$      COMMON /S_MASS1/ AM(99), AM2(99)
c$$$      COMMON /S_CUTOFF/ STR_mass_val, STR_mass_val_hyp, STR_mass_sea

      include 'sib_indx_cmmn.inc'
c      PARAMETER (NS_max = 20, NH_max = 100)
c      COMMON /S_INDX/ IBMRDX(3),ITGRDX(NW_max,3),IHMJDX(NH_max),
c     &     ISMJDX(NS_max),ICSTDX(2*NW_max,3)
      DIMENSION NNSOF(NW_max),NNJET(NW_max),X2JET(NW_max)
      INCLUDE 'sib_utl_cmmn.inc'

      if(Ndebug.gt.1) WRITE(LUN,*)
     &     ' SAMPLE_minijets: (L,NW,NNJET,NNSOF):',
     &     L,NW,(NNJET(ii),ii=1,nw),(NNSOF(ii),ii=1,nw)

      NJET = 0
      NSOF = 0
      Nall = 0
      X1JET = ZERO
      DO JW=1,NW
C...hard sea-sea interactions
         X2JET(JW) = ZERO
         DO JJ=1,NNJET(JW)
           Nall = Nall+1
           NJET = NJET+1
           CALL SAMPLE_hard (L,X1Jj,X2Jj,PTJET)
           X1JET = X1JET + X1Jj           
           X2JET(JW) = X2JET(JW)+X2Jj
           if(Ndebug.gt.2) THEN
              WRITE(LUN,*)
     &             ' SAMPLE_minijets: hard JJ,JW,X1JET,X2JET(JW):',
     &             JJ,JW,X1JET,X2JET(JW)
              WRITE(LUN,*)
     &             '  X1,X2,PT:',X1JJ,X2JJ,PTJET
           ENDIF
           IF ((X2JET(JW).GT.0.9D0).OR.(X1JET.GT.0.9D0)) then
              if(Ndebug.gt.2) WRITE(LUN,*)
     &        ' SAMPLE_minijets: not enough phase space',
     &             ' (Ncall,Njet,lbad):',Ncall,Njet,lBAD
              return
           ENDIF
           FI = TWOPI*S_RNDM(L)
           XM = SQS*sqrt(X1jj*X2jj)
           SQSHALF = HALF*SQS
c           TH = ASIN(MIN((ONE-EPS8),TWO*PTJET/XM))
c     add gluon-gluon string to stack
           call add_prtn
     &     (ZERO,ZERO,SQSHALF*(X1jj-X2jj),SQSHALF*(X1jj+X2jj),
     &          XM,100,0,0,Iref)
           call add_int_ref(Iref,IINTDX(JW))
c     add gluons to stack
           call add_prtn(PTJET*COS(FI),PTJET*SIN(FI),
     &          SQSHALF*X1jj,SQSHALF*X1jj,ZERO,0,1,0,Irefg1)
           call add_prtn(-PTJET*COS(FI),-PTJET*SIN(FI),
     &          -SQSHALF*X2jj,SQSHALF*X2jj,ZERO,0,1,Iref,Irefg2)
c     set up references
c     minijet --> gluon1 --> gluon2 --> minijet
           call add_ref(Irefg1,Irefg2)
           call add_ref(Iref,Irefg1)

         ENDDO

C...soft sea-sea interactions 
         NSOF_JW = 0
         DO JJ=1,NNSOF(JW)-1
c     different soft distributions
            SELECT CASE(IPAR(28))
            CASE(1)
               CALL SAMPLE_soft2 (STR_mass_sea,X1S,X2S,PTSOF)
            CASE(2)
               CALL SAMPLE_soft3 (STR_mass_sea,X1S,X2S,PTSOF)
            CASE(3)
               CALL SAMPLE_soft5 (STR_mass_sea,X1S,X2S,PTSOF)
            CASE(4)
               CALL SAMPLE_soft6 (STR_mass_sea,X1S,X2S,PTSOF)
            CASE default             
               CALL SAMPLE_soft (STR_mass_sea,X1S,X2S,PTSOF)
            END SELECT
            IF ((X2JET(JW)+X2S.LT.0.9D0).AND.(X1JET+X1S.LT.0.9D0)) THEN
               NSOF = NSOF+1
               Nall = Nall+1
               NSOF_JW = NSOF_JW+1
               X1JET = X1JET + X1S
               X2JET(JW) = X2JET(JW)+X2S
c     add to stack
c     add gluon-gluon string to stack
               XM = SQS*SQRT(X1S*X2S)
               SQSHALF = HALF*SQS
               PZ = SQSHALF*(X1S-X2S)
               EN = SQSHALF*(X1S+X2S)
               FI = TWOPI*S_RNDM(L)
               call add_prtn(ZERO,ZERO,PZ,EN,XM,10,0,0,Iref)
               call add_int_ref(Iref,IINTDX(JW))
c     add gluons to stack
               call add_prtn(PTSOF*COS(FI),PTSOF*SIN(FI),
     &              SQSHALF*X1S,SQSHALF*X1S,ZERO,0,1,0,Irefg1)
               call add_prtn(-PTSOF*COS(FI),-PTSOF*SIN(FI),
     &              -SQSHALF*X2S,SQSHALF*X2S,ZERO,0,1,Iref,Irefg2)
c     set up references
c     minijet --> gluon1 --> gluon2 --> minijet
               call add_ref(Irefg1,Irefg2)
               call add_ref(Iref,Irefg1)
               IF(Ndebug.gt.2)THEN
                  WRITE(LUN,*)
     &                 ' SAMPLE_minijets: soft JJ,JW,X1JET,X2JET(JW):',
     &                 JJ,JW,X1JET,X2JET(JW)
                  WRITE(LUN,*)
     &                 '  X1,X2,PT:',X1s,X2s,PTSOF
               ENDIF
            ELSE
               IF(Ndebug.gt.1) WRITE(LUN,*)
     &        ' SAMPLE_minijets: not enough phase space',
     &             ' (Ncall,Nsof,lbad):',Ncall,Njet,lBAD
               RETURN
            ENDIF
         ENDDO
         NNSOF(JW) = NSOF_JW+1
      ENDDO
      lbad = 0

      END
