      SUBROUTINE SAMPLE_RMNT(Kbeam,NW,X1JET,X2JET,Irec,LBAD)
C-----------------------------------------------------------------------
C     routine to sample remnants
C-----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
      include 'sib_nw_prm.inc'
      INCLUDE 'sib_cnt_cmmn.inc'
      
c     external type declarations
      DOUBLE PRECISION X1JET,X2JET
      DIMENSION X2JET(NW_max)
      INTEGER KBEAM,NW,IREC,LBAD

C     COMMONs
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      INCLUDE 'sib_chist_cmmn.inc'
      INCLUDE 'sib_rmnt_cmmn.inc'
      INCLUDE 'sib_chp_cmmn.inc'     
      INCLUDE 'sib_utl_cmmn.inc'     

c     internals
      DOUBLE PRECISION PREM,PREM_NUC,R,R2,S_RNDM,FLVXCHG,ALPH
      INTEGER ITGRMNT,IBMRMNT,I,j,jj,K,NPPLD,NPP0LD,IBMRMNT_OLD,
     &     IBAD,IKBAD,KBM
      DATA PREM /0.D0/ , PREM_NUC /0.D0/
      DIMENSION ITGRMNT(NW_max)

      IF(Ndebug.gt.1) 
     &  WRITE(LUN,*)' sample_RMNT: called with (Kbeam,NW,X1JET,',
     &     'X2JET,JDIF,Irec):',Kbeam,NW,X1JET,(X2JET(JJ),JJ=1,NW),
     &     (JDIF(JJ),JJ=1,NW),Irec

      IF(Irec.eq.0.and.NW.ne.1)then
         WRITE(LUN,*)' sample_RMNT: recursive call inconsistent!'
         call sib_reject
      endif

c     default return point for remnant excitation routine:
c     beam and target sampling
      IBAD = 1

c     set trial counter
      ITRY(2) = 0
c     remember position on parton stack
      call GET_NPP(NPPLD,NPP0LD)

C...  sample no. of remnants
c     ibmrmnt: 0,1..NW : number of excitations on beamside
C     itgrmnt: 0,1 : target side excitation

c     prob. of remnant excitation
      IF(IPAR(78).ne.0)THEN
         PREM = PAR(23)
         PREM_NUC = PAR(23)
         IF(IPAR(84).eq.2.and.IBAR(IABS(KBeam)).eq.0)
     &        PREM = PAR(140)
      ENDIF           

c     define Prem as probablility for remnant survival
c     switch to sampling of remnant de-excitation
      IF(IPAR(79).ne.0) PREM = ONE-PREM     

c     prob. of remnant excitation target side
      IF(IPAR(79).ne.0) PREM_NUC = ONE-PAR(23)
      IF(IPAR(63).eq.1) PREM_NUC = PREM_NUC/dble(NW)

c     turn of remnant for Nw>1
      select case(ipar(77))
      case(1)
c     only beamside
         IF(NW.gt.1) PREM = 0
      case(2)
c     target and beam-side
         IF(NW.gt.1) then
            PREM = zero
            PREM_NUC = zero
         endif
      case default
         continue
      end select

C...  remnant mass dis. exponents
      XRMEX(1) = PAR(98)        ! baryons
      IF(IPAR(84).gt.0)THEN
         XRMEX(2) = PAR(141)    ! mesons
      else
         XRMEX(2) = PAR(98)     ! mesons same as baryons
      endif
      
      IBMRMNT = 0
      DO K=1, NW
c     additionally penalize remnant survival for multiple nucleon interactions
         IF(IPAR(79).eq.2.and.K.gt.1) PREM=ONE-PAR(23)*PAR(128)
c     penalize remnant survival for multiple parton interactions
         IF(IPAR(80).ne.0) THEN
c     multiple interaction penalty for remnant survival, individual interaction
            ALPH = ONE+PAR(129)*(NNSOF(K)+NNJET(K)-1)
            PREM = ONE-(ONE-PREM)**ALPH
            PREM_NUC = ONE-(ONE-PREM_NUC)**ALPH
         ENDIF
         IF(JDIF(K).eq.0)THEN
            R = S_RNDM(k)
            R2 = S_RNDM(k)
            IF(R.LT.PREM) IBMRMNT = IBMRMNT + 1
c     no target side excitation if recursive call (irec=0)!
            IF(R2.LT.PREM_NUC*Irec) THEN
               ITGRMNT(K) = 1
            ELSE
               ITGRMNT(K) = 0
            ENDIF
         ELSE
            ITGRMNT(K) = 0
         ENDIF
         IF(Ndebug.gt.1) 
     &        WRITE(LUN,'(2X,A,1X,I2,1X,F5.3,1X,I2,1X,I2,1X,I2,1X,I2)')
     &        'sample_RMNT: (JW,PREM,NS,NH,IBMRMNT,LTGRMNT):',
     &        K,PREM,NNSOF(k),NNJET(k),IBMRMNT,ITGRMNT(k)
      ENDDO
      IF(IPAR(79).ne.0)THEN
c     Prem was redefined as probablility for remnant destruction
c     therefore invert configuration..
         DO K=1, NW
            IF(JDIF(K).eq.0)THEN
               ITGRMNT(K)=IABS(ITGRMNT(K)-1)
            ENDIF
         ENDDO
c     multiple de-excitations not possible..
         IBMRMNT=MIN(IBMRMNT,1)
         IBMRMNT=IABS(IBMRMNT-1)
      ENDIF
      IF(Ndebug.gt.1) 
     &     WRITE(LUN,*)
     &     ' sample_RMNT: remnant sampling (PREM,NW,LBMRMNT,LTGRMNT): ',
     &     PREM,NW,IBMRMNT,(ITGRMNT(j),j=1,NW)

      IBMRMNT_OLD = IBMRMNT

C...  Sample flavor and momentum fractions
 20   ITRY(2) = ITRY(2) + 1
c     reset parton stack
      call ini_prtn_stck(NPPLD,NPP0LD)
      IBMRMNT = IBMRMNT_OLD

c     retry without counting
c 22   CONTINUE
      IF(ITRY(2).gt.NREJ(2))THEN
         LBAD = 2
         IF(ndebug.gt.1)then 
            WRITE(LUN,*)' sample_RMNT: number of trials exceeded'
            WRITE(LUN,*)' resample minijets...(IREJ,NW,NCALL)',
     &           LBAD, NW, NCALL
         endif
c     raise event call counter
c         NCALL = NCALL + 1
         RETURN
      ENDIF

      Kbm = Kbeam

C..   sample central strings and remnant flavor
      flvXchg = PAR(80)     ! prob. of flv exchange between strgs and rmnt
c     remnant and sea on beam side
      CALL sample_beam(Kbm,NW,flvXchg,IBMRMNT,X1JET,IKBAD)
      select case (IKBAD)
      case(1)
c     resample minijets event
         LBAD = 3
         RETURN
      case(2)
c     too many partons, reject NW, i.e. entire event
         LBAD = 1
         RETURN
      end select

c     remnants and sea on target side
      CALL sample_target(NW,flvXchg,ITGRMNT,X2JET,Irec,IKBAD)
      select case (IKBAD)
      case(1)
c     resample minijets event
         LBAD = 3
         RETURN
      case(2)
c     too many partons, reject NW, i.e. entire event
         LBAD = 1
         RETURN
      end select

C...  sample remnant excitation masses and add to parton stack
c     beam-side (one remnant, formed by several interactions)
c     target-side (possibly NW remnants)

      DO I=1,NW
c     default return point
         IBAD = 1
         IF(IPAR(78).EQ.1)THEN
c$$$            write(lun,*) 
c$$$     &           'SIB_RMNT: multiple excitation model',
c$$$     &           ' not implemented yet!'
c$$$            stop
c     model where beam side remnant can receive mass from multiple target nucleons
            IF(IBMRMNT.gt.0)THEN
c     beam side remnant excited               
               if(ITGRMNT(I).eq.0)then
                  call exct_rmnt(I,1,IBAD)
               else
                  call exct_rmnt(I,3,IBAD)
               endif
               IBMRMNT = IBMRMNT - 1
            ELSE
c     beam side remnant not excited
               if(ITGRMNT(I).ne.0)then
                  call exct_rmnt(I,2,IBAD)
               else
                  call exct_rmnt(I,0,IBAD)
               endif
            ENDIF

         ELSEIF(IPAR(78).eq.2)then
            IF(IBMRMNT.gt.0)then
c     beam side remnant excited, only once!
               IF(ITGRMNT(I).eq.0)then
                  call exct_rmnt(I,1,IBAD)
               else
                  call exct_rmnt(I,3,IBAD)
               endif
               IBMRMNT = 0
            ELSE
c     beam side remnant not excited
               IF(ITGRMNT(I).ne.0)then
                  call exct_rmnt(I,2,IBAD)
               else
                  call exct_rmnt(I,0,IBAD)
               endif
            ENDIF
         ELSE
c     no remnant model
            call exct_rmnt(I,0,IBAD)
         ENDIF
c     catch remant excitation exception, redo sea kinematics..
         IF(IBAD.eq.1) GOTO 20
c     catch severe exception, resample minijet kinematics..
         IF(IBAD.eq.2) THEN
            LBAD = 3
            RETURN              ! resample event
         ENDIF
      ENDDO
      LBAD = 0
      
      END
