      SUBROUTINE SIB_REJECT
c     subroutine dumps state of random number generator 
c     at beginning of event to file then produces fpe/stops
C----------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_rand_cmmn.inc'
      INCLUDE 'sib_nw_prm.inc'
      INCLUDE 'sib_run_cmmn.inc'
      DOUBLE PRECISION XDM
      CHARACTER *(13) FILENA
      DATA FILENA /'sib_rjctn.rnd'/
      WRITE(LUN,*)'SIB_reject: (ncall,KB,iat,ECM)',ncall,kb,iat,sqs
c     restore state before event
!       all pho_rndsi(U2,C2,CD2,CM2,II2,JJ2)
c     dump state to file
!       all PHO_RNDST(2,FILENA)
c     produce floating point error
      XDM = -1.D0
      XDM = LOG(XDM)
      STOP
      END
