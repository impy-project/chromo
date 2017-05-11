      SUBROUTINE SIB_SIGMA_EXT
     &     (L0,SQS,SIGT,SIGEL,SIGINEL,SLOPE,RHO)
C-----------------------------------------------------------------------
C     Hadron-proton cross sections
C     taken from EXTERNAL(!) interpolation table (calculated elsewhere)
C     can be used to run NUCLIB with alternative cross section/int.length
C
C     input:       L     1,2,3      proton-,pion-,kaon-proton
C                  SQS   sqrt(s)
C
C     output:      SIGT       total cross section (mb)
C                  SIGEL      elastic cross section (mb)
C                  SIGINEL    inelastic cross section (mb)
C                  SLOPE      elastic slope parameter (GeV^-2)
C                  RHO        real/imaginary part of forward amplitude
C-----------------------------------------------------------------------
Cf2py integer, intent(in) :: L0
Cf2py double precision, intent(in) :: SQS
Cf2py double precision, intent(out) :: SIGT,SIGEL,SIGINEL,SLOPE,RHO
      IMPLICIT NONE
      SAVE

c     external types
      DOUBLE PRECISION SQS,SIGT,SIGEL,SIGINEL,SLOPE,RHO
      INTEGER L0

      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'

c     external cross section tables
      INCLUDE 'sib_xsctn_21_data.inc'
      
c     internal type declarations
      INTEGER LL,L,J1,NSQS
      DIMENSION LL(39)
      DATA LL /5*0,3*2,4*3,2*1,19*0,6*1/
      DOUBLE PRECISION T,AL,ASQSMIN,ASQSMAX,DASQS

      L = L0
      NSQS = 61
      ASQSMIN = ONE
      ASQSMAX = 7.D0
      DASQS = (ASQSMAX-ASQSMIN)/DBLE(NSQS-1)

      IF(NSQS.LE.0) THEN
        WRITE(LUN,'(//,1X,A)') 
     &    'SIB_SIGMA_EXT: interpolation table not initialized.'
        STOP
      ENDIF
      IF(IABS(L).gt.39)THEN
         WRITE(LUN,*)     
     &        'SIB_SIGMA_EXT: unknown beam particle!',L
         STOP
      ENDIF
      IF(L.GT.3) L=LL(IABS(L))
      IF(L.EQ.0)THEN
         WRITE(LUN,*)     
     &        'SIB_SIGMA_EXT: unknown beam particle!', L
         STOP
      ENDIF
        
      AL = LOG10(SQS)
      J1 = INT((AL-ONE)*10.D0 + 1)
      if((j1.lt.1).or.(j1.ge.NSQS)) then
        write (LUN,'(1x,a,i3,1p,e12.3)') 
     &    'SIB_SIGMA_HP: energy out of range ',L,sqs
        J1 = min(J1,NSQS-1)
        J1 = max(J1,1)
      endif
      T = (AL-ONE)*10.D0 - DBLE(J1-1)
      SIGT    = SSIG_TOT(J1,L)*(ONE-T) + SSIG_TOT(J1+1,L)*T
      SIGINEL = SSIG(J1,L)*(ONE-T) + SSIG(J1+1,L)*T
      SIGEL   = SIGT-SIGINEL
      SLOPE   = SSIG_B(J1,L) *(ONE-T) + SSIG_B(J1+1,L)*T
      RHO     = SSIG_RHO(J1,L) *(ONE-T) + SSIG_RHO(J1+1,L)*T

      END
