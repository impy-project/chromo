      SUBROUTINE SIB_SIGMA_HP
     &     (L0,SQS,SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO)
C-----------------------------------------------------------------------
C     Hadron-proton cross sections, taken from interpolation table
C     calculated by SIBYLL_INI
C
C     input:       L     1      proton-proton
C                        2      pi-proton
C                        3      K-proton
C                  SQS   sqrt(s)
C
C     output:      SIGT       total cross section (mb)
C                  SIGEL      elastic cross section (mb)
C                  SIGINEL    inelastic cross section (mb)
C                  SIGDIF     diffraction dissociation CS (mb)
C                  SLOPE      elastic slope parameter (GeV^-2)
C                  RHO        real/imaginary part of forward amplitude
C-----------------------------------------------------------------------
Cf2py integer, intent(in) :: L0
Cf2py double precision, intent(in) :: SQS
Cf2py double precision, intent(out) :: SIGT,SIGEL,SIGINEL,SLOPE,RHO
Cf2py double precision(3), intent(out) :: SIGDIF
      IMPLICIT NONE
      SAVE

c     external types
      DOUBLE PRECISION SQS,SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO      
      DIMENSION SIGDIF(3)
      INTEGER L0

      INCLUDE 'sib_debug_cmmn.inc'
      include 'sib_int_prm.inc'
      INCLUDE 'sib_ccsig_cmmn.inc'
      INCLUDE 'sib_ccsig2_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'

c     internal type declarations
      INTEGER LL,L,J1
      DIMENSION LL(39)
      DATA LL /5*0,3*2,4*3,2*1,19*0,6*1/

      DOUBLE PRECISION T,AL

      L = L0
      IF(NSQS.LE.0) THEN
        WRITE(LUN,'(//,1X,A)') 
     &    'SIB_SIGMA_HP: interpolation table not initialized.'
        STOP
      ENDIF
      IF(IABS(L).gt.39)THEN
         WRITE(LUN,*)     
     &        'SIB_SIGMA_HP: unknown beam particle!',L
         STOP
      ENDIF
      IF(L.GT.3) L=LL(IABS(L))
      IF(L.EQ.0)THEN
         WRITE(LUN,*)     
     &        'SIB_SIGMA_HP: unknown beam particle!', L
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
      SIGDIF(1) = SSIG_SD1(J1,L)*(ONE-T) + SSIG_SD1(J1+1,L)*T
      SIGDIF(2) = SSIG_SD2(J1,L)*(ONE-T) + SSIG_SD2(J1+1,L)*T
      SIGDIF(3) = SSIG_DD(J1,L)*(ONE-T) + SSIG_DD(J1+1,L)*T
      SLOPE   = SSIG_B(J1,L) *(ONE-T) + SSIG_B(J1+1,L)*T
      RHO     = SSIG_RHO(J1,L) *(ONE-T) + SSIG_RHO(J1+1,L)*T

      END


      SUBROUTINE SIB_SIGMA_HP2
     +     (L,SQS,SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO)
C-----------------------------------------------------------------------
C     Hadron-proton cross sections, taken from interpolation table
C     calculated by SIBYLL_INI
C
C     input:       L     1      proton-proton
C                        2      pi-proton
C                        3      K-proton
C                  SQS   sqrt(s)
C
C     output:      SIGT       total cross section (mb)
C                  SIGEL      elastic cross section (mb)
C                  SIGINEL    inelastic cross section (mb)
C                  SIGDIF     diffraction dissociation CS (mb)
C                             split in high and low mass !!
C                             ( taken from S_CCSIG3 )
C                  SLOPE      elastic slope parameter (GeV^-2)
C                  RHO        real/imaginary part of forward amplitude
C-----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
c     external types      
      DOUBLE PRECISION SQS,SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO
      DIMENSION SIGDIF(3,2)
      INTEGER L

      INCLUDE 'sib_debug_cmmn.inc'
      include 'sib_int_prm.inc'
      INCLUDE 'sib_ccsig_cmmn.inc'
      INCLUDE 'sib_ccsig2_cmmn.inc'
      INCLUDE 'sib_ccsig3_cmmn.inc'

c     internal types
      INTEGER J1
      DOUBLE PRECISION T,AL

      IF(NSQS.LE.0) THEN
        WRITE(LUN,'(//,1X,A)') 
     &    'SIB_SIGMA_HP2: interpolation table not initialized.'
        STOP
      ENDIF
        
      AL = dLOG10(SQS)
      J1 = INT((AL - 1.D0)*10.D0 + 1)
      if((j1.lt.1).or.(j1.ge.NSQS)) then
        write (LUN,'(1x,a,i3,1p,e12.3)') 
     &    'SIB_SIGMA_HP2: energy out of range ',L,sqs
        J1 = min(J1,NSQS-1)
        J1 = max(J1,1)
      endif
      T = (AL-1.D0)*10.D0 - DBLE(J1-1)
      SIGT    = SSIG_TOT(J1,L)*(1.D0-T) + SSIG_TOT(J1+1,L)*T
      SIGINEL = SSIG(J1,L)*(1.D0-T) + SSIG(J1+1,L)*T
      SIGEL   = SIGT-SIGINEL
      SIGDIF(1,1) = SSIG_SD1LM(J1,L)*(1.D0-T) + SSIG_SD1LM(J1+1,L)*T
      SIGDIF(1,2) = SSIG_SD1HM(J1,L)*(1.D0-T) + SSIG_SD1HM(J1+1,L)*T
      SIGDIF(2,1) = SSIG_SD2LM(J1,L)*(1.D0-T) + SSIG_SD2LM(J1+1,L)*T
      SIGDIF(2,2) = SSIG_SD2HM(J1,L)*(1.D0-T) + SSIG_SD2HM(J1+1,L)*T
      SIGDIF(3,1) = SSIG_DDLM(J1,L)*(1.D0-T) + SSIG_DDLM(J1+1,L)*T
      SIGDIF(3,2) = SSIG_DDHM(J1,L)*(1.D0-T) + SSIG_DDHM(J1+1,L)*T
      SLOPE   = SSIG_B(J1,L) *(1.D0-T) + SSIG_B(J1+1,L)*T
      RHO     = SSIG_RHO(J1,L) *(1.D0-T) + SSIG_RHO(J1+1,L)*T

      END


      SUBROUTINE SIB_SIGMA_HAIR (L,SQS,SIGprod,SIGbdif) 
C-----------------------------------------------------------------------
C     Hadron-air cross sections, taken from interpolation table
C     calculated by SIBYLL_INI
C
C     input:       L     1      proton-air
C                        2      pi-air
C                        3      K-air
C                  SQS   sqrt(s)
C
C     output:      SIGprod    particle production cross section (mb)
C                  SIGbdif    q.ela and ela beam diff. cross section
C-----------------------------------------------------------------------
Cf2py integer, intent(in) :: L
Cf2py double precision, intent(in) :: SQS
Cf2py double precision, intent(out) :: SIGprod,SIGbdif
      IMPLICIT NONE
      SAVE
      INCLUDE 'sib_debug_cmmn.inc'
      include 'sib_int_prm.inc'
      INCLUDE 'sib_ccsig_cmmn.inc'

c     external
      DOUBLE PRECISION SQS,SIGPROD,SIGBDIF
      INTEGER L

c     internal
      DOUBLE PRECISION AL,T
      INTEGER J1

      IF(NSQS.LE.0) THEN
        WRITE(LUN,'(//,1X,A)') 
     &    'SIB_SIGMA_HAIR: interpolation table not initialized.'
        STOP
      ENDIF
        
      AL = LOG10(SQS)
      J1 = INT((AL - 1.D0)*10.D0 + 1)
      if((j1.lt.1).or.(j1.ge.NSQS)) then
        write (LUN,'(1x,a,i3,1p,e12.3)') 
     &    'SIB_SIGMA_HAIR: energy out of range ',L,sqs
        J1 = min(J1,NSQS-1)
        J1 = max(J1,1)
      endif
      T = (AL-1.D0)*10.D0 - DBLE(J1-1)
      SIGprod = SSIGN(J1,L)*(1.D0-T) + SSIGN(J1+1,L)*T
      SIGbdif = SSIGNSD(J1,L)*(1.D0-T) + SSIGNSD(J1+1,L)*T

      END

      SUBROUTINE SIB_SIGMA_HNUC (L,IAT,SQS,SIGprod,SIGbdif) 
C-----------------------------------------------------------------------
C     calculate Hadron-nucleus cross sections
C
C     input:       L     1      proton-nuc
C                        2      pi-nuc
C                        3      K-nuc
C                  IAT   0-18   mass number of target nucleus
C                        (beyond A=18 nuclear profiles are inaccurate)
C                  SQS   sqrt(s)
C
C     output:      SIGprod    particle production cross section (mb)
C                  SIGbdif    q.ela and ela beam diff. cross section
C-----------------------------------------------------------------------
Cf2py integer, intent(in) :: L,IAT
Cf2py double precision, intent(in) :: SQS
Cf2py double precision, intent(out) :: SIGprod,SIGbdif
      IMPLICIT NONE
      SAVE

      include 'sib_int_prm.inc'
      INCLUDE 'sib_ccsig_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'
      DOUBLE PRECISION SIGT,SIGEL,SIGINEL,SIGQE,SIGSD,SIGQSD,SIGPPT,
     &     SIGPPEL,SIGPPSD
      INTEGER ITG
      COMMON /NUCSIG/ SIGT,SIGEL,SIGINEL,SIGQE,SIGSD,
     +     SIGQSD,SIGPPT,SIGPPEL,SIGPPSD,ITG

c     external
      DOUBLE PRECISION SQS,SIGPROD,SIGBDIF
      INTEGER L,IAT

c     internal
      DOUBLE PRECISION ALAM
      INTEGER IPARM,ICSMOD

      IF(NSQS.LE.0) THEN
        WRITE(6,'(//,1X,A)') 
     &    'SIB_SIGMA_HNUC: interpolation table not initialized.'
        STOP
      ENDIF

      IF(IAT.ge.0.and.IAT.lt.19)THEN
c     calculate hadron - nucleus cross section
c     dummy arg, coupling derived from dif xsctn
         ALAM = ONE              
c     use Sibyll p-p cross section as input
         ICSMOD = 1             
c     use Goulianos param. for inel. coupling param.
         IPARM = 2 
         CALL SIG_HAD_NUC(L,IAT,SQS,ALAM,ICSMOD,IPARM)
C     particle production cross section        
         SIGprod = SIGT-SIGQE
C     quasi elastic + elastic singl. diff cross section
         SIGbdif = SIGQSD
      ELSE
         WRITE(6,'(//,1X,A)') 
     &     'SIB_SIGMA_HNUC: number of target nucleons too large!',
     &     ' (0<=IAT<=18)'
         SIGprod = -ONE
      ENDIF
      RETURN
      END
