C This file are aux functions for stand-alone use
C of CORSIKA interaction models

      SUBROUTINE CQGSINI( DATDIR, LUN, DEBUGNUM )
C-----------------------------------------------------------------------
C  C(ompact) Q(UARK) G(LUON) S(TRING JET MODEL) INI(TIALZATION)
C
C  INITIALIZES QGSJET MODEL quitely and fully.
C  THIS SUBROUTINE IS CALLED FROM START.
C-----------------------------------------------------------------------
         IMPLICIT NONE
         CHARACTER*256 DATDIR

         COMMON /AREA8/   WWM,BE(4),DC(5),DETA,ALMPT
         DOUBLE PRECISION WWM,BE,DC,DETA,ALMPT

         COMMON /AREA40/  JDIFR
         INTEGER          JDIFR
         INTEGER    LUN, MONIOU
         INTEGER DEBUGNUM, DEBUG
         COMMON /AREA43/ MONIOU
         COMMON /DEBUG/  DEBUG

C  Set output LUN
         MONIOU = LUN

C  Set debug level
         DEBUG = DEBUGNUM

C  TO DISTINGISH QGSJET01D FROM QGSJET01 WE USE DIFFERENT NAMES FOR PSASET
         CALL PSASETC
C  RR   = 0.35D0 ! TO ADJUST SIGMA PP TO 76 MBARN AT TEVATRON/CDF
C  IS SET BY DEFAULT IN XXASET
C  OTHERWISE     RR = 0.53     TO ADJUST SIGMA=80MBARN

C  PARTICULAR MODEL PARAMETERS SETTING
         CALL XXASET

C  SET DIFFRACTION FLAG BY DEFAULT
         JDIFR = 1

C  Call global initialization subroutine
         CALL PSAINI( DATDIR )
C  Call cross-section initialization subroutine
C        CALL QGSSIGINI

         DC(3) = .0003D0   ! To switch off charmed particles set to 0.000
         DC(5) = .01D0     ! To switch off charmed particles set to 0.000

      END

      SUBROUTINE CHEPEVT
C-----------------------------------------------------------------------
C  Convert to HEPEVT common block
C
C-----------------------------------------------------------------------
         IMPLICIT NONE

         INTEGER NPTMAX, ICH, NSP
         DOUBLE PRECISION ESP
         PARAMETER(NPTMAX=95000)
         COMMON /AREA12/ NSP
         COMMON /AREA14/ ESP(4,NPTMAX),ICH(NPTMAX)


         INTEGER NEVHEP,NMXHEP,NHEP,ISTHEP,IDHEP,JMOHEP,JDAHEP
         DOUBLE PRECISION PHEP,VHEP
         PARAMETER (NMXHEP=NPTMAX)
         COMMON /HEPEVT/ NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &           JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),
     &           VHEP(4,NMXHEP)
         INTEGER ICHG
         COMMON /QGCHG/  ICHG(NMXHEP)
C     Particle tables start with the ID -10(rho0) going through 0 (pi0).
         character*12 NAME(-10:10)
         DATA NAME /
     &   'rho0        ','Lambda_cbar-','Dbar0       ','D-          ',
     &   'Lambdabar0  ','K_L0        ','K-          ','nbar0       ',
     &   'pbar-       ','pi-         ','pi0         ','pi+         ',
     &   'p+          ','n0          ','K+          ','K_S0        ',
     &   'Lambda0     ','D+          ','D0          ','Lambda_c+   ',
     &   'eta         '/

         INTEGER IPDGID(-10:10)
         DATA IPDGID /
     &      113, -4122,  -421,  -411, -3122,   130,  -321, -2112, -2212,
     &     -211,   111,   211,  2212,  2112,   321,   310,  3122,   411,
     &      421,  4122,   221/

         DOUBLE PRECISION QMASS(-10:10)
         DATA QMASS /
     &   .548d0,2.27d0,1.868d0,1.868d0,1.116d0,.496d0,.496d0,
     &   0.93827999,0.93827999,.14d0,.14d0,.14d0,0.93827999,
     &   0.93827999,.496d0,.496d0,1.116d0,1.868d0,1.868d0,2.27d0
     &   ,.548d0/
         INTEGER ICHRG(-10:10)
         DATA ICHRG /
     &        0,    -1,     0,    -1,     0,     0,    -1,     0,
     &       -1,    -1,     0,     1,     1,     0,     1,     0,
     &        0,     1,     0,     1,     0/

         INTEGER I

         NHEP = nsp

         DO I=1,nsp
C         WRITE(6,*) I, ich(I), esp(:,I)
            ISTHEP(I) = 1
            NHEP = NSP
            IDHEP(I) = IPDGID(ich(I))
            PHEP(1,I) = esp(3,I)
            PHEP(2,I) = esp(4,I)
            PHEP(3,I) = esp(2,I)
            PHEP(4,I) = esp(1,I)
            PHEP(5,I) = QMASS(ich(I))
            ICHG(I) = ICHRG(ich(I))
         END DO


      END


*-- Author :    D. HECK IK FZK KARLSRUHE       12/01/1996
C=======================================================================

      BLOCK DATA QGSDAT

C-----------------------------------------------------------------------
C  Q(UARK) G(LUON) S(TRING JET MODEL) DAT(A INITIALIZATION)
C
C  Tables for the conversion of QGSJET PIDs to CORSIKA and
C  vise versa
C-----------------------------------------------------------------------

         IMPLICIT NONE

         COMMON /CRQGSLIN/ICTABL,IQTABL
         INTEGER          ICTABL(200),IQTABL(-10:10)
         SAVE
C  FOLLOWING NOTATIONS FOR PARTICLES TYPES ARE USED WITHIN QGSJET:
C             0 - PI0,
C             1 - PI+,
C            -1 - PI-,
C             2 - P,
C            -2 - P-BAR,
C             3 - N,
C            -3 - N-BAR,
C             4 - K+,
C            -4 - K-,
C             5 - K0S,
C            -5 - K0L
C             6 - LAMBDA
C            -6 - LAMBDA-BAR
C             7 - D+
C            -7 - D-
C             8 - D0
C            -8 - D0-BAR
C             9 - LAMBDA_C
C            -9 - LAMBDA_C-BAR
C            10 - ETA

C  ICTABL CONVERTS CORSIKA PARTICLES INTO QGSJET PARTICLES
C  NO CHARMED PARTICLES POSSIBLE AS PROJECTILES
         DATA ICTABL/
     *      0,   0,   0,   0,   0,   0,   1,   1,  -1,  -5,   ! 10
     *      4,  -4,   3,   2,  -2,   5,  10,   6,   0,   0,   ! 20
     *      0,   0,   0,   0,  -3,  -6,   0,   0,   0,   0,   ! 30
     *      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   ! 40
     *      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   ! 50
     *      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   ! 60
     *      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   ! 70
     *     10,  10,  10,  10,   26*0,
C  CHARMED MESONS
C  CONVERT UNKNOWN CHARMED MESONS TO KNOWN D-MESONS
     *     10*0,                                              !110
     *      0,   0,   0,   0,   0,   8,   7,  -7,  -8,   7,   !120
     *     -7,   0,   8,   7,  -7,  -8,   7,  -7,   0,   0,   !130
C  CHARMED BARYONS
     *      0,   0,   0,   0,   0,   0,   9,   0,   0,   0 ,  !140
     *      0,   0,   0,   0,   0,   0,   0,   0,  -9,   0 ,  !150
     *      50*0 /

C  IQTABL CONVERTS QGSJET PARTICLES INTO CORSIKA PARTICLES
C  INCLUDES CHARMED PARTICLES
C  IQTABL RUNS FROM -10:10
         DATA IQTABL/
     *      0, 149, 119, 118,  26,  10,  12,  25,  15,   9,  ! -10 .... -1
     *      7,                                               !   0
     *      8,  14,  13,  11,  16,  18, 117, 116, 137,  17/  !   1 .... 10

      END

*-- Author :    A. Fedynitch ASIoP       14/03/2024
C=======================================================================

      subroutine CQGSHH_HA_CS(gtot,gin,gel,gdp,gdt,gdd)

C-----------------------------------------------------------------------
C Calculate QGSJET cross sections using native functions
C call xxini before calling this function
C-----------------------------------------------------------------------

         implicit double precision (a-h,o-z)
         INTEGER DEBUG
         COMMON /AREA1/  IA(2),ICZ,ICP
         COMMON /AREA2/  S,Y0,WP0,WM0
         COMMON /AREA5/  RD(2),CR1(2),CR2(2),CR3(2)
         COMMON /AREA6/  PI,BM,AM
         COMMON /AREA15/ FP(5),RQ(5),CD(5)
         COMMON /AREA16/ CC(5)
         COMMON /AREA17/ DEL,RS,RS0,FS,ALFP,RR,SH,DELH
         COMMON /AREA43/ MONIOU
         COMMON /AR1/    ANORM
         COMMON /DEBUG/  DEBUG

         dimension gz0(5),gz1(3)
Cf2py    intent(out) :: gtot,gin,gel,gdp,gdt,gdd

         if (ia(1).ne.1) then
            write(*,*) 'Error: CQGSHH_HA_CS called with ia(1) ne 1'
            return
         endif
         ! icz = hadron class (N, pi, K)
c Energy dependent factors:
c Type of the incident hadron (icz = 1: pion, 2: nucleon, 3: kaon, etc
c RS - soft pomeron elastic scattering slope (lambda_ab)
         RS=RQ(ICZ)+ALFP*Y0
c RS0 - initial slope (sum of the pomeron-hadron vertices slopes squared - R_ab)
         RS0=RQ(ICZ)
c RP1 - factor for the impact parameter dependence of the eikonal ( in fm^2 )
         RP1=RS*4.D0*.0391D0/AM**2
c Factor for cross-sections calculation ( in mb )
         G0=PI*RP1/CD(ICZ)*AM**2*10.D0
c SJV - valence-valence cross-section (divided by 8*pi*lambda_ab)
         Z=.2D0*I
         IF(IA(2).EQ.1)THEN
c Hadron-proton interaction
c BM - impact parameter cutoff value
            BM=2.D0*DSQRT(RP1)
c XXFZ - impact parameter integration for the hadron-nucleon interaction eikonal;
c GZ0 - total and absorptive cross-sections (up to a factor); first parameter is
c used only in case of hadron-nucleus interaction (to make convolution with target
c nucleus profile function)
            CALL XXFZ(0.D0,GZ0)
            if (debug .ge.1) write (moniou,*) gz0, g0, Y0, RS0, RP1,
     *         CD(icz), icz, CD
c GTOT - total cross-section
            GTOT=G0*GZ0(1)
c GABS - cut pomerons cross-section
            GABS=G0*GZ0(2)*.5D0
c GD0 - cross-section for the cut between pomerons
            GD0=GTOT-GABS
c GDP - projectile diffraction cross section
            GDP=(1.D0-CC(ICZ))*CC(2)*GD0
c GDT - target diffraction cross section
            GDT=(1.D0-CC(2))*CC(ICZ)*GD0
c  GDD - double diffractive cross section
            GDD=(1.D0-CC(ICZ))*(1.D0-CC(2))*GD0
c GIN - inelastic cross section
            GIN=GABS+GDP+GDT+GDD
            GEL=GD0*CC(ICZ)*CC(2)
c
            IF(DEBUG.GE.1)WRITE (MONIOU,225)GTOT,GIN,GEL,GDP,GDT,GDD
c
  225       FORMAT(2X,'PSAINI: HADRON-PROTON CROSS SECTIONS:',/,
     *      4X,'GTOT=',E10.3,2X,'GIN=',E10.3,2X,'GEL=',E10.3,/,4X,
     *      'GDIFR_PROJ=',E10.3,2X,'GDIFR_TARG=',E10.3,2X,
     *      'G_DOUBLE_DIFR',E10.3)

         ELSE

c Hadron-nucleus interaction
c BM - impact parameter cutoff value
            BM=RD(2)+DLOG(29.D0)
c RRR - Wood-Saxon radius for the target nucleus
            RRR=RD(2)
c RRRM - auxiliary parameter for numerical integration
            RRRM=RRR+DLOG(9.D0)
c ANORM - nuclear density normalization factor multiplied by RP1
            ANORM=1.5D0/PI/RRR**3/(1.D0+(PI/RRR)**2)*RP1

c GAU(GZ) - cross sections calculation ( integration over impact parameters less than
c BM )
            CALL XXGAU(GZ1)
c GAU1(GZ) - cross sections calculation ( integration over impact
c parameters greater than BM )
            CALL XXGAU1(GZ1)
c GIN - total inelastic cross section
            GIN=GZ1(1)+GZ1(2)+GZ1(3)
c
            IF(DEBUG.GE.1)WRITE (MONIOU,224)
     *      GIN*10.D0,GZ1(1)*10.D0,GZ1(2)*10.D0
c
  224       FORMAT(2X,'PSAINI: HADRON-NUCLEUS CROSS SECTIONS:',/,
     *      4X,'GIN=',E10.3,2X,'GDIFR_TARG=',E10.3,2X,
     *      'GDIFR_PROJ=',E10.3)
c GZ - probability to have target diffraction
            gdt=GZ1(1)*10.d0
            gdp=GZ1(2)*10.d0
            GIN=GIN*10.d0
         ENDIF
      end subroutine
