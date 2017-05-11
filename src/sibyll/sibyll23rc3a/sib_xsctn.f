
      SUBROUTINE CUT_PRO (L, SQS, PTmin, NSOFR, NJETR)
C-----------------------------------------------------------------------
C...  Generate a number of soft/hard (jet-)pairs for a 'projectile'
C     (K=1:p),(K=2:pi) interacting with a nucleon at sqrt(s)=SQS(GeV)
C     the interaction structure is only destinguished between nucleons
C     (L=1) and mesons (L=2), for cross sections there is a 
C     distinction between pions and kaons as well (L=2 or 3).
C     For Hyperons the same cross section and interaction structure
C     as for nucleons is used (L=1).
C
C     requires initialization by JET_INI                         /FR'14
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      include 'sib_debug_cmmn.inc'
c      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      include 'sib_cflafr_cmmn.inc'
      include 'sib_int_prm.inc'
      include 'sib_ccsig_cmmn.inc'
c      PARAMETER (NS_max = 20, NH_max = 80)
c      COMMON /S_CCSIG/ SSIG(61,3), PJETC(0:NS_max,0:NH_max,61,2),
c     &     SSIGN(61,3), SSIGNSD(61,3), ALINT(61,3),
c     &     ASQSMIN, ASQSMAX, DASQS, NSQS
c      COMMON /S_CUTOFF/ STR_mass_val, STR_mass_val_hyp, STR_mass_sea
      INCLUDE 'sib_cutoff_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'      

C     check if tables initialized
      IF(NSQS.eq.0) THEN
         WRITE(LUN,*) 'CUT_PRO: tables not initialized! aborting...'
         xa = -ONE
         xa = log(xa)
         stop
      ENDIF
      IF(NDEBUG.GT.1) 
     &     WRITE(LUN,*) ' CUT_PRO: input: L, SQS, PTmin',L, SQS, PTmin

c     choose nucleon or meson table
      K = L
      if(K.eq.3) K = 2

      AL = dLOG10 (SQS)
      IF (AL .LT. ASQSMIN)  THEN
          WRITE(LUN,*)  ' CUT_PRO:  low sqrt(s) ', SQS
          NSOFR = 1
          NJETR = 0
          RETURN
      ENDIF
      IF (AL .GT. ASQSMAX)  THEN
          WRITE(LUN,*)  ' CUT_PRO:  sqrt(s) out of bounds ', SQS
          NJETR = 0
          RETURN
      ENDIF

      J1 = INT((AL - ASQSMIN)/DASQS + 1)
      J1 = MIN(J1,60)
      J1 = MAX(J1,1)
      J2 = J1+1
      T = (AL-ASQSMIN)/DASQS - DBLE(J1-1)

      R = (ONE-EPS8)*S_RNDM(0)
      DO I=0,NS_max
        DO J=0,NH_max
          IF (R.LT.(ONE-T)*PJETC(I,J,J1,K)+T*PJETC(I,J,J2,K)) GOTO 100
        ENDDO
      ENDDO
100   CONTINUE

C...phase space limitation

 120  CONTINUE
      XM = DBLE(2*I)*STR_mass_sea + DBLE(2*J)*PTmin
      PACC = EXP(PAR(9)*(TWO-XM)/SQS)
      IF(S_RNDM(0).GT.PACC) THEN
        IF(I+J.GT.1) THEN
          IF(I.GT.0) THEN
            I = I-1
            GOTO 120
          ELSE IF(J.GT.0) THEN
            J = J-1
            GOTO 120
          ENDIF
        ENDIF
      ENDIF

      NSOFR = I
      NJETR = J

      if(Ndebug.gt.1) 
     &  write(lun,*)' CUT_PRO: (L,SQS,PTmin,Ns,Nh)',K,SQS,PTmin,I,J

      END


      SUBROUTINE JET_INI
C-----------------------------------------------------------------------
C...Compute table of cross sections, and table of probability
C.  for the production of multiple soft and hard interactions
C.
C.  The output of this routine  is the COMMON block /S_CCSIG/
C.  that contains  the cross sections h-p, h-Air, and the 
C.  cumulative probability of NS soft and NH hard interactions
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      include 'sib_int_prm.inc'      
      include 'sib_ccsig_cmmn.inc'
      include 'sib_ccsig2_cmmn.inc'
      include 'sib_ccsig3_cmmn.inc'
c$$$      PARAMETER (NS_max = 20, NH_max = 80)
c$$$      COMMON /S_CCSIG/ SSIG(61,3), PJETC(0:NS_max,0:NH_max,61,2),
c$$$     &     SSIGN(61,3), SSIGNSD(61,3), ALINT(61,3),
c$$$     &     ASQSMIN, ASQSMAX, DASQS, NSQS
c$$$      COMMON /S_CCSIG2/ SSIG_TOT(61,3),SSIG_SD1(61,3),SSIG_SD2(61,3),
c$$$     &    SSIG_DD(61,3),SSIG_B(61,3),SSIG_RHO(61,3)
c$$$
c$$$      COMMON /S_CCSIG3/ SSIG_SD1LM(61,3),SSIG_SD1HM(61,3),
c$$$     &     SSIG_SD2LM(61,3),SSIG_SD2HM(61,3),
c$$$     &     SSIG_DDLM(61,3),SSIG_DDHM(61,3)

c      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      INCLUDE 'sib_debug_cmmn.inc'

      DIMENSION Pjet(0:NS_max,0:NH_max)
      DIMENSION SIG_df(3),SIG_df2(3,2),SIGDIF(3),SIGDIF_pi(3),
     &          PS_tab(61),PH_tab(61),PT_tab(61)
      INCLUDE 'sib_utl_cmmn.inc'

C...spacing in energy for table of cross sections.

      NSQS = 61
      ASQSMIN = ONE
      ASQSMAX = 7.D0
      DASQS = (ASQSMAX-ASQSMIN)/DBLE(NSQS-1)

C...initialization of proton and pion tables
      
      IF(LUN.ne.6) WRITE(6,*)'Calculating cross section tables...'
      DO KK=1,2

         IF(NDEBUG.ge.0)
     &    WRITE(LUN,'(2(/,1X,A,A))') 
     &     'Table: J, sqs,  PT_cut,  SIG_tot,  SIG_inel,  B_el,  ',
     &     'rho,  <n_s>,  <n_h>,  SIG_SD,  SD1_lm,  SD1_hm',
     &     '-----------------------------------------------------',
     &     '----------------------------------------------'

         JINT = KK
         DO J=1, NSQS
           ASQS = ASQSMIN + DASQS*DBLE(J-1)
           SQS = 10.D0**ASQS

           CALL SIB_SIG (JINT, SQS, PTmin,
     &                   SIG_tot, SIG_inel, SIG_df, SIG_df2, B_el, Pjet)

C...low-energy interpolation with data-parametrizations
           call SIB_HADCSL(JINT,SQS,
     &                     SIGTOT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO)
           if(SQS.le.100.D0) then
             SIG_TOT  = SIGTOT
             SIG_inel = SIGINEL
             B_EL     = SLOPE
           else if(SQS.le.1000.D0) then
             Xi = dlog(SQS/100.D0)/2.30258509299405D0
             SIG_TOT  = Xi*SIG_TOT+(ONE-Xi)*SIGTOT
             SIG_inel = Xi*SIG_inel+(ONE-Xi)*SIGINEL
             B_EL     = Xi*B_EL+(ONE-Xi)*SLOPE
           endif

           SSIG_TOT(J,KK) = SIG_TOT
           SSIG(J,KK)     = SIG_inel
           SSIG_SD1(J,KK) = SIGDIF(1)
           SSIG_SD2(J,KK) = SIGDIF(2)
           SSIG_DD(J,KK)  = SIG_df(3)
           SSIG_B(J,KK)   = B_EL
           SSIG_RHO(J,KK) = RHO

           SSIG_SD1LM(J,KK) = SIG_df2(1,1)
           SSIG_SD1HM(J,KK) = SIG_df2(1,2)
           SSIG_SD2LM(J,KK) = SIG_df2(2,1)
           SSIG_SD2HM(J,KK) = SIG_df2(2,2)
           SSIG_DDLM(J,KK) = SIG_df2(3,1)
           SSIG_DDHM(J,KK) = SIG_df2(3,2)

           PSUM = ZERO
           PH = ZERO
           PS = ZERO
           DO NS=0,NS_max
             DO NJ=0,NH_max

               PS = PS+DBLE(NS)*Pjet(NS,NJ)
               PH = PH+DBLE(NJ)*Pjet(NS,NJ)

               PSUM = PSUM+Pjet(NS,NJ)
               PJETC(NS,NJ,J,KK) = PSUM

             ENDDO
           ENDDO
           PS_tab(J) = PS
           PH_tab(J) = PH
           PT_tab(J) = PTmin

           IF(NDEBUG.ge.0)
     &      WRITE(LUN,'(3X,I2,1P,E12.3,0P,4F8.2,6F8.3)') 
     &       JINT,SQS,PTmin,SIG_tot,SIG_inel,B_el,RHO,PS,PH
     &          ,SIGDIF(1)+SIGDIF(2),SIG_df2(1,1),SIG_df2(1,2)

         ENDDO
      ENDDO

C...initialization of kaon tables

      JINT = 3

      IF(NDEBUG.ge.0)
     & WRITE(LUN,'(2(/,1X,A,A))') 
     &  'Table: J, sqs,  PT_cut,  SIG_tot,  SIG_inel,  B_el,  ',
     &  'rho,  <n_s>,  <n_h>',
     &  '-----------------------------------------------------',
     &  '-------------------'
      DO J=1, NSQS
        ASQS = ASQSMIN + DASQS*DBLE(J-1)
        SQS = 10.D0**ASQS
C...use pion cross section rescaled for high-energy extrapolation
        SIG_tot   = SSIG_TOT(J,2)
        SIG_inel  = SSIG(J,2)
        SIG_df(1) = SSIG_SD1(J,2)
        SIG_df(2) = SSIG_SD2(J,2)
        SIG_df(3) = SSIG_DD(J,2)
        B_el = SSIG_B(J,2)
        PTmin = PT_tab(J)
        PS = PS_tab(J)
        PH = PH_tab(J)

C...low-energy interpolation with data-parametrizations
        call SIB_HADCSL(2,SQS,
     &                  SIGTOT_pi,SIGEL_pi,SIGINEL,SIGDIF_pi,SLOPE,RHO)
        call SIB_HADCSL(3,SQS,
     &                  SIGTOT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO)
        SIG_el    = (SIGEL/SIGEL_pi)*(SIG_TOT-SIG_inel)
        SIG_TOT   = (SIGTOT/SIGTOT_pi)*SIG_TOT
        SIG_inel  = SIG_TOT-SIG_el
        SIG_df(3) = (SIGDIF(3)/SIGDIF_pi(3))*SIG_df(3)
        if(SQS.le.100.D0) then
          SIG_TOT  = SIGTOT
          SIG_inel = SIGINEL
          B_EL     = SLOPE
        else if(SQS.le.1000.D0) then
          Xi = dlog(SQS/100.D0)/2.30258509299405D0
          SIG_TOT  = Xi*SIG_TOT+(ONE-Xi)*SIGTOT
          SIG_inel = Xi*SIG_inel+(ONE-Xi)*SIGINEL
          B_EL     = Xi*B_EL+(ONE-Xi)*SLOPE
        endif

        SSIG_TOT(J,3) = SIG_TOT
        SSIG(J,3)     = SIG_inel
        SSIG_SD1(J,3) = SIGDIF(1)
        SSIG_SD2(J,3) = SIGDIF(2)
        SSIG_DD(J,3)  = SIG_df(3)
        SSIG_B(J,3)   = B_EL
        SSIG_RHO(J,3) = RHO

        IF(NDEBUG.ge.0)
     &   WRITE(LUN,'(3X,I2,1P,E12.3,0P,4F8.2,3F8.3)') 
     &    JINT,SQS,PTmin,SIG_tot,SIG_inel,B_el,RHO,PS,PH

      ENDDO


      END


      SUBROUTINE INI_WRITE (LUN)
C-----------------------------------------------------------------------
C   This subroutine prints on unit LUN
C   a table of the cross sections  used in the program
C   and of the average number of hard interactions, and the average
C   number of wounded nucleons in a hadron-air interaction
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      include 'sib_int_prm.inc'
      include 'sib_ccsig_cmmn.inc'
c$$$      PARAMETER (NS_max = 20, NH_max = 80)
c$$$      COMMON /S_CCSIG/ SSIG(61,3), PJETC(0:NS_max,0:NH_max,61,2),
c$$$     &     SSIGN(61,3), SSIGNSD(61,3), ALINT(61,3),
c$$$     &     ASQSMIN, ASQSMAX, DASQS, NSQS
      DIMENSION PJ(2),PS(2),PW(2)

      DATA ATARG /14.514D0/

*      CALL PARAM_PRINT(LUN)
      WRITE (LUN, 10)
      WRITE (LUN, 15)
      WRITE (LUN, 16)
      WRITE (LUN, 18)
10    FORMAT(//,' Table of cross sections, and average number',
     &         ' of minijets and wounded nucleons ')
15    FORMAT('        [sqrt(s) in GeV, cross sections in mbarn]. ')
16    FORMAT(' sqrt(s) sig(pp) sig(pA) <n_s> <n_j> <n_w>',
     &    ' sig(pip) sig(piA) <n_s> <n_j> <n_w>')
18    FORMAT(1X,77('-') )
      DO J=1,61,1
         SQS = 10.D0**(ASQSMIN + DASQS*DBLE(J-1))

         DO K=1,2

           PW(K) = ATARG*SSIG(J,K)/SSIGN(J,K)

           PJ(K) = 0.D0
           PS(K) = 0.D0
           DO NS=0,NS_max
             DO NJ=0,NH_max
               IF(NJ.GT.0) THEN
                 PROB = PJETC(NS,NJ,J,K) - PJETC(NS,NJ-1,J,K)
               ELSE IF(NS.GT.0) THEN
                 PROB = PJETC(NS,NJ,J,K) - PJETC(NS-1,NH_max,J,K)
               ELSE
                 PROB = 0.D0
               ENDIF
               PJ(K) = PJ(K)+DBLE(NJ)*PROB
               PS(K) = PS(K)+DBLE(NS)*PROB
             ENDDO
           ENDDO

         ENDDO

         WRITE(LUN,20) SQS,SSIG(J,1),SSIGN(J,1),PS(1),PJ(1),PW(1)
     &                      ,SSIG(J,2),SSIGN(J,2),PS(2),PJ(2),PW(2)

      ENDDO
      WRITE (LUN, 18)
20    FORMAT (1X,E8.2, 2(2F7.1,1X,3F6.2,1X))

      END


      SUBROUTINE SIG_AIR_INI 
C-----------------------------------------------------------------------
C...Initialize the cross section and interaction lengths on air
C.  (this version initializes p-air, pi-air, and K-air cross sections)
C.
C.  also calculates the low mass beam diffraction cross section in hAir \FR
C.  using the same lambda for all hadrons
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      include 'sib_int_prm.inc'
      include 'sib_ccsig_cmmn.inc'
c$$$      PARAMETER (NS_max = 20, NH_max = 80)
c$$$      COMMON /S_CCSIG/ SSIG(61,3), PJETC(0:NS_max,0:NH_max,61,2),
c$$$     &     SSIGN(61,3), SSIGNSD(61,3), ALINT(61,3),
c$$$     &     ASQSMIN, ASQSMAX, DASQS, NSQS
      COMMON /GLAUB_SCR/ XI_MAX , ALAM(61)
      include 'sib_debug_cmmn.inc'      
      include 'sib_cflafr_cmmn.inc'
c$$$      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      DIMENSION SIGDIF(3)

      DATA AVOG /6.0221367D-04/

      ATARGET = 14.514D0

c      PRINT *,'inel. screening in hadron - nucleus interactions'
c      ALAM = 0.5
c      PRINT *,'const. coupling: ', ALAM
      IF ( IPAR(12).GT.0 ) THEN
         WRITE(LUN,*) 'SIG_AIR_INI:'
         WRITE(LUN,*)'using Goulianos param. for res.coupling..'
         XI_MAX = 0.02D0
         WRITE(LUN,*)'low mass Xi_max: ' , XI_MAX
      ENDIF

C...particle loop (p, pi, K)
      DO K=1,3
         
         WRITE(LUN,'(2(/,1X,A,A))') 
     &        'Table: J, sqs,  SIGtot,   SIGprod,  SIG_SD,  Lambda'
         WRITE(LUN,*) 
     &      '---------------------------------------------------'

        DO J=1,NSQS

           ASQS = ASQSMIN + DASQS*DBLE(J-1)
           SQS = 10.D0**ASQS

           IF (K.EQ.1) THEN
c     Goulianos param. from GAP-2012-056, Mx**2s = 0.02
c     against PDG elastic cross section
              CALL SIB_HADCS1(K,SQS,SIGT1,SIGEL1,SIGINEL1,SLOPE1,RHO1)
              SIGEFF = 0.68D0*(1.D0+36.D0/SQS**2)*
     &             dlog(0.6D0+XI_MAX/1.5D0*SQS**2)
              ALAM(J) = dSQRT(SIGEFF/SIGEL1)
           ENDIF

c           call SIB_HADCSL(k,SQS,
c     &          SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO)

           CALL SIB_SIGMA_HP(K,SQS,SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO)
           CALL SIG_H_AIR2
     &          (SIGT, SLOPE, RHO, ALAM(J),
     &          SSIGT, SSIGEL, SSIGQE, SSIGSD, SSIGQSD)
           IF(IPAR(12).EQ.0)
     &          CALL SIG_H_AIR(SIGT, SLOPE, RHO, SSIGT, SSIGEL, SSIGQE)

           WRITE(LUN,'(1X,1P,5E12.3)') 
     &          SQS,SSIGT,SSIGT-SSIGQE,SSIGQSD,ALAM(J)
C  particle production cross section
           SSIGN(J,K) = SSIGT-SSIGQE
           SSIGNSD(J,K) = SSIGQSD
           ALINT(J,K) = 1.D0/(AVOG*SSIGn(j,K)/ATARGET)
        ENDDO
      ENDDO

      WRITE(LUN,'(1X,A)') 
     &  'SIG_AIR_INI: NUCLIB interaction lengths (p-air, pi-air, K-air)'
      DO K=1,3
         DO J=1,NSQS
            ASQS = ASQSMIN + DASQS*DBLE(J-1)
            SQS = 10.D0**ASQS
          WRITE(LUN,'(1X,1P,4E12.3)') 
     &           SQS,ALINT(J,1),ALINT(J,2),ALINT(J,3)
         ENDDO
      ENDDO
      END
