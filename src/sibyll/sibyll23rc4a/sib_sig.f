      SUBROUTINE SIB_SIG(Jint,SIB_SQS,SIB_PTmin,SIB_SIG_tot,
     &                 SIB_SIG_ine,SIB_diff,SIB_diff2,SIB_B_el,SIB_PJET)
C-----------------------------------------------------------------------
C
C...SIBYLL 2.1 cross sections 
C
C   input parameter: SIB_SQS   c.m.s. energy (GeV)
C                    Jint      1 p-p cross sections
C                              2 pi-p cross sections
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE

      include 'sib_int_prm.inc'
      DOUBLE PRECISION SIB_PJET(0:NS_max,0:NH_max)
      DOUBLE PRECISION SIB_SQS,SIB_PTmin,
     &     SIB_SIG_ine,SIB_SIG_tot,SIB_diff(3),SIB_diff2(3,2),SIB_B_el


      COMMON /SIGMAS/SQS,SIGTOT,SIGEL,SIGINE,
     &               SIGSD1(2),SIGSD2(2),SIGDD(2),
     &               SLOPE,SLOPEc,RHO,PROB(0:NS_max,0:NH_max),SIGSUM


      COMMON /PROFILE/XNUS2,XMUS2,XNUSPI2,
     &                XNUH2,XMUH2,XNUHPI2,
     &                ENHPP,ENHPIP,al1,be1,al2,be2

      COMMON /S_CHDCNV/ABR(2,400),ABP(2,400),ABH(2,400),DB,NB

      DIMENSION XI(50)

      DIMENSION SIG_BRN(3)
      DIMENSION SIG_dif_1(2),SIG_dif_2(2),SIG_dd(2)

      DIMENSION IHAR(2)

      INCLUDE 'sib_qcd_cmmn.inc'
      DATA INIT /0/
      INCLUDE 'sib_qcd_p_data.inc'
      INCLUDE 'sib_qcd_pi_data.inc'
      INCLUDE 'sib_xsctn_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'

      IF(INIT.EQ.0) THEN
*        CALL HAR_INI
        CALL FACT_INI
        IHAR(1) = 0
        IHAR(2) = 0
        INIT = 1
      ENDIF

      ECM = SIB_SQS

      IF(JINT.EQ.1) THEN
c     K = 1 , proton
         DO K=1,NPARFIT
            XI(K) = PARS(K,1)
         ENDDO

      ELSE IF(JINT.EQ.2) THEN
c     K = 2 , pion
         DO K=1,NPARFIT
            XI(K) = PARS(K,2)
         ENDDO

      ENDIF

      XNUS2   = XI(12)
      XMUS2   = XI(13)
      XNUSPI2 = XI(14)

      XNUH2   = XI(15)
      XMUH2   = XI(16)
      XNUHPI2 = XI(17)

      CALL HAD_CONV(IABS(JINT))

      PTCUT = XI(10)+XI(21)*dEXP(XI(22)*dSQRT(2.D0*dLOG(ECM)))
      INDX = abs(JINT)
      IHAR(INDX) = IHAR(INDX)+1
      SIGHAR = SIGQCD(IHAR(INDX),INDX)

      S = ECM**2

      BREG =  ABS(XI(18)) + XI(19)*dLOG(S)
      BPOM =  ABS(XI(12)) + XI(13)*dLOG(S)
      IK = ABS(JINT)
      DO JB=1,NB
        B = DB*DBLE(JB-1)
        ABR(IK,JB) = 2.D0/(8.D0*PI*BREG)*dEXP(-B**2/(4.D0*BREG))
        ABP(IK,JB) = 2.D0/(8.D0*PI*BPOM)*dEXP(-B**2/(4.D0*BPOM))
      ENDDO

C  reggeon
      SIGSR = ABS(XI(2))*S**(-ABS(XI(4)))
      SIG_BRN(1) = SIGSR/CMBARN
C  pomeron (soft part)
      SIGSP = ABS(XI(1))*S**ABS(XI(3))
      SIG_BRN(2) = SIGSP/CMBARN
C  pomeron (hard part)
      SIG_BRN(3) = SIGHAR/CMBARN

C  2x2 channel low-mass model and separate high-mass diffraction
        
      al1 = XI(5)
      be1 = XI(6)
      al2 = al1
      be2 = be1
      EnhPP  = XI(9)
      EnhPiP = EnhPP

      CALL SIG_JET_3 (SIG_brn,JINT,SIG_tot,SIG_ela,SIG_ine,SIG_sum,
     &                SIG_dif_1,SIG_dif_2,SIG_dd,B_el,PROB)

      SIGTOT = SIG_tot*CMBARN
      SIGINE = SIG_ine*CMBARN
      SIGSUM = SIG_sum*CMBARN
      SIGELc = SIGTOT-SIGINE
      SIGEL  = SIG_ela*CMBARN
      SIGSD1(1) = SIG_dif_1(1)*CMBARN
      SIGSD1(2) = SIG_dif_1(2)*CMBARN
      SIGSD2(1) = SIG_dif_2(1)*CMBARN
      SIGSD2(2) = SIG_dif_2(2)*CMBARN
      SIGDD(1)  = SIG_dd(1)*CMBARN
      SIGDD(2)  = SIG_dd(2)*CMBARN
      SLOPE  = B_EL
      SLOPEc = SIG_tot**2/(16.D0*Pi*SIG_ela)

      DE = ABS(SIGEL+SIGINE-SIGTOT)/SIGTOT
      IF(DE.GT.0.01) THEN
        print *,'SIBSIG:      Ecm: ',ECM
        print *,'          SIGTOT: ',SIGTOT
        print *,'        SIGEL1/2: ',SIGEL,SIGELc
        print *,'        SLOPE1/2: ',SLOPE,SLOPEc
        print *,'        SIGDIF 1: ',SIGSD1
        print *,'        SIGDIF 2: ',SIGSD2
        print *,'         SIGDDIF: ',SIGDD
        print *,'      SUM-SIGTOT: ',SIGEL+SIGINE-SIGTOT
      ENDIF

C  SIBYLL interface to single precision

      SIB_PTmin   = PTCUT
      SIB_SIG_tot = SIGTOT
      SIB_SIG_ine = SIGINE
      SIB_diff(1) = SIGSD1(1)+SIGSD1(2)
      SIB_diff(2) = SIGSD2(1)+SIGSD2(2)
      SIB_diff(3) = SIGDD(1)+SIGDD(2)
      SIB_B_el    = SLOPE
      DO I=0,NS_max
        DO K=0,NH_max
          SIB_PJET(I,K) = PROB(I,K)
        ENDDO
      ENDDO
c     full diff. cross section 
c     ( ( b.single , t.single , double ) , ( low mass , high mass ) ) 
      SIB_diff2(1,1) = SIGSD1(1)
      SIB_diff2(1,2) = SIGSD1(2)
      SIB_diff2(2,1) = SIGSD2(1)
      SIB_diff2(2,2) = SIGSD2(2)
      SIB_diff2(3,1) = SIGDD(1)
      SIB_diff2(3,2) = SIGDD(2)
      END


      SUBROUTINE SIG_JET_3 (SIG_brn, JINT, SIG_TOT, SIG_ELA, 
     &        SIG_INE, SIG_sum, SIG_DIF1, SIG_DIF2, SIG_DD, B_EL, P_int)
C-----------------------------------------------------------------------
C
C...This subroutine  receives in INPUT:
C.       SIG_brn (GeV-2)  Born graph cross sections
C.       JINT (1 = pp interaction)    (2 pi-p interaction)
C.       neg. value: without calculation of interaction probabilities
C.
C.  and returns as output:
C.       SIG_???  , B_el
C.       and P_int(0:NS_max,0:NH_max)   interaction probabilities
C
C   two x two -channel approximation for diffraction
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE

      DIMENSION SIG_brn(3)
      PARAMETER (NS_max = 20, NH_max = 80)

      COMMON /S_CFACT/ FACT(0:NH_max), CO_BIN(0:NH_max,0:NH_max)
      COMMON /S_CHDCNV/ABR(2,400),ABP(2,400),ABH(2,400),DB,NB

      COMMON /PROFILE/XNUS2,XMUS2,XNUSPI2,
     &                XNUH2,XMUH2,XNUHPI2,
     &                EnhPP,EnhPiP,al1,be1,al2,be2

      DIMENSION SIG_DIF1(2),SIG_DIF2(2),SIG_DD(2),
     &          P_int(0:NS_max,0:NH_max)
      INCLUDE 'sib_utl_cmmn.inc'
 
      DO J=0,NH_max
        DO I=0,NS_max
          P_int(I,J) = ZERO
        ENDDO
      ENDDO

      ga1 = dsqrt(al1*al1+be1*be1)
      ga2 = dsqrt(al2*al2+be2*be2)

      fe_a_1  = (ONE+al1/ga1)/TWO
      fe_a_2  = (ONE-al1/ga1)/TWO
      fd_a_1  = sqrt(ONE-(al1/ga1)**2)/TWO
      fd_a_2  = -fd_a_1

      fe_b_1  = (ONE+al2/ga2)/TWO
      fe_b_2  = (ONE-al2/ga2)/TWO
      fd_b_1  = dsqrt(ONE-(al2/ga2)**2)/TWO
      fd_b_2  = -fd_b_1

      fe_11 = fe_a_1*fe_b_1
      fe_22 = fe_a_2*fe_b_2
      fe_12 = fe_a_1*fe_b_2
      fe_21 = fe_a_2*fe_b_1

      fd_a_11 = fd_a_1*fe_b_1
      fd_a_22 = fd_a_2*fe_b_2
      fd_a_12 = fd_a_1*fe_b_2
      fd_a_21 = fd_a_2*fe_b_1

      fd_b_11 = fe_a_1*fd_b_1
      fd_b_22 = fe_a_2*fd_b_2
      fd_b_12 = fe_a_1*fd_b_2
      fd_b_21 = fe_a_2*fd_b_1

      fdd_11 = fd_a_1*fd_b_1
      fdd_22 = fd_a_2*fd_b_2
      fdd_12 = fd_a_1*fd_b_2
      fdd_21 = fd_a_2*fd_b_1


      sum_abs = ZERO
      sum_tot = ZERO
      sum_ela = ZERO
      sum_sd_a = ZERO
      sum_sd_b = ZERO
      sum_dd  = ZERO
      sum_B   = ZERO

      IK = ABS(JINT)
      if(JINT.GT.0) then
        I0MAX = NS_max
        J0MAX = NH_max
      ELSE
        I0MAX = 1
        J0MAX = 1
      ENDIF
      SIG_REG = SIG_BRN(1)
      SIG_POM = SIG_BRN(2)
      SIG_HAR = SIG_BRN(3)

      DO JB=1,NB

         B = DB*DBLE(JB-1)

         ABREG = ABR(IK,JB)
         ABPOM = ABP(IK,JB)
         ABHAR = ABH(IK,JB)

         chi2_soft = ABREG*SIG_REG+ABPOM*SIG_POM
         chi2_soft_11 = (ONE-al1+ga1)*(ONE-al2+ga2)*chi2_soft
         chi2_soft_22 = (ONE-al1-ga1)*(ONE-al2-ga2)*chi2_soft
         chi2_soft_12 = (ONE-al1+ga1)*(ONE-al2-ga2)*chi2_soft
         chi2_soft_21 = (ONE-al1-ga1)*(ONE-al2+ga2)*chi2_soft

         chi2_hard = ABHAR*SIG_HAR
         chi2_hard_11 = (ONE-al1+ga1)*(ONE-al2+ga2)*chi2_hard
         chi2_hard_22 = (ONE-al1-ga1)*(ONE-al2-ga2)*chi2_hard
         chi2_hard_12 = (ONE-al1+ga1)*(ONE-al2-ga2)*chi2_hard
         chi2_hard_21 = (ONE-al1-ga1)*(ONE-al2+ga2)*chi2_hard
          

         ef_11  = dexp(-HALF*(chi2_soft_11+chi2_hard_11))
         ef_22  = dexp(-HALF*(chi2_soft_22+chi2_hard_22))
         ef_12 = dexp(-HALF*(chi2_soft_12+chi2_hard_12))
         ef_21 = dexp(-HALF*(chi2_soft_21+chi2_hard_21))

         esf_11  = ef_11**2
         esf_22  = ef_22**2
         esf_12  = ef_12**2
         esf_21  = ef_21**2

         F_ine = B*(ONE - fe_11*esf_11 - fe_12*esf_12 
     &                 - fe_21*esf_21 - fe_22*esf_22)
         F_tot = ONE - fe_11*ef_11 - fe_12*ef_12
     &              - fe_21*ef_21 - fe_22*ef_22
         F_ela = B*F_tot**2
         F_tot = B*F_tot

         F_sd_a = B*(fd_a_11*ef_11 + fd_a_12*ef_12
     &              + fd_a_21*ef_21 + fd_a_22*ef_22)**2
         F_sd_b = B*(fd_b_11*ef_11 + fd_b_12*ef_12
     &              + fd_b_21*ef_21 + fd_b_22*ef_22)**2
         F_dd  = B*(fdd_11*ef_11 + fdd_12*ef_12
     &              + fdd_21*ef_21 + fdd_22*ef_22)**2

         sum_abs = sum_abs+F_ine
         sum_tot = sum_tot+F_tot
         sum_ela = sum_ela+F_ela

         sum_sd_a = sum_sd_a+F_sd_a
         sum_sd_b = sum_sd_b+F_sd_b
         sum_dd  = sum_dd +F_dd 

         sum_B   = sum_b+B**2*F_tot

         fac_11 = B*esf_11
         fac_22 = B*esf_22
         fac_12 = B*esf_12
         fac_21 = B*esf_21
         soft_rec_11 = ONE/chi2_soft_11
         soft_rec_22 = ONE/chi2_soft_22
         soft_rec_12 = ONE/chi2_soft_12
         soft_rec_21 = ONE/chi2_soft_21
         chi2_hard_11 = max(chi2_hard_11,EPS10)
         chi2_hard_22 = max(chi2_hard_22,EPS10)
         chi2_hard_12 = max(chi2_hard_12,EPS10)
         chi2_hard_21 = max(chi2_hard_21,EPS10)
         DO I=0,I0MAX
           soft_rec_11 = soft_rec_11*chi2_soft_11
           soft_rec_22 = soft_rec_22*chi2_soft_22
           soft_rec_12 = soft_rec_12*chi2_soft_12
           soft_rec_21 = soft_rec_21*chi2_soft_21
           hard_rec_11 = ONE/chi2_hard_11
           hard_rec_22 = ONE/chi2_hard_22
           hard_rec_12 = ONE/chi2_hard_12
           hard_rec_21 = ONE/chi2_hard_21
           DO J=0,J0MAX
             hard_rec_11 = hard_rec_11*chi2_hard_11
             hard_rec_22 = hard_rec_22*chi2_hard_22
             hard_rec_12 = hard_rec_12*chi2_hard_12
             hard_rec_21 = hard_rec_21*chi2_hard_21
             P_int(I,J) = P_int(I,J) 
     &                + fe_11*soft_rec_11*hard_rec_11*fac_11
     &                + fe_22*soft_rec_22*hard_rec_22*fac_22
     &                + fe_12*soft_rec_12*hard_rec_12*fac_12
     &                + fe_21*soft_rec_21*hard_rec_21*fac_21
           ENDDO
         ENDDO

      ENDDO

      SIG_abs  = SUM_abs*TWO*PI*DB
      SIG_tot  = SUM_tot*FOUR*PI*DB
      SIG_ela  = SUM_ela*TWO*PI*DB
      SIG_dif1(1) = SUM_sd_a*TWO*PI*DB
      SIG_dif2(1) = SUM_sd_b*TWO*PI*DB
      SIG_dd(1)   = SUM_dd*TWO*PI*DB
      SIG_ine  = SIG_abs + SIG_dif1(1) + SIG_dif2(1) + SIG_dd(1)
      B_EL     = sum_B/SUM_tot/TWO

      SA = ZERO
      P_int(0,0) = ZERO
      DO I=0,I0MAX
        DO J=0,J0MAX
          fac = FACT(I)*FACT(J)
          P_int(I,J) = P_int(I,J)/fac
          SA = SA + P_int(I,J)
        ENDDO
      ENDDO

      SIG_hmsd = EnhPP*(P_int(1,0)+P_int(0,1))*TWO*PI*DB
      SIG_hmdd = be1**2*SIG_hmsd + be2**2*SIG_hmsd
     &          + EnhPP**2*P_int(1,1)*TWO*PI*DB

      SIG_dif1(2) = SIG_hmsd
      SIG_dif2(2) = SIG_hmsd
      SIG_dd(2)   = SIG_hmdd

      SIG_sum = SA*TWO*PI*DB

      DO I=0,I0MAX
        DO J=0,J0MAX
          P_int(I,J) = P_int(I,J)/SA
        ENDDO
      ENDDO

      END
C
C
      SUBROUTINE HAD_CONV(JINT)
C-----------------------------------------------------------------------
C
C...Convolution of hadrons profile
C.  [function A(b) of Durand and Pi]
C.  precalculate and put  in COMMON block
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
C
      COMMON /S_CHDCNV/ABR(2,400),ABP(2,400),ABH(2,400),DB,NB

      DOUBLE PRECISION NU2, MU2, NUPI2, NU, MU, NUPI
      COMMON /S_CH0CNV/ NU2, MU2, NUPI2, NU, MU, NUPI

C
      COMMON /PROFILE/XNUS2,XMUS2,XNUSPI2,
     &                XNUH2,XMUH2,XNUHPI2,
     &                ENHPP,ENHPIP,al1,be1,al2,be2

C...integration constants
      BMAX = 50.
      NB  = 400
      DB = BMAX/FLOAT(NB)

C  soft reggeon interactions

      NU2   = XNUS2
      MU2   = XMUS2
      NUPI2 = XNUSPI2

      NU = SQRT(NU2)
      MU = SQRT(ABS(MU2))
      NUPI = SQRT(NUPI2) 

      DO JB=1,NB
         B = DB*FLOAT(JB-1)
         IF(JINT.EQ.1) THEN
           ABR(JINT,JB) = A_PP(B)
         ELSE
           ABR(JINT,JB) = A_PIP(B)
         ENDIF
      ENDDO 

C  soft pomeron interactions

      NU2   = XNUS2
      MU2   = XMUS2
      NUPI2 = XNUSPI2

      NU = SQRT(NU2)
      MU = SQRT(ABS(MU2))
      NUPI = SQRT(NUPI2)

      DO JB=1,NB
         B = DB*FLOAT(JB-1)
         IF(JINT.EQ.1) THEN
           ABP(JINT,JB) = A_PP(B)
         ELSE
           ABP(JINT,JB) = A_PIP(B)
         ENDIF
      ENDDO

C  hard pomeron interactions

      NU2   = XNUH2
      MU2   = XMUH2
      NUPI2 = XNUHPI2

      NU = SQRT(NU2)
      MU = SQRT(ABS(MU2))
      NUPI = SQRT(NUPI2)

      DB = BMAX/FLOAT(NB)
      DO JB=1,NB
         B = DB*FLOAT(JB-1)
         IF(JINT.EQ.1) THEN
           ABH(JINT,JB) = A_PP(B)
         ELSE
           ABH(JINT,JB) = A_PIP(B)
         ENDIF
      ENDDO

      END
C
C
      DOUBLE PRECISION FUNCTION A_pp (b)
C...Convolution of parton distribution for pp interaction
      IMPLICIT DOUBLE PRECISION (A-Z)
      SAVE
C
      COMMON /S_CH0CNV/ NU2, MU2, NUPI2, NU, MU, NUPI
      data pi / 3.1415926/

      ETA = NU2/MU2
 
      IF(ETA.LT.0.D0) THEN
   
        c = nu**5/(96.*pi)
        if (b .gt. 0.0001D0)  then
           A_pp = c*b**3 * bessk (3, b*nu)
        else
           A_pp = nu**2/(12.*pi)
        endif

      ELSE

        X = B*NU
        Y = B*MU
        C = NU2/(12.*PI)/(1.-ETA)**2
        IF(X.GT.0.0001D0) THEN
          A_PP = C*(1./8.*X**3*BESSK(3,X)
     &          -3./2.*ETA/(1.-ETA)*X**2*BESSK(2,X)
     &          +9*ETA**2/(1.-ETA)**2*X*BESSK1(X)
     &          -24*ETA**3/(1.-ETA)**3*(BESSK0(X)-BESSK0(Y))
     &          +3.*ETA**3/(1.-ETA)**2*Y*BESSK1(Y))
        ELSE
          A_PP = C*(1./8.*8.
     &          -3./2.*ETA/(1.-ETA)*2.
     &          +9*ETA**2/(1.-ETA)**2*1.
     &          -24*ETA**3/(1.-ETA)**3*LOG(MU/NU)
     &          +3.*ETA**3/(1.-ETA)**2*1.)
        ENDIF

      ENDIF

      END
C
C
      DOUBLE PRECISION FUNCTION A_pip (b)
C...Convolution of parton distribution for pip interaction
      IMPLICIT DOUBLE PRECISION (A-Z)
      SAVE
C
      COMMON /S_CH0CNV/ NU2, MU2, NUPI2, NU, MU, NUPI
      data pi / 3.1415926/

      eta = nu2/nupi2
      c = nu2/(2.*pi) * 1./(1.-eta)

      if (b .gt. 0.0001D0)  then
         b1 = b*nu
         b2 = b*nupi
         f1 = 0.5*b1 * bessk1(b1)
         f2 = eta/(1.-eta)*(bessk0(b2)- bessk0(b1))
         A_pip = c*(f1+f2)
      else
         A_pip = c*(0.5 + eta/(1.-eta)*log(nu/nupi))
      endif
      return
      end
C
C
C----------------------------------------------------------------------------
C  Bessel functions
C----------------------------------------------------------------------------
C
      FUNCTION BESSK0(X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
C
*     REAL*8 Y,P1,P2,P3,P4,P5,P6,P7,
*    *    Q1,Q2,Q3,Q4,Q5,Q6,Q7
      DATA P1,P2,P3,P4,P5,P6,P7/-0.57721566D0,0.42278420D0,
     *    0.23069756D0,0.3488590D-1,0.262698D-2,0.10750D-3,0.74D-5/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,-0.7832358D-1,
     * 0.2189568D-1,-0.1062446D-1,0.587872D-2,-0.251540D-2,0.53208D-3/
      IF (X.LE.2.0) THEN
        Y=X*X/4.0
        BESSK0=(-LOG(X/2.0)*BESSI0(X))+(P1+Y*(P2+Y*(P3+
     *        Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        Y=(2.0/X)
        BESSK0=(EXP(-X)/SQRT(X))*(Q1+Y*(Q2+Y*(Q3+
     *        Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END
C
C
      FUNCTION BESSK1(X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
C
*     REAL*8 Y,P1,P2,P3,P4,P5,P6,P7,
*    *    Q1,Q2,Q3,Q4,Q5,Q6,Q7
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,0.15443144D0,-0.67278579D0,
     *    -0.18156897D0,-0.1919402D-1,-0.110404D-2,-0.4686D-4/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,0.23498619D0,
     *    -0.3655620D-1,0.1504268D-1,-0.780353D-2,0.325614D-2,
     *    -0.68245D-3/
      IF (X.LE.2.0) THEN
        Y=X*X/4.0
        BESSK1=(LOG(X/2.0)*BESSI1(X))+(1.0/X)*(P1+Y*(P2+
     *      Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        Y=2.0/X
        BESSK1=(EXP(-X)/SQRT(X))*(Q1+Y*(Q2+Y*(Q3+
     *      Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END
C
C
      FUNCTION BESSK(N,X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
C
      IF (N.LT.2) stop 'bad argument N in BESSK'
      TOX=2.0/X
      BKM=BESSK0(X)
      BK=BESSK1(X)
      DO 11 J=1,N-1
        BKP=BKM+J*TOX*BK
        BKM=BK
        BK=BKP
11    CONTINUE
      BESSK=BK
      RETURN
      END
C
C
      FUNCTION BESSI0(X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
C
*     REAL*8 Y,P1,P2,P3,P4,P5,P6,P7,
*    *    Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,3.5156229D0,3.0899424D0,
     *    1.2067492D0,
     *    0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,
     *    0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,
     *    0.2635537D-1,-0.1647633D-1,0.392377D-2/
      IF (ABS(X).LT.3.75) THEN
        Y=(X/3.75)**2
        BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
        AX=ABS(X)
        Y=3.75/AX
        BESSI0=(EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4
     *      +Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END
C
C
      FUNCTION BESSI1(X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE
C
*     REAL*8 Y,P1,P2,P3,P4,P5,P6,P7,
*    *    Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
      DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,
     *    0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1,
     *    -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1,
     *    -0.2895312D-1,0.1787654D-1,-0.420059D-2/
      IF (ABS(X).LT.3.75) THEN
        Y=(X/3.75)**2
        BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        AX=ABS(X)
        Y=3.75/AX
        BESSI1=(EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+
     *      Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END


      SUBROUTINE FACT_INI
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE

      include 'sib_int_prm.inc'
      COMMON /S_CFACT/ FACT(0:NH_max), CO_BIN(0:NH_max,0:NH_max)
      INCLUDE 'sib_utl_cmmn.inc'
      
      FACT(0) = ONE
      DO J=1,NH_max
         FACT(J) = FACT(J-1)*FLOAT(J)
      ENDDO
      DO J=0,NH_max
         DO K=0,J
            CO_BIN(J,K) = FACT(J)/(FACT(K)*FACT(J-K))
         ENDDO
      ENDDO

      END
