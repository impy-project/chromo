      PROGRAM MAIN
C-----------------------------------------------------------------------
C     
C     Program to study the SIBYLL predictions at a fixed energy 
C     in c.m. frame
C     
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /S_DEBUG/ Ncall, Ndebug
      CHARACTER*6 NAMP
      COMMON /S_CNAM/ NAMP (0:99)
      DOUBLE PRECISION CBR
      INTEGER KDEC,LBARP,IDB
      COMMON /S_CSYDEC/ CBR(223+16+12+8), KDEC(1338+6*(16+12+8)),
     &     LBARP(99), IDB(99)
      COMMON /CCPART/ AA,SQS, XTOT, XXSTAR(-99:99,0:3),
     &    L0,IA,NEVT,NNPARTF(-99:99),NNPARTI(-99:99)
      COMMON /S_MASS1/ AM(99), AM2(99)
      COMMON /S_CLDIF/ LDIFF
      INTEGER NW_max
      PARAMETER (NW_max = 20)
      INTEGER NNSOF,NNJET,JDIF,NWD,NJET,NSOF
      COMMON /S_CHIST/ NNSOF(NW_max),NNJET(NW_max),
     &     JDIF(NW_max),NWD,NJET,NSOF
      SAVE

C...  Initialize the SIBYLL event generator

      CALL SIBYLL_INI

c     intitialize random number generator with default sequence
c     for more options see rndm_dbl.f
      CALL RND_INI

      
C***************************************************************
*     print *,' Particles available in SIBYLL:'
*     print *,' =============================='
*     do I=0,99
*     print *,'   code: ',I,'    name: ',NAMP(I)
*     enddo
C***************************************************************
      
C...  Decide  which events are to be generated

      WRITE (6, 2) ' Possible labels of projectile particle : '
      DO J=6,14
         WRITE (6, 15) J, NAMP(J)
      ENDDO
      WRITE (6, *) ' Give label of projectile particle : '
      READ (5, *) L0
      WRITE (6, *) ' Label of projectile particle : ',L0
      WRITE (6, *) ' Give mass of target nucleus [(1=p),(0:air)] : '
      READ (5, *) IA
      WRITE (6, *) ' Target nucleus : ',IA
      IF (IA .NE. 1) THEN
         WRITE (6, *) ' Give Log10 [E0 eV] : '
         READ (5, *) AA
         ELAB = 10.D0**(AA-9.D0)
         SQS = SQRT(AM2(IABS(L0)) + AM2(13) + 2.D0*AM(13)*ELAB)
      ELSE 
         WRITE (6, *) ' Give  sqrt(s) (GeV)  : '
         READ (5, *) SQS
         ELAB = (SQS**2-AM2(13)-AM2(IABS(L0)))/(2.D0*AM(13))
         WRITE (6,*) 
     &   ' Diffraction codes(-2,-1=nd,0=all,1=nsd,2,3,4=sd1,sd2,dd) :'
         WRITE (6,*) ' Give  diffraction code  : '
         READ (5, *) LDIFF
         WRITE (6,*) ' Diffraction code  : ',LDIFF
      ENDIF


C     debug output
      WRITE (6,*) ' sqrt(s)  = ', SQS
      WRITE (6,*) ' Elab     = ', Elab
      AA = LOG10(ELAB)+9.D0
      WRITE (6,*) ' Log10 [E_lab/eV]  = ', AA

      WRITE (6, *) ' Give number of events to be generated : '
      READ (5, *)  NEVT
C     NEVT= 10000
      WRITE (6, *) ' Number of events to be generated : ',NEVT


C...  Definition of stable particles
      DO J=4,12
         IDB(J) = -abs(IDB(J))
      ENDDO
C     K0s decay
      IDB(12) = abs(IDB(12))
C     Lambda/Anti-lambda stable
      IDB(39) = -abs(IDB(39))
C     Sigmas stable
      do i=34,36
         IDB(i) = -abs(IDB(i))
      enddo
C     Eta stable
*     IDB(23) = -abs(IDB(23))

c     pass event config to sibyll commons, which are read by SIHIST
      call INI_EVENT(SQS, L0, IA, 1)
      
      write(6,*) 'initializing histograms '
      call SIHIST(-1,1.D0)
    
      WRITE (6, *)  ' Starting the event loop.. '

      DO J=1,NEVT

         CALL SIBYLL (L0, IA, SQS)   
         CALL DECSIB
         
         CALL SIHIST(1,1.D0)
        
      ENDDO

      WRITE (6, *) ' End of event loop '

      Xsec = 1.D0
      CALL SIHIST(-2,Xsec)


 1    FORMAT (1X, A, $)
 2    FORMAT (1X, A)
 15   FORMAT (2X, ' L  = ', I2, 2X, A6 )

      END


      SUBROUTINE SIHIST(IMODE,WEIGHT)
C-----------------------------------------------------------------------
C
C   some p-p histograms 
C
C-----------------------------------------------------------------------

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP
      COMMON /S_PLIST1/ LLIST1(NP_max)
      COMMON /S_MASS1/ AM(99), AM2(99)
      INTEGER NW_max
      PARAMETER (NW_max = 20)
      COMMON /S_RUN/ SQSA, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT
      COMMON /S_DEBUG/ Ncall, Ndebug
      COMMON /S_CHP/ ICHP(99), ISTR(99), IBAR(99)
      COMMON /CCPART/ AA,SQS, XTOT, XXSTAR(-99:99,0:3),
     &    L0,IA,NEVT,NNPARTF(-99:99),NNPARTI(-99:99)
      
      INTEGER NNSOF,NNJET,JDIF,NWD,NJET,NSOF
      COMMON /S_CHIST/ NNSOF(NW_max),NNJET(NW_max),
     &     JDIF(NW_max),NWD,NJET,NSOF

      CHARACTER NAMP*6, CODE*18
      COMMON /S_CNAM/ NAMP (0:99)
      SAVE    
      DIMENSION XN(0:3), XL(0:3), XAB(0:3),XPI(0:3),XK(0:3),XG(0:3)

      CHARACTER*72 TITLE
      DIMENSION XLIMB(40),ETRANS(0:5)
      DIMENSION ICHMUL(0:5),IMULC(0:5),IETRA(0:5),
     &          IFEYN1(0:8),IFEYN2(0:8)

      character*4 cintno
      dimension Iimulch(0:20)

C  initialization
***********************************************
      if(imode.eq.-1) then

        IHCALL = 0
     
*        INELA = 0
*        ISDF1 = 0
*        ISDF2 = 0
*        IDDF  = 0
*        IJETP = 0

        XLIM1 = -1.8D0*LOG(SQSA)
        XLIM2 = -XLIM1
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,50,IYYALL)
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,50,IETALL)
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,50,IEALL)

        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,50,IYYCHR)

        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,50,IEAALL)
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,50,IEACHR)
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,50,IEACHR1)
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,50,IEACHR2)
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,50,IEACHR3)
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,50,IETCHR)

C  rapidity distribution of selected particles
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,100,IYYPiP)
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,100,IYYPiM)
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,100,IYYKP)
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,100,IYYKM)

        XLIM1 = 0.D0
        XLIM2 = 5.D0
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,50,IPTALL)
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,50,IPTCHR)

C     energy distribution of selected particles
        XLIM1 = 0.D0
        XLIM2 = SQSA/2.D0
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,100,IEEPiP)
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,100,IEEPiM)
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,100,IEEPi)

C     log. energy distribution of selected particles
        XLIM1 = LOG(0.01D0)
        XLIM2 = LOG(SQSA/2.D0)
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,100,IlEEPiP)
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,100,IlEEPiM)
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,100,IlEEPi)


C  seagull effect
        XLIM1 = -1.D0
        XLIM2 =  1.D0
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,40,ISEAG1)
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,40,ISEAG2)
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,40,ISEAG3)
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,40,ISEAG4)

C  transverse energy distribution
        XLIM1 = 0.D0
        XLIM2 = 50.D0
        DO I=0,5
          CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,50,IETRA(I))
        ENDDO

C  multiplicity distribution
        FACHIS = 2.D0
        XLIM1 = 1.D0
        XLIM2 = 201.D0
        DO 15 I=0,5
          CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,100,IMULC(I))
 15     CONTINUE

C  multiplicity per interaction
        xlim1 = 1.D0
        xlim2 = 301.D0
        do i=0,20
          call newhis(xlim1,xlim2,zero,xlimb,150,Iimulch(i))
        enddo

C  Feynman-x distributions
        XLIM1 = -1.D0
        XLIM2 = 1.D0
        DO I=0,8
          CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,100,IFEYN1(I))
          CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,100,IFEYN2(I))
        ENDDO
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,100,IXFCHR)

        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,200,IFEYN2PiP)
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,200,IFEYN2PiM)
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,200,IFEYN2KP)
        CALL NEWHIS(XLIM1,XLIM2,ZERO,XLIMB,200,IFEYN2KM)


C  elasticity
        CALL NEWHIS(0.D0,1.01D0,ZERO,XLIMB,50,IELAS)

C  particle statistics
*        CALL PARSTA(2,ZERO,ZERO,ZERO,-1)


C  histogram filling
***********************************************
      else if(imode.eq.1) then


        IHCALL = IHCALL+1

*        do i=1,NW
*          IF(JDIF(i).EQ.0) THEN
*            INELA = INELA+1
*          ELSE IF(JDIF(i).EQ.1) THEN
*            ISDF1 = ISDF1+1
*          ELSE IF(JDIF(i).EQ.2) THEN
*            ISDF2 = ISDF2+1
*          ELSE IF(JDIF(i).EQ.3) THEN
*            IDDF  = IDDF+1
*          ENDIF
*        enddo

        IJETP = IJETP+NJET

        DO I=0,5
          ICHMUL(I) = 0
          ETRANS(I) = 0.D0
        ENDDO
 
C...Lorentz boost 
      GAM = SQSA/(AM(13)+AM2(IABS(L0)))

      ARG = 1.D0-1.D0/(GAM**2)
      IF (ARG .GT. 0.D0)  THEN
         BET = SQRT(1.D0-1.D0/(GAM**2))
      ELSE
         BET = 1.D0
      ENDIF
      E0 = (S-AM2(13)-AM2(IABS(L0)))/(2.*AM(13))
*      E0 = (S-2*AM2(13))/(2.D0*AM(13))

      N = 0
      NC = 0
      NPIC = 0
      NPI0 = 0
      NKAON= 0
      NNUC = 0
      NANTI= 0
      YR0 = LOG(SQSA/AM(13))
      XELAM = 0.D0

      DO J=1,NP
         L=LLIST(J)
         LA = IABS(L)
         IF (LLIST1(J) .LE. 0)  THEN
            L1 = MOD(L,10000)
            NNPARTI(L1) = NNPARTI(L1)+1
         ENDIF
         IF (LA .LT. 10000) THEN                 ! stable particles
            NNPARTF(L) = NNPARTF(L) + 1
            N = N+1

            ELAB = GAM*(P(J,4)+BET*P(J,3))
            PLAB = GAM*(P(J,3)+BET*P(J,4))
            IB = IBAR(LA)
            Z = ELAB/E0
*            Z = PLAB/E0
            PT = SQRT(P(J,1)**2 + P(J,2)**2 + 1.D-09)
cjok
cc            YRAP = FRAP(J) + YR0
cc            ETA  = FETA(J) + YR0
            YRAP = FRAP(J)
            ETA  = FETA(J)
cjok
            IF (Z .LE. 0.D0)  THEN
              Z = 0.D0
              AZ = -30.D0
              XXSTAR(L,0) = XXSTAR(L,0) + 1.D0
            ELSE
              AZ = LOG10(Z)
              XXSTAR(L,0) = XXSTAR(L,0) + 1.D0              
              XXSTAR(L,1) = XXSTAR(L,1) + Z                
              XXSTAR(L,2) = XXSTAR(L,2) + Z*Z                
              XXSTAR(L,3) = XXSTAR(L,3) + Z*Z*Z               
            ENDIF
            XTOT = XTOT + Z
*           print *,'SIBYLL: Y,ETA:',YRAP,ETA
C  kinematics
*****************************************
            PT2 = P(J,1)**2+P(J,2)**2
            IF(PT2.LT.1.D-15) THEN
              print *,'WARNING: particle entry skipped (NEV,J,PT2)',
     &          IHCALL,J,PT2
              GOTO 99
            ENDIF
            PABS = SQRT(PT2+P(J,3)**2)
            IF(PABS.LT.1.D-8) THEN
              print *,'WARNING: particle entry skipped (NEV,J,PABS)',
     &          IHCALL,J,PABS
              GOTO 99
            ENDIF
            PT = SQRT(PT2)
*************************************************************
*           IF(ABS(P(J,3)).GE.PABS) THEN
*             print *,'WARNING: particle entry skipped (NEV,J,PABS,PZ)',
*    &          IHCALL,J,PABS,P(J,3)
*             GOTO 99
*           ENDIF
C  pseudorapidity
*           THETA = ACOS(P(J,3)/PABS)
*           ETA = -LOG(TAN(THETA/2.))
*************************************************************
C  high energy approximation for pseudorapidity
            IF(P(J,3).GT.0.D0) THEN
              ETA = LOG((PABS+P(J,3))/PT)
            ELSE
              ETA = LOG(PT/(PABS-P(J,3)))
            ENDIF
C  transverse energy
            IF(ABS(P(J,3)).GE.PABS) THEN
              THETA = 0.D0
            ELSE
              THETA = ACOS(P(J,3)/PABS)
            ENDIF
            ET = P(J,4)*SIN(THETA)
C  rapidity
            PPLUS  = P(J,4)+P(J,3)
            PMINUS = P(J,4)-P(J,3)
            YY = -99999.D0
            IF((PPLUS*PMINUS).GT.1.D-5) YY = 0.5D0*LOG(PPLUS/PMINUS)

*************************************************************
C  check precision of rapidity calculation
*           XMT = SQRT(PT2 + AM(LA)**2)
*           Y2 = SIGN(1.,P(J,3))*LOG((P(J,4)+abs(P(J,3)))/XMT)
*           if(abs(YY-Y2).gt.1.e-2) then
*             write(6,'(a,i3,F8.4,2x,2e12.4)') 
*    &         'WARNING: ID, mass, yy, y2: ',L,AM(LA),yy,y2
*             xm = (P(J,4)-P(J,3))*(P(J,4)+P(J,3))-PT2
*             xm = P(J,4)**2-P(J,3)**2-PT2
*             xm = sqrt(xm)
*             write(6,'(a,2f8.4)') 'mass: (expected, true): ',AM(LA),XM
*           endif
*************************************************************

C  Feynman x
            XF = 2.D0*P(J,3)/SQSA
            XFW = 2.D0*P(J,4)/SQSA
            EE = P(J,4)
****************************************************
*           if(jdif.ne.0) then
*             print *,'J, en, pz, xf ',J,P(J,4),P(J,3),XF 
*           endif
****************************************************

*******************************************

C...all particles

*            CALL PARSTA(L,P(J,4),PT,ETA,2)

C  pseudorapidity distribution of all charged hadrons
            CALL FILHIS(ETA,WEIGHT,IEAALL)
C  rapidity distribution of pi-
            CALL FILHIS(YY,WEIGHT,IYYALL)
C  pt distribution of all charged hadrons
            CALL FILHIS(PT,WEIGHT,IPTALL)
C  transverse energy density
            CALL FILHIS(ETA,ET,IETALL)
C  energy density
            CALL FILHIS(ETA,EE,IEALL)
C  transverse energy distributions
            ABSETA = ABS(ETA)
            ETRANS(0) = ETRANS(0)+ET
            IF(ABSETA.LT.7.D0) THEN
              ETRANS(1) = ETRANS(1)+ET
              IF(ABSETA.LT.5.D0) THEN
                ETRANS(2) = ETRANS(2)+ET
                IF(ABSETA.LT.3.D0) THEN
                  ETRANS(3) = ETRANS(3)+ET
                  IF(ABSETA.LT.2.D0) THEN
                    ETRANS(4) = ETRANS(4)+ET
                    IF(ABSETA.LT.1.5D0) THEN
                      ETRANS(5) = ETRANS(5)+ET
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
C  Feynman-x all
            CALL FILHIS(XF,XFW,IFEYN1(0))
            CALL FILHIS(XF,1.D0,IFEYN2(0))

C...charged

            IF (ICHP(LA)  .NE. 0)  THEN
              NC = NC+1

C  pseudorapidity distribution of all charged hadrons
              CALL FILHIS(ETA,WEIGHT,IEACHR)
C  rapidity distribution of all charged
              CALL FILHIS(YY,WEIGHT,IYYCHR)
C  rapidity distribution of all charged with cuts
              if(abs(xf).lt.0.1D0) CALL FILHIS(ETA,WEIGHT,IEACHR1)
              if(abs(xf).ge.0.1D0) CALL FILHIS(ETA,WEIGHT,IEACHR2)
              if(abs(xf).ge.0.6D0) CALL FILHIS(ETA,WEIGHT,IEACHR3)
C  Feynman-x distribution of all charged
              CALL FILHIS(XF,WEIGHT,IXFCHR)
C  pt distribution of all charged hadrons
              CALL FILHIS(PT,WEIGHT,IPTCHR)
C  transverse energy density
              CALL FILHIS(ETA,ET,IETCHR)

C  rapidity distribution of selected particles
              if (LA.eq.7)  CALL FILHIS(YY,WEIGHT,IYYPiP)
              if (LA.eq.8)  CALL FILHIS(YY,WEIGHT,IYYPiM)
              if (LA.eq.9)  CALL FILHIS(YY,WEIGHT,IYYKP)
              if (LA.eq.10) CALL FILHIS(YY,WEIGHT,IYYKM)

C  energy distribution of selected particles
              if (LA.eq.7)  CALL FILHIS(EE,WEIGHT,IEEPiP)
              if (LA.eq.8)  CALL FILHIS(EE,WEIGHT,IEEPiM)
              if (LA.eq.7.or.LA.eq.8)  CALL FILHIS(EE,WEIGHT,IEEPi)

C  log. energy distribution of selected particles
              if (LA.eq.7)  CALL FILHIS(LOG(EE),WEIGHT,IlEEPiP)
              if (LA.eq.8)  CALL FILHIS(LOG(EE),WEIGHT,IlEEPiM)
              if (LA.eq.7.or.LA.eq.8) CALL FILHIS(LOG(EE),WEIGHT,IlEEPi)

C  seagull
              CALL FILHIS(XF,PT,ISEAG1)
              IF(ICHP(LA).GT.0) CALL FILHIS(XF,PT,ISEAG2)
              IF(ICHP(LA).LT.0) CALL FILHIS(XF,PT,ISEAG3)
              PT2 = PT*PT
              CALL FILHIS(XF,PT2,ISEAG4)

C  charged multiplicity
              ABSETA = ABS(ETA)
              ICHMUL(0) = ICHMUL(0)+1
              IF(ABSETA.LT.7.D0) THEN
                ICHMUL(1) = ICHMUL(1)+1
                IF(ABSETA.LT.5.D0) THEN
                  ICHMUL(2) = ICHMUL(2)+1
                  IF(ABSETA.LT.3.D0) THEN
                    ICHMUL(3) = ICHMUL(3)+1
                    IF(ABSETA.LT.2.D0) THEN
                      ICHMUL(4) = ICHMUL(4)+1
                      IF(ABSETA.LT.1.5D0) THEN
                        ICHMUL(5) = ICHMUL(5)+1
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF

C  Feynman-x charged
              CALL FILHIS(XF,XFW,IFEYN1(1))
              CALL FILHIS(XF,WEIGHT,IFEYN2(1))


            ENDIF

C...pi0's
            IF (LA .EQ. 6)  THEN
               NPI0 = NPI0 + 1

C  Feynman-x pi0
              CALL FILHIS(XF,XFW,IFEYN1(2))
              CALL FILHIS(XF,WEIGHT,IFEYN2(2))

             ENDIF

C...pi+-
            IF (LA .EQ. 7 .OR. LA .EQ. 8)  THEN
               NPIC = NPIC + 1

C  Feynman-x pi+-
              CALL FILHIS(XF,XFW,IFEYN1(3))
              CALL FILHIS(XF,WEIGHT,IFEYN2(3))
              if (LA.eq.7) CALL FILHIS(XF,WEIGHT,IFEYN2PiP)
              if (LA.eq.8) CALL FILHIS(XF,WEIGHT,IFEYN2PiM)


            ENDIF

C...Kaons
            IF (LA .GE. 9 .AND. LA .LE. 12)  THEN
               NKAON = NKAON + 1

C  Feynman-x kaons
              CALL FILHIS(XF,XFW,IFEYN1(4))
              CALL FILHIS(XF,WEIGHT,IFEYN2(4))
              if (LA.eq.9)  CALL FILHIS(XF,WEIGHT,IFEYN2KP)
              if (LA.eq.10) CALL FILHIS(XF,WEIGHT,IFEYN2KM)

            ENDIF

C...Baryons
            IF (IBAR(LA) .EQ. 1 .AND. L .GT. 0)  THEN
               NNUC = NNUC + 1

C  Feynman-x baryons
              CALL FILHIS(XF,XFW,IFEYN1(5))
              CALL FILHIS(XF,1.D0,IFEYN2(5))
              XELA = ELAB/E0
              XELAM = MAX(XELA,XELAM)

            ENDIF

C...Antibaryons
            IF (IBAR(LA) .EQ. 1 .AND. L .LT. 0)  THEN
               NANTI = NANTI + 1

C  Feynman-x antibaryons
              CALL FILHIS(XF,XFW,IFEYN1(6))
              CALL FILHIS(XF,1.D0,IFEYN2(6))

            ENDIF

            if (LLIST(J).eq.13) then
              CALL FILHIS(XF,XFW,IFEYN1(7))
              CALL FILHIS(XF,1.D0,IFEYN2(7))
            endif

            if (LLIST(J).eq.14) then
              CALL FILHIS(XF,XFW,IFEYN1(8))
              CALL FILHIS(XF,1.D0,IFEYN2(8))
            endif
            
            goto 98

 99         continue
            print *,'skipped: ',IHCALL,J,P(J,1),P(J,2),P(J,3),P(J,4)

 98         continue

         ENDIF
      ENDDO


      CALL ADDHIS(IEAALL)
      CALL ADDHIS(IYYALL)
      CALL ADDHIS(IPTALL)
      CALL ADDHIS(IETALL)
      CALL ADDHIS(IEALL)

      CALL ADDHIS(IEACHR)
      CALL ADDHIS(IEACHR1)
      CALL ADDHIS(IEACHR2)
      CALL ADDHIS(IEACHR3)
      CALL ADDHIS(IYYCHR)
      CALL ADDHIS(IXFCHR)
      CALL ADDHIS(IPTCHR)
      CALL ADDHIS(IETCHR)

C  rapidity distribution of selected particles
      CALL ADDHIS(IYYPiP)
      CALL ADDHIS(IYYPiM)
      CALL ADDHIS(IYYKP)
      CALL ADDHIS(IYYKM)

C  energy distribution of selected particles
      CALL ADDHIS(IEEPiP)
      CALL ADDHIS(IEEPiM)
      CALL ADDHIS(IEEPi)

C  log. energy distribution of selected particles
      CALL ADDHIS(IlEEPiP)
      CALL ADDHIS(IlEEPiM)
      CALL ADDHIS(IlEEPi)

C  seagull plots
      CALL ADDHIS(ISEAG1)
      CALL ADDHIS(ISEAG2)
      CALL ADDHIS(ISEAG3)
      CALL ADDHIS(ISEAG4)

C  transverse energy distribution
      DO I=0,5
        CALL FILHIS(REAL(ETRANS(I)),1.D0,IETRA(I))
        CALL ADDHIS(IETRA(I))
      ENDDO

C  charged multiplicity
      DO 12 I=0,5
        CALL FILHIS(REAL(ICHMUL(I)),1.D0,IMULC(I))
        CALL ADDHIS(IMULC(I))
 12   CONTINUE

C  charged multiplicity per interaction
      call filhis(real(ICHMUL(0)),1.D0,Iimulch(0))
      call addhis(Iimulch(0))
      do i=1,20
*        if(Nsof+Njet.eq.i) then
        if(Njet.eq.i) then
          call filhis(real(ICHMUL(0)),1.D0,Iimulch(i))
        endif
        call addhis(Iimulch(i))
      enddo

C  Feynman-x distributions
      DO I=0,8
        CALL ADDHIS(IFEYN1(I))
        CALL ADDHIS(IFEYN2(I))
      ENDDO
      CALL ADDHIS(IFEYN2PiP)
      CALL ADDHIS(IFEYN2PiM)
      CALL ADDHIS(IFEYN2KP)
      CALL ADDHIS(IFEYN2KM)

C  elasticity
      CALL FILHIS(XELAM,WEIGHT,IELAS)
      CALL ADDHIS(IELAS)


C  output of histograms
***********************************************
      else if(imode.eq.-2) then


C  write all to stdout
       lun = 6

       WRITE (LUN,*)   'Sibyll hadronic event generator '
       WRITE (LUN, 41)  L0, NAMP(IABS(L0))
       WRITE (LUN, 42)  IA
       WRITE (LUN, 43)  SQSA
       WRITE (LUN, 44)  AA
       WRITE (LUN,45)   NEVT
 41    FORMAT (' Primary particle = ', I4, 3X, A6)
 42    FORMAT (' Mass of target nucleus = ', I2)
 43    FORMAT (' c.m. energy of the nucleon-nucleon interactions = ', 
     +    F12.2, ' GeV ')
 44    FORMAT( ' Laboratory energy of primary particle = 10**( ',
     +   F5.2,') eV ')   
 45    FORMAT (' Number of events generated = ', I7)
 
C...Composition of final particles
      WRITE(LUN, *) 
      WRITE(LUN, *) ' Composition of final particles '
      CODE = '                  '
      AN = 1.D0/DBLE(NEVT)
      DO J=1,49
         IF (NNPARTF(J) .GT. 0)  THEN
            CODE(1:6) = NAMP(J)
            F = DBLE(NNPARTF(J))
            WRITE(LUN, 430) CODE,NNPARTF(J),F*AN, SQRT(F)*AN
         ENDIF
      ENDDO
      CODE(7:9) = 'bar'
      DO J=-13,-49,-1
         IF (NNPARTF(J) .GT. 0)  THEN
            CODE(1:6) = NAMP(-J)
            F = DBLE(NNPARTF(J))
            WRITE(LUN, 430) CODE,NNPARTF(J),F*AN, SQRT(F)*AN
         ENDIF
      ENDDO 
          
C...Composition of primary particles
      WRITE(LUN, *) 
      WRITE(LUN, *) ' Composition of primary particles '
      CODE = '                  '
      AN = 1.D0/DBLE(NEVT)
      DO J=1,49
         IF (NNPARTI(J) .GT. 0)  THEN
            CODE(1:6) = NAMP(J)
            F = DBLE(NNPARTI(J))
            WRITE(LUN, 430) CODE,NNPARTI(J),F*AN, SQRT(F)*AN
         ENDIF
      ENDDO
      CODE(7:9) = 'bar'
      DO J=-13,-49,-1
         IF (NNPARTI(J) .GT. 0)  THEN
            CODE(1:6) = NAMP(-J)
            F = DBLE(NNPARTI(J))
            WRITE(LUN, 430) CODE,NNPARTI(J),F*AN, SQRT(F)*AN
         ENDIF
      ENDDO
  
C...Study how the energy is diveded among different particles
      WRITE(LUN,*) 
      WRITE(LUN,*) ' Energy distribution : z = E(lab)/E(primary-lab) '
      WRITE(LUN,*) ' Particle           <z**0>      <z**1>   ', 
     +   '   <z**2>      <z**3>'
      YTOT = 0.D0
      DO J=-49,49
        if (XXSTAR(J,0).gt.0) then
          X0 = XXSTAR(J,0)/DBLE(NEVT)
          X1 = XXSTAR(J,1)/DBLE(NEVT)
          X2 = XXSTAR(J,2)/DBLE(NEVT)
          X3 = XXSTAR(J,3)/DBLE(NEVT)
        else
          X1 = 0.D0
          X2 = 0.D0
          X3 = 0.D0
        endif
        IF (J .LT. 0)   THEN
            CODE (1:6) = NAMP(-J)
            CODE(7:9) = 'bar'
        ELSE
            CODE(1:6) = NAMP(J)
            CODE(7:9) = '   '
        ENDIF
        IF (X0 .GT. 1.D-10) WRITE (LUN ,35) CODE, X0, X1, X2, X3
        YTOT  = YTOT + X1
      ENDDO
      WRITE (LUN, * ) 'Check Energy conservation :  ytot = ' , ytot
 
      DO K=0,3
        XN(K) = XXSTAR(13,K)+XXSTAR(14,K)+XXSTAR(-13,K)+XXSTAR(-14,K)
        XN(K) = XN(K)/DBLE(NEVT)
        XL(K) = XXSTAR(13,K)+XXSTAR(14,K)-XXSTAR(-13,K)-XXSTAR(-14,K)
        XL(K) = XL(K)/DBLE(NEVT)
        XAB(K) = +XXSTAR(-13,K)+XXSTAR(-14,K)
        XAB(K) = 2.D0*XAB(K)/DBLE(NEVT)
        XPI(K) = XXSTAR(6,K) +XXSTAR(7,K)+XXSTAR(8,K)
        XPI(K) = XPI(K)/DBLE(NEVT)
        XK(K) = XXSTAR(9,K) +XXSTAR(10,K)+XXSTAR(11,K)+XXSTAR(12,K)
        XK(K) = XK(K)/DBLE(NEVT)
        XG(K) = XXSTAR(1,K)/DBLE(NEVT)
      ENDDO
      WRITE (LUN, *) 
      WRITE(LUN,*) ' Particle           <z**0>       <z**1>   ', 
     +   '     <z**2>       <z**3>'
      WRITE (LUN, 81)  (xn(k) ,k=0,3)
      WRITE (LUN, 82)  (xl(k) ,k=0,3)
      WRITE (LUN, 83)  (xab(k),k=0,3)
      WRITE (LUN, 84)  (xpi(k),k=0,3)
      WRITE (LUN, 85)  (xk(k) ,k=0,3)
      WRITE (LUN, 86)  (xG(k) ,k=0,3)

 81   FORMAT (3X, ' nucleon = ', F12.2, 2X,F12.3, 2X, F12.4, 2X, F12.5) 
 82   FORMAT (3X, ' leading = ', F12.2, 2X,F12.3, 2X, F12.4, 2X, F12.5) 
 83   FORMAT (3X, ' BBbar   = ', F12.2, 2X,F12.3, 2X, F12.4, 2X, F12.5) 
 84   FORMAT (3X, ' pions   = ', F12.2, 2X,F12.3, 2X, F12.4, 2X, F12.5) 
 85   FORMAT (3X, ' kaons   = ', F12.2, 2X,F12.3, 2X, F12.4, 2X, F12.5) 
 86   FORMAT (3X, ' photons = ', F12.2, 2X,F12.3, 2X, F12.4, 2X, F12.5) 
 
35    FORMAT (3X, A9, 3X, 4F12.4)
55    FORMAT (1X, F5.1,2X, 6E11.3)
430   FORMAT (3X,A18, 3X, I7, 3X, F10.4, ' +- ', F8.4)

C  general statistics
*      write(6,*) 
*      write(6,*) '===================================================='
*      write(6,*) 'fraction of non-diff events',DBLE(INELA)/DBLE(NEVT)
*      write(6,*) 'fraction of s-diff-1 events',DBLE(ISDF1)/DBLE(NEVT)
*      write(6,*) 'fraction of s-diff-2 events',DBLE(ISDF2)/DBLE(NEVT)
*      write(6,*) 'fraction of dd-diff events ',DBLE(IDDF)/DBLE(NEVT)
*      write(6,*) 'average number of jet-pairs (single int.,total)',
*     &  DBLE(IJETP)/DBLE(MAX(INELA,1)),DBLE(IJETP)/DBLE(NEVT)
*      write(6,*) '===================================================='

C  particle statistics
        EVENTS = DBLE(NEVT)
*        CALL PARSTA(2,EVENTS,0.D0,0.D0,-2)

        FACMC = 1.D0
        FAC = FACMC*WEIGHT

        TITLE = 'pseudorapidity all hadrons DN/DETA (CMS)'
        CALL OUTHIS(IEAALL,0,TITLE,FACMC,0)
        TITLE = 'rapidity all hadrons DN/DY (CMS)'
        CALL OUTHIS(IYYALL,0,TITLE,FACMC,0)
        TITLE = 'PT dist. all hadrons DSIG/DPT (1/GeV)'
        CALL OUTHIS(IPTALL,1,TITLE,FAC,0)
        TITLE = 'transverse energy all hadrons DET/DETA (GeV)'
        CALL OUTHIS(IETALL,0,TITLE,FACMC,0)
        TITLE = 'energy all hadrons DE/DETA (GeV)'
        CALL OUTHIS(IEALL,0,TITLE,FACMC,0)

        TITLE = 'pseudorapidity charged hadrons DN/DETA (CMS)'
        CALL OUTHIS(IEACHR,0,TITLE,FACMC,0)
        TITLE = 'pseudorapidity charged hadrons DN/DETA (cut 1)'
        CALL OUTHIS(IEACHR1,0,TITLE,FACMC,0)
        TITLE = 'pseudorapidity charged hadrons DN/DETA (cut 2)'
        CALL OUTHIS(IEACHR2,0,TITLE,FACMC,0)
        TITLE = 'pseudorapidity charged hadrons DN/DETA (cut 3)'
        CALL OUTHIS(IEACHR3,0,TITLE,FACMC,0)
        TITLE = 'rapidity charged hadrons DN/DY (CMS)'
        CALL OUTHIS(IYYCHR,0,TITLE,FACMC,0)
        TITLE = 'Feynman-x charged hadrons DN/DXF (CMS)'
        CALL OUTHIS(IXFCHR,0,TITLE,FACMC,0)
        TITLE = 'PT dist. charged hadrons DSIG/DPT (1/GeV)'
        CALL OUTHIS(IPTCHR,1,TITLE,FAC,0)
        TITLE = 'transverse energy charged hadrons DET/DETA (GeV)'
        CALL OUTHIS(IETCHR,0,TITLE,FACMC,0)

C  rapidity distribution of selected particles
        TITLE = 'rapidity distribution DN/DY of pi+'
        CALL OUTHIS(IYYPiP,0,TITLE,FACMC,0)
        TITLE = 'rapidity distribution DN/DY of pi-'
        CALL OUTHIS(IYYPiM,0,TITLE,FACMC,0)
        TITLE = 'rapidity distribution DN/DY of K+'
        CALL OUTHIS(IYYKP,0,TITLE,FACMC,0)
        TITLE = 'rapidity distribution DN/DY of K-'
        CALL OUTHIS(IYYKM,0,TITLE,FACMC,0)

C  energy distribution of selected particles
        TITLE = 'energy distribution DN/DE of pi+'
        CALL OUTHIS(IEEPiP,1,TITLE,FACMC,0)
        TITLE = 'energy distribution DN/DE of pi-'
        CALL OUTHIS(IEEPiM,1,TITLE,FACMC,0)
        TITLE = 'energy distribution DN/DE of pi'
        CALL OUTHIS(IEEPi,1,TITLE,FACMC,0)

C log. energy distribution of selected particles
        TITLE = 'log. energy distribution DN/DLogE of pi+'
        CALL OUTHIS(IlEEPiP,3,TITLE,FACMC,0)
        TITLE = 'energy distribution DN/DLogE of pi-'
        CALL OUTHIS(IlEEPiM,3,TITLE,FACMC,0)
        TITLE = 'energy distribution DN/DLogE of pi'
        CALL OUTHIS(IlEEPi,3,TITLE,FACMC,0)

C  seagull effect
        TITLE = 'SEAGULL PLOT ALL CHARGED'
        CALL OUTHIS(ISEAG1,0,TITLE,FACMC,2)
        TITLE = 'SEAGULL PLOT ALL POSITIVELY CHARGED'
        CALL OUTHIS(ISEAG2,0,TITLE,FACMC,2)
        TITLE = 'SEAGULL PLOT ALL NEGATIVELY CHARGED'
        CALL OUTHIS(ISEAG3,0,TITLE,FACMC,2)
        TITLE = 'PT**2-SEAGULL PLOT ALL CHARGED'
        CALL OUTHIS(ISEAG4,0,TITLE,FACMC,2)

C  transverse energy
        TITLE = 'TRANSVERSE ENERGY 1/N DN/DET (GeV**-1)'
        CALL OUTHIS(IETRA(0),1,TITLE,FACMC,0)
        TITLE = 'TRANSVERSE ENERGY 1/N DN/DET (GeV**-1) ABS(ETA)<7'
        CALL OUTHIS(IETRA(1),1,TITLE,FACMC,0)
        TITLE = 'TRANSVERSE ENERGY 1/N DN/DET (GeV**-1) ABS(ETA)<5'
        CALL OUTHIS(IETRA(2),1,TITLE,FACMC,0)
        TITLE = 'TRANSVERSE ENERGY 1/N DN/DET (GeV**-1) ABS(ETA)<3'
        CALL OUTHIS(IETRA(3),1,TITLE,FACMC,0)
        TITLE = 'TRANSVERSE ENERGY 1/N DN/DET (GeV**-1) ABS(ETA)<2'
        CALL OUTHIS(IETRA(4),1,TITLE,FACMC,0)
        TITLE = 'TRANSVERSE ENERGY 1/N DN/DET (GeV**-1) ABS(ETA)<1.5'
        CALL OUTHIS(IETRA(5),1,TITLE,FACMC,0)

C  charged multiplicity
        TITLE = 'CH. MULTIPLICITY DISTRIBUTION DPN/DN'
        CALL OUTHIS(IMULC(0),1,TITLE,FACHIS,0)
        TITLE = 'CH. MULTIPLICITY DISTRIBUTION DPN/DN ABS(ETA)<7'
        CALL OUTHIS(IMULC(1),1,TITLE,FACHIS,0)
        TITLE = 'CH. MULTIPLICITY DISTRIBUTION DPN/DN ABS(ETA)<5'
        CALL OUTHIS(IMULC(2),1,TITLE,FACHIS,0)
        TITLE = 'CH. MULTIPLICITY DISTRIBUTION DPN/DN ABS(ETA)<3'
        CALL OUTHIS(IMULC(3),1,TITLE,FACHIS,0)
        TITLE = 'CH. MULTIPLICITY DISTRIBUTION DPN/DN ABS(ETA)<2'
        CALL OUTHIS(IMULC(4),1,TITLE,FACHIS,0)
        TITLE = 'CH. MULTIPLICITY DISTRIBUTION DPN/DN ABS(ETA)<1.5'
        CALL OUTHIS(IMULC(5),1,TITLE,FACHIS,0)

        do i=0,20
          write(cintno,'(1x,i3)') i
          title = ' ch. multiplicity for interaction number ' // cintno
          call outhis(Iimulch(i),1,title,fachis,0)
        enddo

C  Feynman-x distribution
        TITLE = 'Feynman-x distribution Xf*dN/dXf (all)'
        CALL OUTHIS(IFEYN1(0),1,TITLE,FACMC,0)
        TITLE = 'Feynman-x distribution dN/dXf (all)'
        CALL OUTHIS(IFEYN2(0),1,TITLE,FACMC,0)

        TITLE = 'Feynman-x distribution Xf*dN/dXf (charged)'
        CALL OUTHIS(IFEYN1(1),1,TITLE,FACMC,0)
        TITLE = 'Feynman-x distribution dN/dXf (charged)'
        CALL OUTHIS(IFEYN2(1),1,TITLE,FACMC,0)

        TITLE = 'Feynman-x distribution Xf*dN/dXf (pi0)'
        CALL OUTHIS(IFEYN1(2),1,TITLE,FACMC,0)
        TITLE = 'Feynman-x distribution dN/dXf (pi0)'
        CALL OUTHIS(IFEYN2(2),1,TITLE,FACMC,0)

        TITLE = 'Feynman-x distribution Xf*dN/dXf (pi+-)'
        CALL OUTHIS(IFEYN1(3),1,TITLE,FACMC,0)
        TITLE = 'Feynman-x distribution dN/dXf (pi+-)'
        CALL OUTHIS(IFEYN2(3),1,TITLE,FACMC,0)

        TITLE = 'Feynman-x distribution dN/dXf (pi+)'
        CALL OUTHIS(IFEYN2PiP,1,TITLE,FACMC,0)
        TITLE = 'Feynman-x distribution dN/dXf (pi-)'
        CALL OUTHIS(IFEYN2PiM,1,TITLE,FACMC,0)

        TITLE = 'Feynman-x distribution Xf*dN/dXf (kaons)'
        CALL OUTHIS(IFEYN1(4),1,TITLE,FACMC,0)
        TITLE = 'Feynman-x distribution dN/dXf (kaons)'
        CALL OUTHIS(IFEYN2(4),1,TITLE,FACMC,0)

        TITLE = 'Feynman-x distribution dN/dXf (K+)'
        CALL OUTHIS(IFEYN2KP,1,TITLE,FACMC,0)
        TITLE = 'Feynman-x distribution dN/dXf (K-)'
        CALL OUTHIS(IFEYN2KM,1,TITLE,FACMC,0)

        TITLE = 'Feynman-x distribution Xf*dN/dXf (baryons)'
        CALL OUTHIS(IFEYN1(5),1,TITLE,FACMC,0)
        TITLE = 'Feynman-x distribution dN/dXf (baryons)'
        CALL OUTHIS(IFEYN2(5),1,TITLE,FACMC,0)

        TITLE = 'Feynman-x distribution Xf*dN/dXf (antibaryons)'
        CALL OUTHIS(IFEYN1(6),1,TITLE,FACMC,0)
        TITLE = 'Feynman-x distribution dN/dXf (antibaryons)'
        CALL OUTHIS(IFEYN2(6),1,TITLE,FACMC,0)

        TITLE = 'Feynman-x distribution Xf*dN/dXf (protons)'
        CALL OUTHIS(IFEYN1(7),1,TITLE,FACMC,0)
        TITLE = 'Feynman-x distribution dN/dXf (protons)'
        CALL OUTHIS(IFEYN2(7),1,TITLE,FACMC,0)

        TITLE = 'Feynman-x distribution Xf*dN/dXf (neutrons)'
        CALL OUTHIS(IFEYN1(8),1,TITLE,FACMC,0)
        TITLE = 'Feynman-x distribution dN/dXf (neutrons)'
        CALL OUTHIS(IFEYN2(8),1,TITLE,FACMC,0)

C  elasticity
        TITLE = 'X_lab of the leading baryon'
        CALL OUTHIS(IELAS,1,TITLE,FACMC,0)


C  unknown imode
      else
        print *,'SIHIST:ERROR:unknown IMODE:',IMODE
        stop
      endif
      
      END

C*****************************************************************

      DOUBLE PRECISION FUNCTION FRAP (I)
C...Compute rapidity for particle I
      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP
      DOUBLE PRECISION AM, AM2
      COMMON /S_MASS1/ AM(99), AM2(99)
      SAVE
      REAL*8 EE,PT,PT2,PZ
      L = IABS(LLIST(I))
      PT2 = P(I,1)**2+P(I,2)**2 + AM2(L) + 1.D-06
      PZ =  P(I,3)
      EE = DSQRT(PT2+PZ**2)
      PT = DSQRT(PT2)
      FRAP = DLOG((EE+PZ)/PT)
      RETURN
      END

      DOUBLE PRECISION FUNCTION FETA (I)
C...Compute pseudorapidity for particle I

      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP
      SAVE
      REAL*8 PP, PT,PT2, PZ
      PT2 = P(I,1)**2+P(I,2)**2+1.D-06
      PZ =  P(I,3)
      PP  = DSQRT(PT2+PZ**2)
      PT  = DSQRT(PT2)
      if(PZ.gt.0.D0) then
        FETA = DLOG((PP+PZ)/PT)      
      else
        FETA = DLOG(PT/(PP-PZ))      
      endif
      RETURN
      END


      SUBROUTINE PARSTA(ID,EE,PT,ETA,IMODE)
C**********************************************************************
C
C     collect particle statistics
C               IMODE  table index
C                      1  collect statistics of resonances
C                      2  collect statistics of final state hadrons
C                         ID     index of particle in /HEPEVS/
C                         EE     particle energy
C                         PT     particle transverse momentum
C                         ETA    pseudorapidity
C                      -1 initialization
C                      -2 output: 
C                         ID = 1 resonances 
C                              2 final state hadrons
C                         EE as the number of events
C
C**********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (DEPS=1.D-7)
C  particle names
      CHARACTER*6 NAMP
      COMMON /S_CNAM/ NAMP (0:99)
C  particle quantum numbers
      COMMON /S_MASS1/ AM(99), AM2(99)
      COMMON /S_CHP/ ICHP(99), ISTR(99), IBAR(99)
C
      SAVE
C  local variables for particle statistics
      PARAMETER(LIMTAB=187)
      DIMENSION ENER(LIMTAB,2),PMUL(LIMTAB,2),PTRA(LIMTAB,2),SUM(2),
     &  CHSUM(2),PTCHR(3,2),PMCHR(3,2),ENCHR(3,2),
     &  CHNALL(2),CHN2P5(2),CHN2P0(2),CHN1P5(2),CHN1P0(2),
     &  PTCALL(2),PTC2P5(2),PTC2P0(2),PTC1P5(2),PTC1P0(2)

      CHARACTER*8 SPCNA(3)
      character*6 name

      IF(IMODE.GT.0) THEN

        if(ID.lt.0) then
C  antibaryons
          IDBAM = abs(ID)
          ICHR  = -ICHP(IDBAM)
          IDBAM = IDBAM+49
        else
C  all others
          IDBAM = ID
          ICHR  = ICHP(IDBAM)
        endif
         
C  collect statistics
        IF(ICHR.NE.0) THEN
          CHSUM(IMODE) = CHSUM(IMODE)+1.D0
          PTCHR(3,IMODE) = PTCHR(3,IMODE)+PT
          PMCHR(3,IMODE) = PMCHR(3,IMODE)+1.D0
          ENCHR(3,IMODE) = ENCHR(3,IMODE)+EE
          IF(ICHR.LT.0) THEN
            K = 1
          ELSE
            K = 2
          ENDIF
          PTCHR(K,IMODE) = PTCHR(K,IMODE)+PT
          PMCHR(K,IMODE) = PMCHR(K,IMODE)+1.D0
          ENCHR(K,IMODE) = ENCHR(K,IMODE)+EE
          CHNALL(IMODE) = CHNALL(IMODE)+1.D0
          PTCALL(IMODE) = PTCALL(IMODE)+PT
          ETAABS = ABS(ETA)
          IF(ETAABS.LT.2.5D0) THEN
            CHN2P5(IMODE) = CHN2P5(IMODE)+1.D0
            PTC2P5(IMODE) = PTC2P5(IMODE)+PT
            IF(ETAABS.LT.2.0D0) THEN
              CHN2P0(IMODE) = CHN2P0(IMODE)+1.D0
              PTC2P0(IMODE) = PTC2P0(IMODE)+PT
              IF(ETAABS.LT.1.5D0) THEN
                CHN1P5(IMODE) = CHN1P5(IMODE)+1.D0
                PTC1P5(IMODE) = PTC1P5(IMODE)+PT
                IF(ETAABS.LT.1.0D0) THEN
                  CHN1P0(IMODE) = CHN1P0(IMODE)+1.D0
                  PTC1P0(IMODE) = PTC1P0(IMODE)+PT
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
        SUM(IMODE) = SUM(IMODE)+1.D0
        IF(IDBAM.LE.LIMTAB) THEN
          PMUL(IDBAM,IMODE) = PMUL(IDBAM,IMODE)+1.D0
          ENER(IDBAM,IMODE) = ENER(IDBAM,IMODE)+EE
          PTRA(IDBAM,IMODE) = PTRA(IDBAM,IMODE)+PT
        ENDIF
      ELSE IF(IMODE.EQ.-1) THEN
C  reset counters for initialization
        II = ID
        CHNALL(II) = 0.
        CHNALL(II) = 0.
        CHN2P5(II) = 0.
        CHN2P0(II) = 0.
        CHN1P5(II) = 0.
        CHN1P0(II) = 0.
        PTC2P5(II) = 0.
        PTC2P0(II) = 0.
        PTC1P5(II) = 0.
        PTC1P0(II) = 0.
        SUM(II)    = 0.
        CHSUM(II) = 0.
        DO 75 I=1,3
          PTCHR(I,II) = 0.
          PMCHR(I,II) = 0.
          ENCHR(I,II) = 0.
 75     CONTINUE
        DO 100 I=1,LIMTAB
          PMUL(I,II) = 0.
          ENER(I,II) = 0.
          PTRA(I,II) = 0.
 100    CONTINUE
        SPCNA(1) = 'neg.chr.'
        SPCNA(2) = 'pos.chr.'
        SPCNA(3) = 'all chr.'
      ELSE IF(IMODE.EQ.-2) THEN
        WRITE(6,'(/1X,A,I10,/1X,A)') 'PARSTA:FINAL OUTPUT:',INT(EE),
     &    '====================='
        IF(EE.LT.0.9D0) RETURN
C  output of collected statistics
        SUM(ID) = SUM(ID)/EE
        CHSUM(ID) = CHSUM(ID)/EE
        DO 150 I=1,3
          PTCHR(I,ID) = PTCHR(I,ID)/EE
          PMCHR(I,ID) = PMCHR(I,ID)/EE
          ENCHR(I,ID) = ENCHR(I,ID)/EE
 150    CONTINUE
        DO 200 I=1,LIMTAB
          PMUL(I,ID) = PMUL(I,ID)/EE
          ENER(I,ID) = ENER(I,ID)/EE
          PTRA(I,ID) = PTRA(I,ID)/EE
 200    CONTINUE
        IF(ID.EQ.1) THEN
          WRITE(6,'(1X,A)') 'STATISTICS FOR PRIMARY RESONANCES'
        ELSE
          WRITE(6,'(1X,A)') 'STATISTICS FOR FINAL STATE PARTICES'
        ENDIF
        WRITE(6,'(5X,A,1p,2E12.4)') 'total,charged multiplicity',
     &    SUM(ID),CHSUM(ID)
        IF(EE*SUM(ID).LE.0.5D0) RETURN
        WRITE(6,'(5X,A,1p,2E12.4)') 'full ETA      N-CH,<PT-CH>:',
     &    CHNALL(ID)/EE,PTCALL(ID)/MAX(CHNALL(ID),1.D0)
        WRITE(6,'(5X,A,1p,2E12.4)') 'abs(ETA)<2.5  N-CH,<PT-CH>:',
     &    CHN2P5(ID)/EE,PTC2P5(ID)/MAX(CHN2P5(ID),1.D0)
        WRITE(6,'(5X,A,1p,2E12.4)') 'abs(ETA)<2.0  N-CH,<PT-CH>:',
     &    CHN2P0(ID)/EE,PTC2P0(ID)/MAX(CHN2P0(ID),1.D0)
        WRITE(6,'(5X,A,1p,2E12.4)') 'abs(ETA)<1.5  N-CH,<PT-CH>:',
     &    CHN1P5(ID)/EE,PTC1P5(ID)/MAX(CHN1P5(ID),1.D0)
        WRITE(6,'(5X,A,1p,2E12.4)') 'abs(ETA)<1.0  N-CH,<PT-CH>:',
     &    CHN1P0(ID)/EE,PTC1P0(ID)/MAX(CHN1P0(ID),1.D0)
        WRITE(6,'(5X,A)')
     &    'PARTICLE    ID    MULT.    REL.MULT.   ENERGY   TRANS.MOM.'
        DO 250 I=1,3
          IF(PMCHR(I,ID).LT.DEPS) GOTO 249
          WRITE(6,'(5X,A8,I6,3X,1p,4E12.4)') SPCNA(I),0,PMCHR(I,ID),
     &      PMCHR(I,ID)/EE,ENCHR(I,ID)/PMCHR(I,ID),
     &      PTCHR(I,ID)/PMCHR(I,ID)
 249      CONTINUE
 250    CONTINUE
        DO 300 I=1,LIMTAB
          IF(PMUL(I,ID).LT.DEPS) GOTO 299
          if(I.le.49) then
            name = NAMP(I)
            IDBAM = I
          else
            name = NAMP(I-49)
            name(4:6) = 'bar'
            IDBAM = 49-I
          endif
          WRITE(6,'(5X,A8,I6,3X,1p,4E12.4)') name,IDBAM,PMUL(I,ID),
     &      PMUL(I,ID)/SUM(ID),ENER(I,ID)/PMUL(I,ID),
     &      PTRA(I,ID)/PMUL(I,ID)
 299      CONTINUE
 300    CONTINUE
      ELSE
        WRITE(6,'(/1X,A,I10)') 'PARSTA:ERROR:INVALID ARGUMENT',ID
        STOP
      ENDIF
      END
