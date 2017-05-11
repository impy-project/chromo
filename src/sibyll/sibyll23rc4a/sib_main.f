C=======================================================================
C          SSSSSS   IIIIIII  BBBBB   YY      YY   L        L
C         S            I     B    B    YY  YY     L        L
C          SSSSS       I     BBBBB       YY       L        L
C               S      I     B    B      YY       L        L
C         SSSSSS    IIIIIII  BBBBB       YY       LLLLLLL  LLLLLLL
C=======================================================================
C  Code for SIBYLL:  hadronic interaction Monte Carlo event generator
C=======================================================================
C   Version 2.3rc4   (Sept-07-2015)
C
C     with CHARM production
C
C       By   Eun-Joo Ahn
C            Ralph Engel
C            R.S. Fletcher
C            T.K. Gaisser
C            Paolo Lipari
C            Felix Riehn
C            Todor Stanev
C
C-----------------------------------------------------------------------
C***  Please  have people who want this code contact one of the authors.
C***  Please report any problems.       ****
C
C      For a correct copy contact:
C                sein@fnal.gov
C                ralph.engel@kit.edu
C                gaisser@bartol.udel.edu
C                stanev@bartol.udel.edu
C                paolo.lipari@roma1.infn.it
C                felix.riehn@kit.edu
C
C-----------------------------------------------------------------------

      SUBROUTINE SIBYLL (K_beam, IATARG, Ecm)
C-----------------------------------------------------------------------
C...Main routine for the production of hadronic events,
C.  generates an inelastic hadronic interaction of 
C.  a `projectile particle' of code K_beam with a 
C.  target nucleus of mass number A = IATARG (integer)
C.  IATARG = 0 is an "air" nucleus  (superposition of oxygen and nitrogen)
C.  with c.m. energy for the hadron-nucleon system Ecm (GeV)
C.  
C.  Allowed values of K_beam: 7,8,9,10,11,12,13,14,-13,-14
C.                            pi+-,K+-,KL,KS,p,n,pbar,nbar
C.                    also:   34,35,36,37,38,39           
C.                            Sig+-,Sig0,Xi0-,Lam0
C.                    --> however, apart from the flavor content
C.                        hadronic interactions are assumed to differ 
C.                        only between mesons and baryons.
C.                        different cross sections are only included in
C.                        the calculation of the interaction length.
C.
C.  The output is contained in COMMON /S_PLIST/ that contains:
C.
C.     NP           number of final particles
C.     P(1:NP, 1:5) 4-momenta + masses of the final particles 
C.     LLIST (1:NP) codes of final particles.
C.  the reaction is studied in the c.m. of  hadron-nucleon system
C.
C.  The COMMON block /S_CHIST/ contains information about the
C.  the structure of the  generated event:
C.    NW   = number of wounded nucleons
C.    NJET = total number of hard interactions
C.    NSOF = total number of soft interactions
C.    NNSOF (1:NW) = number of soft pomeron cuts in each interaction
C.    NNJET (1:NW) = number of minijets produced in each interaction 
C.    XJ1 (1:Index) = x1  for each string
C.    XJ2 (1:Index) = x2   "   "     "
C.    PTJET (1:Index) = pT   "   "     "
C.    NNPJET (1:Index) = total number of particles in each string
C.    NNPSTR (1:2*NW) = number of particles in each `beam string'
C.    JDIF(1:NW) = diffraction code
C----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
c     external type declarations
      DOUBLE PRECISION ECM
      INTEGER K_beam, IATARG

c     COMMONs
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_cldif_cmmn.inc'
      INCLUDE 'sib_plist_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      include 'sib_nw_prm.inc'
      INCLUDE 'sib_chist_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'

c     internal type declarations
      DOUBLE PRECISION Esum,PXsum,PYsum,PZsum,xchgRate
      INTEGER LL,IATARGET,IAIR,KBM,L,NW,IREJ,NF,J
      DIMENSION LL(39)
      DATA LL /5*0,7*2,2*1,19*0,6*1/


      if(Ndebug.gt.0)
     &     WRITE(LUN,'(A41,I3,I3,F8.2)')
     &     'SIBYLL: called with (K_beam,IATARG,Ecm):',
     &     K_beam,IATARG,Ecm

 100  CONTINUE
      
      Ncall = Ncall+1

      IATARGET = IATARG
      IAIR = IABS(MIN(IATARG-1,0))
      KBM = K_beam

      CALL INI_EVENT(ECM,KBM,IATARGET,1)

      L = LL(IABS(K_beam))

C...Generate number NW wounded nucleons, and diffraction code.

1000  CALL SIB_START_EV (Ecm, L, IATARGET, IAIR, NWD, JDIF)
      NW = NWD
C...limits on simulation of pure diffraction dissociation
      IF((LDIFF.NE.0).and.(NW.EQ.1)) THEN 
         IF((LDIFF.EQ.-1) .AND. (JDIF(1).NE.0) ) GOTO 1000
         IF((LDIFF.EQ. 1) .AND. ((JDIF(1).NE.0).AND.(JDIF(1).NE.3)))
     +     GOTO 1000
         IF((LDIFF.EQ. 5) .AND. (JDIF(1).EQ.2)) GOTO 1000
         IF((LDIFF.GE. 2) .AND. (LDIFF.LE.4)) THEN
           JDIF(1) = LDIFF-1
         ENDIF
      ENDIF

C...Diffractive/non-diffractive interactions

      IF((NW.EQ.1).and.(JDIF(1).NE.0)) THEN
        CALL SIB_DIFF (KBM, JDIF(1), Ecm, 1, IREJ)
      ELSE
        CALL SIB_NDIFF (KBM, NW, Ecm, 1, IREJ)
      ENDIF

      IF (IREJ.NE.0) THEN
        if(Ndebug.gt.0) WRITE(LUN,'(A38,F8.2,I3,I3,I3)')
     &    'SIBYLL: rejection (Ecm,Ncall,Nw,JDIF):',Ecm,Ncall,NW,JDIF(1)
c        stop
        GOTO 100
      ENDIF

      do J=1,NP
         if (P(J,4).lt.ZERO ) then
            if(Ndebug.gt.0)
     &           WRITE(LUN,*)'negative energy particle!' , P(J,4)
            CALL SIB_LIST(6)
c            stop
            goto 100
         endif
      enddo

C...Check energy-momentum conservation
      
      CALL PFsum(1,NP,Esum,PXsum,PYsum,PZsum,NF)
      IF (ABS(Esum/(HALF*Ecm*DBLE(NW+1)) - ONE) .GT. EPS3)  THEN
         WRITE(LUN,*) ' SIBYLL: energy not conserved (L,call): ',L,Ncall
         WRITE(LUN,*) ' sqs_inp = ', Ecm, ' sqs_out = ', Esum
         call prnt_prtn_stck
         CALL SIB_LIST(LUN)
         WRITE(LUN,*) ' SIBYLL: event rejected'
c         a = -1.D0
c         a = log(a)
c         stop
         goto 100
      ENDIF
      IF (ABS(PZsum+HALF*Ecm*DBLE(NW-1)) .GT. 0.1D0)  THEN
         if(Ndebug.gt.0)THEN
            WRITE(LUN,*) ' SIBYLL: momentum not conserved (L,call): ',
     &           L,Ncall
            WRITE(LUN,*) ' pz_inp = ', 0., ' pz_out = ', pzsum
         ENDIF
         IF(ndebug.gt.0)then
            call prnt_prtn_stck
            CALL SIB_LIST(LUN)
         endif
         WRITE(LUN,*) ' SIBYLL: event rejected'
c         a = -1.D0
c         a = log(a)
c         stop
         goto 100
      ENDIF

c     exchange pions with vector mesons
      IF(IPAR(45).ne.0) then
         xchgRate = PAR(75)
         call force_vectors(xchgRate,1,NP)
      endif

c     exchange pi0 with charged pions for meson projectiles
      IF(IPAR(50).ne.0.and.IABS(KBM).lt.13) then
         xchgrate = PAR(136)
         call remove_pi0(xchgRate,1,NP)
      endif
      
      
C...list final state particles
      if(Ndebug.gt.10) call sib_list(LUN)

      END



      SUBROUTINE SIBNUC (IAB, IATG, ECM)
C-----------------------------------------------------------------------
C.  Routine that generates the interaction of a nucleus of
C.  mass number IAB with a  target nucleus  of mass IAT
C.  (IAT=0 : air).
C.  SQS (GeV) is the  center of mass energy of each
C.  nucleon - nucleon cross section
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      INCLUDE 'sib_plist_cmmn.inc'
      INCLUDE 'sib_mass_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'
      COMMON /S_PLNUC/ PA(5,40000), LLA(40000), NPA
      COMMON /CKFRAG/ KODFRAG
      PARAMETER (IAMAX=56)
      COMMON /CNUCMS/ B, BMAX, NTRY, NA, NB, NI, NAEL, NBEL
     +         ,JJA(IAMAX), JJB(IAMAX), JJINT(IAMAX,IAMAX)
     +         ,JJAEL(IAMAX), JJBEL(IAMAX)            
      COMMON /FRAGMENTS/ PPP(3,60)
      DIMENSION SIGDIF(3)
      DIMENSION IAF(60)
      DATA RPOX /0.3624/

C...Target mass
      IF (IATG .EQ. 0) THEN
         IATARGET = 14 + 2*INT((ONE+RPOX)*S_RNDM(0))
      ELSE
         IATARGET = IATG
      ENDIF
       
C...Single nucleon (proton) case

      IF (IAB .EQ. 1)  THEN
         NPA = 0
         CALL SIBYLL (13,IATARGET, ECM)
         CALL DECSIB
         DO J=1,NP
            LA = IABS(LLIST(J))
            IF (LA .LT. 10000)  THEN
               NPA = NPA + 1
               LLA(NPA) = LLIST(J)
               DO K=1,5
                  PA(K,NPA) = P(J,K)
               ENDDO
            ENDIF
         ENDDO
         RETURN
      ENDIF


C...Nuclei

      CALL SIB_SIGMA_HP(1,ECM,SIGT,SIGEL,SIG0,SIGDIF,SLOPE,RHO)
      CALL INT_NUC (IATARGET, IAB, SIG0, SIGEL) 

C...fragment spectator nucleons
      NBT = NB + NBEL
      IF (KODFRAG .EQ. 1)  THEN
          CALL FRAGM1(IAB,NBT, NF, IAF)
      ELSE IF(KODFRAG .EQ. 2)  THEN
          CALL FRAGM2(IAB,NBT, NF, IAF)
      ELSE 
          CALL FRAGM (IATARGET, IAB, NBT,B, NF, IAF)
      ENDIF
     
C...Spectator fragments
      NPA = 0
      DO J=1,NF
         NPA = NPA+1
         if(NPA.gt.40000) then
           write(6,'(1x,a,2i8)') 
     &       'SIBNUC: no space left in S_PLNUC (NPA,NF)',NPA,NF
           NPA = NPA-1
           return
         endif
         LLA(NPA) = 1000+IAF(J)
         PA(1,NPA) = ZERO
         PA(2,NPA) = ZERO
         PA(3,NPA) = ECM/TWO
         PA(4,NPA) = ECM/TWO
         PA(5,NPA) = DBLE(IAF(J))*HALF*(AM(13)+AM(14))
      ENDDO

C...Elastically scattered fragments
      DO J=1,NBEL
         NPA = NPA+1
         if(NPA.gt.40000) then
           write(6,'(1x,a,2i8)') 
     &       'SIBNUC: no space left in S_PLNUC (NPA,NBEL)',NPA,NBEL
           NPA = NPA-1
           return
         endif
         LLA(NPA) = 1001
         PA(1,NPA) = ZERO
         PA(2,NPA) = ZERO
         PA(3,NPA) = ECM/TWO
         PA(4,NPA) = ECM/TWO
         PA(5,NPA) = HALF*(AM(13)+AM(14))
      ENDDO

C...Superimpose NB  nucleon interactions
      DO JJ=1,NB
          CALL SIBYLL (13,IATARGET, ECM)
          CALL DECSIB
          DO J=1,NP
             LA = IABS(LLIST(J))
             IF (LA .LT. 10000)   THEN
                NPA = NPA + 1
                if(NPA.gt.40000) then
                  write(6,'(1x,a,2i8)') 
     &              'SIBNUC: no space left in S_PLNUC (NPA,NP)',NPA,NP
                  NPA = NPA-1
                  return
                endif
                LLA(NPA) = LLIST(J)
                DO K=1,5
                    PA(K,NPA) = P(J,K)
                ENDDO
            ENDIF
         ENDDO
      ENDDO

      END
