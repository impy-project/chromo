C=======================================================================
C          SSSSSS   IIIIIII  BBBBB   YY      YY   L        L
C         S            I     B    B    YY  YY     L        L
C          SSSSS       I     BBBBB       YY       L        L
C               S      I     B    B      YY       L        L
C         SSSSSS    IIIIIII  BBBBB       YY       LLLLLLL  LLLLLLL
C=======================================================================
C  Code for SIBYLL:  hadronic interaction Monte Carlo event generator
C=======================================================================
C     Version 2.3d-star-p05 (Jun-01-2017, modified Feb-13-2025)
C
C     ===>
C          Sibyll 2.3d-star
C           ==
C          Sibyll 2.3d with artificially enhanced muon production
C     ===>
C
C     variants:  ad-hoc increase number of muons achieved by replacing
C     pairs/triples of pions with:      
C
C         * rho-mesons
C         * p-pbar/n-nbar pairs
C         * kaons      
C
C     to make this work ALL particle decays are treated internally!
C
C     To switch between variants set IMOD in COMMON S_STAR
C     values are:
C     0    : no enhancement (but internal decays so this is only approx. equal to sib2.3d)
C     1    : rho-meson enhancement
C     2    : baryon pair enhancement
C     3    : kaon enhancement       
C     4    : mixed enhancement (rho&baryon), default
C     5    : rho-mix (rho component of mixed model)
C     6    : baryon-mix (baryon component of mixed model)  
C
C       By   Eun-Joo Ahn
C            Ralph Engel
C            A. Fedynitch      
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
C                afedynitch@gmail.com 
C                gaisser@bartol.udel.edu
C                paolo.lipari@roma1.infn.it
C                friehn@lip.pt
C                stanev@bartol.udel.edu
C     
C     last changes relative to Sibyll 2.3c:
C     * correct energy of nuclear fragments
C     * fix phi asymmetry
C     * fix off-shell particles in remnant
C     * no pi0 suppression in minijets
C     * added cross section tables for hadron-nitrogen and hadron-oxygen
C       (changed S_CCSIG common)
C     * no remnant in high mass diff. events (pi0-had scattering)
C     * repaired had-nuc. cross section routine for kaon beams
C       routine remains inactive in ordinary calls.
C      
C=======================================================================

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
C.                 also:
C.                  hyperons: 34,35,36,37,38,39
C.                            Sig+-,Sig0,Xi0-,Lam0
C.                        
C.                  charmed:  59,60,71,72,74,75
C.                            D+,D-,D0,D0b,Ds+,Ds-
C.                            87,88,89,99
C.                            Xic+,Xic0,LamC+,OmC0      
C.                  rho0:27 is allowed as well to emulate photons!
C.
C.  The output is contained in COMMON /S_PLIST/ that contains:
C.
C.     NP           number of final particles
C.     P(1:NP, 1:5) 4-momenta + masses of the final particles 
C.     LLIST (1:NP) codes of final particles.
C.  the reaction is studied in the c.m. of  hadron-nucleon system
C.
C.  The COMMON block /S_CHIST/ contains information about 
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
c     external type declarations
      DOUBLE PRECISION ECM
      INTEGER K_beam, IATARG

c     COMMONs
C**anfe adding S_STAR common here to wrap it in f2py
      INTEGER IMOD
      COMMON /S_STAR/ IMOD
      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      INTEGER LDIFF
      COMMON /S_CLDIF/ LDIFF
C     The final particle output is contained in COMMON /S_PLIST/    
C     NP           : number of final particles
C     P(1:NP, 1:5) : 4-momenta + masses of the final particles 
C     LLIST (1:NP) : codes of final particles
      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      INTEGER NW_max
      PARAMETER (NW_max = 20)
C     parameters that represent: NW: max. number of wounded nucleons,
C     NS,NH: max. number of soft and hard interactions
c      PARAMETER (NW_max = 20)
C     The COMMON block /S_CHIST/ contains information about the
C     the structure of the  generated event:
C     NWD   = number of wounded nucleons
C     NJET = total number of hard interactions
C     NSOF = total number of soft interactions
C     NNSOF (1:NW) = number of soft pomeron cuts in each interaction
C     NNJET (1:NW) = number of minijets produced in each interaction 
C     JDIF(1:NW) = diffraction code 
C                  0 : non-diff,
C                  1 : beam-diff
C                  2 : target-diff
C                  3 : double-diff
      INTEGER NNSOF,NNJET,JDIF,NWD,NJET,NSOF
      COMMON /S_CHIST/ NNSOF(NW_max),NNJET(NW_max),
     &     JDIF(NW_max),NWD,NJET,NSOF

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN

c     internal type declarations
      DOUBLE PRECISION Esum,PXsum,PYsum,PZsum,xchgRate
      INTEGER LL,IATARGET,IAIR,KBM,L,NW,IREJ,NF,J
      DIMENSION LL(99)
      SAVE
      DATA LL /5*0,7*2,2*1,12*0,2,6*0,6*1,19*0,2,2,10*0,
     &     2,2,0,2,2,11*0,1,1,1,9*0,1/

      if(Ndebug.gt.0)then
        WRITE(LUN,'(A42,I3,I3,1X,F10.2)')
     &     '  SIBYLL: called with (K_beam,IATARG,Ecm):',
     &     K_beam,IATARG,Ecm
        WRITE(LUN,*)'Event type selection LDIFF: ',LDIFF
      endif

 100  CONTINUE
      
      Ncall = Ncall+1

      IATARGET = IATARG
      IAIR = IABS(MIN(IATARG-1,0))
      KBM = K_beam
      
      CALL INI_EVENT(ECM,KBM,IATARGET,1)

      L = LL(IABS(K_beam))
      IF(L.eq.0) THEN
         WRITE(LUN,*)'SIB_MAIN: unknown beam particle! kbeam=',k_beam
         WRITE(6,*)'SIB_MAIN: unknown beam particle! kbeam=',k_beam
         CALL SIB_REJECT('SIB_MAIN        ')
      endif

C...Generate number NW wounded nucleons, and diffraction code.

1000  continue  
      CALL SIB_START_EV (Ecm, L, IATARGET, IAIR, NWD, JDIF)
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
        if(Ndebug.gt.0) WRITE(LUN,'(A38,F10.2,I3,I3,I3)')
     &   '  SIBYLL: rejection (Ecm,Ncall,Nw,JDIF):',Ecm,Ncall,NW,JDIF(1)
        GOTO 100
      ENDIF

      do J=1,NP
         if (P(J,4).lt.0.D0 ) then
            if(Ndebug.gt.0)then
               WRITE(LUN,*)' negative energy particle!' , P(J,4)
               CALL SIB_LIST(LUN)
            endif
            goto 100
         endif
      enddo

C...Check energy-momentum conservation
      
      CALL PFSUM(1,NP,Esum,PXsum,PYsum,PZsum,NF)
      IF (ABS(Esum/(0.5D0*Ecm*DBLE(NW+1)) - 1.D0) .GT. EPS3)  THEN
         WRITE(LUN,*) ' SIBYLL: energy not conserved (L,call): ',L,Ncall
         WRITE(LUN,*) ' sqs_inp = ', Ecm, ' sqs_out = ', Esum
         CALL PRNT_PRTN_STCK
         CALL SIB_LIST(LUN)
         WRITE(LUN,*) ' SIBYLL: event rejected'
c         a = -1.D0
c         a = log(a)
c         stop
         goto 100
      ENDIF
      IF (ABS(PZsum+0.5D0*Ecm*DBLE(NW-1)) .GT. 0.1D0)  THEN
         if(Ndebug.gt.0)THEN
            WRITE(LUN,*) ' SIBYLL: momentum not conserved (L,call): ',
     &           L,Ncall
            WRITE(LUN,*) ' pz_inp = ', 0., ' pz_out = ', pzsum
         ENDIF
         IF(ndebug.gt.0)then
            CALL PRNT_PRTN_STCK
            CALL SIB_LIST(LUN)
            WRITE(LUN,*) ' SIBYLL: event rejected'
         endif
c         a = -1.D0
c         a = log(a)
c         stop
         goto 100
      ENDIF

c     exchange pions with vector mesons
      IF(IPAR(45).ne.0.and.IPAR(94).eq.0) then
         xchgRate = PAR(75)
         CALL FORCE_VECTORS(xchgRate,0,1.D0,0.D0,1,NP)
      endif

c     exchange pi0 with charged pions for meson projectiles
      IF(IPAR(50).ne.0.and.IABS(KBM).lt.13) then
         xchgrate = PAR(136)
         CALL REMOVE_PI0(xchgRate,1,NP)
      endif      

c     more generic muon enhancements (UHECR 2022)
c     this routine modifies the final state! in particular all resonance decays are enacted!
      IF(IPAR(94).ne.0)THEN
         call MUON_ENHANCEMENTS
      ENDIF
c     remove azimuthal asymmetry
      IF(IPAR(99).eq.1)THEN
         CALL RESAMPLE_FPHI(1,NP)
      ENDIF
      
C...list final state particles
      if(Ndebug.gt.10) CALL SIB_LIST(LUN)

      END


C======================================================================

      SUBROUTINE SIBNUC (IAB, IATG, ECM)

C-----------------------------------------------------------------------
C.  Routine that generates the interaction of a nucleus of
C.  mass number IAB with a  target nucleus  of mass IATG
C.  (IATG=0 : air).
C.  SQS (GeV) is the  center of mass energy of each
C.  nucleon - nucleon cross section
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
C     The final particle output is contained in COMMON /S_PLIST/    
C     NP           : number of final particles
C     P(1:NP, 1:5) : 4-momenta + masses of the final particles 
C     LLIST (1:NP) : codes of final particles
      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      COMMON /S_PLNUC/ PA(5,40000), LLA(40000), NPA
      COMMON /CKFRAG/ KODFRAG
      PARAMETER (IAMAX=56)
      COMMON /CNUCMS/ B, BMAX, NTRY, NA, NB, NI, NAEL, NBEL
     +         ,JJA(IAMAX), JJB(IAMAX), JJINT(IAMAX,IAMAX)
     +         ,JJAEL(IAMAX), JJBEL(IAMAX)            
      COMMON /FRAGMENTS/ PPP(3,60)
      DIMENSION SIGDIF(3)
      DIMENSION IAF(60)
      DOUBLE PRECISION FOX
      SAVE
      DATA FOX /0.21522D0/  !atomic percentage of 'non-nitrogen' in air

C...Target mass
      IF (IATG .EQ. 0) THEN
c  select target IATARGET from air composition
         R = S_RNDM(0)
         IATARGET = 14
         IF (R .LT. FOX)  IATARGET = 16
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
     &       ' SIBNUC: no space left in S_PLNUC (NPA,NF)',NPA,NF
           NPA = NPA-1
           return
         endif
         LLA(NPA) = 1000+IAF(J)
         PA(1,NPA) = 0.D0
         PA(2,NPA) = 0.D0
         PA(3,NPA) = IAF(J)*ECM/2.D0
         PA(4,NPA) = IAF(J)*ECM/2.D0
         PA(5,NPA) = DBLE(IAF(J))*0.5D0*(AM(13)+AM(14))
      ENDDO

C...Elastically scattered fragments
      DO J=1,NBEL
         NPA = NPA+1
         if(NPA.gt.40000) then
           write(6,'(1x,a,2i8)') 
     &       ' SIBNUC: no space left in S_PLNUC (NPA,NBEL)',NPA,NBEL
           NPA = NPA-1
           return
         endif
         LLA(NPA) = 1001
         PA(1,NPA) = 0.D0
         PA(2,NPA) = 0.D0
         PA(3,NPA) = ECM/2.D0
         PA(4,NPA) = ECM/2.D0
         PA(5,NPA) = 0.5D0*(AM(13)+AM(14))
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
     &              ' SIBNUC: no space left in S_PLNUC (NPA,NP)',NPA,NP
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
C=======================================================================

      SUBROUTINE SIBYLL_INI

C-----------------------------------------------------------------------
C...Initialization routine for SYBILL 
C.  
C.  the routine fills the COMMON block /CCSIG/ that contains
C.  important information for the generation of events
C.
C     PARAMETER (NS_max = 20, NH_max = 80)
C     COMMON /S_CCSIG/ SSIG(61,3), PJETC(0:NS_max,0:NH_max,61,2),
C    &    SSIGN(61,3,3),SSIGNSD(61,3,3) ALINT(61,3,3), ASQSMIN, ASQSMAX, DASQS, NSQS
C.
C.  NSQS = number of energy points  (61 is current version)
C.  ASQSMIN = log_10 [sqrt(s) GeV]   minimum value
C.  ASQSMIN = log_10 [sqrt(s) GeV]   maximum value
C.  DASQS   = step  in log_10[sqrt(s)]
C.            DASQS = (ASQSMAX - ASQSMIN)/(NSQS-1)
C.
C.  SSIG(J,1) inelastic cross section for pp interaction
C.            at energy: sqrt(s)(GeV) = 10**[ASQSMIN+DASQS*(J-1)]
C.  SSIG(J,2)  inelastic cross section for pi-p interaction
C.  SSIGN(J,1,1) inelastic cross section for p-Air interaction
C.  SSIGN(J,2,1) inelastic cross section for pi-Air interaction
C.  SSIGN(J,1,2) inelastic cross section for p-Nitrogen interaction
C.  SSIGN(J,2,2) inelastic cross section for pi-Nitrogen interaction
C.
C.  PJETC(n_s,n_j,J,1) Cumulative  probability distribution
C.                 for the production of n_s soft interactions and
C.                 n_j (n_j=0:30) jet pairs at sqrt(s) labeled 
C.                 by J, for p-p interaction
C.  PJETC(n_s,n_j,J,2) Same as above for pi-p interaction
C.  ALINT(J,1)   proton-air  interaction length (g cm-2)
C.  ALINT(J,2)   pi-air  interaction length (g cm-2)
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER IMOD
      COMMON /S_STAR/ IMOD
      SAVE

      WRITE(*,100)
 100  FORMAT(' ','====================================================',
     *     /,' ','|                                                  |',
     *     /,' ','|                 S I B Y L L  2.3e-*              |',
     *     /,' ','|                                                  |')
      if(IMOD.ne.0)then
      WRITE(*,110)
 110  FORMAT(' ','|           ==> MUON ENHANCED VERSION <==          |',
     *     /,' ','|                                                  |',      
     *     /,' ','|  ad-hoc increase of muons by replacing pairs of  |',
     *     /,' ','|  pions with baryons,rho0 or kaons !              |',
     *     /,' ','|                                                  |')
      endif
      WRITE(*,120)
 120  FORMAT(' ','|         HADRONIC INTERACTION MONTE CARLO         |',
     *     /,' ','|                        BY                        |',
     *     /,' ','|            Eun-Joo AHN, Felix RIEHN              |',
     *     /,' ','|      R. ENGEL, A. FEDYNITCH, R.S. FLETCHER,      |',
     *     /,' ','|       T.K. GAISSER, P. LIPARI, T. STANEV         |',
     *     /,' ','|                                                  |',
     *     /,' ','| Publication to be cited when using this program: |',
     *     /,' ','| Eun-Joo AHN et al., Phys.Rev. D80 (2009) 094003  |',
     *     /,' ','| F. RIEHN et al., Phys.Rev. D102 (2020) 063002    |',
     *     /,' ','| last modifications: F. Riehn (10/02/2025)        |',
     *     /,' ','====================================================',
     *     /)

      CALL PAR_INI
      CALL SIBYLL_STAR_INI
      CALL DIFF_INI
      CALL JET_INI
      CALL PDF_INI
      CALL BLOCK_INI
      CALL NUC_GEOM_INI
      CALL SIG_AIR_INI
      CALL DEC_INI
c...  charm frag. normalisation
      CALL ZNORMAL
      
      END

C=======================================================================

      SUBROUTINE NO_CHARM
      IMPLICIT NONE
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
c     turn off charm production
c     global charm rate
      PAR(24) = 0.D0
c     minijet string charm rate
      PAR(156) = 0.D0
c     remnant string charm rate
      PAR(107) = 0.D0
c     soft sea charm rate
      PAR(97) = 0.D0
c     valence string charm rate
      PAR(25) = 0.D0
c     minijet charm rate
      PAR(27) = 0.D0
      END

C=======================================================================      
      SUBROUTINE SIBYLL_STAR_INI
      IMPLICIT NONE
      INTEGER IMOD
      COMMON /S_STAR/ IMOD
      IF(IMOD.eq.0)THEN
         CALL STD_INI
      ELSEIF(IMOD.EQ.1)THEN
         CALL VECTOR_INI
      ELSEIF(IMOD.eq.2)THEN
         CALL BARYON_INI
      ELSEIF(IMOD.eq.3)THEN
         CALL STRANGE_INI
      ELSEIF(IMOD.eq.4)THEN
         CALL VECTOR_BARYON_INI
      ELSEIF(IMOD.EQ.5)THEN
         CALL VECTORMIX_INI
      ELSEIF(IMOD.eq.6)THEN
         CALL BARYONMIX_INI         
      ELSE
         WRITE(*,*) 'SIBYLL STAR. Wrong initialization!'
         STOP
      ENDIF

      END
      
C=======================================================================
      SUBROUTINE STD_INI

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)

c     enhancement model (1: vector, 2: strangeness, 3: baryons, -1: no enhancement but do decays internally)
      IPAR(94) = -1

      WRITE(*,*),'===================================================='
      WRITE(*,*),'= NO ENHANCEMENT, BUT DECAYS ARE HANDLED IN SIBYLL ='
      WRITE(*,*),'===================================================='
      END
      
C=======================================================================

      SUBROUTINE BARYON_INI

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
c     muon extension via baryon enhancement
c     
c     enhancement model (1: vector, 2: strangeness, 3: baryons)
      IPAR(94) = 30
c     energy dependence (0:none, 1: logarithmic)
      IPAR(96) = 1
c     which projectiles (0:all, 1: mesons only)
      IPAR(97) = 0
c     LE 
c     exchange rate
      PAR(75) = 0.5
c     power of xf weight (0: no weighting (all weight=1), 1: weight = xf)
      PAR(159) = 0.7            
c     energy threshold (GeV center-of-mass)
      PAR(160) = 5.
c     HE
c     exchange rate
      PAR(161) = 0.25
c     power of xf weight (0: no weighting (all weight=1), 1: weight = xf)
      PAR(163) = 0.0            
c     energy threshold (GeV center-of-mass)
      PAR(162) = 13000.      

      WRITE(*,*),'===================================================='
      WRITE(*,*),'=  BARYON ENHANCEMENT 4  LE & HE                   ='
      WRITE(*,*),'===================================================='
      END
      
C=======================================================================

      SUBROUTINE BARYONMIX_INI

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
c     muon extension via baryon enhancement
c     
c     enhancement model (1: vector, 2: strangeness, 3: baryons)
      IPAR(94) = 3
c     energy dependence (0:none, 1: logarithmic)
      IPAR(96) = 1
c     which projectiles (0:all, 1: mesons only)
      IPAR(97) = 0
c     LE 
c     exchange rate
      PAR(75) = 0.5
c     power of xf weight (0: no weighting (all weight=1), 1: weight = xf)
      PAR(159) = 0.7            
c     energy threshold (GeV center-of-mass)
      PAR(160) = 5.

      WRITE(*,*),'===================================================='
      WRITE(*,*),'=  BARYON ENHANCEMENT 1                            ='
      WRITE(*,*),'===================================================='
      END
      
C=======================================================================

      SUBROUTINE VECTOR_INI

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
c     muon extension via vector (rho0) enhancement
c     
c     enhancement model (1: vector, 2: strangeness, 3: baryons)
      IPAR(94) = 1
c     exchange rate
      PAR(75) = 0.9
c     power of xf weight (0: no weighting (all weight=1), 1: weight = xf)
      PAR(159) = 0.3            
c     energy threshold (GeV center-of-mass)
      PAR(160) = 5.      
c     energy dependence (0:none, 1: logarithmic)
      IPAR(96) = 1
c     which projectiles (0:all, 1: mesons only)
      IPAR(97) = 1

      WRITE(*,*),'===================================================='
      WRITE(*,*),'=  RHO0 ENHANCEMENT 3                              ='
      WRITE(*,*),'===================================================='

      END
      
C=======================================================================

      SUBROUTINE VECTORMIX_INI

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
c     muon extension via vector (rho0) enhancement
c     
c     enhancement model (1: vector, 2: strangeness, 3: baryons)
      IPAR(94) = 1
c     exchange rate
      PAR(75) = 0.95
c     power of xf weight (0: no weighting (all weight=1), 1: weight = xf)
      PAR(159) = 0.4
c     energy threshold (GeV center-of-mass)
      PAR(160) = 5.      
c     energy dependence (0:none, 1: logarithmic)
      IPAR(96) = 1
c     which projectiles (0:all, 1: mesons only)
      IPAR(97) = 1

      WRITE(*,*),'===================================================='
      WRITE(*,*),'=  RHO0 ENHANCEMENT 2                              ='
      WRITE(*,*),'===================================================='

      END
      
C=======================================================================
      
      SUBROUTINE VECTOR_BARYON_INI

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
c     muon extension via vector (rho0) enhancement
c     
c     enhancement model (1: vector, 2: strangeness, 3: baryons)
      IPAR(94) = 4
c     exchange rate
      PAR(75) = 0.95
c     power of xf weight (0: no weighting (all weight=1), 1: weight = xf)
      PAR(159) = 0.4            
c     energy threshold (GeV center-of-mass)
      PAR(160) = 5.      
c     energy dependence (0:none, 1: logarithmic)
      IPAR(96) = 1
c     which projectiles (0:all, 1: mesons only)
      IPAR(97) = 1

c     baryon part
c     exchange rate
      PAR(161) = 0.5
c     power of xf weight (0: no weighting (all weight=1), 1: weight = xf)
      PAR(163) = 0.7            
c     energy threshold (GeV center-of-mass)
      PAR(162) = 5.      
c     which projectiles (0:all, 1: mesons only)
      IPAR(98) = 0      

      WRITE(*,*),'===================================================='
      WRITE(*,*),'=  RHO0 & BARYON ENHANCEMENT 2                     ='
      WRITE(*,*),'===================================================='

      END
      
C=======================================================================

      SUBROUTINE STRANGE_INI

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
c     muon extension via strangeness enhancement
c     
c     enhancement model (1: vector, 2: strangeness, 3: baryons)
      IPAR(94) = 20
c     energy dependence (0:none, 1: logarithmic)
      IPAR(96) = 1
c     which projectiles (0:all, 1: mesons only)
      IPAR(97) = 0
C     LE
c     exchange rate
      PAR(75) = 0.5
c     power of xf weight (0: no weighting (all weight=1), 1: weight = xf)
      PAR(159) = 0.8            
c     energy threshold (GeV center-of-mass)
      PAR(160) = 5.
C     HE
c     exchange rate
      PAR(161) = 0.3
c     power of xf weight (0: no weighting (all weight=1), 1: weight = xf)
      PAR(163) = 0.            
c     energy threshold (GeV center-of-mass)
      PAR(162) = 13000.      

      WRITE(*,*),'===================================================='
      WRITE(*,*),'=  STRANGENESS ENHANCEMENT 3  LE & HE              ='
      WRITE(*,*),'===================================================='
      END
      
C=======================================================================
      
      SUBROUTINE PAR_INI

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      DOUBLE PRECISION FAin, FB0in
      COMMON /S_CZDIS/ FAin, FB0in

      DOUBLE PRECISION FAs1, fAs2
      COMMON /S_CZDISs/ FAs1, fAs2
      DOUBLE PRECISION ZDMAX, EPSI
      COMMON /S_CZDISc/ ZDMAX, EPSI

      DOUBLE PRECISION CLEAD, FLEAD
      COMMON /S_CZLEAD/ CLEAD, FLEAD
      DOUBLE PRECISION CCHIK
      COMMON /S_CPSPL/ CCHIK(4,99)

      PARAMETER ( NPARFIT = 22 )
      DOUBLE PRECISION PARS
      COMMON /XSCTN_FIT/ PARS( 50 , 2 )
      DOUBLE PRECISION STR_mass_val, STR_mass_val_hyp, STR_mass_sea
      COMMON /S_CUTOFF/ STR_mass_val, STR_mass_val_hyp, STR_mass_sea
      COMMON /S_STAR/ IMOD
      SAVE
      DATA (PARS(K,1),K=    1,NPARFIT) /
     &3.9223D+01,4.2055D+01,5.0913D-02,-4.0000D-01,2.0000D-01,
     &5.0000D-01,0.0000D+00,6.0000D-01,9.0000D-02,1.0000D+00,
     &2.0000D+00,3.2327D+00,2.5000D-01,5.4000D-01,1.0000D+00,
     &-8.8000D-01,5.4000D-01,5.0000D-01,9.0000D-01,5.4000D-01,
     &6.5000D-02,9.0000D-01/
      DATA (PARS(K,2),K=    1,NPARFIT) /
     &2.0590D+01,9.6579D+01,5.6069D-02,-7.6393D-01,2.0000D-01,
     &5.0000D-01,0.0000D+00,6.0000D-01,9.0000D-02,1.0000D+00,
     &2.0000D+00,2.9191D+00,2.5000D-01,5.4000D-01,1.0000D+00,
     &-8.8000D-01,5.4000D-01,5.4895D-01,9.0000D-01,5.4000D-01,
     &6.5000D-02,9.0000D-01/
      DATA IMOD /4/
c     
c     adjusted central particle production
c     23rc5.4frgB1 aka retune5 aka Sibyll 2.3.5
c     including muon boost. parameters: ipar45, ipar94, par75, par159, par160
      PAR(1) = 4.0000D-02
      PAR(2) = 2.5000D-01
      PAR(3) = 5.0000D-01
      PAR(4) = 1.4000D-01
      PAR(5) = 3.0000D-01
      PAR(6) = 3.0000D-01
      PAR(7) = 1.5000D-01
      PAR(8) = 1.3903D-02
      PAR(9) = 7.0000D+00
      PAR(10) = 1.0000D+00
      PAR(11) = 6.5000D-02
      PAR(12) = 9.0000D-01
      PAR(13) = 1.0000D-01
      PAR(14) = 6.0000D-02
      PAR(15) = 1.3000D-01
      PAR(16) = 4.0000D-02
      PAR(17) = 4.0000D-02
      PAR(18) = 5.0000D-01
      PAR(19) = 8.0000D-01
      PAR(20) = 8.0000D-01
      PAR(21) = 6.0000D-01
      PAR(22) = 4.0000D+00
      PAR(23) = 7.0000D-01
      PAR(24) = 4.0000D-03
      PAR(25) = 4.0000D-03
      PAR(26) = 2.0000D+01
      PAR(27) = 2.0000D-02
      PAR(28) = 2.0000D+01
      PAR(29) = 0.0000D+00
      PAR(30) = 2.0000D+00
      PAR(31) = 3.3000D-01
      PAR(32) = 0.0000D+00
      PAR(33) = 1.0000D-01
      PAR(34) = 0.0000D+00
      PAR(35) = 0.0000D+00
      PAR(36) = 7.0000D-01
      PAR(37) = 0.0000D+00
      PAR(38) = 5.0000D-01
      PAR(39) = 8.0000D-01
      PAR(40) = 0.0000D+00
      PAR(41) = 1.0000D+00
      PAR(42) = 0.0000D+00
      PAR(43) = 2.3564D-01
      PAR(44) = 9.9000D-01
      PAR(45) = 1.0000D+00
      PAR(46) = 1.8000D-01
      PAR(47) = 2.8000D-01
      PAR(48) = 2.7000D-01
      PAR(49) = 1.0000D-01
      PAR(50) = 6.0000D-01
      PAR(51) = 6.0000D-03
      PAR(52) = 6.0000D-03
      PAR(53) = 6.0000D+00
      PAR(54) = 2.0000D-01
      PAR(55) = 0.0000D+00
      PAR(56) = 0.0000D+00
      PAR(57) = 0.0000D+00
      PAR(58) = 0.0000D+00
      PAR(59) = 6.8345D-01
      PAR(60) = 8.0000D-01
      PAR(61) = 6.6000D-01
      PAR(62) = 0.0000D+00
      PAR(63) = 1.0000D+00
      PAR(64) = 2.5000D-01
      PAR(65) = 3.0000D-01
      PAR(66) = 3.0000D-01
      PAR(67) = 6.0000D-01
      PAR(68) = 6.0000D-03
      PAR(69) = 5.0000D-02
      PAR(70) = 7.0000D-03
      PAR(71) = 1.0000D+00
      PAR(72) = 3.8000D-01
      PAR(73) = 5.0000D-01
      PAR(74) = 6.0000D-01
      PAR(75) = 0.0000D+00      ! rate for ad-hoc exchange of pion w. rho, strange particles or baryons
      PAR(76) = 3.5298D-01
      PAR(77) = 7.0000D-01
      PAR(78) = 2.0000D+00
      PAR(79) = 1.0000D+01
      PAR(80) = 5.0816D-01
      PAR(81) = 1.0000D+04
      PAR(82) = 1.0000D-01
      PAR(83) = 0.0000D+00
      PAR(84) = 6.0000D+00
      PAR(85) = 1.0000D+00
      PAR(86) = 1.0000D+00
      PAR(87) = 3.0000D-01
      PAR(88) = 8.0000D-01
      PAR(89) = 6.0000D-01
      PAR(90) = 1.1000D+01
      PAR(91) = -7.2000D+00
      PAR(92) = 3.5000D+00
      PAR(93) = 1.0000D+00
      PAR(94) = 4.0000D+00
      PAR(95) = 0.0000D+00
      PAR(96) = 1.0000D+00
      PAR(97) = 2.0000D-03
      PAR(98) = 1.5000D+00
      PAR(99) = 5.0000D-01
      PAR(100) = 2.0000D+00
      PAR(101) = 1.0000D+00
      PAR(102) = 0.0000D+00
      PAR(103) = 2.0000D+00
      PAR(104) = 4.0000D-01
      PAR(105) = 1.0000D-01
      PAR(106) = 0.0000D+00
      PAR(107) = 0.0000D+00
      PAR(108) = 0.0000D+00
      PAR(109) = 2.0000D+01
      PAR(110) = 1.5000D+00
      PAR(111) = 0.0000D+00
      PAR(112) = 7.0000D-01
      PAR(113) = 8.0000D-01
      PAR(114) = 2.0000D+00
      PAR(115) = 0.0000D+00
      PAR(116) = 1.0000D+00
      PAR(117) = 0.0000D+00
      PAR(118) = 5.0000D-03
      PAR(119) = 0.0000D+00
      PAR(120) = 1.0000D+00
      PAR(121) = 3.0000D-01
      PAR(122) = 0.0000D+00
      PAR(123) = 3.0000D-01
      PAR(124) = 1.0000D+00
      PAR(125) = 1.0000D+00
      PAR(126) = 1.0000D+00
      PAR(127) = 6.0000D+00
      PAR(128) = 1.0000D+00
      PAR(129) = 8.0000D-02
      PAR(130) = 1.2000D+01
      PAR(131) = 5.0000D-01
      PAR(132) = 5.0000D-01
      PAR(133) = 1.0000D+01
      PAR(134) = -5.0000D+00
      PAR(135) = 6.0000D+00
      PAR(136) = 0.0000D+00
      PAR(137) = 1.2000D+00
      PAR(138) = 0.0000D+00
      PAR(139) = 5.0000D-01
      PAR(140) = 4.5000D-01
      PAR(141) = 1.5000D+00
      PAR(142) = 0.0000D+00
      PAR(143) = 5.0000D-01
      PAR(144) = 9.5000D-01
      PAR(145) = 8.5000D-01
      PAR(146) = 0.0000D+00
      PAR(147) = 3.0000D-01
      PAR(148) = 5.0000D-01
      PAR(149) = 3.0000D-01
      PAR(150) = 4.0000D-03
      PAR(151) = 2.0000D+00
      PAR(152) = 4.0000D+00
      PAR(153) = 1.0000D+01
      PAR(154) = 3.0000D-01
      PAR(155) = 0.0000D+00
      PAR(156) = 5.0000D-01
      PAR(157) = 8.0000D-01
      PAR(158) = 0.0000D+00
      PAR(159) = 0.0000D-00     ! exponent of xf weighting for enhancement
      PAR(160) = 0.0000D+00     ! start energy for enhancements in CoM and GeV
      PAR(161) = 0.0000D+00     ! 2nd HE rate
      PAR(162) = 0.0000D+00     ! start energy for HE enhancements in CoM and GeV
      PAR(163) = 0.0000D+00     ! exponent of xf weighting for enhancement at HE
      PAR(164) = 0.0000D+00
      PAR(165) = 0.0000D+00
      PAR(166) = 0.0000D+00
      PAR(167) = 0.0000D+00
      PAR(168) = 0.0000D+00
      PAR(169) = 0.0000D+00
      PAR(170) = 0.0000D+00
      PAR(171) = 0.0000D+00
      PAR(172) = 0.0000D+00
      PAR(173) = 0.0000D+00
      PAR(174) = 0.0000D+00
      PAR(175) = 0.0000D+00
      PAR(176) = 0.0000D+00
      PAR(177) = 0.0000D+00
      PAR(178) = 0.0000D+00
      PAR(179) = 0.0000D+00
      PAR(180) = 0.0000D+00
      PAR(181) = 0.0000D+00
      PAR(182) = 0.0000D+00
      PAR(183) = 0.0000D+00
      PAR(184) = 0.0000D+00
      PAR(185) = 0.0000D+00
      PAR(186) = 0.0000D+00
      PAR(187) = 0.0000D+00
      PAR(188) = 0.0000D+00
      PAR(189) = 0.0000D+00
      PAR(190) = 0.0000D+00
      PAR(191) = 0.0000D+00
      PAR(192) = 0.0000D+00
      PAR(193) = 0.0000D+00
      PAR(194) = 0.0000D+00
      PAR(195) = 0.0000D+00
      PAR(196) = 0.0000D+00
      PAR(197) = 0.0000D+00
      PAR(198) = 0.0000D+00
      PAR(199) = 0.0000D+00
      PAR(200) = 0.0000D+00
      IPAR(1) = 1
      IPAR(2) = 0
      IPAR(3) = 8
      IPAR(4) = 0
      IPAR(5) = 1
      IPAR(6) = 0
      IPAR(7) = 0
      IPAR(8) = 1
      IPAR(9) = 1
      IPAR(10) = 1
      IPAR(11) = 0
      IPAR(12) = 3
      IPAR(13) = 0
      IPAR(14) = -2
      IPAR(15) = 9
      IPAR(16) = 8
      IPAR(17) = 1
      IPAR(18) = 4
      IPAR(19) = 1
      IPAR(20) = 0
      IPAR(21) = 0
      IPAR(22) = 0
      IPAR(23) = 0
      IPAR(24) = 0
      IPAR(25) = 1
      IPAR(26) = 0
      IPAR(27) = 0
      IPAR(28) = 4
      IPAR(29) = 1
      IPAR(30) = 0
      IPAR(31) = 1
      IPAR(32) = 0
      IPAR(33) = 0
      IPAR(34) = 0
      IPAR(35) = 0
      IPAR(36) = 1
      IPAR(37) = 0
      IPAR(38) = 1
      IPAR(39) = 0
      IPAR(40) = 0
      IPAR(41) = 0
      IPAR(42) = 3
      IPAR(43) = 1
      IPAR(44) = 0
      IPAR(45) = 0  ! model for exchange pi0 with rho0, rate is par(75), switch between before or after decay is ipar(94)
      IPAR(46) = 2
      IPAR(47) = 6
      IPAR(48) = 1
      IPAR(49) = 4
      IPAR(50) = 0
      IPAR(51) = 2
      IPAR(52) = 0
      IPAR(53) = 1
      IPAR(54) = 0
      IPAR(55) = 0
      IPAR(56) = 0
      IPAR(57) = 1
      IPAR(58) = 3
      IPAR(59) = 1
      IPAR(60) = 0
      IPAR(61) = 100
      IPAR(62) = 1
      IPAR(63) = 0
      IPAR(64) = 0
      IPAR(65) = 1
      IPAR(66) = 3
      IPAR(67) = 0
      IPAR(68) = 0
      IPAR(69) = 1
      IPAR(70) = 1
      IPAR(71) = 0
      IPAR(72) = 0
      IPAR(73) = 0
      IPAR(74) = 1
      IPAR(75) = 0
      IPAR(76) = 0
      IPAR(77) = 0
      IPAR(78) = 2
      IPAR(79) = 1
      IPAR(80) = 1
      IPAR(81) = 5
      IPAR(82) = 2
      IPAR(83) = 0
      IPAR(84) = 2
      IPAR(85) = 1
      IPAR(86) = 0
      IPAR(87) = 3
      IPAR(88) = 1
      IPAR(89) = 0
      IPAR(90) = 1
      IPAR(91) = 0
      IPAR(92) = 1
      IPAR(93) = 1
c     ! force vectors (1) or strangeness (2) or baryons (3) after decays have been processed, rate is par75      
      IPAR(94) = 0              
      IPAR(95) = 1
      IPAR(96) = 0              ! energy dependence of enhancements (0: none, 1: logarithmic) threshold in par160
      IPAR(97) = 0              ! enable enhancement for meson projectiles only
      IPAR(98) = 0
      IPAR(99) = 1
      IPAR(100) = 0

C...  valence quark distribution function
c     large x suppression
      do i=1,3                  ! quark flavors
         CCHIK(i,13)=PAR(62)
         CCHIK(i,14)=PAR(62)
      enddo
C...string fragmentation parameters
c     effective quark mass
      STR_mass_val = PAR(36) 
      STR_mass_sea = PAR(41)

C...energy dependence of PTmin
c     pt_cut offset
      PAR(10) = PARS(10 , 1)
c     lambda
      PAR(11) = PARS(21 , 1)
c     c parameter
      PAR(12) = PARS(22 , 1)

C...fragmentation function
      FAin = PAR(20)
      FB0in = PAR(21)

C...Strange fragmentation function
      FAs1 = PAR(35)
      FAs2 = PAR(35)

C...leading baryon fragmentation function
c     hard proton mixing
      CLEAD = PAR(50)

      END
C=======================================================================

      SUBROUTINE PAR_INI_FROM_FILE
      IMPLICIT NONE
c     locals
      CHARACTER*10 FILENA
      CHARACTER*6 CNAME
      CHARACTER*70 NUMBER

      INTEGER ISTAT,J,IVAL,I
      DOUBLE PRECISION VAL      
c     commons
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      
      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      
      DOUBLE PRECISION STR_mass_val, STR_mass_val_hyp, STR_mass_sea
      COMMON /S_CUTOFF/ STR_mass_val, STR_mass_val_hyp, STR_mass_sea

      DOUBLE PRECISION FAin, FB0in
      COMMON /S_CZDIS/ FAin, FB0in

      DOUBLE PRECISION FAs1, fAs2
      COMMON /S_CZDISs/ FAs1, fAs2
      DOUBLE PRECISION ZDMAX, EPSI
      COMMON /S_CZDISc/ ZDMAX, EPSI

      DOUBLE PRECISION CLEAD, FLEAD
      COMMON /S_CZLEAD/ CLEAD, FLEAD
      DOUBLE PRECISION CCHIK
      COMMON /S_CPSPL/ CCHIK(4,99)

      SAVE
      DATA FILENA /'sibyll.par'/
 14   FORMAT(A6,A70)
 15   FORMAT(A5,I3,A2,I8)
 16   FORMAT(A5,I3,A2,F8.2)      
      OPEN(unit=4,file=filena,status='OLD')
      istat = 1
c     set standard parameters (full set)
      CALL PAR_INI
c     read new parameters from file
      IF(ndebug.gt.0)WRITE(LUN,*)'reading parameter file: sibyll.par'
      DO WHILE (istat.ge.0) 
         READ(4,14,iostat=ISTAT) CNAME,NUMBER
         IF(CNAME.eq.'IPAR  ')THEN
            READ(NUMBER,*) j, ival
            IF(ndebug.gt.1)write(LUN,15) 'IPAR(',j,')=', ival
            IPAR(J) = iVAL
         ELSEif(CNAME.eq.'PAR   ')THEN
            READ(NUMBER,*) j, val
            PAR(J) = VAL
            IF(ndebug.gt.1)write(LUN,16) ' PAR(',j,')=', val
         ELSE
            WRITE(LUN,*)'wrong format in parameter file!'
            WRITE(6,*)'wrong format in parameter file!'
            WRITE(LUN,*) CNAME, NUMBER
            stop
         ENDIF
      ENDDO
C     copy parameter values to their respective COMMONs
C...  valence quark distribution function
c     large x suppression
      do i=1,3                  ! quark flavors
         CCHIK(i,13)=PAR(62)
         CCHIK(i,14)=PAR(62)
      enddo
C...string fragmentation parameters
c     effective quark mass
      STR_mass_val = PAR(36) 
      STR_mass_sea = PAR(41)
C...fragmentation function
      FAin = PAR(20)
      FB0in = PAR(21)
C...Strange fragmentation function
      FAs1 = PAR(35)
      FAs2 = PAR(35)
C...leading baryon fragmentation function
c     hard proton mixing
      CLEAD = PAR(50)
      END
      
C=======================================================================
      
      SUBROUTINE MESON_FLV_MRG_INI

C-----------------------------------------------------------------------
c     change flavor merging for pions (favor spin)
C-----------------------------------------------------------------------
      INTEGER KFLV
      COMMON /S_KFLV/ KFLV(4,43)

c     pi+ --> rho+
      KFLV(2,1) = 25
c     pi- --> rho-
      KFLV(1,2) = 26
c     pi0 --> rho0
      KFLV(1,1) = 27
      KFLV(2,2) = 27     
      END
C=======================================================================

      BLOCK DATA PARAM_INI

C-----------------------------------------------------------------------
C....This block data contains default values
C.   of the parameters used in fragmentation
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      
      DOUBLE PRECISION STR_mass_val, STR_mass_val_hyp, STR_mass_sea
      COMMON /S_CUTOFF/ STR_mass_val, STR_mass_val_hyp, STR_mass_sea
      DOUBLE PRECISION FAin, FB0in
      COMMON /S_CZDIS/ FAin, FB0in

      DOUBLE PRECISION FAs1, fAs2
      COMMON /S_CZDISs/ FAs1, fAs2
      DOUBLE PRECISION ZDMAX, EPSI
      COMMON /S_CZDISc/ ZDMAX, EPSI

      DOUBLE PRECISION CLEAD, FLEAD
      COMMON /S_CZLEAD/ CLEAD, FLEAD
      DOUBLE PRECISION CCHIK
      COMMON /S_CPSPL/ CCHIK(4,99)

      INTEGER ITRY, NREJ
      COMMON /S_CNT/ ITRY(20), NREJ(20)

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN

      COMMON /CKFRAG/ KODFRAG
      SAVE

      DATA ITRY /20*0/
      DATA NREJ /5,5,20,10,20,20,14*0/ 
      DATA EPS3,EPS5,EPS8,EPS10 /1.D-3,1.D-5,1.D-8,1.D-10/
      DATA PI,TWOPI,CMBARN/3.14159265358979D0,6.283185308D0,0.389385D0/
      DATA FACN /2.D0,5.D0,15.D0,60.D0,250.D0,1500.D0,12000.D0,
     &     120000.D0/

C...  default output unit
      DATA LUN /7/
c...new fragmentation for charmed particles
      DATA EPSI /2.D0/
C...mass cutoff for soft strings
      data STR_mass_val /0.35D0/ 
      data STR_mass_val_hyp /0.4D0/ 
      data STR_mass_sea /1.D0/ 
C...Longitudinal Fragmentation function
      DATA FAin /0.5D0/, FB0in /0.8D0/
C...Longitudinal Fragmentation function for leading baryons
      DATA CLEAD  /0.6D0/, FLEAD  /0.6D0/
c     strange fragmentation
      data FAs1 /3.D0/, fAs2 /3.D0/
C...  Splitting parameters
      DATA CCHIK /20*0.D0,28*2.D0,8*3.D0,48*0.D0,4*2.D0,24*0.D0,
     &     24*3.D0,76*0.D0,8*2.D0,40*0.D0,12*2.D0,40*0.D0,24*3.D0,
     &     40*0.D0/
C...Parameters of flavor formation 
c     last in use: 160
      DATA PAR/0.04D0,0.3D0,0.3D0,0.14D0,0.3D0,0.3D0,0.15D0,0.D0,7.D0, ! 10
     &     2*0.D0,0.9D0,0.2D0,4*0.04D0,0.5D0,0.8D0,0.5D0,              ! 20 
     &     0.8D0,6.D0,0.5D0,0.004D0,5*0.D0,0.7D0,                      ! 30
     &     2*0.D0,0.1D0,0.D0,3.D0,0.35D0,0.D0,0.5D0,2*0.D0,            ! 40
     &     1.D0,2.D0,0.D0,0.99D0,0.D0,0.3D0,0.45D0,0.6D0,0.6D0,0.6D0,  ! 50
     &     .03D0,.03D0,6.D0,0.2D0,4*0.D0,1.1D0,0.8D0,                  ! 60
     &     .33D0,3.D0,1.D0,.25D0,.3D0,0.3D0,0.6D0,.007D0,.03D0,.007D0, ! 70
     &     1.D0,0.3D0,0.D0,0.3D0,0.0D0,0.2D0,0.5D0,1.0D0,10.D0,0.D0,   ! 80
     &     1000.D0,1000.D0,1.D0,6.D0,1.D0,0.D0,0.3D0,0.8D0,0.3D0,31.D0,! 90
     &     1.D0,6.5D0,1.D0,1.D0,0.D0,1.0D0,0.004D0,1.D0,0.33D0,1.D0,   ! 100
     &     1.D0,0.D0,2.D0,0.3D0,0.15D0,3*0.D0,20.D0,0.25D0,            ! 110
     &     0.D0,0.7D0,0.3D0,0.D0,0.D0,1.D0,3*0.D0,1.D0,                ! 120
     &     0.3D0,0.D0,0.3D0,1.D0,1.D0,1.D0,6.D0,1.D0,1.D0,6.D0,        ! 130
     &     0.0001D0,0.5D0,31.10362D0,-15.29012D0,6.5D0,                ! 135
     &     0.D0,4 *0.D0,                                               ! 140
     &     1.D0,0.D0,0.5D0,0.D0,0.5D0,0.D0,0.3D0,0.8D0,0.08D0,0.004D0, ! 150
     &     2.D0,1.D0,1.D0,1.D0,1.D0,1.D0,0.D0,1.D0,1.D0,0.D0,          ! 160      
     &     40*0.D0/                                                    ! 200
c     last in use:98
      DATA IPAR /9*0,1,0,1,8*0,20*0,    ! 40
     &     9*0,0,2,9*0,                 ! 60
     &     100,25*0,2,1,0,0,0,1,0,0,0,0,0,0,2*0/  ! 100   

C...Fragmentation of nuclei
      DATA KODFRAG /0/
C...Debug label and event counter
      DATA Ndebug /0/
      DATA Ncall /0/

      END

C=======================================================================

      SUBROUTINE PARAM_PRINT(LUN)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      COMMON /S_CZLEAD/ CLEAD, FLEAD
      COMMON /S_CPSPL/ CCHIK(4,99)
      COMMON /S_DEBUG/ Ncall, Ndebug, Lunn
      INTEGER NW_max
      PARAMETER (NW_max = 20)
C--------------------------------------------------------------------
C     SIBYLL common blocks containing event information       \FR'14
C--------------------------------------------------------------------
C     EVENT INFO COMMON
C     contains overall interaction properties, like
C     SQS : center-of-mass energy
C     S   :         "       "     squared
C     PTmin : low pt cut of QCD cross section, 
C             i.e. minimal pt of hard minijets
C     Xmin : low-x bound for PDFs, 
C            i.e. minimal momentum fraction of hard partons
C     Zmin : logarithm of that
C     KB : PID of beam hadron
C     KT() : PID of target
C     IAT : mass number of target
      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
      INTEGER KB,IAT,KT
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT
      DOUBLE PRECISION STR_mass_val, STR_mass_val_hyp, STR_mass_sea
      COMMON /S_CUTOFF/ STR_mass_val, STR_mass_val_hyp, STR_mass_sea
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      DOUBLE PRECISION FAin, FB0in
      COMMON /S_CZDIS/ FAin, FB0in

      DOUBLE PRECISION PPT02
      COMMON /S_CQDIS2/ PPT02(44)
      DOUBLE PRECISION PPT0,ptflag
      COMMON /S_CQDIS/ PPT0(35),ptflag
      SAVE

      WRITE (LUN, 25)
25      FORMAT( /,1x,40('-'), /
     +   ' SIBYLL MONTE CARLO PROGRAM. Version 2.3.f',/,
     +    1x,40('-'),/'  List of parameters: ' )

      WRITE (LUN, 31) FAin, FB0in
31      FORMAT ('  Parameters of longitudinal fragmentation: ', /,
     +          '   f(z) = (1-z)**a * exp(-b * mt**2/z) ', /,
     +          '   a = ', f9.3, 3x, ' b = ', f9.3, ' GeV**-2' )
      WRITE (LUN, 32) CLEAD, 1.D0/FLEAD-1.D0
32      FORMAT ('  Parameters of leading fragmentation: ', /,
     +   '   f(z) = c + (1-z)**a ', /,
     +   '   c = ',f9.3,3x,' a = ',f9.3) 

        WRITE (LUN, 33) str_mass_val, str_mass_sea
 33     FORMAT ('  Mass cuts ', /,
     +          '   val = ', f9.3, 3x, ' sea = ', f9.3, ' GeV' )

      WRITE (LUN, 35) PPT02(1), PPT02(3), PPT02(11),ppt02(10),ppt02(20)
35      FORMAT ('   <pT> of sea partons ', /,
     +   3x,'<pT>(u/d) ',F8.3,2x,'<pT>(s) ',f8.3,2x,'<pT>(qq) ',f8.3,
     +     2x,'<pT>(val) ',f8.3,2x,'<pT>(sea) ',f8.3)

      WRITE (LUN, 120) (PAR(K),K=1,24)
120      FORMAT (1x, ' Parameters of flavor formation: ',/,
     +   3x,'PAR(1) = Prob(qq)/Prob(q)              = ',F10.2,/,
     +   3x,'PAR(2) = Prob(s)/Prob(u)               = ',F10.2,/,
     +   3x,'PAR(3) = Prob(us)/Prob(ud)             = ',F10.2,/,
     +   3x,'PAR(4) = Prob(ud_0)/Prob(ud_1)         = ',F10.2,/,
     +   3x,'PAR(5) = Prob(Vector)/Prob(Scalar)     = ',F10.2,/,
     +   3x,'PAR(6) = Prob(K*)/Prob(K)              = ',F10.2,/,
     +   3x,'PAR(7) = Prob(spin 3/2)/Prob(spin=1/2) = ',F10.2,/,
     +   3x,'PAR(8) = Prob(B-M-Bbar)/Prob(B-Bbar)   = ',F10.2,/,
     +   3x,'PAR(9) = Phase space suppression of MI = ',F10.2,/,
     +   3x,'PAR(10)= Low-energy limit for pt cutoff= ',F10.2,/,
     +   3x,'PAR(11)= Pt cutoff factor for exp      = ',F10.2,/,
     +   3x,'PAR(12)= Pt cutoff factor for exp      = ',F10.2,/,
     +   3x,'PAR(13)= max. mass in diffraction      = ',F10.2,/,
     +   3x,'PAR(14)= Prob(qq)/Prob(q) std. value   = ',F10.2,/,
     +   3x,'PAR(15)= Prob(qq)/Prob(q) in hard jets = ',F10.2,/,
     +   3x,'PAR(16)= Prob(qq)/Prob(q) in diff.     = ',F10.2,/,
     +   3x,'PAR(17)= not used                      = ',F10.2,/,
     +   3x,'PAR(18)= not used                      = ',F10.2,/,
     +   3x,'PAR(19)= not used                      = ',F10.2,/,
     +   3x,'PAR(20)= not used                      = ',F10.2,/,
     +   3x,'PAR(21)= not used                      = ',F10.2,/,
     +   3x,'PAR(22)= effective scale in PDF (Q2)   = ',F10.2,/,
     +   3x,'PAR(23)= not used                      = ',F10.2,/,
     +   3x,'PAR(24)= Prob(s->c)                    = ',F10.2  )

      WRITE (LUN, 130) (IPAR(K),K=1,17)
130      FORMAT (1x, ' Model switches: ',/,
     +   3x,'IPAR(1) = not used                      =',I4,/,
     +   3x,'IPAR(2) = not used                      =',I4,/,
     +   3x,'IPAR(3) = exponential pt                =',I4,/,
     +   3x,'IPAR(4) = decouple qq/q in val. strings =',I4,/,
     +   3x,'IPAR(5) = decouple qq/q in hm. diff.    =',I4,/,
     +   3x,'IPAR(6) = decouple qq/q in hard strings =',I4,/,
     +   3x,'IPAR(7) = remnant (not implemented yet) =',I4,/,
     +   3x,'IPAR(8) = jet kinematic pdf set (DO/GRV)=',I4,/,
     +   3x,'IPAR(9) = smear lowest diff. mass       =',I4,/,
     +   3x,'IPAR(10)= high mass diff. mode (d:ON)   =',I4,/,
     +   3x,'IPAR(11)= leading vec. meson prod. model=',I4,/,
     +   3x,'IPAR(12)= inel. screening in pAir       =',I4,/,
     +   3x,'IPAR(13)= decouple qq/q in val. strings =',I4,/,
     +   3x,'IPAR(14)= fireball model                =',I4,/,
     +   3x,'IPAR(15)= charm production              =',I4,/,
     +   3x,'IPAR(16)= charmed transverse momentum   =',I4,/,
     +   3x,'IPAR(17)= full charm model              =',I4 )

      WRITE (LUN, 40)
      WRITE (LUN, 41) CCHIK (1,13), CCHIK(2,13)
 40   FORMAT(' Parameters of hadron splitting ' )
 41   FORMAT('   p -> [(ud) u] splitting: alpha = ', F10.3, /,
     +       '   p -> [(uu) d] splitting: alpha = ', F10.3 )
c     print rho0 splitting
      WRITE (LUN, 42) CCHIK (1,27), CCHIK(2,27)        
 42   FORMAT(' rho0 -> [u ubar] splitting: alpha = ', F10.3, /,
     +       ' rho0 -> [d dbar] splitting: alpha = ', F10.3 )
c     print d+ splitting
      WRITE (LUN, 43) CCHIK (4,59), CCHIK(2,59)        
 43   FORMAT('  dp -> [c ubar] splitting: alpha = ', F10.3, /,
     +       '  dp -> [dbar c] splitting: alpha = ', F10.3 )
      END

C=======================================================================

      SUBROUTINE SIB_LIST(LUN)

C-----------------------------------------------------------------------
C...This routine prints the event record for the
C.  current event on unit LUN
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      COMMON /S_DEBUG/ Ncall, Ndebug, Lunn
      INTEGER NW_max
      PARAMETER (NW_max = 20)
C--------------------------------------------------------------------
C     SIBYLL common blocks containing event information       \FR'14
C--------------------------------------------------------------------
C     EVENT INFO COMMON
C     contains overall interaction properties, like
C     SQS : center-of-mass energy
C     S   :         "       "     squared
C     PTmin : low pt cut of QCD cross section, 
C             i.e. minimal pt of hard minijets
C     Xmin : low-x bound for PDFs, 
C            i.e. minimal momentum fraction of hard partons
C     Zmin : logarithm of that
C     KB : PID of beam hadron
C     KT() : PID of target
C     IAT : mass number of target
      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
      INTEGER KB,IAT,KT
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT
C     The final particle output is contained in COMMON /S_PLIST/    
C     NP           : number of final particles
C     P(1:NP, 1:5) : 4-momenta + masses of the final particles 
C     LLIST (1:NP) : codes of final particles
      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP

      INTEGER NFORIG,NPORIG,NIORIG,IPFLAG,IIFLAG,KINT
      COMMON /S_PARTO/ NFORIG(NP_max),NPORIG(NP_max),NIORIG(NP_max),
     &IPFLAG,IIFLAG,KINT
      INTEGER LLIST1
      COMMON /S_PLIST1/ LLIST1(8000)
C     parameters that represent: NW: max. number of wounded nucleons,
C     NS,NH: max. number of soft and hard interactions
c      PARAMETER (NW_max = 20)
C     The COMMON block /S_CHIST/ contains information about the
C     the structure of the  generated event:
C     NWD   = number of wounded nucleons
C     NJET = total number of hard interactions
C     NSOF = total number of soft interactions
C     NNSOF (1:NW) = number of soft pomeron cuts in each interaction
C     NNJET (1:NW) = number of minijets produced in each interaction 
C     JDIF(1:NW) = diffraction code 
C                  0 : non-diff,
C                  1 : beam-diff
C                  2 : target-diff
C                  3 : double-diff
      INTEGER NNSOF,NNJET,JDIF,NWD,NJET,NSOF
      COMMON /S_CHIST/ NNSOF(NW_max),NNJET(NW_max),
     &     JDIF(NW_max),NWD,NJET,NSOF
      INTEGER ICHP,ISTR,IBAR
      COMMON /S_CHP/ ICHP(99), ISTR(99), IBAR(99)
      INTEGER IISO,ISPN
      COMMON /S_SPN/ IISO(99), ISPN(99)
      INTEGER ICHM
      COMMON /S_CHM/ ICHM(99)
      CHARACTER*6 NAMP
      COMMON /S_CNAM/ NAMP (0:99)
      INTEGER IRMNT,KRB,KRT
      DOUBLE PRECISION XRMASS,XRMEX
      COMMON /S_RMNT/ XRMASS(2),XRMEX(2),IRMNT(NW_max),KRB,KRT(NW_max)

      CHARACTER*7 CTGT(0:20)
      CHARACTER CODE*18
      CHARACTER*18 NAMDIF(0:3)
      CHARACTER*18 NAMRMNT(0:3)
      SAVE
      DATA CTGT /'Air    ','Proton ',19*'Nucleus'/
      DATA NAMDIF /'Non-diff. event   ',
     &  'Beam diffraction  ','Target diffraction','Double diffraction'/
      DATA NAMRMNT /'No resolvd remnant',
     &  'Beam remnant     ','Target remnant    ','Double remnant    '/

 50   FORMAT(3X,88('-'),/,25X,'SIBYLL EVENT SUMMARY',25X,
     &     /,3X,88('-'))
 52   FORMAT( 3X,'Beam + Target @ Energy:',2X,A6,2X,'+',2X,A7,2X,
     &     '@',1p,E11.3,' GeV')
 53   FORMAT( 3X,'Beam + Target @ Energy:',2X,'Anti-',A6,2X,'+',2X,A7,
     &     2X,'@',1p,E11.3,' GeV')

      WRITE (LUN,50)
      IF (KB .GT. 0 ) THEN
        WRITE (LUN,52)
     &     NAMP(IABS(KB)),CTGT(IAT),SQS
      ELSE 
        WRITE (LUN,53)
     &     NAMP(IABS(KB)),CTGT(IAT),SQS
      ENDIF
      if(NWD.eq.1)THEN
         WRITE (LUN,*) '  ',NAMDIF(JDIF(1))
         IF(jdif(1).eq.0)
     &    WRITE (LUN,*) '  ',NAMRMNT(abs(IRMNT(1)))
      else
         WRITE (LUN,*) '  ',NAMDIF(0)
      endif

      WRITE (LUN,*) '  A/N_w/N_s/N_j = ', IAT , NWD, NSOF, NJET
      WRITE (LUN,100)

C...Print particle list
      kchar = 0
      ibary = 0
      ichmd = 0
      istrg = 0
      DO J=1,NP
        L = MOD(LLIST(J),10000)
        CODE = '                  '
        CODE(1:6) = NAMP(IABS(L))
        IF (L .LT. 0) CODE(7:9) = 'bar'
        IF(IABS(LLIST(J)) .GT. 10000)   CODE(10:10) = '*'
        WRITE (LUN,120) J, CODE, NIORIG(J),JDIF(NIORIG(J)),LLIST1(J), 
     &       NPORIG(J), (P(J,K),K=1,4)
        if(abs(LLIST(J)).LT.10000) then
          kchar = kchar+sign(1,l)*ICHP(iabs(l))
          ibary = ibary+sign(1,l)*IBAR(iabs(l))
          ichmd = ichmd+sign(1,l)*ICHM(iabs(l))
          istrg = istrg+sign(1,l)*ISTR(iabs(l))
        endif
      ENDDO
      CALL PFSUM(1,NP,Esum,PXsum,PYsum,PZsum,NF)
      WRITE(LUN,140) PXsum,PYsum,PZsum,Esum
100      FORMAT(3X,'N  Particle',12X,'Int',2x,'Jdif',2x,'Prnt',2x,'Proc'
     +         ,6x,'PX',9x,'PY',9x,'PZ',9x,'E', /, 3X,88('-'))
120      FORMAT(I6,1X,A18,3I5,I8,2F10.3,1p,2E11.3)
140      FORMAT(3X,88('-'),/,'  Tot =',41X,2F10.3,1p,2E11.3)
         write(LUN,'(1x,a,i3,3x,a,i3)') ' Total charge:     ',kchar,
     &        'total baryon number:',ibary
         write(LUN,'(1x,a,i3,3x,a,i3)') ' Total strangeness:',istrg,
     &        'total charm number: ',ichmd

      RETURN
      END

C=======================================================================

      SUBROUTINE KCODE (J,CODE,NC)

C-----------------------------------------------------------------------
C...Produce the code for parton J
C.  Input K, Output CODE, NC=number of characters
C..................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      CHARACTER*5 CODE
      CHARACTER*1 NAMQ(4)
      SAVE
      DATA NAMQ /'u','d','s','c'/

      CODE = '     '
      IF(J.EQ.0)  THEN
         CODE(1:3) = 'glu'
         NC = 3
         RETURN
      ENDIF
      JA = IABS(J)
      J1 = MOD(JA,10)
      J2 = (JA-J1)/10
      IF(JA .GT. 10) THEN
         CODE(1:1) = NAMQ(J2)
         CODE(2:2) = NAMQ(J1)
         NC = 2
      ELSE
         CODE(1:1) = NAMQ(J1)
         NC = 1      
      ENDIF
      IF (J .LT. 0)  THEN
         CODE(NC+1:NC+3) = 'bar'
         NC = NC+3
      ENDIF
      RETURN
      END
C=======================================================================

      SUBROUTINE SIB_PARTPR(LUN)

C----------------------------------------------------------------
C     prints the particles known to SIBYLL with their internal
C     and PDG labels \FR'13
C----------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)

      INTEGER ICHP,ISTR,IBAR
      COMMON /S_CHP/ ICHP(99), ISTR(99), IBAR(99)

      INTEGER ICHM
      COMMON /S_CHM/ ICHM(99)

      CHARACTER*6 NAMP
      COMMON /S_CNAM/ NAMP (0:99)
      SAVE

      WRITE(LUN,50)
 50   FORMAT(/,2X,16X,'SIBYLL PARTICLE TABLE:',/,2x,80('-'))
      WRITE(LUN,100)
 100  FORMAT(2X,'Particle',4X,'SIB PID',6x,'SIB2PDG',6x,'SIB2PDG^-1', 
     &     4x,'MASS',4x,'STRG',4x,'CHRM',4x,'BRYN'/, 2X,80('-'))

      DO J=1,99
         IA = ISIB_PID2PDG( j )         
         IF(IA.ne.0)THEN
            ISIBPDG2PIDIA=ISIB_PDG2PID( IA )
         ELSE
            WRITE(LUN,'(1X,A,I2)') 'PDG conversion not found! pid=', j
         ENDIF         
         WRITE (LUN,120)  NAMP(J), J, IA, ISIBPDG2PIDIA, AM(J), ISTR(J),
     &        ICHM(J), IBAR(J)
      ENDDO
 120  FORMAT(4X,A6,4X,I4,7X,I7,8X,I4,5X,F9.3,3(6X,I2))

      END

C=======================================================================

      INTEGER FUNCTION ISIB_PID2PDG(Npid)

C----------------------------------------------------------------
C     conversion of SIBYLL internal particle code to PDG standard
C
C     input:     Npid        internal particle number
C     output:    sib_pid2pdg  PDG particle number
C
C     based on similar phojet function \FR'13
C----------------------------------------------------------------
      COMMON /S_PDG2PID/ ID_PDG_LIST(99),ID_LIST(577)
      INTEGER NPIDA,NPID
      SAVE

      Npida = iabs(Npid)
      ISIB_PID2PDG = ID_PDG_LIST(Npida)
      IF(NPID.lt.0)ISIB_PID2PDG = isign(ISIB_PID2PDG,Npid)
      RETURN
      END

C=======================================================================

      INTEGER FUNCTION ISIB_PDG2PID(Npdg)

C-----------------------------------------------------------------------
C     conversion of PDG standard particle code to SIBYLL internal
C
C     input:     Npdg        PDG particle number
C     output:    sib_pdg2pid internal particle id
C
C     based on similar phojet function \FR'13
C----------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON /S_PDG2PID/ IPID_PDG_LIST(99),ID_LIST(577)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN

      DOUBLE PRECISION CBR
      INTEGER KDEC,LBARP,IDB
      COMMON /S_CSYDEC/ CBR(223+16+12+8), KDEC(1338+6*(16+12+8)),
     &     LBARP(99), IDB(99)
      SAVE

      Nin = abs(Npdg)
      if((Nin.gt.999999).or.(Nin.eq.0)) then
C  invalid particle number
        if(ndebug.gt.5) write(6,'(1x,A,I10)')
     &    ' ISIB_PDG2PID: invalid PDG ID number ',Npdg
        ISIB_PDG2PID = 0
        return
      else If(Nin.le.577) then
C  simple case
        Nout = Nin
      else
C  use hash algorithm
        Nout = mod(Nin,577)
      endif

 100  continue

C  particle not in table
      if(ID_list(Nout).Eq.0) then
         if(ndebug.gt.0) write(6,'(1x,A,I10)')
     &    ' ISIB_PDG2PID: particle not in table ',Npdg
        ISIB_PDG2PID = 0
        return
      endif
      ID_out = ID_list(Nout)
      IF(abs(ID_out).gt.99)then
         ISIB_PDG2PID = 0
         return
      else

         if(IPID_PDG_LIST(ID_list(Nout)).eq.Nin) then
C     particle ID found
            ISIB_PDG2PID = ID_list(Nout)
            if (NPDG.lt.0) ISIB_PDG2PID = lbarp( ISIB_PDG2PID )
            return
         else
C     increment and try again
            Nout = Nout + 5
            If(Nout.gt.577) Nout = Mod(Nout,577)
            goto 100
         endif
      endif
      END

C=======================================================================

      SUBROUTINE PDG_INI

C-----------------------------------------------------------------------
C     PDG conversion blocks \FR'13
C----------------------------------------------------------------
      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      PARAMETER ( ID_PDG_MAX = 260 )
      COMMON /S_PDG2PID/ ID_PDG_LIST(99),ID_LIST(577)
      SAVE
      DATA ID_PDG_LIST /22,-11,11,-13,13,111,211,-211,321,-321, !10
     &     130,310,2212,2112,12,-12,14,-14,-2212,-2112,         !20
     &     311,-311,221,331,213,-213,113,323,-323,313,          !30
     &     -313,223,333,3222,3212,3112,3322,3312,3122,2224,     !40
     &     2214,2114,1114,3224,3214,3114,3324,3314,3334,0,      !50
     &     202212,202112,212212,212112,4*0,411,-411,            !60
     &     900111,900211,-900211,7*0,                           !70
     &     421,-421,441,431,-431,433,-433,413,-413,423,         !80
     &     -423,0,443,4222,4212,4112,4232,4132,4122,-15,        !90
     &     15,-16,16,4224,4214,4114,4324,4314,4332/

      IF(Ndebug.gt.2)
     & WRITE(LUN,*) ' INITIALIZING PDG TABLES..'
      CALL SIB_CPCINI(ID_pdg_max,ID_pdg_list,ID_list)
      
      END

C=======================================================================

      SUBROUTINE SIB_CPCINI(Nrows,Number,List)

C-----------------------------------------------------------------------
C     initialization of particle hash table
C
C     input:   Number     vector with Nrows entries according to PDG
C                         convention
C
C     output:  List       vector with hash table
C
C     (this code is based on the function initpns written by
C      Gerry Lynch, LBL, January 1990)
C
C***********************************************************************
      IMPLICIT NONE

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      integer Number(*),List(*),Nrows
      Integer Nin,Nout,Ip,I
      SAVE

      do I = 1,577
        List(I) = 0
      enddo

C    Loop over all of the elements in the Number vector

        Do 500 Ip = 1,Nrows
            Nin = Number(Ip)

C    Calculate a list number for this particle id number
            If(Nin.Gt.999999.or.Nin.Le.0) Then
                 Nout = -1
            Else If(Nin.Le.577) Then
                 Nout = Nin
            Else
               Nout = Mod(Nin,577)
            End If

 200        continue

            If(Nout.Lt.0) Then
C    Count the bad entries
               IF(Ndebug.gt.3) Write(LUN,'(1x,a,i10)')
     &            ' SIB_CPCINI: invalid particle ID',Nin
               Go to 500
            End If
            If(List(Nout).eq.0) Then
                List(Nout) = Ip
            Else
                If(Nin.eq.Number(List(Nout))) Then
                  IF(Ndebug.gt.3)Write(LUN,'(1x,a,i10)')
     &              ' SIB_CPCINI: double particle  ID',Nin
                End If
                Nout = Nout + 5
                If(Nout.Gt.577) Nout = Mod(Nout, 577)

                Go to 200
            End If
 500    Continue

      END
C=======================================================================

      SUBROUTINE PFSUM(N1,N2,ETOT,PXT,PYT,PZT,NF)

C-----------------------------------------------------------------------
C...Return the energy,px,py,pz and the number of stable
C.  particles in the list between N1 and N2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
c      COMMON /S_PLIST/ P(8000,5), LLIST(8000), NP
C     The final particle output is contained in COMMON /S_PLIST/    
C     NP           : number of final particles
C     P(1:NP, 1:5) : 4-momenta + masses of the final particles 
C     LLIST (1:NP) : codes of final particles
      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP
      SAVE

      NF=0
      ETOT=0.D0
      PXT=0.D0
      PYT=0.D0
      PZT=0.D0
      DO J=N1,N2
         L = LLIST(J)     
         IF (IABS(L) .LT. 10000)  THEN
           NF = NF+1
           ETOT = ETOT + ABS( P(J,4) )
           PXT = PXT + P(J,1)
           PYT = PYT + P(J,2)
           PZT = PZT + P(J,3)
         ENDIF
      ENDDO
      RETURN
      END

C=======================================================================

      SUBROUTINE RESAMPLE_FPHI(N1,N2)

C-----------------------------------------------------------------------
C...Resample the azimuth angle, resetting px,py
C.  particles in the list between N1 and N2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
c      COMMON /S_PLIST/ P(8000,5), LLIST(8000), NP
C     The final particle output is contained in COMMON /S_PLIST/    
C     NP           : number of final particles
C     P(1:NP, 1:5) : 4-momenta + masses of the final particles 
C     LLIST (1:NP) : codes of final particles
      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN
      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      SAVE

      PXT=0.D0
      PYT=0.D0
      DO J=N1,N2
         L = LLIST(J)     
         IF (IABS(L) .LT. 10000)  THEN
           NF = NF+1
           PHI = TWOPI*S_RNDM(J)
           PT = dsqrt(P(J,1)**2 + P(J,2)**2)
c     new px and py
           P(J,1) = PT*dCOS(PHI)
           P(J,2) = PT*dSIN(PHI)
           PXT = PXT + P(J,1)
           PYT = PYT + P(J,2)
         ENDIF
      ENDDO
      IF(NDEBUG.gt.0)THEN
         WRITE(LUN,*) 'resampled azimuthal angle'
         WRITE(LUN,*) 'total shift in px, py:', PXT, PYT
      ENDIF
      RETURN
      END

C=======================================================================

      SUBROUTINE QNUM (JQ,JS,JC,JB,JBA, NC, NF)

C-----------------------------------------------------------------------
C...Return the quantum numbers of one event
C.  JQ = charge, JB = baryon number, JS = strangeness, JC = charmedness
C.  JBA = (number of baryons+antibaryons)
C.  NC  = number of charged particles
C.  NF  = number of final particles
C..................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
C     The final particle output is contained in COMMON /S_PLIST/    
C     NP           : number of final particles
C     P(1:NP, 1:5) : 4-momenta + masses of the final particles 
C     LLIST (1:NP) : codes of final particles
      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP

      INTEGER ICHP,ISTR,IBAR
      COMMON /S_CHP/ ICHP(99), ISTR(99), IBAR(99)

      INTEGER IISO,ISPN
      COMMON /S_SPN/ IISO(99), ISPN(99)

      INTEGER ICHM
      COMMON /S_CHM/ ICHM(99)
      SAVE

      JQ = 0
      JB = 0
      JS = 0
      JC = 0
      JBA= 0
      NC = 0
      NF = 0
      DO J=1,NP
          L = LLIST(J)
          LL = IABS(L)
          IF (LL .LT. 10000)  THEN
              IF(ICHP(LL) .NE. 0) NC = NC + 1
              NF = NF + 1
              JQ = JQ + ICHP(LL)*ISIGN(1,L)
              JB = JB + IBAR(LL)*ISIGN(1,L)
              JBA= JBA+ IBAR(LL)
              JS = JS + ISTR(LL)*ISIGN(1,L)
              JC = JC + ICHM(LL)*ISIGN(1,L)
          ENDIF
      ENDDO
      RETURN
      END

C=======================================================================

      BLOCK DATA KFLV_INI

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER KFLV
      COMMON /S_KFLV/ KFLV(4,43)
      SAVE
      DATA (KFLV(1,i),i=1,4) /6,8,10,71/
      DATA (KFLV(1,i),i=5,43) /6*0,40,13,34,84,6*0,13,14,39,89,6*0,
     &     34,39,37,87,6*0,84,85,87/      
      DATA (KFLV(2,i),i=1,4) /7,6,21,59/
      DATA (KFLV(2,i),i=5,43) /6*0,13,14,39,89,6*0,14,43,36,86,6*0,
     &     39,36,38,88,6*0,84,85,87/     
      DATA (KFLV(3,i),i=1,4) /9,22,33,74/
      DATA (KFLV(3,i),i=5,43) /6*0,34,39,35,87,6*0,39,36,38,88,6*0,
     &     35,36,49,99,6*0,84,85,87/
      DATA (KFLV(4,i),i=1,4) /72,60,75,83/
      DATA (KFLV(4,i),i=5,43) /6*0,84,85,87,0,6*0,85,86,88,0,6*0,
     &     87,88,99,0,6*0,0,0,0/

      END
C=======================================================================

      SUBROUTINE SIB_I4FLAV (IFL1, IFL2_A, IRNK, IFL2, KF)

C-----------------------------------------------------------------------
C.  This subroutine receives as input IFL1 the flavor code
C.  of a quark (antiquark) and  generates the antiquark (quark)
C.  of flavor code IFL2 that combine with the original parton
C.  to compose an hadron of code KF.
C.
C.  updated to 4 FLAVORS \FR'13
C.  Baryon sector is from jetset code
C.  assuming D*_s+- are J=1, only Charm=1 baryons
C.
C.  If (IFL2_A.NE.0) returns an hadron KF composed of IFL1 and IFL2_A
c
c     Input: IFL1 - flavor of first quark
c            IFL2_A - flavor of second quark ( if 0 randomly chosen ) 
c            IRNK - position in hadron chain
c     Output: IFL2 - flavor of second quark partner to be passed on
c             KF - final hadron
C-----------------------------------------------------------------------
Cf2py integer,intent(out) :: ifl2
Cf2py integer,intent(out) :: kf
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      IMPLICIT INTEGER(I-N)

      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN

      DIMENSION KFLA(4,4,2), CDIAG(16), KDIAG(8)
      DIMENSION KBAR(40), CFR(28), KFR(80)
      SAVE
      DATA KFLA /0,8,10,71,7,0,22,59,9,21,0,74,72,60,75,0, ! spin-zero mesons
     +     0,26,29,80,25,0,31,78,28,30,0,76,81,79,77,0/ ! spin-one mesons
      DATA CDIAG/.5D0,.25D0,.5D0,.25D0,1.D0,.5D0,2.D0,1.D0, !spin-zero diagonal mesons
     +     .5D0,0.D0,.5D0,0.D0,1.D0,1.D0,2.D0,1.D0/ ! spin-one diagonal mesons
      DATA KDIAG /6,23,24,73,27,32,33,83/
      DATA KBAR /13,14,34,35,36,37,38,84,85,86, !jetset -> sibyll part. code map
     +     87,88,99,3*0,39,89,87,88, 
     +     40,41,42,43,44,45,46,47,48,49,      
     +     94,95,96,97,98,99,4*0/ ! spin-3/2 css baryon added to 1/2 css
      DATA CFR /0.75D0,0.D0,0.5D0,0.D0,0.D0,1.D0,0.1667D0,0.3333D0,
     $          0.0833D0,0.6667D0,0.1667D0,0.3333D0,-3.D0,1.D0,-2.D0,
     $          -2.D0,1.D0,0.D0,0.D0,-3.D0,1.D0,1.D0,1.D0,5*0.D0/
      DATA KFR/0,16,17,19,100,104,109,115,0,26,27,29,122,126,131,137
     +  ,0,40,42,47,144,158,178,205,0,1,3,6,10,15,21,28,0,0,56,57,240,
     +  246,256,271,0,0,1,3,6,10,15,21,60,61,64,70,292,307,328,356,
     +  0,1,3,6,10,15,21,28,16*0/

      IF(NDEBUG.gt.6)
     &     WRITE(LUN,*)' SIB_FLAV: input:',IFL1, IFL2_A, IRNK, IFL2, KF

c     set rho0 / ( omega, phi ) ratio, i.e. I=1 to I=0
c     default: 0.5, 0.0 ( phi only created from s-sbar)
      CDIAG(8+1) =  1.D0-PAR(143) ! u-flavor, Prob. I=1 vs 0
      CDIAG(8+3) =  1.D0-PAR(143) ! d-flavor, Prob. I=1 vs 0

      XDIQ = 1.D0
      
      IARNK = IABS(IRNK)
      IFLA = IABS(IFL1)
c     check if diq production allowed?
c     for strings with leading diquarks the immediate formation of another diquark may be forbidden
      if(ifla.gt.100.and.mod(ifla,100).lt.10)then
         XDIQ = PAR(158)
         ifl1 = mod(ifl1,100)
         IFLA = IABS(IFL1)
      endif
      
      IFL2A = IFL2_A
      IF (IFL2A .NE. 0)  THEN
c     combine existing flavors to hadron
c     three cases: input diquark (MB=2): need to sample additional quark,
c                  input quark (MB=0,1): sample quark (0) or diquark (1)?         
         IFL2A = MOD(IFL2A,100)
         IFL2 = IFL2A
         IFLB = IABS(IFL2A)
         MB = 0
         IF (IFLB .GT. 10)   MB=1
         IF (IFLA .GT. 10)   MB=2
      ELSE
c     sample new flavor
         MB = 2
         IF (IFLA .LT. 10)   THEN
             MB = 1
             IF ((1.D0+PAR(1))*S_RNDM(0).LT. 1.D0)  MB=0
             XDIQ = 1.D0             
c     suppress baryons close to the string end
c     IPAR(55) defines largest forbidden rank
c     PAR(101) is the rejection probability
             IF (IPAR(54).eq.1)THEN
                IF(IARNK.le.IPAR(55).and.S_RNDM(1).lt.PAR(101)) MB=0
             ENDIF
         ENDIF
      ENDIF

 50   IF (MB .EQ. 0)  THEN
c     flavor open, sample from u,d,s,c
         IF (IFL2A.EQ.0)THEN
            IF(IPAR(69).eq.2)THEN
c     asymmetric between u,d
               IFL2 = MIN(2,1+INT((2.D0+PAR(115))*S_RNDM(0)))
               IFLS = 3*INT(INT((2.D0+PAR(2))*S_RNDM(1))*0.5D0)
               IFL2 = MAX(IFL2,IFLS)
               IFL2 = ISIGN(IFL2,-IFL1)
            ELSE
c     symmetric in u,d
               IFL2=ISIGN(1+INT((2.D0+PAR(2))*S_RNDM(0)),-IFL1)
            ENDIF
            IF(IABS(IFL2).eq.3) THEN
               IF(S_RNDM(1).lt.PAR(24)*PAR(125))
     +              IFL2=ISIGN(4,-IFL1)
            ENDIF
         ENDIF
         IFLD = MAX(IFL1,IFL2)
         IFLE = MIN(IFL1,IFL2)
         GOTO 100
      ENDIF

C...  Decide if the diquark must be split
c     if diquark is from previous splitting (popcorn) do NOT split diquark
c     jump to sample quark and form baryon      
      IF (MB .EQ. 2 .AND. IFLA .GT. 100)   THEN
         IFLA = MOD(IFLA,100)
           GOTO 200
      ENDIF
c     split diquark? if yes sample single flavor and form meson
c     diquark with any flavor combination is passed on with id+100        
      IF (MB .EQ. 2 .AND. IFL2A .EQ. 0)   THEN
          IF (S_RNDM(0) .LT. PAR(8))  THEN
             MB = 0
             IFLG = MOD(IFL1,10)
             IFLH =(IFL1-IFLG)/10
             IF (S_RNDM(1) .GT. 0.5D0)  THEN
                IFLDUM = IFLG
                IFLG = IFLH
                IFLH = IFLDUM
             ENDIF
             IFL11=IFLG
             IFL22=ISIGN(1+INT((2.D0+PAR(2))*S_RNDM(2)),-IFL1)
             IFLD = MAX(IFL11,IFL22)
             IFLE = MIN(IFL11,IFL22)
             IFL2 = -IFLH*10+IFL22
             IF (S_RNDM(3) .GT. 0.5D0)  IFL2 = IFL22*10-IFLH
c     limit diquark splitting to B-M-B (default: yes)             
             IF(IPAR(92).eq.1) IFL2 = IFL2+ISIGN(100,IFL2)
          ENDIF
      ENDIF
       
C...Form a meson: consider spin and flavor mixing for the diagonal states
 100  IF (MB .EQ. 0)  THEN
         IF1 = IABS(IFLD)
         IF2 = IABS(IFLE)
         IFLC = MAX(IF1,IF2)
         KSP = INT(PAR(5)+S_RNDM(0))
         KSP = MIN(KSP,1)
         IF (IFLC.EQ.3)  KSP = INT(PAR(6)+S_RNDM(1))
         IF (IFLC.EQ.4)  KSP = INT(PAR(6)+S_RNDM(2))
         IF (IF1 .NE. IF2)   THEN
            KF = KFLA(IF1,IF2,KSP+1)
         ELSE
            R = S_RNDM(0)
            JF=1+INT(R+CDIAG(8*KSP+2*IF1-1))+
     +             INT(R+CDIAG(8*KSP+2*IF1))
            JF = MIN(JF,4)
            KF=KDIAG(JF+4*KSP)
c     suppress neutral pions
            IF(KF.eq.6)THEN
               IF(IPAR(82).eq.1.and.
     +              S_RNDM(kf).lt.PAR(137))then
                  IF(IFL2A.ne.0) goto 100
                  IF(IFLA.gt.10) mb = 2                  
                  GOTO 50
               endif
c     suppress neutral pions, depending on rank
               IF(IPAR(82).eq.2.and.S_RNDM(3).lt.PAR(137).and.
     +              irnk.gt.0.and.irnk.lt.2) then
                  IF(IFL2A.ne.0) goto 100
                  IF(IFLA.gt.10) mb = 2
                  GOTO 50
               endif
            ENDIF
c     suppress rank1 (leading) omega
            IF(KF.eq.32)THEN
               IF(IPAR(83).ne.0.and.
     +              S_RNDM(kf).lt.PAR(138))then
                  IF(IFL2A.ne.0) goto 100
                  IF(IFLA.gt.10) mb = 2                  
                  GOTO 50
               endif
            ENDIF
         ENDIF
c         PRINT*,' I4FLAV returns :(IFL1,IFL2,LL)',IFL1,IFL2,KF
        IF(NDEBUG.gt.6)
     &     WRITE(LUN,*)' SIB_FLAV: output:',IFL1, IFL2_A, IRNK, IFL2, KF
         RETURN
      ENDIF

C...Form a baryon
 200  IF (IFL2A .NE. 0)   THEN
         IF (MB .EQ. 1)  THEN
            IFLD = IFLA
            IFLE = IFLB/10
            IFLF = MOD(IFLB,10)
         ELSE
            IFLD = IFLB
            IFLE = IFLA/10
            IFLF = MOD(IFLA,10)
         ENDIF
         LFR = 3+2*((2*(IFLE-IFLF))/(1+IABS(IFLE-IFLF)))
         IF(IFLD.NE.IFLE.AND.IFLD.NE.IFLF)  LFR=LFR+1
      ELSE
 110     CONTINUE
         IF(MB.EQ.1)   THEN     ! generate diquark
            IFLD = IFLA
 120        IFLE = 1+INT((2.D0+PAR(2)*PAR(3))*S_RNDM(0))
            IFLF = 1+INT((2.D0+PAR(2)*PAR(3))*S_RNDM(1))          
            IF(IFLD.NE.4)THEN
               IF(IFLE.EQ.3)THEN 
                  IF(S_RNDM(2).lt.PAR(24)*PAR(125))
     +                 IFLE=4
               ENDIF
               IF(IFLF.EQ.3.and.IFLE.NE.4)THEN 
                  IF(S_RNDM(3).lt.PAR(24)*PAR(125))
     +                 IFLF=4
               ENDIF
            ENDIF
            IF(IFLE.GE.IFLF.AND.PAR(4).LT.S_RNDM(4))    GOTO 120
            IF(IFLE.LT.IFLF.AND.PAR(4)*S_RNDM(5).GT.1.D0) GOTO 120     
            IFL2=ISIGN(10*IFLE+IFLF,IFL1)
         ELSE                   ! generate quark
            IF(IPAR(69).eq.2)THEN
c     asymmetric between u,d
               IFL2 = MIN(2,1+INT((2.D0+PAR(115))*S_RNDM(6)))
               IFLS = 3*(INT((2.D0+PAR(2))*S_RNDM(7))/2)
               IFL2 = MAX(IFL2,IFLS)
               IFL2 = ISIGN(IFL2,IFL1)
            ELSE
c     symmetric in u,d
               IFL2=ISIGN(1+INT((2.D0+PAR(2))*S_RNDM(8)),IFL1)
            ENDIF
            IFLE=IFLA/10
            IFLF=MOD(IFLA,10)
            IF(IABS(IFL2).EQ.3.and.IFLF.ne.4.and.IFLE.ne.4) THEN
               IF(S_RNDM(9).lt.PAR(24)*PAR(125))
     +              IFL2=ISIGN(4,IFL1)
            ENDIF
            IFLD=IABS(IFL2)
         ENDIF
C...SU(6) factors for baryon formation
         LFR=3+2*((2*(IFLE-IFLF))/(1+IABS(IFLE-IFLF)))
         IF(IFLD.NE.IFLE.AND.IFLD.NE.IFLF)  LFR=LFR+1
         WT = CFR(2*LFR-1)+PAR(7)*CFR(2*LFR)
         IF(IFLE.LT.IFLF)  WT=WT/3.D0
         IF (WT.LT.S_RNDM(0)) GOTO 110
      ENDIF

C...Form Baryon
      IFLG=MAX(IFLD,IFLE,IFLF)
      IFLI=MIN(IFLD,IFLE,IFLF)
      IFLH=IFLD+IFLE+IFLF-IFLG-IFLI
c      IF(IFLG+IFLH.gt.7) GOTO 200 ! forbid double charmed
      KSP=2+2*INT(1.D0-CFR(2*LFR-1)+(CFR(2*LFR-1)+PAR(7)*
     1       CFR(2*LFR))*S_RNDM(0))

C...Distinguish Lambda- and Sigma- like particles
      IF (KSP.EQ.2.AND.IFLG.GT.IFLH.AND.IFLH.GT.IFLI)  THEN
         IF(IFLE.GT.IFLF.AND.IFLD.NE.IFLG) KSP=2+INT(0.75D0+S_RNDM(1))
         IF(IFLE.LT.IFLF.AND.IFLD.EQ.IFLG) KSP=3
         IF(IFLE.LT.IFLF.AND.IFLD.NE.IFLG) KSP=2+INT(0.25D0+S_RNDM(2))
      ENDIF
      KF=KFR(16*KSP-16+IFLG)+KFR(16*KSP-8+IFLH)+IFLI
      IF(KBAR(KF-40).eq.0)THEN
         WRITE(LUN,*)' jetset code missing,flvs:',kf,IFLG,IFLH,IFLI
         CALL SIB_REJECT('SIB_I4FLAV      ')
      ENDIF
      KF=KBAR(KF-40)
      IF(KF.le.14)THEN
         IF(PAR(106).gt.S_RNDM(3).and.IARNK.le.IPAR(61)) KF=KF-13+51
     &        +2*INT(PAR(108)+S_RNDM(4))
      ENDIF
      KF=ISIGN(KF,IFL1)
c     if leading baryon, mark quark to supress baryon production in the next iteration
c     i.e. forbid: Blead-Bbar-B combination
      if(iarnk.eq.1.and.IPAR(93).eq.1.and.iabs(mod(ifl1,100)).gt.10)then
         IFL2 = IFL2 + ISIGN(100,IFL2)
      endif

      IF(NDEBUG.gt.6)
     &     WRITE(LUN,*)' SIB_FLAV: output:',IFL1, IFL2_A, IRNK, IFL2, KF
      
      RETURN
      END
C=======================================================================

      SUBROUTINE SIB_ICFLAV( Q2, IS0, IS, IFL )

C-----------------------------------------------------------------------
C     Routine that samples symmetric between the available flavors
C     Input: Q2 - mass scale, usually filled with s_hat
C           IS0 - input flavor sign
C     
C     Output: IFL - flavor code: u,d,s,c or anti-quarks
C              IS - flavor sign: quark or anti-quark, if 0 passed then
C                   a new value is sampled      
C     Parameters:   kT_s and kT_c i.e. width of the fermi function
C-----------------------------------------------------------------------
C     f2py double precision,intent(in) :: q2
Cf2py integer,intent(in) :: is0
Cf2py integer,intent(out) :: is
Cf2py integer,intent(out) :: ifl
      IMPLICIT NONE
      DOUBLE PRECISION Q2
      INTEGER IFL,IS
      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      
      DOUBLE PRECISION XMS2,XMC2,P_S,P_C,S_RNDM,QMASS,FERMI
      INTEGER IFL1,IS0
      
      IF( NDEBUG.gt.6 )
     &     WRITE(LUN,*)'  SIB_ICFLAV: input (Q2,IFL,IS):',Q2,IFL,IS
      
c     quark or antiquark, sampled if input is zero
      IF(IS0.eq.0) THEN
         IS = -1 + 2*INT((2.D0-EPS8)*S_RNDM(IS0))
      ELSE
         IS = IS0
      ENDIF
      
c     strange and charm quark masses
      XMS2 = 4*QMASS(3)**2
      XMC2 = 4*QMASS(4)**2 * PAR(153)
      
c     strange and charm parameters
      IF(IPAR(89).eq.1)THEN
c     exponential thresholds
         P_S = PAR(154) *  EXP(-PAR(151)/Q2)
         P_C = PAR(156) * EXP(-PAR(152)/Q2)
      ELSE
c     fermi func. threshold      
c     P_s: 0 (u,d only) --> 1 (u,d,s equal) --> 2 (u+d,s+c equal)
c     P_c: 0 (s only) --> 1 (s,c) equal
         P_S = PAR(154) * FERMI( Q2, XMS2, -PAR(151) )
     &        + PAR(155) * FERMI( Q2, XMC2, -PAR(152) )
         P_C = PAR(156) * 0.5D0*FERMI( Q2, XMC2, -PAR(152) )
      ENDIF
      IF(NDEBUG.gt.6)THEN
         WRITE(LUN,*)'  SIB_ICFLAV: (4*M_S**2, P_S, kT):',
     &        xms2, P_s, PAR(151)
         WRITE(LUN,*)'  SIB_ICFLAV: (4*M_C**2, P_C, kT):',
     &        xmc2, P_c, PAR(152)
      ENDIF

c     sample u,d,s
      IFL1 = MIN(INT((2.D0+P_S)*S_RNDM(IS0))+1,3)

c     replace s with c
      IFL1 = IFL1 + IFL1/3*MIN(INT(P_C+S_RNDM(IS0)),1)

      IFL = IS*IFL1

      IF(NDEBUG.gt.6)
     &     WRITE(LUN,*)'  SIB_ICFLAV: output (Q2,IFL,IS):',Q2,IFL,IS

      END
C=======================================================================      
      
      SUBROUTINE SIB_DIFF (L0, JDIF1, Ecm, Irec, IREJ)

C-----------------------------------------------------------------------
C...diffraction dissociation
C.  INPUT L0 = index of "beam particle"
C.             the target is assumed to be a proton.
C.    JDIF1 = 1  "beam diffraction"
C.          = 2  "target diffraction"
C.          = 3  "double diffraction"
C     Irec  flag to avoid recursive calls of SIB_DIFF and SIB_NDIFF
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
C     The final particle output is contained in COMMON /S_PLIST/    
C     NP           : number of final particles
C     P(1:NP, 1:5) : 4-momenta + masses of the final particles 
C     LLIST (1:NP) : codes of final particles
      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP
      DOUBLE PRECISION XM2MIN,ALXMIN,SLOP0,ASLOP,BSLOP,XMASS
      COMMON /S_DIFMAss/ XM2MIN(6),ALXMIN(6),SLOP0,ASLOP,BSLOP,XMASS(2)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)


      INTEGER NFORIG,NPORIG,NIORIG,IPFLAG,IIFLAG,KINT
      COMMON /S_PARTO/ NFORIG(NP_max),NPORIG(NP_max),NIORIG(NP_max),
     &IPFLAG,IIFLAG,KINT

      INTEGER LRNK
      COMMON /SIB_RNK/ LRNK(8000)
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)

      INTEGER ICHP,ISTR,IBAR
      COMMON /S_CHP/ ICHP(99), ISTR(99), IBAR(99)

      INTEGER IISO,ISPN
      COMMON /S_SPN/ IISO(99), ISPN(99)

      INTEGER ICHM
      COMMON /S_CHM/ ICHM(99)

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      
      DIMENSION P0(5),P1(5),P2(5)
      
C     mapping array from particle space to diff. mass
c     six groups: proton, pion, kaons, hyperons,
c                 charmed mesons, charmed baryons
      INTEGER KK,I
      DIMENSION KK(99)      
      SAVE
      DATA (KK(I), I= 1,39) /5*0,3*2,4*3,2*1,6*0,6*2,3,6*2,6*4/      
      DATA (KK(I), I=40,99) /19*0,5,5,10*0,5,5,0,5,5,11*0,6,6,6,9*0,6/

      if(Ndebug.gt.1) 
     &  WRITE(LUN,*)' SIB_DIFF: called with (L0,JDIF1,Ecm):',
     &  L0,JDIF1,Ecm

      if(Irec.eq.1) THEN
         Ipflag= -1
         IIFLAG = 1
c     add incoming target particles
         PZ = PAWT(ECM,AM(IABS(L0)),AM(13))
         E2 = SQRT(PZ**2+AM2(13))
         CALL ADD_PRTN(0.D0,0.D0,-PZ,E2,AM(13),13,-2,0,IREFout)

c     add interactions
        xjdif = dble(jdif1)
        CALL ADD_PRTN(0.D0,0.D0,xjdif,ecm,0.D0,1,-1,IREFout,IREF)
      ENDIF
      CALL GET_NPP(NPP_0,NPP0_0)

      IDBAD = 0
      NTRY = 0
 20   IREJ = 1
      CALL INI_PRTN_STCK(NPP_0,NPP0_0)

      IF(NTRY.gt.20*Irec) RETURN ! zero tolerance for recursive calls 
      NTRY = NTRY + 1
     
      LL = L0
      LA = IABS(L0)
      XM2MAX = PAR(13)*Ecm*Ecm
      if(Ndebug.gt.1) 
     &   WRITE(LUN,*)' SIB_DIFF: max diff. mass (M,Xi):',XM2MAX,PAR(13)
      
C...Double diffraction
      IF (JDIF1 .EQ. 3)   THEN
         K = MAX(1,2-IBAR(LA)-ISTR(LA)-ICHM(LA))
         IF(Irec.eq.1) K = KK(LA)
c     minimal mass if larger than particle mass plus one pion
         XMMIN = XM2MIN(K)
         IF(Irec.eq.0) XMMIN = MAX(XMMIN,(AM(LA)+AM(7)+0.02D0)**2)
         XMB2 = XM2DIS(XMMIN,XM2MAX,1.D0)
         XMB = SQRT (XMB2)
         XMT2 = XM2DIS(XM2MIN(1),XM2MAX,1.D0)
         XMT = SQRT (XMT2)
         CALL TRANSFONSHELL(ECM,XMB,XMT,XM2MAX,0,P1,P2,IBAD)
         IF(IBAD.ne.0) goto 20
         XMASS(1) = XMB
         IF(Irec.eq.1)THEN
c     add diffractive system to parton stack
            CALL ADD_PRTN_4VEC(P1,3,0,0,Iref)
            CALL ADD_INT_REF(Iref,1)
            CALL ADD_PRTN_4VEC(P2,-3,0,0,Iref)
            CALL ADD_INT_REF(Iref,1)
         ENDIF
         if(Ndebug.gt.1) 
     &        write(lun,*)' double-diff.: (kb,xmb,kt,xmt)',LL,xmb,13,xmt
         CALL DIFDEC (LL, Irec, IDBAD, P1)
         IF(IDBAD.eq.1)goto 20
         Ipflag= -2
         XMASS(2) = XMT
         CALL DIFDEC (13, Irec, IDBAD, P2)
         IF(IDBAD.eq.1)goto 20
         IREJ = 0
         RETURN
      ENDIF

C...Single diffraction
      IF (JDIF1.EQ. 1)  THEN
         K = MAX(1,2-IBAR(LA))
         IF(Irec.eq.1) K = KK(LA)
         EM  = AM(13)
         EM2 = AM2(13)
         L = 13
         ZD = -1.D0
         if(Ndebug.gt.1) 
     &        write(lun,*)' single-diff. (beam): (kb)',LL
      ELSE
         K = 1
         EM  = AM(LA)
         EM2 = AM2(LA)
         L = LL
         LL = 13
         ZD = +1.D0
         if(Ndebug.gt.1) 
     &        write(lun,*)' single-diff. (target): (kt)', LL

      ENDIF
C...Generate the mass of the diffracted system Mx (1/Mx**2 distribution)
      XMMIN = XM2MIN(K)
      IF(Irec.eq.0) XMMIN = MAX(XMMIN,(AM(LA)+AM(7)+0.02D0)**2)
      XM2 = XM2DIS(XMMIN,XM2MAX,1.D0)
      ALX = log(XM2)
c... added part
      X = XM2/XM2MAX*PAR(13)
      IF (X.GT.PAR(13)-0.05D0) THEN
        PRO = 0.5D0*(1.D0+(X-PAR(13))/0.05D0)
        IF (S_RNDM(0).LT.PRO) X = 2.D0*PAR(13)-X
        XM2 = XM2MAX*X/PAR(13)
      ENDIF
c...

      XM = SQRT (XM2)
      XMB = XM
      XMT = XM
      XMASS(1) = XMB
      XMASS(2) = XMT

C..   kinematics
      CALL TRANSFONSHELL(ECM,XMB,EM,XM2MAX,0,P1,P2,IBAD)
      IF(IBAD.ne.0) goto 20

C...Generate the Kinematics of the pseudoelastic hadron
      NP = NP+1
      P(NP,4) = P2(4)
      P(NP,3) = abs(P2(3))*ZD
      P(NP,1) = p2(1)
      P(NP,2) = p2(2)
      P(NP,5) = EM
      LLIST(NP) = L
      NPORIG(NP) = IPFLAG
      niorig(NP) = iiflag
      LRNK(NP) = 0
      
C...Generating the hadronic system recoiling against the produced particle
      P0(5) = SQRT(XM2)
      P0(4) = P1(4)
      DO J=1,3
         P0(J) = -P(NP,J)
      ENDDO
      IF(Irec.eq.1)THEN
c     add diffractive system to parton stack
         CALL ADD_PRTN_4VEC(P1,JDIF1,0,0,Iref)
         CALL ADD_INT_REF(Iref,1)
         CALL ADD_PRTN_4VEC(P2,int(zd),0,0,Iref)
         CALL ADD_INT_REF(Iref,1)
      ENDIF
      CALL DIFDEC (LL, Irec, IDBAD, P0)
      IF(IDBAD.eq.1)goto 20
      IREJ = 0

      END
      
C=======================================================================
      
      SUBROUTINE DIFF_INI

C----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER I,NPION
      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      DOUBLE PRECISION XM2MIN,ALXMIN,SLOP0,ASLOP,BSLOP,XMASS
      COMMON /S_DIFMAss/ XM2MIN(6),ALXMIN(6),SLOP0,ASLOP,BSLOP,XMASS(2)
      INTEGER NIPAR_max,NPAR_max      
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)
      SAVE      
C...  Diffractive mass parameters from Sibyll 2.1
c     minimal mass
      DATA (XM2MIN(I), I=1,3) /1.5D0, 0.2D0, 0.6D0/ ! M_x**2(min) GeV**2
      DATA (ALXMIN(I), I=1,3)   ! log[M_x**2(min)]
     &     /0.405465D0,-1.6094379D0,-0.5108256D0/ 
C...  pt spectrum
      DATA SLOP0 /6.5D0/                 ! b (slope_ for Mx**2 > 5 GeV**2
      DATA ASLOP /31.10362D0/            ! fit to the slope parameter.
      DATA BSLOP /-15.29012D0/

C     minimal mass for strange and charmed hadrons: 
C     m_beam + n_pi * m_pi
      NPION = IPAR(86)
      
C     hyperons (4), lowest mass: lambda
      XM2MIN(4) = AM2(39) + NPION * AM2(7)
      ALXMIN(4) = log(XM2MIN(4))
      
C     charmed mesons (5), lowest mass: Dmeson
      XM2MIN(5) = AM2(59) + NPION * AM2(7)
      ALXMIN(5) = log(XM2MIN(5))
      
C     charmed baryons (6), lowest mass: lambda_c
      XM2MIN(6) = AM2(89) + NPION * AM2(7)
      ALXMIN(6) = log(XM2MIN(6))

c     debug output
      IF(NDEBUG.gt.1)THEN
         WRITE(LUN,*)'DIFF_INI: setting diff. mass parameters'
         WRITE(LUN,*)' min mass: ', (XM2MIN(I), I=1,6)
         WRITE(LUN,*)' log min mass: ', (ALXMIN(I), I=1,6)
      ENDIF
      
      END

C=======================================================================      

      DOUBLE PRECISION FUNCTION SIGELA_PN(plab)

C-----------------------------------------------------------------------
C
C     low-energy pn/np elastic cross section
C     (based on spline interpolations)
C
C                                              (R.Engel 02/01)
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      COMMON / SPAPPR / SECD1 , SECDN , VOFINT , IERR , NXY

      dimension PTPP(100),STPP(100),DERIV(100,2),Z(10),FV(10),FD(10,2)
      logical init
      SAVE
      data init /.true./

C  pn elastic cross section
      DATA (PTPP(K),K=    1,   18) /
     &  -1.0128D+00,-8.8365D-01,-7.8000D-01,-6.8973D-01,-5.7462D-01,
     &  -4.2138D-01,-2.9384D-01,-1.1581D-01, 1.1309D-01, 5.3273D-01,
     &   9.6497D-01, 1.4860D+00, 2.0449D+00, 2.6798D+00, 3.5939D+00,
     &   4.9903D+00, 6.2215D+00, 6.8942D+00/
      DATA (STPP(K),K=    1,   18) /
     &1.0001D+02,8.2414D+01,6.5819D+01,5.4660D+01,4.7794D+01,4.0500D+01,
     &3.5781D+01,3.3208D+01,2.9921D+01,2.3919D+01,1.8633D+01,1.4206D+01,
     &1.1068D+01,9.0752D+00,7.5167D+00,6.6817D+00,6.8455D+00,6.8568D+00/


C  initialize cross section tables

      if(init) then
        N = 18
        M = 0
        CALL SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,-1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)') 
     &      ' SIGELA_PN: spline initialization failed: ',IERR
          stop
        endif
        NXY_save = NXY
        init = .false.
      endif

C  spline interpolation 

      sigela_pn = 0.D0
      Z(1) = log(plab)

      if((Z(1).gt.PTPP(1)).and.(Z(1).lt.PTPP(N))) then
        M = 1
        NXY = NXY_save
        CALL SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)')
     &      ' SIGELA_PN: spline interpolation failed: ',IERR
          return
        endif
      else
        return
      endif

      sigela_pn = FV(1) 

      END

C=======================================================================

      DOUBLE PRECISION FUNCTION SIGELA_PP(plab)

C-----------------------------------------------------------------------
C
C     low-energy pp elastic cross section
C     (based on spline interpolations)
C
C                                              (R.Engel 02/01)
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      COMMON / SPAPPR / SECD1 , SECDN , VOFINT , IERR , NXY

      dimension PTPP(100),STPP(100),DERIV(100,2),Z(10),FV(10),FD(10,2)
      logical init
      SAVE
      data init /.true./

C  pp elastic cross section
      DATA (PTPP(K),K=    1,   20) /
     &  -1.0548D+00,-9.9070D-01,-8.2516D-01,-6.8608D-01,-4.7199D-01,
     &  -2.7085D-01,-1.0784D-01,-7.6152D-03, 1.6806D-01, 3.3154D-01,
     &   5.4551D-01, 8.2275D-01, 1.3768D+00, 2.0058D+00, 2.9862D+00,
     &   3.7151D+00, 4.3182D+00, 5.1348D+00, 5.6750D+00, 6.2152D+00/
      DATA (STPP(K),K=    1,   20) /
     &4.2555D+01,3.7310D+01,2.8426D+01,2.4873D+01,2.2758D+01,2.2166D+01,
     &2.3350D+01,2.4450D+01,2.5212D+01,2.4535D+01,2.2927D+01,1.9459D+01,
     &1.4213D+01,1.0745D+01,8.4602D+00,7.3604D+00,6.8528D+00,6.6836D+00,
     &6.6836D+00,6.6836D+00/


C  initialize cross section tables

      if(init) then
        N = 20
        M = 0
        CALL SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,-1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)') 
     &      ' SIGELA_PP: spline initialization failed: ',IERR
          stop
        endif
        NXY_save = NXY
        init = .false.
      endif

C  spline interpolation 

      sigela_pp = 0.D0
      Z(1) = log(plab)

      if((Z(1).gt.PTPP(1)).and.(Z(1).lt.PTPP(N))) then
        M = 1
        NXY = NXY_save
        CALL SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)')
     &      ' SIGELA_PP: spline interpolation failed: ',IERR
          return
        endif
      else
        return
      endif

      sigela_pp = FV(1) 

      END

C=======================================================================

      DOUBLE PRECISION FUNCTION SIGTOT_PN(plab)

C-----------------------------------------------------------------------
C
C     low-energy pn and np total cross section
C     (based on spline interpolations)
C
C                                              (R.Engel 02/01)
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      COMMON / SPAPPR / SECD1 , SECDN , VOFINT , IERR , NXY

      dimension PTPP(100),STPP(100),DERIV(100,2),Z(10),FV(10),FD(10,2)
      logical init
      SAVE
      data init /.true./

C  pn total cross section
      DATA (PTPP(K),K=    1,   17) /
     &  -1.0129D+00,-8.4520D-01,-7.4136D-01,-5.3626D-01,-3.3210D-01,
     &  -1.2859D-01, 8.7237D-02, 3.1519D-01, 6.7022D-01, 1.0889D+00,
     &   1.5714D+00, 2.0792D+00, 2.6760D+00, 3.9453D+00, 4.9226D+00,
     &   5.6207D+00, 6.7629D+00/
      DATA (STPP(K),K=    1,   17) /
     &1.0000D+02,7.9053D+01,6.0976D+01,4.5194D+01,3.6729D+01,3.3429D+01,
     &3.3142D+01,3.7303D+01,4.0316D+01,4.1607D+01,4.0746D+01,3.9885D+01,
     &3.8594D+01,3.8307D+01,3.8881D+01,3.9168D+01,4.1320D+01/


C  initialize cross section tables

      if(init) then
        N = 17
        M = 0
        CALL SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,-1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)') 
     &      ' SIGTOT_PN: spline initialization failed: ',IERR
          stop
        endif
        NXY_save = NXY
        init = .false.
      endif

C  spline interpolation 

      sigtot_pn = 0.D0
      Z(1) = log(plab)

      if((Z(1).gt.PTPP(1)).and.(Z(1).lt.PTPP(N))) then
        M = 1
        NXY = NXY_save
        CALL SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)')
     &      ' SIGTOT_PN: spline interpolation failed: ',IERR
          return
        endif
      else
        return
      endif

      sigtot_pn = FV(1) 

      END

C=======================================================================

      DOUBLE PRECISION FUNCTION SIGTOT_PP(plab)

C-----------------------------------------------------------------------
C
C     low-energy pp 
C     (based on spline interpolations)
C
C                                              (R.Engel 02/01)
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      COMMON / SPAPPR / SECD1 , SECDN , VOFINT , IERR , NXY

      dimension PTPP(100),STPP(100),DERIV(100,2),Z(10),FV(10),FD(10,2)
      logical init
      SAVE
      data init /.true./

C  pp total cross section
      DATA (PTPP(K),K=    1,   23) /
     &  -1.4202D+00,-1.2583D+00,-1.0464D+00,-8.3253D-01,-6.0471D-01,
     &  -3.6376D-01,-8.4289D-02, 6.8739D-02, 1.9666D-01, 3.2471D-01,
     &   4.2673D-01, 5.5375D-01, 7.5675D-01, 1.0737D+00, 1.5176D+00,
     &   2.1393D+00, 2.7230D+00, 3.5353D+00, 4.3223D+00, 5.1728D+00,
     &   5.7949D+00, 6.2392D+00, 6.9122D+00/
      DATA (STPP(K),K=    1,   23) /
     &9.2081D+01,7.0000D+01,4.2437D+01,2.8579D+01,2.3858D+01,2.2335D+01,
     &2.3858D+01,2.8883D+01,3.5888D+01,4.3807D+01,4.7157D+01,4.7766D+01,
     &4.7157D+01,4.4569D+01,4.1523D+01,3.9695D+01,3.8782D+01,3.8173D+01,
     &3.8173D+01,3.8477D+01,3.9391D+01,4.0000D+01,4.1523D+01/


C  initialize cross section tables

      if(init) then
        N = 23
        M = 0
        CALL SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,-1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)') 
     &      ' SIGTOT_PP: spline initialization failed: ',IERR
          stop
        endif
        NXY_save = NXY
        init = .false.
      endif

C  spline interpolation 

      sigtot_pp = 0.D0
      Z(1) = log(plab)

      if((Z(1).gt.PTPP(1)).and.(Z(1).lt.PTPP(N))) then
        M = 1
        NXY = NXY_save
        CALL SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)')
     &      ' SIGTOT_PP: spline interpolation failed: ',IERR
          return
        endif
      else
        return
      endif

      sigtot_pp = FV(1) 

      END

C=======================================================================

      DOUBLE PRECISION FUNCTION SIGELA_PIPP(plab)

C-----------------------------------------------------------------------
C
C     low-energy pi+p elastic cross section
C     (based on spline interpolations)
C
C                                              (R.Engel 05/01)
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      COMMON / SPAPPR / SECD1 , SECDN , VOFINT , IERR , NXY

      dimension PTPP(100),STPP(100),DERIV(100,2),Z(10),FV(10),FD(10,2)
      logical init
      SAVE
      data init /.true./

C  pi+p elastic cross section
      DATA (PTPP(K),K=    1,   24) /
     &  -9.1117D-01,-8.4887D-01,-7.8656D-01,-6.6196D-01,-5.3736D-01,
     &  -4.4390D-01,-3.6083D-01,-2.6738D-01,-1.8431D-01,-5.9706D-02,
     &   5.4515D-02, 1.3758D-01, 2.4142D-01, 3.5564D-01, 4.0756D-01,
     &   5.1140D-01, 6.9830D-01, 1.0410D+00, 1.6225D+00, 2.2455D+00,
     &   2.9620D+00, 3.7407D+00, 4.6026D+00, 5.5163D+00/
      DATA (STPP(K),K=    1,   24) /
     &7.3812D+01,5.8453D+01,4.5967D+01,3.1602D+01,2.2652D+01,1.6133D+01,
     &1.2044D+01,9.2818D+00,8.3978D+00,9.9448D+00,1.2818D+01,1.4144D+01,
     &1.6354D+01,1.8011D+01,1.7238D+01,1.2928D+01,1.0055D+01,7.1823D+00,
     &5.5249D+00,4.6409D+00,3.6464D+00,2.9834D+00,3.2044D+00,3.0939D+00/


C  initialize cross section tables

      if(init) then
        N = 24
        M = 0
        CALL SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,-1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)') 
     &      ' SIGELA_PIPP: spline initialization failed: ',IERR
          stop
        endif
        NXY_save = NXY
        init = .false.
      endif

C  spline interpolation 

      sigela_pipp = 0.D0
      Z(1) = log(plab)

      if((Z(1).gt.PTPP(1)).and.(Z(1).lt.PTPP(N))) then
        M = 1
        NXY = NXY_save
        CALL SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)')
     &      ' SIGELA_PIPP: spline interpolation failed: ',IERR
          return
        endif
      else
        return
      endif

      sigela_pipp = FV(1) 

      END

C=======================================================================

      DOUBLE PRECISION FUNCTION SIGTOT_PIPP(plab)

C-----------------------------------------------------------------------
C
C     low-energy pi+p total cross section
C     (based on spline interpolations)
C
C                                              (R.Engel 05/01)
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      COMMON / SPAPPR / SECD1 , SECDN , VOFINT , IERR , NXY

      dimension PTPP(100),STPP(100),DERIV(100,2),Z(10),FV(10),FD(10,2)
      logical init
      SAVE
      data init /.true./

C  pi+p total cross section
      DATA (PTPP(K),K=    1,   37) /
     &  -9.2155D-01,-8.6963D-01,-8.0733D-01,-7.2426D-01,-5.4774D-01,
     &  -4.7505D-01,-4.1275D-01,-3.6083D-01,-3.0891D-01,-2.2585D-01,
     &  -1.7393D-01,-8.0473D-02, 2.3363D-02, 1.5835D-01, 2.3104D-01,
     &   2.9334D-01, 3.1411D-01, 3.5564D-01, 4.1794D-01, 4.2833D-01,
     &   4.9063D-01, 5.7370D-01, 6.7754D-01, 7.2945D-01, 8.1252D-01,
     &   8.8521D-01, 9.9943D-01, 1.1033D+00, 1.4044D+00, 1.7782D+00,
     &   2.1313D+00, 2.6712D+00, 3.2942D+00, 3.8342D+00, 4.6441D+00,
     &   5.4748D+00, 5.8382D+00/
      DATA (STPP(K),K=    1,   37) /
     &7.3812D+01,6.4420D+01,5.0939D+01,3.7790D+01,2.3867D+01,1.8674D+01,
     &1.6022D+01,1.5138D+01,1.4365D+01,1.5138D+01,1.7127D+01,2.0773D+01,
     &2.4420D+01,2.7845D+01,3.3591D+01,3.9116D+01,4.0773D+01,4.1215D+01,
     &4.0000D+01,3.8232D+01,3.3370D+01,3.0608D+01,2.9061D+01,2.8619D+01,
     &2.9834D+01,3.0829D+01,3.0497D+01,2.9061D+01,2.7514D+01,2.5746D+01,
     &2.4862D+01,2.3646D+01,2.3094D+01,2.2873D+01,2.3204D+01,2.3978D+01,
     &2.4420D+01/


C  initialize cross section tables

      if(init) then
        N = 37
        M = 0
        CALL SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,-1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)') 
     &      ' SIGTOT_PIPP: spline initialization failed: ',IERR
          stop
        endif
        NXY_save = NXY
        init = .false.
      endif

C  spline interpolation 

      sigtot_pipp = 0.D0
      Z(1) = log(plab)

      if((Z(1).gt.PTPP(1)).and.(Z(1).lt.PTPP(N))) then
        M = 1
        NXY = NXY_save
        CALL SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)')
     &      ' SIGTOT_PIPP: spline interpolation failed: ',IERR
          return
        endif
      else
        return
      endif

      sigtot_pipp = FV(1) 

      END

C=======================================================================

      DOUBLE PRECISION FUNCTION SIGELA_PIMP(plab)

C-----------------------------------------------------------------------
C
C     low-energy pi-p elastic cross section
C     (based on spline interpolations)
C
C                                              (R.Engel 05/01)
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      COMMON / SPAPPR / SECD1 , SECDN , VOFINT , IERR , NXY

      dimension PTPP(100),STPP(100),DERIV(100,2),Z(10),FV(10),FD(10,2)
      logical init
      SAVE
      data init /.true./

C  pi-p elastic cross section
      DATA (PTPP(K),K=    1,   56) /
     &  -1.8980D+00,-1.5458D+00,-1.4323D+00,-1.3602D+00,-1.2880D+00,
     &  -1.2571D+00,-1.1845D+00,-1.1531D+00,-1.1112D+00,-1.0691D+00,
     &  -1.0063D+00,-9.1252D-01,-8.2935D-01,-7.0477D-01,-6.0118D-01,
     &  -4.6652D-01,-4.1489D-01,-3.9435D-01,-3.6334D-01,-3.4267D-01,
     &  -3.0100D-01,-2.6966D-01,-2.4866D-01,-2.1741D-01,-1.6542D-01,
     &  -1.1357D-01,-9.2992D-02,-8.2923D-02,-4.1875D-02,-1.1054D-02,
     &   3.0281D-02, 7.2145D-02, 8.2958D-02, 1.1458D-01, 1.5645D-01,
     &   2.6051D-01, 3.4368D-01, 3.8539D-01, 4.7900D-01, 5.3080D-01,
     &   6.3455D-01, 7.4898D-01, 9.1527D-01, 1.1023D+00, 1.3412D+00,
     &   1.5594D+00, 1.9541D+00, 2.4007D+00, 2.7122D+00, 3.0653D+00,
     &   3.4392D+00, 3.8130D+00, 4.2387D+00, 5.0175D+00, 5.3602D+00,
     &   5.8897D+00/
      DATA (STPP(K),K=    1,   56) /
     &2.9793D+00,9.7103D+00,1.5007D+01,1.9862D+01,2.3393D+01,2.5269D+01,
     &2.6041D+01,2.4276D+01,2.1076D+01,1.6772D+01,1.3021D+01,1.0372D+01,
     &9.6000D+00,9.8207D+00,1.1697D+01,1.4234D+01,1.6441D+01,1.8207D+01,
     &1.9310D+01,2.0083D+01,1.8979D+01,1.7545D+01,1.5779D+01,1.5007D+01,
     &1.4455D+01,1.5007D+01,1.6441D+01,1.8869D+01,2.2621D+01,2.5159D+01,
     &2.6703D+01,2.4166D+01,2.0855D+01,1.7214D+01,1.4676D+01,1.2910D+01,
     &1.2138D+01,1.0814D+01,9.6000D+00,1.0483D+01,1.1145D+01,9.6000D+00,
     &8.3862D+00,7.5034D+00,6.6207D+00,6.0690D+00,4.9655D+00,4.4138D+00,
     &4.4138D+00,3.7517D+00,3.3103D+00,3.2000D+00,3.3103D+00,3.3103D+00,
     &3.3103D+00,3.5310D+00/


C  initialize cross section tables

      if(init) then
        N = 56
        M = 0
        CALL SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,-1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)') 
     &      ' SIGELA_PIMP: spline initialization failed: ',IERR
          stop
        endif
        NXY_save = NXY
        init = .false.
      endif

C  spline interpolation 

      sigela_pimp = 0.D0
      Z(1) = log(plab)

      if((Z(1).gt.PTPP(1)).and.(Z(1).lt.PTPP(N))) then
        M = 1
        NXY = NXY_save
        CALL SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)')
     &      ' SIGELA_PIMP: spline interpolation failed: ',IERR
          return
        endif
      else
        return
      endif

      sigela_pimp = FV(1) 

      END

C=======================================================================

      DOUBLE PRECISION FUNCTION SIGTOT_PIMP(plab)

C-----------------------------------------------------------------------
C
C     low-energy pi-p total cross section
C     (based on spline interpolations)
C
C                                              (R.Engel 05/01)
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      COMMON / SPAPPR / SECD1 , SECDN , VOFINT , IERR , NXY

      dimension PTPP(100),STPP(100),DERIV(100,2),Z(10),FV(10),FD(10,2)
      logical init
      SAVE
      data init /.true./

C  pi-p total cross section
      DATA (PTPP(K),K=    1,   53) /
     &  -1.9302D+00,-1.8269D+00,-1.6617D+00,-1.5490D+00,-1.4577D+00,
     &  -1.3146D+00,-1.2630D+00,-1.2211D+00,-1.1686D+00,-1.1364D+00,
     &  -1.0937D+00,-1.0305D+00,-9.4645D-01,-8.5245D-01,-7.6915D-01,
     &  -6.7584D-01,-5.2057D-01,-4.3813D-01,-4.0781D-01,-3.6669D-01,
     &  -3.1507D-01,-2.8372D-01,-2.6240D-01,-2.0995D-01,-1.7861D-01,
     &  -1.1661D-01,-9.6329D-02,-7.6149D-02,-3.5817D-02,-5.0811D-03,
     &   1.5958D-02, 5.8095D-02, 1.1175D-01, 1.7444D-01, 1.9540D-01,
     &   2.8868D-01, 3.7173D-01, 4.5500D-01, 5.4845D-01, 6.4176D-01,
     &   7.1436D-01, 8.3919D-01, 9.6397D-01, 1.3069D+00, 1.7018D+00,
     &   2.0447D+00, 2.5952D+00, 3.1249D+00, 3.6130D+00, 4.1426D+00,
     &   4.8175D+00, 5.3159D+00, 5.9284D+00/
      DATA (STPP(K),K=    1,   53) /
     &1.1145D+01,1.5007D+01,2.2179D+01,3.4428D+01,5.0428D+01,6.7862D+01,
     &7.0952D+01,6.7972D+01,6.3007D+01,5.5393D+01,4.6566D+01,3.9614D+01,
     &3.1779D+01,2.7586D+01,2.5821D+01,2.6924D+01,3.0676D+01,3.5531D+01,
     &4.1931D+01,4.5131D+01,4.7448D+01,4.5903D+01,4.1600D+01,3.7517D+01,
     &3.6083D+01,3.8400D+01,4.2152D+01,4.6676D+01,5.5945D+01,5.9145D+01,
     &5.7048D+01,5.2414D+01,3.9062D+01,3.6083D+01,3.4538D+01,3.5862D+01,
     &3.6083D+01,3.4538D+01,3.4538D+01,3.5641D+01,3.6303D+01,3.4538D+01,
     &3.3214D+01,3.1117D+01,2.8690D+01,2.7145D+01,2.5600D+01,2.4717D+01,
     &2.4166D+01,2.4166D+01,2.3945D+01,2.4055D+01,2.5159D+01/


C  initialize cross section tables

      if(init) then
        N = 53
        M = 0
        CALL SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,-1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)') 
     &      ' SIGTOT_PIMP: spline initialization failed: ',IERR
          stop
        endif
        NXY_save = NXY
        init = .false.
      endif

C  spline interpolation 

      sigtot_pimp = 0.D0
      Z(1) = log(plab)

      if((Z(1).gt.PTPP(1)).and.(Z(1).lt.PTPP(N))) then
        M = 1
        NXY = NXY_save
        CALL SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)')
     &      ' SIGTOT_PIMP: spline interpolation failed: ',IERR
          return
        endif
      else
        return
      endif

      sigtot_pimp = FV(1) 

      END

C=======================================================================

      DOUBLE PRECISION FUNCTION SIGELA_KPP(plab)

C-----------------------------------------------------------------------
C
C     low-energy K+p elastic cross section
C     (based on spline interpolations)
C
C                                              (R.Engel 05/01)
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      COMMON / SPAPPR / SECD1 , SECDN , VOFINT , IERR , NXY

      dimension PTPP(100),STPP(100),DERIV(100,2),Z(10),FV(10),FD(10,2)
      logical init
      SAVE
      data init /.true./

C  K+p elastic cross section
      DATA (PTPP(K),K=    1,   22) /
     &  -1.1500D+00,-8.0733D-01,-5.4774D-01,-4.1275D-01,-2.5700D-01,
     &  -8.0474D-02, 7.5281D-02, 2.5180D-01, 3.7641D-01, 5.3216D-01,
     &   6.8792D-01, 8.4368D-01, 1.0929D+00, 1.5913D+00, 1.9340D+00,
     &   2.3182D+00, 2.8166D+00, 3.2215D+00, 3.4708D+00, 3.9276D+00,
     &   4.6233D+00, 5.5475D+00/
      DATA (STPP(K),K=    1,   22) /
     &1.2227D+01,1.2570D+01,1.2499D+01,1.2498D+01,1.2428D+01,1.2012D+01,
     &1.1183D+01,1.0284D+01,9.4544D+00,8.2796D+00,6.8977D+00,5.9300D+00,
     &4.6854D+00,3.6461D+00,3.2293D+00,3.0193D+00,2.6704D+00,2.4602D+00,
     &2.3203D+00,2.0407D+00,2.2426D+00,2.5809D+00/


C  initialize cross section tables

      if(init) then
        N = 22
        M = 0
        CALL SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,-1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)') 
     &      ' SIGELA_KPP: spline initialization failed: ',IERR
          stop
        endif
        NXY_save = NXY
        init = .false.
      endif

C  spline interpolation 

      sigela_kpp = 0.D0
      Z(1) = log(plab)

      if((Z(1).gt.PTPP(1)).and.(Z(1).lt.PTPP(N))) then
        M = 1
        NXY = NXY_save
        CALL SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)')
     &      ' SIGELA_KPP: spline interpolation failed: ',IERR
          return
        endif
      else
        return
      endif

      sigela_kpp = FV(1) 

      END

C=======================================================================

      DOUBLE PRECISION FUNCTION SIGTOT_KPP(plab)

C-----------------------------------------------------------------------
C
C     low-energy K+p total cross section
C     (based on spline interpolations)
C
C                                              (R.Engel 05/01)
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      COMMON / SPAPPR / SECD1 , SECDN , VOFINT , IERR , NXY

      dimension PTPP(100),STPP(100),DERIV(100,2),Z(10),FV(10),FD(10,2)
      logical init
      SAVE
      data init /.true./

C  K+p total cross section
      DATA (PTPP(K),K=    1,   20) /
     &  -1.0981D+00,-7.1388D-01,-4.7505D-01,-3.1930D-01,-1.7393D-01,
     &  -8.0474D-02, 2.3363D-02, 9.6049D-02, 1.9989D-01, 3.2449D-01,
     &   4.6986D-01, 6.2562D-01, 8.3329D-01, 1.0825D+00, 1.4355D+00,
     &   2.1001D+00, 2.6920D+00, 3.5434D+00, 4.6337D+00, 5.7448D+00/
      DATA (STPP(K),K=    1,   20) /
     &1.2158D+01,1.2362D+01,1.2429D+01,1.2428D+01,1.3187D+01,1.4429D+01,
     &1.5809D+01,1.7327D+01,1.8224D+01,1.8430D+01,1.7945D+01,1.7806D+01,
     &1.7459D+01,1.7250D+01,1.7041D+01,1.7381D+01,1.7446D+01,1.7853D+01,
     &1.8881D+01,2.0529D+01/


C  initialize cross section tables

      if(init) then
        N = 20
        M = 0
        CALL SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,-1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)') 
     &      ' SIGTOT_KPP: spline initialization failed: ',IERR
          stop
        endif
        NXY_save = NXY
        init = .false.
      endif

C  spline interpolation 

      sigtot_kpp = 0.D0
      Z(1) = log(plab)

      if((Z(1).gt.PTPP(1)).and.(Z(1).lt.PTPP(N))) then
        M = 1
        NXY = NXY_save
        CALL SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)')
     &      ' SIGTOT_KPP: spline interpolation failed: ',IERR
          return
        endif
      else
        return
      endif

      sigtot_kpp = FV(1) 

      END


C=======================================================================

      DOUBLE PRECISION FUNCTION SIGELA_KMP(plab)

C-----------------------------------------------------------------------
C
C     low-energy K-p elastic cross section
C     (based on spline interpolations)
C
C                                              (R.Engel 05/01)
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      COMMON / SPAPPR / SECD1 , SECDN , VOFINT , IERR , NXY

      dimension PTPP(100),STPP(100),DERIV(100,2),Z(10),FV(10),FD(10,2)
      logical init
      SAVE
      data init /.true./

C  K-p elastic cross section
      DATA (PTPP(K),K=    1,   36) /
     &  -1.7871D+00,-1.4709D+00,-1.2813D+00,-1.1867D+00,-1.0179D+00,
     &  -8.8055D-01,-8.0666D-01,-7.9648D-01,-7.7560D-01,-6.5951D-01,
     &  -5.6450D-01,-4.7995D-01,-3.9539D-01,-3.4256D-01,-2.7894D-01,
     &  -2.4691D-01,-2.0439D-01,-1.1952D-01,-1.3598D-02, 6.0479D-02,
     &   1.1311D-01, 1.4462D-01, 2.0784D-01, 2.6053D-01, 3.2387D-01,
     &   4.4022D-01, 5.5672D-01, 6.9424D-01, 8.6348D-01, 1.2127D+00,
     &   1.6678D+00, 2.3770D+00, 3.2133D+00, 3.9226D+00, 4.6425D+00,
     &   5.1612D+00/
      DATA (STPP(K),K=    1,   36) /
     &6.8962D+01,5.6135D+01,4.7307D+01,4.0271D+01,3.5582D+01,3.2549D+01,
     &3.0480D+01,2.6617D+01,2.3858D+01,2.0410D+01,1.7927D+01,1.6549D+01,
     &1.5308D+01,1.4343D+01,1.5310D+01,1.7794D+01,1.9451D+01,2.1108D+01,
     &2.1661D+01,2.1386D+01,1.8490D+01,1.6144D+01,1.3386D+01,1.1041D+01,
     &9.3860D+00,8.4219D+00,8.8376D+00,7.8738D+00,6.4965D+00,4.7080D+00,
     &3.8869D+00,3.3456D+00,2.6682D+00,2.5409D+00,2.6896D+00,2.6974D+00/


C  initialize cross section tables

      if(init) then
        N = 36
        M = 0
        CALL SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,-1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)') 
     &      ' SIGELA_KMP: spline initialization failed: ',IERR
          stop
        endif
        NXY_save = NXY
        init = .false.
      endif

C  spline interpolation 

      sigela_kmp = 0.D0
      Z(1) = log(plab)

      if((Z(1).gt.PTPP(1)).and.(Z(1).lt.PTPP(N))) then
        M = 1
        NXY = NXY_save
        CALL SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)')
     &      ' SIGELA_KMP: spline interpolation failed: ',IERR
          return
        endif
      else
        return
      endif

      sigela_kmp = FV(1) 

      END

C=======================================================================

      DOUBLE PRECISION FUNCTION SIGTOT_KMP(plab)

C-----------------------------------------------------------------------
C
C     low-energy K-p total cross section
C     (based on spline interpolations)
C
C                                              (R.Engel 05/01)
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      COMMON / SPAPPR / SECD1 , SECDN , VOFINT , IERR , NXY

      dimension PTPP(100),STPP(100),DERIV(100,2),Z(10),FV(10),FD(10,2)
      logical init
      SAVE
      data init /.true./

C  K-p total cross section
      DATA (PTPP(K),K=    1,   43) /
     &  -1.3500D+00,-1.2345D+00,-9.8216D-01,-8.2491D-01,-7.4143D-01,
     &  -6.1508D-01,-4.5679D-01,-3.7223D-01,-2.9802D-01,-2.6595D-01,
     &  -1.7037D-01,-1.0660D-01,-2.1599D-02,-2.5037D-04, 6.3445D-02,
     &   8.4428D-02, 1.3703D-01, 1.5769D-01, 1.8898D-01, 2.4156D-01,
     &   3.3667D-01, 3.5796D-01, 4.1106D-01, 5.1700D-01, 5.9099D-01,
     &   6.5431D-01, 6.9651D-01, 7.7067D-01, 8.5538D-01, 9.6104D-01,
     &   1.1303D+00, 1.3209D+00, 1.4266D+00, 1.5853D+00, 1.8075D+00,
     &   1.9769D+00, 2.4743D+00, 3.0353D+00, 3.5222D+00, 4.0515D+00,
     &   4.6550D+00, 5.1949D+00, 5.7455D+00/
      DATA (STPP(K),K=    1,   43) /
     &9.7669D+01,8.8840D+01,7.2700D+01,5.8076D+01,4.6625D+01,4.0142D+01,
     &3.5315D+01,3.4074D+01,3.5041D+01,3.7939D+01,4.0838D+01,4.3185D+01,
     &4.6084D+01,4.7740D+01,4.9397D+01,4.7603D+01,4.4430D+01,3.9601D+01,
     &3.5186D+01,3.1876D+01,3.0221D+01,3.1325D+01,3.2982D+01,3.3674D+01,
     &3.2571D+01,3.0640D+01,2.9261D+01,2.9814D+01,2.9953D+01,2.8023D+01,
     &2.6922D+01,2.6924D+01,2.5684D+01,2.4859D+01,2.4034D+01,2.3761D+01,
     &2.2112D+01,2.1155D+01,2.0472D+01,2.0480D+01,2.0627D+01,2.0773D+01,
     &2.1472D+01/


C  initialize cross section tables

      if(init) then
        N = 43
        M = 0
        CALL SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,-1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)') 
     &      ' SIGTOT_KMP: spline initialization failed: ',IERR
          stop
        endif
        NXY_save = NXY
        init = .false.
      endif

C  spline interpolation 

      sigtot_kmp = 0.D0
      Z(1) = log(plab)

      if((Z(1).gt.PTPP(1)).and.(Z(1).lt.PTPP(N))) then
        M = 1
        NXY = NXY_save
        CALL SPLIN3(PTPP,STPP,DERIV,N,100,Z,FV,FD,M,10,1)
        if(IERR.ne.0) then
          write(6,'(1x,a,i6)')
     &      ' SIGTOT_KMP: spline interpolation failed: ',IERR
          return
        endif
      else
        return
      endif

      sigtot_kmp = FV(1) 

      END


C=======================================================================

      SUBROUTINE SPLIN3(X,Y,DERIV,N,NC,Z,FVALUE,FDERIV,M,MC,IOP)

C-----------------------------------------------------------------------
C
C     CERN LIBRARY PROGRAM NO E-209.
C
C     REVISED VERSION JULY 1973.
C
C     CHANGED BY R.ENGEL (10/10/93) TO CONFORM WITH F77 STANDARD
C
C     PURPOSE = TO COMPUTE A NATURAL SPLINE APPROXIMATION OF THIRD ORDER
C               FOR A FUNCTION Y(X) GIVEN IN THE N POINTS (X(I),Y(I)) ,
C               I=1(1)N.
C
C     PARAMETERS (IN LIST).
C
C     X       = AN ARRAY STORING THE INPUT ARGUMENTS.DIMENSION X(N).
C     Y       = AN ARRAY STORING THE INPUT FUNCTION VALUES.THE ELEMENT
C               Y(I) REPRESENT THE FUNCTION VALUE Y(X) FOR X=X(I).
C     DERIV   = AN ARRAY USED FOR STORING THE COMPUTED DERIVATIVES OF
C               THE FUNCTION Y(X).IN DERIV(I,1) AND DERIV(I,2) ARE STOR-
C               ED THE FIRST-AND SECOND ORDER DERIVATIVES OF Y(X) FOR
C               X=X(I) RESPECTIVELY.
C     N       = NUMBER OF INPUT FUNCTION VALUES.
C     NC      = ARRAY DERIV IS DIMENSIONED DERIV(NC,2) IN CALLING
C               PROGRAM.
C     Z       = AN ARRAY STORING THE ARGUMENTS FOR THE INTERPOLATED
C               VALUES TO BE COMPUTED.
C     FVALUE  = AN ARRAY STORING THE COMPUTED INTERPOLATED VALUES.
C               FVALUE(J) REPRESENT THE FUNCTION VALUE FVALUE(Z) FOR
C               Z=Z(J).
C     FDERIV    = AN ARRAY USED FOR STORING THE DERIVATIVES OF THE COM-
C               PUTED INTERPOLATED VALUES.EXPLANATION AS FOR DERIV.
C     M       = NUMBER OF INTERPOLATED VALUES TO BE COMPUTED.
C     MC      = ARRAY FDERIV IS DIMENSIONED FDERIV(MC,2) IN CALLING
C               PROGRAM.
C     IOP     = OPTION PARAMETER.FOR IOP.LE.0 THE DERIVATIVES FOR EACH
C               SUB-INTERVAL IN THE SPLINE APPROXIMATION ARE COMPUTED.
C                                  IOP=-1, THE SECOND ORDER END-POINT
C                                          DERIVATIVES ARE COMPUTED BY
C                                          LINEAR EXTRAPOLATION.
C                                  IOP=0 , THE SECOND ORDER END-POINT
C                                          DERIVATIVES ASSUMED TO BE GI-
C                                          VEN (SEE COMMON /SPAPPR/).
C                                  IOP=1 , COMPUTE SPLINE APPROXIMATIONS
C                                          FOR THE ARGUMENTS GIVEN IN
C                                          THE ARRAY Z,THE DERIVATIVES
C                                          BEEING ASSUMED TO HAVE BEEN
C                                          CALCULATED IN A PREVIOUS CALL
C                                          ON THE ROUTINE.
C
C     PARAMETERS (IN COMMON BLOCK / SPAPPR /).
C
C     SECD1   = VALUE OF THE SECOND DERIVATIVE D2Y(X)/DX2 FOR THE INPUT
C               ARGUMENT X=X(1).
C     SECDN   = VALUE OF THE SECOND DERIVATIVE D2Y(X)/DX2 FOR THE INPUT
C               ARGUMENT X=X(N).
C               NB. VALUES HAVE TO BE ASSIGNED TO SECD1 AND SECDN IN THE
C               CALLING PROGRAM.IF A NATURAL SPLINE FIT IS WANTED PUT
C               SECD1=SECDN=0.
C     VOFINT  = COMPUTED APPROXIMATION FOR THE INTEGRAL OF Y(X) TAKEN
C               FROM X(1) TO X(N).
C     IERR    = ERROR PARAMETER.IERR=0,NO ERRORS OCCURED.
C                               IERR=1,THE NUMBER OF POINTS TOO SMALL
C                                      I.E.N LESS THAN 4.
C                               IERR=2,THE ARGUMENTS X(I) NOT IN INCREA-
C                                      SING ORDER.
C                               IERR=3,ARGUMENT TO BE USED IN INTERPOLA-
C                                      TION ABOVE RANGE.
C                               IERR=4,ARGUMENT TO BE USED IN INTERPOLA-
C                                      TION BELOW RANGE.
C     NXY     = N (SEE ABOVE),HAS TO BE STORED FOR ENTRIES CORRESPONDING
C               TO IOP=1.
C
C**********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
C
      DIMENSION X(NC) , Y(NC) , DERIV(NC,2) , Z(MC) , FVALUE(MC) ,
     1          FDERIV(MC,2)
C
      COMMON / SPAPPR / SECD1 , SECDN , VOFINT , IERR , NXY
      SAVE
      DATA THIRD , SIXTH / .333333333333333D0 , .166666666666667D0 /
C
C 1000 
      IF (IOP.GT.0) GO TO 1110
C
      IERR=0
C
C     CHECK IF ENOUGH DATA-POINTS ARE AVAILABLEI.E. IF N LESS THAN 4 NO
C     THIRD ORDER SPLINE APPROXIMATION IS POSSIBLE.
C
      IF (N.GE.4) GO TO 1010
C
      IERR=1
      GO TO 2000
C
C     START CALCULATION OF COEFFICIENTS TO BE USED IN THE SYSTEM OF EQU-
C     ATIONS FOR THE SECOND ORDER DERIVATIVES OF Y(X).
C
 1010 IF (IOP.NE.-1) GO TO 1015
      SECD1=0.D0
      SECDN = 0.D0
      BET1=1.D0/(1.D0+0.5D0*(X(2)-X(1))/(X(3)-X(2)))
      ALF1=BET1*(1.D0- ((X(2)-X(1))/(X(3)-X(2)))**2)
      BETN=1.D0/(1.D0+0.5D0*(X(N)-X(N-1))/(X(N-1)-X(N-2)))
      ALFN=BETN*(1.D0- ((X(N)-X(N-1))/(X(N-1)-X(N-2)))**2)
C
 1015 DERIV(1,2)=SECD1
      DERIV(N,2)=SECDN
      DERIV(1,1)=0.D0
      DXPLUS=X(2)-X(1)
C
C     CHECK IF ARGUMENTS ARE IN INCREASING ORDER.IF NOT PRINT ERROR
C     MESSAGE AND STOP.
C
      IF ( DXPLUS.GT.0.D0) GO TO 1020
      IN=1
      IERR=2
      GO TO 2000
C
 1020 DYPLUS=(Y(2)-Y(1))/DXPLUS
      IU=N-1
      DO 1040 I=2,IU
      DXMIN =DXPLUS
      DYMIN =DYPLUS
      DXPLUS=X(I+1)-X(I)
C
C     CHECK IF ARGUMENTS ARE IN INCREASING ORDER.IF NOT PRINT ERROR
C     MESSAGE AND STOP.
C
      IF (DXPLUS.GT.0.D0) GO TO 1030
C
      IN=I
      IERR=2
      GO TO 2000
C
 1030 DXINV =1.D0/(DXPLUS+DXMIN)
      DYPLUS=(Y(I+1)-Y(I))/DXPLUS
      DIVDIF=DXINV*(DYPLUS-DYMIN)
      ALF   =0.5D0*DXINV*DXMIN
      BET   =0.5D0-ALF
C
      IF (I.EQ.2)  DIVDIF=DIVDIF-THIRD*ALF*DERIV(1,2)
      IF (I.EQ.IU) DIVDIF=DIVDIF-THIRD*BET*DERIV(N,2)
      IF (I.EQ.2) ALF=0.D0
C
      IF (IOP.NE.-1) GO TO 1035
      IF (I.NE.2) GO TO 1032
      BET=BET*ALF1
      DIVDIF=DIVDIF*BET1
      GO TO 1035
 1032 IF (I.NE.IU) GO TO 1035
      ALF=ALF*ALFN
      DIVDIF=DIVDIF*BETN
C
 1035 DXINV =1.D0/(1.D0+ALF*DERIV(I-1,1))
      DERIV(I,1)=-DXINV*BET
      DERIV(I,2)= DXINV*(3.D0*DIVDIF-ALF*DERIV(I-1,2))
 1040 CONTINUE
C
C     COMPUTE THE SECOND DERIVATIVES BY BACKWARDS RECURRENCE RELATION.
C     THE SECOND ORDER DERIVATIVES FOR X=X(N-1) ALREADY COMPUTED.
C
C 1050 
      DO 1060 I=2,IU
      J=N-I
      DERIV(J,2)=DERIV(J,1)*DERIV(J+1,2)+DERIV(J,2)
 1060 CONTINUE
C
      IF (IOP.NE.-1) GO TO 1070
      DERIV(1,2)=((X(3)-X(1))/(X(3)-X(2)))*DERIV(2,2)-((X(2)-X(1))/(X(3)
     $-X(2)))*DERIV(3,2)
      DERIV(N,2)=-((X(N)-X(N-1))/(X(N-1)-X(N-2)))*DERIV(N-2,2)+((X(N)-X(
     $N-2))/(X(N-1)-X(N-2)))*DERIV(N-1,2)
C
C     CALCULATION OF THE SECOND ORDER DERIVATIVES FINISHED.START CAL-
C     CULATION OF THE FIRST ORDER DERIVATIVES AND OF THE INTEGRAL.
C
 1070 VOFINT=0.D0
      DO 1080 I=1,IU
      DXPLUS=X(I+1)-X(I)
      DYPLUS=Y(I+1)-Y(I)
      DIVDIF=DYPLUS/DXPLUS
      DERIV(I,1)=DIVDIF-DXPLUS*(THIRD*DERIV(I,2)+SIXTH*DERIV(I+1,2))
      DXPLUS=0.5D0*DXPLUS
      VOFINT=VOFINT+DXPLUS*(Y(I+1)+Y(I)-THIRD*(DERIV(I+1,2)+DERIV(I,2))*
     $DXPLUS**2)
 1080 CONTINUE
C
C     COMPUTE THE LAST FIRST ORDER DERIVATIVE.
C
      DXPLUS=X(N)-X(N-1)
      DYPLUS=Y(N)-Y(N-1)
      DIVDIF=DYPLUS/DXPLUS
      DERIV(N,1)=DIVDIF+DXPLUS*(SIXTH*DERIV(N-1,2)+THIRD*DERIV(N,2))
C
C     CALCULATION OF FIRST ORDER DERIVATIVES AND INTEGRAL FINISHED.
C
C     SET VALUE OF N IN COMMON BLOCK / SPAPPR /.
C
      NXY=N
C
C     COMPUTE INTERPOLATED VALUES IF ANY.
C
 1110 IF (M.LT.1) RETURN
C
      XL=X(1)
      XU=X(2)
      IP=3
      IL=0
C
C 1120 
      DO 1160 J=1,M
      ARG=Z(J)
      IF (ARG.GT.XU) GO TO 1170
      IF (ARG.LT.XL) GO TO 1190
C
C     ARGUMENT IN CORRECT RANGE.CHECK IF POLYNOMIAL COEFFICIENTS HAVE
C     TO BE CALCULATED.
C
C 1130 
      IF (IL.GT.0) GO TO 1150
C
C     COMPUTE POLYNOMIAL COEFFICIENTS.
C
 1140 II=IP-2
      A0=Y(II)
      A1=DERIV(II,1)
      A4=DERIV(II,2)
      A6=(DERIV(II+1,2)-A4)/(XU-XL)
      A2=0.5D0*A4
      A3=SIXTH*A6
      A5=0.5D0*A6
      IL=1
C
C     CALCULATION OF POLYNOMIAL COEFFICIENTS FINISHED.COMPUTE VALUES.
C
 1150 ARG=ARG-XL
      FVALUE(J)=((A3*ARG+A2)*ARG+A1)*ARG+A0
      FDERIV(J,1)=(A5*ARG+A4)*ARG+A1
      FDERIV(J,2)=A6*ARG+A4
C
 1155 CONTINUE
      GOTO 1160
C
C     RANGE MOVING
C
C
C     ARGUMENT ABOVE PRESENT RANGE.SHIFT RANGE UPWARDS.
C
 1170 IF(IP.GT.NXY) GO TO 1185
      IPP=IP
      DO 1180 I=IPP,NXY
      IF (ARG.GT.X(I)) GO TO 1180
      XL=X(I-1)
      XU=X(I)
      IP=I+1
      IL=0
      GO TO 1140
C
 1180 CONTINUE
C
C     ARGUMENT  OUT OF RANGE,I.E. ARG GREATER THAN X(N).
C
 1185 IERR=3
      IP=NXY+1
      GO TO 2010
C
C     ARGUMENT BELOW PRESENT RANGE.SHIFT DOWNWARDS.
C
 1190 IPP=IP
      DO 1200 I=1,IPP
      II=IP-I-2
      IF (II.EQ.0) GO TO 1210
      IF (ARG.LT.X(II)) GO TO 1200
      XL=X(II)
      XU=X(II+1)
      IP=II+2
      IL=0
      GO TO 1140
C
 1200 CONTINUE
C
C     ARGUMENT OUT OF RANGE,I.E. ARG LESS THAN X(1).
C
 1210 IERR=4
      IP=3
      GO TO 2010
C
 2010 WRITE(6,3000)  IERR , ARG
C
      FVALUE(J)=0.D0
      FDERIV(J,1)=0.D0
      FDERIV(J,2)=0.D0
C
      II=IP-2
      XL=X(II)
      XU=X(II+1)
      IL=0
      GO TO 1155
C
C
C     END OF INTERPOLATION LOOP
C
 1160 CONTINUE
C
C     CALCULATION OF INTERPOLATED VALUES FINISHED.
C
      RETURN
C
C     PRINT ERROR MESSAGES.
C
 2000 IF (IERR.EQ.1) WRITE(6,3000)  IERR
      IF (IERR.EQ.2) WRITE(6,3000)  IERR , X(IN) , X(IN+1)
      RETURN
C
 3000 FORMAT(//5X,'*** SUBROUTINE SPLIN3 ERROR NO ',I2,' ***',
     $       2(4X,E21.14))
C
      END
C=======================================================================

      SUBROUTINE FRAG_VLNCE(IDX,LBAD)

C-----------------------------------------------------------------------
C     routine that fragments a quark - quark system               \FR'14
C
C     INPUT: IDX : parton stack index of central string
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER IDX,LBAD

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
C     The final particle output is contained in COMMON /S_PLIST/    
C     NP           : number of final particles
C     P(1:NP, 1:5) : 4-momenta + masses of the final particles 
C     LLIST (1:NP) : codes of final particles
      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP
      INTEGER NW_max
      PARAMETER (NW_max = 20)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
C     parameters that represent: NW: max. number of wounded nucleons,
C     NS,NH: max. number of soft and hard interactions
c      PARAMETER (NW_max = 20)
C     The COMMON block /S_CHIST/ contains information about the
C     the structure of the  generated event:
C     NWD   = number of wounded nucleons
C     NJET = total number of hard interactions
C     NSOF = total number of soft interactions
C     NNSOF (1:NW) = number of soft pomeron cuts in each interaction
C     NNJET (1:NW) = number of minijets produced in each interaction 
C     JDIF(1:NW) = diffraction code 
C                  0 : non-diff,
C                  1 : beam-diff
C                  2 : target-diff
C                  3 : double-diff
      INTEGER NNSOF,NNJET,JDIF,NWD,NJET,NSOF
      COMMON /S_CHIST/ NNSOF(NW_max),NNJET(NW_max),
     &     JDIF(NW_max),NWD,NJET,NSOF


      INTEGER NFORIG,NPORIG,NIORIG,IPFLAG,IIFLAG,KINT
      COMMON /S_PARTO/ NFORIG(NP_max),NPORIG(NP_max),NIORIG(NP_max),
     &IPFLAG,IIFLAG,KINT
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN

      DOUBLE PRECISION PST,PBM,PTG,PSTH,P1,P2,GABE,EE,
     &     PAR1_def,PAR24_def,PX1,PY1,PX2,PY2,GAM,BET,P1TOT,P2TOT,
     &     SIF,COF,COD,SID,ANORF,PZ
      DIMENSION PST(5),PBM(5),PTG(5),PSTH(5),P1(4),P2(4),GABE(4)
      INTEGER LSTH,IPID,IBMST,ITGST,ISTH,IFLB,IFLT,IST,I,IFBAD,JJ,
     &     NOLD,II,K,J
      SAVE
      
      LBAD = 2
      LSTH = 0
      
c     references are:
c     string --> bm-parton --> tg-parton (--> merged string/hadron)
c     read string 4momentum from stack
      CALL RD_PRTN_4VEC(IDX,PST,IPID,IBMST)
      CALL RD_PRTN_4VEC(IBMST,PBM,IFLB,ITGST)
      CALL RD_PRTN_4VEC(ITGST,PTG,IFLT,ISTH)

C     kinematic variables
      EE = PST(5)            ! string mass
      
      IF(NDEBUG.gt.1) WRITE(LUN,*)' FRAG_VLNCE: IDX,EE,IFLB,IFLT',
     &     IDX,EE,IFLB,IFLT

      IF(IDX.ne.ISTH) then
c     read merged string and add hadron to final particle stack..
         CALL RD_PRTN_4VEC(ISTH,PstH,LSTH,IST)
         IF(NDEBUG.gt.1) WRITE(LUN,*)' FRAG_VLNCE: found merged string',
     &        LSTH,(PSTH(I),I=1,5)         
         IF(IDX.ne.IST) then
            write(lun,*) ' FRAG_VLNCE: reference loop broken!' , IDX
            CALL SIB_REJECT('FRAG_VLNCE      ')
         endif
         NP = NP + 1
         DO I=1,4
            P(NP,I) = PST(I)
         ENDDO
         P(NP,5) = AM(IABS(LSTH))
         LLIST(NP) = LSTH
         NPORIG(NP) = IPFLAG*2+KINT
         niorig(NP) = iiflag
         LBAD = 0
         RETURN
      ENDIF

c     baryon production setup
      PAR1_def = PAR(1)
      if( NSOF+NJET.gt.0) then
         PAR(1)= PAR(15)
      else
         PAR(1)= PAR(14)
      endif

c     charm fractions in different parameterizations
      PAR24_def = PAR(24)
      IF(IPAR(15).gt.2.and.IPAR(15).ne.7.and.IPAR(15).lt.12)THEN
         PAR(24) = PAR(25)*EXP(-PAR(26)/EE)
      ENDIF

      IF(NDEBUG.gt.1)
     &     WRITE(LUN,*)' FRAG_VLNCE: parameters (CHM,DIQ,STR,VEC,POP)',
     &     PAR(24),PAR(1),PAR(2),PAR(5),PAR(8)

      NOLD=NP
      IF(IPAR(38).eq.1.or.IPAR(38).eq.2)THEN
C...  rotate strings instead of attaching all pt to string end hadrons
         PX1 = 0.D0
         PY1 = 0.D0
         PX2 = 0.D0
         PY2 = 0.D0
      ELSEIF(IPAR(38).eq.0.or.IPAR(38).eq.3)THEN
c     assign pt to hadrons at string end (old model)
         PX1 = PBM(1)
         PY1 = PBM(2)
         PX2 = PTG(1)
         PY2 = PTG(2)
         GAM = PST(4)/EE
         BET = PST(3)/PST(4)
      ENDIF         

C...  fragment strings in string restframe
      CALL STRING_FRAG_4FLV
     &     (EE,IFLB,IFLT,PX1,PY1,PX2,PY2,IFBAD,1)

      PAR(24) = PAR24_def
      PAR(1) = PAR1_def
      KINT= 0
      IF (IFBAD .EQ. 1) then
         if(Ndebug.gt.1) 
     &        WRITE(LUN,*)' STRING_FRAG: rejection (Ncall):',Ncall
         RETURN
      ENDIF

C...  rotate and boost string
      IF(IPAR(38).eq.1.or.IPAR(38).eq.2)THEN
C     boost quark momentum to string center-of-mass 
c     to calculate rotation angles in string center-of-mass
         do jj=1,3
            gabe(jj) = PST(jj)/PST(5)
         enddo
         GABE(4) = PST(4)/PST(5)
         CALL SIB_ALTRA(gabe(4),-gabe(1),-gabe(2),-gabe(3),
     &        PBM(1),pbm(2),pbm(3),pbm(4),
     &        P1TOT,p1(1),p1(2),p1(3),p1(4))
         CALL SIB_ALTRA(gabe(4),-gabe(1),-gabe(2),-gabe(3),
     &        PTG(1),pTG(2),ptg(3),ptg(4),
     &        P2TOT,p2(1),p2(2),p2(3),p2(4))

c     should be back-to-back...
         IF(ndebug.gt.1)THEN
            write(lun,*)
     &      ' FRAG_VLNCE: string c.m. momentum, parton 1 (Pabs,P(i)):' ,
     &           P1TOT, (P1(j),j=1,4)
            write(lun,*)
     &      ' FRAG_VLNCE: string c.m. momentum, parton 2 (Pabs,P(i)):' ,
     &           P2TOT, (P2(j),j=1,4)
            write(lun,*) '  partons should be back to back...'
         ENDIF
c     rotation factors
         COD= P1(3)/P1TOT
         SID= DSQRT(P1(1)**2+P1(2)**2)/P1TOT
         COF=1.D0
         SIF=0.D0
         IF(P1TOT*SID.GT.EPS5) THEN
            COF=P1(1)/(SID*P1TOT)
            SIF=P1(2)/(SID*P1TOT)
            ANORF=DSQRT(COF*COF+SIF*SIF)
            COF=COF/ANORF
            SIF=SIF/ANORF
         ENDIF
c     rotate string final state
         DO K=NOLD+1,NP
            CALL SIB_TRANI(P(K,1),P(k,2),P(k,3),cod,sid,cof,sif
     &           ,P2(1),P2(2),P2(3))
            do ii=1,3
               P(K,ii)=P2(ii)
            enddo
         ENDDO
c     boost to hadron - hadron center-of-mass
         DO K=NOLD+1,NP
            CALL SIB_ALTRA(gabe(4),gabe(1),gabe(2),
     &           gabe(3),P(k,1),p(k,2),p(k,3),p(k,4),
     &           P1TOT,p2(1),p2(2),p2(3),p2(4))
            do ii=1,4
               P(K,ii)=P2(ii)
            enddo
         ENDDO
      ELSEIF(IPAR(38).eq.0.or.IPAR(38).eq.3)THEN
C...  boost string
         DO K=NOLD+1,NP
            PZ = P(K,3)
            P(K,3) = GAM*(PZ+BET*P(K,4))
            P(K,4) = GAM*(P(K,4)+BET*PZ)
         ENDDO
      ENDIF
      LBAD = 0
      END


C-----------------------------------------------------------------------
C     fragmentation functions in SIBYLL                        \FR'14
C=======================================================================

      FUNCTION ZDIS_4FLV (IFL1,IFL2, XMT2)

C-----------------------------------------------------------------------
C...z distribution
c     includes charmed fragmentation (Peterson/SLAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      DOUBLE PRECISION FAin, FB0in
      COMMON /S_CZDIS/ FAin, FB0in

      DOUBLE PRECISION FAs1, fAs2
      COMMON /S_CZDISs/ FAs1, fAs2
      DOUBLE PRECISION ZDMAX, EPSI
      COMMON /S_CZDISc/ ZDMAX, EPSI
      INTEGER NW_max
      PARAMETER (NW_max = 20)
C--------------------------------------------------------------------
C     SIBYLL common blocks containing event information       \FR'14
C--------------------------------------------------------------------

C     EVENT INFO COMMON
C     contains overall interaction properties, like
C     SQS : center-of-mass energy
C     S   :         "       "     squared
C     PTmin : low pt cut of QCD cross section, 
C             i.e. minimal pt of hard minijets
C     Xmin : low-x bound for PDFs, 
C            i.e. minimal momentum fraction of hard partons
C     Zmin : logarithm of that
C     KB : PID of beam hadron
C     KT() : PID of target
C     IAT : mass number of target
      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
      INTEGER KB,IAT,KT
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      SAVE

      IAFL1 = IABS(mod(IFL1,100))
      IAFL2 = IABS(mod(IFL2,100))
c     SLAC-Peterson fragmentation function for charm
      IF ((IAFL1/10.eq.4.or.mod(IAFL1,10).eq.4)
     +     .or.(IAFL2/10.eq.4.or.mod(IAFL2,10).eq.4))THEN
 90      z = max(S_RNDM(0),1.e-8)
         tcp = zmefn(z,epsi)/zdmax
         if (tcp .lt. S_RNDM(1)) goto 90
         zdis_4flv = z
      else                      
c     original lund function, non charm
         fa=fain                ! lund parameter a
         fb0=fb0in              ! lund parameter b
c     parameters for hard scattering (gluon) fragmentation
         IF(IPAR(6).eq.2)THEN
            fa= PAR(18)
            fb0= PAR(19)
         ENDIF   
c     special parameters for strange fragmentation
c     only active for baryon beams (or K0,K0bar)
C     DH   correction  may 10-1996
         if (iabs(kb).ge.13) then ! baryons only
            if (iafl2.eq.3)  fa=fain+fas2
            if (iafl1.eq.3)  fa=fain+fas1
         endif
c     special parameters for baryon fragmentation
c     similar to pythia
         IF((IAFL1+IAFL2).gt.10.and.
     &        (IPAR(36).eq.1.or.IPAR(20).eq.3))then
            fa = fain + PAR(45)
            fb0 = PAR(60)
         ENDIF        
         FB = FB0*XMT2
         IF(FA.GT.0.01D0.AND.ABS(FA-1.D0)/FB.LE.0.01D0)
     +         ZMAX=FB/(1.D0+FB)+(1.D0-FA)*FB**2/(1.D0+FB)**3
         IF(FA.GT.0.01D0.AND.ABS(FA-1.D0)/FB.GT.0.01D0)
     +     ZMAX=0.5D0*(1.D0+FB-DSQRT((1.D0-FB)**2+4.D0*FA*FB))/(1.D0-FA)
         IF(ZMAX.LT.0.1D0)  ZDIV=2.75D0*ZMAX
         IF(ZMAX.GT.0.85D0) 
     +        ZDIV=ZMAX-0.6D0/FB**2+(FA/FB)*dLOG((0.01D0+FA)/FB)
C...  Choice if z, preweighted for peaks at low or high z
 100     Z=max(S_RNDM(0),1.e-8)
         IDIV=1
         FPRE=1.D0
         IF (ZMAX.LT.0.1D0)  THEN
            IF(1.D0.LT.S_RNDM(1)*(1.D0-dLOG(ZDIV)))  IDIV=2
            IF (IDIV.EQ.1)  Z=ZDIV*Z
            IF (IDIV.EQ.2)  Z=ZDIV**Z
            IF (IDIV.EQ.2)  FPRE=ZDIV/Z
         ELSEIF (ZMAX.GT.0.85D0)  THEN
            IF(1.D0.LT.S_RNDM(2)*(FB*(1.D0-ZDIV)+1.D0)) IDIV=2
            IF (IDIV.EQ.1)  Z=ZDIV+dLOG(Z)/FB
            IF (IDIV.EQ.1)  FPRE=dEXP(FB*(Z-ZDIV))
            IF (IDIV.EQ.2)  Z=ZDIV+Z*(1.D0-ZDIV)
         ENDIF
C...weighting according to the correct formula
         IF (Z.LE.FB/(50.D0+FB).OR.Z.GE.1.D0)  GOTO 100
         FVAL=(ZMAX/Z)*dEXP(FB*(1.D0/ZMAX-1.D0/Z))
         IF(FA.GT.0.01D0)  FVAL=((1.D0-Z)/(1.D0-ZMAX))**FA*FVAL
         IF(FVAL.LT.S_RNDM(3)*FPRE)  GOTO 100
         ZDIS_4FLV=Z
         
      ENDIF
      
      RETURN
      END
C=======================================================================
      
      SUBROUTINE ZNORMAL

C-----------------------------------------------------------------------
C...  normalisation for Peterson/SLAC frag. func
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      DOUBLE PRECISION ZDMAX, EPSI
      COMMON /S_CZDISc/ ZDMAX, EPSI

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      SAVE

c     get the maximum zmefn value first for normalisation
      jmax = 1000
      zdmax = 1.D-10

      DO j = 1, jmax
         z = dble(j)/dble(jmax+1)
         zdmax = max(zdmax, zmefn(z,epsi))
      enddo
      if (ndebug .gt. 0) WRITE(LUN,*)' ZDMAX,EPS:',zdmax, epsi
      RETURN
      END
C=======================================================================

      FUNCTION ZMEFN(z,eps)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      SAVE

C...  Peterson/SLAC frag. func
cdh   zmefn = (z*(1.D0-z**(-1)-eps/(1.D0-z))**2)**(-1)
      zmefn = 1.D0/(z*(1.D0-z**(-1)-eps/(1.D0-z))**2)
      RETURN
      END

C=======================================================================

      FUNCTION ZBLEAD (LB)

C-----------------------------------------------------------------------
C...fragmentation function for leading baryon
C.  simple form:  f(z) = a + x**b
C   INPUT : LB = particle code.
C..................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      DOUBLE PRECISION CLEAD, FLEAD
      COMMON /S_CZLEAD/ CLEAD, FLEAD
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)

      INTEGER ICHP,ISTR,IBAR
      COMMON /S_CHP/ ICHP(99), ISTR(99), IBAR(99)

      INTEGER IISO,ISPN
      COMMON /S_SPN/ IISO(99), ISPN(99)

      INTEGER ICHM
      COMMON /S_CHM/ ICHM(99)

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      SAVE

c      ncall = ncall + 1
c      print*,'leading baryon frag. called:',lb,ncall

C...  leading z lower bound
c     used for protons only in Sib21 (if ..)
c     used for all baryons alike in Sib22 (else..)    
      ZLMIN = PAR(55)
      ZSMR = PAR(56)

      IF(IPAR(30).ne.0)THEN
C     Sibyll 2.1 hard fragmentation function

        IC = ICHP(LB)*ISIGN(1,LB)
      
        if (LB.ge.34.and.LB.le.39)  then ! Lambda's and Sigma's
           IF(IPAR(35).eq.1)then
              zblead=zdisn(1)     ! zblead**2   !soft
           ELSE
 665          ZBLEAD = S_RNDM(LB)
              if (zblead.le.0.01D0) goto 665
           ENDIF
c     zblead=zdisn(1) ! blead**2   ! soft
        elseif (ic.eq.0)     then
           if(IPAR(30).eq.2)then
 555          zblead = S_RNDM(1)
              if (zblead .le. 0.01D0) goto 555     
           else
              zblead=zdisn(1)     ! blead**2   !soft
           endif
        elseif (ic.eq.1)  then   ! fast protons only
           if (abs(lb).eq.13) then
 661          IF (S_RNDM(2) .LT. CLEAD)  THEN
 666             ZBLEAD = S_RNDM(0)
                 if (zblead.le.0.01D0) goto 666
              ELSE
                 zblead=1.D0-zdisn(1) ! zblead**2   !hard
              ENDIF
c     truncated zblead to fix antiprotons
              if (zblead.le.ZLMIN+ZSMR*(1.D0-2.D0*S_RNDM(LB))) goto 661
           else
              zblead=zdisn(1)     ! zblead**2   !hard
           endif   
        else if (ic.eq.2)  then   ! fast delta++
           zblead=1.D0- zdisn(1)    ! (zblead)**.3333
        else
           zblead=S_RNDM(0)       ! zdisn(1)     !hard
        endif
        RETURN
      ELSE
C...  Sein's flat baryon fragmentation function a.k.a. Sibyll 2.2
 999     zblead = S_RNDM(0)
         if (zblead .le. 0.01D0) goto 999
c     truncated zblead to fix instring pair production (antiprotons)
         if (zblead.le.ZLMIN+ZSMR*(1.D0-2.D0*S_RNDM(LB))) goto 999
         RETURN
      ENDIF
      END

C=======================================================================

      FUNCTION ZDISN (n)

C-----------------------------------------------------------------------
C...Generate (1-x)**n
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      SAVE

666   rmin=1.1D0
      do i=1,n+1
         R1=S_RNDM(i)
         IF (R1.LE.RMIN) RMIN=R1
      ENDDO
      ZDISn=RMIN
      if (zdisn.le.0.01D0) goto 666
      if (zdisn.ge.0.99D0) goto 666
      END
C=======================================================================

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
      IMPLICIT INTEGER(I-N)

      INTEGER NS_max, NH_max
      PARAMETER (NS_max = 20, NH_max = 80)
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

      PARAMETER ( NPARFIT = 22 )
      DOUBLE PRECISION PARS
      COMMON /XSCTN_FIT/ PARS( 50 , 2 )

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN

      COMMON /QCD_XSCTN/SIGQCD(61,2),INIT
      DOUBLE PRECISION SIGQCD
      SAVE
      DATA INIT /0/
      DATA (SIGQCD(K,1),K=    1,   61) /
     &8.4663D-02,1.8246D-01,3.3880D-01,5.6845D-01,8.8686D-01,1.3116D+00,
     &1.8626D+00,2.5645D+00,3.4445D+00,4.5343D+00,5.8715D+00,7.4962D+00,
     &9.4579D+00,1.1811D+01,1.4620D+01,1.7955D+01,2.1890D+01,2.6522D+01,
     &3.1952D+01,3.8303D+01,4.5704D+01,5.4307D+01,6.4284D+01,7.5818D+01,
     &8.9121D+01,1.0447D+02,1.2213D+02,1.4240D+02,1.6562D+02,1.9221D+02,
     &2.2260D+02,2.5733D+02,2.9694D+02,3.4207D+02,3.9348D+02,4.5194D+02,
     &5.1838D+02,5.9376D+02,6.7921D+02,7.7609D+02,8.8578D+02,1.0099D+03,
     &1.1504D+03,1.3090D+03,1.4882D+03,1.6903D+03,1.9183D+03,2.1754D+03,
     &2.4650D+03,2.7912D+03,3.1582D+03,3.5707D+03,4.0341D+03,4.5538D+03,
     &5.1360D+03,5.7883D+03,6.5193D+03,7.3358D+03,8.2428D+03,9.2498D+03,
     &1.0369D+04/
      DATA (SIGQCD(K,2),K=    1,   61) /
     &1.5665D-01,2.8800D-01,4.7863D-01,7.4235D-01,1.0949D+00,1.5547D+00,
     &2.1433D+00,2.8859D+00,3.8118D+00,4.9547D+00,6.3534D+00,8.0525D+00,
     &1.0103D+01,1.2563D+01,1.5498D+01,1.8986D+01,2.3111D+01,2.7971D+01,
     &3.3678D+01,4.0358D+01,4.8154D+01,5.7228D+01,6.7762D+01,7.9965D+01,
     &9.4071D+01,1.1034D+02,1.2909D+02,1.5063D+02,1.7536D+02,2.0370D+02,
     &2.3613D+02,2.7321D+02,3.1553D+02,3.6379D+02,4.1875D+02,4.8129D+02,
     &5.5238D+02,6.3311D+02,7.2470D+02,8.2854D+02,9.4614D+02,1.0792D+03,
     &1.2298D+03,1.3999D+03,1.5920D+03,1.8089D+03,2.0534D+03,2.3291D+03,
     &2.6396D+03,2.9892D+03,3.3825D+03,3.8248D+03,4.3219D+03,4.8803D+03,
     &5.5072D+03,6.2109D+03,7.0001D+03,7.8849D+03,8.8764D+03,9.9871D+03,
     &1.1231D+04/


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

      PTCUT = XI(10)+XI(21)*dEXP(XI(22)*DSQRT(2.D0*dLOG(ECM)))
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
      SLOPEc = SIG_tot**2/(16.D0*PI*SIG_ela)

      DE = ABS(SIGEL+SIGINE-SIGTOT)/SIGTOT
      IF(DE.GT.0.01D0) THEN
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

C=======================================================================

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
      IMPLICIT INTEGER(I-N)

      DIMENSION SIG_brn(3)
      PARAMETER (NS_max = 20, NH_max = 80)

      COMMON /S_CFACT/ FACT(0:NH_max), CO_BIN(0:NH_max,0:NH_max)
      COMMON /S_CHDCNV/ABR(2,400),ABP(2,400),ABH(2,400),DB,NB

      COMMON /PROFILE/XNUS2,XMUS2,XNUSPI2,
     &                XNUH2,XMUH2,XNUHPI2,
     &                EnhPP,EnhPiP,al1,be1,al2,be2

      DIMENSION SIG_DIF1(2),SIG_DIF2(2),SIG_DD(2),
     &          P_int(0:NS_max,0:NH_max)

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN

      DOUBLE PRECISION EPS100
      SAVE
      DATA EPS100 /1.D-100/
      
      DO J=0,NH_max
        DO I=0,NS_max
          P_int(I,J) = 0.D0
        ENDDO
      ENDDO

      ga1 = dsqrt(al1*al1+be1*be1)
      ga2 = dsqrt(al2*al2+be2*be2)

      fe_a_1  = (1.D0+al1/ga1)/2.D0
      fe_a_2  = (1.D0-al1/ga1)/2.D0
      fd_a_1  = sqrt(1.D0-(al1/ga1)**2)/2.D0
      fd_a_2  = -fd_a_1

      fe_b_1  = (1.D0+al2/ga2)/2.D0
      fe_b_2  = (1.D0-al2/ga2)/2.D0
      fd_b_1  = dsqrt(1.D0-(al2/ga2)**2)/2.D0
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


      sum_abs = 0.D0
      sum_tot = 0.D0
      sum_ela = 0.D0
      sum_sd_a = 0.D0
      sum_sd_b = 0.D0
      sum_dd  = 0.D0
      sum_B   = 0.D0

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
         chi2_soft_11 = (1.D0-al1+ga1)*(1.D0-al2+ga2)*chi2_soft
         chi2_soft_22 = (1.D0-al1-ga1)*(1.D0-al2-ga2)*chi2_soft
         chi2_soft_12 = (1.D0-al1+ga1)*(1.D0-al2-ga2)*chi2_soft
         chi2_soft_21 = (1.D0-al1-ga1)*(1.D0-al2+ga2)*chi2_soft

         chi2_hard = ABHAR*SIG_HAR
         chi2_hard_11 = (1.D0-al1+ga1)*(1.D0-al2+ga2)*chi2_hard
         chi2_hard_22 = (1.D0-al1-ga1)*(1.D0-al2-ga2)*chi2_hard
         chi2_hard_12 = (1.D0-al1+ga1)*(1.D0-al2-ga2)*chi2_hard
         chi2_hard_21 = (1.D0-al1-ga1)*(1.D0-al2+ga2)*chi2_hard
          
         ef_11 = max(-0.5D0*(chi2_soft_11+chi2_hard_11),log(EPS100))
         ef_22 = max(-0.5D0*(chi2_soft_22+chi2_hard_22),log(EPS100))
         ef_12 = max(-0.5D0*(chi2_soft_12+chi2_hard_12),log(EPS100))
         ef_21 = max(-0.5D0*(chi2_soft_21+chi2_hard_21),log(EPS100))
         
         ef_11 = dexp(ef_11)
         ef_22 = dexp(ef_22)
         ef_12 = dexp(ef_12)
         ef_21 = dexp(ef_21)

         esf_11  = max(ef_11,EPS100)**2
         esf_22  = max(ef_22,EPS100)**2
         esf_12  = max(ef_12,EPS100)**2
         esf_21  = max(ef_21,EPS100)**2

         F_ine = B*(1.D0 - fe_11*esf_11 - fe_12*esf_12 
     &                   - fe_21*esf_21 - fe_22*esf_22)
         F_tot = 1.D0 - fe_11*ef_11 - fe_12*ef_12
     &                - fe_21*ef_21 - fe_22*ef_22
         F_ela = B*F_tot**2
         F_tot = B*F_tot

         F_sd_a = B*(fd_a_11*ef_11 + fd_a_12*ef_12
     &             + fd_a_21*ef_21 + fd_a_22*ef_22)**2
         F_sd_b = B*(fd_b_11*ef_11 + fd_b_12*ef_12
     &             + fd_b_21*ef_21 + fd_b_22*ef_22)**2
         F_dd  = B*(fdd_11*ef_11 + fdd_12*ef_12
     &            + fdd_21*ef_21 + fdd_22*ef_22)**2

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
         soft_rec_11 = 1.D0/chi2_soft_11
         soft_rec_22 = 1.D0/chi2_soft_22
         soft_rec_12 = 1.D0/chi2_soft_12
         soft_rec_21 = 1.D0/chi2_soft_21
         chi2_hard_11 = max(chi2_hard_11,EPS100)
         chi2_hard_22 = max(chi2_hard_22,EPS100)
         chi2_hard_12 = max(chi2_hard_12,EPS100)
         chi2_hard_21 = max(chi2_hard_21,EPS100)
         DO I=0,I0MAX
           soft_rec_11 = max(soft_rec_11*chi2_soft_11,EPS100)
           soft_rec_22 = max(soft_rec_22*chi2_soft_22,EPS100)
           soft_rec_12 = max(soft_rec_12*chi2_soft_12,EPS100)
           soft_rec_21 = max(soft_rec_21*chi2_soft_21,EPS100)
           hard_rec_11 = max(1.D0/chi2_hard_11,EPS100)
           hard_rec_22 = max(1.D0/chi2_hard_22,EPS100)
           hard_rec_12 = max(1.D0/chi2_hard_12,EPS100)
           hard_rec_21 = max(1.D0/chi2_hard_21,EPS100)
           DO J=0,J0MAX
             hard_rec_11 = max(hard_rec_11*chi2_hard_11,EPS100)
             hard_rec_22 = max(hard_rec_22*chi2_hard_22,EPS100)
             hard_rec_12 = max(hard_rec_12*chi2_hard_12,EPS100)
             hard_rec_21 = max(hard_rec_21*chi2_hard_21,EPS100)
             P_int(I,J) = P_int(I,J) 
     &                + fe_11*soft_rec_11*hard_rec_11*fac_11
     &                + fe_22*soft_rec_22*hard_rec_22*fac_22
     &                + fe_12*soft_rec_12*hard_rec_12*fac_12
     &                + fe_21*soft_rec_21*hard_rec_21*fac_21
           ENDDO
         ENDDO

      ENDDO

      SIG_abs  = SUM_abs*TWOPI*DB
      SIG_tot  = SUM_tot*4.D0*PI*DB
      SIG_ela  = SUM_ela*TWOPI*DB
      SIG_dif1(1) = SUM_sd_a*TWOPI*DB
      SIG_dif2(1) = SUM_sd_b*TWOPI*DB
      SIG_dd(1)   = SUM_dd*TWOPI*DB
      SIG_ine  = SIG_abs + SIG_dif1(1) + SIG_dif2(1) + SIG_dd(1)
      B_EL     = sum_B/SUM_tot/2.D0

      SA = 0.D0
      P_int(0,0) = 0.D0
      DO I=0,I0MAX
        DO J=0,J0MAX
          fac = FACT(I)*FACT(J)
          P_int(I,J) = P_int(I,J)/fac
          SA = SA + P_int(I,J)
        ENDDO
      ENDDO

      SIG_hmsd = EnhPP*(P_int(1,0)+P_int(0,1))*TWOPI*DB
      SIG_hmdd = be1**2*SIG_hmsd + be2**2*SIG_hmsd
     &          + EnhPP**2*P_int(1,1)*TWOPI*DB

      SIG_dif1(2) = SIG_hmsd
      SIG_dif2(2) = SIG_hmsd
      SIG_dd(2)   = SIG_hmdd

      SIG_sum = SA*TWOPI*DB

      DO I=0,I0MAX
        DO J=0,J0MAX
          P_int(I,J) = P_int(I,J)/SA
        ENDDO
      ENDDO

      END

C=======================================================================

      SUBROUTINE HAD_CONV(JINT)

C-----------------------------------------------------------------------
C
C...Convolution of hadrons profile
C.  [function A(b) of Durand and Pi]
C.  precalculate and put  in COMMON block
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER(I-N)
C
      COMMON /S_CHDCNV/ABR(2,400),ABP(2,400),ABH(2,400),DB,NB

      DOUBLE PRECISION NU2, MU2, NUPI2, NU, MU, NUPI
      COMMON /S_CH0CNV/ NU2, MU2, NUPI2, NU, MU, NUPI

C
      COMMON /PROFILE/XNUS2,XMUS2,XNUSPI2,
     &                XNUH2,XMUH2,XNUHPI2,
     &                ENHPP,ENHPIP,al1,be1,al2,be2
      SAVE

C...integration constants
      BMAX = 50.D0
      NB  = 400
      DB = BMAX/DBLE(NB)

C  soft reggeon interactions

      NU2   = XNUS2
      MU2   = XMUS2
      NUPI2 = XNUSPI2

      NU = SQRT(NU2)
      MU = SQRT(ABS(MU2))
      NUPI = SQRT(NUPI2) 

      DO JB=1,NB
         B = DB*DBLE(JB-1)
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
         B = DB*DBLE(JB-1)
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

      DB = BMAX/DBLE(NB)
      DO JB=1,NB
         B = DB*DBLE(JB-1)
         IF(JINT.EQ.1) THEN
           ABH(JINT,JB) = A_PP(B)
         ELSE
           ABH(JINT,JB) = A_PIP(B)
         ENDIF
      ENDDO

      END

C=======================================================================

      DOUBLE PRECISION FUNCTION A_PP (b)

C-----------------------------------------------------------------------
C...Convolution of parton distribution for pp interaction
      IMPLICIT DOUBLE PRECISION (A-Z)
C
      DOUBLE PRECISION NU2, MU2, NUPI2, NU, MU, NUPI
      COMMON /S_CH0CNV/ NU2, MU2, NUPI2, NU, MU, NUPI
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN
      SAVE

      ETA = NU2/MU2
 
      IF(ETA.LT.0.D0) THEN
   
        c = nu**5/(96.D0*PI)
        if (b .gt. 0.0001D0)  then
           A_pp = c*b**3 * bessk (3, b*nu)
        else
           A_pp = nu**2/(12.D0*PI)
        endif

      ELSE

        X = B*NU
        Y = B*MU
        C = NU2/(12.D0*PI)/(1.D0-ETA)**2
        IF(X.GT.0.0001D0) THEN
          A_PP = C*(1.D0/8.D0*X**3*BESSK(3,X)
     &          -3.D0/2.D0*ETA/(1.D0-ETA)*X**2*BESSK(2,X)
     &          + 9.D0*ETA**2/(1.D0-ETA)**2*X*BESSK1(X)
     &          -24.D0*ETA**3/(1.D0-ETA)**3*(BESSK0(X)-BESSK0(Y))
     &          + 3.D0*ETA**3/(1.D0-ETA)**2*Y*BESSK1(Y))
        ELSE
          A_PP = C*(1.D0 /8.D0*8.D0
     &          -3.D0/2.D0*ETA/(1.D0-ETA)*2.D0
     &          +9.D0*ETA**2/(1.D0-ETA)**2*1.D0
     &          -24.D0*ETA**3/(1.D0-ETA)**3*LOG(MU/NU)
     &          +3.D0*ETA**3/(1.D0-ETA)**2*1.D0)
        ENDIF

      ENDIF

      END

C=======================================================================

      DOUBLE PRECISION FUNCTION A_PIP (b)

C-----------------------------------------------------------------------
C...Convolution of parton distribution for pip interaction
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-Z)
C
      DOUBLE PRECISION NU2, MU2, NUPI2, NU, MU, NUPI
      COMMON /S_CH0CNV/ NU2, MU2, NUPI2, NU, MU, NUPI
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN
      SAVE

      eta = nu2/nupi2
      c = nu2/(2.D0*PI) * 1.D0/(1.D0-eta)

      if (b .gt. 0.0001D0)  then
         b1 = b*nu
         b2 = b*nupi
         f1 = 0.5D0*b1 * bessk1(b1)
         f2 = eta/(1.D0-eta)*(bessk0(b2)- bessk0(b1))
         A_pip = c*(f1+f2)
      else
         A_pip = c*(0.5D0 + eta/(1.D0-eta)*log(nu/nupi))
      endif
      return
      end
C
C
C-----------------------------------------------------------------------
C  Bessel functions
C=======================================================================

      FUNCTION BESSK0(X)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER(I-N)
C
      DOUBLE PRECISION P1,P2,P3,P4,P5,P6,P7,
     *                 Q1,Q2,Q3,Q4,Q5,Q6,Q7
      SAVE
      DATA P1,P2,P3,P4,P5,P6,P7/-0.57721566D0,0.42278420D0,
     *    0.23069756D0,0.3488590D-1,0.262698D-2,0.10750D-3,0.74D-5/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,-0.7832358D-1,
     * 0.2189568D-1,-0.1062446D-1,0.587872D-2,-0.251540D-2,0.53208D-3/

      IF (X.LE.2.0D0) THEN
        Y=X*X/4.D0
        BESSK0=(-LOG(X/2.D0)*BESSI0(X))+(P1+Y*(P2+Y*(P3+
     *        Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        Y=(2.D0/X)
        BESSK0=(EXP(-X)/SQRT(X))*(Q1+Y*(Q2+Y*(Q3+
     *        Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END
C
C=======================================================================

      FUNCTION BESSK1(X)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER(I-N)
C
      DOUBLE PRECISION P1,P2,P3,P4,P5,P6,P7,
     *                 Q1,Q2,Q3,Q4,Q5,Q6,Q7
      SAVE
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,0.15443144D0,-0.67278579D0,
     *    -0.18156897D0,-0.1919402D-1,-0.110404D-2,-0.4686D-4/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,0.23498619D0,
     *    -0.3655620D-1,0.1504268D-1,-0.780353D-2,0.325614D-2,
     *    -0.68245D-3/

      IF (X.LE.2.D0) THEN
        Y=X*X/4.D0
        BESSK1=(LOG(X/2.D0)*BESSI1(X))+(1.D0/X)*(P1+Y*(P2+
     *      Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        Y=2.D0/X
        BESSK1=(EXP(-X)/SQRT(X))*(Q1+Y*(Q2+Y*(Q3+
     *      Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END
C
C=======================================================================

      FUNCTION BESSK(N,X)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      SAVE
C
      IF (N.LT.2) stop 'bad argument N in BESSK'
      TOX=2.D0/X
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
C=======================================================================

      FUNCTION BESSI0(X)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER(I-N)
C
      DOUBLE PRECISION P1,P2,P3,P4,P5,P6,P7,
     *                 Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
      SAVE
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,3.5156229D0,3.0899424D0,
     *    1.2067492D0, 0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,
     *    0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,
     *    0.2635537D-1,-0.1647633D-1,0.392377D-2/

      IF (ABS(X).LT.3.75D0) THEN
        Y=(X/3.75D0)**2
        BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
        AX=ABS(X)
        Y=3.75D0/AX
        BESSI0=(EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4
     *      +Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END
C
C=======================================================================

      FUNCTION BESSI1(X)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER(I-N)
C
      DOUBLE PRECISION P1,P2,P3,P4,P5,P6,P7,
     *                 Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
      SAVE
      DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,
     *    0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1,
     *    -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1,
     *    -0.2895312D-1,0.1787654D-1,-0.420059D-2/

      IF (ABS(X).LT.3.75D0) THEN
        Y=(X/3.75D0)**2
        BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        AX=ABS(X)
        Y=3.75D0/AX
        BESSI1=(EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+
     *      Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END

C=======================================================================

      SUBROUTINE FACT_INI

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NS_max, NH_max
      PARAMETER (NS_max = 20, NH_max = 80)
      COMMON /S_CFACT/ FACT(0:NH_max), CO_BIN(0:NH_max,0:NH_max)

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      SAVE
      
      FACT(0) = 1.D0
      FACT(NS_max) = 1.D0  ! avoid unused warning and keep parameter block
      DO J=1,NH_max
         FACT(J) = FACT(J-1)*DBLE(J)
      ENDDO
      DO J=0,NH_max
         DO K=0,J
            CO_BIN(J,K) = FACT(J)/(FACT(K)*FACT(J-K))
         ENDDO
      ENDDO

      END
cC=======================================================================
c
c      SUBROUTINE SAMPLE_SOFT (STR_mass_min, X1,X2,PT)
c
C-----------------------------------------------------------------------
C...  Routine for the sampling the kinematical variables of sea quarks
C.     according to (1-x)**b / x**2
C.  INPUT:  STR_mass_min : minimal string mass ** 2 = x1 * x2 * s
C.          SLOPE : large x suppression exponent
C.  OUTPUT:  gluon 4momenta on parton stack (GeV)                /FR'14
C-----------------------------------------------------------------------
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      IMPLICIT INTEGER(I-N)
c
c      INTEGER NW_max
c      PARAMETER (NW_max = 20)
C--------------------------------------------------------------------
C     SIBYLL common blocks containing event information       \FR'14
C--------------------------------------------------------------------
c
c     EVENT INFO COMMON
c     contains overall interaction properties, like
c     SQS : center-of-mass energy
C     S   :         "       "     squared
C     PTmin : low pt cut of QCD cross section, 
C             i.e. minimal pt of hard minijets
C     Xmin : low-x bound for PDFs, 
C            i.e. minimal momentum fraction of hard partons
C     Zmin : logarithm of that
C     KB : PID of beam hadron
C     KT() : PID of target
C     IAT : mass number of target
c      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
c      INTEGER KB,IAT,KT
c      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT
c
c      INTEGER NCALL, NDEBUG, LUN
c      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
c
c      DOUBLE PRECISION PPT02
c      COMMON /S_CQDIS2/ PPT02(44)
c      INTEGER NIPAR_max,NPAR_max
c      PARAMETER (NPAR_max=200,NIPAR_max=100)
c      DOUBLE PRECISION PAR
c      INTEGER IPAR
c      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
c
C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
c      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
c      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10
c
c      DOUBLE PRECISION PI,TWOPI,CMBARN
c      COMMON /SIB_CST/ PI,TWOPI,CMBARN
c
c      DOUBLE PRECISION FACN
c      DIMENSION FACN(3:10)
c      COMMON /SIB_FAC/ FACN
c      SAVE
c
c      SLOPE = max(1.D0,PAR(42))
c      ZSOF = 2.D0*dLOG(STR_mass_min/SQS) ! minim. mass ~ x1 * x2
c 50   XMIN = dEXP(ZSOF)
c      axmin = 1.D0/xmin
c 100  Z1 = -1.D0*dLOG(axmin-(axmin-1.D0)*S_RNDM(0))
c      x1 = dexp(z1)
c      XR = dlog(1.D0-X1) - dlog(1.D0-xmin)
c      if(SLOPE*XR.le.log(S_RNDM(0))) goto 100
c
c 200  Z2 = -1.D0*dLOG(axmin-(axmin-1.D0)*S_RNDM(0))
c      X2 = dEXP(Z2)
c      XR = dlog(1.D0-X2) - dlog(1.D0-dEXP(ZSOF))
c      if(SLOPE*XR.le.dlog(S_RNDM(0))) goto 200     
c
c      IF(Z1+Z2.LE.ZSOF) GOTO 50
c      STR_mass2 = dsqrt(X1*X2*S)/2.D0
c      PPTT = PPT02(10)
c 150  PT = PPTT*dSQRT(-dLOG(MAX(EPS10,S_RNDM(0))))
c      IF(IPAR(3).eq.6)THEN
c         XM = 0.D0
c         XM2 = XM**2
c         RNDM = MAX(EPS10,S_RNDM(IFL))
c         XMT = PPTT * dLOG(RNDM) - XM
c         XMT2 = XMT**2
c         PT = dSQRT(XMT2-XM2)
c      ENDIF
c      IF(PT.GT.PTmin) GOTO 150
c      IF(PT.GE.STR_mass2) GOTO 150
c      END
c
cC=======================================================================
c
c      SUBROUTINE SAMPLE_SOFT2 (STR_mass_min, X1,X2,PT)
c
C-----------------------------------------------------------------------
C...Routine for sampling the kinematical variables
C.  that characterize a soft cut pomeron (x1,x2, pT)
C.  from the differential cross section:
C.     d3sigma/(dx1 dx2 dpT)
C.      ~ 1/x_i**a .*. exp(-mT)
C.  INPUT: STR_mass_min : minimal string mass defined by kinematic limits
C.                        of the string fragmentation
C.  PAR:   PAR(42) : exponent a
C.  OUTPUT:  X1, X2, PT (GeV)
C-----------------------------------------------------------------------
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      IMPLICIT INTEGER(I-N)
c      INTEGER NW_max
c      PARAMETER (NW_max = 20)
C--------------------------------------------------------------------
C     SIBYLL common blocks containing event information       \FR'14
C--------------------------------------------------------------------
c
C     EVENT INFO COMMON
C     contains overall interaction properties, like
C     SQS : center-of-mass energy
C     S   :         "       "     squared
C     PTmin : low pt cut of QCD cross section, 
C             i.e. minimal pt of hard minijets
C     Xmin : low-x bound for PDFs, 
C            i.e. minimal momentum fraction of hard partons
C     Zmin : logarithm of that
C     KB : PID of beam hadron
C     KT() : PID of target
C     IAT : mass number of target
c      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
c      INTEGER KB,IAT,KT
c      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT
c
c      INTEGER NCALL, NDEBUG, LUN
c      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
c
c      DOUBLE PRECISION PPT02
c      COMMON /S_CQDIS2/ PPT02(44)
c      INTEGER NIPAR_max,NPAR_max
c      PARAMETER (NPAR_max=200,NIPAR_max=100)
c      DOUBLE PRECISION PAR
c      INTEGER IPAR
c      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
c
cC--------------------------------------------------------------------
cC     SIBYLL utility common blocks containing constants       \FR'14
cC--------------------------------------------------------------------
c      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
c      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10
c
c      DOUBLE PRECISION PI,TWOPI,CMBARN
c      COMMON /SIB_CST/ PI,TWOPI,CMBARN
c
c      DOUBLE PRECISION FACN
c      DIMENSION FACN(3:10)
c      COMMON /SIB_FAC/ FACN
c      SAVE
c
c      SLOPE = PAR(42)
c      ZSOF = 2.D0*dLOG(STR_mass_min/SQS) ! zmin
c      zsof = zsof * slope
c 100  Z1=1.D0/SLOPE*(-zsof*S_RNDM(0)+zsof)
c      Z2=1.D0/SLOPE*(-zsof*S_RNDM(0)+zsof)
cc      print *,'zsof,z1,z2',zsof,z1,z2
c      IF(Z1+Z2.LE.ZSOF) GOTO 100
c      X1=dEXP(Z1)
c      X2=dEXP(Z2)
c      STR_mass2 = sqrt(X1*X2*S)/2.D0
c      if(str_mass2.lt.0.9D0)goto 100
c      PPTT = PPT02(10)
cc      print *,'ptmin,str_mass:',ptmin,str_mass2
c 150  PT = PPTT*dSQRT(-dLOG(MAX(EPS10,S_RNDM(0))))
c      IF(IPAR(3).eq.6)THEN
c         XM = 0.D0
c         XM2 = XM**2
c         RNDM = MAX(EPS10,S_RNDM(IFL))
c         XMT = PPTT * dLOG(RNDM) - XM
c         XMT2 = XMT**2
c         PT = dSQRT(XMT2-XM2)
c      ENDIF
c      IF(PT.GT.PTmin) GOTO 150
c      IF(PT.GE.STR_mass2) GOTO 150
c      PHI = TWOPI*S_RNDM(L)
c      END
cC=======================================================================
cc
c      SUBROUTINE SAMPLE_SOFT3 (STR_mass_min, X1,X2,PT)
c
cC-----------------------------------------------------------------------
cC...Routine for the sampling the kinematical variables
cC.  that characterize a soft cut pomeron (x1,x2, pT)
cC.  from the differential cross section:
cC.     d3sigma/(dx1 dx2 dpT)
cC.  INPUT:  L=1 incident proton, L=2  incident pi
cC.          (soft strings identical for pi and p interactions)
cC.  OUTPUT:  X1, X2, PT (GeV)
cC-----------------------------------------------------------------------
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      IMPLICIT INTEGER(I-N)
c      INTEGER NW_max
c      PARAMETER (NW_max = 20)
cC--------------------------------------------------------------------
cC     SIBYLL common blocks containing event information       \FR'14
cC--------------------------------------------------------------------
c
cC     EVENT INFO COMMON
cC     contains overall interaction properties, like
cC     SQS : center-of-mass energy
cC     S   :         "       "     squared
cC     PTmin : low pt cut of QCD cross section, 
cC             i.e. minimal pt of hard minijets
cC     Xmin : low-x bound for PDFs, 
cC            i.e. minimal momentum fraction of hard partons
cC     Zmin : logarithm of that
cC     KB : PID of beam hadron
cC     KT() : PID of target
cC     IAT : mass number of target
c      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
c      INTEGER KB,IAT,KT
c      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT
c
c      INTEGER NCALL, NDEBUG, LUN
c      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
c
c      DOUBLE PRECISION PPT02
c      COMMON /S_CQDIS2/ PPT02(44)
c      INTEGER NIPAR_max,NPAR_max
c      PARAMETER (NPAR_max=200,NIPAR_max=100)
c      DOUBLE PRECISION PAR
c      INTEGER IPAR
c      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
c
cC--------------------------------------------------------------------
cC     SIBYLL utility common blocks containing constants       \FR'14
cC--------------------------------------------------------------------
c      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
c      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10
c
c      DOUBLE PRECISION PI,TWOPI,CMBARN
c      COMMON /SIB_CST/ PI,TWOPI,CMBARN
c
c      DOUBLE PRECISION FACN
c      DIMENSION FACN(3:10)
c      COMMON /SIB_FAC/ FACN
c      SAVE
c
c      SLOPE = max(1.D0,PAR(42))
c      ZSOF = 2.D0*dLOG(STR_mass_min/SQS) ! minim. mass ~ x1 * x2
c 100  Z1=-ZSOF*S_RNDM(0)+ZSOF   ! sample envelope 1/x
c      X1 = dEXP(Z1)
cc      print *,'z1,x1:',z1,x1
c      XR = dlog(1.D0-X1) - dlog(1.D0-dEXP(ZSOF))
cc      print *,'ratio:',(1.-X1)/(1.-EXP(ZSOF)),(1.-X1),1.-EXP(ZSOF)
cc      print *,'log ratio:',xr,log(1.-X1),log(1.-EXP(ZSOF))
c      if(SLOPE*XR.le.dlog(S_RNDM(0))) goto 100
c
c 200  Z2=-ZSOF*S_RNDM(0)+ZSOF   ! sample envelope 1/x
c      X2 = dEXP(Z2)
c      XR = dlog(1.D0-X2) - dlog(1.D0-dEXP(ZSOF))
c      if(SLOPE*XR.le.dlog(S_RNDM(0))) goto 200     
cc      print *,'zsof,z1,z2',zsof,z1,z2
c      IF(Z1+Z2.LE.ZSOF) GOTO 100
c      STR_mass2 = sqrt(X1*X2*S)/2.D0
c      PPTT = PPT02(10)
c      IF(IPAR(3).eq.8) PPTT = PPT02(20)
c 150  PT = PPTT*dSQRT(-dLOG(MAX(EPS10,S_RNDM(0))))
c      IF(IPAR(3).ge.6)THEN
c         XM = 0.D0
c         XM2 = XM**2
c         RNDM = MAX(EPS10,S_RNDM(IFL))
c         XMT = PPTT * dLOG(RNDM) - XM
c         XMT2 = XMT**2
c         PT = dSQRT(XMT2-XM2)
c      ENDIF
c      IF(PT.GT.PTmin) GOTO 150
c      IF(PT.GE.STR_mass2) GOTO 150
c      PHI = TWOPI*S_RNDM(L)
c      END
cC=======================================================================
c
c      SUBROUTINE SAMPLE_SOFT5 (STR_mass_min, X1,X2,PT)
c
cC-----------------------------------------------------------------------
cC...Routine for the sampling the kinematical variables of sea quarks
cC.     according to (1-x)**b / x**2
cC.  INPUT:  STR_mass_min : minimal string mass ** 2 = x1 * x2 * s
cC.          SLOPE : large x suppression exponent
cC.  OUTPUT:  X1, X2, PT (GeV)                                   /FR'14
cC-----------------------------------------------------------------------
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      IMPLICIT INTEGER(I-N)
c      INTEGER NW_max
c      PARAMETER (NW_max = 20)
cc      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, kb ,kt
cC--------------------------------------------------------------------
cC     SIBYLL common blocks containing event information       \FR'14
cC--------------------------------------------------------------------
c
cC     EVENT INFO COMMON
cC     contains overall interaction properties, like
cC     SQS : center-of-mass energy
cC     S   :         "       "     squared
cC     PTmin : low pt cut of QCD cross section, 
cC             i.e. minimal pt of hard minijets
cC     Xmin : low-x bound for PDFs, 
cC            i.e. minimal momentum fraction of hard partons
cC     Zmin : logarithm of that
cC     KB : PID of beam hadron
cC     KT() : PID of target
cC     IAT : mass number of target
c      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
c      INTEGER KB,IAT,KT
c      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT
c
c      INTEGER NCALL, NDEBUG, LUN
c      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
c
c      DOUBLE PRECISION PPT02
c      COMMON /S_CQDIS2/ PPT02(44)
c      INTEGER NIPAR_max,NPAR_max
c      PARAMETER (NPAR_max=200,NIPAR_max=100)
c      DOUBLE PRECISION PAR
c      INTEGER IPAR
c      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
c
cC--------------------------------------------------------------------
cC     SIBYLL utility common blocks containing constants       \FR'14
cC--------------------------------------------------------------------
c      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
c      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10
c
c      DOUBLE PRECISION PI,TWOPI,CMBARN
c      COMMON /SIB_CST/ PI,TWOPI,CMBARN
c
c      DOUBLE PRECISION FACN
c      DIMENSION FACN(3:10)
c      COMMON /SIB_FAC/ FACN
c      SAVE
c
c      SLOPE = max(1.D0,PAR(42))
c      ZSOF = 2.D0*dLOG(STR_mass_min/SQS) ! minim. mass ~ x1 * x2
c 50   XMIN = dEXP(ZSOF)
c      axmin = 1.D0/xmin
c 100  Z1 = -1.D0*dLOG(axmin-(axmin-1.D0)*S_RNDM(0))
c      x1 = dexp(z1)
c      XR = dlog(1.D0-X1) - dlog(1.D0-xmin)
c      if(SLOPE*XR.le.log(S_RNDM(0))) goto 100
c
c 200  Z2 = -1.D0*dLOG(axmin-(axmin-1.D0)*S_RNDM(0))
c      X2 = dEXP(Z2)
c      XR = dlog(1.D0-X2) - dlog(1.D0-dEXP(ZSOF))
c      if(SLOPE*XR.le.dlog(S_RNDM(0))) goto 200     
c
c      IF(Z1+Z2.LE.ZSOF) GOTO 50
c      STR_mass2 = dsqrt(X1*X2*S)/2.D0
c      PPTT = PPT02(10)
c      IF(IPAR(3).eq.8) PPTT = PPT02(20)
c 150  PT = PPTT*dSQRT(-dLOG(MAX(EPS10,S_RNDM(0))))
c      IF(IPAR(3).ge.6)THEN
c         XM = 0.D0
c         XM2 = XM**2
c         RNDM = MAX(EPS10,S_RNDM(IFL))
c         XMT = PPTT * dLOG(RNDM) - XM
c         XMT2 = XMT**2
c         PT = dSQRT(XMT2-XM2)
c      ENDIF
c      IF(PT.GT.PTmin) GOTO 150
c      IF(PT.GE.STR_mass2) GOTO 150
c      END
c
C=======================================================================

      SUBROUTINE SAMPLE_SOFT6 (STR_mass_min, X1,X2,PT)

C-----------------------------------------------------------------------
C...Routine for the sampling the kinematical variables of sea quarks
C.     according to (1-x)**b / x
C.  INPUT:  STR_mass_min : minimal string mass ** 2 = x1 * x2 * s
C.          SLOPE : large x suppression exponent
C.  OUTPUT:  X1, X2, PT (GeV)                                   /FR'14
C-----------------------------------------------------------------------
Cf2py double precision, intent(in) :: STR_mass_min
Cf2py double precision, intent(out) :: X1
Cf2py double precision, intent(out) :: X2
Cf2py double precision, intent(out) :: PT

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      INTEGER NW_max
      PARAMETER (NW_max = 20)
C--------------------------------------------------------------------
C     SIBYLL common blocks containing event information       \FR'14
C--------------------------------------------------------------------

C     EVENT INFO COMMON
C     contains overall interaction properties, like
C     SQS : center-of-mass energy
C     S   :         "       "     squared
C     PTmin : low pt cut of QCD cross section, 
C             i.e. minimal pt of hard minijets
C     Xmin : low-x bound for PDFs, 
C            i.e. minimal momentum fraction of hard partons
C     Zmin : logarithm of that
C     KB : PID of beam hadron
C     KT() : PID of target
C     IAT : mass number of target
      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
      INTEGER KB,IAT,KT
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)

      DOUBLE PRECISION PPT02
      COMMON /S_CQDIS2/ PPT02(44)

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      SAVE
      
      NOSLOPE = 0
      SLOPE = PAR(42)
      IF(SLOPE.lt.0.5D0) NOSLOPE = 1
      XMAX = 0.8D0
      ZSOF = 2.D0*LOG(STR_mass_min/SQS)       ! minim. mass ~ x1 * x2
      XMINA = MAX(EXP(ZSOF),EPS10)
      AXMINA = 1.D0/XMINA
      IF(ndebug.gt.2)
     &    write(lun,*) ' SAMPLE_SOFT6: Mmin,ZSOF,XMINA,XMAX,SLOPE:',
     &     STR_mass_min,ZSOF,XMINA,XMAX,SLOPE
      
 100  X1 = XM2DIS(XMINA,XMAX,1.D0)            ! ~(1/x)**alpha
      IF(NOSLOPE.eq.1) goto 200
      XRNDM = S_RNDM(0)
      XR = LOG(1.D0-X1)-LOG(1.D0-XMINA)
      IF(ndebug.gt.5)
     &     write(lun,*) '  X1,XR,SLOPE*XR:',X1,XR,SLOPE*XR
      if(SLOPE*XR.le.LOG(max(xrndm,eps10))) goto 100

 200  X2 = XM2DIS(XMINA,XMAX,1.D0)            ! ~(1/x)**alpha
      IF(NOSLOPE.eq.1) goto 300
      XRNDM = S_RNDM(1)
      XR = log(1.D0-X2) - log(1.D0-XMINA)
      IF(ndebug.gt.5)
     &    write(lun,*) '  X2,XR,SLOPE*XR:',X2,XR,SLOPE*XR
      if(SLOPE*XR.le.log(max(xrndm,eps10))) goto 200

 300  Z1 = log(X1)
      Z2 = log(X2)
      IF(Z1+Z2.LE.ZSOF) GOTO 100     
      STR_mass2 = sqrt(X1*X2*S)/2.D0
      PPTT = PPT02(10)
      IF(IPAR(3).eq.8) PPTT = PPT02(20)
      IF(ndebug.gt.2)
     &    write(lun,*) ' SAMPLE_SOFT6: PPTT,Mmin2,PTmin:',
     &PPTT,STR_mass2,PTmin
 150  PT = PPTT*SQRT(-LOG(MAX(EPS10,S_RNDM(0))))
      IF(IPAR(3).ge.6)THEN
         XM = 0.D0
         XM2 = XM**2
         RNDM = MAX(EPS10,S_RNDM(1))
         XMT = PPTT * LOG(RNDM) - XM
         XMT2 = XMT**2
         PT = SQRT(XMT2-XM2)
      ENDIF
      IF(ndebug.gt.2)
     &    write(lun,*) '  XM,XMT2,PT:',XM,XMT2,PT
      IF(PT.GT.PTmin) GOTO 150
      IF(PT.GE.STR_mass2) GOTO 150
      END
C=======================================================================

      SUBROUTINE SIB_START_EV (SQS, L, IA, IAFLG, NW, JDIF)

C-----------------------------------------------------------------------
C...Beginning of a SIBYLL interaction 
C.
C.  add l.m. Glauber SD cross section for pAir  13/FR
C.
C.  INPUT : SQS = c.m.s. energy (GeV)
C.          L = 1:proton, 2:charged pion
C.          IA = mass of target nucleon
C.          IAFLG = target is air
C. 
C.  OUTPUT: NW    = number of wounded nucleons
C.          JDIF(JW)  = diffraction code    !!!! changed to field !!!!
C.                  (0 : non-diffractive interaction)
C.                  (1 : forward diffraction)
C.                  (2 : backward diffraction)
C.                  (3 : double diffraction)
C.
C-----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
c     external type declarations
      INTEGER NW_max,JDIF,IA,L,IAFLG,NW
      DOUBLE PRECISION SQS
      PARAMETER (NW_max = 20)
      DIMENSION JDIF(NW_max)
      
      DOUBLE PRECISION B, BMAX
      INTEGER NTRY, NA
      COMMON /S_CNCM0/ B, BMAX, NTRY, NA
      DOUBLE PRECISION XM2MIN,ALXMIN,SLOP0,ASLOP,BSLOP,XMASS
      COMMON /S_DIFMAss/ XM2MIN(6),ALXMIN(6),SLOP0,ASLOP,BSLOP,XMASS(2)
      DOUBLE PRECISION XI_MAX, ALAM
      COMMON /GLAUB_SCR/ XI_MAX, ALAM(61)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)

c     local type declarations
      DOUBLE PRECISION SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO,
     &     SIGPROD,SIGBDIF,SIGELA,S_RNDM,S,PF,PB,PD,P0,P1,P2,R
      DIMENSION SIGDIF(3)
      INTEGER K
      SAVE

      IF(NDEBUG.gt.0)
     &WRITE(LUN,*)'SIB_START_EV:', SQS, L, IA, IAFLG, NW, JDIF
      
C...sample number of wounded nucleons
c     read hadron-nucleon cross section from table
      CALL SIB_SIGMA_HP(L,SQS,SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO) 

      IF (IA .GT. 1)  THEN
         IF(IPAR(12).NE.0)THEN
            IF(IPAR(12).eq.3)THEN
c     distinguish between nuclear cross sections..
               IF(IAFLG.eq.0)THEN
c     if target is nucleus calc. hadron-nucleus cross section
                  CALL SIB_SIGMA_HNUC(L,IA,SQS,SIGprod,SIGbdif,SIGela)
               ELSE
c     if target is air read hadron-air cross section from table
                  CALL SIB_SIGMA_HAIR(L,SQS,SIGprod,SIGbdif)
               ENDIF
            ELSE
c     always use air cross section...
               CALL SIB_SIGMA_HAIR(L,SQS,SIGprod,SIGbdif)
            ENDIF
C     2channel low-mass (coherent) diffraction?
            IF(S_RNDM(L).LT.SIGbdif/SIGprod)THEN
               NW = 1
               JDIF(1) = 1
               RETURN
            ENDIF
         ENDIF
c     sample number of wounded nucleons
         CALL INT_H_NUC (IA, SIGT, SLOPE, RHO) 
      ELSE
         NA = 1
      ENDIF      
      NW = NA

      IF(NDEBUG.gt.0) 
     &   WRITE(LUN,'(A50,2I3,1P,3E10.3)')
     &   '  START_EVT: IA, NW, SIGT, SLOPE, RHO:',IA,NW,SIGT,SLOPE,RHO
C...new treatment of diffraction 
      IF(IA.GT.1) THEN
c     hadron-nucleus case
         IF(NW.eq.1)THEN
            IF(IPAR(12).NE.0)THEN
c     high mass (incoherent) diffraction?
               S = SQS ** 2
               PF =(1.D0-dLOG(S*XI_MAX/XM2MIN(L))/
     &              dLOG(S*PAR(13)/XM2MIN(L)))*SIGDIF(1)/SIGINEL
               PB = SIGDIF(2)/SIGINEL
               PD = SIGDIF(3)/SIGINEL
            ELSE
               PF = SIGDIF(1)/SIGINEL
               PB = SIGDIF(2)/SIGINEL
               PD = SIGDIF(3)/SIGINEL
            ENDIF
         ELSE
c     Nw>1:
            IF(IPAR(12).EQ.1)THEN
c     all interactions with Nw>1 are non-diff.
               DO K=1, NW
                  JDIF(K) = 0
               ENDDO
               RETURN
            ELSE
c     some Nw>1 are attached by diff. 
               PF = PAR(124)*SIGDIF(1)/SIGINEL
               PB = PAR(124)*SIGDIF(2)/SIGINEL
               PD = PAR(124)*SIGDIF(3)/SIGINEL
            ENDIF
         ENDIF
      ELSE
c     hadron-nucleon case
         PF = SIGDIF(1)/SIGINEL
         PB = SIGDIF(2)/SIGINEL
         PD = SIGDIF(3)/SIGINEL
      ENDIF
      P0 = 1.D0-PF-PB-PD
      P1 = P0 + PF
      P2 = P1 + PB
      DO K=1, NW
         R = S_RNDM(0)
         IF (R .LT. P0)  THEN
            JDIF(K) = 0
         ELSE IF (R .LT. P1)  THEN
            JDIF(K) = 1
         ELSE IF (R .LT. P2)  THEN
            JDIF(K) = 2
         ELSE 
            JDIF(K) = 3
         ENDIF
      ENDDO
      
      END
C=======================================================================

      SUBROUTINE INI_EVENT(ECM,KBEAM,IATARG,IMOD)

C-----------------------------------------------------------------------
C     initializes the stacks and event info common
c     if Imod : 0 - initiate subevent in recursive call
c                  ( keeps the final hadron stack intact )
C             : 1 - initiate entire new event
C-----------------------------------------------------------------------
      IMPLICIT NONE
c     external type declarations
      DOUBLE PRECISION ECM
      INTEGER KBEAM,IATARG,IMOD

c     COMMONs
      INTEGER NW_max
      PARAMETER (NW_max = 20)
      INTEGER NS_max, NH_max
      PARAMETER (NS_max = 20, NH_max = 80)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
C     The final particle output is contained in COMMON /S_PLIST/    
C     NP           : number of final particles
C     P(1:NP, 1:5) : 4-momenta + masses of the final particles 
C     LLIST (1:NP) : codes of final particles
      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP
C     parameters that represent: NW: max. number of wounded nucleons,
C     NS,NH: max. number of soft and hard interactions
c      PARAMETER (NW_max = 20)
C     The COMMON block /S_CHIST/ contains information about the
C     the structure of the  generated event:
C     NWD   = number of wounded nucleons
C     NJET = total number of hard interactions
C     NSOF = total number of soft interactions
C     NNSOF (1:NW) = number of soft pomeron cuts in each interaction
C     NNJET (1:NW) = number of minijets produced in each interaction 
C     JDIF(1:NW) = diffraction code 
C                  0 : non-diff,
C                  1 : beam-diff
C                  2 : target-diff
C                  3 : double-diff
      INTEGER NNSOF,NNJET,JDIF,NWD,NJET,NSOF
      COMMON /S_CHIST/ NNSOF(NW_max),NNJET(NW_max),
     &     JDIF(NW_max),NWD,NJET,NSOF

      INTEGER IBMRDX,ITGRDX,IHMJDX,ISMJDX,ICSTDX,IINTDX
      COMMON /S_INDX/ IBMRDX(3),ITGRDX(NW_max,3),
     &     IHMJDX(NW_max*NH_max),IINTDX(NW_max),
     &     ISMJDX(NW_max*NS_max),ICSTDX(2*NW_max,3)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
C--------------------------------------------------------------------
C     SIBYLL common blocks containing event information       \FR'14
C--------------------------------------------------------------------

C     EVENT INFO COMMON
C     contains overall interaction properties, like
C     SQS : center-of-mass energy
C     S   :         "       "     squared
C     PTmin : low pt cut of QCD cross section, 
C             i.e. minimal pt of hard minijets
C     Xmin : low-x bound for PDFs, 
C            i.e. minimal momentum fraction of hard partons
C     Zmin : logarithm of that
C     KB : PID of beam hadron
C     KT() : PID of target
C     IAT : mass number of target
      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
      INTEGER KB,IAT,KT
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT
      INTEGER II2,JJ2
      DOUBLE PRECISION U2,C2,CD2,CM2
      COMMON /SIB_RAND/ U2(97),C2,CD2,CM2,II2,JJ2

c     local types
      DOUBLE PRECISION PZ,E1,PAWT,S_RNDM,R,FOX
      INTEGER KK,JJ,II,KBA,IREFout,JN
      SAVE
      DATA FOX /0.21522D0/  !atomic percentage of 'non-nitrogen' in air
 
      
      IF(NDEBUG.gt.0.and.IMOD.eq.1) 
     &     WRITE(LUN,'(A50,F10.2,I4,I3,I3)')
     &     '  INI_EVENT: called with (ECM,KBEAM,IATARG,NCALL):',
     &     ECM,KBEAM,IATARG,NCALL

c     set final particle stack to zero
      IF(IMOD.eq.1)then
         NP = 0
         NWD = 0
         NJET = 0
         NSOF = 0
      endif

      CALL INI_PRTN_STCK(0,0)

c     clear index cache
      do kk=1,3
         IBMRDX(kk) = 0
      ENDDO
      do jj=1,NW_max
         do kk=1,3
            ICSTDX(jj,kk) = 0
            ICSTDX(jj+1,kk) = 0
            ITGRDX(jj,kk) = 0
            IINTDX(jj) = 0
         ENDDO
         do ii=1,NH_max
            IHMJDX(NH_max*(JJ-1)+II) = 0
         enddo
         do ii=1,NS_max
            ISMJDX(NS_max*(JJ-1)+II) = 0
         enddo
      ENDDO

      SQS   = Ecm
      S     = SQS*SQS
      
      KB = KBEAM
      KBA = IABS(KBEAM)
c     add beam particles to parton stack, lvl -2
      PZ = PAWT(SQS,AM(KBA),AM(13))
      E1 = SQRT(PZ**2+AM2(KBA))
      CALL ADD_PRTN(0.D0,0.D0,PZ,E1,AM(KBA),KB,-2,0,IREFout)
      IF(IMOD.eq.1)THEN
         IAT = IATARG
         IF(IATARG.EQ.1)THEN
            KT(1) = 13
         ELSE
            IF(IATARG.eq.0)THEN
C...  Generate an 'air' interaction by choosing Nitrogen or Oxygen
               R = S_RNDM(0)
               IATARG = 14
               IF (R .LT. FOX)  IATARG = 16
               if (NDEBUG.gt.0) 
     *           WRITE(lun,*)'fox,rndm,iatarg,eps:',fox,r,iatarg,eps8
            ENDIF
            DO JN=1,IATARG
c     for nuclear target: proton (13) or neutron (14)
               KT(JN) = 13 + INT((2.D0-EPS8)*S_RNDM(JN))
            ENDDO
         ENDIF
      ELSE
         KT(1) = IATARG
      ENDIF

C...energy-dependent transverse momentum cutoff
c...EJA correction 2007.03.27
      IF(IPAR(27).eq.1)THEN
         PTmin = PAR(10)+PAR(11)*EXP(PAR(12)*SQRT(LOG(SQS)))
      else
         PTmin = PAR(10)+PAR(11)*EXP(PAR(12)*SQRT(LOG(S)))
      endif
      XMIN = 4.D0*PTmin**2/S
      ZMIN = LOG(XMIN)
      IF(ndebug.gt.0)then
         write(lun,*) ' INI_EVENT: ncall:', ncall
         write(lun,'(2X,A33,F10.2,1X,F16.2,F8.5,E10.3,F10.5)')
     &        'INI_EVENT: (SQS,S,PTmin,Xmin,Zmin)',
     &        SQS,S,PTmin,Xmin,Zmin
         write(lun,*) ' INI_EVENT: KB,IAT,IATARG,KT',KB,IAT,IATARG
         write(lun,*) '         ',(KT(jj),jj=1,IATARG)
      endif

      CALL PTSETUP_4FLV(ECM)

      return
      END
C-----------------------------------------------------------------------
C     parton level administration tools for SIBYLL                \FR'14
C-----------------------------------------------------------------------

C...  COMMON /S_PRTNS/ : parton stack
c     PP: 4momentum of parton, px,py,pz,energy,mass
c     LPID(1): parton id, i.e. flavor (u:1,d:2,s:3,c:4) for quarks
c     LPID(2): level of parton
c              fragmenting systems (strings,remnants) are marked as level0
c              partons that make up these systems are marked as level1
c     LPID(3): 'downward' reference 
c               pointer from level1 partons to their level0 parent
c     LPID(4): 'upward' reference 
c               pointer from level0 partons to their level-1 parent
c     LVL0IDX: index cache for level0 partons
c     NPP: total number of partons on stack
c     NPP0: number of level0 partons on stack

C=======================================================================

      SUBROUTINE ADD_PRTN(PX,PY,PZ,E,XMS,IPID,LVL,IREFin,IREFout)

C-----------------------------------------------------------------------
C     routine to add a parton to the stack \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      PARAMETER (NPP_max = 1000, NPP0_max = 500)
      COMMON /S_PRTNS/ PP(NPP_max,5), LPID(NPP_max,4), LVL0IDX(NPP0_max)
     &     ,NPP,NPP0
      SAVE

      PP(NPP+1,1) = PX
      PP(NPP+1,2) = PY
      PP(NPP+1,3) = PZ
      PP(NPP+1,4) = E
      PP(NPP+1,5) = XMS
      LPID(NPP+1,1) = IPID
      LPID(NPP+1,2) = LVL     
      LPID(NPP+1,3) = IREFin
      NPP = NPP + 1
c     level0 index
      IF(LVL.eq.0)THEN
         LVL0IDX(NPP0+1) = NPP
         NPP0 = NPP0 + 1
      ENDIF
      IREFout = NPP
      IF(NDEBUG.gt.6)THEN
         WRITE(LUN,*) ' ADD_PRTN: (#,PID,LEVEL,REF)',
     &        NPP,LPID(NPP,1),LPID(NPP,2),LPID(NPP,3)
         WRITE(LUN,*) '  4momentum:        ',(PP(NPP,JJ),JJ=1,5)
      ENDIF
      END

C=======================================================================

      SUBROUTINE ADD_PRTN_4VEC(PIN,IPID,LVL,IREFin,IREFout)

C-----------------------------------------------------------------------
C     wrapper for ADD_PRTN to add 4momentum directly \FR'14
C----------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      DIMENSION PIN(5)
      SAVE

      CALL ADD_PRTN
     &     (PIN(1),PIN(2),PIN(3),PIN(4),PIN(5),IPID,LVL,IREFin,IRF)
      IREFout = IRF
      END

C=======================================================================

      SUBROUTINE ADD_REF(IDX,Irefin)

C-----------------------------------------------------------------------
C     routine to add a reference label to a particle
C     after it has been added to the stack        \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      PARAMETER (NPP_max = 1000, NPP0_max = 500)
      COMMON /S_PRTNS/ PP(NPP_max,5), LPID(NPP_max,4), LVL0IDX(NPP0_max)
     &     ,NPP,NPP0
      SAVE
      
c      IF(LPID(IDX,3).ne.0)  WRITE(LUN,*)
c     &     ' ADD_REF: warning particle already has defined reference,',
c     &     IDX,' overwritting..'
      IF(NDEBUG.gt.6)
     &WRITE(LUN,*) ' ADD_REF: (IDX,REFin)',IDX,Irefin
      LPID(IDX,3) = Irefin
      END

C=======================================================================

      SUBROUTINE RD_REF(IDX,Irefout)

C-----------------------------------------------------------------------
C     routine to add a reference label to a particle
C     after it has been added to the stack        \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      PARAMETER (NPP_max = 1000, NPP0_max = 500)
      COMMON /S_PRTNS/ PP(NPP_max,5), LPID(NPP_max,4), LVL0IDX(NPP0_max)
     &     ,NPP,NPP0
      SAVE

      Irefout = LPID(IDX,3)
      IF(NDEBUG.gt.6)
     &  WRITE(LUN,*) ' RD_ref: (IDX,REFout)',IDX,Irefout
      END

C=======================================================================

      SUBROUTINE ADD_INT_REF(IDX,Irefin)

C-----------------------------------------------------------------------
C     routine to add a reference label to an interaction
C     after it has been added to the stack        \FR'15
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      PARAMETER (NPP_max = 1000, NPP0_max = 500)
      COMMON /S_PRTNS/ PP(NPP_max,5), LPID(NPP_max,4), LVL0IDX(NPP0_max)
     &     ,NPP,NPP0
      SAVE

      IF(NDEBUG.gt.6)
     &    WRITE(LUN,*) ' ADD_INT_REF: (IDX,REFin)',IDX,Irefin
      LPID(IDX,4) = Irefin
      END

C=======================================================================

      SUBROUTINE RD_INT(IDX,Irefout,Iout)

C-----------------------------------------------------------------------
C     routine to add a reference label to an interaction
C     after it has been added to the stack        \FR'15
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      PARAMETER (NPP_max = 1000, NPP0_max = 500)
      COMMON /S_PRTNS/ PP(NPP_max,5), LPID(NPP_max,4), LVL0IDX(NPP0_max)
     &     ,NPP,NPP0
      SAVE

      Irefout = LPID(IDX,4)
      IF(Irefout.ne.0) Iout = LPID(Irefout,1)
      IF(NDEBUG.gt.6)
     &  WRITE(LUN,*) ' RD_INT: (IDX,REFout,Iint)',IDX,Irefout,Iout
      END

C=======================================================================

      SUBROUTINE EDT_PRTN(IDX,PX,PY,PZ,EN,XMS,IREFout)

C-----------------------------------------------------------------------
C     routine to edit a parton already on stack   \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      PARAMETER (NPP_max = 1000, NPP0_max = 500)
      COMMON /S_PRTNS/ PP(NPP_max,5), LPID(NPP_max,4), LVL0IDX(NPP0_max)
     &     ,NPP,NPP0
      SAVE

      IF(NDEBUG.gt.6)THEN
         WRITE(LUN,*) ' EDT_PRTN: (#,PID,LEVEL,REF)',
     &        IDX,LPID(IDX,1),LPID(IDX,2),LPID(IDX,3)
         WRITE(LUN,*) '  initial 4momentum:',(PP(IDX,JJ),JJ=1,5)
      ENDIF
      PP(IDX,1) = PX
      PP(IDX,2) = PY
      PP(IDX,3) = PZ
      PP(IDX,4) = EN
      PP(IDX,5) = XMS
c     return reference to other partons
      IREFout = LPID(IDX,3)
      IF(NDEBUG.gt.6)
     &     WRITE(LUN,*) '  final 4momentum:  ',(PP(IDX,JJ),JJ=1,5)
      END

C=======================================================================

      SUBROUTINE RD_PRTN(IDX,PX,PY,PZ,EN,XMS,IFL,IREFout)

C-----------------------------------------------------------------------
C     routine to read a parton from the stack     \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      PARAMETER (NPP_max = 1000, NPP0_max = 500)
      COMMON /S_PRTNS/ PP(NPP_max,5), LPID(NPP_max,4), LVL0IDX(NPP0_max)
     &     ,NPP,NPP0
      SAVE

      IF(NDEBUG.gt.6)THEN
         WRITE(LUN,*) ' RD_PRTN: (#,PID,LEVEL,REF)',
     &        IDX,LPID(IDX,1),LPID(IDX,2),LPID(IDX,3)
         WRITE(LUN,*) '  4momentum:        ',(PP(IDX,JJ),JJ=1,5)
      ENDIF
      PX = PP(IDX,1)
      PY = PP(IDX,2)
      PZ = PP(IDX,3)
      EN = PP(IDX,4)
      XMS = PP(IDX,5)
      IFL = LPID(IDX,1)
c     return reference to other partons
      IREFout = LPID(IDX,3)
      END

C=======================================================================

      SUBROUTINE RD_PRTN_4VEC(IDX,Pin,IFL,IREFout)

C-----------------------------------------------------------------------
C     routine to read a parton from the stack     \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      PARAMETER (NPP_max = 1000, NPP0_max = 500)
      COMMON /S_PRTNS/ PP(NPP_max,5), LPID(NPP_max,4), LVL0IDX(NPP0_max)
     &     ,NPP,NPP0

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      DIMENSION Pin(5)
      SAVE

      IF(IDX.EQ.0) THEN
         WRITE(LUN,*) ' RD_PRTN_4VEC: invalid index!',IDX
         xa = -1.D0
         xa = log(xa)
         RETURN
      ELSE         
         do ii = 1,5
            PIN(ii) = PP(IDX,ii)
         enddo
         IFL = LPID(IDX,1)
c     return reference to other partons
         IREFout = LPID(IDX,3)
         IF(NDEBUG.gt.6)THEN
            WRITE(LUN,*) ' RD_PRTN: (#,PID,LEVEL,REF)',
     &           IDX,IFL,LPID(IDX,2),IREFout
            WRITE(LUN,*) '  4momentum:        ',(PIN(JJ),JJ=1,5)
         ENDIF

      ENDIF
      END

C=======================================================================

      SUBROUTINE ITR_LVL0_PRTN(JJ,IDX,LID)

C-----------------------------------------------------------------------
C     routine that serves as iterator over the level0
C     partons on the stack                        \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      PARAMETER (NPP_max = 1000, NPP0_max = 500)
      COMMON /S_PRTNS/ PP(NPP_max,5), LPID(NPP_max,4), LVL0IDX(NPP0_max)
     &     ,NPP,NPP0
      SAVE
 
      IDX = LVL0IDX(JJ)
      IF(ndebug.gt.6)
     &  WRITE(LUN,*) ' ITR_LVL0_PRTN: JJ,IDX',JJ,IDX
      LID = LPID(IDX,1)
      IF(JJ+1.gt.NPP0) THEN
         JJ = -1
         RETURN
      ELSE
         JJ = JJ + 1
      ENDIF      
      END

C=======================================================================

      SUBROUTINE INI_PRTN_STCK(NOLD,N0OLD)

C-----------------------------------------------------------------------
C     reset parton stack                          \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      PARAMETER (NPP_max = 1000, NPP0_max = 500)
      COMMON /S_PRTNS/ PP(NPP_max,5), LPID(NPP_max,4), LVL0IDX(NPP0_max)
     &     ,NPP,NPP0
      SAVE

      IF(NDEBUG.gt.6) WRITE(LUN,*) ' PRTN_STCK: reset .. '
      IF(NDEBUG.gt.6) WRITE(LUN,*) '  old state: NPP,NPP0',NPP,NPP0
      
      NPP = NOLD
      NPP0 = N0OLD

      IF(NDEBUG.gt.6) WRITE(LUN,*) '  new state: NPP,NPP0',NPP,NPP0

      END

C=======================================================================

      SUBROUTINE GET_NPP(NPPLD,NPP0LD)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      PARAMETER (NPP_max = 1000, NPP0_max = 500)
      COMMON /S_PRTNS/ PP(NPP_max,5), LPID(NPP_max,4), LVL0IDX(NPP0_max)
     &     ,NPP,NPP0
      SAVE

      NPPLD = NPP
      NPP0LD = NPP0
      END

C=======================================================================

      SUBROUTINE GET_LVL0(NPP0LD,IDXLIST)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      PARAMETER (NPP_max = 1000, NPP0_max = 500)
      COMMON /S_PRTNS/ PP(NPP_max,5), LPID(NPP_max,4), LVL0IDX(NPP0_max)
     &     ,NPP,NPP0
      DIMENSION IDXLIST(NPP0_max)
      INTEGER   N
      SAVE

      NPP0LD = NPP0
      DO N = 1, NPP0_max
        IDXLIST(N) = LVL0IDX(N)
      ENDDO

      END

C=======================================================================

      SUBROUTINE PRNT_PRTN_STCK

C-----------------------------------------------------------------------
C     as the name suggests, prints the current state
C     of the parton stack                         
C     print unit is defined in S_DEBUG:LUN        \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      PARAMETER (NPP_max = 1000, NPP0_max = 500)
      COMMON /S_PRTNS/ PP(NPP_max,5), LPID(NPP_max,4), LVL0IDX(NPP0_max)
     &     ,NPP,NPP0
      CHARACTER*5 CDE
      CHARACTER*9 CODE

      CHARACTER*6 NAMP
      COMMON /S_CNAM/ NAMP (0:99)
      SAVE

      WRITE (LUN,50) 
 50   FORMAT(3X,88('-'),/,21X,'SIBYLL PARTON LEVEL EVENT SUMMARY',21X,
     &     /,3X,75('-'),13('-'))

c     beam particles
      WRITE(LUN,*) '  BEAM PARTICLES'
 52   FORMAT(4X,'#',3X,'PID',2x,'LVL',2x,'REF',20x,'PX',9x,'PY',7x,
     +     'PZ',9x,'E',11X,'Mass', /, 3X,75('-'),13('-'))
      WRITE (LUN,52)
      DO J=1,NPP
         IF(LPID(J,2).eq.-2)then
            WRITE (LUN,60) J, (LPID(J,KK),KK=1,3), (PP(J,K),K=1,5)
         ENDIF
      ENDDO
c     level -2 format
 60   FORMAT(4I5,14X,2F11.3,1p,2E11.3,0p,F9.3)
      WRITE(LUN,61)
 61   FORMAT(3X,75('-'),13('-'))

c     interactions
      WRITE(LUN,*) '  INTERACTIONS'
 62   FORMAT(4X,'#',3X,'PID',2x,'LVL',2x,'REF',20x,'NSOF',8x,'NJET',7x,
     +     'JDIF',7x,'E',11X,'Mass', /, 3X,75('-'),13('-'))
      WRITE (LUN,62)
      DO J=1,NPP
         IF(LPID(J,2).eq.-1)then
            WRITE (LUN,63) J, (LPID(J,KK),KK=1,3), (PP(J,K),K=1,5)
         ENDIF
      ENDDO
c     level -1 format
 63   FORMAT(4I5,12X,4F12.0,F11.3)
 64   FORMAT(3X,75('-'),13('-'))
      WRITE(LUN,64)
      
c     partons
      WRITE (LUN,100)
      DO J=1,NPP
         IF(LPID(J,2).eq.0)then
            WRITE (LUN,120) J, (LPID(J,KK),KK=1,3), (PP(J,K),K=1,5)
         elseif(LPID(J,2).eq.1)then
            CALL KCODE(LPID(J,1),cde,nc)
            WRITE (LUN,121) J, CDE(1:nc),(LPID(J,KK),KK=2,3), 
     &           (PP(J,K),K=1,5)
         elseif(LPID(J,2).eq.2)then
            CODE = '        '
            L = LPID(J,1)
            CODE(1:6) = NAMP(IABS(L))
            IF (L .LT. 0) CODE(7:9) = 'bar'
            WRITE (LUN,122) J,CODE,(LPID(J,KK),KK=2,3), (PP(J,K),K=1,5)
         endif
      ENDDO
      CALL PPSUM(1,NPP,Esum,PXsum,PYsum,PZsum,NF)
      WRITE(LUN,140) PXsum,PYsum,PZsum,Esum
 
 100  FORMAT(4X,'#',3X,'PID',2x,'LVL',2x,'REF',20x,'PX',9x,'PY',7x,
     +     'PZ',9x,'E',11X,'Mass', /, 3X,75('-'),13('-'))
c     level 0 format
 120  FORMAT(4I5,14X,2F11.3,1p,2E11.3,0p,F11.3)
c     level 1 format cjoe
 121  FORMAT(I7,1X,A5,2I5,14X,2F11.3,1p,2E11.3,0p,F11.3)
c     level 2 format
 122  FORMAT(I10,1X,A9,2I5,10X,2F11.3,1p,2E11.3,0p,F11.3)
 140  FORMAT(3X,75('-'),13('-'),/,'  Tot = ',26X,2F11.3,1p,2e11.3)

      END

C=======================================================================

      SUBROUTINE PPSUM(N1,N2,ETOT,PXT,PYT,PZT,NF)

C-----------------------------------------------------------------------
C     Return the energy,px,py,pz of level0 partons 
C     in the list between N1 and N2
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      PARAMETER (NPP_max = 1000, NPP0_max = 500)
      COMMON /S_PRTNS/ PP(NPP_max,5), LPID(NPP_max,4), LVL0IDX(NPP0_max)
     &     ,NPP,NPP0

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      SAVE

      NF=0
      ETOT=0.D0
      PXT=0.D0
      PYT=0.D0
      PZT=0.D0
      DO J=N1,N2
         IF (LPID(J,2) .EQ. 0)  THEN
           NF = NF+1
           ETOT = ETOT + ABS( PP(J,4) )
           PXT = PXT + PP(J,1)
           PYT = PYT + PP(J,2)
           PZT = PZT + PP(J,3)
         ENDIF
      ENDDO
      RETURN
      END
C=======================================================================

      SUBROUTINE FOUR_LENGTH(XP,XM2)

C-----------------------------------------------------------------------
C     Calculate the length of a 4vector (+---)    \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      DIMENSION XP(5)
      SAVE

      XM2 = XP(4)**2 - XP(1)**2 - XP(2)**2 - XP(3)**2
      END
C=======================================================================

      DOUBLE PRECISION FUNCTION CALC_INVM(XP1,XP2)

C-----------------------------------------------------------------------
C     Calculate the invariant mass of two 4vectors FR'15
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      DIMENSION XP1(5),XP2(5)
      SAVE

      CALC_INVM = (XP1(4)+ XP2(4))**2
      DO I=1,3
         CALC_INVM = CALC_INVM-(XP1(I)+XP2(I))**2         
      ENDDO
      CALC_INVM = SQRT(CALC_INVM)
      END

C=======================================================================
      
      SUBROUTINE GET_XMT2(IDX,XM2)

C-----------------------------------------------------------------------
C     Calculate the transverse mass of a parton 
C     on the stack                                \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      PARAMETER (NPP_max = 1000, NPP0_max = 500)
      COMMON /S_PRTNS/ PP(NPP_max,5), LPID(NPP_max,4), LVL0IDX(NPP0_max)
     &     ,NPP,NPP0
      SAVE

      XM2 = PP(IDX,1)**2 + PP(IDX,2)**2 + PP(IDX,5)**2
      END
C=======================================================================

      SUBROUTINE GET_IMASS2(IDX,XM2)

C-----------------------------------------------------------------------
C     Calculate the invariant mass squared of a parton 
C     on the stack  (+---)                        \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      PARAMETER (NPP_max = 1000, NPP0_max = 500)
      COMMON /S_PRTNS/ PP(NPP_max,5), LPID(NPP_max,4), LVL0IDX(NPP0_max)
     &     ,NPP,NPP0
      SAVE

      XM2 = PP(IDX,1)**2 + PP(IDX,2)**2 + PP(IDX,3)**2
      XM2 = PP(IDX,4)**2 - XM2
      END

C=======================================================================

      SUBROUTINE GET_MASS(IDX,XM)

C-----------------------------------------------------------------------
C     read mass of parton on stack                \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      PARAMETER (NPP_max = 1000, NPP0_max = 500)
      COMMON /S_PRTNS/ PP(NPP_max,5), LPID(NPP_max,4), LVL0IDX(NPP0_max)
     &     ,NPP,NPP0

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      SAVE

      IF(IDX.EQ.0) THEN 
         XM2 = 0.D0
      else
         XM = PP(IDX,5)
      ENDIF
      END
C=======================================================================

      SUBROUTINE GET_MASS2(IDX,XM2)

C-----------------------------------------------------------------------
C     read mass of parton on stack                \FR'14
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      PARAMETER (NPP_max = 1000, NPP0_max = 500)
      COMMON /S_PRTNS/ PP(NPP_max,5), LPID(NPP_max,4), LVL0IDX(NPP0_max)
     &     ,NPP,NPP0

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      SAVE

      IF(IDX.EQ.0) THEN 
         XM2 = 0.D0
      else
         XM2 = PP(IDX,5)**2
      ENDIF
      END

C=======================================================================

      SUBROUTINE GET_VRTLTY(IDX,XX)

C-----------------------------------------------------------------------
C     calculate virtuality of parton on stack     \FR'14
C     = on-shell mass - inv. mass
C-------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      PARAMETER (NPP_max = 1000, NPP0_max = 500)
      COMMON /S_PRTNS/ PP(NPP_max,5), LPID(NPP_max,4), LVL0IDX(NPP0_max)
     &     ,NPP,NPP0
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10
      SAVE

      IF(IDX.EQ.0) XM2 = 0.D0
      CALL GET_IMASS2(IDX,xm2)
      XX = PP(IDX,5)**2-xm2
      END

C=======================================================================

      SUBROUTINE ADD_4VECS(P1,P2,POUT)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      DIMENSION P1(5),P2(5),POUT(5)
      SAVE

      DO II=1,4
         POUT(II) = P1(II) + P2(II)
      ENDDO
      CALL FOUR_LENGTH(POUT,XM2)
      IF(XM2.LT.0)THEN
c     virtual particle
         POUT(5) = -1.D0
         IF(NDEBUG.gt.6)then
            WRITE(LUN,*)
     &           ' ADD_4VECS: resulting particle virtual!! (m**2):',XM2
            WRITE(LUN,*) ' p**2' , POUT(1)**2+POUT(2)**2+POUT(3)**2
            WRITE(LUN,*) ' E**2: ', POUT(4)**2
         ENDIF
      ELSE
         POUT(5) = sqrt(xm2)
      ENDIF
      END
C=======================================================================

      SUBROUTINE DECPAR (LA,P0,ND,LL,P)

C-----------------------------------------------------------------------
C...This subroutine generates the decay of a particle
C.  with ID = LA, and 5-momentum P0(1:5)
C.  into ND particles of 5-momenta P(j,1:5) (j=1:ND)
C.
C.  If the initial particle code is LA=0
C.  then ND and LL(1:ND) are considered as  input and
C.  the routine generates a phase space decay into ND
C.  particles of codes LL(1:nd)
C.
C.  june 1992
C.  This version  contains the decay of polarized muons
C.  The muon codes are  L =  4 : mu+ R
C.                          -4 : mu+ L
C.                           5 : mu- L
C.                          -5 : mu- R
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      DOUBLE PRECISION CBR
      INTEGER KDEC,LBARP,IDB
      COMMON /S_CSYDEC/ CBR(223+16+12+8), KDEC(1338+6*(16+12+8)),
     &     LBARP(99), IDB(99)
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)
      DIMENSION P0(5), LL(10), P(10,5)
      DIMENSION PV(10,5), RORD(10), UE(3),BE(3)

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      SAVE
      
C...Phase space decay into the particles in the list
      IF (LA .EQ. 0)  THEN
          MAT = 0
          MBST = 0
          PS = 0.D0
          DO J=1,ND
CDH          following statements corrected by D.H. dec 20.,1995
             P (J,5) = AM(IABS(LL(J)))
             PV(J,5) = AM(IABS(LL(J)))
             PS = PS+P(J,5)
          ENDDO
          DO J=1,4
             PV(1,J) = P0(J)
          ENDDO
          PV(1,5) = P0(5)
          GOTO 140
      ENDIF
         
C...Choose decay channel
      L = IABS(LA)
      ND=0
      IDC = IDB(L)-1
      IF (IDC+1 .LE.0)  RETURN
      RBR = S_RNDM(0)
110   IDC=IDC+1
      IF(RBR.GT.CBR(IDC))  GOTO 110

      KD =6*(IDC-1)+1
      ND = KDEC(KD)
      MAT= KDEC(KD+1)
      MBST=0
      IF (MAT .GT.0 .AND. P0(4) .GT. 20.D0*P0(5)) MBST=1
      IF (MAT .GT.0 .AND. MBST .EQ. 0) 
     +        BETA = DSQRT(P0(1)**2+P0(2)**2+P0(3)**2)/P0(4)
      PS = 0.D0
c     reduce omega mass by 50MeV to allow on-shell N(1710) decay
      Xmomega = am(32)
      IF(L.eq.53.or.L.eq.54) AM(32) = AM(32)-0.05D0
      DO J=1,ND
         LL(J) = KDEC(KD+1+J)
         P(J,5)  = AM(LL(J))
         PV(J,5) = AM(LL(J))
         PS = PS + P(J,5)
      ENDDO
      AM(32) = Xmomega
      DO J=1,4
         PV(1,J) = 0.D0
         IF (MBST .EQ. 0)  PV(1,J) = P0(J)
      ENDDO
      IF (MBST .EQ. 1)  PV(1,4) = P0(5)
      PV(1,5) = P0(5)

140   IF (ND .EQ. 2) GOTO 280

      IF (ND .EQ. 1)  THEN
         DO J=1,4
            P(1,J) = P0(J)
         ENDDO
         RETURN
      ENDIF

C...Calculate maximum weight for ND-particle decay
      WWTMAX = 1.D0/FACN(ND)      
      PMAX=PV(1,5)-PS+P(ND,5)
      PMIN=0.D0
      DO IL=ND-1,1,-1
         PMAX = PMAX+P(IL,5)
         PMIN = PMIN+P(IL+1,5)
         WWTMAX = WWTMAX*PAWT(PMAX,PMIN,P(IL,5))
      ENDDO

C...generation of the masses, compute weight, if rejected try again
240   RORD(1) = 1.D0
      DO 260 IL1=2,ND-1
        RSAV = S_RNDM(0)
        DO 250 IL2=IL1-1,1,-1
          IF(RSAV.LE.RORD(IL2))   GOTO 260
250     RORD(IL2+1)=RORD(IL2)
260   RORD(IL2+1)=RSAV
      RORD(ND) = 0.D0
      WT = 1.D0      
      DO 270 IL=ND-1,1,-1
        PV(IL,5)=PV(IL+1,5)+P(IL,5)+(RORD(IL)-RORD(IL+1))*(PV(1,5)-PS)
270   WT=WT*PAWT(PV(IL,5),PV(IL+1,5),P(IL,5))
      IF (WT.LT.S_RNDM(1)*WWTMAX)   GOTO 240

C...Perform two particle decays in respective cm frame
280   DO 300 IL=1,ND-1
        PA=PAWT(PV(IL,5),PV(IL+1,5),P(IL,5))
        UE(3)=2.D0*S_RNDM(IL)-1.D0
        PHI=TWOPI*S_RNDM(3)
        UT = DSQRT(1.D0-UE(3)**2)
        UE(1) = UT*dCOS(PHI)
        UE(2) = UT*dSIN(PHI)
        DO 290 J=1,3
          P(IL,J)=PA*UE(J)
290     PV(IL+1,J)=-PA*UE(J)
        P(IL,4)=DSQRT(PA**2+P(IL,5)**2)
300   PV(IL+1,4)=DSQRT(PA**2+PV(IL+1,5)**2)

C...Lorentz transform decay products to lab frame
      DO 310 J=1,4
310   P(ND,J)=PV(ND,J)
      DO 340 IL=ND-1,1,-1
        DO 320 J=1,3
320     BE(J)=PV(IL,J)/PV(IL,4)
        GA=PV(IL,4)/PV(IL,5)
        DO 340 I=IL,ND
          BEP = BE(1)*P(I,1)+BE(2)*P(I,2)+BE(3)*P(I,3)
          DO 330 J=1,3
330       P(I,J)=P(I,J)+GA*(GA*BEP/(1.D0+GA)+P(I,4))*BE(J)
340   P(I,4)=GA*(P(I,4)+BEP)
      
C...Weak decays
      IF (MAT .EQ. 1)  THEN
         F1=P(2,4)*P(3,4)-P(2,1)*P(3,1)-P(2,2)*P(3,2)-P(2,3)*P(3,3)     
         IF (MBST.EQ.1)  THEN
C           WT = P0(5)*P(1,4)*F1
            WT = P0(5)*(P(1,4)+DBLE(LA/L)*P(1,3))*F1
         ENDIF
         IF (MBST.EQ.0)  THEN  
            WT=F1*(P(1,4)*P0(4)-P(1,1)*P0(1)-P(1,2)*P0(2)-P(1,3)*P0(3))
            IF(L.lt.50)
     +           WT= WT-DBLE(LA/L)*(P0(4)*BETA*P(1,4)-P0(4)*P(1,3))*F1
         ENDIF
         WTMAX = P0(5)**4/8.D0
         IF(WT.LT.S_RNDM(0)*WTMAX)   GOTO 240
      ENDIF

C...Boost back for rapidly moving particle
      IF (MBST .EQ. 1)   THEN
         DO 440 J=1,3
440      BE(J)=P0(J)/P0(4)
         GA= P0(4)/P0(5)
         DO 460 I=1,ND
           BEP=BE(1)*P(I,1)+BE(2)*P(I,2)+BE(3)*P(I,3)
           DO 450 J=1,3
450        P(I,J)=P(I,J)+GA*(GA*BEP/(1.D0+GA)+P(I,4))*BE(J)
460      P(I,4)=GA*(P(I,4)+BEP)
      ENDIF

C...labels for antiparticle decay
      IF (LA .LT. 0 .AND. L .GT. 18)  THEN
           DO J=1,ND
            LL(J) = LBARP(LL(J))
         ENDDO
      ENDIF

      RETURN
      END

C=======================================================================

      BLOCK DATA DATDEC

C-----------------------------------------------------------------------
C...initialization of SIBYLL particle data
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      DOUBLE PRECISION CBR
      INTEGER KDEC,LBARP,IDB
      COMMON /S_CSYDEC/ CBR(223+16+12+8), KDEC(1338+6*(16+12+8)),
     &     LBARP(99), IDB(99)
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)

      DOUBLE PRECISION AW,AW2
      COMMON /S_WIDTH1/ AW(99), AW2(99)

      INTEGER ICHP,ISTR,IBAR
      COMMON /S_CHP/ ICHP(99), ISTR(99), IBAR(99)

      INTEGER IISO,ISPN
      COMMON /S_SPN/ IISO(99), ISPN(99)

      INTEGER ICHM
      COMMON /S_CHM/ ICHM(99)

      CHARACTER*6 NAMP
      COMMON /S_CNAM/ NAMP (0:99)
      SAVE
c     CBR contains the normed sum of the branching ratios of the decay channels
c     indexed by IDB, i.e. a particle with 4 decay channels will have the entries
c     [B1/Btot, (B1+B2)/Btot, (B1+B2+B3)/Btot, 1.]
      DATA CBR /3*1.D0,0.D0,1.D0,1.D0,0.6354D0,0.8422D0,0.8981D0,
     + 0.9157D0,0.9492D0,1.D0,0.6354D0,0.8422D0,0.8981D0,0.9157D0,
     + 0.9492D0,1.D0,0.1965D0,0.3224D0,0.4579D0,0.5934D0,0.7967D0,1.D0,
     + 0.6925D0,1.D0,3*0.D0,0.5D0,1.D0,0.5D0,1.D0,
     + 0.3941D0,0.7197D0,0.9470D0,0.9930D0,1.D0,                     ! eta
     + 0.4285D0,0.7193D0,0.9487D0,0.9750D0,0.9973D0,0.9999D0,1.D0,   ! eta'
     + 3*1.D0,                                                       ! rho-mesons
     + 0.6670D0,1.D0,                                                ! K*+
     + 0.4894D0,0.8317D0,0.9850D0,0.9981D0,0.9994D0,0.9997D0,1.D0,   ! phi(1020)
     + 2*0.D0,                                                       ! (empty)      
     + 0.6670D0,1.D0,                                                ! K*-
     + 0.6670D0,1.D0,                                                ! K*0
     + 0.6670D0,1.D0,                                                ! K*0 bar
     + 0.8940D0,0.9830D0,1.D0,                                       ! omega
     + 4*0.D0,                                                       ! (empty)
     + 0.5160D0,5*1.D0,0.6410D0,2*1.D0,0.67D0,1.D0,0.33D0,2*1.D0,
     + 0.88D0,0.94D0,1.D0,0.88D0,0.94D0,1.D0,0.88D0,0.94D0,1.D0,0.33D0,
     + 1.D0,0.67D0,1.D0,0.678D0,0.914D0,1.D0,0.217D0,0.398D0,0.506D0,
     + 0.595D0,0.684D0,0.768D0,0.852D0,0.923D0,0.976D0,1.D0,0.217D0,
     + 0.398D0,0.506D0,0.595D0,0.684D0,0.768D0,0.852D0,0.923D0,0.976D0,
     + 1.D0,0.2490D0,0.4604D0,0.5338D0,0.5703D0,0.7440D0,0.7840D0,
     + 0.8460D0,0.8880D0,0.9230D0,0.9650D0,1.D0,0.2490D0,0.4604D0,
     + 0.5338D0,0.5703D0,0.7440D0,0.7840D0,0.8460D0,0.8880D0,0.9230D0,
     + 0.9650D0,1.D0,0.1666D0,0.3332D0,0.4998D0,0.6664D0,0.8330D0,1.D0,
     + 0.6770D0,0.9840D0,1.D0,
     + 0.6770D0,0.9840D0,1.D0,0.6190D0,1.D0,0.6190D0,1.D0,0.0602D0,
     + 0.1203D0,1.D0,3*1.D0,0.06D0,0.08D0,0.14D0,0.16D0,0.73D0,0.855D0,
     + 0.98D0,1.D0,0.08D0,0.16D0,0.92D0,1.D0,0.2335D0,0.4283D0,0.6446D0,
     + 0.7099D0,0.8080D0,0.9080D0,0.9380D0,0.9540D0,0.9840D0,1.D0,
     + 3*1.D0,0.5D0,1.D0,0.5D0,1.D0,0.08D0,0.16D0,0.92D0,1.D0,0.942D0,
     + 1.D0,0.942D0,1.D0,0.2493D0,0.4061D0,0.5602D0,0.6860D0,0.7608D0,
     + 0.8305D0,0.8818D0,0.9277D0,0.9691D0,1.D0,0.2493D0,0.4061D0,
     + 0.5602D0,0.6860D0,0.7608D0,0.8305D0,0.8818D0,0.9277D0,0.9691D0,
     + 1.D0,
     & 0.466D0,0.7D0,0.899D0,1.D0,0.466D0,0.7D0,0.899D0,1.D0, ! N1440+-
     & 0.3334D0,0.5D0,0.6334D0,0.7634D0,0.8734D0,0.9394D0,1.D0, ! N1710+
     & 0.3334D0,0.5D0,0.6334D0,0.7634D0,0.8734D0,0.9394D0,1.D0, ! N1710-
     & 0.5D0, 1.D0, 0.5D0, 1.D0, 0.5D0, 1.0D0, ! pi1+-0      
     & 0.6666D0,1.D0, 0.6666D0,1.D0,0.6666D0,1.D0,0.6666D0,1.D0/ ! K0*
      DATA AM / 0.0,2*0.511D-3, 2*0.10566, 0.13497, 2*0.13957,
     +   2*0.49368, 2*0.49761, 0.93827, 0.93957, 4*0.0,0.93827,
     +   0.93957, 2*0.49761, 0.54785,0.95766,2*0.76690,0.76850,
     +   2*0.89166D0,2*0.89600,0.78265,1.01946D0,1.18937D0,1.19264D0,
     +   1.19745,1.31486,1.32171,1.11568,1.23100,1.23500,
     +   1.23400,1.23300,1.38280,1.38370,1.38720,
     +   1.53180,1.53500,1.67245,0.,1.44,1.44,1.71,1.71,4*0.0,
     +   2*1.86926,1.30,1.30,1.30,4*1.430, 3*0.0,      
     +   2*1.86484,2.9803,2*1.9685,2*2.1123,2*2.01027,2*2.00697,
     +   0.0,3.09692,2.45402,2.4529,2.45376,2.4679,2.4710,
     +   2.28646, 2*1.777, 2*0.0, 2.5184,2.5175, 2.5180, 2.6466,
     +   2.6461, 2.6975 /
      DATA AM2 /0.0,2*2.61121D-07,2*0.011164,0.018217,0.019480,
     + 0.019480,0.243720,0.243720,0.247616,0.247616,0.880351,
     + 0.882792,0.000000,0.000000,0.000000,0.000000,0.880351,
     + 0.882792,0.247616,0.247616,0.300140,0.917113,0.588136,
     + 0.588136,0.590592,0.795058,0.795058,0.802816,0.802816,
     + 0.612541,1.039299,1.414601,1.422390,1.433887,1.728857,
     + 1.746917,1.244742,1.515361,1.525225,1.522765,1.520289,
     + 1.912136,1.914626,1.924324,2.346411,2.356225,2.797022,
     + 0.,2.0736,2.0736,2.9241,2.9241,4*0.0, 2*3.49414,
     + 1.690, 1.690, 1.690, 4*2.0449, 3*0.0, 2*3.477628, 8.882188,       
     + 2*3.8750,2*4.4618,2*4.041186,2*4.027928, 0.0, 9.590914, 6.022214,
     + 6.016718, 6.020938,6.09053, 6.105841, 5.227899, 2*3.158, 2*0.0,
     + 6.342339, 6.337806, 6.340323,7.004492, 7.001845, 7.276506/
      DATA AW /24*0.D0,0.022231D0,0.022231D0,0.022231D0,0.002581D0,
     &     0.002581D0,0.D0,0.D0,7.20801D-05,1.81476D-05,6*0.D0,
     &     0.013689D0,0.013689D0,0.013689D0,0.013689D0,0.001296D0,
     &     0.001295D0,0.00155D0,8.281D-05,9.801D-05,0.D0,0.D0,0.09D0,
     &     0.01D0,0.09D0,0.01D0,6*0.D0,0.1D0,0.1D0,0.1D0,4*0.27D0,
     &     32*0.D0/
      DATA AW2 /24*0.D0,0.022231D0,0.022231D0,0.022231D0,0.002581D0,
     &     0.002581D0,0.D0,0.D0,7.20801D-05,1.81476D-05,6*0.D0,
     &     0.013689D0,0.013689D0,0.013689D0,0.013689D0,0.001296D0,
     &     0.001295D0,0.00155D0,8.281D-05,9.801D-05,0.D0,0.D0,0.09D0,
     &     0.01D0,0.09D0,0.01D0,6*0.D0,0.01D0,0.01D0,0.01D0,4*0.0729D0,
     &     32*0.D0/
c     IDB is the index to the branching ratios (CBR) and decay channels (KDEC).
c     always indicates the first decay channel
      DATA IDB /
     +     0,0,0,1,2,                                        ! leptons
     +     3,5,6,7,13,19,25,                                 ! pions and kaons
     +     8*0,30,32,34,39,46,47,48,49,60,62,64,66,51, !69,       ! meson resonances
     +     73,75,76,77,78,79,81,82,84,86,87,90,93,96,98,100, ! baryons : Sibyll 2.1
     +     0,224,228,232,239,4*0,                            ! Nucleon resonaces
     +     103,113,246,248,250, 252,254,256,258,3*0,      
     +     123,134,145,204,214,200,202,151,154,157,159,0,
     +     161,164,165,166,167,175,179,4*0,189,190,191,192,194,196 /
c     KDEC contains decay channels, format is [ND, MAT, LL(1:4)]
c     where ND is the number of particles in the final state (max 4)
C     MAT is 0, 1 for semi-leptonic (weak decay) or not
c     (adds primitive matrix element)
c     LL(1:4) are the particle ids of the final state particles
      DATA KDEC /
     + 3,1,15,2,18,0,3,1,16,3,17,0,2,0,1,1,8*0,2,0,4,17,0,0,2,0,5,18,0,
     + 0,2,0,4,17,0,0,2,0,7,6,0,0,3,0,7,7,8,0,3,0,7,6,6,0,3,1,17,4,6,0,
     + 3,1,15,2,6,0,2,0,5,18,0,0,2,0,8,6,0,0,3,0,8,8,7,0,3,0,8,6,6,0,3,
     + 1,18,5,6,0,3,1,16,3,6,0,3,0,6,6,6,0,3,0,7,8,6,0,3,1,18,5,7,0,3,
     + 1,17,4,8,0,3,1,16,3,7,0,3,1,15,2,8,0,2,0,7,8,0,0,2,0,6,6,20*0,1,
     + 0,11,3*0,1,0,12,0,0,0,1,0,11,0,0,0,1,0,12,0,0,0,2,0,1,1,0,0,3,0,
     + 6,6,6,0,3,0,7,8,6,0,3,0,1,7,8,0,3,0,1,3,2,0,
     + 3,0,7,8,23,0, 3,0,6,6,23,0, 2,0,1,27,0,0, 2,0,1,32,0,0,           ! eta'
     + 2,0,1,1,0,0, 3,0,6,6,6,0, 3,0,1,4,5,0,                            ! eta'
     + 2,0,7,6,0,0,                                                      ! rho+
     + 2,0,8,6,0,0,                                                      ! rho-
     + 2,0,7,8,0,0,                                                      ! rho0
     + 2,0,21,7,0,0, 2,0,9,6,0,0,                                        ! K*+
     + 2,0,9,10,0,0, 2,0,11,12,0,0, 3,0,7,8,6,0, 2,0,1,23,0,0,           ! phi(1020)
     + 2,0,1,6,0,0, 2,0,2,3,0,0, 2,0,4,5,0,0,                            ! phi(1020)                  
     + 12*0,
     + 2,0,22,8,0,0, 2,0,10,6,0,0,                                       ! K*-
     + 2,0,9,8,0,0, 2,0,21,6,0,0,                                        ! K*0
     + 2,0,10,7,0,0, 2,0,22,6,0,0,                                       ! K*0 bar
     + 3,0,7,8,6,0, 2,0,1,6,0,0, 2,0,7,8,0,0,                            ! omega
     + 24*0,
     + 2,0,13,6,0,0,2,0,14,7,0,0,2,0,39,1,0,0,2,                         ! baryons
     + 0,14,8,0,0,2,0,39,6,0,0,2,0,39,8,0,0,2,0,13,8,0,0,2,0,
     + 14,6,0,0,2,0,13,7,0,0,2,0,13,6,
     + 0,0,2,0,14,7,0,0,2,0,13,8,0,0,2,0,14,6,0,0,2,0,14,8,0,0,2,0,
     + 39,7,0,0,2,0,34,6,0,0,2,0,35,7,0,0,2,0,39,6,0,0,2,0,34,8,0,0,
     + 2,0,36,7,0,0,2,0,39,8,0,0,2,
     + 0,35,8,0,0,2,0,36,6,0,0,2,0,37,6,0,0,2,0,38,7,0,0,2,0,
     + 37,8,0,0,2,0,38,6,0,0,2,0,39,10,0,0,2,0,37,8,0,0,2,0,38,6,0,0,
     + 3,0,22,7,6,0,3,0,22,9,22,0,2,0,22,7,0,0,3,1,2,15,22,0,3,1,4,17,
     + 22,0,3,1,2,15,31,0,3,1,4,17,31,0,2,0,31,25,0,0,3,0,33,7,6,0,
     + 3,0,10,7,7,0,
     + 3,0,21,8,6,0,3,0,21,10,21,0,2,0,21,8,0,0,3,1,3,16,21,0,3,1,5,18,
     + 21,0,3,1,3,16,30,0,3,1,5,18,30,0,2,0,30,26,0,0,3,0,33,8,6,0,
     + 3,0,9,8,8,0,
     + 2,0,29,7,0,0,2,0,31,6,0,0,2,0,22,6,0,0,2,0,10,7,0,0,2,0,31,27,0,
     + 0,2,0,30,27,0,0,2,0,29,25,0,0,3,1,2,15,10,0,3,1,2,15,29,0,
     + 3,1,4,17,10,0,3,1,4,17,29,0,
     + 2,0,28,8,0,0,2,0,30,6,0,0,2,0,21,6,0,0,2,0,9,8,0,0,2,0,30,27,0,
     + 0,2,0,31,27,0,0,2,0,28,26,0,0,3,1,3,16,9,0,3,1,3,16,28,0,
     + 3,1,5,18,9,0,3,1,5,18,28,0,
     + 3,0,6,21,22,0,3,0,6,9,10,0,3,0,23,6,6,0,3,0,23,7,8,0,3,0,24,6,6,
     + 0,3,0,24,7,8,0,
     + 2,0,71,7,0,0,2,0,59,6,0,0,2,0,59,1,0,0,
     + 2,0,72,8,0,0,2,0,60,6,0,0,2,0,60,1,0,0,
     + 2,0,71,6,0,0,2,0,71,1,0,0,2,0,72,6,0,0,2,0,72,1,0,0,
     + 2,0,2,3,0,0,2,0,4,5,0,0,3,0,6,7,8,0,
     + 2,0,89,7,0,0,2,0,89,6,0,0,2,0,89,8,0,0,
     + 3,1,2,15,22,0,3,1,2,15,33,0,3,1,4,17,22,0,3,1,4,17,33,0,2,0,7,22,
     + 0,0,2,0,9,22,0,0,2,0,7,33,0,0,2,0,9,33,0,0,
     + 3,1,2,15,10,0,3,1,4,17,10,0,2,0,7,10,0,0,2,0,9,10,0,0,
     + 3,0,7,10,13,0,3,0,7,22,14,0,3,0,7,8,13,0,3,0,9,10,13,0,3,0,9,22,
     + 14,0,3,0,22,8,40,0,3,1,2,15,39,0,3,1,2,15,14,0,3,1,4,17,39,0,3,
     + 1,4,17,14,0,
     + 2,0,89,7,0,0,2,0,89,6,0,0,2,0,89,8,0,0,
     + 2,0,87,6,0,0,2,0,87,1,0,0,2,0,88,6,0,0,2,0,88,1,0,0,
     + 3,1,2,15,10,0,3,1,4,17,10,0,2,0,7,10,0,0,2,0,9,10,0,0 ,
     + 2,0,74,1,0,0 ,2,0,74,6,0,0 , 2,0,75,1,0,0 ,2,0,75,6,0,0, !C=1,S=1 mesons
     + 2,0,23,25,0,0, 4,0,9,10,7,6, 3,0,9,10,7,0, 2,0,33,7,0,0, 
     + 3,1,23,2,15,0, 3,1,33,2,15,0, 2,0,23,7,0,0, 4,0,12,10,7,7,
     + 2,0,9,12,0,0, 4,0,7,8,7,8, 2,0,23,26,0,0, 4,0,10,9,8,6, !  | D*(_s)
     + 3,0,10,9,8,0, 2,0,33,8,0,0, 3,1,23,3,16,0, 3,1,33,3,16,0,! v
     + 2,0,23,8,0,0, 4,0,12,9,8,8, 2,0,10,12,0,0, 4,0,7,8,7,8, ! ----
     & 2,0,14,7,0,0, 2,0,13,6,0,0, 3,0,14,7,6,0, 3,0,13,7,8,0, ! N-res
     & 2,0,13,8,0,0, 2,0,14,6,0,0, 3,0,13,6,8,0, 3,0,14,7,8,0, 
     & 3,0,14,7,6,0, 3,0,13,7,8,0, 2,0,14,7,0,0, 2,0,13,32,0,0,  
     & 2,0,39,9,0,0, 2,0,13,6,0,0, 2,0,13,23,0,0,
     & 3,0,13,8,6,0, 3,0,14,7,8,0, 2,0,13,8,0,0, 2,0,14,32,0,0,  
     & 2,0,39,21,0,0, 2,0,14,6,0,0, 2,0,14,23,0,0, ! ---
     & 2,0,25,8,0,0, 2,0,26,7,0,0,    ! pi10 |  
     & 2,0,25,6,0,0, 2,0,27,7,0,0,    !  +   v
     & 2,0,27,8,0,0, 2,0,26,6,0,0,    !  -  ---
     & 2,0,21,7,0,0, 2,0,9,6,0,0, 2,0,22,8,0,0, 2,0,10,6,0,0,    ! k0* |  
     & 2,0,9,8,0,0, 2,0,21,6,0,0, 2,0,10,7,0,0, 2,0,22,6,0,0/ !        v
      DATA LBARP/1,3,2,5,4,6,8,7,10,9,11,12,-13,-14,16,15,18,17,13,14,
     +  22,21,23,24,26,25,27,29,28,31,30,32,33,-34,-35,-36,-37,-38,-39,
     +  -40,-41,-42,-43,-44,-45,-46,-47,-48,-49,0,-51,-52,-53,-54,4*0,
     +  60,59,61,63,62,65,64,67,66,3*0,72,71,
     +  73,75,74,77,76,79,78,81,80,0,83,-84,-85,-86,-87,-88,-89,
     +  91,90,93,92,-94,-95,-96,-97,-98,-99 /
      DATA ICHP /0,1,-1,1,-1,0,1,-1,1,-1,0,0,1,0,4*0,-1,0,4*0, !24
     + 1,-1,0,1,-1,4*0,1,0,-1,0,-1,0,2,1,0,-1,1,0,-1,0,-1,-1,  !49
     + 0,1,0,1,0,4*0,1,-1,0,1,-1,1,-1,0,0,3*0,                 !70
     + 0,0,0,1,-1,1,-1,1,-1,0,0,0,                             !82
     + 0,2,1,0,1,0,1,1,-1,2*0,2,1,0,1,0,0 / ! charmed baryons + tau

      DATA ISTR /8*0,-1,+1,-1,-1,8*0,-1,+1,5*0,-1,+1,-1,+1,2*0, ! mesons
     +     3*1,2*2,1,4*0,3*1,2*2,3,0,4*0, ! 54
     +     4*0,2*0,3*0,-1,1,-1,1,3*0,2*0,0,-1,1,-1,1,2*0,2*0,0,0, ! 83
     +     3*0,2*1,0,4*0,3*0,2*1,2 / ! charmed baryons
      DATA IBAR /12*0,2*1,4*0,2*-1,13*0,16*1,0,4*1,4*0,
     +     2*0,10*0,2*0,0,4*0,2*0,2*0,0,0,6*1,4*0,6*1 /
      DATA ICHM /58*0,1,-1,10*0,1,-1,0,1,-1,+1,-1,1,-1,1,-1,0,0,
     +     6*1,4*0,6*1/
      DATA NAMP /
     +     '     ','gam   ','e+','e-','mu+','mu-','pi0',
     +     'pi+','pi-','k+', 'k-', 'k0l','k0s',
     +     'p', 'n', 'nue', 'nueb', 'num', 'numb', 'pbar', 'nbar',
     +     'k0', 'k0b', 'eta', 'etap', 'rho+', 'rho-','rho0',
     +     'k*+','k*-','k*0','k*0b','omeg', 'phi', 'SIG+', 'SIG0',
     +     'SIG-','XI0','XI-','LAM','DELT++','DELT+','DELT0','DELT-',
     +     'SIG*+','SIG*0','SIG*-', 'XI*0', 'XI*-', 'OME-',
     +     '     ','N144_+','N144_0','N171_+','N171_0',
     +     4*'     ', 'D+', 'D-','pi1_0 ','pi1_+ ','pi1_- ',
     +     'k0*_+','k0*_-','k0*_0','k0*_0b',
     +     3*'     ', 'D0', 'D0b', 'eta_c', 
     +     'D_s+','D_s-','D*_s+','D*_s-','D*+', 'D*-', 'D*0', 'D*0b',
     +     '     ', 'J/psi',
     +     'SIGc++', 'SIGc+', 'SIGc0','XI_c+','XI_c0','LAM_c+',
     +     'tau+  ','tau-  ','nut   ','nutb  ',
     +     'SIc*++','SIGc*+','SIGc*0', 'XI_c*+', 'XI_c*0',
     +     'OME_c0'  /
      END
C->
C=======================================================================

      SUBROUTINE DECPR (LUN)

C-----------------------------------------------------------------------
C...Print on unit LUN the list of particles and decay channels
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      DOUBLE PRECISION CBR
      INTEGER KDEC,LBARP,IDB
      COMMON /S_CSYDEC/ CBR(223+16+12+8), KDEC(1338+6*(16+12+8)),
     &     LBARP(99), IDB(99)
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)

      DOUBLE PRECISION AW,AW2
      COMMON /S_WIDTH1/ AW(99), AW2(99)

      CHARACTER*6 NAMP
      COMMON /S_CNAM/ NAMP (0:99)
      DIMENSION LL(4)
      SAVE
      
 100  FORMAT(/,1X,75('-'),/,28X,'SIBYLL DECAY TABLE')
      WRITE(LUN,100)
 101  FORMAT(1X,75('-'),/,2X,'PID',1X,'Particle',6X,'Mass',9X,'Width',/,
     +     4X,'Channel',1X,'Br.frac.',1X,'Nf',2X,'MAT',1X,
     +     'Final Particles',/,1X,75('-'))
      WRITE(LUN,101)
      DO L=1,99
         IF(MOD(L,10).EQ.0)WRITE(LUN,101)
         IDC = IDB(L)-1
         NC = 0
         WRITE (LUN,10) L,NAMP(L), AM(L), AW(L)
         IF(IDC+1 .GT. 0)  THEN
            CB = 0.D0
110         IDC=IDC+1
            NC = NC+1
            CBOLD = CB
            CB = CBR(IDC)
            BR = CB-CBOLD
            KD = 6*(IDC-1)+1
            ND = KDEC(KD)
            MAT= KDEC(KD+1)
            DO J=1,ND
               LL(J) = KDEC(KD+1+J)
            ENDDO
            WRITE (LUN,15) NC,BR,ND,MAT, (NAMP(LL(J)),J=1,ND)
            IF (CB .LT. 1.D0)  GOTO 110
         ENDIF
      ENDDO
      RETURN
10    FORMAT(2X,I3,2X,A6,3X,F10.4,3X,F10.4)
15    FORMAT(5X,I2,2X,F9.4,I4,I4,2X,3(A6,2X))
      END

C=======================================================================

      SUBROUTINE DEC_DEBUG (L,P0, ND, LL, PD)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      CHARACTER*6 NAMP
      COMMON /S_CNAM/ NAMP (0:99)
      DIMENSION P0(5), LL(10), PD(10,5)
      SAVE

      ETOT = 0.D0
      DO J=1,ND
         ETOT = ETOT + PD(J,4)
      ENDDO
      WRITE(*,*)  NAMP(IABS(L)),' -> ', (NAMP(IABS(LL(J))),J=1,ND)
      WRITE(*,*)  ' Ei, Ef = ', P0(4), ETOT, ' L = ', L
      RETURN
      END
C=======================================================================

      SUBROUTINE DEC_INI

C-----------------------------------------------------------------------
C     decay initialization routine
C     sets which particles should decay and wich should be stable
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      DOUBLE PRECISION CBR
      INTEGER KDEC,LBARP,IDB
      COMMON /S_CSYDEC/ CBR(223+16+12+8), KDEC(1338+6*(16+12+8)),
     &     LBARP(99), IDB(99)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      SAVE
      
      if ( ndebug .gt. 0 ) then
        write(lun,*)' -----------------------------------------'
        write(lun,*)' SIBYLL DEC_INI: setting particle decays!'
        write(lun,*)'  to be used in stand-alone SIBYLL only ! '
        write(lun,*)' -----------------------------------------'
      endif
      
C...  Definition of stable particles
      DO J=4,12
         IDB(J) = -abs(IDB(J))
      ENDDO
c----------------------------------------------------------
c     if the folowing is commented out then all particles
c     except leptons, protons and neutrons are UNSTABLE
c----------------------------------------------------------
c     all particles with t<0.3e-10s are considered unstable
c     i.e. all the mesons from K0s onwards(K0l is stable)
c----------------------------------------------------------
C     K0s stable
      if (ndebug .gt. 0 ) write(lun,*)' making K0s stable..'
      IDB(12) = -abs(IDB(12))

C     Lambda/Anti-lambda stable
      if (ndebug .gt. 0 ) write(lun,*)' making LAMBDA stable..'
      IDB(39) = -abs(IDB(39))

c     Sigmas stable
      if (ndebug .gt. 0 ) write(lun,*)' making SIGMAs stable..'
      do i=34,36
         IDB(i) = -abs(IDB(i))
      enddo
      IDB(35) = -abs(IDB(35))

      if (ndebug .gt. 0 ) 
     *       write(lun,*)' ------------------------------------------'
      end
C=======================================================================

      SUBROUTINE DEC_INI_EXT

C-----------------------------------------------------------------------
C     decay initialization routine for muon enhancements
C     sets which particles should decay and wich should be stable
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      DOUBLE PRECISION CBR
      INTEGER KDEC,LBARP,IDB
      COMMON /S_CSYDEC/ CBR(223+16+12+8), KDEC(1338+6*(16+12+8)),
     &     LBARP(99), IDB(99)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      SAVE
      
      if ( ndebug .gt. 0 ) then
        write(lun,*)' -----------------------------------------'
        write(lun,*)' SIBYLL DEC_INI: setting particle decays!'
        write(lun,*)'  to be used in stand-alone SIBYLL only ! '
        write(lun,*)' -----------------------------------------'
      endif

C...  Make all particles unstable
      DO J=1,99
         IDB(J) = abs(IDB(J))
      ENDDO
      
C...  Definition of stable particles (muon, pions, kaons)
      DO J=4,12
         IDB(J) = -abs(IDB(J))
      ENDDO

C     Lambda/Anti-lambda stable
      if (ndebug .gt. 0 ) write(lun,*)' making LAMBDA stable..'
      IDB(39) = -abs(IDB(39))

c     Charged sigmas stable
      if (ndebug .gt. 0 ) write(lun,*)' making SIGMAs stable..'
      IDB(34) = -abs(IDB(34))
      IDB(36) = -abs(IDB(36))

      if (ndebug .gt. 0 ) 
     *       write(lun,*)' ------------------------------------------'
      end
C=======================================================================
      
      SUBROUTINE STRING_FRAG_4FLV
     +     (E0,IFL1,IFL2,PX1,PY1,PX2,PY2,IFBAD,IFQRK)

C-----------------------------------------------------------------------
C.  This routine fragments a string of energy E0
C.  the ends of the strings  have flavors IFL1 and IFL2
C.  the particles produced are in the  jet-jet frame
C.  with IFL1 going in the +z direction
C.     E0 = total energy in jet-jet system
C.  This version consider also a primordial pT attached
C.  to the ends of the string PX1,PY1,  PX2,PY2
C.  OUTPUT:  IFBAD =1  kinematically impossible decay
c       2010.03.11 ifqrk - leading quark flag
c       1 in valence quark, 0 in others
c
c      Modified Nov. 91.  RSF and TSS to fragment symmetrically
c      ie forward and backward are fragmented as leading.
c      Change- Dec. 92  RSF.  call to ptdis moved- to use flavor
c      of NEW quark in fragmentation.
c
c     includes 4 FLAVORS \FR'13
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
C     The final particle output is contained in COMMON /S_PLIST/    
C     NP           : number of final particles
C     P(1:NP, 1:5) : 4-momenta + masses of the final particles 
C     LLIST (1:NP) : codes of final particles
      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP
      DOUBLE PRECISION ZLIST
      COMMON /S_ZLIST/ ZLIST(8000)
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)
      INTEGER NW_max
      PARAMETER (NW_max = 20)
C--------------------------------------------------------------------
C     SIBYLL common blocks containing event information       \FR'14
C--------------------------------------------------------------------

C     EVENT INFO COMMON
C     contains overall interaction properties, like
C     SQS : center-of-mass energy
C     S   :         "       "     squared
C     PTmin : low pt cut of QCD cross section, 
C             i.e. minimal pt of hard minijets
C     Xmin : low-x bound for PDFs, 
C            i.e. minimal momentum fraction of hard partons
C     Zmin : logarithm of that
C     KB : PID of beam hadron
C     KT() : PID of target
C     IAT : mass number of target
      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
      INTEGER KB,IAT,KT
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)

      INTEGER NFORIG,NPORIG,NIORIG,IPFLAG,IIFLAG,KINT
      COMMON /S_PARTO/ NFORIG(NP_max),NPORIG(NP_max),NIORIG(NP_max),
     &IPFLAG,IIFLAG,KINT
      DOUBLE PRECISION FAin, FB0in
      COMMON /S_CZDIS/ FAin, FB0in

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN

      INTEGER LRNK
      COMMON /SIB_RNK/ LRNK(8000)

      INTEGER ICHP,ISTR,IBAR
      COMMON /S_CHP/ ICHP(99), ISTR(99), IBAR(99)

      INTEGER IISO,ISPN
      COMMON /S_SPN/ IISO(99), ISPN(99)

      INTEGER ICHM
      COMMON /S_CHM/ ICHM(99)

      DIMENSION WW(2,2), PTOT(4), PX(3),PY(3),IFL(3),ILEAD(2)
      DIMENSION LPOINT(8000), PMQ(3), IRNK(2), LRES(6:99)
      LOGICAL LRANK
      SAVE
      DATA LRANK/.true./

      DATA (LRES(I),I=6, 39)
     &     /27,25,26,28,29,9,9,41,42,19*0,44,45,46,47,48,39/
      DATA (LRES(I),I=40, 49) /40,41,42,43,44,45,46,47,48,49/
      DATA (LRES(I),I=50, 83) 
     &     /0,51,52,53,54,4*0,78,79,10*0,71,72,73,76,77,76,
     &     77,78,79,80,81,0,83/
      DATA (LRES(I),I=84, 99) /94,95,96,97,98,89,4*0,94,95,96,97,98,99/
      
      IF(Ndebug.gt.3) THEN
        WRITE(LUN,*)
     &        ' STRING_FRAG_4FLV: called with ',
     &        '(E0,IFL1,IFL2,PX1,PY1,PX2,PY2,IVAL)',
     &        E0,IFL1,IFL2,PX1,PY1,PX2,PY2,IFQRK
        WRITE(LUN,*)' STRING_FRAG_4FLV: NP before fragmentation:',NP
      ENDIF

c...  remember initial values
c     strange fraction
      par2_def = PAR(2)
c     vector model
      IPAR11_def = IPAR(11)
c     vector fraction
      PAR5_def = PAR(5)      
c     charm fraction
      PAR24_def = PAR(24)
c     popcorn fraction
      PAR8_def = PAR(8)

C...initialise
      NTRY = 0
      IFBAD = 0
 200  NTRY = NTRY + 1

c     reset parameters after rejection
      PAR(2) = PAR2_def
      PAR(5) = PAR5_def
      PAR(24) = PAR24_def
      IPAR(11) = IPAR11_def
      PAR(8) = PAR8_def

      IF (NTRY .GT. 50)  THEN
         IFBAD = 1
         RETURN
      ENDIF
      I = NP
      DO K=1,2
         WW(K,1) = 1.D0
         WW(K,2) = 0.D0
         IRNK(K) = 0
      ENDDO
      PX(1) = PX1
      PY(1) = PY1
      PX(2) = PX2
      PY(2) = PY2
      PX(3) = 0.D0
      PY(3) = 0.D0
      PTOT (1) = PX1+PX2
      PTOT (2) = PY1+PY2
      PTOT (3) = 0.D0
      PTOT (4) = E0
c     turn on/off splitting of leading diquark
c     (1: no splitting, 0: diq may be split, producing leading meson)
      IFL(1) = IFL1+ISIGN(100,IFL1)*MIN(1,IABS(IFL1)/10)*IPAR(90)
      IFL(2) = IFL2+ISIGN(100,IFL2)*MIN(1,IABS(IFL2)/10)*IPAR(90)     
      PMQ(1) = QMASS(IFL(1))
      PMQ(2) = QMASS(IFL(2))

      ILEAD(1) = 0
      ILEAD(2) = 0
      IBLEAD = 0
      IF(IABS(IFQRK).eq.1) THEN
         ILEAD(1) = 1
         ILEAD(2) = 1
      ENDIF
c     switch leading baryon fragmentation function on/off
      IF(IPAR(20).eq.0) GOTO 300
c     set flags for leading baryon
C
C      SET FLAG FOR GENERATION OF LEADING PARTICLES. 
C      "AND" IS FOR PPBAR ( DIQUARK AT BOTH ENDS)
C      "OR" IS FOR PP, PPI, ( DIQUARK AT ONE END.)
C
      IF (IABS(IFL1) .GT. 10 .AND. IABS(IFL2) .GT. 10)  THEN
         IBLEAD = 2
         I = I+1
         JT = INT(1.5D0+S_RNDM(0))
         GOTO 350
      ENDIF         
      IF (IABS(IFL1) .GT. 10 .OR. IABS(IFL2) .GT. 10)  THEN
         IBLEAD = 1
         I = I+1
         JT = 2
         IF (IABS(IFL2) .GT. 10) JT = 1
         GOTO 350
      ENDIF         

C...produce new particle: side, pT
 300  continue
      I=I+1
      if(i.gt.8000) then
        write(LUN,'(1x,a,i8)') 
     &        ' STRING_FRAG_4FLV: no space left in S_PLIST:',I
        CALL SIB_REJECT('STRING_FRAG_4FLV')
      endif
      IF (IBLEAD .GT. 0)  THEN
         JT = 3 - JT   
         GO TO 350              
      ENDIF
c     
c 349  continue
c     choose side (1 or 2)
      JT=INT(1.5D0+S_RNDM(0))
c     set 'other' side
 350  JR=3-JT
c     remember side particle was produced
      LPOINT(I) = JT
c     increase rank counter
      IRNK(JT) = ISIGN(ABS(IRNK(JT))+1,1-JT)
c     set particle rank
      LRNK(I) = IRNK(JT)

      nporig(I)= Ipflag*2 + KINT
      niorig(I)= iiflag
      IF(ILEAD(JT).eq.1) nporig(I)= -1 * nporig(I)
      nforig(I) = 0

 555  CONTINUE
c
c.... CHARM config
c
      charmPARdef=PAR(24)
      IF(IPAR(15).lt.9)THEN
c     no s->c
         PAR(24) = 0.D0
         IF (IFQRK.EQ.1) THEN
c     ifqrk = 1 (valence quark attatched) 
            IF(IPAR(15).ge.1) THEN
c     enforce s->c at string end
               IF(ILEAD(JT).eq.1) PAR(24)=charmPARdef
c     produce charm in all strings
               IF(IPAR(15).eq.8) PAR(24)=charmPARdef
            ELSE
c     compatibility to broken version
               PAR(24)=charmPARdef
            ENDIF
         ELSE
c     no val. quark at string end or diff
            PAR(24)=charmPARdef
         ENDIF
      ENDIF
c
C.... Vector meson config
c
c     increase vec.meson ratio for leading particle in str. diff.
      IF(IFQRK.eq.-1)THEN
         IF(IPAR(66).eq.1)THEN
            IF(ILEAD(JT).EQ.1)THEN
              IF(IBAR(IABS(kb)).eq.0.or.IPAR(70).eq.1) PAR(5) = PAR(113)
           ENDIF
        ELSEIF(IPAR(66).eq.2)THEN
            IF(IBAR(IABS(kb)).eq.0.or.IPAR(70).eq.1) PAR(5) = PAR(113)

        ELSEIF(IPAR(66).eq.3)THEN
c     increase vector meson rate for meson beam
c     on beam side (rank+) only!            
            IF(ILEAD(JT).EQ.1)THEN               
               IF(IBAR(IABS(kb)).eq.0.and.IRNK(JT).gt.0)
     &              PAR(5) = PAR(113)
c     always incr. vector rate for diff. strings independent of beam type
               IF(IPAR(70).eq.1) PAR(5) = PAR(113)               
            ENDIF         
            
         ENDIF           
      endif
      
c...  switch off for proton beam
      IF(IPAR(31).eq.1)then
c         print*,'ipar11,ipar11def,1-kb/13,kb',IPAR(11),ipar11_def,
c     +        max((1-iabs(kb)/13),0),kb
         IPAR(11) = IPAR(11)*max((1-iabs(kb)/13),0) ! meson beam only
      endif
c     increase vec.meson ratio for leading quarks
      IF(IABS(IFQRK).eq.1)THEN
         IF(IPAR(11).le.-5.and.IPAR(11).ge.-7
     &        .and.ilead(jt).eq.1)
     &        PAR(5) = 9.D0
         
c     increase vec.meson ratio for diff.
         IF(IFQRK.eq.-1.and.IPAR(11).le.-4.and.IPAR(11).ge.-7)
     &        PAR(5) = 9.D0

c     increase vec.meson ratio for leading particle in str. diff. (lvec16)
         IF(IFQRK.eq.-1.and.IPAR(11).le.-11.and.ILEAD(JT).EQ.1)
     &        PAR(5) = 99.D0
      ENDIF

c...  suppress leading charm for pion and kaon beams
      IF(IPAR(15).eq.11)then
         IF((1-IABS(KB)/13)*ILEAD(JT).gt.0) PAR(24)=0.D0
      ENDIF

C...  suppress rank-1 baryon through popcorn
      IF(IBLEAD .GT. 0.and.abs(ifl(jt)).gt.10
     &     .and.abs(ifl(3)).lt.10) PAR(8)=PAR(63)*PAR(8)

C...  leading strange/charm
      IF(ILEAD(JT).eq.1.and.IPAR(39).gt.0) PAR(2) = PAR(65)
      
c     scale valence string end charm for assoc. prod.      
      IF(IPAR(41).eq.1)THEN
         IF(ILEAD(JT).eq.1.and.IFQRK.eq.1) PAR(24) = PAR(71)*PAR(24)
      ENDIF

c     suppress direct pi0 for meson projectiles
c     rate set by par( 137 )
      ipar82_def = IPAR(82)
c     skip if baryon projectile or minijet (i.e no flavor attached)
      if(ibar(iabs(kb)).ne.0.or.ifqrk.eq.0) IPAR(82) = 0

c     suppress direct omega for meson projectiles
c     rate set by par( 138 )
      ipar83_def = IPAR(83)
c     skip if baryon projectile or central string
      if(ibar(iabs(kb)).ne.0.or.(ifqrk.gt.0.and.IPAR(83).eq.2))
     &     IPAR(83) = 0

c     change rho0 / omega ratio
      PAR143_def = PAR(143)
      IF(IPAR(81).eq.1)THEN
c     change if beam is meson
         if(ibar(iabs(kb)).eq.0) PAR(143) = PAR(144)
      ELSEIF(IPAR(81).eq.2)THEN         
c     change if beam is meson, on meson side only         
         if(ibar(iabs(kb)).eq.0.and.IRNK(JT).gt.0) PAR(143) = PAR(144) 
      ELSEIF(IPAR(81).eq.3)THEN         
c     change if beam is meson, on meson side only, for leading only
         if(ibar(iabs(kb)).eq.0.and.ISIGN(ILEAD(JT),IRNK(JT)).eq.1)
     &        PAR(143) = PAR(144)
      ELSEIF(IPAR(81).eq.4)THEN         
c     change if beam is meson, on meson side only, for diff. strings only
         if(ibar(iabs(kb)).eq.0.and.IFQRK.eq.-1)
     &        PAR(143) = PAR(144)
      ELSEIF(IPAR(81).eq.5)THEN         
c     change if beam is meson, for leading on meson side only and
c     for diff. strings only
         if(ibar(iabs(kb)).eq.0.and.IFQRK.eq.-1.and.
     &        ISIGN(ILEAD(JT),IRNK(JT)).eq.1) PAR(143) = PAR(144)
      ENDIF
      
C...particle ID and pt.

      CALL SIB_I4FLAV (IFL(JT), 0, IRNK(JT), IFL(3), LLIST(I))

c     reset strange fraction
      PAR(2) = PAR2_def
c     reset vec.meson production
      PAR(5) = PAR5_def
c     reset charm fraction
      PAR(24) = PAR24_def
c     reset popcorn
      PAR(8) = par8_def

c     reset pi0 suppr.
      IPAR(82) = ipar82_def

c     reset omega suppr.
      IPAR(83) = ipar83_def

c     reset rho0 / omega ratio
      PAR(143) = PAR143_def
      
c     reject iso 0 spin 1 for meson projectiles
      IF(IBAR(IABS(KB)).eq.0)THEN
c     reject leading spin1,isospin singlett
         IF(ILEAD(JT).EQ.1.and.LLIST(I).eq.32.and.
     +        PAR(136).gt.S_RNDM(I)) LLIST(I) = 27
      endif
      
c     replace leading or all pi0 with rho0
      IF(IFQRK.eq.-1) THEN 
         IF(IPAR(67).eq.1)THEN
            IF(ILEAD(JT).EQ.1) THEN 
c     replace leading pi0 with rho0
               IF(ABS(LLIST(I)).EQ.6) LLIST(I) = 27*isign(1,LLIST(I))  
            ENDIF
         ELSEIF(IPAR(67).eq.2)THEN
c     replace all pi0 with rho0 for all beams
            IF(ABS(LLIST(I)).EQ.6) LLIST(I) = 27*isign(1,LLIST(I))
         ELSEIF(IPAR(67).eq.3)THEN
c     replace all pi0 with rho0 for meson beam only
            IF(IBAR(IABS(KB)).eq.0)THEN
               IF(ABS(LLIST(I)).EQ.6) LLIST(I) = 27*isign(1,LLIST(I))
            ENDIF
         ELSEIF(IPAR(67).eq.4)THEN
c     replace all pi0 with rho0 for meson beam only
c     replace some beam mesons with their vector partner
            IF(IBAR(IABS(KB)).eq.0)THEN
               IF(ABS(LLIST(I)).EQ.6) LLIST(I) = 27*isign(1,LLIST(I))
c     reject leading spin1,isospin singlett
               IF(ILEAD(JT).EQ.1.and.LLIST(I).eq.32.and.
     +              PAR(136).gt.S_RNDM(I)) LLIST(I) = 27
               IF(S_RNDM(0).lt.PAR(120).and.LLIST(I).eq.KB) 
     &              LLIST(I) = LRES(LLIST(I))
            ENDIF
         ENDIF
      ENDIF

c     replace leading pi0 by rho0's
      IF(IABS(IFQRK).eq.1)THEN
         IF(ABS(IPAR(11)).ge.2.and.IPAR(11).ge.-3)THEN
            IF(ilead(jt).EQ.1) then 
               IF(ABS(LLIST(I)).EQ.6) THEN
                  LLIST(I) = 27*isign(1,LLIST(I))
               endif
            endif
        
c     replace leading pi0 in string diff by rho0's (lvec15)
         ELSEIF(IFQRK.eq.-1.and.IPAR(11).eq.-10)THEN
            IF(ILEAD(JT).EQ.1) THEN 
               IF(ABS(LLIST(I)).EQ.6) THEN
                  LLIST(I) = 27*isign(1,LLIST(I))
               ENDIF
            ENDIF
c     replace leading pi0 in string diff by rho0's 
c     in addition to increased leading vec.meson ratio (lvec20)
         ELSEIF(IFQRK.eq.-1.and.IPAR(11).eq.-15)THEN
            IF(ILEAD(JT).EQ.1) THEN 
               IF(ABS(LLIST(I)).EQ.6) THEN
                  LLIST(I) = 27*isign(1,LLIST(I))
               ENDIF
            ENDIF     
c     replace leading omega in string diff by rho0's 
c     in addition to increased leading vec.meson ratio (lvec21)
         ELSEIF(IFQRK.eq.-1.and.IPAR(11).eq.-16)THEN
            IF(ILEAD(JT).EQ.1) THEN 
               IF(ABS(LLIST(I)).EQ.32) 
     &              LLIST(I) = 27*isign(1,LLIST(I))
            ENDIF     
c     replace leading omega in string diff by rho0's 
c     suppress pi0 in diff. strings
c     in addition to increased leading vec.meson ratio (lvec22)
         ELSEIF(IFQRK.eq.-1.and.IPAR(11).eq.-17)THEN
            IF(ILEAD(JT).EQ.1) THEN 
c     print*,'replacing leading omega with rho0'
               IF(ABS(LLIST(I)).EQ.32)
     &              LLIST(I) = 27*isign(1,LLIST(I))
            ENDIF
            IF(LLIST(I).EQ.6) then
c     print*,'pi0 found! start again.. '
               GOTO 555
            endif

c     replace all for diff.
         ELSEIF(IFQRK.eq.-1.and.IPAR(11).lt.0.and.
     &           IPAR(11).ge.-3) then
            IF(ABS(LLIST(I)).EQ.6)  LLIST(I) = 27*isign(1,LLIST(I))

c     increased vec.meson ratio and replace pi0 with rho0 in str.diff
         ELSEIF(IFQRK.eq.-1.and.IPAR(11).eq.-7) then
            IF(ABS(LLIST(I)).EQ.6)  LLIST(I) = 27*isign(1,LLIST(I))  

c     replace leading pi's by vec.mesons, iso-spin conserving
         ELSEIF(IPAR(11).eq.-8.and.IPAR(11).lt.0)THEN
            PAR(5) = 9.D0
            IF(ilead(jt).EQ.1.and.
     $                       INT((PAR(5)+1.D0)*S_RNDM(0)).gt.1) then 
               IF(ABS(LLIST(I)).EQ.6) LLIST(I) = 27*isign(1,LLIST(I))
               IF(ABS(LLIST(I)).EQ.7) LLIST(I) = 25*isign(1,LLIST(I))
c     IF(ABS(LLIST(I)).EQ.8) LLIST(I) = 26*isign(1,LLIST(I))
            endif

c     replace almost all for diff.
         ELSEIF(IFQRK.eq.-1.and.IPAR(11).eq.-8.and.IPAR(11).lt.0) then
            PAR(5) = 9.D0
            if( INT((PAR(5)+1.D0)*S_RNDM(0)).gt.1 ) then
               IF(ABS(LLIST(I)).EQ.6)  LLIST(I) = 27*isign(1,LLIST(I))
               IF(ABS(LLIST(I)).EQ.7) LLIST(I) = 25*isign(1,LLIST(I))
            endif
      
c     replace leading pi0's by vec.mesons
         ELSEIF(IPAR(11).eq.-9.and.IPAR(11).lt.0)THEN
            PCHF = 0.1D0
            IF(ilead(jt).EQ.1.and.ABS(LLIST(I)).EQ.6) 
     &           LLIST(I) = 27*isign(1,LLIST(I))
            if(ilead(jt).EQ.1.and.ABS(LLIST(I)).EQ.7)then
               if(S_RNDM(0).lt.PCHF) LLIST(I) = 25*isign(1,LLIST(I))
            endif        

c     replace for string diff.
         ELSEIF(IFQRK.eq.-1.and.IPAR(11).eq.-9) then
            IF(ABS(LLIST(I)).EQ.6) 
     &           LLIST(I) = 27*isign(1,LLIST(I))
            if(ABS(LLIST(I)).EQ.7)then
               if(S_RNDM(0).lt.PCHF) 
     &              LLIST(I) = 25*isign(1,LLIST(I))
            endif
         ELSE
            CONTINUE
         ENDIF
      ENDIF

c     reset vec.meson ratio
      PAR(5) = 0.3D0
      IF(IABS(IFQRK).eq.1) ILEAD(JT) = 0
      
      PMQ(3) = QMASS(IFL(3))
      P(I,5) = AM(IABS(LLIST(I)))
      CALL PTDIS_4FLV (IFL(3), PX(3),PY(3))

C...fill transverse momentum
      P(I,1) = PX(JT) + PX(3)
      P(I,2) = PY(JT) + PY(3)
      XMT2 = P(I,5)**2+P(I,1)**2+P(I,2)**2

C...test end of fragmentation

      WREM2 = PTOT(4)**2-PTOT(1)**2-PTOT(2)**2-PTOT(3)**2
c      IF (WREM2 .LT. 0.1)  GOTO 200
      IF (WREM2 .LT. 0.1D0)  GOTO 200
c      WMIN = PMQ(1)+PMQ(2)+2.*PMQ(3)+ 1.1 + (2.*S_RNDM(0)-1.)*0.2
      WMIN=PMQ(1)+PMQ(2)+2.D0*PMQ(3)+PAR(59)+(2.D0*S_RNDM(0)-1.D0)*0.2D0
      IF (WREM2 .LT. WMIN**2) Then
         if (IABS(ifl(3)).ne.3.and.IABS(IFL(3)).ne.4) GOTO 400
         goto 200
      endif

C...Choose z
      IF(IABS(IFQRK).eq.1) THEN
c     valence strings: ( str.diff and non diff. )
         IF(IPAR(11).EQ.1) THEN
c     use hard distribution for leading quarks ( no exchange )
            IF(ILEAD(JT).eq.1) THEN
               Z = ZBLEAD (IABS(LLIST(I)))
            ELSE
               IF(IBLEAD .GT. 0.and.iabs(ifl(jt)).gt.10
     &              .and.iabs(ifl(3)).lt.10)  THEN
                  Z = ZBLEAD (IABS(LLIST(I)))   
                  IBLEAD = IBLEAD - 1
               ELSE
                  Z = ZDIS_4FLV (IFL(3),ifl(jt),XMT2)
               ENDIF
            ENDIF
c     use hard frag. for leading particles
         ELSEIF(IPAR(11).ge.3.or.IPAR(11).eq.-3.or.IPAR(11).eq.-6
     &           .or.IPAR(11).eq.-7) THEN
            IF(ILEAD(jt).eq.1) THEN
               Z = ZBLEAD (IABS(LLIST(I)))
            ELSE
               IF(IBLEAD .GT. 0.and.iabs(ifl(jt)).gt.10
     &              .and.iabs(ifl(3)).lt.10)  THEN
                  Z = ZBLEAD (IABS(LLIST(I)))   
                  IBLEAD = IBLEAD - 1
               ELSE
                  Z = ZDIS_4FLV (IFL(3),ifl(jt),XMT2)
               ENDIF
            ENDIF
         ELSEIF(IPAR(11).EQ.-11) THEN
c     very hard leading frag. for diff and non. diff val. strings (lvec16)
            IF(IBLEAD .GT. 0.and.iabs(ifl(jt)).gt.10) THEN
               Z = ZBLEAD (IABS(LLIST(I)))
            ELSEIF(ILEAD(jt).eq.1)THEN
               Z = 1.D0 - ZDISN(1)
            ELSE
               Z = ZDIS_4FLV (IFL(3),ifl(jt),XMT2)
            ENDIF

         ELSEIF(IPAR(11).EQ.-12.OR.IPAR(11).LE.-15.or.IPAR(68).eq.1)THEN
c     very hard leading frag. for diff. val. strings only (lvec17)
            IF(IBLEAD .GT. 0.and.iabs(ifl(jt)).gt.10) THEN
               Z = ZBLEAD (IABS(LLIST(I)))
            ELSEIF(ILEAD(jt).eq.1.and.IFQRK.eq.-1)THEN
               Z = 1.D0 - ZDISN(1)
            ELSE
               Z = ZDIS_4FLV (IFL(3),ifl(jt),XMT2)
            ENDIF

         ELSEIF(IPAR(11).EQ.-13.AND.IFQRK.eq.-1) THEN
c     hard leading frag. for diff. val. strings only (lvec18)
            IF(IBLEAD .GT. 0.and.iabs(ifl(jt)).gt.10) THEN
               Z = ZBLEAD (IABS(LLIST(I)))
            ELSEIF(ILEAD(jt).eq.1)THEN
               Z = S_RNDM(JT)
            ELSE
               Z = ZDIS_4FLV (IFL(3),ifl(jt),XMT2)
            ENDIF
         ELSEIF(IPAR(11).EQ.-14.AND.IFQRK.eq.-1) THEN
c     hard leading frag. for diff. AND ndiff. val. strings (lvec19)
            IF(IBLEAD .GT. 0.and.iabs(ifl(jt)).gt.10) THEN
               Z = ZBLEAD (IABS(LLIST(I)))
            ELSEIF(ILEAD(jt).eq.1)THEN
               Z = S_RNDM(JT)
            ELSE
               Z = ZDIS_4FLV (IFL(3),ifl(jt),XMT2)
            ENDIF
            
         ELSE

c     hard leading baryons only ( standard )
            IF(IBLEAD .GT. 0.and.iabs(ifl(jt)).gt.10
     &           .and.abs(ifl(3)).lt.10)  THEN
c           print*,'calling zblead: i,id,jt,ncall', i,llist(i),jt,ncall
               IF(IPAR(20).eq.3)THEN
c     use lund function with different parameters for leading baryon
                  fa_def = FAin
                  fb_def = FB0in
                  FAin = PAR(57)
                  FB0in = PAR(58)
                  z = zdis_4flv(IFL(3),ifl(jt),xmt2)
c     set parameters to initial values again
                  FAin = fa_def
                  FB0in = fb_def
               ELSE
                  Z = ZBLEAD (IABS(LLIST(I)))
               ENDIF
               IBLEAD = IBLEAD - 1
            ELSE
               Z = ZDIS_4FLV (IFL(3),ifl(jt),XMT2)
            ENDIF
         ENDIF
      ELSE
c     non valence string
         IF (IBLEAD .GT. 0.and.iabs(ifl(jt)).gt.10
     &        .and.iabs(ifl(3)).lt.10)  THEN
C     Special frag. for leading Baryon only
c            print*,'calling zblead: i,id,jt,ncall', i,llist(i),jt,ncall
            Z = ZBLEAD (IABS(LLIST(I)))   
            IBLEAD = IBLEAD - 1
         ELSE
            Z = ZDIS_4FLV (IFL(3),ifl(jt),XMT2)
         ENDIF
      ENDIF
      IF(IPAR(20).eq.2)IBLEAD = 2
      IF(IFQRK.eq.1) ILEAD(JT) = 0

      ZLIST(I) = Z
      WW(JT,2) = Z*WW(JT,1)
      WW(JR,2) = XMT2/(WW(JT,2)*E0**2)

      P(I,3) = WW(1,2)*0.5D0*E0 - WW(2,2)*0.5D0*E0
      P(I,4) = WW(1,2)*0.5D0*E0 + WW(2,2)*0.5D0*E0

      DO J=1,4
         PTOT (J) = PTOT(J) - P(I,J)
      ENDDO
      DO K=1,2
         WW(K,1) = WW(K,1) - WW(K,2)
      ENDDO

C...Reset pT and flavor at ends of the string
      PX(JT) = -PX(3)
      PY(JT) = -PY(3)
      IFL(JT) =-IFL(3)
      PMQ(JT) = PMQ(3)

      GOTO 300

C...Final two hadrons
 400  IAFL1 = IABS(mod(IFL(JR),100))
      IAFL2 = IABS(mod(IFL(3),100))
      IF(NDEBUG.gt.5)
     &     write(lun,*)'STRING_FRAG: final flavors:', IFL(JR), -IFL(3)
      
C..   check if flavor combination is allowed..
      
c     reject anti-baryon next to leading baryon
c     remaining anti-quark from leading baryon is marked by id+100
      IF((IABS(IFL(JR)).gt.100.and.IAFL2.gt.10).or.
     & (IABS(IFL(3)).gt.100.and.IAFL1.gt.10)) GOTO 200
      
      IF(IPAR(40).eq.0)THEN
c     reject two diquarks, two anti-diquarks AND diquark anti-diquark pairs
         IF (IAFL1*IAFL2 .GT. 100)  GOTO 200 
      ELSE
c     ONLY reject two diquarks or two anti-diquarks (unphysical) 
c     AND KEEP diquark anti-diquark pairs 
         IF (mod(IFL(JR),100)*mod(IFL(3),100).GT.100) GOTO 200 
      ENDIF

      IF ((IAFL1/10.eq.4.or.mod(IAFL1,10).eq.4)
     +     .and.(IAFL2/10.eq.4.or.mod(IAFL2,10).eq.4))
     +     GOTO 200             ! reject two charm quarks

C.... Vector meson config
c     increase vec.meson ration for diff.
      IF(IFQRK.eq.-1.and.IPAR(11).le.-4.and.IPAR(11).gt.-8) PAR(5) =9.D0
c     increase vec.meson ration for leading quarks in valence interactions
      IF(IABS(IFQRK).eq.1.and.IPAR(11).le.-5.and.ilead(jr).eq.1
     &     .and.IPAR(11).gt.-8) PAR(5) = 9.D0

c     suppress direct pi0 for meson projectiles
c     rate set by par( 137 )
 666  ipar82_def = IPAR(82)
c     skip if baryon projectile
      if(ibar(iabs(kb)).ne.0.or.ifqrk.eq.0) IPAR(82) = 0

c     suppress direct omega for meson projectiles
c     rate set by par( 138 )
      ipar83_def = IPAR(83)
c     skip if baryon projectile or central string     
      if(ibar(iabs(kb)).ne.0.or.(ifqrk.gt.0.and.IPAR(83).eq.2))
     &     IPAR(83) = 0

c     set current rank
      IRNK(JR)=ISIGN(IABS(IRNK(JR))+1,1-JR)
      
c     change rho0 / omega ratio
      IF(IPAR(81).eq.1)THEN
c     change if beam is meson
         if(ibar(iabs(kb)).eq.0) PAR(143) = PAR(144)
      ELSEIF(IPAR(81).eq.2)THEN         
c     change if beam is meson, on meson side only         
         if(ibar(iabs(kb)).eq.0.and.IRNK(JR).gt.0) PAR(143) = PAR(144)  
      ELSEIF(IPAR(81).eq.3)THEN         
c     change if beam is meson, on meson side only, for leading only
         if(ibar(iabs(kb)).eq.0.and.ISIGN(ILEAD(JR),IRNK(JR)).eq.1)
     &        PAR(143) = PAR(144)
      ELSEIF(IPAR(81).eq.4)THEN         
c     change if beam is meson, on meson side only, for diff. strings only
         if(ibar(iabs(kb)).eq.0.and.IFQRK.eq.-1)
     &        PAR(143) = PAR(144)
      ELSEIF(IPAR(81).eq.5)THEN         
c     change if beam is meson, for leading on meson side only and
c     for diff. strings only
         if(ibar(iabs(kb)).eq.0.and.IFQRK.eq.-1.and.
     &        ISIGN(ILEAD(JR),IRNK(JR)).eq.1) PAR(143) = PAR(144)
      ENDIF

c     increase vec.meson ratio for leading particle in str. diff.
      IF(IPAR(66).eq.1)THEN
         IF(ILEAD(JT).EQ.1.and.IFQRK.eq.-1)THEN
            IF(IBAR(IABS(kb)).eq.0.or.IPAR(70).eq.1) PAR(5) = PAR(113)
         ENDIF

      ELSEIF(IPAR(66).eq.2)THEN
         IF(IFQRK.eq.-1)THEN
            IF(IBAR(IABS(kb)).eq.0.or.IPAR(70).eq.1) PAR(5) = PAR(113)
         ENDIF

      ELSEIF(IPAR(66).eq.3)THEN
c     increase vector meson rate for meson beam
c     on beam side (rank+) only!
         IF(IFQRK.eq.-1)THEN
            IF(ILEAD(JR).EQ.1)THEN               
               IF(IBAR(IABS(kb)).eq.0.and.IRNK(JR).gt.0)
     &              PAR(5) = PAR(113)
c     always incr. vector rate for diff. strings independent of beam type
               IF(IPAR(70).eq.1) PAR(5) = PAR(113)               
            ENDIF
         ENDIF
      ENDIF

      CALL SIB_I4FLAV (IFL(JR), -IFL(3), IRNK(JR), IFLA, LLIST(I+1))

      IPAR(82) = ipar82_def
      IPAR(83) = ipar83_def
      PAR(143) = PAR143_def
      
      nporig(I+1)= Ipflag*2 + KINT
      niorig(I+1)= iiflag
      IF(ILEAD(1).eq.1.or.ILEAD(2).eq.1) nporig(I+1)= -1 * nporig(I+1)

c     replace leading or all pi0 with rho0
      IF(IFQRK.eq.-1) THEN 
         IF(IPAR(67).eq.1)THEN
            IF(ILEAD(JR).EQ.1) THEN 
               IF(IABS(LLIST(I+1)).EQ.6) 
     &              LLIST(I+1) = 27*isign(1,LLIST(I+1))           
            ENDIF
         ELSEIF(IPAR(67).eq.2)THEN
            IF(IABS(LLIST(I+1)).EQ.6) LLIST(I+1) =27*isign(1,LLIST(I+1))
         ELSEIF(IPAR(67).eq.3)THEN
            IF(IBAR(IABS(KB)).eq.0)THEN
               IF(ABS(LLIST(I+1)).EQ.6)LLIST(I+1)=27*isign(1,LLIST(I+1))
            ENDIF
         ENDIF
      ENDIF
      
c     replace all for diff.
      IF(IABS(IFQRK).EQ.1)THEN
         IF(IFQRK.eq.-1.and.IPAR(11).lt.0
     &        .and.IPAR(11).ge.-3) then
            IF(ABS(LLIST(I+1)).EQ.6) LLIST(I+1) = 27*isign(1,LLIST(I+1))
         endif
c     replace all for leading val.
         IF(IPAR(11).le.-2.and.IPAR(11).ge.-3) then
            if( ilead(jr).eq.1 ) then
               IF(IABS(LLIST(I+1)).EQ.6)
     &              LLIST(I+1) = 27*isign(1,LLIST(I+1))
            endif
         endif

c     increased vec.meson ratio and replace pi0 with rho0
         IF(IFQRK.eq.-1.and.IPAR(11).eq.-7) then
           IF(IABS(LLIST(I+1)).EQ.6) LLIST(I+1) = 27*isign(1,LLIST(I+1))
c     IF(ABS(LLIST(I+1)).EQ.7)  LLIST(I+1) = 25*isign(1,LLIST(I+1))
         endif
         
c     replace all for diff. ( same as lvec6 but for rhop as well )
c     reset vec.meson ratio
         IF(IFQRK.eq.-1.and.IPAR(11).eq.-8) then
            PAR(5) = 9.D0
            if( INT((PAR(5)+1.D0)*S_RNDM(0)).gt.1 ) then
               IF(IABS(LLIST(I+1)).EQ.6)
     &              LLIST(I+1) = 27*isign(1,LLIST(I+1))
               IF(IABS(LLIST(I+1)).EQ.7)
     &              LLIST(I+1) = 25*isign(1,LLIST(I+1))
            endif
         endif
c     replace leading pseudoscalar by vector
         IF(IPAR(11).eq.-8.and.ilead(jr).eq.1) then
            PAR(5) = 9.D0
            if( INT((PAR(5)+1.D0)*S_RNDM(0)).gt.1 ) then
               IF(IABS(LLIST(I+1)).EQ.6) 
     &              LLIST(I+1) = 27*isign(1,LLIST(I+1))
               IF(IABS(LLIST(I+1)).EQ.7)
     &              LLIST(I+1) = 25*isign(1,LLIST(I+1))
            endif
         endif
         
c     replace all pi0 for string diff.( same as lvec7 but for rhop as well )
         IF(IFQRK.eq.-1.and.IPAR(11).eq.-9) then
            if(IABS(LLIST(I+1)).EQ.6) LLIST(I+1) =27*isign(1,LLIST(I+1))
         endif
c     replace leading pi0 by vector
         IF(IPAR(11).eq.-9.and.ILEAD(JR).eq.1) then
            if(IABS(LLIST(I+1)).EQ.6) LLIST(I+1) =27*isign(1,LLIST(I+1))
         endif

c     replace leading omega in string diff by rho0's 
c     suppress pi0 in diff. strings
c     in addition to increased leading vec.meson ratio (lvec22)
         IF(IFQRK.eq.-1.and.IPAR(11).eq.-17)THEN
            IF(IABS(LLIST(I+1)).EQ.6)THEN
c     print*,'found pi0, restarting..'
               GOTO 666
            ENDIF
         ENDIF
         ILEAD(JR)= 0
      ENDIF

c     reject iso 0 spin 1 (omega) for meson projectiles
      IF(IBAR(IABS(KB)).eq.0)THEN
c     reject leading spin1,isospin singlett
         IF(ILEAD(JR).EQ.1.and.LLIST(I+1).eq.32.and.
     +        PAR(136).gt.S_RNDM(I)) LLIST(I+1) = 27
      endif
      
c     reset vec.mes. ratio
      PAR(5) = PAR5_def
      PAR(24) = charmPARdef
      IPAR(11) = IPAR11_def

      P(I,1)   = PX(JT)+PX(3)      
      P(I,2)   = PY(JT)+PY(3)
      LPOINT(I) = JT
      I1 = I+1
      nforig(I1) = 0      
      P(I1,5) = AM(IABS(LLIST(I1)))
      P(I1,1) = PX(JR)-PX(3)      
      P(I1,2) = PY(JR)-PY(3)   
      LPOINT(I1) = JR 
      LRNK(I1) = IRNK(JR)
      XM1 = P(I,5)**2+P(I,1)**2+P(I,2)**2
      XM2 = P(I1,5)**2+P(I1,1)**2+P(I1,2)**2
      IF (DSQRT(XM1)+dSQRT(XM2) .GT. dSQRT(WREM2)) GOTO 200

c...RE & EJA fix
      PT2 = (P(I,1)+P(I1,1))**2+(P(I,2)+P(I1,2))**2
      WREMPT = dsqrt(WREM2+PT2)
      EA1 = (WREM2+XM1-XM2+PT2)/(2.D0*WREMPT)

      PA2 = (EA1**2-XM1)
      if (pa2.gt.0.D0)  then
            PA = dSQRT(PA2)
      else
            goto 200
      endif
      BA = PTOT(3)/PTOT(4)
      GA = PTOT(4)/WREMPT
      SGN = DBLE(3-2*JT)
      P(I,3) = GA*(BA*EA1+SGN*PA)
      P(I,4) = GA*(EA1+BA*SGN*PA)
      P(I+1,3) = PTOT(3)-P(I,3)
      P(I+1,4) = PTOT(4)-P(I,4)

c     mark as final hadrons
      ZLIST(I) = 0.D0
      ZLIST(I+1) = 0.D0

      NA= NP+1
      NP=I+1
         
C...reorder  particles along chain (in rank)
      IF (LRANK)  THEN
      N1 = NA-1
      N2 = 0
      DO J=NA,NP
         IF(P(J,4).lt.0) THEN
            NP=NA-1
            GOTO 200            ! negative energy bug 'fix'
         ENDIF
         IF(LPOINT(J) .EQ. 2)  THEN
            N2=N2+1
            LLIST (NP+N2) = LLIST(J)
            LRNK(NP+N2) = LRNK(J)
            ZLIST (NP+N2) = ZLIST(J)
            nporig(NP+N2) = nporig(J)
            niorig(NP+N2) = niorig(J)
            nforig(NP+N2) = 0
            DO K=1,5
               P(NP+N2,K)=P(J,K)
            ENDDO
         ELSE
            N1= N1+1
            IF (N1.LT.J)   THEN
               LLIST(N1) = LLIST(J)
               LRNK(N1) = LRNK(J)
               ZLIST(N1) = ZLIST(J)
               nporig(N1) = nporig(J)
               niorig(N1) = niorig(J)
               nforig(N1) = nforig(J)
               DO K=1,5
                  P(N1,K) = P(J,K)
               ENDDO
            ENDIF
         ENDIF
      ENDDO
      JJ=N1
      DO J=NP+N2,NP+1,-1
         JJ= JJ+1
         LLIST(JJ) = LLIST(J)
         LRNK(JJ) = LRNK(J)
         ZLIST(JJ) = ZLIST(J)
         nporig(JJ) = nporig(J)
         niorig(JJ) = niorig(J)
         nforig(JJ) = nforig(J)
         DO K=1,5
            P(JJ,K) = P(J,K)
         ENDDO
      ENDDO
      ENDIF

      if(Ndebug.gt.3)
     &  WRITE(LUN,*)' STRING_FRAG_4FLV: NP after fragmentation:',NP

      END


C=======================================================================

      SUBROUTINE GG_FRAG_4FLV (E0)

C-----------------------------------------------------------------------
C...This routine fragments a  gluon-gluon system
C.  of mass E0 (GeV)
C.  the particles produced are in the  jet-jet frame
C.  oriented along the z axis
C...........................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
C     The final particle output is contained in COMMON /S_PLIST/    
C     NP           : number of final particles
C     P(1:NP, 1:5) : 4-momenta + masses of the final particles 
C     LLIST (1:NP) : codes of final particles
      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP
      DOUBLE PRECISION ZLIST
      COMMON /S_ZLIST/ ZLIST(8000)
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)

      INTEGER NFORIG,NPORIG,NIORIG,IPFLAG,IIFLAG,KINT
      COMMON /S_PARTO/ NFORIG(NP_max),NPORIG(NP_max),NIORIG(NP_max),
     &IPFLAG,IIFLAG,KINT

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN

      DIMENSION WW(2,2),PTOT(4),PX(3),PY(3),IFL(3),PMQ(3)
      SAVE

      if(Ndebug.gt.3) then
        WRITE(LUN,*)
     &    ' GG_FRAG_4FLV: called with (E0)',
     &    E0
        WRITE(LUN,*)' GG_FRAG_4FLV: NP before fragmentation:',NP
      endif

C...  'leading' strange fraction
      PAR2_def = PAR(2)
      IF(IPAR(39).eq.2) PAR(2) = PAR(66)

      PAR24_def = PAR(24)
C     leading charm fraction
      IF(IPAR(87).eq.1) PAR(24) = PAR(150)
      IF(IPAR(87).eq.2) PAR(24) = PAR(150)*PAR(24)

      E0S = E0**2
      
C...Generate the 'forward' leading particle.
100   I = NP+1
c     dummy rank argument
      IDM = 5
c     sample new flavor, i.e. split gluon into quark-antiquark, quark or antiquark
      if( IPAR(87).eq.3 )THEN
C     flavor threshold model            
c     u,d -> u,d,s -> u,d,s,c
         CALL SIB_ICFLAV(E0S,0,I0,IFL1)
      ELSE
c     default u,d,s model, same rates as in hadronization (string frag.)
         I0 = INT(-1 + 2.D0*INT((2.D0-EPS8)*S_RNDM(I)))               
         CALL SIB_I4FLAV(I0,0,IDM,IFL1, LDUM)
      ENDIF
c     form first hadron from new flavor
      CALL SIB_I4FLAV(IFL1,0,IDM,IFL2, LLIST(I))
      CALL PTDIS_4FLV(IFL1,PX1,PY1)
      CALL PTDIS_4FLV(IFL2,PX2,PY2)
      P(I,1) = PX1+PX2
      P(I,2) = PY1+PY2
      P(I,5) = AM(IABS(LLIST(I)))
      XM1 = P(I,5)**2+P(I,1)**2+P(I,2)**2
      Z1 = ZDIS_4FLV (IFL1,1,0.25D0*XM1)
      Z2 = ZDIS_4FLV (IFL2,1,0.25D0*XM1)
      T1  = 4.D0*XM1/(E0S*(Z1+Z2))
      P(I,4) = 0.25D0*E0*(Z1+Z2 + T1)
      P(I,3) = 0.25D0*E0*(Z1+Z2 - T1)

      nforig(I)= 0
      nporig(I)= Ipflag*3 + KINT
      niorig(I)= iiflag
      ZLIST(I) = Z1 + Z2

C...Generate the 'backward' leading particle.
      I = I+1
      IF( IPAR(87).eq.3 )THEN
         CALL SIB_ICFLAV(E0S,-I0,IDM,IFL3)
      ELSE
         CALL SIB_I4FLAV(-I0,0,IDM,IFL3, LDUM)
      ENDIF
      CALL SIB_I4FLAV(IFL3,0,IDM,IFL4, LLIST(I))
      CALL PTDIS_4FLV(IFL3,PX3,PY3)
      CALL PTDIS_4FLV(IFL4,PX4,PY4)
      P(I,1) = PX3+PX4
      P(I,2) = PY3+PY4
      P(I,5) = AM(IABS(LLIST(I)))
      XM2 = P(I,5)**2+P(I,1)**2+P(I,2)**2
      Z3 = ZDIS_4FLV (IFL3,1,0.25D0*XM2)
      Z4 = ZDIS_4FLV (IFL4,1,0.25D0*XM2)
      T2  = 4.D0*XM2/(E0S*(Z3+Z4))
      P(I,4) = 0.25D0*E0*( Z3+Z4 + T2)
      P(I,3) = 0.25D0*E0*(-Z3-Z4 + T2)

      nforig(I)= 0
      nporig(I)= Ipflag*3 + KINT
      niorig(I)= iiflag
      ZLIST(I) = Z3 + Z4
c      PAR24def = PAR(24)
c     charm QCD fusion
c      IF(IPAR(17).eq.2) PAR(24) = 0.

c     reset strange fraction
      PAR(2) = PAR2_def

c     reset charm fraction
      PAR(24) = PAR24_def
      
C...Fragment the two remaning strings
      N0 = 0
      DO KS=1,2
      
      NTRY = 0
200   NTRY = NTRY+1
      I = NP+2+N0
      IF (NTRY .GT. 30)  GOTO 100

      IF (KS .EQ. 1)  THEN
         WW(1,1) = 0.5D0 * (1.D0 - Z1 - 0.5D0*T2) 
         WW(2,1) = 0.5D0 * (1.D0 - Z3 - 0.5D0*T1)
         PX(1) = -PX1
         PY(1) = -PY1
         PX(2) = -PX3
         PY(2) = -PY3
         IFL(1) = -IFL1
         IFL(2) = -IFL3
      ELSE
         WW(1,1) = 0.5D0 * (1.D0 - Z2 - 0.5D0*T2) 
         WW(2,1) = 0.5D0 * (1.D0 - Z4 - 0.5D0*T1)
         PX(1) = -PX2
         PY(1) = -PY2
         PX(2) = -PX4
         PY(2) = -PY4
         IFL(1) = -IFL2
         IFL(2) = -IFL4
      ENDIF
      PX(3) = 0.D0
      PY(3) = 0.D0
      PTOT (1) = PX(1)+PX(2)
      PTOT (2) = PY(1)+PY(2)
      PTOT (3) = 0.5D0*E0*(WW(1,1)-WW(2,1))
      PTOT (4) = 0.5D0*E0*(WW(1,1)+WW(2,1))

      PMQ(1) = QMASS(IFL(1))
      PMQ(2) = QMASS(IFL(2))

C...produce new particle: side, pT
300   I=I+1
      if(i.gt.8000) then
        write(LUN,'(1x,a,i8)') 
     &    ' GG_FRAG: no space left in S_PLIST:',I
        CALL SIB_REJECT ('GG_FRAG         ')
      endif
      nforig(I)= 0
      nporig(I)= Ipflag*2 + KINT
      niorig(I)= iiflag

      JT=INT(1.5D0+S_RNDM(0))
      JR=3-JT
c      CALL PTDIS (IFL(JT), PX(3),PY(3))

C...particle ID
      CALL SIB_I4FLAV (IFL(JT), 0, IDM, IFL(3), LLIST(I))
      PMQ(3) = QMASS(IFL(3))
      P(I,5) = AM(IABS(LLIST(I)))

      CALL PTDIS_4FLV (IFL(3), PX(3),PY(3))
      
C...test end of fragmentation
      WREM2 = PTOT(4)**2-PTOT(1)**2-PTOT(2)**2-PTOT(3)**2
      IF (WREM2 .LT. 0.1D0)  GOTO 200
      WMIN = PMQ(1)+PMQ(2)+2.D0*PMQ(3)+1.1D0+(2.D0*S_RNDM(0)-1.D0)*0.2D0
      IF (WREM2 .LT. WMIN**2)THEN
         GOTO 400
      ENDIF

C...fill transverse momentum
      P(I,1) = PX(JT) + PX(3)
      P(I,2) = PY(JT) + PY(3)

C...Choose z
      XMT2 = P(I,5)**2+P(I,1)**2+P(I,2)**2
      Z = ZDIS_4FLV (ifl(3),IFL(JT), XMT2)

      ZLIST(I) = Z      
      WW(JT,2) = Z*WW(JT,1)
      WW(JR,2) = XMT2/(WW(JT,2)*E0S)

      P(I,3) = WW(1,2)*0.5D0*E0 - WW(2,2)*0.5D0*E0
      P(I,4) = WW(1,2)*0.5D0*E0 + WW(2,2)*0.5D0*E0

      DO J=1,4
         PTOT (J) = PTOT(J) - P(I,J)
      ENDDO
      DO K=1,2
         WW(K,1) = WW(K,1) - WW(K,2)
      ENDDO

C...Reset pT and flavor at ends of the string
      PX(JT) = -PX(3)
      PY(JT) = -PY(3)
      IFL(JT) =-IFL(3)
      PMQ(JT) = PMQ(3)
      GOTO 300

C...Final two hadrons
 400  IAFL1 = mod(IABS(IFL(JR)),100)
      IAFL2 = mod(IABS(IFL(3)),100)
      IF (IAFL1*IAFL2 .GT. 100)  GOTO 200 ! reject two diquarks
      IF ((IAFL1/10.eq.4.or.mod(IAFL1,10).eq.4)
     +     .and.(IAFL2/10.eq.4.or.mod(IAFL2,10).eq.4))
     +     GOTO 200             ! reject two charm quarks

      CALL SIB_I4FLAV (IFL(JR), -IFL(3), IDM, IFLA, LLIST(I+1))
      P(I+1,5) = AM(IABS(LLIST(I+1)))
      P(I,1)   = PX(JT)+PX(3)      
      P(I,2)   = PY(JT)+PY(3)      
      nporig(I)= Ipflag*2 + KINT
      niorig(I)= iiflag
      I1 = I+1
      nporig(I1)= Ipflag*2 + KINT
      niorig(I1)= iiflag
      P(I1,1) = PX(JR)-PX(3)      
      P(I1,2) = PY(JR)-PY(3)      
      XM1 = P(I,5)**2+P(I,1)**2+P(I,2)**2
      XM2 = P(I1,5)**2+P(I1,1)**2+P(I1,2)**2
      IF (dSQRT(XM1)+dSQRT(XM2) .GT. dSQRT(WREM2)) GOTO 200
      if (ptot(4).le.0.D0) goto 200
      PT2 = (P(I,1)+P(I1,1))**2+(P(I,2)+P(I1,2))**2
      WREMPT = dsqrt(WREM2+PT2)
      EA1 = (WREM2+XM1-XM2+PT2)/(2.D0*WREMPT)
      PA2 = (EA1**2-XM1)
      if (PA2.ge.0.D0) then
        PA = dSQRT(PA2)
      else
         goto 200
      endif
      BA = PTOT(3)/PTOT(4)
      GA = PTOT(4)/WREMPT
      SGN = DBLE(3-2*JT)
      P(I,3) = GA*(BA*EA1+SGN*PA)
      P(I,4) = GA*(EA1+BA*SGN*PA)
      P(I+1,3) = PTOT(3)-P(I,3)
      P(I+1,4) = PTOT(4)-P(I,4)
      ZLIST(I) = 0.D0
      ZLIST(I+1) = 0.D0
      N0 = I-NP-1
      ENDDO                  ! loop on two `remaining strings'

      NP = I+1
c      PAR(24) = PAR24def 
      IF(Ndebug.gt.3) then
        WRITE(LUN,*)' GG_FRAG_4FLV: NP after fragmentation:',NP
      ENDIF
      RETURN
      END

C=======================================================================

      SUBROUTINE DIFDEC (L0, Irec, IBAD, P0)

C-----------------------------------------------------------------------
C..."decay" of an excited state with the quantum numbers
C.   of particle L0 and the 5-momentum P0
C.   - low energy: phase space decay (fire ball model)
C.   - intermediate energy: one-string decay (longitudinal phase space)
C.   - high energy: pomeron-hadron scattering (multi-string model) 
C-----------------------------------------------------------------------
      IMPLICIT NONE

c     external types
      INTEGER L0, Irec, IBAD
      DOUBLE PRECISION P0
      DIMENSION P0(5)
      
      INTEGER NW_max
      PARAMETER (NW_max = 20)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
C     The final particle output is contained in COMMON /S_PLIST/    
C     NP           : number of final particles
C     P(1:NP, 1:5) : 4-momenta + masses of the final particles 
C     LLIST (1:NP) : codes of final particles
      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP

      INTEGER NFORIG,NPORIG,NIORIG,IPFLAG,IIFLAG,KINT
      COMMON /S_PARTO/ NFORIG(NP_max),NPORIG(NP_max),NIORIG(NP_max),
     &IPFLAG,IIFLAG,KINT

      INTEGER LRNK
      COMMON /SIB_RNK/ LRNK(8000)

      INTEGER ICHP,ISTR,IBAR
      COMMON /S_CHP/ ICHP(99), ISTR(99), IBAR(99)

      INTEGER IISO,ISPN
      COMMON /S_SPN/ IISO(99), ISPN(99)

      INTEGER ICHM
      COMMON /S_CHM/ ICHM(99)
C--------------------------------------------------------------------
C     SIBYLL common blocks containing event information       \FR'14
C--------------------------------------------------------------------

C     EVENT INFO COMMON
C     contains overall interaction properties, like
C     SQS : center-of-mass energy
C     S   :         "       "     squared
C     PTmin : low pt cut of QCD cross section, 
C             i.e. minimal pt of hard minijets
C     Xmin : low-x bound for PDFs, 
C            i.e. minimal momentum fraction of hard partons
C     Zmin : logarithm of that
C     KB : PID of beam hadron
C     KT() : PID of target
C     IAT : mass number of target
      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
      INTEGER KB,IAT,KT
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN

c     internal types
      INTEGER LL,LCON,LRES,LRES1,NTRYS,NRJECT,LA,N1,IREJ,I,J,IFLA,
     &     IFL1,IFL2,IFBAD,NPI,IRES,LA1,JQQ,JQTOT,K,JQR,
     &     KB_0,IAT_0
      DOUBLE PRECISION PD,BE,EMIN,EMIN2,PCHEX,PRES,DELTAE,
     &     SQS_0,S_0,PTmin_0,XMIN_0,ZMIN_0,
     &     PAR1_def,PAR24_def,PAR53_def,GA,BEP,S_RNDM,AV,GASDEV,PCXG,
     &     XI1,XI2,XSMR         !,FERMI
      DIMENSION LL(10), PD(10,5), BE(3), LCON(6:99),LRES1(6:99)
      DIMENSION LRES(6:99)
      SAVE
      EXTERNAL GASDEV
      DATA (LRES(k),k=6,22)  /27,25,26,28,29,0,0,51,52,6*0,30,31/
      DATA (LRES(k),k=23,33) /23,24,25,26,27,28,29,30,31,27,27/
      DATA (LRES(k),k=34,49) /34,35,36,37,38,39,40,41,42,43,34,35,36,
     &     37,38,49/
      DATA (LRES(k),k=50,83) /0,51,52,53,54,4*0,78,79,10*0,80,81,73,
     &     74,75,76,77,78,79,80,81,0,83/
      DATA (LRES(k),k=84,99) /94,95,96,97,98,89,4*0,94,95,96,97,98,99/
      
      DATA EMIN /0.7D0/
      DATA EMIN2 /10.D0/
      DATA LCON /7,6,6,11,11,9,9,14,13,19*0,35,34,35,38,37,39,
     &     19*0,71,72,10*0,59,60,73,10*0,85,86,85,88,87,89,10*0/
      DATA LRES1 /27,25,26,11,11,9,9,14,13,19*0,35,34,35,38,37,39,
     &     19*0,78,79,10*0,80,81,83,10*0,94,95,96,97,98,89,10*0/      
      DATA PCHEX /0.33D0/            ! probability of charge exchange
      DATA PRES /0.7D0/         ! probability of forming a resonance
      DATA NRJECT /0/

      IF(NDEBUG.gt.2)
     &     WRITE(LUN,'(2X,A,1x,I2,1x,I2,/,5(2x,F10.3))')
     &     'DIFDEC: (L0,Irec,P0):',L0,Irec,(P0(i),i=1,5)
      
      
      NTRYS = 0

      LA = IABS(L0)
      DELTAE = P0(5) - AM(LA)
      IF(IBAR(LA).ne.0.or.IPAR(65).eq.0)THEN
c     baryons
         EMIN = PAR(30)
      ELSE
c     mesons
         EMIN = PAR(112)
      ENDIF
c      IBAD = 0
      PAR1_def= PAR(1)
      if(Irec.gt.0) PAR(1)= PAR(16)
c      XSMR = 0.5D0
c     XI2=FERMI(DELTAE,EMIN2,XSMR)
c     XI1=FERMI(DELTAE,EMIN,XSMR)
      XSMR=PAR(131)*EMIN
      XI1=MAX((EMIN-DELTAE)/XSMR,0.D0)      
      XSMR=PAR(131)*EMIN2
      XI2=MAX((EMIN2-DELTAE)/XSMR,0.D0)
      if(Ndebug.gt.2) 
     &     WRITE(LUN,'(1x,A29,2(2x,F5.2),2(2x,F8.3))')
     &     '  DIFDEC: EMIN1,EMIN2,XI1,XI2',
     &     EMIN,EMIN2,Xi1,Xi2
      
C...  pomeron-hadron scattering (pi0 is used instead of pomeron)      
      IF ((IPAR(10).gt.0).and.(Irec.gt.0).and.
     &     (DELTAE.gt.EMIN2.or.S_RNDM(LA).gt.XI2))  THEN
         if(Ndebug.gt.2) 
     &        WRITE(LUN,*)' DIFDEC: central (L0,DELTAE,NP,XI):',
     &        L0,DELTAE,NP,XI2
         N1 = NP+1
         if(irec.gt.0.and.IPAR(5).eq.1) par(1)= par(15)
 50      CONTINUE
         IPFLAG= IPFLAG*100
c     create subevent
         SQS_0   = SQS
         S_0     = S
         PTmin_0 = PTmin 
         XMIN_0  = XMIN
         ZMIN_0  = ZMIN
         KB_0    = KB
         IAT_0   = IAT
         CALL INI_EVENT(P0(5),L0,6,0)
c     create L0 - pi0 interaction, pi0(pid=6) target
         CALL SIB_NDIFF(L0, 1, P0(5), 0, IREJ) ! ori
c     restore main event
         SQS   = SQS_0
         S     = S_0
         PTmin = PTmin_0         
         XMIN  = XMIN_0
         ZMIN  = ZMIN_0
         KB    = KB_0
         IAT   = IAT_0
         IF(IREJ.NE.0) THEN
            NP = N1-1
            GOTO 50
         ENDIF
         PAR(1) = PAR1_def
         DO J=1,3
            BE(J)=P0(J)/P0(4)
         ENDDO
         GA=P0(4)/P0(5)
         if(P0(3).lt.0.D0) then
           do i=N1,NP
             P(I,3) = -P(I,3)
           enddo
         endif
         DO I=N1,NP
            BEP=BE(1)*P(I,1)+BE(2)*P(I,2)+BE(3)*P(I,3)
            DO J=1,3
               P(I,J)=P(I,J)+GA*(GA*BEP/(1.D0+GA)+P(I,4))*BE(J)
            ENDDO
            P(I,4)=GA*(P(I,4)+BEP)
         ENDDO

C..."string-like" decay
      ELSE IF (DELTAE .GT. EMIN .or. S_RNDM(LA).gt.XI1)  THEN          
         IF(NDEBUG.gt.2) 
     &        WRITE(LUN,'(2X,A,3(2x,F8.3))')
     &        'DIFDEC: string-like, (DELTAE,E0,central prob.):',
     &        DELTAE,P0(5),1.D0-XI2
c     set charge exchange probability, i.e. prob for p* -> n + pip
         PAR53_def = PAR(53)
         PAR(53) = PAR(130)
         N1 = NP+1
         CALL HSPLI(L0,IFL1,IFL2)
         PAR(53) = PAR53_def
         IF (P0(3) .GT. 0.D0.and.L0.gt.0)  THEN
            IFLA = IFL2
            IFL2 = IFL1
            IFL1 = IFLA
         ENDIF
c     randomize flavor orientation in string
         IF(IPAR(25).eq.1.and.S_RNDM(L0).gt.PAR(39))THEN
            IFLA = IFL2
            IFL2 = IFL1
            IFL1 = IFLA
         ENDIF
         PAR24_def = PAR(24)
         IF(IPAR(15).eq.2.and.IPAR(15).eq.3.and.IPAR(15).ne.7.and.
     &        IPAR(15).lt.12)THEN
         PAR(24) = PAR(25)*dEXP(-PAR(26)/P0(5))
         ELSEIF(IPAR(15).eq.7)THEN
            PAR(24) = PAR(25)
         ENDIF
 10      CONTINUE
         IPFLAG = IPFLAG*10
         CALL STRING_FRAG_4FLV 
     +        (P0(5), IFL1, IFL2, 0.D0,0.D0,0.D0,0.D0,IFBAD,-1)
         IF (IFBAD .EQ. 1)then
            if(ndebug.gt.1)
     &           WRITE(lun,*)' SIB_DIFF: string-frag rejection! ',
     &           '(M,NCALL)',P0(5),NCALL
            NTRYS = NTRYS + 1
            NP = N1-1
            IFBAD = 0
            IF(NTRYS.gt.5)then ! resample diff. mass              
               NP = 0
               IBAD = 1
               PAR(24) = PAR24_def
               RETURN
            endif
            GOTO 10
         ENDIF
         DO J=1,3
            BE(J)=P0(J)/P0(4)
         ENDDO
         GA=P0(4)/P0(5)
         DO I=N1,NP
            BEP=BE(1)*P(I,1)+BE(2)*P(I,2)+BE(3)*P(I,3)
            DO J=1,3
               P(I,J)=P(I,J)+GA*(GA*BEP/(1.D0+GA)+P(I,4))*BE(J)
            ENDDO
            P(I,4)=GA*(P(I,4)+BEP)
         ENDDO
         PAR(24) = PAR24_def

C...Phase space decay of the excited state
      ELSEIF(DELTAE.GT.AM(7)+0.02D0)THEN
         if(Ndebug.gt.2) 
     &        WRITE(LUN,*)' DIFDEC: fireball, (DELTAE,string prob.):',
     &        DELTAE,1.D0-XI1
         IF(IPAR(14).GT.0.and.IPAR(14).NE.7)THEN
            IF(IPAR(14).eq.5) PCHEX = 0.D0
            NPI=0
            IRES = 0
            IF (S_RNDM(0).LT.PRES) THEN
               IF (LA.LT.9) THEN
c     if kinematically possible produce rho0 in charge exchange
                  LL(1) = LRES(LA)
                  DELTAE = P0(5) -  AM(LRES(LA))
                  IF (DELTAE.GT.AM(7)+0.02D0) GOTO 100
               ENDIF
            ENDIF
c     switch charge exchange on/off
            IF( S_RNDM(1).LT.PCHEX)THEN
               LL(1) = LCON(LA)*ISIGN(1,L0)
               IF( (L0 .EQ. 6) .OR. (L0 .EQ. 11) )
     .              LL(1) = LL(1)+INT((2.D0-EPS8)*S_RNDM(2))
            ELSE
               LL(1) = L0
            ENDIF
            
            DELTAE = P0(5) - AM(LA)
 100        AV = 2.D0*dSQRT(DELTAE)
            LA1 = IABS(LL(1))
            NPI = INT(AV*(2.D0+0.5D0*GASDEV(LA)))
            IF (IPAR(14).EQ.6)THEN
               IF(NPI.LT.1.OR.NPI.GT.9.OR.AM(LA1)+NPI*AM(7)+0.02D0
     .              .GT.P0(5))  GOTO 100
            ELSE
               IF(NPI.LT.0.OR.NPI.GT.9.OR.AM(LA1)+NPI*AM(7)+0.02D0
     .              .GT.P0(5))  GOTO 100
            ENDIF
c     create resonances inside fireball..
            IF(IPAR(14).ge.2
     +           .and.DELTAE.GE.AM(LA1)+AM(27)+(NPI-1)*AM(7)+0.02D0)
     +           IRES = 1
            IF(IPAR(14).ge.3.and.DELTAE.GE.AM(LA1)+NPI*AM(27)+0.02D0) 
     +           IRES=3
            JQQ = ICHP(LA)*ISIGN(1,L0)-
     .           ICHP(IABS(LL(1)))*ISIGN(1,LL(1))  
 120        JQTOT = 0
            DO K=2,NPI
               LL(K) = 6+INT(S_RNDM(K)*(3.D0-EPS8))
c     suppress pi0 in fireball
               IF(IPAR(14).ge.4)
     +              LL(K) = 7+INT(S_RNDM(0)*(2.D0-EPS8))
c     IF(IRES.EQ.1.and.S_RNDM(LA).LT.0.5D0)
               IF(IRES.EQ.1) THEN
                  LL(K) = 27-INT(S_RNDM(1)*(3.D0-EPS8))
                  IRES = 2
               ENDIF
               IF(IRES.EQ.3)
     +              LL(K) = 27-INT(S_RNDM(2)*(3.D0-EPS8))
               JQTOT = JQTOT + ICHP(LL(K))
            ENDDO
            JQR = JQQ-JQTOT
            IF (JQR.LT.-1.OR.JQR.GT.1)  GOTO 120
            LL(NPI+1) = 6+JQR
            IF (LL(NPI+1) .EQ. 5)  LL(NPI+1)=8
            CALL DECPAR (0,P0,NPI+1,LL, PD)
            DO J=1,NPI+1
               NP = NP+1
               LLIST(NP) = LL(J)
               nporig(NP)= Ipflag*2
               lrnk(Np) = 0
               niorig(NP)= iiflag
               DO K=1,5
                  P(NP,K) = PD(J,K)
               ENDDO
            ENDDO

         ELSEIF (IPAR(14).EQ.7.AND.LA.LT.9) THEN
c     all diff states go to resonances for pi beam ..
            NPI=0
            IRES = 0
            LL(1) = LRES1(LA)
            DELTAE = P0(5) -  AM(LL(1))
cdh         IF( DELTAE.LT.AM(7)+0.02D0) GOTO 222
            IF( DELTAE.LT.AM(7)+0.02D0) THEN
              IF(IPAR(14).EQ.7)  DELTAE = P0(5) - AM(LA)
              AV = 2.D0*DSQRT(DELTAE)
 201          NPI = INT(AV*(1.D0+0.5D0*GASDEV(LA)))
c              print *,'npi:',npi,'av',av,'p05',p0(5),am(la),deltae
              IF(NPI.LE.0.OR.NPI.GT.9.OR.AM(LA)+NPI*AM(7)+0.02D0
     .             .GT.P0(5))  GOTO 201
              IF (S_RNDM(0).LT.PCHEX)  THEN
                 LL(NPI+1) = LCON(LA)*ISIGN(1,L0)
                 IF( (L0 .EQ. 6) .OR. (L0 .EQ. 11) )
     .                LL(NPI+1) = LL(NPI+1)+INT((2.D0-EPS8)*S_RNDM(1))
              ELSE
                 LL(NPI+1) = L0
              ENDIF
              JQQ = ICHP(LA)*ISIGN(1,L0)-
     .             ICHP(IABS(LL(NPI+1)))*ISIGN(1,LL(NPI+1))
 221          JQTOT = 0
              DO K=1,NPI-1
                 LL(K) = 6+INT(S_RNDM(K)*(3.D0-EPS8))
                 JQTOT = JQTOT + ICHP(LL(K))
              ENDDO
              JQR = JQQ-JQTOT
              IF (JQR.LT.-1.OR.JQR.GT.1)  GOTO 221
              LL(NPI) = 6+JQR
              IF (LL(NPI) .EQ. 5)  LL(NPI)=8
              CALL DECPAR (0,P0,NPI+1,LL, PD)
              DO J=1,NPI+1
                 NP = NP+1
                 LLIST(NP) = LL(J)
                 NPORIG(NP) = IPFLAG*2
                 lrnk(Np) = 0
                 niorig(NP)= iiflag
                 DO K=1,5
                    P(NP,K) = PD(J,K)
                 ENDDO
              ENDDO

            ELSE
              IF( S_RNDM(0).LT.PAR(31))THEN
                 LL(1) = LRES1(LCON(LA))
                 IF( (L0 .EQ. 6) .OR. (L0 .EQ. 11) )
     .                LL(1) = LRES1(IABS(L0)+INT((2.D0-EPS8)*S_RNDM(1)))
              ENDIF
 300          AV = 2.D0*dSQRT(DELTAE)
              LA1 = IABS(LL(1))
              NPI = INT(AV*(2.D0+0.5D0*GASDEV(LA)))
              IF(ABS(PAR(32)).gt.0.D0)
     &             NPI = INT(AV*(PAR(32)+0.5D0*GASDEV(LA)))
              IF(NPI.LT.0.OR.NPI.GT.9.OR.AM(LA1)+NPI*AM(7)+0.02D0
     .             .GT.P0(5))  GOTO 300
c     create resonances inside fireball..
c              IRES=3
              JQQ = ICHP(LA)*ISIGN(1,L0)-
     .             ICHP(IABS(LL(1)))*ISIGN(1,LL(1))
 320          JQTOT = 0
              DO K=2,NPI
                 LL(K) = 6+INT(S_RNDM(K)*(3.D0-EPS8))
c     suppress pi0 in fireball
c                 IF(IPAR(14).ge.4)
c     +                LL(K) = 7+INT(S_RNDM(0)*1.99999D0)
c       IF(IRES.EQ.1.and.S_RNDM(LA).LT.0.5D0)
c                 LL(K) = 27-INT(S_RNDM(0)*2.99999D0)
                 JQTOT = JQTOT + ICHP(LL(K))
              ENDDO
              JQR = JQQ-JQTOT
              IF (JQR.LT.-1.OR.JQR.GT.1)  GOTO 320
              LL(NPI+1) = 6+JQR
              IF (LL(NPI+1) .EQ. 5)  LL(NPI+1)=8
              CALL DECPAR (0,P0,NPI+1,LL, PD)
              DO J=1,NPI+1
                 NP = NP+1
                 LLIST(NP) = LL(J)
                 nporig(NP)= Ipflag*2
                 lrnk(Np) = 0
                 niorig(NP)= iiflag
                 DO K=1,5
                    P(NP,K) = PD(J,K)
                 ENDDO
              ENDDO
            ENDIF

         ELSEIF (IPAR(14).LE.-1) THEN
C...  generalized fireball model
            IF(Ndebug.gt.2) 
     &           WRITE(LUN,*)' DIFDEC: using generalized fireball!'
c     set charge exchange probability, 
c     i.e. prob for p* -> n + pip
            PCXG = PAR(61)
            CALL FIREBALL_4FLV(L0,P0,PCXG,IFBAD)
            IF(IFBAD.eq.1)THEN
               IF(ndebug.gt.0)THEN
                  IF(NRJECT.le.10)THEN
                     WRITE(LUN,*)
     &                    ' DIFDEC: warning: fireball rejection! ',
     &                    'diff. mass to low to dissociate beam!'
                     WRITE(LUN,*)
     &                ' DIFDEC: m_Beam, DELTAE ,AM(7)+0.02, NCALL: ', 
     &                AM(LA),DELTAE,'>',AM(7)+0.02D0,NCALL
                  ENDIF
                  IF(NRJECT.eq.10) 
     &             write(lun,*)' this was the last warning.. good luck!'
               ENDIF
               NRJECT = NRJECT + 1
               NP = 0
               IBAD = 1
               RETURN
            ENDIF

         ELSE
cdh 222       IF(IPAR(14).EQ.7)  DELTAE = P0(5) - AM(LA)
           IF(IPAR(14).EQ.7)  DELTAE = P0(5) - AM(LA)
            AV = 2.D0*dSQRT(DELTAE)
 200        NPI = INT(AV*(1.D0+0.5D0*GASDEV(0)))
c            print *,'npi:',npi,'av',av,'p05',p0(5),am(la),deltae
            IF(NPI.LE.0.OR.NPI.GT.9.OR.AM(LA)+NPI*AM(7)+0.02D0
     .           .GT.P0(5))  GOTO 200
            IF (S_RNDM(0).LT.PCHEX)  THEN
               LL(NPI+1) = LCON(LA)*ISIGN(1,L0)
               IF( (L0 .EQ. 6) .OR. (L0 .EQ. 11) )
     .              LL(NPI+1) = LL(NPI+1)+INT((2.D0-EPS8)*S_RNDM(1))
            ELSE
               LL(NPI+1) = L0
            ENDIF
            JQQ = ICHP(LA)*ISIGN(1,L0)-
     .           ICHP(IABS(LL(NPI+1)))*ISIGN(1,LL(NPI+1))  
 220        JQTOT = 0
            DO K=1,NPI-1
               LL(K) = 6+INT(S_RNDM(K)*(3.D0-EPS8))
               JQTOT = JQTOT + ICHP(LL(K))
            ENDDO
            JQR = JQQ-JQTOT
            IF (JQR.LT.-1.OR.JQR.GT.1)  GOTO 220
            LL(NPI) = 6+JQR
            IF (LL(NPI) .EQ. 5)  LL(NPI)=8
            CALL DECPAR (0,P0,NPI+1,LL, PD)
            DO J=1,NPI+1
               NP = NP+1
               LLIST(NP) = LL(J)
               NPORIG(NP) = IPFLAG*2
               lrnk(Np) = 0
               niorig(NP)= iiflag
               DO K=1,5
                  P(NP,K) = PD(J,K)
               ENDDO
            ENDDO
         ENDIF
      ELSE
         if (ndebug .gt. 0) then
           IF(NRJECT.le.10)THEN
            WRITE(LUN,*) ' DIFDEC rejection! ',
     &           'diff. mass to low to dissociate beam!'
            WRITE(LUN,*) ' DIFDEC: LA, m_Beam, DELTAE, NCALL : ', 
     &           LA, AM(LA),DELTAE,'>',AM(7)+0.02D0,NCALL
            IF(Irec.ne.1) 
     &           WRITE(LUN,*) '   was recursive call! (ECM):',P0(5)
           ENDIF
           IF(NRJECT.eq.10) 
     &        write(lun,*)' this was the last warning.. good luck!'
         endif
         NRJECT = NRJECT + 1            
         NP = 0
         IBAD = 1
         RETURN
      ENDIF
      PAR(1) = PAR1_def
      END
C=======================================================================

      SUBROUTINE EXCT_RMNT(JW,KRMNT,IREJ)

C-----------------------------------------------------------------------
C     routine to produce massive excitations of beam and/or target \FR'14
C-----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)
      INTEGER NW_max
      PARAMETER (NW_max = 20)
      INTEGER NS_max, NH_max
      PARAMETER (NS_max = 20, NH_max = 80)
C--------------------------------------------------------------------
C     SIBYLL common blocks containing event information       \FR'14
C--------------------------------------------------------------------

C     EVENT INFO COMMON
C     contains overall interaction properties, like
C     SQS : center-of-mass energy
C     S   :         "       "     squared
C     PTmin : low pt cut of QCD cross section, 
C             i.e. minimal pt of hard minijets
C     Xmin : low-x bound for PDFs, 
C            i.e. minimal momentum fraction of hard partons
C     Zmin : logarithm of that
C     KB : PID of beam hadron
C     KT() : PID of target
C     IAT : mass number of target
      DOUBLE PRECISION SQS,S,PTMIN,XMIN,ZMIN
      INTEGER KB,IAT,KT
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT

      INTEGER IBMRDX,ITGRDX,IHMJDX,ISMJDX,ICSTDX,IINTDX
      COMMON /S_INDX/ IBMRDX(3),ITGRDX(NW_max,3),
     &     IHMJDX(NW_max*NH_max),IINTDX(NW_max),
     &     ISMJDX(NW_max*NS_max),ICSTDX(2*NW_max,3)

      INTEGER IRMNT,KRB,KRT
      DOUBLE PRECISION XRMASS,XRMEX
      COMMON /S_RMNT/ XRMASS(2),XRMEX(2),IRMNT(NW_max),KRB,KRT(NW_max)

      INTEGER ICHP,ISTR,IBAR
      COMMON /S_CHP/ ICHP(99), ISTR(99), IBAR(99)

      INTEGER IISO,ISPN
      COMMON /S_SPN/ IISO(99), ISPN(99)

      INTEGER ICHM
      COMMON /S_CHM/ ICHM(99)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
C     parameters that represent: NW: max. number of wounded nucleons,
C     NS,NH: max. number of soft and hard interactions
c      PARAMETER (NW_max = 20)
C     The COMMON block /S_CHIST/ contains information about the
C     the structure of the  generated event:
C     NWD   = number of wounded nucleons
C     NJET = total number of hard interactions
C     NSOF = total number of soft interactions
C     NNSOF (1:NW) = number of soft pomeron cuts in each interaction
C     NNJET (1:NW) = number of minijets produced in each interaction 
C     JDIF(1:NW) = diffraction code 
C                  0 : non-diff,
C                  1 : beam-diff
C                  2 : target-diff
C                  3 : double-diff
      INTEGER NNSOF,NNJET,JDIF,NWD,NJET,NSOF
      COMMON /S_CHIST/ NNSOF(NW_max),NNJET(NW_max),
     &     JDIF(NW_max),NWD,NJET,NSOF

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      INTEGER ITRY, NREJ
      COMMON /S_CNT/ ITRY(20), NREJ(20)
      DOUBLE PRECISION XM2MIN,ALXMIN,SLOP0,ASLOP,BSLOP,XMASS
      COMMON /S_DIFMAss/ XM2MIN(6),ALXMIN(6),SLOP0,ASLOP,BSLOP,XMASS(2)

      DOUBLE PRECISION P1(5),P2(5),P1N(5),P2N(5),PBM1(5),PBM2(5),PBM(5),
     &     PTG1(5),PTG2(5),PTG(5),PTT(5),GABE(4)
      DOUBLE PRECISION XMB,XMB2,ALPHA,XMSQMIN,XM2MAX,XM2,SHAT,ECM,EE,EE2
      DOUBLE PRECISION XMFRAC,XSFRAC,XMT,XMT2,XMT12,XMT22,P1TOT,P2TOT
      DOUBLE PRECISION DELTAE,XMMIN,COD,COF,SID,SIF,ANORF,PX,PY,PZ
      DOUBLE PRECISION XM1,ETOT,XI,XM2DIS,S_RNDM
c     DOUBLE PRECISION XDUMMY

      INTEGER   IMRG2HAD,LL(99)
      INTEGER   IBM1,IBM2,IBMST1,IBMST2,ITG1,ITG2,ITGST1,ITGST2,ITGH
      INTEGER   IDM,IFL,IBMH, IREF,I, II,K,J,JJ,L01,L02,NP0LD,NPLD
      INTEGER   JW,IREJ,KRMNT,LREJ,IBD,ICST11,ICST21
      INTEGER   IFLB1,IFLB2,IFLT1,IFLT2,L0,IDHAD,ISTH,IBMST,ITGST
      INTEGER   IFL1,IFL2,IMRG,IMST,IMST1,IMRGBAR,ICST2,LBD
      INTEGER   IMST11,IMST2,IMST21,ISTH1,ISTH2,IAFL1,IAFL2!,IMST22

      SAVE
      DATA LL /5*0,7*2,2*1,12*0,2,6*0,6*1,19*0,2,2,10*0,
     &     2,2,0,2,2,11*0,1,1,1,9*0,1/     

      
c     default return point, beam and target sampling
c      IREJ = 1

      IF(NDEBUG.gt.2)
     &     WRITE(LUN,*) ' EXCT_RMNT: input (JW,KRMNT,IREJ)', 
     &     JW,KRMNT,IREJ

      IF(NDEBUG.gt.3)THEN
         write(LUN,*)'  beam remnant index: (lvl0,flv1,flv2)  ',IBMRDX
         write(LUN,*)'  1st central string index: (lvl0,bm,tg)', 
     &        (ICSTDX(2*(JW-1)+1,ii),ii=1,3)
         write(LUN,*)'  2nd central string index: (lvl0,bm,tg)',
     &        (ICSTDX(2*(JW-1)+2,ii),ii=1,3)
         write(LUN,*)'  target remnant index: (lvl0,flv1,flv2)',
     &        (ITGRDX(JW,ii),ii=1,3)
      ENDIF

      ITRY(5) = 0

C...  select indices depending on configuration
C     krmnt = 0 : no excitation on either side
c           = 1 : beam side excited remnant
c           = 2 : target side
c           = 3 : both sides

c     write remnant configuration to remnant common
      IRMNT(JW) = KRMNT
      IF(KRMNT.eq.1)THEN
c     beam side remnant only
c     proto-remnant position: IBMRDX(1)
c     partons in : IBMRDX(2:3)
         IBM1 = IBMRDX(2)
         IBM2 = IBMRDX(3)
c     target side to transfer energy from: 
c     (sofar always choose valence pair)
         ITG1 = ICSTDX(2*(JW-1)+1,3)
         ITG2 = ICSTDX(2*(JW-1)+2,3)
c     beam-side partons to go into central strings
         IBMST1 = ICSTDX(2*(JW-1)+1,2)
         IBMST2 = ICSTDX(2*(JW-1)+2,2)
c     target-side partons to go into central strings
         ITGST1 = ITG1
         ITGST2 = ITG2

      ELSEIF(KRMNT.eq.2)THEN
c     target side remnant only
c     proto-remnant in ITGRDX(JW,1)
         ITG1 = ITGRDX(JW,2)
         ITG2 = ITGRDX(JW,3)
c     transfer energy from beam remnant or 
c     central strings with valence quarks
c     in ICSTDX(JW+0:1,2)
c     means no beam remnant --> get from valence strings
         IBM1 = ICSTDX(2*(JW-1)+1,2)
         IBM2 = ICSTDX(2*(JW-1)+2,2)
c     beam-side partons to go into central strings
         IBMST1 = IBM1
         IBMST2 = IBM2
c     target-side partons to go into central strings
         ITGST1 = ICSTDX(2*(JW-1)+1,3)
         ITGST2 = ICSTDX(2*(JW-1)+2,3)

      ELSEIF(KRMNT.eq.3)THEN
c     beam and target side remnant
c     transfer energy from pairs in rmnt or central strings
c     listed in I?RDX and ICSTDX()
         IBM1 = IBMRDX(2)
         IBM2 = IBMRDX(3)
         ITG1 = ITGRDX(JW,2)
         ITG2 = ITGRDX(JW,3)

      ELSEIF(KRMNT.eq.0)THEN
c     no excited remnant case, jump straight to central strings..
         GOTO 100

      ENDIF

      IF(NDEBUG.gt.3)then
         write(lun,*) '  beam parton1:  ',IBM1
         write(lun,*) '  beam parton2:  ',IBM2
         write(lun,*) '  target parton1:',ITG1
         write(lun,*) '  target parton2:',ITG2
      endif

c     save status of parton stack
      CALL GET_NPP(NPLD,NP0LD)

 10   ITRY(5) = ITRY(5) + 1
      IF(ITRY(5).GT.NREJ(5))THEN
         IF(NDEBUG.gt.2) 
     &        WRITE(LUN,*) ' EXCT_RMNT: no. of trials exceeded, ',
     &        NREJ(5), 'resample minijets ...' , IREJ
         RETURN
      ENDIF
c     reset parton stack after rmnt mass rejection
      CALL INI_PRTN_STCK(NPLD,NP0LD)

C..   construct 4momenta of proto-remnants
c     index of beam remnant on stack: IBMRDX(1)

C..   center-of-mass energy of parton system (s hat)
c     calculated in hadron-hadron frame
c     for first interaction (jw=1) partons are massless and collinear (sum pt=0)
c     in this case ecm = SQS*SQRT(XB*XT), xb,t=x1+x2 
c      for jw>1 beam partons may have already acquired mass and additional pt
c     therefore ecm = sqs*sqrt(xb*xt) + corr.
c     IRDX: index of remnant on parton stack
c      SHAT = S*XB*XT+XM2+(XT/XB)*XMT2

c     with 4momenta of partons on stack, momentum fractions are obsolete
c     center-of-mass energy is simply: shat = (pbm+ptg)**2

c     construct total 4momentum
c     add beam-side parton momenta, in had.-had. frame
      CALL RD_PRTN_4VEC(IBM1,PBM1,IFL,IDM)
      CALL RD_PRTN_4VEC(IBM2,PBM2,IFL,IBMH)
      CALL ADD_4VECS(PBM1,PBM2,PBM)
      
c     target-side parton momenta, in had.-had. frame
      CALL RD_PRTN_4VEC(ITG1,PTG1,IFL,IDM)
      CALL RD_PRTN_4VEC(ITG2,PTG2,IFL,IDM)
      CALL ADD_4VECS(PTG1,PTG2,PTG)
      
c     add beam and target side to get total 4momentum
      CALL ADD_4VECS(PBM,PTG,PTT)
      SHAT = PTT(5)**2
      ECM = PTT(5)
c     catch virtual remnants
      IF(PTT(5).LT.0.D0) THEN
         IF(NDEBUG.GT.2)THEN
            WRITE(LUN,*) ' EXCT_RMNT: too little mass left (Shat):',
     &           SHAT
            WRITE(LUN,*) '        resample minijets...'
         ENDIF
         LREJ = 2
         RETURN                 ! resample minijets
      ENDIF


      IF(NDEBUG.GT.2) WRITE(LUN,*) ' EXCT_RMNT: try no.',ITRY(5)
      IF(NDEBUG.GT.3)THEN
         write(LUN,*) '  4momenta before scattering:'
         write(LUN,*) '  PBM1:' , (PBM1(jj),jj=1,5)
         write(LUN,*) '  PBM2:' , (PBM2(jj),jj=1,5)
         write(LUN,*) '  PBM:' , (PBM(jj),jj=1,5)

         write(LUN,*) '  PTG1:' , (PTG1(jj),jj=1,5)
         write(LUN,*) '  PTG2:' , (PTG2(jj),jj=1,5)
         write(LUN,*) '  PTG:' , (PTG(jj),jj=1,5)

         write(LUN,*) '  PTT:' , (PTT(jj),jj=1,5)
      ENDIF

      IF(NDEBUG.gt.2)
     &     WRITE(LUN,*)' EXCT_RMNT: SHAT:',SHAT

      XMFRAC = PAR(81)
      XSFRAC = PAR(82)

c     exponent of remnant mass distribution (1/Mx**2)**alpha
c     by default: alpha = 1
c     different for baryons and mesons
c      ALPHA = PAR(98)

C..   Sample masses
      IF(KRMNT.eq.1)THEN
         XM2MAX = MIN(XSFRAC*S,XMFRAC*AM2(IABS(KB)))
         XM2MAX = MAX(XM2MAX,1.D0)

c     mass of target-side: 0
         XMT = 0.D0
         XMT2 = 0.D0
c     get remnant mass
c     (might have received mass from prior interaction)
         CALL GET_MASS2(IBMRDX(1),XM2)
c     allowing excitation to fallback to beam means min.
c     mass is beam mass, or more exact smallest mass of hadrons 
c     with flavors in remnant
         IF(IPAR(64).eq.1)THEN
c     remnant mass can also decrease through interactions
            XMSQMIN = AM2(IABS(KB))
         ELSE
c     remnant mass only increased by multiple interactions..
            XMSQMIN = MAX(AM2(IABS(KB)),XM2)
         ENDIF
C     select exponent from COMMON
         ALPHA = XRMEX(LL(IABS(KB)))
c     sample beam mass
         XMB2 = XM2DIS(XMSQMIN,XM2MAX,ALPHA)
         IF(NDEBUG.gt.2)
     &        WRITE(LUN,*)' EXCT_RMNT: XM2min,XM2max,ALPHA,XM2:',
     &        XMSQMIN,XM2MAX,ALPHA,XMB2
c     check if resonance or massive hadron has to be formed
         CALL SEL_RES(XMB2,KRB,IBMRDX(1),IBMH)
         XMB = dsqrt(XMB2)

      ELSEIF(KRMNT.eq.2)THEN
c     target side mass
         XM2MAX = MIN(XSFRAC*S,XMFRAC*AM2(IABS(KT(JW))))
         XM2MAX = MAX(XM2MAX,1.D0)

         XMB = 0.D0
         XMB2 = 0.D0
         XMSQMIN = AM2(KT(JW))
C     select exponent from COMMON
         ALPHA = XRMEX(LL(IABS(KT(JW))))
         XMT2 = XM2DIS(XMSQMIN,XM2MAX,ALPHA)
         IF(NDEBUG.gt.2)
     &        WRITE(LUN,*)' EXCT_RMNT: XM2min,XM2max,ALPHA,XM2:',
     &        XMSQMIN,XM2MAX,ALPHA,XMT2

c     check if resonance or massive hadron has to be formed
         CALL SEL_RES(XMT2,KRT(JW),ITGRDX(JW,1),ITGH)
         XMT = dsqrt(XMT2)

      ELSEIF(KRMNT.eq.3)THEN
         XM2MAX = MIN(XSFRAC*S,XMFRAC*AM2(IABS(KB)))
         XM2MAX = MAX(XM2MAX,1.D0)

         CALL GET_MASS2(IBMRDX(1),XM2)
         IF(IPAR(64).eq.1)THEN
c     remnant mass can also decrease through interactions
            XMSQMIN = AM2(IABS(KB))
         ELSE
c     remnant mass only increased by multiple interactions..
            XMSQMIN = MAX(AM2(IABS(KB)),XM2)
         ENDIF
C     select exponent from COMMON
         ALPHA = XRMEX(LL(IABS(KB)))        
         XMB2 = XM2DIS(XMSQMIN,XM2MAX,ALPHA)
         IF(NDEBUG.gt.2)
     &        WRITE(LUN,*)' EXCT_RMNT: XM2min,XM2max,ALPHA,xm2:',
     &        XMSQMIN,XM2MAX,ALPHA,XMB2

c     check if resonance or massive hadron has to be formed
         CALL SEL_RES(XMB2,KRB,IBMRDX(1),IBMH)
         XMB = SQRT(XMB2)
         
c     target always nucleon
         XM2MAX = MIN(XSFRAC*S,XMFRAC*AM2(IABS(KT(JW))))
         XM2MAX = MAX(XM2MAX,1.D0)

         XMSQMIN = AM2(IABS(KT(JW)))
C     select exponent from COMMON
         ALPHA = XRMEX(LL(IABS(KT(JW))))        
         XMT2 = XM2DIS(XMSQMIN,XM2MAX,ALPHA)
         IF(NDEBUG.gt.2)
     &        WRITE(LUN,*)' EXCT_RMNT: XM2min,XM2max,alpha,XM2:',
     &        XMSQMIN,XM2MAX,ALPHA,XMT2

c     check if resonance or massive hadron has to be formed
         CALL SEL_RES(XMT2,KRT(JW),ITGRDX(JW,1),ITGH)
         XMT = SQRT(XMT2)

      ENDIF
c     write excitation mass to output common
      XRMASS(1) = XMB
      XRMASS(2) = XMT

c     minimal mass requirement
c      IF(SHAT.lt.XMB2+XMT2+0.3) GOTO 10
      IF(SHAT.lt.XMB2+XMT2+2.D0*XMB*XMT+0.3D0) GOTO 10

C     transfer cm energy to mass of particle in parton-parton cm
      CALL TRANSFONSHELL(ECM,XMB,XMT,XM2MAX,1,P1,P2,IBD)
      IF(IBD.eq.1) THEN
         IF(NDEBUG.gt.2) WRITE(LUN,*) ' EXCT_RMNT: excitation rejected!'
         RETURN
      ENDIF

C...  Boost 4momenta to hadron-hadron center-of-mass
c     along z only if initial partons do not carry transverse momentum
c     (cancels between val1 and val2)
c     with multiple nucleons interacting beam val partons can aquire 
c     transverse momentum from the target. in this case need arbitrary boost
      DO K = 1,4
         GABE(k) = PTT(k)/PTT(5)
      ENDDO
      CALL SIB_ALTRA(GABE(4), GABE(1), GABE(2), GABE(3),
     &     P1(1),P1(2),P1(3),P1(4),
     &     P1TOT,P1N(1),P1N(2),P1N(3),P1N(4))
      P1N(5)=P1(5)
      CALL SIB_ALTRA(GABE(4), GABE(1), GABE(2), GABE(3),
     &     P2(1),P2(2),P2(3),P2(4),
     &     P2TOT,P2N(1),P2N(2),P2N(3),P2N(4))
      P2N(5)=P2(5)

C...  Calculate new 4momentum of partons in had.-had. frame
c     P1,P2: momenta after scattering in parton-parton cm.
c     P1n,P2n: momenta after scattering in had.-had. cm
c     PBM1,2: momenta of beam partons in had.-had. before scattering
c     PTG1,2: momenta of target partons in had.-had. before scattering
c     PBM: combined momentum of all beam partons before scattering
c     PTG: combined momentum of all target partons before scattering

c     energy and z component
      DO II=3,4
         PBM1(II) = PBM1(II)*P1n(II)/PBM(II)
         PBM2(II) = PBM2(II)*P1n(II)/PBM(II)

         PTG1(II) = PTG1(II)*abs(P2N(II)/PTG(II))
         PTG2(ii) = PTG2(ii)*abs(P2N(II)/PTG(II))
      ENDDO

c     if transverse momentum prior to interaction zero then
c     assign transverse momentum of partons according to random fraction
      IF(ABS(PBM(1)).LT.EPS10.or.ABS(PBM(2)).LT.EPS10)THEN
         DO II = 1,2
            XI = S_RNDM(II)
            PBM1(II) = XI*P1N(Ii)
            PBM2(II) = (1.D0-XI)*P1N(II)
         ENDDO
      ELSE
         DO II=1,2
            PBM1(II) = PBM1(II)*P1N(II)/PBM(II)
            PBM2(II) = PBM2(II)*P1N(II)/PBM(II)
         ENDDO         
      ENDIF

      IF(ABS(PTG(1)).LT.EPS10.or.ABS(PTG(2)).LT.EPS10)THEN
         DO II=1,2
            XI = S_RNDM(II)
            PTG1(II) = XI*P2N(II)
            PTG2(II) = (1.D0-XI)*P2N(II)
         ENDDO
      ELSE
         DO II=1,2
            PTG1(II) = PTG1(II)*P2N(II)/PTG(II)
            PTG2(II) = PTG2(II)*P2N(II)/PTG(II)
         ENDDO                  
      ENDIF

      IF(NDEBUG.GT.3)THEN
         write(LUN,*) '  parton 4momenta after scattering:'
         write(LUN,*) '   PBM1:' , (PBM1(jj),jj=1,5)
         write(LUN,*) '   PBM2:' , (PBM2(jj),jj=1,5)
         write(LUN,*) '   sum: ' , (PBM2(jj)+PBM1(jj),jj=1,5)
         write(LUN,*) '   PTG1:' , (PTG1(jj),jj=1,5)
         write(LUN,*) '   PTG2:' , (PTG2(jj),jj=1,5)
         write(LUN,*) '   sum: ' , (PTG2(jj)+PTG1(jj),jj=1,5)
      ENDIF
      
C...  change parton 4momenta on stack
      CALL EDT_PRTN(IBM1,PBM1(1),PBM1(2),PBM1(3),PBM1(4),PBM1(5),IDM)
      CALL EDT_PRTN(IBM2,PBM2(1),PBM2(2),PBM2(3),PBM2(4),PBM2(5),IDM)

      CALL EDT_PRTN(ITG1,PTG1(1),PTG1(2),PTG1(3),PTG1(4),PTG1(5),IDM)
      CALL EDT_PRTN(ITG2,PTG2(1),PTG2(2),PTG2(3),PTG2(4),PTG2(5),IDM)
         
C...  add remnants
c     references are circular: 
c     rmnt --> parton1 --> parton2 --> lvl2 rmnt (hadron) --> rmnt
      IF(KRMNT.eq.1)THEN
c     beam side remnant, add only if does not exist yet otherwise edit
         IF(IBMRDX(1).eq.0)THEN
            CALL ADD_PRTN
     &           (P1N(1),P1N(2),P1N(3),P1N(4),P1N(5),2,0,IBM1,IBMRDX(1))
         ELSE
            CALL EDT_PRTN
     &           (IBMRDX(1),P1N(1),P1N(2),P1N(3),P1N(4),P1N(5),IREF)
         ENDIF
c     add beam hadron as hypothetical final state
         IF(IBMH.eq.0)THEN
            CALL ADD_PRTN
     &         (P1N(1),P1N(2),P1N(3),P1N(4),P1N(5),KRB,2,IBMRDX(1),IBMH)
         ELSE
            CALL EDT_PRTN
     &           (IBMH,P1N(1),P1N(2),P1N(3),P1N(4),P1N(5),IREF)
         ENDIF
c     add references rmnt --> parton1 etc
         CALL ADD_REF(IBMRDX(1),IBM1)
         CALL ADD_REF(IBM1,IBM2)
         CALL ADD_REF(IBM2,IBMH)

      ELSEIF(KRMNT.eq.2)THEN
c     add target side remnant
         IF(ITGRDX(JW,1).eq.0)THEN
            CALL ADD_PRTN
     &           (P2N(1),P2N(2),P2N(3),P2N(4),P2N(5),
     &           -2,0,0,ITGRDX(JW,1))
         ELSE
            CALL EDT_PRTN
     &           (ITGRDX(JW,1),P2N(1),P2N(2),P2N(3),P2N(4),P2N(5),IREF)
         ENDIF
         IF(ITGH.eq.0)THEN
c     add target hadron as hypothetical final state, always nucleon
            CALL ADD_PRTN
     &           (P2N(1),P2N(2),P2N(3),P2N(4),P2N(5),
     &           KRT(JW),2,ITGRDX(JW,1),ITGH)
         ELSE
            CALL EDT_PRTN
     &           (ITGH,P2N(1),P2N(2),P2N(3),P2N(4),P2N(5),IREF)
         ENDIF

c     add references rmnt --> parton1 etc
         CALL ADD_REF(ITGRDX(JW,1),ITG1)
         CALL ADD_REF(ITG1,ITG2)
         CALL ADD_REF(ITG2,ITGH)

      ELSEIF(KRMNT.eq.3)THEN
c     beam side remnant, add only if does not exist yet, otherwise edit
         IF(IBMRDX(1).EQ.0)THEN
            CALL ADD_PRTN
     &           (P1N(1),P1N(2),P1N(3),P1N(4),P1N(5),2,0,0,IBMRDX(1))
         ELSE
            CALL EDT_PRTN
     &           (IBMRDX(1),P1N(1),P1N(2),P1N(3),P1N(4),P1N(5),IREF)
         ENDIF
c     add beam hadron as hypothetical final state
         IF(IBMH.EQ.0)THEN
            CALL ADD_PRTN
     &         (P1N(1),P1N(2),P1N(3),P1N(4),P1N(5),KRB,2,IBMRDX(1),IBMH)
         ELSE
            CALL EDT_PRTN
     &           (IBMH,P1N(1),P1N(2),P1N(3),P1N(4),P1N(5),IREF)
         ENDIF
         CALL ADD_REF(IBMRDX(1),IBM1)
         CALL ADD_REF(IBM1,IBM2)
         CALL ADD_REF(IBM2,IBMH)

c     add target side remnant
         IF(ITGRDX(JW,1).eq.0)THEN
            CALL ADD_PRTN
     &           (P2N(1),P2N(2),P2N(3),P2N(4),P2N(5),-2,0,0,IREF)
            ITGRDX(JW,1) = IREF
         ELSE
            CALL EDT_PRTN
     &           (ITGRDX(JW,1),P2N(1),P2N(2),P2N(3),P2N(4),P2N(5),IREF)
         ENDIF
         IF(ITGH.eq.0)THEN
c     add target hadron as hypothetical final state
            CALL ADD_PRTN
     &           (P2N(1),P2N(2),P2N(3),P2N(4),P2N(5),
     &           KRT(JW),2,ITGRDX(JW,1),ITGH)
         ELSE
            CALL EDT_PRTN
     &           (ITGH,P2N(1),P2N(2),P2N(3),P2N(4),P2N(5),IREF)
         ENDIF
c     add references rmnt --> parton1 etc
         CALL ADD_REF(ITGRDX(JW,1),ITG1)
         CALL ADD_REF(ITG1,ITG2)
         CALL ADD_REF(ITG2,ITGH)

      ENDIF

 100  IF(JDIF(JW).ne.0.and.NWD.ne.1)THEN
c     incoherent diffraction case
c     add parton 4momenta to obtain c.m energy
         
c     beam side
         IBMST1 = ICSTDX(2*(JW-1)+1,2)
         IBMST2 = ICSTDX(2*(JW-1)+2,2)

c     target side
         ITGST1 = ICSTDX(2*(JW-1)+1,3)
         ITGST2 = ICSTDX(2*(JW-1)+2,3)
         
         CALL RD_PRTN_4VEC(IBMST1,PBM1,IFLB1,IDM)
         CALL RD_PRTN_4VEC(IBMST2,PBM2,IFLB2,IDM)
         CALL ADD_4VECS(PBM1,PBM2,PBM)
         CALL RD_PRTN_4VEC(ITGST1,PTG1,IFLT1,IDM)
         CALL RD_PRTN_4VEC(ITGST2,PTG2,IFLT2,IDM)
         CALL ADD_4VECS(PTG1,PTG2,PTG)
c     total 4momentum
         CALL ADD_4VECS(PBM,PTG,PTT)
c     add diffractive system to parton stack
c     references are: diff --> diff. hadron 
c     --> beam parton1 --> beam parton2 --> target parton1 etc
         CALL ADD_PRTN_4VEC(PTT,-10*JDIF(JW),0,IBMST1,IREF)
         CALL ADD_INT_REF(IREF,IINTDX(JW))
c     both string indices point to diff. system
         ICSTDX(2*(JW-1)+1,1) = IREF
         ICSTDX(2*(JW-1)+2,1) = IREF
c     add diff. beam hadron to stack
c     model assumes remnant always excited in first interaction
         L0 = KB
c     if not first interaction or remnant excited, merge sea pair to hadron
         IF(KRMNT.ne.0.or.JW.ne.1) THEN       
            L0 = IMRG2HAD(IFLB1,IFLB2)
c     CALL SIB_I4FLAV(IFLB1,IFLB2,IDM,IDM1,L0)
         ENDIF
c     check kinematic limits
c     m2_max should be smaller than m2_min
         IREJ = 1
         EE = PTT(5)
         EE2 = PTT(5)**2
         K = 2-IBAR(IABS(L0))
         IF(JDIF(jw).gt.1)THEN
            DELTAE = EE-AM(13)
            XMMIN=max(XM2MIN(1),(AM(IABS(l0))+AM(7)+0.02D0)**2)
         ELSE
            DELTAE = EE-AM(IABS(L0))
            XMMIN=max(XM2MIN(K),(AM(IABS(l0))+AM(7)+0.02D0)**2)
         ENDIF
c         print *,'jw,jdif,nwd,l0,ifl1,ifl2,deltae,xmin,ee,xmax',
c     &        jw,jdif(jw),nwd,l0,ifl1,ifl2,deltae,xmmin,ee,par(13)*ee2
         IF(DELTAE.lt.AM(7)+0.02D0) THEN
            IF(ndebug.gt.2) 
     &           WRITE(lun,*) ' EXCT_RMNT: inchoherent diff. :',
     &           ' not enough mass left for excitation! (DELTAE,PION,',
     &           'IREJ,NCALL)',DELTAE,AM(7)+0.02D0,IREJ,NCALL
            RETURN
         ENDIF
         IF(PAR(13)*EE2.lt.XMMIN)THEN
            IF(ndebug.gt.2)
     &           WRITE(lun,*) ' EXCT_RMNT: inchoherent diff. :',
     &           ' not enough mass left for excitation! (min,max,',
     &           'IREJ,NCALL)',PAR(13)*EE2,XMMIN,IREJ,NCALL
            RETURN
         ENDIF
         CALL ADD_PRTN_4VEC(PTT,L0,2,IBMST1,IDHAD)
         CALL ADD_REF(IREF,IDHAD)
c     reset references of partons
         CALL ADD_REF(IBMST1,IBMST2)
         CALL ADD_REF(IBMST2,ITGST1)
         CALL ADD_REF(ITGST1,ITGST2)
         CALL ADD_REF(ITGST2,IREF)
         IF(ndebug.gt.2) THEN
            WRITE(LUN,*) ' EXCT_RMNT: incoherent diff. ',
     &           '(IDX,IDX2,JDIF,ECM,L0)',IREF,IDHAD,JDIF(JW),PTT(5),L0
            WRITE(LUN,*) ' EXCT_RMNT: DELTAE,XM2MAX:',DELTAE,PAR(13)*EE2
         ENDIF
         IREJ = 0
         RETURN
      ENDIF

C...  add central strings to stack
c     partons designated for central strings 
c     are indexed in ICSTDX(JW,2:3)
c     pstr_j = p_j_bm + p_j_tg
c     string mass ** 2 = pstr_j ** 2
c     --> read momenta from stack, add beam and target side, 
c     references are set in a loop:
c     string --> beam-parton --> target-parton --> string
c     then write string 4momentum on stack
      IMRG = 0
      DO JJ=1,2
         ISTH = 0
         IBMST = ICSTDX(2*(JW-1)+JJ,2)
         ITGST = ICSTDX(2*(JW-1)+JJ,3)
         CALL RD_PRTN_4VEC(IBMST,PBM1,IFL1,IDM)
         CALL RD_PRTN_4VEC(ITGST,PTG1,IFL2,IDM)
         CALL ADD_4VECS(PBM1,PTG1,PTT)
c     transverse mass of string end partons (pt**2)
         CALL GET_XMT2(IBMST,XMT12)
         CALL GET_XMT2(ITGST,XMT22)
c     available mass for string
         EE = SQRT(PTT(4)**2-PTT(3)**2)
c     catch virtual strings
         IF(PTT(5).lt.0.D0) THEN
            IREJ = 1
            IF(ndebug.gt.2)
     &           write(LUN,*)' EXCT_RMNT: virt. string (M):',EE
            IF(ndebug.gt.3)then
               CALL GET_IMASS2(IBMST,XM2)
               write(LUN,*) '  PBM1:', (PBM1(j),j=1,5),XM2
               CALL GET_IMASS2(ITGST,xm2)
               write(LUN,*) '  PTG1:', (PTG1(j),j=1,5),XM2
               write(LUN,*) '  Ptot:', (PTT(j),j=1,5)
            ENDIF
c               stop
            RETURN
         ENDIF
c     minimal string mass requirement
         IF(EE.lt.sqrt(XMT12)+sqrt(XMT22)+PAR(123))THEN
            IAFL1 = IABS(IFL1)
            IAFL2 = IABS(IFL2)
            IF(IPAR(74).eq.1)THEN
c     try to form single meson, set merge flag
               IF(IAFL1.gt.10.and.IAFL2.gt.10) THEN
c     skip if two diquarks need merging..                  
                  IREJ = 1
                  RETURN
               ENDIF
               IF((IAFL1/10.eq.4.or.mod(IAFL1,10).eq.4)
     +              .and.(IAFL2/10.eq.4.or.mod(IAFL2,10).eq.4)) THEN
c     skip if two charm quarks need merging..                  
                  IREJ = 1
                  RETURN
               ENDIF                             
               L0 = IMRG2HAD(IFL1,IFL2)
               IF(EE.gt.AM(IABS(L0))) then
                  IMRG = IMRG + JJ
                  CALL ADD_PRTN_4VEC(PTT,L0,2,IBMST,ISTH)
                  IF(ndebug.gt.2)then
                     write(lun,*)
     &                    ' EXCT_RMNT: c.string mass too low! ',
     &                    'merge into hadron..',l0
                  ENDIF
               ENDIF
            ELSE
               IF(ndebug.gt.2)then
                  write(lun,*)
     &                 ' EXCT_RMNT: c.string kinematic rejection!'
                  write(lun,*) ' EE,limit,XMT1,XMT2:',
     &                 EE,sqrt(XMT12)+sqrt(XMT22)+0.3D0,sqrt(XMT12),
     &                 sqrt(XMT22)
                  write(lun,*) ' return to momentum sampling..'
               endif
               IREJ = 1
               RETURN
            ENDIF
         ENDIF
c     add central string to stack, refering to beam-end parton
         CALL ADD_PRTN_4VEC(PTT,1,0,IBMST,IREF)
         ICSTDX(2*(JW-1)+JJ,1) = IREF
         CALL ADD_INT_REF(Iref,IINTDX(JW))
c     add reference to target parton to beam parton
         CALL ADD_REF(IBMST,ITGST)
         IF(ISTH.ne.0) THEN
c     if string merged to hadron add reference corresponding reference            
            CALL ADD_REF(ITGST,ISTH)
            CALL ADD_REF(ISTH,IREF)
         ELSE
c     add reference to corresponding central string to target parton
            CALL ADD_REF(ITGST,IREF)
         ENDIF
      ENDDO
      
c     form single hadron from string if mass was too low ..
c     need to put hadron on shell by exchanging energy with other string            
      IF(IMRG.eq.1.or.IMRG.eq.2)THEN
         IF(ndebug.gt.2)
     &        WRITE(LUN,*)' EXCT_RMNT: merging one string..',IMRG
c     one string merged
c     index of merged string and its last parton
         IMST = ICSTDX(2*(JW-1)+IMRG,1)
         IMST1 = ICSTDX(2*(JW-1)+IMRG,3)
c     index of ordinary string
         IMRGBAR = 3-IMRG
         ICST2 = ICSTDX(2*(JW-1)+IMRGBAR,1)
c     read 4momenta
         CALL RD_REF(IMST1,ISTH)
         CALL RD_PRTN_4VEC(ISTH,P1,L0,IREF)
c     string two
         CALL RD_PRTN_4VEC(ICST2,P2,IFL2,IDM)
c     cm energy
         CALL ADD_4VECS(P1,P2,PTT)
         IF(ndebug.gt.2)THEN
            write(lun,*)' EXCT_RMNT: string A :',(P1(i),i=1,5)
            write(lun,*)' EXCT_RMNT: string B :',(P2(i),i=1,5)
            write(lun,*)' EXCT_RMNT: total :',(PTT(i),i=1,5)
         ENDIF
         ECM = PTT(5)
         XM1 = AM(IABS(L0))
         XM2 = P2(5)
         CALL TRANSFONSHELL(ECM,XM1,XM2,1.D0,3,P1N,P2N,LBD)
         IF(LBD.eq.1) THEN
            IF(NDEBUG.gt.2)
     &           WRITE(LUN,*)' EXCT_RMNT: mass transfer failed!'
            RETURN
         ENDIF
c     by definition p1n is along +z in string cm, need to invert if pzA < pzB
c         IF(P2(3).gt.P1(3)) CALL SWTCH_LMNTS(P1N(3),P2N(3))

C..   rotate parton-parton axis onto string-string axis
c     therefore boost to parton-parton cm
c     to calc. rotation angles BEFORE interaction !
         DO K = 1,4
            GABE(K) = PTT(K)/PTT(5)
         enddo         
         CALL SIB_ALTRA(GABE(4),-GABE(1),-GABE(2),-GABE(3),
     &        P1(1),P1(2),P1(3),P1(4),
     &        P1TOT,PBM1(1),PBM1(2),PBM1(3),PBM1(4))
c     rotation factors
         COD= PBM1(3)/P1TOT
         SID= DSQRT(PBM1(1)**2+PBM1(2)**2)/P1TOT
         COF=1.D0
         SIF=0.D0
         IF(P1TOT*SID.GT.EPS5) THEN
            COF=PBM1(1)/(SID*P1TOT)
            SIF=PBM1(2)/(SID*P1TOT)
            ANORF=DSQRT(COF*COF+SIF*SIF)
            COF=COF/ANORF
            SIF=SIF/ANORF
         ENDIF
         IF(ndebug.gt.2)THEN
            write(lun,*)' EXCT_RMNT: momentum in cm:',(PBM1(i),i=1,5)
            write(lun,*)' EXCT_RMNT: rotation factors:',COD,SID,COF,SIF
            write(lun,*)' EXCT_RMNT: rotation angles (theta,phi):',
     &           ACOS(COD),ACOS(COF),ASIN(SID),ASIN(SIF)
            write(lun,*)' EXCT_RMNT: momentum:',
     &           sqrt(P1N(1)**2+P1N(2)**2+P1N(3)**2)
         ENDIF
c     rotate parton momenta after interaction, still in parton-parton frame
         CALL SIB_TRANI(P1N(1),P1N(2),P1N(3),COD,SID,COF,SIF
     &        ,PX,PY,PZ)
         P1N(1)=PX
         P1N(2)=PY
         P1N(3)=PZ
         CALL SIB_TRANI(P2N(1),P2N(2),P2N(3),COD,SID,COF,SIF
     &        ,PX,PY,PZ)
         P2N(1)=PX
         P2N(2)=PY
         P2N(3)=PZ
         IF(ndebug.gt.2) write(lun,*)' EXCT_RMNT: momentum*:',
     &        sqrt(P1N(1)**2+P1N(2)**2+P1N(3)**2)

c     boost back to hadron-hadron
         DO K = 1,4
            GABE(K) = PTT(K)/PTT(5)
         ENDDO
         CALL SIB_ALTRA(GABE(4), GABE(1), GABE(2), GABE(3),
     &        P1N(1),P1N(2),P1N(3),P1N(4),
     &        P1TOT,P1(1),P1(2),P1(3),P1(4))
         P1(5)=P1N(5)
         CALL SIB_ALTRA(GABE(4), GABE(1), GABE(2), GABE(3),
     &        P2N(1),P2N(2),P2N(3),P2N(4),
     &        P2TOT,P2(1),P2(2),P2(3),P2(4))
         p2(5)=p2n(5)
         IF(ndebug.gt.2)THEN
            write(lun,*)' EXCT_RMNT: momenta after scattering:'
            write(lun,*)' EXCT_RMNT: hadron A :',(P1(i),i=1,5)
            write(lun,*)' EXCT_RMNT: string B :',(P2(i),i=1,5)
         ENDIF

c     edit partons on stack
         CALL EDT_PRTN
     &        (ISTH,P1(1),P1(2),P1(3),P1(4),P1(5),IREF)
         ICST11 = ICSTDX(2*(JW-1)+IMRG,2)
         CALL EDT_PRTN
     &        (IMST,P1(1),P1(2),P1(3),P1(4),P1(5),ICST11)
         ICST21 = ICSTDX(2*(JW-1)+IMRGBAR,2)
         CALL EDT_PRTN
     &        (ICST2,P2(1),P2(2),P2(3),P2(4),P2(5),ICST21)
         
      ELSEIF(IMRG.eq.3)THEN
         IF(ndebug.gt.2)
     &        WRITE(LUN,*)' EXCT_RMNT: merge both strings..'

c     both strings merged
c     index of merged string and its last parton
         IMST1 = ICSTDX(2*(JW-1)+1,1)
         IMST11 = ICSTDX(2*(JW-1)+1,3)
c     index of ordinary string
         IMST2 = ICSTDX(2*(JW-1)+2,1)
         IMST21 = ICSTDX(2*(JW-1)+2,3)
c     read 4momenta
         CALL RD_REF(IMST11,ISTH1)
         CALL RD_PRTN_4VEC(ISTH1,P1,L01,IREF)
c     string two
         CALL RD_REF(IMST21,ISTH2)
         CALL RD_PRTN_4VEC(ISTH2,P2,L02,IREF)
         XM1 = AM(IABS(L01))
         XM2 = AM(IABS(L02))
c     cm energy
         CALL ADD_4VECS(P1,P2,PTT)
         ECM = PTT(5)
         ETOT = PTT(4)
         IF(ndebug.gt.2)THEN
            write(lun,*)' EXCT_RMNT: string A :',(P1(i),i=1,5)
            write(lun,*)' EXCT_RMNT: string B :',(P2(i),i=1,5)
            write(lun,*)' EXCT_RMNT: total :',(PTT(i),i=1,5)
         ENDIF

         CALL TRANSFONSHELL(ecm,xm1,xm2,1.D0,3,P1n,P2n,LBD)
         IF(LBD.eq.1) THEN
            IF(NDEBUG.gt.2)
     &           WRITE(LUN,*)' EXCT_RMNT: mass transfer failed!'
            RETURN
         ENDIF
c     by definition p1n is along +z in string cm, need to invert if pzA < pzB
c         IF(P2(3).gt.P1(3)) CALL SWTCH_LMNTS(P1N(3),P2N(3))
c     rotate parton-parton axis onto string-string axis
c     boost to parton-parton cm to calc. rotation angles BEFORE interaction!
         DO K = 1,4
            GABE(K) = PTT(K)/PTT(5)
         ENDDO
         CALL SIB_ALTRA(GABE(4),-GABE(1),-GABE(2),-GABE(3),
     &        P1(1),P1(2),P1(3),P1(4),
     &        P1TOT,PBM1(1),PBM1(2),PBM1(3),PBM1(4))
c     rotation factors
         COD= PBM1(3)/P1TOT
         SID= DSQRT(PBM1(1)**2+PBM1(2)**2)/P1TOT
         COF=1.D0
         SIF=0.D0
         IF(P1TOT*SID.GT.EPS5) THEN
            COF=PBM1(1)/(SID*P1TOT)
            SIF=PBM1(2)/(SID*P1TOT)
            ANORF=DSQRT(COF*COF+SIF*SIF)
            COF=COF/ANORF
            SIF=SIF/ANORF
         ENDIF
c     rotate parton momenta after interaction
         CALL SIB_TRANI(P1N(1),P1N(2),P1N(3),COD,SID,COF,SIF
     &        ,PX,PY,PZ)
         P1N(1)=PX
         P1N(2)=PY
         P1N(3)=PZ
         CALL SIB_TRANI(P2N(1),P2N(2),P2N(3),COD,SID,COF,SIF
     &        ,PX,PY,PZ)
         P2N(1)=PX
         P2N(2)=PY
         P2N(3)=PZ

c     boost massive hadrons back to hadron-hadron
         CALL SIB_ALTRA(GABE(4), GABE(1), GABE(2), GABE(3),
     &        P1N(1),P1N(2),P1N(3),P1N(4),
     &        P1TOT,P1(1),P1(2),P1(3),P1(4))
         P1(5)=P1N(5)
         CALL SIB_ALTRA(GABE(4), GABE(1), GABE(2), GABE(3),
     &        P2N(1),P2N(2),P2N(3),P2N(4),
     &        P2TOT,P2(1),P2(2),P2(3),P2(4))
         P2(5)=P2N(5)
         IF(ndebug.gt.2)THEN
            write(lun,*)' EXCT_RMNT: hadron A :',(P1(i),i=1,5)
            write(lun,*)' EXCT_RMNT: hadron B :',(P2(i),i=1,5)
         ENDIF

c     edit partons on stack
         CALL EDT_PRTN
     &        (ISTH1,P1(1),P1(2),P1(3),P1(4),P1(5),IREF)
         ICST11 = ICSTDX(2*(JW-1)+1,2)
         CALL EDT_PRTN
     &        (IMST1,P1(1),P1(2),P1(3),P1(4),P1(5),ICST11)

         CALL EDT_PRTN
     &        (ISTH2,P2(1),P2(2),P2(3),P2(4),P2(5),IREF)
         ICST21 = ICSTDX(2*(JW-1)+2,2)
         CALL EDT_PRTN
     &        (IMST2,P2(1),P2(2),P2(3),P2(4),P2(5),ICST21)
         
      ENDIF

      IREJ = 0

      RETURN
      END
C=======================================================================

      SUBROUTINE FIREBALL_4FLV(L0,P0,PCHEXin,IREJ)

C-----------------------------------------------------------------------
C... "decay" of an excited state with the quantum numbers
C.   of particle L0 and the 5-momentum P0
C.   4 flavor generalization /FR'13
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
C     The final particle output is contained in COMMON /S_PLIST/    
C     NP           : number of final particles
C     P(1:NP, 1:5) : 4-momenta + masses of the final particles 
C     LLIST (1:NP) : codes of final particles
      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP
      INTEGER NW_max
      PARAMETER (NW_max = 20)
C--------------------------------------------------------------------
C     SIBYLL common blocks containing event information       \FR'14
C--------------------------------------------------------------------

C     EVENT INFO COMMON
C     contains overall interaction properties, like
C     SQS : center-of-mass energy
C     S   :         "       "     squared
C     PTmin : low pt cut of QCD cross section, 
C             i.e. minimal pt of hard minijets
C     Xmin : low-x bound for PDFs, 
C            i.e. minimal momentum fraction of hard partons
C     Zmin : logarithm of that
C     KB : PID of beam hadron
C     KT() : PID of target
C     IAT : mass number of target
      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
      INTEGER KB,IAT,KT
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)

      INTEGER ICHP,ISTR,IBAR
      COMMON /S_CHP/ ICHP(99), ISTR(99), IBAR(99)

      INTEGER IISO,ISPN
      COMMON /S_SPN/ IISO(99), ISPN(99)

      INTEGER ICHM
      COMMON /S_CHM/ ICHM(99)

      CHARACTER*6 NAMP
      COMMON /S_CNAM/ NAMP (0:99)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)

      INTEGER NFORIG,NPORIG,NIORIG,IPFLAG,IIFLAG,KINT
      COMMON /S_PARTO/ NFORIG(NP_max),NPORIG(NP_max),NIORIG(NP_max),
     &IPFLAG,IIFLAG,KINT

      DIMENSION P0(5), LL(10), PD(10,5), IFL(3), INONLEAD(2)
      DIMENSION LRESCHEX(6:99), LRES(6:99), LCON(6:99), LPIC(-1:1)
      DIMENSION LSTR(6:99), LPICS(-2:2)
      
C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      SAVE
c     charge exchange map
      DATA (LCON(I),I=6, 33)  /7,6,6,21,22,9,9,14,13,4*0,20,19,9,10,23,
     &     24,27,27,25,30,31,28,29,32,33/
      DATA (LCON(I),I=34, 49) 
     &     /35,34,35,38,37,39,41,42,41,42,45,44,45,48,47,49/
      DATA (LCON(I),I=50, 83) /0,52,51,54,53,4*0,71,72,10*0,
     &     59,60,73,74,75,76,77,80,81,78,79,0,83/
      DATA (LCON(I),I=84, 99) /84,85,86,87,88,89,4*0,94,95,96,97,98,99/
c     pion charge conversion map
      DATA LPIC /8,6,7/
c     kaon charge conversion map      
      DATA LPICS /9,21,0,22,10/     
c     charge exchange to resonances map
      DATA (LRESCHEX(I),I=6, 33) /26,27,27,30,31,9,9,42,41,19*0/
      DATA (LRESCHEX(I),I=34, 39) /45,44,45,48,47,39/ 
      DATA (LRESCHEX(I),I=40, 49) /41,42,43,42,45,46,45,48,47,49/
      DATA (LRESCHEX(I),I=50, 83) 
     &     /0,52,51,54,53,4*0,60,59,10*0,71,72,73,75,74,
     &     77,76,79,78,80,81,0,83/
      DATA (LRESCHEX(I),I=84, 99) 
     &     /84,85,86,87,88,89,4*0,94,95,96,97,98,99/
c     resonance excitation map
      DATA (LRES(I),I=6, 39) 
     &     /27,25,26,28,29,9,9,41,42,19*0,44,45,46,47,48,39/
      DATA (LRES(I),I=40, 49) /40,41,42,43,44,45,46,47,48,49/
      DATA (LRES(I),I=50, 83) 
     &     /0,51,52,53,54,4*0,78,79,10*0,71,72,73,76,77,76,
     &     77,78,79,80,81,0,83/
      DATA (LRES(I),I=84, 99) /94,95,96,97,98,89,4*0,94,95,96,97,98,99/
c     strangeness excitation map
      DATA LSTR(6:27) /6,7,8,9,10,11,12,34,39,6*0,21,22,23,24,25,26,27/
      DATA LSTR(28:39) /28,29,30,31,32,33,44,45,46,47,48,39/
      DATA LSTR(40:49) /40,41,42,43,44,45,46,47,48,49/
      DATA LSTR(50:83) /0,51,52,53,54,4*0,78,79,10*0,71,72,73,76,77,76,
     &     77,78,79,80,81,0,83/
      DATA LSTR(84:99) /94,95,96,97,98,89,4*0,94,95,96,97,98,99/
      
c...  charge exchange reaction rate
c      DATA PCHEX /0.33/

c     default parameter: PAR(61)
      PCHEX = PCHEXin

c     split charge exchange between 2 and 3+ fireballs
      IF(IPAR(91).eq.1.and.NPI.gt.2)THEN
         PCHEX = 1.D0-PCHEX
      ENDIF

c     hyperon production rate
      PLAM = PAR(157)
      
c...  suppression of high mass particles in fireball
c     xmpsuppr = prob. accepting additional proton
      XMPSUPPR=PAR(33)
      IF(ABS(XMPSUPPR).lt.EPS3) THEN
         WRITE(LUN,*)
     &        ' Error: too low mass suppression in 4 flv fireball!'
         WRITE(LUN,*)
     &        ' Probably PAR(33)/IPAR(14) not properly set, aborting..'
         STOP
      ENDIF
      XTEMPH=(AM(6)-AM(13))/dLOG(XMPSUPPR)

      IF(Ndebug.gt.3) THEN
         WRITE(LUN,*)' FIRBALL_4FLV: called with (L0,P0):',
     &        L0,P0
         WRITE(LUN,*)' 2nd Proton rejection prob.:',XMPSUPPR
         WRITE(LUN,*)' fireball temperature:',XTEMPH
         WRITE(LUN,*)' charge exchange prob.:',PCHEX
         WRITE(LUN,*)' multiplicity width:',PAR(38)
      ENDIF

c...  special vector resonance treatment for meson projectiles
c     i.e. spin exchange probability
      PAR5def = PAR(5)
      IF(IPAR(14).eq.-2.and.abs(kb).lt.13)THEN
         PAR(5)=PAR(34)
      ENDIF

      NTRY=0
 100  NTRY=NTRY+1
      IF(NTRY.GT.20)THEN
         WRITE(LUN,*)' FIRBALL_4FLV: unable to sample 4flv fireball!'
         WRITE(LUN,*)' lacking rejection mechanism, abort..'
         CALL SIB_REJECT ('FIRBALL_4FLV    ')
c         RETURN
      ENDIF

      LA = ABS(L0)
      ISGN = ISIGN(1,L0)
      DELTAE = P0(5) - AM(LA)
      IF(DELTAE.lt.AM(6)+0.02D0)THEN
         IREJ = 1
         IF(ndebug.gt.3)
     &    WRITE(LUN,*)' FIRBALL_4FLV:  too low mass!! aborting...',IREJ
c         xa=-1.
c         xa=log(xa)
c         stop        
         RETURN
      ENDIF
      AV = 2.D0*SQRT(DELTAE)

c...  select number of particles in fireball
c     at least two
 200  XRNDM = GASDEV(LA)
      NPI = INT(AV*(1.D0+PAR(38)*XRNDM))
      XMMIN = AM(LA)+DBLE(NPI-1)*AM(6)+0.02D0
      IF(Ndebug.gt.3)
     &     WRITE(LUN,*)'  NPI,av,rndm,xmin,delta',
     &     NPI,av,XRNDM,xmmin,P0(5)-XMMIN

      IF((NPI.LE.1).OR.(NPI.GT.9).OR.(P0(5).LT.XMMIN))THEN
         GOTO 200
      ENDIF
      IF(Ndebug.gt.3) 
     &  WRITE(LUN,*)' FIRBALL_4FLV: No. of particles sampled. ',
     &  '(NPI,DELTAE,NTRY):',NPI,DELTAE,NTRY

c...  sample particle list      
      NTRYL=0
 210  CONTINUE
c...  special vector resonance treatment with meson projectile
      IF(IPAR(14).eq.-3.and.LA.lt.13)THEN
c     form resonance from meson beam
cdh      IF(NTRY.GT.5) GOTO 211
         IF(NTRY.GT.5) THEN
c     split last hadron again to start hadron chain
           CALL HSPLI (LL(I+1),IFL(1),IFL(2))

           IF(Ndebug.gt.3)
     &          WRITE(LUN,*)' FIRBALL_4FLV: Input hadron split. ',
     &          '(L0,IFL1,IFL2):',LL(I+1),IFL(1),IFL(2)
           WREM = P0(5)
           WREM2 = AM2(ABS(LL(1)))
           INONLEAD(1)=0
           INONLEAD(2)=0
         ELSE
           I=1
           IF(PCHEX.gt.S_RNDM(LA))THEN
              LL(I)=LRESCHEX(LA)
              CALL HSPLI(LCON(LA),IFL1,IFL2)
              IFL(1)=IFL1
              IFL(2)=IFL2
           ELSE
              LL(I)=LRES(LA)
              CALL HSPLI(L0,IFL1,IFL2)
              IFL(1)=-IFL1
              IFL(2)=-IFL2
           ENDIF
           WREM = P0(5)-AM(ABS(LL(1)))
           WREM2 = AM2(ABS(LL(1)))
           INONLEAD(1)=1
           INONLEAD(2)=1
         ENDIF

      ELSE
c...  baryon projectile
c     first two particles defined by charge exchange
         I=1
         LA1=LA
c     add strangeness
         XLIMLAM=sqrt(AM2(35)+AM2(9)+0.4)
         IF(S_RNDM(LA1).lt.PLAM*(1-IABS(ISTR(LA))).and.
     &        DELTAE.gt.XLIMLAM)THEN
            LA1 = LSTR(LA)
c            print *,'xlim<deltae?: ',xlimlam,deltae
            IF(Ndebug.gt.3)
     &write(lun,*)' FIRBALL_4FLV: producing hyperon:',namp(LA),namp(LA1)
         endif        
         IF(PCHEX.gt.S_RNDM(LA1))THEN
            L1=LCON(LA1)
            if(la.eq.42) l1 = l1 + 2 * int(2.D0*S_RNDM(L1))
            LL(I)=L1*ISGN
c            WRITE(LUN,*)' charge exchange!',ISGN*LA,'->',L1
         ELSE
            L1=LA1
            LL(I)=LA1*ISGN
         ENDIF
c     determine remaining charge and strangeness         
         IDQ=ICHP(LA1)*ISGN-ICHP(L1)*ISIGN(1,LL(I))
         IDS=ISTR(LA)*ISGN-ISTR(L1)*ISIGN(1,LL(I))
         IF(ABS(IDQ).gt.1) write(lun,*) 'LA,LA1,L1',LA,LA1,L1
         IF(IABS(IDS).gt.1)
     &        write(lun,*) 'too much strangeness,LA,LA1,L1:'
     &        ,namp(LA),namp(LA1),namp(L1)
         IF(IDS.ne.0)THEN
            IDX = IDS-IDQ
            LL(I+1)=LPICS(IDX)  ! compensate with strange meson if 
         ELSE
            LL(I+1)=LPIC(IDQ)   ! compensate with meson
         ENDIF         
         IF(NPI.eq.2) GOTO 300
c     split last hadron again to start hadron chain
cdh 211     CALL HSPLI (LL(I+1),IFL(1),IFL(2))
         CALL HSPLI (LL(I+1),IFL(1),IFL(2))

         IF(Ndebug.gt.3) 
     &        WRITE(LUN,*)' FIRBALL_4FLV: Input hadron split. ',
     &        '(L0,IFL1,IFL2):',LL(I+1),IFL(1),IFL(2)
         WREM = P0(5)
         WREM2 = AM2(ABS(LL(1)))
         INONLEAD(1)=0
         INONLEAD(2)=0
      ENDIF

      IF(NTRYL.gt.20) GOTO 100
      NTRYL=NTRYL+1

 230  I=I+1    
      JT=INT(1.5D0+S_RNDM(I))
      JR=3-JT
      NTRYS=0
      IFLB=IFL(JT)
      IDM = 5
 240  CALL SIB_I4FLAV (IFL(JT), 0, IDM, IFL(3), LL(I))
      IF(NTRYS.gt.50) GOTO 210    
      NTRYS=NTRYS+1
      W=dEXP(-AM(ABS(LL(I)))/XTEMPH)
      IF(Ndebug.gt.4) 
     &  WRITE(LUN,*)' FIRBALL_4FLV: flavor added: ',
     &  '(I,NTRYS,LL(I),IFL3,W):',I,NTRYS,LL(I),IFL(3),W
      IF(W.LT.S_RNDM(I).and.INONLEAD(JT).eq.1) GOTO 240

c...  kinematic limits...     
      WREM = WREM-AM(IABS(LL(I)))
      WREM2_2=WREM2+2.D0*dSQRT(WREM2)*AM(IABS(LL(I)))+AM2(IABS(LL(I)))
      IF(Ndebug.gt.4) 
     &  WRITE(LUN,*)' FIRBALL_4FLV: kinematic limits: ',
     &  '(I,NTRYS,P05**2,WREM2):',I,NTRYS,P0(5)**2,WREM2_2
      IF(WREM2_2+0.2D0*S_RNDM(I+1).ge.P0(5)**2) GOTO 240
      WREM2=WREM2_2
      IF(Ndebug.gt.3) 
     & WRITE(LUN,*)
     & ' FIRBALL_4FLV: Hadron added: (KF,NAMP,I,NONlead,WRME2)',
     & LL(I),NAMP(ABS(LL(I))),I,INONLEAD(JT),WREM2

      IFL(JT)=-IFL(3)
      INONLEAD(JT)=1
      IF(I.lt.NPI-1) GOTO 230
      IF(ABS(IFL(JT)).gt.3.and.ABS(IFL(JR)).gt.3) THEN
         IFL(JT)=IFLB
         GOTO 240
      ENDIF

c...  close list
      I=I+1
      NTRYC=0
c$$$      IAFL1 = IABS(mod(IFL(JR),100))
c$$$      IAFL2 = IABS(mod(IFL(jt),100))
c$$$      IF ((IAFL1/10.eq.4.or.mod(IAFL1,10).eq.4)
c$$$     +     .and.(IAFL2/10.eq.4.or.mod(IAFL2,10).eq.4))
c$$$     +     GOTO 100             ! reject two charm quarks
c$$$      IF(IAFL1*IAFL2.GT.100)  GOTO 100
 250  CALL SIB_I4FLAV (IFL(JT), IFL(JR), IDM, IFL(3), LL(I))
      IF(NTRYC.gt.10) GOTO 210
      NTRYC=NTRYC+1
      WREM2_2=WREM2+2.D0*dSQRT(WREM2)*AM(ABS(LL(I)))+AM2(ABS(LL(I)))
      IF(Ndebug.gt.5) 
     &  WRITE(LUN,*)' FIRBALL_4FLV: closing List: (IFL1,IFL2,KF,',
     &             'NAMP,I,NTRYC,WREM2)',
     &  IFL(JT),IFL(JR),LL(I),NAMP(ABS(LL(I))),I,NTRYC,WREM2_2

      IF(WREM2_2+0.2D0*S_RNDM(I).ge.P0(5)**2) GOTO 250

 300  IF(Ndebug.gt.3) 
     &     WRITE(LUN,*)
     &     ' FIRBALL_4FLV: flavors sampled. (NPI,LL,WREM,NTRYL):',
     &     NPI,(LL(ii),ii=1,NPI),WREM,NTRYL

c...  fill phasespace
      CALL DECPAR (0,P0,NPI,LL,PD)
      DO J=1,NPI
         NP = NP+1
         LLIST(NP) = LL(J)
         NPORIG(NP) = IPFLAG*2
         niorig(NP)= iiflag
         DO K=1,5
            P(NP,K) = PD(J,K)
         ENDDO
      ENDDO
      PAR(5)=PAR5def
      IREJ = 0
      RETURN
      END
C=======================================================================

      SUBROUTINE SIG_RPP2014(L,KT,SQS,SLOPE,SIGT,SIGEL,SIGINEL,RHO)

C-----------------------------------------------------------------------
C     implementation of the PDG RPP 2014 cross section fit
C     proton-, pion-, kaon-nucleon interactions
C      
c     projectile dependent parameters are stored in amp array
c     dimensions are: (beam,target,exchange mode)
c     cross section is used for interaction length in AIR
c     therefore proton and neutron cross sections are averaged.
c
C     Input:
c     L : beam id (1: proton, 2: pion, 3: kaon)
c     KT: target id (0: Nucleon, 1: proton, 2: neutron)
c     SQS: c.m. energy in GeV
c     SLOPE: fit does not include elastic slope, need input to calc
c            elastic and inelastic cross section
c     Output:      
c     SIGT,SIGEL,SIGINEL,RHO
c     cross sections and ratio of real and imaginary part of ela. amp.
C-----------------------------------------------------------------------
      IMPLICIT NONE
c     external types
      DOUBLE PRECISION SQS,SIGT,SIGEL,SIGINEL,SLOPE,RHO!,SIGDIF
      integer l,kt
c     commons
      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)
      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
c     internal types
      DOUBLE PRECISION S,S0,SIG,RHO1,XI
      INTEGER k,i,INIT
C     universal constants and parameters
      DOUBLE PRECISION M0,ETA1,ETA2,H
      DOUBLE PRECISION AMP(3,2,3)
      DOUBLE PRECISION XMA(3),XMB(2)
      SAVE
      DATA M0,ETA1,ETA2,H /2.076D0,0.412D0,0.5626D0,0.2838D0/
c     hadron-proton
      DATA (AMP(1,1,i),i=1,3) /33.73D0, 13.67D0, 7.77D0 /
      DATA (AMP(2,1,i),i=1,3) /18.08D0, 10.44D0, 1.977D0 /
      DATA (AMP(3,1,i),i=1,3) /15.84D0, 5.12D0, 3.538D0 /
c     hadron-neutron
      DATA (AMP(1,2,i),i=1,3) /33.77D0, 14.05D0, 6.93D0 /
      DATA (AMP(2,2,i),i=1,3) /18.08D0, 10.44D0, 1.977D0 /
      DATA (AMP(3,2,i),i=1,3) /15.73D0, 4.81D0, 1.86D0 /
      DATA INIT/0/
c     particle masses
c     DATA XMA /0.93827D0,0.13957D0,0.493667D0/
c     DATA XMB /0.93827D0,0.939565D0/

      IF(INIT.EQ.0) THEN
c  use the masses from the mass table
        XMA(1) = AM(13)     ! proton
        XMA(2) = AM(7)      ! pi+
        XMA(3) = AM(9)      ! K+
        XMB(1) = AM(13)     ! proton
        XMB(2) = AM(14)     ! neutron
        INIT = 1
      ENDIF

      s = SQS**2
      sigt = 0.D0
      rho = 0.D0
      k = kt
 100  if(kt.eq.0.and.k.lt.2) k = k + 1
      s0=XMA(l)+XMB(k)+M0
      s0=s0**2
      xi=s/s0
c     print *,'s,s0,xi',s,s0,xi
c     print *,'eta1,eta2,h,M0',eta1,eta2,h,M0
c     print *,'P,R1,R2',amp(l,k,1),amp(l,k,2),amp(l,k,3)
c     print *,H*log(xi)**2,amp(l,k,1),amp(l,k,2)*(1.D0/xi)**eta1,
c     &        amp(l,k,3)*(1.D0/xi)**eta2
      sig = H*log(xi)**2+amp(l,k,1)+amp(l,k,2)*(1.D0/xi)**eta1
     &     +amp(l,k,3)*(1.D0/xi)**eta2
c     print *,'sig',sig
c     print *,'pi,0.5D0,0.D0',pi,0.5D0,0.D0
c     print *,pi*h*log(xi),amp(l,k,2)*xi**(-eta1),tan(eta1*pi*0.5D0),
c     &        amp(l,k,3)*xi**(-eta2),(tan(pi*eta2*0.5D0)+EPS5)
      rho1 = PI*h*log(xi)-amp(l,k,2)*xi**(-eta1)*tan(eta1*PI*0.5D0)
     &     +amp(l,k,3)*xi**(-eta2)/(tan(PI*eta2*0.5D0)+EPS5)
c     print *,'rho:',rho1
      rho = rho + rho1/sig
      sigt = sigt + sig
c     write(LUN,*) ' l,k,sig,rho:',l,k,sig,rho
      if(kt.eq.0.and.k.lt.2) goto 100
      if(kt.eq.0) then
         sigt = sigt*0.5D0
         rho = rho*0.5D0
      endif
c     derive elastic and inelastic cross section
      sigel = sigt**2*(1.D0+rho**2)/(16.D0*PI*slope*cmbarn)
      siginel = sigt-sigel
      IF(ndebug.gt.2)
     &  write(LUN,*)
     &  ' SIG_RPP2014: L,KT,SQS,SIGT,SIGEL,SIGINEL,SLOPE,RHO',
     &     L,KT,SQS,SIGT,SIGEL,SIGINEL,SLOPE,RHO
      end
C=======================================================================

      DOUBLE PRECISION FUNCTION FERMI(XARG,X0,XALPH)

C-----------------------------------------------------------------------
C     fermi function, used to smoothen samplings
C     f = 1/(1+exp((x-x0)/alpha))
C-----------------------------------------------------------------------
      IMPLICIT NONE
c     externals
      DOUBLE PRECISION XARG,X0,XALPH,XE
c     COMMONs

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      SAVE

c     internals
      xe = max((xarg-x0)/xalph,-10.D0)
      fermi=1.D0+exp(xe)
      fermi=1.D0/fermi
      END
C=======================================================================
      
      SUBROUTINE SEL_RES(XM2in,KDin,IRDX,IKDH)
      
C--------------------------------------------------------------------
C     routine that checks if excitation should go into resonant state
C     or rather should fallback to on-shell beam hadron
C     Input: XM2in : squared excitation mass
C            KDin : projectile hadron code
C            IRDX : reference to remnant on stack
C     Output: adds hadron to stack
C             IKDH : parton stack index of final hadron
C--------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)

      DOUBLE PRECISION AW,AW2
      COMMON /S_WIDTH1/ AW(99), AW2(99)
      INTEGER          MRES(6:99,2)
      DOUBLE PRECISION XM2,XM1,DELTAE,EMIN1,EMIN2
      INTEGER          KD
      SAVE

      DATA (MRES(k,1),k=6,22)  /27,25,26,28,29,0,0,51,52,6*0,30,31/
      DATA (MRES(k,1),k=23,33) /23,24,25,26,27,28,29,30,31,27,27/
      DATA (MRES(k,1),k=34,49) /34,35,36,37,38,39,40,41,42,43,34,35,36,
     &     37,38,49/
      DATA (MRES(k,1),k=50,83) /0,51,52,53,54,4*0,78,79,10*0,80,81,73,
     &     74,75,76,77,78,79,80,81,0,83/
      DATA (MRES(k,1),k=84,99) /94,95,96,97,98,89,4*0,94,95,96,97,98,99/

      DATA (MRES(k,2),k=6,22)  /61,62,63,64,65,0,0,53,54,6*0,66,67/
      DATA (MRES(k,2),k=23,33) /61,61,62,63,61,64,65,66,67,61,61/
      DATA (MRES(k,2),k=34,49) /34,35,36,37,38,39,40,41,42,43,44,45,46,
     &     47,48,49/
      DATA (MRES(k,2),k=50,83) /0,51,52,53,54,4*0,78,79,10*0,80,81,73,74
     &     ,75,76,77,78,79,80,81,0,83/
      DATA (MRES(k,2),k=84,99) /94,95,96,97,98,89,4*0,94,95,96,97,98,99/

      XM2 = XM2in
      XM1 = sqrt(XM2)
      KD = KDin

C     thresholds
c     fallback threshold
      EMIN1 = PAR(76)
      
c     resonance threshold
      EMIN2 = PAR(77)
      
c     parton stack index of incoming hadron
      IKDH = 0
      
c     if too low, fallback on beam
      IF(ndebug.gt.2)
     &     write(lun,*)' SEL_RES: input (XM2in,KDin,IRDX):',XM2,KD,IRDX
      DELTAE = XM1-AM(ABS(KD))
      IF(ndebug.gt.1)then
         write(lun,*)' SEL_RES: DELTAE,EMIN1,EMIN2',deltae,emin1,emin2
         write(lun,*)' SEL_RES: XM,XM1,XM2',
     &        XM1,emin1+AM(ABS(KD)),emin2+AM(ABS(KD))
      endif
      IF(DELTAE.LT.EMIN1)THEN
c     fallback to beam region
         KDH = kd
         XM1 = AM(abs(kd))
         XM2 = AM2(abs(kd))

      ELSEIF(DELTAE.LT.EMIN2)THEN
c     form resonance
         II = 1
         KDH = KD
         DO WHILE (II.le.2.and.KDH.eq.KD)
            KDD = IABS(KD)
            
c     K0s and K0l projection on K0 and K0bar
cdh         IF(KDD.eq.11.or.KDD.eq.12)KDD=21
cdh  &                              +INT((2.D0-EPS10)*S_RNDM(KD))
            IF(KDD.eq.11.or.KDD.eq.12)KDD=21
     &                              +INT(0.5D0+S_RNDM(KD))
            IL = MRES(KDD,II)
            IF(ndebug.gt.2) then
               write(lun,*) ' SEL_RES: res. select (KD,II,IL):',
     &         KD,II,IL
            ENDif
cdh   to prevent  index of array AW2 out of range
            IF(IL.eq.0) write(lun,*) ' SEL_RES: KD,KDD:' , KD,KDD
            IF(IL.eq.0) CALL SIB_REJECT('SEL_RES         ')
c     sample probability for resonance to occur at this mass
c     from the relativistic breit-wigner dist.
c     scale widths to artificially increase or decrease resonance occurence
            XWDTH = PAR(94)*AW2(IL)
            PRES = BREIT_WIGNER(XM2,AM2(IL),XWDTH)
            IF(ndebug.gt.2)
     &           write(lun,*)
     &           ' SEL_RES: res. proposal (AM2,AW2,Prob.):',
     &           AM2(IL),XWDTH,PRES
            IF(S_RNDM(ii).lt.PRES) KDH = ISIGN(IL,KD)            
            II = II + 1
         ENDDO
c     no resonance selected, fallback to beam or phasespace decay?
         IF(IPAR(59).eq.1.and.KDH.eq.KD)THEN
c     distinguish regions in deltaE
            IF(DELTAE.LT.EMIN1)THEN
c     fallback to beam
               XM1 = AM(abs(kdh))
               XM2 = AM2(abs(kdh))           
            ELSE
               KDH = 0
            ENDIF
         ELSE
c     case where resonance has been selected 
c     or no overlap between resonance and phasespace region exists
c     set mass to pole masses of selected particles
            XM1 = AM(abs(kdh))
            XM2 = AM2(abs(kdh))
         ENDIF
      ELSE
c     neither resonance nor fallback
         KDH = 0
      ENDIF
      IF(KDH.ne.0)THEN
c     add new beam hadron to stack
         XM2in = XM2
         CALL ADD_PRTN
     &        (0.D0,0.D0,0.D0,0.D0,XM1,KDH,2,IRDX,IKDH)
      endif
      IF(ndebug.gt.2)
     &     write(lun,*)' SEL_RES: output (XM2in,KDin,KDH):',XM2,KD,KDH

      RETURN
      END

C=======================================================================

      DOUBLE PRECISION FUNCTION BREIT_WIGNER(S,XM2,XWDTH2)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      SAVE

C     peak set to one
      x1 = (s-xm2)**2+xm2*xwdth2
      breit_wigner = xm2*xwdth2/x1
      end
C=======================================================================

      DOUBLE PRECISION FUNCTION TBREIT_WIGNER(S,XM2,XWDTH2)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      SAVE
C     breit-wigner truncated at 2*gamma from peak
C     peak set to one
      DATA N /10/

      XMLOW = MAX(XM2-N*XWDTH2,0.D0)
      XCUT = SIGN(1.D0,S-XMLOW)
      XCUT = MAX(XCUT,0.D0)
      x1 = (S-xm2)**2+xm2*xwdth2
      TBREIT_WIGNER = xcut * xm2*xwdth2/x1
      
      end
C=======================================================================

      SUBROUTINE FRAG_MINIJET(IDX,IBAD)

C-----------------------------------------------------------------------
C     routine that fragments a gluon - gluon system \FR'14
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER IDX,IBAD

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
C     The final particle output is contained in COMMON /S_PLIST/    
C     NP           : number of final particles
C     P(1:NP, 1:5) : 4-momenta + masses of the final particles 
C     LLIST (1:NP) : codes of final particles
      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP
      INTEGER NW_max
      PARAMETER (NW_max = 20)
C--------------------------------------------------------------------
C     SIBYLL common blocks containing event information       \FR'14
C--------------------------------------------------------------------

C     EVENT INFO COMMON
C     contains overall interaction properties, like
C     SQS : center-of-mass energy
C     S   :         "       "     squared
C     PTmin : low pt cut of QCD cross section, 
C             i.e. minimal pt of hard minijets
C     Xmin : low-x bound for PDFs, 
C            i.e. minimal momentum fraction of hard partons
C     Zmin : logarithm of that
C     KB : PID of beam hadron
C     KT() : PID of target
C     IAT : mass number of target
      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
      INTEGER KB,IAT,KT
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
C     parameters that represent: NW: max. number of wounded nucleons,
C     NS,NH: max. number of soft and hard interactions
c      PARAMETER (NW_max = 20)
C     The COMMON block /S_CHIST/ contains information about the
C     the structure of the  generated event:
C     NWD   = number of wounded nucleons
C     NJET = total number of hard interactions
C     NSOF = total number of soft interactions
C     NNSOF (1:NW) = number of soft pomeron cuts in each interaction
C     NNJET (1:NW) = number of minijets produced in each interaction 
C     JDIF(1:NW) = diffraction code 
C                  0 : non-diff,
C                  1 : beam-diff
C                  2 : target-diff
C                  3 : double-diff
      INTEGER NNSOF,NNJET,JDIF,NWD,NJET,NSOF
      COMMON /S_CHIST/ NNSOF(NW_max),NNJET(NW_max),
     &     JDIF(NW_max),NWD,NJET,NSOF

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      
      DOUBLE PRECISION PGG,PST,PBM,PTG,E0,PT2JET,PTJET,TH,FI,S_RNDM,
     &     PAR1_def,PAR24_def,PAR3_def,PAR2_1_def,PAR2_2_def,PAR5_def,
     &     PAR6_def,PAR24_2_def,XM,QMASS,DBETJ      
      DIMENSION PST(5),PBM(5),PTG(5)
      INTEGER IST,ITGST,IBMST,IPID,IFLB,IFLT,NOLD,IS,IFL1,IFBAD,IDM,
     &     ipar82_def
      SAVE
      DATA PGG /1.D0/

C     read partons from stack
c     references are string --> bm-parton --> tg-parton
c     read string 4momentum from stack
      CALL RD_PRTN_4VEC(IDX,PST,IPID,IBMST)
      CALL RD_PRTN_4VEC(IBMST,PBM,IFLB,ITGST)
      CALL RD_PRTN_4VEC(ITGST,PTG,IFLT,IST)
      IF(IDX.ne.IST) then
         write(lun,*) ' FRAG_MINIJET: reference loop broken!' , IDX
         CALL SIB_REJECT('FRAG_MINIJET    ')
      endif

C..   kinematic variables
      E0 = PST(5)            ! string mass
      PT2JET = PBM(1)**2 + PBM(2)**2
      PTJET = sqrt(PT2JET)
      TH = ASIN(MIN((1.D0-EPS8),2.D0*PTJET/E0))
c      FI = ASIN(MIN((1.D0-EPS8),PBM(2)/PTJET))
      FI = TWOPI*S_RNDM(IDX)
c      TH = PST(1)
c      FI = PST(2)

      IF(NDEBUG.gt.1) WRITE(LUN,*)' FRAG_MINIJET: IDX,EE,IFLB,IFLT,PT',
     &     IDX,E0,IFLB,IFLT,PTJET,IBAD
      IF(NDEBUG.gt.1) WRITE(LUN,*)' FRAG_MINIJET: PTJET,TH,FI:',
     &     PTJET,TH,FI

C...  parameter setup (string fragmentation)

c     baryon production setup
      PAR1_def = PAR(1)
      if( NSOF+NJET.gt.0) then
         PAR(1)= PAR(15)
      else
         PAR(1)= PAR(14)
      endif
            
C...  charm setup
      PAR24_def = PAR(24)
      IF(IPAR(15).eq.2.or.IPAR(15).eq.3)THEN
         PAR(24) = PAR(25)*EXP(-PAR(26)/E0)
      ELSEIF(IPAR(15).eq.4)THEN
         PAR(24) = PAR(27)*EXP(-PAR(26)/E0)
      ELSEIF(IPAR(15).eq.5)THEN
         PAR(24) = PAR(27)*EXP(-PAR(26)/E0)
         PAR(29) = PAR(27)*EXP(-PAR(28)/E0)
      ELSEIF(IPAR(15).eq.6.or.IPAR(15).eq.8.or.IPAR(15).eq.9.or.
     &        IPAR(15).eq.11)THEN
         PAR(24) = PAR(27)*EXP(-PAR(28)/E0)
      ELSEIF(IPAR(15).eq.7)THEN
         PAR(24) = PAR(27)
      ELSEIF(IPAR(15).eq.10)THEN
         WRITE(LUN,*)' FRAG_minijet: charm model not implemented!'
         CALL SIB_REJECT('FRAG_minijet    ')
      ENDIF

C...  strange setup
      PAR2_1_def = PAR(2)
      PAR3_def = PAR(3)
      IF(IPAR(42).eq.1)THEN
c     change to constant value 
         PAR(2) = PAR(72)
      ELSEIF(IPAR(42).eq.2)THEN
c     change according to string mass, saturating
         PAR(2) = PAR(72)*EXP(-PAR(73)/E0)
      ELSEIF(IPAR(42).eq.3)THEN
c     change strange diq fraction as well
         PAR(2) = PAR(72)       ! P_s / P_ud
         PAR(3) = PAR(73)       ! P_us / P_ud
      ENDIF

C...  vector setup
      PAR5_def = PAR(5)
      PAR6_def = PAR(6)
      IF(IPAR(43).eq.1)THEN
c     change vector rate and kaon vector rate
         PAR(5) = PAR(74)       ! P_vec
         PAR(6) = PAR(74)       ! P_K* from K
         
      ENDIF
      
C...  switch off pi0 suppression
c     should only be applied for remnant, diff and valence
c     in case of meson projectile
      ipar82_def = IPAR(82)
      IF(IPAR(95).eq.1)THEN
         IPAR(82) = 0
      ENDIF
      
      NOLD = NP
      IF ( (E0.LT.8.D0) .OR. (S_RNDM(0).GT.PGG)) THEN
C...  one string case, q - qbar
         
C     sample flavor for q-qbar minijet        
         IF( IPAR(87).eq.3 )THEN
C     flavor threshold model            
c     u,d -> u,d,s -> u,d,s,c
c     s and transition from massive to massless at m_s and m_c thresholds
c     beyond the charm mass all flavors are equally likely
            CALL SIB_ICFLAV(E0**2,0,IDM,IFL1)
            
         ELSE
C     default u,d,s model, same rates as in hadronization (string frag.)
            PAR2_2_def = PAR(2)
            PAR24_2_def = PAR(24)
C     set 'leading' strange fraction         
            IF(IPAR(39).eq.2) PAR(2) = PAR(66)         
c     leading charm fraction
            IF( IPAR(87).eq.1 )THEN
               PAR(24) = PAR(150)
            ELSEIF( IPAR(87).eq.2 )THEN
               PAR(24) = PAR(150)*PAR(24)
            ENDIF

            IS = -1 + 2*INT((2.D0-EPS8)*S_RNDM(0))
 100        IFL1 = IS*(INT((2.D0+PAR(2))*S_RNDM(0))+1)
            XM = 2.D0*QMASS(IFL1)+0.3D0
            if(E0.LE.XM) GOTO 100
            IF(IABS(IFL1).eq.3)THEN
               IF(S_RNDM(IFL1).lt.PAR(24)*PAR(125))IFL1 = IS*4
               XM = 2.D0*QMASS(IFL1)+0.3D0
               if(E0.LE.XM) GOTO 100
            ENDIF
            PAR(2) = PAR2_2_def
            PAR(24) = PAR24_2_def                        
         ENDIF
      
         CALL STRING_FRAG_4FLV 
     &        (E0,IFL1,-IFL1,0.D0,0.D0,0.D0,0.D0,IFBAD,0)
         if(IFBAD.gt.0) then
            IF(ndebug.gt.1)
     &       WRITE(LUN,*)
     &           ' JET_FRAG: rejection in STRING_FRAG (IFL,E0,NCALL):',
     &           IFL1,E0,NCALL
            PAR(24) = PAR24_def
            PAR(1) = PAR1_def
            PAR(2) = PAR2_1_def
            PAR(5) = PAR5_def
            PAR(6) = PAR6_def
            PAR(3) = PAR3_def
            IPAR(82) = ipar82_def       
            RETURN
         ENDIF
      ELSE
C...  two string case, gluon - gluon
         CALL GG_FRAG_4FLV(E0)
      ENDIF

c      DBETJ = (DX1J-DX2J)/(DX1J+DX2J)
      DBETJ = PST(3)/PST(4)
      CALL SIROBO (NOLD+1,NP,TH,FI,0.D0,0.D0,DBETJ)

      if(Ndebug.gt.1) WRITE(LUN,*)
     &     ' JET_FRAG: particles produced:',NP-NOLD
      PAR(24) = PAR24_def
      PAR(1) = PAR1_def
      PAR(2) = PAR2_1_def
      PAR(5) = PAR5_def
      PAR(6) = PAR6_def
      PAR(3) = PAR3_def
      IPAR(82) = ipar82_def  
      IBAD = 0
      END
C=======================================================================

      SUBROUTINE INT_H_NUC (IA, SIGT, SLOPE, RHO) 

C-----------------------------------------------------------------------
C...Compute with a montecarlo method the "multiple interaction structure"
C.  of an hadron-nucleus collision.
C.  
C.
C.  INPUT : IA               = mass of target nucleus
C.          SIGT (mbarn)     = total hp cross section
C.          SLOPE (GeV**-2)  = slope of hp elastic scattering
C.          RHO              = real/imaginary part of forward elastic
C.                             scattering amplitude
C.
C.  OUTPUT : in COMMON block /CNCMS0/
C.           B = impact parameter (fm)
C.           BMAX = maximum impact parameter for generation
C.           NTRY = number of "trials" before one interaction
C.           NA = number of wounded nucleons in A
C. Author : P.Lipari  (may 1993)
C---------------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      PARAMETER (IAMAX=56)
      COMMON /S_CNCM0/ B, BMAX, NTRY, NA
      DIMENSION XA(IAMAX), YA(IAMAX)

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      SAVE


      PI=4.d0*atan(1.d0)

      CC = SIGT/(4.D0*PI*SLOPE*CMBARN)         
      DEN = 2.D0*SLOPE*CMBARN*0.1D0
      BMAX = 1.D0*10.D0             ! fm
      NTRY = 0
      CALL NUC_CONF (IA, XA, YA)
1000  CONTINUE
      B = BMAX*dSQRT(S_RNDM(0))
      PHI = 2.D0*PI*S_RNDM(NTRY)
      BX = B*DCOS(PHI)
      BY = B*DSIN(PHI)
      NTRY = NTRY + 1
      NA = 0
      DO JA=1,IA
         S = (XA(JA)-BX)**2 + (YA(JA)-BY)**2
         F = dEXP(-S/DEN)
         PEL = CC*CC*(1.D0+RHO*RHO)*F*F
         PINEL  = 2.D0*CC*F-PEL
         R = S_RNDM(JA)
         IF (R .LT. PINEL)  THEN
            NA = NA + 1
         ENDIF
      ENDDO
      IF (NA .EQ. 0 .and. NTRY .lt. 1000)  GOTO 1000

      RETURN
      END
C=======================================================================

      SUBROUTINE SIB_REJECT(text)

C-----------------------------------------------------------------------
c     subroutine dumps state of random number generator 
c     at beginning of event to file then produces fpe/stops
C----------------------------------------------------------
      IMPLICIT NONE

      character*16  text
      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      INTEGER II2,JJ2
      DOUBLE PRECISION U2,C2,CD2,CM2
      COMMON /SIB_RAND/ U2(97),C2,CD2,CM2,II2,JJ2
      INTEGER NW_max
      PARAMETER (NW_max = 20)
C--------------------------------------------------------------------
C     SIBYLL common blocks containing event information       \FR'14
C--------------------------------------------------------------------

C     EVENT INFO COMMON
C     contains overall interaction properties, like
C     SQS : center-of-mass energy
C     S   :         "       "     squared
C     PTmin : low pt cut of QCD cross section, 
C             i.e. minimal pt of hard minijets
C     Xmin : low-x bound for PDFs, 
C            i.e. minimal momentum fraction of hard partons
C     Zmin : logarithm of that
C     KB : PID of beam hadron
C     KT() : PID of target
C     IAT : mass number of target
      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
      INTEGER KB,IAT,KT
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT
      DOUBLE PRECISION XDM
c      CHARACTER*13 FILENA
      SAVE
c      DATA FILENA /'sib_rjctn.rnd'/

      WRITE(LUN,*)
     &     ' SIB_REJECT:(from,ncall,KB,iat,ECM) ',
     &                   text,ncall,kb,iat,sqs
c     produce floating point error
      XDM = -1.D0
      XDM = LOG(XDM)
      STOP
      END
C=======================================================================

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
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
c      COMMON /S_DEBUG/ Ncall, Ndebug, Lun
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      INTEGER NS_max, NH_max
      PARAMETER (NS_max = 20, NH_max = 80)
      
      DOUBLE PRECISION SSIG,PJETC,SSIGN,SSIGNSD,SSIGNEL,ALINT,ASQSMIN,
     &     ASQSMAX,DASQS
      INTEGER NSQS
      COMMON /S_CCSIG/ SSIG(61,3), PJETC(0:NS_max,0:NH_max,61,2),
     &     SSIGN(61,3,3), SSIGNSD(61,3,3), SSIGNEL(61,3,3), 
     &     ALINT(61,3,3), ASQSMIN, ASQSMAX, DASQS, NSQS
      DOUBLE PRECISION STR_mass_val, STR_mass_val_hyp, STR_mass_sea
      COMMON /S_CUTOFF/ STR_mass_val, STR_mass_val_hyp, STR_mass_sea

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      SAVE

C     check if tables initialized
      IF(NSQS.eq.0) THEN
         WRITE(LUN,*) ' CUT_PRO: tables not initialized! aborting...'
         xa = -1.D0
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

      R = (1.D0-EPS8)*S_RNDM(0)
      DO I=0,NS_max
        DO J=0,NH_max
          IF (R.LT.(1.D0-T)*PJETC(I,J,J1,K)+T*PJETC(I,J,J2,K)) GOTO 100
        ENDDO
      ENDDO
100   CONTINUE

C...phase space limitation

 120  CONTINUE
      XM = DBLE(2*I)*STR_mass_sea + DBLE(2*J)*PTmin
      PACC = EXP(PAR(9)*(2.D0-XM)/SQS)
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
     &  write(lun,*)' CUT_PRO: (L,SQS,PTmin,Ns,Nh) ',K,SQS,PTmin,I,J

      END

C=======================================================================

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
      IMPLICIT INTEGER(I-N)

      INTEGER NS_max, NH_max
      PARAMETER (NS_max = 20, NH_max = 80)
      
      DOUBLE PRECISION SSIG,PJETC,SSIGN,SSIGNSD,SSIGNEL,ALINT,ASQSMIN,
     &     ASQSMAX,DASQS
      INTEGER NSQS
      COMMON /S_CCSIG/ SSIG(61,3), PJETC(0:NS_max,0:NH_max,61,2),
     &     SSIGN(61,3,3), SSIGNSD(61,3,3), SSIGNEL(61,3,3), 
     &     ALINT(61,3,3), ASQSMIN, ASQSMAX, DASQS, NSQS
      DOUBLE PRECISION SSIG_TOT,SSIG_SD1,SSIG_SD2,SSIG_DD,SSIG_B,
     &     SSIG_RHO
      COMMON /S_CCSIG2/ SSIG_TOT(61,3),SSIG_SD1(61,3),SSIG_SD2(61,3),
     &    SSIG_DD(61,3),SSIG_B(61,3),SSIG_RHO(61,3)
      DOUBLE PRECISION SSIG_SD1LM,SSIG_SD1HM,SSIG_SD2LM,SSIG_SD2HM,
     &     SSIG_DDLM,SSIG_DDHM
      COMMON /S_CCSIG3/ SSIG_SD1LM(61,3),SSIG_SD1HM(61,3),
     &     SSIG_SD2LM(61,3),SSIG_SD2HM(61,3),
     &     SSIG_DDLM(61,3),SSIG_DDHM(61,3)
      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN

      DIMENSION Pjet(0:NS_max,0:NH_max)
      DIMENSION SIG_df(3),SIG_df2(3,2),SIGDIF(3),SIGDIF_pi(3),
     &          PS_tab(61),PH_tab(61),PT_tab(61)

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      SAVE

C...spacing in energy for table of cross sections.

      NSQS = 61
      ASQSMIN = 1.D0
      ASQSMAX = 7.D0
      DASQS = (ASQSMAX-ASQSMIN)/DBLE(NSQS-1)

C...initialization of proton and pion tables
      
      IF(LUN.ne.6) WRITE(6,*)' Calculating cross section tables...'
      DO KK=1,2

         IF(NDEBUG.gt.0)
     &    WRITE(LUN,'(2(/,1X,A,A))') 
     &     'Table: J, sqs,  PT_cut,  SIG_tot, SIG_inel, B_el,  ',
     &     'rho,    <n_s>,  <n_h>, SIG_SD, SD1_lm, SD1_hm',
     &     '---------------------------------------------------',
     &     '----------------------------------------------'

         JINT = KK
         DO J=1, NSQS
           ASQS = ASQSMIN + DASQS*DBLE(J-1)
           SQS = 10.D0**ASQS

           CALL SIB_SIG (JINT, SQS, PTmin,
     &                   SIG_tot, SIG_inel, SIG_df, SIG_df2, B_el, Pjet)

C...low-energy interpolation with data-parametrizations
           CALL SIB_HADCSL(JINT,SQS,
     &                     SIGTOT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO)
           if(SQS.le.100.D0) then
             SIG_TOT  = SIGTOT
             SIG_inel = SIGINEL
             B_EL     = SLOPE
           else if(SQS.le.1000.D0) then
             Xi = dlog(SQS/100.D0)/2.30258509299405D0
             SIG_TOT  = Xi*SIG_TOT+(1.D0-Xi)*SIGTOT
             SIG_inel = Xi*SIG_inel+(1.D0-Xi)*SIGINEL
             B_EL     = Xi*B_EL+(1.D0-Xi)*SLOPE
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

           PSUM = 0.D0
           PH = 0.D0
           PS = 0.D0
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

           IF(NDEBUG.gt.0)
     &      WRITE(LUN,'(3X,I2,1P,E12.3,0P,4F8.2,6F8.3)') 
     &       JINT,SQS,PTmin,SIG_tot,SIG_inel,B_el,RHO,PS,PH
     &          ,SIGDIF(1)+SIGDIF(2),SIG_df2(1,1),SIG_df2(1,2)

         ENDDO
      ENDDO

C...initialization of kaon tables

      JINT = 3

      IF(NDEBUG.gt.0)
     & WRITE(LUN,'(2(/,1X,A,A))') 
     &  'Table: J, sqs,  PT_cut,  SIG_tot, SIG_inel, B_el,  ',
     &  'rho,    <n_s>,  <n_h>',
     &  '---------------------------------------------------',
     &  '---------------------'
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
        CALL SIB_HADCSL(2,SQS,
     &                  SIGTOT_pi,SIGEL_pi,SIGINEL,SIGDIF_pi,SLOPE,RHO)
        CALL SIB_HADCSL(3,SQS,
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
          SIG_TOT  = Xi*SIG_TOT+(1.D0-Xi)*SIGTOT
          SIG_inel = Xi*SIG_inel+(1.D0-Xi)*SIGINEL
          B_EL     = Xi*B_EL+(1.D0-Xi)*SLOPE
        endif

        SSIG_TOT(J,3) = SIG_TOT
        SSIG(J,3)     = SIG_inel
        SSIG_SD1(J,3) = SIGDIF(1)
        SSIG_SD2(J,3) = SIGDIF(2)
        SSIG_DD(J,3)  = SIG_df(3)
        SSIG_B(J,3)   = B_EL
        SSIG_RHO(J,3) = RHO

        IF(NDEBUG.gt.0)
     &   WRITE(LUN,'(3X,I2,1P,E12.3,0P,4F8.2,3F8.3)') 
     &    JINT,SQS,PTmin,SIG_tot,SIG_inel,B_el,RHO,PS,PH

      ENDDO

      END

C=======================================================================

      SUBROUTINE INI_WRITE (LUN)

C-----------------------------------------------------------------------
C   This subroutine prints on unit LUN
C   a table of the cross sections  used in the program
C   and of the average number of hard interactions, and the average
C   number of wounded nucleons in a hadron-air interaction
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NS_max, NH_max
      PARAMETER (NS_max = 20, NH_max = 80)
      
      DOUBLE PRECISION SSIG,PJETC,SSIGN,SSIGNSD,SSIGNEL,ALINT,ASQSMIN,
     &     ASQSMAX,DASQS
      INTEGER NSQS
      COMMON /S_CCSIG/ SSIG(61,3), PJETC(0:NS_max,0:NH_max,61,2),
     &     SSIGN(61,3,3), SSIGNSD(61,3,3), SSIGNEL(61,3,3), 
     &     ALINT(61,3,3), ASQSMIN, ASQSMAX, DASQS, NSQS
      DIMENSION PJ(2),PS(2),PW(2)

      SAVE
      DATA ATARGET /14.514D0/

      if ( ndebug .gt. 3 ) CALL PARAM_PRINT(LUN)
      if ( ndebug .gt. 0 ) THEN
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

           PW(K) = ATARGET*SSIG(J,K)/SSIGN(J,K,1)

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

         WRITE(LUN,20) SQS,SSIG(J,1),SSIGN(J,1,1),PS(1),PJ(1),PW(1)
     &                      ,SSIG(J,2),SSIGN(J,2,1),PS(2),PJ(2),PW(2)

        ENDDO

        WRITE(LUN, 18)
      endif
20    FORMAT(1p,E10.2,2(2F7.1,1X,3F6.2,1X))

      return
      END

C=======================================================================

            SUBROUTINE SIG_AIR_INI 

C-----------------------------------------------------------------------
C...  Initialize the cross section and interaction lengths on air, 
C.    nitrogen and oxygen      
C.  (this version initializes p-air, pi-air, and K-air cross sections)
C.
C.  also calculates the low mass beam diffraction cross section in hAir \FR
C.  using the same lambda for all hadrons
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NS_max, NH_max
      PARAMETER (NS_max = 20, NH_max = 80)
      
      DOUBLE PRECISION SSIG,PJETC,SSIGN,SSIGNSD,SSIGNEL,ALINT,ASQSMIN,
     &     ASQSMAX,DASQS
      INTEGER NSQS
      COMMON /S_CCSIG/ SSIG(61,3), PJETC(0:NS_max,0:NH_max,61,2),
     &     SSIGN(61,3,3), SSIGNSD(61,3,3), SSIGNEL(61,3,3), 
     &     ALINT(61,3,3), ASQSMIN, ASQSMAX, DASQS, NSQS
      COMMON /GLAUB_SCR/ XI_MAX , ALAM(61)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      DIMENSION SIGDIF(3)
      DIMENSION ITARGC(3)
      CHARACTER*3 TARGN
      DIMENSION TARGN(3)
      SAVE
      DATA AVOG /6.0221367D-04/
      DATA ATARGET /14.514D0/
      DATA ITARGC /0,14,16/
      DATA TARGN /'air','nit','oxy'/

      IF ( IPAR(12).GT.0 ) THEN
         if (ndebug.gt.0) then
           WRITE(LUN,*) ' SIG_AIR_INI:'
           WRITE(LUN,*)' using Goulianos param. for res.coupling..'
         endif
         XI_MAX = 0.02D0
         if (ndebug.gt.0)WRITE(LUN,*)' low mass Xi_max: ' , XI_MAX
      ENDIF

C...target loop (air, N, O)      
      DO IK=1,3
         IAT = ITARGC(IK)
         WRITE(6,*) 'SIG_AIR_INI: initializing target: (i,A)',
     &    ik, IAT, TARGN(IK) , '..'
C...particle loop (p, pi, K)
         DO K=1,3
         
            if (NDEBUG .gt. 0 ) then
               WRITE(6,'(/,1X,A,A)') 
     &          'Table: J, IK, sqs,    SIGtot,     SIGprod,    SIG_SD,',
     &'     Lambda  '
               WRITE(6,*) 
     &              '-------------------------------------------------',
     &              '-------------'
            endif
            DO J=1,NSQS

               ASQS = ASQSMIN + DASQS*DBLE(J-1)
               SQS = 10.D0**ASQS

               IF (K.EQ.1) THEN
c     Goulianos param. from GAP-2012-056, Mx**2s = 0.02
c     against PDG elastic cross section
                  CALL SIB_HADCS1
     &             (K,SQS,SIGT1,SIGEL1,SIGINEL1,SLOPE1,RHO1)
                  SIGEFF = 0.68D0*(1.D0+36.D0/SQS**2)*
     &                 dlog(0.6D0+XI_MAX/1.5D0*SQS**2)
                  ALAM(J) = dSQRT(SIGEFF/SIGEL1)
               ENDIF
               CALL SIB_SIGMA_HP(K,SQS,
     &              SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO)
               IF(IK.eq.1)THEN
c     fixed O-N mixture
                  CALL SIG_H_AIR
     &                 (SIGT, SLOPE, RHO, ALAM(J),
     &                 SSIGT, SSIGEL, SSIGQE, SSIGSD, SSIGQSD)
               ELSE
                  CALL SIG_H_NUC
     &                 (IAT, SIGT, SLOPE, RHO, ALAM(J),
     &                 SSIGT, SSIGEL, SSIGQE, SSIGSD, SSIGQSD)
               ENDIF
               if (ndebug .gt. 0 ) WRITE(6,'(1X,I2,1P,5E12.3)') 
     &              K,SQS,SSIGT,SSIGT-SSIGQE,SSIGQSD,ALAM(J)
C     particle production cross section
               SSIGN(J,K,IK) = SSIGT-SSIGQE
c     diffractive cross section               
               SSIGNSD(J,K,IK) = SSIGQSD
c     elastic cross section
               SSIGNEL(J,K,IK) = SSIGEL
c     interaction length
               IF(IK.eq.1)then
                  ALINT(J,K,IK) = 1.D0/(AVOG*SSIGn(j,K,IK)/ATARGET)
               else
                  ALINT(J,K,IK) = 1.D0/(AVOG*SSIGn(j,K,IK)/IAT)
               endif
            ENDDO
         ENDDO
     
         if (ndebug .gt. 0 ) then
            WRITE(6,'(/,1X,A)') 
     &           ' SIG_AIR_INI: NUCLIB interaction lengths [g/cm**2]'
            WRITE(6,*) 'target:', TARGN(IK)
            WRITE(6,'(1X,A)') 
     &           '     sqs,       p-targ,      pi-targ,     K-targ'
            DO J=1,NSQS
               ASQS = ASQSMIN + DASQS*DBLE(J-1)
               SQS = 10.D0**ASQS
               WRITE(6,'(1X,1P,4E12.3)') 
     &              SQS,ALINT(J,1,IK),ALINT(J,2,IK),ALINT(J,3,IK)
            ENDDO
         endif
      ENDDO
      END
C=======================================================================

      SUBROUTINE SAMPLE_TARGET(NW,XCHG,KRMNT,XJET,Irec,IREJ)

C-----------------------------------------------------------------------/
C...Subroutine to sample valence and sea quark kinematic variables
C     on the target side
C.    fills IFLT,X2 and PXT,PYT
C.    1,2 are valence quarks, 3,4 are additional sea quarks
C.    transverse momentum is shared between the val. and sea pairs
C.    X and flv are exchanged occasionally, not pt so far
C-------------------------------------------------------------------      
      IMPLICIT NONE

      INTEGER NW_max
      PARAMETER (NW_max = 20)
c     external types
      DOUBLE PRECISION XJET,XCHG
      DIMENSION XJET(NW_max)
      INTEGER KRMNT,NW,IREC,IREJ
      DIMENSION KRMNT(NW_max)


      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      
      DOUBLE PRECISION STR_mass_val, STR_mass_val_hyp, STR_mass_sea
      COMMON /S_CUTOFF/ STR_mass_val, STR_mass_val_hyp, STR_mass_sea
C--------------------------------------------------------------------
C     SIBYLL common blocks containing event information       \FR'14
C--------------------------------------------------------------------

C     EVENT INFO COMMON
C     contains overall interaction properties, like
C     SQS : center-of-mass energy
C     S   :         "       "     squared
C     PTmin : low pt cut of QCD cross section, 
C             i.e. minimal pt of hard minijets
C     Xmin : low-x bound for PDFs, 
C            i.e. minimal momentum fraction of hard partons
C     Zmin : logarithm of that
C     KB : PID of beam hadron
C     KT() : PID of target
C     IAT : mass number of target
      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
      INTEGER KB,IAT,KT
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT
      INTEGER NS_max, NH_max
      PARAMETER (NS_max = 20, NH_max = 80)

      INTEGER IBMRDX,ITGRDX,IHMJDX,ISMJDX,ICSTDX,IINTDX
      COMMON /S_INDX/ IBMRDX(3),ITGRDX(NW_max,3),
     &     IHMJDX(NW_max*NH_max),IINTDX(NW_max),
     &     ISMJDX(NW_max*NS_max),ICSTDX(2*NW_max,3)

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN

      INTEGER IRMNT,KRB,KRT
      DOUBLE PRECISION XRMASS,XRMEX
      COMMON /S_RMNT/ XRMASS(2),XRMEX(2),IRMNT(NW_max),KRB,KRT(NW_max)

c     internal types
      DOUBLE PRECISION XX,X2,PX,PXT,PY,PYT,PZ,PZ1,PZ2
      DIMENSION XX(2*NW_max+2),PX(2*NW_max+2),PY(2*NW_max+2)
      DIMENSION X2(4*NW_max),PXT(4*NW_max),PYT(4*NW_max)
      INTEGER IFL,IFLT,IREJ1,J,J1,J2,J3,J4,JJ,JJJ,JI,I,KID,Iref1,
     &     Iref,KID1
      DIMENSION IFL(2*NW_max+2),IFLT(4*NW_max)
      SAVE

      IREJ1 = 1

      IF(ndebug.gt.2) 
     +     WRITE(LUN,*)
     +     ' SAMPLE_TARGET: NW,XCHG,LRMNT,XJET,IREC,IREJ',
     +     NW,XCHG,(KRMNT(j),j=1,NW),(XJET(j),j=1,NW),IREC,IREJ

      DO J=1,NW ! zero arrays
         j1 = 1+4*(j-1)
         j2 = j1 + 1
         j3 = j2 + 1
         j4 = j3 + 1
         X2(j1) = 0.D0
         X2(j2) = 0.D0
         X2(j3) = 0.D0
         X2(j4) = 0.D0
         PXT(j1) = 0.D0
         PXT(j2) = 0.D0
         PXT(j3) = 0.D0
         PXT(j4) = 0.D0
         PyT(j1) = 0.D0
         PyT(j2) = 0.D0
         PyT(j3) = 0.D0
         PyT(j4) = 0.D0
      ENDDO

      DO j=1,NW
c     read target id from event info 
         KID = KT(J)
c     reset rejection
         IREJ = IREJ1
c     always fills remnant partons into 1,2 and c.strings into 3,4
c     so far only one interaction possible (beam is always a single hadron!)        
         CALL SAMPLE_PROJECTILE
     +        (KID,1,KRMNT(j),XCHG,XJET(j),XX,PX,PY,IFL,KID1,IREJ)
         IF(IREJ.ne.0) RETURN

c     write to target variables
         do jj=3-2*KRMNT(j),4
            ji = jj+4*(j-1)
            IFLT(ji) = IFL(jj)
            X2(ji) = XX(jj)
            PXT(ji) = PX(jj)
            PYT(ji) = PY(jj)
         enddo

         IF(KRMNT(j).ne.0)THEN
c     by convention hadron is split such that diq is 2nd flv
c     for string frag routine argument flv1 is along +z, flv2 -z
c     by convention again flv2 in the remnant is passed to +z and flv1 to -z
c     therefor on the target side the flavors need to be switched such that
c     the diq is along -z
            j1 = 1+4*(j-1)
            j2 = j1 + 1
            CALL ISWTCH_LMNTS(IFLT(j1),IFLT(j2))
         ENDIF

c     central strings
c     flavors need to be switched as well (strictly speaking color)
c     in dual-parton model: q : color , diq : anticolor
c     need to combine q with diq for color neutral system..
         j3 = 3+4*(j-1)
         j4 = j3 + 1
         CALL ISWTCH_LMNTS(IFLT(j3),IFLT(j4))
         CALL SWTCH_LMNTS(X2(j3),X2(j4))
         
c     reset remnant id 
c     might have changed in flavor exchange (actually color)...
         KRT(J) = KID1
      ENDDO

C..   write target partons to stack
      DO I=1,NW
         IF(KRMNT(I).ne.0)THEN
c     add proto-remnant
            j1 = 1+4*(i-1)
            j2 = j1 + 1
            CALL ADD_PRTN(PXT(J1)+PXT(J2),PYT(J1)+PYT(J2),
     &           -0.5D0*SQS*(X2(J1)+X2(j2)),0.5D0*SQS*(X2(J1)+X2(j2)),
     &           0.D0,-2,0,0,Iref1)
            ITGRDX(I,1) = Iref1
            CALL ADD_INT_REF(Iref1,IINTDX(I))
c     add quarks to stack
            do j = 1,2
               jj = 4*(i-1)+j
               jjj = 4*(i-1)+j + 2
               pz1 = (0.5D0*SQS*X2(JJ))**2
c               PZ1 = (0.5D0*SQS*X2(JJ))**2-PXT(JJ)**2-PYT(JJ)**2
               CALL ADD_PRTN(PXT(JJ),PYT(JJ),-sqrt(pz1),
     &              0.5D0*SQS*X2(JJ),0.D0,IFLT(JJ),1,Iref1,Iref)
               ITGRDX(I,j+1) = Iref
               pz2 = (0.5D0*SQS*X2(JJj))**2
c               pz2 = (0.5D0*SQS*X2(JJj))**2-PXT(JJj)**2-PYT(JJj)**2
               CALL ADD_PRTN(PXT(JJj),PYT(JJj),-sqrt(pz2),
     &              0.5D0*SQS*X2(JJj),0.D0,IFLT(JJj),1,0,Iref)
               ICSTDX(2*(I-1)+j,3) = Iref
            enddo
         else
            do j = 3,4
               jj = 4*(i-1)+j
               pz = (0.5D0*SQS*X2(JJ))**2
c               pz = (0.5D0*SQS*X2(JJ))**2-PXT(JJ)**2-PYT(JJ)**2
               CALL ADD_PRTN(PXT(JJ),PYT(JJ),-sqrt(pz),
     &              0.5D0*SQS*X2(JJ),0.D0,IFLT(JJ),1,0,Iref)
               ICSTDX(2*(I-1)+(J-2),3) = Iref
            enddo
         ENDIF
      ENDDO
      IF(NDEBUG.GT.3) CALL PRNT_PRTN_STCK

      IREJ = 0
      END
C=======================================================================

      SUBROUTINE SIGMA_NUC_AIR(IA,ECM,KINT)

C-----------------------------------------------------------------------
C.  wrapping for SIGMA_NUC in NUCLIB
C...Compute with a montecarlo method the "production"
C.  and "quasi-elastic" cross section for  
C.  a nucleus-nucleus interaction
C.  nucleon - nucleon cross section is taken from 
C.  the table calculated by SIBYLL_INI
C.
C.  INPUT : IA            = mass of target nucleus
C.          ECM          = c.m. energy
C.          KINT            = number  of interactions to generate
C.  OUTPUT : SIGMA (mbarn) = "production" cross section
C.           DSIGMA   "    = error
C.           SIGQE    "    = "quasi-elastic" cross section
C.           DSIGQE   "    = error
C.      in COMMON /NUCNUCSIG/ 
C.           additional output is in the common block  /CPROBAB/
C.           Prob(n_A), Prob(n_B), Prob(n_int)
C..........................................................................
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON /NUCNUCSIG/ SIGPROD,DSIGPROD,SIGQE,DSIGQE,IBE,ITG
      DIMENSION SIGDIF(3)
      SAVE
      DATA NDB /0/
      
      DSIGPROD = 0.D0
      DSIGQE = 0.D0

      CALL SIB_SIGMA_HP(1,ECM,SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO)
      CALL SIGMA_AIR(IA,SIGINEL,SIGEL,KINT,SIGPROD,DSIGPROD,
     +     SIGQE,DSIGQE)
      IBE = IA
      ITG = 0
      IF(DSIGPROD/SIGPROD.gt.0.1D0)THEN
         IF( NDB.EQ.0 ) 
     +     PRINT*,'SIG_NUC_AIR: warning! : large error in cross section'
         NDB = 1
      ENDIF
      RETURN
      END

C=======================================================================

      SUBROUTINE SIG_NUC_AIR(IA,SIGPP,SIGPPEL,KINT)

C-----------------------------------------------------------------------
C.  wrapping for SIGMA_NUC in NUCLIB
C...Compute with a montecarlo method the "production"
C.  and "quasi-elastic" cross section for  
C.  a nucleus-nucleus interaction
C.
C.  INPUT : IA            = mass of target nucleus
C.          IB            = mass of projectile nucleus
C.          SIGPP (mbarn)  = inelastic pp cross section
C.          SIGPPEL        = elastic pp cross section
C.          KINT            = number  of interactions to generate
C.  OUTPUT : SIGMA (mbarn) = "production" cross section
C.           DSIGMA   "    = error
C.           SIGQE    "    = "quasi-elastic" cross section
C.           DSIGQE   "    = error
C.      in COMMON /NUCNUCSIG/ 
C.           additional output is in the common block  /CPROBAB/
C.           Prob(n_A), Prob(n_B), Prob(n_int)
C..........................................................................
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON /NUCNUCSIG/ SIGPROD,DSIGPROD,SIGQE,DSIGQE,IBE,ITG
      SAVE

      DSIGPROD = 0.D0
      DSIGQE = 0.D0
      CALL SIGMA_AIR(IA,SIGPP,SIGPPEL,KINT,SIGPROD,DSIGPROD,
     +     SIGQE,DSIGQE)
      IBE = IA
      ITG = 0
      IF(DSIGPROD/SIGPROD.gt.0.1D0)THEN
         IF( NDB.EQ.0 ) 
     +     PRINT*,'SIG_NUC_AIR: warning! : large error in cross section'
         NDB = 1
      ENDIF
      RETURN
      END

C=======================================================================

      SUBROUTINE SIG_NUC_NUC(IA,IB,SIGPP,SIGPPEL,KINT)

C-----------------------------------------------------------------------
C.  wrapping for SIGMA_NUC in NUCLIB
C...Compute with a montecarlo method the "production"
C.  and "quasi-elastic" cross section for  
C.  a nucleus-nucleus interaction
C.
C.  INPUT : IA            = mass of target nucleus
C.          IB            = mass of projectile nucleus
C.          SIGPP (mbarn)  = inelastic pp cross section
C.          SIGPPEL        = elastic pp cross section
C.          KINT            = number  of interactions to generate
C.  OUTPUT : SIGMA (mbarn) = "production" cross section
C.           DSIGMA   "    = error
C.           SIGQE    "    = "quasi-elastic" cross section
C.           DSIGQE   "    = error
C.      in COMMON /NUCNUCSIG/ 
C.           additional output is in the common block  /CPROBAB/
C.           Prob(n_A), Prob(n_B), Prob(n_int)
C..........................................................................
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON /NUCNUCSIG/ SIGPROD,DSIGPROD,SIGQE,DSIGQE,IBE,ITG
      SAVE

      DSIGPROD = 0.D0
      DSIGQE = 0.D0
      CALL SIGMA_MC(IA,IB,SIGPP,SIGPPEL,KINT,SIGPROD,DSIGPROD,
     +     SIGQE,DSIGQE)
      IBE = IA
      ITG = IB
      IF(DSIGPROD/SIGPROD.gt.0.1D0)THEN
         IF( NDB.EQ.0 ) 
     +     PRINT*,'SIG_NUC_NUC: warning! : large error in cross section'
         NDB = 1
      ENDIF
      RETURN
      END
C=======================================================================

      SUBROUTINE SIG_HAD_NUC(L,IA,ECM,ALAM,ICSMOD,IPARM)

C-----------------------------------------------------------------------
C**********************************************************************
C...Subroutine to compute hadron-nucleus cross sections
C.  according to:
C.  R.J. Glauber and G.Matthiae  Nucl.Phys. B21, 135, (1970)
C.
C.
C.  INPUT :  L projectile particle (1:p , 2:pi, 3:K )
C.           IA mass-number of target nucleus
C.           SSIG  (mbarn) total pp cross section
C.           SLOPE (GeV**-2)  elastic scattering slope for pp
C.           ALPHA    real/imaginary part of the forward pp elastic
C.                                               scattering amplitude
C.           ALAM: inel. screening coupling
C.
C.  OUTPUT : ( in COMMON block /NUCSIG/ )
C.           SIGT  = Total cross section
C.           SIGEL = Elastic cross section
C.           SIGQEL  = Elastic + Quasi elastic cross section
C.           SIGSD  = beam single diff. cross section
C......................................................................
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON /NUCSIG/ SIGT,SIGEL,SIGINEL,SIGQE,SIGSD,
     +     SIGQSD,SIGPPT,SIGPPEL,SIGPPSD,ITG
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)
c      double precision dplab
c      double precision DSSIG,DSLOPE,DALPHA,DALAM
c      DOUBLE PRECISION SG1,SGEL1,SGQE1,SGSD1,SGQSD1
      DIMENSION SSIGDIF(3),XM(4)     
      SAVE
c     DATA XM / 0.93956563, 0.13956995, 0.493677, 0.93956563 /
      DATA GEV2MB /0.3893D0/
      DATA INIT/0/

      IF(INIT.EQ.0) THEN
c  use the masses from the mass table
cdh     XM(1) = AM(14)     ! neutron
        XM(1) = AM(13)     ! proton
        XM(2) = AM(7)      ! pi+
        XM(3) = AM(9)      ! K+
        XM(4) = AM(14)     ! neutron
        INIT = 1
      ENDIF

      xma = XM(L)
      xmb = (XM(1)+XM(4))/2.D0

      Plab = dsqrt(((ecm**2-xma**2-xmb**2)/(2.D0*xmb))**2-xma**2)

C     hadron proton cross section to be used for calculation

      IF( ICSMOD.EQ.1 ) THEN
c     sibyll 2.1 cross section

         CALL SIB_SIGMA_HP(L,ECM,SSIG,SSIGEL,SSIGINEL,SSIGDIF,SLOPE,RHO)

      ELSEIF( ICSMOD.EQ.0 ) THEN
c     cross section parametrizations

         if(Ecm.gt.12.D0) then

           CALL SIB_HADCSL(L,ECM,SSIG,SSIGEL,SSIGINEL,SSIGDIF,SLOPE,RHO)

         else
c     low energy parametrization
            SSIG = (sigtot_pp(Plab)+sigtot_pn(plab))/2.D0
            SSIGEL  = (sigela_pp(Plab)+sigela_pn(plab))/2.D0
C     parametrization from U. Dersch et al. Nucl Phys. B579 (2000) 277
            RHO = 6.8D0/plab**0.742D0-6.6D0/plab**0.599D0+0.124D0
            SLOPE = (1.D0+RHO**2)*SIGTOT**2/(16.D0*PI*SIGEL)/GEV2MB
            SSIGDIF(1) = 0.D0
            SSIGDIF(2) = 0.D0
            SSIGDIF(3) = 0.D0
         endif
      ENDIF
      SSIGSD = SSIGDIF(1) + SSIGDIF(2)

c     energy dependence of lambda parameter
      if( IPARM.eq.1 ) then

c     empirical parametrization
         SIGEFF = 0.25D0*Ecm**2/(Ecm**2+10.D0**2)*dLOG(1000.D0*Ecm**2)
     &        -1.5D0/2.D0
         SIGEFF = MAX(0.D0,SIGEFF)
         
         ALAM = dsqrt(SIGEFF/SSIGEL)

         SSIGSD = 2.D0 * SIGEFF
         
      elseif( IPARM.EQ.2 ) then
         
c     lambda derived from proton interactions
         CALL SIB_HADCS1(1,ECM,SIGT1,SSIGEL1,SIGINEL1,SLOPE1,RHO1)
C     parametrization by Goulianos for diff. interaction
         SIGEFF = 0.68D0*(1.D0+36.D0/Ecm**2)
     &        *LOG(0.6D0+0.02D0/1.5D0*Ecm**2)
         SIGEFF = MAX(0.D0,SIGEFF)
         ALAM = sqrt(SIGEFF/SSIGEL1)
         
         SSIGSD = 2.D0 * SIGEFF
         
      elseif( IPARM.eq.3)then

C     data from Paolo Lipari's note
         SIGTOT = 129.D0
         SIGEL  = 0.3D0*SIGTOT
         SIGEFF = ECM*0.01D0*SIGTOT
         RHO    = 0.D0
         SLOPE  = (1.D0+RHO**2)*SIGTOT**2/(16.D0*PI*SIGEL)/GEV2MB
         ALAM   = dsqrt(SIGEFF/SIGEL)
         
         SSIG = SIGTOT
         SSIGEL = SIGEL
         SSIGSD = 2.D0 * SIGEFF
      endif 

      ALPHA = RHO

C     hadron - nucleon cross section
      
      IF( IA.EQ.0 ) THEN
         CALL SIG_H_AIR
     +        (SSIG,SLOPE,ALPHA,ALAM,SG1,SGEL1,SGQE1,SGSD1,SGQSD1)
      else
         CALL GLAUBER2
     +        (IA,SSIG,SLOPE,ALPHA,ALAM,SG1,SGEL1,SGQE1,SGSD1,SGQSD1)
      endif

      ITG = IA

      SIGPPT = SSIG
      SIGPPEL = SSIGEL
      SIGPPSD = SSIGSD
      SIGT  = SG1
      SIGEL = SGEL1
      SIGQE = SGQE1
      SIGSD = SGSD1
      SIGQSD = SGQSD1
      SIGINEL = SIGT - SIGEL

      RETURN
      END
C=======================================================================

      SUBROUTINE SIG_H_AIR
     +     (SSIG,SLOPE,ALPHA,ALAM,SIGT,SIGEL,SIGQE,SIGSD,SIGQSD)

C-----------------------------------------------------------------------
C**********************************************************************
C...Subroutine to compute hadron-air cross sections
C.  according to:
C.  R.J. Glauber and G.Matthiae  Nucl.Phys. B21, 135, (1970)
C.
C.  Air is a linear combination of Nitrogen and oxygen
C.
C.  INPUT :  SSIG  (mbarn) total pp cross section
C.           SLOPE (GeV**-2)  elastic scattering slope for pp
C.           ALPHA    real/imaginary part of the forward pp elastic
C.                                               scattering amplitude
C.           ALAM  coupling to inel. intermediat states
C.  OUTPUT : SIGT  = Total cross section
C.           SIGEL = Elastic cross section
C.           SIGQEL  = Elastic + Quasi elastic cross section
C.           SIGSD   = single diff. cross section (beam) 
C.           SIGQSD  = Elastic + Quasi elastic SD cross section (beam)
C.
C.  ALSO including interface from single precision in SIBYLL to
C.       double precision in GLAUBER2
C......................................................................
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      SAVE
      DATA FOX /0.21522D0/  !atomic percentage of 'non-nitrogen' in air

      CALL GLAUBER2
     +  (14,SSIG,SLOPE,ALPHA,ALAM,SIG1,SIGEL1,SIGQE1,SIGSD1,SIGQSD1)
      CALL GLAUBER2
     +  (16,SSIG,SLOPE,ALPHA,ALAM,SIG2,SIGEL2,SIGQE2,SIGSD2,SIGQSD2)

      SIGT  = (1.D0-FOX)*SIG1   + FOX*SIG2
      SIGEL = (1.D0-FOX)*SIGEL1 + FOX*SIGEL2
      SIGQE = (1.D0-FOX)*SIGQE1 + FOX*SIGQE2
      SIGSD = (1.D0-FOX)*SIGSD1 + FOX*SIGSD2
      SIGQSD = (1.D0-FOX)*SIGQSD1 + FOX*SIGQSD2
      RETURN
      END

C=======================================================================

      SUBROUTINE SIG_H_NUC
     +     (IAT,SSIG,SLOPE,ALPHA,ALAM,SIGT,SIGEL,SIGQE,SIGSD,SIGQSD)

C-----------------------------------------------------------------------
C**********************************************************************
C...Subroutine to compute hadron-nucleus cross sections
C.  according to:
C.  R.J. Glauber and G.Matthiae  Nucl.Phys. B21, 135, (1970)
C.
C.  INPUT :  IAT   nucleon number in target nucleus
C.           SSIG  (mbarn) total pp cross section
C.           SLOPE (GeV**-2)  elastic scattering slope for pp
C.           ALPHA    real/imaginary part of the forward pp elastic
C.                                               scattering amplitude
C.  OUTPUT : SIGT  = Total cross section
C.           SIGEL = Elastic cross section
C.           SIGQEL  = Elastic + Quasi elastic cross section
C.           SIGSD   = single diff. cross section (beam) 
C.           SIGQSD  = Elastic + Quasi elastic SD cross section (beam)
C.
C......................................................................
      IMPLICIT NONE
      INTEGER IAT
      DOUBLE PRECISION SSIG,SLOPE,ALPHA,ALAM      
      DOUBLE PRECISION SIG1,SIGEL1,SIGQE1,SIGSD1,SIGQSD1
      DOUBLE PRECISION SIGT,SIGEL,SIGQE,SIGSD,SIGQSD
      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      SAVE
      IF(IAT.eq.0.or.IAT.gt.18) THEN
         WRITE(LUN,'(//,1X,A)') 
     &        ' SIG_H_NUC: number of target nucleons too large!',
     &        ' (1<=IAT<=18)'
         SIGT = -1.D0
         STOP
      ENDIF

      CALL GLAUBER2
     +  (IAT,SSIG,SLOPE,ALPHA,ALAM,SIG1,SIGEL1,SIGQE1,SIGSD1,SIGQSD1)
      SIGT  = SIG1   
      SIGEL = SIGEL1 
      SIGQE = SIGQE1 
      SIGSD = SIGSD1 
      SIGQSD = SIGQSD1
      RETURN
      END

C=======================================================================
      
      SUBROUTINE GLAUBER2
     +     (JA,SSIG,SLOPE,ALPHA,ALAM,SIGT,SIGEL,SIGQEL,SIGSD,SIGQSD)

C-----------------------------------------------------------------------
C...Subroutine to compute hadron-Nucleus cross sections
C.  according to:
C.  R.J. Glauber and G.Matthiae  Nucl.Phys. B21, 135, (1970)
C.
C.  This formulas assume that the target nucleus  density is
C.  modeled by a shell-model form.  A reasonable range of models
C.  is  4 < JA < 18
C.
C.  This is a modified version with a two-channel model for inelastic
C.  intermediate states of low mass (R. Engel 2012/03/26)
C.
C.  INPUT :  A = mass number of the nucleus
C.           SSIG  (mbarn) total pp cross section
C.           SLOPE (GeV**-2)  elastic scattering slope for pp
C.           ALAM  enhancement factor (sqrt of sigma_sd1/sigma_ela)
C.           ALPHA    real/imaginary part of the forward pp elastic
C.                                               scattering amplitude
C.  OUTPUT : SIGT  = Total cross section
C.           SIGEL = Elastic cross section
C.           SIGQEL  = Elastic + Quasi elastic cross section
C.           SIGSD = single diff. cross section
C.           SIGQSD = Quasi single diff. cross section
C.
C. Internally  everything is computed in GeV (length = GeV**-1)
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON /CA0SH/ R0, R02
      COMPLEX*16  ZS1, ZS2, ZP1, ZP2, Z1, Z2, OM12
      DIMENSION RR(18)
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN
      SAVE
      DATA BMAX /100.D0/            ! GeV**-1
      DATA NB /500/
C...data on Sqrt[<r**2>] (fm). (A=5,8 are not correct).
C   From Barett and Jackson
      DATA RR /0.81,2.095,1.88,1.674, 2.56,2.56,2.41,2.5,2.519,2.45
     +          ,2.37, 2.460, 2.440, 2.54, 2.58, 2.718, 2.662,2.789 /

      A = DBLE(JA)
C...Parameter of shell model density
      R0 = RR(JA)/0.197D0/dSQRT(5.D0/2.D0 - 4.D0/A)    ! GeV**-1
      R02 = R0*R0

      SIG1 = (1.D0+ALAM) * SSIG/CMBARN            ! GeV**-2
      SIG2 = (1.D0-ALAM) * SSIG/CMBARN
      SIG12 = dSQRT((1.D0+ALAM)*(1.D0-ALAM)) * SSIG/CMBARN
      DB = BMAX/DBLE(NB)
      SUM0 = 0.D0
      SUM1 = 0.D0
      SUM2 = 0.D0
      SUM3 = 0.D0
      SUM4 = 0.D0
      DO JB=1,NB

        B = DB*(DBLE(JB)-0.5D0)

        GS1 = GLAUBGS_D (B,SLOPE, SIG1)
        XS1 = (1.D0- GS1)
        YS1 = GS1*ALPHA
        ZS1 = DCMPLX(XS1,YS1)

        GP1 = GLAUBGP_D (B,SLOPE, SIG1)
        XP1 = (1.D0- GP1)
        YP1 = GP1*ALPHA
        ZP1 = DCMPLX(XP1,YP1)

        Z1 = ZS1**4 * ZP1**(A-4.D0)

        GS2 = GLAUBGS_D (B,SLOPE, SIG2)
        XS2 = (1.D0- GS2)
        YS2 = GS2*ALPHA
        ZS2 = DCMPLX(XS2,YS2)

        GP2 = GLAUBGP_D (B,SLOPE, SIG2)
        XP2 = (1.D0- GP2)
        YP2 = GP2*ALPHA
        ZP2 = DCMPLX(XP2,YP2)

        Z2 = ZS2**4 * ZP2**(A-4.D0)

        XZ = 0.5D0 * DREAL(Z1+Z2)
        YZ = 0.5D0 * DIMAG(Z1+Z2)

        XZ2 = 0.5D0 * DREAL(Z2-Z1)
        YZ2 = 0.5D0 * DIMAG(Z2-Z1)

        SUM0 = SUM0 + (1.D0-XZ)*B

        SUM1 = SUM1 + ((1.D0-XZ)**2 + YZ**2)*B

        SUM3 = SUM3 + (XZ2**2 + YZ2**2)*B

        OMS1 = OMEGAS_D(B,SIG1,SLOPE,ALPHA)
        OMS2 = OMEGAS_D(B,SIG2,SLOPE,ALPHA)
        OMS12 = OMEGAS_D(B,SIG12,SLOPE,ALPHA)

        OMP1 = OMEGAP_D(B,SIG1,SLOPE,ALPHA)
        OMP2 = OMEGAP_D(B,SIG2,SLOPE,ALPHA)
        OMP12 = OMEGAP_D(B,SIG12,SLOPE,ALPHA)

        OM1 = (1.D0 - 2.D0*GS1 + OMS1)**4
     &      * (1.D0 - 2.D0*GP1 + OMP1)**(A-4.D0)
        OM2 = (1.D0 - 2.D0*GS2 + OMS2)**4
     &      * (1.D0 - 2.D0*GP2 + OMP2)**(A-4.D0)
        OM12 = (1.D0 - GS1*DCMPLX(1.D0,ALPHA)-GS2*DCMPLX(1.D0,-ALPHA)
     &               + OMS12)**4
     &       * (1.D0 - GP1*DCMPLX(1.D0,ALPHA)-GP2*DCMPLX(1.D0,-ALPHA)
     &               + OMP12)**(A-4.D0)
        SUM2 = SUM2 + (1.D0-2.D0*XZ + (OM1+OM2)/4.D0
     &                 + DREAL(OM12)/2.D0)*B
        SUM4 = SUM4 + ((OM1+OM2)/4.D0
     &                 - DREAL(OM12)/2.D0)*B

      ENDDO

      SIGT =   SUM0 * DB * 4.D0*PI * CMBARN
      SIGEL =  SUM1 * DB * TWOPI * CMBARN
      SIGQEL = SUM2 * DB * TWOPI * CMBARN
      SIGSD =  SUM3 * DB * TWOPI * CMBARN
      SIGQSD = SUM4 * DB * TWOPI * CMBARN
      END

C=======================================================================

      FUNCTION GLAUBGS_D (B,SLOPE, SIG)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON /CA0SH/ A0, A02
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN
      SAVE

      GAMMA2 = A02/4.D0 + 0.5D0*SLOPE
      ARG = B**2/(4.D0*GAMMA2)
      GLAUBGS_D = SIG/(8.D0*PI*GAMMA2) * EXP(-ARG)
      RETURN
      END

C=======================================================================

      FUNCTION GLAUBGP_D (B,SLOPE, SIG)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON /CA0SH/ A0, A02
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN
      SAVE

      GAMMA2 = A02/4.D0 + 0.5D0*SLOPE
      ARG = B**2/(4.D0*GAMMA2)
      C1 = 1.D0- A02/(6.D0*GAMMA2)*(1.D0-ARG)
      GLAUBGP_D = SIG/(8.D0*PI*GAMMA2) *  C1 * EXP(-ARG)
      RETURN
      END

C=======================================================================

      FUNCTION OMEGAS_D (B, SIG, SLOPE, RHO)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON /CA0SH/ A0, A02
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN
      SAVE

      ETA2 = 0.25D0*(A02 + SLOPE)
      F02 = SIG*SIG*(1.D0+RHO*RHO)/(16.D0*PI**2)
      ARG = -B*B/(4.D0*ETA2)
      OMEGAS_D = F02/(4.D0*ETA2*SLOPE) *EXP(ARG)
      RETURN
      END

C=======================================================================

      FUNCTION OMEGAP_D (B, SIG, SLOPE, RHO)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON /CA0SH/ A0, A02
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN
      SAVE

      ETA2 = 0.25D0*(A02 + SLOPE)
      F02 = SIG*SIG*(1.D0+RHO*RHO)/(16.D0*PI**2)
      ARG = -B*B/(4.D0*ETA2)
      OMEGAP_D=F02/(4.D0*ETA2*SLOPE)*(1.D0-A02/(6.D0*ETA2)*(1.D0+ARG))
     $                                         *EXP(ARG)
      RETURN
      END
C=======================================================================

      SUBROUTINE REMOVE_PI0(XRATE,N1,N2)

C-----------------------------------------------------------------------
C     routine to exchange pi0 on stack with charged pions
C     violating charge conservation.
C     final pions will be off-shell
C      
C     Input: exchange rate and stack positions inbetween
C     which pions shall be exchanged.
C     
C---------------------------------------------------------     
      IMPLICIT NONE
c     Commons
      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
C     The final particle output is contained in COMMON /S_PLIST/    
C     NP           : number of final particles
C     P(1:NP, 1:5) : 4-momenta + masses of the final particles 
C     LLIST (1:NP) : codes of final particles
      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
C     external types
      DOUBLE PRECISION XRATE
      INTEGER N1,N2
C     internals
      INTEGER I,LL,LA,IFPI0
      DOUBLE PRECISION S_RNDM
      SAVE

      IF(NDEBUG.gt.0)write(lun,*)
     &               ' REMOVE_PI0: Rate,Mode:',xrate,IPAR(50)
C     select exchange model      
      IF(IPAR(50).eq.1)THEN
C     stack loop     
         DO I=N1,N2
            LL = MOD(LLIST(I),10000)
            LA = IABS(LL)
c     IF(LA.eq.6)THEN
            IFPI0=(1-MIN(IABS(1-LA/6),1))*MAX(1-MOD(LA,6),0)
c     replace with pi+ or pi-
            LL=LL+IFPI0*(2-INT(MIN((2.D0+XRATE)*S_RNDM(LA),
     &                                                 3.D0-EPS10)))
            LLIST(I) = LL
            IF(NDEBUG.gt.1)
     &           WRITE(LUN,*) ' REMOVE_PI0: LA,IFPI0,LNEW:',LA,IFPI0,LL
         ENDDO
      ENDIF         
      END
C=======================================================================

      SUBROUTINE SAMPLE_SEA_INDV(KRMNT,XMINA,XMINA_SEA,NSEA,
     &     XREM0,ALPHA,ASUP,XQMASS,XMAX,XX,IREJ)

C-----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      INTEGER NW_max
      PARAMETER (NW_max = 20)

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      INTEGER ITRY, NREJ
      COMMON /S_CNT/ ITRY(20), NREJ(20)

      DOUBLE PRECISION XMINA,XMINA_SEA,XREM0,ALPHA,ASUP,XQMASS,XMAX
      INTEGER NSEA,KRMNT
      DOUBLE PRECISION XX
      DIMENSION XX(2*NW_max+2)
      INTEGER IREJ

      DOUBLE PRECISION XREM,XKIN,X1,X2,pt,S_RNDM,XQM
      INTEGER ICNT2,J,jj1,jj2
      SAVE
      DATA ICNT2 /0/
      
      IF(ndebug.gt.2)
     &    write(lun,*)' SAMPLE_SEA_INDV: called with ',
     &     '(KRMNT,XMINA,XMINA_SEA,NSEA,XREM0,ALPHA,ASUP,XQMASS,XMAX):',
     &     KRMNT,XMINA,XMINA_SEA,NSEA,XREM0,ALPHA,ASUP,XQMASS,XMAX
      XREM = 0.D0
      XKIN = 0.1D0
      XQM = XQMASS
      ITRY(4) = 0
      DO WHILE ( XREM .lt. XMINA )
         XREM = XREM0
         IF ( XREM .LT. 2.D0*XMINA + Nsea*XMINA_SEA
     &        +XKIN*(1.5D0-S_RNDM(ICNT2)) ) THEN
            IREJ = 2            ! resample event
            RETURN
         ENDIF
         IF(ITRY(4).gt.Nsea/2*NREJ(4))THEN
            ICNT2 = ICNT2 + 1
            IF(ndebug.gt.2)THEN
               IF(ICNT2.le.5)THEN
                  write(lun,*)' SAMPLE_SEA_INDV: rejection!' 
                  write(lun,*)' reached max. no. of trials!', NREJ(4)
                  write(lun,*)' XREM0,N,XMIN:' ,XREM0,Nsea,XMINA_SEA
               ENDIF
               IF(ICNT2.eq.5) 
     &              write(lun,*)' last warning of this type..'
            ENDIF
            IREJ = IPAR(51)
            RETURN
         ENDIF
         DO j=1,Nsea/2
c     scale for interactions other than first if Nw>1
            IF(IPAR(75).eq.1.and.J.gt.1) XQM = XQM*PAR(118)
            CALL SAMPLE_SEA(ALPHA,ASUP,XQM,XMAX,x1,x2,pt)
            jj1 = 2 + 2*(j-1) + 1
            IF(KRMNT.eq.0) jj1 = 4+2*(j-1) + 1
            jj2 = jj1 + 1
            XX(jj1) = x1
            XX(jj2) = x2
            XREM = XREM - XX(jj1) - XX(jj2)
            IF(NDEBUG.gt.2) 
     &           WRITE(LUN,*) '  x-frac: JW,X3,X4,XREM',
     &           J,XX(jj1),XX(jj2),XREM
         ENDDO
         ITRY(4) = ITRY(4) + 1
         IF(NDEBUG.gt.1) WRITE(LUN,*) 
     &        ' SAMPLE_SEA_INDV: ISMPL,XREM0,XREM,XMINA,XMINSEA',
     &        ITRY(4),XREM0,XREM,XMINA,XMINA_SEA
      ENDDO
      XREM0 = XREM
      IREJ = 0
      END
C=======================================================================

      SUBROUTINE FORCE_BARYONS(XRATE1,IEDEP,ETHR,XFEXP,IPROJ,N1,N2)

C-----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
C     The final particle output is contained in COMMON /S_PLIST/    
C     NP           : number of final particles
C     P(1:NP, 1:5) : 4-momenta + masses of the final particles 
C     LLIST (1:NP) : codes of final particles
      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP
      INTEGER NW_max
      PARAMETER (NW_max = 20)
      INTEGER LLIST1
      COMMON /S_PLIST1/ LLIST1(8000)
C--------------------------------------------------------------------
C     SIBYLL common blocks containing event information       \FR'14
C--------------------------------------------------------------------

C     EVENT INFO COMMON
C     contains overall interaction properties, like
C     SQS : center-of-mass energy
C     S   :         "       "     squared
C     PTmin : low pt cut of QCD cross section, 
C             i.e. minimal pt of hard minijets
C     Xmin : low-x bound for PDFs, 
C            i.e. minimal momentum fraction of hard partons
C     Zmin : logarithm of that
C     KB : PID of beam hadron
C     KT() : PID of target
C     IAT : mass number of target
      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
      INTEGER KB,IAT,KT
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)

      INTEGER NFORIG,NPORIG,NIORIG,IPFLAG,IIFLAG,KINT
      COMMON /S_PARTO/ NFORIG(NP_max),NPORIG(NP_max),NIORIG(NP_max),
     &IPFLAG,IIFLAG,KINT

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN

c     external types
      double precision xrate1, ethr, xfexp
      integer n1,n2,iedep,iproj

c     internal types
      integer ncomb
      parameter (ncomb=7)
      integer kk,i1,i2,i3,kba,ineutron,incrmax,incr1,incr2,id1,id2,
     &     iflip,icount,ipart,icomb,ismeson
      DIMENSION ipart(ncomb,3),ismeson(99)
      double precision S_RNDM,xrate,ecm_max,pmin,xmass2
     &     ,ptot,ptmp,p1,p2,p3,gamma,gambet,xm,enew,pnew,ph,pn,p1new
     &     ,p2new,xf,xff,pts,pz2,px,py,p1newt,p2newt
c     &     ,sthp,cthp,pt1p,p2pt,p1p,p1pt,p3p,p3pt,pt2p,p2p
      
      DIMENSION ptot(5),ptmp(5),p1(5),p2(5),p3(5),
     &    gambet(4),p1new(5),p2new(5),pn(4)!,p2p(4),p3p(4),p1p(4)
      EQUIVALENCE (gambet(4),gamma)
      SAVE
      DATA IPART /6,7,8,6,6,7,8, 6,8,7,6,7,8,7, 0,0,0,6,8,6,6/
      DATA ISMESON /5*0,7*1,87*0/
      DATA ECM_MAX /1.36987D+05/ ! 10^19eV in ecm and gev for protons
      DATA PMIN /0.05D0/          ! minimal momentum in GeV
      
      KBA = IABS(KB)
      
c     energy dependence
      if(iedep.eq.0) then
         xrate = xrate1
      elseif(iedep.eq.1)then
         xrate=xrate1 * max(0.D0,log(SQS/ethr)/log(ecm_max/ethr))
      endif
      if(iproj.eq.1) then
c     meson projectiles only
         if(ismeson(kba).ne.1) return         
      endif      

c      CALL SIB_LIST(6)
c      print*,'xrate', xrate
      icount = 0
c     check if candidate (pi0, pi+, pi-)      
c      icomb = 1
c     loop over particle stack
      DO I1=N1,N2
c         print*,'checking index:',I1
         DO icomb=1,ncomb
c            print*,'**** checking index:',I1
c            print*,'combination:',icomb
            IF(llist(I1).eq.ipart(icomb,1)) then
c     next neighbours         
               incrmax = max(0,min(5,MIN(N2,NP)-I1))
c               print*,'max increment:',incrmax
               do incr1=1,incrmax
                  I2 = I1 + incr1
                  if(llist(I2).ne.ipart(icomb,2)) goto 300 ! next particle2
                  if(ipart(icomb,3).eq.0)then
c     two-particle combination
c                     print*,' 2 particle combination found!'
c                     print*,' indices:' ,i1, i2
c                     print*,' ids:', llist(i1), llist(i2)
c     sufficient mass?
                     do kk=1,5
                        p1(kk) = p(i1,kk)
                        p2(kk) = p(i2,kk)
                        p3(kk) = 0.D0
                     enddo
                     call ADD_4VECS(P1,P2,ptot)
                     call four_length(ptot,xmass2)
                     if(xmass2.le.4*(pmin**2+am2(13))) then
c                       print*,'insufficient mass!'
                        goto 300 ! next particle2
                     else
                        goto 500 ! implement exchange
                     endif
                     
                  endif
                  do incr2=1,incrmax
                     if(incr1.eq.incr2) goto 400 ! next particle3
                     I3 = I1 + incr2
c     print*,'checking..',i1,i2,i3
                     if(llist(I3).ne.ipart(icomb,3)) goto 400 ! next particle3
c     pi0,pi+,pi- combination found!
c                     print*,' combination found!'
c                     print*,' indices:' ,i1, i2, i3
c                     print*,' ids:', llist(i1), llist(i2), llist(i3)
c     sufficient mass?
                     do kk=1,5
                        p1(kk) = p(i1,kk)
                        p2(kk) = p(i2,kk)
                        p3(kk) = p(i3,kk)
                     enddo
                     call ADD_4VECS(P1,P2,ptmp)
                     call ADD_4VECS(ptmp,P3,ptot)
                     call four_length(ptot,xmass2)
                     if(xmass2.le.4*(pmin**2+am2(13))) then
c                        print*,'insufficient mass!'
                        goto 400 ! next particle3
                     endif                     
 500                 continue                     
                     xm = sqrt(xmass2)
c                     print*,'-------'
c                     print*,' sufficient mass! mass', xm                     
c     calc. initial momenta in 2-particle center-of-mass (eg. p1'--> p1p)
c     use these to calculate the angles which then define pt after the exchange
                     do kk=1,4
                        gambet(kk) = ptot(kk)/xm
                     enddo
c                     print*,'ptot:', ptot
c                     print*,'gambet',gambet
CC                     call sib_altra(gambet(4),-gambet(1),-gambet(2)
CC     &                    ,-gambet(3),p1(1),p1(2),p1(3),p1(4),
CC     &                    p1pt,p1p(1),p1p(2),p1p(3),p1p(4))
CC                     call sib_altra(gambet(4),-gambet(1),-gambet(2)
CC     &                    ,-gambet(3),p2(1),p2(2),p2(3),p2(4),
CC     &                    p2pt,p2p(1),p2p(2),p2p(3),p2p(4))
CC                     call sib_altra(gambet(4),-gambet(1),-gambet(2)
CC     &                    ,-gambet(3),p3(1),p3(2),p3(3),p3(4),
CC     &                    p3pt,p3p(1),p3p(2),p3p(3),p3p(4))                     
c                     print*,'p1:',p1
CC                     pt1p = sqrt(p1p(1)**2+p1p(2)**2)
CC                     pt2p = sqrt(p2p(1)**2+p2p(2)**2)
c                     print*,'p1p:',p1p,sqrt(p1p(4)**2-p1p(3)**2-pt1p**2)
c                     print*,'p2:',p2
c                     print*,'p2p:',p2p
c                     print*,'p3:',p3
c                     print*,'p3p:',p3p
c                     print*,'pt prime, p prime:',pt1p,p1pt
c                     print*,'pt 2 prime, p 2 prime:',
c     &                    pt2p,p2pt
c                     print*,'pt 3 prime, p 2 prime:',
c     &                    sqrt(p3p(1)**2+p3p(2)**2),p3pt                     
CC                     sthp = min(pt1p/p1pt,pt2p/p2pt) ! sin(theta')
CC                     cthp = max(p1p(3)/p1pt,p2p(3)/p2pt) ! cos(theta')
c                     print*,'sthp,cthp:', sthp, cthp, acos(cthp)*180./PI       
c     realize exchange?
                     xf = 2.D0*p(i1,3)/SQS
                     xff = abs(xf)**xfexp
                     IF(S_RNDM(I1).lt.xff.and.S_RNDM(I1).lt.xrate)then
                        icount = icount + 1
c                        print*,' replacing with baryon pair'
                        ineutron = int(0.5D0+S_RNDM(I2))
c                        print*,' ineutron:', ineutron
                        if(xmass2.le.4*(pmin**2+am2(14))) ineutron = 0 ! neutron not possible
                        iflip = -1 + 2*int(0.5D0 + S_RNDM(I3))
                        id1 = (13+ineutron)*iflip
                        id2 = (-13-ineutron)*iflip
                        if(ndebug.gt.0)then
                           write(lun,*) 'force baryons:'
                           if(ipart(icomb,3).eq.0)then
                              write(lun,*) 'replacing indices:' ,i1, i2
                              write(lun,*) 'ids:', llist(i1), llist(i2)
                           else
                              write(lun,*) 'replacing indices:',i1,i2,i3
                              write(lun,*) 'ids:', llist(i1), llist(i2),
     &                             llist(i3)
                           endif
                           write(lun,*)' with baryon pair!'
                        endif
c                        print*,' new ids: ', id1, id2
c     determine momenta of new particles and write to stack
                        enew = xm/2.D0
                        pnew = sqrt(xmass2/4.D0-am2(abs(id1)))
                        ph = S_RNDM(I2)*TWOPI
c     diq pt
 80                     call ptdis_4flv(12,pn(1),pn(2))
                        call ptdis_4flv(1,px,py)
                        pn(1) = 1.8*pn(1) + px
                        pn(2) = 1.8*pn(2) + py
                        pts = pn(1)**2+pn(2)**2
                        pz2 = pnew**2-pts
                        if(pz2.lt.0.D0) goto 80                        
                        pn(1) = cos(ph) * sqrt(pts)
                        pn(2) = sin(ph) * sqrt(pts)
                        pn(3) = sqrt(pz2)
                        pn(4) = enew
c                        print*,'enew,pnew:', enew,pnew
c                        print*,'p1n',pn
c                        print*,'ptn',pnew*sthp,sqrt(pn(1)**2+pn(2)**2)
c                        print*,'m2',enew**2-pnew**2
c                        print*,'m2 vec:',pn(4)**2-pn(1)**2-pn(2)**2
c     &                       -pn(3)**2
c                        print*,'gamma:', gamma
                        p1new(5) = am(abs(id1))
                        p2new(5) = am(abs(id2))
                        call sib_altra(gamma,
     &                       gambet(1),gambet(2),gambet(3),
     &                       pn(1),pn(2),pn(3),pn(4),
     &                       p1newt,p1new(1),p1new(2),p1new(3),p1new(4))
                        call sib_altra(gamma,
     &                       gambet(1),gambet(2),gambet(3),
     &                       -pn(1),-pn(2),-pn(3),pn(4),
     &                       p2newt,p2new(1),p2new(2),p2new(3),p2new(4))          
c                        print*,'p1new:' , p1new
c                        print*,'p2new:' , p2new
c                        call four_length(p1new,xm)
c                        print*,'p1mass_2:', xm, am2(abs(id1))
c                        sqrt(p1new(4)**2-p1new(1)**2
c     &                       -p1new(2)**2-p1new(3)**2)
c                        call four_length(p2new,xm)
c                        print*,'p2mass_2:', xm, am2(abs(id2))

c                        print*,'pt1 :', sqrt(p1(1)**2+p1(2)**2)
c                        print*,'pt2 :', sqrt(p2(1)**2+p2(2)**2)
                        
c                        print*,'pt1 new:', sqrt(p1new(1)**2+p1new(2)**2)
c                        print*,'pt2 new:', sqrt(p2new(1)**2+p2new(2)**2)
c                        if(sqrt(p2new(1)**2+p2new(2)**2).gt.80)stop
c                        if(sqrt(p1new(1)**2+p1new(2)**2).gt.80)stop
c     write to stack
                        llist(i1) = id1
                        llist(i2) = id2
                        do kk=1,5
                           p(i1,kk) = p1new(kk)
                           p(i2,kk) = p2new(kk)
                        enddo
c     3-->2 exchange! reduce stack by moving last particle to position 3
                        if(ipart(icomb,3).ne.0)then
                           llist(i3) = llist(NP)
                           llist1(i3) = llist1(NP)
                           nporig(i3) = nporig(NP)
                           nforig(i3) = nforig(NP)
                           do kk=1,5
                              p(i3,kk) = p(NP,kk)
                           enddo                     
                           NP = NP - 1
                        endif
                     endif      ! end realize exchange?
                     goto 100   ! next index I1
 400                 continue
                  enddo         ! end I3
 300              continue
               enddo            ! end I2
            ENDIF               ! check I1?
c 200        continue            
         ENDDO                  ! end combinations
 100     continue
      ENDDO ! end i1
c      print*,icount, 'replacements done'
      if(ndebug.ge.5) CALL SIB_LIST(6)
      END
C=======================================================================

      SUBROUTINE FORCE_STRANGE(XRATE1,IEDEP,ETHR,XFEXP,IPROJ,N1,N2)

C-----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
C     The final particle output is contained in COMMON /S_PLIST/    
C     NP           : number of final particles
C     P(1:NP, 1:5) : 4-momenta + masses of the final particles 
C     LLIST (1:NP) : codes of final particles
      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP
      INTEGER NW_max
      PARAMETER (NW_max = 20)
      INTEGER LLIST1
      COMMON /S_PLIST1/ LLIST1(8000)
C--------------------------------------------------------------------
C     SIBYLL common blocks containing event information       \FR'14
C--------------------------------------------------------------------

C     EVENT INFO COMMON
C     contains overall interaction properties, like
C     SQS : center-of-mass energy
C     S   :         "       "     squared
C     PTmin : low pt cut of QCD cross section, 
C             i.e. minimal pt of hard minijets
C     Xmin : low-x bound for PDFs, 
C            i.e. minimal momentum fraction of hard partons
C     Zmin : logarithm of that
C     KB : PID of beam hadron
C     KT() : PID of target
C     IAT : mass number of target
      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
      INTEGER KB,IAT,KT
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)

      INTEGER NFORIG,NPORIG,NIORIG,IPFLAG,IIFLAG,KINT
      COMMON /S_PARTO/ NFORIG(NP_max),NPORIG(NP_max),NIORIG(NP_max),
     &IPFLAG,IIFLAG,KINT

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN

c     external types
      double precision xrate1, ethr, xfexp
      integer n1,n2,iedep,iproj

c     internal types
      integer ncomb
      parameter (ncomb=7)
      integer kk,i1,i2,i3,kba,ineutral,incrmax,incr1,incr2,id1,id2,
     &     iflip,icount,ipart,icomb,ismeson
      DIMENSION ipart(ncomb,3),ismeson(99)
      double precision S_RNDM,xrate,ecm_max,pmin,xmass2,pts
     &     ,ptot,ptmp,p1,p2,p3,gamma,gambet,xm,enew,pnew,ph,pn,p1new,
     &     p2new,xf,xff,p1newt,p2newt,px,py,pz2
c     &     ,sthp,cthp,p1p,p1pt,ptp
      
      DIMENSION ptot(5),ptmp(5),p1(5),p2(5),
     &     p3(5),gambet(4),p1new(5),p2new(5),pn(4)!,p1p(4)
      EQUIVALENCE (gambet(4),gamma)
      SAVE
      DATA IPART /6,7,8,6,6,7,8, 6,8,7,6,7,8,7, 0,0,0,6,8,6,6/
      DATA ISMESON /5*0,7*1,87*0/
      DATA ECM_MAX /1.36987D+05/ ! 10^19eV in ecm and gev for protons
      DATA PMIN /0.05D0/          ! minimal momentum in GeV
      
      KBA = IABS(KB)
      
c     energy dependence
      if(iedep.eq.0) then
         xrate = xrate1
      elseif(iedep.eq.1)then
         xrate=xrate1 * max(0.D0,log(SQS/ethr)/log(ecm_max/ethr))
      endif
      if(IPROJ.eq.1) then
c     meson projectiles only
         if(ismeson(kba).ne.1) return         
      endif      

c      CALL SIB_LIST(6)
c      print*,'xrate', xrate
      icount = 0
c     check if candidate (pi0, pi+, pi-)      
c      icomb = 1
c     loop over particle stack
      DO I1=N1,N2
c         print*,'checking index:',I1
         DO icomb=1,ncomb
c            print*,'**** checking index:',I1
c            print*,'combination:',icomb
            IF(llist(I1).eq.ipart(icomb,1)) then
c     next neighbours         
               incrmax = max(0,min(5,MIN(N2,NP)-I1))
c               print*,'max increment:',incrmax
               do incr1=1,incrmax
                  I2 = I1 + incr1
                  if(llist(I2).ne.ipart(icomb,2)) goto 300 ! next particle2
                  if(ipart(icomb,3).eq.0)then
c     two-particle combination
c                     print*,' 2 particle combination found!'
c                     print*,' indices:' ,i1, i2
c                     print*,' ids:', llist(i1), llist(i2)
c     sufficient mass?
                     do kk=1,5
                        p1(kk) = p(i1,kk)
                        p2(kk) = p(i2,kk)
                     enddo
                     call ADD_4VECS(P1,P2,ptot)
                     call four_length(ptot,xmass2)
                     if(xmass2.le.4*(pmin**2+am2(9))) then
c                       print*,'insufficient mass!'
                        goto 300 ! next particle2
                     else
                        goto 500 ! implement exchange
                     endif
                     
                  endif
                  do incr2=1,incrmax
                     if(incr1.eq.incr2) goto 400 ! next particle3
                     I3 = I1 + incr2
c                     print*,'checking..',i1,i2,i3
                     if(llist(I3).ne.ipart(icomb,3)) goto 400 ! next particle3
c     pi0,pi+,pi- combination found!
c                     print*,' combination found!'
c                     print*,' indices:' ,i1, i2, i3
c                     print*,' ids:', llist(i1), llist(i2), llist(i3)
c     sufficient mass?
                     do kk=1,5
                        p1(kk) = p(i1,kk)
                        p2(kk) = p(i2,kk)
                        p3(kk) = p(i3,kk)
                     enddo
                     call ADD_4VECS(P1,P2,ptmp)
                     call ADD_4VECS(ptmp,P3,ptot)
                     call four_length(ptot,xmass2)
                     if(xmass2.le.4*(pmin**2+am2(13))) then
c                        print*,'insufficient mass!'
                        goto 400 ! next particle3
                     endif
 500                 continue                     
                     xm = sqrt(xmass2)                    
c     calc. initial momenta in 2-particle center-of-mass (eg. p1'--> p1p)
c     use these to calculate the angles which then define pt after the exchange
                     do kk=1,4
                        gambet(kk) = ptot(kk)/xm
                     enddo
CC                     call sib_altra(gambet(4),-gambet(1),-gambet(2)
CC     &                    ,-gambet(3),p1(1),p1(2),p1(3),p1(4),
CC     &                    p1pt,p1p(1),p1p(2),p1p(3),p1p(4))
CC                     ptp = sqrt(p1p(1)**2+p1p(2)**2)
CC                     sthp = ptp/p1pt ! sin(theta')
CC                     cthp = p1p(3)/p1pt ! cos(theta')                     
c     realize exchange?
                     xf = 2.D0*p(i1,3)/SQS
                     xff = abs(xf)**xfexp
                     IF(S_RNDM(I1).lt.xff.and.S_RNDM(I1).lt.xrate)then
                        icount = icount + 1
c                        print*,' replacing with new pair'
                        ineutral = int(0.5D0+S_RNDM(I2)) ! 0 or 1
c                        print*,' ineutral:', ineutral
                        if(xmass2.le.4*(pmin**2+am2(21))) ineutral = 0 ! neutron not possible
                        iflip = int(0.5D0+S_RNDM(I3)) ! 0 or 1
                        if(0.5D0.lt.S_RNDM(I1).and.
     &                       xmass2.gt.4*(pmin**2+am2(28)))then
                           if(xmass2.le.4*(pmin**2+am2(30)))
     &                          ineutral = 0 ! neutron not possible
                           id1 = 28+iflip + 2*ineutral ! k*+:28 k*-: 29 k*0: 30 k*0bar: 31
                           id2 = 29-iflip + 2*ineutral
                        else
                           id1 = 9+iflip + 12*ineutral ! k+:9 k-: 10 k0: 21 k0bar: 22
                           id2 = 10-iflip + 12*ineutral                           
                        endif                           
c                        print*,' new ids: ', id1, id2
c     determine momenta of new particles and write to stack
                        enew = xm/2.D0
                        pnew = sqrt(xmass2/4.D0-am2(abs(id1)))                        
                        ph = S_RNDM(I2)*TWOPI
c     diq pt
 80                     call ptdis_4flv(3,pn(1),pn(2))
                        call ptdis_4flv(1,px,py)
                        pn(1) = 1.8*pn(1) + px
                        pn(2) = 1.8*pn(2) + py
                        pts = pn(1)**2+pn(2)**2
                        pz2 = pnew**2-pts
                        if(pz2.lt.0.D0) goto 80                        
                        pn(1) = cos(ph) * sqrt(pts)
                        pn(2) = sin(ph) * sqrt(pts)
                        pn(3) = sqrt(pz2)
                        pn(4) = enew                        
                       call sib_altra(gamma,
     &                       gambet(1),gambet(2),gambet(3),
     &                       pn(1),pn(2),pn(3),pn(4),
     &                       p1newt,p1new(1),p1new(2),p1new(3),p1new(4))
                        call sib_altra(gamma,
     &                       gambet(1),gambet(2),gambet(3),
     &                       -pn(1),-pn(2),-pn(3),pn(4),
     &                       p2newt,p2new(1),p2new(2),p2new(3),p2new(4))          
                        p1new(5) = am(abs(id1))
                        p2new(5) = am(abs(id2))
c     write to stack
                        llist(i1) = id1
                        llist(i2) = id2
                        do kk=1,5
                           p(i1,kk) = p1new(kk)
                           p(i2,kk) = p2new(kk)
                        enddo
c     3-->2 exchange! reduce stack by moving last particle to position 3
                        if(ipart(icomb,3).ne.0)then
                           llist(i3) = llist(NP)
                           llist1(i3) = llist1(NP)
                           nporig(i3) = nporig(NP)
                           nforig(i3) = nforig(NP)
                           do kk=1,5
                              p(i3,kk) = p(NP,kk)
                           enddo                     
                           NP = NP - 1
                        endif
                     endif      ! end realize exchange?
                     goto 100   ! next index I1
 400                 continue
                  enddo         ! end I3
 300              continue
               enddo            ! end I2
            ENDIF               ! check I1?
c 200        continue            
         ENDDO                  ! end combinations
 100     continue
      ENDDO ! end i1
c      print*,icount, 'replacements done'
      if(ndebug.ge.5) CALL SIB_LIST(6)
      END
      
C=======================================================================

      SUBROUTINE FORCE_VECTORS1(XRATE1,IEDEP,ETHR,XFEXP,IPROJ,N1,N2)

C-----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
C     The final particle output is contained in COMMON /S_PLIST/    
C     NP           : number of final particles
C     P(1:NP, 1:5) : 4-momenta + masses of the final particles 
C     LLIST (1:NP) : codes of final particles
      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP
      INTEGER NW_max
      PARAMETER (NW_max = 20)
      INTEGER LLIST1
      COMMON /S_PLIST1/ LLIST1(8000)
C--------------------------------------------------------------------
C     SIBYLL common blocks containing event information       \FR'14
C--------------------------------------------------------------------

C     EVENT INFO COMMON
C     contains overall interaction properties, like
C     SQS : center-of-mass energy
C     S   :         "       "     squared
C     PTmin : low pt cut of QCD cross section, 
C             i.e. minimal pt of hard minijets
C     Xmin : low-x bound for PDFs, 
C            i.e. minimal momentum fraction of hard partons
C     Zmin : logarithm of that
C     KB : PID of beam hadron
C     KT() : PID of target
C     IAT : mass number of target
      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
      INTEGER KB,IAT,KT
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)

      INTEGER NFORIG,NPORIG,NIORIG,IPFLAG,IIFLAG,KINT
      COMMON /S_PARTO/ NFORIG(NP_max),NPORIG(NP_max),NIORIG(NP_max),
     &IPFLAG,IIFLAG,KINT

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN

c     external types
      double precision xrate1, ethr, xfexp
      integer n1,n2,iedep,iproj

c     internal types
      integer kk,i1,i2,kba,incrmax,incr1,id1,id2,
     &     icount,ismeson,aid2
      
      DIMENSION ismeson(99)
      double precision S_RNDM,xrate,ecm_max,pmin,xmass2,xmin2,
     &     ptot,p1,p2,gamma,gambet,xm,pnew,ph,pn,p1new,
     &     p2new,pawt,enew1,enew2,xf,xff,p2newt,p1newt,px,py,pts,pz2
c     &     ,ptp,p1p,p1pt,sthp, cthp
      
      DIMENSION ptot(5),p1(5),p2(5),
     &     gambet(4),p1new(5),p2new(5),pn(4)!,p1p(4)
      EQUIVALENCE (gambet(4),gamma)
      SAVE
      DATA ISMESON /5*0,7*1,87*0/
      DATA ECM_MAX /1.36987D+05/ ! 10^19eV in ecm and gev for protons
      DATA PMIN /0.05D0/          ! minimal momentum in GeV
      
      KBA = IABS(KB)
      
c     energy dependence
      if(iedep.eq.0) then
         xrate = xrate1
      elseif(iedep.eq.1)then
         xrate=xrate1 * max(0.D0,log(SQS/ethr)/log(ecm_max/ethr))
      endif
      if(IPROJ.eq.1) then
c     meson projectiles only
         if(ismeson(kba).ne.1) return         
      endif      

c      CALL SIB_LIST(6)
c      print*,'xrate', xrate
      icount = 0
c     loop over particle stack
      DO I1=N1,N2
c         print*,'**** checking index:',I1
         IF(llist(I1).eq.6) then
c     next neighbours         
            incrmax = max(0,min(5,MIN(N2,NP)-I1))
c            print*,'max increment:',incrmax
            do incr1=1,incrmax
               I2 = I1 + incr1
               if(abs(llist(i2)).gt.10000) goto 300
c     two-particle combination
c               print*,' checking two particle combination..'
c               print*,' indices:' ,i1, i2
c               print*,' ids:', llist(i1), llist(i2)
c     sufficient mass?
               do kk=1,5
                  p1(kk) = p(i1,kk)
                  p2(kk) = p(i2,kk)
               enddo
               call ADD_4VECS(P1,P2,ptot)
               call four_length(ptot,xmass2)
               aid2 = mod(abs(llist(i2)),10000)
               xmin2 = am2(27)+am2(aid2)
               xmin2 = xmin2+2*(pmin**2+am(27)*am(aid2))
               if(xmass2.le.xmin2) then
c                  print*,'insufficient mass!', xmass2, xmin2
                  goto 300      ! next particle2
               endif
               xm = sqrt(xmass2)
c     print*,' sufficient mass! mass', xm, xmin2
               do kk=1,4
                  gambet(kk) = ptot(kk)/xm
               enddo
CC               call sib_altra(gambet(4),-gambet(1),-gambet(2)
CC     &              ,-gambet(3),p1(1),p1(2),p1(3),p1(4),
CC     &              p1pt,p1p(1),p1p(2),p1p(3),p1p(4))
CC               ptp = sqrt(p1p(1)**2+p1p(2)**2)
CC               sthp = ptp/p1pt  ! sin(theta')
CC               cthp = p1p(3)/p1pt ! cos(theta')
c     realize exchange?
               xf = 2.D0*p(i1,3)/SQS
               xff = abs(xf)**xfexp
               IF(S_RNDM(I1).lt.xff.and.S_RNDM(I1).lt.xrate)then
                  icount = icount + 1
                  id1 = 27
                  id2 = llist(i2)
                  aid2 = mod(abs(id2),10000)
                  if(ndebug.gt.0)then
                     print*,' indices:' ,i1, i2
                     print*,' new ids: ', id1, id2
                  endif
c     determine momenta of new particles and write to stack
                  do kk=1,4
                     gambet(kk) = ptot(kk)/xm
                  enddo
                  enew1 = (xmass2-am2(aid2)+am2(id1))/(2.D0*xm)
                  enew2 = (xmass2+am2(aid2)-am2(id1))/(2.D0*xm)
                  pnew = pawt(xm,am(aid2),am(id1)) 
                  ph = S_RNDM(I2)*TWOPI
CCc     use angles from initial state
CC                  pn(1) = pnew * sthp * cos(ph)
CC                  pn(2) = pnew * sthp * sin(ph)
CC                  pn(3) = pnew * cthp
 80               call ptdis_4flv(2,pn(1),pn(2))
                  call ptdis_4flv(1,px,py)
                  pn(1) = 1.*pn(1) + px
                  pn(2) = 1.*pn(2) + py
                  pts = pn(1)**2+pn(2)**2
                  pz2 = pnew**2-pts
                  if(pz2.lt.0.D0) goto 80                        
                  pn(1) = cos(ph) * sqrt(pts)
                  pn(2) = sin(ph) * sqrt(pts)
                  pn(3) = sqrt(pz2)
                  call sib_altra(gamma,
     &                 gambet(1),gambet(2),gambet(3),
     &                 pn(1),pn(2),pn(3),enew1,
     &                 p1newt,p1new(1),p1new(2),p1new(3),p1new(4))
                  call sib_altra(gamma,
     &                 gambet(1),gambet(2),gambet(3),
     &                 -pn(1),-pn(2),-pn(3),enew2,
     &                 p2newt,p2new(1),p2new(2),p2new(3),p2new(4))          
                  p1new(5) = am(abs(id1))
                  p2new(5) = am(aid2)
c     write to stack
                  llist(i1) = id1
                  llist(i2) = id2
                  do kk=1,5
                     p(i1,kk) = p1new(kk)
                     p(i2,kk) = p2new(kk)
                  enddo
               endif            ! end realize exchange?
               goto 100         ! next index I1
 300           continue
            enddo               ! end I2
         ENDIF                  ! check I1?
 100     continue
      ENDDO ! end i1
c      print*,icount, 'replacements done'
      if(ndebug.ge.5) CALL SIB_LIST(6)
      END
      
C=======================================================================
      
      SUBROUTINE FORCE_VECTORS(XRATE1,IEDEP,ETHR,XFEXP,N1,N2)

C-----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
C     The final particle output is contained in COMMON /S_PLIST/    
C     NP           : number of final particles
C     P(1:NP, 1:5) : 4-momenta + masses of the final particles 
C     LLIST (1:NP) : codes of final particles
      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP
      INTEGER NW_max
      PARAMETER (NW_max = 20)
C--------------------------------------------------------------------
C     SIBYLL common blocks containing event information       \FR'14
C--------------------------------------------------------------------

C     EVENT INFO COMMON
C     contains overall interaction properties, like
C     SQS : center-of-mass energy
C     S   :         "       "     squared
C     PTmin : low pt cut of QCD cross section, 
C             i.e. minimal pt of hard minijets
C     Xmin : low-x bound for PDFs, 
C            i.e. minimal momentum fraction of hard partons
C     Zmin : logarithm of that
C     KB : PID of beam hadron
C     KT() : PID of target
C     IAT : mass number of target
      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
      INTEGER KB,IAT,KT
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)

      INTEGER NFORIG,NPORIG,NIORIG,IPFLAG,IIFLAG,KINT
      COMMON /S_PARTO/ NFORIG(NP_max),NPORIG(NP_max),NIORIG(NP_max),
     &IPFLAG,IIFLAG,KINT

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN

c     external types
      double precision xrate1, ethr, xfexp
      integer n1,n2,iedep

c     internal types
      integer ipi2vec,ismeson,lcon,lreschex,ll,la,la_new,i,j,kba
      DIMENSION IPI2VEC(99),ISMESON(99)
      double precision pz2,xmts,xf,xfs,S_RNDM,xff,xrate,ecm_max !,pts
      
      DIMENSION LCON(6:43),LRESCHEX(6:39)
      INTEGER IFIRST
      SAVE
c     charge exchange map, i.e. pip -> pi0 ...
      DATA LCON /7,6,6,22,21,9,9,14,13,4*0,20,19,10,9,23,24,27,27,25,
     &     31,30,29,28,32,33,35,34,35,38,37,39,41,42,41,42/
c     charge and spin exchange map, i.e. pip -> rho0
c     approximate, proton and neutron should go to N(1520) not Delta
      DATA LRESCHEX /26,27,27,31,30,9,9,42,41,19*0,45,44,45,48,47,39/ 
      DATA IFIRST /0/
      DATA ECM_MAX /1.36987D+05/  ! 10^19eV in ecm and gev for protons

      if(ifirst.eq.0)then
         print *,'initializing..'
         do j=1,99
            IPI2VEC(J) = J
            ismeson(j) = 0
         enddo
         IPI2VEC(6) = 27        ! pi(0) ---> rho(0)
         IPI2VEC(7) = 25        ! pi+   ---> rho+
         IPI2VEC(8) = 26        ! pi-   ---> rho-

         ISMESON(7) = 1         ! pi+
         ISMESON(8) = 1         ! pi-
         ISMESON(9) = 1         ! k+
         ISMESON(10) = 1        ! k-
         ISMESON(11) = 1        ! k0s
         ISMESON(12) = 1        ! k0l

         ifirst = 1
      endif

      KBA = IABS(KB)
      
c     energy dependence
      if(iedep.eq.0) then
         xrate = xrate1
      elseif(iedep.eq.1)then
         xrate=xrate1 * max(0.D0,log(SQS/ethr)/log(ecm_max/ethr))
      endif
      if(IPAR(97).eq.1) then
c     meson projectiles only
         if(ismeson(kba).ne.1) return         
      endif           
      
      IF(IPAR(45).eq.1)THEN
c     trivial exchange model      
         do I=N1,N2
c     replace pions with vector mesons
            LL = mod(llist(I),10000)
            LA = abs(LL)
            IF(S_RNDM(I).lt.xrate)then
c     put back on mass shell
               la_new = IPI2VEC(LA)
               xmts = p(i,1)**2 + p(i,2)**2 + am2(la_new)
               pz2 = p(i,4)**2 - xmts
               if(pz2.gt.EPS8)then
                  p(i,3) = sign(sqrt(pz2),p(i,3))
                  p(i,5) = am(la_new)
                  LLIST(I) = ISIGN(la_new,ll)
               endif
            endif
         enddo

      ELSEIF(IPAR(45).eq.2)THEN
c     large xf only, neutral pions only
         do I=N1,N2
            LL = mod(llist(I),10000)
            LA = abs(LL)
            IF(LA.eq.6)then
               xf = 2.D0*p(i,3)/SQS
               IF(S_RNDM(I).lt.xrate*xf)then
c     exhcange and put back on mass shell
                  la_new = IPI2VEC(la)
                  xmts = p(i,1)**2 + p(i,2)**2 + am2(la_new)
                  pz2 = p(i,4)**2 - xmts
                  if(pz2.gt.EPS8)then
                     p(i,3) = sign(dsqrt(pz2),p(i,3))
                     p(i,5) = am(la_new)
                     LLIST(I) = ISIGN(la_new,ll)
                  endif
               endif
            endif
         enddo

      ELSEIF(IPAR(45).eq.3)THEN
c     large xf only, charge and spin exchange
         do I=N1,N2
            LL = mod(llist(I),10000)
            LA = abs(LL)
            IF(ll.eq.LCON(KBA))then
               xf = 2.D0*p(i,3)/sqs
               IF(S_RNDM(I).lt.xrate*xf)then
c     replace charge exchange product of beam with
c     charge and spin exchange product, i.e.
c     pip-beam -> rho0 instead of pip-beam -> pi0
c     so replace pi0 with rho0 in final state
                  la_new = LRESCHEX(KBA)
c     put back on mass shell
                  xmts = p(i,1)**2 + p(i,2)**2 + am2(la_new)
                  pz2 = p(i,4)**2 - xmts
                  if(pz2.gt.EPS8)then
                     p(i,3) = sign(dsqrt(pz2),p(i,3))
                     p(i,5) = am(la_new)
                     LLIST(I) = ISIGN(la_new,ll)
                  endif
               endif
            endif
         enddo
         
      ELSEIF(IPAR(45).eq.4)THEN
c     large xf only, charge and spin exchange
         do I=N1,N2
            LL = mod(llist(I),10000)
            LA = abs(ll)
            IF(LL.eq.LCON(KBA))then
               xf = 2.D0*p(i,3)/sqs
               xfs = xf ** 2
               IF(S_RNDM(I).lt.xrate*xfs)then
c     replace charge exchange product of beam with
c     charge and spin exchange product, i.e.
c     pip-beam -> rho0 instead of pip-beam -> pi0
c     so replace pi0 with rho0 in final state
                  la_new = LRESCHEX(KBA)
c     put back on mass shell
                  xmts = p(i,1)**2 + p(i,2)**2 + am2(la_new)
                  pz2 = p(i,4)**2 - xmts
                  if(pz2.gt.EPS8)then
                     p(i,3) = sign(dsqrt(pz2),p(i,3))
                     p(i,5) = am(la_new)
                     LLIST(I) = ISIGN(la_new,ll)
                  endif
               endif
            endif
         enddo
         
      ELSEIF(IPAR(45).eq.5)THEN
c     trivial exchange model, neutral pions only      
         do I=N1,N2
c     replace pions with vector mesons
            LL = mod(llist(I),10000)
            LA = abs(LL)
            IF(LA.eq.6.and.S_RNDM(I).lt.xrate)then
c     put back on mass shell
               la_new = IPI2VEC(LA)
               xmts = p(i,1)**2 + p(i,2)**2 + am2(la_new)
               pz2 = p(i,4)**2 - xmts
               if(pz2.gt.EPS8)then
                  p(i,3) = sign(sqrt(pz2),p(i,3))
                  p(i,5) = am(la_new)
                  LLIST(I) = ISIGN(la_new,ll)
               endif
            endif
         enddo

      ELSEIF(IPAR(45).eq.6)THEN
c     trivial exchange model, neutral pions only, weighted by xF
         do I=N1,N2
c     replace pions with vector mesons
            LL = mod(llist(I),10000)
            LA = abs(LL)
            IF(LA.eq.6)then
               xf = 2.D0*p(i,3)/SQS
               xff = abs(xf)**xfexp
               IF(S_RNDM(I).lt.xff.and.S_RNDM(LA).lt.xrate)then
c     put back on mass shell
                  la_new = IPI2VEC(LA)
                  xmts = p(i,1)**2 + p(i,2)**2 + am2(la_new)
                  pz2 = p(i,4)**2 - xmts
                  if(pz2.gt.EPS8)then
                     p(i,3) = sign(sqrt(pz2),p(i,3))
                     p(i,5) = am(la_new)
                     LLIST(I) = ISIGN(la_new,ll)
                  endif
               endif
            endif
         enddo
         
      ENDIF
      if(ndebug.ge.5) CALL SIB_LIST(6)
      END
C=======================================================================

      SUBROUTINE SAMPLE_BEAM(KID,NW,XCHG,KRMNT,XJET,IREJ)

C-----------------------------------------------------------------------
C...Subroutine to sample valence and sea quark kinematics
C.    fills IFL?,X? and PX?,PY?
C.    1,2 are valence quarks, 3,4 are additional sea quarks
C.    transverse momentum is shared between the val. and sea pairs
C.    X and flv are exchanged occasionally
C-------------------------------------------------------------------      
      IMPLICIT NONE

      DOUBLE PRECISION XCHG,XJET
      INTEGER KID,NW,KRMNT,IREJ

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      INTEGER NW_max
      PARAMETER (NW_max = 20)
      INTEGER NS_max, NH_max
      PARAMETER (NS_max = 20, NH_max = 80)
C--------------------------------------------------------------------
C     SIBYLL common blocks containing event information       \FR'14
C--------------------------------------------------------------------

C     EVENT INFO COMMON
C     contains overall interaction properties, like
C     SQS : center-of-mass energy
C     S   :         "       "     squared
C     PTmin : low pt cut of QCD cross section, 
C             i.e. minimal pt of hard minijets
C     Xmin : low-x bound for PDFs, 
C            i.e. minimal momentum fraction of hard partons
C     Zmin : logarithm of that
C     KB : PID of beam hadron
C     KT() : PID of target
C     IAT : mass number of target
      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
      INTEGER KB,IAT,KT
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT

      INTEGER IBMRDX,ITGRDX,IHMJDX,ISMJDX,ICSTDX,IINTDX
      COMMON /S_INDX/ IBMRDX(3),ITGRDX(NW_max,3),
     &     IHMJDX(NW_max*NH_max),IINTDX(NW_max),
     &     ISMJDX(NW_max*NS_max),ICSTDX(2*NW_max,3)

      INTEGER IRMNT,KRB,KRT
      DOUBLE PRECISION XRMASS,XRMEX
      COMMON /S_RMNT/ XRMASS(2),XRMEX(2),IRMNT(NW_max),KRB,KRT(NW_max)

      DOUBLE PRECISION X1,PXB,PYB
      DIMENSION X1(2*NW_max+2),PXB(2*NW_max+2),PYB(2*NW_max+2)
      INTEGER IFLB,KID1,J,J1,J2,J3,J4,Iref1,Iref,Idm
      DIMENSION IFLB(2*NW_max+2)

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      SAVE

c     default rejection
c     options are: 1: resample minijets (Xjet)..
c                  2: resample non-diff event (Ns,Nh)..
c                  3: resample event (Nw,diff,ndiff)..
      IREJ = 1

      IF(ndebug.gt.2) 
     +     WRITE(LUN,*)
     +     ' SAMPLE_BEAM: KID,NW,XCHG,KRMNT,XJET,IREJ',
     +     KID,NW,XCHG,KRMNT,XJET,IREJ

      CALL SAMPLE_PROJECTILE
     +     (KID,NW,KRMNT,XCHG,XJET,X1,PXB,PYB,IFLB,KID1,IREJ)    
      IF(IREJ.ne.0) RETURN

c     set remnant id to beam
c     will be changed if flavor is exchanged between central strings and remnant
      KRB = KID1        

C..   write beam partons to stack
c     order is: val1, val2, q, qbar etc
      IF(KRMNT.ne.0)THEN
         j1 = 1
         j2 = 2
c     add proto-remnant (still massless)
         CALL ADD_PRTN(PXB(J1)+PXB(J2),PYB(J1)+PYB(J2),
     &        0.5D0*SQS*(X1(J1)+X1(J2)),
     &        0.5D0*SQS*(X1(J1)+X1(J2)),0.D0,2,0,0,Iref1)
         IBMRDX(1) = Iref1
c     beam remnant always associated with first interaction
         CALL ADD_INT_REF(Iref1,IINTDX(1))
c     add quarks designated for remnant
         IF(KID.lt.0)THEN
c     if beam is antibaryon then hspli puts diq into 1st flv
c     need to switch to fit call to string frag routine 
c     such that diq is along +z
            CALL ISWTCH_LMNTS(IFLB(j1),IFLB(j2))
         ENDIF
         CALL ADD_PRTN(PXB(J1),PYB(J1),0.5D0*SQS*X1(J1),
     &        0.5D0*SQS*X1(J1),0.D0,IFLB(J1),1,Iref1,Iref)
         IBMRDX(2) = Iref
         CALL ADD_PRTN(PXB(J2),PYB(J2),0.5D0*SQS*X1(J2),
     &        0.5D0*SQS*X1(J2),0.D0,IFLB(J2),1,Idm,Iref)
         IBMRDX(3) = Iref
      ENDIF
      DO j=1,NW
         j3 = 3+(j-1)*2
         j4 = j3+1
c     add sea quarks
         CALL ADD_PRTN(PXB(J3),PYB(J3),0.5D0*SQS*X1(J3),
     &        0.5D0*SQS*X1(J3),0.D0,IFLB(J3),1,0,Iref)
         ICSTDX(2*(J-1)+1,2) = Iref
         CALL ADD_PRTN(PXB(J4),PYB(J4),0.5D0*SQS*X1(J4),
     &        0.5D0*SQS*X1(J4),0.D0,IFLB(J4),1,0,Iref)
         ICSTDX(2*(J-1)+2,2) = Iref
c     add parton index to cache
      ENDDO
      IF(NDEBUG.GT.3) CALL PRNT_PRTN_STCK

      IREJ = 0

      END
C=======================================================================

      SUBROUTINE FRAG_INCHRNT_DIFF(IDX,LBAD)

C-----------------------------------------------------------------------
C     routine that fragments a diffractive system               \FR'15
C
C     INPUT: IDX : parton stack index of 4momentum
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER IDX,LBAD

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
C     The final particle output is contained in COMMON /S_PLIST/    
C     NP           : number of final particles
C     P(1:NP, 1:5) : 4-momenta + masses of the final particles 
C     LLIST (1:NP) : codes of final particles
      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP

      INTEGER NFORIG,NPORIG,NIORIG,IPFLAG,IIFLAG,KINT
      COMMON /S_PARTO/ NFORIG(NP_max),NPORIG(NP_max),NIORIG(NP_max),
     &IPFLAG,IIFLAG,KINT

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN

      DOUBLE PRECISION PST,PDIFF,GABE,P2,EE,P1TOT
      DIMENSION PST(5),PDIFF(5),GABE(4),P2(4)
      INTEGER IDIFF1,IDIFF,IPID,L0,JDIFF,NOLD,LXBAD,K,II
      SAVE

      LBAD = 2

c     references are diff --> diff.hadron --> bm-partons --> tg-partons
c     only diff and diff. hadron are read out
c     read diff 4momentum from stack
      CALL RD_PRTN_4VEC(IDX,PST,IPID,IDIFF1)
      CALL RD_PRTN_4VEC(IDIFF1,PDIFF,L0,IDIFF)
      
C     kinematic variables
      EE = PDIFF(5)             ! center of mass energy in diff. system
      
c     set diffraction code of system (1:beam,2:target,3:double)
      JDIFF = ABS(IPID)/10

      IF(NDEBUG.gt.1) WRITE(LUN,*)' FRAG_INCHRNT_DIFF: IDX,EE,L0',
     &     IDX,EE,L0

      IPFLAG = -1

      NOLD = NP

c     diffractive interaction in center-of-mass system of (sea,rmnt)-nuc
      CALL SIB_DIFF(L0,JDIFF,EE,0,LXBAD)
      IF(LXBAD.ne.0) THEN
         IF(NDEBUG.gt.1) 
     &        WRITE(LUN,*)' FRAG_INCHRNT_DIFF: fragmentation rejection' 
         RETURN
      ENDIF
      IF(NDEBUG.gt.1) 
     &     WRITE(LUN,*)' FRAG_INCHRNT_DIFF: particles before/after :',
     &     NOLD,NP

c     boost to hadron - hadron center-of-mass
      do ii=1,4
         gabe(ii) = PDIFF(ii)/PDIFF(5)
      enddo
      DO K=NOLD+1,NP
         CALL SIB_ALTRA(gabe(4),gabe(1),gabe(2),
     &        gabe(3),P(k,1),p(k,2),p(k,3),p(k,4),
     &        P1TOT,p2(1),p2(2),p2(3),p2(4))
         do ii=1,4
            P(K,ii)=P2(ii)
         enddo
      ENDDO

      LBAD = 0
      END
C=======================================================================

      SUBROUTINE SAMPLE_MINIJET
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
      IMPLICIT NONE
      
c     external types
      INTEGER L,NW,NNJET,NNSOF,NJET,NSOF,LBAD
      DOUBLE PRECISION X1JET,X2JET
      
      INTEGER NS_max, NH_max
      PARAMETER (NS_max = 20, NH_max = 80)
      INTEGER NW_max
      PARAMETER (NW_max = 20)

C--------------------------------------------------------------------
C     SIBYLL common blocks containing event information       \FR'14
C--------------------------------------------------------------------

C     EVENT INFO COMMON
C     contains overall interaction properties, like
C     SQS : center-of-mass energy
C     S   :         "       "     squared
C     PTmin : low pt cut of QCD cross section, 
C             i.e. minimal pt of hard minijets
C     Xmin : low-x bound for PDFs, 
C            i.e. minimal momentum fraction of hard partons
C     Zmin : logarithm of that
C     KB : PID of beam hadron
C     KT() : PID of target
C     IAT : mass number of target
      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
      INTEGER KB,IAT,KT
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      
      DOUBLE PRECISION STR_mass_val, STR_mass_val_hyp, STR_mass_sea
      COMMON /S_CUTOFF/ STR_mass_val, STR_mass_val_hyp, STR_mass_sea
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)

      INTEGER IBMRDX,ITGRDX,IHMJDX,ISMJDX,ICSTDX,IINTDX
      COMMON /S_INDX/ IBMRDX(3),ITGRDX(NW_max,3),
     &     IHMJDX(NW_max*NH_max),IINTDX(NW_max),
     &     ISMJDX(NW_max*NS_max),ICSTDX(2*NW_max,3)
      DIMENSION NNSOF(NW_max),NNJET(NW_max),X2JET(NW_max)

c     internal types
      INTEGER NALL,JW,JJ,IREF,IREFG1,IREFG2,NSOF_JW,II
      DOUBLE PRECISION X1JJ,X2JJ,PTJET,FI,S_RNDM,SQSHALF,XM,
     &     X1S,X2S,PTSOF,PZ,EN     

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      SAVE

      if(Ndebug.gt.1) WRITE(LUN,*)
     &     ' SAMPLE_MINIJETS: (L,NW,NNJET,NNSOF):',
     &     L,NW,(NNJET(ii),ii=1,nw),(NNSOF(ii),ii=1,nw)

      IF(L.eq.0) THEN
         WRITE(LUN,*) 'SAMPLE_minijets: unknown particle? L=',L
         CALL SIB_REJECT('SAMPLE_minijets ')
      ENDIF
      
      NJET = 0
      NSOF = 0
      Nall = 0
      X1JET = 0.D0
      DO JW=1,NW
C...hard sea-sea interactions
         X2JET(JW) = 0.D0
         DO JJ=1,NNJET(JW)
           Nall = Nall+1
           NJET = NJET+1
           CALL SAMPLE_HARD (L,X1Jj,X2Jj,PTJET)
           X1JET = X1JET + X1Jj           
           X2JET(JW) = X2JET(JW)+X2Jj
           if(Ndebug.gt.2) THEN
              WRITE(LUN,*)
     &             ' SAMPLE_MINIJETS: hard JJ,JW,X1JET,X2JET(JW):',
     &             JJ,JW,X1JET,X2JET(JW)
              WRITE(LUN,*)
     &             '  X1,X2,PT:',X1JJ,X2JJ,PTJET
           ENDIF
           IF ((X2JET(JW).GT.0.9D0).OR.(X1JET.GT.0.9D0)) then
              if(Ndebug.gt.2) WRITE(LUN,*)
     &        ' SAMPLE_MINIJETS: not enough phase space',
     &             ' (Ncall,Njet,lbad):',Ncall,Njet,lBAD
              return
           ENDIF
           FI = TWOPI*S_RNDM(JJ)
           XM = SQS*sqrt(X1jj*X2jj)
           SQSHALF = 0.5D0*SQS
c           TH = ASIN(MIN((1.D0-EPS8),2.D0*PTJET/XM))
c     add gluon-gluon string to stack
           CALL ADD_PRTN
     &          (0.D0,0.D0,SQSHALF*(X1jj-X2jj),SQSHALF*(X1jj+X2jj),
     &          XM,100,0,0,Iref)
           CALL ADD_INT_REF(Iref,IINTDX(JW))
c     add gluon-gluon system to hard minijet index
           IHMJDX(NJET) = Iref
c     add gluons to stack
           CALL ADD_PRTN(PTJET*COS(FI),PTJET*SIN(FI),
     &          SQSHALF*X1jj,SQSHALF*X1jj,0.D0,0,1,0,Irefg1)
           CALL ADD_PRTN(-PTJET*COS(FI),-PTJET*SIN(FI),
     &          -SQSHALF*X2jj,SQSHALF*X2jj,0.D0,0,1,Iref,Irefg2)
c     set up references
c     minijet --> gluon1 --> gluon2 --> minijet
           CALL ADD_REF(Irefg1,Irefg2)
           CALL ADD_REF(Iref,Irefg1)

         ENDDO

C...soft sea-sea interactions 
         NSOF_JW = 0
         DO JJ=1,NNSOF(JW)-1
c     different soft distributions
            CALL SAMPLE_SOFT6 (STR_mass_sea,X1S,X2S,PTSOF)
            IF ((X2JET(JW)+X2S.LT.0.9D0).AND.(X1JET+X1S.LT.0.9D0)) THEN
               NSOF = NSOF+1
               Nall = Nall+1
               NSOF_JW = NSOF_JW+1
               X1JET = X1JET + X1S
               X2JET(JW) = X2JET(JW)+X2S
c     add to stack
c     add gluon-gluon string to stack
               XM = SQS*SQRT(X1S*X2S)
               SQSHALF = 0.5D0*SQS
               PZ = SQSHALF*(X1S-X2S)
               EN = SQSHALF*(X1S+X2S)
               FI = TWOPI*S_RNDM(JJ)
               CALL ADD_PRTN(0.D0,0.D0,PZ,EN,XM,10,0,0,Iref)
               CALL ADD_INT_REF(Iref,IINTDX(JW))
c     add gluons to stack
c     add gluon-gluon system to soft minijet index
               ISMJDX(NSOF) = Iref                              
               CALL ADD_PRTN(PTSOF*COS(FI),PTSOF*SIN(FI),
     &              SQSHALF*X1S,SQSHALF*X1S,0.D0,0,1,0,Irefg1)
               CALL ADD_PRTN(-PTSOF*COS(FI),-PTSOF*SIN(FI),
     &              -SQSHALF*X2S,SQSHALF*X2S,0.D0,0,1,Iref,Irefg2)
c     set up references
c     minijet --> gluon1 --> gluon2 --> minijet
               CALL ADD_REF(Irefg1,Irefg2)
               CALL ADD_REF(Iref,Irefg1)
               IF(Ndebug.gt.2)THEN
                  WRITE(LUN,*)
     &                 ' SAMPLE_MINIJETS: soft JJ,JW,X1JET,X2JET(JW):',
     &                 JJ,JW,X1JET,X2JET(JW)
                  WRITE(LUN,*)
     &                 '  X1,X2,PT:',X1s,X2s,PTSOF
               ENDIF
            ELSE
               IF(Ndebug.gt.1) WRITE(LUN,*)
     &        ' SAMPLE_MINIJETS: not enough phase space',
     &             ' (Ncall,Nsof,lbad):',Ncall,Njet,lBAD
               RETURN
            ENDIF
         ENDDO
         NNSOF(JW) = NSOF_JW+1
      ENDDO
      lbad = 0

      END
C=======================================================================

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
      IMPLICIT NONE

c     external types
      DOUBLE PRECISION SQS,SIGT,SIGEL,SIGINEL,SLOPE,RHO
      INTEGER L0

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN

c     external cross section tables
C     cross sections in model: 23rc1_sib23
      INTEGER K
      DOUBLE PRECISION SSIG_TOT(61,3)
      DOUBLE PRECISION SSIG(61,3)
      DOUBLE PRECISION SSIG_B(61,3)
      DOUBLE PRECISION SSIG_RHO(61,3)
c     internal type declarations
      DOUBLE PRECISION T,AL,ASQSMIN,ASQSMAX,DASQS
      INTEGER LL,L,J1,NSQS
      DIMENSION LL(39)
      SAVE
C     proton-proton:
C     total cross section
      DATA (SSIG_TOT(K,1),K=    1,   61) /
     &3.8328D+01,3.8267D+01,3.8435D+01,3.8838D+01,3.9463D+01,
     &4.0288D+01,4.1277D+01,4.2391D+01,4.3586D+01,4.4918D+01,
     &4.6354D+01,4.7836D+01,4.9394D+01,5.1050D+01,5.2835D+01,
     &5.4789D+01,5.6957D+01,5.9392D+01,6.2151D+01,6.5294D+01,
     &6.8883D+01,7.2529D+01,7.6458D+01,8.0673D+01,8.5172D+01,
     &8.9955D+01,9.5017D+01,1.0035D+02,1.0595D+02,1.1181D+02,
     &1.1790D+02,1.2423D+02,1.3077D+02,1.3751D+02,1.4444D+02,
     &1.5156D+02,1.5885D+02,1.6631D+02,1.7392D+02,1.8169D+02,
     &1.8960D+02,1.9766D+02,2.0584D+02,2.1416D+02,2.2260D+02,
     &2.3115D+02,2.3982D+02,2.4860D+02,2.5749D+02,2.6648D+02,
     &2.7556D+02,2.8475D+02,2.9403D+02,3.0340D+02,3.1287D+02,
     &3.2242D+02,3.3206D+02,3.4179D+02,3.5159D+02,3.6149D+02,
     &3.7146D+02/
C     inel. cross section
      DATA (SSIG(K,1),K=    1,   61) /
     &3.0881D+01,3.1156D+01,3.1540D+01,3.2046D+01,3.2673D+01,
     &3.3410D+01,3.4236D+01,3.5126D+01,3.6050D+01,3.7062D+01,
     &3.8139D+01,3.9280D+01,4.0476D+01,4.1740D+01,4.3092D+01,
     &4.4556D+01,4.6161D+01,4.7937D+01,4.9918D+01,5.2137D+01,
     &5.4629D+01,5.7057D+01,5.9635D+01,6.2361D+01,6.5230D+01,
     &6.8236D+01,7.1376D+01,7.4643D+01,7.8029D+01,8.1529D+01,
     &8.5138D+01,8.8847D+01,9.2654D+01,9.6552D+01,1.0054D+02,
     &1.0461D+02,1.0875D+02,1.1298D+02,1.1727D+02,1.2164D+02,
     &1.2607D+02,1.3057D+02,1.3512D+02,1.3974D+02,1.4441D+02,
     &1.4914D+02,1.5393D+02,1.5877D+02,1.6365D+02,1.6859D+02,
     &1.7357D+02,1.7860D+02,1.8368D+02,1.8880D+02,1.9397D+02,
     &1.9918D+02,2.0443D+02,2.0972D+02,2.1505D+02,2.2042D+02,
     &2.2583D+02/
C     slope parameter
      DATA (SSIG_B(K,1),K=    1,   61) /
     &1.0828D+01,1.1096D+01,1.1363D+01,1.1629D+01,1.1894D+01,
     &1.2159D+01,1.2424D+01,1.2688D+01,1.2953D+01,1.3217D+01,
     &1.3482D+01,1.3728D+01,1.3980D+01,1.4237D+01,1.4500D+01,
     &1.4770D+01,1.5047D+01,1.5333D+01,1.5632D+01,1.5945D+01,
     &1.6278D+01,1.6613D+01,1.6961D+01,1.7324D+01,1.7703D+01,
     &1.8100D+01,1.8515D+01,1.8949D+01,1.9404D+01,1.9880D+01,
     &2.0378D+01,2.0899D+01,2.1443D+01,2.2010D+01,2.2600D+01,
     &2.3212D+01,2.3845D+01,2.4499D+01,2.5173D+01,2.5867D+01,
     &2.6579D+01,2.7309D+01,2.8055D+01,2.8819D+01,2.9599D+01,
     &3.0394D+01,3.1205D+01,3.2031D+01,3.2870D+01,3.3724D+01,
     &3.4590D+01,3.5470D+01,3.6362D+01,3.7266D+01,3.8181D+01,
     &3.9109D+01,4.0047D+01,4.0995D+01,4.1955D+01,4.2924D+01,
     &4.3903D+01/
C     
      DATA (SSIG_RHO(K,1),K=    1,   61) /
     &-1.8490D-01,-1.2654D-01,-7.7648D-02,-3.7250D-02,-4.2495D-03,
     &2.2457D-02,4.3908D-02,6.1032D-02,7.4637D-02,8.5403D-02,
     &9.3897D-02,1.0058D-01,1.0583D-01,1.0995D-01,1.1318D-01,
     &1.1571D-01,1.1768D-01,1.1923D-01,1.2044D-01,1.2138D-01,
     &1.2212D-01,1.2269D-01,1.2314D-01,1.2349D-01,1.2376D-01,
     &1.2398D-01,1.2414D-01,1.2427D-01,1.2437D-01,1.2445D-01,
     &1.2451D-01,1.2456D-01,1.2460D-01,1.2463D-01,1.2465D-01,
     &1.2467D-01,1.2468D-01,1.2470D-01,1.2470D-01,1.2471D-01,
     &1.2472D-01,1.2472D-01,1.2472D-01,1.2473D-01,1.2473D-01,
     &1.2473D-01,1.2473D-01,1.2473D-01,1.2473D-01,1.2473D-01,
     &1.2473D-01,1.2473D-01,1.2473D-01,1.2473D-01,1.2473D-01,
     &1.2473D-01,1.2473D-01,1.2473D-01,1.2473D-01,1.2473D-01,
     &1.2473D-01/
C     pion-proton:
C     total cross section
      DATA (SSIG_TOT(K,2),K=    1,   61) /
     &2.3119D+01,2.3225D+01,2.3487D+01,2.3867D+01,2.4328D+01,
     &2.4886D+01,2.5529D+01,2.6249D+01,2.7038D+01,2.7890D+01,
     &2.8802D+01,2.9725D+01,3.0766D+01,3.1961D+01,3.3355D+01,
     &3.4994D+01,3.6931D+01,3.9223D+01,4.1928D+01,4.5104D+01,
     &4.8811D+01,5.2129D+01,5.5692D+01,5.9498D+01,6.3545D+01,
     &6.7832D+01,7.2350D+01,7.7094D+01,8.2059D+01,8.7235D+01,
     &9.2612D+01,9.8183D+01,1.0394D+02,1.0987D+02,1.1596D+02,
     &1.2221D+02,1.2862D+02,1.3518D+02,1.4188D+02,1.4871D+02,
     &1.5568D+02,1.6278D+02,1.7001D+02,1.7735D+02,1.8481D+02,
     &1.9239D+02,2.0008D+02,2.0788D+02,2.1578D+02,2.2378D+02,
     &2.3189D+02,2.4009D+02,2.4839D+02,2.5679D+02,2.6528D+02,
     &2.7386D+02,2.8253D+02,2.9129D+02,3.0014D+02,3.0908D+02,
     &3.1810D+02/
C     inel. cross section
      DATA (SSIG(K,2),K=    1,   61) /
     &1.9941D+01,2.0212D+01,2.0566D+01,2.0995D+01,2.1492D+01,
     &2.1955D+01,2.2477D+01,2.3056D+01,2.3685D+01,2.4360D+01,
     &2.5076D+01,2.5721D+01,2.6455D+01,2.7304D+01,2.8298D+01,
     &2.9466D+01,3.0844D+01,3.2465D+01,3.4364D+01,3.6574D+01,
     &3.9128D+01,4.1429D+01,4.3864D+01,4.6428D+01,4.9117D+01,
     &5.1926D+01,5.4847D+01,5.7875D+01,6.1006D+01,6.4233D+01,
     &6.7551D+01,7.0956D+01,7.4444D+01,7.8010D+01,8.1651D+01,
     &8.5363D+01,8.9145D+01,9.2994D+01,9.6906D+01,1.0088D+02,
     &1.0491D+02,1.0901D+02,1.1315D+02,1.1736D+02,1.2161D+02,
     &1.2592D+02,1.3028D+02,1.3469D+02,1.3915D+02,1.4366D+02,
     &1.4821D+02,1.5281D+02,1.5746D+02,1.6215D+02,1.6688D+02,
     &1.7166D+02,1.7648D+02,1.8134D+02,1.8625D+02,1.9119D+02,
     &1.9618D+02/
C     slope parameter
      DATA (SSIG_B(K,2),K=    1,   61) /
     &1.0120D+01,1.0270D+01,1.0416D+01,1.0559D+01,1.0698D+01,
     &1.0836D+01,1.0971D+01,1.1105D+01,1.1238D+01,1.1371D+01,
     &1.1502D+01,1.1435D+01,1.1392D+01,1.1377D+01,1.1395D+01,
     &1.1452D+01,1.1549D+01,1.1690D+01,1.1878D+01,1.2118D+01,
     &1.2413D+01,1.2781D+01,1.3163D+01,1.3558D+01,1.3967D+01,
     &1.4391D+01,1.4829D+01,1.5282D+01,1.5751D+01,1.6236D+01,
     &1.6738D+01,1.7256D+01,1.7791D+01,1.8343D+01,1.8912D+01,
     &1.9498D+01,2.0100D+01,2.0718D+01,2.1351D+01,2.1999D+01,
     &2.2661D+01,2.3338D+01,2.4029D+01,2.4733D+01,2.5451D+01,
     &2.6182D+01,2.6926D+01,2.7682D+01,2.8450D+01,2.9231D+01,
     &3.0023D+01,3.0827D+01,3.1642D+01,3.2468D+01,3.3305D+01,
     &3.4152D+01,3.5010D+01,3.5878D+01,3.6757D+01,3.7645D+01,
     &3.8543D+01/
C     
      DATA (SSIG_RHO(K,2),K=    1,   61) /
     &-6.7332D-02,-3.0879D-02,-5.4256D-04,2.4410D-02,4.4739D-02,
     &6.1172D-02,7.4371D-02,8.4920D-02,9.3315D-02,9.9976D-02,
     &1.0525D-01,1.0941D-01,1.1269D-01,1.1528D-01,1.1731D-01,
     &1.1891D-01,1.2016D-01,1.2115D-01,1.2192D-01,1.2253D-01,
     &1.2300D-01,1.2338D-01,1.2367D-01,1.2390D-01,1.2408D-01,
     &1.2422D-01,1.2433D-01,1.2442D-01,1.2449D-01,1.2454D-01,
     &1.2458D-01,1.2462D-01,1.2464D-01,1.2466D-01,1.2468D-01,
     &1.2469D-01,1.2470D-01,1.2471D-01,1.2471D-01,1.2472D-01,
     &1.2472D-01,1.2472D-01,1.2473D-01,1.2473D-01,1.2473D-01,
     &1.2473D-01,1.2473D-01,1.2473D-01,1.2473D-01,1.2473D-01,
     &1.2473D-01,1.2473D-01,1.2473D-01,1.2473D-01,1.2473D-01,
     &1.2473D-01,1.2473D-01,1.2473D-01,1.2473D-01,1.2473D-01,
     &1.2473D-01/
C     kaon-proton:
C     total cross section
      DATA (SSIG_TOT(K,3),K=    1,   61) /
     &1.8299D+01,1.8827D+01,1.9408D+01,2.0016D+01,2.0633D+01,
     &2.1318D+01,2.2044D+01,2.2810D+01,2.3615D+01,2.4458D+01,
     &2.5339D+01,2.6253D+01,2.7209D+01,2.8235D+01,2.9372D+01,
     &3.0683D+01,3.2250D+01,3.4173D+01,3.6576D+01,3.9602D+01,
     &4.3417D+01,4.6380D+01,4.9560D+01,5.2954D+01,5.6563D+01,
     &6.0384D+01,6.4411D+01,6.8639D+01,7.3062D+01,7.7674D+01,
     &8.2464D+01,8.7426D+01,9.2551D+01,9.7831D+01,1.0326D+02,
     &1.0883D+02,1.1454D+02,1.2037D+02,1.2634D+02,1.3243D+02,
     &1.3864D+02,1.4496D+02,1.5139D+02,1.5793D+02,1.6458D+02,
     &1.7133D+02,1.7817D+02,1.8512D+02,1.9215D+02,1.9928D+02,
     &2.0650D+02,2.1380D+02,2.2119D+02,2.2867D+02,2.3623D+02,
     &2.4387D+02,2.5160D+02,2.5940D+02,2.6728D+02,2.7524D+02,
     &2.8328D+02/
C     inel. cross section
      DATA (SSIG(K,3),K=    1,   61) /
     &1.6131D+01,1.6687D+01,1.7256D+01,1.7835D+01,1.8414D+01,
     &1.8990D+01,1.9596D+01,2.0228D+01,2.0887D+01,2.1572D+01,
     &2.2282D+01,2.3007D+01,2.3748D+01,2.4525D+01,2.5373D+01,
     &2.6337D+01,2.7475D+01,2.8859D+01,3.0574D+01,3.2718D+01,
     &3.5399D+01,3.7521D+01,3.9768D+01,4.2138D+01,4.4626D+01,
     &4.7228D+01,4.9939D+01,5.2752D+01,5.5666D+01,5.8673D+01,
     &6.1770D+01,6.4952D+01,6.8215D+01,7.1555D+01,7.4969D+01,
     &7.8453D+01,8.2007D+01,8.5626D+01,8.9308D+01,9.3052D+01,
     &9.6855D+01,1.0072D+02,1.0463D+02,1.0861D+02,1.1263D+02,
     &1.1671D+02,1.2084D+02,1.2501D+02,1.2924D+02,1.3352D+02,
     &1.3784D+02,1.4220D+02,1.4662D+02,1.5107D+02,1.5558D+02,
     &1.6012D+02,1.6471D+02,1.6934D+02,1.7401D+02,1.7872D+02,
     &1.8348D+02/
C     slope parameter
      DATA (SSIG_B(K,3),K=    1,   61) /
     &8.8352D+00,9.1363D+00,9.4011D+00,9.6374D+00,9.8515D+00,
     &1.0048D+01,1.0230D+01,1.0402D+01,1.0564D+01,1.0720D+01,
     &1.0870D+01,1.1058D+01,1.1205D+01,1.1322D+01,1.1419D+01,
     &1.1511D+01,1.1611D+01,1.1734D+01,1.1897D+01,1.2116D+01,
     &1.2413D+01,1.2781D+01,1.3163D+01,1.3558D+01,1.3967D+01,
     &1.4391D+01,1.4829D+01,1.5282D+01,1.5751D+01,1.6236D+01,
     &1.6738D+01,1.7256D+01,1.7791D+01,1.8343D+01,1.8912D+01,
     &1.9498D+01,2.0100D+01,2.0718D+01,2.1351D+01,2.1999D+01,
     &2.2661D+01,2.3338D+01,2.4029D+01,2.4733D+01,2.5451D+01,
     &2.6182D+01,2.6926D+01,2.7682D+01,2.8450D+01,2.9231D+01,
     &3.0023D+01,3.0827D+01,3.1642D+01,3.2468D+01,3.3305D+01,
     &3.4152D+01,3.5010D+01,3.5878D+01,3.6757D+01,3.7645D+01,
     &3.8543D+01/
C     
      DATA (SSIG_RHO(K,3),K=    1,   61) /
     &-2.4506D-02,9.2028D-03,3.5513D-02,5.5961D-02,7.1799D-02,
     &8.4036D-02,9.3471D-02,1.0074D-01,1.0632D-01,1.1061D-01,
     &1.1391D-01,1.1643D-01,1.1837D-01,1.1986D-01,1.2100D-01,
     &1.2187D-01,1.2254D-01,1.2305D-01,1.2345D-01,1.2375D-01,
     &1.2398D-01,1.2416D-01,1.2429D-01,1.2439D-01,1.2447D-01,
     &1.2453D-01,1.2458D-01,1.2462D-01,1.2464D-01,1.2467D-01,
     &1.2468D-01,1.2469D-01,1.2470D-01,1.2471D-01,1.2472D-01,
     &1.2472D-01,1.2472D-01,1.2473D-01,1.2473D-01,1.2473D-01,
     &1.2473D-01,1.2473D-01,1.2473D-01,1.2473D-01,1.2473D-01,
     &1.2473D-01,1.2473D-01,1.2473D-01,1.2473D-01,1.2473D-01,
     &1.2473D-01,1.2473D-01,1.2473D-01,1.2473D-01,1.2473D-01,
     &1.2473D-01,1.2473D-01,1.2473D-01,1.2473D-01,1.2473D-01,
     &1.2473D-01/
      
      DATA LL /5*0,3*2,4*3,2*1,19*0,6*1/


      L = L0
      NSQS = 61
      ASQSMIN = 1.D0
      ASQSMAX = 7.D0
      DASQS = (ASQSMAX-ASQSMIN)/DBLE(NSQS-1)

      IF(NSQS.LE.0) THEN
        WRITE(LUN,'(//,1X,A)') 
     &    'SIB_SIGMA_EXT: interpolation table not initialized.'
        STOP
      ENDIF
      IF(IABS(L).gt.39)THEN
         WRITE(LUN,*)     
     &        ' SIB_SIGMA_EXT: unknown beam particle!',L
         STOP
      ENDIF
      IF(L.GT.3) L=LL(IABS(L))
      IF(L.EQ.0)THEN
         WRITE(LUN,*)     
     &        ' SIB_SIGMA_EXT: unknown beam particle!', L
         STOP
      ENDIF
        
      AL = LOG10(SQS)
      J1 = INT((AL-1.D0)*10.D0 + 1)
      if((j1.lt.1).or.(j1.gt.NSQS)) then
        if (ndebug .gt. 0) 
     *    write (LUN,'(1x,a,i3,1p,e12.3)') 
     &      ' SIB_SIGMA_EXT: energy out of range ',L,sqs
      endif
      if((j1.lt.1).or.(j1.ge.NSQS)) then
        J1 = min(J1,NSQS-1)
        J1 = max(J1,1)
      endif
      T = (AL-1.D0)*10.D0 - DBLE(J1-1)
      SIGT    = SSIG_TOT(J1,L)*(1.D0-T) + SSIG_TOT(J1+1,L)*T
      SIGINEL = SSIG(J1,L)*(1.D0-T) + SSIG(J1+1,L)*T
      SIGEL   = SIGT-SIGINEL
      SLOPE   = SSIG_B(J1,L) *(1.D0-T) + SSIG_B(J1+1,L)*T
      RHO     = SSIG_RHO(J1,L) *(1.D0-T) + SSIG_RHO(J1+1,L)*T

      END
C=======================================================================

      SUBROUTINE SAMPLE_PROJECTILE
     +     (KID,KINT,LRMNT,XCHG,XJET,XX,PX,PY,IFL,KID1,IREJ)

C-----------------------------------------------------------------------
C...  Subroutine to sample sea and valence quarks in a hadron.
C.    variables are stored in xx,px,py and ifl arrays.
C.    for each interaction the hadron undergoes there is one 
C.    pair of partons attached to the ends of two strings
C.    (one cut pomeron)
C.    In addition flavor and momentum may be set aside for the remnant
C.    arrays are filled: rmnt1,rmnt2, c.str1,c.str2, etc..
C.    i.e. positions 1 and 2 are reserved for remnant.
C.
C.    Input: KINT  : number of interactions the hadron takes part in
C.           KID   : particle id of hadron
C.           LRMNT : remnant excitation flag,
C.                   defines if valence quarks need to be sampled
C.           XCHG  : flavor exchange prob. between remnant and 
C.                   central strings
C.           XJET  : momentum fraction already asigned to minijets
C.           IREJ  : rejection flag, default set in calling routine
C.
C.    Output: XX,IFL,PX,PY  : arrays of momentum fractions, flavor 
C.                            and transverse momentum
C.            KID1 : new hadron id (in case of flavor exchange)
C-------------------------------------------------------------------      
      IMPLICIT NONE

C     include COMMONs
      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      
      DOUBLE PRECISION STR_mass_val, STR_mass_val_hyp, STR_mass_sea
      COMMON /S_CUTOFF/ STR_mass_val, STR_mass_val_hyp, STR_mass_sea
      INTEGER NW_max
      PARAMETER (NW_max = 20)
C--------------------------------------------------------------------
C     SIBYLL common blocks containing event information       \FR'14
C--------------------------------------------------------------------

C     EVENT INFO COMMON
C     contains overall interaction properties, like
C     SQS : center-of-mass energy
C     S   :         "       "     squared
C     PTmin : low pt cut of QCD cross section, 
C             i.e. minimal pt of hard minijets
C     Xmin : low-x bound for PDFs, 
C            i.e. minimal momentum fraction of hard partons
C     Zmin : logarithm of that
C     KB : PID of beam hadron
C     KT() : PID of target
C     IAT : mass number of target
      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
      INTEGER KB,IAT,KT
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)

      INTEGER ICHP,ISTR,IBAR
      COMMON /S_CHP/ ICHP(99), ISTR(99), IBAR(99)

      INTEGER IISO,ISPN
      COMMON /S_SPN/ IISO(99), ISPN(99)

      INTEGER ICHM
      COMMON /S_CHM/ ICHM(99)

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      INTEGER ITRY, NREJ
      COMMON /S_CNT/ ITRY(20), NREJ(20)

C     input type declarations
      INTEGER KID,KINT,LRMNT
      DOUBLE PRECISION XCHG,XJET
      
C     output type declarations
      DOUBLE PRECISION XX,PX,PY
      INTEGER IFL,KID1,IREJ
      DIMENSION XX(2*NW_max+2),PX(2*NW_max+2),PY(2*NW_max+2),
     &     IFL(2*NW_max+2)

c     local type declarations
      INTEGER ICNT1,ICNT2,J,JJ,j1,j2,j3,j4,KRMNT,IRNK,
     &     IDXVAL,IDX,ISWTD,i,IFLS,NVAL,NSEA,IR,IDUM,IDUM2,KIDA,IMRG2HAD
      DOUBLE PRECISION XSEAJET,XVAL,XMINA,XMINA_SEA,GAMMA,XREM,XMINA2,
     &     XMAX2,ALPHA,XM2DIS,ASUP,XMAX,XQM,S_RNDM,
     &     CHIDIS,CHI,GAMDIQ,XSUPP,XSUPP1,PAR53_def,PAR5_def,PAR6_def,
     &     PAR7_def,PAR143_def,XSUM,STR_mass,PTS,XSCL
      SAVE
      DATA ICNT1,ICNT2 /0,0/
      
C..   initialization
      ITRY(3) = 0
      XVAL = 0.D0
      XSCL = 1.D0
      XSEAJET = 0.D0
      XSUM = 0.D0
      DO J=1,KINT               ! zero arrays
         j1 = 1+2*(j-1)
         j2 = j1 + 1
         j3 = 3+2*(j-1)
         j4 = j3 + 1
         XX(j1) = 0.D0
         XX(j2) = 0.D0
         XX(j3) = 0.D0
         XX(j4) = 0.D0
         PX(j1) = 0.D0
         PX(j2) = 0.D0
         PX(j3) = 0.D0
         PX(j4) = 0.D0
      ENDDO

      KRMNT = MIN(LRMNT,1)

      IF(ndebug.gt.3) 
     +     WRITE(LUN,*)
     +     ' SAMPLE_PROJECTILE: KID,KINT,KRMNT,XCHG,XJET,IREJ',
     +     KID,KINT,KRMNT,XCHG,XJET,IREJ

      KID1 = KID
      KIDA = IABS(KID)
      
c     number of valence quarks to sample
c     if remnant is resolved (krmnt=1) no valence pair needed
      Nval = 2*(1-KRMNT)

c     number of sea quarks to sample (one pair per interaction)
c     if remnant is not resolved then on pair less is needed 
c     (valence pair takes role of one sea pair)
      Nsea = 2*(KINT-(1-KRMNT))

      IF(ndebug.gt.3) 
     +     WRITE(LUN,*)
     +     ' SAMPLE_PROJECTILE: number of partons to sample ',
     +     '(tot,val,sea):',Nval+Nsea,Nval,Nsea

c     change proton splitting to enhance charge exchange by allowing
c     ud more often than uu, default scenario is ud,du,uu: 3:1:2
      PAR53_def = PAR(53)
      PAR(53) = PAR(84)
c     change proton splitting in case no remnant is present
      IF(LRMNT.eq.0) PAR(53) = PAR(127)

 20   ITRY(3) = ITRY(3) + 1
      IF(ITRY(3).gt.NREJ(3)) THEN
         ICNT1 = ICNT1 + 1
         IF(ICNT1.lt.10)THEN
          if (NDEBUG.gt.0) then
            WRITE(LUN,*)' SAMPLE_PROJECTILE: trials exceeded! return..'
            WRITE(LUN,*)
     +           '  KID,KINT,KRMNT,XCHG,XJET,XVAL,IREJ,NCALL',
     +           KID,KINT,KRMNT,XCHG,XJET,XVAL,IREJ,NCALL           
          endif
         ENDIF
         PAR(53) = PAR53_def
         RETURN 
      ENDIF

C...  kinematic limits
 22   XSEAJET = XJET
      IF(KRMNT.eq.0)THEN
c     minimal momentum fraction for valences
         XMINA = 2.D0*STR_mass_val/SQS
c     default for valence quarks: 0.35 , xmin@10GeV = 0.07
c     taken from COMMON s_cutoff
         IF(ISTR(KIDA)*IBAR(KIDA).ne.0)
     &        XMINA = 2.D0*STR_mass_val_hyp/SQS
      ELSE
         IF(IPAR(47).eq.4.or.IPAR(47).eq.4.or.IPAR(47).eq.6)then
c     no valence sampling model
c     if remnant present then the minimal remnant mass has to be provided
            XMINA = PAR(96)*AM(IABS(KID))/SQS            
         ELSEIF(IPAR(47).lt.4)THEN
c     valences sampled, even if combined again in remnant
            XMINA = 2.D0*STR_mass_val/SQS
         ELSEIF(IPAR(47).eq.7)THEN
c     minimal remnant mass not requiered,
c     mass is taken from central strings anyway..
            XMINA = AM(IABS(KID))/SQS
         ENDIF
      ENDIF
         
c     minimal momentum fraction for sea partons
      IF(IPAR(47).eq.0.or.IPAR(47).eq.3)THEN
c     same as valence quarks
         STR_mass = STR_mass_val
      ELSEIF(IPAR(47).eq.1.or.IPAR(47).eq.2.or.IPAR(47).gt.4)THEN
c     set by parameter
         STR_mass = PAR(87)
      ELSEIF(IPAR(47).eq.4)THEN
c     same as soft minijets
         STR_mass = STR_mass_sea
      ENDIF
      IF(IPAR(72).eq.2.and.KINT.gt.1)THEN
         STR_mass = STR_mass * PAR(118)
      ENDIF
      XMINA_SEA = 2.D0*STR_mass/SQS
c     default for sea quarks: 1.0 , xmin@10GeV = 0.2
c     taken from COMMON s_cutoff or s_cflafr
c     should be the same as min. string mass in SAMPLE_SOFT !

c     dependence on number of interactions
      IF(IPAR(72).eq.1.and.KINT.gt.1)THEN
         XMINA_SEA = XMINA_SEA * PAR(118)
      ENDIF

C..   check if enough energy left to sample all partons
      IF (1.D0-XJET.LT.(Nsea*XMINA_SEA+2.D0*XMINA))THEN
         ICNT2 = ICNT2 + 1
         IF(ICNT2.le.10)THEN
            IF(NDEBUG.gt.3)THEN
               write(lun,*)' SAMPLE_PROJECTILE: rejection!' 
               write(lun,*)'  too little energy to sample all partons!'
               write(lun,*)'  (NW,Ntot,Nval,Nsea,XMIN,XMIN*N',
     &              'XREM,XALL,NCALL,ICNT:)',KINT,nval+nsea,Nval,nsea,
     &              2*xmina,nsea*xmina_sea,1.D0-xjet,
     &              Nsea*XMINA_SEA+2*XMINA,NCALL,ICNT2
               IF(ICNT2.eq.10) write(lun,*)' last warning ! good luck..'
            ENDIF
         ENDIF

         IREJ = 2
         PAR(53) = PAR53_def
         RETURN
      ENDIF


C...  sample sea partons
c     if no additional partons need to be sampled 
C     jump straight to valence quarks
      IF(Nsea.EQ.0) GOTO 100

C     select sea quark model
      IF(IPAR(47).eq.0.or.IPAR(47).eq.3.or.IPAR(47).eq.4.or.
     &     IPAR(47).eq.5.or.IPAR(47).eq.7)THEN
         GAMMA = PAR(103)
         IF(IPAR(73).eq.1.and.KINT.gt.1) GAMMA = PAR(119)
         CALL SAMPLE_SEA_TOT
     &        (KRMNT,KINT,NSEA,GAMMA,XJET,STR_MASS,XSEAJET,XX)

      ELSEIF(IPAR(47).eq.1)THEN
c     sample from 1/x individually then reject if too large
         XREM = 0.D0
         XMINA2 = XMINA_SEA ** 2
         XMAX2 = 0.8D0**2
         ALPHA = 1.D0
         DO WHILE ( XREM .lt. 2*XMINA )
            XREM = 1.D0-XJET
            IF(NDEBUG.gt.3)
     &           WRITE(LUN,*) '  N,XREM,XMINA,XMIN2,XMAX2,ALPHA',
     &           Nsea,XREM,XMINA_SEA,XMINA2,XMAX2,ALPHA
            DO j=1,Nsea
               jj = 2 + j
               IF(KRMNT.eq.0) jj = 4+j
               XX(jj) = XM2DIS(XMINA2,XMAX2,ALPHA)
               IF(NDEBUG.gt.3) 
     &           WRITE(LUN,*) ' J,X,XREM',JJ,XX(JJ),XREM
               XREM = XREM - XX(jj)
            ENDDO
         ENDDO
         XSEAJET = 1.D0-XREM

      ELSEIF(IPAR(47).eq.2.or.IPAR(47).eq.6)THEN
c     sample from (1-x)**b / x with common mass constraint
         XREM = 1.D0-XJET
         XMAX = PAR(88)
         ALPHA = PAR(85)
         ASUP = PAR(86)
         XQM = STR_mass
         CALL SAMPLE_SEA_INDV(KRMNT,XMINA,XMINA_SEA,NSEA,
     &        XREM,ALPHA,ASUP,XQM,XMAX,XX,IR)
         IF(IR.ne.0)THEN
            IREJ = IR
            PAR(53) = PAR53_def
            RETURN
         ENDIF

         XSEAJET = 1.D0-XREM
      ENDIF

C...  sample sea flavor: u,d,s,c
c     write to ifl after valences..
      DO J=1+Nval/2,KINT
         j3 = 3+2*(j-1)
         j4 = j3 + 1
c     turn on strange sea..
         IF(IPAR(29).eq.1)THEN
            IF(IPAR(69).ne.0)THEN
c     sample asymmetric u,d
               IFL(j3) = MIN(2,1+INT((2.D0+PAR(114))*S_RNDM(KID)))
c     sample strange
               IFLS = 3*(INT((2+PAR(43))*S_RNDM(j3))/2)
               IFL(j3) = MAX(IFL(j3),IFLS)
            else
               IFL(j3) = 1+INT((2.D0+PAR(43))*S_RNDM(j4))
            endif
c     sample charm
c     scale up for mesons
            IF(IPAR(76).eq.1) XSCL=XSCL+(1-IABS(IBAR(KIDA)))*PAR(126)
            IF(IFL(j3).eq.3.and.S_RNDM(kid).lt.PAR(97)*PAR(125)*XSCL)
     &           IFL(j3) = 4
         ELSE
            IFL(j3) = INT(1.5D0+S_RNDM(KID))
         ENDIF
         IFL(j4) = -IFL(j3)
         IF(NDEBUG.gt.3) 
     &        WRITE(LUN,*) '  flavor: JW,FLV1,FLV2',J,IFL(j3),IFL(j4)

C...  sample sea pt
 33      IF(IPAR(49).eq.1)THEN
c     in-string pt for sea partons
c     flavor and cm energy dependent avg, exponential dist.
c     avg pt (defined in subroutine ptsetup ):
c     u,d : PAR(46)+PAR(68)*log10(sqs/20.D0)**2
c     s:    PAR(47)+PAR(70)*log10(sqs/20.D0)**2
c     diq:  PAR(48)+PAR(69)*log10(sqs/20.D0)**2
            CALL PTDIS_4FLV (IFL(j3),PX(j3),PY(j3))
            PX(j4) = -PX(j3)
            PY(j4) = -PY(j3)

         ELSEIF(IPAR(49).eq.2)THEN
c     'primordial' pt
c     c.m. energy dependent avg, exponential
c     same for all flavors
c     avg: PAR(49)+PAR(69)*log10(sqs/20.)**2
            CALL PTDIS_4FLV (10,PX(j3),PY(j3))
            PX(j4) = -PX(j3)
            PY(j4) = -PY(j3)

         ELSEIF(IPAR(49).eq.3)THEN
c     constant pt
            PX(j3) = EPS5
            PY(j3) = EPS5
            PX(j4) = -PX(j3)
            PY(j4) = -PY(j3)

         ELSEIF(IPAR(49).eq.4)THEN
c     sea pt, same as primordial but different params..
c     c.m. energy dependent avg, exponential
c     same for all flavors
c     avg: PAR(132)
            CALL PTDIS_4FLV (30,PX(j3),PY(j3))
            PX(j4) = -PX(j3)
            PY(j4) = -PY(j3)
         ENDIF
c     limit parton virtuality         
         PTS = MAX(PX(j3)**2+PY(j3)**2,EPS10)
         IF((XX(j3)**2+XX(J4)**2)/PTS.lt.8.D0*PAR(122)/S) GOTO 33
         IF(NDEBUG.gt.3) 
     &        WRITE(LUN,*)'  pt: JW,PX,PY,pt',J,Px(j3),Py(j3),sqrt(pts)
      ENDDO     

C...  Prepare the valence partons
 100  XVAL=1.D0-XSEAJET
      IF(ndebug.gt.3)
     &     write(lun,*) ' SAMPLE_PROJECTILE: val. fraction remaining:',
     &     XVAL

      IF(IPAR(47).eq.7)THEN
c     no remnant, sample valence quarks
         IF(KRMNT.eq.0) THEN
c     enough momentum left?
            IF (XVAL.LT.XMINA) goto 20 ! reject sea kinematics
         ELSE
c     sample remnant
            IF(IPAR(53).eq.1)THEN
c     momentum dis: x**alpha
               IF(S_RNDM(KID).gt.XVAL**(PAR(100)+1)) GOTO 22
            ENDIF
c     split remnant momentum into partons, just to fill slots
            
         ENDIF            
      ELSE
         IF(KRMNT.eq.0.or.IPAR(47).lt.4)THEN
            IF (XVAL.LT.XMINA) goto 20 ! reject sea kinematics
         ENDIF
c     remnant momentum fraction
         IF(KRMNT.ne.0.and.IPAR(53).eq.1)THEN
            IF(S_RNDM(KID).gt.XVAL**(PAR(100)+1)) GOTO 22
         ENDIF
      ENDIF
c     valence quarks are in 1,2 of IFL,XX etc.
      IDXVAL = 3
      IF(KRMNT.ne.0) IDXVAL = 1
      CALL HSPLI (KID,IFL(IDXVAL),IFL(IDXVAL+1))
 110  CHI = CHIDIS(KID,IFL(IDXVAL),IFL(IDXVAL+1))
      XX(IDXVAL) = MAX(CHI*XVAL,XMINA)
      XX(IDXVAL) = MIN(XX(IDXVAL),XVAL-XMINA)
C     FOR MESONS, SPLIT ENERGY SYMMETRICALLY.
      IF (IABS(KID).LT.13.AND.S_RNDM(0).LE.0.5D0) 
     &     XX(IDXVAL)=XVAL-XX(IDXVAL)
      XX(IDXVAL+1)=XVAL-XX(IDXVAL)
      IF(ndebug.gt.3)
     &     write(lun,*) ' SAMPLE_PROJECTILE: val. sampled (x1,x2):',
     &     XX(IDXVAL),XX(IDXVAL+1)
c     for baryons force diq distribution
      IF(IBAR(IABS(KID)).ne.0.and.IPAR(47).ne.7)THEN
         IF(IPAR(52).eq.1)THEN
            GAMDIQ=PAR(95)
            IF(S_RNDM(KID).gt.XX(IDXVAL+1)**(GAMDIQ+1)) GOTO 110
         ELSE
            IF(KRMNT.eq.0.or.IPAR(47).lt.4.and.IPAR(53).eq.0)THEN
c     force diquark distribution
               GAMDIQ=PAR(95)
               IF(S_RNDM(KID).gt.XX(IDXVAL+1)**(GAMDIQ+1)) GOTO 20
            ENDIF
         ENDIF
      ENDIF
C...  val. quark transverse momentum
      CALL PTDIS_4FLV (10,PX(IDXVAL),PY(IDXVAL))
      PX(IDXVAL+1) = -PX(IDXVAL)
      PY(IDXVAL+1) = -PY(IDXVAL)     
      IF(ndebug.gt.3)
     &     write(lun,*) ' SAMPLE_PROJECTILE: val. pt (px,py):',
     &     PX(IDXVAL),PY(IDXVAL)

C...  exchange flavor between central strings and remnant
c     there is one pair of strings for each interaction with another hadron
c     in general allowed for both flavors but diquarks usually strongly suppressed
c     Xchg : prob. of flv exchange between strgs and rmnt
      IF(KRMNT.ne.0)THEN
         do idx=1,2
            iswtd = 0
            i = 1
            XSUPP = 1.D0
            IF(iabs(ifl(idx)).gt.10)THEN
c     suppress exchange of diq: prob_exchange = prob0 * xsupp
               XSUPP = PAR(83)
            ELSEIF(IPAR(46).eq.2)THEN
c     suppress exchange for fast quark ( usually in mesons )
               IF(xx(idx).gt.xx(3-idx)) XSUPP = PAR(139)
            ENDIF
            DO WHILE (ISWTD.eq.0.and.i.le.KINT)
c     sea flavor index
               jj = idx+2*i
c     forbid exchange for charmed hadrons if sea pair is charmed too
c     needed to avoid double charmed particles
               XSUPP1 = XSUPP
               IF(IABS(KID).gt.50.and.IABS(IFL(JJ)).eq.4) XSUPP1 = 0.D0
               if(S_RNDM(I).lt.XCHG*XSUPP1) THEN               
c     exchange flavor between remnant and sea
                  CALL ISWTCH_LMNTS(ifl(jj),ifl(idx))
c     also exchange momentum fraction
                  IF(IPAR(46).ne.0) CALL SWTCH_LMNTS(xx(jj),xx(idx))
c     change flavor id accordingly, i.e. reassamble remnant from new flavor
                  IF(IPAR(58).eq.0)THEN
c     combine to any hadron that matches flavor wise, ignoring (iso)spin
                     CALL SIB_I4FLAV(ifl(1),ifl(2),idum,idum2,KID1)
                  ELSEIF(IPAR(58).eq.1)THEN
c     combine to lightest hadron
                     KID1 = IMRG2HAD(IFL(1),IFL(2))
                  ELSEIF(IPAR(58).eq.2.or.IPAR(58).eq.3)THEN
c     combine to any hadron that matches flavor wise, ignoring (iso)spin
c     set vector meson rate
                     PAR5_def = PAR(5)
                     PAR(5) = PAR(104)
c     set strange vector rate
                     PAR6_def = PAR(6)
                     PAR(6) = PAR(121)
c     set spin3/2 vs spin1/2 rate
                     PAR7_def = PAR(7)
                     PAR(7) = PAR(105)
c     set rho / omega-phi rate
                     PAR143_def = PAR(143)
                     if(ibar(iabs(kb)).eq.0.and.IPAR(85).eq.1)
     &                    PAR(143) = PAR(145)
                     irnk = 0                     
                     IF(IPAR(58).eq.3) irnk = 1                     
                     CALL SIB_I4FLAV(ifl(1),ifl(2),irnk,idum2,KID1)
                     PAR(5) = PAR5_def
                     PAR(6) = PAR6_def
                     PAR(7) = PAR7_def
                     PAR(143) = PAR143_def

c     reject spin1,isospin singlett
                     IF(KID1.eq.32.and.PAR(111).gt.S_RNDM(KID1))
     &                    KID1 = 27
                  ENDIF
                  iswtd = 1
               endif
               i = i + 1
            ENDDO
         enddo
      ENDIF
      IF(ndebug.gt.3)THEN        
         WRITE(LUN,*)' SAMPLE_PROJECTILE: rmnt PID,NTRY: ',KID1,ITRY(3)
         WRITE(LUN,*)' SAMPLE_PROJECTILE: output: I,FLV,PX,PY,X,XSUM'
      ENDIF
      XSUM = XJET
      DO j=IDXVAL,2*(KINT+Krmnt)+2*(1-Krmnt)
         XSUM = XSUM + XX(j)
         IF(NDEBUG.gt.3) WRITE(LUN,*) j,IFL(j),PX(J),PY(J),XX(j),XSUM
      ENDDO
      IF(ABS(XSUM-1.D0).gt.EPS3) THEN
         WRITE(LUN,*)' SAMPLE_PROJECTILE: parton sum incomplete!',
     &        '(ID,XSUM,NCALL):' , KID1,XSUM, NCALL,' aborting..'
         CALL SIB_REJECT('SAMPLE_PROJECTIL')
      ENDIF
      IREJ = 0

      END
C=======================================================================

      SUBROUTINE DECSIB 

C-----------------------------------------------------------------------
C...Decay all unstable particle in Sibyll
C.  decayed particle have the code increased by 10000
C
C   changed to allow for multiple calls to DECSIB in one event
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      DOUBLE PRECISION CBR
      INTEGER KDEC,LBARP,IDB
      COMMON /S_CSYDEC/ CBR(223+16+12+8), KDEC(1338+6*(16+12+8)),
     &     LBARP(99), IDB(99)
C     The final particle output is contained in COMMON /S_PLIST/    
C     NP           : number of final particles
C     P(1:NP, 1:5) : 4-momenta + masses of the final particles 
C     LLIST (1:NP) : codes of final particles
      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP
      INTEGER LLIST1
      COMMON /S_PLIST1/ LLIST1(8000)

      INTEGER NFORIG,NPORIG,NIORIG,IPFLAG,IIFLAG,KINT
      COMMON /S_PARTO/ NFORIG(NP_max),NPORIG(NP_max),NIORIG(NP_max),
     &IPFLAG,IIFLAG,KINT
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)

      INTEGER LRNK
      COMMON /SIB_RNK/ LRNK(8000)
      DIMENSION P0(5), LL(10), PD(10,5)
      SAVE

c     call pythia decay routine      
c      IF(IPAR(44).eq.1) CALL PYDEC

c     decay with sibyll
      NN = 1
      IF(IPAR(44).ne.1)THEN
         DO J=1,NP
            LLIST1(J) = 0
         ENDDO
      ENDIF
      DO WHILE (NN .LE. NP)
         L= LLIST(NN)
         LA = IABS(L)
         if(LA.lt.100) then
           IF (IDB(LA) .GT. 0)  THEN
              DO K=1,5
                P0(K) = P(NN,K)
              ENDDO
              CALL DECPAR (L,P0,ND,LL,PD)
              LLIST(NN) = LLIST(NN)+ISIGN(10000,LLIST(NN))
              DO J=1,ND
                NP = NP+1
                if(NP.gt.8000) then
                  write(LUN,'(1x,a,2i8)') 
     &              ' DECSIB: no space left in S_PLIST (NP,ND):',NP,ND
                  NP = NP-1
                  return
                endif
                DO K=1,5
                  P(NP,K) = PD(J,K)
                ENDDO
                LLIST(NP)=LL(J)
                LLIST1(NP)=NN
                LRNK(NP)=LRNK(NN)
                NPORIG(NP)= NPORIG(NN)
                niorig(NP)= NIORIG(NN)
                NFORIG(NP) = L
              ENDDO
           ENDIF
         endif
         NN = NN+1
      ENDDO
      
c      CALL SIB_LIST(20)
      END

C=======================================================================
      SUBROUTINE DEC_SET(Imode)
C-----------------------------------------------------------------------
C...  Routine to store/write decay configuration
C-----------------------------------------------------------------------     
      IMPLICIT NONE
c     decay common
      DOUBLE PRECISION CBR
      INTEGER KDEC,LBARP,IDB
      COMMON /S_CSYDEC/ CBR(223+16+12+8), KDEC(1338+6*(16+12+8)),
     &     LBARP(99), IDB(99)

      INTEGER Imode,IDBST,Iset,I
      DIMENSION IDBST(99)      
      SAVE
      DATA Iset/0/      
      IF(IMODE.EQ.1)THEN
         if(Iset.eq.0) then
            print*,'no decay configuration to restore known!',
     &           ' run dec_set with -1 first!'
            stop
         endif
         DO I=1,99
            IDB(I)=IDBST(i)
         ENDDO
      ELSEIF(IMODE.EQ.-1)THEN
         DO I=1,99
            IDBST(I)=IDB(i)
         ENDDO
         Iset = 1
      ELSE
         PRINT*,'WRONG USE OF DEC_SET!'
         STOP
      ENDIF         
      END
C=======================================================================
      
      SUBROUTINE MUON_ENHANCEMENTS 
C-----------------------------------------------------------------------
C...  routine that enacts different enhancements of muons in eas by
c     modifying the final state
c      
c     Parameters:
c      
c     model: IPAR(94)
c      
c       -1: none but execute decays in sibyll event generation
c        0: none/do not call
c        1: rho0
c        2: strangeness (single threshold)
c        20: strangeness (dual threshold, allows different enhancement at LE and HE )
c        3: baryons (single threshold)
c        30: baryons (dual threshold, allows different enhancement at LE and HE )
c        4: mixed enhancement: rho0 & baryons
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
C     The final particle output is contained in COMMON /S_PLIST/    
C     NP           : number of final particles
C     P(1:NP, 1:5) : 4-momenta + masses of the final particles 
C     LLIST (1:NP) : codes of final particles
      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      SAVE

c     store initial decay configuration
      call dec_set(-1)      
c     set decay configuration for enhancements
      call dec_ini_ext
c     run decay routine
      call decsib
c     reset decay configuration
      call dec_set(1)
c     post-process
      xchgRate = PAR(75)        ! rate
      iedep = IPAR(96)          ! energy dependence (0: none, 1: log)
      eThr = PAR(160)           ! energy threshold
      xfexp = PAR(159)          ! xf-weighting
      iproj = IPAR(97)          ! projectiles (0: all, 1: mesons only)
      if(IPAR(94).eq.1)then      
         CALL FORCE_VECTORS1(xchgRate,iedep,ethr,xfexp,iproj,1,NP)
         
      elseif(IPAR(94).eq.2)then
         call FORCE_STRANGE(xchgRate,iedep,ethr,xfexp,iproj,1,NP)

      elseif(IPAR(94).eq.20)then
c     low-energy enhancement
         call FORCE_STRANGE(xchgRate,iedep,ethr,xfexp,iproj,1,NP)
c     high-energy enhancement         
         xchgRateHE = PAR(161)
         ethrHE = PAR(162)
         xfHE = PAR(163)
         call FORCE_STRANGE(xchgRateHE,iedep,ethrHE,xfHE,iproj,1,NP)

      elseif(IPAR(94).eq.3)then         
         call FORCE_BARYONS(xchgRate,iedep,ethr,xfexp,iproj,1,NP)

      elseif(IPAR(94).eq.30)then
c     low-energy enhancement
         call FORCE_BARYONS(xchgRate,iedep,ethr,xfexp,iproj,1,NP)
c     high-energy enhancement         
         xchgRateHE = PAR(161)
         ethrHE = PAR(162)
         xfHE = PAR(163)
         call FORCE_BARYONS(xchgRateHE,iedep,ethrHE,xfHE,iproj,1,NP)

      elseif(IPAR(94).eq.4)then ! mixed enhancement rho0 & baryons
c     parameters of vector enhancement are default ones, set by vector_ini
c     (rate: par75, edep: ipar96, Ethr: par160, xfexponent: par159)
         CALL FORCE_VECTORS1(xchgRate,iedep,ethr,xfexp,iproj,1,NP)
c     parameters of baryon enhancement are rate: par161, Ethr: par162, xfexponent: par163, projectile: ipar98
c     energy dependence is same as vector enhancement
         xRateBAR = PAR(161)
         ethrBAR = PAR(162)
         xfBAR = PAR(163)
         iprojBAR = IPAR(98)
         call FORCE_BARYONS(xRateBAR,iedep,ethrBAR,xfBAR,iprojBAR,1,NP)
         
      endif
      END
      
C=======================================================================

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

c     external types
      DOUBLE PRECISION SQS,SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO      
      DIMENSION SIGDIF(3)
      INTEGER L0

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      INTEGER NS_max, NH_max
      PARAMETER (NS_max = 20, NH_max = 80)
      
      DOUBLE PRECISION SSIG,PJETC,SSIGN,SSIGNSD,SSIGNEL,ALINT,ASQSMIN,
     &     ASQSMAX,DASQS
      INTEGER NSQS
      COMMON /S_CCSIG/ SSIG(61,3), PJETC(0:NS_max,0:NH_max,61,2),
     &     SSIGN(61,3,3), SSIGNSD(61,3,3), SSIGNEL(61,3,3), 
     &     ALINT(61,3,3), ASQSMIN, ASQSMAX, DASQS, NSQS
      DOUBLE PRECISION SSIG_TOT,SSIG_SD1,SSIG_SD2,SSIG_DD,SSIG_B,
     &     SSIG_RHO
      COMMON /S_CCSIG2/ SSIG_TOT(61,3),SSIG_SD1(61,3),SSIG_SD2(61,3),
     &    SSIG_DD(61,3),SSIG_B(61,3),SSIG_RHO(61,3)

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN

c     internal type declarations
      DOUBLE PRECISION T,AL
      INTEGER LL,L,J1
      DIMENSION LL(39)
      SAVE
      DATA LL /5*0,3*2,4*3,2*1,19*0,6*1/


      L = L0
      IF(NSQS.LE.0) THEN
        WRITE(LUN,'(//,1X,A)') 
     &    ' SIB_SIGMA_HP: interpolation table not initialized.'
        STOP
      ENDIF
      IF(IABS(L).gt.39)THEN
         WRITE(LUN,*)     
     &        ' SIB_SIGMA_HP: unknown beam particle!',L
         STOP
      ENDIF
      IF(L.GT.3) L=LL(IABS(L))
      IF(L.EQ.0)THEN
         WRITE(LUN,*)     
     &        ' SIB_SIGMA_HP: unknown beam particle!', L
         STOP
      ENDIF
        
      AL = LOG10(SQS)
      J1 = INT((AL-1.D0)*10.D0 + 1)
      if((j1.lt.1).or.(j1.gt.NSQS)) then
        if(ndebug.gt.0)
     &         write (LUN,'(1x,a,i3,1p,e12.3)') 
     &         ' SIB_SIGMA_HP: energy out of range ',L,sqs
      endif
      if((j1.lt.1).or.(j1.ge.NSQS)) then
        J1 = min(J1,NSQS-1)
        J1 = max(J1,1)
      endif
      T = (AL-1.D0)*10.D0 - DBLE(J1-1)
      SIGT    = SSIG_TOT(J1,L)*(1.D0-T) + SSIG_TOT(J1+1,L)*T
      SIGINEL = SSIG(J1,L)*(1.D0-T) + SSIG(J1+1,L)*T
      SIGEL   = SIGT-SIGINEL
      SIGDIF(1) = SSIG_SD1(J1,L)*(1.D0-T) + SSIG_SD1(J1+1,L)*T
      SIGDIF(2) = SSIG_SD2(J1,L)*(1.D0-T) + SSIG_SD2(J1+1,L)*T
      SIGDIF(3) = SSIG_DD(J1,L)*(1.D0-T) + SSIG_DD(J1+1,L)*T
      SLOPE   = SSIG_B(J1,L) *(1.D0-T) + SSIG_B(J1+1,L)*T
      RHO     = SSIG_RHO(J1,L) *(1.D0-T) + SSIG_RHO(J1+1,L)*T

      END

C=======================================================================

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
c     external types      
      DOUBLE PRECISION SQS,SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO
      DIMENSION SIGDIF(3,2)
      INTEGER L

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      INTEGER NS_max, NH_max
      PARAMETER (NS_max = 20, NH_max = 80)
      
      DOUBLE PRECISION SSIG,PJETC,SSIGN,SSIGNSD,SSIGNEL,ALINT,ASQSMIN,
     &     ASQSMAX,DASQS
      INTEGER NSQS
      COMMON /S_CCSIG/ SSIG(61,3), PJETC(0:NS_max,0:NH_max,61,2),
     &     SSIGN(61,3,3), SSIGNSD(61,3,3), SSIGNEL(61,3,3), 
     &     ALINT(61,3,3), ASQSMIN, ASQSMAX, DASQS, NSQS
      DOUBLE PRECISION SSIG_TOT,SSIG_SD1,SSIG_SD2,SSIG_DD,SSIG_B,
     &     SSIG_RHO
      COMMON /S_CCSIG2/ SSIG_TOT(61,3),SSIG_SD1(61,3),SSIG_SD2(61,3),
     &    SSIG_DD(61,3),SSIG_B(61,3),SSIG_RHO(61,3)
      DOUBLE PRECISION SSIG_SD1LM,SSIG_SD1HM,SSIG_SD2LM,SSIG_SD2HM,
     &     SSIG_DDLM,SSIG_DDHM
      COMMON /S_CCSIG3/ SSIG_SD1LM(61,3),SSIG_SD1HM(61,3),
     &     SSIG_SD2LM(61,3),SSIG_SD2HM(61,3),
     &     SSIG_DDLM(61,3),SSIG_DDHM(61,3)

c     internal types
      INTEGER J1
      DOUBLE PRECISION T,AL
      SAVE

      IF(NSQS.LE.0) THEN
        WRITE(LUN,'(//,1X,A)') 
     &    ' SIB_SIGMA_HP2: interpolation table not initialized.'
        STOP
      ENDIF
        
      AL = dLOG10(SQS)
      J1 = INT((AL - 1.D0)*10.D0 + 1)
      if((j1.lt.1).or.(j1.gt.NSQS)) then
        if(ndebug.gt.0)write(LUN,'(1x,a,i3,1p,e12.3)') 
     &    ' SIB_SIGMA_HP2: energy out of range ',L,sqs
      endif
      if((j1.lt.1).or.(j1.ge.NSQS)) then
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

C=======================================================================

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

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      INTEGER NS_max, NH_max
      PARAMETER (NS_max = 20, NH_max = 80)
      
      DOUBLE PRECISION SSIG,PJETC,SSIGN,SSIGNSD,SSIGNEL,ALINT,ASQSMIN,
     &     ASQSMAX,DASQS
      INTEGER NSQS
      COMMON /S_CCSIG/ SSIG(61,3), PJETC(0:NS_max,0:NH_max,61,2),
     &     SSIGN(61,3,3), SSIGNSD(61,3,3), SSIGNEL(61,3,3), 
     &     ALINT(61,3,3), ASQSMIN, ASQSMAX, DASQS, NSQS

c     external
      DOUBLE PRECISION SQS,SIGPROD,SIGBDIF
      INTEGER L

c     internal
      DOUBLE PRECISION AL,T
      INTEGER J1
      SAVE

      IF(NSQS.LE.0) THEN
        WRITE(LUN,'(//,1X,A)') 
     &    ' SIB_SIGMA_HAIR: interpolation table not initialized.'
        STOP
      ENDIF
        
      AL = LOG10(SQS)
      J1 = INT((AL - 1.D0)*10.D0 + 1)
      if((j1.lt.1).or.(j1.gt.NSQS)) then
        if (ndebug .gt. 0) 
     &          write (LUN,'(1x,a,i3,1p,e12.3)') 
     &          ' SIB_SIGMA_HAIR: energy out of range ',L,sqs
      endif
      if((j1.lt.1).or.(j1.ge.NSQS)) then
        J1 = min(J1,NSQS-1)
        J1 = max(J1,1)
      endif
      T = (AL-1.D0)*10.D0 - DBLE(J1-1)
      SIGprod = SSIGN(J1,L,1)*(1.D0-T) + SSIGN(J1+1,L,1)*T
      SIGbdif = SSIGNSD(J1,L,1)*(1.D0-T) + SSIGNSD(J1+1,L,1)*T

      END
C=======================================================================
      SUBROUTINE SIB_SIGMA_HNUC (L,IAT,SQS,SIGprod,SIGbdif,SIGela) 

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
C                  SIGela     elastic cross section
C-----------------------------------------------------------------------
Cf2py integer, intent(in) :: L,IAT
Cf2py double precision, intent(in) :: SQS
Cf2py double precision, intent(out) :: SIGprod,SIGbdif,SIGela
      IMPLICIT NONE

      INTEGER NS_max, NH_max
      PARAMETER (NS_max = 20, NH_max = 80)
      
      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      DOUBLE PRECISION SSIG,PJETC,SSIGN,SSIGNSD,SSIGNEL,ALINT,ASQSMIN,
     &     ASQSMAX,DASQS
      INTEGER NSQS
      COMMON /S_CCSIG/ SSIG(61,3), PJETC(0:NS_max,0:NH_max,61,2),
     &     SSIGN(61,3,3), SSIGNSD(61,3,3), SSIGNEL(61,3,3), 
     &     ALINT(61,3,3), ASQSMIN, ASQSMAX, DASQS, NSQS

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      DOUBLE PRECISION SIGT,SIGEL,SIGINEL,SIGQE,SIGSD,SIGQSD,SIGPPT,
     &     SIGPPEL,SIGPPSD
      INTEGER ITG
      COMMON /NUCSIG/ SIGT,SIGEL,SIGINEL,SIGQE,SIGSD,
     +     SIGQSD,SIGPPT,SIGPPEL,SIGPPSD,ITG

c     external
      DOUBLE PRECISION SQS,SIGPROD,SIGBDIF,SIGELA
      INTEGER L,IAT

c     internal
      DOUBLE PRECISION ALAM,T,AL
      INTEGER IPARM,ICSMOD,IK,J1
      SAVE
      IF(NSQS.LE.0) THEN
        WRITE(LUN,'(//,1X,A)') 
     &    ' SIB_SIGMA_HNUC: interpolation table not initialized.'
        STOP
      ENDIF
      IF(IAT.eq.0.or.IAT.eq.14.or.IAT.eq.16)THEN
c     read cross section from table
         IF(IAT.eq.0) THEN
            IK=1
         ELSEIF(IAT.eq.14)THEN
            IK=2
         ELSE
            IK=3
         ENDIF            
         AL = LOG10(SQS)
         J1 = INT((AL - 1.D0)*10.D0 + 1)
         if((j1.lt.1).or.(j1.gt.NSQS)) then
            if (ndebug .gt. 0) 
     &           write (LUN,'(1x,a,i3,1p,e12.3)') 
     &           ' SIB_SIGMA_HNUC: energy out of range ',L,sqs
         endif
         if((j1.lt.1).or.(j1.ge.NSQS)) then
            J1 = min(J1,NSQS-1)
            J1 = max(J1,1)
         endif
         T = (AL-1.D0)*10.D0 - DBLE(J1-1)
         SIGprod = SSIGN(J1,L,IK)*(1.D0-T) + SSIGN(J1+1,L,IK)*T
         SIGbdif = SSIGNSD(J1,L,IK)*(1.D0-T) + SSIGNSD(J1+1,L,IK)*T
         SIGela  = SSIGNEL(J1,L,IK)*(1.D0-T) + SSIGNEL(J1+1,L,IK)*T
      ELSEIF(IAT.lt.19)THEN
c     calculate cross section         
         IF(ndebug.gt.0)THEN
            WRITE(LUN,'(1X,A,2I4,F8.2)')
     &           'SIB_SIGMA_HNUC: L,IAT,SQS:',L,IAT,SQS
         ENDIF
c     calculate hadron - nucleus cross section
c     dummy arg, coupling derived from dif xsctn
         ALAM = 1.D0              
c     use Sibyll p-p cross section as input
         ICSMOD = 1             
c     use Goulianos param. for inel. coupling param.
         IPARM = 2 
         CALL SIG_HAD_NUC(L,IAT,SQS,ALAM,ICSMOD,IPARM)
C     particle production cross section        
         SIGprod = SIGT-SIGQE
C     quasi elastic + elastic singl. diff cross section
         SIGbdif = SIGQSD
c     elastic cross section
         SIGela = SIGel
         if(ndebug.gt.0)THEN
            WRITE(LUN,'(1X,A,3F8.2)')
     &           'SIB_SIGMA_HNUC: SIGprod, SIGbdif, ALAM:',
     &           SIGprod, SIGbdif, ALAM
         ENDIF
      ELSE
         WRITE(LUN,'(//,1X,A)') 
     &     ' SIB_SIGMA_HNUC: number of target nucleons too large!',
     &     ' (0<=IAT<=18)'
         SIGprod = -1.D0
      ENDIF
      RETURN
      END

C----------------------------------------------------------------------
C     sampling routines for hard partons in SIBYLL
C     includes GRV98 pdf table and initialization routines
C----------------------------------------------------------------------
C=======================================================================

      SUBROUTINE SAMPLE_HARD (L,X1,X2,PT)

C-----------------------------------------------------------------------
C...Routine for sampling the kinematical variables 
C.  that determine a jet-jet (gluon - gluon) system (x1,x2, pT) 
C.  from the differential cross section:
C.     d3sigma/(dx1 dx2 dpT)
C.  This version assumes the `single parton approximation'
C.  INPUT:  L=1 incident proton, L=2  incident pi
C.          NPLD: position on parton stack
C.  OUTPUT:  gluon 4momenta
C-----------------------------------------------------------------------
      IMPLICIT NONE

c     external types
      INTEGER L
      DOUBLE PRECISION X1,X2,PT
      
      INTEGER NW_max
      PARAMETER (NW_max = 20)
C--------------------------------------------------------------------
C     SIBYLL common blocks containing event information       \FR'14
C--------------------------------------------------------------------

C     EVENT INFO COMMON
C     contains overall interaction properties, like
C     SQS : center-of-mass energy
C     S   :         "       "     squared
C     PTmin : low pt cut of QCD cross section, 
C             i.e. minimal pt of hard minijets
C     Xmin : low-x bound for PDFs, 
C            i.e. minimal momentum fraction of hard partons
C     Zmin : logarithm of that
C     KB : PID of beam hadron
C     KT() : PID of target
C     IAT : mass number of target
      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
      INTEGER KB,IAT,KT
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT
      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
c     internal types
      DOUBLE PRECISION Z1,Z2,SIG,S_RNDM,Q2,ZSAMPLE      
      SAVE
 
      IF(ndebug.gt.2)then
         write(lun,*) ' SAMPLE_HARD: (SQS,S,PTmin,Xmin,Zmin)',
     &        SQS,S,PTmin,Xmin,Zmin
      endif

 100  Z1=ZSAMPLE (ZMIN,L)       ! beam L=1,2 for proton or pion
      Z2=ZSAMPLE (ZMIN,1)       ! target always a nucleon
      SIG=1.D0-XMIN*dEXP(-Z1-Z2)
      IF (SIG .LT. S_RNDM(0))  GOTO 100      
      X1=dEXP(Z1)
      X2=dEXP(Z2)
      IF (X1.gt.0.9D0.or.X2.gt.0.9D0) GOTO 100
      Q2=PTmin**2/(1.D0-S_RNDM(L)*SIG)
      IF(Q2.gt.S*X1*X2) goto 100
      PT=dSQRT(Q2*(1.D0-Q2/(S*X1*X2)))

      IF(ndebug.gt.2)then
         write(lun,*) ' SAMPLE_HARD: (X1,X2,PT)',
     &        X1,X2,PT
      endif

      RETURN
      END

C=======================================================================

      FUNCTION ZSAMPLE (ZMIN,L)

C-----------------------------------------------------------------------
C...This function returns as output a value z=log(x)
C.  distributed as f(x) = g(x) + 4/9 *(q(x) + qbar(x))
C.  from a minimum value ZMIN to 0,
C.  for a proton (L=1) or a pi (L=2)
C.  needs to be initialised with: CALL ZSAMPLE_INI
C.....................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      COMMON /S_CZGEN/ XA(2),XB(2),XMAX,ZA(2),ZB(2),ZMAX,
     +     DX(2),DZ(2),APART(2),FFA(2),FFB(2),
     +     DFX(2),DFZ(2),XX(200,2),ZZ(200,2),FFX(200,2),FFZ(200,2),
     +     NX,NZ
      PARAMETER (b=0.268D0)
      PARAMETER (bpi=3.7D0)
      PARAMETER (cpi=0.698D0)
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10
      SAVE

      F = PART_INT(ZMIN,L)*S_RNDM(0)
      IF (F .GE. FFA(L))  THEN
         IF(IPAR(8).EQ.0)THEN
            ZSAMPLE = ZA(L) - (F-FFA(L))/APART(L)
         ELSE
            if(L.eq.1) then
               ZSAMPLE = -1.D0/b * dLOG( 1.D0 - F / APART(L) ) 
            else
               ZSAMPLE = -1.D0 * ( (F - cpi)/APART(L) )**(1.D0/bpi)
            endif
         ENDIF
      ELSE IF (F .GE. FFB(L))  THEN
         JF = INT((F-FFB(L))/DFZ(L) + 1.D0)
         JF = min(JF,199)
         F0 = FFB(L) + DFZ(L)*DBLE(JF-1)
         T = (F-F0)/DFZ(L)
         ZSAMPLE = ZZ(JF,L)*(1.D0-T)+ZZ(JF+1,L)*T
      ELSE
         JF = INT(F/DFX(L)+1.D0)
         JF = min(JF,199)
         F0 = DFX(L)*DBLE(JF-1)
         T = (F-F0)/DFX(L)
         X = XX(JF,L)*(1.D0-T)+XX(JF+1,L)*T
         ZSAMPLE = dLOG(X)
      ENDIF

      RETURN
      END

C=======================================================================

      FUNCTION PART_INT (ZMIN,L)

C-----------------------------------------------------------------------
C...This function returns as output the integral of
C.  the parton structure function:
C.     f(x) = g(x) + 4/9 *(q(x) + qbar(x))
C.  from xmin = exp(zmin) to 1 
C.  for a proton (L=1) or a pi (L=2)
C.  needs to be initialised with: CALL ZSAMPLE_INI
C.....................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      COMMON /S_CZGEN/ XA(2),XB(2),XMAX,ZA(2),ZB(2),ZMAX,
     +     DX(2),DZ(2),APART(2),FFA(2),FFB(2),
     +     DFX(2),DFZ(2),XX(200,2),ZZ(200,2),FFX(200,2),FFZ(200,2),
     +     NX,NZ
      DOUBLE PRECISION b,bpi,cpi
      PARAMETER (b=0.268D0)
      PARAMETER (bpi=3.7D0)
      PARAMETER (cpi=0.698D0)
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10
      SAVE

      IF (ZMIN .LT. ZA(L))  THEN
         IF(IPAR(8).EQ.0)THEN
            PART_INT = FFA(L) + APART(L) * (ZA(L) - ZMIN)
         ELSE
            if(L.eq.1) then
               PART_INT = APART(L) * ( 1.D0 - dEXP(-b*ZMIN) ) 
            else
               PART_INT = APART(L) * ( -ZMIN )**bpi + cpi
            endif
         ENDIF
      ELSE IF (ZMIN .LT. ZB(L)) THEN
         JZ = INT((ZB(L)-ZMIN)/DZ(L)+1.D0)
         JZ = min(JZ,199)
         Z0 = ZB(L)-DZ(L)*DBLE(JZ-1)
         T = (Z0-ZMIN)/DZ(L)
         PART_INT = FFZ(JZ,L)*(1.D0-T) + FFZ(JZ+1,L)*T

      ELSE
         X = EXP(ZMIN)
         JX = INT((XMAX-X)/DX(L)+1.D0)
         JX = min(JX,199)
         X0 = XMAX-DX(L)*DBLE(JX-1)
         T = (X0-X)/DX(L)
         PART_INT = FFX(JX,L)*(1.D0-T) + FFX(JX+1,L)*T
      
      ENDIF
      RETURN
      END

C=======================================================================

      SUBROUTINE GRV_INI

C-----------------------------------------------------------------------
C...This subroutine initializes the COMMON block
C   used for sampling z, according to the GRV98LO
C   pdf set
C..................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      COMMON /S_CZGEN/ XA(2),XB(2),XMAX,ZA(2),ZB(2),ZMAX,
     +     DX(2),DZ(2),APART(2),FFA(2),FFB(2),
     +     DFX(2),DFZ(2),XX(200,2),ZZ(200,2),FFX(200,2),FFZ(200,2),
     +     NX,NZ
      DOUBLE PRECISION b,bpi,cpi
      PARAMETER (b=0.268D0)
      PARAMETER (bpi=3.7D0)
      PARAMETER (cpi=0.698D0)
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10
      SAVE

      IPAR(8) = 1

      XA(1) = 1.D-06
      XB(1) = 1.D-01

      XA(2) = 1.D-04
      XB(2) = 1.D-01

      XMAX = 0.8D0
      ZMAX = dLOG(XMAX)
      NX = 200
      NZ = 200

      DO L=1,2

         ZA(L) = dLOG(XA(L))
         ZB(L) = dLOG(XB(L))
         DX(L) = (XMAX-XB(L))/DBLE(NX)
         DZ(L) = (ZB(L)-ZA(L))/DBLE(NZ)

C     large x: interpolation in x
         FFX(1,L) = 0.D0
         DO J=2,NX
            X = XMAX - DX(L)*(DBLE(J)-1.D0)
            G = PARTON(X,L)/X
            FFX(J,L) = FFX(J-1,L)+G*DX(L)
         ENDDO
         CALL INVERT_ARRAY (FFX(1,L),XMAX,-DX(L),NX,XX(1,L),FMIN,DFX(L))

C     small x: interpolation in log(x)
         FFZ(1,L) = FFX(NX,L)
         DO J=2,NZ
            Z = ZB(L) - DZ(L)*(DBLE(J)-1.D0)
            X = dEXP(Z)
            G = PARTON(X,L)
            FFZ(J,L) = FFZ(J-1,L)+G*DZ(L)
         ENDDO
         CALL INVERT_ARRAY(FFZ(1,L),ZB(L),-DZ(L),NZ,ZZ(1,L),FMIN,DFZ(L))
         FFA(L) = FFZ(NZ,L)
         FFB(L) = FFX(NX,L)
         
C     very small x:  f(x) = A/x**b b=1.268
         IF(L.eq.1) THEN
            APART(L) = FFA(L) / ( 1.D0 - dEXP(-b*ZA(L)) )
         ELSE
            APART(L) = ( FFA(L) - cpi ) / ( -ZA(L) )**bpi
         ENDIF
      ENDDO

      RETURN
      END

C=======================================================================

      SUBROUTINE ZSAMPLE_INI

C-----------------------------------------------------------------------
C...This subroutine initialise the generation of
C.  z = log(x)  for the generation  of z according
C.  to the structure functions
C..................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      COMMON /S_CZGEN/ XA(2),XB(2),XMAX,ZA(2),ZB(2),ZMAX,
     +     DX(2),DZ(2),APART(2),FFA(2),FFB(2),
     +     DFX(2),DFZ(2),XX(200,2),ZZ(200,2),FFX(200,2),FFZ(200,2),
     +     NX,NZ
      SAVE

      IPAR(8) = 0

      XA(1) = 1.D-04
      XB(1) = 1.D-01
      XMAX = 0.8D0
      ZA(1) = dLOG(XA(1))
      ZB(1) = dLOG(XB(1))
      ZMAX = dLOG(XMAX)
      NX = 200
      NZ = 200
      DX(1) = (XMAX-XB(1))/DBLE(NX-1)
      DZ(1) = (ZB(1)-ZA(1))/DBLE(NZ-1)

      XA(2) = 1.D-04
      XB(2) = 1.D-01
      XMAX = 0.8D0
      ZA(2) = dLOG(XA(2))
      ZB(2) = dLOG(XB(2))
      ZMAX = dLOG(XMAX)
      NX = 200
      NZ = 200
      DX(2) = (XMAX-XB(2))/DBLE(NX-1)
      DZ(2) = (ZB(2)-ZA(2))/DBLE(NZ-1)

      DO L=1,2
            
C     very small x:  f(x) = A/x
         APART(L) = PARTON(0.D0,L)
            
C         large x: interpolation in x
         FFX(1,L) = 0.D0
         DO J=2,NX
            X = XMAX - DX(L)*(DBLE(J)-0.5D0)
            G = PARTON(X,L)/X
            FFX(J,L) = FFX(J-1,L)+G*DX(L)
         ENDDO
         CALL INVERT_ARRAY (FFX(1,L),XMAX,-DX(L),NX,XX(1,L),FMIN,DFX(L))
            
C     small x: interpolation in log(x)
         FFZ(1,L) = FFX(NX,L)
         DO J=2,NZ
            Z = ZB(L) - DZ(L)*(DBLE(J)-0.5D0)
            X = dEXP(Z)
            G = PARTON(X,L)
            FFZ(J,L) = FFZ(J-1,L)+G*DZ(L)
         ENDDO
         CALL INVERT_ARRAY(FFZ(1,L),ZB(L),-DZ(L),NZ,ZZ(1,L),FMIN,DFZ(L))
         FFA(L) = FFZ(NZ,L)
         FFB(L) = FFX(NX,L)
         
      ENDDO
      RETURN
      END

C=======================================================================

      FUNCTION PARTON(X,L)

C-----------------------------------------------------------------------
C...This function returns the structure function
C.   f(x) = x * [ g(x) + 4/9 *(q(x) + qbar(x)) ]
C.  for a proton. 
C................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      PARAMETER (beta=1.925978D0)
      SAVE
      DATA NOUTP /0/

c     effective scale
      Q2inp = PAR(22)
      IF (L .EQ. 2)  GOTO 1000

      IF(IPAR(8).eq.0) THEN
C...  Eichten et al.  (set 1)
c     tp060203 100      uv = 1.78 * x**0.5 * (1.-x**1.51)**3.5
         uv = 1.78D0 * x**0.5D0 * (1.D0-x**1.51D0)**3.5D0
         dv = 0.67D0 * x**0.4D0 * (1.D0-x**1.51D0)**4.5D0
         us = 0.182D0 * (1.D0-x)**8.54D0
         ss = 0.081D0 * (1.D0-x)**8.54D0
         qq0 = uv + dv + 4.D0*us + 2.D0*ss
         glu0 = (2.62D0 + 9.17D0*x)* (1.D0-x)**5.9D0
      ELSE
         IF( NOUTP.eq.0 ) print *,' using GRV pdf set'
         IF( NOUTP.eq.0 ) print *,' Q2 scale in pdf:',Q2INP
         NOUTP = 1

         CALL SIB_DOR98LO (X, Q2inp, UV, DV, US, DS, SS, GL)
         qq0 = uv + dv + 4.D0* (us + ds) + 2.D0*ss
         glu0 = gl
      ENDIF
      parton = glu0 + 4.D0/9.D0*qq0
      RETURN

 1000 CONTINUE
      IF(IPAR(8).eq.0) THEN
C...Owens set 1   from STRF from Wisc. Pheno. group. for q2=q2_min
         AV=0.4D0
         BV=0.7D0
c      BETA=GGAMMA(AV)*GGAMMA(BV+1.)/GGAMMA(AV+BV+1.)  =1.925978
         uv=X**(AV)*(1.D0-X)**BV/BETA
         dv=uv

         A=0.9D0
         BET=5.D0
         us=(A*(1.D0-X)**BET)/6.D0
         
         A=0.888D0
         BET=3.11D0
         GA1=6.D0
         glu0=A*(1.D0-X)**BET*(1.D0+GA1*X)
c   Bug Fix thanks to Sue Kashahara- correct factor in front of 
c   sea quarks for Owens S.F.  5-94
         qq0 = uv + dv + 6.D0*us
         parton = (glu0 + 4.D0/9.D0*qq0)
         RETURN
      ELSE

c     duv = valence quark distribution
c     dgl = gluon distribution
c     dus = sea quark distribution (u,d,s)
c     dds = sea charm quark ( neglected )
c     dss = sea bottom quark ( neglected )
         CALL DORPLO (X, Q2inp, uv, gl, us, ds, ss)
         qq0 = uv + dv + 4.D0*us
         glu0 = gl
         parton = (glu0 + 4.D0/9.D0*qq0)
         RETURN
      ENDIF
      END
C=======================================================================

      SUBROUTINE PDF_INI

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      COMMON /S_CZGEN/ XA(2),XB(2),XMAX,ZA(2),ZB(2),ZMAX,
     +     DX(2),DZ(2),APART(2),FFA(2),FFB(2),
     +     DFX(2),DFZ(2),XX(200,2),ZZ(200,2),FFX(200,2),FFZ(200,2),
     +     NX,NZ
      SAVE

      IF(IPAR(8).eq.0) THEN
         if (ndebug .gt. 0 ) WRITE(LUN,*)
     *       ' PDF_INI: calcuLating pdf table using Eichten param..'
         CALL ZSAMPLE_INI
      ELSEIF(IPAR(8).eq.2) THEN
         if (ndebug .gt. 0 ) then
           WRITE(LUN,*)' PDF_INI: calculating pdf table using GRV',
     *                  '  param..'
           WRITE(LUN,*)' does not work with -fbounds-check !!'
         endif
         CALL GRV_INI
      ELSE
         if (ndebug .gt. 0 ) WRITE(LUN,*)
     *        ' PDF_INI: using common table of GRV parametrization..'
      ENDIF
      if (ndebug .gt. 0 )  THEN
           WRITE(LUN,*)APART(1),FFA(1),FFB(1),DX(1),DZ(1)
           WRITE(LUN,*)APART(2),FFA(2),FFB(2),DX(2),DZ(2)
      ENDIF
      END

C=======================================================================

      BLOCK DATA PDFINI

C-----------------------------------------------------------------------
C..   tabled parton distribution function
c     Proton: GRV98LO , Eur.Phys.J. C5(1998) 461-470
c     Pion:   GRV91 , Z. Phys. C53, 651-655 (1992)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      IMPLICIT INTEGER(I-N)
      COMMON /S_CZGEN/ XA(2),XB(2),XMAX,ZA(2),ZB(2),ZMAX,
     +     DX(2),DZ(2),APART(2),FFA(2),FFB(2),
     +     DFX(2),DFZ(2),XX(200,2),ZZ(200,2),FFX(200,2),FFZ(200,2),
     +     NX,NZ
      SAVE
      DATA XA /1.D-06,0.0001D0/
      DATA XB /0.1D0,0.1D0/
      DATA XMAX /0.800000011921D0/
      DATA ZMAX /-0.223143532872D0/
      DATA NX /200/
      DATA NZ /200/
      DATA ZA /-13.8155D0,-9.21034D0/
      DATA ZB /-2.30259D0,-2.30259D0/
      DATA DX /0.00351759D0,0.00351759D0/
      DATA DZ /0.0578539D0,0.0347123D0/
      DATA DFX /0.00952501D0,0.00847474D0/
      DATA DFZ /1.93863D0,0.326082D0/
      DATA APART /-9.80215D0,0.0178207D0/
      DATA FFA /387.684D0,66.5767D0/
      DATA FFB /1.89548D0,1.68647D0/
      
      DATA (FFX(K,1),K=1,200 ) /
     &0.000D+00,6.380D-05,1.315D-04,2.034D-04,2.795D-04,
     &3.601D-04,4.454D-04,5.356D-04,6.309D-04,7.315D-04,
     &8.377D-04,9.497D-04,1.068D-03,1.192D-03,1.323D-03,
     &1.460D-03,1.605D-03,1.756D-03,1.916D-03,2.083D-03,
     &2.258D-03,2.441D-03,2.633D-03,2.835D-03,3.045D-03,
     &3.265D-03,3.496D-03,3.736D-03,3.988D-03,4.250D-03,
     &4.524D-03,4.810D-03,5.108D-03,5.418D-03,5.742D-03,
     &6.078D-03,6.429D-03,6.794D-03,7.174D-03,7.570D-03,
     &7.981D-03,8.408D-03,8.852D-03,9.313D-03,9.793D-03,
     &1.029D-02,1.081D-02,1.134D-02,1.190D-02,1.247D-02,
     &1.307D-02,1.369D-02,1.433D-02,1.500D-02,1.568D-02,
     &1.640D-02,1.714D-02,1.790D-02,1.869D-02,1.951D-02,
     &2.035D-02,2.123D-02,2.213D-02,2.307D-02,2.403D-02,
     &2.503D-02,2.607D-02,2.713D-02,2.823D-02,2.937D-02,
     &3.054D-02,3.176D-02,3.301D-02,3.430D-02,3.563D-02,
     &3.701D-02,3.842D-02,3.989D-02,4.139D-02,4.295D-02,
     &4.455D-02,4.620D-02,4.791D-02,4.966D-02,5.147D-02,
     &5.334D-02,5.526D-02,5.724D-02,5.927D-02,6.137D-02,
     &6.353D-02,6.576D-02,6.805D-02,7.041D-02,7.284D-02,
     &7.534D-02,7.791D-02,8.056D-02,8.329D-02,8.609D-02,
     &8.898D-02,9.195D-02,9.500D-02,9.814D-02,1.014D-01,
     &1.047D-01,1.081D-01,1.116D-01,1.153D-01,1.190D-01,
     &1.228D-01,1.267D-01,1.308D-01,1.350D-01,1.392D-01,
     &1.436D-01,1.481D-01,1.528D-01,1.575D-01,1.624D-01,
     &1.674D-01,1.725D-01,1.778D-01,1.832D-01,1.888D-01,
     &1.946D-01,2.005D-01,2.066D-01,2.128D-01,2.193D-01,
     &2.259D-01,2.327D-01,2.397D-01,2.469D-01,2.543D-01,
     &2.619D-01,2.698D-01,2.778D-01,2.862D-01,2.947D-01,
     &3.035D-01,3.125D-01,3.218D-01,3.314D-01,3.413D-01,
     &3.514D-01,3.618D-01,3.726D-01,3.836D-01,3.950D-01,
     &4.067D-01,4.188D-01,4.312D-01,4.440D-01,4.572D-01,
     &4.708D-01,4.848D-01,4.992D-01,5.141D-01,5.294D-01,
     &5.452D-01,5.615D-01,5.783D-01,5.956D-01,6.134D-01,
     &6.319D-01,6.509D-01,6.706D-01,6.909D-01,7.118D-01,
     &7.334D-01,7.558D-01,7.789D-01,8.029D-01,8.276D-01,
     &8.532D-01,8.797D-01,9.072D-01,9.356D-01,9.650D-01,
     &9.956D-01,1.027D+00,1.060D+00,1.094D+00,1.130D+00,
     &1.167D+00,1.205D+00,1.245D+00,1.287D+00,1.331D+00,
     &1.376D+00,1.423D+00,1.473D+00,1.525D+00,1.579D+00,
     &1.636D+00,1.696D+00,1.759D+00,1.826D+00,1.895D+00/
      
      DATA (FFX(K,2),K=1,200 ) /
     &0.000D+00,7.266D-04,1.470D-03,2.231D-03,3.009D-03,
     &3.805D-03,4.619D-03,5.450D-03,6.300D-03,7.168D-03,
     &8.055D-03,8.961D-03,9.886D-03,1.083D-02,1.179D-02,
     &1.278D-02,1.378D-02,1.481D-02,1.585D-02,1.692D-02,
     &1.800D-02,1.911D-02,2.024D-02,2.139D-02,2.256D-02,
     &2.376D-02,2.498D-02,2.622D-02,2.748D-02,2.877D-02,
     &3.008D-02,3.142D-02,3.278D-02,3.416D-02,3.557D-02,
     &3.701D-02,3.847D-02,3.996D-02,4.147D-02,4.301D-02,
     &4.458D-02,4.617D-02,4.779D-02,4.945D-02,5.112D-02,
     &5.283D-02,5.457D-02,5.634D-02,5.813D-02,5.996D-02,
     &6.182D-02,6.371D-02,6.563D-02,6.759D-02,6.957D-02,
     &7.159D-02,7.365D-02,7.573D-02,7.786D-02,8.001D-02,
     &8.221D-02,8.443D-02,8.670D-02,8.900D-02,9.134D-02,
     &9.372D-02,9.614D-02,9.860D-02,1.011D-01,1.036D-01,
     &1.062D-01,1.088D-01,1.115D-01,1.142D-01,1.170D-01,
     &1.197D-01,1.226D-01,1.255D-01,1.284D-01,1.314D-01,
     &1.344D-01,1.375D-01,1.406D-01,1.438D-01,1.470D-01,
     &1.503D-01,1.536D-01,1.570D-01,1.605D-01,1.640D-01,
     &1.675D-01,1.712D-01,1.748D-01,1.786D-01,1.824D-01,
     &1.862D-01,1.901D-01,1.941D-01,1.982D-01,2.023D-01,
     &2.065D-01,2.107D-01,2.151D-01,2.195D-01,2.239D-01,
     &2.285D-01,2.331D-01,2.378D-01,2.426D-01,2.474D-01,
     &2.524D-01,2.574D-01,2.625D-01,2.677D-01,2.730D-01,
     &2.784D-01,2.839D-01,2.895D-01,2.951D-01,3.009D-01,
     &3.068D-01,3.128D-01,3.189D-01,3.251D-01,3.314D-01,
     &3.378D-01,3.443D-01,3.510D-01,3.578D-01,3.647D-01,
     &3.717D-01,3.789D-01,3.862D-01,3.937D-01,4.012D-01,
     &4.090D-01,4.169D-01,4.249D-01,4.331D-01,4.415D-01,
     &4.500D-01,4.587D-01,4.676D-01,4.767D-01,4.859D-01,
     &4.954D-01,5.050D-01,5.148D-01,5.249D-01,5.352D-01,
     &5.457D-01,5.564D-01,5.674D-01,5.786D-01,5.901D-01,
     &6.019D-01,6.139D-01,6.262D-01,6.388D-01,6.517D-01,
     &6.649D-01,6.785D-01,6.923D-01,7.066D-01,7.212D-01,
     &7.362D-01,7.516D-01,7.673D-01,7.836D-01,8.002D-01,
     &8.174D-01,8.350D-01,8.532D-01,8.718D-01,8.911D-01,
     &9.109D-01,9.313D-01,9.524D-01,9.742D-01,9.966D-01,
     &1.020D+00,1.044D+00,1.069D+00,1.094D+00,1.121D+00,
     &1.149D+00,1.177D+00,1.207D+00,1.238D+00,1.271D+00,
     &1.304D+00,1.339D+00,1.376D+00,1.414D+00,1.454D+00,
     &1.496D+00,1.540D+00,1.586D+00,1.635D+00,1.686D+00/
      
      DATA (FFZ(K,1),K=1,200 ) /
     &1.895D+00,2.014D+00,2.137D+00,2.263D+00,2.393D+00,
     &2.527D+00,2.665D+00,2.807D+00,2.953D+00,3.103D+00,
     &3.257D+00,3.417D+00,3.580D+00,3.748D+00,3.921D+00,
     &4.098D+00,4.281D+00,4.469D+00,4.663D+00,4.861D+00,
     &5.065D+00,5.274D+00,5.489D+00,5.710D+00,5.937D+00,
     &6.170D+00,6.409D+00,6.654D+00,6.906D+00,7.164D+00,
     &7.430D+00,7.702D+00,7.981D+00,8.267D+00,8.561D+00,
     &8.862D+00,9.171D+00,9.487D+00,9.811D+00,1.014D+01,
     &1.048D+01,1.083D+01,1.119D+01,1.156D+01,1.193D+01,
     &1.232D+01,1.271D+01,1.311D+01,1.352D+01,1.395D+01,
     &1.438D+01,1.482D+01,1.527D+01,1.573D+01,1.621D+01,
     &1.669D+01,1.718D+01,1.769D+01,1.821D+01,1.874D+01,
     &1.928D+01,1.983D+01,2.040D+01,2.097D+01,2.156D+01,
     &2.217D+01,2.278D+01,2.341D+01,2.406D+01,2.471D+01,
     &2.539D+01,2.607D+01,2.677D+01,2.749D+01,2.822D+01,
     &2.896D+01,2.973D+01,3.050D+01,3.130D+01,3.211D+01,
     &3.293D+01,3.378D+01,3.464D+01,3.552D+01,3.642D+01,
     &3.733D+01,3.827D+01,3.922D+01,4.020D+01,4.119D+01,
     &4.220D+01,4.323D+01,4.429D+01,4.536D+01,4.646D+01,
     &4.758D+01,4.872D+01,4.988D+01,5.106D+01,5.227D+01,
     &5.350D+01,5.476D+01,5.604D+01,5.735D+01,5.868D+01,
     &6.003D+01,6.142D+01,6.282D+01,6.426D+01,6.572D+01,
     &6.721D+01,6.873D+01,7.028D+01,7.186D+01,7.346D+01,
     &7.510D+01,7.677D+01,7.847D+01,8.020D+01,8.196D+01,
     &8.375D+01,8.558D+01,8.744D+01,8.934D+01,9.127D+01,
     &9.324D+01,9.524D+01,9.728D+01,9.936D+01,1.015D+02,
     &1.036D+02,1.058D+02,1.080D+02,1.103D+02,1.126D+02,
     &1.150D+02,1.174D+02,1.198D+02,1.223D+02,1.248D+02,
     &1.274D+02,1.300D+02,1.327D+02,1.354D+02,1.381D+02,
     &1.409D+02,1.438D+02,1.467D+02,1.496D+02,1.526D+02,
     &1.557D+02,1.588D+02,1.619D+02,1.652D+02,1.684D+02,
     &1.718D+02,1.751D+02,1.786D+02,1.821D+02,1.856D+02,
     &1.892D+02,1.929D+02,1.967D+02,2.005D+02,2.043D+02,
     &2.083D+02,2.122D+02,2.163D+02,2.204D+02,2.246D+02,
     &2.289D+02,2.332D+02,2.376D+02,2.421D+02,2.467D+02,
     &2.513D+02,2.560D+02,2.608D+02,2.656D+02,2.706D+02,
     &2.756D+02,2.807D+02,2.859D+02,2.911D+02,2.965D+02,
     &3.019D+02,3.074D+02,3.130D+02,3.187D+02,3.245D+02,
     &3.304D+02,3.364D+02,3.425D+02,3.486D+02,3.549D+02,
     &3.612D+02,3.677D+02,3.743D+02,3.809D+02,3.877D+02/
      
      DATA (FFZ(K,2),K=1,200 ) /
     &1.686D+00,1.738D+00,1.791D+00,1.844D+00,1.899D+00,
     &1.955D+00,2.011D+00,2.069D+00,2.128D+00,2.188D+00,
     &2.249D+00,2.311D+00,2.374D+00,2.438D+00,2.504D+00,
     &2.570D+00,2.638D+00,2.708D+00,2.778D+00,2.850D+00,
     &2.923D+00,2.997D+00,3.072D+00,3.149D+00,3.228D+00,
     &3.307D+00,3.388D+00,3.471D+00,3.555D+00,3.640D+00,
     &3.727D+00,3.815D+00,3.905D+00,3.997D+00,4.090D+00,
     &4.184D+00,4.281D+00,4.378D+00,4.478D+00,4.579D+00,
     &4.682D+00,4.787D+00,4.893D+00,5.002D+00,5.112D+00,
     &5.224D+00,5.337D+00,5.453D+00,5.571D+00,5.690D+00,
     &5.811D+00,5.935D+00,6.060D+00,6.188D+00,6.317D+00,
     &6.449D+00,6.583D+00,6.719D+00,6.857D+00,6.997D+00,
     &7.139D+00,7.284D+00,7.431D+00,7.580D+00,7.732D+00,
     &7.886D+00,8.042D+00,8.201D+00,8.363D+00,8.526D+00,
     &8.693D+00,8.862D+00,9.033D+00,9.207D+00,9.384D+00,
     &9.563D+00,9.746D+00,9.930D+00,1.012D+01,1.031D+01,
     &1.050D+01,1.070D+01,1.090D+01,1.110D+01,1.130D+01,
     &1.151D+01,1.172D+01,1.194D+01,1.215D+01,1.237D+01,
     &1.260D+01,1.283D+01,1.306D+01,1.329D+01,1.353D+01,
     &1.377D+01,1.401D+01,1.426D+01,1.451D+01,1.476D+01,
     &1.502D+01,1.528D+01,1.554D+01,1.581D+01,1.608D+01,
     &1.636D+01,1.664D+01,1.692D+01,1.721D+01,1.750D+01,
     &1.780D+01,1.810D+01,1.840D+01,1.871D+01,1.902D+01,
     &1.934D+01,1.966D+01,1.998D+01,2.031D+01,2.065D+01,
     &2.098D+01,2.133D+01,2.167D+01,2.203D+01,2.238D+01,
     &2.274D+01,2.311D+01,2.348D+01,2.385D+01,2.423D+01,
     &2.462D+01,2.501D+01,2.541D+01,2.581D+01,2.621D+01,
     &2.662D+01,2.704D+01,2.746D+01,2.789D+01,2.832D+01,
     &2.875D+01,2.920D+01,2.965D+01,3.010D+01,3.056D+01,
     &3.103D+01,3.150D+01,3.198D+01,3.246D+01,3.295D+01,
     &3.344D+01,3.395D+01,3.445D+01,3.497D+01,3.549D+01,
     &3.601D+01,3.655D+01,3.709D+01,3.763D+01,3.819D+01,
     &3.875D+01,3.931D+01,3.989D+01,4.047D+01,4.105D+01,
     &4.165D+01,4.225D+01,4.286D+01,4.347D+01,4.410D+01,
     &4.473D+01,4.537D+01,4.601D+01,4.666D+01,4.732D+01,
     &4.799D+01,4.867D+01,4.935D+01,5.005D+01,5.075D+01,
     &5.146D+01,5.217D+01,5.290D+01,5.363D+01,5.437D+01,
     &5.512D+01,5.588D+01,5.665D+01,5.743D+01,5.821D+01,
     &5.901D+01,5.981D+01,6.062D+01,6.145D+01,6.228D+01,
     &6.312D+01,6.397D+01,6.483D+01,6.570D+01,6.658D+01/
      
      DATA (XX(K,1),K=1,200 ) /
     &8.000D-01,6.472D-01,5.944D-01,5.597D-01,5.335D-01,
     &5.121D-01,4.941D-01,4.785D-01,4.647D-01,4.522D-01,
     &4.409D-01,4.306D-01,4.210D-01,4.122D-01,4.039D-01,
     &3.961D-01,3.887D-01,3.817D-01,3.751D-01,3.688D-01,
     &3.628D-01,3.571D-01,3.516D-01,3.463D-01,3.413D-01,
     &3.365D-01,3.318D-01,3.273D-01,3.230D-01,3.188D-01,
     &3.147D-01,3.108D-01,3.070D-01,3.033D-01,2.998D-01,
     &2.963D-01,2.929D-01,2.896D-01,2.864D-01,2.833D-01,
     &2.802D-01,2.773D-01,2.744D-01,2.715D-01,2.688D-01,
     &2.661D-01,2.634D-01,2.608D-01,2.583D-01,2.558D-01,
     &2.534D-01,2.510D-01,2.487D-01,2.464D-01,2.442D-01,
     &2.420D-01,2.398D-01,2.377D-01,2.356D-01,2.336D-01,
     &2.316D-01,2.296D-01,2.277D-01,2.257D-01,2.239D-01,
     &2.220D-01,2.202D-01,2.184D-01,2.167D-01,2.150D-01,
     &2.132D-01,2.116D-01,2.099D-01,2.083D-01,2.067D-01,
     &2.051D-01,2.036D-01,2.020D-01,2.005D-01,1.990D-01,
     &1.976D-01,1.961D-01,1.947D-01,1.933D-01,1.919D-01,
     &1.905D-01,1.891D-01,1.878D-01,1.865D-01,1.852D-01,
     &1.839D-01,1.826D-01,1.814D-01,1.801D-01,1.789D-01,
     &1.777D-01,1.765D-01,1.753D-01,1.741D-01,1.730D-01,
     &1.718D-01,1.707D-01,1.696D-01,1.685D-01,1.674D-01,
     &1.663D-01,1.653D-01,1.642D-01,1.632D-01,1.622D-01,
     &1.611D-01,1.601D-01,1.591D-01,1.581D-01,1.572D-01,
     &1.562D-01,1.552D-01,1.543D-01,1.534D-01,1.524D-01,
     &1.515D-01,1.506D-01,1.497D-01,1.488D-01,1.479D-01,
     &1.471D-01,1.462D-01,1.453D-01,1.445D-01,1.437D-01,
     &1.428D-01,1.420D-01,1.412D-01,1.404D-01,1.396D-01,
     &1.388D-01,1.380D-01,1.372D-01,1.365D-01,1.357D-01,
     &1.349D-01,1.342D-01,1.335D-01,1.327D-01,1.320D-01,
     &1.313D-01,1.306D-01,1.299D-01,1.292D-01,1.284D-01,
     &1.278D-01,1.271D-01,1.264D-01,1.257D-01,1.251D-01,
     &1.244D-01,1.237D-01,1.231D-01,1.224D-01,1.218D-01,
     &1.212D-01,1.205D-01,1.199D-01,1.193D-01,1.187D-01,
     &1.181D-01,1.175D-01,1.169D-01,1.163D-01,1.157D-01,
     &1.151D-01,1.145D-01,1.139D-01,1.134D-01,1.128D-01,
     &1.123D-01,1.117D-01,1.112D-01,1.106D-01,1.101D-01,
     &1.095D-01,1.090D-01,1.085D-01,1.079D-01,1.074D-01,
     &1.069D-01,1.064D-01,1.059D-01,1.054D-01,1.049D-01,
     &1.044D-01,1.039D-01,1.034D-01,1.029D-01,1.024D-01,
     &1.019D-01,1.014D-01,1.010D-01,1.005D-01,1.000D-01/
      
      DATA (XX(K,2),K=1,200 ) /
     &8.000D-01,7.632D-01,7.331D-01,7.073D-01,6.846D-01,
     &6.643D-01,6.458D-01,6.289D-01,6.132D-01,5.986D-01,
     &5.849D-01,5.721D-01,5.600D-01,5.485D-01,5.376D-01,
     &5.272D-01,5.172D-01,5.077D-01,4.986D-01,4.899D-01,
     &4.815D-01,4.734D-01,4.656D-01,4.581D-01,4.508D-01,
     &4.438D-01,4.370D-01,4.304D-01,4.240D-01,4.178D-01,
     &4.118D-01,4.059D-01,4.002D-01,3.947D-01,3.893D-01,
     &3.840D-01,3.789D-01,3.739D-01,3.690D-01,3.643D-01,
     &3.597D-01,3.551D-01,3.507D-01,3.464D-01,3.421D-01,
     &3.380D-01,3.340D-01,3.300D-01,3.261D-01,3.223D-01,
     &3.186D-01,3.150D-01,3.114D-01,3.079D-01,3.045D-01,
     &3.011D-01,2.978D-01,2.945D-01,2.914D-01,2.883D-01,
     &2.852D-01,2.822D-01,2.792D-01,2.763D-01,2.735D-01,
     &2.707D-01,2.679D-01,2.652D-01,2.625D-01,2.599D-01,
     &2.574D-01,2.548D-01,2.523D-01,2.499D-01,2.475D-01,
     &2.451D-01,2.428D-01,2.405D-01,2.382D-01,2.360D-01,
     &2.338D-01,2.316D-01,2.295D-01,2.274D-01,2.254D-01,
     &2.233D-01,2.213D-01,2.193D-01,2.174D-01,2.155D-01,
     &2.136D-01,2.117D-01,2.099D-01,2.081D-01,2.063D-01,
     &2.045D-01,2.028D-01,2.011D-01,1.994D-01,1.977D-01,
     &1.961D-01,1.944D-01,1.929D-01,1.913D-01,1.897D-01,
     &1.882D-01,1.867D-01,1.851D-01,1.837D-01,1.822D-01,
     &1.808D-01,1.793D-01,1.779D-01,1.765D-01,1.752D-01,
     &1.738D-01,1.725D-01,1.711D-01,1.698D-01,1.686D-01,
     &1.673D-01,1.660D-01,1.648D-01,1.635D-01,1.623D-01,
     &1.611D-01,1.599D-01,1.588D-01,1.576D-01,1.564D-01,
     &1.553D-01,1.542D-01,1.531D-01,1.520D-01,1.509D-01,
     &1.498D-01,1.488D-01,1.477D-01,1.467D-01,1.457D-01,
     &1.447D-01,1.437D-01,1.427D-01,1.417D-01,1.407D-01,
     &1.398D-01,1.388D-01,1.379D-01,1.369D-01,1.360D-01,
     &1.351D-01,1.342D-01,1.333D-01,1.324D-01,1.316D-01,
     &1.307D-01,1.299D-01,1.290D-01,1.282D-01,1.273D-01,
     &1.265D-01,1.257D-01,1.249D-01,1.241D-01,1.233D-01,
     &1.225D-01,1.218D-01,1.210D-01,1.203D-01,1.195D-01,
     &1.188D-01,1.180D-01,1.173D-01,1.166D-01,1.159D-01,
     &1.152D-01,1.144D-01,1.138D-01,1.131D-01,1.124D-01,
     &1.117D-01,1.110D-01,1.104D-01,1.097D-01,1.091D-01,
     &1.084D-01,1.078D-01,1.072D-01,1.065D-01,1.059D-01,
     &1.053D-01,1.047D-01,1.041D-01,1.035D-01,1.029D-01,
     &1.023D-01,1.017D-01,1.012D-01,1.006D-01,1.000D-01/
      
      DATA (ZZ(K,1),K=1,200 ) /
     &-2.303D+00,-3.084D+00,-3.649D+00,-4.098D+00,
     &-4.472D+00,-4.795D+00,-5.080D+00,-5.335D+00,
     &-5.568D+00,-5.781D+00,-5.978D+00,-6.161D+00,
     &-6.333D+00,-6.494D+00,-6.647D+00,-6.792D+00,
     &-6.929D+00,-7.060D+00,-7.186D+00,-7.306D+00,
     &-7.421D+00,-7.532D+00,-7.639D+00,-7.742D+00,
     &-7.842D+00,-7.938D+00,-8.031D+00,-8.122D+00,
     &-8.210D+00,-8.295D+00,-8.378D+00,-8.459D+00,
     &-8.538D+00,-8.614D+00,-8.689D+00,-8.762D+00,
     &-8.834D+00,-8.904D+00,-8.972D+00,-9.039D+00,
     &-9.104D+00,-9.168D+00,-9.231D+00,-9.293D+00,
     &-9.353D+00,-9.412D+00,-9.470D+00,-9.528D+00,
     &-9.584D+00,-9.639D+00,-9.693D+00,-9.746D+00,
     &-9.799D+00,-9.851D+00,-9.901D+00,-9.951D+00,
     &-1.000D+01,-1.005D+01,-1.010D+01,-1.014D+01,
     &-1.019D+01,-1.024D+01,-1.028D+01,-1.033D+01,
     &-1.037D+01,-1.041D+01,-1.046D+01,-1.050D+01,
     &-1.054D+01,-1.058D+01,-1.062D+01,-1.066D+01,
     &-1.070D+01,-1.074D+01,-1.078D+01,-1.082D+01,
     &-1.086D+01,-1.089D+01,-1.093D+01,-1.097D+01,
     &-1.101D+01,-1.104D+01,-1.108D+01,-1.111D+01,
     &-1.115D+01,-1.118D+01,-1.122D+01,-1.125D+01,
     &-1.128D+01,-1.132D+01,-1.135D+01,-1.138D+01,
     &-1.141D+01,-1.145D+01,-1.148D+01,-1.151D+01,
     &-1.154D+01,-1.157D+01,-1.160D+01,-1.163D+01,
     &-1.166D+01,-1.169D+01,-1.172D+01,-1.175D+01,
     &-1.178D+01,-1.181D+01,-1.184D+01,-1.186D+01,
     &-1.189D+01,-1.192D+01,-1.195D+01,-1.198D+01,
     &-1.200D+01,-1.203D+01,-1.206D+01,-1.208D+01,
     &-1.211D+01,-1.214D+01,-1.216D+01,-1.219D+01,
     &-1.221D+01,-1.224D+01,-1.226D+01,-1.229D+01,
     &-1.231D+01,-1.234D+01,-1.236D+01,-1.239D+01,
     &-1.241D+01,-1.244D+01,-1.246D+01,-1.248D+01,
     &-1.251D+01,-1.253D+01,-1.255D+01,-1.258D+01,
     &-1.260D+01,-1.262D+01,-1.264D+01,-1.267D+01,
     &-1.269D+01,-1.271D+01,-1.273D+01,-1.276D+01,
     &-1.278D+01,-1.280D+01,-1.282D+01,-1.284D+01,
     &-1.286D+01,-1.289D+01,-1.291D+01,-1.293D+01,
     &-1.295D+01,-1.297D+01,-1.299D+01,-1.301D+01,
     &-1.303D+01,-1.305D+01,-1.307D+01,-1.309D+01,
     &-1.311D+01,-1.313D+01,-1.315D+01,-1.317D+01,
     &-1.319D+01,-1.321D+01,-1.323D+01,-1.325D+01,
     &-1.327D+01,-1.329D+01,-1.330D+01,-1.332D+01,
     &-1.334D+01,-1.336D+01,-1.338D+01,-1.340D+01,
     &-1.342D+01,-1.343D+01,-1.345D+01,-1.347D+01,
     &-1.349D+01,-1.351D+01,-1.352D+01,-1.354D+01,
     &-1.356D+01,-1.358D+01,-1.360D+01,-1.361D+01,
     &-1.363D+01,-1.365D+01,-1.366D+01,-1.368D+01,
     &-1.370D+01,-1.372D+01,-1.373D+01,-1.375D+01,
     &-1.377D+01,-1.378D+01,-1.380D+01,-1.382D+01/
      
      DATA (ZZ(K,2),K=1,200 ) /
     &-2.303D+00,-2.512D+00,-2.700D+00,-2.871D+00,
     &-3.029D+00,-3.175D+00,-3.310D+00,-3.438D+00,
     &-3.557D+00,-3.670D+00,-3.778D+00,-3.880D+00,
     &-3.977D+00,-4.070D+00,-4.159D+00,-4.245D+00,
     &-4.328D+00,-4.407D+00,-4.484D+00,-4.558D+00,
     &-4.630D+00,-4.699D+00,-4.767D+00,-4.832D+00,
     &-4.896D+00,-4.958D+00,-5.019D+00,-5.078D+00,
     &-5.135D+00,-5.191D+00,-5.246D+00,-5.300D+00,
     &-5.352D+00,-5.403D+00,-5.453D+00,-5.503D+00,
     &-5.551D+00,-5.598D+00,-5.645D+00,-5.690D+00,
     &-5.735D+00,-5.779D+00,-5.822D+00,-5.864D+00,
     &-5.906D+00,-5.947D+00,-5.988D+00,-6.027D+00,
     &-6.067D+00,-6.105D+00,-6.143D+00,-6.181D+00,
     &-6.217D+00,-6.254D+00,-6.290D+00,-6.325D+00,
     &-6.360D+00,-6.394D+00,-6.428D+00,-6.462D+00,
     &-6.495D+00,-6.528D+00,-6.560D+00,-6.592D+00,
     &-6.624D+00,-6.655D+00,-6.686D+00,-6.716D+00,
     &-6.746D+00,-6.776D+00,-6.805D+00,-6.835D+00,
     &-6.863D+00,-6.892D+00,-6.920D+00,-6.948D+00,
     &-6.976D+00,-7.003D+00,-7.030D+00,-7.057D+00,
     &-7.084D+00,-7.110D+00,-7.136D+00,-7.162D+00,
     &-7.188D+00,-7.213D+00,-7.238D+00,-7.263D+00,
     &-7.288D+00,-7.312D+00,-7.336D+00,-7.360D+00,
     &-7.384D+00,-7.408D+00,-7.431D+00,-7.455D+00,
     &-7.478D+00,-7.501D+00,-7.523D+00,-7.546D+00,
     &-7.568D+00,-7.590D+00,-7.612D+00,-7.634D+00,
     &-7.656D+00,-7.677D+00,-7.698D+00,-7.720D+00,
     &-7.741D+00,-7.761D+00,-7.782D+00,-7.803D+00,
     &-7.823D+00,-7.843D+00,-7.863D+00,-7.883D+00,
     &-7.903D+00,-7.923D+00,-7.943D+00,-7.962D+00,
     &-7.981D+00,-8.001D+00,-8.020D+00,-8.039D+00,
     &-8.057D+00,-8.076D+00,-8.095D+00,-8.113D+00,
     &-8.132D+00,-8.150D+00,-8.168D+00,-8.186D+00,
     &-8.204D+00,-8.222D+00,-8.239D+00,-8.257D+00,
     &-8.274D+00,-8.292D+00,-8.309D+00,-8.326D+00,
     &-8.343D+00,-8.360D+00,-8.377D+00,-8.394D+00,
     &-8.411D+00,-8.427D+00,-8.444D+00,-8.460D+00,
     &-8.476D+00,-8.493D+00,-8.509D+00,-8.525D+00,
     &-8.541D+00,-8.557D+00,-8.572D+00,-8.588D+00,
     &-8.604D+00,-8.619D+00,-8.635D+00,-8.650D+00,
     &-8.666D+00,-8.681D+00,-8.696D+00,-8.711D+00,
     &-8.726D+00,-8.741D+00,-8.756D+00,-8.771D+00,
     &-8.786D+00,-8.800D+00,-8.815D+00,-8.829D+00,
     &-8.844D+00,-8.858D+00,-8.872D+00,-8.887D+00,
     &-8.901D+00,-8.915D+00,-8.929D+00,-8.943D+00,
     &-8.957D+00,-8.971D+00,-8.985D+00,-8.998D+00,
     &-9.012D+00,-9.026D+00,-9.039D+00,-9.053D+00,
     &-9.066D+00,-9.080D+00,-9.093D+00,-9.106D+00,
     &-9.119D+00,-9.133D+00,-9.146D+00,-9.159D+00,
     &-9.172D+00,-9.185D+00,-9.197D+00,-9.210D+00/
      END
C=======================================================================

      DOUBLE PRECISION FUNCTION CHIDIS (KPARTin, IFL1, IFL2)

C-----------------------------------------------------------------------
C...Generate CHI (fraction of energy of a hadron carried by 
C.                the valence quark, or diquark, as specified by IFL1)
C.  INPUT KPART = code of particle
C.        IFL1, IFL2 = codes of partons (3, 3bar of color)
C.........................................................
      IMPLICIT NONE
c     external types
      INTEGER KPARTIN,IFL1,IFL2
c     COMMONs
      INTEGER NW_max
      PARAMETER (NW_max = 20)
C--------------------------------------------------------------------
C     SIBYLL common blocks containing event information       \FR'14
C--------------------------------------------------------------------

C     EVENT INFO COMMON
C     contains overall interaction properties, like
C     SQS : center-of-mass energy
C     S   :         "       "     squared
C     PTmin : low pt cut of QCD cross section, 
C             i.e. minimal pt of hard minijets
C     Xmin : low-x bound for PDFs, 
C            i.e. minimal momentum fraction of hard partons
C     Zmin : logarithm of that
C     KB : PID of beam hadron
C     KT() : PID of target
C     IAT : mass number of target
      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
      INTEGER KB,IAT,KT
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      DOUBLE PRECISION CCHIK
      COMMON /S_CPSPL/ CCHIK(4,99)
      
      DOUBLE PRECISION STR_mass_val, STR_mass_val_hyp, STR_mass_sea
      COMMON /S_CUTOFF/ STR_mass_val, STR_mass_val_hyp, STR_mass_sea
c     internal types
      DOUBLE PRECISION CUT,S_RNDM
      INTEGER KPART,IFQ
      SAVE

      kpart=IABS(kpartin)
      IFQ=IABS(IFL1)
      IF (IFQ.GT.10) IFQ=IABS(IFL2)
      CUT=2.D0*STR_mass_val/SQS
c     hyperon beam cut
      IF(kpart.gt.14) CUT=2.D0*STR_mass_val_hyp/SQS
100   CHIDIS=S_RNDM(0)**2
      if (chidis.lt.cut) goto 100
      if (chidis.gt.(1.D0-cut)) goto 100
      IF((CHIDIS**2/(CHIDIS**2+CUT**2))**0.5D0
     +   *(1.D0-CHIDIS)**CCHIK(IFQ,KPART).LT.S_RNDM(1)) GOTO 100
      CHIDIS = MAX(0.5D0*CUT,CHIDIS)
      CHIDIS = MIN(1.D0-CUT,CHIDIS)
c     diquarks or charm quarks      
      IF (IABS(IFL1).GT.3)  CHIDIS=1.D0-CHIDIS
      RETURN
      END
C=======================================================================

      FUNCTION QMASS(IFL)

C-----------------------------------------------------------------------
C...Return quark or diquark constituent masses
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      DIMENSION QMAS(4)
      SAVE
      DATA QMAS /0.325D0,0.325D0,0.5D0,1.5D0/

      IFLA = IABS(IFL)
      IFLA = MOD(IFLA,100)
      IF (IFLA .LE. 4)       THEN
         QMASS = QMAS(IFLA)
      ELSE
         QMA = QMAS(IFLA/10)
         QMB = QMAS(MOD(IFLA,10))
         QMASS = QMA+QMB
      ENDIF
      RETURN
      END
C=======================================================================

      FUNCTION XM2DIS(XM2MIN,XM2MAX,ALPHA)

C-----------------------------------------------------------------------
C     function that samples mass**2 from (1/M**2)**alpha
C     with alpha <= 1                                             
C     INPUT: Mmin**2 : minimal mass
C            Mmax**2 : maximal mass
C            alpha   : slope                                      \FR'14
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      SAVE
      
c     reduced alpha
      ALPHArdc = 2.d0*(ALPHA-1.d0)
      AMIN = LOG(XM2MIN)
      AMAX = LOG(XM2MAX)
      ADLT = AMAX-AMIN
      IF(ABS(ALPHArdc).LT.1.d-3)THEN
c     alpha = 1
         XRNDM = MAX(S_RNDM(0),1.D-10)
         AX = AMIN+ADLT*XRNDM
         XM2DIS = EXP(AX)
      ELSEIF(ALPHArdc.LT.0.D0.and.ALPHA.gt.0.D0)THEN
c     0 < alpha < 1
         XRNDM = MAX(S_RNDM(0),1.D-10)
c     AX = AMAX-LOG(XRNDM)*ALPHArdc
         DX = XM2MAX**(1.D0-ALPHA)*XRNDM +
     +        XM2MIN**(1.D0-ALPHA)*(1.D0-XRNDM)
         AX = LOG(DX)/(1.D0-ALPHA)
         XM2DIS = EXP(AX)
      ELSEIF(ALPHArdc.GE.1.D0)THEN
c     alpha >= 1
         ALPHAr = 1.D0-ALPHA
         XMINA = XM2MIN**ALPHAr
         XMAXA = XM2MAX**ALPHAr
         XDLT = XMAXA-XMINA
         XRNDM = MAX(S_RNDM(0),1.D-10)
         Z = LOG(XMINA+XDLT*XRNDM)/ALPHAR
         XM2DIS = EXP(Z)
      ELSE
         WRITE(6,*) 'M2DIS: undefined exponent in mass distribution!',
     &        ALPHA
         XM2DIS = 0.D0
         CALL SIB_REJECT('M2DIS           ')
      ENDIF
      END
C=======================================================================

      SUBROUTINE EXCTDEC( IDX, LBAD)

C-----------------------------------------------------------------------
C     routine to fragment an excited system with known flavor via
C     resonance decay
C-----------------------------------------------------------------------
      IMPLICIT NONE
c     external variables
      INTEGER IDX,LBAD

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
C     The final particle output is contained in COMMON /S_PLIST/    
C     NP           : number of final particles
C     P(1:NP, 1:5) : 4-momenta + masses of the final particles 
C     LLIST (1:NP) : codes of final particles
      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP

      INTEGER NFORIG,NPORIG,NIORIG,IPFLAG,IIFLAG,KINT
      COMMON /S_PARTO/ NFORIG(NP_max),NPORIG(NP_max),NIORIG(NP_max),
     &IPFLAG,IIFLAG,KINT

      INTEGER LRNK
      COMMON /SIB_RNK/ LRNK(8000)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)

      INTEGER ICHP,ISTR,IBAR
      COMMON /S_CHP/ ICHP(99), ISTR(99), IBAR(99)

      INTEGER IISO,ISPN
      COMMON /S_SPN/ IISO(99), ISPN(99)

      INTEGER ICHM
      COMMON /S_CHM/ ICHM(99)

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN

c     local variables
      DOUBLE PRECISION P0,BE,PR1,PR2,PRH,GABE,P2,
     &     PAR2_def,PAR8_def,PAR24_def,DELTAE,PCXG,
     &     EMIN1,EMIN2,EMIN3,EMIN4,S_RNDM,GA,PTR,PTOT,P1TOT,PX,PY,
     &     COD,SID,COF,SIF,ANORF,BEP
      DIMENSION P0(5),BE(3),PR1(5),PR2(5),PRH(5),GABE(4),
     &     P2(5)
      INTEGER IPID,IR1DX,IFLR1,IR2DX,IFLR2,IRH,IRHPID,IR,
     &     KK,KD,IFAIL,N1,IFBAD,J,K,I
      SAVE
      
c      LBAD = 1

c     initial parameters
      PAR2_def = PAR(2)         ! ud/s rate
      PAR8_def = PAR(8)         ! popcorn rate
      PAR24_def = PAR(24)       ! c/s rate
      if(ndebug.gt.1)
     &     WRITE(LUN,*) ' EXCTDEC: IDX,IREJ',IDX,LBAD
      
c     read remnant 4momentum from stack
      CALL RD_PRTN_4VEC(IDX,P0,IPID,IR1DX)
      CALL RD_PRTN_4VEC(IR1DX,PR1,IFLR1,IR2DX)
      CALL RD_PRTN_4VEC(IR2DX,PR2,IFLR2,IRH)
      CALL RD_PRTN_4VEC(IRH,PRH,IRHPID,IR)
      IPFLAG = IPID
      IF(IDX.ne.IR)then
         write(lun,*) ' EXCTDEC: reference loop broken!',IDX,IR
         CALL SIB_REJECT('EXCTDEC         ')
      endif
      IF(NDEBUG.GT.2)THEN
         WRITE(LUN,*) ' EXCTDEC: P0:' , (P0(kk),kk=1,5)
         WRITE(LUN,*) ' EXCTDEC: PR1:' , (PR1(kk),kk=1,5)
         WRITE(LUN,*) ' EXCTDEC: PR2:' , (PR2(kk),kk=1,5)
         WRITE(LUN,*) ' EXCTDEC: PH:' , (PRH(kk),kk=1,5)
      ENDIF
      
C     identity of remnant
c     form hadron from flavors in remnant
c     (not preserving spin or isospin!)
c      CALL SIB_I4FLAV(iflr1,iflr2,Idm, KD )
      KD = IRHPID

c     available kinetic energy
      DELTAE = P0(5)-AM(ABS(KD))
c     fallback region: 0 < DELTAE < EMIN1
      EMIN1 = PAR(76)
c     resonance region: EMIN1 < DELTAE < EMIN2
      EMIN2 = PAR(77)
c     phasespace decay region: EMIN2 < DELTAE < EMIN3
      EMIN3 = PAR(78)
c     string decay region: EMIN3 < DELTAE < EMIN4
      EMIN4 = PAR(79)

      IF(NDEBUG.gt.2)THEN
         WRITE(LUN,*) 
     &        ' EXCTDEC: MASS,IFL1,IFL2,PID',P0(5),IFLR1,IFLR2,KD
         WRITE(LUN,*) ' EXCTDEC: DELTAE,EMIN1,EMIN2,EMIN3',
     &        DELTAE,EMIN1,EMIN2,EMIN3
      ENDIF
      
c     strange quark rate
      IF(IPAR(48).eq.1)THEN
         PAR(2) = PAR(89)
      ENDIF

c     charm quark rate
      IF(IPAR(62).eq.1)THEN
         PAR(24) = PAR(107)
      ENDIF
     
c     popcorn rate in remnant
      IF(IPAR(56).eq.1)THEN
         PAR(8) = PAR(102)
      ENDIF     

      IF(DELTAE.lt.EMIN2.and.abs(P0(5)-am(abs(kd))).lt.EPS5)THEN
c     beam or resonance region
         IF(NDEBUG.gt.1) then 
            if(DELTAE.lt.EMIN1)then
               WRITE(LUN,*)' EXCTDEC: fallback to beam..'
            else
               WRITE(LUN,*)' EXCTDEC: forming resonance..'
            endif
         endif
         NP = NP + 1
         LLIST(NP) = KD
         NPORIG(NP) = IPFLAG
         LRNK(NP) = 0
         niorig(NP) = iiflag
         DO kk=1,5
            P(NP,KK) = P0(KK)
         ENDDO
         LBAD = 0
         PAR(2) = PAR2_def
         PAR(8) = PAR8_def
         PAR(24) = PAR24_def
         RETURN         

      ELSEIF(DELTAE.lt.EMIN3)THEN
c     phasespace decay region
         IF(NDEBUG.gt.1) WRITE(LUN,*)' EXCTDEC: phasespace decay ..'
         IPFLAG = IPID/iabs(IPID) + ISIGN(1000,IPID)
c     set charge exchange probability, 
c     i.e. prob for p* -> n + pip
         PCXG = PAR(99)
         CALL FIREBALL_4FLV(KD,P0,PCXG,IFAIL)
         PAR(2) = PAR2_def
         PAR(8) = PAR8_def
         PAR(24) = PAR24_def
         IF(IFAIL.eq.1) THEN
            IF(ndebug.gt.0)
     &           WRITE(LUN,*) ' EXCTDEC: remnant frag. rejection!'
            LBAD = 1
            RETURN
         ENDIF
         LBAD = 0
         RETURN

c      ELSEIF(DELTAE.lt.EMIN4)THEN
      ELSE
C     string fragmentation region
         IF(NDEBUG.gt.1) WRITE(LUN,*)' EXCTDEC: string decay ..'
         N1 = NP+1
         IPFLAG = IPFLAG + ISIGN(3000,IPID)
c     for meson remnant quark and anti-quark should be treated equally
c     therefor switch randomly
         IF(IBAR(ABS(KD)).eq.0.and.S_RNDM(KD).lt.0.5D0)
     &        CALL ISWTCH_LMNTS(IFLR1,IFLR2)

c     turn remnant string around
         IF(IPAR(23).eq.1)THEN
            IF(S_RNDM(0).gt.PAR(39))
     &           CALL ISWTCH_LMNTS(IFLR1,IFLR2)
         ENDIF

         CALL STRING_FRAG_4FLV 
     +        (P0(5), IFLR2, IFLR1, 0.D0,0.D0,0.D0,0.D0,IFBAD,1)
         IF (IFBAD .EQ. 1)THEN
            IF(ndebug.gt.0)
     &           WRITE(LUN,*) ' EXCTDEC: remnant frag. rejection!'
            LBAD = 1
            PAR(2) = PAR2_def
            PAR(8) = PAR8_def
            PAR(24) = PAR24_def
            RETURN
         ENDIF
         DO J=1,3
            BE(J)=P0(J)/P0(4)
            GABE(J)=P0(J)/P0(5)
         ENDDO
         GA=P0(4)/P0(5)
         GABE(4)=P0(4)/P0(5)
C...  rotate and boost string
         IF(IPAR(38).eq.1.or.IPAR(38).eq.3)THEN
c     sample additional soft pt for remnant partons
            CALL PTDIS_4FLV(0,PX,PY)
            PTR = SQRT(PX**2+PY**2)
            PTOT = SQRT(4.D0*PTR**2+P0(5)**2)*0.5D0
c     rotation factors
            COD = 0.5D0*P0(5)/PTOT
            SID = PTR/PTOT
c            COD= 1.D0/SQRT(1.D0+4.D0*PTR**2/P0(5))
c            SID= 2.D0*PTR/P0(5)*COD
            COF=1.D0
            SIF=0.D0
            IF(PTOT*SID.GT.EPS5) THEN
               COF=PX/(SID*PTOT)
               SIF=PY/(SID*PTOT)
               ANORF=DSQRT(COF*COF+SIF*SIF)
               COF=COF/ANORF
               SIF=SIF/ANORF
            ENDIF
            IF(ndebug.gt.3)THEN
            write(lun,*)' EXCTDEC: rotation factors (cod,sid,cof,sif):',
     &           cod,sid,cof,sif
            write(lun,*)' EXCTDEC: rotation angles (theta,phi):',
     &           ACOS(cod),ACOS(cof)
            ENDIF
c     rotate string final state
            DO K=N1,NP
               CALL SIB_TRANI(P(K,1),P(k,2),P(k,3),cod,sid,cof,sif
     &              ,P2(1),P2(2),P2(3))
               do j=1,3
                  P(K,j)=P2(j)
               enddo
            ENDDO
c     boost to hadron-hadron center-of-mass
            IF(ndebug.gt.3)
     &        write(lun,*) ' EXCTDEC: boost to had-had (gabe,gam):',
     &           (gabe(j),j=1,4)
            DO K=N1,NP
               NPORIG(K) = IPFLAG
               niorig(K) = iiflag
               CALL SIB_ALTRA(gabe(4),gabe(1),gabe(2),
     &              gabe(3),P(k,1),p(k,2),p(k,3),p(k,4),
     &              P1TOT,p2(1),p2(2),p2(3),p2(4))
               do j=1,4
                  P(K,j)=P2(j)
               enddo
            ENDDO
         ELSEIF(IPAR(38).eq.2.or.IPAR(38).eq.0)THEN            
C...  boost string
            DO I=N1,NP
               NPORIG(I) = IPFLAG
               niorig(I) = iiflag
               BEP=BE(1)*P(I,1)+BE(2)*P(I,2)+BE(3)*P(I,3)
               DO J=1,3
                  P(I,J)=P(I,J)+GA*(GA*BEP/(1.D0+GA)+P(I,4))*BE(J)
               ENDDO
               P(I,4)=GA*(P(I,4)+BEP)
            ENDDO
         ENDIF
      ENDIF
      LBAD = 0
      PAR(2) = PAR2_def
      PAR(8) = PAR8_def
      PAR(24) = PAR24_def
      return
      END
C=======================================================================

      SUBROUTINE PTDIS_4FLV (IFL,PX,PY)

C-----------------------------------------------------------------------
C...Generate pT
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      DOUBLE PRECISION PPT02
      COMMON /S_CQDIS2/ PPT02(44)
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      SAVE

      IF(IFL.eq.0)THEN
c     quark confinement pt
         PPTT = PAR(110)
         XM = 0.325D0
         XM2 = XM**2
         RNDM = MAX(EPS10,S_RNDM(IFL))
         XMT = PPTT * LOG(RNDM) - XM
         XMT2 = XMT**2
         PT = SQRT(XMT2-XM2)         
      ELSE
         IFLA = IABS(IFL)
         IFLA = MOD(IFLA,100)
         PPTT = PPT02(IFLA)
c     Gaussian distribution
         PT = PPTT*SQRT(-LOG(MAX(EPS10,S_RNDM(IFL))))
         IF (IPAR(3).GE.1) THEN
            IF(MOD(IFLA,10).NE.0) THEN
               XM = QMASS(IFL)
            ELSE
               XM = 0.5D0        ! pomeron mass
               IF(IPAR(3).ge.6) XM = 0.D0
            ENDIF
c     exponential transverse mass
            XM2 = XM**2
            RNDM = MAX(EPS10,S_RNDM(IFL))
            XMT = PPTT * LOG(RNDM) - XM
            XMT2 = XMT**2
            PT = SQRT(XMT2-XM2)
         ENDIF      
      ENDIF
      PHI= TWOPI*S_RNDM(IFL)
      PX=PT*COS(PHI)
      PY=PT*SIN(PHI)
      RETURN
      END

C=======================================================================

      SUBROUTINE PTSETUP_4FLV(ECM)

C-----------------------------------------------------------------------
C     moved from sib_ndiff to seperate subroutine 
c     so that changes will affect diff. /FR'13
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)

      DOUBLE PRECISION PPT02
      COMMON /S_CQDIS2/ PPT02(44)
      SAVE

      SQS = ECM

c     NA22 piC retune
      PTU=.3D0+.08D0*dlog10(sqs/30.D0)
      PTS=.45D0+.08D0*dlog10(sqs/30.D0)
      PTQQ=.6D0+.08D0*dlog10(sqs/30.D0)
      PTPOM= .6D0+.08D0*dlog10(sqs/30.D0)
      if ( IPAR(3).eq.1 ) then 
c     pt0
         ptu=.15D0+.007D0*dlog10(sqs/20.D0)**2
         pts=.3D0+.007D0*dlog10(sqs/20.D0)**2
         ptqq=.3D0+.03D0*dlog10(sqs/20.D0)**2
         ptpom= .6D0+.08D0*dlog10(sqs/30.D0)
      elseif ( IPAR(3).eq.2 ) then
C     pt1
         ptu=.15D0+.007D0*dlog10(sqs/20.D0)**2
         pts=.32D0+.007D0*dlog10(sqs/20.D0)**2
         ptqq=.4D0+.007D0*dlog10(sqs/20.D0)**2
         ptpom= .6D0+.08D0*dlog10(sqs/30.D0)
c     pt2
      elseif ( IPAR(3).eq.3 ) then
         ptu=.17D0+.007D0*dlog10(sqs/20.D0)**2
         pts=.3D0+.007D0*dlog10(sqs/20.D0)**2
         ptqq=.3D0+.03D0*dlog10(sqs/20.D0)**2
         ptpom = .6D0+.08D0*dlog10(sqs/30.D0)
      elseif ( IPAR(3).eq.5 ) then
         PTU=.16D0+.007D0*dlog10(sqs/20.D0)**2
         PTS=.28D0+.007D0*dlog10(sqs/20.D0)**2
         PTQQ= .3D0+.03D0*dlog10(sqs/20.D0)**2
         PTPOM = .23D0+.03D0*dlog10(sqs/20.D0)**2
      elseif ( IPAR(3).eq.6 ) then
         PTU=.16D0+.007D0*dlog10(sqs/20.D0)**2
         PTS=.28D0+.007D0*dlog10(sqs/20.D0)**2
         PTQQ= .3D0+.03D0*dlog10(sqs/20.D0)**2
         PTPOM = .23D0+.03D0*dlog10(sqs/20.D0)**2
      elseif ( IPAR(3).eq.7 ) then
         PTU= PAR(46) + .007D0*dlog10(sqs/20.D0)**2
         PTS= PAR(47) + .007D0*dlog10(sqs/20.D0)**2
         PTQQ= PAR(48) + .03D0*dlog10(sqs/20.D0)**2
         PTPOM = PAR(49) + .03D0*dlog10(sqs/20.D0)**2
      elseif ( IPAR(3).eq.8 ) then
         ASQS = MAX(log10(SQS/PAR(109)),0.D0)
         PTU= PAR(46) + PAR(68)*ASQS**2
         PTS= PAR(47) + PAR(70)*ASQS**2
         PTQQ= PAR(48) + PAR(69)*ASQS**2
         PTPOM = PAR(49) + PAR(51)*ASQS**2
         PTSEA = PAR(67) + PAR(52)*ASQS**2
      endif
      PPT02 (1) = PTU
      PPT02 (2) = PTU
      PPT02 (3) = PTS
c     valence pt
      PPT02 (10) = PTPOM
      DO J=11,33
         PPT02(J) = PTQQ
      ENDDO
c     soft minijet pt
      PPT02 (20) = PTSEA
c     sea quark pt
      PPT02 (30) = PAR(132)
c     charm pt
      ASQS = MAX(log10(SQS/30.D0),0.D0)
      IF(IPAR(16).eq.8)THEN
         PTCHM= PAR(147) + PAR(149)*ASQS
         PTCHB= PAR(148) + PAR(149)*ASQS         
      ELSE
c     rc4a charm pt
         PTCHM=0.308D0 + .165D0*ASQS
         PTCHB=0.5D0 + .165D0*ASQS         
      ENDIF
      PPT02(4) = PTCHM
      PPT02(14) = PTCHB
      PPT02(24) = PTCHB
      DO J=34,44
         PPT02(J) = PTCHB
      ENDDO
     
      IF(ndebug.gt.2)THEN
         WRITE(LUN,*)' PTSETUP_4FLV: (sqs,(u,d),s,diq,pom,cm,cb)',sqs
     +      ,ppt02(1),ppt02(3),ppt02(11), ppt02(10),ppt02(4),ppt02(34)
      ENDIF

      RETURN
      END
C=======================================================================

      INTEGER FUNCTION IMRG2HAD(IFLB1,IFLB2)

C-----------------------------------------------------------------------
C     -----------------------------------------------------
C     function that merges two flavors into lightest hadron
C     -----------------------------------------------------
      IMPLICIT NONE
c     flavor merging array
      INTEGER KFLV
      COMMON /S_KFLV/ KFLV(4,43)
      INTEGER IFLB1,IFLB2,IFLA,IFLB,IFL1,IFL2
      SAVE

      IFLA = IFLB1
      IFLB = IFLB2
c     order by flavor, meson: antiquark-quark, baryon: quark-diquark
      IF(IFLB.lt.IFLA) CALL ISWTCH_LMNTS(ifla,iflb)
c     if antibaryon switch again..
      IF(IFLB.lt.0) CALL ISWTCH_LMNTS(ifla,iflb)
      IFL1 = IABS(IFLA)
      IFL2 = IABS(IFLB)
      IMRG2HAD = ISIGN(KFLV(IFL1,IFL2),IFLB)
      END

C=======================================================================

      SUBROUTINE SAMPLE_SEA_TOT
     &     (KRMNT,KINT,NSEA,XGAM,XJET,STR_MASS,XSJ,XX)

C-----------------------------------------------------------------------
C   input parameter: xgam,xjet,str_mass,  Nsea,KINT,krmnt
c   outpt parameter: xsj,xx
C-----------------------------------------------------------------------
      IMPLICIT NONE

c     include COMMON blocks
      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      INTEGER NW_max
      PARAMETER (NW_max = 20)
C--------------------------------------------------------------------
C     SIBYLL common blocks containing event information       \FR'14
C--------------------------------------------------------------------
C     EVENT INFO COMMON
C     contains overall interaction properties, like
C     SQS : center-of-mass energy
C     S   :         "       "     squared
C     PTmin : low pt cut of QCD cross section, 
C             i.e. minimal pt of hard minijets
C     Xmin : low-x bound for PDFs, 
C            i.e. minimal momentum fraction of hard partons
C     Zmin : logarithm of that
C     KB : PID of beam hadron
C     KT() : PID of target
C     IAT : mass number of target
      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
      INTEGER KB,IAT,KT
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      
c     input/output type definitions
      DOUBLE PRECISION XGAM,XJET,STR_MASS,XSEA,XX,XSJ
      DIMENSION XX(2*NW_max+2)

      INTEGER NSEA,KINT,KRMNT

c     local type definitions
      DOUBLE PRECISION AC,GAMMAX,S_RNDM,XA,XREM,R,Z,Z1,Z2,XMINA
      INTEGER j,jj,ilast
      SAVE
      DATA AC /-0.2761856692D0/ ! log(2) - gamma(Eulero)

      GAMMAX = xgam
      XMINA = 2.D0*STR_mass/SQS
      IF(IPAR(73).eq.1.and.KINT.gt.1) GAMMAX = PAR(119)
      IF(ndebug.gt.3) THEN
         WRITE(LUN,*)' IMRG2HAD: called with ',
     &        '(KRMNT,KINT,NSEA,XGAM,XJET,STR_MASS):', 
     &        KRMNT,KINT,NSEA,XGAM,XJET,STR_MASS
         
         WRITE(LUN,*)' IMRG2HAD: XMIN,XMIN*N,XREM:',
     &        XMINA,NSEA*XMINA,1.D0-XJET
      ENDIF
c     sample total fraction for sea partons..
      Z1 = LOG(DBLE(NSEA))
 50   Z2 = LOG(0.5D0*SQS*(1.D0-XJET)/STR_MASS-2.D0)
      R = S_RNDM(0)
      Z=(Z1+AC)*(1.D0+R*(((Z2+AC)/(Z1+AC))**NSEA-1.D0))
     &     **(1.D0/DBLE(NSEA))-AC
      XSEA = XMINA*EXP(Z)
      IF(ndebug.gt.3) WRITE(LUN,*) '  total SEA fraction:' , xsea
      IF ( (1.D0-XSEA)**GAMMAX .LT. S_RNDM(1)) GOTO 50
c     maximal fraction remaining for valence..
 60   XREM = XSEA - DBLE(Nsea)*XMINA
      IF(ndebug.gt.3) 
     &     WRITE(LUN,*) '  Xsea,xval,xjet:',
     &     xsea,1.D0-XSEA-XJET,xjet
      
C...  Split the energy of sea partons among the different partons
      DO j=1,Nsea-1
         jj = 2+j
         IF(KRMNT.eq.0) jj = 4+j
c     fraction for first parton
         XA = XREM*S_RNDM(J)
c     for interactions other than first decrease energy fraction
c     (beam side hadron can participate in multiple binary collisions)
c     IF(KINT.gt.1.and.j.gt.2*KRMNT) XA=SIGN(ABS(XA)**PAR(116),XA)
         XX(jj) = XMINA + XA
c     new remainder
         XREM = XREM - XA
         IF(ndebug.gt.3) write(lun,*)'  x1,j,rem,xa',xX(jj),jj,xrem,xa
      ENDDO
c     last parton..
      ilast = 2+Nsea
      IF(KRMNT.eq.0) ilast = 4+Nsea
      XX(ILAST) = XMINA + XREM

c     break symmetry between nucleon interactions
c     first interaction takes most energy
      IF(KINT.gt.1.and.IPAR(71).eq.1)THEN
         JJ = 3
         IF(KRMNT.eq.0) JJ = 5
         if(ndebug.gt.4) write(lun,*) '  x1+x2,p*xeq:',
     &        XX(JJ)+XX(JJ+1),PAR(117)*XSEA/KINT
         IF(XX(JJ)+XX(JJ+1).lt.PAR(117)*XSEA/KINT) GOTO 60
      ENDIF

      XSJ = XSJ + XSEA
      IF(ndebug.gt.3)THEN  
         write(lun,*)'  x1,N,rem',xx(ilast),ilast,xrem
         write(lun,*) '  xseajet',xsj
      endif

      END
C-----------------------------------------------------------------------
C
C     dummy subroutines, remove to link PDFLIB
C
C=======================================================================
c
c     SUBROUTINE PDFSET(PARAM,VALUE)
c
c-----------------------------------------------------------------------
c     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     IMPLICIT INTEGER(I-N)
c     DIMENSION PARAM(20),VALUE(20)
c     CHARACTER*20 PARAM
c     END
c
c=======================================================================
c
c     SUBROUTINE STRUCTM(XI,SCALE2,UV,DV,US,DS,SS,CS,BS,TS,GL)
c
c-----------------------------------------------------------------------
c     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     IMPLICIT INTEGER(I-N)
c     END
c
c=======================================================================
c
c     SUBROUTINE STRUCTP(XI,SCALE2,P2,IP2,UV,DV,US,DS,SS,CS,BS,TS,GL)
c
c-----------------------------------------------------------------------
c     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     IMPLICIT INTEGER(I-N)
c     END
c
C-----------------------------------------------------------------------
C
C=======================================================================

      SUBROUTINE SIB_NDIFF(K_beam, NW, Ecm, Irec, IREJ)

C-----------------------------------------------------------------------
C     routine that samples and fragments a non-diffractive interaction
C
C     3 stages: 0: setup
C               1: sampling of event structure (number of parton interactions)
C                  (labeled as 2000)
C               2: sampling of kinematics
C                  (labeled as 3000)
C               3: fragmentation
C-----------------------------------------------------------------------
      IMPLICIT NONE

c     external types
      DOUBLE PRECISION ECM
      INTEGER K_beam, NW, Irec, IREJ

c     COMMONs
      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      INTEGER NW_max
      PARAMETER (NW_max = 20)
C--------------------------------------------------------------------
C     SIBYLL common blocks containing event information       \FR'14
C--------------------------------------------------------------------
C     EVENT INFO COMMON
C     contains overall interaction properties, like
C     SQS : center-of-mass energy
C     S   :         "       "     squared
C     PTmin : low pt cut of QCD cross section, 
C             i.e. minimal pt of hard minijets
C     Xmin : low-x bound for PDFs, 
C            i.e. minimal momentum fraction of hard partons
C     Zmin : logarithm of that
C     KB : PID of beam hadron
C     KT() : PID of target
C     IAT : mass number of target
      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
      INTEGER KB,IAT,KT
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT
C     The final particle output is contained in COMMON /S_PLIST/    
C     NP           : number of final particles
C     P(1:NP, 1:5) : 4-momenta + masses of the final particles 
C     LLIST (1:NP) : codes of final particles
      DOUBLE PRECISION P
      INTEGER NP,LLIST,NP_max
      PARAMETER (NP_max=8000)
      COMMON /S_PLIST/ P(NP_max,5), LLIST(NP_max), NP

      INTEGER NFORIG,NPORIG,NIORIG,IPFLAG,IIFLAG,KINT
      COMMON /S_PARTO/ NFORIG(NP_max),NPORIG(NP_max),NIORIG(NP_max),
     &IPFLAG,IIFLAG,KINT
C     parameters that represent: NW: max. number of wounded nucleons,
C     NS,NH: max. number of soft and hard interactions
c      PARAMETER (NW_max = 20)
C     The COMMON block /S_CHIST/ contains information about the
C     the structure of the  generated event:
C     NWD   = number of wounded nucleons
C     NJET = total number of hard interactions
C     NSOF = total number of soft interactions
C     NNSOF (1:NW) = number of soft pomeron cuts in each interaction
C     NNJET (1:NW) = number of minijets produced in each interaction 
C     JDIF(1:NW) = diffraction code 
C                  0 : non-diff,
C                  1 : beam-diff
C                  2 : target-diff
C                  3 : double-diff
      INTEGER NNSOF,NNJET,JDIF,NWD,NJET,NSOF
      COMMON /S_CHIST/ NNSOF(NW_max),NNJET(NW_max),
     &     JDIF(NW_max),NWD,NJET,NSOF
      INTEGER NS_max, NH_max
      PARAMETER (NS_max = 20, NH_max = 80)

      INTEGER IBMRDX,ITGRDX,IHMJDX,ISMJDX,ICSTDX,IINTDX
      COMMON /S_INDX/ IBMRDX(3),ITGRDX(NW_max,3),
     &     IHMJDX(NW_max*NH_max),IINTDX(NW_max),
     &     ISMJDX(NW_max*NS_max),ICSTDX(2*NW_max,3)
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)
      INTEGER ITRY, NREJ
      COMMON /S_CNT/ ITRY(20), NREJ(20)

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN

c     internal type declarations
      DOUBLE PRECISION X2JET,SQS_0,PZ,E2,PAWT,xnsof,xnjet,xjdif,x1jet,
     &     Esum,PXsum,PYsum,PZsum
      DIMENSION X2JET(NW_max)
      INTEGER LL,LXBAD,NP_0,NPP_0,NPP0_0,J,JJ,I,KBA,L,NPP_1,NPP0_1,
     &     IREFout,IREF,nj,ns,nv,II,Idm,LPID,NF,NPP,NPP0
      DIMENSION LL(99)      
      SAVE
      DATA LL /5*0,7*2,2*1,12*0,2,6*0,6*1,19*0,2,2,10*0,
     &     2,2,0,2,2,11*0,1,1,1,9*0,1/


C..   setup stage
      IREJ = 1
c     default return point is kinematic sampling stage
      LXBAD = 3

c     remember initial setup
      NP_0    = NP
      SQS_0   = SQS
c     remember position on parton stack
      CALL GET_NPP(NPP_0,NPP0_0)

c     set interaction properties
c      IF(Irec.ne.1) CALL INI_EVENT(ECM,K_beam,Idm,Irec)

      IF(ndebug.gt.0)then
         IF(Irec.eq.0)THEN
            WRITE(LUN,*) 
     &           ' SIB_NDIFF: recursive call with (ecm,kb,kt,np,jdif):',
     &           ecm,k_beam,kt(1),(jdif(j),j=1,NW),NP
         ELSE
            WRITE(LUN,*)' SIB_NDIFF: regular call with (ECM,KB,NW,KT,',
     &           'JDIF,NP):',ecm,k_beam,NW,(kt(ii),ii=1,NW),
     &           (jdif(j),j=1,NW),NP
         ENDIF
      ENDIF
      
 2000 CONTINUE

c     reset parton stack
      CALL INI_PRTN_STCK(NPP_0,NPP0_0)

C...  sample multiple interaction configuration
      KBA = IABS(K_beam)
      L = LL(KBA)
      DO I=1,NW
        if(JDIF(I).eq.0) then
           CALL CUT_PRO(L, SQS, PTmin, NNSOF(I), NNJET(I))
        else
          NNSOF(I) = 1
          NNJET(I) = 0
        endif
c     add incoming target particles
        PZ = PAWT(SQS,AM(KBA),AM(KT(I)))
        E2 = SQRT(PZ**2+AM2(KT(I)))
        CALL ADD_PRTN(0.D0,0.D0,-PZ,E2,AM(KT(I)),KT(I),-2,0,IREFout)

c     add interactions
        xjdif = dble(jdif(I))
        xnjet = dble(nnjet(I))
        xnsof = dble(nnsof(I))
        CALL ADD_PRTN(xnsof,xnjet,xjdif,sqs,0.D0,I,-1,IREFout,IREF)
c     write parton stack index to interaction index
        IINTDX(I) = IREF
      ENDDO
c     remember state of parton stack
      CALL GET_NPP(NPP_1,NPP0_1)

C...  kinematic sampling stage

C...  sample x values
      ITRY(1) = 0
 3000 CONTINUE
      ITRY(1) = ITRY(1)+1
      IF(ITRY(1).GT.NREJ(1)) THEN 
c         NCALL = NCALL + 1
         GOTO 2000
      ENDIF
      NP = NP_0
      CALL INI_PRTN_STCK(NPP_1,NPP0_1)

      CALL SAMPLE_MINIJET(L,NW,NNJET,NNSOF,NJET,NSOF,x1jet,x2jet,lxbad)
      IF(LXBAD.eq.3)THEN
c     reject kinematics
         GOTO 3000
      ELSEIF(LXBAD.eq.2)THEN
c     reject kinematics and event structure
c         NCALL = NCALL + 1
         GOTO 2000
      ELSEIF(LXBAD.eq.1)THEN
c     reject entire event
         if(Ndebug.gt.0) 
     &        WRITE(LUN,*)' SIB_NDIFF: minijet rejection (Ncall):',Ncall
c     restore initial state
         NP    = NP_0
         CALL INI_PRTN_STCK(NPP_0,NPP0_0)
         SQS   = SQS_0
         S     = SQS*SQS
         RETURN
      ENDIF

C...  Prepare 2*NW valence/sea color strings and/or remnant.

c     default return point, jump back to sampling interaction structure
c      LXBAD = 2
      CALl SAMPLE_RMNT(K_beam,NW,X1Jet,X2JET,Irec,LXBAD)
      IF(LXBAD.eq.3)THEN
c     reject kinematics
         GOTO 3000
      ELSEIF(LXBAD.eq.2)THEN         
c     reject kinematics and event structure
c         NCALL = NCALL + 1
         GOTO 2000
      ELSEIF(LXBAD.eq.1)THEN         
c     reject entire event
         if(Ndebug.gt.0) 
     &   WRITE(LUN,*)' SIB_NDIFF: rmnt rejection (Ncall,NW):',Ncall,NW
c     restore initial state
         NP    = NP_0
         CALL INI_PRTN_STCK(NPP_0,NPP0_0)
         SQS   = SQS_0
         S     = SQS*SQS
         RETURN
      ENDIF

C     Check parton final state..
      CALL GET_NPP(NPP,NPP0)
      CALL PPSUM(1,NPP,Esum,PXsum,PYsum,PZsum,NF)
      IF(ABS(Esum/(0.5D0*Ecm*DBLE(NW+1))-1.D0).GT.EPS3)THEN
         WRITE(LUN,*) ' SIB_NDIFF: energy not conserved! : ',Ncall
         WRITE(LUN,*) '  sqs_inp = ', Ecm, ' sqs_out = ', Esum
         CALL PRNT_PRTN_STCK
         WRITE(LUN,*) ' SIB_NDIFF: event rejected! ',
     &        'partons do not conserve energy'
         WRITE(LUN,*)' (Ncall,NW,NPP,NJET,NSOF):',Ncall,NW,NPP,NJET,NSOF
c     CALL SIB_REJECT('SIB_NDIFF       ')
c     restore initial state
         NP    = NP_0
         CALL INI_PRTN_STCK(NPP_0,NPP0_0)
         SQS   = SQS_0
         S     = SQS*SQS
         RETURN
      ENDIF
      IF(NDEBUG.gt.0) THEN
         IF(NDEBUG.gt.1) CALL PRNT_PRTN_STCK
         WRITE(LUN,*) ' SIB_NDIFF: entering fragmentation stage...'
      ENDIF

C...  Fragmentation stage
      nj = 0
      ns = 0
      nv = 0
      II = NPP0_0+1
      DO WHILE (II.gt.0)
c     default return point: reject event if fragmentation fails
         LXBAD = 1         
c     loop over level0 partons
         CALL ITR_LVL0_PRTN(II,JJ,LPID)
c     read interaction
         CALL RD_INT(jj,Idm,iiflag)

C...  Fragmentation of soft/hard sea color strings
         IF(LPID.eq.100)THEN
            nj = nj + 1
            ipflag = 100
            KINT = nj
            CALL FRAG_MINIJET(jj,LXBAD)
            IF(LXBAD.ne.0) RETURN

         ELSEIF(LPID.eq.10)THEN
            ns = ns + 1
            ipflag = 10
            KINT = ns
            CALL FRAG_MINIJET(jj,LXBAD)
            IF(LXBAD.ne.0) RETURN

C...  fragment 'valence' strings
         ELSEIF(LPID.eq.1)THEN
            nv = nv + 1
            KINT = nv
            ipflag = 1
            CALL FRAG_VLNCE(jj,LXBAD)
            IF(LXBAD.ne.0) RETURN

C...  fragment remnants
         ELSEIF(IABS(LPID).eq.2)THEN
            CALL EXCTDEC(JJ,LXBAD)
            IF(LXBAD.ne.0) RETURN

C...  fragment incoherent diffraction
         ELSEIF(LPID.eq.-10.or.LPID.eq.-20.or.LPID.eq.-30)THEN
            CALL FRAG_INCHRNT_DIFF(jj,lxbad)
            IF(LXBAD.ne.0) RETURN

         ENDIF
      ENDDO
      IREJ = 0
      
      END
C=======================================================================

      SUBROUTINE SAMPLE_RMNT(Kbeam,NW,X1JET,X2JET,Irec,LBAD)

C-----------------------------------------------------------------------
C     routine to sample remnants
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NW_max
      PARAMETER (NW_max = 20)
      INTEGER ITRY, NREJ
      COMMON /S_CNT/ ITRY(20), NREJ(20)
      
c     external type declarations
      DOUBLE PRECISION X1JET,X2JET
      DIMENSION X2JET(NW_max)
      INTEGER KBEAM,NW,IREC,LBAD

C     COMMONs
      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
C     parameters that represent: NW: max. number of wounded nucleons,
C     NS,NH: max. number of soft and hard interactions
c      PARAMETER (NW_max = 20)
C     The COMMON block /S_CHIST/ contains information about the
C     the structure of the  generated event:
C     NWD   = number of wounded nucleons
C     NJET = total number of hard interactions
C     NSOF = total number of soft interactions
C     NNSOF (1:NW) = number of soft pomeron cuts in each interaction
C     NNJET (1:NW) = number of minijets produced in each interaction 
C     JDIF(1:NW) = diffraction code 
C                  0 : non-diff,
C                  1 : beam-diff
C                  2 : target-diff
C                  3 : double-diff
      INTEGER NNSOF,NNJET,JDIF,NWD,NJET,NSOF
      COMMON /S_CHIST/ NNSOF(NW_max),NNJET(NW_max),
     &     JDIF(NW_max),NWD,NJET,NSOF

      INTEGER IRMNT,KRB,KRT
      DOUBLE PRECISION XRMASS,XRMEX
      COMMON /S_RMNT/ XRMASS(2),XRMEX(2),IRMNT(NW_max),KRB,KRT(NW_max)

      INTEGER ICHP,ISTR,IBAR
      COMMON /S_CHP/ ICHP(99), ISTR(99), IBAR(99)

      INTEGER IISO,ISPN
      COMMON /S_SPN/ IISO(99), ISPN(99)

      INTEGER ICHM
      COMMON /S_CHM/ ICHM(99)

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN

c     internals
      DOUBLE PRECISION PREM,PREM_NUC,R,R2,S_RNDM,FLVXCHG,ALPH
      INTEGER ITGRMNT,IBMRMNT,I,j,jj,K,NPPLD,NPP0LD,IBMRMNT_OLD,
     &     IBAD,IKBAD,KBM
      DIMENSION ITGRMNT(NW_max)
      SAVE
      DATA PREM /0.D0/ , PREM_NUC /0.D0/

      IF(Ndebug.gt.1) 
     &  WRITE(LUN,*)' SAMPLE_RMNT: called with (Kbeam,NW,X1JET,',
     &     'X2JET,JDIF,Irec):',Kbeam,NW,X1JET,(X2JET(JJ),JJ=1,NW),
     &     (JDIF(JJ),JJ=1,NW),Irec

      IF(Irec.eq.0.and.NW.ne.1)then
         WRITE(LUN,*)' SAMPLE_RMNT: recursive call inconsistent!'
         CALL SIB_REJECT('SAMPLE_RMNT     ')
      endif

c     default return point for remnant excitation routine:
c     beam and target sampling
      IBAD = 1

c     set trial counter
      ITRY(2) = 0
c     remember position on parton stack
      CALL GET_NPP(NPPLD,NPP0LD)

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
      IF(IPAR(79).ne.0) PREM = 1.D0-PREM     

c     prob. of remnant excitation target side
      IF(IPAR(79).ne.0) PREM_NUC = 1.D0-PAR(23)
      IF(IPAR(63).eq.1) PREM_NUC = PREM_NUC/dble(NW)

c     turn of remnant for Nw>1
      IF(IPAR(77).eq.1)THEN
c     only beamside
         IF(NW.gt.1) PREM = 0
      ELSEIF(IPAR(77).eq.2)THEN
c     target and beam-side
         IF(NW.gt.1) then
            PREM = 0.D0
            PREM_NUC = 0.D0
         endif
      ELSE
         CONTINUE
      ENDIF

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
         IF(IPAR(79).eq.2.and.K.gt.1) PREM=1.D0-PAR(23)*PAR(128)
c     penalize remnant survival for multiple parton interactions
         IF(IPAR(80).ne.0) THEN
c     multiple interaction penalty for remnant survival, individual interaction
            ALPH = 1.D0+PAR(129)*DBLE(NNSOF(K)+NNJET(K)-1)
            PREM = 1.D0-(1.D0-PREM)**ALPH
            PREM_NUC = 1.D0-(1.D0-PREM_NUC)**ALPH
         ENDIF
         IF(JDIF(K).eq.0)THEN
            R = S_RNDM(k)
            R2 = S_RNDM(0)
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
     &        'SAMPLE_RMNT: (JW,PREM,NS,NH,IBMRMNT,LTGRMNT):',
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
         IBMRMNT=IABS(IBMRMNT-1)*Irec
      ENDIF
      IF(Ndebug.gt.1) 
     &     WRITE(LUN,*)
     &     ' SAMPLE_RMNT: remnant sampling (PREM,NW,LBMRMNT,LTGRMNT): ',
     &     PREM,NW,IBMRMNT,(ITGRMNT(j),j=1,NW)

      IBMRMNT_OLD = IBMRMNT

C...  Sample flavor and momentum fractions
 20   ITRY(2) = ITRY(2) + 1
c     reset parton stack
      CALL INI_PRTN_STCK(NPPLD,NPP0LD)
      IBMRMNT = IBMRMNT_OLD

c     retry without counting
c 22   CONTINUE
      IF(ITRY(2).gt.NREJ(2))THEN
         LBAD = 2
         IF(ndebug.gt.1)then 
            WRITE(LUN,*)' SAMPLE_RMNT: number of trials exceeded'
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
      CALL SAMPLE_BEAM(Kbm,NW,flvXchg,IBMRMNT,X1JET,IKBAD)
      IF(IKBAD.eq.1)THEN
c     resample minijets event
         LBAD = 3
         RETURN
      ELSEIF(IKBAD.eq.2)THEN
c     too many partons, reject NW, i.e. entire event
         LBAD = 1
         RETURN
      ENDIF

c     remnants and sea on target side
      CALL SAMPLE_TARGET(NW,flvXchg,ITGRMNT,X2JET,Irec,IKBAD)
      IF(IKBAD.eq.1)THEN
c     resample minijets event
         LBAD = 3
         RETURN
      ELSEIF(IKBAD.eq.2)THEN
c     too many partons, reject NW, i.e. entire event
         LBAD = 1
         RETURN
      ENDIF

C...  sample remnant excitation masses and add to parton stack
c     beam-side (one remnant, formed by several interactions)
c     target-side (possibly NW remnants)

      DO I=1,NW
c     default return point
         IBAD = 1
         IF(IPAR(78).EQ.1)THEN
c$$$            write(lun,*) 
c$$$     &           ' SIB_RMNT: multiple excitation model',
c$$$     &           ' not implemented yet!'
c$$$            stop
c     model where beam side remnant can receive mass from multiple target nucleons
            IF(IBMRMNT.gt.0)THEN
c     beam side remnant excited               
               if(ITGRMNT(I).eq.0)then
                  CALL EXCT_RMNT(I,1,IBAD)
               else
                  CALL EXCT_RMNT(I,3,IBAD)
               endif
               IBMRMNT = IBMRMNT - 1
            ELSE
c     beam side remnant not excited
               if(ITGRMNT(I).ne.0)then
                  CALL EXCT_RMNT(I,2,IBAD)
               else
                  CALL EXCT_RMNT(I,0,IBAD)
               endif
            ENDIF

         ELSEIF(IPAR(78).eq.2)then
            IF(IBMRMNT.gt.0)then
c     beam side remnant excited, only once!
               IF(ITGRMNT(I).eq.0)then
                  CALL EXCT_RMNT(I,1,IBAD)
               else
                  CALL EXCT_RMNT(I,3,IBAD)
               endif
               IBMRMNT = 0
            ELSE
c     beam side remnant not excited
               IF(ITGRMNT(I).ne.0)then
                  CALL EXCT_RMNT(I,2,IBAD)
               else
                  CALL EXCT_RMNT(I,0,IBAD)
               endif
            ENDIF
         ELSE
c     no remnant model
            CALL EXCT_RMNT(I,0,IBAD)
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
C=======================================================================

      SUBROUTINE SIB_HADCSL(L,ECM,SIGTOT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO)

C-----------------------------------------------------------------------
C     low-energy cross section parametrizations (target always proton)
C
C     input:   L           beam particle: (1 - proton,
C                                          2 - pion,
C                                          3 - kaon)
C                          target is always proton
C              ECM         c.m. energy (GeV)
C
C     output:  SIGTOT      total cross section (mb)
C              SIGEL       elastic cross section (mb)
C              SIGDIF      diffractive cross section (sd-1,sd-2,dd, mb)
C              SLOPE       forward elastic slope (GeV**-2)
C              RHO         real/imaginary part of elastic amplitude
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      DIMENSION SIGDIF(3)

      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      SAVE

C  proton-proton cross section as reference
      CALL SIB_HADCS1(1,ECM,SIGTOT,SIGEL,SIGINEL,SLOPE,RHO)

C  parametrization for diffraction
      Xi_min = 1.5D0/(ECM*ECM)
      Xi_max = PAR(13)
      SIGeff = SIGEL
      CALL SIB_HADCS2(ECM,Xi_min,Xi_max,SIGeff,SIGDIF)

      if(L.eq.1) return

C  regge motivated rescaling of diffraction dissociation
      sigtot_pp = SIGTOT
      sigel_pp  = SIGEL
      slope_pp  = SLOPE
      CALL SIB_HADCS1(L,ECM,SIGTOT,SIGEL,SIGINEL,SLOPE,RHO)
      SIGDIF(1) = slope_pp/SLOPE*SIGTOT/sigtot_pp*SIGDIF(1)
      SIGDIF(2) = slope_pp/SLOPE*SIGEL/sigel_pp*SIGDIF(2)
      SIGDIF(3) = SIGTOT/sigtot_pp*SIGDIF(3)

      END

C=======================================================================

      SUBROUTINE SIB_HADCS1(L,ECM,SIGTOT,SIGEL,SIGINEL,SLOPE,RHO)

C-----------------------------------------------------------------------
C     low-energy cross section parametrizations
C
C     input:   L           beam particle: (1 - proton,
C                                          2 - pion,
C                                          3 - kaon)
C                          target is always proton
C              ECM         c.m. energy (GeV)
C
C     output:  SIGTOT      total cross section (mb)
C              SIGEL       elastic cross section (mb)
C              SIGDIF      diffractive cross section (sd-1,sd-2,dd, mb)
C              SLOPE       forward elastic slope (GeV**-2)
C              RHO         real/imaginary part of elastic amplitude
C
C     comments:
C     - low-energy data interpolation uses PDG fits from 1992
C     - slopes from ???, new fit to pp data
C     - high-energy extrapolation by Donnachie-Landshoff like fit made
C       by PDG 1996
C     - analytic extension of amplitude to calculate rho
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)
      DIMENSION TPDG92(7,2,6),TPDG96(9,6),BURQ83(3,6),XMA(6)
      SAVE

      DATA TPDG92  /
     &  3.D0, 2100.D0, 48.D0, 0.D0, 1.D0, 0.522D0, -4.51D0,
     &  3.D0, 2100.D0, 11.9D0, 26.9D0, -1.21D0, 0.169D0, -1.85D0,
     &  5.D0, 2100.D0, 38.4D0, 77.6D0, -0.64D0, 0.26D0, -1.2D0,
     &  5.D0, 2100.D0, 10.2D0, 52.7D0, -1.16D0, 0.125D0, -1.28D0,
     &  4.D0, 340.D0,  16.4D0, 19.3D0, -0.42D0, 0.19D0, 0.D0,
     &  4.D0, 340.D0,  0.D0, 11.4D0, -0.4D0, 0.079D0, 0.D0,
     &  2.5D0, 370.D0, 33.D0, 14.D0, -1.36D0, 0.456D0, -4.03D0,
     &  2.5D0, 370.D0, 1.76D0, 11.2D0, -0.64D0, 0.043D0, 0.D0,
     &  2.D0, 310.D0,  18.1D0, 0.D0, 1.D0, 0.26D0, -1.D0,
     &  2.D0, 310.D0,  5.D0, 8.1D0, -1.8D0, 0.16D0, -1.3D0,
     &  3.D0, 310.D0,  32.1D0, 0.D0, 1.D0, 0.66D0, -5.6D0,
     &  3.D0, 310.D0,  7.3D0, 0.D0, 1.D0, 0.29D0, -2.4D0  /

      DATA TPDG96  /
     &  50.D0, 22.D0,0.079D0,0.25D0,0.D0,
     &         77.15D0,-21.05D0,0.46D0,0.9D0,
     &  50.D0, 22.D0,0.079D0,0.25D0,0.D0,
     &         77.15D0,21.05D0,0.46D0,0.9D0,
     &  10.D0, 13.70D0,0.079D0,0.25D0,0.D0,
     &         31.85D0,-4.05D0,0.45D0,0.9D0,
     &  10.D0, 13.70D0,0.079D0,0.25D0,0.D0,
     &         31.85D0,4.05D0,0.45D0,0.9D0,
     &  10.D0, 12.20D0,0.079D0,0.25D0,0.D0,
     &         17.35D0,-9.05D0,0.50D0,0.9D0,
     &  10.D0, 12.20D0,0.079D0,0.25D0,0.D0,
     &         17.35D0,9.05D0,0.50D0,0.9D0  /

      DATA BURQ83 /
     &  8.557D0,  0.00D0, 0.574D0,
     &  11.13D0,  7.23D0, 0.30D0,
     &  9.11D0,  -0.73D0, 0.28D0,
     &  9.11D0,   0.65D0, 0.28D0,
     &  8.55D0,  -5.98D0, 0.28D0,
     &  8.55D0,   1.60D0, 0.28D0  /

c     DATA XMA / 2*0.93956563D0, 2*0.13956995D0, 2*0.493677D0 /
      DATA GEV2MB /0.389365D0/
      DATA INIT/0/

      IF(INIT.EQ.0) THEN
c  use the internal masses 
        XMA(1) =  AM(13)   ! proton
        XMA(2) =  AM(14)   ! neutron
        XMA(3) =  AM(7)    ! pi+
        XMA(4) =  AM(8)    ! pi-
        XMA(5) =  AM(9)    ! K+
        XMA(6) =  AM(10)   ! K-
        INIT = 1
      ENDIF

C  find index
      IF    (L.eq.1) THEN
        K = 1                            ! p p
      ELSEIF(L.eq.2) THEN
        K = 3                            ! pi+ p
*       K = 4                            ! pi- p
      ELSEIF(L.eq.3) THEN
        K = 5                            ! K+ p
*       K = 6                            ! K- p
      ELSE
        GOTO 100
      ENDIF

C  calculate lab momentum
      SS = ECM**2
      E1 = (SS-XMA(1)**2-XMA(K)**2)/(2.D0*XMA(1))
      PL = dSQRT((E1-XMA(K))*(E1+XMA(K)))
      PLL = dLOG(PL)

C  check against lower limit
      IF(ECM.LE.XMA(1)+XMA(K)) GOTO 200

      XP  = TPDG96(2,K)*SS**TPDG96(3,K)
      YP  = TPDG96(6,K)/SS**TPDG96(8,K)
      YM  = TPDG96(7,K)/SS**TPDG96(8,K)

      PHR = dTAN(PI/2.D0*(1.D0-TPDG96(8,K)))
      PHP = dTAN(PI/2.D0*(1.D0+TPDG96(3,K)))
      RHO = (-YP/PHR + YM*PHR - XP/PHP)/(YP+YM+XP)

      SLOPE = BURQ83(1,K)+BURQ83(2,K)/dSQRT(PL)+BURQ83(3,K)*PLL

C  select energy range and interpolation method
      IF(PL.LT.TPDG96(1,K)) THEN
        SIGTOT = TPDG92(3,1,K)+TPDG92(4,1,K)*PL**TPDG92(5,1,K)
     &          + TPDG92(6,1,K)*PLL**2+TPDG92(7,1,K)*PLL
        SIGEL  = TPDG92(3,2,K)+TPDG92(4,2,K)*PL**TPDG92(5,2,K)
     &          + TPDG92(6,2,K)*PLL**2+TPDG92(7,2,K)*PLL
      ELSE IF(PL.LT.TPDG92(2,1,K)) THEN
        SIGTO1 = TPDG92(3,1,K)+TPDG92(4,1,K)*PL**TPDG92(5,1,K)
     &          + TPDG92(6,1,K)*PLL**2+TPDG92(7,1,K)*PLL
        SIGEL1 = TPDG92(3,2,K)+TPDG92(4,2,K)*PL**TPDG92(5,2,K)
     &          + TPDG92(6,2,K)*PLL**2+TPDG92(7,2,K)*PLL
        SIGTO2 = YP+YM+XP
        SIGEL2 = SIGTO2**2/(16.D0*PI*SLOPE*GEV2MB)*(1.D0+RHO**2)
        X2 = dLOG(PL/TPDG96(1,K))/dLOG(TPDG92(2,1,K)/TPDG96(1,K))
        X1 = 1.D0 - X2
        SIGTOT = SIGTO2*X2 + SIGTO1*X1
        SIGEL  = SIGEL2*X2 + SIGEL1*X1
      ELSE
        SIGTOT = YP+YM+XP
        SIGEL  = SIGTOT**2/(16.D0*PI*SLOPE*GEV2MB)*(1.D0+RHO**2)
      ENDIF
      SIGINEL = SIGTOT-SIGEL

      RETURN

 100  CONTINUE
        WRITE(LUN,'(1X,2A,2I7)') ' SIB_HADCS1: ',
     &    'invalid beam particle: ',L
        RETURN

 200  CONTINUE
        WRITE(LUN,'(1X,2A,1P,E12.4)') ' SIB_HADCS1: ',
     &    'energy too small (Ecm): ',ECM

      END
C=======================================================================

      SUBROUTINE SIB_HADCS2(SQS,Xi_min,Xi_max,SIGeff,SIGDIF)

C-----------------------------------------------------------------------
C   cross section for diffraction dissociation 
C
C   - single diffraction dissociation:
C     Goulianos' parametrization (Ref: PL B358 (1995) 379)
C   - double diffration dissociation: simple scaling model using 
C     single diff. cross section
C
C     in addition rescaling for different particles is applied using
C     internal rescaling tables (not implemented yet)
C
C     input:     SQS         c.m. energy (GeV)
C                Xi_min      min. diff mass (squared) = Xi_min*SQS**2
C                Xi_max      max. diff mass (squared) = Xi_max*SQS**2
C                SIGeff      effective cross section for DD scaling
C
C     output:    sig_sd1     cross section for diss. of particle 1 (mb)
C                sig_sd2     cross section for diss. of particle 2 (mb)
C                sig_dd      cross section for diss. of both particles
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      DIMENSION SIGDIF(3)
      DIMENSION Xpos1(96),Xwgh1(96),Xpos2(96),Xwgh2(96)
      DOUBLE PRECISION xil,xiu,tl,tu
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN
      SAVE

C  model parameters
      DATA delta    / 0.104D0 /
      DATA alphap   / 0.25D0 /
      DATA beta0    / 6.56D0 /
      DATA gpom0    / 1.21D0 /
      DATA xm_p     / 0.938D0 /
      DATA x_rad2   / 0.71D0 /

C  integration precision
      DATA Ngau1    / 32 /
      DATA Ngau2    / 32 /

      DATA GEV2MB /0.389365D0/

      SIGDIF(1) = 0.D0
      SIGDIF(2) = 0.D0
      SIGDIF(3) = 0.D0

      XIL = dLOG(Xi_min)
      XIU = dLOG(Xi_max)

      if(XIL.ge.XIU) return

      SS = SQS*SQS
      xm4_p2 = 4.D0*xm_p**2
      fac = beta0**2/(16.D0*PI)

      t1 = -5.D0
      t2 = 0.D0
      tl = x_rad2/3.D0/(1.D0-t1/x_rad2)**3
      tu = x_rad2/3.D0/(1.D0-t2/x_rad2)**3

C  flux renormalization and cross section for pp/ppbar case

      Xnorm  = 0.D0

      xil = dlog(1.5D0/SS)
      xiu = dlog(0.1D0)

      IF(xiu.LE.xil) goto 1000

      CALL SIB_GAUSET(xil,xiu,Ngau1,xpos1,xwgh1)
      CALL SIB_GAUSET(tl,tu,Ngau2,xpos2,xwgh2)

      do i1=1,Ngau1

        xi = dexp(xpos1(i1))
        w_xi = Xwgh1(i1)

        do i2=1,Ngau2

          tt = x_rad2-x_rad2*(x_rad2/(3.D0*xpos2(i2)))**(1.D0/3.D0)

          alpha_t =  1.D0+delta+alphap*tt
          f2_t = ((xm4_p2-2.8D0*tt)/(xm4_p2-tt))**2
            
          Xnorm = Xnorm
     &      + f2_t*xi**(2.D0-2.D0*alpha_t)*Xwgh2(i2)*w_xi

        enddo
      enddo   

      Xnorm = Xnorm*fac

 1000 continue

      XIL = dLOG(Xi_min)
      XIU = dLOG(Xi_max)

      T1 = -5.D0
      T2 = 0.D0

      TL = x_rad2/3.D0/(1.D0-t1/x_rad2)**3
      TU = x_rad2/3.D0/(1.D0-t2/x_rad2)**3

C  single diffraction diss. cross section 

      CSdiff = 0.D0

      CALL SIB_GAUSET(XIL,XIU,NGAU1,XPOS1,XWGH1)
      CALL SIB_GAUSET(TL,TU,NGAU2,XPOS2,XWGH2)

      do i1=1,Ngau1

        xi = dexp(xpos1(i1))
        w_xi = Xwgh1(i1)*beta0*gpom0*(xi*ss)**delta

        do i2=1,Ngau2

          tt = x_rad2-x_rad2*(x_rad2/(3.D0*xpos2(i2)))**(1.D0/3.D0)

          alpha_t =  1.D0+delta+alphap*tt
          f2_t = ((xm4_p2-2.8D0*tt)/(xm4_p2-tt))**2

          CSdiff = CSdiff 
     &      + f2_t*xi**(2.D0-2.D0*alpha_t)*Xwgh2(i2)*w_xi

        enddo
      enddo

      CSdiff = CSdiff*fac*GEV2MB/MAX(1.D0,Xnorm)

*     write(LUN,'(1x,1p,4e14.3)') 
*    &  sqrt(SS),Xnorm,2.d0*CSdiff*MAX(1.d0,Xnorm),2.d0*CSdiff

      SIGDIF(1) = CSdiff
      SIGDIF(2) = CSdiff

C  double diff. dissociation from simple probability consideration
*     Pdiff = 0.5d0-sqrt(0.25d0-CSdiff/SIGeff)
      Pdiff = CSdiff/SIGeff
      SIGDIF(3) = Pdiff*Pdiff*SIGeff

      END
C=======================================================================

      SUBROUTINE SIB_GAUSET(AX,BX,NX,Z,W)

C-----------------------------------------------------------------------
C
C     N-point gauss zeros and weights for the interval (AX,BX) are
C           stored in  arrays Z and W respectively.
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
C
      COMMON /GQCOM/A(273),X(273),KTAB(96)
      DIMENSION Z(NX),W(NX)
      SAVE
      DATA INIT/0/
C
      ALPHA=0.5D0*(BX+AX)
      BETA=0.5D0*(BX-AX)
      N=NX
*
*  the N=1 case:
      IF(N.NE.1) GO TO 1
      Z(1)=ALPHA
      W(1)=BX-AX
      RETURN
*
*  the Gauss cases:
    1 IF((N.LE.16).AND.(N.GT.1)) GO TO 2
      IF(N.EQ.20) GO TO 2
      IF(N.EQ.24) GO TO 2
      IF(N.EQ.32) GO TO 2
      IF(N.EQ.40) GO TO 2
      IF(N.EQ.48) GO TO 2
      IF(N.EQ.64) GO TO 2
      IF(N.EQ.80) GO TO 2
      IF(N.EQ.96) GO TO 2
*
*  the extended Gauss cases:
      IF((N/96)*96.EQ.N) GO TO 3
*
C  jump to center of intervall intrgration:
      GO TO 100
*
C  get Gauss point array
*
    2 CALL PO106BD
C     -print out message
*     IF(INIT.LE.20)THEN
*       INIT=init+1
*       WRITE (6,*) ' initialization of Gauss int. N=',N
*     ENDIF
C  extract real points
      K=KTAB(N)
      M=N/2
      DO 21 J=1,M
C       extract values from big array
        JTAB=K-1+J
        WTEMP=BETA*A(JTAB)
        DELTA=BETA*X(JTAB)
C       store them backward
        Z(J)=ALPHA-DELTA
        W(J)=WTEMP
C       store them forward
        JP=N+1-J
        Z(JP)=ALPHA+DELTA
        W(JP)=WTEMP
   21 CONTINUE
C     store central point (odd N)
      IF((N-M-M).EQ.0) RETURN
      Z(M+1)=ALPHA
      JMID=K+M
      W(M+1)=BETA*A(JMID)
      RETURN
C
C  get ND96 times chained 96 Gauss point array
C
    3 CALL PO106BD
C  print out message
      IF(INIT.LE.20)THEN
        INIT=init+1
        WRITE (6,*) ' initialization of extended Gauss int. N=',N
      ENDIF
C     -extract real points
      K=KTAB(96)
      ND96=N/96
      DO 31 J=1,48
C       extract values from big array
        JTAB=K-1+J
        WTEMP=BETA*A(JTAB)
        DELTA=BETA*X(JTAB)
        WTeMP=WTEMP/ND96
        DeLTA=DELTA/ND96
        DO 32 JD96=0,ND96-1
          ZCNTR= (ALPHA-BETA)+ BETA*DBLE(2*JD96+1)/DBLE(ND96)
C         store them backward
          Z(J+JD96*96)=ZCNTR-DELTA
          W(J+JD96*96)=WTEMP
C         store them forward
          JP=96+1-J
          Z(JP+JD96*96)=ZCNTR+DELTA
          W(JP+JD96*96)=WTEMP
   32   CONTINUE
   31 CONTINUE
      RETURN
*
C  the center of intervall cases:
  100 CONTINUE
C  print out message
      IF(INIT.LE.20)THEN
        INIT=init+1
        WRITE (6,*) ' init. of center of intervall int. N=',N
      ENDIF
C  put in constant weight and equally spaced central points
      N=IABS(N)
      DO 111 IN=1,N
        WIN=(BX-AX)/DBLE(N)
        Z(IN)=AX  + (DBLE(IN)-.5D0)*WIN
  111 W(IN)=WIN
      RETURN
      END

C=======================================================================

      SUBROUTINE PO106BD

C-----------------------------------------------------------------------
C
C     store big arrays needed for Gauss integral, CERNLIB D106BD
C     (arrays A,X,ITAB copied on B,Y,LTAB)
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
C
      COMMON /GQCOM/ B(273),Y(273),LTAB(96)
      DIMENSION      A(273),X(273),KTAB(96)
      SAVE
C
C-----TABLE OF INITIAL SUBSCRIPTS FOR N=2(1)16(4)96
      DATA KTAB(2)/1/
      DATA KTAB(3)/2/
      DATA KTAB(4)/4/
      DATA KTAB(5)/6/
      DATA KTAB(6)/9/
      DATA KTAB(7)/12/
      DATA KTAB(8)/16/
      DATA KTAB(9)/20/
      DATA KTAB(10)/25/
      DATA KTAB(11)/30/
      DATA KTAB(12)/36/
      DATA KTAB(13)/42/
      DATA KTAB(14)/49/
      DATA KTAB(15)/56/
      DATA KTAB(16)/64/
      DATA KTAB(20)/72/
      DATA KTAB(24)/82/
      DATA KTAB(28)/82/
      DATA KTAB(32)/94/
      DATA KTAB(36)/94/
      DATA KTAB(40)/110/
      DATA KTAB(44)/110/
      DATA KTAB(48)/130/
      DATA KTAB(52)/130/
      DATA KTAB(56)/130/
      DATA KTAB(60)/130/
      DATA KTAB(64)/154/
      DATA KTAB(68)/154/
      DATA KTAB(72)/154/
      DATA KTAB(76)/154/
      DATA KTAB(80)/186/
      DATA KTAB(84)/186/
      DATA KTAB(88)/186/
      DATA KTAB(92)/186/
      DATA KTAB(96)/226/
C
C-----TABLE OF ABSCISSAE (X) AND WEIGHTS (A) FOR INTERVAL (-1,+1).
C
C-----N=2
      DATA X(1)/0.577350269189626D0  /, A(1)/1.000000000000000D0  /
C-----N=3
      DATA X(2)/0.774596669241483D0  /, A(2)/0.555555555555556D0  /
      DATA X(3)/0.000000000000000D0  /, A(3)/0.888888888888889D0  /
C-----N=4
      DATA X(4)/0.861136311594053D0  /, A(4)/0.347854845137454D0  /
      DATA X(5)/0.339981043584856D0  /, A(5)/0.652145154862546D0  /
C-----N=5
      DATA X(6)/0.906179845938664D0  /, A(6)/0.236926885056189D0  /
      DATA X(7)/0.538469310105683D0  /, A(7)/0.478628670499366D0  /
      DATA X(8)/0.000000000000000D0  /, A(8)/0.568888888888889D0  /
C-----N=6
      DATA X(9)/0.932469514203152D0  /, A(9)/0.171324492379170D0  /
      DATA X(10)/0.661209386466265D0 /, A(10)/0.360761573048139D0 /
      DATA X(11)/0.238619186083197D0 /, A(11)/0.467913934572691D0 /
C-----N=7
      DATA X(12)/0.949107912342759D0 /, A(12)/0.129484966168870D0 /
      DATA X(13)/0.741531185599394D0 /, A(13)/0.279705391489277D0 /
      DATA X(14)/0.405845151377397D0 /, A(14)/0.381830050505119D0 /
      DATA X(15)/0.000000000000000D0 /, A(15)/0.417959183673469D0 /
C-----N=8
      DATA X(16)/0.960289856497536D0 /, A(16)/0.101228536290376D0 /
      DATA X(17)/0.796666477413627D0 /, A(17)/0.222381034453374D0 /
      DATA X(18)/0.525532409916329D0 /, A(18)/0.313706645877887D0 /
      DATA X(19)/0.183434642495650D0 /, A(19)/0.362683783378362D0 /
C-----N=9
      DATA X(20)/0.968160239507626D0 /, A(20)/0.081274388361574D0 /
      DATA X(21)/0.836031107326636D0 /, A(21)/0.180648160694857D0 /
      DATA X(22)/0.613371432700590D0 /, A(22)/0.260610696402935D0 /
      DATA X(23)/0.324253423403809D0 /, A(23)/0.312347077040003D0 /
      DATA X(24)/0.000000000000000D0 /, A(24)/0.330239355001260D0 /
C-----N=10
      DATA X(25)/0.973906528517172D0 /, A(25)/0.066671344308688D0 /
      DATA X(26)/0.865063366688985D0 /, A(26)/0.149451349150581D0 /
      DATA X(27)/0.679409568299024D0 /, A(27)/0.219086362515982D0 /
      DATA X(28)/0.433395394129247D0 /, A(28)/0.269266719309996D0 /
      DATA X(29)/0.148874338981631D0 /, A(29)/0.295524224714753D0 /
C-----N=11
      DATA X(30)/0.978228658146057D0 /, A(30)/0.055668567116174D0 /
      DATA X(31)/0.887062599768095D0 /, A(31)/0.125580369464905D0 /
      DATA X(32)/0.730152005574049D0 /, A(32)/0.186290210927734D0 /
      DATA X(33)/0.519096129206812D0 /, A(33)/0.233193764591990D0 /
      DATA X(34)/0.269543155952345D0 /, A(34)/0.262804544510247D0 /
      DATA X(35)/0.000000000000000D0 /, A(35)/0.272925086777901D0 /
C-----N=12
      DATA X(36)/0.981560634246719D0 /, A(36)/0.047175336386512D0 /
      DATA X(37)/0.904117256370475D0 /, A(37)/0.106939325995318D0 /
      DATA X(38)/0.769902674194305D0 /, A(38)/0.160078328543346D0 /
      DATA X(39)/0.587317954286617D0 /, A(39)/0.203167426723066D0 /
      DATA X(40)/0.367831498998180D0 /, A(40)/0.233492536538355D0 /
      DATA X(41)/0.125233408511469D0 /, A(41)/0.249147045813403D0 /
C-----N=13
      DATA X(42)/0.984183054718588D0 /, A(42)/0.040484004765316D0 /
      DATA X(43)/0.917598399222978D0 /, A(43)/0.092121499837728D0 /
      DATA X(44)/0.801578090733310D0 /, A(44)/0.138873510219787D0 /
      DATA X(45)/0.642349339440340D0 /, A(45)/0.178145980761946D0 /
      DATA X(46)/0.448492751036447D0 /, A(46)/0.207816047536889D0 /
      DATA X(47)/0.230458315955135D0 /, A(47)/0.226283180262897D0 /
      DATA X(48)/0.000000000000000D0 /, A(48)/0.232551553230874D0 /
C-----N=14
      DATA X(49)/0.986283808696812D0 /, A(49)/0.035119460331752D0 /
      DATA X(50)/0.928434883663574D0 /, A(50)/0.080158087159760D0 /
      DATA X(51)/0.827201315069765D0 /, A(51)/0.121518570687903D0 /
      DATA X(52)/0.687292904811685D0 /, A(52)/0.157203167158194D0 /
      DATA X(53)/0.515248636358154D0 /, A(53)/0.185538397477938D0 /
      DATA X(54)/0.319112368927890D0 /, A(54)/0.205198463721296D0 /
      DATA X(55)/0.108054948707344D0 /, A(55)/0.215263853463158D0 /
C-----N=15
      DATA X(56)/0.987992518020485D0 /, A(56)/0.030753241996117D0 /
      DATA X(57)/0.937273392400706D0 /, A(57)/0.070366047488108D0 /
      DATA X(58)/0.848206583410427D0 /, A(58)/0.107159220467172D0 /
      DATA X(59)/0.724417731360170D0 /, A(59)/0.139570677926154D0 /
      DATA X(60)/0.570972172608539D0 /, A(60)/0.166269205816994D0 /
      DATA X(61)/0.394151347077563D0 /, A(61)/0.186161000015562D0 /
      DATA X(62)/0.201194093997435D0 /, A(62)/0.198431485327111D0 /
      DATA X(63)/0.000000000000000D0 /, A(63)/0.202578241925561D0 /
C-----N=16
      DATA X(64)/0.989400934991650D0 /, A(64)/0.027152459411754D0 /
      DATA X(65)/0.944575023073233D0 /, A(65)/0.062253523938648D0 /
      DATA X(66)/0.865631202387832D0 /, A(66)/0.095158511682493D0 /
      DATA X(67)/0.755404408355003D0 /, A(67)/0.124628971255534D0 /
      DATA X(68)/0.617876244402644D0 /, A(68)/0.149595988816577D0 /
      DATA X(69)/0.458016777657227D0 /, A(69)/0.169156519395003D0 /
      DATA X(70)/0.281603550779259D0 /, A(70)/0.182603415044924D0 /
      DATA X(71)/0.095012509837637D0 /, A(71)/0.189450610455069D0 /
C-----N=20
      DATA X(72)/0.993128599185094D0 /, A(72)/0.017614007139152D0 /
      DATA X(73)/0.963971927277913D0 /, A(73)/0.040601429800386D0 /
      DATA X(74)/0.912234428251325D0 /, A(74)/0.062672048334109D0 /
      DATA X(75)/0.839116971822218D0 /, A(75)/0.083276741576704D0 /
      DATA X(76)/0.746331906460150D0 /, A(76)/0.101930119817240D0 /
      DATA X(77)/0.636053680726515D0 /, A(77)/0.118194531961518D0 /
      DATA X(78)/0.510867001950827D0 /, A(78)/0.131688638449176D0 /
      DATA X(79)/0.373706088715419D0 /, A(79)/0.142096109318382D0 /
      DATA X(80)/0.227785851141645D0 /, A(80)/0.149172986472603D0 /
      DATA X(81)/0.076526521133497D0 /, A(81)/0.152753387130725D0 /
C-----N=24
      DATA X(82)/0.995187219997021D0 /, A(82)/0.012341229799987D0 /
      DATA X(83)/0.974728555971309D0 /, A(83)/0.028531388628933D0 /
      DATA X(84)/0.938274552002732D0 /, A(84)/0.044277438817419D0 /
      DATA X(85)/0.886415527004401D0 /, A(85)/0.059298584915436D0 /
      DATA X(86)/0.820001985973902D0 /, A(86)/0.073346481411080D0 /
      DATA X(87)/0.740124191578554D0 /, A(87)/0.086190161531953D0 /
      DATA X(88)/0.648093651936975D0 /, A(88)/0.097618652104113D0 /
      DATA X(89)/0.545421471388839D0 /, A(89)/0.107444270115965D0 /
      DATA X(90)/0.433793507626045D0 /, A(90)/0.115505668053725D0 /
      DATA X(91)/0.315042679696163D0 /, A(91)/0.121670472927803D0 /
      DATA X(92)/0.191118867473616D0 /, A(92)/0.125837456346828D0 /
      DATA X(93)/0.064056892862605D0 /, A(93)/0.127938195346752D0 /
C-----N=32
      DATA X(94)/0.997263861849481D0 /, A(94)/0.007018610009470D0 /
      DATA X(95)/0.985611511545268D0 /, A(95)/0.016274394730905D0 /
      DATA X(96)/0.964762255587506D0 /, A(96)/0.025392065309262D0 /
      DATA X(97)/0.934906075937739D0 /, A(97)/0.034273862913021D0 /
      DATA X(98)/0.896321155766052D0 /, A(98)/0.042835898022226D0 /
      DATA X(99)/0.849367613732569D0 /, A(99)/0.050998059262376D0 /
      DATA X(100)/0.794483795967942D0/, A(100)/0.058684093478535D0/
      DATA X(101)/0.732182118740289D0/, A(101)/0.065822222776361D0/
      DATA X(102)/0.663044266930215D0/, A(102)/0.072345794108848D0/
      DATA X(103)/0.587715757240762D0/, A(103)/0.078193895787070D0/
      DATA X(104)/0.506899908932229D0/, A(104)/0.083311924226946D0/
      DATA X(105)/0.421351276130635D0/, A(105)/0.087652093004403D0/
      DATA X(106)/0.331868602282127D0/, A(106)/0.091173878695763D0/
      DATA X(107)/0.239287362252137D0/, A(107)/0.093844399080804D0/
      DATA X(108)/0.144471961582796D0/, A(108)/0.095638720079274D0/
      DATA X(109)/0.048307665687738D0/, A(109)/0.096540088514727D0/
C-----N=40
      DATA X(110)/0.998237709710559D0/, A(110)/0.004521277098533D0/
      DATA X(111)/0.990726238699457D0/, A(111)/0.010498284531152D0/
      DATA X(112)/0.977259949983774D0/, A(112)/0.016421058381907D0/
      DATA X(113)/0.957916819213791D0/, A(113)/0.022245849194166D0/
      DATA X(114)/0.932812808278676D0/, A(114)/0.027937006980023D0/
      DATA X(115)/0.902098806968874D0/, A(115)/0.033460195282547D0/
      DATA X(116)/0.865959503212259D0/, A(116)/0.038782167974472D0/
      DATA X(117)/0.824612230833311D0/, A(117)/0.043870908185673D0/
      DATA X(118)/0.778305651426519D0/, A(118)/0.048695807635072D0/
      DATA X(119)/0.727318255189927D0/, A(119)/0.053227846983936D0/
      DATA X(120)/0.671956684614179D0/, A(120)/0.057439769099391D0/
      DATA X(121)/0.612553889667980D0/, A(121)/0.061306242492928D0/
      DATA X(122)/0.549467125095128D0/, A(122)/0.064804013456601D0/
      DATA X(123)/0.483075801686178D0/, A(123)/0.067912045815233D0/
      DATA X(124)/0.413779204371605D0/, A(124)/0.070611647391286D0/
      DATA X(125)/0.341994090825758D0/, A(125)/0.072886582395804D0/
      DATA X(126)/0.268152185007253D0/, A(126)/0.074723169057968D0/
      DATA X(127)/0.192697580701371D0/, A(127)/0.076110361900626D0/
      DATA X(128)/0.116084070675255D0/, A(128)/0.077039818164247D0/
      DATA X(129)/0.038772417506050D0/, A(129)/0.077505947978424D0/
C-----N=48
      DATA X(130)/0.998771007252426D0/, A(130)/0.003153346052305D0/
      DATA X(131)/0.993530172266350D0/, A(131)/0.007327553901276D0/
      DATA X(132)/0.984124583722826D0/, A(132)/0.011477234579234D0/
      DATA X(133)/0.970591592546247D0/, A(133)/0.015579315722943D0/
      DATA X(134)/0.952987703160430D0/, A(134)/0.019616160457355D0/
      DATA X(135)/0.931386690706554D0/, A(135)/0.023570760839324D0/
      DATA X(136)/0.905879136715569D0/, A(136)/0.027426509708356D0/
      DATA X(137)/0.876572020274247D0/, A(137)/0.031167227832798D0/
      DATA X(138)/0.843588261624393D0/, A(138)/0.034777222564770D0/
      DATA X(139)/0.807066204029442D0/, A(139)/0.038241351065830D0/
      DATA X(140)/0.767159032515740D0/, A(140)/0.041545082943464D0/
      DATA X(141)/0.724034130923814D0/, A(141)/0.044674560856694D0/
      DATA X(142)/0.677872379632663D0/, A(142)/0.047616658492490D0/
      DATA X(143)/0.628867396776513D0/, A(143)/0.050359035553854D0/
      DATA X(144)/0.577224726083972D0/, A(144)/0.052890189485193D0/
      DATA X(145)/0.523160974722233D0/, A(145)/0.055199503699984D0/
      DATA X(146)/0.466902904750958D0/, A(146)/0.057277292100403D0/
      DATA X(147)/0.408686481990716D0/, A(147)/0.059114839698395D0/
      DATA X(148)/0.348755886292160D0/, A(148)/0.060704439165893D0/
      DATA X(149)/0.287362487355455D0/, A(149)/0.062039423159892D0/
      DATA X(150)/0.224763790394689D0/, A(150)/0.063114192286254D0/
      DATA X(151)/0.161222356068891D0/, A(151)/0.063924238584648D0/
      DATA X(152)/0.097004699209462D0/, A(152)/0.064466164435950D0/
      DATA X(153)/0.032380170962869D0/, A(153)/0.064737696812683D0/
C-----N=64
      DATA X(154)/0.999305041735772D0/, A(154)/0.001783280721696D0/
      DATA X(155)/0.996340116771955D0/, A(155)/0.004147033260562D0/
      DATA X(156)/0.991013371476744D0/, A(156)/0.006504457968978D0/
      DATA X(157)/0.983336253884625D0/, A(157)/0.008846759826363D0/
      DATA X(158)/0.973326827789910D0/, A(158)/0.011168139460131D0/
      DATA X(159)/0.961008799652053D0/, A(159)/0.013463047896718D0/
      DATA X(160)/0.946411374858402D0/, A(160)/0.015726030476024D0/
      DATA X(161)/0.929569172131939D0/, A(161)/0.017951715775697D0/
      DATA X(162)/0.910522137078502D0/, A(162)/0.020134823153530D0/
      DATA X(163)/0.889315445995114D0/, A(163)/0.022270173808383D0/
      DATA X(164)/0.865999398154092D0/, A(164)/0.024352702568710D0/
      DATA X(165)/0.840629296252580D0/, A(165)/0.026377469715054D0/
      DATA X(166)/0.813265315122797D0/, A(166)/0.028339672614259D0/
      DATA X(167)/0.783972358943341D0/, A(167)/0.030234657072402D0/
      DATA X(168)/0.752819907260531D0/, A(168)/0.032057928354851D0/
      DATA X(169)/0.719881850171610D0/, A(169)/0.033805161837141D0/
      DATA X(170)/0.685236313054233D0/, A(170)/0.035472213256882D0/
      DATA X(171)/0.648965471254657D0/, A(171)/0.037055128540240D0/
      DATA X(172)/0.611155355172393D0/, A(172)/0.038550153178615D0/
      DATA X(173)/0.571895646202634D0/, A(173)/0.039953741132720D0/
      DATA X(174)/0.531279464019894D0/, A(174)/0.041262563242623D0/
      DATA X(175)/0.489403145707052D0/, A(175)/0.042473515123653D0/
      DATA X(176)/0.446366017253464D0/, A(176)/0.043583724529323D0/
      DATA X(177)/0.402270157963991D0/, A(177)/0.044590558163756D0/
      DATA X(178)/0.357220158337668D0/, A(178)/0.045491627927418D0/
      DATA X(179)/0.311322871990210D0/, A(179)/0.046284796581314D0/
      DATA X(180)/0.264687162208767D0/, A(180)/0.046968182816210D0/
      DATA X(181)/0.217423643740007D0/, A(181)/0.047540165714830D0/
      DATA X(182)/0.169644420423992D0/, A(182)/0.047999388596458D0/
      DATA X(183)/0.121462819296120D0/, A(183)/0.048344762234802D0/
      DATA X(184)/0.072993121787799D0/, A(184)/0.048575467441503D0/
      DATA X(185)/0.024350292663424D0/, A(185)/0.048690957009139D0/
C-----N=80
      DATA X(186)/0.999553822651630D0/, A(186)/0.001144950003186D0/
      DATA X(187)/0.997649864398237D0/, A(187)/0.002663533589512D0/
      DATA X(188)/0.994227540965688D0/, A(188)/0.004180313124694D0/
      DATA X(189)/0.989291302499755D0/, A(189)/0.005690922451403D0/
      DATA X(190)/0.982848572738629D0/, A(190)/0.007192904768117D0/
      DATA X(191)/0.974909140585727D0/, A(191)/0.008683945269260D0/
      DATA X(192)/0.965485089043799D0/, A(192)/0.010161766041103D0/
      DATA X(193)/0.954590766343634D0/, A(193)/0.011624114120797D0/
      DATA X(194)/0.942242761309872D0/, A(194)/0.013068761592401D0/
      DATA X(195)/0.928459877172445D0/, A(195)/0.014493508040509D0/
      DATA X(196)/0.913263102571757D0/, A(196)/0.015896183583725D0/
      DATA X(197)/0.896675579438770D0/, A(197)/0.017274652056269D0/
      DATA X(198)/0.878722567678213D0/, A(198)/0.018626814208299D0/
      DATA X(199)/0.859431406663111D0/, A(199)/0.019950610878141D0/
      DATA X(200)/0.838831473580255D0/, A(200)/0.021244026115782D0/
      DATA X(201)/0.816954138681463D0/, A(201)/0.022505090246332D0/
      DATA X(202)/0.793832717504605D0/, A(202)/0.023731882865930D0/
      DATA X(203)/0.769502420135041D0/, A(203)/0.024922535764115D0/
      DATA X(204)/0.744000297583597D0/, A(204)/0.026075235767565D0/
      DATA X(205)/0.717365185362099D0/, A(205)/0.027188227500486D0/
      DATA X(206)/0.689637644342027D0/, A(206)/0.028259816057276D0/
      DATA X(207)/0.660859898986119D0/, A(207)/0.029288369583267D0/
      DATA X(208)/0.631075773046871D0/, A(208)/0.030272321759557D0/
      DATA X(209)/0.600330622829751D0/, A(209)/0.031210174188114D0/
      DATA X(210)/0.568671268122709D0/, A(210)/0.032100498673487D0/
      DATA X(211)/0.536145920897131D0/, A(211)/0.032941939397645D0/
      DATA X(212)/0.502804111888784D0/, A(212)/0.033733214984611D0/
      DATA X(213)/0.468696615170544D0/, A(213)/0.034473120451753D0/
      DATA X(214)/0.433875370831756D0/, A(214)/0.035160529044747D0/
      DATA X(215)/0.398393405881969D0/, A(215)/0.035794393953416D0/
      DATA X(216)/0.362304753499487D0/, A(216)/0.036373749905835D0/
      DATA X(217)/0.325664370747701D0/, A(217)/0.036897714638276D0/
      DATA X(218)/0.288528054884511D0/, A(218)/0.037365490238730D0/
      DATA X(219)/0.250952358392272D0/, A(219)/0.037776364362001D0/
      DATA X(220)/0.212994502857666D0/, A(220)/0.038129711314477D0/
      DATA X(221)/0.174712291832646D0/, A(221)/0.038424993006959D0/
      DATA X(222)/0.136164022809143D0/, A(222)/0.038661759774076D0/
      DATA X(223)/0.097408398441584D0/, A(223)/0.038839651059051D0/
      DATA X(224)/0.058504437152420D0/, A(224)/0.038958395962769D0/
      DATA X(225)/0.019511383256793D0/, A(225)/0.039017813656306D0/
C-----N=96
      DATA X(226)/0.999689503883230D0/, A(226)/0.000796792065552D0/
      DATA X(227)/0.998364375863181D0/, A(227)/0.001853960788946D0/
      DATA X(228)/0.995981842987209D0/, A(228)/0.002910731817934D0/
      DATA X(229)/0.992543900323762D0/, A(229)/0.003964554338444D0/
      DATA X(230)/0.988054126329623D0/, A(230)/0.005014202742927D0/
      DATA X(231)/0.982517263563014D0/, A(231)/0.006058545504235D0/
      DATA X(232)/0.975939174585136D0/, A(232)/0.007096470791153D0/
      DATA X(233)/0.968326828463264D0/, A(233)/0.008126876925698D0/
      DATA X(234)/0.959688291448742D0/, A(234)/0.009148671230783D0/
      DATA X(235)/0.950032717784437D0/, A(235)/0.010160770535008D0/
      DATA X(236)/0.939370339752755D0/, A(236)/0.011162102099838D0/
      DATA X(237)/0.927712456722308D0/, A(237)/0.012151604671088D0/
      DATA X(238)/0.915071423120898D0/, A(238)/0.013128229566961D0/
      DATA X(239)/0.901460635315852D0/, A(239)/0.014090941772314D0/
      DATA X(240)/0.886894517402420D0/, A(240)/0.015038721026994D0/
      DATA X(241)/0.871388505909296D0/, A(241)/0.015970562902562D0/
      DATA X(242)/0.854959033434601D0/, A(242)/0.016885479864245D0/
      DATA X(243)/0.837623511228187D0/, A(243)/0.017782502316045D0/
      DATA X(244)/0.819400310737931D0/, A(244)/0.018660679627411D0/
      DATA X(245)/0.800308744139140D0/, A(245)/0.019519081140145D0/
      DATA X(246)/0.780369043867433D0/, A(246)/0.020356797154333D0/
      DATA X(247)/0.759602341176647D0/, A(247)/0.021172939892191D0/
      DATA X(248)/0.738030643744400D0/, A(248)/0.021966644438744D0/
      DATA X(249)/0.715676812348967D0/, A(249)/0.022737069658329D0/
      DATA X(250)/0.692564536642171D0/, A(250)/0.023483399085926D0/
      DATA X(251)/0.668718310043916D0/, A(251)/0.024204841792364D0/
      DATA X(252)/0.644163403784967D0/, A(252)/0.024900633222483D0/
      DATA X(253)/0.618925840125468D0/, A(253)/0.025570036005349D0/
      DATA X(254)/0.593032364777572D0/, A(254)/0.026212340735672D0/
      DATA X(255)/0.566510418561397D0/, A(255)/0.026826866725591D0/
      DATA X(256)/0.539388108324357D0/, A(256)/0.027412962726029D0/
      DATA X(257)/0.511694177154667D0/, A(257)/0.027970007616848D0/
      DATA X(258)/0.483457973920596D0/, A(258)/0.028497411065085D0/
      DATA X(259)/0.454709422167743D0/, A(259)/0.028994614150555D0/
      DATA X(260)/0.425478988407300D0/, A(260)/0.029461089958167D0/
      DATA X(261)/0.395797649828908D0/, A(261)/0.029896344136328D0/
      DATA X(262)/0.365696861472313D0/, A(262)/0.030299915420827D0/
      DATA X(263)/0.335208522892625D0/, A(263)/0.030671376123669D0/
      DATA X(264)/0.304364944354496D0/, A(264)/0.031010332586313D0/
      DATA X(265)/0.273198812591049D0/, A(265)/0.031316425596861D0/
      DATA X(266)/0.241743156163840D0/, A(266)/0.031589330770727D0/
      DATA X(267)/0.210031310460567D0/, A(267)/0.031828758894411D0/
      DATA X(268)/0.178096882367618D0/, A(268)/0.032034456231992D0/
      DATA X(269)/0.145973714654896D0/, A(269)/0.032206204794030D0/
      DATA X(270)/0.113695850110665D0/, A(270)/0.032343822568575D0/
      DATA X(271)/0.081297495464425D0/, A(271)/0.032447163714064D0/
      DATA X(272)/0.048812985136049D0/, A(272)/0.032516118713868D0/
      DATA X(273)/0.016276744849602D0/, A(273)/0.032550614492363D0/
      DATA IBD/0/

      IF(IBD.NE.0) RETURN
      IBD=1
      DO 10 I=1,273
        B(I) = A(I)
10      Y(I) = X(I)
      DO 20 I=1,96
20      LTAB(I) = KTAB(I)
      RETURN
      END
C=======================================================================

      SUBROUTINE SIB_ALTRA(GA,BGX,BGY,BGZ,PCX,PCY,PCZ,EC,P,PX,PY,PZ,E)

C-----------------------------------------------------------------------
C
C    arbitrary Lorentz transformation
C
C     Input: GA : gamma factor
C            BG? : components of gamma * beta
C            PC?,EC : components of initial 4 vector
C
C     Output: P?,E : components of 4vector in final frame
C             P : 3-norm in final frame, a.k.a momentum
C
C     PHO_ALTRA taken from PHOJET /FR'14
C*********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      DOUBLE PRECISION P,E
      SAVE

c     consistency check: (gamma*beta)**2 = gamma**2 - 1
      BETGAM2 = BGX**2+BGY**2+BGZ**2
      xtst = 1.D0-BETGAM2/GA**2 - 1.D0/GA**2
      IF(abs(xtst).gt.1.D-5) THEN
         WRITE(LUN,*) ' SIB_ALTRA: transf. inconsistent!'
         WRITE(LUN,*) ' SIB_ALTRA: input (GA,GABE):',GA,BGX,BGY,BGZ
      ENDIF
      IF(GA.LT.1.D0) THEN
         WRITE(LUN,*) ' SIB_ALTRA: you are joking right? GAMMA=',GA
         CALL SIB_REJECT('SIB_ALTRA       ')
      ENDIF
      EP=PCX*BGX+PCY*BGY+PCZ*BGZ
      PE=EP/(GA+1.D0)+EC
      PX=PCX+BGX*PE
      PY=PCY+BGY*PE
      PZ=PCZ+BGZ*PE
      P=DSQRT(PX*PX+PY*PY+PZ*PZ)
      E=GA*EC+EP
      END

C=======================================================================

      SUBROUTINE SIB_TRANS(XO,YO,ZO,CDE,SDE,CFE,SFE,X,Y,Z)

C-----------------------------------------------------------------------
C
C  rotation of coordinate frame (1) de rotation around y axis
C                               (2) fe rotation around z axis
C  (inverse rotation to SIB_TRANI)
C
C     Input: ?0 : vector components in initial frame
C            C? : cosine of rotation angle
C            S? : sine of rotation angle
C            DE : angle of rotation around y axis 
C                 (polar angle in spherical coord.)
C            FE : angle of rotation around z axis 
C                 (azimuthal angle in spherical coord.)
C
C     Output: X,Y,Z: components of vector in rotated frame
C
C     PHO_TRANS taken from PHOJET \FR'14
C**********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      SAVE
 
      X= CDE*CFE*XO-SFE*YO+SDE*CFE*ZO
      Y= CDE*SFE*XO+CFE*YO+SDE*SFE*ZO
      Z=-SDE    *XO       +CDE    *ZO
      END

C=======================================================================

      SUBROUTINE SIB_TRANI(XO,YO,ZO,CDE,SDE,CFE,SFE,X,Y,Z)

C-----------------------------------------------------------------------
C
C  rotation of coordinate frame (1) -fe rotation around z axis
C                               (2) -de rotation around y axis
C  (inverse rotation to SIB_TRANS)
C
C     Input: ?0 : vector components in initial frame
C            C? : cosine of rotation angle
C            S? : sine of rotation angle
C            DE : angle of rotation around y axis 
C                 (polar angle in spherical coord.)
C            FE : angle of rotation around z axis 
C                 (azimuthal angle in spherical coord.)
C
C     Output: X,Y,Z: components of vector in rotated frame
C
C     PHO_TRANS taken from PHOJET \FR'14
C**********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      SAVE

      X= CDE*CFE*XO+CDE*SFE*YO-SDE*ZO
      Y=-SFE    *XO+CFE*    YO
      Z= SDE*CFE*XO+SDE*SFE*YO+CDE*ZO
      END

C=======================================================================

      SUBROUTINE SIROBO( NBEG, NEND, THE, PHI, DBEX, DBEY, DBEZ)

C-----------------------------------------------------------------------
C   THIS IS A SLIGHTLY ALTERED VERSION OF "LUROBO" [JETSET63.PYTHIA]   *
C SET TO WORK IN THE SIBYL ENVIROMENT. THE TRANSFORMATION IS PERFORMED *
C ON PARTICLES NUMBER FROM NBEG TO NEND. COMMON BLOCKS CHANGED.        *
C                                      TSS,   Oct '87                  *
C  modification  use directly BETA in double precision in input (PL)   *
C **********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON /S_PLIST/ P(8000,5), LLIST(8000), NP
      DIMENSION ROT(3,3),PV(3),DP(4)
      SAVE

      IF(THE**2+PHI**2 .LE. 1.D-20) GO TO 131
C...ROTATE (TYPICALLY FROM Z AXIS TO DIRECTION THETA,PHI)
       ROT(1,1)=dCOS(THE)*dCOS(PHI)
       ROT(1,2)=-dSIN(PHI)
       ROT(1,3)=dSIN(THE)*dCOS(PHI)
       ROT(2,1)=dCOS(THE)*dSIN(PHI)
       ROT(2,2)=dCOS(PHI)
       ROT(2,3)=dSIN(THE)*dSIN(PHI)
       ROT(3,1)=-dSIN(THE)
       ROT(3,2)=0.D0
       ROT(3,3)=dCOS(THE)
       DO 120 I=NBEG,NEND
       DO 100 J=1,3
 100   PV(J)=P(I,J)
       DO 110 J=1,3
 110   P(I,J)=ROT(J,1)*PV(1)+ROT(J,2)*PV(2)+ROT(J,3)*PV(3)
 120   CONTINUE
 131    IF(DBEX**2+DBEY**2+DBEZ**2 .LE. 1.D-20) GO TO 151
C...LORENTZ BOOST (TYPICALLY FROM REST TO MOMENTUM/ENERGY=BETA)
       DGA=1.D0/DSQRT(1D0-DBEX**2-DBEY**2-DBEZ**2)
       DO 140 I=NBEG, NEND
       DO 130 J=1,4
 130   DP(J)=P(I,J)
       DBEP=DBEX*DP(1)+DBEY*DP(2)+DBEZ*DP(3)
       DGABEP=DGA*(DGA*DBEP/(1.D0+DGA)+DP(4))
       P(I,1)=DP(1)+DGABEP*DBEX
       P(I,2)=DP(2)+DGABEP*DBEY
       P(I,3)=DP(3)+DGABEP*DBEZ
       P(I,4)=DGA*(DP(4)+DBEP)
 140   CONTINUE
 151   RETURN
      END


C=======================================================================

      SUBROUTINE ISWTCH_LMNTS(ia,ib)

C-----------------------------------------------------------------------
      IMPLICIT INTEGER(I-N)
      SAVE

      itmp = ia
      ia = ib
      ib = itmp
      end
C=======================================================================

      SUBROUTINE SWTCH_LMNTS(a,b)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      SAVE

      tmp = a
      a = b
      b = tmp
      end
C=======================================================================

      DOUBLE PRECISION FUNCTION PAWT(A,B,C)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10
      SAVE

C...  c.m.s. Momentum in two particle decays
      PAWT = SQRT((A**2-(B+C)**2+EPS10)*(A**2-(B-C)**2))/(2.D0*A)
      END

C=======================================================================

      SUBROUTINE HSPLI (KF,KP1,KP2)

C-----------------------------------------------------------------------
C...This subroutine splits one hadron of code KF
C.  into 2 partons of code KP1 and KP2
C.  KP1 refers to a color triplet [q or (qq)bar]         
C.  KP2 to a a color anti-triplet [qbar or (qq)]         
C.  allowed inputs:
C.  KF = 6:14 pi0,pi+-,k+-,k0L,k0s, p,n
C.     = -13,-14  pbar,nbar
C.     = 34:39 Sig+, Sig0, Sig-, Xi0, Xi-, Lam0 
C.     = 49: Omega-
C.   \FR'16
C------------------------------------------------
      IMPLICIT NONE

c     external types
      INTEGER KF,KP1,KP2

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)

c     internal types
      INTEGER KPP
      DOUBLE PRECISION R,XBUG,S_RNDM
      SAVE

      IF(IABS(KF).eq.6.or.IABS(KF).eq.27)THEN     ! pi0, rho0
         R = S_RNDM(0)              
         XBUG = 0.D0
         IF(IPAR(19).eq.1) XBUG = 0.5D0
         IF (R.LE.XBUG)  THEN
            KP1 = 1                  
            KP2 = -1
         ELSE
            KP1 = 2
            KP2 = -2
         ENDIF

      ELSEIF(IABS(KF).eq.7)THEN ! pi+
         KP1 = 1                  
         KP2 = -2

      ELSEIF(IABS(KF).eq.8)THEN ! pi-
         KP1 = 2                  
         KP2 = -1

      ELSEIF(IABS(KF).eq.9)THEN                 ! K+
         KP1 = 1                  
         KP2 = -3
      ELSEIF(IABS(KF).eq.10)THEN                ! K-
         KP1 = 3                  
         KP2 = -1
      ELSEIF(IABS(KF).eq.11.or.IABS(KF).eq.12)THEN             ! K0S/K0L
         KP1 = 2
         KP2 = -3
         IF (S_RNDM(1).GT. 0.5D0)  THEN
            KP1 = 3
            KP2 = -2
         ENDIF
      ELSEIF(IABS(KF).eq.21)THEN                ! K0
         KP1 = 2
         KP2 = -3
      ELSEIF(IABS(KF).eq.22)THEN                ! K0bar
         KP1 = 3
         KP2 = -2
      ELSEIF(IABS(KF).eq.33)THEN                 ! phi
         KP1 = 3
         KP2 = -3
      ELSEIF(IABS(KF).eq.13.or.IABS(KF).eq.41)THEN             ! p/pbar,delta+
         R = PAR(53)*S_RNDM(KF)
         IF (R .LT.3.D0)       THEN
            KP1 = 1
            KP2 = 12
         ELSEIF (R .LT. 4.D0)  THEN
            KP1 = 1
            KP2 = 21
         ELSE
            KP1 = 2
            KP2 = 11
         ENDIF
      ELSEIF(IABS(KF).eq.14.or.IABS(KF).eq.42)THEN             ! n/nbar,delta0
         R = 6.D0*S_RNDM(0)                  
         IF (R .LT.3.D0)       THEN
            KP1 = 2
            KP2 = 12
         ELSEIF (R .LT. 4.D0)  THEN
            KP1 = 2
            KP2 = 21
         ELSE
            KP1 = 1
            KP2 = 22
         ENDIF
      ELSEIF(IABS(KF).eq.40)THEN                ! delta++
         KP1 = 1
         KP2 = 11
      ELSEIF(IABS(KF).eq.43)THEN                ! delta-
         KP1 = 2
         KP2 = 22
      ELSEIF(IABS(KF).eq.34)THEN                !Sigma+
         R = 6.D0*S_RNDM(0)                  
         IF (R .LT.3.D0)       THEN
            KP1 = 3
            KP2 = 11
         ELSEIF (R .LT. 4.D0)  THEN
            KP1 = 1
            KP2 = 31
         ELSE
            KP1 = 1
            KP2 = 13
         ENDIF
      ELSEIF(IABS(KF).eq.35.or.IABS(KF).eq.39)THEN             !Sigma0/Lambda0     
c     all configurations equally likely --> Knuth shuffle
c     setup quarks: u,d,s
         CALL SHFFL_QRKS(1,2,3,KP1,KP2)
         
      ELSEIF(IABS(KF).eq.36)THEN                !Sigma-
         R = 6.D0*S_RNDM(0)                  
         IF (R .LT.3.D0)       THEN
            KP1 = 3
            KP2 = 22
         ELSEIF (R .LT. 4.D0)  THEN
            KP1 = 2
            KP2 = 32
         ELSE
            KP1 = 2
            KP2 = 23
         ENDIF
      ELSEIF(IABS(KF).eq.37)THEN                !Xi0
         R = 6.D0*S_RNDM(0)                  
         IF (R .LT.3.D0)       THEN
            KP1 = 1
            KP2 = 33
         ELSEIF (R .LT. 4.D0)  THEN
            KP1 = 3
            KP2 = 13
         ELSE
            KP1 = 1
            KP2 = 33
         ENDIF
      ELSEIF(IABS(KF).eq.38)THEN                !Xi-
         R = 6.D0*S_RNDM(0)                  
         IF (R .LT.3.D0)       THEN
            KP1 = 2
            KP2 = 33
         ELSEIF (R .LT. 4.D0)  THEN
            KP1 = 3
            KP2 = 23
         ELSE
            KP1 = 2
            KP2 = 33
         ENDIF
      ELSEIF(IABS(KF).eq.49)THEN                 ! Omega-
         KP1 = 3
         KP2 = 33

      ELSEIF(IABS(KF).eq.59)THEN                 ! D+
         KP1 = 4
         KP2 = -2

      ELSEIF(IABS(KF).eq.60)THEN                 ! D-
         KP1 = 2
         KP2 = -4

      ELSEIF(IABS(KF).eq.71)THEN                 ! D0
         KP1 = 4
         KP2 = -1

      ELSEIF(IABS(KF).eq.72)THEN                 ! D0bar
         KP1 = 1
         KP2 = -4

      ELSEIF(IABS(KF).eq.73)THEN                 ! eta_c
         KP1 = 4
         KP2 = -4

      ELSEIF(IABS(KF).eq.74)THEN                 ! Ds+
         KP1 = 4
         KP2 = -3

      ELSEIF(IABS(KF).eq.75)THEN                 ! Ds-
         KP1 = 3
         KP2 = -4

      ELSEIF(IABS(KF).eq.76)THEN                 ! Ds*+
         KP1 = 4
         KP2 = -3

      ELSEIF(IABS(KF).eq.77)THEN                 ! Ds*-
         KP1 = 3
         KP2 = -4

      ELSEIF(IABS(KF).eq.78)THEN                 ! D*+
         KP1 = 4
         KP2 = -2

      ELSEIF(IABS(KF).eq.79)THEN                 ! D*-
         KP1 = 2
         KP2 = -4

      ELSEIF(IABS(KF).eq.80)THEN                 ! D*0
         KP1 = 4
         KP2 = -1

      ELSEIF(IABS(KF).eq.81)THEN                 ! D*0bar
         KP1 = 1
         KP2 = -4

      ELSEIF(IABS(KF).eq.83)THEN                 ! J/psi
         KP1 = 4
         KP2 = -4
         
      ELSEIF(IABS(KF).eq.84)THEN                  ! Sigma_c++
         R = 6.D0*S_RNDM(0)                  
         IF (R .LT.3.D0)       THEN
            KP1 = 4
            KP2 = 11
         ELSEIF (R .LT. 4.D0)  THEN
            KP1 = 1
            KP2 = 41
         ELSE
            KP1 = 1
            KP2 = 14
         ENDIF

      ELSEIF(IABS(KF).eq.85.or.IABS(KF).eq.89)THEN               ! Sigma_c+ / Lambda_c+
c     setup quarks: u,d,c
         CALL SHFFL_QRKS(1,2,4,KP1,KP2)

      ELSEIF(IABS(KF).eq.86)THEN                  ! Sigma_c0
         R = 6.D0*S_RNDM(0)                  
         IF (R .LT.3.D0)       THEN
            KP1 = 4
            KP2 = 22
         ELSEIF (R .LT. 4.D0)  THEN
            KP1 = 2
            KP2 = 42
         ELSE
            KP1 = 2
            KP2 = 24
         ENDIF

      ELSEIF(IABS(KF).eq.87)THEN               ! Xi_c+
c     setup quarks: u,s,c
         CALL SHFFL_QRKS(1,3,4,KP1,KP2)

      ELSEIF(IABS(KF).eq.88)THEN                  ! Xi_c0
         CALL SHFFL_QRKS(2,3,4,KP1,KP2)

      ELSEIF(IABS(KF).eq.99)THEN                  ! Omega_c0
         R = 6.D0*S_RNDM(0)                  
         IF (R .LT.3.D0)       THEN
            KP1 = 4
            KP2 = 33
         ELSEIF (R .LT. 4.D0)  THEN
            KP1 = 3
            KP2 = 43
         ELSE
            KP1 = 3
            KP2 = 34
         ENDIF         

      ELSE
C...  Test for good input
         WRITE(LUN,*)
     &        ' HSPLI : Routine entered with illegal particle code ',KF
         CALL SIB_REJECT('HSPLI           ')
      ENDIF

C     if anti-baryon, invert valences
      IF (KF .LT. 0) THEN
         KPP = KP1
         KP1 = -KP2
         KP2 = -KPP
      ENDIF
      RETURN
      END

C=======================================================================      

      SUBROUTINE SHFFL_QRKS(IQF1,IQF2,IQF3,KF1,KF2)

C-----------------------------------------------------------------------
c     routine to shuffle 3 quark flavors
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER IQF1,IQF2,IQF3,KF1,KF2
      INTEGER KPL,JJ,II,IFL
      DOUBLE PRECISION S_RNDM
      DIMENSION KPL(3)
c     quark flavors to shuffle
      KPL(1) = IQF1
      KPL(2) = IQF2
      KPL(3) = IQF3
c     Knuth shuffle..
      DO II=3,2,-1
         JJ=1+INT(II*S_RNDM(II))
         IFL=KPL(jj)
         KPL(jj)=KPL(ii)
         KPL(ii)=IFL
      ENDDO
      KF1=KPL(1)
      KF2=KPL(2)*10+KPL(3)     
      END
            
C.=========================================================================
C.  Library of programs for the generation of nucleus-nucleus interactions
C.  and the study of nucleus-induced cosmic ray showers
C.
C.  September 2001  changes in FPNI, and SIGMA_INI,
C.                  new SIGMA_PP, SIGMA_PPI, SIGMA_KP  (R. Engel)
C.
C.  may  1996       small bug  corrected by Dieter Heck in NUC_CONF 
C.
C.  march 1996      small modification to the superposition code
C.
C.  February 1996   change to FPNI to give an interaction length
C.                   also  at very low energy  
C.
C.  Version 1.01  september 1995 
C.       (small corrections P.L.)
C.       the random number generator is called as S_RNDM(0)
C.  ------------------------------------------------------
C.  Version 1.00  April 1992
C.
C.  Authors:
C.
C.     J. Engel
C.     T.K Gaisser
C.     P.Lipari
C.     T. Stanev
C. 
C.  This set of routines  when used in  the simulation of cosmic ray
C.  showers have only three  "contact points" with the "external world"
C.
C.    (i) SUBROUTINE NUC_NUC_INI
C.        (no  calling arguments)         
C.         to be called once during general initialization
C.    (ii) SUBROUTINE HEAVY (IA, E0)
C.         where IA (integer) is the mass number of the projectile
C.         nucleus  and E0 (TeV) is the energy per nucleon
C.         The output (positions of first interaction for the IA
C.         nucleons of the projectile) is  contained in  the common block:
C.           COMMON /C1STNC/ XX0(60),XX(60),YY(60),AX(60),AY(60)
C.         In detail:
C.             XX0(j)   (g cm-2) =  position of interaction
C.             XX(j) (mm)    x-distance from shower axis
C.             YY(j) (mm)    y-distance from shower axis
C.             AX(j) (radiants)  Theta_x with respect to original direction
C.             AY(j) (radiants)  Theta_y with respect to original direction
C.      
C.    (iii)  FUNCTION FPNI (E,L)
C.           Interaction length in air.
C.           E (TeV) is the energy of the particle, L is the particle
C.           code (NOTE: "Sibyll" codes are used : L =1-18) 
C.           WANRNING : The nucleus-nucleus cross section
C.           tabulated in the program are "matched" to the p-Air
C.           cross section calculated  with this FUNCTION, in other words 
C.           they are both calculated with the same input pp cross section
C==========================================================================

      SUBROUTINE NUC_NUC_INI

C-----------------------------------------------------------------------
C...Initialization for the generation of nucleus-nucleus interactions
C.  INPUT : E0 (TeV) Energy per nucleon of the beam nucleus
C........................................................................
      SAVE

      CALL NUC_GEOM_INI                       ! nucleus profiles
      CALL SIGMA_INI                          ! initialize pp cross sections

      RETURN
      END
C=======================================================================

      SUBROUTINE FRAGM1 (IA,NW, NF, IAF)

C-----------------------------------------------------------------------
C...Nuclear Fragmentation 
C.  total dissolution of nucleus
C.......................................................................
      SAVE

      DIMENSION IAF(60)
      NF = IA-NW
      DO J=1,NF
         IAF(J) = 1
      ENDDO
      RETURN
      END
C->
C=======================================================================

      SUBROUTINE FRAGM2 (IA,NW, NF, IAF)

C-----------------------------------------------------------------------
C...Nuclear Fragmentation 
C.  Spectator in one single fragment 
C.......................................................................
      SAVE

      DIMENSION IAF(60)
      IF (IA-NW .GT. 0)  THEN
         NF = 1
         IAF(1) = IA-NW
      ELSE
         NF = 0
      ENDIF
      RETURN
      END

C-----------------------------------------------------------------------
C...Code of fragmentation  of spectator nucleons
C.  based on Jon Engel  abrasion-ablation algorithms
C=======================================================================

      BLOCK DATA FRAG_DATA

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

C...Data for the fragmentation of  nucleus  projectiles
      COMMON /FRAGMOD/A(10,10,20),AE(10,10,20),ERES(10,10),NFLAGG(10,10)
      SAVE
      DATA (NFLAGG(I, 1),I=1,10)  / 
     +    0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
      DATA (NFLAGG(I, 2),I=1,10)  / 
     +    0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
      DATA (NFLAGG(I, 3),I=1,10)  / 
     +    0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
      DATA (NFLAGG(I, 4),I=1,10)  / 
     +    0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
      DATA (NFLAGG(I, 5),I=1,10)  / 
     +    0,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
      DATA (NFLAGG(I, 6),I=1,10)  / 
     +    0,  0,  0,  0,  0,  0,  0,  1,  1,  1 /
      DATA (NFLAGG(I, 7),I=1,10)  / 
     +    1,  1,  1,  1,  1,  1,  1,  1,  1,  1 /
      DATA (NFLAGG(I, 8),I=1,10)  / 
     +    1,  1,  1,  1,  1,  1,  1,  1,  1,  1 /
      DATA (NFLAGG(I, 9),I=1,10)  / 
     +    1,  1,  1,  1,  1,  1,  1,  1,  1,  1 /
      DATA (NFLAGG(I,10),I=1,10)  / 
     +    1,  1,  1,  1,  1,  1,  1,  1,  1,  1 /
      DATA (A(I, 1, 1),I=1,10)  / 
     +  .438D-01,.172D0  ,.283D0  ,.511D0  ,.715D0  ,.920D0  ,1.19D0  ,
     +  1.37D0  ,1.65D0  ,2.14D0   /
      DATA (A(I, 1, 2),I=1,10)  / 
     +  .147D-01,.249D-01,.439D-01,.592D-01,.776D-01,.886D-01,.108D0  ,
     +  .117D0  ,.126D0  ,.128D0   /
      DATA (A(I, 1, 3),I=1,10)  / 
     +  .216D-02,.627D-02,.834D-02,.108D-01,.144D-01,.152D-01,.196D-01,
     +  .200D-01,.210D-01,.224D-01 /
      DATA (A(I, 1, 4),I=1,10)  / 
     +  .593D-01,.653D-01,.116D0  ,.145D0  ,.184D0  ,.204D0  ,.234D0  ,
     +  .257D0  ,.271D0  ,.248D0   /
      DATA (A(I, 1, 5),I=1,10)  / 
     +  .000D+00,.918D-02,.362D-02,.805D-02,.436D-02,.728D-02,.466D-02,
     +  .707D-02,.932D-02,.130D-01 /
      DATA (A(I, 1, 6),I=1,10)  / 
     +  .000D+00,.180D-02,.247D-02,.208D-02,.224D-02,.214D-02,.226D-02,
     +  .233D-02,.230D-02,.194D-02 /
      DATA (A(I, 1, 7),I=1,10)  / 
     +  .000D+00,.106D-02,.703D-03,.687D-03,.739D-03,.674D-03,.819D-03,
     +  .768D-03,.756D-03,.720D-03 /
      DATA (A(I, 1, 8),I=1,10)  / 
     +  .000D+00,.000D+00,.188D-02,.130D-02,.138D-02,.117D-02,.124D-02,
     +  .119D-02,.111D-02,.829D-03 /
      DATA (A(I, 1, 9),I=1,10)  / 
     +  .000D+00,.000D+00,.302D-03,.258D-03,.249D-03,.208D-03,.248D-03,
     +  .222D-03,.210D-03,.187D-03 /
      DATA (A(I, 1,10),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.235D-03,.222D-03,.172D-03,.181D-03,
     +  .166D-03,.152D-03,.124D-03 /
      DATA (A(I, 1,11),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.238D-03,.179D-03,.145D-03,.156D-03,
     +  .138D-03,.129D-03,.111D-03 /
      DATA (A(I, 1,12),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.368D-03,.400D-03,.255D-03,.262D-03,
     +  .221D-03,.182D-03,.112D-03 /
      DATA (A(I, 1,13),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.753D-04,.712D-04,.527D-04,
     +  .537D-04,.538D-04,.487D-04 /
      DATA (A(I, 1,14),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.103D-03,.589D-04,.578D-04,
     +  .468D-04,.385D-04,.269D-04 /
      DATA (A(I, 1,15),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.444D-04,.372D-04,
     +  .318D-04,.284D-04,.218D-04 /
      DATA (A(I, 1,16),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.487D-04,.473D-04,
     +  .338D-04,.243D-04,.122D-04 /
      DATA (A(I, 1,17),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.121D-04,.117D-04,
     +  .932D-05,.792D-05,.583D-05 /
      DATA (A(I, 1,18),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.147D-04,
     +  .101D-04,.756D-05,.496D-05 /
      DATA (A(I, 1,19),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.755D-05,
     +  .612D-05,.505D-05,.341D-05 /
      DATA (A(I, 1,20),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,
     +  .630D-05,.444D-05,.282D-05 /
      DATA (A(I, 2, 1),I=1,10)  / 
     +  .269D0  ,.510D0  ,.738D0  ,1.12D0  ,1.46D0  ,1.83D0  ,2.22D0  ,
     +  2.57D0  ,3.00D0  ,3.67D0   /
      DATA (A(I, 2, 2),I=1,10)  / 
     +  .121D0  ,.133D0  ,.190D0  ,.234D0  ,.293D0  ,.332D0  ,.395D0  ,
     +  .431D0  ,.468D0  ,.502D0   /
      DATA (A(I, 2, 3),I=1,10)  / 
     +  .227D-01,.374D-01,.474D-01,.578D-01,.722D-01,.794D-01,.960D-01,
     +  .102D0  ,.110D0  ,.120D0   /
      DATA (A(I, 2, 4),I=1,10)  / 
     +  .287D0  ,.196D0  ,.270D0  ,.314D0  ,.373D0  ,.408D0  ,.462D0  ,
     +  .498D0  ,.529D0  ,.523D0   /
      DATA (A(I, 2, 5),I=1,10)  / 
     +  .000D+00,.433D-01,.218D-01,.384D-01,.263D-01,.385D-01,.298D-01,
     +  .405D-01,.504D-01,.671D-01 /
      DATA (A(I, 2, 6),I=1,10)  / 
     +  .000D+00,.151D-01,.177D-01,.159D-01,.173D-01,.173D-01,.187D-01,
     +  .196D-01,.201D-01,.191D-01 /
      DATA (A(I, 2, 7),I=1,10)  / 
     +  .000D+00,.457D-02,.607D-02,.610D-02,.677D-02,.670D-02,.784D-02,
     +  .787D-02,.806D-02,.803D-02 /
      DATA (A(I, 2, 8),I=1,10)  / 
     +  .000D+00,.000D+00,.702D-02,.536D-02,.558D-02,.510D-02,.554D-02,
     +  .546D-02,.538D-02,.489D-02 /
      DATA (A(I, 2, 9),I=1,10)  / 
     +  .000D+00,.000D+00,.190D-02,.199D-02,.205D-02,.191D-02,.221D-02,
     +  .214D-02,.213D-02,.204D-02 /
      DATA (A(I, 2,10),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.226D-02,.219D-02,.195D-02,.208D-02,
     +  .204D-02,.203D-02,.194D-02 /
      DATA (A(I, 2,11),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.213D-02,.195D-02,.175D-02,.191D-02,
     +  .183D-02,.179D-02,.166D-02 /
      DATA (A(I, 2,12),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.588D-03,.186D-02,.137D-02,.141D-02,
     +  .128D-02,.117D-02,.947D-03 /
      DATA (A(I, 2,13),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.554D-03,.562D-03,.454D-03,
     +  .485D-03,.505D-03,.509D-03 /
      DATA (A(I, 2,14),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.490D-03,.533D-03,.531D-03,
     +  .476D-03,.437D-03,.369D-03 /
      DATA (A(I, 2,15),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.427D-03,.382D-03,
     +  .358D-03,.340D-03,.294D-03 /
      DATA (A(I, 2,16),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.239D-03,.298D-03,
     +  .238D-03,.196D-03,.134D-03 /
      DATA (A(I, 2,17),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.299D-04,.893D-04,
     +  .796D-04,.744D-04,.683D-04 /
      DATA (A(I, 2,18),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.127D-03,
     +  .107D-03,.916D-04,.720D-04 /
      DATA (A(I, 2,19),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.397D-04,
     +  .630D-04,.565D-04,.461D-04 /
      DATA (A(I, 2,20),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,
     +  .511D-04,.459D-04,.402D-04 /
      DATA (A(I, 3, 1),I=1,10)  / 
     +  .708D0  ,1.02D0  ,1.41D0  ,1.91D0  ,2.42D0  ,3.00D0  ,3.53D0  ,
     +  4.09D0  ,4.71D0  ,5.57D0   /
      DATA (A(I, 3, 2),I=1,10)  / 
     +  .397D0  ,.410D0  ,.539D0  ,.648D0  ,.795D0  ,.910D0  ,1.06D0  ,
     +  1.17D0  ,1.29D0  ,1.42D0   /
      DATA (A(I, 3, 3),I=1,10)  / 
     +  .845D-01,.122D0  ,.157D0  ,.190D0  ,.232D0  ,.262D0  ,.307D0  ,
     +  .335D0  ,.366D0  ,.402D0   /
      DATA (A(I, 3, 4),I=1,10)  / 
     +  .210D0  ,.379D0  ,.450D0  ,.490D0  ,.574D0  ,.636D0  ,.709D0  ,
     +  .769D0  ,.820D0  ,.849D0   /
      DATA (A(I, 3, 5),I=1,10)  / 
     +  .000D+00,.102D0  ,.675D-01,.104D0  ,.858D-01,.115D0  ,.102D0  ,
     +  .129D0  ,.154D0  ,.194D0   /
      DATA (A(I, 3, 6),I=1,10)  / 
     +  .000D+00,.392D-01,.615D-01,.593D-01,.649D-01,.674D-01,.735D-01,
     +  .779D-01,.817D-01,.828D-01 /
      DATA (A(I, 3, 7),I=1,10)  / 
     +  .000D+00,.539D-02,.222D-01,.238D-01,.269D-01,.280D-01,.320D-01,
     +  .334D-01,.350D-01,.361D-01 /
      DATA (A(I, 3, 8),I=1,10)  / 
     +  .000D+00,.000D+00,.838D-02,.130D-01,.133D-01,.131D-01,.141D-01,
     +  .144D-01,.149D-01,.152D-01 /
      DATA (A(I, 3, 9),I=1,10)  / 
     +  .000D+00,.000D+00,.228D-02,.647D-02,.688D-02,.687D-02,.772D-02,
     +  .786D-02,.811D-02,.824D-02 /
      DATA (A(I, 3,10),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.664D-02,.828D-02,.802D-02,.845D-02,
     +  .869D-02,.902D-02,.930D-02 /
      DATA (A(I, 3,11),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.338D-02,.735D-02,.710D-02,.767D-02,
     +  .767D-02,.776D-02,.756D-02 /
      DATA (A(I, 3,12),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.280D-03,.262D-02,.349D-02,.342D-02,
     +  .322D-02,.312D-02,.291D-02 /
      DATA (A(I, 3,13),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.618D-03,.161D-02,.138D-02,
     +  .148D-02,.155D-02,.166D-02 /
      DATA (A(I, 3,14),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.313D-03,.128D-02,.161D-02,
     +  .150D-02,.144D-02,.134D-02 /
      DATA (A(I, 3,15),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.645D-03,.118D-02,
     +  .115D-02,.111D-02,.103D-02 /
      DATA (A(I, 3,16),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.117D-03,.497D-03,
     +  .581D-03,.501D-03,.401D-03 /
      DATA (A(I, 3,17),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.115D-04,.997D-04,
     +  .202D-03,.203D-03,.206D-03 /
      DATA (A(I, 3,18),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.877D-04,
     +  .242D-03,.263D-03,.226D-03 /
      DATA (A(I, 3,19),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.158D-04,
     +  .881D-04,.152D-03,.136D-03 /
      DATA (A(I, 3,20),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,
     +  .358D-04,.997D-04,.117D-03 /
      DATA (A(I, 4, 1),I=1,10)  / 
     +  .945D0  ,1.29D0  ,1.40D0  ,1.98D0  ,2.73D0  ,3.17D0  ,3.77D0  ,
     +  4.29D0  ,4.78D0  ,5.54D0   /
      DATA (A(I, 4, 2),I=1,10)  / 
     +  .581D0  ,.599D0  ,.645D0  ,.839D0  ,1.10D0  ,1.25D0  ,1.47D0  ,
     +  1.64D0  ,1.78D0  ,1.99D0   /
      DATA (A(I, 4, 3),I=1,10)  / 
     +  .127D0  ,.182D0  ,.202D0  ,.264D0  ,.344D0  ,.387D0  ,.455D0  ,
     +  .504D0  ,.549D0  ,.611D0   /
      DATA (A(I, 4, 4),I=1,10)  / 
     +  .183D0  ,.464D0  ,.351D0  ,.444D0  ,.642D0  ,.659D0  ,.772D0  ,
     +  .830D0  ,.882D0  ,.930D0   /
      DATA (A(I, 4, 5),I=1,10)  / 
     +  .000D+00,.122D0  ,.803D-01,.136D0  ,.134D0  ,.173D0  ,.164D0  ,
     +  .203D0  ,.239D0  ,.300D0   /
      DATA (A(I, 4, 6),I=1,10)  / 
     +  .000D+00,.393D-01,.766D-01,.872D-01,.108D0  ,.111D0  ,.123D0  ,
     +  .132D0  ,.139D0  ,.145D0   /
      DATA (A(I, 4, 7),I=1,10)  / 
     +  .000D+00,.416D-02,.289D-01,.360D-01,.454D-01,.477D-01,.549D-01,
     +  .583D-01,.618D-01,.654D-01 /
      DATA (A(I, 4, 8),I=1,10)  / 
     +  .000D+00,.000D+00,.761D-02,.157D-01,.214D-01,.205D-01,.233D-01,
     +  .241D-01,.255D-01,.271D-01 /
      DATA (A(I, 4, 9),I=1,10)  / 
     +  .000D+00,.000D+00,.238D-02,.803D-02,.123D-01,.123D-01,.140D-01,
     +  .145D-01,.153D-01,.160D-01 /
      DATA (A(I, 4,10),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.695D-02,.150D-01,.154D-01,.166D-01,
     +  .172D-01,.181D-01,.192D-01 /
      DATA (A(I, 4,11),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.355D-02,.104D-01,.143D-01,.156D-01,
     +  .158D-01,.164D-01,.165D-01 /
      DATA (A(I, 4,12),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.112D-03,.276D-02,.568D-02,.736D-02,
     +  .684D-02,.691D-02,.661D-02 /
      DATA (A(I, 4,13),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.740D-03,.222D-02,.339D-02,
     +  .352D-02,.382D-02,.409D-02 /
      DATA (A(I, 4,14),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.369D-03,.160D-02,.322D-02,
     +  .375D-02,.375D-02,.355D-02 /
      DATA (A(I, 4,15),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.750D-03,.190D-02,
     +  .298D-02,.319D-02,.299D-02 /
      DATA (A(I, 4,16),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.260D-03,.673D-03,
     +  .117D-02,.156D-02,.126D-02 /
      DATA (A(I, 4,17),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.283D-05,.131D-03,
     +  .363D-03,.618D-03,.690D-03 /
      DATA (A(I, 4,18),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.205D-03,
     +  .378D-03,.709D-03,.844D-03 /
      DATA (A(I, 4,19),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.654D-05,
     +  .150D-03,.341D-03,.527D-03 /
      DATA (A(I, 4,20),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,
     +  .957D-04,.197D-03,.406D-03 /
      DATA (A(I, 5, 1),I=1,10)  / 
     +  1.16D0  ,1.70D0  ,2.19D0  ,2.79D0  ,3.33D0  ,3.90D0  ,4.49D0  ,
     +  5.07D0  ,5.66D0  ,6.38D0   /
      DATA (A(I, 5, 2),I=1,10)  / 
     +  .779D0  ,.899D0  ,1.09D0  ,1.28D0  ,1.51D0  ,1.71D0  ,1.96D0  ,
     +  2.18D0  ,2.39D0  ,2.62D0   /
      DATA (A(I, 5, 3),I=1,10)  / 
     +  .167D0  ,.263D0  ,.334D0  ,.408D0  ,.482D0  ,.548D0  ,.632D0  ,
     +  .700D0  ,.767D0  ,.840D0   /
      DATA (A(I, 5, 4),I=1,10)  / 
     +  .203D0  ,.565D0  ,.845D0  ,.867D0  ,.906D0  ,.961D0  ,1.08D0  ,
     +  1.13D0  ,1.21D0  ,1.25D0   /
      DATA (A(I, 5, 5),I=1,10)  / 
     +  .000D+00,.129D0  ,.152D0  ,.237D0  ,.208D0  ,.268D0  ,.258D0  ,
     +  .312D0  ,.368D0  ,.450D0   /
      DATA (A(I, 5, 6),I=1,10)  / 
     +  .000D+00,.460D-01,.126D0  ,.174D0  ,.182D0  ,.188D0  ,.208D0  ,
     +  .219D0  ,.233D0  ,.239D0   /
      DATA (A(I, 5, 7),I=1,10)  / 
     +  .000D+00,.289D-02,.380D-01,.611D-01,.788D-01,.845D-01,.974D-01,
     +  .103D0  ,.111D0  ,.117D0   /
      DATA (A(I, 5, 8),I=1,10)  / 
     +  .000D+00,.000D+00,.137D-01,.223D-01,.374D-01,.436D-01,.488D-01,
     +  .488D-01,.524D-01,.547D-01 /
      DATA (A(I, 5, 9),I=1,10)  / 
     +  .000D+00,.000D+00,.162D-02,.114D-01,.198D-01,.263D-01,.315D-01,
     +  .323D-01,.348D-01,.364D-01 /
      DATA (A(I, 5,10),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.149D-01,.240D-01,.320D-01,.428D-01,
     +  .436D-01,.469D-01,.493D-01 /
      DATA (A(I, 5,11),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.562D-02,.194D-01,.290D-01,.408D-01,
     +  .460D-01,.492D-01,.500D-01 /
      DATA (A(I, 5,12),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.476D-04,.106D-01,.134D-01,.191D-01,
     +  .227D-01,.264D-01,.253D-01 /
      DATA (A(I, 5,13),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.281D-02,.679D-02,.879D-02,
     +  .123D-01,.165D-01,.190D-01 /
      DATA (A(I, 5,14),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.542D-04,.847D-02,.125D-01,
     +  .144D-01,.173D-01,.192D-01 /
      DATA (A(I, 5,15),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.652D-02,.982D-02,
     +  .129D-01,.159D-01,.192D-01 /
      DATA (A(I, 5,16),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.109D-03,.688D-02,
     +  .751D-02,.845D-02,.905D-02 /
      DATA (A(I, 5,17),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.823D-06,.237D-02,
     +  .318D-02,.446D-02,.569D-02 /
      DATA (A(I, 5,18),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.604D-03,
     +  .610D-02,.673D-02,.827D-02 /
      DATA (A(I, 5,19),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.716D-06,
     +  .412D-02,.519D-02,.617D-02 /
      DATA (A(I, 5,20),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,
     +  .710D-03,.543D-02,.674D-02 /
      DATA (A(I, 6, 1),I=1,10)  / 
     +  1.36D0  ,2.08D0  ,2.67D0  ,3.30D0  ,3.94D0  ,4.62D0  ,5.18D0  ,
     +  3.60D0  ,3.64D0  ,3.95D0   /
      DATA (A(I, 6, 2),I=1,10)  / 
     +  1.07D0  ,1.33D0  ,1.58D0  ,1.82D0  ,2.10D0  ,2.44D0  ,2.74D0  ,
     +  1.78D0  ,1.73D0  ,1.80D0   /
      DATA (A(I, 6, 3),I=1,10)  / 
     +  .158D0  ,.276D0  ,.402D0  ,.506D0  ,.609D0  ,.700D0  ,.802D0  ,
     +  .638D0  ,.629D0  ,.658D0   /
      DATA (A(I, 6, 4),I=1,10)  / 
     +  .308D0  ,.739D0  ,1.02D0  ,1.12D0  ,1.26D0  ,1.35D0  ,1.57D0  ,
     +  1.94D0  ,1.71D0  ,1.55D0   /
      DATA (A(I, 6, 5),I=1,10)  / 
     +  .000D+00,.217D0  ,.183D0  ,.324D0  ,.276D0  ,.395D0  ,.393D0  ,
     +  .558D0  ,.602D0  ,.681D0   /
      DATA (A(I, 6, 6),I=1,10)  / 
     +  .000D+00,.658D-01,.251D0  ,.267D0  ,.299D0  ,.326D0  ,.386D0  ,
     +  .452D0  ,.475D0  ,.409D0   /
      DATA (A(I, 6, 7),I=1,10)  / 
     +  .000D+00,.198D-02,.774D-01,.136D0  ,.149D0  ,.164D0  ,.187D0  ,
     +  .210D0  ,.238D0  ,.256D0   /
      DATA (A(I, 6, 8),I=1,10)  / 
     +  .000D+00,.000D+00,.290D-01,.122D0  ,.139D0  ,.128D0  ,.129D0  ,
     +  .137D0  ,.147D0  ,.167D0   /
      DATA (A(I, 6, 9),I=1,10)  / 
     +  .000D+00,.000D+00,.699D-03,.617D-01,.750D-01,.801D-01,.905D-01,
     +  .974D-01,.105D0  ,.122D0   /
      DATA (A(I, 6,10),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.310D-01,.112D0  ,.127D0  ,.140D0  ,
     +  .143D0  ,.155D0  ,.176D0   /
      DATA (A(I, 6,11),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.277D-02,.889D-01,.143D0  ,.150D0  ,
     +  .175D0  ,.184D0  ,.208D0   /
      DATA (A(I, 6,12),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.202D-04,.343D-01,.959D-01,.109D0  ,
     +  .115D0  ,.112D0  ,.116D0   /
      DATA (A(I, 6,13),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.186D-02,.435D-01,.512D-01,
     +  .744D-01,.856D-01,.103D0   /
      DATA (A(I, 6,14),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.144D-04,.427D-01,.786D-01,
     +  .911D-01,.993D-01,.108D0   /
      DATA (A(I, 6,15),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.466D-02,.518D-01,
     +  .848D-01,.109D0  ,.119D0   /
      DATA (A(I, 6,16),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.655D-05,.330D-01,
     +  .586D-01,.617D-01,.594D-01 /
      DATA (A(I, 6,17),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.228D-06,.328D-02,
     +  .190D-01,.301D-01,.454D-01 /
      DATA (A(I, 6,18),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.218D-04,
     +  .272D-01,.501D-01,.707D-01 /
      DATA (A(I, 6,19),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.146D-06,
     +  .441D-02,.378D-01,.556D-01 /
      DATA (A(I, 6,20),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,
     +  .160D-03,.204D-01,.679D-01 /
      DATA (A(I, 7, 1),I=1,10)  / 
     +  .522D0  ,.862D0  ,1.14D0  ,1.40D0  ,1.70D0  ,1.94D0  ,2.26D0  ,
     +  2.48D0  ,2.72D0  ,3.95D0   /
      DATA (A(I, 7, 2),I=1,10)  / 
     +  .314D0  ,.450D0  ,.588D0  ,.692D0  ,.834D0  ,.936D0  ,1.09D0  ,
     +  1.18D0  ,1.28D0  ,1.80D0   /
      DATA (A(I, 7, 3),I=1,10)  / 
     +  .814D-01,.147D0  ,.189D0  ,.226D0  ,.272D0  ,.302D0  ,.351D0  ,
     +  .378D0  ,.406D0  ,.658D0   /
      DATA (A(I, 7, 4),I=1,10)  / 
     +  .252D0  ,.864D0  ,1.01D0  ,.851D0  ,.837D0  ,.774D0  ,.763D0  ,
     +  .757D0  ,.748D0  ,1.55D0   /
      DATA (A(I, 7, 5),I=1,10)  / 
     +  .000D+00,.225D0  ,.180D0  ,.276D0  ,.193D0  ,.240D0  ,.190D0  ,
     +  .228D0  ,.259D0  ,.681D0   /
      DATA (A(I, 7, 6),I=1,10)  / 
     +  .000D+00,.485D-01,.272D0  ,.273D0  ,.253D0  ,.216D0  ,.206D0  ,
     +  .197D0  ,.191D0  ,.409D0   /
      DATA (A(I, 7, 7),I=1,10)  / 
     +  .000D+00,.137D-02,.752D-01,.137D0  ,.152D0  ,.134D0  ,.125D0  ,
     +  .119D0  ,.116D0  ,.256D0   /
      DATA (A(I, 7, 8),I=1,10)  / 
     +  .000D+00,.000D+00,.220D-01,.155D0  ,.175D0  ,.155D0  ,.116D0  ,
     +  .977D-01,.858D-01,.167D0   /
      DATA (A(I, 7, 9),I=1,10)  / 
     +  .000D+00,.000D+00,.326D-03,.695D-01,.881D-01,.106D0  ,.897D-01,
     +  .782D-01,.706D-01,.122D0   /
      DATA (A(I, 7,10),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.261D-01,.124D0  ,.131D0  ,.156D0  ,
     +  .141D0  ,.121D0  ,.176D0   /
      DATA (A(I, 7,11),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.785D-03,.864D-01,.130D0  ,.170D0  ,
     +  .182D0  ,.172D0  ,.208     /
      DATA (A(I, 7,12),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.896D-05,.225D-01,.105D0  ,.126D0  ,
     +  .126D0  ,.135D0  ,.116D0   /
      DATA (A(I, 7,13),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.542D-03,.427D-01,.553D-01,
     +  .744D-01,.980D-01,.103D0   /
      DATA (A(I, 7,14),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.515D-05,.377D-01,.831D-01,
     +  .985D-01,.104D0  ,.108D0   /
      DATA (A(I, 7,15),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.285D-02,.495D-01,
     +  .871D-01,.106D0  ,.119D0   /
      DATA (A(I, 7,16),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.110D-05,.284D-01,
     +  .588D-01,.657D-01,.594D-01 /
      DATA (A(I, 7,17),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.722D-07,.176D-02,
     +  .170D-01,.305D-01,.454D-01 /
      DATA (A(I, 7,18),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.148D-05,
     +  .213D-01,.492D-01,.707D-01 /
      DATA (A(I, 7,19),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.323D-07,
     +  .722D-02,.359D-01,.556D-01 /
      DATA (A(I, 7,20),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,
     +  .461D-05,.155D-01,.679D-01 /
      DATA (A(I, 8, 1),I=1,10)  / 
     +  .630D0  ,.974D0  ,1.29D0  ,1.58D0  ,1.89D0  ,2.16D0  ,2.49D0  ,
     +  2.75D0  ,3.02D0  ,3.95D0   /
      DATA (A(I, 8, 2),I=1,10)  / 
     +  .328D0  ,.459D0  ,.613D0  ,.735D0  ,.879D0  ,.994D0  ,1.15D0  ,
     +  1.27D0  ,1.38D0  ,1.80D0   /
      DATA (A(I, 8, 3),I=1,10)  / 
     +  .748D-01,.121D0  ,.164D0  ,.197D0  ,.235D0  ,.265D0  ,.310D0  ,
     +  .339D0  ,.370D0  ,.658D0   /
      DATA (A(I, 8, 4),I=1,10)  / 
     +  .194D0  ,.211D0  ,.337D0  ,.344D0  ,.339D0  ,.351D0  ,.390    ,
     +  .419D0  ,.442D0  ,1.55D0   /
      DATA (A(I, 8, 5),I=1,10)  / 
     +  .000D+00,.869D-01,.725D-01,.113D0  ,.810D-01,.106D0  ,.951D-01,
     +  .120D0  ,.143D0  ,.681D0   /
      DATA (A(I, 8, 6),I=1,10)  / 
     +  .000D+00,.288D-01,.102D0  ,.922D-01,.857D-01,.845D-01,.932D-01,
     +  .983D-01,.102D0  ,.409D0   /
      DATA (A(I, 8, 7),I=1,10)  / 
     +  .000D+00,.668D-03,.533D-01,.575D-01,.493D-01,.482D-01,.539D-01,
     +  .558D-01,.582D-01,.256D0   /
      DATA (A(I, 8, 8),I=1,10)  / 
     +  .000D+00,.000D+00,.205D-01,.808D-01,.510D-01,.409D-01,.406D-01,
     +  .394D-01,.389D-01,.167D0   /
      DATA (A(I, 8, 9),I=1,10)  / 
     +  .000D+00,.000D+00,.999D-04,.647D-01,.385D-01,.325D-01,.325D-01,
     +  .316D-01,.314D-01,.122D0   /
      DATA (A(I, 8,10),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.169D-01,.834D-01,.611D-01,.565D-01,
     +  .533D-01,.519D-01,.176D0   /
      DATA (A(I, 8,11),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.107D-03,.769D-01,.922D-01,.805D-01,
     +  .745D-01,.711D-01,.208D0   /
      DATA (A(I, 8,12),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.180D-05,.143D-01,.983D-01,.775D-01,
     +  .627D-01,.541D-01,.116D0   /
      DATA (A(I, 8,13),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.157D-04,.346D-01,.507D-01,
     +  .479D-01,.455D-01,.103D0   /
      DATA (A(I, 8,14),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.752D-06,.248D-01,.721D-01,
     +  .728D-01,.611D-01,.108D0   /
      DATA (A(I, 8,15),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.686D-04,.356D-01,
     +  .731D-01,.791D-01,.119D0   /
      DATA (A(I, 8,16),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.838D-07,.151D-01,
     +  .470D-01,.567D-01,.594D-01 /
      DATA (A(I, 8,17),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.759D-08,.400D-04,
     +  .193D-01,.313D-01,.454D-01 /
      DATA (A(I, 8,18),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.385D-07,
     +  .921D-02,.353D-01,.707D-01 /
      DATA (A(I, 8,19),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.219D-08,
     +  .348D-03,.226D-01,.556D-01 /
      DATA (A(I, 8,20),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,
     +  .212D-07,.149D-01,.679D-01 /
      DATA (A(I, 9, 1),I=1,10)  / 
     +  .736D0  ,1.13D0  ,1.49D0  ,1.82D0  ,2.20D0  ,2.49D0  ,2.86D0  ,
     +  3.17D0  ,3.49D0  ,3.95D0   /
      DATA (A(I, 9, 2),I=1,10)  / 
     +  .339D0  ,.492D0  ,.658D0  ,.789D0  ,.958D0  ,1.08D0  ,1.25D0  ,
     +  1.37D0  ,1.50D0  ,1.80D0   /
      DATA (A(I, 9, 3),I=1,10)  / 
     +  .680D-01,.110D0  ,.150D0  ,.180D0  ,.222D0  ,.247D0  ,.289    ,
     +  .318D0  ,.349D0  ,.658D0   /
      DATA (A(I, 9, 4),I=1,10)  / 
     +  .110D0  ,.104D0  ,.157D0  ,.156D0  ,.210D0  ,.205D0  ,.246D0  ,
     +  .274D0  ,.300D0  ,1.55D0   /
      DATA (A(I, 9, 5),I=1,10)  / 
     +  .000D+00,.379D-01,.347D-01,.477D-01,.486D-01,.576D-01,.569D-01,
     +  .732D-01,.893D-01,.681D0   /
      DATA (A(I, 9, 6),I=1,10)  / 
     +  .000D+00,.223D-01,.354D-01,.312D-01,.436D-01,.400D-01,.489D-01,
     +  .548D-01,.600D-01,.409D0   /
      DATA (A(I, 9, 7),I=1,10)  / 
     +  .000D+00,.338D-03,.149D-01,.142D-01,.215D-01,.188D-01,.248D-01,
     +  .278D-01,.307D-01,.256D0   /
      DATA (A(I, 9, 8),I=1,10)  / 
     +  .000D+00,.000D+00,.553D-02,.862D-02,.150D-01,.106D-01,.145D-01,
     +  .165D-01,.181D-01,.167D0   /
      DATA (A(I, 9, 9),I=1,10)  / 
     +  .000D+00,.000D+00,.375D-04,.641D-02,.111D-01,.792D-02,.112D-01,
     +  .127D-01,.140D-01,.122D0   /
      DATA (A(I, 9,10),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.112D-01,.200D-01,.127D-01,.176D-01,
     +  .200D-01,.220D-01,.176D0   /
      DATA (A(I, 9,11),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.244D-04,.261D-01,.162D-01,.232D-01,
     +  .263D-01,.287D-01,.208D0   /
      DATA (A(I, 9,12),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.455D-06,.635D-02,.121D-01,.186D-01,
     +  .201D-01,.207D-01,.116D0   /
      DATA (A(I, 9,13),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.146D-05,.922D-02,.116D-01,
     +  .145D-01,.165D-01,.103D0   /
      DATA (A(I, 9,14),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.135D-06,.128D-01,.202D-01,
     +  .215D-01,.220D-01,.108D0   /
      DATA (A(I, 9,15),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.237D-05,.229D-01,
     +  .259D-01,.271D-01,.119D0   /
      DATA (A(I, 9,16),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.100D-07,.534D-02,
     +  .210D-01,.193D-01,.594D-01 /
      DATA (A(I, 9,17),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.915D-09,.847D-06,
     +  .119D-01,.125D-01,.454D-01 /
      DATA (A(I, 9,18),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.298D-08,
     +  .101D-01,.242D-01,.707D-01 /
      DATA (A(I, 9,19),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.196D-09,
     +  .243D-05,.234D-01,.556D-01 /
      DATA (A(I, 9,20),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,
     +  .575D-09,.364D-02,.679D-01 /
      DATA (A(I,10, 1),I=1,10)  / 
     +  .959D0  ,1.46D0  ,1.92D0  ,2.34D0  ,2.80D0  ,3.24D0  ,3.64D0  ,
     +  4.05D0  ,4.48D0  ,3.95D0   /
      DATA (A(I,10, 2),I=1,10)  / 
     +  .343D0  ,.516D0  ,.692D0  ,.836D0  ,1.01D0  ,1.16D0  ,1.31D0  ,
     +  1.46D0  ,1.61D0  ,1.80D0   /
      DATA (A(I,10, 3),I=1,10)  / 
     +  .512D-01,.837D-01,.115D0  ,.138D0  ,.169D0  ,.195D0  ,.220D0  ,
     +  .245D0  ,.270D0  ,.658D0   /
      DATA (A(I,10, 4),I=1,10)  / 
     +  .274D-01,.361D-01,.510D-01,.562D-01,.703D-01,.828D-01,.877D-01,
     +  .996D-01,.111D0  ,1.55D0   /
      DATA (A(I,10, 5),I=1,10)  / 
     +  .000D+00,.850D-02,.875D-02,.118D-01,.124D-01,.170D-01,.154D-01,
     +  .194D-01,.237D-01,.681D0   /
      DATA (A(I,10, 6),I=1,10)  / 
     +  .000D+00,.345D-02,.519D-02,.533D-02,.691D-02,.842D-02,.844D-02,
     +  .987D-02,.113D-01,.409D0   /
      DATA (A(I,10, 7),I=1,10)  / 
     +  .000D+00,.722D-04,.130D-02,.135D-02,.189D-02,.240D-02,.235D-02,
     +  .281D-02,.331D-02,.256D0   /
      DATA (A(I,10, 8),I=1,10)  / 
     +  .000D+00,.000D+00,.283D-03,.272D-03,.394D-03,.557D-03,.480D-03,
     +  .616D-03,.775D-03,.167D0   /
      DATA (A(I,10, 9),I=1,10)  / 
     +  .000D+00,.000D+00,.457D-05,.122D-03,.192D-03,.275D-03,.225D-03,
     +  .292D-03,.373D-03,.122D0   /
      DATA (A(I,10,10),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.119D-03,.185D-03,.278D-03,.201D-03,
     +  .274D-03,.364D-03,.176D0   /
      DATA (A(I,10,11),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.140D-05,.129D-03,.200D-03,.137D-03,
     +  .188D-03,.252D-03,.208D0   /
      DATA (A(I,10,12),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.207D-07,.307D-04,.518D-04,.278D-04,
     +  .421D-04,.608D-04,.116D0   /
      DATA (A(I,10,13),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.306D-07,.252D-04,.111D-04,
     +  .188D-04,.295D-04,.103D0   /
      DATA (A(I,10,14),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.321D-08,.220D-04,.104D-04,
     +  .162D-04,.243D-04,.108D0   /
      DATA (A(I,10,15),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.770D-08,.632D-05,
     +  .105D-04,.162D-04,.119D0   /
      DATA (A(I,10,16),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.117D-09,.199D-05,
     +  .321D-05,.492D-05,.594D-01 /
      DATA (A(I,10,17),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.888E-11,.323D-09,
     +  .106D-05,.192D-05,.454D-01 /
      DATA (A(I,10,18),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.174E-10,
     +  .131D-05,.218D-05,.707D-01 /
      DATA (A(I,10,19),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.994E-12,
     +  .233D-09,.104D-05,.556D-01 /
      DATA (A(I,10,20),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,
     +  .144E-11,.724D-06,.679D-01 /
      DATA (AE(I, 1, 1),I=1,10)  / 
     +  7.27D0  ,6.29D0  ,7.76D0  ,6.70D0  ,8.17D0  ,7.34D0  ,8.70D0  ,
     +  8.02D0  ,7.37D0  ,6.18D0   /
      DATA (AE(I, 1, 2),I=1,10)  / 
     +  7.41D0  ,7.52D0  ,8.14D0  ,8.20D0  ,8.96D0  ,9.05D0  ,9.96D0  ,
     +  10.0D0  ,10.1D0  ,9.86D0   /
      DATA (AE(I, 1, 3),I=1,10)  / 
     +  7.72D0  ,7.69D0  ,9.17D0  ,8.99D0  ,10.6D0  ,10.5D0  ,12.1D0  ,
     +  12.1D0  ,12.0D0  ,11.5D0   /
      DATA (AE(I, 1, 4),I=1,10)  / 
     +  7.90D0  ,8.48D0  ,9.50D0  ,9.94D0  ,10.8D0  ,11.4D0  ,12.2D0  ,
     +  12.8D0  ,13.3D0  ,13.8D0   /
      DATA (AE(I, 1, 5),I=1,10)  / 
     +  .000D+00,8.52D0  ,9.59D0  ,10.1D0  ,11.1D0  ,11.8D0  ,12.7D0  ,
     +  13.3D0  ,13.8D0  ,14.4D0   /
      DATA (AE(I, 1, 6),I=1,10)  / 
     +  .000D+00,9.00D0  ,10.7D0  ,11.7D0  ,13.2D0  ,14.2D0  ,15.6D0  ,
     +  16.5D0  ,17.3D0  ,18.0D0   /
      DATA (AE(I, 1, 7),I=1,10)  / 
     +  .000D+00,9.01D0  ,11.1D0  ,11.9D0  ,14.3D0  ,15.0D0  ,17.4D0  ,
     +  18.0D0  ,18.6D0  ,18.8D0   /
      DATA (AE(I, 1, 8),I=1,10)  / 
     +  .000D+00,.000D+00,11.2D0  ,12.4D0  ,14.5D0  ,15.7D0  ,17.6D0  ,
     +  18.8D0  ,19.9D0  ,20.9D0   /
      DATA (AE(I, 1, 9),I=1,10)  / 
     +  .000D+00,.000D+00,11.4D0  ,12.7D0  ,15.5D0  ,16.6D0  ,19.3D0  ,
     +  20.2D0  ,21.1D0  ,21.7D0   /
      DATA (AE(I, 1,10),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,13.2D0  ,15.8D0  ,17.3D0  ,19.9D0  ,
     +  21.2D0  ,22.4D0  ,23.2D0   /
      DATA (AE(I, 1,11),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,13.2D0  ,16.3D0  ,17.8D0  ,20.8D0  ,
     +  22.1D0  ,23.3D0  ,24.2D0   /
      DATA (AE(I, 1,12),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,13.4D0  ,16.2D0  ,18.2D0  ,21.0D0  ,
     +  22.8D0  ,24.4D0  ,25.9D0   /
      DATA (AE(I, 1,13),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,16.5D0  ,18.4D0  ,21.6D0  ,
     +  23.2D0  ,24.8D0  ,26.2D0   /
      DATA (AE(I, 1,14),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,16.7D0  ,19.0D0  ,22.3D0  ,
     +  24.3D0  ,26.1D0  ,27.4D0   /
      DATA (AE(I, 1,15),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,19.1D0  ,22.8D0  ,
     +  24.7D0  ,26.6D0  ,28.2D0   /
      DATA (AE(I, 1,16),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,19.2D0  ,23.0D0  ,
     +  25.3D0  ,27.5D0  ,29.5D0   /
      DATA (AE(I, 1,17),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,19.6D0  ,23.3D0  ,
     +  25.6D0  ,27.8D0  ,29.6D0   /
      DATA (AE(I, 1,18),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,23.6D0  ,
     +  26.2D0  ,28.5D0  ,30.4D0   /
      DATA (AE(I, 1,19),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,23.7D0  ,
     +  26.3D0  ,28.8D0  ,31.0D0   /
      DATA (AE(I, 1,20),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,
     +  26.5D0  ,29.2D0  ,31.5D0   /
      DATA (AE(I, 2, 1),I=1,10)  / 
     +  8.74D0  ,8.16D0  ,9.25D0  ,8.45D0  ,9.46D0  ,8.90D0  ,9.83D0  ,
     +  9.38D0  ,8.96D0  ,8.15D0   /
      DATA (AE(I, 2, 2),I=1,10)  / 
     +  8.96D0  ,9.30D0  ,9.95D0  ,10.0D0  ,10.8D0  ,10.9D0  ,11.7D0  ,
     +  11.8D0  ,11.9D0  ,11.8D0   /
      DATA (AE(I, 2, 3),I=1,10)  / 
     +  9.44D0  ,9.66D0  ,11.0D0  ,11.0D0  ,12.3D0  ,12.5D0  ,13.7D0  ,
     +  13.9D0  ,14.0D0  ,13.8D0   /
      DATA (AE(I, 2, 4),I=1,10)  / 
     +  8.86D0  ,9.81D0  ,10.8D0  ,11.2D0  ,12.0D0  ,12.6D0  ,13.4D0  ,
     +  14.0D0  ,14.5D0  ,15.1D0   /
      DATA (AE(I, 2, 5),I=1,10)  / 
     +  .000D+00,10.2D0  ,11.4D0  ,12.0D0  ,12.9D0  ,13.6D0  ,14.5D0  ,
     +  15.1D0  ,15.7D0  ,16.3D0   /
      DATA (AE(I, 2, 6),I=1,10)  / 
     +  .000D+00,10.7D0  ,12.5D0  ,13.5D0  ,15.1D0  ,16.0D0  ,17.5D0  ,
     +  18.3D0  ,19.2D0  ,19.9D0   /
      DATA (AE(I, 2, 7),I=1,10)  / 
     +  .000D+00,11.5D0  ,12.9D0  ,13.9D0  ,16.1D0  ,17.0D0  ,19.1D0  ,
     +  19.8D0  ,20.6D0  ,21.0D0   /
      DATA (AE(I, 2, 8),I=1,10)  / 
     +  .000D+00,.000D+00,12.4D0  ,13.8D0  ,15.9D0  ,17.2D0  ,19.1D0  ,
     +  20.3D0  ,21.4D0  ,22.3D0   /
      DATA (AE(I, 2, 9),I=1,10)  / 
     +  .000D+00,.000D+00,13.4D0  ,14.5D0  ,17.1D0  ,18.3D0  ,20.9D0  ,
     +  21.9D0  ,23.0D0  ,23.7D0   /
      DATA (AE(I, 2,10),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,14.9D0  ,17.5D0  ,19.1D0  ,21.6D0  ,
     +  22.9D0  ,24.1D0  ,25.0D0   /
      DATA (AE(I, 2,11),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,15.0D0  ,18.0D0  ,19.6D0  ,22.4D0  ,
     +  23.8D0  ,25.2D0  ,26.2D0   /
      DATA (AE(I, 2,12),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,16.2D0  ,17.3D0  ,19.4D0  ,22.2D0  ,
     +  24.0D0  ,25.7D0  ,27.2D0   /
      DATA (AE(I, 2,13),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,17.8D0  ,19.8D0  ,22.9D0  ,
     +  24.6D0  ,26.2D0  ,27.7D0   /
      DATA (AE(I, 2,14),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,19.1D0  ,20.4D0  ,23.7D0  ,
     +  25.7D0  ,27.6D0  ,29.1D0   /
      DATA (AE(I, 2,15),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,20.5D0  ,24.1D0  ,
     +  26.1D0  ,28.1D0  ,29.9D0   /
      DATA (AE(I, 2,16),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,20.9D0  ,23.9D0  ,
     +  26.4D0  ,28.7D0  ,30.7D0   /
      DATA (AE(I, 2,17),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,22.4D0  ,24.2D0  ,
     +  26.7D0  ,29.0D0  ,30.9D0   /
      DATA (AE(I, 2,18),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,24.8D0  ,
     +  27.3D0  ,29.7D0  ,31.8D0   /
      DATA (AE(I, 2,19),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,26.1D0  ,
     +  27.3D0  ,29.9D0  ,32.3D0   /
      DATA (AE(I, 2,20),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,
     +  27.4D0  ,30.1D0  ,32.6D0   /
      DATA (AE(I, 3, 1),I=1,10)  / 
     +  11.0D0  ,11.0D0  ,11.7D0  ,11.3D0  ,11.9D0  ,11.4D0  ,12.1D0  ,
     +  11.7D0  ,11.5D0  ,11.0D0   /
      DATA (AE(I, 3, 2),I=1,10)  / 
     +  11.2D0  ,12.0D0  ,12.7D0  ,12.9D0  ,13.6D0  ,13.7D0  ,14.4D0  ,
     +  14.6D0  ,14.7D0  ,14.6D0   /
      DATA (AE(I, 3, 3),I=1,10)  / 
     +  12.1D0  ,12.6D0  ,13.7D0  ,13.9D0  ,15.0D0  ,15.2D0  ,16.3D0  ,
     +  16.5D0  ,16.7D0  ,16.7D0   /
      DATA (AE(I, 3, 4),I=1,10)  / 
     +  12.6D0  ,11.3D0  ,12.4D0  ,13.0D0  ,13.8D0  ,14.2D0  ,15.0D0  ,
     +  15.6D0  ,16.1D0  ,16.6D0   /
      DATA (AE(I, 3, 5),I=1,10)  / 
     +  .000D+00,12.6D0  ,13.7D0  ,14.4D0  ,15.3D0  ,16.0D0  ,16.8D0  ,
     +  17.5D0  ,18.1D0  ,18.6D0   /
      DATA (AE(I, 3, 6),I=1,10)  / 
     +  .000D+00,14.0D0  ,14.6D0  ,15.8D0  ,17.4D0  ,18.4D0  ,19.8D0  ,
     +  20.6D0  ,21.5D0  ,22.2D0   /
      DATA (AE(I, 3, 7),I=1,10)  / 
     +  .000D+00,16.0D0  ,15.2D0  ,16.3D0  ,18.3D0  ,19.3D0  ,21.1D0  ,
     +  22.0D0  ,22.8D0  ,23.5D0   /
      DATA (AE(I, 3, 8),I=1,10)  / 
     +  .000D+00,.000D+00,15.6D0  ,15.1D0  ,17.2D0  ,18.6D0  ,20.6D0  ,
     +  21.8D0  ,22.9D0  ,23.8D0   /
      DATA (AE(I, 3, 9),I=1,10)  / 
     +  .000D+00,.000D+00,17.8D0  ,16.3D0  ,18.8D0  ,20.1D0  ,22.5D0  ,
     +  23.6D0  ,24.7D0  ,25.6D0   /
      DATA (AE(I, 3,10),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,17.5D0  ,19.0D0  ,20.7D0  ,23.1D0  ,
     +  24.5D0  ,25.8D0  ,26.8D0   /
      DATA (AE(I, 3,11),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,19.2D0  ,19.4D0  ,21.1D0  ,23.8D0  ,
     +  25.4D0  ,26.8D0  ,28.0D0   /
      DATA (AE(I, 3,12),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,20.7D0  ,19.6D0  ,19.7D0  ,22.4D0  ,
     +  24.4D0  ,26.2D0  ,27.9D0   /
      DATA (AE(I, 3,13),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,21.6D0  ,20.4D0  ,23.2D0  ,
     +  25.1D0  ,26.9D0  ,28.5D0   /
      DATA (AE(I, 3,14),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,23.5D0  ,22.0D0  ,23.8D0  ,
     +  26.1D0  ,28.1D0  ,29.9D0   /
      DATA (AE(I, 3,15),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,23.7D0  ,24.2D0  ,
     +  26.3D0  ,28.5D0  ,30.4D0   /
      DATA (AE(I, 3,16),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,25.4D0  ,24.8D0  ,
     +  25.6D0  ,28.1D0  ,30.5D0   /
      DATA (AE(I, 3,17),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,26.9D0  ,26.8D0  ,
     +  26.1D0  ,28.4D0  ,30.8D0   /
      DATA (AE(I, 3,18),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,28.8D0  ,
     +  27.6D0  ,29.0D0  ,31.5D0   /
      DATA (AE(I, 3,19),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,30.5D0  ,
     +  29.2D0  ,28.9D0  ,31.5D0   /
      DATA (AE(I, 3,20),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,
     +  31.0D0  ,30.0D0  ,31.7D0   /
      DATA (AE(I, 4, 1),I=1,10)  / 
     +  13.0D0  ,13.2D0  ,14.8D0  ,14.2D0  ,14.2D0  ,14.1D0  ,14.5D0  ,
     +  14.4D0  ,14.3D0  ,14.0D0   /
      DATA (AE(I, 4, 2),I=1,10)  / 
     +  13.5D0  ,14.5D0  ,16.1D0  ,15.9D0  ,16.0D0  ,16.3D0  ,16.8D0  ,
     +  17.0D0  ,17.1D0  ,17.2D0   /
      DATA (AE(I, 4, 3),I=1,10)  / 
     +  14.9D0  ,15.3D0  ,17.2D0  ,17.1D0  ,17.5D0  ,17.8D0  ,18.6D0  ,
     +  18.9D0  ,19.1D0  ,19.3D0   /
      DATA (AE(I, 4, 4),I=1,10)  / 
     +  15.1D0  ,13.5D0  ,16.4D0  ,16.7D0  ,16.4D0  ,17.3D0  ,17.8D0  ,
     +  18.5D0  ,19.0D0  ,19.6D0   /
      DATA (AE(I, 4, 5),I=1,10)  / 
     +  .000D+00,15.6D0  ,17.5D0  ,17.7D0  ,17.8D0  ,18.6D0  ,19.2D0  ,
     +  19.9D0  ,20.3D0  ,21.1D0   /
      DATA (AE(I, 4, 6),I=1,10)  / 
     +  .000D+00,18.0D0  ,18.4D0  ,19.2D0  ,19.8D0  ,20.9D0  ,22.0D0  ,
     +  23.1D0  ,23.6D0  ,24.7D0   /
      DATA (AE(I, 4, 7),I=1,10)  / 
     +  .000D+00,27.4D0  ,19.1D0  ,19.8D0  ,20.7D0  ,21.8D0  ,23.2D0  ,
     +  24.4D0  ,24.9D0  ,25.9D0   /
      DATA (AE(I, 4, 8),I=1,10)  / 
     +  .000D+00,.000D+00,18.9D0  ,18.9D0  ,19.3D0  ,21.1D0  ,22.5D0  ,
     +  24.0D0  ,24.7D0  ,26.0D0   /
      DATA (AE(I, 4, 9),I=1,10)  / 
     +  .000D+00,.000D+00,21.1D0  ,19.7D0  ,20.7D0  ,22.3D0  ,24.0D0  ,
     +  25.6D0  ,26.3D0  ,27.7D0   /
      DATA (AE(I, 4,10),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,21.0D0  ,21.1D0  ,22.9D0  ,24.6D0  ,
     +  26.5D0  ,27.3D0  ,29.0D0   /
      DATA (AE(I, 4,11),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,21.3D0  ,22.4D0  ,23.1D0  ,25.0D0  ,
     +  27.1D0  ,27.9D0  ,29.8D0   /
      DATA (AE(I, 4,12),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,36.6D0  ,21.5D0  ,22.2D0  ,23.1D0  ,
     +  25.6D0  ,26.8D0  ,29.1D0   /
      DATA (AE(I, 4,13),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,22.9D0  ,23.1D0  ,23.7D0  ,
     +  26.2D0  ,27.3D0  ,29.6D0   /
      DATA (AE(I, 4,14),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,30.5D0  ,23.6D0  ,25.0D0  ,
     +  26.9D0  ,28.2D0  ,30.7D0   /
      DATA (AE(I, 4,15),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,25.4D0  ,26.2D0  ,
     +  27.2D0  ,28.3D0  ,31.0D0   /
      DATA (AE(I, 4,16),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,24.5D0  ,25.9D0  ,
     +  27.4D0  ,27.6D0  ,30.7D0   /
      DATA (AE(I, 4,17),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,43.3D0  ,28.4D0  ,
     +  27.5D0  ,27.9D0  ,30.9D0   /
      DATA (AE(I, 4,18),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,27.2D0  ,
     +  29.1D0  ,29.0D0  ,31.4D0   /
      DATA (AE(I, 4,19),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,51.3D0  ,
     +  30.6D0  ,29.5D0  ,31.4D0   /
      DATA (AE(I, 4,20),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,
     +  28.8D0  ,30.6D0  ,32.4D0   /
      DATA (AE(I, 5, 1),I=1,10)  / 
     +  15.0D0  ,14.9D0  ,15.5D0  ,15.4D0  ,15.9D0  ,15.8D0  ,16.2D0  ,
     +  16.2D0  ,16.1D0  ,15.9D0   /
      DATA (AE(I, 5, 2),I=1,10)  / 
     +  15.4D0  ,16.1D0  ,17.0D0  ,17.4D0  ,18.0D0  ,18.2D0  ,18.7D0  ,
     +  18.9D0  ,19.0D0  ,19.1D0   /
      DATA (AE(I, 5, 3),I=1,10)  / 
     +  17.1D0  ,17.2D0  ,18.3D0  ,18.7D0  ,19.3D0  ,19.6D0  ,20.3D0  ,
     +  20.6D0  ,20.8D0  ,20.9D0   /
      DATA (AE(I, 5, 4),I=1,10)  / 
     +  14.7D0  ,14.8D0  ,15.0D0  ,16.0D0  ,17.0D0  ,17.7D0  ,18.1D0  ,
     +  19.0D0  ,19.4D0  ,20.0D0   /
      DATA (AE(I, 5, 5),I=1,10)  / 
     +  .000D+00,16.7D0  ,17.6D0  ,18.1D0  ,18.6D0  ,19.2D0  ,19.7D0  ,
     +  20.4D0  ,20.8D0  ,21.2D0   /
      DATA (AE(I, 5, 6),I=1,10)  / 
     +  .000D+00,17.8D0  ,18.2D0  ,19.2D0  ,20.0D0  ,21.0D0  ,21.9D0  ,
     +  23.0D0  ,23.6D0  ,24.3D0   /
      DATA (AE(I, 5, 7),I=1,10)  / 
     +  .000D+00,35.2D0  ,18.9D0  ,20.3D0  ,20.6D0  ,21.5D0  ,22.6D0  ,
     +  23.7D0  ,24.2D0  ,24.7D0   /
      DATA (AE(I, 5, 8),I=1,10)  / 
     +  .000D+00,.000D+00,16.4D0  ,18.9D0  ,18.8D0  ,19.6D0  ,20.7D0  ,
     +  22.3D0  ,23.1D0  ,23.9D0   /
      DATA (AE(I, 5, 9),I=1,10)  / 
     +  .000D+00,.000D+00,33.9D0  ,19.8D0  ,20.3D0  ,20.7D0  ,21.9D0  ,
     +  23.4D0  ,24.1D0  ,24.8D0   /
      DATA (AE(I, 5,10),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,18.0D0  ,20.0D0  ,21.4D0  ,22.0D0  ,
     +  23.8D0  ,24.6D0  ,25.4D0   /
      DATA (AE(I, 5,11),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,26.4D0  ,20.4D0  ,21.2D0  ,22.3D0  ,
     +  23.8D0  ,24.7D0  ,25.5D0   /
      DATA (AE(I, 5,12),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,41.7D0  ,18.2D0  ,19.8D0  ,21.1D0  ,
     +  22.6D0  ,23.4D0  ,24.6D0   /
      DATA (AE(I, 5,13),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,22.5D0  ,20.0D0  ,21.7D0  ,
     +  22.8D0  ,23.7D0  ,24.7D0   /
      DATA (AE(I, 5,14),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,54.1D0  ,19.9D0  ,21.9D0  ,
     +  23.2D0  ,24.3D0  ,25.3D0   /
      DATA (AE(I, 5,15),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,21.2D0  ,22.2D0  ,
     +  23.6D0  ,24.9D0  ,25.5D0   /
      DATA (AE(I, 5,16),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,44.9D0  ,21.9D0  ,
     +  23.8D0  ,25.2D0  ,25.6D0   /
      DATA (AE(I, 5,17),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,47.8D0  ,22.7D0  ,
     +  23.8D0  ,24.9D0  ,26.3D0   /
      DATA (AE(I, 5,18),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,35.5D0  ,
     +  23.9D0  ,25.9D0  ,26.6D0   /
      DATA (AE(I, 5,19),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,64.3D0  ,
     +  24.1D0  ,25.7D0  ,27.1D0   /
      DATA (AE(I, 5,20),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,
     +  34.0D0  ,25.7D0  ,27.7D0   /
      DATA (AE(I, 6, 1),I=1,10)  / 
     +  16.6D0  ,16.5D0  ,16.8D0  ,16.7D0  ,17.0D0  ,16.5D0  ,16.7D0  ,
     +  18.3D0  ,18.9D0  ,19.0D0   /
      DATA (AE(I, 6, 2),I=1,10)  / 
     +  16.2D0  ,16.6D0  ,17.2D0  ,17.4D0  ,17.9D0  ,17.4D0  ,17.7D0  ,
     +  20.7D0  ,22.0D0  ,22.6D0   /
      DATA (AE(I, 6, 3),I=1,10)  / 
     +  18.9D0  ,18.7D0  ,18.8D0  ,18.6D0  ,18.9D0  ,18.6D0  ,18.9D0  ,
     +  21.0D0  ,22.3D0  ,22.9D0   /
      DATA (AE(I, 6, 4),I=1,10)  / 
     +  18.3D0  ,12.7D0  ,14.2D0  ,15.0D0  ,15.7D0  ,16.1D0  ,16.3D0  ,
     +  16.5D0  ,17.9D0  ,19.0D0   /
      DATA (AE(I, 6, 5),I=1,10)  / 
     +  .000D+00,15.7D0  ,15.1D0  ,15.3D0  ,16.5D0  ,16.4D0  ,16.4D0  ,
     +  17.0D0  ,18.3D0  ,19.4D0   /
      DATA (AE(I, 6, 6),I=1,10)  / 
     +  .000D+00,22.9D0  ,14.9D0  ,15.2D0  ,16.2D0  ,16.9D0  ,17.4D0  ,
     +  18.2D0  ,19.5D0  ,21.1D0   /
      DATA (AE(I, 6, 7),I=1,10)  / 
     +  .000D+00,40.7D0  ,18.4D0  ,15.9D0  ,17.1D0  ,17.7D0  ,18.9D0  ,
     +  19.5D0  ,20.3D0  ,21.1D0   /
      DATA (AE(I, 6, 8),I=1,10)  / 
     +  .000D+00,.000D+00,23.3D0  ,16.2D0  ,16.3D0  ,17.3D0  ,18.7D0  ,
     +  19.5D0  ,20.3D0  ,21.1D0   /
      DATA (AE(I, 6, 9),I=1,10)  / 
     +  .000D+00,.000D+00,49.2D0  ,19.0D0  ,19.1D0  ,19.4D0  ,20.2D0  ,
     +  20.8D0  ,21.6D0  ,22.0D0   /
      DATA (AE(I, 6,10),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,27.2D0  ,21.2D0  ,20.8D0  ,21.4D0  ,
     +  22.3D0  ,22.8D0  ,23.3D0   /
      DATA (AE(I, 6,11),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,45.6D0  ,25.0D0  ,22.8D0  ,23.9D0  ,
     +  23.6D0  ,24.3D0  ,24.4D0   /
      DATA (AE(I, 6,12),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,45.8D0  ,29.7D0  ,25.1D0  ,25.3D0  ,
     +  25.3D0  ,26.0D0  ,26.3D0   /
      DATA (AE(I, 6,13),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,42.7D0  ,29.0D0  ,28.0D0  ,
     +  27.0D0  ,27.2D0  ,27.6D0   /
      DATA (AE(I, 6,14),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,62.0D0  ,32.0D0  ,30.0D0  ,
     +  29.8D0  ,29.5D0  ,29.6D0   /
      DATA (AE(I, 6,15),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,44.5D0  ,34.4D0  ,
     +  32.7D0  ,31.5D0  ,31.8D0   /
      DATA (AE(I, 6,16),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,75.6D0  ,37.1D0  ,
     +  34.6D0  ,34.4D0  ,34.4D0   /
      DATA (AE(I, 6,17),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,51.2D0  ,45.2D0  ,
     +  39.0D0  ,37.5D0  ,36.4D0   /
      DATA (AE(I, 6,18),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,74.9D0  ,
     +  42.3D0  ,39.9D0  ,38.3D0   /
      DATA (AE(I, 6,19),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,69.5D0  ,
     +  50.7D0  ,42.3D0  ,41.4D0   /
      DATA (AE(I, 6,20),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,
     +  66.3D0  ,48.0D0  ,43.4D0   /
      DATA (AE(I, 7, 1),I=1,10)  / 
     +  27.0D0  ,25.8D0  ,26.3D0  ,26.2D0  ,26.7D0  ,26.7D0  ,27.1D0  ,
     +  27.1D0  ,27.2D0  ,19.0D0   /
      DATA (AE(I, 7, 2),I=1,10)  / 
     +  29.1D0  ,28.9D0  ,29.7D0  ,30.3D0  ,31.0D0  ,31.4D0  ,32.0D0  ,
     +  32.3D0  ,32.7D0  ,22.6D0   /
      DATA (AE(I, 7, 3),I=1,10)  / 
     +  31.6D0  ,29.7D0  ,30.9D0  ,31.4D0  ,32.5D0  ,33.1D0  ,34.0D0  ,
     +  34.6D0  ,35.1D0  ,22.9D0   /
      DATA (AE(I, 7, 4),I=1,10)  / 
     +  27.4D0  ,19.9D0  ,20.8D0  ,22.8D0  ,24.6D0  ,26.4D0  ,28.2D0  ,
     +  29.6D0  ,30.8D0  ,19.0D0   /
      DATA (AE(I, 7, 5),I=1,10)  / 
     +  .000D+00,24.6D0  ,24.1D0  ,25.0D0  ,27.2D0  ,28.7D0  ,30.7D0  ,
     +  31.8D0  ,32.9D0  ,19.4D0   /
      DATA (AE(I, 7, 6),I=1,10)  / 
     +  .000D+00,35.6D0  ,25.2D0  ,25.6D0  ,27.9D0  ,30.4D0  ,32.7D0  ,
     +  34.6D0  ,36.3D0  ,21.1D0   /
      DATA (AE(I, 7, 7),I=1,10)  / 
     +  .000D+00,45.4D0  ,30.9D0  ,28.2D0  ,29.0D0  ,31.2D0  ,34.0D0  ,
     +  35.8D0  ,37.4D0  ,21.1D0   /
      DATA (AE(I, 7, 8),I=1,10)  / 
     +  .000D+00,.000D+00,38.2D0  ,29.6D0  ,29.4D0  ,30.3D0  ,33.2D0  ,
     +  35.5D0  ,37.6D0  ,21.1D0   /
      DATA (AE(I, 7, 9),I=1,10)  / 
     +  .000D+00,.000D+00,59.3D0  ,34.5D0  ,33.7D0  ,32.9D0  ,35.4D0  ,
     +  37.6D0  ,39.6D0  ,22.0D0   /
      DATA (AE(I, 7,10),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,44.5D0  ,37.8D0  ,37.5D0  ,37.2D0  ,
     +  39.0D0  ,41.4D0  ,23.3D0   /
      DATA (AE(I, 7,11),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,67.0D0  ,43.6D0  ,42.0D0  ,40.8D0  ,
     +  41.4D0  ,43.0D0  ,24.4D0   /
      DATA (AE(I, 7,12),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,49.9D0  ,50.9D0  ,44.6D0  ,43.9D0  ,
     +  44.2D0  ,44.2D0  ,26.3D0   /
      DATA (AE(I, 7,13),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,67.2D0  ,50.5D0  ,48.7D0  ,
     +  48.1D0  ,47.2D0  ,27.6D0   /
      DATA (AE(I, 7,14),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,68.1D0  ,55.2D0  ,52.3D0  ,
     +  51.5D0  ,51.6D0  ,29.6D0   /
      DATA (AE(I, 7,15),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,68.7D0  ,58.6D0  ,
     +  56.5D0  ,55.7D0  ,31.8D0   /
      DATA (AE(I, 7,16),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,89.3D0  ,62.9D0  ,
     +  60.0D0  ,59.1D0  ,34.4D0   /
      DATA (AE(I, 7,17),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,56.0D0  ,72.9D0  ,
     +  66.3D0  ,64.2D0  ,36.4D0   /
      DATA (AE(I, 7,18),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,105.D0  ,
     +  71.3D0  ,68.3D0  ,38.3D0   /
      DATA (AE(I, 7,19),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,73.4D0  ,
     +  76.8D0  ,72.4D0  ,41.4D0   /
      DATA (AE(I, 7,20),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,
     +  107.D0  ,79.9D0  ,43.4D0   /
      DATA (AE(I, 8, 1),I=1,10)  / 
     +  35.5D0  ,35.3D0  ,35.7D0  ,35.7D0  ,36.3D0  ,36.3D0  ,36.7D0  ,
     +  36.7D0  ,36.7D0  ,19.0D0   /
      DATA (AE(I, 8, 2),I=1,10)  / 
     +  40.6D0  ,41.4D0  ,41.9D0  ,42.3D0  ,43.2D0  ,43.5D0  ,44.0D0  ,
     +  44.3D0  ,44.5D0  ,22.6D0   /
      DATA (AE(I, 8, 3),I=1,10)  / 
     +  45.4D0  ,45.7D0  ,46.4D0  ,47.0D0  ,48.1D0  ,48.7D0  ,49.4D0  ,
     +  49.8D0  ,50.2D0  ,22.9D0   /
      DATA (AE(I, 8, 4),I=1,10)  / 
     +  43.9D0  ,44.3D0  ,43.4D0  ,45.1D0  ,47.3D0  ,48.7D0  ,49.6D0  ,
     +  50.5D0  ,51.3D0  ,19.0D0   /
      DATA (AE(I, 8, 5),I=1,10)  / 
     +  .000D+00,49.3D0  ,49.6D0  ,50.5D0  ,53.2D0  ,54.2D0  ,55.4D0  ,
     +  56.1D0  ,56.8D0  ,19.4D0   /
      DATA (AE(I, 8, 6),I=1,10)  / 
     +  .000D+00,59.1D0  ,53.0D0  ,55.4D0  ,58.0D0  ,60.0D0  ,61.2D0  ,
     +  62.5D0  ,63.6D0  ,21.1D0   /
      DATA (AE(I, 8, 7),I=1,10)  / 
     +  .000D+00,54.5D0  ,57.1D0  ,59.2D0  ,62.3D0  ,64.4D0  ,66.0D0  ,
     +  67.3D0  ,68.5D0  ,21.1D0   /
      DATA (AE(I, 8, 8),I=1,10)  / 
     +  .000D+00,.000D+00,65.9D0  ,62.1D0  ,65.1D0  ,67.6D0  ,69.4D0  ,
     +  71.1D0  ,72.6D0  ,21.1D0   /
      DATA (AE(I, 8, 9),I=1,10)  / 
     +  .000D+00,.000D+00,72.2D0  ,67.1D0  ,70.5D0  ,73.1D0  ,75.1D0  ,
     +  76.8D0  ,78.4D0  ,22.0D0   /
      DATA (AE(I, 8,10),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,80.1D0  ,75.0D0  ,78.0D0  ,80.0D0  ,
     +  82.1D0  ,83.9D0  ,23.3D0   /
      DATA (AE(I, 8,11),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,94.5D0  ,82.2D0  ,82.8D0  ,85.1D0  ,
     +  87.3D0  ,89.2D0  ,24.4D0   /
      DATA (AE(I, 8,12),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,56.8D0  ,92.5D0  ,87.2D0  ,89.4D0  ,
     +  91.9D0  ,94.1D0  ,26.3D0   /
      DATA (AE(I, 8,13),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,116.D0  ,96.2D0  ,94.4D0  ,
     +  97.0D0  ,99.2D0  ,27.6D0   /
      DATA (AE(I, 8,14),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,78.1D0  ,104.D0  ,102.D0  ,
     +  102.D0  ,105.D0  ,29.6D0   /
      DATA (AE(I, 8,15),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,128.D0  ,111.D0  ,
     +  109.D0  ,110.D0  ,31.8D0   /
      DATA (AE(I, 8,16),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,104.D0  ,118.D0  ,
     +  117.D0  ,115.D0  ,34.4D0   /
      DATA (AE(I, 8,17),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,64.4D0  ,138.D0  ,
     +  124.D0  ,122.D0  ,36.4D0   /
      DATA (AE(I, 8,18),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,133.D0  ,
     +  133.D0  ,132.D0  ,38.3D0   /
      DATA (AE(I, 8,19),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,83.6D0  ,
     +  146.D0  ,139.D0  ,41.4D0   /
      DATA (AE(I, 8,20),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,
     +  166.D0  ,147.D0  ,43.4D0   /
      DATA (AE(I, 9, 1),I=1,10)  / 
     +  43.3D0  ,43.2D0  ,43.6D0  ,43.8D0  ,44.1D0  ,44.3D0  ,44.7D0  ,
     +  44.8D0  ,44.8D0  ,19.0D0   /
      DATA (AE(I, 9, 2),I=1,10)  / 
     +  50.9D0  ,51.4D0  ,52.0D0  ,52.6D0  ,53.1D0  ,53.6D0  ,54.2D0  ,
     +  54.5D0  ,54.7D0  ,22.6D0   /
      DATA (AE(I, 9, 3),I=1,10)  / 
     +  58.0D0  ,58.4D0  ,59.3D0  ,60.1D0  ,60.7D0  ,61.5D0  ,62.3D0  ,
     +  62.7D0  ,63.1D0  ,22.9D0   /
      DATA (AE(I, 9, 4),I=1,10)  / 
     +  62.0D0  ,63.9D0  ,63.7D0  ,65.7D0  ,65.5D0  ,67.5D0  ,68.2D0  ,
     +  68.9D0  ,69.7D0  ,19.0D0   /
      DATA (AE(I, 9, 5),I=1,10)  / 
     +  .000D+00,72.2D0  ,72.5D0  ,74.2D0  ,74.2D0  ,76.1D0  ,77.0D0  ,
     +  77.8D0  ,78.6D0  ,19.4D0   /
      DATA (AE(I, 9, 6),I=1,10)  / 
     +  .000D+00,80.4D0  ,80.5D0  ,83.1D0  ,83.0D0  ,85.5D0  ,86.8D0  ,
     +  88.1D0  ,89.2D0  ,21.1D0   /
      DATA (AE(I, 9, 7),I=1,10)  / 
     +  .000D+00,63.4D0  ,88.5D0  ,91.3D0  ,91.1D0  ,94.0D0  ,95.8D0  ,
     +  97.3D0  ,98.6D0  ,21.1D0   /
      DATA (AE(I, 9, 8),I=1,10)  / 
     +  .000D+00,.000D+00,98.8D0  ,98.6D0  ,97.8D0  ,102.D0  ,104.D0  ,
     +  106.D0  ,108.D0  ,21.1D0   /
      DATA (AE(I, 9, 9),I=1,10)  / 
     +  .000D+00,.000D+00,84.1D0  ,107.D0  ,107.D0  ,111.D0  ,113.D0  ,
     +  116.D0  ,117.D0  ,22.0D0   /
      DATA (AE(I, 9,10),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,116.D0  ,115.D0  ,119.D0  ,122.D0  ,
     +  125.D0  ,127.D0  ,23.3D0   /
      DATA (AE(I, 9,11),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,111.D0  ,123.D0  ,127.D0  ,131.D0  ,
     +  134.D0  ,137.D0  ,24.4D0   /
      DATA (AE(I, 9,12),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,65.6D0  ,136.D0  ,135.D0  ,140.D0  ,
     +  143.D0  ,146.D0  ,26.3D0   /
      DATA (AE(I, 9,13),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,146.D0  ,144.D0  ,149.D0  ,
     +  152.D0  ,155.D0  ,27.6D0   /
      DATA (AE(I, 9,14),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,88.7D0  ,152.D0  ,158.D0  ,
     +  162.D0  ,165.D0  ,29.6D0   /
      DATA (AE(I, 9,15),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,181.D0  ,167.D0  ,
     +  171.D0  ,174.D0  ,31.8D0   /
      DATA (AE(I, 9,16),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,117.D0  ,174.D0  ,
     +  180.D0  ,183.D0  ,34.4D0   /
      DATA (AE(I, 9,17),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,72.0D0  ,201.D0  ,
     +  189.D0  ,192.D0  ,36.4D0   /
      DATA (AE(I, 9,18),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,151.D0  ,
     +  198.D0  ,201.D0  ,38.3D0   /
      DATA (AE(I, 9,19),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,95.2D0  ,
     +  220.D0  ,210.D0  ,41.4D0   /
      DATA (AE(I, 9,20),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,
     +  192.D0  ,217.D0  ,43.4D0   /
      DATA (AE(I,10, 1),I=1,10)  / 
     +  62.1D0  ,62.1D0  ,62.6D0  ,62.9D0  ,63.3D0  ,63.3D0  ,64.0D0  ,
     +  64.0D0  ,64.0D0  ,19.0D0   /
      DATA (AE(I,10, 2),I=1,10)  / 
     +  75.1D0  ,75.4D0  ,76.3D0  ,76.8D0  ,77.6D0  ,77.9D0  ,78.8D0  ,
     +  79.0D0  ,79.3D0  ,22.6D0   /
      DATA (AE(I,10, 3),I=1,10)  / 
     +  87.5D0  ,88.3D0  ,89.4D0  ,90.2D0  ,91.3D0  ,91.9D0  ,93.0D0  ,
     +  93.5D0  ,93.9D0  ,22.9D0   /
      DATA (AE(I,10, 4),I=1,10)  / 
     +  104.D0  ,104.D0  ,105.D0  ,106.D0  ,107.D0  ,108.D0  ,109.D0  ,
     +  110.D0  ,110.D0  ,19.0D0   /
      DATA (AE(I,10, 5),I=1,10)  / 
     +  .000D+00,122.D0  ,122.D0  ,123.D0  ,124.D0  ,125.D0  ,126.D0  ,
     +  127.D0  ,128.D0  ,19.4D0   /
      DATA (AE(I,10, 6),I=1,10)  / 
     +  .000D+00,138.D0  ,139.D0  ,140.D0  ,142.D0  ,143.D0  ,144.D0  ,
     +  146.D0  ,147.D0  ,21.1D0   /
      DATA (AE(I,10, 7),I=1,10)  / 
     +  .000D+00,85.3D0  ,158.D0  ,159.D0  ,161.D0  ,162.D0  ,164.D0  ,
     +  166.D0  ,167.D0  ,21.1D0   /
      DATA (AE(I,10, 8),I=1,10)  / 
     +  .000D+00,.000D+00,176.D0  ,177.D0  ,179.D0  ,181.D0  ,183.D0  ,
     +  184.D0  ,186.D0  ,21.1D0   /
      DATA (AE(I,10, 9),I=1,10)  / 
     +  .000D+00,.000D+00,114.D0  ,199.D0  ,201.D0  ,202.D0  ,205.D0  ,
     +  206.D0  ,207.D0  ,22.0D0   /
      DATA (AE(I,10,10),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,218.D0  ,219.D0  ,220.D0  ,224.D0  ,
     +  225.D0  ,226.D0  ,23.3D0   /
      DATA (AE(I,10,11),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,150.D0  ,238.D0  ,238.D0  ,243.D0  ,
     +  244.D0  ,245.D0  ,24.4D0   /
      DATA (AE(I,10,12),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,85.8D0  ,255.D0  ,255.D0  ,261.D0  ,
     +  262.D0  ,263.D0  ,26.3D0   /
      DATA (AE(I,10,13),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,195.D0  ,272.D0  ,279.D0  ,
     +  279.D0  ,280.D0  ,27.6D0   /
      DATA (AE(I,10,14),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,115.D0  ,290.D0  ,296.D0  ,
     +  297.D0  ,298.D0  ,29.6D0   /
      DATA (AE(I,10,15),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,263.D0  ,313.D0  ,
     +  314.D0  ,315.D0  ,31.8D0   /
      DATA (AE(I,10,16),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,150.D0  ,330.D0  ,
     +  331.D0  ,332.D0  ,34.4D0   /
      DATA (AE(I,10,17),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,90.0D0  ,319.D0  ,
     +  349.D0  ,349.D0  ,36.4D0   /
      DATA (AE(I,10,18),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,196.D0  ,
     +  366.D0  ,367.D0  ,38.3D0   /
      DATA (AE(I,10,19),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,122.D0  ,
     +  387.D0  ,384.D0  ,41.4D0   /
      DATA (AE(I,10,20),I=1,10)  / 
     +  .000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,.000D+00,
     +  247.D0  ,401.D0  ,43.4D0     /
      DATA (ERES(I, 1),I=1,10)  / 10*0.D0/
      DATA (ERES(I, 2),I=1,10)  / 10*0.D0/
      DATA (ERES(I, 3),I=1,10)  / 10*0.D0/
      DATA (ERES(I, 4),I=1,10)  / 10*0.D0/
      DATA (ERES(I, 5),I=1,10)  / 10*0.D0/
      DATA (ERES(I, 6),I=1,10)  / 
     +    0.000D0, 0.000D0, 0.000D0, 0.000D0, 0.000D0, 0.000D0, 0.000D0,
     +    2.780D0, 2.880D0, 2.890D0 /
      DATA (ERES(I, 7),I=1,10)  / 
     +    1.500D0, 2.460D0, 2.510D0, 2.610D0, 2.700D0, 2.920D0, 3.070D0,
     +    3.200D0, 3.330D0, 2.890D0 /
      DATA (ERES(I, 8),I=1,10)  / 
     +    4.470D0, 4.350D0, 4.390D0, 4.550D0, 4.660D0, 4.890D0, 4.980D0,
     +    5.100D0, 5.220D0, 2.890D0 /
      DATA (ERES(I, 9),I=1,10)  / 
     +    7.480D0, 7.380D0, 7.370D0, 7.480D0, 7.510D0, 7.630D0, 7.660D0,
     +    7.750D0, 7.820D0, 2.890D0 /
      DATA (ERES(I,10),I=1,10)  / 
     +   15.270D0,15.190D0,15.200D0,15.370D0,15.380D0,15.430D0,15.540D0,
     +   15.590D0,15.630D0, 2.890D0 /
      END
C->
C=======================================================================

      SUBROUTINE FRAGM (IAT,IAP, NW,B, NF, IAF)

C-----------------------------------------------------------------------
C...Nuclear Fragmentation, Abrasion-ablation model, 
C...Based on Jon Engel's routines ABRABL 
C...This most recent version adds for all prefragment
C...masses > 10 the model calculation for the fragment
C...mass distribution and the energy carried by the fragment
C...of W. Friedmann
C...The average values are used to implement the model
C...in the montecarlo fashion / TSS, Dec '91
C.
C.  INPUT: IAP = mass of incident nucleus
C.         IAT = mass of target   nucleus
C.         NW = number of wounded nucleons in the beam nucleus
C.         B  = impact parameter in the interaction
C.     
C.  OUTPUT : NF = number of fragments  of the spectator nucleus
C.           IAF(1:NF) = mass number of each fragment
C.           PF(3,60) in common block /FRAGMENTS/ contains
C.           the three momentum components (MeV/c) of each
C.           fragment in the projectile frame
C..............................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      COMMON /FRAGMENTS/ PPP(3,60)
      COMMON /FRAGMOD/A(10,10,20),AE(10,10,20),ERES(10,10),NFLAGG(10,10)
      DIMENSION IAF(60)
      DIMENSION AA(10), EAA(10) 
      SAVE
      EXTERNAL GASDEV
      DATA AA/10.D0,15.D0,20.D0,25.D0,30.D0,35.D0,40.D0,45.D0,50.D0,
     $        56.D0/
      DATA EAA/1.D0,2.D0,4.D0,6.D0,8.D0,10.D0,12.D0,16.D0,20.D0,30.D0/

      AP=IAP
      AT=IAT
      NPF = IAP - NW
      IF (NPF .EQ. 0) THEN
         NF = 0
         RETURN
      ENDIF

      EB = ESTAR(AP,AT, B)
      EBP = ESTARP (NPF, NW)
C CONTRIBUTION TO E* FROM ENERGY DEPOSITED BY SECONDARIES
      EB = EB + EBP
C TOTAL E* IS THE SUM OF THE TWO COMPONENTS

C.....Prefragment transverse momentum (MeV/nucleon)...
            FK = FERMK(AP)
C FERMI MOMENTUM OF THE PROJECTILE NUCLEUS
            IF (NW .LT. IAP) THEN
            SIG = FK*DSQRT(NW*NPF/(AP-1.D0))/3.162D0
C GAUSSIAN SIGMA IN ALL THREE DIRECTION
            ELSE
            SIG = FK/3.162D0
C THIS IS NOT CORRECT, TOO LARGE !!!!!!!!!!!!!!
            ENDIF
             PPFX = SIG*GASDEV(0)/NPF
             PPFY = SIG*GASDEV(1)/NPF
C THREE MOMENTUM COMPONENTS PER NUCLEON FOR THE PREFRAGMENT

C.............Crude model for small prefragment mass .......
            IF (NPF .LT. 10) THEN
                 CALL EVAP(NPF, EB, EPS, NNUC, NALP)
C   EPS IS THE KINETIC ENERGY CARRIED BY THE EVAPORATED NUCLEONS
               ETOT = 938.D0 + EPS
                 PP = SQRT((ETOT*ETOT - 8.79844D5)/3.D0)
C   AVERAGE MOMENTUM OF EVAPORATED NUCLEONS IN EACH DIRECTION
                 NUC = NPF - NNUC - 4*NALP
                 NF = 0
                 IF (NUC .GT. 0) THEN
                    NF = NF + 1
                    IAF(NF) = NUC
                    PPP(1,NF) = NUC*PPFX
                    PPP(2,NF) = NUC*PPFY
                 ENDIF
                 IF (NALP .NE. 0) THEN
                 DO I=1,NALP
                   NF = NF + 1
                    IAF(NF) = 4
                   CALL SINCO(S1,C1)
                   CALL SINCO(S2,C2)
                   PXE = 4.D0*PP*S1*S2
                   PYE = 4.D0*PP*S1*C2
                   PPP(1,NF) = 4.D0*PPFX + PXE
                   PPP(2,NF) = 4.D0*PPFY + PYE
                   PPP(1,1) = PPP(1,1) - PXE
                   PPP(2,1) = PPP(2,1) - PYE
                 ENDDO
                 ENDIF
                 IF (NNUC .NE. 0) THEN
                 DO I=1,NNUC
                    NF = NF + 1
                    IAF(NF) = 1
                    CALL SINCO(S1,C1)
                    CALL SINCO(S2,C2)
                    PXE = PP*S1*S2
                    PYE = PP*S1*C2
                    PPP(1,NF) = 4.D0*PPFX + PXE
                    PPP(2,NF) = 4.D0*PPFY + PYE
                    PPP(1,1) = PPP(1,1) - PXE
                    PPP(2,1) = PPP(2,1) - PYE
                 ENDDO
                 ENDIF
                 RETURN
            ENDIF

C.........More refined model calculation .............
      JA = NPF/5 -1
      IF (JA .LT. 10) THEN
      IF ((NPF - AA(JA)) .GT. (AA(JA+1)-NPF)) JA = JA + 1
      ENDIF
      ARAT = DBLE(NPF)/AA(JA)
      DO J=1,10
      IF (EB .LT. EAA(J)) GO TO 29
      ENDDO
      JE = 10
      GO TO 39
   29      JE = J
   39      IF (JE .GT. 1 .AND. JE .NE. 10) THEN
      IF ((EB - EAA(J-1)) .LT. (EAA(J)-EB)) JE = J - 1
      ENDIF
      ERAT = EB/EAA(JE)
        IF (EB .LT. 1.D0) THEN
        ERAT = EB
        ENDIF
C INTERPOLATE BETWEEN EB=0. (NOTHING HAPPENS) AND EB = 1. MeV
         IF (JA .EQ. 10 .AND. JE .GT. 6) THEN
            WRITE(*,*)' JA=',JA,',   JE=',JE
         ENDIF
   43    ESUM = 0.D0
      NSUM = 0
      JF = 0
      DO J=20,1,-1
        FR =  A(JA, JE, J)*ARAT*ERAT
        N1 = INT(1.D0 + FR)
        FR1 = FR/DBLE(N1)
        DO K=1, N1
          IF (S_RNDM(0) .LT. FR1) THEN
            JF = JF + 1
            IAF(JF) = J
            NSUM = NSUM + J
            EKIN = ERAT*AE(JA,JE, J)
            IF (EKIN .GT. 0.D0) THEN
              ESUM = ESUM + EKIN
              ETOT = 938.D0*IAF(JF) + EKIN
              PP = DSQRT(2.D0*(ETOT*ETOT - IAF(JF)**2*8.79844D5)/3.D0)
              CALL SINCO(S1,C1)
              CALL SINCO(S2,C2)
              PPP(1,JF) = PP*S1*S2 + IAF(JF)*PPFX
              PPP(2,JF) = PP*S1*C2 + IAF(JF)*PPFY
            ENDIF
            IF (NSUM .GT. NPF) THEN
C           WRITE(*,*)' WARNING, NSUM=', NSUM,',  NPF=',NPF
C           WRITE(*,*)'  ARAT =', ARAT
              GO TO 43
            ELSE
              IF (NSUM .EQ. NPF) THEN
                GO TO 44
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      IF (NFLAGG(JA,JE) .EQ. 0) THEN
C 'THE RESIDUE' IS A NUCLEAR FRAGMENT
        JF = JF + 1
        IAF(JF) = NPF - NSUM
        F1 = NPF*EB - ESUM
        IF (F1 .LT. 0.D0) F1 = 0.D0
C GIVE THE REST OF EB TO THE FRAGMENT
        EKIN = F1
        IF (EKIN .GT. 0.D0) THEN
          ETOT = 938.D0*IAF(JF) + EKIN
          PP = DSQRT(2.D0*(ETOT*ETOT - IAF(JF)**2*8.79844D5)/3.D0)
          CALL SINCO(S1,C1)
          CALL SINCO(S2,C2)
          PPP(1,JF) = PP*S1*S2 + IAF(JF)*PPFX
          PPP(2,JF) = PP*S1*C2 + IAF(JF)*PPFY
        ENDIF
      ELSE
C 'THE RESIDUE' CONSISTS OF SPECTATOR NUCLEONS
        N1 = NPF - NSUM
        DO K=1,N1
          JF = JF + 1
          IAF(JF) = 1
          EKIN = ERAT*ERES(JA,JE)
          IF (EKIN .GT. 0.D0) THEN
            ETOT = 938.D0*IAF(JF) + EKIN
            PP = DSQRT(2.D0*(ETOT*ETOT - IAF(JF)**2*8.79844D5)/3.D0)
            CALL SINCO(S1,C1)
            CALL SINCO(S2,C2)
            PPP(1,JF) = PP*S1*S2 + PPFX
            PPP(2,JF) = PP*S1*C2 + PPFY
          ENDIF
        ENDDO
      ENDIF
  44  NF = JF
      RETURN
      END
C->
C=======================================================================

      FUNCTION ESTARP (NPF, NW)

C-----------------------------------------------------------------------
C CONTRIBUTION TO E* FROM ENERGY DEPOSITED BY SECONDARIES
C VERY NAIVE VERSION INCORPORATING HUEFFNER'S IDEAS
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      SAVE

      APF = NPF
      F1 = 15.3D0/APF**0.666666666D0
C AVERAGE KINETIC ENERGY/NUCLEON IN PREFRAGMENT (MeV)
C PER PATHLENGTH EQUAL TO THE PREFRAGMENT RADIUS
      ESTARP = 0.D0
      DO I=1,NW
        IF (S_RNDM(0) .GT. 0.5D0) THEN
          F2 = F1*RDIS(0)
          ESTARP = ESTARP + F2
        ENDIF
      ENDDO
C SAMPLE RANDOMLY PER WOUNDED NUCLEON, x NW
      RETURN
      END
C=======================================================================
      
      FUNCTION RDIS(Idum)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      dimension probr(20)
      SAVE
      data probr/
     *      0.10000D0, 0.15748D0, 0.21778D0, 0.28605D0, 0.36060D0,
     *      0.43815D0, 0.51892D0, 0.60631D0, 0.70002D0, 0.79325D0,
     *      0.88863D0, 0.98686D0, 1.10129D0, 1.21202D0, 1.32932D0,
     *      1.44890D0, 1.57048D0, 1.70139D0, 1.83417D0, 2.00000D0/

      rdis = idum
      nr = INT(20.D0*S_RNDM(0) + 1.D0)
      if (nr .eq. 1) then
        f1 = 0.D0
      else
        f1 = probr(nr-1)
      endif
      dr = probr(nr) - f1
      rdis = f1 + dr*S_RNDM(1)
      return
      end

C=======================================================================

      FUNCTION ESTAR(ap,at,b)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      SAVE

c      real*4 ap,at,b,estar
      sigma=4.5D0  !total n-n cross section in fm**2
      rt=.82d0*at**.33333333D0 !target radius
      rp=.82d0*ap**.33333333D0 !projectile radius
      alpha=rt**2/rp**2
      beta=b**2/rt**2
      f=at*sigma/(PI*rt**2)
      alf = log(f)
      alalf = log(alpha)
      gfac=0.d0
      gfac1=0.d0
      s1=0.D0
      s2=0.D0
      s3=0.D0
      ii=1
      do n=0,10 ! This limit may not need to be so high.
         if(n.ge.2) then
            gfac1=gfac
            gfac=gfac+log(float(n)) 
         endif
         g0=n*alf -n*beta*alpha/(n+alpha)+alalf
         g1=g0-log(alpha+n)-gfac
         g2=(n+2)*log(f)-(n+2)*beta*alpha/(n+2+alpha) 
     >      +log(n+2+alpha+beta*alpha**2)-3.d0*log(n+2.d0+alpha)-gfac
         g3=g0-2.d0*log(n+alpha)-gfac1
         ii=-ii
         s1=s1+ii*exp(g1)
         s2=s2+ii*exp(g2)
         if(n.ge.1) s3=s3+ii*exp(g3)
      enddo

      pb=s1
      e1b=197.D0**2/(2.D0*938.d0*rp**2*pb) *s2
c      a=b*(s3/pb-1)
c      a=-b*s3/pb
c      e2b=-.5* 938. * (41./(ap**.333))**2 * a**2 /(197.**2)
c      estar=e1b+e2b
      estar = e1b
      return
      end
C=======================================================================

      SUBROUTINE EVAP(npf,eb,eps,nnuc,nalp)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      SAVE

      eps=7.5D0+sqrt(8.D0*eb)
      n=min(npf*int(eb/eps),npf)
      nalp=n/5
      nnuc=n-4*nalp
      return
      end
C->
C=======================================================================

      FUNCTION FERMK(A)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      DIMENSION AA(6), FK(6)
      SAVE
      DATA AA/4.D0, 6.D0, 12.D0, 24.D0, 40.D0, 57.D0/
      DATA FK/130.D0,169.D0,221.D0,235.D0,251.D0,260.D0/

      DO I=2,4
      IF (A .LT. AA(I)) GO TO 25
      ENDDO
      I = 5
   25      F11 = AA(I-1)
      F12 = AA(I)
      F13 = AA(I+1)
      F21 = FK(I-1)
      F22 = FK(I)
      F23 = FK(I+1)
      FERMK = QUAD_INT(A,F11,F12,F13, F21,F22,F23)
      RETURN
      END

C*=======================================================================
C. Multiple interaction structure
C========================================================================

      SUBROUTINE INT_NUC (IA, IB, SIG0, SIGEL) 

C-----------------------------------------------------------------------
C...Compute with a montecarlo code  the  "multiple interaction structure"
C.  of a nucleus-nucleus interaction
C.
C.  INPUT : IA            = mass of target nucleus
C.          IB            = mass of projectile nucleus
C.          SIG0 (mbarn)  = inelastic pp cross section
C.          SIGEL(mbarn)  = elastic pp cross section
C.
C.  OUTPUT : in common block /CNUCMS/
C.           B = impact parameter (fm)
C.           BMAX = maximum impact parameter for generation
C.           NTRY = number of "trials" before one interaction
C.           NA = number of wounded nucleons in A
C.           NB =    "        "        "     in B
C.           NI = number of nucleon-nucleon inelastic interactions 
C.           NAEL = number of elastically scattered nucleons in  A 
C.           NBEL =    "         "           "          "    in  B
C.           JJA(J)  [J=1:IA]   = number of inelastic interactions 
C.                                of J-th nucleon of nucleus A
C.           JJB(J)  [J=1:IB]   = number of inelastic interactions 
C.                                of J-th nucleon of nucleus B
C.           JJAEL(J)  [J=1:IA]   = number of elastic interactions 
C.                                of J-th nucleon of nucleus A
C.           JJBEL(J)  [J=1:IB]   = number of elastic interactions 
C.                                of J-th nucleon of nucleus B
C.           JJINT(J,K)  [J=1:NB, K=1:NA]  (0 = no interaction) 
C.                                         (1 = interaction )
C.                                         between nucleon J of A and K of B
C-----------------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      PARAMETER (IAMAX=56)
      COMMON /CNUCMS/ B, BMAX, NTRY, NA, NB, NI, NAEL, NBEL
     +         ,JJA(IAMAX), JJB(IAMAX), JJINT(IAMAX,IAMAX)
     +         ,JJAEL(IAMAX), JJBEL(IAMAX)
      DIMENSION XA(IAMAX), YA(IAMAX), XB(IAMAX), YB(IAMAX)
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN
      SAVE

      SIGT = SIG0 + SIGEL
      R2  = 0.1D0 * SIG0/PI
      R2T = 0.1D0 * SIGT/PI
      BMAX = 15.D0                             ! fm
      NTRY = 0
      CALL NUC_CONF (IA, XA, YA)
      CALL NUC_CONF (IB, XB, YB)
      NI = 0
      NIEL = 0
      DO JA=1,IA
         JJA(JA) = 0
         JJAEL(JA) = 0
      ENDDO
      DO JB=1,IB
         JJB(JB) = 0
         JJBEL(JB) = 0
         DO JA=1,IA
            JJINT(JB,JA) = 0
         ENDDO
      ENDDO
1000  B = BMAX*SQRT(S_RNDM(0))
      PHI = TWOPI*S_RNDM(1)
      BX = B*COS(PHI)
      BY = B*SIN(PHI)
      NTRY = NTRY+1
      DO JA=1,IA
         DO JB=1,IB
            S = (XA(JA)-XB(JB)-BX)**2 + (YA(JA)-YB(JB)-BY)**2
            IF (S .LT. R2)  THEN
               NI = NI + 1
               JJA(JA) = JJA(JA)+1
               JJB(JB) = JJB(JB)+1
               JJINT(JB,JA) = 1
            ELSE IF (S .LT. R2T)  THEN
               NIEL = NIEL + 1
               JJAEL(JA) = JJAEL(JA)+1
               JJBEL(JB) = JJBEL(JB)+1
            ENDIF
         ENDDO
      ENDDO
      IF (NI + NIEL .EQ. 0)  GOTO 1000
      NA = 0
      NB = 0
      NAEL = 0
      NBEL = 0
      DO JA=1,IA
         IF (JJA(JA) .GT. 0)  THEN
            NA = NA + 1
         ELSE
            IF (JJAEL(JA) .GT. 0)  NAEL = NAEL+1
         ENDIF
      ENDDO
      DO JB=1,IB
         IF (JJB(JB) .GT. 0)  THEN
            NB = NB + 1
         ELSE
            IF (JJBEL(JB) .GT. 0)  NBEL = NBEL+1
         ENDIF
      ENDDO
      RETURN
      END
C=======================================================================

       SUBROUTINE NUC_CONF (IA, XX, YY)

C-----------------------------------------------------------------------
C...This routine generates the configuration  of a nucleus 
C.  need an initialization call to NUC_GEOM_INI
C.
C.  INPUT  : IA = mass number of the nucleus
C.  OUTPUT : XX(1:IA), YY(1:IA) (fm) = position in impact parameter
C.                                     space of the IA nucleons
C...................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      PARAMETER (IAMAX=56)
      DIMENSION XX(IAMAX), YY(IAMAX)
      PARAMETER (NB=401)
      COMMON /CPROFA/ ZMIN, DZ, BBZ(NB,IAMAX)
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN
      SAVE

      DO J=1,IA
         Z = S_RNDM(J)
         JZ = INT((Z-ZMIN)/DZ)+1
         JZ = MIN(JZ,400)
         T = (Z-ZMIN)/DZ - DBLE(JZ-1)
         B = BBZ(JZ,IA)*(1.D0-T) + BBZ(JZ+1,IA)*T
         PHI = TWOPI*S_RNDM(J+1)
         XX(J) = B*COS(PHI)
         YY(J) = B*SIN(PHI)
      ENDDO
      RETURN
      END
C=======================================================================

      SUBROUTINE NUC_GEOM_INI

C-----------------------------------------------------------------------
C...Initialize all nucleus profiles
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      PARAMETER (NB=401)
      PARAMETER (IAMAX=56)
      COMMON /CPROF/ DB, BMAX, BB(NB), TB(NB), A
      COMMON /CPROFA/ ZMIN, DZ, BBZ(NB,IAMAX)
      DIMENSION FFB(NB), GGB(NB)
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN
      SAVE

      CALL SHELL_INI
      CALL WOOD_SAXON_INI
      DO IA= 2,IAMAX
           JA = IA
         CALL NUC_PROFIL(JA)
         DO K=1,NB
           FFB(K) = BB(K)*TB(K) * TWOPI
         ENDDO            
         GGB(1) = 0.D0
         GGB(NB) = 1.D0
         DO K=2,NB-1
           GGB(K) = GGB(K-1) + FFB(K-1)*DB
         ENDDO            
         CALL INVERT_ARRAY(GGB,0.D0,DB,NB, BBZ(1,IA), ZMIN, DZ)
      ENDDO
      RETURN
      END
C=======================================================================

      SUBROUTINE NUC_PROFIL (JA)

C-----------------------------------------------------------------------
C...Compute the profile function T(b)
C.  normalised as INT[d2b T(b) = 1]
C.  INPUT : JA = integer mass number of nucleus
C...............................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      PARAMETER (NB=401)
      EXTERNAL DENSA
      DOUBLE PRECISION DENSA
      COMMON /CC01/  B
      COMMON /CCDA/ JJA
      COMMON /CPROF/ DB, BMAX, BB(NB), TB(NB), A
      SAVE

      BMAX = 7.5D0
      DB = BMAX/DBLE(NB-1)
      JJA = JA
      A = JA
      DO JB=1,NB
        B = DB*DBLE(JB-1)
        BB(JB) = B
        IF (JA .LE. 18)  THEN
            TB(JB) = PROFNUC (B, JA)
         ELSE
            TB(JB) = 2.D0*GAUSS (DENSA,0.D0,BMAX)
         ENDIF
      ENDDO
      RETURN
      END
C=======================================================================

      SUBROUTINE NUC1_PROFIL (AA)

C-----------------------------------------------------------------------
C...Compute the profile function T(b)
C.  normalised as INT[d2b T(b) = 1]
C.  INPUT : AA = mass number of nucleus
C...............................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      PARAMETER (NB=401)
      EXTERNAL DENSA
      DOUBLE PRECISION DENSA
      COMMON /CC01/  B
      COMMON /CPROF/ DB, BMAX, BB(NB), TB(NB), A
      SAVE

      A = AA
      IA1 = INT(AA)
      IA2 = IA1 + 1
      U = AA - DBLE(IA1)
      BMAX = 7.5D0
      DB = BMAX/DBLE(NB-1)
      DO JB=1,NB
         B = DB*DBLE(JB-1)
         BB(JB) = B
         IF (A .LE. 18.D0)  THEN
             T1 = PROFNUC (B, IA1)
             T2 = PROFNUC (B, IA2)
          ELSE
             JJA = IA1
             T1 = 2.D0*GAUSS (DENSA,0.D0,BMAX)
             JJA = IA2
             T2 = 2.D0*GAUSS (DENSA,0.D0,BMAX)
          ENDIF
          TB(JB) = (1.D0-U)*T1  + U*T2
      ENDDO
      RETURN
      END

C*======================================================================
C.   Code about nuclear densities
C=======================================================================

      FUNCTION DENS_NUC (R, JA)

C-----------------------------------------------------------------------
C....Nuclear density (normalised to 1)
C.   for a nucleus of mass number JA
C.   INPUT R = radial coordinate  (fm)
C.         JA = integer mass number
C.  OUTPUT (fm**-3)
C--------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON /CWOOD/ RR0(19:56), AA0(19:56), CC0(19:56)
      SAVE

      IF (JA .GT. 18)  THEN
         DENS_NUC = WOOD_SAXON(R,JA)
      ELSE IF (JA .NE. 4)  THEN
         DENS_NUC = HELIUM(R)
      ELSE
         DENS_NUC = SHELL(R,JA)
      ENDIF
      RETURN
      END
C=======================================================================

      FUNCTION WOOD_SAXON (R, JA) 

C-----------------------------------------------------------------------
C....Wood-Saxon nuclear density (normalised to 1)
C.   for a nucleus of mass number A.
C.   INPUT R =  (fm)
C.         JA = mass number
C.   OUTPUT (fm**-3)
C------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON /CWOOD/ RR0(19:56), AA0(19:56), CC0(19:56)
      SAVE

      WOOD_SAXON = CC0(JA)/(1.D0+EXP((R-RR0(JA))/AA0(JA)))
      RETURN
      END      
C=======================================================================

      FUNCTION HELIUM (R)

C-----------------------------------------------------------------------
C... Helium density from Barrett and Jackson
C.   INPUT R = r coordinate (fm)
C.   OUTPUT (fm**-3)
C........................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      SAVE
      DATA R0 /0.964D0/, CA /0.322D0/   ! fm
      DATA W /0.517D0/, CC /5.993224D-02/

      HELIUM = CC*(1.D0+W*(R/R0)**2)/(1.D0 + EXP((R-R0)/CA))
      RETURN
      END
C=======================================================================

      FUNCTION SHELL (R,JA)

C-----------------------------------------------------------------------
C...Density in the shell model
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON /CSHELL/ RR0(18), RR02(18)
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN
      SAVE

      R0 = RR0(JA)
      C1 = MIN(1.D0,4.D0/DBLE(JA))
      CS = 1.D0/(R0**3 * PI**1.5D0)
      CP = 2.D0*CS/3.D0
      FS = EXP(-(R/R0)**2)
      FP = (R/R0)**2 * FS
      SHELL = C1*CS*FS + (1.D0-C1)*CP*FP
      RETURN
      END
C=======================================================================

      FUNCTION PROFNUC (B, JA)

C-----------------------------------------------------------------------
C...This function return
C.  the profile T(b) for a nucleus of mass number A
C.  INPUT B = impact parameter (GeV**-1)
C.        JA = integer mass number
C.  OUTPUT  (fm**-2)
C.
C.  The  density of the nucleus is the `shell model density'
C.  the parameter r0 must beinitialized in the common block
C.............................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON /CSHELL/ RR0(18), RR02(18)
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN
      SAVE

      B2 = B*B
      ARG = B2/RR02(JA)
      TS = EXP(-ARG)
      TP = TS*(2.D0*B2+RR02(JA))/(3.D0*RR02(JA))
      CS = MIN(1.D0,4.D0/DBLE(JA))
      PROFNUC = (CS*TS + (1.D0-CS)*TP)/(PI*RR02(JA))
      RETURN
      END
C=======================================================================

      SUBROUTINE SHELL_INI

C-----------------------------------------------------------------------
C...Initialize the parameter  of the shell model
C.  for the nuclei with    6 < A < 18
C..............................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      COMMON /CSHELL/ RR0(18), RR02(18)
      DIMENSION RR(18)
      SAVE
C...Data on Sqrt[<r**2>]  in fermi
      DATA RR /0.81D0, 2.095D0,  1.88D0, 1.674D0,  -1.D0,
     +         2.56D0, 2.41D0,    -1.D0, 2.519D0, 2.45D0,
     +         2.37D0, 2.460D0, 2.440D0,  2.54D0, 2.58D0, 
     +         2.718D0,2.662D0, 2.789D0/

      DO JA=1,18
         A = DBLE(JA)
         RMED = RR(JA)
         IF (RMED .LE. 0.D0)   RMED = 0.5D0*(RR(JA-1) + RR(JA+1))
         C = MAX(1.5D0,(5.D0/2.D0 - 4.D0/A) )
         R0 = RMED/SQRT(C)
         RR0 (JA) = R0
         RR02(JA) = R0*R0
      ENDDO
      RETURN
      END
C->
C=======================================================================

      SUBROUTINE WOOD_SAXON_INI

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON /CWOOD/ RR0(19:56), AA0(19:56), CC0(19:56)
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN
      SAVE

C...Wood-Saxon parameters from  table 6.2   of Barrett and Jackson
      RR0 (19) = 2.59D0
      AA0 (19) = 0.564D0
      RR0 (20) = 2.74D0
      AA0 (20) = 0.569D0
      RR0 (22) = 2.782D0
      AA0 (22) = 0.549D0
      RR0 (24) = 2.99D0
      AA0 (24) = 0.548D0
      RR0 (27) = 2.84D0
      AA0 (27) = 0.569D0
      RR0 (28) = 3.14D0
      AA0 (28) = 0.537D0
      RR0 (29) = 3.77D0
      AA0 (29) = 0.52D0
      RR0 (48) = 3.912D0
      AA0 (48) = 0.5234D0
      RR0 (56) = 3.98D0
      AA0 (56) = 0.569D0
      DO J=19, 56
         IF (RR0(J) .LE. 0.D0)  THEN
            RR0(J) = 1.05D0*DBLE(J)**0.333333333333D0
            AA0(J) = 0.545D0
         ENDIF
         CC0(J)=3.D0/(4.D0*PI*RR0(J)**3)/(1.D0+((AA0(J)*PI)/RR0(J))**2)
      ENDDO
      RETURN
      END
C=======================================================================

      FUNCTION DENSA (Z)

C-----------------------------------------------------------------------
C....Woods Saxon nuclear density (normalised to 1)
C.   for a nucleus of mass number A.
C.   INPUT z = z coordinate (fm)
C.         JA = integer mass number
C.         B (in common /CC01/)  impact parameter  (fm)
C.  OUTPUT (fm**-3)
C--------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON /CC01/  B
      COMMON /CCDA/ JA
      COMMON /CWOOD/ RR0(19:56), AA0(19:56), CC0(19:56)
      SAVE

      R = SQRT (Z*Z + B*B)
      DENSA = CC0(JA)/(1.D0+EXP((R-RR0(JA))/AA0(JA)))
      RETURN
      END

C*=====================================================================
C. Cross sections
C======================================================================

      SUBROUTINE SIGMA_AIR (IB,SIG0,SIGEL,KINT,
     +                            SIGMA,DSIGMA,SIGQE,DSIGQE)

C-----------------------------------------------------------------------
C...Compute with a montecarlo method the "production"
C.  and "quasi-elastic" cross section for  
C.  a nucleus-air  interaction 
C.
C.  INPUT : IB            = mass of projectile nucleus
C.          SIG0 (mbarn)  = inelastic pp cross section
C.          KINT            = number  of interactions to generate
C.  OUTPUT : SIGMA (mbarn) = "production" cross section
C.           DSIGMA   "    = error
C.           SIGQE    "    = "quasi-elastic" cross section
C.           DSIGQE   "    = error
C.           additional output is in the common block  /CPROBAB/
C..........................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      PARAMETER (IAMAX=56)
      PARAMETER (IAMAX2=3136)          ! IAMAX*IAMAX
      COMMON  /CPROBAB/ PROBA(IAMAX), DPROBA(IAMAX), 
     +   PROBB(IAMAX), DPROBB(IAMAX), PROBI(IAMAX2), DPROBI(IAMAX2),
     +   P1AEL(0:IAMAX),DP1AEL(0:IAMAX),P1BEL(0:IAMAX), DP1BEL(0:IAMAX),
     +   P2AEL(0:IAMAX),DP2AEL(0:IAMAX),P2BEL(0:IAMAX), DP2BEL(0:IAMAX)
      COMMON /CNUCMS/ B, BMAX, NTRY, NA, NB, NI, NAEL, NBEL
     +         ,JJA(IAMAX), JJB(IAMAX), JJINT(IAMAX,IAMAX)
     +         ,JJAEL(IAMAX), JJBEL(IAMAX)
      DIMENSION  MMA(0:IAMAX), MMB(0:IAMAX), MMI(0:IAMAX2)
      DIMENSION  M1AEL(0:IAMAX), M1BEL(0:IAMAX)
      DIMENSION  M2AEL(0:IAMAX), M2BEL(0:IAMAX)
      DOUBLE PRECISION FOX
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN
      SAVE
      DATA FOX /0.21522D0/  !atomic percentage of 'non-nitrogen' in air

      R2 = 0.1D0 * SIG0/PI
      BMAX = 15.D0                           ! fm
      SIGMA0 = PI*BMAX*BMAX*10.              ! mbarn
      IA = 16
      DO J=1,IA
         MMA(J) = 0
         M1AEL(J) = 0
         M2AEL(J) = 0
      ENDDO
      DO J=1,IB
         MMB(J) = 0
         M1BEL(J) = 0
         M2BEL(J) = 0
      ENDDO
      DO J=1,IA*IB
         MMI(J) = 0
      ENDDO
      NN = 0
      M = 0
      DO KK=1,KINT
c  select target IA from air composition
         R = S_RNDM(KK)
         IA = 14
         IF (R .LT. FOX)  IA = 16

         CALL INT_NUC (IA, IB, SIG0, SIGEL) 
         NN = NN + NTRY
         MMI(NI) = MMI(NI) + 1
         MMA(NA) = MMA(NA)+1
         MMB(NB) = MMB(NB)+1
         IF (NI .GT. 0)  THEN
            M = M+1
            M1AEL(NAEL) = M1AEL(NAEL)+1
            M1BEL(NBEL) = M1BEL(NBEL)+1
         ELSE
            M2AEL(NAEL) = M2AEL(NAEL)+1
            M2BEL(NBEL) = M2BEL(NBEL)+1
         ENDIF
      ENDDO
      MQE = KINT - M
      SIGMA  = SIGMA0 * DBLE(M)/DBLE(NN)
      DSIGMA = SIGMA0 * SQRT(DBLE(M))/DBLE(NN)
      SIGQE  = SIGMA0 * DBLE(MQE)/DBLE(NN)
      DSIGQE = SIGMA0 * SQRT(DBLE(MQE))/DBLE(NN)
      DO J=1,IA
         PROBA(J) = DBLE(MMA(J))/DBLE(M)
         DPROBA(J) = SQRT(DBLE(MMA(J)))/DBLE(M)
      ENDDO
      DO J=1,IB
         PROBB(J) = DBLE(MMB(J))/DBLE(M)
         DPROBB(J) = SQRT(DBLE(MMB(J)))/DBLE(M)
      ENDDO
      DO J=1,IA*IB
         PROBI(J) = DBLE(MMI(J))/DBLE(M)
         DPROBI(J) = SQRT(DBLE(MMI(J)))/DBLE(M)
      ENDDO
      DO J=0,IA
         P1AEL(J) = DBLE(M1AEL(J))/DBLE(M)
         DP1AEL(J) = SQRT(DBLE(M1AEL(J)))/DBLE(M)
         P2AEL(J) = DBLE(M2AEL(J))/DBLE(MQE)
         DP2AEL(J) = SQRT(DBLE(M2AEL(J)))/DBLE(MQE)
      ENDDO
      DO J=0,IB
         P1BEL(J) = DBLE(M1BEL(J))/DBLE(M)
         DP1BEL(J) = SQRT(DBLE(M1BEL(J)))/DBLE(M)
         P2BEL(J) = DBLE(M2BEL(J))/DBLE(MQE)
         DP2BEL(J) = SQRT(DBLE(M2BEL(J)))/DBLE(MQE)
      ENDDO
      RETURN
      END
C->
C=======================================================================

      SUBROUTINE SIGMA_MC (IA,IB,SIG0,SIGEL,KINT,
     +                            SIGMA,DSIGMA,SIGQE,DSIGQE)

C-----------------------------------------------------------------------
C...Compute with a montecarlo method the "production"
C.  and "quasi-elastic" cross section for  
C.  a nucleus-nucleus interaction
C.
C.  INPUT : IA            = mass of target nucleus
C.          IB            = mass of projectile nucleus
C.          SIG0 (mbarn)  = inelastic pp cross section
C.          KINT            = number  of interactions to generate
C.  OUTPUT : SIGMA (mbarn) = "production" cross section
C.           DSIGMA   "    = error
C.           SIGQE    "    = "quasi-elastic" cross section
C.           DSIGQE   "    = error
C.           additional output is in the common block  /CPROBAB/
C.           Prob(n_A), Prob(n_B), Prob(n_int)
C..........................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      PARAMETER (IAMAX=56)
      PARAMETER (IAMAX2=3136)          ! IAMAX*IAMAX
      COMMON  /CPROBAB/ PROBA(IAMAX), DPROBA(IAMAX), 
     +   PROBB(IAMAX), DPROBB(IAMAX), PROBI(IAMAX2), DPROBI(IAMAX2),
     +   P1AEL(0:IAMAX),DP1AEL(0:IAMAX),P1BEL(0:IAMAX), DP1BEL(0:IAMAX),
     +   P2AEL(0:IAMAX),DP2AEL(0:IAMAX),P2BEL(0:IAMAX), DP2BEL(0:IAMAX)
      COMMON /CNUCMS/ B, BMAX, NTRY, NA, NB, NI, NAEL, NBEL
     +         ,JJA(IAMAX), JJB(IAMAX), JJINT(IAMAX,IAMAX)
     +         ,JJAEL(IAMAX), JJBEL(IAMAX)
      DIMENSION  MMA(0:IAMAX), MMB(0:IAMAX), MMI(0:IAMAX2)
      DIMENSION  M1AEL(0:IAMAX), M1BEL(0:IAMAX)
      DIMENSION  M2AEL(0:IAMAX), M2BEL(0:IAMAX)
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN
      SAVE

      R2 = 0.1D0 * SIG0/PI
      BMAX = 15.D0                           ! fm
      SIGMA0 = PI*BMAX*BMAX*10.D0              ! mbarn
      DO J=1,IA
         MMA(J) = 0
         M1AEL(J) = 0
         M2AEL(J) = 0
      ENDDO
      DO J=1,IB
         MMB(J) = 0
         M1BEL(J) = 0
         M2BEL(J) = 0
      ENDDO
      DO J=1,IA*IB
         MMI(J) = 0
      ENDDO
      NN = 0
      M = 0
      DO KK=1,KINT
         CALL INT_NUC (IA, IB, SIG0, SIGEL) 
         NN = NN + NTRY
         MMI(NI) = MMI(NI) + 1
         MMA(NA) = MMA(NA)+1
         MMB(NB) = MMB(NB)+1
         IF (NI .GT. 0)  THEN
            M = M+1
            M1AEL(NAEL) = M1AEL(NAEL)+1
            M1BEL(NBEL) = M1BEL(NBEL)+1
         ELSE
            M2AEL(NAEL) = M2AEL(NAEL)+1
            M2BEL(NBEL) = M2BEL(NBEL)+1
         ENDIF
      ENDDO
      MQE = KINT - M
      SIGMA  = SIGMA0 * DBLE(M)/DBLE(NN)
      DSIGMA = SIGMA0 * SQRT(DBLE(M))/DBLE(NN)
      SIGQE  = SIGMA0 * DBLE(MQE)/DBLE(NN)
      DSIGQE = SIGMA0 * SQRT(DBLE(MQE))/DBLE(NN)
      DO J=1,IA
         PROBA(J) = DBLE(MMA(J))/DBLE(M)
         DPROBA(J) = SQRT(DBLE(MMA(J)))/DBLE(M)
      ENDDO
      DO J=1,IB
         PROBB(J) = DBLE(MMB(J))/DBLE(M)
         DPROBB(J) = SQRT(DBLE(MMB(J)))/DBLE(M)
      ENDDO
      DO J=1,IA*IB
         PROBI(J) = DBLE(MMI(J))/DBLE(M)
         DPROBI(J) = SQRT(DBLE(MMI(J)))/DBLE(M)
      ENDDO
      DO J=0,IA
         P1AEL(J) = DBLE(M1AEL(J))/DBLE(M)
         DP1AEL(J) = SQRT(DBLE(M1AEL(J)))/DBLE(M)
         P2AEL(J) = DBLE(M2AEL(J))/DBLE(MQE)
         DP2AEL(J) = SQRT(DBLE(M2AEL(J)))/DBLE(MQE)
      ENDDO
      DO J=0,IB
         P1BEL(J) = DBLE(M1BEL(J))/DBLE(M)
         DP1BEL(J) = SQRT(DBLE(M1BEL(J)))/DBLE(M)
         P2BEL(J) = DBLE(M2BEL(J))/DBLE(MQE)
         DP2BEL(J) = SQRT(DBLE(M2BEL(J)))/DBLE(MQE)
      ENDDO
      RETURN
      END

C*=============================================================
C.  Cross sections
C*=============================================================

C Glauber h-air cross section calculation moved to inelScreen src file..

C-----------------------------------------------------------------------
C.  Fit of Block and Cahn to pp and pbar-p cross sections
C-----------------------------------------------------------------------
C=======================================================================

      SUBROUTINE BLOCK(SQS,SIG1,SIG2,SLOP1,SLOP2,
     +                 RHO1,RHO2,SIGEL1,SIGEL2)

C-----------------------------------------------------------------------
C...p-p and pbar-p cross sections
C.  Parametrization of  Block and Cahn
C
C.  INPUT  : SQS   (GeV)  = c.m. energy
C.  
C.  OUPUT : SIG1 (mbarn)    = pp  total  cross section 
C.          SLOP1 (GeV**2)  = slope of elastic scattering
C.          RHO1            = Real/Imaginary part of the amplitude
C.                            for forward elastic  scattering (pp)
C.          SIGEL1 (mbarn)  = pp  elastic scattering  cross section
C.          [1 -> 2   : pp -> pbar p]
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN
      SAVE

      S = SQS*SQS
      CALL FPLUS  (S, FR, FI)
      CALL FMINUS (S, GR, GI)
      SIG1 = FI-GI
      SIG2 = FI+GI
      RHO1 = (FR-GR)/(FI-GI)
      RHO2 = (FR+GR)/(FI+GI)
      CALL SSLOPE (S, BP, BM)
      SLOP1 = BP - GI/FI*(BM-BP)
      SLOP2 = BP + GI/FI*(BM-BP)
      SIGEL1 = SIG1**2*(1.D0+RHO1**2)/(16.D0*PI*SLOP1)/CMBARN
      SIGEL2 = SIG2**2*(1.D0+RHO2**2)/(16.D0*PI*SLOP2)/CMBARN
      RETURN
      END
C=======================================================================

      SUBROUTINE FPLUS (S, FR, FI)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON /BLOCKC/ AA, BETA, S0, CC, AMU, DD, ALPHA, A0
      COMPLEX*16 Z1, Z2, Z3
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN
      SAVE

      F1 = LOG(S/S0)
      Z1 = DCMPLX(F1,-PI/2.D0)
      Z1 = Z1*Z1
      Z2 = 1.D0 + A0*Z1
      Z3 = Z1/Z2
      F2 = CC*S**(AMU-1.D0)
      F3 = 0.5D0*PI*(1.-AMU)
      FI = AA + F2*COS(F3) + BETA*DREAL(Z3)
      FR = -BETA*DIMAG(Z3)+F2*SIN(F3)
      RETURN
      END
C=======================================================================

      SUBROUTINE FMINUS (S, FR, FI)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON /BLOCKC/ AA, BETA, S0, CC, AMU, DD, ALPHA, A0
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN
      SAVE

      F1 = S**(ALPHA-1.D0)
      F2 = 0.5D0*PI*(1.D0-ALPHA)
      FR = -DD*F1*COS(F2)
      FI = -DD*F1*SIN(F2)
      RETURN
      END
C=======================================================================

      SUBROUTINE SSLOPE (S, BP, BM)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON /BLOCKD/ CP, DP, EP, CM, DM
      SAVE

      AL = LOG(S)
      BP = CP + DP*AL + EP*AL*AL
      BM = CM + DM*AL
      RETURN
      END
C=======================================================================

      SUBROUTINE BLOCK_INI

C-----------------------------------------------------------------------
C...Parameters of fit IFIT=1 of Block and Cahn
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON /BLOCKC/ AA, BETA, S0, CC, AMU, DD, ALPHA, A0
      COMMON /BLOCKD/ CP, DP, EP, CM, DM
      SAVE

      AA = 41.74D0
      BETA = 0.66D0
      S0 = 338.5D0
      CC = 0.D0
      AMU = 0.D0
      DD = -39.37D0
      ALPHA = 0.48D0
      A0 = 0.D0
      CP = 10.90D0
      DP = -0.08D0
      EP = 0.043D0
      CM = 23.27D0
      DM = 0.93D0
      RETURN
      END

C*=============================================================
C.  Nucleus-nucleus cross sections
C=======================================================================

      SUBROUTINE SIGNUC_INI (IA,E0)

C-----------------------------------------------------------------------
C...This subroutine receives in INPUT E0 (TeV)
C.  energy per nucleon and computes the cross sections
C.  and interactions lengths for  all nuclei
C.  with A  between 2 and IA
C.  The output is contained in common block /CLENNN/
C.
C.  Attention: the tabulated cross sections are obtained with
C.  new p-p cross sections as used in SIBYLL 2x,
C.  in addition field dimensions changed (RE 04/2000)
C.
C........................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON /CLENNN/ SSIGNUC(60), ALNUC(60)
      DIMENSION SIGMA(6,56), SIGQE(6,56)
      DIMENSION AA(6)
      SAVE
      DATA NE /6/, AMIN /1.D0/, DA /1.D0/
      DATA AA /1.D0,2.D0,3.D0,4.D0,5.D0,6.D0/
      DATA AVOG /6.0221367D-04/
      DATA ATARGET /14.514D0/            ! effective masss of air
C...Data on `inelastic-production' nucleus-air cross section
      DATA (SIGMA(J, 2),J=1,6) /
     &3.842D+02,4.287D+02,4.940D+02,5.887D+02,6.922D+02,7.767D+02/
      DATA (SIGMA(J, 3),J=1,6) /
     &4.601D+02,5.149D+02,5.595D+02,6.663D+02,7.641D+02,8.446D+02/
      DATA (SIGMA(J, 4),J=1,6) /
     &4.881D+02,5.373D+02,6.005D+02,6.895D+02,7.716D+02,8.967D+02/
      DATA (SIGMA(J, 5),J=1,6) /
     &5.874D+02,6.176D+02,7.181D+02,7.993D+02,9.089D+02,1.031D+03/
      DATA (SIGMA(J, 6),J=1,6) /
     &7.054D+02,7.399D+02,8.388D+02,9.463D+02,1.080D+03,1.197D+03/
      DATA (SIGMA(J, 7),J=1,6) /
     &7.192D+02,7.611D+02,8.449D+02,9.539D+02,1.061D+03,1.176D+03/
      DATA (SIGMA(J, 8),J=1,6) /
     &7.550D+02,7.975D+02,9.153D+02,9.944D+02,1.126D+03,1.236D+03/
      DATA (SIGMA(J, 9),J=1,6) /
     &7.929D+02,8.392D+02,9.265D+02,1.059D+03,1.167D+03,1.262D+03/
      DATA (SIGMA(J, 10),J=1,6) /
     &8.157D+02,8.644D+02,9.512D+02,1.058D+03,1.182D+03,1.298D+03/
      DATA (SIGMA(J, 11),J=1,6) /
     &8.039D+02,8.587D+02,9.534D+02,1.055D+03,1.182D+03,1.298D+03/
      DATA (SIGMA(J, 12),J=1,6) /
     &8.515D+02,8.957D+02,9.869D+02,1.122D+03,1.253D+03,1.366D+03/
      DATA (SIGMA(J, 13),J=1,6) /
     &8.769D+02,9.100D+02,1.018D+03,1.119D+03,1.252D+03,1.341D+03/
      DATA (SIGMA(J, 14),J=1,6) /
     &9.058D+02,9.532D+02,1.057D+03,1.171D+03,1.302D+03,1.391D+03/
      DATA (SIGMA(J, 15),J=1,6) /
     &9.555D+02,9.799D+02,1.098D+03,1.201D+03,1.342D+03,1.444D+03/
      DATA (SIGMA(J, 16),J=1,6) /
     &1.009D+03,1.058D+03,1.149D+03,1.290D+03,1.414D+03,1.520D+03/
      DATA (SIGMA(J, 17),J=1,6) /
     &9.907D+02,1.045D+03,1.166D+03,1.290D+03,1.384D+03,1.516D+03/
      DATA (SIGMA(J, 18),J=1,6) /
     &1.036D+03,1.121D+03,1.198D+03,1.328D+03,1.470D+03,1.592D+03/
      DATA (SIGMA(J, 19),J=1,6) /
     &1.083D+03,1.162D+03,1.250D+03,1.371D+03,1.516D+03,1.661D+03/
      DATA (SIGMA(J, 20),J=1,6) /
     &1.146D+03,1.215D+03,1.295D+03,1.443D+03,1.544D+03,1.744D+03/
      DATA (SIGMA(J, 21),J=1,6) /
     &1.158D+03,1.234D+03,1.292D+03,1.467D+03,1.618D+03,1.750D+03/
      DATA (SIGMA(J, 22),J=1,6) /
     &1.153D+03,1.205D+03,1.329D+03,1.451D+03,1.596D+03,1.734D+03/
      DATA (SIGMA(J, 23),J=1,6) /
     &1.210D+03,1.274D+03,1.356D+03,1.493D+03,1.655D+03,1.803D+03/
      DATA (SIGMA(J, 24),J=1,6) /
     &1.212D+03,1.273D+03,1.398D+03,1.489D+03,1.641D+03,1.800D+03/
      DATA (SIGMA(J, 25),J=1,6) /
     &1.236D+03,1.315D+03,1.423D+03,1.561D+03,1.669D+03,1.855D+03/
      DATA (SIGMA(J, 26),J=1,6) /
     &1.279D+03,1.345D+03,1.431D+03,1.595D+03,1.734D+03,1.889D+03/
      DATA (SIGMA(J, 27),J=1,6) /
     &1.228D+03,1.304D+03,1.438D+03,1.546D+03,1.714D+03,1.836D+03/
      DATA (SIGMA(J, 28),J=1,6) /
     &1.289D+03,1.370D+03,1.451D+03,1.597D+03,1.754D+03,1.913D+03/
      DATA (SIGMA(J, 29),J=1,6) /
     &1.411D+03,1.469D+03,1.613D+03,1.777D+03,1.910D+03,2.075D+03/
      DATA (SIGMA(J, 30),J=1,6) /
     &1.347D+03,1.401D+03,1.498D+03,1.642D+03,1.816D+03,1.975D+03/
      DATA (SIGMA(J, 31),J=1,6) /
     &1.359D+03,1.448D+03,1.551D+03,1.694D+03,1.858D+03,2.007D+03/
      DATA (SIGMA(J, 32),J=1,6) /
     &1.358D+03,1.460D+03,1.559D+03,1.698D+03,1.842D+03,1.974D+03/
      DATA (SIGMA(J, 33),J=1,6) /
     &1.418D+03,1.448D+03,1.578D+03,1.727D+03,1.872D+03,2.047D+03/
      DATA (SIGMA(J, 34),J=1,6) /
     &1.433D+03,1.466D+03,1.605D+03,1.738D+03,1.892D+03,2.019D+03/
      DATA (SIGMA(J, 35),J=1,6) /
     &1.430D+03,1.511D+03,1.602D+03,1.752D+03,1.935D+03,2.060D+03/
      DATA (SIGMA(J, 36),J=1,6) /
     &1.462D+03,1.499D+03,1.653D+03,1.805D+03,1.920D+03,2.057D+03/
      DATA (SIGMA(J, 37),J=1,6) /
     &1.470D+03,1.520D+03,1.656D+03,1.818D+03,1.946D+03,2.131D+03/
      DATA (SIGMA(J, 38),J=1,6) /
     &1.470D+03,1.542D+03,1.691D+03,1.800D+03,1.968D+03,2.133D+03/
      DATA (SIGMA(J, 39),J=1,6) /
     &1.495D+03,1.588D+03,1.676D+03,1.834D+03,1.969D+03,2.163D+03/
      DATA (SIGMA(J, 40),J=1,6) /
     &1.525D+03,1.551D+03,1.722D+03,1.833D+03,2.020D+03,2.192D+03/
      DATA (SIGMA(J, 41),J=1,6) /
     &1.526D+03,1.615D+03,1.709D+03,1.899D+03,2.040D+03,2.181D+03/
      DATA (SIGMA(J, 42),J=1,6) /
     &1.510D+03,1.567D+03,1.716D+03,1.892D+03,2.056D+03,2.197D+03/
      DATA (SIGMA(J, 43),J=1,6) /
     &1.557D+03,1.658D+03,1.776D+03,1.898D+03,2.092D+03,2.200D+03/
      DATA (SIGMA(J, 44),J=1,6) /
     &1.556D+03,1.645D+03,1.752D+03,1.920D+03,2.091D+03,2.243D+03/
      DATA (SIGMA(J, 45),J=1,6) /
     &1.583D+03,1.663D+03,1.798D+03,1.940D+03,2.051D+03,2.263D+03/
      DATA (SIGMA(J, 46),J=1,6) /
     &1.599D+03,1.642D+03,1.799D+03,1.941D+03,2.107D+03,2.268D+03/
      DATA (SIGMA(J, 47),J=1,6) /
     &1.611D+03,1.692D+03,1.811D+03,1.956D+03,2.107D+03,2.264D+03/
      DATA (SIGMA(J, 48),J=1,6) /
     &1.625D+03,1.706D+03,1.819D+03,1.986D+03,2.139D+03,2.354D+03/
      DATA (SIGMA(J, 49),J=1,6) /
     &1.666D+03,1.737D+03,1.854D+03,1.971D+03,2.160D+03,2.318D+03/
      DATA (SIGMA(J, 50),J=1,6) /
     &1.648D+03,1.747D+03,1.856D+03,2.023D+03,2.181D+03,2.352D+03/
      DATA (SIGMA(J, 51),J=1,6) /
     &1.653D+03,1.763D+03,1.868D+03,2.015D+03,2.203D+03,2.386D+03/
      DATA (SIGMA(J, 52),J=1,6) /
     &1.690D+03,1.720D+03,1.902D+03,2.027D+03,2.189D+03,2.357D+03/
      DATA (SIGMA(J, 53),J=1,6) /
     &1.690D+03,1.750D+03,1.921D+03,2.059D+03,2.208D+03,2.417D+03/
      DATA (SIGMA(J, 54),J=1,6) /
     &1.705D+03,1.781D+03,1.911D+03,2.073D+03,2.242D+03,2.411D+03/
      DATA (SIGMA(J, 55),J=1,6) /
     &1.714D+03,1.806D+03,1.896D+03,2.100D+03,2.253D+03,2.411D+03/
      DATA (SIGMA(J, 56),J=1,6) /
     &1.774D+03,1.813D+03,1.954D+03,2.098D+03,2.280D+03,2.482D+03/
 
      DATA (SIGQE(J, 2),J=1,6) /
     &4.141D+01,3.708D+01,5.428D+01,8.696D+01,1.403D+02,1.885D+02/
      DATA (SIGQE(J, 3),J=1,6) /
     &4.357D+01,3.894D+01,5.177D+01,9.675D+01,1.447D+02,2.029D+02/
      DATA (SIGQE(J, 4),J=1,6) /
     &4.123D+01,3.933D+01,6.070D+01,9.482D+01,1.474D+02,2.023D+02/
      DATA (SIGQE(J, 5),J=1,6) /
     &4.681D+01,4.287D+01,6.381D+01,1.050D+02,1.519D+02,2.198D+02/
      DATA (SIGQE(J, 6),J=1,6) /
     &5.407D+01,5.195D+01,6.723D+01,1.108D+02,1.750D+02,2.368D+02/
      DATA (SIGQE(J, 7),J=1,6) /
     &4.975D+01,4.936D+01,6.880D+01,1.162D+02,1.689D+02,2.329D+02/
      DATA (SIGQE(J, 8),J=1,6) /
     &5.361D+01,5.027D+01,6.858D+01,1.177D+02,1.759D+02,2.412D+02/
      DATA (SIGQE(J, 9),J=1,6) /
     &4.980D+01,5.063D+01,7.210D+01,1.196D+02,1.806D+02,2.299D+02/
      DATA (SIGQE(J, 10),J=1,6) /
     &5.170D+01,5.070D+01,7.105D+01,1.182D+02,1.679D+02,2.411D+02/
      DATA (SIGQE(J, 11),J=1,6) /
     &4.950D+01,4.950D+01,7.286D+01,1.137D+02,1.769D+02,2.477D+02/
      DATA (SIGQE(J, 12),J=1,6) /
     &5.262D+01,5.133D+01,7.110D+01,1.204D+02,1.789D+02,2.501D+02/
      DATA (SIGQE(J, 13),J=1,6) /
     &5.320D+01,5.378D+01,6.847D+01,1.200D+02,1.805D+02,2.442D+02/
      DATA (SIGQE(J, 14),J=1,6) /
     &5.638D+01,5.271D+01,6.985D+01,1.209D+02,1.867D+02,2.610D+02/
      DATA (SIGQE(J, 15),J=1,6) /
     &5.294D+01,5.353D+01,7.435D+01,1.211D+02,1.899D+02,2.612D+02/
      DATA (SIGQE(J, 16),J=1,6) /
     &5.668D+01,5.254D+01,7.557D+01,1.269D+02,1.917D+02,2.707D+02/
      DATA (SIGQE(J, 17),J=1,6) /
     &5.456D+01,5.721D+01,7.481D+01,1.208D+02,1.859D+02,2.658D+02/
      DATA (SIGQE(J, 18),J=1,6) /
     &5.901D+01,5.382D+01,7.591D+01,1.246D+02,1.872D+02,2.874D+02/
      DATA (SIGQE(J, 19),J=1,6) /
     &6.328D+01,6.116D+01,8.451D+01,1.318D+02,2.088D+02,2.749D+02/
      DATA (SIGQE(J, 20),J=1,6) /
     &5.779D+01,5.924D+01,8.382D+01,1.370D+02,2.062D+02,2.837D+02/
      DATA (SIGQE(J, 21),J=1,6) /
     &7.155D+01,5.732D+01,8.231D+01,1.363D+02,2.047D+02,2.820D+02/
      DATA (SIGQE(J, 22),J=1,6) /
     &6.699D+01,5.651D+01,8.511D+01,1.477D+02,2.031D+02,2.921D+02/
      DATA (SIGQE(J, 23),J=1,6) /
     &6.179D+01,6.269D+01,9.395D+01,1.437D+02,2.195D+02,2.964D+02/
      DATA (SIGQE(J, 24),J=1,6) /
     &6.784D+01,6.028D+01,8.622D+01,1.279D+02,2.214D+02,2.867D+02/
      DATA (SIGQE(J, 25),J=1,6) /
     &6.589D+01,5.795D+01,8.890D+01,1.385D+02,2.055D+02,2.988D+02/
      DATA (SIGQE(J, 26),J=1,6) /
     &6.364D+01,6.325D+01,8.942D+01,1.421D+02,2.128D+02,3.083D+02/
      DATA (SIGQE(J, 27),J=1,6) /
     &6.449D+01,6.664D+01,8.986D+01,1.453D+02,2.140D+02,2.932D+02/
      DATA (SIGQE(J, 28),J=1,6) /
     &7.284D+01,6.139D+01,8.867D+01,1.425D+02,2.179D+02,2.978D+02/
      DATA (SIGQE(J, 29),J=1,6) /
     &7.221D+01,7.085D+01,9.079D+01,1.482D+02,2.277D+02,2.913D+02/
      DATA (SIGQE(J, 30),J=1,6) /
     &6.928D+01,6.294D+01,8.935D+01,1.463D+02,2.265D+02,2.834D+02/
      DATA (SIGQE(J, 31),J=1,6) /
     &6.611D+01,6.586D+01,9.133D+01,1.461D+02,2.201D+02,2.959D+02/
      DATA (SIGQE(J, 32),J=1,6) /
     &6.401D+01,6.177D+01,8.971D+01,1.480D+02,2.155D+02,3.152D+02/
      DATA (SIGQE(J, 33),J=1,6) /
     &7.057D+01,6.918D+01,8.410D+01,1.465D+02,2.288D+02,3.088D+02/
      DATA (SIGQE(J, 34),J=1,6) /
     &6.453D+01,7.020D+01,9.272D+01,1.517D+02,2.189D+02,2.999D+02/
      DATA (SIGQE(J, 35),J=1,6) /
     &6.741D+01,6.295D+01,9.323D+01,1.536D+02,2.190D+02,2.930D+02/
      DATA (SIGQE(J, 36),J=1,6) /
     &6.807D+01,7.046D+01,1.025D+02,1.565D+02,2.315D+02,3.090D+02/
      DATA (SIGQE(J, 37),J=1,6) /
     &8.082D+01,6.565D+01,9.160D+01,1.572D+02,2.229D+02,3.125D+02/
      DATA (SIGQE(J, 38),J=1,6) /
     &6.494D+01,6.964D+01,9.089D+01,1.653D+02,2.336D+02,3.120D+02/
      DATA (SIGQE(J, 39),J=1,6) /
     &6.833D+01,6.860D+01,8.933D+01,1.601D+02,2.261D+02,3.167D+02/
      DATA (SIGQE(J, 40),J=1,6) /
     &7.021D+01,6.866D+01,8.437D+01,1.588D+02,2.249D+02,2.941D+02/
      DATA (SIGQE(J, 41),J=1,6) /
     &7.122D+01,6.205D+01,9.545D+01,1.582D+02,2.335D+02,3.395D+02/
      DATA (SIGQE(J, 42),J=1,6) /
     &7.265D+01,6.936D+01,9.486D+01,1.505D+02,2.379D+02,3.248D+02/
      DATA (SIGQE(J, 43),J=1,6) /
     &7.048D+01,7.539D+01,9.192D+01,1.566D+02,2.532D+02,3.182D+02/
      DATA (SIGQE(J, 44),J=1,6) /
     &6.650D+01,7.139D+01,9.862D+01,1.602D+02,2.289D+02,3.077D+02/
      DATA (SIGQE(J, 45),J=1,6) /
     &7.511D+01,6.893D+01,9.245D+01,1.641D+02,2.519D+02,3.381D+02/
      DATA (SIGQE(J, 46),J=1,6) /
     &6.437D+01,6.894D+01,8.697D+01,1.544D+02,2.391D+02,3.213D+02/
      DATA (SIGQE(J, 47),J=1,6) /
     &7.980D+01,6.958D+01,1.022D+02,1.609D+02,2.408D+02,3.246D+02/
      DATA (SIGQE(J, 48),J=1,6) /
     &7.265D+01,7.313D+01,8.989D+01,1.578D+02,2.387D+02,3.235D+02/
      DATA (SIGQE(J, 49),J=1,6) /
     &6.959D+01,6.337D+01,9.084D+01,1.656D+02,2.331D+02,3.226D+02/
      DATA (SIGQE(J, 50),J=1,6) /
     &7.371D+01,6.807D+01,9.726D+01,1.535D+02,2.445D+02,3.189D+02/
      DATA (SIGQE(J, 51),J=1,6) /
     &7.882D+01,6.680D+01,9.377D+01,1.629D+02,2.448D+02,3.297D+02/
      DATA (SIGQE(J, 52),J=1,6) /
     &7.223D+01,6.794D+01,9.925D+01,1.738D+02,2.446D+02,3.162D+02/
      DATA (SIGQE(J, 53),J=1,6) /
     &7.703D+01,6.971D+01,9.601D+01,1.595D+02,2.484D+02,3.265D+02/
      DATA (SIGQE(J, 54),J=1,6) /
     &7.549D+01,7.459D+01,8.984D+01,1.645D+02,2.348D+02,3.201D+02/
      DATA (SIGQE(J, 55),J=1,6) /
     &7.891D+01,6.840D+01,1.017D+02,1.698D+02,2.501D+02,3.429D+02/
      DATA (SIGQE(J, 56),J=1,6) /
     &7.545D+01,6.673D+01,1.057D+02,1.684D+02,2.424D+02,3.181D+02/

      ASQS = 0.5D0*LOG10(1.876D+03*E0)
      JE = MIN(INT((ASQS-AMIN)/DA)+1,NE-2)
      DO JA=2,IA
         ABEAM = DBLE(JA)
         S1 = QUAD_INT(ASQS, AA(JE),AA(JE+1),AA(JE+2),
     +                   SIGMA(JE,JA),SIGMA(JE+1,JA),SIGMA(JE+2,JA))
         S2 = QUAD_INT(ASQS, AA(JE),AA(JE+1),AA(JE+2),
     +                   SIGQE(JE,JA),SIGQE(JE+1,JA),SIGQE(JE+2,JA))
         SSIGNUC(JA) = S1 + S2
         ALNUC(JA) = ATARGET/(AVOG*SSIGNUC(JA))
      ENDDO
      ALNUC(1) = FPNI(E0, 13)
      SSIGNUC(1) = ATARGET/(AVOG*ALNUC(1))

      RETURN
      END


C*=======================================================================
C.  General utilities
C=======================================================================

      FUNCTION QUAD_INT (R,X0,X1,X2,V0,V1,V2)

C-----------------------------------------------------------------------
C...Quadratic interpolation
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      SAVE

      R0=R-X0
      R1=R-X1
      R2=R-X2
      S0=X0-X1
      S1=X0-X2
      S2=X1-X2
      QUAD_INT = V0*R1*R2/(S0*S1)-V1*R0*R2/(S0*S2)+V2*R0*R1/(S1*S2)
      RETURN
      END
C=======================================================================

      FUNCTION GAUSS (FUN, A,B)

C-----------------------------------------------------------------------
C...Returns the  8 points Gauss-Legendre integral
C.  of function FUN from A to B
C...........................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      DIMENSION X(8), W(8)
      SAVE
      DATA X/.0950125098D0, .2816035507D0, .4580167776D0, .6178762444D0,
     1       .7554044083D0, .8656312023D0, .9445750230D0, .9894009349D0/
      DATA W/.1894506104D0, .1826034150D0, .1691565193D0, .1495959888D0,
     1       .1246289712D0, .0951585116D0, .0622535239D0, .0271524594D0/

      XM = 0.5D0*(B+A)
      XR = 0.5D0*(B-A)
      SS = 0.D0
      DO J=1,8
        DX = XR*X(J)
        SS = SS + W(J) * (FUN(XM+DX) + FUN(XM-DX))
      ENDDO
      GAUSS = XR*SS
      RETURN
      END
C=======================================================================

      SUBROUTINE INVERT_ARRAY (yy, xmin, dx, n, xnew, ymin, dy)

C-----------------------------------------------------------------------
C..    This subroutine receives one   array
C      of n y values in input yy(1:n)
C      that correspond to  equispaced values of x_j = xmin + dx*(j-1)
C
C      and "reverse" the array returning an array of  x values
C      xnew (1:n) that  corresponds to equispaced values of y
C      The relation is assumed monotonous but can be 
C      increasing or decreasing
C..............................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      dimension  yy(n), xnew (n)
      SAVE

      ymin = yy(1)
      ymax = yy(n)
      dy = (ymax - ymin)/float(n-1)
      xnew (1) = xmin
      xnew (n) = xmin + dx*float(n-1)
      k0 = 1
      do j=2,n-1
         y = ymin + float(j-1)*dy 
         do k=k0,n
            if((yy(k) .gt. y) .eqv. (yy(n) .gt. yy(1))) goto 100
         enddo
100      y2 = yy(k)
         y1 = yy(k-1)
         k0 = k-1
         x1 = xmin + dx*float(k-2)
         x2 = x1+dx
         xnew (j)  = x1 + dx* (y-y1)/(y2-y1)
      enddo
      return
      end
C->
C=======================================================================

      SUBROUTINE SINCO(S,C)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN
      SAVE

      F = TWOPI*S_RNDM(0)
      C = COS (F)
      S = SIN (F)
      RETURN
      END

C***********************************************************************
C.  Cross sections for cascade calculations (FPNI)
C=======================================================================
      
      SUBROUTINE SIGMA_PP (E0, SIGT, SIGEL, SIGINEL, SLOPE, RHO) 

C-----------------------------------------------------------------------
C...p-p cross sections
C.
C.  this routine serves the purpose to calculate cascades with different 
C.  cross sections
C.
C. INPUT: E0 = Laboratory Energy  (TeV)
C. 
C. OUTPUT: SIGT = total cross section
C.         SIGEL = elastic cross section
C.         SIGINEL = inelastic cross section
C.         SLOPE = slope of elastic scattering (GeV**-2)
C.         RHO = Imaginary/Real part of forward elastic amplitude
C.   
C.  (old cross section tables end at 10^6 GeV)
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      DIMENSION SSIG0(51)
      DIMENSION SIGDIF(3)
      COMMON /CSPA/ ICSPA2(3)
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN
      SAVE

C...p-p inelastic cross sections (mbarn)
      DATA (SSIG0(J),J=1,51) /
     +      32.05D0,  32.06D0,  32.08D0,  32.13D0,  32.22D0,  32.36D0,
     +      32.56D0,  32.85D0,  33.24D0,  33.75D0,  34.37D0,  35.14D0,
     +      36.05D0,  37.12D0,  38.37D0,  39.78D0,  41.36D0,  43.13D0,
     +      45.07D0,  47.18D0,  49.47D0,  51.91D0,  54.54D0,  57.28D0,
     +      60.15D0,  63.15D0,  66.28D0,  69.48D0,  72.80D0,  76.22D0,
     +      79.71D0,  83.27D0,  86.87D0,  90.55D0,  94.26D0,  98.05D0,
     +     101.89D0, 105.75D0, 109.71D0, 113.65D0, 117.60D0, 121.55D0,
     +     125.53D0, 129.56D0, 133.60D0, 137.70D0, 141.77D0, 145.84D0,
     +     149.92D0, 154.02D0, 158.15D0/

      ICSPA = ICSPA2(1)

      SQS = SQRT(2000.D0*0.938D0*E0)

*  pre-LHC SIBYLL2.1 model
      
      IF(ICSPA.EQ.-2) THEN

         CALL SIB_SIGMA_EXT(3,SQS,SIGT,SIGEL,SIGINEL,SLOPE,RHO)      

*  old standard NUCLIB/SIBYLL model

      ELSE IF(ICSPA.EQ.-1) THEN

        AL = LOG10(SQS)
        if(AL.le.1.D0) then
          SIGINEL = SSIG0(1)
        else
          J1 = INT((AL - 1.D0)*10.D0) + 1
          J1 = min(J1,50)
          T = (AL-1.D0)*10.D0 - DBLE(J1-1)
          SIGINEL = SSIG0(J1)*(1.D0-T) + SSIG0(J1+1)*T
        endif
        CALL BLOCK(SQS,SIGT1,SIGT2,SLOP1,SLOP2,RHO1,RHO2,SIGEL1,SIGEL2)
        R = SIGEL1/SIGT1
        RHO = RHO1
        SIGT  = SIGINEL/(1.D0-R)
        SIGEL = SIGINEL*R/(1.D0-R)
        SLOPE = SIGT**2/(SIGEL * 16.D0*PI) * (1.D0+RHO1**2) /CMBARN

*  cross section as calculated in SIBYLL

      ELSE IF(ICSPA.EQ.0) THEN

        CALL SIB_SIGMA_HP(1,SQS,SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO)

*  Donnachie-Landshoff  (sig-tot)

      ELSE IF(ICSPA.EQ.1) THEN

        CALL BLOCK(SQS,SIGT1,SIGT2,SLOP1,SLOP2,RHO1,RHO2,
     +             SIGEL1,SIGEL2)
        R = SIGEL1/SIGT1
        RHO = RHO1

        DELDL = 0.0808D0
        EPSDL = -0.4525D0
        S = SQS*SQS
        SIGT = 21.7D0*S**DELDL+56.08D0*S**EPSDL
        SIGEL = R*SIGT
        SIGINEL = SIGT-SIGEL
        SLOPE = SIGT**2/(SIGEL * 16.D0*PI) * (1.D0+RHO**2) /CMBARN

*  Donnachie-Landshoff (sig-tot and sig-el)

      ELSE IF(ICSPA.EQ.2) THEN

        DELDL = 0.0808D0
        EPSDL = -0.4525D0
        S = SQS*SQS
        SIGT = 21.7D0*S**DELDL+56.08D0*S**EPSDL
        IMODEL = 1
        IF(IMODEL.EQ.1) THEN
          ALPHAP = 0.25D0
          SLOPE = 8.5D0+2.D0*ALPHAP*LOG(S)
        ELSE IF(IMODEL.EQ.2) THEN
          ALPHAP = 0.3D0
          SLOPE = 8.D0+2.D0*ALPHAP*LOG(S)
        ENDIF
        SIGEL = SIGT**2/(16.D0*PI*SLOPE*CMBARN)
        SIGINEL = SIGT-SIGEL
        RHO = 0.D0

*  geometrical scaling with Donnachie-Landshoff sig-tot

      ELSE IF(ICSPA.EQ.3) THEN

        R = 0.17D0

        DELDL = 0.0808D0
        EPSDL = -0.4525D0
        S = SQS*SQS
        SIGT = 21.7D0*S**DELDL+56.08D0*S**EPSDL

        SIGEL = R*SIGT
        SIGINEL = SIGT-SIGEL
        SLOPE = SIGT**2/(16.D0*PI*SIGEL)/CMBARN
        RHO = 0.D0

c ICSPA=4 reserved for CONEX_EXTENSION
c      ELSE IF(ICSPA.EQ.4) THEN

*  cross section from 2014 Review of Particle Physics
        
      ELSE IF(ICSPA.EQ.5) THEN
         
c     elastic slope not included in fit
c     taking slope parameterization from sigma_pp Donnie.-Landshoff
         ALPHAP = 0.25D0
         SLOPE = 8.5D0+4.D0*ALPHAP*LOG(SQS)
         
         CALL SIG_RPP2014(1,1,SQS,SLOPE,SIGT,SIGEL,SIGINEL,RHO)

      ENDIF

      RETURN
      END

C=======================================================================

      SUBROUTINE SIGMA_PIP (E0, SIGT, SIGEL, SIGINEL, SLOPE, RHO) 

C-----------------------------------------------------------------------
C...pi-p cross sections
C.
C.  this routine serves the purpose to calculate cascades with different 
C.  cross sections
C.
C. INPUT: E0 = Laboratory Energy  (TeV)
C. 
C. OUTPUT: SIGT = total cross section
C.         SIGEL = elastic cross section
C.         SIGINEL = inelastic cross section
C.         SLOPE = slope of elastic scattering (GeV**-2)
C.         RHO = Imaginary/Real part of forward elastic amplitude
C.
C.  (old cross section tables end at 10^6 GeV)
C-----------------------------------------------------------------------
Cf2py double precision,intent(in) :: e0
Cf2py double precision,intent(out) :: sigt, sigel, siginel, slope, rho
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      DIMENSION SSIG0(51)
      DIMENSION SIGDIF(3)
      COMMON /CSPA/ ICSPA2(3)
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN
      SAVE

C...pi-p inelastic cross sections (mbarn)
      DATA (SSIG0(J),J=1,51) /
     +      20.76D0,  20.78D0,  20.81D0,  20.88D0,  20.98D0,  21.13D0,
     +      21.33D0,  21.61D0,  21.96D0,  22.39D0,  22.92D0,  23.56D0,
     +      24.31D0,  25.18D0,  26.18D0,  27.32D0,  28.60D0,  30.04D0,
     +      31.64D0,  33.40D0,  35.34D0,  37.43D0,  39.72D0,  42.16D0,
     +      44.77D0,  47.56D0,  50.53D0,  53.66D0,  56.99D0,  60.50D0,
     +      64.17D0,  68.03D0,  72.05D0,  76.27D0,  80.67D0,  85.27D0,
     +      90.08D0,  95.04D0, 100.27D0, 105.65D0, 111.21D0, 116.94D0,
     +     122.87D0, 129.03D0, 135.37D0, 141.93D0, 148.62D0, 155.49D0,
     +     162.48D0, 169.60D0, 176.94D0/

      ICSPA = ICSPA2(2)

      SQS = SQRT(2000.D0*0.938D0*E0)
      
*  pre-LHC SIBYLL2.1 model
      
      IF(ICSPA.EQ.-2) THEN

         CALL SIB_SIGMA_EXT(2,SQS,SIGT,SIGEL,SIGINEL,SLOPE,RHO)
      
*  old standard NUCLIB/SIBYLL model

      ELSE IF(ICSPA.EQ.-1) THEN

        AL = LOG10(SQS)
        if(AL.le.1.D0) then
          SIGINEL = SSIG0(1)
        else
          J1 = INT((AL - 1.D0)*10.D0) + 1
          J1 = min(J1,50)
          T = (AL-1.D0)*10.D0 - DBLE(J1-1)
          SIGINEL = SSIG0(J1)*(1.D0-T) + SSIG0(J1+1)*T
        endif
        CALL BLOCK(SQS,SIGT1,SIGT2,SLOP1,SLOP2,RHO1,RHO2,SIGEL1,SIGEL2)
        R = SIGEL1/SIGT1
        RHO = RHO1
        SIGT  = SIGINEL/(1.D0-R)
        SIGEL = SIGINEL*R/(1.D0-R)
        SLOPE = SIGT**2/(SIGEL * 16.D0*PI) * (1.D0+RHO1**2) /CMBARN

*  cross section as calculated in SIBYLL

      ELSE IF(ICSPA.EQ.0) THEN

        CALL SIB_SIGMA_HP(2,SQS,SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO)

*  Donnachie-Landshoff  (sig-tot)

      ELSE IF(ICSPA.EQ.1) THEN

        CALL BLOCK(SQS,SIGT1,SIGT2,SLOP1,SLOP2,RHO1,RHO2,
     +             SIGEL1,SIGEL2)
        R = SIGEL1/SIGT1
        RHO = RHO1

        DELDL = 0.0808D0
        EPSDL = -0.4525D0
        S = SQS*SQS
        SIGT = 13.63D0*S**DELDL+(36.02D0+27.56D0)/2.D0*S**EPSDL
        SIGEL = R*SIGT
        SIGINEL = SIGT-SIGEL
        SLOPE = SIGT**2/(SIGEL * 16.D0*PI) * (1.D0+RHO**2) /CMBARN

*  Donnachie-Landshoff (sig-tot and sig-el)

      ELSE IF(ICSPA.EQ.2) THEN

        DELDL = 0.0808D0
        EPSDL = -0.4525D0
        S = SQS*SQS
        SIGT = 13.63D0*S**DELDL+(36.02D0+27.56D0)/2.D0*S**EPSDL
        IMODEL = 1
        IF(IMODEL.EQ.1) THEN
          ALPHAP = 0.25D0
          SLOPE = 8.5D0+2.D0*ALPHAP*LOG(S)
        ELSE IF(IMODEL.EQ.2) THEN
          ALPHAP = 0.3D0
          SLOPE = 8.D0+2.D0*ALPHAP*LOG(S)
        ENDIF
        SIGEL = SIGT**2/(16.D0*PI*SLOPE*CMBARN)
        SIGINEL = SIGT-SIGEL
        RHO = 0.

*  geometrical scaling with Donnachie-Landshoff sig-tot

      ELSE IF(ICSPA.EQ.3) THEN

        R = 0.17D0

        DELDL = 0.0808D0
        EPSDL = -0.4525D0
        S = SQS*SQS
        SIGT = 13.63D0*S**DELDL+(36.02D0+27.56D0)/2.D0*S**EPSDL

        SIGEL = R*SIGT
        SIGINEL = SIGT-SIGEL
        SLOPE = SIGT**2/(16.D0*PI*SIGEL)/CMBARN
        RHO = 0.D0

c ICSPA=4 reserved for CONEX_EXTENSION
c      ELSE IF(ICSPA.EQ.4) THEN

*  cross section from 2014 Review of Particle Physics
        
      ELSE IF(ICSPA.EQ.5) THEN
         
c     elastic slope not included in fit
c     taking slope parameterization from sigma_pp Donnie.-Landshoff
         ALPHAP = 0.25D0
         SLOPE = 8.5D0+4.D0*ALPHAP*LOG(SQS)
         
         CALL SIG_RPP2014(2,1,SQS,SLOPE,SIGT,SIGEL,SIGINEL,RHO)
         
      ENDIF


      RETURN
      END

C=======================================================================

      SUBROUTINE SIGMA_KP (E0, SIGT, SIGEL, SIGINEL, SLOPE, RHO) 

C-----------------------------------------------------------------------
C...K-p cross sections
C.
C.  this routine serves the purpose to calculate cascades with different 
C.  cross sections
C.
C.  if old cross sections are selected then sigma_pi = sigma_K
C.
C. INPUT: E0 = Laboratory Energy  (TeV)
C. 
C. OUTPUT: SIGT = total cross section
C.         SIGEL = elastic cross section
C.         SIGINEL = inelastic cross section
C.         SLOPE = slope of elastic scattering (GeV**-2)
C.         RHO = Imaginary/Real part of forward elastic amplitude
C.
C.  (old cross section tables end at 10^6 GeV)
C-----------------------------------------------------------------------
Cf2py double precision,intent(in) :: e0
Cf2py double precision,intent(out) :: sigt, sigel, siginel, slope, rho
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      DIMENSION SSIG0(51)
      DIMENSION SIGDIF(3)
      COMMON /CSPA/ ICSPA2(3)
      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN
      SAVE

C...pi-p inelastic cross sections (mbarn)
      DATA (SSIG0(J),J=1,51) /
     +      20.76D0,  20.78D0,  20.81D0,  20.88D0,  20.98D0,  21.13D0,
     +      21.33D0,  21.61D0,  21.96D0,  22.39D0,  22.92D0,  23.56D0,
     +      24.31D0,  25.18D0,  26.18D0,  27.32D0,  28.60D0,  30.04D0,
     +      31.64D0,  33.40D0,  35.34D0,  37.43D0,  39.72D0,  42.16D0,
     +      44.77D0,  47.56D0,  50.53D0,  53.66D0,  56.99D0,  60.50D0,
     +      64.17D0,  68.03D0,  72.05D0,  76.27D0,  80.67D0,  85.27D0,
     +      90.08D0,  95.04D0, 100.27D0, 105.65D0, 111.21D0, 116.94D0,
     +     122.87D0, 129.03D0, 135.37D0, 141.93D0, 148.62D0, 155.49D0,
     +     162.48D0, 169.60D0, 176.94D0/

      ICSPA = ICSPA2(3)
      
      SQS = SQRT(2000.D0*0.938D0*E0)      

*  pre-LHC SIBYLL2.1 model
      
      IF(ICSPA.EQ.-2) THEN

         CALL SIB_SIGMA_EXT(3,SQS,SIGT,SIGEL,SIGINEL,SLOPE,RHO)
      
*  old standard NUCLIB/SIBYLL model

      ELSE IF(ICSPA.EQ.-1) THEN

        AL = LOG10(SQS)
        if(AL.le.1.D0) then
          SIGINEL = SSIG0(1)
        else
          J1 = INT((AL - 1.D0)*10.D0) + 1
          J1 = min(J1,50)
          T = (AL-1.D0)*10.D0 - DBLE(J1-1)
          SIGINEL = SSIG0(J1)*(1.D0-T) + SSIG0(J1+1)*T
        endif
        CALL BLOCK(SQS,SIGT1,SIGT2,SLOP1,SLOP2,RHO1,RHO2,SIGEL1,SIGEL2)
        R = SIGEL1/SIGT1
        RHO = RHO1
        SIGT  = SIGINEL/(1.D0-R)
        SIGEL = SIGINEL*R/(1.D0-R)
        SLOPE = SIGT**2/(SIGEL * 16.D0*PI) * (1.D0+RHO1**2) /CMBARN

*  cross section as calculated in SIBYLL

      ELSE IF(ICSPA.EQ.0) THEN

        CALL SIB_SIGMA_HP(3,SQS,SIGT,SIGEL,SIGINEL,SIGDIF,SLOPE,RHO)

*  Donnachie-Landshoff  (sig-tot)

      ELSE IF(ICSPA.EQ.1) THEN

        CALL BLOCK(SQS,SIGT1,SIGT2,SLOP1,SLOP2,RHO1,RHO2,
     +             SIGEL1,SIGEL2)
        R = SIGEL1/SIGT1
        RHO = RHO1

        DELDL = 0.0808D0
        EPSDL = -0.4525D0
        S = SQS*SQS
        SIGT = 11.82D0*S**DELDL+(26.36D0+ 8.15D0)/2.D0*S**EPSDL
        SIGEL = R*SIGT
        SIGINEL = SIGT-SIGEL
        SLOPE = SIGT**2/(SIGEL * 16.D0*PI) * (1.D0+RHO**2) /CMBARN

*  Donnachie-Landshoff (sig-tot and sig-el)

      ELSE IF(ICSPA.EQ.2) THEN

        DELDL = 0.0808D0
        EPSDL = -0.4525D0
        S = SQS*SQS
        SIGT = 11.82D0*S**DELDL+(26.36D0+ 8.15D0)/2.D0*S**EPSDL
        IMODEL = 1
        IF(IMODEL.EQ.1) THEN
          ALPHAP = 0.25D0
          SLOPE = 8.5D0+2.D0*ALPHAP*LOG(S)
        ELSE IF(IMODEL.EQ.2) THEN
          ALPHAP = 0.3D0
          SLOPE = 8.D0+2.D0*ALPHAP*LOG(S)
        ENDIF
        SIGEL = SIGT**2/(16.D0*PI*SLOPE*CMBARN)
        SIGINEL = SIGT-SIGEL
        RHO = 0.D0

*  geometrical scaling with Donnachie-Landshoff sig-tot

      ELSE IF(ICSPA.EQ.3) THEN

        R = 0.17D0

        DELDL = 0.0808D0
        EPSDL = -0.4525D0
        S = SQS*SQS
        SIGT = 11.82D0*S**DELDL+(26.36D0+ 8.15D0)/2.D0*S**EPSDL

        SIGEL = R*SIGT
        SIGINEL = SIGT-SIGEL
        SLOPE = SIGT**2/(16.D0*PI*SIGEL)/CMBARN
        RHO = 0.D0
        
c ICSPA=4 reserved for CONEX_EXTENSION
c      ELSE IF(ICSPA.EQ.4) THEN


*  cross section from 2014 Review of Particle Physics
        
      ELSE IF(ICSPA.EQ.5) THEN
         
c     elastic slope not included in fit
c     taking slope parameterization from sigma_pp Donnie.-Landshoff
         ALPHAP = 0.25D0
         SLOPE = 8.5D0+4.D0*ALPHAP*LOG(SQS)
         
         CALL SIG_RPP2014(3,1,SQS,SLOPE,SIGT,SIGEL,SIGINEL,RHO)

      ENDIF

      RETURN
      END

C=======================================================================

      SUBROUTINE SIGMA_INI 

C-----------------------------------------------------------------------
C.  Initialize the cross section and interaction lengths in air
C.  cross section model can be chosen, per particle, by setting ICSPA2()
C.  default is Sibyll cross section (0,0,0)      
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON /CSAIR/ ASQSMIN, ASQSMAX, DASQS,
     &     SSIG0(61,3),SSIGA(61,3),ALINT(61,3),NSQS
      
      COMMON /CSPA/ ICSPA2(3)

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      SAVE
      DATA ICSPA2 /0,0,0/
      DATA AVOG /6.0221367D-04/
      DATA ATARGET /14.514D0/            ! effective masss of air

      IF(NDEBUG.gt.0)
     &     write(lun,*) ' SIGMA_INI: using cross section model no.',
     &     (ICSPA2(i),i=1,3)

      CALL BLOCK_INI

C...Loop on c.m. energy 
      NSQS = 61
      SQSMIN = 10.D0
      SQSMAX = 1.d+07
      ASQSMIN = LOG10(SQSMIN)
      ASQSMAX = LOG10(SQSMAX)
      DASQS = (ASQSMAX-ASQSMIN)/DBLE(NSQS-1)
      DO J=1,NSQS
         ASQS = ASQSMIN + DASQS*DBLE(J-1)
         SQS = 10.D0**ASQS
         E0 = SQS*SQS/(2.D0*0.938D0) * 1.D-03       ! TeV
C...p-air
         CALL SIGMA_PP (E0, SIGT, SIGEL, SIGINEL, SLOPE, RHO)
C     using parametrization by Goulianos for diff. cross section
c     (depends on elastic cross section)
c     used to determine coupling to intermediate resonances in Glauber calc (ALAM)
c     assumed to be universal, i.e. same coupling used for proton, pion and kaons
         CALL SIB_HADCS1(1,SQS,SIGT1,SIGEL1,SIGINEL1,SLOPE1,RHO1)
         SIGEFF = 0.68D0*(1.D0+36.D0/SQS**2)
     &        *LOG(0.6D0+0.02D0/1.5D0*SQS**2)
         SIGEFF = MAX(0.D0,SIGEFF)
         ALAM = sqrt(SIGEFF/SIGEL1)
         SSIGSD = 2.D0 * SIGEFF        
         CALL SIG_H_AIR (SIGT, SLOPE, RHO, ALAM,
     &        SSIGT, SSIGEL, SSIGQE, SIGSD, SIGQSD )
         SSIGA(J,1) = SSIGT-SSIGQE ! had-air production cross section
         SSIG0(J,1) = SIGINEL   ! had-nucleon inel. cross section
         ALINT(J,1) = 1.D0/(AVOG*SSIGA(J,1)/ATARGET) ! interaction length in air
C...pi-air
         CALL SIGMA_PIP (E0, SIGT, SIGEL, SIGINEL, SLOPE, RHO) 
         CALL  SIG_H_AIR (SIGT, SLOPE, RHO, ALAM,
     &        SSIGT, SSIGEL, SSIGQE, SIGSD, SIGQSD )
         SSIGA(J,2) = SSIGT-SSIGQE
         SSIG0(J,2) = SIGINEL
         ALINT(J,2) = 1.D0/(AVOG*SSIGA(J,2)/ATARGET)
C...K-air
         CALL SIGMA_KP (E0, SIGT, SIGEL, SIGINEL, SLOPE, RHO) 
         CALL  SIG_H_AIR (SIGT, SLOPE, RHO, ALAM,
     &        SSIGT, SSIGEL, SSIGQE, SIGSD, SIGQSD )
         SSIGA(J,3) = SSIGT-SSIGQE
         SSIG0(J,3) = SIGINEL
         ALINT(J,3) = 1.D0/(AVOG*SSIGA(J,3)/ATARGET)
      ENDDO

      if (ndebug .gt. 0 ) THEN
        WRITE(LUN,'(1X,A)') 
     &  ' SIGMA_INI: NUCLIB interaction lengths [g/cm**2]'
        WRITE(LUN,'(1X,A)') 
     &  '     sqs,       p-air,      pi-air,     K-air'
      DO J=1,NSQS
         SQS = 10.D0**(ASQSMIN + DASQS*DBLE(J-1))
         WRITE(LUN,'(1X,1P,4E12.3)') 
     &        SQS,ALINT(J,1),ALINT(J,2),ALINT(J,3)
        ENDDO
      endif

      RETURN
      END

C=======================================================================

      FUNCTION FPNI (E,Linp)

C-----------------------------------------------------------------------
C...This function  returns the interaction length 
C.  of an hadronic particle travelling in air
C.
C.  INPUT:   E (TeV)   particle energy
C.           Linp      particle code
C.  OUTPUT:  FPNI      (g cm-2)
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
            
      COMMON /CSAIR/ ASQSMIN, ASQSMAX, DASQS,
     &     SSIG0(61,3),SSIGA(61,3),ALINT(61,3),NSQS

      DIMENSION KK(6:14)
      SAVE
      DATA KK /3*2, 4*3, 2*1/

      SQS = SQRT(2000.D0*E*0.937D0)                        ! GeV
      AL = LOG10 (SQS)
      L = abs(Linp)
      IF (AL .LE. ASQSMIN)  THEN
         FPNI = ALINT(1,KK(L))
      ELSE
         T = (AL-ASQSMIN)/DASQS
         J = INT(T)
         J = MIN(J,NSQS-2)
         T = T-DBLE(J)
         FPNI = ((1.D0-T)*ALINT(J+1,KK(L)) + T*ALINT(J+2,KK(L)))
      ENDIF
      RETURN
      END

C=======================================================================
      
      FUNCTION FSIGHAIR (E,Linp)

C-----------------------------------------------------------------------
C...This function returns the production cross section
C.  of an hadronic particle with air calculated in NUCLIB (SIGMA_INI)     
C.
C.  INPUT:   E (TeV)   particle energy
C.           Linp      particle code
C.  OUTPUT:  SIG_PROD  (mb)
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
            
      COMMON /CSAIR/ ASQSMIN, ASQSMAX, DASQS,
     &     SSIG0(61,3),SSIGA(61,3),ALINT(61,3),NSQS

      DIMENSION KK(6:14)
      SAVE
      DATA KK /3*2, 4*3, 2*1/

      SQS = SQRT(2000.D0*E*0.937D0)                        ! GeV
      AL = LOG10 (SQS)
      L = abs(Linp)
      IF (AL .LE. ASQSMIN)  THEN
         FSIGHAIR = SSIGA(1,KK(L))
      ELSE
         T = (AL-ASQSMIN)/DASQS
         J = INT(T)
         J = MIN(J,NSQS-2)
         T = T-DBLE(J)
         FSIGHAIR = ((1.D0-T)*SSIGA(J+1,KK(L)) + T*SSIGA(J+2,KK(L)))
      ENDIF     
      RETURN
      END

C=======================================================================

      SUBROUTINE INT_LEN_INI

C-----------------------------------------------------------------------
C...Initialize the interaction lengths from NUCLIB
C-----------------------------------------------------------------------
      SAVE
      
      CALL NUC_GEOM_INI                 ! nucleus profiles
      CALL SIGMA_INI                    ! initialize cross sections

      RETURN
      END
C=======================================================================

      SUBROUTINE TRANSFONSHELL(ECM,XM1in,XM2in,XMAX,IMOD,P1,P2,LBAD)

C-----------------------------------------------------------------------
C     samples 2 --> 2 scattering that puts a particle on its mass shell
C
C     particle1 is along +z, always receives mass
C     particle2 is along -z, mass only sampled if both aquire mass
C
C     DEPENDS: slope-parameter in s_difmass
C
C     INPUT: ECM : center-of-mass energy of scattering particles
C            M1in  : mass of first particle
C            M2in  : mass of second particle
C            XMAX  : maximal mass that can be obtained
C            IMOD  : remnant or diffraction mode
C     
C     OUTPUT: P1,P2 : final state 4vectors in two-particle c.m.   \FR'14
C-----------------------------------------------------------------------
      IMPLICIT NONE
      
c     external types
      DOUBLE PRECISION ECM,XM1in,XM2in,XMAX,P1,P2
      DIMENSION P1(5),P2(5)
      INTEGER IMOD,LBAD

      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      INTEGER NIPAR_max,NPAR_max
      PARAMETER (NPAR_max=200,NIPAR_max=100)
      DOUBLE PRECISION PAR
      INTEGER IPAR
      COMMON /S_CFLAFR/ PAR(NPAR_max), IPAR(NIPAR_max)
      INTEGER ITRY, NREJ
      COMMON /S_CNT/ ITRY(20), NREJ(20)
      DOUBLE PRECISION AM,AM2
      COMMON /S_MASS1/ AM(99), AM2(99)

c     internal types
      DOUBLE PRECISION XMB2,XMT2,AXMX,S,X1,X2,ALX,SLOP0,SLOPE,DB,
     &     T,PTS,PZB2,PZT2,PT,PHI,XMB,XMT,S_RNDM,PTSWTCH

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN
      DOUBLE PRECISION SLOP0_0,ASLOP,BSLOP
      INTEGER II
      SAVE
      DATA SLOP0_0 /6.5D0/        ! b (slope_ for Mx**2 > 5 GeV**2
      DATA ASLOP /31.10362D0/     ! fit to the slope parameter.
      DATA BSLOP /-15.29012D0/

      IF(NDEBUG.gt.3)
     &     WRITE(LUN,*) ' TRANSFONSHELL: called with (Ecm,M1,M2,XMAX):',
     &     ECM,XM1in,XM2in,XMAX
      
      XMB2 = XM1in**2
      XMT2 = XM2in**2     

      AXMX = LOG(XMAX)
      
      ITRY(6) = 0
      LBAD = 1

C     remnant pt parameters
c     distribution is: exp(-slope*t)
c     slope = aslop + bslop * log(Mx**2)
c     (by default same as in diff.
c      scale with paramterers 90 and 91)

c     diff. pt paramters
      ASLOP = PAR(133)
      BSLOP = PAR(134)
      SLOP0_0 = PAR(135)
      
      S = ECM*ECM
      X1 = 1.D0-(XMT2-XMB2)/S
      X2 = 2.D0-X1
      IF(X2.LT.EPS5) RETURN

 60   ITRY(6) = ITRY(6) + 1
      IF(ITRY(6).GT.NREJ(6)) RETURN
c     sample transverse momentum
      ALX = LOG(MAX(XMT2,XMB2))
c     set slope of pt distribution
      IF(IMOD.eq.0)THEN
c     diffraction dissociation
         SLOP0 = SLOP0_0*PAR(93)
         SLOPE = MAX(SLOP0,ASLOP+BSLOP*ALX)
         PTSWTCH = 1.D0

      ELSEIF(IMOD.eq.1)THEN
c     remnant excitation
         IF(IPAR(57).eq.0)THEN
            ALX = ALX-LOG(AM2(13))
            SLOP0 = PAR(92)
            DB = (SLOP0-PAR(90))/AXMX
            SLOPE = MAX(SLOP0,PAR(90)+DB*PAR(91)*ALX)
         ELSE
            ALX = ALX-LOG(AM2(13))
            SLOP0 = PAR(92)
            SLOPE = MAX(SLOP0,PAR(90)+PAR(91)*ALX)
         ENDIF
         PTSWTCH = 1.D0

      ELSEIF(IMOD.eq.3)THEN
c     no pt
         PTSWTCH = 0.D0
         SLOPE = 1.D0
      ENDIF
      IF(ndebug.gt.3)
     &     WRITE(LUN,*) ' TRANSFONSHELL: (SLOP0,SLOPE,log(M**2)):',
     &     SLOP0,SLOPE,ALX
      T = -DLOG(MAX(EPS10,S_RNDM(0)))/SLOPE
      PTS = T*X1*PTSWTCH
      PZB2 = S*0.5D0*0.5D0*X1*X1-XMB2-PTS
      PZT2 = S*0.5D0*0.5D0*X2*X2-XMT2-PTS
      IF(NDEBUG.gt.3) 
     &     WRITE(LUN,*) ' TRANSFONSHELL: (PTS,PZB2,PZT2):',PTS,PZB2,PZT2
c      IF (ABS(PZB2)-PZT2.GT.EPS10) GOTO 60
      IF (PZB2.lt.0.D0.or.PZT2.LT.0.D0) GOTO 60
      PT = DSQRT(PTS)
      PHI = TWOPI*S_RNDM(1)
      XMB = sqrt(XMB2)
      XMT = sqrt(XMT2)
      P2(4) = 0.5D0*ECM*X2
      P2(3) = -DSQRT(PZT2)
      P2(1) = PT*dCOS(PHI)
      P2(2) = PT*dSIN(PHI)
      P2(5) = XMT

      P1(4) = 0.5D0*ECM*X1
      P1(3) = DSQRT(PZB2)
      do ii = 1,2
         P1(ii) = -P2(ii)
      enddo
      P1(5) = XMB
      IF(NDEBUG.gt.3) THEN
          WRITE(LUN,*) ' TRANSFONSHELL: (P1):',(p1(ii),ii=1,5)
          WRITE(LUN,*) ' TRANSFONSHELL: (P2):',(p2(ii),ii=1,5)
       ENDIF
      LBAD = 0
      END
C=======================================================================

      SUBROUTINE SAMPLE_SEA (ALPHA,ASUP,XMASS,XMAX,X1,X2,PT)

C-----------------------------------------------------------------------
C.    Routine that samples the kinematical variables of a sea quark pair.
C.  INPUT:  STR_mass_min : minimal string mass ** 2 = x1 * x2 * s
C.          ASUP : large x suppression exponent
C.  OUTPUT:  X1, X2, PT (GeV)                                   /FR'14
C-----------------------------------------------------------------------
Cf2py double precision, intent(in) :: ALPHA,ASUP,XMASS,XMAX
Cf2py double precision, intent(out) :: X1,X2,PT
      IMPLICIT NONE

c     include COMMONs
      INTEGER NCALL, NDEBUG, LUN
      COMMON /S_DEBUG/ NCALL, NDEBUG, LUN
      INTEGER NW_max
      PARAMETER (NW_max = 20)
C--------------------------------------------------------------------
C     SIBYLL common blocks containing event information       \FR'14
C--------------------------------------------------------------------

C     EVENT INFO COMMON
C     contains overall interaction properties, like
C     SQS : center-of-mass energy
C     S   :         "       "     squared
C     PTmin : low pt cut of QCD cross section, 
C             i.e. minimal pt of hard minijets
C     Xmin : low-x bound for PDFs, 
C            i.e. minimal momentum fraction of hard partons
C     Zmin : logarithm of that
C     KB : PID of beam hadron
C     KT() : PID of target
C     IAT : mass number of target
      DOUBLE PRECISION SQS,S,PTmin,XMIN,ZMIN
      INTEGER KB,IAT,KT
      COMMON /S_RUN/ SQS, S, PTmin, XMIN, ZMIN, KB, KT(NW_max), IAT

C--------------------------------------------------------------------
C     SIBYLL utility common blocks containing constants       \FR'14
C--------------------------------------------------------------------
      DOUBLE PRECISION EPS3,EPS5,EPS8,EPS10
      COMMON /SIB_EPS/ EPS3,EPS5,EPS8,EPS10

      DOUBLE PRECISION PI,TWOPI,CMBARN
      COMMON /SIB_CST/ PI,TWOPI,CMBARN

      DOUBLE PRECISION FACN
      DIMENSION FACN(3:10)
      COMMON /SIB_FAC/ FACN

c     external type declarations
      DOUBLE PRECISION ALPHA,ASUP,XMASS,XMAX,X1,X2,PT

c     internal types
      DOUBLE PRECISION XMINA,XM2DIS,XR,SLOPE,S_RNDM,XRNDM
      SAVE
      
      IF(ndebug.gt.3)
     &    write(lun,*) ' SAMPLE_SEA: alpha,asup,qmass,xmax',
     &    ALPHA,ASUP,XMASS,XMAX

c     min. momentum fraction for massive quarks
c     i.e. sample from 1/(x+x_min)
      XMINA = 2.D0*XMASS/SQS
      IF(ndebug.gt.3)
     &     write(lun,*) ' SAMPLE_SEA: xmina:',XMINA
c     exponent of large x suppression: (1-x)**b, b=0 or b>1
      IF(ABS(ASUP).lt.EPS3)THEN
c     b = 0 , no suppression, sample bare 1/(x+xmin)       
         X1 = XM2DIS(XMINA,XMAX,ALPHA) ! ~(1/x)**alpha
         X2 = XM2DIS(XMINA,XMAX,ALPHA) ! ~(1/x)**alpha
         
      ELSEIF(ASUP.ge.EPS3)THEN
c     b >= 1 , sample bare (1-x)**b/(x+xmin)
         SLOPE = MAX(ASUP,EPS3)
c     quark
 100     X1 = XM2DIS(XMINA,XMAX,ALPHA) ! ~(1/x)**alpha
         XR = LOG(1.D0-X1)-LOG(1.D0-XMINA)
         XRNDM = S_RNDM(1)
         IF(ndebug.gt.4)
     &        write(lun,*) '  X1,XR,SLOPE*XR:',X1,XR,SLOPE*XR
         if(SLOPE*XR.le.LOG(max(XRNDM,eps10))) goto 100

c     anti-quark
 200     X2 = XM2DIS(XMINA,XMAX,ALPHA) ! ~(1/x)**alpha
         XR = log(1.D0-X2)-log(1.D0-XMINA)
         XRNDM = S_RNDM(2)
         IF(ndebug.gt.4)
     &        write(lun,*) '  X2,XR,SLOPE*XR,XRNDM:',
     &    X2,XR,SLOPE*XR,XRNDM
         if(SLOPE*XR.le.log(max(XRNDM,eps10))) goto 200     
      ELSE
         WRITE(LUN,*) ' SAMPLE_SEA: suppression exponent out of range.'
         WRITE(LUN,*) ' SAMPLE_SEA: ASUP:',ASUP
         STOP
      ENDIF

c     sample pt
c     not yet implemented... to avoid problem with virtual partons
      pt = 0.D0
      IF(ndebug.gt.3)
     &     write(lun,*) ' SAMPLE_SEA: X1,X2,PT:',X1,X2,PT

      END
C**********************************************
C
C     contains the src for pion and proton pdf
C     parametrizations according to GRV
C     ( see function head for refs. )
C
C     1 pion pdf
C     2 proton pdf GRV98LO
C
C**********************************************

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
*         G R V - P I O N - P A R A M E T R I Z A T I O N S       *
*                                                                 *
*                 FOR A DETAILED EXPLANATION SEE :                *
*              M. GLUECK, E.REYA, A.VOGT: DO-TH 91/16             *
*                                                                 *
*   THE PARAMETRIZATIONS ARE FITTED TO THE PARTON DISTRIBUTIONS   *
*   FOR Q ** 2 BETWEEN MU ** 2 (=  0.25 / 0.30  GEV ** 2  IN LO   *
*   / HO) AND  1.E8 GEV ** 2  AND FOR X BETWEEN  1.E-5  AND  1.   *
*   REGIONS, WHERE THE DISTRIBUTION UNDER CONSIDERATION IS NEG-   *
*   LIGIBLE, I.E. BELOW ABOUT 1.E-4, WERE EXCLUDED FROM THE FIT.  *
*                                                                 *
*              HEAVY QUARK THRESHOLDS  Q(H) = M(H) :              *
*         M(C)  =  1.5,  M(B)  =  4.5,  M(T)  =  100  GEV         *
*                                                                 *
*      CORRESPONDING LAMBDA(F) VALUES FOR F ACTIVE FLAVOURS :     *
*      LO :   LAMBDA(3)  =  0.232,   LAMBDA(4)  =  0.200,         *
*             LAMBDA(5)  =  0.153,   LAMBDA(6)  =  0.082  GEV     *
*      HO :   LAMBDA(3)  =  0.248,   LAMBDA(4)  =  0.200,         *
*             LAMBDA(5)  =  0.131,   LAMBDA(6)  =  0.053  GEV     *
*                                                                 *
*   HO DISTRIBUTION REFER TO THE MS-BAR SCHEME OF BARDEEN ET AL.  *
*                                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C=======================================================================

       SUBROUTINE DORPLO (X, Q2, VAP, GLP, QBP, CBP, BBP)

C-----------------------------------------------------------------------
       IMPLICIT DOUBLE PRECISION (A - Z)
       SAVE

       MU2  = 0.25D0
       LAM2 = 0.232D0 * 0.232D0
       S  = LOG (LOG(Q2/LAM2) / LOG(MU2/LAM2))
       DS = SQRT (S)
       S2 = S * S
C...X * VALENCE :
       NV  =  0.519D0 + 0.180D0 * S - 0.011D0 * S2
       AKV =  0.499D0 - 0.027D0 * S
       AGV =  0.381D0 - 0.419D0 * S
       DV  =  0.367D0 + 0.563D0 * S
       VAP =  DORFVP (X, NV, AKV, AGV, DV)
C...X * GLUON :
       ALG =  0.599D0
       BEG =  1.263D0
       AKG =  0.482D0 + 0.341D0 * DS
       BKG =   0.0D0
       AGG =  0.678D0 + 0.877D0 * S  - 0.175D0 * S2
       BGG =  0.338D0 - 1.597D0 * S
       CG  =   0.0D0  - 0.233D0 * S  + 0.406D0 * S2
       DG  =  0.390D0 + 1.053D0 * S
       EG  =  0.618D0 + 2.070D0 * S
       ESG =  3.676D0
       GLP =DORFGP(X, S, ALG, BEG, AKG, BKG, AGG, BGG, CG, DG, EG, ESG)
C...X * QBAR (SU(3)-SYMMETRIC SEA) :
       SL  =   0.0D0
       ALS =   0.55D0
       BES =   0.56D0
       AKS =  2.538D0 - 0.763D0 * S
       AGS = -0.748D0
       BS  =  0.313D0 + 0.935D0 * S
       DS  =  3.359D0
       EST =  4.433D0 + 1.301D0 * S
       ESS =   9.30D0 - 0.887D0 * S
       QBP =  DORFQP (X, S, SL, ALS, BES, AKS, AGS, BS, DS, EST, ESS)
C...X * CBAR = X * C :
       SC  =  0.888D0
       ALC =   1.02D0
       BEC =   0.39D0
       AKC =   0.0D0
       AGC =   0.0D0
       BC  =  1.008D0
       DC  =  1.208D0 + 0.771D0 * S
       EC  =   4.40D0 + 1.493D0 * S
       ESC =  2.032D0 + 1.901D0 * S
       CBP =  DORFQP (X, S, SC, ALC, BEC, AKC, AGC, BC, DC, EC, ESC)
C...X * BBAR = X * B :
       SBO =  1.351D0
       ALB =   1.03D0
       BEB =   0.39D0
       AKB =   0.0D0
       AGB =   0.0D0
       BBO =   0.0D0
       DB  =  0.697D0 + 0.855D0 * S
       EB  =   4.51D0 + 1.490D0 * S
       ESB =  3.056D0 + 1.694D0 * S
       BBP =  DORFQP (X, S, SBO, ALB, BEB, AKB, AGB, BBO, DB, EB, ESB)
       RETURN
       END
C
C=======================================================================

       SUBROUTINE DORPHO (X, Q2, VAP, GLP, QBP, CBP, BBP)

C-----------------------------------------------------------------------
       IMPLICIT DOUBLE PRECISION (A - Z)
       SAVE

       MU2  = 0.3D0
       LAM2 = 0.248D0 * 0.248D0
       S  = LOG (LOG(Q2/LAM2) / LOG(MU2/LAM2))
       DS = SQRT (S)
       S2 = S * S
C...X * VALENCE :
       NV  =  0.456D0 + 0.150D0 * DS + 0.112D0 * S - 0.019D0 * S2
       AKV =  0.505D0 - 0.033D0 * S
       AGV =  0.748D0 - 0.669D0 * DS - 0.133D0 * S
       DV  =  0.365D0 + 0.197D0 * DS + 0.394D0 * S
       VAP =  DORFVP (X, NV, AKV, AGV, DV)
C...X * GLUON :
       ALG =  1.096D0
       BEG =  1.371D0
       AKG =  0.437D0 - 0.689D0 * DS
       BKG = -0.631D0
       AGG =  1.324D0 - 0.441D0 * DS - 0.130D0 * S
       BGG = -0.955D0 + 0.259D0 * S
       CG  =  1.075D0 - 0.302D0 * S
       DG  =  1.158D0 + 1.229D0 * S
       EG  =   0.0D0  + 2.510D0 * S
       ESG =  2.604D0 + 0.165D0 * S
       GLP =DORFGP(X, S, ALG, BEG, AKG, BKG, AGG, BGG, CG, DG, EG, ESG)
C...X * QBAR (SU(3)-SYMMETRIC SEA) :
       SL  =   0.0D0
       ALS =   0.85D0
       BES =   0.96D0
       AKS = -0.350D0 + 0.806D0 * S
       AGS = -1.663D0
       BS  =  3.148D0
       DS  =  2.273D0 + 1.438D0 * S
       EST =  3.214D0 + 1.545D0 * S
       ESS =  1.341D0 + 1.938D0 * S
       QBP =  DORFQP (X, S, SL, ALS, BES, AKS, AGS, BS, DS, EST, ESS)
C...X * CBAR = X * C :
       SC  =  0.820D0
       ALC =   0.98D0
       BEC =   0.0D0
       AKC =   0.0D0  - 0.457D0 * S
       AGC =   0.0D0
       BC  =  -1.00D0 +  1.40 D0* S
       DC  =  1.318D0 + 0.584D0 * S
       EC  =   4.45D0 + 1.235D0 * S
       ESC =  1.496D0 + 1.010D0 * S
       CBP =  DORFQP (X, S, SC, ALC, BEC, AKC, AGC, BC, DC, EC, ESC)
C...X * BBAR = X * B :
       SBO =  1.297D0
       ALB =   0.99D0
       BEB =   0.0D0
       AKB =   0.0D0  - 0.172D0 * S
       AGB =   0.0D0
       BBO =   0.0D0
       DB  =  1.447D0 + 0.485D0 * S
       EB  =   4.79D0 + 1.164D0 * S
       ESB =  1.724D0 + 2.121D0 * S
       BBP =  DORFQP (X, S, SBO, ALB, BEB, AKB, AGB, BBO, DB, EB, ESB)
       RETURN
       END
C
C=======================================================================

       FUNCTION DORFVP (X, N, AK, AG, D)

C-----------------------------------------------------------------------
       IMPLICIT DOUBLE PRECISION (A - Z)
       SAVE

       DX = SQRT (X)
       DORFVP = N * X**AK * (1.D0+ AG*DX) * (1.D0- X)**D
       RETURN
       END
C
C=======================================================================

       FUNCTION DORFGP (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)

C-----------------------------------------------------------------------
       IMPLICIT DOUBLE PRECISION (A - Z)
       SAVE

       DX = SQRT (X)
       LX = LOG (1.D0/X)
       DORFGP = (X**AK * (AG + BG*DX + C*X) * LX**BK + S**AL
     1       * EXP (-E + SQRT (ES * S**BE * LX))) * (1.D0- X)**D
       RETURN
       END
C
C=======================================================================

       FUNCTION DORFQP (X, S, ST, AL, BE, AK, AG, B, D, E, ES)

C-----------------------------------------------------------------------
       IMPLICIT DOUBLE PRECISION (A - Z)
       SAVE

       DX = SQRT (X)
       LX = LOG (1./X)
       IF (S .LE. ST) THEN
          DORFQP = 0.0D0
       ELSE
          DORFQP = (S-ST)**AL / LX**AK * (1.D0+AG*DX+B*X) * (1.D0- X)**D
     1           * EXP(-E + SQRT(ES * S**BE * LX))
       END IF
       RETURN
       END
C=======================================================================

      DOUBLE PRECISION FUNCTION SIB_DOR92FS(X,S,ST,AL,BE,AK,AG,B,D,E,ES)

C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A - Z)
      SAVE

       DX = SQRT (X)
       LX = LOG (1.D0/X)
       IF (S .LE. ST) THEN
         SIB_DOR92FS = 0.D0
       ELSE
         SIB_DOR92FS = (S-ST)**AL/LX**AK*(1.D0+AG*DX+B*X)*(1.D0-X)**D
     1          * EXP (-E + SQRT (ES * S**BE * LX))
       END IF

      END

C=======================================================================

      DOUBLE PRECISION FUNCTION SIB_DBFINT(NARG,ARG,NA,ENT,TABLE)

C-----------------------------------------------------------------------
C
C     routine based on CERN library E104
C
C     multi-dimensional interpolation routine, needed for PHOJET
C     internal cross section tables and several PDF sets (GRV98 and AGL)
C
C     changed to avoid recursive function calls (R.Engel, 09/98)
C
C***********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      INTEGER NA(NARG), INDE(32)
      DOUBLE PRECISION ARG(NARG),ENT(NARG),TABLE(*),WEIGHT(32)
      SAVE


      DBFINT     =  0.D0
      SIB_DBFINT =  0.D0
      IF(NARG .LT. 1  .OR.  NARG .GT. 5)  RETURN

           LMAX      =  0
           ISTEP     =  1
           KNOTS     =  1
           INDE(1)   =  1
           WEIGHT(1) =  1.D0
           DO 100    N  =  1, NARG
              X     =  ARG(N)
              NDIM  =  NA(N)
              LOCA  =  LMAX
              LMIN  =  LMAX + 1
              LMAX  =  LMAX + NDIM
              IF(NDIM .GT. 2)  GOTO 10
              IF(NDIM .EQ. 1)  GOTO 100
              H  =  X - ENT(LMIN)
              IF(ABS(H) .LT. 0.D-8)  GOTO 90
              ISHIFT  =  ISTEP
              IF(ABS(X-ENT(LMIN+1)) .LT. 0.D-8)  GOTO 21
              ISHIFT  =  0
              ETA     =  H / (ENT(LMIN+1) - ENT(LMIN))
              GOTO 30
   10         LOCB  =  LMAX + 1
   11         LOCC  =  (LOCA+LOCB) / 2
              IF(X-ENT(LOCC))  12, 20, 13
   12         LOCB  =  LOCC
              GOTO 14
   13         LOCA  =  LOCC
   14         IF(LOCB-LOCA .GT. 1)  GOTO 11
              LOCA    =  MIN ( MAX (LOCA,LMIN), LMAX-1 )
              ISHIFT  =  (LOCA - LMIN) * ISTEP
              ETA     =  (X - ENT(LOCA)) / (ENT(LOCA+1) - ENT(LOCA))
              GOTO 30
   20         ISHIFT  =  (LOCC - LMIN) * ISTEP
   21         DO 22  K  =  1, KNOTS
                 INDE(K)  =  INDE(K) + ISHIFT
   22         CONTINUE
              GOTO 90
   30         DO 31  K  =  1, KNOTS
                 INDE(K)         =  INDE(K) + ISHIFT
                 INDE(K+KNOTS)   =  INDE(K) + ISTEP
                 WEIGHT(K+KNOTS) =  WEIGHT(K) * ETA
                 WEIGHT(K)       =  WEIGHT(K) - WEIGHT(K+KNOTS)
   31         CONTINUE
              KNOTS  =  2*KNOTS
   90         ISTEP  =  ISTEP * NDIM
  100      CONTINUE
           DO 200    K  =  1, KNOTS
              I  =  INDE(K)
              DBFINT =  DBFINT + WEIGHT(K) * TABLE(I)
  200      CONTINUE

      SIB_DBFINT = DBFINT

      END

C=======================================================================

      SUBROUTINE SIB_DOR98LO (Xinp, Q2inp, UV, DV, US, DS, SS, GL)

C-----------------------------------------------------------------------
C***********************************************************************
C
C   GRV98 parton densities, leading order set
C
C                  For a detailed explanation see
C                   M. Glueck, E. Reya, A. Vogt :
C        hep-ph/9806404  =  DO-TH 98/07  =  WUE-ITP-98-019
C                  (To appear in Eur. Phys. J. C)
C
C   interpolation routine based on the original GRV98PA routine,
C   adapted to define interpolation table as DATA statements
C
C                                                   (R.Engel, 09/98)
C
C
C   INPUT:   X  =  Bjorken-x        (between  1.E-9 and 1.)
C            Q2 =  scale in GeV**2  (between  0.8 and 1.E6)
C
C   OUTPUT:  UV = u - u(bar),  DV = d - d(bar),  US = u(bar),
C            DS = d(bar),  SS = s = s(bar),  GL = gluon.
C            Always x times the distribution is returned.
C
C******************************************************i****************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER(I-N)

      PARAMETER (NX=68, NQ=27, NARG=2)
      DIMENSION XUVF(NX,NQ), XDVF(NX,NQ), XDEF(NX,NQ), XUDF(NX,NQ),
     1          XSF(NX,NQ), XGF(NX,NQ),
     2          XT(NARG), NA(NARG), ARRF(NX+NQ)

      DIMENSION XUVF_L(NX*NQ), XDVF_L(NX*NQ), XDEF_L(NX*NQ),
     &  XUDF_L(NX*NQ), XSF_L(NX*NQ), XGF_L(NX*NQ)

      EQUIVALENCE (XUVF(1,1),XUVF_L(1))
      EQUIVALENCE (XDVF(1,1),XDVF_L(1))
      EQUIVALENCE (XDEF(1,1),XDEF_L(1))
      EQUIVALENCE (XUDF(1,1),XUDF_L(1))
      EQUIVALENCE (XSF(1,1),XSF_L(1))
      EQUIVALENCE (XGF(1,1),XGF_L(1))
      SAVE

      DATA (ARRF(K),K=    1,   95) /
     &  -2.0723D+01,-2.0135D+01,-1.9560D+01,-1.8983D+01,-1.8421D+01,
     &  -1.7833D+01,-1.7258D+01,-1.6680D+01,-1.6118D+01,-1.5530D+01,
     &  -1.4955D+01,-1.4378D+01,-1.3816D+01,-1.3479D+01,-1.3122D+01,
     &  -1.2717D+01,-1.2311D+01,-1.1913D+01,-1.1513D+01,-1.1176D+01,
     &  -1.0820D+01,-1.0414D+01,-1.0009D+01,-9.6108D+00,-9.2103D+00,
     &  -8.8739D+00,-8.5172D+00,-8.1117D+00,-7.7063D+00,-7.3082D+00,
     &  -6.9078D+00,-6.5713D+00,-6.2146D+00,-5.8091D+00,-5.4037D+00,
     &  -5.0056D+00,-4.6052D+00,-4.2687D+00,-3.9120D+00,-3.5066D+00,
     &  -3.1011D+00,-2.8134D+00,-2.5257D+00,-2.3026D+00,-2.0794D+00,
     &  -1.8971D+00,-1.7430D+00,-1.6094D+00,-1.4917D+00,-1.3863D+00,
     &  -1.2910D+00,-1.2040D+00,-1.1239D+00,-1.0498D+00,-9.8083D-01,
     &  -9.1629D-01,-7.9851D-01,-6.9315D-01,-5.9784D-01,-5.1083D-01,
     &  -4.3078D-01,-3.5667D-01,-2.8768D-01,-2.2314D-01,-1.6252D-01,
     &  -1.0536D-01,-5.1293D-02, 0.0000D+00,-2.2314D-01, 0.0000D+00,
     &   2.6236D-01, 5.8779D-01, 9.9325D-01, 1.3863D+00, 1.8563D+00,
     &   2.3026D+00, 2.7726D+00, 3.2189D+00, 3.6889D+00, 4.1589D+00,
     &   4.6052D+00, 5.1930D+00, 5.7683D+00, 6.3456D+00, 6.9078D+00,
     &   7.4955D+00, 8.0709D+00, 8.6482D+00, 9.2103D+00, 9.9988D+00,
     &   1.0736D+01, 1.1513D+01, 1.2301D+01, 1.3039D+01, 1.3816D+01/
      DATA (XUVF_L(K),K=    1,  114) /
     &2.3186D+00,2.2915D+00,2.2645D+00,2.2385D+00,2.2140D+00,2.1876D+00,
     &2.1623D+00,2.1366D+00,2.1121D+00,2.0862D+00,2.0612D+00,2.0358D+00,
     &2.0110D+00,1.9963D+00,1.9806D+00,1.9624D+00,1.9446D+00,1.9263D+00,
     &1.9072D+00,1.8904D+00,1.8724D+00,1.8515D+00,1.8294D+00,1.8085D+00,
     &1.7865D+00,1.7680D+00,1.7483D+00,1.7249D+00,1.6993D+00,1.6715D+00,
     &1.6385D+00,1.6141D+00,1.5884D+00,1.5597D+00,1.5337D+00,1.5121D+00,
     &1.4985D+00,1.4980D+00,1.5116D+00,1.5555D+00,1.6432D+00,1.7434D+00,
     &1.8861D+00,2.0327D+00,2.2174D+00,2.4015D+00,2.5849D+00,2.7671D+00,
     &2.9488D+00,3.1308D+00,3.3142D+00,3.4998D+00,3.6885D+00,3.8826D+00,
     &4.0815D+00,4.2069D+00,4.5481D+00,4.8830D+00,5.2116D+00,5.5351D+00,
     &5.8553D+00,6.1665D+00,6.4745D+00,6.7767D+00,7.0735D+00,7.3628D+00,
     &7.6283D+00,0.0000D+00,2.3948D+00,2.3665D+00,2.3388D+00,2.3126D+00,
     &2.2860D+00,2.2592D+00,2.2327D+00,2.2065D+00,2.1810D+00,2.1541D+00,
     &2.1284D+00,2.1020D+00,2.0760D+00,2.0605D+00,2.0443D+00,2.0259D+00,
     &2.0068D+00,1.9873D+00,1.9676D+00,1.9500D+00,1.9312D+00,1.9081D+00,
     &1.8860D+00,1.8635D+00,1.8406D+00,1.8221D+00,1.8007D+00,1.7764D+00,
     &1.7489D+00,1.7195D+00,1.6855D+00,1.6600D+00,1.6332D+00,1.6031D+00,
     &1.5760D+00,1.5532D+00,1.5397D+00,1.5376D+00,1.5507D+00,1.5929D+00,
     &1.6784D+00,1.7759D+00,1.9129D+00,2.0531D+00,2.2292D+00,2.4032D+00/
      DATA (XUVF_L(K),K=  115,  228) /
     &2.5752D+00,2.7449D+00,2.9135D+00,3.0810D+00,3.2491D+00,3.4183D+00,
     &3.5898D+00,3.7650D+00,3.9437D+00,4.0443D+00,4.3402D+00,4.6262D+00,
     &4.9009D+00,5.1640D+00,5.4156D+00,5.6530D+00,5.8759D+00,6.0779D+00,
     &6.2540D+00,6.3836D+00,6.4062D+00,0.0000D+00,2.4808D+00,2.4513D+00,
     &2.4236D+00,2.3948D+00,2.3680D+00,2.3397D+00,2.3127D+00,2.2853D+00,
     &2.2585D+00,2.2307D+00,2.2026D+00,2.1762D+00,2.1490D+00,2.1332D+00,
     &2.1164D+00,2.0964D+00,2.0766D+00,2.0565D+00,2.0353D+00,2.0171D+00,
     &1.9969D+00,1.9738D+00,1.9501D+00,1.9258D+00,1.9026D+00,1.8821D+00,
     &1.8594D+00,1.8330D+00,1.8046D+00,1.7734D+00,1.7378D+00,1.7112D+00,
     &1.6829D+00,1.6514D+00,1.6228D+00,1.5994D+00,1.5840D+00,1.5808D+00,
     &1.5927D+00,1.6334D+00,1.7157D+00,1.8093D+00,1.9406D+00,2.0735D+00,
     &2.2394D+00,2.4019D+00,2.5615D+00,2.7178D+00,2.8718D+00,3.0246D+00,
     &3.1766D+00,3.3284D+00,3.4820D+00,3.6370D+00,3.7952D+00,3.8716D+00,
     &4.1225D+00,4.3580D+00,4.5798D+00,4.7847D+00,4.9730D+00,5.1395D+00,
     &5.2832D+00,5.3945D+00,5.4634D+00,5.4612D+00,5.2940D+00,0.0000D+00,
     &2.5823D+00,2.5527D+00,2.5226D+00,2.4928D+00,2.4650D+00,2.4358D+00,
     &2.4071D+00,2.3783D+00,2.3505D+00,2.3212D+00,2.2928D+00,2.2636D+00,
     &2.2360D+00,2.2185D+00,2.2005D+00,2.1801D+00,2.1591D+00,2.1376D+00,
     &2.1153D+00,2.0960D+00,2.0747D+00,2.0505D+00,2.0247D+00,1.9991D+00/
      DATA (XUVF_L(K),K=  229,  342) /
     &1.9746D+00,1.9523D+00,1.9287D+00,1.9000D+00,1.8693D+00,1.8361D+00,
     &1.7994D+00,1.7711D+00,1.7409D+00,1.7076D+00,1.6772D+00,1.6517D+00,
     &1.6345D+00,1.6302D+00,1.6408D+00,1.6789D+00,1.7574D+00,1.8457D+00,
     &1.9692D+00,2.0939D+00,2.2474D+00,2.3969D+00,2.5419D+00,2.6837D+00,
     &2.8216D+00,2.9573D+00,3.0915D+00,3.2246D+00,3.3583D+00,3.4917D+00,
     &3.6273D+00,3.6791D+00,3.8823D+00,4.0673D+00,4.2350D+00,4.3813D+00,
     &4.5072D+00,4.6083D+00,4.6757D+00,4.7055D+00,4.6825D+00,4.5674D+00,
     &4.2566D+00,0.0000D+00,2.7025D+00,2.6705D+00,2.6393D+00,2.6093D+00,
     &2.5790D+00,2.5484D+00,2.5184D+00,2.4880D+00,2.4590D+00,2.4277D+00,
     &2.3971D+00,2.3669D+00,2.3380D+00,2.3200D+00,2.3002D+00,2.2782D+00,
     &2.2557D+00,2.2331D+00,2.2092D+00,2.1887D+00,2.1660D+00,2.1400D+00,
     &2.1126D+00,2.0859D+00,2.0586D+00,2.0351D+00,2.0094D+00,1.9786D+00,
     &1.9453D+00,1.9096D+00,1.8707D+00,1.8406D+00,1.8084D+00,1.7728D+00,
     &1.7392D+00,1.7128D+00,1.6933D+00,1.6875D+00,1.6949D+00,1.7295D+00,
     &1.8023D+00,1.8845D+00,1.9991D+00,2.1134D+00,2.2525D+00,2.3868D+00,
     &2.5160D+00,2.6405D+00,2.7609D+00,2.8781D+00,2.9929D+00,3.1059D+00,
     &3.2180D+00,3.3292D+00,3.4407D+00,3.4675D+00,3.6225D+00,3.7573D+00,
     &3.8710D+00,3.9617D+00,4.0270D+00,4.0642D+00,4.0675D+00,4.0263D+00,
     &3.9240D+00,3.7262D+00,3.3217D+00,0.0000D+00,2.8135D+00,2.7813D+00/
      DATA (XUVF_L(K),K=  343,  456) /
     &2.7489D+00,2.7166D+00,2.6850D+00,2.6527D+00,2.6212D+00,2.5898D+00,
     &2.5592D+00,2.5267D+00,2.4943D+00,2.4636D+00,2.4320D+00,2.4129D+00,
     &2.3929D+00,2.3695D+00,2.3453D+00,2.3211D+00,2.2959D+00,2.2740D+00,
     &2.2496D+00,2.2221D+00,2.1931D+00,2.1653D+00,2.1356D+00,2.1112D+00,
     &2.0830D+00,2.0503D+00,2.0147D+00,1.9766D+00,1.9361D+00,1.9037D+00,
     &1.8696D+00,1.8318D+00,1.7966D+00,1.7677D+00,1.7459D+00,1.7378D+00,
     &1.7430D+00,1.7738D+00,1.8407D+00,1.9169D+00,2.0223D+00,2.1273D+00,
     &2.2537D+00,2.3742D+00,2.4892D+00,2.5990D+00,2.7043D+00,2.8056D+00,
     &2.9038D+00,3.0000D+00,3.0936D+00,3.1864D+00,3.2782D+00,3.2867D+00,
     &3.4021D+00,3.4971D+00,3.5691D+00,3.6188D+00,3.6422D+00,3.6335D+00,
     &3.5908D+00,3.5036D+00,3.3552D+00,3.1085D+00,2.6634D+00,0.0000D+00,
     &2.9406D+00,2.9062D+00,2.8726D+00,2.8385D+00,2.8060D+00,2.7720D+00,
     &2.7392D+00,2.7058D+00,2.6734D+00,2.6399D+00,2.6057D+00,2.5722D+00,
     &2.5390D+00,2.5194D+00,2.4975D+00,2.4728D+00,2.4471D+00,2.4216D+00,
     &2.3945D+00,2.3712D+00,2.3458D+00,2.3152D+00,2.2856D+00,2.2545D+00,
     &2.2237D+00,2.1966D+00,2.1672D+00,2.1312D+00,2.0926D+00,2.0521D+00,
     &2.0093D+00,1.9748D+00,1.9384D+00,1.8975D+00,1.8601D+00,1.8275D+00,
     &1.8036D+00,1.7924D+00,1.7948D+00,1.8206D+00,1.8808D+00,1.9499D+00,
     &2.0450D+00,2.1390D+00,2.2512D+00,2.3570D+00,2.4564D+00,2.5501D+00/
      DATA (XUVF_L(K),K=  457,  570) /
     &2.6391D+00,2.7240D+00,2.8053D+00,2.8834D+00,2.9590D+00,3.0326D+00,
     &3.1042D+00,3.0942D+00,3.1727D+00,3.2289D+00,3.2628D+00,3.2739D+00,
     &3.2574D+00,3.2103D+00,3.1297D+00,3.0047D+00,2.8211D+00,2.5467D+00,
     &2.0897D+00,0.0000D+00,3.0557D+00,3.0193D+00,2.9840D+00,2.9497D+00,
     &2.9150D+00,2.8801D+00,2.8454D+00,2.8109D+00,2.7771D+00,2.7412D+00,
     &2.7065D+00,2.6716D+00,2.6360D+00,2.6149D+00,2.5923D+00,2.5663D+00,
     &2.5395D+00,2.5120D+00,2.4834D+00,2.4589D+00,2.4330D+00,2.4011D+00,
     &2.3676D+00,2.3363D+00,2.3027D+00,2.2736D+00,2.2422D+00,2.2040D+00,
     &2.1629D+00,2.1194D+00,2.0750D+00,2.0384D+00,1.9996D+00,1.9565D+00,
     &1.9160D+00,1.8811D+00,1.8541D+00,1.8409D+00,1.8399D+00,1.8611D+00,
     &1.9143D+00,1.9764D+00,2.0622D+00,2.1459D+00,2.2457D+00,2.3385D+00,
     &2.4249D+00,2.5051D+00,2.5806D+00,2.6515D+00,2.7182D+00,2.7823D+00,
     &2.8427D+00,2.9008D+00,2.9564D+00,2.9332D+00,2.9828D+00,3.0094D+00,
     &3.0142D+00,2.9955D+00,2.9537D+00,2.8796D+00,2.7735D+00,2.6260D+00,
     &2.4242D+00,2.1388D+00,1.6900D+00,0.0000D+00,3.1718D+00,3.1348D+00,
     &3.0971D+00,3.0610D+00,3.0260D+00,2.9896D+00,2.9533D+00,2.9173D+00,
     &2.8818D+00,2.8449D+00,2.8072D+00,2.7709D+00,2.7340D+00,2.7121D+00,
     &2.6877D+00,2.6605D+00,2.6319D+00,2.6032D+00,2.5732D+00,2.5471D+00,
     &2.5180D+00,2.4851D+00,2.4511D+00,2.4170D+00,2.3817D+00,2.3505D+00/
      DATA (XUVF_L(K),K=  571,  684) /
     &2.3172D+00,2.2762D+00,2.2328D+00,2.1868D+00,2.1400D+00,2.1012D+00,
     &2.0601D+00,2.0136D+00,1.9704D+00,1.9335D+00,1.9035D+00,1.8868D+00,
     &1.8827D+00,1.8990D+00,1.9452D+00,2.0005D+00,2.0763D+00,2.1507D+00,
     &2.2377D+00,2.3179D+00,2.3917D+00,2.4592D+00,2.5218D+00,2.5799D+00,
     &2.6336D+00,2.6843D+00,2.7314D+00,2.7753D+00,2.8166D+00,2.7824D+00,
     &2.8054D+00,2.8081D+00,2.7893D+00,2.7474D+00,2.6818D+00,2.5888D+00,
     &2.4646D+00,2.3032D+00,2.0902D+00,1.8025D+00,1.3740D+00,0.0000D+00,
     &3.2793D+00,3.2385D+00,3.2014D+00,3.1643D+00,3.1270D+00,3.0888D+00,
     &3.0517D+00,3.0141D+00,2.9773D+00,2.9392D+00,2.9009D+00,2.8610D+00,
     &2.8230D+00,2.8000D+00,2.7754D+00,2.7459D+00,2.7163D+00,2.6858D+00,
     &2.6545D+00,2.6270D+00,2.5962D+00,2.5617D+00,2.5271D+00,2.4903D+00,
     &2.4527D+00,2.4207D+00,2.3851D+00,2.3421D+00,2.2960D+00,2.2476D+00,
     &2.1987D+00,2.1578D+00,2.1146D+00,2.0670D+00,2.0202D+00,1.9796D+00,
     &1.9468D+00,1.9282D+00,1.9203D+00,1.9319D+00,1.9712D+00,2.0197D+00,
     &2.0872D+00,2.1524D+00,2.2288D+00,2.2981D+00,2.3606D+00,2.4177D+00,
     &2.4692D+00,2.5159D+00,2.5591D+00,2.5981D+00,2.6339D+00,2.6669D+00,
     &2.6962D+00,2.6528D+00,2.6566D+00,2.6395D+00,2.6028D+00,2.5437D+00,
     &2.4622D+00,2.3555D+00,2.2200D+00,2.0488D+00,1.8335D+00,1.5506D+00,
     &1.1442D+00,0.0000D+00,3.3868D+00,3.3470D+00,3.3075D+00,3.2689D+00/
      DATA (XUVF_L(K),K=  685,  798) /
     &3.2300D+00,3.1909D+00,3.1517D+00,3.1129D+00,3.0747D+00,3.0335D+00,
     &2.9946D+00,2.9537D+00,2.9140D+00,2.8896D+00,2.8638D+00,2.8337D+00,
     &2.8021D+00,2.7705D+00,2.7373D+00,2.7075D+00,2.6767D+00,2.6403D+00,
     &2.6031D+00,2.5649D+00,2.5258D+00,2.4917D+00,2.4537D+00,2.4080D+00,
     &2.3597D+00,2.3091D+00,2.2580D+00,2.2150D+00,2.1692D+00,2.1186D+00,
     &2.0701D+00,2.0257D+00,1.9901D+00,1.9679D+00,1.9571D+00,1.9629D+00,
     &1.9955D+00,2.0378D+00,2.0963D+00,2.1529D+00,2.2178D+00,2.2766D+00,
     &2.3287D+00,2.3749D+00,2.4162D+00,2.4529D+00,2.4850D+00,2.5140D+00,
     &2.5392D+00,2.5617D+00,2.5798D+00,2.5298D+00,2.5151D+00,2.4811D+00,
     &2.4282D+00,2.3561D+00,2.2611D+00,2.1439D+00,2.0005D+00,1.8252D+00,
     &1.6091D+00,1.3345D+00,9.5375D-01,0.0000D+00,3.4912D+00,3.4507D+00,
     &3.4100D+00,3.3696D+00,3.3310D+00,3.2893D+00,3.2496D+00,3.2088D+00,
     &3.1686D+00,3.1278D+00,3.0865D+00,3.0438D+00,3.0020D+00,2.9766D+00,
     &2.9494D+00,2.9180D+00,2.8850D+00,2.8520D+00,2.8174D+00,2.7877D+00,
     &2.7550D+00,2.7169D+00,2.6762D+00,2.6369D+00,2.5958D+00,2.5594D+00,
     &2.5195D+00,2.4721D+00,2.4211D+00,2.3680D+00,2.3145D+00,2.2695D+00,
     &2.2214D+00,2.1684D+00,2.1154D+00,2.0706D+00,2.0303D+00,2.0058D+00,
     &1.9909D+00,1.9920D+00,2.0177D+00,2.0531D+00,2.1031D+00,2.1511D+00,
     &2.2060D+00,2.2548D+00,2.2972D+00,2.3339D+00,2.3655D+00,2.3927D+00/
      DATA (XUVF_L(K),K=  799,  912) /
     &2.4159D+00,2.4357D+00,2.4520D+00,2.4644D+00,2.4735D+00,2.4171D+00,
     &2.3878D+00,2.3397D+00,2.2743D+00,2.1907D+00,2.0861D+00,1.9611D+00,
     &1.8128D+00,1.6351D+00,1.4227D+00,1.1584D+00,8.0371D-01,0.0000D+00,
     &3.5892D+00,3.5473D+00,3.5055D+00,3.4637D+00,3.4230D+00,3.3809D+00,
     &3.3396D+00,3.2976D+00,3.2571D+00,3.2126D+00,3.1696D+00,3.1272D+00,
     &3.0840D+00,3.0569D+00,3.0286D+00,2.9959D+00,2.9619D+00,2.9273D+00,
     &2.8910D+00,2.8598D+00,2.8266D+00,2.7863D+00,2.7448D+00,2.7029D+00,
     &2.6598D+00,2.6219D+00,2.5804D+00,2.5305D+00,2.4773D+00,2.4214D+00,
     &2.3662D+00,2.3191D+00,2.2698D+00,2.2126D+00,2.1577D+00,2.1092D+00,
     &2.0674D+00,2.0393D+00,2.0210D+00,2.0173D+00,2.0367D+00,2.0654D+00,
     &2.1076D+00,2.1485D+00,2.1942D+00,2.2338D+00,2.2678D+00,2.2959D+00,
     &2.3193D+00,2.3386D+00,2.3539D+00,2.3660D+00,2.3738D+00,2.3789D+00,
     &2.3799D+00,2.3197D+00,2.2776D+00,2.2186D+00,2.1426D+00,2.0495D+00,
     &1.9397D+00,1.8097D+00,1.6583D+00,1.4814D+00,1.2736D+00,1.0200D+00,
     &6.8880D-01,0.0000D+00,3.7157D+00,3.6699D+00,3.6275D+00,3.5842D+00,
     &3.5420D+00,3.4972D+00,3.4542D+00,3.4107D+00,3.3678D+00,3.3234D+00,
     &3.2774D+00,3.2332D+00,3.1870D+00,3.1600D+00,3.1297D+00,3.0952D+00,
     &3.0595D+00,3.0231D+00,2.9850D+00,2.9534D+00,2.9160D+00,2.8740D+00,
     &2.8312D+00,2.7872D+00,2.7408D+00,2.7014D+00,2.6568D+00,2.6045D+00/
      DATA (XUVF_L(K),K=  913, 1026) /
     &2.5481D+00,2.4895D+00,2.4315D+00,2.3817D+00,2.3283D+00,2.2697D+00,
     &2.2106D+00,2.1591D+00,2.1128D+00,2.0807D+00,2.0578D+00,2.0477D+00,
     &2.0583D+00,2.0796D+00,2.1122D+00,2.1433D+00,2.1777D+00,2.2069D+00,
     &2.2299D+00,2.2483D+00,2.2618D+00,2.2718D+00,2.2778D+00,2.2803D+00,
     &2.2797D+00,2.2749D+00,2.2668D+00,2.2019D+00,2.1468D+00,2.0761D+00,
     &1.9902D+00,1.8883D+00,1.7711D+00,1.6370D+00,1.4847D+00,1.3103D+00,
     &1.1091D+00,8.7047D-01,5.6856D-01,0.0000D+00,3.8327D+00,3.7877D+00,
     &3.7424D+00,3.6981D+00,3.6540D+00,3.6083D+00,3.5637D+00,3.5184D+00,
     &3.4753D+00,3.4271D+00,3.3800D+00,3.3325D+00,3.2860D+00,3.2564D+00,
     &3.2258D+00,3.1893D+00,3.1519D+00,3.1135D+00,3.0738D+00,3.0389D+00,
     &3.0010D+00,2.9580D+00,2.9118D+00,2.8654D+00,2.8178D+00,2.7758D+00,
     &2.7289D+00,2.6738D+00,2.6146D+00,2.5530D+00,2.4924D+00,2.4399D+00,
     &2.3845D+00,2.3213D+00,2.2605D+00,2.2040D+00,2.1540D+00,2.1186D+00,
     &2.0908D+00,2.0749D+00,2.0772D+00,2.0914D+00,2.1145D+00,2.1368D+00,
     &2.1613D+00,2.1804D+00,2.1941D+00,2.2037D+00,2.2088D+00,2.2101D+00,
     &2.2083D+00,2.2031D+00,2.1942D+00,2.1826D+00,2.1665D+00,2.0987D+00,
     &2.0321D+00,1.9516D+00,1.8571D+00,1.7497D+00,1.6281D+00,1.4923D+00,
     &1.3406D+00,1.1697D+00,9.7635D-01,7.5209D-01,4.7638D-01,0.0000D+00,
     &3.9497D+00,3.9009D+00,3.8555D+00,3.8080D+00,3.7630D+00,3.7163D+00/
      DATA (XUVF_L(K),K= 1027, 1140) /
     &3.6699D+00,3.6231D+00,3.5765D+00,3.5285D+00,3.4807D+00,3.4305D+00,
     &3.3810D+00,3.3511D+00,3.3185D+00,3.2805D+00,3.2414D+00,3.2016D+00,
     &3.1598D+00,3.1244D+00,3.0837D+00,3.0383D+00,2.9908D+00,2.9424D+00,
     &2.8919D+00,2.8477D+00,2.7990D+00,2.7403D+00,2.6784D+00,2.6142D+00,
     &2.5507D+00,2.4960D+00,2.4362D+00,2.3710D+00,2.3058D+00,2.2463D+00,
     &2.1931D+00,2.1539D+00,2.1216D+00,2.0996D+00,2.0940D+00,2.1012D+00,
     &2.1154D+00,2.1294D+00,2.1444D+00,2.1543D+00,2.1597D+00,2.1610D+00,
     &2.1585D+00,2.1523D+00,2.1432D+00,2.1307D+00,2.1155D+00,2.0964D+00,
     &2.0742D+00,2.0035D+00,1.9273D+00,1.8396D+00,1.7387D+00,1.6273D+00,
     &1.5032D+00,1.3665D+00,1.2164D+00,1.0501D+00,8.6515D-01,6.5470D-01,
     &4.0284D-01,0.0000D+00,4.0572D+00,4.0093D+00,3.9616D+00,3.9140D+00,
     &3.8670D+00,3.8185D+00,3.7706D+00,3.7224D+00,3.6746D+00,3.6251D+00,
     &3.5744D+00,3.5233D+00,3.4720D+00,3.4406D+00,3.4062D+00,3.3671D+00,
     &3.3263D+00,3.2847D+00,3.2414D+00,3.2046D+00,3.1620D+00,3.1150D+00,
     &3.0653D+00,3.0145D+00,2.9619D+00,2.9153D+00,2.8641D+00,2.8032D+00,
     &2.7388D+00,2.6715D+00,2.6056D+00,2.5481D+00,2.4880D+00,2.4171D+00,
     &2.3496D+00,2.2862D+00,2.2282D+00,2.1865D+00,2.1502D+00,2.1217D+00,
     &2.1086D+00,2.1086D+00,2.1149D+00,2.1216D+00,2.1275D+00,2.1295D+00,
     &2.1273D+00,2.1212D+00,2.1119D+00,2.0992D+00,2.0837D+00,2.0653D+00/
      DATA (XUVF_L(K),K= 1141, 1254) /
     &2.0442D+00,2.0194D+00,1.9912D+00,1.9193D+00,1.8359D+00,1.7412D+00,
     &1.6366D+00,1.5214D+00,1.3956D+00,1.2594D+00,1.1115D+00,9.5033D-01,
     &7.7356D-01,5.7585D-01,3.4506D-01,0.0000D+00,4.1710D+00,4.1201D+00,
     &4.0712D+00,4.0213D+00,3.9730D+00,3.9228D+00,3.8734D+00,3.8233D+00,
     &3.7726D+00,3.7217D+00,3.6699D+00,3.6160D+00,3.5640D+00,3.5311D+00,
     &3.4960D+00,3.4549D+00,3.4121D+00,3.3689D+00,3.3237D+00,3.2848D+00,
     &3.2425D+00,3.1917D+00,3.1399D+00,3.0866D+00,3.0319D+00,2.9838D+00,
     &2.9306D+00,2.8668D+00,2.7992D+00,2.7291D+00,2.6605D+00,2.6007D+00,
     &2.5375D+00,2.4631D+00,2.3919D+00,2.3261D+00,2.2643D+00,2.2183D+00,
     &2.1772D+00,2.1426D+00,2.1222D+00,2.1155D+00,2.1135D+00,2.1130D+00,
     &2.1102D+00,2.1039D+00,2.0941D+00,2.0815D+00,2.0652D+00,2.0466D+00,
     &2.0251D+00,2.0014D+00,1.9746D+00,1.9450D+00,1.9116D+00,1.8381D+00,
     &1.7481D+00,1.6484D+00,1.5404D+00,1.4225D+00,1.2963D+00,1.1611D+00,
     &1.0161D+00,8.6047D-01,6.9193D-01,5.0691D-01,2.9581D-01,0.0000D+00,
     &4.2754D+00,4.2238D+00,4.1737D+00,4.1233D+00,4.0740D+00,4.0219D+00,
     &3.9713D+00,3.9196D+00,3.8675D+00,3.8160D+00,3.7618D+00,3.7060D+00,
     &3.6510D+00,3.6173D+00,3.5808D+00,3.5380D+00,3.4941D+00,3.4493D+00,
     &3.4027D+00,3.3623D+00,3.3163D+00,3.2647D+00,3.2114D+00,3.1563D+00,
     &3.0989D+00,3.0489D+00,2.9929D+00,2.9263D+00,2.8563D+00,2.7837D+00/
      DATA (XUVF_L(K),K= 1255, 1368) /
     &2.7122D+00,2.6501D+00,2.5825D+00,2.5073D+00,2.4327D+00,2.3623D+00,
     &2.2962D+00,2.2474D+00,2.2020D+00,2.1616D+00,2.1335D+00,2.1209D+00,
     &2.1113D+00,2.1034D+00,2.0929D+00,2.0795D+00,2.0634D+00,2.0439D+00,
     &2.0222D+00,1.9982D+00,1.9716D+00,1.9428D+00,1.9113D+00,1.8773D+00,
     &1.8394D+00,1.7649D+00,1.6692D+00,1.5658D+00,1.4547D+00,1.3360D+00,
     &1.2095D+00,1.0761D+00,9.3485D-01,7.8430D-01,6.2380D-01,4.5010D-01,
     &2.5625D-01,0.0000D+00,4.3798D+00,4.3275D+00,4.2762D+00,4.2239D+00,
     &4.1730D+00,4.1196D+00,4.0674D+00,4.0143D+00,3.9623D+00,3.9056D+00,
     &3.8502D+00,3.7935D+00,3.7370D+00,3.7018D+00,3.6642D+00,3.6200D+00,
     &3.5742D+00,3.5277D+00,3.4786D+00,3.4371D+00,3.3901D+00,3.3359D+00,
     &3.2800D+00,3.2235D+00,3.1639D+00,3.1115D+00,3.0537D+00,2.9847D+00,
     &2.9116D+00,2.8364D+00,2.7623D+00,2.6973D+00,2.6275D+00,2.5497D+00,
     &2.4705D+00,2.3972D+00,2.3281D+00,2.2747D+00,2.2253D+00,2.1793D+00,
     &2.1444D+00,2.1253D+00,2.1081D+00,2.0939D+00,2.0755D+00,2.0555D+00,
     &2.0332D+00,2.0081D+00,1.9814D+00,1.9522D+00,1.9205D+00,1.8875D+00,
     &1.8520D+00,1.8139D+00,1.7725D+00,1.6968D+00,1.5976D+00,1.4911D+00,
     &1.3772D+00,1.2577D+00,1.1320D+00,1.0005D+00,8.6242D-01,7.1750D-01,
     &5.6466D-01,4.0150D-01,2.2333D-01,0.0000D+00,4.4809D+00,4.4265D+00,
     &4.3735D+00,4.3193D+00,4.2670D+00,4.2128D+00,4.1585D+00,4.1039D+00/
      DATA (XUVF_L(K),K= 1369, 1482) /
     &4.0509D+00,3.9928D+00,3.9351D+00,3.8769D+00,3.8180D+00,3.7821D+00,
     &3.7434D+00,3.6974D+00,3.6501D+00,3.6019D+00,3.5513D+00,3.5093D+00,
     &3.4594D+00,3.4035D+00,3.3456D+00,3.2870D+00,3.2250D+00,3.1715D+00,
     &3.1110D+00,3.0396D+00,2.9639D+00,2.8863D+00,2.8096D+00,2.7429D+00,
     &2.6702D+00,2.5884D+00,2.5068D+00,2.4296D+00,2.3560D+00,2.3003D+00,
     &2.2464D+00,2.1951D+00,2.1530D+00,2.1283D+00,2.1045D+00,2.0843D+00,
     &2.0591D+00,2.0328D+00,2.0047D+00,1.9749D+00,1.9429D+00,1.9096D+00,
     &1.8740D+00,1.8369D+00,1.7978D+00,1.7560D+00,1.7116D+00,1.6360D+00,
     &1.5322D+00,1.4233D+00,1.3084D+00,1.1885D+00,1.0637D+00,9.3449D-01,
     &7.9961D-01,6.6020D-01,5.1453D-01,3.6103D-01,1.9641D-01,0.0000D+00,
     &4.6169D+00,4.5608D+00,4.5060D+00,4.4504D+00,4.3960D+00,4.3395D+00,
     &4.2837D+00,4.2262D+00,4.1710D+00,4.1106D+00,4.0517D+00,3.9908D+00,
     &3.9300D+00,3.8920D+00,3.8509D+00,3.8030D+00,3.7538D+00,3.7035D+00,
     &3.6494D+00,3.6055D+00,3.5556D+00,3.4966D+00,3.4351D+00,3.3738D+00,
     &3.3090D+00,3.2518D+00,3.1888D+00,3.1141D+00,3.0348D+00,2.9533D+00,
     &2.8730D+00,2.8020D+00,2.7264D+00,2.6400D+00,2.5551D+00,2.4732D+00,
     &2.3941D+00,2.3329D+00,2.2742D+00,2.2147D+00,2.1644D+00,2.1317D+00,
     &2.0986D+00,2.0700D+00,2.0363D+00,2.0021D+00,1.9668D+00,1.9299D+00,
     &1.8922D+00,1.8532D+00,1.8125D+00,1.7704D+00,1.7270D+00,1.6809D+00/
      DATA (XUVF_L(K),K= 1483, 1596) /
     &1.6327D+00,1.5570D+00,1.4497D+00,1.3373D+00,1.2215D+00,1.1020D+00,
     &9.7897D-01,8.5304D-01,7.2349D-01,5.9074D-01,4.5411D-01,3.1307D-01,
     &1.6547D-01,0.0000D+00,4.7403D+00,4.6834D+00,4.6262D+00,4.5696D+00,
     &4.5140D+00,4.4557D+00,4.3978D+00,4.3393D+00,4.2817D+00,4.2191D+00,
     &4.1578D+00,4.0941D+00,4.0310D+00,3.9917D+00,3.9492D+00,3.8995D+00,
     &3.8481D+00,3.7958D+00,3.7411D+00,3.6937D+00,3.6405D+00,3.5806D+00,
     &3.5171D+00,3.4520D+00,3.3840D+00,3.3254D+00,3.2596D+00,3.1812D+00,
     &3.0985D+00,3.0137D+00,2.9301D+00,2.8556D+00,2.7782D+00,2.6879D+00,
     &2.5974D+00,2.5119D+00,2.4281D+00,2.3629D+00,2.2982D+00,2.2324D+00,
     &2.1730D+00,2.1332D+00,2.0922D+00,2.0570D+00,2.0152D+00,1.9739D+00,
     &1.9323D+00,1.8902D+00,1.8474D+00,1.8039D+00,1.7589D+00,1.7129D+00,
     &1.6654D+00,1.6163D+00,1.5652D+00,1.4896D+00,1.3789D+00,1.2649D+00,
     &1.1487D+00,1.0300D+00,9.0896D-01,7.8619D-01,6.6149D-01,5.3498D-01,
     &4.0654D-01,2.7586D-01,1.4208D-01,0.0000D+00,4.8699D+00,4.8107D+00,
     &4.7518D+00,4.6928D+00,4.6350D+00,4.5750D+00,4.5152D+00,4.4524D+00,
     &4.3956D+00,4.3299D+00,4.2674D+00,4.2014D+00,4.1350D+00,4.0939D+00,
     &4.0503D+00,3.9982D+00,3.9448D+00,3.8905D+00,3.8328D+00,3.7846D+00,
     &3.7300D+00,3.6664D+00,3.5991D+00,3.5326D+00,3.4620D+00,3.3998D+00,
     &3.3311D+00,3.2494D+00,3.1632D+00,3.0752D+00,2.9881D+00,2.9120D+00/
      DATA (XUVF_L(K),K= 1597, 1710) /
     &2.8299D+00,2.7339D+00,2.6398D+00,2.5493D+00,2.4611D+00,2.3911D+00,
     &2.3215D+00,2.2482D+00,2.1812D+00,2.1342D+00,2.0854D+00,2.0427D+00,
     &1.9932D+00,1.9453D+00,1.8978D+00,1.8504D+00,1.8030D+00,1.7545D+00,
     &1.7059D+00,1.6565D+00,1.6056D+00,1.5535D+00,1.4989D+00,1.4245D+00,
     &1.3108D+00,1.1959D+00,1.0798D+00,9.6219D-01,8.4358D-01,7.2422D-01,
     &6.0451D-01,4.8425D-01,3.6380D-01,2.4286D-01,1.2189D-01,0.0000D+00,
     &4.9964D+00,4.9356D+00,4.8755D+00,4.8147D+00,4.7550D+00,4.6935D+00,
     &4.6315D+00,4.5697D+00,4.5062D+00,4.4406D+00,4.3752D+00,4.3061D+00,
     &4.2380D+00,4.1962D+00,4.1500D+00,4.0963D+00,4.0405D+00,3.9832D+00,
     &3.9245D+00,3.8728D+00,3.8172D+00,3.7504D+00,3.6811D+00,3.6108D+00,
     &3.5381D+00,3.4734D+00,3.4018D+00,3.3164D+00,3.2269D+00,3.1352D+00,
     &3.0446D+00,2.9657D+00,2.8794D+00,2.7800D+00,2.6821D+00,2.5867D+00,
     &2.4930D+00,2.4184D+00,2.3433D+00,2.2634D+00,2.1877D+00,2.1342D+00,
     &2.0772D+00,2.0279D+00,1.9713D+00,1.9172D+00,1.8642D+00,1.8120D+00,
     &1.7600D+00,1.7076D+00,1.6553D+00,1.6027D+00,1.5491D+00,1.4938D+00,
     &1.4374D+00,1.3637D+00,1.2481D+00,1.1325D+00,1.0166D+00,9.0047D-01,
     &7.8428D-01,6.6889D-01,5.5381D-01,4.3953D-01,3.2652D-01,2.1461D-01,
     &1.0498D-01,0.0000D+00,5.1134D+00,5.0511D+00,4.9886D+00,4.9273D+00,
     &4.8660D+00,4.8016D+00,4.7382D+00,4.6744D+00,4.6106D+00,4.5420D+00/
      DATA (XUVF_L(K),K= 1711, 1824) /
     &4.4742D+00,4.4028D+00,4.3320D+00,4.2892D+00,4.2413D+00,4.1858D+00,
     &4.1281D+00,4.0682D+00,4.0067D+00,3.9556D+00,3.8955D+00,3.8271D+00,
     &3.7556D+00,3.6829D+00,3.6071D+00,3.5401D+00,3.4662D+00,3.3777D+00,
     &3.2849D+00,3.1898D+00,3.0960D+00,3.0140D+00,2.9244D+00,2.8224D+00,
     &2.7183D+00,2.6191D+00,2.5219D+00,2.4431D+00,2.3628D+00,2.2767D+00,
     &2.1931D+00,2.1332D+00,2.0695D+00,2.0145D+00,1.9514D+00,1.8920D+00,
     &1.8340D+00,1.7775D+00,1.7215D+00,1.6664D+00,1.6108D+00,1.5553D+00,
     &1.4995D+00,1.4421D+00,1.3839D+00,1.3103D+00,1.1944D+00,1.0782D+00,
     &9.6271D-01,8.4822D-01,7.3481D-01,6.2240D-01,5.1184D-01,4.0291D-01,
     &2.9618D-01,1.9206D-01,9.1846D-02,0.0000D+00,5.2367D+00,5.1713D+00,
     &5.1071D+00,5.0425D+00,4.9800D+00,4.9141D+00,4.8489D+00,4.7833D+00,
     &4.7181D+00,4.6457D+00,4.5768D+00,4.5034D+00,4.4300D+00,4.3847D+00,
     &4.3353D+00,4.2782D+00,4.2182D+00,4.1570D+00,4.0921D+00,4.0385D+00,
     &3.9782D+00,3.9074D+00,3.8331D+00,3.7575D+00,3.6781D+00,3.6086D+00,
     &3.5313D+00,3.4401D+00,3.3439D+00,3.2455D+00,3.1483D+00,3.0623D+00,
     &2.9694D+00,2.8629D+00,2.7561D+00,2.6527D+00,2.5508D+00,2.4669D+00,
     &2.3816D+00,2.2887D+00,2.1979D+00,2.1317D+00,2.0613D+00,2.0002D+00,
     &1.9307D+00,1.8659D+00,1.8033D+00,1.7426D+00,1.6834D+00,1.6247D+00,
     &1.5668D+00,1.5085D+00,1.4504D+00,1.3916D+00,1.3311D+00,1.2591D+00/
      DATA (XUVF_L(K),K= 1825, 1836) /
     &1.1415D+00,1.0256D+00,9.1107D-01,7.9840D-01,6.8736D-01,5.7902D-01,
     &4.7260D-01,3.6895D-01,2.6838D-01,1.7161D-01,8.0264D-02,0.0000D+00/
      DATA (XDVF_L(K),K=    1,  114) /
     &1.4230D+00,1.4064D+00,1.3903D+00,1.3749D+00,1.3590D+00,1.3424D+00,
     &1.3271D+00,1.3114D+00,1.2962D+00,1.2803D+00,1.2647D+00,1.2492D+00,
     &1.2340D+00,1.2246D+00,1.2155D+00,1.2044D+00,1.1927D+00,1.1814D+00,
     &1.1695D+00,1.1589D+00,1.1479D+00,1.1347D+00,1.1214D+00,1.1080D+00,
     &1.0944D+00,1.0824D+00,1.0700D+00,1.0544D+00,1.0371D+00,1.0188D+00,
     &9.9884D-01,9.8287D-01,9.6563D-01,9.4645D-01,9.2847D-01,9.1313D-01,
     &9.0246D-01,8.9955D-01,9.0461D-01,9.2737D-01,9.7648D-01,1.0343D+00,
     &1.1168D+00,1.2030D+00,1.3129D+00,1.4240D+00,1.5357D+00,1.6492D+00,
     &1.7643D+00,1.8818D+00,2.0016D+00,2.1253D+00,2.2535D+00,2.3853D+00,
     &2.5225D+00,2.5620D+00,2.7906D+00,3.0230D+00,3.2574D+00,3.4983D+00,
     &3.7459D+00,4.0062D+00,4.2803D+00,4.5790D+00,4.9150D+00,5.3263D+00,
     &5.9228D+00,0.0000D+00,1.4698D+00,1.4526D+00,1.4360D+00,1.4199D+00,
     &1.4030D+00,1.3864D+00,1.3702D+00,1.3542D+00,1.3386D+00,1.3221D+00,
     &1.3059D+00,1.2896D+00,1.2740D+00,1.2644D+00,1.2544D+00,1.2425D+00,
     &1.2309D+00,1.2185D+00,1.2061D+00,1.1953D+00,1.1836D+00,1.1697D+00,
     &1.1558D+00,1.1417D+00,1.1275D+00,1.1154D+00,1.1011D+00,1.0844D+00,
     &1.0663D+00,1.0471D+00,1.0261D+00,1.0092D+00,9.9133D-01,9.7103D-01,
     &9.5184D-01,9.3560D-01,9.2380D-01,9.1922D-01,9.2378D-01,9.4563D-01,
     &9.9235D-01,1.0474D+00,1.1262D+00,1.2078D+00,1.3110D+00,1.4146D+00/
      DATA (XDVF_L(K),K=  115,  228) /
     &1.5192D+00,1.6241D+00,1.7298D+00,1.8375D+00,1.9471D+00,2.0592D+00,
     &2.1741D+00,2.2925D+00,2.4144D+00,2.4425D+00,2.6407D+00,2.8375D+00,
     &3.0361D+00,3.2345D+00,3.4343D+00,3.6388D+00,3.8488D+00,4.0682D+00,
     &4.3043D+00,4.5737D+00,4.9280D+00,0.0000D+00,1.5226D+00,1.5047D+00,
     &1.4874D+00,1.4702D+00,1.4530D+00,1.4363D+00,1.4193D+00,1.4023D+00,
     &1.3860D+00,1.3690D+00,1.3520D+00,1.3351D+00,1.3190D+00,1.3083D+00,
     &1.2983D+00,1.2858D+00,1.2733D+00,1.2606D+00,1.2476D+00,1.2362D+00,
     &1.2237D+00,1.2092D+00,1.1943D+00,1.1795D+00,1.1645D+00,1.1509D+00,
     &1.1365D+00,1.1185D+00,1.0994D+00,1.0784D+00,1.0566D+00,1.0388D+00,
     &1.0195D+00,9.9801D-01,9.7765D-01,9.6019D-01,9.4712D-01,9.4158D-01,
     &9.4524D-01,9.6454D-01,1.0088D+00,1.0604D+00,1.1346D+00,1.2112D+00,
     &1.3076D+00,1.4038D+00,1.4995D+00,1.5957D+00,1.6918D+00,1.7888D+00,
     &1.8877D+00,1.9877D+00,2.0896D+00,2.1940D+00,2.2999D+00,2.3168D+00,
     &2.4844D+00,2.6497D+00,2.8098D+00,2.9678D+00,3.1219D+00,3.2743D+00,
     &3.4260D+00,3.5742D+00,3.7237D+00,3.8717D+00,4.0300D+00,0.0000D+00,
     &1.5849D+00,1.5662D+00,1.5482D+00,1.5298D+00,1.5130D+00,1.4944D+00,
     &1.4769D+00,1.4593D+00,1.4423D+00,1.4243D+00,1.4066D+00,1.3894D+00,
     &1.3720D+00,1.3607D+00,1.3499D+00,1.3366D+00,1.3237D+00,1.3101D+00,
     &1.2963D+00,1.2840D+00,1.2709D+00,1.2553D+00,1.2396D+00,1.2232D+00/
      DATA (XDVF_L(K),K=  229,  342) /
     &1.2075D+00,1.1932D+00,1.1776D+00,1.1584D+00,1.1377D+00,1.1152D+00,
     &1.0922D+00,1.0729D+00,1.0524D+00,1.0294D+00,1.0074D+00,9.8843D-01,
     &9.7377D-01,9.6751D-01,9.6901D-01,9.8606D-01,1.0264D+00,1.0745D+00,
     &1.1435D+00,1.2136D+00,1.3018D+00,1.3894D+00,1.4758D+00,1.5619D+00,
     &1.6474D+00,1.7332D+00,1.8194D+00,1.9063D+00,1.9941D+00,2.0832D+00,
     &2.1725D+00,2.1789D+00,2.3166D+00,2.4460D+00,2.5708D+00,2.6884D+00,
     &2.7987D+00,2.9025D+00,2.9974D+00,3.0823D+00,3.1538D+00,3.2013D+00,
     &3.2043D+00,0.0000D+00,1.6586D+00,1.6391D+00,1.6202D+00,1.6014D+00,
     &1.5830D+00,1.5638D+00,1.5457D+00,1.5267D+00,1.5087D+00,1.4899D+00,
     &1.4711D+00,1.4517D+00,1.4340D+00,1.4224D+00,1.4107D+00,1.3972D+00,
     &1.3827D+00,1.3684D+00,1.3535D+00,1.3404D+00,1.3263D+00,1.3096D+00,
     &1.2927D+00,1.2758D+00,1.2575D+00,1.2422D+00,1.2250D+00,1.2046D+00,
     &1.1821D+00,1.1579D+00,1.1331D+00,1.1127D+00,1.0905D+00,1.0655D+00,
     &1.0415D+00,1.0207D+00,1.0042D+00,9.9612D-01,9.9507D-01,1.0089D+00,
     &1.0451D+00,1.0887D+00,1.1514D+00,1.2146D+00,1.2936D+00,1.3711D+00,
     &1.4469D+00,1.5220D+00,1.5960D+00,1.6694D+00,1.7428D+00,1.8159D+00,
     &1.8894D+00,1.9620D+00,2.0344D+00,2.0313D+00,2.1357D+00,2.2333D+00,
     &2.3215D+00,2.4009D+00,2.4706D+00,2.5292D+00,2.5750D+00,2.6036D+00,
     &2.6096D+00,2.5783D+00,2.4673D+00,0.0000D+00,1.7269D+00,1.7065D+00/
      DATA (XDVF_L(K),K=  343,  456) /
     &1.6866D+00,1.6676D+00,1.6480D+00,1.6279D+00,1.6089D+00,1.5891D+00,
     &1.5701D+00,1.5502D+00,1.5307D+00,1.5113D+00,1.4910D+00,1.4799D+00,
     &1.4673D+00,1.4526D+00,1.4373D+00,1.4221D+00,1.4060D+00,1.3922D+00,
     &1.3771D+00,1.3596D+00,1.3414D+00,1.3234D+00,1.3045D+00,1.2879D+00,
     &1.2689D+00,1.2468D+00,1.2227D+00,1.1966D+00,1.1706D+00,1.1487D+00,
     &1.1248D+00,1.0980D+00,1.0724D+00,1.0495D+00,1.0310D+00,1.0212D+00,
     &1.0181D+00,1.0291D+00,1.0609D+00,1.1002D+00,1.1563D+00,1.2136D+00,
     &1.2840D+00,1.3528D+00,1.4201D+00,1.4854D+00,1.5492D+00,1.6125D+00,
     &1.6751D+00,1.7368D+00,1.7981D+00,1.8579D+00,1.9157D+00,1.9057D+00,
     &1.9875D+00,2.0577D+00,2.1190D+00,2.1700D+00,2.2094D+00,2.2370D+00,
     &2.2484D+00,2.2403D+00,2.2047D+00,2.1261D+00,1.9567D+00,0.0000D+00,
     &1.8047D+00,1.7833D+00,1.7626D+00,1.7418D+00,1.7220D+00,1.7009D+00,
     &1.6810D+00,1.6603D+00,1.6403D+00,1.6193D+00,1.5986D+00,1.5775D+00,
     &1.5570D+00,1.5441D+00,1.5309D+00,1.5156D+00,1.4991D+00,1.4828D+00,
     &1.4658D+00,1.4510D+00,1.4350D+00,1.4160D+00,1.3966D+00,1.3772D+00,
     &1.3565D+00,1.3386D+00,1.3184D+00,1.2942D+00,1.2680D+00,1.2404D+00,
     &1.2125D+00,1.1887D+00,1.1631D+00,1.1342D+00,1.1064D+00,1.0813D+00,
     &1.0608D+00,1.0480D+00,1.0426D+00,1.0500D+00,1.0774D+00,1.1111D+00,
     &1.1608D+00,1.2107D+00,1.2719D+00,1.3315D+00,1.3886D+00,1.4445D+00/
      DATA (XDVF_L(K),K=  457,  570) /
     &1.4984D+00,1.5505D+00,1.6020D+00,1.6524D+00,1.7009D+00,1.7480D+00,
     &1.7926D+00,1.7763D+00,1.8327D+00,1.8794D+00,1.9154D+00,1.9405D+00,
     &1.9531D+00,1.9537D+00,1.9362D+00,1.8986D+00,1.8325D+00,1.7203D+00,
     &1.5163D+00,0.0000D+00,1.8755D+00,1.8533D+00,1.8314D+00,1.8106D+00,
     &1.7890D+00,1.7672D+00,1.7464D+00,1.7248D+00,1.7038D+00,1.6817D+00,
     &1.6601D+00,1.6385D+00,1.6160D+00,1.6033D+00,1.5889D+00,1.5721D+00,
     &1.5552D+00,1.5380D+00,1.5199D+00,1.5042D+00,1.4871D+00,1.4670D+00,
     &1.4463D+00,1.4249D+00,1.4036D+00,1.3843D+00,1.3630D+00,1.3364D+00,
     &1.3086D+00,1.2791D+00,1.2500D+00,1.2245D+00,1.1971D+00,1.1662D+00,
     &1.1361D+00,1.1090D+00,1.0858D+00,1.0721D+00,1.0641D+00,1.0676D+00,
     &1.0898D+00,1.1195D+00,1.1627D+00,1.2069D+00,1.2603D+00,1.3118D+00,
     &1.3607D+00,1.4079D+00,1.4534D+00,1.4968D+00,1.5392D+00,1.5794D+00,
     &1.6181D+00,1.6552D+00,1.6888D+00,1.6690D+00,1.7073D+00,1.7353D+00,
     &1.7530D+00,1.7595D+00,1.7531D+00,1.7338D+00,1.6988D+00,1.6428D+00,
     &1.5583D+00,1.4293D+00,1.2136D+00,0.0000D+00,1.9470D+00,1.9238D+00,
     &1.9021D+00,1.8782D+00,1.8570D+00,1.8343D+00,1.8123D+00,1.7898D+00,
     &1.7680D+00,1.7449D+00,1.7222D+00,1.6994D+00,1.6760D+00,1.6624D+00,
     &1.6469D+00,1.6299D+00,1.6118D+00,1.5933D+00,1.5742D+00,1.5574D+00,
     &1.5392D+00,1.5179D+00,1.4955D+00,1.4738D+00,1.4506D+00,1.4300D+00/
      DATA (XDVF_L(K),K=  571,  684) /
     &1.4069D+00,1.3792D+00,1.3492D+00,1.3178D+00,1.2868D+00,1.2597D+00,
     &1.2307D+00,1.1976D+00,1.1654D+00,1.1363D+00,1.1108D+00,1.0945D+00,
     &1.0840D+00,1.0845D+00,1.1017D+00,1.1268D+00,1.1637D+00,1.2016D+00,
     &1.2473D+00,1.2910D+00,1.3324D+00,1.3719D+00,1.4090D+00,1.4450D+00,
     &1.4784D+00,1.5109D+00,1.5404D+00,1.5681D+00,1.5925D+00,1.5689D+00,
     &1.5916D+00,1.6043D+00,1.6067D+00,1.5981D+00,1.5779D+00,1.5449D+00,
     &1.4949D+00,1.4262D+00,1.3303D+00,1.1932D+00,9.7657D-01,0.0000D+00,
     &2.0122D+00,1.9881D+00,1.9640D+00,1.9418D+00,1.9190D+00,1.8954D+00,
     &1.8721D+00,1.8492D+00,1.8262D+00,1.8024D+00,1.7784D+00,1.7550D+00,
     &1.7300D+00,1.7157D+00,1.6999D+00,1.6818D+00,1.6627D+00,1.6435D+00,
     &1.6233D+00,1.6058D+00,1.5866D+00,1.5643D+00,1.5417D+00,1.5178D+00,
     &1.4926D+00,1.4705D+00,1.4465D+00,1.4174D+00,1.3856D+00,1.3527D+00,
     &1.3198D+00,1.2914D+00,1.2605D+00,1.2257D+00,1.1915D+00,1.1601D+00,
     &1.1326D+00,1.1142D+00,1.1016D+00,1.0982D+00,1.1114D+00,1.1321D+00,
     &1.1637D+00,1.1958D+00,1.2352D+00,1.2722D+00,1.3071D+00,1.3397D+00,
     &1.3704D+00,1.3995D+00,1.4267D+00,1.4516D+00,1.4736D+00,1.4942D+00,
     &1.5100D+00,1.4848D+00,1.4955D+00,1.4964D+00,1.4873D+00,1.4675D+00,
     &1.4366D+00,1.3933D+00,1.3349D+00,1.2585D+00,1.1565D+00,1.0171D+00,
     &8.0601D-01,0.0000D+00,2.0789D+00,2.0539D+00,2.0294D+00,2.0053D+00/
      DATA (XDVF_L(K),K=  685,  798) /
     &1.9820D+00,1.9581D+00,1.9336D+00,1.9096D+00,1.8860D+00,1.8609D+00,
     &1.8367D+00,1.8106D+00,1.7860D+00,1.7706D+00,1.7543D+00,1.7350D+00,
     &1.7150D+00,1.6945D+00,1.6735D+00,1.6550D+00,1.6349D+00,1.6112D+00,
     &1.5864D+00,1.5617D+00,1.5356D+00,1.5128D+00,1.4868D+00,1.4555D+00,
     &1.4224D+00,1.3876D+00,1.3532D+00,1.3231D+00,1.2904D+00,1.2536D+00,
     &1.2173D+00,1.1838D+00,1.1545D+00,1.1338D+00,1.1185D+00,1.1113D+00,
     &1.1199D+00,1.1362D+00,1.1627D+00,1.1895D+00,1.2222D+00,1.2529D+00,
     &1.2813D+00,1.3080D+00,1.3324D+00,1.3546D+00,1.3756D+00,1.3938D+00,
     &1.4103D+00,1.4232D+00,1.4319D+00,1.4055D+00,1.4052D+00,1.3959D+00,
     &1.3768D+00,1.3480D+00,1.3084D+00,1.2576D+00,1.1928D+00,1.1110D+00,
     &1.0066D+00,8.6804D-01,6.6615D-01,0.0000D+00,2.1434D+00,2.1178D+00,
     &2.0930D+00,2.0676D+00,2.0440D+00,2.0184D+00,1.9935D+00,1.9686D+00,
     &1.9439D+00,1.9179D+00,1.8915D+00,1.8663D+00,1.8400D+00,1.8239D+00,
     &1.8067D+00,1.7863D+00,1.7654D+00,1.7440D+00,1.7219D+00,1.7025D+00,
     &1.6814D+00,1.6565D+00,1.6311D+00,1.6045D+00,1.5766D+00,1.5526D+00,
     &1.5250D+00,1.4925D+00,1.4574D+00,1.4213D+00,1.3849D+00,1.3532D+00,
     &1.3191D+00,1.2800D+00,1.2418D+00,1.2062D+00,1.1743D+00,1.1517D+00,
     &1.1338D+00,1.1237D+00,1.1272D+00,1.1399D+00,1.1608D+00,1.1828D+00,
     &1.2092D+00,1.2341D+00,1.2570D+00,1.2774D+00,1.2962D+00,1.3135D+00/
      DATA (XDVF_L(K),K=  799,  912) /
     &1.3280D+00,1.3406D+00,1.3511D+00,1.3588D+00,1.3613D+00,1.3335D+00,
     &1.3246D+00,1.3067D+00,1.2801D+00,1.2441D+00,1.1985D+00,1.1418D+00,
     &1.0724D+00,9.8806D-01,8.8293D-01,7.4746D-01,5.5665D-01,0.0000D+00,
     &2.2035D+00,2.1769D+00,2.1514D+00,2.1259D+00,2.1000D+00,2.0743D+00,
     &2.0488D+00,2.0226D+00,1.9973D+00,1.9702D+00,1.9428D+00,1.9166D+00,
     &1.8890D+00,1.8729D+00,1.8548D+00,1.8337D+00,1.8116D+00,1.7895D+00,
     &1.7662D+00,1.7461D+00,1.7239D+00,1.6980D+00,1.6714D+00,1.6436D+00,
     &1.6146D+00,1.5889D+00,1.5604D+00,1.5266D+00,1.4895D+00,1.4515D+00,
     &1.4138D+00,1.3806D+00,1.3448D+00,1.3040D+00,1.2638D+00,1.2261D+00,
     &1.1920D+00,1.1669D+00,1.1469D+00,1.1341D+00,1.1335D+00,1.1420D+00,
     &1.1583D+00,1.1760D+00,1.1971D+00,1.2168D+00,1.2343D+00,1.2501D+00,
     &1.2640D+00,1.2762D+00,1.2866D+00,1.2942D+00,1.2996D+00,1.3020D+00,
     &1.3003D+00,1.2725D+00,1.2557D+00,1.2312D+00,1.1982D+00,1.1569D+00,
     &1.1068D+00,1.0465D+00,9.7460D-01,8.8884D-01,7.8459D-01,6.5333D-01,
     &4.7359D-01,0.0000D+00,2.2800D+00,2.2524D+00,2.2256D+00,2.1987D+00,
     &2.1730D+00,2.1459D+00,2.1192D+00,2.0922D+00,2.0656D+00,2.0374D+00,
     &2.0100D+00,1.9802D+00,1.9520D+00,1.9346D+00,1.9156D+00,1.8937D+00,
     &1.8706D+00,1.8475D+00,1.8228D+00,1.8017D+00,1.7783D+00,1.7509D+00,
     &1.7221D+00,1.6937D+00,1.6627D+00,1.6354D+00,1.6050D+00,1.5688D+00/
      DATA (XDVF_L(K),K=  913, 1026) /
     &1.5301D+00,1.4898D+00,1.4503D+00,1.4150D+00,1.3772D+00,1.3339D+00,
     &1.2911D+00,1.2510D+00,1.2138D+00,1.1866D+00,1.1637D+00,1.1458D+00,
     &1.1403D+00,1.1441D+00,1.1548D+00,1.1669D+00,1.1817D+00,1.1950D+00,
     &1.2065D+00,1.2163D+00,1.2249D+00,1.2313D+00,1.2355D+00,1.2379D+00,
     &1.2379D+00,1.2348D+00,1.2275D+00,1.1987D+00,1.1744D+00,1.1427D+00,
     &1.1035D+00,1.0570D+00,1.0018D+00,9.3862D-01,8.6494D-01,7.7913D-01,
     &6.7747D-01,5.5266D-01,3.8741D-01,0.0000D+00,2.3524D+00,2.3243D+00,
     &2.2963D+00,2.2689D+00,2.2420D+00,2.2137D+00,2.1858D+00,2.1579D+00,
     &2.1301D+00,2.1011D+00,2.0718D+00,2.0424D+00,2.0120D+00,1.9937D+00,
     &1.9743D+00,1.9509D+00,1.9267D+00,1.9020D+00,1.8763D+00,1.8541D+00,
     &1.8295D+00,1.8006D+00,1.7713D+00,1.7402D+00,1.7077D+00,1.6794D+00,
     &1.6475D+00,1.6087D+00,1.5679D+00,1.5259D+00,1.4840D+00,1.4470D+00,
     &1.4072D+00,1.3615D+00,1.3163D+00,1.2738D+00,1.2336D+00,1.2045D+00,
     &1.1783D+00,1.1563D+00,1.1459D+00,1.1457D+00,1.1504D+00,1.1577D+00,
     &1.1662D+00,1.1742D+00,1.1807D+00,1.1857D+00,1.1886D+00,1.1902D+00,
     &1.1899D+00,1.1878D+00,1.1830D+00,1.1751D+00,1.1633D+00,1.1345D+00,
     &1.1039D+00,1.0667D+00,1.0230D+00,9.7228D-01,9.1417D-01,8.4905D-01,
     &7.7478D-01,6.9004D-01,5.9155D-01,4.7371D-01,3.2191D-01,0.0000D+00,
     &2.4233D+00,2.3947D+00,2.3653D+00,2.3365D+00,2.3090D+00,2.2800D+00/
      DATA (XDVF_L(K),K= 1027, 1140) /
     &2.2512D+00,2.2220D+00,2.1934D+00,2.1628D+00,2.1319D+00,2.1007D+00,
     &2.0700D+00,2.0512D+00,2.0301D+00,2.0057D+00,1.9809D+00,1.9549D+00,
     &1.9281D+00,1.9049D+00,1.8791D+00,1.8497D+00,1.8175D+00,1.7854D+00,
     &1.7507D+00,1.7209D+00,1.6878D+00,1.6474D+00,1.6047D+00,1.5603D+00,
     &1.5164D+00,1.4777D+00,1.4358D+00,1.3879D+00,1.3403D+00,1.2952D+00,
     &1.2523D+00,1.2206D+00,1.1913D+00,1.1661D+00,1.1505D+00,1.1462D+00,
     &1.1460D+00,1.1481D+00,1.1518D+00,1.1545D+00,1.1559D+00,1.1562D+00,
     &1.1548D+00,1.1523D+00,1.1478D+00,1.1414D+00,1.1331D+00,1.1212D+00,
     &1.1055D+00,1.0763D+00,1.0405D+00,9.9877D-01,9.5130D-01,8.9815D-01,
     &8.3813D-01,7.7188D-01,6.9792D-01,6.1492D-01,5.2020D-01,4.0920D-01,
     &2.7020D-01,0.0000D+00,2.4906D+00,2.4607D+00,2.4307D+00,2.4014D+00,
     &2.3730D+00,2.3427D+00,2.3127D+00,2.2828D+00,2.2528D+00,2.2213D+00,
     &2.1903D+00,2.1577D+00,2.1250D+00,2.1053D+00,2.0839D+00,2.0583D+00,
     &2.0318D+00,2.0051D+00,1.9771D+00,1.9527D+00,1.9259D+00,1.8935D+00,
     &1.8607D+00,1.8269D+00,1.7917D+00,1.7606D+00,1.7253D+00,1.6833D+00,
     &1.6387D+00,1.5925D+00,1.5465D+00,1.5061D+00,1.4624D+00,1.4121D+00,
     &1.3623D+00,1.3152D+00,1.2700D+00,1.2349D+00,1.2036D+00,1.1745D+00,
     &1.1544D+00,1.1457D+00,1.1410D+00,1.1389D+00,1.1378D+00,1.1357D+00,
     &1.1332D+00,1.1290D+00,1.1244D+00,1.1176D+00,1.1099D+00,1.0996D+00/
      DATA (XDVF_L(K),K= 1141, 1254) /
     &1.0875D+00,1.0729D+00,1.0538D+00,1.0249D+00,9.8511D-01,9.3994D-01,
     &8.8948D-01,8.3410D-01,7.7332D-01,7.0681D-01,6.3377D-01,5.5280D-01,
     &4.6214D-01,3.5755D-01,2.2965D-01,0.0000D+00,2.5589D+00,2.5291D+00,
     &2.4979D+00,2.4676D+00,2.4370D+00,2.4060D+00,2.3753D+00,2.3443D+00,
     &2.3135D+00,2.2809D+00,2.2486D+00,2.2146D+00,2.1810D+00,2.1602D+00,
     &2.1376D+00,2.1114D+00,2.0841D+00,2.0557D+00,2.0265D+00,2.0011D+00,
     &1.9730D+00,1.9392D+00,1.9055D+00,1.8697D+00,1.8327D+00,1.8003D+00,
     &1.7635D+00,1.7197D+00,1.6727D+00,1.6246D+00,1.5770D+00,1.5346D+00,
     &1.4890D+00,1.4363D+00,1.3841D+00,1.3341D+00,1.2867D+00,1.2492D+00,
     &1.2151D+00,1.1824D+00,1.1578D+00,1.1451D+00,1.1356D+00,1.1298D+00,
     &1.1233D+00,1.1169D+00,1.1105D+00,1.1027D+00,1.0940D+00,1.0840D+00,
     &1.0726D+00,1.0592D+00,1.0444D+00,1.0265D+00,1.0045D+00,9.7613D-01,
     &9.3249D-01,8.8451D-01,8.3193D-01,7.7510D-01,7.1373D-01,6.4749D-01,
     &5.7554D-01,4.9725D-01,4.1072D-01,3.1254D-01,1.9551D-01,0.0000D+00,
     &2.6244D+00,2.5927D+00,2.5615D+00,2.5299D+00,2.4990D+00,2.4671D+00,
     &2.4356D+00,2.4034D+00,2.3717D+00,2.3377D+00,2.3034D+00,2.2689D+00,
     &2.2340D+00,2.2126D+00,2.1892D+00,2.1616D+00,2.1331D+00,2.1040D+00,
     &2.0736D+00,2.0471D+00,2.0180D+00,1.9830D+00,1.9472D+00,1.9112D+00,
     &1.8717D+00,1.8375D+00,1.7996D+00,1.7538D+00,1.7053D+00,1.6548D+00/
      DATA (XDVF_L(K),K= 1255, 1368) /
     &1.6053D+00,1.5612D+00,1.5138D+00,1.4590D+00,1.4045D+00,1.3516D+00,
     &1.3023D+00,1.2626D+00,1.2251D+00,1.1889D+00,1.1601D+00,1.1441D+00,
     &1.1302D+00,1.1201D+00,1.1098D+00,1.0996D+00,1.0888D+00,1.0782D+00,
     &1.0659D+00,1.0531D+00,1.0388D+00,1.0228D+00,1.0047D+00,9.8480D-01,
     &9.6040D-01,9.3234D-01,8.8589D-01,8.3563D-01,7.8162D-01,7.2366D-01,
     &6.6215D-01,5.9658D-01,5.2617D-01,4.5043D-01,3.6787D-01,2.7575D-01,
     &1.6826D-01,0.0000D+00,2.6886D+00,2.6564D+00,2.6234D+00,2.5908D+00,
     &2.5600D+00,2.5268D+00,2.4943D+00,2.4612D+00,2.4283D+00,2.3924D+00,
     &2.3582D+00,2.3219D+00,2.2860D+00,2.2642D+00,2.2394D+00,2.2113D+00,
     &2.1817D+00,2.1512D+00,2.1198D+00,2.0920D+00,2.0618D+00,2.0268D+00,
     &1.9890D+00,1.9503D+00,1.9098D+00,1.8739D+00,1.8343D+00,1.7867D+00,
     &1.7365D+00,1.6843D+00,1.6329D+00,1.5870D+00,1.5377D+00,1.4807D+00,
     &1.4239D+00,1.3692D+00,1.3169D+00,1.2751D+00,1.2350D+00,1.1954D+00,
     &1.1624D+00,1.1425D+00,1.1247D+00,1.1110D+00,1.0963D+00,1.0827D+00,
     &1.0687D+00,1.0547D+00,1.0396D+00,1.0240D+00,1.0070D+00,9.8853D-01,
     &9.6834D-01,9.4569D-01,9.1962D-01,8.9220D-01,8.4321D-01,7.9105D-01,
     &7.3592D-01,6.7777D-01,6.1620D-01,5.5143D-01,4.8272D-01,4.0962D-01,
     &3.3102D-01,2.4455D-01,1.4574D-01,0.0000D+00,2.7496D+00,2.7153D+00,
     &2.6835D+00,2.6504D+00,2.6180D+00,2.5834D+00,2.5502D+00,2.5161D+00/
      DATA (XDVF_L(K),K= 1369, 1482) /
     &2.4824D+00,2.4466D+00,2.4095D+00,2.3736D+00,2.3360D+00,2.3124D+00,
     &2.2875D+00,2.2580D+00,2.2274D+00,2.1960D+00,2.1631D+00,2.1347D+00,
     &2.1032D+00,2.0670D+00,2.0277D+00,1.9882D+00,1.9458D+00,1.9086D+00,
     &1.8675D+00,1.8179D+00,1.7658D+00,1.7122D+00,1.6586D+00,1.6112D+00,
     &1.5600D+00,1.5010D+00,1.4420D+00,1.3855D+00,1.3294D+00,1.2858D+00,
     &1.2435D+00,1.2006D+00,1.1641D+00,1.1410D+00,1.1193D+00,1.1023D+00,
     &1.0837D+00,1.0664D+00,1.0496D+00,1.0329D+00,1.0157D+00,9.9745D-01,
     &9.7803D-01,9.5735D-01,9.3539D-01,9.1075D-01,8.8302D-01,8.5608D-01,
     &8.0509D-01,7.5168D-01,6.9580D-01,6.3743D-01,5.7619D-01,5.1233D-01,
     &4.4547D-01,3.7496D-01,2.9995D-01,2.1862D-01,1.2745D-01,0.0000D+00,
     &2.8331D+00,2.7978D+00,2.7648D+00,2.7299D+00,2.6960D+00,2.6609D+00,
     &2.6263D+00,2.5910D+00,2.5561D+00,2.5197D+00,2.4802D+00,2.4424D+00,
     &2.4030D+00,2.3791D+00,2.3526D+00,2.3216D+00,2.2897D+00,2.2570D+00,
     &2.2225D+00,2.1925D+00,2.1595D+00,2.1199D+00,2.0799D+00,2.0383D+00,
     &1.9938D+00,1.9551D+00,1.9121D+00,1.8601D+00,1.8054D+00,1.7494D+00,
     &1.6932D+00,1.6435D+00,1.5898D+00,1.5280D+00,1.4659D+00,1.4056D+00,
     &1.3471D+00,1.3010D+00,1.2550D+00,1.2078D+00,1.1652D+00,1.1383D+00,
     &1.1114D+00,1.0902D+00,1.0668D+00,1.0451D+00,1.0248D+00,1.0039D+00,
     &9.8353D-01,9.6205D-01,9.4076D-01,9.1705D-01,8.9229D-01,8.6577D-01/
      DATA (XDVF_L(K),K= 1483, 1596) /
     &8.3604D-01,8.0985D-01,7.5687D-01,7.0190D-01,6.4516D-01,5.8700D-01,
     &5.2660D-01,4.6452D-01,3.9995D-01,3.3310D-01,2.6289D-01,1.8826D-01,
     &1.0655D-01,0.0000D+00,2.9096D+00,2.8732D+00,2.8390D+00,2.8027D+00,
     &2.7690D+00,2.7325D+00,2.6961D+00,2.6597D+00,2.6231D+00,2.5833D+00,
     &2.5456D+00,2.5047D+00,2.4650D+00,2.4391D+00,2.4120D+00,2.3799D+00,
     &2.3462D+00,2.3123D+00,2.2763D+00,2.2451D+00,2.2108D+00,2.1692D+00,
     &2.1276D+00,2.0835D+00,2.0378D+00,1.9974D+00,1.9525D+00,1.8983D+00,
     &1.8413D+00,1.7827D+00,1.7243D+00,1.6725D+00,1.6166D+00,1.5520D+00,
     &1.4872D+00,1.4244D+00,1.3627D+00,1.3136D+00,1.2649D+00,1.2130D+00,
     &1.1663D+00,1.1352D+00,1.1040D+00,1.0787D+00,1.0514D+00,1.0264D+00,
     &1.0021D+00,9.7883D-01,9.5548D-01,9.3171D-01,9.0763D-01,8.8283D-01,
     &8.5596D-01,8.2732D-01,7.9601D-01,7.7056D-01,7.1598D-01,6.6027D-01,
     &6.0340D-01,5.4514D-01,4.8601D-01,4.2556D-01,3.6359D-01,2.9984D-01,
     &2.3396D-01,1.6486D-01,9.0844D-02,0.0000D+00,2.9880D+00,2.9510D+00,
     &2.9150D+00,2.8782D+00,2.8430D+00,2.8048D+00,2.7677D+00,2.7301D+00,
     &2.6924D+00,2.6517D+00,2.6110D+00,2.5696D+00,2.5280D+00,2.5017D+00,
     &2.4728D+00,2.4393D+00,2.4042D+00,2.3687D+00,2.3313D+00,2.2988D+00,
     &2.2631D+00,2.2204D+00,2.1768D+00,2.1312D+00,2.0828D+00,2.0405D+00,
     &1.9928D+00,1.9364D+00,1.8772D+00,1.8164D+00,1.7558D+00,1.7018D+00/
      DATA (XDVF_L(K),K= 1597, 1710) /
     &1.6434D+00,1.5762D+00,1.5084D+00,1.4432D+00,1.3783D+00,1.3261D+00,
     &1.2741D+00,1.2182D+00,1.1669D+00,1.1315D+00,1.0961D+00,1.0671D+00,
     &1.0360D+00,1.0071D+00,9.7992D-01,9.5371D-01,9.2801D-01,9.0200D-01,
     &8.7588D-01,8.4862D-01,8.2038D-01,7.9020D-01,7.5770D-01,7.3298D-01,
     &6.7721D-01,6.2090D-01,5.6394D-01,5.0631D-01,4.4841D-01,3.8970D-01,
     &3.3019D-01,2.6973D-01,2.0791D-01,1.4420D-01,7.7416D-02,0.0000D+00,
     &3.0661D+00,3.0288D+00,2.9911D+00,2.9537D+00,2.9160D+00,2.8778D+00,
     &2.8392D+00,2.8000D+00,2.7610D+00,2.7200D+00,2.6782D+00,2.6345D+00,
     &2.5900D+00,2.5625D+00,2.5329D+00,2.4982D+00,2.4617D+00,2.4247D+00,
     &2.3857D+00,2.3518D+00,2.3145D+00,2.2697D+00,2.2245D+00,2.1764D+00,
     &2.1269D+00,2.0819D+00,2.0331D+00,1.9746D+00,1.9126D+00,1.8497D+00,
     &1.7862D+00,1.7303D+00,1.6696D+00,1.5995D+00,1.5285D+00,1.4608D+00,
     &1.3929D+00,1.3377D+00,1.2826D+00,1.2228D+00,1.1669D+00,1.1279D+00,
     &1.0882D+00,1.0555D+00,1.0205D+00,9.8876D-01,9.5876D-01,9.2969D-01,
     &9.0171D-01,8.7356D-01,8.4551D-01,8.1668D-01,7.8701D-01,7.5564D-01,
     &7.2196D-01,6.9797D-01,6.4121D-01,5.8469D-01,5.2810D-01,4.7131D-01,
     &4.1460D-01,3.5783D-01,3.0063D-01,2.4338D-01,1.8544D-01,1.2660D-01,
     &6.6270D-02,0.0000D+00,3.1379D+00,3.0995D+00,3.0600D+00,3.0213D+00,
     &2.9840D+00,2.9442D+00,2.9047D+00,2.8641D+00,2.8239D+00,2.7813D+00/
      DATA (XDVF_L(K),K= 1711, 1824) /
     &2.7383D+00,2.6928D+00,2.6470D+00,2.6191D+00,2.5880D+00,2.5519D+00,
     &2.5145D+00,2.4761D+00,2.4357D+00,2.4004D+00,2.3615D+00,2.3153D+00,
     &2.2678D+00,2.2180D+00,2.1669D+00,2.1208D+00,2.0699D+00,2.0087D+00,
     &1.9447D+00,1.8795D+00,1.8139D+00,1.7558D+00,1.6930D+00,1.6205D+00,
     &1.5467D+00,1.4759D+00,1.4054D+00,1.3484D+00,1.2895D+00,1.2267D+00,
     &1.1663D+00,1.1242D+00,1.0808D+00,1.0449D+00,1.0065D+00,9.7194D-01,
     &9.3967D-01,9.0840D-01,8.7834D-01,8.4891D-01,8.1928D-01,7.8930D-01,
     &7.5803D-01,7.2562D-01,6.9124D-01,6.6796D-01,6.1058D-01,5.5392D-01,
     &4.9752D-01,4.4176D-01,3.8633D-01,3.3127D-01,2.7648D-01,2.2186D-01,
     &1.6735D-01,1.1268D-01,5.7652D-02,0.0000D+00,3.2129D+00,3.1726D+00,
     &3.1325D+00,3.0928D+00,3.0540D+00,3.0127D+00,2.9717D+00,2.9303D+00,
     &2.8887D+00,2.8449D+00,2.8001D+00,2.7537D+00,2.7060D+00,2.6766D+00,
     &2.6453D+00,2.6073D+00,2.5683D+00,2.5286D+00,2.4866D+00,2.4501D+00,
     &2.4107D+00,2.3628D+00,2.3125D+00,2.2620D+00,2.2079D+00,2.1597D+00,
     &2.1067D+00,2.0440D+00,1.9778D+00,1.9097D+00,1.8421D+00,1.7819D+00,
     &1.7169D+00,1.6416D+00,1.5664D+00,1.4922D+00,1.4189D+00,1.3583D+00,
     &1.2971D+00,1.2300D+00,1.1652D+00,1.1200D+00,1.0729D+00,1.0343D+00,
     &9.9254D-01,9.5513D-01,9.2006D-01,8.8711D-01,8.5555D-01,8.2426D-01,
     &7.9305D-01,7.6193D-01,7.2963D-01,6.9636D-01,6.6128D-01,6.3868D-01/
      DATA (XDVF_L(K),K= 1825, 1836) /
     &5.8093D-01,5.2428D-01,4.6858D-01,4.1372D-01,3.5972D-01,3.0648D-01,
     &2.5392D-01,2.0208D-01,1.5083D-01,1.0018D-01,5.0068D-02,0.0000D+00/
      DATA (XDEF_L(K),K=    1,  114) /
     &4.3007D-01,4.2474D-01,4.1967D-01,4.1458D-01,4.0970D-01,4.0443D-01,
     &3.9925D-01,3.9397D-01,3.8864D-01,3.8302D-01,3.7707D-01,3.7100D-01,
     &3.6470D-01,3.6080D-01,3.5639D-01,3.5109D-01,3.4531D-01,3.3914D-01,
     &3.3238D-01,3.2609D-01,3.1913D-01,3.1062D-01,3.0152D-01,2.9176D-01,
     &2.8100D-01,2.7114D-01,2.5952D-01,2.4467D-01,2.2784D-01,2.0937D-01,
     &1.9117D-01,1.7470D-01,1.5685D-01,1.3678D-01,1.1825D-01,1.0349D-01,
     &9.4854D-02,9.5054D-02,1.0589D-01,1.3527D-01,1.8584D-01,2.3426D-01,
     &2.9021D-01,3.3527D-01,3.7670D-01,4.0255D-01,4.1326D-01,4.0880D-01,
     &3.8831D-01,3.5045D-01,2.9287D-01,2.1298D-01,1.0773D-01,0.0000D+00,
     &0.0000D+00,2.0644D-01,1.5422D-01,1.0950D-01,7.3614D-02,4.6726D-02,
     &2.7433D-02,1.4144D-02,6.5080D-03,2.4719D-03,0.0000D+00,0.0000D+00,
     &0.0000D+00,0.0000D+00,4.4398D-01,4.3864D-01,4.3346D-01,4.2809D-01,
     &4.2290D-01,4.1747D-01,4.1205D-01,4.0650D-01,4.0098D-01,3.9480D-01,
     &3.8873D-01,3.8226D-01,3.7560D-01,3.7145D-01,3.6678D-01,3.6108D-01,
     &3.5488D-01,3.4833D-01,3.4123D-01,3.3464D-01,3.2718D-01,3.1811D-01,
     &3.0838D-01,2.9811D-01,2.8670D-01,2.7630D-01,2.6412D-01,2.4861D-01,
     &2.3110D-01,2.1209D-01,1.9355D-01,1.7681D-01,1.5878D-01,1.3870D-01,
     &1.2044D-01,1.0620D-01,9.8341D-02,9.9345D-02,1.1086D-01,1.4055D-01,
     &1.9033D-01,2.3696D-01,2.8983D-01,3.3137D-01,3.6834D-01,3.8982D-01/
      DATA (XDEF_L(K),K=  115,  228) /
     &3.9672D-01,3.8896D-01,3.6609D-01,3.2678D-01,2.6933D-01,1.9181D-01,
     &9.1683D-02,0.0000D+00,0.0000D+00,1.8955D-01,1.4041D-01,9.8873D-02,
     &6.5928D-02,4.1462D-02,2.3905D-02,1.2324D-02,5.6113D-03,2.1050D-03,
     &0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,4.5980D-01,4.5420D-01,
     &4.4884D-01,4.4319D-01,4.3780D-01,4.3208D-01,4.2642D-01,4.2053D-01,
     &4.1457D-01,4.0824D-01,4.0181D-01,3.9484D-01,3.8780D-01,3.8328D-01,
     &3.7831D-01,3.7223D-01,3.6559D-01,3.5853D-01,3.5072D-01,3.4400D-01,
     &3.3590D-01,3.2633D-01,3.1598D-01,3.0508D-01,2.9301D-01,2.8197D-01,
     &2.6915D-01,2.5289D-01,2.3470D-01,2.1511D-01,1.9623D-01,1.7918D-01,
     &1.6098D-01,1.4092D-01,1.2294D-01,1.0928D-01,1.0224D-01,1.0401D-01,
     &1.1623D-01,1.4620D-01,1.9488D-01,2.3948D-01,2.8894D-01,3.2681D-01,
     &3.5905D-01,3.7613D-01,3.7908D-01,3.6817D-01,3.4299D-01,3.0266D-01,
     &2.4596D-01,1.7115D-01,7.6792D-02,0.0000D+00,0.0000D+00,1.7267D-01,
     &1.2670D-01,8.8446D-02,5.8458D-02,3.6380D-02,2.0551D-02,1.0608D-02,
     &4.7732D-03,1.7670D-03,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &4.7845D-01,4.7258D-01,4.6687D-01,4.6107D-01,4.5540D-01,4.4938D-01,
     &4.4336D-01,4.3728D-01,4.3070D-01,4.2403D-01,4.1702D-01,4.0968D-01,
     &4.0210D-01,3.9723D-01,3.9181D-01,3.8522D-01,3.7808D-01,3.7047D-01,
     &3.6211D-01,3.5469D-01,3.4619D-01,3.3582D-01,3.2478D-01,3.1314D-01/
      DATA (XDEF_L(K),K=  229,  342) /
     &3.0021D-01,2.8848D-01,2.7488D-01,2.5781D-01,2.3886D-01,2.1865D-01,
     &1.9932D-01,1.8196D-01,1.6359D-01,1.4359D-01,1.2596D-01,1.1295D-01,
     &1.0678D-01,1.0933D-01,1.2234D-01,1.5242D-01,1.9969D-01,2.4187D-01,
     &2.8742D-01,3.2112D-01,3.4825D-01,3.6067D-01,3.5959D-01,3.4546D-01,
     &3.1813D-01,2.7719D-01,2.2151D-01,1.5037D-01,6.2862D-02,0.0000D+00,
     &0.0000D+00,1.5516D-01,1.1270D-01,7.7856D-02,5.0916D-02,3.1337D-02,
     &1.7279D-02,8.9355D-03,3.9672D-03,1.4465D-03,0.0000D+00,0.0000D+00,
     &0.0000D+00,0.0000D+00,5.0059D-01,4.9450D-01,4.8826D-01,4.8213D-01,
     &4.7610D-01,4.6972D-01,4.6326D-01,4.5655D-01,4.4999D-01,4.4265D-01,
     &4.3505D-01,4.2703D-01,4.1870D-01,4.1345D-01,4.0758D-01,4.0034D-01,
     &3.9260D-01,3.8434D-01,3.7539D-01,3.6725D-01,3.5804D-01,3.4696D-01,
     &3.3492D-01,3.2231D-01,3.0852D-01,2.9601D-01,2.8154D-01,2.6348D-01,
     &2.4363D-01,2.2272D-01,2.0295D-01,1.8526D-01,1.6669D-01,1.4678D-01,
     &1.2956D-01,1.1726D-01,1.1212D-01,1.1548D-01,1.2910D-01,1.5906D-01,
     &2.0458D-01,2.4395D-01,2.8508D-01,3.1418D-01,3.3593D-01,3.4343D-01,
     &3.3827D-01,3.2104D-01,2.9189D-01,2.5067D-01,1.9688D-01,1.3016D-01,
     &5.0498D-02,0.0000D+00,0.0000D+00,1.3742D-01,9.8602D-02,6.7357D-02,
     &4.3555D-02,2.6444D-02,1.4175D-02,7.3561D-03,3.2181D-03,1.1530D-03,
     &0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,5.2114D-01,5.1454D-01/
      DATA (XDEF_L(K),K=  343,  456) /
     &5.0806D-01,5.0160D-01,4.9520D-01,4.8843D-01,4.8165D-01,4.7456D-01,
     &4.6738D-01,4.5962D-01,4.5149D-01,4.4293D-01,4.3400D-01,4.2833D-01,
     &4.2194D-01,4.1420D-01,4.0580D-01,3.9678D-01,3.8741D-01,3.7848D-01,
     &3.6878D-01,3.5682D-01,3.4416D-01,3.3062D-01,3.1602D-01,3.0269D-01,
     &2.8749D-01,2.6857D-01,2.4798D-01,2.2641D-01,2.0626D-01,1.8828D-01,
     &1.6960D-01,1.4976D-01,1.3293D-01,1.2126D-01,1.1684D-01,1.2099D-01,
     &1.3505D-01,1.6471D-01,2.0841D-01,2.4521D-01,2.8248D-01,3.0770D-01,
     &3.2484D-01,3.2845D-01,3.1999D-01,3.0047D-01,2.7030D-01,2.2924D-01,
     &1.7739D-01,1.1482D-01,4.2174D-02,0.0000D+00,0.0000D+00,1.2330D-01,
     &8.7586D-02,5.9211D-02,3.7890D-02,2.2733D-02,1.1877D-02,6.1865D-03,
     &2.6713D-03,9.4247D-04,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &5.4423D-01,5.3740D-01,5.3068D-01,5.2385D-01,5.1700D-01,5.0982D-01,
     &5.0256D-01,4.9509D-01,4.8731D-01,4.7895D-01,4.7023D-01,4.6094D-01,
     &4.5130D-01,4.4506D-01,4.3820D-01,4.2973D-01,4.2069D-01,4.1108D-01,
     &4.0069D-01,3.9131D-01,3.8063D-01,3.6796D-01,3.5430D-01,3.3991D-01,
     &3.2433D-01,3.1014D-01,2.9407D-01,2.7418D-01,2.5281D-01,2.3056D-01,
     &2.0999D-01,1.9171D-01,1.7291D-01,1.5321D-01,1.3677D-01,1.2578D-01,
     &1.2220D-01,1.2696D-01,1.4132D-01,1.7056D-01,2.1212D-01,2.4603D-01,
     &2.7912D-01,3.0023D-01,3.1274D-01,3.1234D-01,3.0087D-01,2.7925D-01/
      DATA (XDEF_L(K),K=  457,  570) /
     &2.4820D-01,2.0782D-01,1.5841D-01,1.0056D-01,3.5470D-02,0.0000D+00,
     &0.0000D+00,1.0941D-01,7.6864D-02,5.1391D-02,3.2506D-02,1.9250D-02,
     &9.7741D-03,5.1192D-03,2.1775D-03,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.0000D+00,0.0000D+00,5.6542D-01,5.5814D-01,5.5101D-01,5.4385D-01,
     &5.3670D-01,5.2913D-01,5.2140D-01,5.1352D-01,5.0533D-01,4.9639D-01,
     &4.8702D-01,4.7710D-01,4.6670D-01,4.6011D-01,4.5270D-01,4.4365D-01,
     &4.3394D-01,4.2383D-01,4.1271D-01,4.0253D-01,3.9137D-01,3.7783D-01,
     &3.6325D-01,3.4810D-01,3.3163D-01,3.1674D-01,2.9988D-01,2.7922D-01,
     &2.5706D-01,2.3429D-01,2.1333D-01,1.9484D-01,1.7592D-01,1.5634D-01,
     &1.4028D-01,1.2985D-01,1.2692D-01,1.3218D-01,1.4678D-01,1.7535D-01,
     &2.1492D-01,2.4628D-01,2.7582D-01,2.9349D-01,3.0215D-01,2.9865D-01,
     &2.8479D-01,2.6176D-01,2.3025D-01,1.9073D-01,1.4372D-01,9.0030D-02,
     &3.1431D-02,0.0000D+00,0.0000D+00,9.8561D-02,6.8571D-02,4.5400D-02,
     &2.8439D-02,1.6650D-02,8.2414D-03,4.3377D-03,1.8226D-03,0.0000D+00,
     &0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,5.8660D-01,5.7912D-01,
     &5.7170D-01,5.6412D-01,5.5660D-01,5.4858D-01,5.4040D-01,5.3194D-01,
     &5.2336D-01,5.1383D-01,5.0381D-01,4.9326D-01,4.8220D-01,4.7515D-01,
     &4.6719D-01,4.5756D-01,4.4719D-01,4.3619D-01,4.2441D-01,4.1376D-01,
     &4.0188D-01,3.8750D-01,3.7220D-01,3.5617D-01,3.3884D-01,3.2317D-01/
      DATA (XDEF_L(K),K=  571,  684) /
     &3.0561D-01,2.8413D-01,2.6132D-01,2.3801D-01,2.1667D-01,1.9794D-01,
     &1.7898D-01,1.5951D-01,1.4381D-01,1.3395D-01,1.3154D-01,1.3722D-01,
     &1.5183D-01,1.7978D-01,2.1726D-01,2.4615D-01,2.7227D-01,2.8668D-01,
     &2.9185D-01,2.8560D-01,2.6981D-01,2.4566D-01,2.1405D-01,1.7560D-01,
     &1.3093D-01,8.1317D-02,2.8821D-02,0.0000D+00,0.0000D+00,8.9016D-02,
     &6.1335D-02,4.0241D-02,2.4960D-02,1.4451D-02,6.9787D-03,3.6912D-03,
     &1.5320D-03,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &6.0621D-01,5.9821D-01,5.9043D-01,5.8253D-01,5.7470D-01,5.6625D-01,
     &5.5768D-01,5.4870D-01,5.3948D-01,5.2962D-01,5.1919D-01,5.0796D-01,
     &4.9620D-01,4.8867D-01,4.8027D-01,4.7003D-01,4.5907D-01,4.4740D-01,
     &4.3484D-01,4.2392D-01,4.1127D-01,3.9627D-01,3.8010D-01,3.6326D-01,
     &3.4524D-01,3.2900D-01,3.1064D-01,2.8853D-01,2.6510D-01,2.4135D-01,
     &2.1970D-01,2.0080D-01,1.8175D-01,1.6242D-01,1.4701D-01,1.3753D-01,
     &1.3572D-01,1.4160D-01,1.5623D-01,1.8343D-01,2.1902D-01,2.4571D-01,
     &2.6885D-01,2.8059D-01,2.8292D-01,2.7441D-01,2.5704D-01,2.3223D-01,
     &2.0062D-01,1.6317D-01,1.2079D-01,7.4733D-02,2.7461D-02,0.0000D+00,
     &0.0000D+00,8.1334D-02,5.5577D-02,3.6150D-02,2.2243D-02,1.2749D-02,
     &6.0264D-03,3.2009D-03,1.3143D-03,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.0000D+00,0.0000D+00,6.2581D-01,6.1778D-01,6.0953D-01,6.0134D-01/
      DATA (XDEF_L(K),K=  685,  798) /
     &5.9310D-01,5.8428D-01,5.7523D-01,5.6587D-01,5.5625D-01,5.4565D-01,
     &5.3457D-01,5.2280D-01,5.1030D-01,5.0236D-01,4.9350D-01,4.8267D-01,
     &4.7104D-01,4.5899D-01,4.4560D-01,4.3381D-01,4.2066D-01,4.0485D-01,
     &3.8801D-01,3.7047D-01,3.5165D-01,3.3476D-01,3.1574D-01,2.9293D-01,
     &2.6889D-01,2.4469D-01,2.2279D-01,2.0369D-01,1.8458D-01,1.6537D-01,
     &1.5025D-01,1.4125D-01,1.3980D-01,1.4589D-01,1.6046D-01,1.8686D-01,
     &2.2052D-01,2.4502D-01,2.6530D-01,2.7444D-01,2.7406D-01,2.6361D-01,
     &2.4491D-01,2.1954D-01,1.8819D-01,1.5193D-01,1.1170D-01,6.9146D-02,
     &2.6829D-02,0.0000D+00,0.0000D+00,7.4387D-02,5.0398D-02,3.2529D-02,
     &1.9840D-02,1.1260D-02,5.2109D-03,2.7796D-03,1.1291D-03,0.0000D+00,
     &0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,6.4510D-01,6.3663D-01,
     &6.2809D-01,6.1948D-01,6.1090D-01,6.0165D-01,5.9256D-01,5.8263D-01,
     &5.7237D-01,5.6121D-01,5.4960D-01,5.3710D-01,5.2390D-01,5.1555D-01,
     &5.0615D-01,4.9474D-01,4.8273D-01,4.6980D-01,4.5603D-01,4.4343D-01,
     &4.2983D-01,4.1325D-01,3.9561D-01,3.7731D-01,3.5765D-01,3.4017D-01,
     &3.2063D-01,2.9709D-01,2.7258D-01,2.4795D-01,2.2572D-01,2.0647D-01,
     &1.8735D-01,1.6824D-01,1.5339D-01,1.4470D-01,1.4366D-01,1.4990D-01,
     &1.6437D-01,1.8986D-01,2.2169D-01,2.4408D-01,2.6175D-01,2.6863D-01,
     &2.6585D-01,2.5363D-01,2.3397D-01,2.0813D-01,1.7714D-01,1.4205D-01/
      DATA (XDEF_L(K),K=  799,  912) /
     &1.0396D-01,6.4602D-02,2.6785D-02,0.0000D+00,0.0000D+00,6.8343D-02,
     &4.5962D-02,2.9434D-02,1.7812D-02,1.0015D-02,4.5458D-03,2.4331D-03,
     &9.7866D-04,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &6.6281D-01,6.5407D-01,6.4523D-01,6.3631D-01,6.2740D-01,6.1775D-01,
     &6.0821D-01,5.9770D-01,5.8724D-01,5.7535D-01,5.6321D-01,5.5021D-01,
     &5.3640D-01,5.2763D-01,5.1775D-01,5.0583D-01,4.9310D-01,4.7946D-01,
     &4.6520D-01,4.5225D-01,4.3811D-01,4.2074D-01,4.0247D-01,3.8355D-01,
     &3.6315D-01,3.4516D-01,3.2502D-01,3.0091D-01,2.7589D-01,2.5090D-01,
     &2.2842D-01,2.0903D-01,1.8987D-01,1.7087D-01,1.5631D-01,1.4790D-01,
     &1.4709D-01,1.5345D-01,1.6771D-01,1.9243D-01,2.2253D-01,2.4307D-01,
     &2.5846D-01,2.6327D-01,2.5857D-01,2.4493D-01,2.2441D-01,1.9832D-01,
     &1.6773D-01,1.3380D-01,9.7606D-02,6.1077D-02,2.7123D-02,4.1687D-04,
     &0.0000D+00,6.3316D-02,4.2290D-02,2.6899D-02,1.6166D-02,9.0143D-03,
     &4.0214D-03,2.1587D-03,8.6042D-04,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.0000D+00,0.0000D+00,6.8558D-01,6.7623D-01,6.6716D-01,6.5776D-01,
     &6.4840D-01,6.3825D-01,6.2778D-01,6.1697D-01,6.0589D-01,5.9350D-01,
     &5.8071D-01,5.6677D-01,5.5220D-01,5.4293D-01,5.3246D-01,5.1980D-01,
     &5.0630D-01,4.9221D-01,4.7690D-01,4.6348D-01,4.4839D-01,4.3024D-01,
     &4.1112D-01,3.9125D-01,3.7016D-01,3.5134D-01,3.3054D-01,3.0571D-01/
      DATA (XDEF_L(K),K=  913, 1026) /
     &2.8005D-01,2.5463D-01,2.3186D-01,2.1230D-01,1.9311D-01,1.7422D-01,
     &1.5985D-01,1.5187D-01,1.5138D-01,1.5783D-01,1.7178D-01,1.9543D-01,
     &2.2331D-01,2.4162D-01,2.5415D-01,2.5666D-01,2.4964D-01,2.3438D-01,
     &2.1293D-01,1.8681D-01,1.5680D-01,1.2430D-01,9.0488D-02,5.7352D-02,
     &2.7942D-02,7.0995D-03,2.4780D-03,5.7612D-02,3.8138D-02,2.4057D-02,
     &1.4329D-02,7.9111D-03,3.4566D-03,1.8603D-03,7.3347D-04,0.0000D+00,
     &0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,7.0709D-01,6.9744D-01,
     &6.8784D-01,6.7803D-01,6.6830D-01,6.5763D-01,6.4678D-01,6.3540D-01,
     &6.2360D-01,6.1071D-01,5.9715D-01,5.8240D-01,5.6710D-01,5.5722D-01,
     &5.4625D-01,5.3291D-01,5.1856D-01,5.0380D-01,4.8797D-01,4.7363D-01,
     &4.5801D-01,4.3900D-01,4.1917D-01,3.9846D-01,3.7656D-01,3.5717D-01,
     &3.3564D-01,3.1017D-01,2.8397D-01,2.5816D-01,2.3508D-01,2.1538D-01,
     &1.9615D-01,1.7737D-01,1.6324D-01,1.5559D-01,1.5535D-01,1.6175D-01,
     &1.7537D-01,1.9793D-01,2.2384D-01,2.4005D-01,2.5009D-01,2.5051D-01,
     &2.4150D-01,2.2495D-01,2.0291D-01,1.7668D-01,1.4739D-01,1.1625D-01,
     &8.4583D-02,5.4470D-02,2.9013D-02,1.3147D-02,1.4553D-02,5.2777D-02,
     &3.4672D-02,2.1686D-02,1.2821D-02,7.0105D-03,3.0093D-03,1.6226D-03,
     &6.3321D-04,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &7.2796D-01,7.1795D-01,7.0799D-01,6.9776D-01,6.8760D-01,6.7649D-01/
      DATA (XDEF_L(K),K= 1027, 1140) /
     &6.6523D-01,6.5299D-01,6.4099D-01,6.2720D-01,6.1289D-01,5.9763D-01,
     &5.8140D-01,5.7108D-01,5.5954D-01,5.4555D-01,5.3082D-01,5.1501D-01,
     &4.9841D-01,4.8352D-01,4.6718D-01,4.4758D-01,4.2678D-01,4.0543D-01,
     &3.8267D-01,3.6267D-01,3.4052D-01,3.1445D-01,2.8771D-01,2.6154D-01,
     &2.3817D-01,2.1835D-01,1.9910D-01,1.8043D-01,1.6662D-01,1.5905D-01,
     &1.5900D-01,1.6548D-01,1.7871D-01,2.0015D-01,2.2403D-01,2.3835D-01,
     &2.4610D-01,2.4469D-01,2.3394D-01,2.1634D-01,1.9372D-01,1.6761D-01,
     &1.3910D-01,1.0920D-01,7.9530D-02,5.2165D-02,3.0250D-02,1.8723D-02,
     &2.5275D-02,4.8575D-02,3.1676D-02,1.9677D-02,1.1540D-02,6.2533D-03,
     &2.6411D-03,1.4253D-03,5.5072D-04,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.0000D+00,0.0000D+00,7.4788D-01,7.3751D-01,7.2708D-01,7.1644D-01,
     &7.0580D-01,6.9430D-01,6.8256D-01,6.6975D-01,6.5712D-01,6.4276D-01,
     &6.2791D-01,6.1180D-01,5.9490D-01,5.8409D-01,5.7199D-01,5.5739D-01,
     &5.4166D-01,5.2544D-01,5.0821D-01,4.9288D-01,4.7590D-01,4.5544D-01,
     &4.3393D-01,4.1178D-01,3.8837D-01,3.6775D-01,3.4513D-01,3.1844D-01,
     &2.9125D-01,2.6472D-01,2.4110D-01,2.2115D-01,2.0189D-01,1.8330D-01,
     &1.6955D-01,1.6237D-01,1.6243D-01,1.6875D-01,1.8164D-01,2.0201D-01,
     &2.2410D-01,2.3665D-01,2.4236D-01,2.3927D-01,2.2710D-01,2.0852D-01,
     &1.8563D-01,1.5962D-01,1.3170D-01,1.0314D-01,7.5292D-02,5.0347D-02/
      DATA (XDEF_L(K),K= 1141, 1254) /
     &3.1513D-02,2.3688D-02,3.4520D-02,4.4988D-02,2.9140D-02,1.7975D-02,
     &1.0472D-02,5.6268D-03,2.3442D-03,1.2646D-03,4.8432D-04,0.0000D+00,
     &0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,7.6812D-01,7.5731D-01,
     &7.4653D-01,7.3551D-01,7.2440D-01,7.1234D-01,6.9989D-01,6.8692D-01,
     &6.7357D-01,6.5855D-01,6.4312D-01,6.2624D-01,6.0850D-01,5.9719D-01,
     &5.8457D-01,5.6934D-01,5.5297D-01,5.3626D-01,5.1802D-01,5.0223D-01,
     &4.8440D-01,4.6329D-01,4.4109D-01,4.1826D-01,3.9408D-01,3.7291D-01,
     &3.4966D-01,3.2243D-01,2.9475D-01,2.6790D-01,2.4406D-01,2.2399D-01,
     &2.0470D-01,1.8621D-01,1.7262D-01,1.6558D-01,1.6576D-01,1.7201D-01,
     &1.8441D-01,2.0372D-01,2.2403D-01,2.3482D-01,2.3856D-01,2.3398D-01,
     &2.2040D-01,2.0103D-01,1.7782D-01,1.5205D-01,1.2492D-01,9.7540D-02,
     &7.1452D-02,4.8817D-02,3.2832D-02,2.8412D-02,4.3068D-02,4.1684D-02,
     &2.6819D-02,1.6431D-02,9.5049D-03,5.0674D-03,2.0840D-03,1.1231D-03,
     &4.2643D-04,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &7.8709D-01,7.7617D-01,7.6509D-01,7.5353D-01,7.4210D-01,7.2955D-01,
     &7.1666D-01,7.0326D-01,6.8906D-01,6.7364D-01,6.5743D-01,6.3988D-01,
     &6.2140D-01,6.0962D-01,5.9645D-01,5.8083D-01,5.6382D-01,5.4630D-01,
     &5.2750D-01,5.1079D-01,4.9267D-01,4.7078D-01,4.4780D-01,4.2425D-01,
     &3.9948D-01,3.7773D-01,3.5398D-01,3.2619D-01,2.9811D-01,2.7093D-01/
      DATA (XDEF_L(K),K= 1255, 1368) /
     &2.4686D-01,2.2668D-01,2.0735D-01,1.8888D-01,1.7555D-01,1.6865D-01,
     &1.6887D-01,1.7500D-01,1.8693D-01,2.0522D-01,2.2377D-01,2.3300D-01,
     &2.3501D-01,2.2902D-01,2.1428D-01,1.9427D-01,1.7084D-01,1.4533D-01,
     &1.1889D-01,9.2655D-02,6.8174D-02,4.7575D-02,3.4123D-02,3.2605D-02,
     &5.0454D-02,3.8820D-02,2.4822D-02,1.5113D-02,8.6857D-03,4.5962D-03,
     &1.8704D-03,1.0050D-03,3.7856D-04,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.0000D+00,0.0000D+00,8.0606D-01,7.9455D-01,7.8312D-01,7.7128D-01,
     &7.5940D-01,7.4610D-01,7.3287D-01,7.1917D-01,7.0456D-01,6.8825D-01,
     &6.7140D-01,6.5313D-01,6.3390D-01,6.2170D-01,6.0798D-01,5.9180D-01,
     &5.7419D-01,5.5596D-01,5.3636D-01,5.1934D-01,5.0050D-01,4.7790D-01,
     &4.5436D-01,4.3012D-01,4.0458D-01,3.8238D-01,3.5808D-01,3.2984D-01,
     &3.0133D-01,2.7388D-01,2.4957D-01,2.2930D-01,2.0996D-01,1.9168D-01,
     &1.7832D-01,1.7159D-01,1.7177D-01,1.7770D-01,1.8921D-01,2.0651D-01,
     &2.2344D-01,2.3117D-01,2.3152D-01,2.2426D-01,2.0844D-01,1.8790D-01,
     &1.6440D-01,1.3914D-01,1.1342D-01,8.8280D-02,6.5276D-02,4.6578D-02,
     &3.5360D-02,3.6411D-02,5.6986D-02,3.6256D-02,2.3040D-02,1.3948D-02,
     &7.9676D-03,4.1856D-03,1.6876D-03,9.0394D-04,3.3789D-04,0.0000D+00,
     &0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,8.2409D-01,8.1223D-01,
     &8.0027D-01,7.8810D-01,7.7580D-01,7.6250D-01,7.4852D-01,7.3383D-01/
      DATA (XDEF_L(K),K= 1369, 1482) /
     &7.1879D-01,7.0216D-01,6.8466D-01,6.6571D-01,6.4580D-01,6.3303D-01,
     &6.1887D-01,6.0161D-01,5.8362D-01,5.6485D-01,5.4490D-01,5.2736D-01,
     &5.0788D-01,4.8465D-01,4.6048D-01,4.3549D-01,4.0949D-01,3.8678D-01,
     &3.6198D-01,3.3325D-01,3.0435D-01,2.7667D-01,2.5212D-01,2.3179D-01,
     &2.1241D-01,1.9410D-01,1.8093D-01,1.7428D-01,1.7445D-01,1.8022D-01,
     &1.9133D-01,2.0758D-01,2.2299D-01,2.2941D-01,2.2823D-01,2.1990D-01,
     &2.0319D-01,1.8211D-01,1.5852D-01,1.3371D-01,1.0856D-01,8.4430D-02,
     &6.2776D-02,4.5758D-02,3.6514D-02,3.9756D-02,6.2597D-02,3.4019D-02,
     &2.1502D-02,1.2943D-02,7.3506D-03,3.8366D-03,1.5351D-03,8.1923D-04,
     &3.0383D-04,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &8.4844D-01,8.3627D-01,8.2378D-01,8.1114D-01,7.9820D-01,7.8411D-01,
     &7.6977D-01,7.5436D-01,7.3871D-01,7.2101D-01,7.0269D-01,6.8280D-01,
     &6.6180D-01,6.4849D-01,6.3365D-01,6.1605D-01,5.9682D-01,5.7721D-01,
     &5.5628D-01,5.3805D-01,5.1772D-01,4.9378D-01,4.6868D-01,4.4295D-01,
     &4.1599D-01,3.9262D-01,3.6722D-01,3.3788D-01,3.0847D-01,2.8040D-01,
     &2.5562D-01,2.3513D-01,2.1572D-01,1.9746D-01,1.8447D-01,1.7787D-01,
     &1.7810D-01,1.8358D-01,1.9394D-01,2.0894D-01,2.2227D-01,2.2689D-01,
     &2.2385D-01,2.1408D-01,1.9620D-01,1.7461D-01,1.5108D-01,1.2667D-01,
     &1.0243D-01,7.9635D-02,5.9715D-02,4.4804D-02,3.7997D-02,4.3894D-02/
      DATA (XDEF_L(K),K= 1483, 1596) /
     &6.9391D-02,3.1240D-02,1.9603D-02,1.1712D-02,6.6036D-03,3.4150D-03,
     &1.3549D-03,7.1812D-04,2.6373D-04,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.0000D+00,0.0000D+00,8.7089D-01,8.5819D-01,8.4535D-01,8.3207D-01,
     &8.1860D-01,8.0424D-01,7.8877D-01,7.7320D-01,7.5642D-01,7.3822D-01,
     &7.1895D-01,6.9816D-01,6.7640D-01,6.6244D-01,6.4701D-01,6.2817D-01,
     &6.0860D-01,5.8841D-01,5.6672D-01,5.4767D-01,5.2667D-01,5.0182D-01,
     &4.7599D-01,4.4955D-01,4.2190D-01,3.9787D-01,3.7196D-01,3.4199D-01,
     &3.1220D-01,2.8382D-01,2.5874D-01,2.3816D-01,2.1874D-01,2.0063D-01,
     &1.8770D-01,1.8107D-01,1.8121D-01,1.8638D-01,1.9622D-01,2.0994D-01,
     &2.2156D-01,2.2456D-01,2.1986D-01,2.0892D-01,1.9015D-01,1.6817D-01,
     &1.4465D-01,1.2070D-01,9.7309D-02,7.5665D-02,5.7234D-02,4.4095D-02,
     &3.9289D-02,4.7307D-02,7.4739D-02,2.8958D-02,1.8046D-02,1.0716D-02,
     &6.0010D-03,3.0801D-03,1.2145D-03,6.3833D-04,2.3251D-04,0.0000D+00,
     &0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,8.9366D-01,8.8058D-01,
     &8.6727D-01,8.5353D-01,8.3950D-01,8.2436D-01,8.0890D-01,7.9205D-01,
     &7.7476D-01,7.5566D-01,7.3557D-01,7.1393D-01,6.9120D-01,6.7672D-01,
     &6.6059D-01,6.4145D-01,6.2086D-01,5.9962D-01,5.7716D-01,5.5756D-01,
     &5.3584D-01,5.1022D-01,4.8344D-01,4.5615D-01,4.2780D-01,4.0320D-01,
     &3.7671D-01,3.4621D-01,3.1594D-01,2.8727D-01,2.6196D-01,2.4126D-01/
      DATA (XDEF_L(K),K= 1597, 1710) /
     &2.2177D-01,2.0361D-01,1.9078D-01,1.8427D-01,1.8432D-01,1.8918D-01,
     &1.9834D-01,2.1079D-01,2.2065D-01,2.2210D-01,2.1587D-01,2.0383D-01,
     &1.8424D-01,1.6197D-01,1.3849D-01,1.1505D-01,9.2463D-02,7.1949D-02,
     &5.4952D-02,4.3474D-02,4.0525D-02,5.0376D-02,7.9517D-02,2.6835D-02,
     &1.6616D-02,9.8004D-03,5.4489D-03,2.7768D-03,1.0900D-03,5.6728D-04,
     &2.0489D-04,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &9.1643D-01,9.0298D-01,8.8901D-01,8.7472D-01,8.6030D-01,8.4449D-01,
     &8.2790D-01,8.1090D-01,7.9278D-01,7.7287D-01,7.5201D-01,7.2942D-01,
     &7.0580D-01,6.9067D-01,6.7395D-01,6.5357D-01,6.3264D-01,6.1082D-01,
     &5.8728D-01,5.6718D-01,5.4478D-01,5.1825D-01,4.9075D-01,4.6263D-01,
     &4.3360D-01,4.0844D-01,3.8138D-01,3.5032D-01,3.1963D-01,2.9065D-01,
     &2.6511D-01,2.4428D-01,2.2479D-01,2.0678D-01,1.9385D-01,1.8735D-01,
     &1.8722D-01,1.9179D-01,2.0029D-01,2.1158D-01,2.1961D-01,2.1971D-01,
     &2.1194D-01,1.9894D-01,1.7862D-01,1.5609D-01,1.3279D-01,1.0972D-01,
     &8.8007D-02,6.8578D-02,5.2905D-02,4.2942D-02,4.1624D-02,5.3065D-02,
     &8.3506D-02,2.4920D-02,1.5334D-02,8.9876D-03,4.9653D-03,2.5112D-03,
     &9.8300D-04,5.0629D-04,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,
     &0.0000D+00,0.0000D+00,9.3762D-01,9.2325D-01,9.0916D-01,8.9432D-01,
     &8.7930D-01,8.6312D-01,8.4579D-01,8.2807D-01,8.0954D-01,7.8866D-01/
      DATA (XDEF_L(K),K= 1711, 1824) /
     &7.6704D-01,7.4360D-01,7.1911D-01,7.0343D-01,6.8612D-01,6.6512D-01,
     &6.4349D-01,6.2048D-01,5.9676D-01,5.7574D-01,5.5261D-01,5.2556D-01,
     &4.9731D-01,4.6862D-01,4.3881D-01,4.1318D-01,3.8556D-01,3.5408D-01,
     &3.2299D-01,2.9375D-01,2.6794D-01,2.4706D-01,2.2744D-01,2.0939D-01,
     &1.9662D-01,1.9016D-01,1.8990D-01,1.9412D-01,2.0192D-01,2.1208D-01,
     &2.1863D-01,2.1745D-01,2.0845D-01,1.9458D-01,1.7365D-01,1.5094D-01,
     &1.2783D-01,1.0526D-01,8.4228D-02,6.5746D-02,5.1203D-02,4.2521D-02,
     &4.2531D-02,5.5238D-02,8.6619D-02,2.3321D-02,1.4266D-02,8.3142D-03,
     &4.5684D-03,2.2945D-03,8.9721D-04,4.5700D-04,0.0000D+00,0.0000D+00,
     &0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,9.5912D-01,9.4446D-01,
     &9.2967D-01,9.1446D-01,8.9890D-01,8.8176D-01,8.6424D-01,8.4567D-01,
     &8.2630D-01,8.0492D-01,7.8242D-01,7.5817D-01,7.3271D-01,7.1653D-01,
     &6.9849D-01,6.7725D-01,6.5433D-01,6.3091D-01,6.0625D-01,5.8456D-01,
     &5.6088D-01,5.3305D-01,5.0402D-01,4.7461D-01,4.4411D-01,4.1800D-01,
     &3.8988D-01,3.5790D-01,3.2644D-01,2.9690D-01,2.7087D-01,2.4987D-01,
     &2.3039D-01,2.1219D-01,1.9955D-01,1.9298D-01,1.9248D-01,1.9636D-01,
     &2.0355D-01,2.1258D-01,2.1752D-01,2.1512D-01,2.0490D-01,1.9021D-01,
     &1.6876D-01,1.4586D-01,1.2296D-01,1.0090D-01,8.0587D-02,6.3034D-02,
     &4.9591D-02,4.2122D-02,4.3355D-02,5.7203D-02,8.9336D-02,2.1802D-02/
      DATA (XDEF_L(K),K= 1825, 1836) /
     &1.3258D-02,7.6843D-03,4.1967D-03,2.0952D-03,8.1932D-04,4.1202D-04,
     &0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00,0.0000D+00/
      DATA (XUDF_L(K),K=    1,  114) /
     &1.8987D-02,1.9947D-02,2.0980D-02,2.2068D-02,2.3225D-02,2.4540D-02,
     &2.5957D-02,2.7526D-02,2.9229D-02,3.1232D-02,3.3453D-02,3.6003D-02,
     &3.8855D-02,4.0763D-02,4.2980D-02,4.5778D-02,4.8895D-02,5.2320D-02,
     &5.6174D-02,5.9765D-02,6.3980D-02,6.9315D-02,7.5299D-02,8.1888D-02,
     &8.9292D-02,9.6162D-02,1.0414D-01,1.1410D-01,1.2505D-01,1.3674D-01,
     &1.4937D-01,1.6060D-01,1.7296D-01,1.8730D-01,2.0166D-01,2.1531D-01,
     &2.2821D-01,2.3833D-01,2.4848D-01,2.6049D-01,2.7586D-01,2.9166D-01,
     &3.1456D-01,3.3942D-01,3.7230D-01,4.0597D-01,4.3921D-01,4.7071D-01,
     &4.9846D-01,5.2057D-01,5.3433D-01,5.3610D-01,5.2141D-01,4.8433D-01,
     &4.1719D-01,6.3794D-01,6.7411D-01,7.2040D-01,7.8812D-01,8.9495D-01,
     &1.0702D+00,1.3629D+00,1.8763D+00,2.8399D+00,4.8968D+00,1.0506D+01,
     &3.7793D+01,0.0000D+00,3.1111D-02,3.2336D-02,3.3580D-02,3.4906D-02,
     &3.6247D-02,3.7773D-02,3.9337D-02,4.1056D-02,4.2876D-02,4.5001D-02,
     &4.7299D-02,4.9897D-02,5.2761D-02,5.4666D-02,5.6867D-02,5.9620D-02,
     &6.2679D-02,6.6018D-02,6.9775D-02,7.3275D-02,7.7353D-02,8.2522D-02,
     &8.8327D-02,9.4694D-02,1.0184D-01,1.0846D-01,1.1615D-01,1.2575D-01,
     &1.3628D-01,1.4752D-01,1.5964D-01,1.7036D-01,1.8215D-01,1.9580D-01,
     &2.0933D-01,2.2213D-01,2.3411D-01,2.4341D-01,2.5275D-01,2.6387D-01,
     &2.7831D-01,2.9333D-01,3.1510D-01,3.3876D-01,3.6995D-01,4.0170D-01/
      DATA (XUDF_L(K),K=  115,  228) /
     &4.3298D-01,4.6172D-01,4.8742D-01,5.0700D-01,5.1856D-01,5.1873D-01,
     &5.0352D-01,4.6746D-01,4.0418D-01,6.1801D-01,6.5339D-01,6.9923D-01,
     &7.6627D-01,8.7125D-01,1.0408D+00,1.3199D+00,1.8020D+00,2.6920D+00,
     &4.5574D+00,9.5310D+00,3.2877D+01,0.0000D+00,5.1176D-02,5.2640D-02,
     &5.4100D-02,5.5603D-02,5.7095D-02,5.8737D-02,6.0416D-02,6.2154D-02,
     &6.4016D-02,6.6046D-02,6.8273D-02,7.0765D-02,7.3444D-02,7.5182D-02,
     &7.7263D-02,7.9781D-02,8.2626D-02,8.5707D-02,8.9176D-02,9.2402D-02,
     &9.6182D-02,1.0098D-01,1.0635D-01,1.1227D-01,1.1893D-01,1.2513D-01,
     &1.3230D-01,1.4128D-01,1.5115D-01,1.6164D-01,1.7300D-01,1.8301D-01,
     &1.9397D-01,2.0660D-01,2.1907D-01,2.3072D-01,2.4154D-01,2.4985D-01,
     &2.5817D-01,2.6810D-01,2.8136D-01,2.9535D-01,3.1585D-01,3.3824D-01,
     &3.6743D-01,3.9701D-01,4.2565D-01,4.5205D-01,4.7460D-01,4.9184D-01,
     &5.0110D-01,4.9954D-01,4.8363D-01,4.4878D-01,3.8940D-01,5.9452D-01,
     &6.2820D-01,6.7181D-01,7.3612D-01,8.3598D-01,9.9560D-01,1.2543D+00,
     &1.6953D+00,2.4947D+00,4.1415D+00,8.4275D+00,2.7797D+01,0.0000D+00,
     &8.6266D-02,8.7847D-02,8.9380D-02,9.0869D-02,9.2337D-02,9.3826D-02,
     &9.5315D-02,9.6842D-02,9.8333D-02,1.0003D-01,1.0178D-01,1.0370D-01,
     &1.0575D-01,1.0710D-01,1.0872D-01,1.1075D-01,1.1295D-01,1.1538D-01,
     &1.1821D-01,1.2088D-01,1.2396D-01,1.2796D-01,1.3252D-01,1.3756D-01/
      DATA (XUDF_L(K),K=  229,  342) /
     &1.4331D-01,1.4870D-01,1.5500D-01,1.6291D-01,1.7166D-01,1.8100D-01,
     &1.9111D-01,2.0002D-01,2.0977D-01,2.2095D-01,2.3189D-01,2.4200D-01,
     &2.5123D-01,2.5821D-01,2.6512D-01,2.7351D-01,2.8514D-01,2.9789D-01,
     &3.1683D-01,3.3731D-01,3.6424D-01,3.9124D-01,4.1697D-01,4.4030D-01,
     &4.6002D-01,4.7419D-01,4.8085D-01,4.7740D-01,4.6086D-01,4.2728D-01,
     &3.7241D-01,5.6656D-01,5.9684D-01,6.3694D-01,6.9622D-01,7.8804D-01,
     &9.3343D-01,1.1653D+00,1.5545D+00,2.2504D+00,3.6537D+00,7.2124D+00,
     &2.2653D+01,0.0000D+00,1.4838D-01,1.4960D-01,1.5068D-01,1.5161D-01,
     &1.5242D-01,1.5316D-01,1.5373D-01,1.5426D-01,1.5470D-01,1.5511D-01,
     &1.5554D-01,1.5602D-01,1.5660D-01,1.5698D-01,1.5750D-01,1.5830D-01,
     &1.5923D-01,1.6034D-01,1.6181D-01,1.6324D-01,1.6509D-01,1.6746D-01,
     &1.7054D-01,1.7402D-01,1.7811D-01,1.8208D-01,1.8687D-01,1.9296D-01,
     &1.9986D-01,2.0734D-01,2.1554D-01,2.2281D-01,2.3075D-01,2.3983D-01,
     &2.4863D-01,2.5660D-01,2.6366D-01,2.6883D-01,2.7387D-01,2.8026D-01,
     &2.8982D-01,3.0088D-01,3.1780D-01,3.3626D-01,3.6021D-01,3.8399D-01,
     &4.0666D-01,4.2682D-01,4.4278D-01,4.5386D-01,4.5774D-01,4.5230D-01,
     &4.3509D-01,4.0314D-01,3.5321D-01,5.3325D-01,5.5916D-01,5.9448D-01,
     &6.4707D-01,7.2797D-01,8.5557D-01,1.0563D+00,1.3882D+00,1.9717D+00,
     &3.1223D+00,5.9601D+00,1.7750D+01,0.0000D+00,2.3139D-01,2.3138D-01/
      DATA (XUDF_L(K),K=  343,  456) /
     &2.3120D-01,2.3076D-01,2.3006D-01,2.2907D-01,2.2788D-01,2.2645D-01,
     &2.2489D-01,2.2308D-01,2.2120D-01,2.1929D-01,2.1743D-01,2.1630D-01,
     &2.1526D-01,2.1411D-01,2.1311D-01,2.1231D-01,2.1171D-01,2.1148D-01,
     &2.1150D-01,2.1182D-01,2.1271D-01,2.1412D-01,2.1601D-01,2.1822D-01,
     &2.2096D-01,2.2496D-01,2.2961D-01,2.3481D-01,2.4086D-01,2.4622D-01,
     &2.5214D-01,2.5891D-01,2.6537D-01,2.7104D-01,2.7588D-01,2.7922D-01,
     &2.8235D-01,2.8664D-01,2.9413D-01,3.0352D-01,3.1845D-01,3.3481D-01,
     &3.5617D-01,3.7737D-01,3.9689D-01,4.1403D-01,4.2736D-01,4.3558D-01,
     &4.3712D-01,4.3016D-01,4.1245D-01,3.8197D-01,3.3645D-01,5.0322D-01,
     &5.2507D-01,5.5559D-01,6.0172D-01,6.7286D-01,7.8413D-01,9.5797D-01,
     &1.2422D+00,1.7341D+00,2.6883D+00,4.9868D+00,1.4177D+01,0.0000D+00,
     &3.6389D-01,3.6098D-01,3.5780D-01,3.5400D-01,3.5016D-01,3.4553D-01,
     &3.4044D-01,3.3521D-01,3.2971D-01,3.2369D-01,3.1755D-01,3.1120D-01,
     &3.0494D-01,3.0120D-01,2.9724D-01,2.9287D-01,2.8855D-01,2.8449D-01,
     &2.8072D-01,2.7770D-01,2.7469D-01,2.7175D-01,2.6933D-01,2.6740D-01,
     &2.6613D-01,2.6556D-01,2.6563D-01,2.6631D-01,2.6763D-01,2.6975D-01,
     &2.7268D-01,2.7539D-01,2.7857D-01,2.8224D-01,2.8565D-01,2.8841D-01,
     &2.9040D-01,2.9139D-01,2.9220D-01,2.9395D-01,2.9888D-01,3.0633D-01,
     &3.1877D-01,3.3296D-01,3.5147D-01,3.6947D-01,3.8604D-01,3.9986D-01/
      DATA (XUDF_L(K),K=  457,  570) /
     &4.1008D-01,4.1548D-01,4.1467D-01,4.0620D-01,3.8830D-01,3.5965D-01,
     &3.1902D-01,4.7020D-01,4.8772D-01,5.1303D-01,5.5185D-01,6.1224D-01,
     &7.0699D-01,8.5323D-01,1.0903D+00,1.4950D+00,2.2640D+00,4.0723D+00,
     &0.0000D+00,0.0000D+00,5.2666D-01,5.1909D-01,5.1100D-01,5.0238D-01,
     &4.9333D-01,4.8312D-01,4.7293D-01,4.6180D-01,4.5066D-01,4.3890D-01,
     &4.2692D-01,4.1467D-01,4.0262D-01,3.9542D-01,3.8784D-01,3.7925D-01,
     &3.7080D-01,3.6267D-01,3.5482D-01,3.4841D-01,3.4190D-01,3.3492D-01,
     &3.2852D-01,3.2287D-01,3.1768D-01,3.1409D-01,3.1066D-01,3.0785D-01,
     &3.0564D-01,3.0446D-01,3.0380D-01,3.0388D-01,3.0402D-01,3.0458D-01,
     &3.0488D-01,3.0475D-01,3.0386D-01,3.0263D-01,3.0116D-01,3.0045D-01,
     &3.0296D-01,3.0852D-01,3.1888D-01,3.3085D-01,3.4677D-01,3.6222D-01,
     &3.7600D-01,3.8707D-01,3.9488D-01,3.9799D-01,3.9530D-01,3.8568D-01,
     &3.6791D-01,3.4080D-01,3.0424D-01,4.4195D-01,4.5570D-01,4.7648D-01,
     &5.0935D-01,5.6099D-01,6.4225D-01,7.6680D-01,9.6736D-01,1.3053D+00,
     &1.9393D+00,3.3976D+00,0.0000D+00,0.0000D+00,7.4015D-01,7.2498D-01,
     &7.0940D-01,6.9297D-01,6.7620D-01,6.5800D-01,6.3935D-01,6.2047D-01,
     &6.0114D-01,5.8076D-01,5.6065D-01,5.4030D-01,5.2035D-01,5.0839D-01,
     &4.9583D-01,4.8167D-01,4.6773D-01,4.5434D-01,4.4113D-01,4.3035D-01,
     &4.1922D-01,4.0719D-01,3.9582D-01,3.8536D-01,3.7557D-01,3.6805D-01/
      DATA (XUDF_L(K),K=  571,  684) /
     &3.6079D-01,3.5336D-01,3.4710D-01,3.4173D-01,3.3719D-01,3.3400D-01,
     &3.3124D-01,3.2819D-01,3.2494D-01,3.2158D-01,3.1765D-01,3.1400D-01,
     &3.1011D-01,3.0684D-01,3.0682D-01,3.1046D-01,3.1856D-01,3.2861D-01,
     &3.4189D-01,3.5475D-01,3.6597D-01,3.7463D-01,3.8003D-01,3.8108D-01,
     &3.7681D-01,3.6631D-01,3.4865D-01,3.2327D-01,2.9078D-01,4.1488D-01,
     &4.2529D-01,4.4193D-01,4.6945D-01,5.1322D-01,5.8236D-01,6.8846D-01,
     &8.5739D-01,1.1394D+00,1.6617D+00,2.8395D+00,0.0000D+00,0.0000D+00,
     &9.8501D-01,9.5975D-01,9.3420D-01,9.0757D-01,8.8092D-01,8.5237D-01,
     &8.2383D-01,7.9445D-01,7.6556D-01,7.3524D-01,7.0484D-01,6.7495D-01,
     &6.4547D-01,6.2798D-01,6.0969D-01,5.8904D-01,5.6882D-01,5.4932D-01,
     &5.3014D-01,5.1443D-01,4.9826D-01,4.8058D-01,4.6380D-01,4.4815D-01,
     &4.3330D-01,4.2167D-01,4.1020D-01,3.9827D-01,3.8748D-01,3.7784D-01,
     &3.6931D-01,3.6303D-01,3.5669D-01,3.4992D-01,3.4358D-01,3.3710D-01,
     &3.3025D-01,3.2429D-01,3.1817D-01,3.1242D-01,3.1001D-01,3.1195D-01,
     &3.1802D-01,3.2610D-01,3.3719D-01,3.4770D-01,3.5674D-01,3.6357D-01,
     &3.6695D-01,3.6631D-01,3.6075D-01,3.4960D-01,3.3214D-01,3.0855D-01,
     &2.7931D-01,3.9198D-01,3.9931D-01,4.1263D-01,4.3550D-01,4.7310D-01,
     &5.3259D-01,6.2375D-01,7.6876D-01,1.0087D+00,1.4464D+00,2.4185D+00,
     &0.0000D+00,0.0000D+00,1.2917D+00,1.2523D+00,1.2128D+00,1.1722D+00/
      DATA (XUDF_L(K),K=  685,  798) /
     &1.1321D+00,1.0894D+00,1.0473D+00,1.0044D+00,9.6262D-01,9.1838D-01,
     &8.7565D-01,8.3283D-01,7.9186D-01,7.6734D-01,7.4146D-01,7.1300D-01,
     &6.8484D-01,6.5787D-01,6.3134D-01,6.0963D-01,5.8730D-01,5.6294D-01,
     &5.3947D-01,5.1767D-01,4.9689D-01,4.8039D-01,4.6398D-01,4.4675D-01,
     &4.3087D-01,4.1650D-01,4.0371D-01,3.9342D-01,3.8361D-01,3.7293D-01,
     &3.6284D-01,3.5305D-01,3.4307D-01,3.3468D-01,3.2613D-01,3.1788D-01,
     &3.1306D-01,3.1309D-01,3.1715D-01,3.2346D-01,3.3232D-01,3.4066D-01,
     &3.4779D-01,3.5251D-01,3.5401D-01,3.5184D-01,3.4519D-01,3.3347D-01,
     &3.1650D-01,2.9433D-01,2.6872D-01,3.6968D-01,3.7446D-01,3.8477D-01,
     &4.0368D-01,4.3551D-01,4.8654D-01,5.6457D-01,6.8832D-01,8.9135D-01,
     &1.2583D+00,2.0601D+00,0.0000D+00,0.0000D+00,1.6499D+00,1.5928D+00,
     &1.5356D+00,1.4773D+00,1.4202D+00,1.3601D+00,1.3009D+00,1.2413D+00,
     &1.1836D+00,1.1235D+00,1.0650D+00,1.0076D+00,9.5212D-01,9.1919D-01,
     &8.8569D-01,8.4733D-01,8.1006D-01,7.7436D-01,7.3955D-01,7.1104D-01,
     &6.8173D-01,6.4966D-01,6.1893D-01,5.9026D-01,5.6287D-01,5.4114D-01,
     &5.1941D-01,4.9621D-01,4.7490D-01,4.5564D-01,4.3786D-01,4.2408D-01,
     &4.1024D-01,3.9562D-01,3.8175D-01,3.6853D-01,3.5541D-01,3.4455D-01,
     &3.3366D-01,3.2286D-01,3.1565D-01,3.1397D-01,3.1618D-01,3.2069D-01,
     &3.2744D-01,3.3383D-01,3.3911D-01,3.4194D-01,3.4194D-01,3.3844D-01/
      DATA (XUDF_L(K),K=  799,  912) /
     &3.3088D-01,3.1887D-01,3.0224D-01,2.8177D-01,2.5901D-01,3.4945D-01,
     &3.5200D-01,3.5959D-01,3.7518D-01,4.0212D-01,4.4590D-01,5.1305D-01,
     &6.1934D-01,7.9273D-01,1.1025D+00,1.7693D+00,0.0000D+00,0.0000D+00,
     &2.0413D+00,1.9626D+00,1.8840D+00,1.8053D+00,1.7284D+00,1.6480D+00,
     &1.5697D+00,1.4911D+00,1.4157D+00,1.3375D+00,1.2620D+00,1.1875D+00,
     &1.1168D+00,1.0751D+00,1.0321D+00,9.8410D-01,9.3682D-01,8.9196D-01,
     &8.4816D-01,8.1245D-01,7.7582D-01,7.3576D-01,6.9745D-01,6.6154D-01,
     &6.2742D-01,6.0036D-01,5.7319D-01,5.4409D-01,5.1721D-01,4.9291D-01,
     &4.7049D-01,4.5284D-01,4.3541D-01,4.1671D-01,3.9926D-01,3.8274D-01,
     &3.6660D-01,3.5348D-01,3.4035D-01,3.2727D-01,3.1788D-01,3.1459D-01,
     &3.1499D-01,3.1792D-01,3.2291D-01,3.2764D-01,3.3124D-01,3.3250D-01,
     &3.3120D-01,3.2663D-01,3.1834D-01,3.0608D-01,2.8998D-01,2.7085D-01,
     &2.5085D-01,3.3191D-01,3.3258D-01,3.3808D-01,3.5072D-01,3.7379D-01,
     &4.1182D-01,4.7005D-01,5.6257D-01,7.1233D-01,9.7788D-01,1.5412D+00,
     &0.0000D+00,0.0000D+00,2.6325D+00,2.5188D+00,2.4060D+00,2.2942D+00,
     &2.1863D+00,2.0740D+00,1.9650D+00,1.8571D+00,1.7537D+00,1.6473D+00,
     &1.5453D+00,1.4458D+00,1.3515D+00,1.2965D+00,1.2394D+00,1.1767D+00,
     &1.1150D+00,1.0560D+00,9.9927D-01,9.5301D-01,9.0565D-01,8.5400D-01,
     &8.0462D-01,7.5858D-01,7.1481D-01,6.7994D-01,6.4502D-01,6.0799D-01/
      DATA (XUDF_L(K),K=  913, 1026) /
     &5.7349D-01,5.4206D-01,5.1299D-01,4.9028D-01,4.6789D-01,4.4387D-01,
     &4.2168D-01,4.0096D-01,3.8070D-01,3.6457D-01,3.4857D-01,3.3249D-01,
     &3.2026D-01,3.1503D-01,3.1326D-01,3.1423D-01,3.1703D-01,3.1974D-01,
     &3.2120D-01,3.2086D-01,3.1799D-01,3.1221D-01,3.0315D-01,2.9072D-01,
     &2.7522D-01,2.5796D-01,2.4114D-01,3.1079D-01,3.0956D-01,3.1267D-01,
     &3.2223D-01,3.4089D-01,3.7246D-01,4.2134D-01,4.9853D-01,6.2305D-01,
     &8.4191D-01,1.2983D+00,0.0000D+00,0.0000D+00,3.2997D+00,3.1427D+00,
     &2.9900D+00,2.8374D+00,2.6927D+00,2.5421D+00,2.3973D+00,2.2549D+00,
     &2.1191D+00,1.9809D+00,1.8488D+00,1.7209D+00,1.6001D+00,1.5300D+00,
     &1.4576D+00,1.3771D+00,1.2999D+00,1.2268D+00,1.1551D+00,1.0975D+00,
     &1.0385D+00,9.7437D-01,9.1327D-01,8.5649D-01,8.0236D-01,7.5952D-01,
     &7.1667D-01,6.7091D-01,6.2847D-01,5.9005D-01,5.5422D-01,5.2636D-01,
     &4.9890D-01,4.6976D-01,4.4269D-01,4.1752D-01,3.9377D-01,3.7477D-01,
     &3.5594D-01,3.3710D-01,3.2226D-01,3.1511D-01,3.1131D-01,3.1067D-01,
     &3.1132D-01,3.1227D-01,3.1198D-01,3.1021D-01,3.0606D-01,2.9926D-01,
     &2.8958D-01,2.7716D-01,2.6233D-01,2.4655D-01,2.3275D-01,2.9229D-01,
     &2.8941D-01,2.9061D-01,2.9753D-01,3.1273D-01,3.3909D-01,3.8034D-01,
     &4.4548D-01,5.5028D-01,7.3256D-01,1.1074D+00,0.0000D+00,0.0000D+00,
     &4.0557D+00,3.8486D+00,3.6460D+00,3.4480D+00,3.2579D+00,3.0626D+00/
      DATA (XUDF_L(K),K= 1027, 1140) /
     &2.8756D+00,2.6929D+00,2.5196D+00,2.3441D+00,2.1778D+00,2.0170D+00,
     &1.8670D+00,1.7797D+00,1.6902D+00,1.5909D+00,1.4960D+00,1.4058D+00,
     &1.3191D+00,1.2484D+00,1.1764D+00,1.0991D+00,1.0253D+00,9.5689D-01,
     &8.9197D-01,8.4046D-01,7.8904D-01,7.3442D-01,6.8367D-01,6.3780D-01,
     &5.9520D-01,5.6218D-01,5.2934D-01,4.9500D-01,4.6300D-01,4.3370D-01,
     &4.0611D-01,3.8431D-01,3.6284D-01,3.4121D-01,3.2389D-01,3.1494D-01,
     &3.0926D-01,3.0697D-01,3.0594D-01,3.0501D-01,3.0330D-01,3.0019D-01,
     &2.9492D-01,2.8734D-01,2.7718D-01,2.6476D-01,2.5057D-01,2.3646D-01,
     &2.2503D-01,2.7558D-01,2.7132D-01,2.7089D-01,2.7569D-01,2.8794D-01,
     &3.1000D-01,3.4491D-01,4.0016D-01,4.8886D-01,6.4191D-01,9.5232D-01,
     &0.0000D+00,0.0000D+00,4.8799D+00,4.6116D+00,4.3560D+00,4.1035D+00,
     &3.8608D+00,3.6163D+00,3.3822D+00,3.1557D+00,2.9412D+00,2.7247D+00,
     &2.5209D+00,2.3248D+00,2.1421D+00,2.0368D+00,1.9287D+00,1.8094D+00,
     &1.6955D+00,1.5877D+00,1.4841D+00,1.4003D+00,1.3154D+00,1.2237D+00,
     &1.1368D+00,1.0563D+00,9.8015D-01,9.2005D-01,8.5978D-01,7.9615D-01,
     &7.3715D-01,6.8369D-01,6.3441D-01,5.9609D-01,5.5830D-01,5.1865D-01,
     &4.8192D-01,4.4872D-01,4.1747D-01,3.9300D-01,3.6895D-01,3.4483D-01,
     &3.2508D-01,3.1459D-01,3.0709D-01,3.0328D-01,3.0056D-01,2.9840D-01,
     &2.9543D-01,2.9107D-01,2.8485D-01,2.7655D-01,2.6610D-01,2.5368D-01/
      DATA (XUDF_L(K),K= 1141, 1254) /
     &2.4019D-01,2.2736D-01,2.1837D-01,2.6080D-01,2.5542D-01,2.5362D-01,
     &2.5693D-01,2.6661D-01,2.8505D-01,3.1490D-01,3.6226D-01,4.3798D-01,
     &5.6769D-01,8.2836D-01,0.0000D+00,0.0000D+00,5.8340D+00,5.4940D+00,
     &5.1700D+00,4.8532D+00,4.5515D+00,4.2463D+00,3.9559D+00,3.6752D+00,
     &3.4138D+00,3.1496D+00,2.9022D+00,2.6648D+00,2.4450D+00,2.3189D+00,
     &2.1896D+00,2.0476D+00,1.9120D+00,1.7843D+00,1.6621D+00,1.5639D+00,
     &1.4648D+00,1.3569D+00,1.2556D+00,1.1618D+00,1.0734D+00,1.0037D+00,
     &9.3416D-01,8.6065D-01,7.9257D-01,7.3145D-01,6.7463D-01,6.3082D-01,
     &5.8786D-01,5.4262D-01,5.0118D-01,4.6374D-01,4.2883D-01,4.0146D-01,
     &3.7490D-01,3.4814D-01,3.2612D-01,3.1397D-01,3.0482D-01,2.9958D-01,
     &2.9536D-01,2.9178D-01,2.8756D-01,2.8208D-01,2.7504D-01,2.6611D-01,
     &2.5539D-01,2.4319D-01,2.3031D-01,2.1877D-01,2.1195D-01,2.4673D-01,
     &2.4036D-01,2.3746D-01,2.3912D-01,2.4677D-01,2.6223D-01,2.8748D-01,
     &3.2792D-01,3.9255D-01,5.0271D-01,7.2095D-01,0.0000D+00,0.0000D+00,
     &6.8578D+00,6.4388D+00,6.0380D+00,5.6501D+00,5.2825D+00,4.9103D+00,
     &4.5613D+00,4.2230D+00,3.9070D+00,3.5911D+00,3.2966D+00,3.0156D+00,
     &2.7567D+00,2.6078D+00,2.4563D+00,2.2905D+00,2.1319D+00,1.9837D+00,
     &1.8421D+00,1.7287D+00,1.6141D+00,1.4902D+00,1.3730D+00,1.2663D+00,
     &1.1652D+00,1.0858D+00,1.0067D+00,9.2337D-01,8.4648D-01,7.7710D-01/
      DATA (XUDF_L(K),K= 1255, 1368) /
     &7.1333D-01,6.6392D-01,6.1566D-01,5.6531D-01,5.1904D-01,4.7761D-01,
     &4.3908D-01,4.0927D-01,3.8022D-01,3.5109D-01,3.2686D-01,3.1318D-01,
     &3.0244D-01,2.9602D-01,2.9031D-01,2.8538D-01,2.8024D-01,2.7382D-01,
     &2.6607D-01,2.5668D-01,2.4571D-01,2.3364D-01,2.2155D-01,2.1116D-01,
     &2.0617D-01,2.3421D-01,2.2704D-01,2.2320D-01,2.2366D-01,2.2952D-01,
     &2.4241D-01,2.6402D-01,2.9884D-01,3.5437D-01,4.4860D-01,6.3331D-01,
     &0.0000D+00,0.0000D+00,7.9784D+00,7.4673D+00,6.9820D+00,6.5121D+00,
     &6.0712D+00,5.6250D+00,5.2080D+00,4.8065D+00,4.4309D+00,4.0590D+00,
     &3.7131D+00,3.3843D+00,3.0816D+00,2.9094D+00,2.7332D+00,2.5420D+00,
     &2.3595D+00,2.1895D+00,2.0271D+00,1.8966D+00,1.7658D+00,1.6248D+00,
     &1.4933D+00,1.3718D+00,1.2579D+00,1.1683D+00,1.0795D+00,9.8589D-01,
     &8.9996D-01,8.2253D-01,7.5153D-01,6.9648D-01,6.4287D-01,5.8736D-01,
     &5.3655D-01,4.9109D-01,4.4891D-01,4.1655D-01,3.8518D-01,3.5367D-01,
     &3.2738D-01,3.1221D-01,3.0006D-01,2.9246D-01,2.8544D-01,2.7940D-01,
     &2.7319D-01,2.6601D-01,2.5763D-01,2.4782D-01,2.3676D-01,2.2486D-01,
     &2.1329D-01,2.0405D-01,2.0083D-01,2.2267D-01,2.1489D-01,2.1027D-01,
     &2.0967D-01,2.1409D-01,2.2473D-01,2.4320D-01,2.7316D-01,3.2113D-01,
     &4.0209D-01,5.5899D-01,0.0000D+00,0.0000D+00,9.1575D+00,8.5458D+00,
     &7.9700D+00,7.4123D+00,6.8876D+00,6.3653D+00,5.8736D+00,5.4042D+00/
      DATA (XUDF_L(K),K= 1369, 1482) /
     &4.9684D+00,4.5359D+00,4.1366D+00,3.7576D+00,3.4110D+00,3.2138D+00,
     &3.0122D+00,2.7943D+00,2.5871D+00,2.3944D+00,2.2102D+00,2.0646D+00,
     &1.9163D+00,1.7581D+00,1.6109D+00,1.4753D+00,1.3483D+00,1.2486D+00,
     &1.1500D+00,1.0462D+00,9.5130D-01,8.6585D-01,7.8770D-01,7.2741D-01,
     &6.6891D-01,6.0781D-01,5.5266D-01,5.0342D-01,4.5788D-01,4.2322D-01,
     &3.8960D-01,3.5594D-01,3.2768D-01,3.1125D-01,2.9779D-01,2.8890D-01,
     &2.8091D-01,2.7385D-01,2.6670D-01,2.5886D-01,2.4989D-01,2.3976D-01,
     &2.2861D-01,2.1703D-01,2.0604D-01,1.9777D-01,1.9598D-01,2.1238D-01,
     &2.0408D-01,1.9879D-01,1.9735D-01,2.0048D-01,2.0933D-01,2.2523D-01,
     &2.5120D-01,2.9296D-01,3.6305D-01,4.9711D-01,0.0000D+00,0.0000D+00,
     &1.0956D+01,1.0188D+01,9.4660D+00,8.7704D+00,8.1209D+00,7.4727D+00,
     &6.8721D+00,6.2972D+00,5.7646D+00,5.2434D+00,4.7595D+00,4.3051D+00,
     &3.8911D+00,3.6559D+00,3.4174D+00,3.1598D+00,2.9153D+00,2.6889D+00,
     &2.4732D+00,2.3031D+00,2.1311D+00,1.9475D+00,1.7771D+00,1.6202D+00,
     &1.4748D+00,1.3609D+00,1.2481D+00,1.1301D+00,1.0222D+00,9.2549D-01,
     &8.3728D-01,7.6947D-01,7.0373D-01,6.3561D-01,5.7438D-01,5.1959D-01,
     &4.6984D-01,4.3187D-01,3.9529D-01,3.5864D-01,3.2783D-01,3.0967D-01,
     &2.9444D-01,2.8428D-01,2.7469D-01,2.6638D-01,2.5813D-01,2.4942D-01,
     &2.3986D-01,2.2937D-01,2.1819D-01,2.0682D-01,1.9665D-01,1.8966D-01/
      DATA (XUDF_L(K),K= 1483, 1596) /
     &1.8971D-01,1.9926D-01,1.9036D-01,1.8442D-01,1.8192D-01,1.8362D-01,
     &1.9037D-01,2.0318D-01,2.2459D-01,2.5904D-01,3.1665D-01,4.2407D-01,
     &0.0000D+00,0.0000D+00,1.2798D+01,1.1861D+01,1.0986D+01,1.0144D+01,
     &9.3643D+00,8.5887D+00,7.8706D+00,7.1866D+00,6.5568D+00,5.9419D+00,
     &5.3754D+00,4.8419D+00,4.3593D+00,4.0864D+00,3.8109D+00,3.5127D+00,
     &3.2315D+00,2.9714D+00,2.7252D+00,2.5309D+00,2.3356D+00,2.1269D+00,
     &1.9338D+00,1.7578D+00,1.5939D+00,1.4656D+00,1.3394D+00,1.2075D+00,
     &1.0875D+00,9.8023D-01,8.8256D-01,8.0772D-01,7.3533D-01,6.6054D-01,
     &5.9364D-01,5.3423D-01,4.8009D-01,4.3930D-01,4.0003D-01,3.6079D-01,
     &3.2768D-01,3.0809D-01,2.9130D-01,2.7993D-01,2.6898D-01,2.5976D-01,
     &2.5062D-01,2.4123D-01,2.3116D-01,2.2040D-01,2.0917D-01,1.9814D-01,
     &1.8865D-01,1.8272D-01,1.8428D-01,1.8820D-01,1.7883D-01,1.7238D-01,
     &1.6914D-01,1.6979D-01,1.7482D-01,1.8534D-01,2.0325D-01,2.3214D-01,
     &2.8022D-01,3.6659D-01,0.0000D+00,0.0000D+00,1.4900D+01,1.3767D+01,
     &1.2708D+01,1.1700D+01,1.0766D+01,9.8403D+00,8.9832D+00,8.1757D+00,
     &7.4366D+00,6.7121D+00,6.0486D+00,5.4300D+00,4.8704D+00,4.5555D+00,
     &4.2371D+00,3.8955D+00,3.5734D+00,3.2760D+00,2.9952D+00,2.7738D+00,
     &2.5528D+00,2.3175D+00,2.1001D+00,1.9012D+00,1.7176D+00,1.5750D+00,
     &1.4344D+00,1.2880D+00,1.1547D+00,1.0364D+00,9.2859D-01,8.4652D-01/
      DATA (XUDF_L(K),K= 1597, 1710) /
     &7.6723D-01,6.8578D-01,6.1255D-01,5.4848D-01,4.9034D-01,4.4649D-01,
     &4.0456D-01,3.6275D-01,3.2738D-01,3.0624D-01,2.8805D-01,2.7544D-01,
     &2.6343D-01,2.5315D-01,2.4318D-01,2.3314D-01,2.2263D-01,2.1166D-01,
     &2.0051D-01,1.8983D-01,1.8102D-01,1.7610D-01,1.7901D-01,1.7764D-01,
     &1.6791D-01,1.6102D-01,1.5715D-01,1.5684D-01,1.6056D-01,1.6899D-01,
     &1.8376D-01,2.0786D-01,2.4776D-01,3.1470D-01,0.0000D+00,0.0000D+00,
     &1.7212D+01,1.5853D+01,1.4590D+01,1.3390D+01,1.2283D+01,1.1191D+01,
     &1.0185D+01,9.2395D+00,8.3762D+00,7.5315D+00,6.7670D+00,6.0503D+00,
     &5.4086D+00,5.0481D+00,4.6843D+00,4.2940D+00,3.9280D+00,3.5917D+00,
     &3.2752D+00,3.0252D+00,2.7768D+00,2.5132D+00,2.2690D+00,2.0490D+00,
     &1.8445D+00,1.6857D+00,1.5301D+00,1.3685D+00,1.2219D+00,1.0920D+00,
     &9.7438D-01,8.8478D-01,7.9825D-01,7.1007D-01,6.3111D-01,5.6196D-01,
     &5.0016D-01,4.5321D-01,4.0867D-01,3.6435D-01,3.2686D-01,3.0431D-01,
     &2.8470D-01,2.7109D-01,2.5789D-01,2.4674D-01,2.3605D-01,2.2547D-01,
     &2.1459D-01,2.0348D-01,1.9237D-01,1.8201D-01,1.7376D-01,1.6982D-01,
     &1.7398D-01,1.6789D-01,1.5795D-01,1.5065D-01,1.4630D-01,1.4521D-01,
     &1.4773D-01,1.5443D-01,1.6659D-01,1.8664D-01,2.1966D-01,2.6878D-01,
     &0.0000D+00,0.0000D+00,1.9526D+01,1.7951D+01,1.6470D+01,1.5074D+01,
     &1.3790D+01,1.2527D+01,1.1370D+01,1.0282D+01,9.2958D+00,8.3330D+00/
      DATA (XUDF_L(K),K= 1711, 1824) /
     &7.4603D+00,6.6536D+00,5.9285D+00,5.5219D+00,5.1141D+00,4.6768D+00,
     &4.2681D+00,3.8926D+00,3.5402D+00,3.2626D+00,2.9882D+00,2.6963D+00,
     &2.4284D+00,2.1851D+00,1.9619D+00,1.7885D+00,1.6187D+00,1.4429D+00,
     &1.2838D+00,1.1431D+00,1.0159D+00,9.1924D-01,8.2663D-01,7.3180D-01,
     &6.4793D-01,5.7429D-01,5.0828D-01,4.5904D-01,4.1215D-01,3.6558D-01,
     &3.2620D-01,3.0238D-01,2.8167D-01,2.6700D-01,2.5302D-01,2.4098D-01,
     &2.2975D-01,2.1873D-01,2.0756D-01,1.9633D-01,1.8532D-01,1.7533D-01,
     &1.6763D-01,1.6450D-01,1.6959D-01,1.5953D-01,1.4943D-01,1.4185D-01,
     &1.3716D-01,1.3545D-01,1.3705D-01,1.4238D-01,1.5258D-01,1.6945D-01,
     &1.9705D-01,2.3049D-01,0.0000D+00,0.0000D+00,2.2141D+01,2.0286D+01,
     &1.8570D+01,1.6948D+01,1.5466D+01,1.4010D+01,1.2679D+01,1.1431D+01,
     &1.0303D+01,9.2106D+00,8.2239D+00,7.3077D+00,6.4926D+00,6.0348D+00,
     &5.5765D+00,5.0879D+00,4.6321D+00,4.2138D+00,3.8233D+00,3.5162D+00,
     &3.2122D+00,2.8907D+00,2.5960D+00,2.3300D+00,2.0856D+00,1.8954D+00,
     &1.7110D+00,1.5199D+00,1.3476D+00,1.1955D+00,1.0584D+00,9.5478D-01,
     &8.5531D-01,7.5417D-01,6.6439D-01,5.8623D-01,5.1682D-01,4.6468D-01,
     &4.1541D-01,3.6662D-01,3.2538D-01,3.0035D-01,2.7843D-01,2.6291D-01,
     &2.4798D-01,2.3522D-01,2.2346D-01,2.1203D-01,2.0062D-01,1.8935D-01,
     &1.7843D-01,1.6874D-01,1.6163D-01,1.5920D-01,1.6520D-01,1.5147D-01/
      DATA (XUDF_L(K),K= 1825, 1836) /
     &1.4120D-01,1.3349D-01,1.2844D-01,1.2620D-01,1.2701D-01,1.3118D-01,
     &1.3954D-01,1.5369D-01,1.7631D-01,1.9416D-01,0.0000D+00,0.0000D+00/
      DATA (XSF_L(K),K=    1,  114) /
     &8.9277D-03,9.2838D-03,9.6380D-03,9.9960D-03,1.0349D-02,1.0719D-02,
     &1.1082D-02,1.1442D-02,1.1792D-02,1.2148D-02,1.2489D-02,1.2817D-02,
     &1.3124D-02,1.3295D-02,1.3474D-02,1.3661D-02,1.3835D-02,1.3985D-02,
     &1.4121D-02,1.4217D-02,1.4303D-02,1.4379D-02,1.4419D-02,1.4434D-02,
     &1.4412D-02,1.4366D-02,1.4286D-02,1.4158D-02,1.3991D-02,1.3790D-02,
     &1.3553D-02,1.3335D-02,1.3094D-02,1.2821D-02,1.2580D-02,1.2410D-02,
     &1.2357D-02,1.2459D-02,1.2790D-02,1.3571D-02,1.5018D-02,1.6665D-02,
     &1.9113D-02,2.1832D-02,2.5587D-02,2.9818D-02,3.4535D-02,3.9813D-02,
     &4.5737D-02,5.2358D-02,5.9765D-02,6.8021D-02,7.7185D-02,8.7258D-02,
     &9.8198D-02,1.1073D-01,1.4216D-01,1.8364D-01,2.3959D-01,3.1758D-01,
     &4.3050D-01,6.0203D-01,8.8214D-01,1.3845D+00,2.4294D+00,5.2463D+00,
     &1.8903D+01,0.0000D+00,1.4987D-02,1.5468D-02,1.5936D-02,1.6403D-02,
     &1.6855D-02,1.7319D-02,1.7760D-02,1.8194D-02,1.8600D-02,1.9008D-02,
     &1.9382D-02,1.9730D-02,2.0033D-02,2.0199D-02,2.0359D-02,2.0523D-02,
     &2.0654D-02,2.0760D-02,2.0831D-02,2.0870D-02,2.0886D-02,2.0858D-02,
     &2.0798D-02,2.0680D-02,2.0523D-02,2.0363D-02,2.0127D-02,1.9825D-02,
     &1.9464D-02,1.9060D-02,1.8607D-02,1.8200D-02,1.7750D-02,1.7240D-02,
     &1.6759D-02,1.6362D-02,1.6103D-02,1.6050D-02,1.6240D-02,1.6916D-02,
     &1.8336D-02,2.0030D-02,2.2586D-02,2.5447D-02,2.9418D-02,3.3874D-02/
      DATA (XSF_L(K),K=  115,  228) /
     &3.8821D-02,4.4375D-02,5.0509D-02,5.7343D-02,6.4974D-02,7.3385D-02,
     &8.2640D-02,9.2732D-02,1.0354D-01,1.1667D-01,1.4809D-01,1.8910D-01,
     &2.4387D-01,3.1940D-01,4.2764D-01,5.9054D-01,8.5228D-01,1.3150D+00,
     &2.2623D+00,4.7596D+00,1.6445D+01,0.0000D+00,2.5010D-02,2.5616D-02,
     &2.6180D-02,2.6758D-02,2.7279D-02,2.7792D-02,2.8274D-02,2.8729D-02,
     &2.9134D-02,2.9513D-02,2.9836D-02,3.0110D-02,3.0324D-02,3.0417D-02,
     &3.0492D-02,3.0537D-02,3.0551D-02,3.0517D-02,3.0432D-02,3.0326D-02,
     &3.0181D-02,2.9954D-02,2.9663D-02,2.9316D-02,2.8913D-02,2.8508D-02,
     &2.8021D-02,2.7422D-02,2.6741D-02,2.5997D-02,2.5204D-02,2.4500D-02,
     &2.3734D-02,2.2858D-02,2.2019D-02,2.1281D-02,2.0698D-02,2.0402D-02,
     &2.0365D-02,2.0844D-02,2.2137D-02,2.3807D-02,2.6404D-02,2.9338D-02,
     &3.3433D-02,3.8036D-02,4.3135D-02,4.8799D-02,5.5061D-02,6.1999D-02,
     &6.9633D-02,7.8024D-02,8.7156D-02,9.6998D-02,1.0742D-01,1.2099D-01,
     &1.5162D-01,1.9121D-01,2.4363D-01,3.1510D-01,4.1638D-01,5.6669D-01,
     &8.0557D-01,1.2216D+00,2.0572D+00,4.2084D+00,1.3911D+01,0.0000D+00,
     &4.2554D-02,4.3210D-02,4.3820D-02,4.4379D-02,4.4862D-02,4.5317D-02,
     &4.5708D-02,4.6037D-02,4.6300D-02,4.6434D-02,4.6540D-02,4.6530D-02,
     &4.6426D-02,4.6317D-02,4.6155D-02,4.5919D-02,4.5622D-02,4.5267D-02,
     &4.4833D-02,4.4425D-02,4.3932D-02,4.3298D-02,4.2582D-02,4.1785D-02/
      DATA (XSF_L(K),K=  229,  342) /
     &4.0903D-02,4.0097D-02,3.9179D-02,3.8047D-02,3.6815D-02,3.5547D-02,
     &3.4199D-02,3.3020D-02,3.1748D-02,3.0298D-02,2.8905D-02,2.7644D-02,
     &2.6563D-02,2.5882D-02,2.5485D-02,2.5614D-02,2.6651D-02,2.8199D-02,
     &3.0731D-02,3.3652D-02,3.7768D-02,4.2390D-02,4.7530D-02,5.3188D-02,
     &5.9436D-02,6.6257D-02,7.3734D-02,8.1918D-02,9.0696D-02,1.0004D-01,
     &1.0978D-01,1.2357D-01,1.5274D-01,1.8999D-01,2.3888D-01,3.0452D-01,
     &3.9656D-01,5.3136D-01,7.4246D-01,1.1043D+00,1.8158D+00,3.6023D+00,
     &0.0000D+00,0.0000D+00,7.3602D-02,7.4085D-02,7.4460D-02,7.4729D-02,
     &7.4904D-02,7.4982D-02,7.4902D-02,7.4713D-02,7.4446D-02,7.3972D-02,
     &7.3397D-02,7.2626D-02,7.1803D-02,7.1200D-02,7.0479D-02,6.9610D-02,
     &6.8654D-02,6.7624D-02,6.6495D-02,6.5467D-02,6.4313D-02,6.2898D-02,
     &6.1380D-02,5.9788D-02,5.8079D-02,5.6557D-02,5.4876D-02,5.2866D-02,
     &5.0733D-02,4.8592D-02,4.6341D-02,4.4415D-02,4.2370D-02,4.0073D-02,
     &3.7825D-02,3.5778D-02,3.3956D-02,3.2702D-02,3.1749D-02,3.1334D-02,
     &3.1922D-02,3.3216D-02,3.5534D-02,3.8322D-02,4.2321D-02,4.6830D-02,
     &5.1816D-02,5.7335D-02,6.3369D-02,6.9947D-02,7.7109D-02,8.4752D-02,
     &9.2948D-02,1.0153D-01,1.1031D-01,1.2405D-01,1.5100D-01,1.8509D-01,
     &2.2905D-01,2.8761D-01,3.6847D-01,4.8537D-01,6.6543D-01,9.6831D-01,
     &1.5524D+00,2.9766D+00,0.0000D+00,0.0000D+00,1.1509D-01,1.1500D-01/
      DATA (XSF_L(K),K=  343,  456) /
     &1.1474D-01,1.1430D-01,1.1371D-01,1.1292D-01,1.1196D-01,1.1079D-01,
     &1.0948D-01,1.0791D-01,1.0620D-01,1.0426D-01,1.0215D-01,1.0076D-01,
     &9.9224D-02,9.7466D-02,9.5472D-02,9.3507D-02,9.1346D-02,8.9460D-02,
     &8.7382D-02,8.4914D-02,8.2326D-02,7.9663D-02,7.6874D-02,7.4459D-02,
     &7.1794D-02,6.8694D-02,6.5489D-02,6.2266D-02,5.8964D-02,5.6164D-02,
     &5.3226D-02,4.9916D-02,4.6721D-02,4.3794D-02,4.1128D-02,3.9225D-02,
     &3.7654D-02,3.6613D-02,3.6666D-02,3.7626D-02,3.9655D-02,4.2227D-02,
     &4.6000D-02,5.0288D-02,5.5044D-02,6.0308D-02,6.6020D-02,7.2218D-02,
     &7.8943D-02,8.6079D-02,9.3611D-02,1.0141D-01,1.0925D-01,1.2274D-01,
     &1.4748D-01,1.7840D-01,2.1791D-01,2.6997D-01,3.4109D-01,4.4280D-01,
     &5.9706D-01,8.5325D-01,1.3371D+00,2.4909D+00,0.0000D+00,0.0000D+00,
     &1.8131D-01,1.7986D-01,1.7802D-01,1.7597D-01,1.7372D-01,1.7110D-01,
     &1.6825D-01,1.6515D-01,1.6187D-01,1.5820D-01,1.5428D-01,1.5016D-01,
     &1.4582D-01,1.4314D-01,1.4017D-01,1.3677D-01,1.3315D-01,1.2951D-01,
     &1.2571D-01,1.2248D-01,1.1891D-01,1.1472D-01,1.1045D-01,1.0615D-01,
     &1.0173D-01,9.7944D-02,9.3854D-02,8.9131D-02,8.4347D-02,7.9597D-02,
     &7.4799D-02,7.0788D-02,6.6599D-02,6.1932D-02,5.7438D-02,5.3307D-02,
     &4.9546D-02,4.6816D-02,4.4417D-02,4.2536D-02,4.1862D-02,4.2361D-02,
     &4.3960D-02,4.6198D-02,4.9612D-02,5.3553D-02,5.7974D-02,6.2830D-02/
      DATA (XSF_L(K),K=  457,  570) /
     &6.8141D-02,7.3865D-02,7.9970D-02,8.6422D-02,9.3160D-02,1.0006D-01,
     &1.0685D-01,1.1989D-01,1.4199D-01,1.6937D-01,2.0407D-01,2.4925D-01,
     &3.1029D-01,3.9635D-01,5.2529D-01,7.3579D-01,1.1263D+00,2.0347D+00,
     &0.0000D+00,0.0000D+00,2.6278D-01,2.5883D-01,2.5460D-01,2.5007D-01,
     &2.4526D-01,2.3995D-01,2.3437D-01,2.2848D-01,2.2242D-01,2.1578D-01,
     &2.0894D-01,2.0181D-01,1.9465D-01,1.9018D-01,1.8540D-01,1.7984D-01,
     &1.7415D-01,1.6846D-01,1.6261D-01,1.5768D-01,1.5234D-01,1.4615D-01,
     &1.3987D-01,1.3368D-01,1.2736D-01,1.2199D-01,1.1628D-01,1.0975D-01,
     &1.0321D-01,9.6788D-02,9.0380D-02,8.5059D-02,7.9532D-02,7.3436D-02,
     &6.7594D-02,6.2243D-02,5.7363D-02,5.3720D-02,5.0502D-02,4.7772D-02,
     &4.6346D-02,4.6358D-02,4.7497D-02,4.9377D-02,5.2401D-02,5.5965D-02,
     &6.0009D-02,6.4489D-02,6.9334D-02,7.4546D-02,8.0117D-02,8.5936D-02,
     &9.1972D-02,9.8056D-02,1.0398D-01,1.1644D-01,1.3628D-01,1.6068D-01,
     &1.9127D-01,2.3085D-01,2.8377D-01,3.5756D-01,4.6698D-01,6.4315D-01,
     &9.6485D-01,1.6969D+00,0.0000D+00,0.0000D+00,3.6944D-01,3.6187D-01,
     &3.5380D-01,3.4525D-01,3.3659D-01,3.2716D-01,3.1761D-01,3.0767D-01,
     &2.9759D-01,2.8675D-01,2.7586D-01,2.6462D-01,2.5339D-01,2.4660D-01,
     &2.3933D-01,2.3101D-01,2.2257D-01,2.1415D-01,2.0571D-01,1.9854D-01,
     &1.9083D-01,1.8216D-01,1.7338D-01,1.6480D-01,1.5613D-01,1.4885D-01/
      DATA (XSF_L(K),K=  571,  684) /
     &1.4115D-01,1.3244D-01,1.2380D-01,1.1542D-01,1.0713D-01,1.0031D-01,
     &9.3226D-02,8.5515D-02,7.8171D-02,7.1449D-02,6.5307D-02,6.0723D-02,
     &5.6523D-02,5.2878D-02,5.0622D-02,5.0109D-02,5.0720D-02,5.2187D-02,
     &5.4770D-02,5.7950D-02,6.1582D-02,6.5595D-02,6.9997D-02,7.4716D-02,
     &7.9677D-02,8.4886D-02,9.0221D-02,9.5543D-02,1.0065D-01,1.1245D-01,
     &1.3012D-01,1.5166D-01,1.7859D-01,2.1305D-01,2.5881D-01,3.2188D-01,
     &4.1454D-01,5.6186D-01,8.2718D-01,1.4188D+00,0.0000D+00,0.0000D+00,
     &4.9195D-01,4.7916D-01,4.6620D-01,4.5277D-01,4.3908D-01,4.2463D-01,
     &4.0985D-01,3.9491D-01,3.7975D-01,3.6377D-01,3.4790D-01,3.3178D-01,
     &3.1592D-01,3.0640D-01,2.9622D-01,2.8462D-01,2.7303D-01,2.6160D-01,
     &2.5012D-01,2.4047D-01,2.3023D-01,2.1867D-01,2.0717D-01,1.9597D-01,
     &1.8477D-01,1.7546D-01,1.6568D-01,1.5468D-01,1.4387D-01,1.3343D-01,
     &1.2319D-01,1.1482D-01,1.0622D-01,9.6828D-02,8.7978D-02,7.9884D-02,
     &7.2526D-02,6.6973D-02,6.1948D-02,5.7359D-02,5.4304D-02,5.3263D-02,
     &5.3381D-02,5.4456D-02,5.6601D-02,5.9380D-02,6.2613D-02,6.6252D-02,
     &7.0174D-02,7.4432D-02,7.8943D-02,8.3559D-02,8.8282D-02,9.2963D-02,
     &9.7382D-02,1.0858D-01,1.2441D-01,1.4363D-01,1.6745D-01,1.9778D-01,
     &2.3771D-01,2.9246D-01,3.7200D-01,4.9738D-01,7.2010D-01,1.2083D+00,
     &0.0000D+00,0.0000D+00,6.4521D-01,6.2534D-01,6.0540D-01,5.8499D-01/
      DATA (XSF_L(K),K=  685,  798) /
     &5.6467D-01,5.4301D-01,5.2143D-01,4.9951D-01,4.7813D-01,4.5538D-01,
     &4.3325D-01,4.1083D-01,3.8899D-01,3.7591D-01,3.6210D-01,3.4648D-01,
     &3.3091D-01,3.1578D-01,3.0062D-01,2.8797D-01,2.7469D-01,2.5979D-01,
     &2.4501D-01,2.3066D-01,2.1649D-01,2.0481D-01,1.9252D-01,1.7884D-01,
     &1.6549D-01,1.5274D-01,1.4029D-01,1.3018D-01,1.1985D-01,1.0865D-01,
     &9.8135D-02,8.8550D-02,7.9829D-02,7.3318D-02,6.7269D-02,6.1748D-02,
     &5.7838D-02,5.6250D-02,5.5826D-02,5.6474D-02,5.8181D-02,6.0533D-02,
     &6.3373D-02,6.6563D-02,7.0085D-02,7.3865D-02,7.7842D-02,8.1937D-02,
     &8.6092D-02,9.0169D-02,9.3962D-02,1.0448D-01,1.1858D-01,1.3561D-01,
     &1.5663D-01,1.8318D-01,2.1803D-01,2.6529D-01,3.3349D-01,4.3985D-01,
     &6.2661D-01,1.0291D+00,0.0000D+00,0.0000D+00,8.2462D-01,7.9558D-01,
     &7.6680D-01,7.3764D-01,7.0860D-01,6.7834D-01,6.4822D-01,6.1798D-01,
     &5.8880D-01,5.5792D-01,5.2800D-01,4.9801D-01,4.6912D-01,4.5197D-01,
     &4.3393D-01,4.1360D-01,3.9348D-01,3.7394D-01,3.5462D-01,3.3856D-01,
     &3.2180D-01,3.0303D-01,2.8460D-01,2.6681D-01,2.4932D-01,2.3502D-01,
     &2.2005D-01,2.0359D-01,1.8747D-01,1.7224D-01,1.5746D-01,1.4551D-01,
     &1.3337D-01,1.2028D-01,1.0805D-01,9.6986D-02,8.6877D-02,7.9334D-02,
     &7.2326D-02,6.5799D-02,6.1060D-02,5.8911D-02,5.7957D-02,5.8189D-02,
     &5.9441D-02,6.1387D-02,6.3834D-02,6.6632D-02,6.9732D-02,7.3070D-02/
      DATA (XSF_L(K),K=  799,  912) /
     &7.6595D-02,8.0190D-02,8.3816D-02,8.7358D-02,9.0631D-02,1.0046D-01,
     &1.1304D-01,1.2815D-01,1.4670D-01,1.7006D-01,2.0049D-01,2.4154D-01,
     &3.0039D-01,3.9121D-01,5.4894D-01,8.8378D-01,0.0000D+00,0.0000D+00,
     &1.0199D+00,9.8025D-01,9.4100D-01,9.0151D-01,8.6283D-01,8.2243D-01,
     &7.8262D-01,7.4321D-01,7.0465D-01,6.6494D-01,6.2647D-01,5.8811D-01,
     &5.5152D-01,5.2985D-01,5.0721D-01,4.8183D-01,4.5681D-01,4.3274D-01,
     &4.0883D-01,3.8916D-01,3.6878D-01,3.4589D-01,3.2366D-01,3.0238D-01,
     &2.8152D-01,2.6437D-01,2.4685D-01,2.2733D-01,2.0858D-01,1.9085D-01,
     &1.7375D-01,1.6000D-01,1.4607D-01,1.3115D-01,1.1722D-01,1.0469D-01,
     &9.3284D-02,8.4739D-02,7.6803D-02,6.9420D-02,6.3844D-02,6.1178D-02,
     &5.9720D-02,5.9561D-02,6.0398D-02,6.1984D-02,6.4051D-02,6.6494D-02,
     &6.9202D-02,7.2161D-02,7.5274D-02,7.8453D-02,8.1651D-02,8.4728D-02,
     &8.7564D-02,9.6777D-02,1.0806D-01,1.2157D-01,1.3806D-01,1.5882D-01,
     &1.8566D-01,2.2170D-01,2.7301D-01,3.5168D-01,4.8696D-01,7.7010D-01,
     &0.0000D+00,0.0000D+00,1.3158D+00,1.2585D+00,1.2024D+00,1.1462D+00,
     &1.0919D+00,1.0352D+00,9.8042D-01,9.2608D-01,8.7345D-01,8.1987D-01,
     &7.6814D-01,7.1724D-01,6.6882D-01,6.4053D-01,6.1093D-01,5.7796D-01,
     &5.4572D-01,5.1470D-01,4.8433D-01,4.5934D-01,4.3358D-01,4.0495D-01,
     &3.7717D-01,3.5082D-01,3.2513D-01,3.0408D-01,2.8258D-01,2.5918D-01/
      DATA (XSF_L(K),K=  913, 1026) /
     &2.3648D-01,2.1538D-01,1.9510D-01,1.7888D-01,1.6255D-01,1.4508D-01,
     &1.2895D-01,1.1443D-01,1.0131D-01,9.1507D-02,8.2387D-02,7.3778D-02,
     &6.7147D-02,6.3813D-02,6.1721D-02,6.1065D-02,6.1373D-02,6.2475D-02,
     &6.4105D-02,6.6079D-02,6.8362D-02,7.0856D-02,7.3440D-02,7.6143D-02,
     &7.8812D-02,8.1388D-02,8.3726D-02,9.2167D-02,1.0190D-01,1.1355D-01,
     &1.2780D-01,1.4554D-01,1.6841D-01,1.9900D-01,2.4223D-01,3.0775D-01,
     &4.1920D-01,6.4849D-01,0.0000D+00,0.0000D+00,1.6483D+00,1.5703D+00,
     &1.4940D+00,1.4180D+00,1.3449D+00,1.2694D+00,1.1966D+00,1.1250D+00,
     &1.0566D+00,9.8644D-01,9.1985D-01,8.5482D-01,7.9312D-01,7.5722D-01,
     &7.1986D-01,6.7849D-01,6.3821D-01,5.9972D-01,5.6214D-01,5.3143D-01,
     &4.9987D-01,4.6500D-01,4.3136D-01,3.9956D-01,3.6875D-01,3.4379D-01,
     &3.1832D-01,2.9044D-01,2.6397D-01,2.3923D-01,2.1580D-01,1.9706D-01,
     &1.7829D-01,1.5838D-01,1.3999D-01,1.2356D-01,1.0875D-01,9.7664D-02,
     &8.7392D-02,7.7645D-02,7.0035D-02,6.6062D-02,6.3365D-02,6.2239D-02,
     &6.2062D-02,6.2731D-02,6.3942D-02,6.5526D-02,6.7390D-02,6.9436D-02,
     &7.1635D-02,7.3891D-02,7.6122D-02,7.8246D-02,8.0196D-02,8.7884D-02,
     &9.6357D-02,1.0648D-01,1.1880D-01,1.3413D-01,1.5386D-01,1.7993D-01,
     &2.1655D-01,2.7189D-01,3.6486D-01,5.5332D-01,0.0000D+00,0.0000D+00,
     &2.0271D+00,1.9234D+00,1.8224D+00,1.7226D+00,1.6272D+00,1.5293D+00/
      DATA (XSF_L(K),K= 1027, 1140) /
     &1.4356D+00,1.3438D+00,1.2568D+00,1.1682D+00,1.0841D+00,1.0026D+00,
     &9.2625D-01,8.8207D-01,8.3568D-01,7.8523D-01,7.3607D-01,6.8926D-01,
     &6.4385D-01,6.0685D-01,5.6892D-01,5.2730D-01,4.8731D-01,4.4961D-01,
     &4.1331D-01,3.8417D-01,3.5441D-01,3.2210D-01,2.9168D-01,2.6323D-01,
     &2.3631D-01,2.1500D-01,1.9374D-01,1.7129D-01,1.5067D-01,1.3231D-01,
     &1.1579D-01,1.0349D-01,9.2080D-02,8.1205D-02,7.2626D-02,6.8039D-02,
     &6.4761D-02,6.3188D-02,6.2549D-02,6.2795D-02,6.3617D-02,6.4835D-02,
     &6.6329D-02,6.8017D-02,6.9809D-02,7.1667D-02,7.3520D-02,7.5270D-02,
     &7.6864D-02,8.3899D-02,9.1206D-02,1.0002D-01,1.1070D-01,1.2399D-01,
     &1.4094D-01,1.6341D-01,1.9474D-01,2.4163D-01,3.1971D-01,4.7587D-01,
     &0.0000D+00,0.0000D+00,2.4392D+00,2.3049D+00,2.1760D+00,2.0502D+00,
     &1.9296D+00,1.8065D+00,1.6895D+00,1.5750D+00,1.4674D+00,1.3585D+00,
     &1.2554D+00,1.1565D+00,1.0638D+00,1.0103D+00,9.5527D-01,8.9449D-01,
     &8.3572D-01,7.8018D-01,7.2635D-01,6.8280D-01,6.3819D-01,5.8948D-01,
     &5.4299D-01,4.9923D-01,4.5740D-01,4.2371D-01,3.8978D-01,3.5296D-01,
     &3.1832D-01,2.8629D-01,2.5599D-01,2.3212D-01,2.0840D-01,1.8346D-01,
     &1.6065D-01,1.4043D-01,1.2229D-01,1.0880D-01,9.6294D-02,8.4335D-02,
     &7.4905D-02,6.9717D-02,6.5897D-02,6.3914D-02,6.2851D-02,6.2731D-02,
     &6.3183D-02,6.4075D-02,6.5225D-02,6.6597D-02,6.8048D-02,6.9577D-02/
      DATA (XSF_L(K),K= 1141, 1254) /
     &7.1093D-02,7.2525D-02,7.3842D-02,8.0241D-02,8.6615D-02,9.4292D-02,
     &1.0360D-01,1.1517D-01,1.2992D-01,1.4936D-01,1.7633D-01,2.1652D-01,
     &2.8294D-01,4.1389D-01,0.0000D+00,0.0000D+00,2.9162D+00,2.7470D+00,
     &2.5840D+00,2.4244D+00,2.2743D+00,2.1215D+00,1.9764D+00,1.8358D+00,
     &1.7035D+00,1.5708D+00,1.4463D+00,1.3268D+00,1.2152D+00,1.1514D+00,
     &1.0857D+00,1.0132D+00,9.4449D-01,8.7867D-01,8.1556D-01,7.6453D-01,
     &7.1252D-01,6.5602D-01,6.0218D-01,5.5192D-01,5.0387D-01,4.6545D-01,
     &4.2679D-01,3.8521D-01,3.4602D-01,3.1005D-01,2.7623D-01,2.4962D-01,
     &2.2332D-01,1.9577D-01,1.7070D-01,1.4856D-01,1.2874D-01,1.1402D-01,
     &1.0040D-01,8.7343D-02,7.6984D-02,7.1254D-02,6.6892D-02,6.4508D-02,
     &6.3019D-02,6.2518D-02,6.2667D-02,6.3211D-02,6.4031D-02,6.5064D-02,
     &6.6243D-02,6.7458D-02,6.8679D-02,6.9830D-02,7.0885D-02,7.6672D-02,
     &8.2192D-02,8.8844D-02,9.6930D-02,1.0696D-01,1.1972D-01,1.3654D-01,
     &1.5978D-01,1.9411D-01,2.5048D-01,3.6023D-01,0.0000D+00,0.0000D+00,
     &3.4281D+00,3.2194D+00,3.0180D+00,2.8239D+00,2.6400D+00,2.4537D+00,
     &2.2781D+00,2.1087D+00,1.9503D+00,1.7915D+00,1.6433D+00,1.5021D+00,
     &1.3711D+00,1.2958D+00,1.2191D+00,1.1350D+00,1.0536D+00,9.7846D-01,
     &9.0526D-01,8.4668D-01,7.8697D-01,7.2243D-01,6.6110D-01,6.0402D-01,
     &5.4971D-01,5.0652D-01,4.6307D-01,4.1647D-01,3.7287D-01,3.3288D-01/
      DATA (XSF_L(K),K= 1255, 1368) /
     &2.9545D-01,2.6636D-01,2.3751D-01,2.0740D-01,1.8012D-01,1.5611D-01,
     &1.3467D-01,1.1881D-01,1.0414D-01,9.0105D-02,7.8839D-02,7.2563D-02,
     &6.7703D-02,6.4930D-02,6.3070D-02,6.2241D-02,6.2071D-02,6.2347D-02,
     &6.2882D-02,6.3645D-02,6.4526D-02,6.5473D-02,6.6427D-02,6.7333D-02,
     &6.8194D-02,7.3430D-02,7.8217D-02,8.3974D-02,9.1017D-02,9.9745D-02,
     &1.1088D-01,1.2552D-01,1.4563D-01,1.7528D-01,2.2351D-01,3.1636D-01,
     &0.0000D+00,0.0000D+00,3.9892D+00,3.7328D+00,3.4900D+00,3.2549D+00,
     &3.0344D+00,2.8108D+00,2.6014D+00,2.4001D+00,2.2123D+00,2.0253D+00,
     &1.8518D+00,1.6860D+00,1.5339D+00,1.4463D+00,1.3575D+00,1.2608D+00,
     &1.1678D+00,1.0809D+00,9.9767D-01,9.3087D-01,8.6314D-01,7.8996D-01,
     &7.2083D-01,6.5671D-01,5.9602D-01,5.4775D-01,4.9935D-01,4.4773D-01,
     &3.9951D-01,3.5571D-01,3.1467D-01,2.8272D-01,2.5135D-01,2.1871D-01,
     &1.8923D-01,1.6331D-01,1.4031D-01,1.2332D-01,1.0762D-01,9.2560D-02,
     &8.0473D-02,7.3714D-02,6.8385D-02,6.5246D-02,6.3019D-02,6.1878D-02,
     &6.1420D-02,6.1413D-02,6.1734D-02,6.2226D-02,6.2861D-02,6.3564D-02,
     &6.4288D-02,6.4985D-02,6.5657D-02,7.0367D-02,7.4522D-02,7.9506D-02,
     &8.5651D-02,9.3297D-02,1.0298D-01,1.1572D-01,1.3323D-01,1.5884D-01,
     &2.0039D-01,2.7925D-01,0.0000D+00,0.0000D+00,4.5788D+00,4.2729D+00,
     &3.9840D+00,3.7039D+00,3.4438D+00,3.1812D+00,2.9349D+00,2.6996D+00/
      DATA (XSF_L(K),K= 1369, 1482) /
     &2.4810D+00,2.2644D+00,2.0633D+00,1.8732D+00,1.6979D+00,1.5988D+00,
     &1.4974D+00,1.3865D+00,1.2812D+00,1.1834D+00,1.0891D+00,1.0143D+00,
     &9.3839D-01,8.5662D-01,7.7948D-01,7.0838D-01,6.4106D-01,5.8780D-01,
     &5.3454D-01,4.7781D-01,4.2528D-01,3.7737D-01,3.3289D-01,2.9818D-01,
     &2.6446D-01,2.2932D-01,1.9770D-01,1.7005D-01,1.4552D-01,1.2746D-01,
     &1.1078D-01,9.4770D-02,8.1957D-02,7.4689D-02,6.8915D-02,6.5457D-02,
     &6.2902D-02,6.1493D-02,6.0768D-02,6.0515D-02,6.0585D-02,6.0863D-02,
     &6.1298D-02,6.1789D-02,6.2311D-02,6.2835D-02,6.3340D-02,6.7601D-02,
     &7.1162D-02,7.5516D-02,8.0878D-02,8.7566D-02,9.6095D-02,1.0725D-01,
     &1.2258D-01,1.4495D-01,1.8090D-01,2.4841D-01,0.0000D+00,0.0000D+00,
     &5.4774D+00,5.0929D+00,4.7320D+00,4.3841D+00,4.0592D+00,3.7350D+00,
     &3.4329D+00,3.1454D+00,2.8799D+00,2.6172D+00,2.3747D+00,2.1466D+00,
     &1.9383D+00,1.8195D+00,1.6996D+00,1.5689D+00,1.4457D+00,1.3301D+00,
     &1.2211D+00,1.1339D+00,1.0456D+00,9.5119D-01,8.6259D-01,7.8097D-01,
     &7.0419D-01,6.4380D-01,5.8358D-01,5.1955D-01,4.6051D-01,4.0719D-01,
     &3.5768D-01,3.1962D-01,2.8220D-01,2.4360D-01,2.0909D-01,1.7895D-01,
     &1.5240D-01,1.3282D-01,1.1484D-01,9.7655D-02,8.3739D-02,7.5857D-02,
     &6.9509D-02,6.5616D-02,6.2633D-02,6.0853D-02,5.9819D-02,5.9271D-02,
     &5.9038D-02,5.9046D-02,5.9192D-02,5.9432D-02,5.9709D-02,6.0008D-02/
      DATA (XSF_L(K),K= 1483, 1596) /
     &6.0340D-02,6.4032D-02,6.6851D-02,7.0446D-02,7.4870D-02,8.0457D-02,
     &8.7554D-02,9.6862D-02,1.0964D-01,1.2821D-01,1.5779D-01,2.1189D-01,
     &0.0000D+00,0.0000D+00,6.3982D+00,5.9307D+00,5.4920D+00,5.0710D+00,
     &4.6822D+00,4.2915D+00,3.9337D+00,3.5898D+00,3.2756D+00,2.9660D+00,
     &2.6817D+00,2.4150D+00,2.1724D+00,2.0348D+00,1.8961D+00,1.7457D+00,
     &1.6034D+00,1.4714D+00,1.3471D+00,1.2473D+00,1.1476D+00,1.0408D+00,
     &9.4083D-01,8.4932D-01,7.6350D-01,6.9606D-01,6.2897D-01,5.5833D-01,
     &4.9315D-01,4.3444D-01,3.8044D-01,3.3861D-01,2.9817D-01,2.5642D-01,
     &2.1917D-01,1.8685D-01,1.5838D-01,1.3752D-01,1.1831D-01,9.9987D-02,
     &8.5224D-02,7.6762D-02,6.9910D-02,6.5655D-02,6.2297D-02,6.0213D-02,
     &5.8897D-02,5.8096D-02,5.7624D-02,5.7400D-02,5.7322D-02,5.7351D-02,
     &5.7432D-02,5.7560D-02,5.7758D-02,6.0939D-02,6.3212D-02,6.6167D-02,
     &6.9884D-02,7.4560D-02,8.0552D-02,8.8432D-02,9.9242D-02,1.1491D-01,
     &1.3966D-01,1.8320D-01,0.0000D+00,0.0000D+00,7.4490D+00,6.8826D+00,
     &6.3540D+00,5.8477D+00,5.3805D+00,4.9187D+00,4.4884D+00,4.0843D+00,
     &3.7147D+00,3.3516D+00,3.0193D+00,2.7088D+00,2.4279D+00,2.2696D+00,
     &2.1091D+00,1.9368D+00,1.7739D+00,1.6237D+00,1.4821D+00,1.3692D+00,
     &1.2557D+00,1.1358D+00,1.0238D+00,9.2133D-01,8.2567D-01,7.5070D-01,
     &6.7656D-01,5.9850D-01,5.2688D-01,4.6263D-01,4.0371D-01,3.5842D-01/
      DATA (XSF_L(K),K= 1597, 1710) /
     &3.1427D-01,2.6933D-01,2.2930D-01,1.9466D-01,1.6427D-01,1.4208D-01,
     &1.2168D-01,1.0226D-01,8.6560D-02,7.7553D-02,7.0202D-02,6.5576D-02,
     &6.1860D-02,5.9487D-02,5.7920D-02,5.6852D-02,5.6166D-02,5.5736D-02,
     &5.5458D-02,5.5289D-02,5.5193D-02,5.5163D-02,5.5243D-02,5.7935D-02,
     &5.9740D-02,6.2111D-02,6.5158D-02,6.9050D-02,7.4078D-02,8.0683D-02,
     &8.9776D-02,1.0288D-01,1.2351D-01,1.5725D-01,0.0000D+00,0.0000D+00,
     &8.6044D+00,7.9255D+00,7.2940D+00,6.6940D+00,6.1391D+00,5.5940D+00,
     &5.0907D+00,4.6180D+00,4.1841D+00,3.7622D+00,3.3775D+00,3.0195D+00,
     &2.6967D+00,2.5153D+00,2.3331D+00,2.1364D+00,1.9521D+00,1.7815D+00,
     &1.6211D+00,1.4944D+00,1.3683D+00,1.2334D+00,1.1084D+00,9.9465D-01,
     &8.8864D-01,8.0585D-01,7.2432D-01,6.3866D-01,5.6038D-01,4.9058D-01,
     &4.2648D-01,3.7768D-01,3.3036D-01,2.8189D-01,2.3907D-01,2.0214D-01,
     &1.6987D-01,1.4635D-01,1.2479D-01,1.0428D-01,8.7748D-02,7.8203D-02,
     &7.0386D-02,6.5431D-02,6.1373D-02,5.8719D-02,5.6916D-02,5.5642D-02,
     &5.4751D-02,5.4118D-02,5.3653D-02,5.3314D-02,5.3067D-02,5.2897D-02,
     &5.2861D-02,5.5140D-02,5.6493D-02,5.8378D-02,6.0860D-02,6.4090D-02,
     &6.8261D-02,7.3828D-02,8.1439D-02,9.2423D-02,1.0952D-01,1.3424D-01,
     &0.0000D+00,0.0000D+00,9.7645D+00,8.9701D+00,8.2340D+00,7.5357D+00,
     &6.8926D+00,6.2607D+00,5.6834D+00,5.1374D+00,4.6459D+00,4.1625D+00/
      DATA (XSF_L(K),K= 1711, 1824) /
     &3.7261D+00,3.3206D+00,2.9567D+00,2.7529D+00,2.5476D+00,2.3274D+00,
     &2.1217D+00,1.9320D+00,1.7541D+00,1.6131D+00,1.4740D+00,1.3257D+00,
     &1.1879D+00,1.0631D+00,9.4732D-01,8.5726D-01,7.6844D-01,6.7586D-01,
     &5.9131D-01,5.1597D-01,4.4748D-01,3.9504D-01,3.4470D-01,2.9317D-01,
     &2.4779D-01,2.0880D-01,1.7478D-01,1.5007D-01,1.2748D-01,1.0600D-01,
     &8.8713D-02,7.8704D-02,7.0472D-02,6.5220D-02,6.0885D-02,5.7993D-02,
     &5.5967D-02,5.4536D-02,5.3470D-02,5.2665D-02,5.2054D-02,5.1577D-02,
     &5.1203D-02,5.0930D-02,5.0809D-02,5.2731D-02,5.3716D-02,5.5192D-02,
     &5.7203D-02,5.9902D-02,6.3412D-02,6.8123D-02,7.4602D-02,8.3905D-02,
     &9.8185D-02,1.1515D-01,0.0000D+00,0.0000D+00,1.1069D+01,1.0141D+01,
     &9.2840D+00,8.4741D+00,7.7316D+00,7.0038D+00,6.3364D+00,5.7137D+00,
     &5.1475D+00,4.6031D+00,4.1059D+00,3.6477D+00,3.2381D+00,3.0086D+00,
     &2.7788D+00,2.5333D+00,2.3033D+00,2.0926D+00,1.8951D+00,1.7404D+00,
     &1.5854D+00,1.4229D+00,1.2715D+00,1.1352D+00,1.0089D+00,9.1089D-01,
     &8.1457D-01,7.1424D-01,6.2332D-01,5.4229D-01,4.6872D-01,4.1295D-01,
     &3.5903D-01,3.0454D-01,2.5654D-01,2.1539D-01,1.7965D-01,1.5373D-01,
     &1.3011D-01,1.0766D-01,8.9530D-02,7.9108D-02,7.0483D-02,6.4943D-02,
     &6.0331D-02,5.7203D-02,5.4990D-02,5.3395D-02,5.2144D-02,5.1206D-02,
     &5.0454D-02,4.9840D-02,4.9351D-02,4.8978D-02,4.8801D-02,5.0351D-02/
      DATA (XSF_L(K),K= 1825, 1836) /
     &5.1017D-02,5.2127D-02,5.3737D-02,5.5934D-02,5.8835D-02,6.2800D-02,
     &6.8260D-02,7.6135D-02,8.7873D-02,0.0000D+00,0.0000D+00,0.0000D+00/
      DATA (XGF_L(K),K=    1,  114) /
     &1.0646D+00,1.0934D+00,1.1214D+00,1.1484D+00,1.1741D+00,1.1999D+00,
     &1.2242D+00,1.2466D+00,1.2676D+00,1.2873D+00,1.3042D+00,1.3194D+00,
     &1.3313D+00,1.3376D+00,1.3430D+00,1.3472D+00,1.3502D+00,1.3504D+00,
     &1.3501D+00,1.3478D+00,1.3430D+00,1.3356D+00,1.3267D+00,1.3149D+00,
     &1.3003D+00,1.2857D+00,1.2680D+00,1.2451D+00,1.2189D+00,1.1899D+00,
     &1.1575D+00,1.1282D+00,1.0947D+00,1.0543D+00,1.0121D+00,9.6983D-01,
     &9.2809D-01,8.9556D-01,8.6663D-01,8.4606D-01,8.4971D-01,8.7714D-01,
     &9.3569D-01,1.0140D+00,1.1325D+00,1.2706D+00,1.4268D+00,1.6005D+00,
     &1.7918D+00,2.0014D+00,2.2301D+00,2.4791D+00,2.7490D+00,3.0404D+00,
     &3.3541D+00,3.5718D+00,4.2579D+00,5.0478D+00,5.9674D+00,7.0458D+00,
     &8.3375D+00,9.9284D+00,1.1949D+01,1.4650D+01,1.8560D+01,2.5096D+01,
     &4.0067D+01,0.0000D+00,1.6404D+00,1.6723D+00,1.7014D+00,1.7287D+00,
     &1.7533D+00,1.7768D+00,1.7973D+00,1.8152D+00,1.8297D+00,1.8417D+00,
     &1.8498D+00,1.8540D+00,1.8544D+00,1.8526D+00,1.8489D+00,1.8424D+00,
     &1.8335D+00,1.8221D+00,1.8091D+00,1.7949D+00,1.7784D+00,1.7555D+00,
     &1.7310D+00,1.7034D+00,1.6713D+00,1.6428D+00,1.6093D+00,1.5680D+00,
     &1.5230D+00,1.4754D+00,1.4241D+00,1.3785D+00,1.3278D+00,1.2681D+00,
     &1.2068D+00,1.1462D+00,1.0867D+00,1.0400D+00,9.9665D-01,9.6041D-01,
     &9.4923D-01,9.6563D-01,1.0117D+00,1.0781D+00,1.1816D+00,1.3028D+00/
      DATA (XGF_L(K),K=  115,  228) /
     &1.4397D+00,1.5912D+00,1.7573D+00,1.9376D+00,2.1326D+00,2.3425D+00,
     &2.5677D+00,2.8078D+00,3.0611D+00,3.2398D+00,3.7904D+00,4.4126D+00,
     &5.1162D+00,5.9322D+00,6.8841D+00,8.0278D+00,9.4403D+00,1.1276D+01,
     &1.3844D+01,1.7948D+01,2.6821D+01,0.0000D+00,2.5295D+00,2.5563D+00,
     &2.5800D+00,2.5995D+00,2.6174D+00,2.6286D+00,2.6363D+00,2.6395D+00,
     &2.6379D+00,2.6306D+00,2.6184D+00,2.6000D+00,2.5768D+00,2.5598D+00,
     &2.5397D+00,2.5137D+00,2.4839D+00,2.4516D+00,2.4161D+00,2.3833D+00,
     &2.3459D+00,2.3000D+00,2.2499D+00,2.1966D+00,2.1407D+00,2.0900D+00,
     &2.0320D+00,1.9647D+00,1.8929D+00,1.8190D+00,1.7411D+00,1.6734D+00,
     &1.5997D+00,1.5142D+00,1.4279D+00,1.3438D+00,1.2617D+00,1.1967D+00,
     &1.1353D+00,1.0800D+00,1.0501D+00,1.0526D+00,1.0849D+00,1.1369D+00,
     &1.2228D+00,1.3250D+00,1.4410D+00,1.5691D+00,1.7085D+00,1.8587D+00,
     &2.0200D+00,2.1915D+00,2.3728D+00,2.5633D+00,2.7603D+00,2.9047D+00,
     &3.3315D+00,3.8026D+00,4.3243D+00,4.9121D+00,5.5828D+00,6.3648D+00,
     &7.3038D+00,8.4817D+00,1.0068D+01,1.2484D+01,1.7398D+01,0.0000D+00,
     &3.9781D+00,3.9859D+00,3.9880D+00,3.9845D+00,3.9763D+00,3.9582D+00,
     &3.9337D+00,3.9028D+00,3.8636D+00,3.8159D+00,3.7613D+00,3.6984D+00,
     &3.6287D+00,3.5836D+00,3.5326D+00,3.4703D+00,3.4046D+00,3.3350D+00,
     &3.2612D+00,3.1962D+00,3.1248D+00,3.0388D+00,2.9485D+00,2.8565D+00/
      DATA (XGF_L(K),K=  229,  342) /
     &2.7591D+00,2.6752D+00,2.5823D+00,2.4756D+00,2.3627D+00,2.2510D+00,
     &2.1352D+00,2.0365D+00,1.9308D+00,1.8097D+00,1.6896D+00,1.5737D+00,
     &1.4618D+00,1.3735D+00,1.2886D+00,1.2087D+00,1.1551D+00,1.1411D+00,
     &1.1545D+00,1.1903D+00,1.2550D+00,1.3356D+00,1.4282D+00,1.5306D+00,
     &1.6419D+00,1.7606D+00,1.8869D+00,2.0194D+00,2.1574D+00,2.2992D+00,
     &2.4432D+00,2.5568D+00,2.8674D+00,3.2008D+00,3.5626D+00,3.9572D+00,
     &4.3932D+00,4.8857D+00,5.4544D+00,6.1386D+00,7.0188D+00,8.2895D+00,
     &1.0709D+01,0.0000D+00,6.3697D+00,6.3265D+00,6.2740D+00,6.2091D+00,
     &6.1391D+00,6.0517D+00,5.9560D+00,5.8525D+00,5.7367D+00,5.6106D+00,
     &5.4709D+00,5.3235D+00,5.1695D+00,5.0724D+00,4.9662D+00,4.8411D+00,
     &4.7105D+00,4.5784D+00,4.4412D+00,4.3226D+00,4.1943D+00,4.0442D+00,
     &3.8903D+00,3.7360D+00,3.5773D+00,3.4420D+00,3.2967D+00,3.1301D+00,
     &2.9593D+00,2.7916D+00,2.6229D+00,2.4802D+00,2.3301D+00,2.1613D+00,
     &1.9957D+00,1.8382D+00,1.6875D+00,1.5691D+00,1.4545D+00,1.3433D+00,
     &1.2614D+00,1.2264D+00,1.2177D+00,1.2342D+00,1.2749D+00,1.3313D+00,
     &1.3987D+00,1.4740D+00,1.5559D+00,1.6431D+00,1.7346D+00,1.8295D+00,
     &1.9260D+00,2.0232D+00,2.1174D+00,2.2034D+00,2.4118D+00,2.6289D+00,
     &2.8563D+00,3.0948D+00,3.3486D+00,3.6231D+00,3.9250D+00,4.2677D+00,
     &4.6847D+00,5.2492D+00,6.2650D+00,0.0000D+00,9.3778D+00,9.2428D+00/
      DATA (XGF_L(K),K=  343,  456) /
     &9.0960D+00,8.9365D+00,8.7665D+00,8.5746D+00,8.3714D+00,8.1544D+00,
     &7.9343D+00,7.6882D+00,7.4352D+00,7.1724D+00,6.9027D+00,6.7360D+00,
     &6.5571D+00,6.3494D+00,6.1374D+00,5.9260D+00,5.7093D+00,5.5249D+00,
     &5.3270D+00,5.0995D+00,4.8700D+00,4.6419D+00,4.4114D+00,4.2173D+00,
     &4.0129D+00,3.7786D+00,3.5451D+00,3.3173D+00,3.0900D+00,2.9004D+00,
     &2.7040D+00,2.4853D+00,2.2734D+00,2.0742D+00,1.8851D+00,1.7372D+00,
     &1.5941D+00,1.4536D+00,1.3433D+00,1.2893D+00,1.2607D+00,1.2587D+00,
     &1.2792D+00,1.3153D+00,1.3616D+00,1.4149D+00,1.4736D+00,1.5361D+00,
     &1.6012D+00,1.6677D+00,1.7344D+00,1.7990D+00,1.8589D+00,1.9261D+00,
     &2.0646D+00,2.2044D+00,2.3456D+00,2.4882D+00,2.6342D+00,2.7823D+00,
     &2.9370D+00,3.1022D+00,3.2902D+00,3.5288D+00,3.9528D+00,0.0000D+00,
     &1.3926D+01,1.3617D+01,1.3298D+01,1.2959D+01,1.2612D+01,1.2230D+01,
     &1.1845D+01,1.1442D+01,1.1036D+01,1.0599D+01,1.0158D+01,9.7041D+00,
     &9.2562D+00,8.9827D+00,8.6974D+00,8.3632D+00,8.0255D+00,7.6946D+00,
     &7.3614D+00,7.0802D+00,6.7814D+00,6.4439D+00,6.1064D+00,5.7775D+00,
     &5.4468D+00,5.1723D+00,4.8858D+00,4.5617D+00,4.2425D+00,3.9361D+00,
     &3.6353D+00,3.3874D+00,3.1301D+00,2.8506D+00,2.5816D+00,2.3318D+00,
     &2.0965D+00,1.9136D+00,1.7368D+00,1.5622D+00,1.4211D+00,1.3452D+00,
     &1.2937D+00,1.2737D+00,1.2719D+00,1.2868D+00,1.3119D+00,1.3437D+00/
      DATA (XGF_L(K),K=  457,  570) /
     &1.3799D+00,1.4189D+00,1.4596D+00,1.5003D+00,1.5401D+00,1.5761D+00,
     &1.6073D+00,1.6574D+00,1.7377D+00,1.8158D+00,1.8902D+00,1.9601D+00,
     &2.0263D+00,2.0884D+00,2.1452D+00,2.1990D+00,2.2512D+00,2.3118D+00,
     &2.4354D+00,0.0000D+00,1.9256D+01,1.8699D+01,1.8142D+01,1.7563D+01,
     &1.6980D+01,1.6355D+01,1.5725D+01,1.5081D+01,1.4443D+01,1.3769D+01,
     &1.3097D+01,1.2422D+01,1.1755D+01,1.1358D+01,1.0937D+01,1.0454D+01,
     &9.9818D+00,9.5167D+00,9.0465D+00,8.6570D+00,8.2473D+00,7.7870D+00,
     &7.3320D+00,6.8911D+00,6.4569D+00,6.0969D+00,5.7223D+00,5.3051D+00,
     &4.8992D+00,4.5131D+00,4.1351D+00,3.8285D+00,3.5148D+00,3.1749D+00,
     &2.8517D+00,2.5534D+00,2.2748D+00,2.0598D+00,1.8527D+00,1.6465D+00,
     &1.4780D+00,1.3832D+00,1.3129D+00,1.2758D+00,1.2566D+00,1.2544D+00,
     &1.2628D+00,1.2778D+00,1.2971D+00,1.3186D+00,1.3412D+00,1.3637D+00,
     &1.3845D+00,1.4021D+00,1.4142D+00,1.4518D+00,1.4945D+00,1.5327D+00,
     &1.5661D+00,1.5941D+00,1.6160D+00,1.6309D+00,1.6386D+00,1.6381D+00,
     &1.6291D+00,1.6176D+00,1.6271D+00,0.0000D+00,2.5945D+01,2.5063D+01,
     &2.4160D+01,2.3234D+01,2.2336D+01,2.1370D+01,2.0417D+01,1.9450D+01,
     &1.8508D+01,1.7517D+01,1.6548D+01,1.5580D+01,1.4645D+01,1.4085D+01,
     &1.3496D+01,1.2836D+01,1.2181D+01,1.1547D+01,1.0921D+01,1.0404D+01,
     &9.8614D+00,9.2547D+00,8.6616D+00,8.0926D+00,7.5352D+00,7.0774D+00/
      DATA (XGF_L(K),K=  571,  684) /
     &6.6043D+00,6.0842D+00,5.5816D+00,5.1040D+00,4.6450D+00,4.2749D+00,
     &3.8995D+00,3.4941D+00,3.1134D+00,2.7651D+00,2.4423D+00,2.1941D+00,
     &1.9564D+00,1.7198D+00,1.5241D+00,1.4112D+00,1.3220D+00,1.2705D+00,
     &1.2348D+00,1.2175D+00,1.2113D+00,1.2119D+00,1.2167D+00,1.2238D+00,
     &1.2321D+00,1.2398D+00,1.2460D+00,1.2491D+00,1.2470D+00,1.2752D+00,
     &1.2894D+00,1.2998D+00,1.3055D+00,1.3049D+00,1.2991D+00,1.2860D+00,
     &1.2655D+00,1.2370D+00,1.1998D+00,1.1564D+00,1.1181D+00,0.0000D+00,
     &3.3362D+01,3.2051D+01,3.0740D+01,2.9429D+01,2.8133D+01,2.6758D+01,
     &2.5422D+01,2.4082D+01,2.2784D+01,2.1435D+01,2.0130D+01,1.8839D+01,
     &1.7597D+01,1.6865D+01,1.6098D+01,1.5241D+01,1.4397D+01,1.3587D+01,
     &1.2791D+01,1.2130D+01,1.1444D+01,1.0687D+01,9.9507D+00,9.2501D+00,
     &8.5659D+00,8.0104D+00,7.4390D+00,6.8118D+00,6.2125D+00,5.6506D+00,
     &5.1096D+00,4.6780D+00,4.2434D+00,3.7769D+00,3.3424D+00,2.9475D+00,
     &2.5842D+00,2.3061D+00,2.0409D+00,1.7770D+00,1.5572D+00,1.4290D+00,
     &1.3248D+00,1.2609D+00,1.2112D+00,1.1814D+00,1.1636D+00,1.1530D+00,
     &1.1469D+00,1.1433D+00,1.1407D+00,1.1378D+00,1.1337D+00,1.1269D+00,
     &1.1152D+00,1.1360D+00,1.1320D+00,1.1243D+00,1.1127D+00,1.0960D+00,
     &1.0739D+00,1.0461D+00,1.0122D+00,9.7100D-01,9.2292D-01,8.6909D-01,
     &8.1432D-01,0.0000D+00,4.2364D+01,4.0483D+01,3.8640D+01,3.6792D+01/
      DATA (XGF_L(K),K=  685,  798) /
     &3.4991D+01,3.3112D+01,3.1295D+01,2.9487D+01,2.7748D+01,2.5953D+01,
     &2.4235D+01,2.2543D+01,2.0935D+01,1.9990D+01,1.9011D+01,1.7921D+01,
     &1.6852D+01,1.5830D+01,1.4831D+01,1.4013D+01,1.3165D+01,1.2236D+01,
     &1.1337D+01,1.0485D+01,9.6616D+00,8.9943D+00,8.3137D+00,7.5711D+00,
     &6.8670D+00,6.2090D+00,5.5842D+00,5.0866D+00,4.5873D+00,4.0564D+00,
     &3.5646D+00,3.1234D+00,2.7185D+00,2.4107D+00,2.1172D+00,1.8273D+00,
     &1.5836D+00,1.4407D+00,1.3211D+00,1.2459D+00,1.1839D+00,1.1433D+00,
     &1.1153D+00,1.0949D+00,1.0794D+00,1.0667D+00,1.0555D+00,1.0443D+00,
     &1.0317D+00,1.0172D+00,9.9883D-01,1.0131D+00,9.9503D-01,9.7446D-01,
     &9.5064D-01,9.2316D-01,8.9156D-01,8.5528D-01,8.1439D-01,7.6837D-01,
     &7.1718D-01,6.6210D-01,6.0243D-01,0.0000D+00,5.2603D+01,5.0038D+01,
     &4.7540D+01,4.5053D+01,4.2652D+01,4.0175D+01,3.7784D+01,3.5407D+01,
     &3.3154D+01,3.0851D+01,2.8651D+01,2.6507D+01,2.4488D+01,2.3310D+01,
     &2.2084D+01,2.0735D+01,1.9418D+01,1.8166D+01,1.6951D+01,1.5960D+01,
     &1.4935D+01,1.3817D+01,1.2742D+01,1.1732D+01,1.0759D+01,9.9749D+00,
     &9.1794D+00,8.3186D+00,7.5044D+00,6.7510D+00,6.0386D+00,5.4762D+00,
     &4.9137D+00,4.3200D+00,3.7728D+00,3.2842D+00,2.8391D+00,2.5026D+00,
     &2.1835D+00,1.8677D+00,1.6033D+00,1.4461D+00,1.3138D+00,1.2277D+00,
     &1.1557D+00,1.1057D+00,1.0689D+00,1.0407D+00,1.0176D+00,9.9768D-01/
      DATA (XGF_L(K),K=  799,  912) /
     &9.7951D-01,9.6199D-01,9.4331D-01,9.2359D-01,9.0058D-01,9.0921D-01,
     &8.8156D-01,8.5244D-01,8.2081D-01,7.8702D-01,7.5025D-01,7.1005D-01,
     &6.6667D-01,6.1984D-01,5.6969D-01,5.1748D-01,4.5895D-01,0.0000D+00,
     &6.3459D+01,6.0127D+01,5.6900D+01,5.3695D+01,5.0615D+01,4.7464D+01,
     &4.4440D+01,4.1483D+01,3.8684D+01,3.5826D+01,3.3122D+01,3.0500D+01,
     &2.8040D+01,2.6617D+01,2.5143D+01,2.3518D+01,2.1950D+01,2.0455D+01,
     &1.9011D+01,1.7842D+01,1.6646D+01,1.5337D+01,1.4094D+01,1.2920D+01,
     &1.1799D+01,1.0903D+01,9.9940D+00,9.0166D+00,8.0967D+00,7.2512D+00,
     &6.4551D+00,5.8279D+00,5.2081D+00,4.5519D+00,3.9568D+00,3.4237D+00,
     &2.9425D+00,2.5798D+00,2.2371D+00,1.8995D+00,1.6161D+00,1.4477D+00,
     &1.3046D+00,1.2096D+00,1.1285D+00,1.0709D+00,1.0274D+00,9.9290D-01,
     &9.6399D-01,9.3860D-01,9.1550D-01,8.9324D-01,8.7036D-01,8.4674D-01,
     &8.2129D-01,8.2506D-01,7.9094D-01,7.5633D-01,7.2031D-01,6.8307D-01,
     &6.4387D-01,6.0237D-01,5.5907D-01,5.1344D-01,4.6618D-01,4.1810D-01,
     &3.6329D-01,0.0000D+00,7.9498D+01,7.4941D+01,7.0580D+01,6.6266D+01,
     &6.2169D+01,5.8002D+01,5.4045D+01,5.0164D+01,4.6539D+01,4.2847D+01,
     &3.9386D+01,3.6065D+01,3.2968D+01,3.1180D+01,2.9347D+01,2.7330D+01,
     &2.5394D+01,2.3566D+01,2.1811D+01,2.0388D+01,1.8944D+01,1.7368D+01,
     &1.5877D+01,1.4488D+01,1.3164D+01,1.2111D+01,1.1051D+01,9.9162D+00/
      DATA (XGF_L(K),K=  913, 1026) /
     &8.8542D+00,7.8839D+00,6.9777D+00,6.2689D+00,5.5695D+00,4.8410D+00,
     &4.1789D+00,3.5909D+00,3.0635D+00,2.6689D+00,2.2973D+00,1.9324D+00,
     &1.6270D+00,1.4446D+00,1.2882D+00,1.1839D+00,1.0926D+00,1.0266D+00,
     &9.7585D-01,9.3473D-01,8.9976D-01,8.6898D-01,8.4068D-01,8.1374D-01,
     &7.8714D-01,7.6011D-01,7.3262D-01,7.3148D-01,6.9170D-01,6.5270D-01,
     &6.1357D-01,5.7426D-01,5.3417D-01,4.9316D-01,4.5166D-01,4.0914D-01,
     &3.6649D-01,3.2429D-01,2.7651D-01,0.0000D+00,9.7091D+01,9.1127D+01,
     &8.5440D+01,7.9869D+01,7.4603D+01,6.9275D+01,6.4220D+01,5.9343D+01,
     &5.4780D+01,5.0195D+01,4.5912D+01,4.1816D+01,3.8028D+01,3.5857D+01,
     &3.3637D+01,3.1205D+01,2.8880D+01,2.6695D+01,2.4601D+01,2.2923D+01,
     &2.1219D+01,1.9374D+01,1.7634D+01,1.6009D+01,1.4488D+01,1.3276D+01,
     &1.2064D+01,1.0772D+01,9.5709D+00,8.4795D+00,7.4649D+00,6.6775D+00,
     &5.9046D+00,5.1015D+00,4.3733D+00,3.7372D+00,3.1677D+00,2.7434D+00,
     &2.3459D+00,1.9566D+00,1.6317D+00,1.4368D+00,1.2699D+00,1.1572D+00,
     &1.0581D+00,9.8558D-01,9.2913D-01,8.8297D-01,8.4349D-01,8.0862D-01,
     &7.7667D-01,7.4686D-01,7.1760D-01,6.8906D-01,6.6005D-01,6.5493D-01,
     &6.1209D-01,5.7080D-01,5.3038D-01,4.9085D-01,4.5137D-01,4.1231D-01,
     &3.7316D-01,3.3442D-01,2.9613D-01,2.5928D-01,2.1912D-01,0.0000D+00,
     &1.1660D+02,1.0899D+02,1.0178D+02,9.4752D+01,8.8142D+01,8.1480D+01/
      DATA (XGF_L(K),K= 1027, 1140) /
     &7.5219D+01,6.9198D+01,6.3578D+01,5.7986D+01,5.2800D+01,4.7867D+01,
     &4.3328D+01,4.0736D+01,3.8088D+01,3.5213D+01,3.2469D+01,2.9907D+01,
     &2.7451D+01,2.5501D+01,2.3516D+01,2.1392D+01,1.9391D+01,1.7546D+01,
     &1.5800D+01,1.4426D+01,1.3057D+01,1.1607D+01,1.0266D+01,9.0517D+00,
     &7.9294D+00,7.0617D+00,6.2165D+00,5.3397D+00,4.5572D+00,3.8687D+00,
     &3.2598D+00,2.8078D+00,2.3859D+00,1.9745D+00,1.6317D+00,1.4267D+00,
     &1.2497D+00,1.1305D+00,1.0247D+00,9.4657D-01,8.8556D-01,8.3542D-01,
     &7.9253D-01,7.5465D-01,7.2037D-01,6.8840D-01,6.5775D-01,6.2793D-01,
     &5.9852D-01,5.9015D-01,5.4553D-01,5.0339D-01,4.6306D-01,4.2411D-01,
     &3.8622D-01,3.4909D-01,3.1294D-01,2.7773D-01,2.4373D-01,2.1150D-01,
     &1.7848D-01,0.0000D+00,1.3738D+02,1.2796D+02,1.1904D+02,1.1042D+02,
     &1.0233D+02,9.4222D+01,8.6662D+01,7.9409D+01,7.2655D+01,6.6001D+01,
     &5.9833D+01,5.4007D+01,4.8672D+01,4.5642D+01,4.2552D+01,3.9214D+01,
     &3.6040D+01,3.3082D+01,3.0272D+01,2.8026D+01,2.5779D+01,2.3361D+01,
     &2.1093D+01,1.9009D+01,1.7062D+01,1.5526D+01,1.4003D+01,1.2396D+01,
     &1.0916D+01,9.5845D+00,8.3611D+00,7.4188D+00,6.5021D+00,5.5589D+00,
     &4.7169D+00,3.9865D+00,3.3389D+00,2.8617D+00,2.4178D+00,1.9872D+00,
     &1.6283D+00,1.4143D+00,1.2296D+00,1.1049D+00,9.9315D-01,9.1079D-01,
     &8.4623D-01,7.9317D-01,7.4768D-01,7.0802D-01,6.7178D-01,6.3836D-01/
      DATA (XGF_L(K),K= 1141, 1254) /
     &6.0703D-01,5.7658D-01,5.4733D-01,5.3630D-01,4.9100D-01,4.4879D-01,
     &4.0920D-01,3.7138D-01,3.3521D-01,3.0054D-01,2.6721D-01,2.3523D-01,
     &2.0485D-01,1.7634D-01,1.4852D-01,0.0000D+00,1.6103D+02,1.4938D+02,
     &1.3848D+02,1.2798D+02,1.1818D+02,1.0840D+02,9.9309D+01,9.0651D+01,
     &8.2647D+01,7.4733D+01,6.7469D+01,6.0672D+01,5.4433D+01,5.0913D+01,
     &4.7343D+01,4.3482D+01,3.9833D+01,3.6452D+01,3.3242D+01,3.0689D+01,
     &2.8134D+01,2.5404D+01,2.2863D+01,2.0531D+01,1.8362D+01,1.6652D+01,
     &1.4967D+01,1.3197D+01,1.1573D+01,1.0120D+01,8.7877D+00,7.7679D+00,
     &6.7819D+00,5.7685D+00,4.8731D+00,4.0967D+00,3.4122D+00,2.9097D+00,
     &2.4451D+00,1.9953D+00,1.6222D+00,1.3995D+00,1.2076D+00,1.0771D+00,
     &9.6151D-01,8.7563D-01,8.0819D-01,7.5269D-01,7.0548D-01,6.6395D-01,
     &6.2666D-01,5.9253D-01,5.6034D-01,5.3005D-01,5.0122D-01,4.8790D-01,
     &4.4273D-01,4.0115D-01,3.6251D-01,3.2632D-01,2.9224D-01,2.5988D-01,
     &2.2931D-01,2.0039D-01,1.7324D-01,1.4805D-01,1.2201D-01,0.0000D+00,
     &1.8591D+02,1.7193D+02,1.5886D+02,1.4632D+02,1.3469D+02,1.2310D+02,
     &1.1237D+02,1.0218D+02,9.2839D+01,8.3643D+01,7.5256D+01,6.7382D+01,
     &6.0231D+01,5.6204D+01,5.2127D+01,4.7743D+01,4.3601D+01,3.9784D+01,
     &3.6172D+01,3.3310D+01,3.0455D+01,2.7410D+01,2.4579D+01,2.2009D+01,
     &1.9599D+01,1.7727D+01,1.5886D+01,1.3956D+01,1.2193D+01,1.0620D+01/
      DATA (XGF_L(K),K= 1255, 1368) /
     &9.1866D+00,8.0925D+00,7.0383D+00,5.9623D+00,5.0119D+00,4.1917D+00,
     &3.4750D+00,2.9503D+00,2.4663D+00,1.9999D+00,1.6141D+00,1.3840D+00,
     &1.1856D+00,1.0518D+00,9.3192D-01,8.4324D-01,7.7348D-01,7.1642D-01,
     &6.6779D-01,6.2531D-01,5.8732D-01,5.5231D-01,5.2039D-01,4.9037D-01,
     &4.6218D-01,4.4711D-01,4.0225D-01,3.6159D-01,3.2438D-01,2.8982D-01,
     &2.5765D-01,2.2765D-01,1.9954D-01,1.7331D-01,1.4889D-01,1.2621D-01,
     &9.6984D-02,0.0000D+00,2.1269D+02,1.9609D+02,1.8060D+02,1.6582D+02,
     &1.5214D+02,1.3863D+02,1.2613D+02,1.1431D+02,1.0351D+02,9.2957D+01,
     &8.3294D+01,7.4318D+01,6.6188D+01,6.1617D+01,5.7019D+01,5.2073D+01,
     &4.7428D+01,4.3153D+01,3.9122D+01,3.5941D+01,3.2764D+01,2.9404D+01,
     &2.6282D+01,2.3458D+01,2.0836D+01,1.8796D+01,1.6786D+01,1.4693D+01,
     &1.2792D+01,1.1101D+01,9.5678D+00,8.4010D+00,7.2773D+00,6.1402D+00,
     &5.1403D+00,4.2791D+00,3.5311D+00,2.9851D+00,2.4835D+00,2.0017D+00,
     &1.6039D+00,1.3677D+00,1.1646D+00,1.0265D+00,9.0375D-01,8.1271D-01,
     &7.4135D-01,6.8280D-01,6.3328D-01,5.9018D-01,5.5184D-01,5.1677D-01,
     &4.8494D-01,4.5537D-01,4.2797D-01,4.1146D-01,3.6736D-01,3.2788D-01,
     &2.9207D-01,2.5923D-01,2.2901D-01,2.0110D-01,1.7527D-01,1.5131D-01,
     &1.2926D-01,1.0839D-01,6.9776D-02,0.0000D+00,2.4043D+02,2.2104D+02,
     &2.0300D+02,1.8582D+02,1.7003D+02,1.5443D+02,1.4007D+02,1.2658D+02/
      DATA (XGF_L(K),K= 1369, 1482) /
     &1.1426D+02,1.0227D+02,9.1332D+01,8.1197D+01,7.2119D+01,6.6989D+01,
     &6.1846D+01,5.6342D+01,5.1188D+01,4.6448D+01,4.2002D+01,3.8498D+01,
     &3.5016D+01,3.1335D+01,2.7931D+01,2.4848D+01,2.2009D+01,1.9797D+01,
     &1.7637D+01,1.5389D+01,1.3354D+01,1.1550D+01,9.9187D+00,8.6824D+00,
     &7.4988D+00,6.3022D+00,5.2549D+00,4.3589D+00,3.5788D+00,3.0139D+00,
     &2.4962D+00,2.0005D+00,1.5931D+00,1.3514D+00,1.1435D+00,1.0028D+00,
     &8.7751D-01,7.8479D-01,7.1218D-01,6.5272D-01,6.0250D-01,5.5920D-01,
     &5.2061D-01,4.8590D-01,4.5422D-01,4.2519D-01,3.9858D-01,3.8094D-01,
     &3.3789D-01,2.9975D-01,2.6524D-01,2.3401D-01,2.0560D-01,1.7956D-01,
     &1.5565D-01,1.3374D-01,1.1354D-01,9.4096D-02,3.9275D-02,0.0000D+00,
     &2.8195D+02,2.5830D+02,2.3640D+02,2.1554D+02,1.9645D+02,1.7774D+02,
     &1.6058D+02,1.4448D+02,1.2990D+02,1.1575D+02,1.0299D+02,9.1121D+01,
     &8.0574D+01,7.4642D+01,6.8724D+01,6.2402D+01,5.6498D+01,5.1101D+01,
     &4.6042D+01,4.2081D+01,3.8152D+01,3.4014D+01,3.0201D+01,2.6780D+01,
     &2.3611D+01,2.1171D+01,1.8789D+01,1.6329D+01,1.4107D+01,1.2148D+01,
     &1.0386D+01,9.0557D+00,7.7874D+00,6.5118D+00,5.4006D+00,4.4539D+00,
     &3.6370D+00,3.0467D+00,2.5088D+00,1.9959D+00,1.5762D+00,1.3274D+00,
     &1.1142D+00,9.7065D-01,8.4265D-01,7.4825D-01,6.7451D-01,6.1445D-01,
     &5.6374D-01,5.2024D-01,4.8166D-01,4.4741D-01,4.1643D-01,3.8830D-01/
      DATA (XGF_L(K),K= 1483, 1596) /
     &3.6282D-01,3.4411D-01,3.0249D-01,2.6607D-01,2.3369D-01,2.0474D-01,
     &1.7852D-01,1.5489D-01,1.3341D-01,1.1384D-01,9.5862D-02,7.7509D-02,
     &0.0000D+00,0.0000D+00,3.2379D+02,2.9556D+02,2.6960D+02,2.4513D+02,
     &2.2265D+02,2.0073D+02,1.8071D+02,1.6202D+02,1.4515D+02,1.2887D+02,
     &1.1419D+02,1.0071D+02,8.8650D+01,8.1931D+01,7.5233D+01,6.8140D+01,
     &6.1510D+01,5.5467D+01,4.9832D+01,4.5419D+01,4.1070D+01,3.6493D+01,
     &3.2295D+01,2.8536D+01,2.5086D+01,2.2426D+01,1.9846D+01,1.7175D+01,
     &1.4781D+01,1.2681D+01,1.0797D+01,9.3831D+00,8.0380D+00,6.6897D+00,
     &5.5221D+00,4.5337D+00,3.6831D+00,3.0714D+00,2.5159D+00,1.9884D+00,
     &1.5586D+00,1.3048D+00,1.0886D+00,9.4191D-01,8.1217D-01,7.1679D-01,
     &6.4238D-01,5.8194D-01,5.3136D-01,4.8766D-01,4.4965D-01,4.1594D-01,
     &3.8570D-01,3.5847D-01,3.3403D-01,3.1456D-01,2.7454D-01,2.3977D-01,
     &2.0922D-01,1.8216D-01,1.5795D-01,1.3622D-01,1.1669D-01,9.9012D-02,
     &8.2668D-02,6.4604D-02,0.0000D+00,0.0000D+00,3.7071D+02,3.3727D+02,
     &3.0660D+02,2.7790D+02,2.5169D+02,2.2608D+02,2.0283D+02,1.8123D+02,
     &1.6179D+02,1.4311D+02,1.2635D+02,1.1097D+02,9.7357D+01,8.9759D+01,
     &8.2263D+01,7.4239D+01,6.6821D+01,6.0073D+01,5.3813D+01,4.8927D+01,
     &4.4114D+01,3.9072D+01,3.4471D+01,3.0351D+01,2.6592D+01,2.3699D+01,
     &2.0903D+01,1.8031D+01,1.5459D+01,1.3211D+01,1.1204D+01,9.7024D+00/
      DATA (XGF_L(K),K= 1597, 1710) /
     &8.2828D+00,6.8644D+00,5.6367D+00,4.6059D+00,3.7241D+00,3.0915D+00,
     &2.5189D+00,1.9786D+00,1.5396D+00,1.2816D+00,1.0611D+00,9.1306D-01,
     &7.8207D-01,6.8594D-01,6.1118D-01,5.5075D-01,5.0031D-01,4.5732D-01,
     &4.1996D-01,3.8671D-01,3.5732D-01,3.3101D-01,3.0775D-01,2.8769D-01,
     &2.4931D-01,2.1637D-01,1.8763D-01,1.6241D-01,1.4002D-01,1.2013D-01,
     &1.0238D-01,8.6311D-02,7.1348D-02,5.2982D-02,0.0000D+00,0.0000D+00,
     &4.2142D+02,3.8237D+02,3.4660D+02,3.1292D+02,2.8259D+02,2.5300D+02,
     &2.2626D+02,2.0148D+02,1.7927D+02,1.5797D+02,1.3896D+02,1.2163D+02,
     &1.0632D+02,9.7858D+01,8.9366D+01,8.0488D+01,7.2234D+01,6.4771D+01,
     &5.7843D+01,5.2468D+01,4.7182D+01,4.1663D+01,3.6633D+01,3.2165D+01,
     &2.8082D+01,2.4971D+01,2.1960D+01,1.8866D+01,1.6118D+01,1.3723D+01,
     &1.1595D+01,1.0008D+01,8.5101D+00,7.0232D+00,5.7443D+00,4.6705D+00,
     &3.7584D+00,3.1066D+00,2.5189D+00,1.9659D+00,1.5193D+00,1.2575D+00,
     &1.0346D+00,8.8517D-01,7.5338D-01,6.5695D-01,5.8219D-01,5.2200D-01,
     &4.7218D-01,4.2954D-01,3.9258D-01,3.6043D-01,3.3190D-01,3.0663D-01,
     &2.8431D-01,2.6413D-01,2.2746D-01,1.9612D-01,1.6912D-01,1.4557D-01,
     &1.2488D-01,1.0660D-01,9.0362D-02,7.5731D-02,6.1890D-02,4.2720D-02,
     &0.0000D+00,0.0000D+00,4.7166D+02,4.2676D+02,3.8580D+02,3.4749D+02,
     &3.1273D+02,2.7927D+02,2.4899D+02,2.2108D+02,1.9611D+02,1.7230D+02/
      DATA (XGF_L(K),K= 1711, 1824) /
     &1.5107D+02,1.3178D+02,1.1483D+02,1.0548D+02,9.6179D+01,8.6383D+01,
     &7.7331D+01,6.9156D+01,6.1613D+01,5.5763D+01,5.0019D+01,4.4056D+01,
     &3.8633D+01,3.3819D+01,2.9446D+01,2.6108D+01,2.2889D+01,1.9617D+01,
     &1.6706D+01,1.4179D+01,1.1938D+01,1.0276D+01,8.7112D+00,7.1630D+00,
     &5.8345D+00,4.7275D+00,3.7856D+00,3.1171D+00,2.5164D+00,1.9532D+00,
     &1.4997D+00,1.2350D+00,1.0108D+00,8.6027D-01,7.2804D-01,6.3166D-01,
     &5.5726D-01,4.9745D-01,4.4802D-01,4.0623D-01,3.7002D-01,3.3850D-01,
     &3.1081D-01,2.8644D-01,2.6509D-01,2.4476D-01,2.0951D-01,1.7979D-01,
     &1.5426D-01,1.3217D-01,1.1290D-01,9.5951D-02,8.0975D-02,6.7483D-02,
     &5.4483D-02,3.4309D-02,0.0000D+00,0.0000D+00,5.2745D+02,4.7595D+02,
     &4.2900D+02,3.8543D+02,3.4589D+02,3.0795D+02,2.7377D+02,2.4235D+02,
     &2.1434D+02,1.8771D+02,1.6408D+02,1.4266D+02,1.2392D+02,1.1358D+02,
     &1.0335D+02,9.2593D+01,8.2702D+01,7.3780D+01,6.5553D+01,5.9207D+01,
     &5.2983D+01,4.6535D+01,4.0700D+01,3.5531D+01,3.0842D+01,2.7278D+01,
     &2.3855D+01,2.0386D+01,1.7301D+01,1.4635D+01,1.2282D+01,1.0538D+01,
     &8.9065D+00,7.2932D+00,5.9178D+00,4.7769D+00,3.8086D+00,3.1240D+00,
     &2.5114D+00,1.9387D+00,1.4794D+00,1.2125D+00,9.8604D-01,8.3538D-01,
     &7.0309D-01,6.0683D-01,5.3289D-01,4.7378D-01,4.2493D-01,3.8387D-01,
     &3.4846D-01,3.1778D-01,2.9097D-01,2.6744D-01,2.4699D-01,2.2688D-01/
      DATA (XGF_L(K),K= 1825, 1836) /
     &1.9308D-01,1.6489D-01,1.4079D-01,1.2009D-01,1.0214D-01,8.6447D-02,
     &7.2603D-02,6.0131D-02,4.7893D-02,2.6613D-02,0.0000D+00,0.0000D+00/

*
      X = Xinp
*...CHECK OF X AND Q2 VALUES :
      IF ( (X.LT.0.99D-9) .OR. (X.GT.1.D0) ) THEN
*        WRITE(6,91) X
*  91     FORMAT (2X,'PHO_DOR98LO: x out of range',1p,E12.4)
         X = 0.99D-9
*        STOP
      ENDIF

      Q2 = Q2inp
      IF ( (Q2.LT.0.799D0) .OR. (Q2.GT.1.E6) ) THEN
*        WRITE(6,92) Q2
*  92     FORMAT (2X,'PHO_DOR98LO: Q2 out of range',1p,E12.4)
         Q2 = 0.99D6
*        STOP
      ENDIF

*
*...INTERPOLATION :
      NA(1) = NX
      NA(2) = NQ
      XT(1) = DLOG(X)
      XT(2) = DLOG(Q2)
      X1 = 1.D0- X
      XV = X**0.5D0
      XS = X**(-0.2D0)
      UV = SIB_DBFINT(NARG,XT,NA,ARRF,XUVF) * X1**3 * XV
      DV = SIB_DBFINT(NARG,XT,NA,ARRF,XDVF) * X1**4 * XV
      DE = SIB_DBFINT(NARG,XT,NA,ARRF,XDEF) * X1**7 * XV
      UD = SIB_DBFINT(NARG,XT,NA,ARRF,XUDF) * X1**7 * XS
      US = 0.5D0 * (UD - DE)
      DS = 0.5D0 * (UD + DE)
      SS = SIB_DBFINT(NARG,XT,NA,ARRF,XSF)  * X1**7 * XS
      GL = SIB_DBFINT(NARG,XT,NA,ARRF,XGF)  * X1**5 * XS

      END
