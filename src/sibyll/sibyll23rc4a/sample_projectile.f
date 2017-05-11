      SUBROUTINE SAMPLE_projectile
     +     (KID,NINT,LRMNT,XCHG,XJET,XX,PX,PY,IFL,KID1,IREJ)
C...  Subroutine to sample sea and valence quarks in a hadron.
C.    variables are stored in xx,px,py and ifl arrays.
C.    for each interaction the hadron undergoes there is one 
C.    pair of partons attached to the ends of two strings
C.    (one cut pomeron)
C.    In addition flavor and momentum may be set aside for the remnant
C.    arrays are filled: rmnt1,rmnt2, c.str1,c.str2, etc..
C.    i.e. positions 1 and 2 are reserved for remnant.
C.
C.    Input: Nint  : number of interactions the hadron takes part in
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
      SAVE

C     include COMMONs
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      INCLUDE 'sib_cutoff_cmmn.inc'
      INCLUDE 'sib_nw_prm.inc'
      INCLUDE 'sib_run_cmmn.inc'
      INCLUDE 'sib_mass_cmmn.inc'
      INCLUDE 'sib_chp_cmmn.inc'
      INCLUDE 'sib_utl_cmmn.inc'
      INCLUDE 'sib_cnt_cmmn.inc'

C     input type declarations
      INTEGER KID,NINT,LRMNT
      DOUBLE PRECISION XCHG,XJET
      
C     output type declarations
      DOUBLE PRECISION XX,PX,PY
      INTEGER IFL,KID1,IREJ
      DIMENSION XX(2*NW_max+2),PX(2*NW_max+2),PY(2*NW_max+2),
     &     IFL(2*NW_max+2)

c     local type declarations
      INTEGER ICNT1,ICNT2,J,JJ,j1,j2,j3,j4,KRMNT,IRNK,
     &     IDXVAL,IDX,ISWTD,i,IFLS,NVAL,NSEA,IR,IDUM,IDUM2,KIDA,IMRG2HAD
      DATA ICNT1,ICNT2 /0,0/
      DOUBLE PRECISION XSEAJET,XVAL,XMINA,XMINA_SEA,GAMMA,XREM,XMINA2,
     &     XMAX2,ALPHA,XM2DIS,ASUP,XMAX,XQM,S_RNDM,
     &     CHIDIS,CHI,GAMDIQ,XSUPP,PAR53_def,PAR5_def,PAR6_def,PAR7_def,
     &     PAR143_def,XSUM,STR_mass,PTS,XSCL
      
C..   initialization
      ITRY(3) = 0
      XVAL = ZERO
      XSCL = ONE
      XSEAJET = ZERO
      XSUM = ZERO
      DO J=1,Nint               ! zero arrays
         j1 = 1+2*(j-1)
         j2 = j1 + 1
         j3 = 3+2*(j-1)
         j4 = j3 + 1
         XX(j1) = ZERO
         XX(j2) = ZERO
         XX(j3) = ZERO
         XX(j4) = ZERO
         PX(j1) = ZERO
         PX(j2) = ZERO
         PX(j3) = ZERO
         PX(j4) = ZERO
      ENDDO

      KRMNT = MIN(LRMNT,1)

      IF(ndebug.gt.3) 
     +     WRITE(LUN,*)
     +     ' SAMPLE_projectile: KID,Nint,KRMNT,XCHG,XJET,IREJ',
     +     KID,Nint,KRMNT,XCHG,XJET,IREJ

      KID1 = KID
      KIDA = IABS(KID)
      
c     number of valence quarks to sample
c     if remnant is resolved (krmnt=1) no valence pair needed
      Nval = 2*(1-KRMNT)

c     number of sea quarks to sample (one pair per interaction)
c     if remnant is not resolved then on pair less is needed 
c     (valence pair takes role of one sea pair)
      Nsea = 2*(NINT-(1-KRMNT))

      IF(ndebug.gt.3) 
     +     WRITE(LUN,*)
     +     ' SAMPLE_projectile: number of partons to sample ',
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
            WRITE(LUN,*) 'SAMPLE_projectile: trials exceeded! return..'
            WRITE(LUN,*)
     +           ' KID,Nint,KRMNT,XCHG,XJET,XVAL,IREJ,NCALL',
     +           KID,Nint,KRMNT,XCHG,XJET,XVAL,IREJ,NCALL           
         ENDIF
         PAR(53) = PAR53_def
         RETURN 
      ENDIF

C...  kinematic limits
 22   XSEAJET = XJET
      IF(KRMNT.eq.0)THEN
c     minimal momentum fraction for valences
         XMINA = TWO*STR_mass_val/SQS
c     default for valence quarks: 0.35 , xmin@10GeV = 0.07
c     taken from COMMON s_cutoff
         IF(ISTR(KIDA)*IBAR(KIDA).ne.0)
     &        XMINA = TWO*STR_mass_val_hyp/SQS
      ELSE
         SELECT CASE(IPAR(47))
         CASE(4,5,6)
c     no valence sampling model
c     if remnant present then the minimal remnant mass has to be provided
            XMINA = PAR(96)*AM(ABS(KID))/SQS
         CASE(0,1,2,3)
c     valences sampled, even if combined again in remnant
            XMINA = TWO*STR_mass_val/SQS
         CASE(7)
c     minimal remnant mass not requiered,
c     mass is taken from central strings anyway..
            XMINA = AM(ABS(KID))/SQS
         END SELECT
      ENDIF
         
c     minimal momentum fraction for sea partons
      SELECT CASE(IPAR(47))
      CASE(0,3)
c     same as valence quarks
         STR_mass = STR_mass_val
      CASE(1,2,5,6,7)
c     set by parameter
         STR_mass = PAR(87)
      CASE(4)
c     same as soft minijets
         STR_mass = STR_mass_sea
      END SELECT
      IF(IPAR(72).eq.2.and.NINT.gt.1)THEN
         STR_mass = STR_mass * PAR(118)
      ENDIF
      XMINA_SEA = TWO*STR_mass/SQS
c     default for sea quarks: 1.0 , xmin@10GeV = 0.2
c     taken from COMMON s_cutoff or s_cflafr
c     should be the same as min. string mass in SAMPLE_soft !

c     dependence on number of interactions
      IF(IPAR(72).eq.1.and.NINT.gt.1)THEN
         XMINA_SEA = XMINA_SEA * PAR(118)
      ENDIF

C..   check if enough energy left to sample all partons
      IF (ONE-XJET.LT.(Nsea*XMINA_SEA+2*XMINA))THEN
         ICNT2 = ICNT2 + 1
         IF(ICNT2.le.10)THEN
            IF(NDEBUG.gt.3)THEN
               write(lun,*)' SAMPLE_projectile: rejection!' 
               write(lun,*)'  too little energy to sample all partons!'
               write(lun,*)' (NW,Ntot,Nval,Nsea,XMIN,XMIN*N',
     &              'XREM,XALL,NCALL,ICNT:)',nint,nval+nsea,Nval,nsea,
     &              2*xmina,nsea*xmina_sea,one-xjet,
     &              Nsea*XMINA_SEA+2*XMINA,NCALL,ICNT2
               IF(ICNT2.eq.10) write(lun,*)'last warning ! good luck..'
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
      SELECT CASE(IPAR(47))
c     sample from combined sea distribution
C     (1-x)**gamma / sqrt( x**2 + 2*m**2/s )
      CASE(0,3,4,5,7)
         GAMMA = PAR(103)
         IF(IPAR(73).eq.1.and.NINT.gt.1) GAMMA = PAR(119)
         CALL SAMPLE_SEA_TOT
     &        (KRMNT,NINT,NSEA,GAMMA,XJET,STR_MASS,XSEAJET,XX)

c     sample from 1/x individually then reject if too large
      CASE(1)
         XREM = ZERO
         XMINA2 = XMINA_SEA ** 2
         XMAX2 = 0.8D0**2
         ALPHA = ONE
         DO WHILE ( XREM .lt. 2*XMINA )
            XREM = ONE-XJET
            IF(NDEBUG.gt.3)
     &           WRITE(LUN,*) 'N,XREM,XMINA,XMIN2,XMAX2,ALPHA',
     &           Nsea,XREM,XMINA_SEA,XMINA2,XMAX2,ALPHA
            DO j=1,Nsea
               jj = 2 + j
               IF(KRMNT.eq.0) jj = 4+j
               XX(jj) = XM2DIS(XMINA2,XMAX2,ALPHA)
               IF(NDEBUG.gt.3) 
     &           WRITE(LUN,*) 'J,X,XREM',JJ,XX(JJ),XREM
               XREM = XREM - XX(jj)
            ENDDO
         ENDDO
         XSEAJET = ONE-XREM

c     sample from (1-x)**b / x with common mass constraint
      CASE(2,6)
         XREM = ONE-XJET
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

         XSEAJET = ONE-XREM

      END SELECT

C...  sample sea flavor: u,d,s,c
c     write to ifl after valences..
      DO J=1+Nval/2,Nint
         j3 = 3+2*(j-1)
         j4 = j3 + 1
c     turn on strange sea..
         IF(IPAR(29).eq.1)THEN
            IF(IPAR(69).ne.0)THEN
c     sample asymmetric u,d
               IFL(j3) = MIN(2,1+INT((TWO+PAR(114))*S_RNDM(KID)))
c     sample strange
               IFLS = 3*(INT((2+PAR(43))*S_RNDM(j3))/2)
               IFL(j3) = MAX(IFL(j3),IFLS)
            else
               IFL(j3) = 1+INT((TWO+PAR(43))*S_RNDM(KID))
            endif
c     sample charm
c     scale up for mesons
            IF(IPAR(76).eq.1) XSCL=XSCL+(1-IABS(IBAR(KIDA)))*PAR(126)
            IF(IFL(j3).eq.3.and.s_rndm(kid).lt.PAR(97)*PAR(125)*XSCL)
     &           IFL(j3) = 4
         ELSE
            IFL(j3) = INT(ONEHALF+S_RNDM(KID))
         ENDIF
         IFL(j4) = -IFL(j3)
         IF(NDEBUG.gt.3) 
     &        WRITE(LUN,*) 'flavor: JW,FLV1,FLV2',J,IFL(j3),IFL(j4)

C...  sample sea pt
 33      SELECT CASE(IPAR(49))
c     in-string pt for sea partons
c     flavor and cm energy dependent avg, exponential dist.
c     avg pt (defined in subroutine ptsetup ):
c     u,d : PAR(46)+PAR(68)*log10(sqs/20.D0)**2
c     s:    PAR(47)+PAR(70)*log10(sqs/20.D0)**2
c     diq:  PAR(48)+PAR(69)*log10(sqs/20.D0)**2
         CASE(1)
            CALL PTDIS_4FLV (IFL(j3),PX(j3),PY(j3))
            PX(j4) = -PX(j3)
            PY(j4) = -PY(j3)
            
c     'primordial' pt
c     c.m. energy dependent avg, exponential
c     same for all flavors
c     avg: PAR(49)+PAR(69)*log10(sqs/20.)**2
         CASE(2)
            CALL PTDIS_4FLV (10,PX(j3),PY(j3))
            PX(j4) = -PX(j3)
            PY(j4) = -PY(j3)

         CASE(3)
c     constant pt
            PX(j3) = EPS5
            PY(j3) = EPS5
            PX(j4) = -PX(j3)
            PY(j4) = -PY(j3)
         CASE(4)
c     sea pt, same as primordial but different params..
c     c.m. energy dependent avg, exponential
c     same for all flavors
c     avg: PAR(132)
            CALL PTDIS_4FLV (30,PX(j3),PY(j3))
            PX(j4) = -PX(j3)
            PY(j4) = -PY(j3)
            
         END SELECT
c     limit parton virtuality         
         PTS = MAX(PX(j3)**2+PY(j3)**2,EPS10)
         IF((XX(j3)**2+XX(J4)**2)/PTS.lt.EIGHT*PAR(122)/S) GOTO 33
         IF(NDEBUG.gt.3) 
     &        WRITE(LUN,*) 'pt: JW,PX,PY,pt',J,Px(j3),Py(j3),sqrt(pts)
      ENDDO     

C...  Prepare the valence partons
 100  XVAL=ONE-XSEAJET
      IF(ndebug.gt.3)
     &     write(lun,*) ' SAMPLE_projectile: val. fraction remaining:',
     &     XVAL

      SELECT CASE(IPAR(47))
      CASE(7)
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
      CASE default
         IF(KRMNT.eq.0.or.IPAR(47).lt.4)THEN
            IF (XVAL.LT.XMINA) goto 20 ! reject sea kinematics
         ENDIF
c     remnant momentum fraction
         IF(KRMNT.ne.0.and.IPAR(53).eq.1)THEN
            IF(S_RNDM(KID).gt.XVAL**(PAR(100)+1)) GOTO 22
         ENDIF
      END SELECT
c     valence quarks are in 1,2 of IFL,XX etc.
      IDXVAL = 3
      IF(KRMNT.ne.0) IDXVAL = 1
      CALL HSPLI (KID,IFL(IDXVAL),IFL(IDXVAL+1))
 110  CHI = CHIDIS(KID,IFL(IDXVAL),IFL(IDXVAL+1))
      XX(IDXVAL) = MAX(CHI*XVAL,XMINA)
      XX(IDXVAL) = MIN(XX(IDXVAL),XVAL-XMINA)
C     FOR MESONS, SPLIT ENERGY SYMETRICALLY.
      IF (ABS(KID).LT.13.AND.S_RNDM(0).LE.HALF) 
     &     XX(IDXVAL)=XVAL-XX(IDXVAL)
      XX(IDXVAL+1)=XVAL-XX(IDXVAL)
      IF(ndebug.gt.3)
     &     write(lun,*) 'SAMPLE_projectile: val. sampled (x1,x2):',
     &     XX(IDXVAL),XX(IDXVAL+1)
c     for baryons force diq distribution
      IF(IBAR(ABS(KID)).ne.0.and.IPAR(47).ne.7)THEN
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
     &     write(lun,*) 'SAMPLE_projectile: val. pt (px,py):',
     &     PX(IDXVAL),PY(IDXVAL)

C...  exchange flavor between central strings and remnant
c     there is one pair of strings for each interaction with another hadron
c     in general allowed for both flavors but diquarks usually strongly suppressed
c     Xchg : prob. of flv exchange between strgs and rmnt
      IF(KRMNT.ne.0)THEN
         do idx=1,2
            iswtd = 0
            i = 1
            XSUPP = ONE
            IF(iabs(ifl(idx)).gt.10)THEN
c     suppress exchange of diq: prob_exchange = prob0 * xsupp
               XSUPP = PAR(83)
            ELSEIF(IPAR(46).eq.2)THEN
c     suppress exchange for fast quark ( usually in mesons )
               IF(xx(idx).gt.xx(3-idx)) XSUPP = PAR(139)
            ENDIF
            DO WHILE (ISWTD.eq.0.and.i.le.Nint)
               if(S_RNDM(KID).lt.XCHG*XSUPP) then 
                  jj = idx+2*i
c     exchange flavor between remnant and sea
                  call iswtch_lmnts(ifl(jj),ifl(idx))
c     also exchange momentum fraction
                  IF(IPAR(46).ne.0) call swtch_lmnts(xx(jj),xx(idx))
c     change flavor id accordingly, i.e. reassamble remnant from new flavor
                  SELECT CASE(IPAR(58))
                  CASE(0)
c     combine to any hadron that matches flavor wise, ignoring (iso)spin
                     call sib_i4flav(ifl(1),ifl(2),idum,idum2,KID1)
                     
                  CASE(1)
c     combine to lightest hadron
                     KID1 = IMRG2HAD(IFL(1),IFL(2))

                  CASE(2,3)
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
                     call sib_i4flav(ifl(1),ifl(2),irnk,idum2,KID1)
                     PAR(5) = PAR5_def
                     PAR(6) = PAR6_def
                     PAR(7) = PAR7_def
                     PAR(143) = PAR143_def

c     reject spin1,isospin singlett
                     IF(KID1.eq.32.and.PAR(111).gt.S_RNDM(KID1))
     &                    KID1 = 27
                  END SELECT
                  iswtd = 1
               endif
               i = i + 1
            ENDDO
         enddo
      ENDIF
      IF(ndebug.gt.3)THEN        
         WRITE(LUN,*)' SAMPLE_projectile: rmnt PID,NTRY: ',KID1,ITRY(3)
         WRITE(LUN,*)' SAMPLE_projectile: output: I,FLV,PX,PY,X,XSUM'
      ENDIF
      XSUM = XJET
      DO j=IDXVAL,2*(Nint+Krmnt)+2*(1-Krmnt)
         XSUM = XSUM + XX(j)
         IF(NDEBUG.gt.3) WRITE(LUN,*) j,IFL(j),PX(J),PY(J),XX(j),XSUM
      ENDDO
      IF(ABS(XSUM-ONE).gt.EPS3) THEN
         WRITE(LUN,*)' SAMPLE_projectile: parton sum incomplete!',
     &        '(ID,XSUM,NCALL):' , KID1,XSUM, NCALL,' aborting..'
         xsum = -1.D0
         xsum = log(xsum)
         call sib_reject
      ENDIF
      IREJ = 0

      END
