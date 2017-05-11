      SUBROUTINE SIBYLL_INI
C-----------------------------------------------------------------------
C...Initialization routine for SYBILL 
C.  
C.  the routine fills the COMMON block /CCSIG/ that contains
C.  important information for the generation of events
C.
C     PARAMETER (NS_max = 20, NH_max = 80)
C     COMMON /S_CCSIG/ SSIG(61,3), PJETC(0:NS_max,0:NH_max,61,2),
C    &    SSIGN(61,3),SSIGNSD(61,3) ALINT(61,3), ASQSMIN, ASQSMAX, DASQS, NSQS
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
C.  SSIGN(J,1) inelastic cross section for p-Air interaction
C.  SSIGN(J,2) inelastic cross section for pi-Air interaction
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
      SAVE

      WRITE(*,100)
 100  FORMAT(' ','====================================================',
     *     /,' ','|                                                  |',
     *     /,' ','|                 S I B Y L L  2.3rc3              |',
     *     /,' ','|                                                  |',
     *     /,' ','|         HADRONIC INTERACTION MONTE CARLO         |',
     *     /,' ','|                        BY                        |',
     *     /,' ','|            Eun-Joo AHN, Felix RIEHN              |',
     *     /,' ','|     R. ENGEL, R.S. FLETCHER, T.K. GAISSER        |',
     *     /,' ','|               P. LIPARI, T. STANEV               |',
     *     /,' ','|                                                  |',
     *     /,' ','| Publication to be cited when using this program: |',
     *     /,' ','| R. Engel et al., Proc. 26th ICRC, 1 (1999) 415   |',
     *     /,' ','|                                                  |',
     *     /,' ','| last tampered with by: F. Riehn (05/12/2015)     |',
     *     /,' ','|                                                  |',
     *     /,' ','|     --> version under heavy development !! <--   |',
     *     /,' ','====================================================',
     *     /)

!       CALL RND_INI
      CALL PAR_INI
      CALL JET_INI
      CALL PDF_INI
      CALL BLOCK_INI
      CALL NUC_GEOM_INI
      CALL SIG_AIR_INI
      CALL DEC_INI
c...  charm frag. normalisation
      call znormal

      END


      SUBROUTINE PAR_INI
C------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      INCLUDE 'sib_cflafr_cmmn.inc'
      INCLUDE 'sib_czdis_cmmn.inc'
      INCLUDE 'sib_czdiss_cmmn.inc'
      INCLUDE 'sib_czdisc_cmmn.inc'
      INCLUDE 'sib_czlead_cmmn.inc'
      INCLUDE 'sib_cpspl_cmmn.inc'
      INCLUDE 'sib_xsctn_cmmn.inc'
      INCLUDE 'sib_xsctn_p_data.inc'
      INCLUDE 'sib_xsctn_pi_data.inc'
      INCLUDE 'sib_cutoff_cmmn.inc'

C..   PAR_INI_23rc1_4_16qmSea1_2qmVal3valB0rMax3_1dMax1iPnlty3_sib23
      PAR(1) = 0.0399999991059
      PAR(2) = 0.3
      PAR(3) = 0.5
      PAR(4) = 0.140000000596
      PAR(5) = 0.300000011921
      PAR(6) = 0.300000011921
      PAR(7) = 0.15000000596
      PAR(8) = 0.5
      PAR(9) = 7.0
      PAR(10) = 1.0
      PAR(11) = 0.065
      PAR(12) = 0.9
      PAR(13) = 0.1
      PAR(14) = 0.08
      PAR(15) = 0.1
      PAR(16) = 0.04
      PAR(17) = 0.0399999991059
      PAR(18) = 0.5
      PAR(19) = 0.8
      PAR(20) = 0.8
      PAR(21) = 0.8
      PAR(22) = 4.0
      PAR(23) = 0.7
      PAR(24) = 0.004
      PAR(25) = 0.004
      PAR(26) = 10.0
      PAR(27) = 0.08
      PAR(28) = 10.0
      PAR(29) = 0.0
      PAR(30) = 2.0
      PAR(31) = 0.33
      PAR(32) = 0.0
      PAR(33) = 0.10000000149
      PAR(34) = 0.0
      PAR(35) = 0.0
      PAR(36) = 0.7
      PAR(37) = 0.0
      PAR(38) = 2.0
      PAR(39) = 1.0
      PAR(40) = 0.0
      PAR(41) = 1.0
      PAR(42) = 3.0
      PAR(43) = 0.2
      PAR(44) = 0.990000009537
      PAR(45) = 1.0
      PAR(46) = 0.18
      PAR(47) = 0.28
      PAR(48) = 0.3
      PAR(49) = 0.1
      PAR(50) = 0.600000023842
      PAR(51) = 0.006
      PAR(52) = 0.006
      PAR(53) = 6.0
      PAR(54) = 0.20000000298
      PAR(55) = 0.0
      PAR(56) = 0.0
      PAR(57) = 0.0
      PAR(58) = 0.0
      PAR(59) = 0.6
      PAR(60) = 0.800000011921
      PAR(61) = 0.66
      PAR(62) = 0.0
      PAR(63) = 1.0
      PAR(64) = 0.25
      PAR(65) = 0.300000011921
      PAR(66) = 0.300000011921
      PAR(67) = 0.6
      PAR(68) = 0.006
      PAR(69) = 0.05
      PAR(70) = 0.007
      PAR(71) = 3.0
      PAR(72) = 0.3
      PAR(73) = 0.5
      PAR(74) = 0.6
      PAR(75) = 0.0
      PAR(76) = 0.2
      PAR(77) = 0.7
      PAR(78) = 0.9
      PAR(79) = 10.0
      PAR(80) = 0.6
      PAR(81) = 10000.0
      PAR(82) = 0.02
      PAR(83) = 0.0
      PAR(84) = 6.0
      PAR(85) = 1.0
      PAR(86) = 1.0
      PAR(87) = 0.25
      PAR(88) = 0.800000011921
      PAR(89) = 0.6
      PAR(90) = 7.0
      PAR(91) = -2.5
      PAR(92) = 0.2
      PAR(93) = 1.0
      PAR(94) = 4.0
      PAR(95) = 0.0
      PAR(96) = 1.0
      PAR(97) = 0.002
      PAR(98) = 2.0
      PAR(99) = 0.330000013113
      PAR(100) = 2.0
      PAR(101) = 1.0
      PAR(102) = 0.0
      PAR(103) = 2.0
      PAR(104) = 0.7
      PAR(105) = 0.1
      PAR(106) = 0.0
      PAR(107) = 0.0
      PAR(108) = 0.0
      PAR(109) = 20.0
      PAR(110) = 1.5
      PAR(111) = 0.0
      PAR(112) = 0.7
      PAR(113) = 0.9
      PAR(114) = 2.0
      PAR(115) = 0.0
      PAR(116) = 1.0
      PAR(117) = 0.0
      PAR(118) = 0.005
      PAR(119) = 0.0
      PAR(120) = 1.0
      PAR(121) = 0.300000011921
      PAR(122) = 0.0
      PAR(123) = 0.300000011921
      PAR(124) = 1.0
      PAR(125) = 1.0
      PAR(126) = 1.0
      PAR(127) = 6.0
      PAR(128) = 1.0
      PAR(129) = 0.08
      PAR(130) = 12.0
      PAR(131) = 0.5
      PAR(132) = 0.5
      PAR(133) = 10.0
      PAR(134) = -5.0
      PAR(135) = 6.0
      PAR(136) = 0.0
      PAR(137) = 0.0
      PAR(138) = 0.0
      PAR(139) = 0.0
      PAR(140) = 0.0
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
      IPAR(16) = 7
      IPAR(17) = 1
      IPAR(18) = 4
      IPAR(19) = 1
      IPAR(20) = 0
      IPAR(21) = 0
      IPAR(22) = 0
      IPAR(23) = 0
      IPAR(24) = 0
      IPAR(25) = 0
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
      IPAR(41) = 1
      IPAR(42) = 3
      IPAR(43) = 1
      IPAR(44) = 0
      IPAR(45) = 0
      IPAR(46) = 1
      IPAR(47) = 6
      IPAR(48) = 1
      IPAR(49) = 4
      IPAR(50) = 10
      IPAR(51) = 2
      IPAR(52) = 0
      IPAR(53) = 1
      IPAR(54) = 0
      IPAR(55) = 0
      IPAR(56) = 0
      IPAR(57) = 1
      IPAR(58) = 2
      IPAR(59) = 1
      IPAR(60) = 0
      IPAR(61) = 100
      IPAR(62) = 1
      IPAR(63) = 0
      IPAR(64) = 0
      IPAR(65) = 1
      IPAR(66) = 1
      IPAR(67) = 3
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

      SUBROUTINE MESON_FLV_MRG_INI
c     change flavor merging for pions (favor spin)
      INCLUDE 'sib_kflv_cmmn.inc'
c     pi+ --> rho+
      KFLV(2,1) = 25
c     pi- --> rho-
      KFLV(1,2) = 26
c     pi0 --> rho0
      KFLV(1,1) = 27
      KFLV(2,2) = 27     
      END


      BLOCK DATA PARAM_INI
C-----------------------------------------------------------------------
C....This block data contains default values
C.   of the parameters used in fragmentation
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      INCLUDE 'sib_cutoff_cmmn.inc'
      INCLUDE 'sib_czdis_cmmn.inc'
      INCLUDE 'sib_czdiss_cmmn.inc'
      INCLUDE 'sib_czdisc_cmmn.inc'
      INCLUDE 'sib_czlead_cmmn.inc'
      INCLUDE 'sib_cpspl_cmmn.inc'
      INCLUDE 'sib_difmass_cmmn.inc'

      INCLUDE 'sib_cnt_cmmn.inc'
      INCLUDE 'sib_cnt_data.inc'

      INCLUDE 'sib_utl_cmmn.inc'
      INCLUDE 'sib_utl_data.inc'
      
      COMMON /CKFRAG/ KODFRAG

C...  default output unit
      DATA LUN /6/

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
C...Splitting parameters
      DATA CCHIK /15*0.D0,21*2.D0,6*3.D0,57*0.D0,18*3.D0/
C...Parameters of flavor formation 
c     last in use: 135
      DATA PAR /0.04,0.3,0.3,0.14,0.3,0.3,0.15,0.,7., ! 10
     &     0.,0.,0.2,0.,0.04,0.04,0.04,0.04, 0.5,0.8, 0.5, ! 20
     &     0.8,6.,0.5, 0.004,5*0.,0.7,                ! 30
     &     2*0.,0.1,0.,3.,0.35,0.,0.5,2*0.,           ! 40
     &     1.,2.,0.,0.99,0.,0.3,0.45,0.6,0.6,0.6,     ! 50
     &     .03,.03,6.,0.2,4*0.,1.1,0.8,               ! 60
     &     0.33,3.,1.,0.25,0.3,0.3,0.6,.007,.03,.007, ! 70
     &     1.,0.3,0.,0.3,0.0,0.2,0.5,1.0,10.,0.,      ! 80
     &     1000.,1000.,1.,6.,1.,0.,0.3,0.8,0.3,31.,   ! 90
     &     1.,6.5,1.,1.,0.,1.0,0.004,1.,0.33,1.,      ! 100
     &     1.,0.,2.,0.3,0.15,3*0.,20.,0.25,           ! 110
     &     0.,0.7,0.3,0.,0.,1.,3*0.,1.,               ! 120
     &     0.3,0.,0.3,1.,1.,1.,6.,1.,1.,6.,           ! 130
     &     0.0001,0.5,31.10362D0,-15.29012D0,6.5D0,   ! 135
     &     5 *0./                                     ! 140
c     last in use:80
      DATA IPAR /9*0,1,0,0,8*0,20*0,9*0,10,2,9*0,100,19*0/

C...Fragmentation of nuclei
      DATA KODFRAG /0/
C...Debug label and event counter
      DATA Ndebug /0/
      DATA Ncall /0/
C...Diffractive mass parameters
      DATA XM2MIN /1.5D0, 0.2D0, 0.6D0/                  ! M_x**2(min) GeV**2
      DATA ALXMIN /0.405465D0,-1.6094379D0,-0.5108256D0/ ! log[M_x**2(min)]
      DATA SLOP0 /6.5D0/                 ! b (slope_ for Mx**2 > 5 GeV**2
      DATA ASLOP /31.10362D0/            ! fit to the slope parameter.
      DATA BSLOP /-15.29012D0/

      END


      SUBROUTINE PARAM_PRINT(LUN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      COMMON /S_CZLEAD/ CLEAD, FLEAD
      COMMON /S_CPSPL/ CCHIK(3,39)
      COMMON /S_DEBUG/ Ncall, Ndebug, Lunn
      include 'sib_nw_prm.inc'
      include 'sib_run_cmmn.inc'
      include 'sib_cutoff_cmmn.inc'
      INCLUDE 'sib_cflafr_cmmn.inc'
      INCLUDE 'sib_czdis_cmmn.inc'
      INCLUDE 'sib_cqdis2_cmmn.inc'
      INCLUDE 'sib_cqdis_cmmn.inc'

      WRITE (LUN, 25)
25      FORMAT( //,1x,40('-'), /
     +   ' SIBYLL MONTE CARLO PROGRAM. Version 2.2f',/,1x,40('-'),/
     +   ' List of parameters: ' )

      WRITE (LUN, 31) FAin, FB0in
31      FORMAT (' Parameters of longitudinal fragmentation: ', /,
     +          '  f(z) = (1-z)**a * exp(-b * mt**2/z) ', /,
     +          '  a = ', f9.3, 3x, ' b = ', f9.3, ' GeV**-2' )
      WRITE (LUN, 32) CLEAD, 1.D0/FLEAD-1.D0
32      FORMAT (' Parameters of leading fragmentation: ', /,
     +   '  f(z) = c + (1-z)**a ', /,
     +   '  c = ',f9.3,3x,' a = ',f9.3) 

        WRITE (LUN, 33) str_mass_val, str_mass_sea
 33     FORMAT (' Mass cuts ', /,
     +          '  val = ', f9.3, 3x, ' sea = ', f9.3, ' GeV' )

      WRITE (LUN, 35) PPT02(1), PPT02(3), PPT02(11),ppt02(10),ppt02(20)
35      FORMAT (' <pT> of sea partons ', /,
     +   2x,'<pT>(u/d) ',F8.3,2x,'<pT>(s) ',f8.3,2x,'<pT>(qq) ',f8.3,
     +     2x,'<pT>(val) ',f8.3,2x,'<pT>(sea) ',f8.3)

      WRITE (LUN, 120) (PAR(K),K=1,24)
120      FORMAT (1x, 'Parameters of flavor formation: ',/,
     +   3x,'PAR(1) = Prob(qq)/Prob(q) =              ',F10.2,/,
     +   3x,'PAR(2) = Prob(s)/Prob(u)  =              ',F10.2,/,
     +   3x,'PAR(3) = Prob(us)/Prob(ud) =             ',F10.2,/,
     +   3x,'PAR(4) = Prob(ud_0)/Prob(ud_1) =         ',F10.2,/,
     +   3x,'PAR(5) = Prob(Vector)/Prob(Scalar) =     ',F10.2,/,
     +   3x,'PAR(6) = Prob(K*)/Prob(K) =              ',F10.2,/,
     +   3x,'PAR(7) = Prob(spin 3/2)/Prob(spin=1/2) = ',F10.2,/,
     +   3x,'PAR(8) = Prob(B-M-Bbar)/Prob(B-Bbar) =   ',F10.2,/,
     +   3x,'PAR(9) = Phase space suppression of MI = ',F10.2,/,
     +   3x,'PAR(10)= Low-energy limit for pt cutoff= ',F10.2,/,
     +   3x,'PAR(11)= Pt cutoff factor for exp      = ',F10.2,/,
     +   3x,'PAR(12)= Pt cutoff factor for exp      = ',F10.2,/,
     +   3x,'PAR(13)= max. mass in diffraction      = ',F10.2,/,
     +   3x,'PAR(14)= Prob(qq)/Prob(q) std. value   = ',F10.2,/,
     +   3x,'PAR(15)= Prob(qq)/Prob(q) in hard jets = ',F10.2,/,
     +   3x,'PAR(16)= Prob(qq)/Prob(q) in diff.     = ',F10.2,/,
     +   3x,'PAR(17)= not used = ',F10.2,/,
     +   3x,'PAR(18)= not used = ',F10.2,/,
     +   3x,'PAR(19)= not used = ',F10.2,/,
     +   3x,'PAR(20)= not used = ',F10.2,/,
     +   3x,'PAR(21)= not used = ',F10.2,/,
     +   3x,'PAR(22)= effective scale in PDF (Q2) = ',F10.2,/,
     +   3x,'PAR(23)= not used = ',F10.2,/,
     +   3x,'PAR(24)= Prob(s->c) = ',F10.2)

      WRITE (LUN, 130) (IPAR(K),K=1,16)
130      FORMAT (1x, 'Model switches: ',/,
     +   3x,'IPAR(1) = not used =                     ',I4,/,
     +   3x,'IPAR(2) = not used =                     ',I4,/,
     +   3x,'IPAR(3) = exponential pt =               ',I4,/,
     +   3x,'IPAR(4) = decouple qq/q in val. strings= ',I4,/,
     +   3x,'IPAR(5) = decouple qq/q in hm. diff. =   ',I4,/,
     +   3x,'IPAR(6) = decouple qq/q in hard strings= ',I4,/,
     +   3x,'IPAR(7) = remnant (not implemented yet)= ',I4,/,
     +   3x,'IPAR(8) = jet kinematic pdf set (DO/GRV)=',I4,/,
     +   3x,'IPAR(9) = smear lowest diff. mass =      ',I4,/,
     +   3x,'IPAR(10)= high mass diff. mode (d:ON)=   ',I4,/,
     +   3x,'IPAR(11)= leading vec. meson prod. model=',I4,/,
     +   3x,'IPAR(12)= inel. screening in pAir=       ',I4,/,
     +   3x,'IPAR(13)= decouple qq/q in val. strings= ',I4,/,
     +   3x,'IPAR(14)= fireball model =               ',I4,/,
     +   3x,'IPAR(15) = charm production =            ',I4,/,
     +   3x,'IPAR(16) = charmed transverse momentum = ',I4,/,
     +   3x,'IPAR(17)= full charm model =             ',I4)

      WRITE (LUN, 40)
      WRITE (LUN, 41) CCHIK (1,13), CCHIK(2,13)
40      FORMAT(' Parameters of hadron splitting ' )
41      FORMAT('   p -> [(ud) u] splitting: alpha = ', F10.3, /,
     +         '   p -> [(uu) d] splitting: alpha = ', F10.3 )

      END


      SUBROUTINE SIB_LIST(LUN)
C-----------------------------------------------------------------------
C...This routine prints the event record for the
C.  current event on unit LUN
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      COMMON /S_DEBUG/ Ncall, Ndebug, Lunn
      include 'sib_nw_prm.inc'
      INCLUDE 'sib_run_cmmn.inc'
      INCLUDE 'sib_plist_cmmn.inc'
      INCLUDE 'sib_parto_cmmn.inc'
      INCLUDE 'sib_plist1_cmmn.inc'
      INCLUDE 'sib_chist_cmmn.inc'
      INCLUDE 'sib_chp_cmmn.inc'
      INCLUDE 'sib_cnam_cmmn.inc'
      INCLUDE 'sib_rmnt_cmmn.inc'

      CHARACTER*7 CTGT(0:20)
      DATA CTGT /'Air    ','Proton ',19*'Nucleus'/

      CHARACTER CODE*18
      CHARACTER*18 NAMDIF(0:3)
      DATA NAMDIF /'Non-diff. event   ',
     &  'Beam diffraction  ','Target diffraction','Double diffraction'/
      CHARACTER*18 NAMRMNT(0:3)
      DATA NAMRMNT /'No resolvd remnant',
     &  'Beam remnant     ','Target remnant    ','Double remnant    '/

 50   FORMAT(3X,88('-'),/,25X,'SIBYLL EVENT SUMMARY',25X,
     &     /,3X,88('-'))

      WRITE (LUN,50)
 52   FORMAT( 3X,'Beam + Target @ Energy:',2X,A6,2X,'+',2X,A7,2X,
     &     '@',2X,E9.3,1X,'GeV')

      WRITE (LUN,52)
     &     NAMP(IABS(KB)),CTGT(IAT),SQS
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
      ichar = 0
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
          ichar = ichar+sign(1,l)*ICHP(iabs(l))
          ibary = ibary+sign(1,l)*IBAR(iabs(l))
          ichmd = ichmd+sign(1,l)*ICHM(iabs(l))
          istrg = istrg+sign(1,l)*ISTR(iabs(l))
        endif
      ENDDO
      CALL PFsum(1,NP,Esum,PXsum,PYsum,PZsum,NF)
      WRITE(LUN,140) PXsum,PYsum,PZsum,Esum
100      FORMAT(3X,'N  Particle',12X,'Int',2x,'Jdif',2x,'Prnt',2x,'Proc'
     +         ,6x,'PX',9x,'PY',9x,'PZ',9x,'E', /, 3X,88('-'))
120      FORMAT(1X,I4,1X,A18,1X,I4,1X,I4,1X,I4,3X,I5,1X,
     +        2(F9.3,2X),1X,2(E9.3,2X))
140      FORMAT(3X,88('-'),/,1X,'Tot = ',29X,2(F9.3,2X),G9.3,2X,E9.3)
         write(LUN,'(1x,a,i3,3x,a,i3)') 'Total charge:',ichar,
     &        'baryon number:',ibary
         write(LUN,'(1x,a,i3,3x,a,i3)') 'Total strangeness:',istrg,
     &        'charm number:',ichmd

      END



      SUBROUTINE KCODE (J,CODE,NC)
C...Produce the code for parton J
C.  Input K, Output CODE, NC=number of characters
C..................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      CHARACTER*5 CODE
      CHARACTER*1 NAMQ(4)
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

      SUBROUTINE SIB_PARTPR(LUN)
C----------------------------------------------------------------
C     prints the particles known to SIBYLL with their internal
C     and PDG labels \FR'13
C----------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'sib_mass_cmmn.inc'
      INCLUDE 'sib_cnam_cmmn.inc'

      WRITE(LUN,50)
 50   FORMAT(2X,'SIBYLL PARTICLE TABLE:',/,2x,56('-'))
      WRITE(LUN,100)
 100  FORMAT(2X,'Particle',4X,'SIB PID',6x,'SIB2PDG',6x,'SIB2PDG^-1', 
     &     4x,'MASS'/, 2X,56('-'))

      DO J=1,99
         IA = ISIB_PID2PDG( j )         
         IF(IA.ne.0)
     &        WRITE (LUN,120)  NAMP(J), J, IA, ISIB_PDG2PID( IA ), AM(J)
      ENDDO
 120  FORMAT(4X,A6,4X,I4,8X,I6,8X,I4,5X,F9.3)

      END


      INTEGER FUNCTION ISIB_PID2PDG(Npid)
C----------------------------------------------------------------
C     conversion of SIBYLL internal particle code to PDG standard
C
C     input:     Npid        internal particle number
C     output:    sib_pid2pdg  PDG particle number
C
C     based on similar phojet function \FR'13
C----------------------------------------------------------------
      SAVE
      COMMON /S_PDG2PID/ ID_PDG_LIST(99),ID_LIST(577)
      INTEGER NPIDA,NPID

      Npida = iabs(Npid)
      isib_pid2pdg = ID_PDG_LIST(Npida)
      IF(NPID.lt.0)isib_pid2pdg = isign(isib_pid2pdg,Npid)
      RETURN
      END



      INTEGER FUNCTION ISIB_PDG2PID(Npdg)
C----------------------------------------------------------------
C     conversion of PDG standard particle code to SIBYLL internal
C
C     input:     Npdg        PDG particle number
C     output:    sib_pdg2pid internal particle id
C
C     based on similar phojet function \FR'13
C----------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      COMMON /S_PDG2PID/ IPID_PDG_LIST(99),ID_LIST(577)
      INCLUDE 'sib_debug_cmmn.inc'
      INCLUDE 'sib_csydec_cmmn.inc'

      Nin = abs(Npdg)
      if((Nin.gt.99999).or.(Nin.eq.0)) then
C  invalid particle number
        if(ndebug.gt.5) write(6,'(1x,A,I10)')
     &    'isib_pdg2pid: invalid PDG ID number ',Npdg
        isib_pdg2pid = 0
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
     &    'isib_pdg2pid: particle not in table ',Npdg
        isib_pdg2pid = 0
        return
      endif
      ID_out = ID_list(Nout)
      IF(abs(ID_out).gt.99)then
         isib_pdg2pid = 0
         return
      else

         if(IPID_pdg_list(ID_list(Nout)).eq.Nin) then
C     particle ID found
            isib_pdg2pid = ID_list(Nout)
            if (NPDG.lt.0) isib_pdg2pid = lbarp( isib_pdg2pid )
            return
         else
C     increment and try again
            Nout = Nout + 5
            If(Nout.gt.577) Nout = Mod(Nout,577)
            goto 100
         endif
      endif
      END



      SUBROUTINE PDG_INI
C----------------------------------------------------------------
C     PDG conversion blocks \FR'13
C----------------------------------------------------------------
      SAVE
      INCLUDE 'sib_debug_cmmn.inc'
      PARAMETER ( ID_PDG_MAX = 260 )
      COMMON /S_PDG2PID/ ID_PDG_LIST(99),ID_LIST(577)
      DATA ID_PDG_LIST /22,-11,11,-13,13,111,211,-211,321,-321, !10
     &     130,310,2212,2112,12,-12,14,-14,-2212,-2112,         !20
     &     311,-311,221,331,213,-213,113,10321,-10321,10311,    !30
     &     -10311,223,333,3222,3212,3112,3322,3312,3122,2224,   !40
     &     2214,2114,1114,3224,3214,3114,3324,3314,3334,0,      !50
     &     202212,202112,212212,212112,4*0,411,-411,            !60
     &     900111,900211,-900211,7*0,                           !70
     &     421,-421,441,431,-431,433,-433,413,-413,423,         !80
     &     -423,0,443,4222,4212,4112,4232,4132,4122,0,          !90
     &     0,0,0,4224,4214,4114,4324,4314,4332/

      IF(Ndebug.gt.2)
     & WRITE(6,*) 'INITIALIZING PDG TABLES..'
      CALL pho_cpcini(ID_pdg_max,ID_pdg_list,ID_list)
      
      END


      SUBROUTINE pho_cpcini(Nrows,Number,List)
C***********************************************************************
C
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
      SAVE
      integer Number(*),List(*),Nrows
      Integer Nin,Nout,Ip,I

      do I = 1,577
        List(I) = 0
      enddo

C    Loop over all of the elements in the Number vector

        Do 500 Ip = 1,Nrows
            Nin = Number(Ip)

C    Calculate a list number for this particle id number
            If(Nin.Gt.99999.or.Nin.Le.0) Then
                 Nout = -1
            Else If(Nin.Le.577) Then
                 Nout = Nin
            Else
                 Nout = Mod(Nin,577)
            End If

 200        continue

            If(Nout.Lt.0) Then
C    Count the bad entries
c                Write(6,'(1x,a,i10)')
c     &            'pho_cpcini: invalid particle ID',Nin
                Go to 500
            End If
            If(List(Nout).eq.0) Then
                List(Nout) = Ip
            Else
                If(Nin.eq.Number(List(Nout))) Then
c                  Write(6,'(1x,a,i10)')
c     &              'pho_cpcini: double particle ID',Nin
                End If
                Nout = Nout + 5
                If(Nout.Gt.577) Nout = Mod(Nout, 577)

                Go to 200
            End If
 500      Continue

      END

      SUBROUTINE PFsum(N1,N2,ETOT,PXT,PYT,PZT,NF)
C...Return the energy,px,py,pz and the number of stable
C.  particles in the list between N1 and N2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
c      COMMON /S_PLIST/ P(8000,5), LLIST(8000), NP
      INCLUDE 'sib_plist_cmmn.inc'
      DATA ZERO /0.D0/
      NF=0
      ETOT=ZERO
      PXT=ZERO
      PYT=ZERO
      PZT=ZERO
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


      SUBROUTINE QNUM (JQ,JS,JC,JB,JBA, NC, NF)
C...Return the quantum numbers of one event
C.  JQ = charge, JB = baryon number, JS = strangeness, JC = charmedness
C.  JBA = (number of baryons+antibaryons)
C.  NC  = number of charged particles
C.  NF  = number of final particles
C..................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      INCLUDE 'sib_plist_cmmn.inc'
      INCLUDE 'sib_chp_cmmn.inc'
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


      BLOCK DATA KFLV_INI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      INCLUDE 'sib_kflv_cmmn.inc'
      DATA (KFLV(1,i),i=1,4) /6,8,10,71/
      DATA (KFLV(1,i),i=5,34) /6*0,40,13,34,84,6*0,13,14,39,89,6*0,
     &     34,39,37,87/
      DATA (KFLV(2,i),i=1,4) /7,6,21,59/
      DATA (KFLV(2,i),i=5,34) /6*0,13,14,39,89,6*0,14,43,36,86,6*0,
     &     39,36,38,88/
      DATA (KFLV(3,i),i=1,4) /9,22,33,74/
      DATA (KFLV(3,i),i=5,34) /6*0,34,39,35,87,6*0,39,36,38,88,6*0,
     &     35,36,49,99/
      DATA (KFLV(4,i),i=1,4) /72,60,75,83/
      DATA (KFLV(4,i),i=5,34) /6*0,84,85,87,0,6*0,85,86,88,0,6*0,
     &     87,88,99,0/

      END
