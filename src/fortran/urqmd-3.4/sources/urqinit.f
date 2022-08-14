*-- Author :    J. Wentz  IK FZK KARLSRUHE     10/05/1999
C=======================================================================

      SUBROUTINE URQINI(MDEBUG, iudebug)

C-----------------------------------------------------------------------
C  URQ(MD) INI(TIALIZATION)
C**anfe Modified to accept the LUN directly via MDEBUG
C
C  FIRST INITIALIZATION
C  THIS SUBROUTINE IS CALLED FROM START.
C-----------------------------------------------------------------------

      IMPLICIT NONE
      include 'boxinc.f'
      include 'coms.f'
      include 'inputs.f'
      include 'options.f'

      INTEGER          i,io,ia,ie,id, MDEBUG, iudebug
      CHARACTER        CTPStrg(numctp)*60, CTOStrng(numcto)*60
      CHARACTER        debug_output*60
C-----------------------------------------------------------------------

c c Open debug output unit
c       if (debug_output.eq.'') then
c           MDEBUG = 6
c       else
c           MDEBUG = 42
c           open(unit=MDEBUG, file = debug_output, status = 'UNKNOWN')
c       endif


C  SET NMAX TO DEFAULT VALUE
      call set0
      call params

C  THIS IS THE SUBSTITUE FOR THE URQMD INPUT ROUTINE
C  INITIALIZE COUNTERS
      boxflag = 0
      mbflag  = 0
      edens   = 0.d0
      para    = 0
      solid   = 0
      mbox    = 0
      io      = 0

C  THE FOLLOWING FLAGS CHECK, WHETHER ALL NECESSARY INPUT IS GIVEN
C  PROJECTILE
      prspflg = 0
C  TARGET
      trspflg = 0
C
      srtflag = 0
      firstev = 0
C  EXCITATION FUNCTION
      nsrt    = 1
      npb     = 1
      efuncflag = 0
C  DEFAULT NUMBER OF EVENTS
      nevents = 1
C  DEFAULT NUMBER OF TIMESTEPS
      nsteps  = 1000

C  SKIP CONDITIONS ON UNIT 13, 14, 15, 16 & 18
C  SUPPRESS ALL OUTPUT
      bf13 = .true.
      bf14 = .true.
      bf15 = .true.
      bf16 = .true.
      bf18 = .true.
      bf19 = .true.
      bf20 = .true.
C  SET DEBUG OUTPUT DEPENDING ON CHOSEN DEBUG LEVEL
C  SET THE OUTPUT OF UNITS 13, 14, 15 TO THE DEBUG OUTPUT UNIT
      IF     ( IUDEBUG .EQ. 1 ) THEN
        bf13 = .true.
        bf14 = .false.
        call uounit(14,MDEBUG)
        bf15 = .true.
      ELSEIF ( IUDEBUG .EQ. 2 ) THEN
        bf13 = .false.
        call uounit(13,MDEBUG)
        bf14 = .true.
        bf15 = .true.
      ELSEIF ( IUDEBUG .GT. 2 ) THEN
        bf13 = .true.
        bf14 = .true.
        bf15 = .false.
        call uounit(15,MDEBUG)
      ENDIF
      do  i = 1, numcto
         CTOdc(i) = '  '
      enddo
      do  i = 1, numctp
         CTPdc(i) = '  '
      enddo

c default settings for CTParam and CTOption cccccccccccccccccccccccccccccc
      CTParam(1)=1.d0
      CTPStrg(1)='scaling factor for decay-width'
      CTParam(2)=0.52d0
      CTPStrg(2)='used for minimal stringmass & el/inel cut in makestr'
      CTParam(3)=2d0
      CTPStrg(3)='velocity exponent for modified AQM'
      CTParam(4)=0.3d0
      CTPStrg(4)='transverse pion mass, used in make22 & strexct'
      CTParam(5)=0d0
      CTPStrg(5)='probabil. for quark rearrangement in cluster'
      CTParam(6)=0.37d0
      CTPstrg(6)='strangeness probability'
      CTParam(7)=0.d0
      CTPStrg(7)='charm probability (not yet implemented in UQMD)'
      CTParam(8)=0.093d0
      CTPStrg(8)='probability to create a diquark'
      CTParam(9)=0.35d0
      CTPStrg(9)='kinetic energy cut off for last string break'
      CTParam(10)=0.25d0
      CTPStrg(10)='min. kinetic energy for hadron in string'
      CTParam(11)=0.d0
      CTPStrg(11)='fraction of non groundstate resonances'
      CTParam(12)=.5d0
      CTPStrg(12)='probability for rho 770 in String'
      CTParam(13)=.27d0
      CTPStrg(13)='probability for rho 1450 (rest->rho1700)'
      CTParam(14)=.49d0
      CTPStrg(14)='probability for omega 782'
      CTParam(15)=.27d0
      CTPStrg(15)='probability for omega 1420(rest->om1600)'
      CTParam(16)=1.0d0
      CTPStrg(16)='mass cut betw. rho770 and rho 1450'
      CTParam(17)=1.6d0
      CTPSTRG(17)='mass cut betw. rho1450 and rho1700'
      CTParam(18)=.85d0
      CTPStrg(18)='mass cut betw. om 782 and om1420'
      CTParam(19)=1.55d0
      CTPStrg(19)='mass cut betw. om1420 and om1600'
      CTParam(20)=0.0d0
      CTPStrg(20)=' distance for second projectile'
      CTParam(21)=0.0d0
      CTPStrg(21)=' deformation parameter'
      CTParam(25)=.9d0
      CTPStrg(25)=' probability for diquark not to break'
      CTParam(26)=50d0
      CTPStrg(26)=' maximum trials to get string masses'
      CTParam(27)=1d0
      CTPStrg(27)=' scaling factor for xmin in string excitation'
      CTParam(28)=1d0
      CTPStrg(28)=' scaling factor for transverse fermi motion'
      CTParam(29)=1d0
      CTPStrg(29)=' double strange di-quark suppression factor '
      CTParam(30)=1.5
      CTPStrg(30)=' radius offset for initialisation  '
      CTParam(31)=1.6d0
      CTPStrg(31)=' sigma of gaussian for tranverse momentum tranfer '
      CTParam(32)=0d0
      CTPStrg(32)=' alpha-1 for valence quark distribution  '
      CTParam(33)=2.5d0
      CTPStrg(33)=' betav for valence quark distribution  (DPM)'
      CTParam(34)=0.1
      CTPStrg(34)=' minimal x multiplied with ecm  '
      CTParam(35)=3.0
      CTPStrg(35)=' offset for cut for the FSM '
      CTParam(36)=0.275d0
      CTPStrg(36)=' fragmentation function parameter a  '
      CTParam(37)=0.42d0
      CTPStrg(37)=' fragmentation function parameter b  '
      CTParam(38)=1.08d0
      CTPStrg(38)=' diquark pt scaling factor '
      CTParam(39)=0.8d0
      CTPStrg(39)=' strange quark pt scaling factor '
      CTParam(40)=0.5d0
      CTPStrg(40)=' betas-1 for valence quark distribution (LEM)'
      CTParam(41)=0.0
      CTPStrg(41)=' distance of initialisation'
      CTParam(42)=0.55d0
      CTPStrg(42)=' width of gaussian -> pt in string-fragmentation '
      CTParam(43)=5.d0
      CTPStrg(43)=' maximum kinetic energy in mesonic clustr '
      CTParam(44)=.8d0
      CTPStrg(44)=' prob. of double vs. single excitation for AQM inel.'
      CTParam(45)=0.5
      CTPStrg(45)=' offset for minimal mass generation of strings'
      CTParam(46)=800000
      CTPStrg(46)=' maximal number of rejections for initialisation'
      CTParam(47)=1.0
      CTPStrg(47)=' field feynman fragmentation funct. param. a'
      CTParam(48)=2.0
      CTPStrg(48)=' field feynman fragmentation funct. param. b'
      CTParam(49)=0.5
      CTPStrg(49)='additional single strange diquark suppression factor'
      CTParam(50)=1d0
      CTPStrg(50)=' enhancement factor for 0- mesons'
      CTParam(51)=1d0
      CTPStrg(51)=' enhancement factor for 1- mesons'
      CTParam(52)=1d0
      CTPStrg(52)=' enhancement factor for 0+ mesons'
      CTParam(53)=1d0
      CTPStrg(53)=' enhancement factor for 1+ mesons'
      CTParam(54)=1d0
      CTPStrg(54)=' enhancement factor for 2+ mesons'
      CTParam(55)=1d0
      CTPStrg(55)=' enhancement factor for 1+-mesons'
      CTParam(56)=1d0
      CTPStrg(56)=' enhancement factor for 1-*mesons'
      CTParam(57)=1d0
      CTPStrg(57)=' enhancement factor for 1-*mesons'
      CTParam(58)=1.d0
      CTPStrg(58)=' scaling factor for DP time-delay'
      CTParam(59)=0.7d0
      CTPStrg(59)='scaling factor for leading hadron x-section (PYTHIA)'
      CTParam(60)=3.0d0
      CTPStrg(60)=' resonance/string transition energy for s-chanel'
      CTParam(61)=0.2d0
      CTPStrg(61)=' cell size for hydro grid in fm/c'
      CTParam(62)=200
      CTPStrg(62)=' total hydro grid size, number of cells'
      CTParam(63)=1.d0
      CTPStrg(63)=' minimal hydro start time'
      CTParam(64)=5.d0
      CTPStrg(64)=' factor for freezeout criterium (x*e0)'
      CTParam(65)=1.d0
      CtPStrg(65)=' factor for variation of thydro_start'
      CTParam(66)=1.d10
      CTPStrg(66)=' Rapidity cut for initial state set to'
      CTParam(67)=1.d0
      CTPStrg(67)=' Number of testparticles per real particle'
      CTParam(68)=1.d0
      CTPStrg(68)=' Width of 3d-Gauss for hydro initial state mapping'
      CTParam(69)=0.0d0
      CTPStrg(69)=' Quark density cut for initial state,units  1/rho0/3'
      CTParam(70)=1.0d10
      CTPStrg(70)='Cut in Paseudorapidity-range for the Core density'
      CTParam(71)=2.d0
      CTPStrg(71)='Hypersurface is determined avery nth timestep'
      CTParam(72)=55d-2
      CTPStrg(72)="Ratio Sigma0/(Sigma0+Lambda0) in s-exchange reaction"

cbb Note: If you add more CTParams, please make sure that all parameters
c   are included in the standard event header output in output.f.
c   Currently, 72 CTPs are written.
cc
      CTOption(1)=1
      CTOStrng(1)=' resonance widths are mass dependent '
      CTOption(2)=0
      CTOStrng(2)=' conservation of scattering plane'
      CTOption(3)=0
      CTOStrng(3)=' use modified detailed balance'
      CTOption(4)=0
      CTOStrng(4)=' no initial conf. output '
      CTOption(5)=0
      CTOStrng(5)=' fixed impact parameter'
      CTOption(6)=0
      CTOStrng(6)=' no first collisions inside proj/target'
      CTOption(7)=1
      CTOStrng(7)=' elastic cross section enabled (<>0:total=inelast)'
      CTOption(8)=0
      CTOStrng(8)=' extrapolate branching ratios '
      CToption(9)=0
      CTOStrng(9)=' use tabulated pp cross sections '
      CTOption(10)=0
      CTOStrng(10)=' enable Pauli Blocker'
      CTOption(11)=0
      CTOStrng(11)=' mass reduction for cascade initialisation'
      CTOption(12)=0
      CTOStrng(12)=' string condition =0 (.ne.0 no strings)'
      CTOption(13)=0
      CTOStrng(13)=' enhanced file16 output '
      CTOption(14)=0
      CTOStrng(14)=' cos(the) is distributet between -1..1 '
      CTOption(15)=0
      CTOStrng(15)=' allow mm&mb-scattering'
      CTOption(16)=0
      CTOStrng(16)=' propagate without collisions'
      CTOption(17)=0
      CTOStrng(17)=' colload after every timestep '
      CTOption(18)=0
      CTOStrng(18)=' final decay of unstable particles'
      CTOption(19)=0
      CTOStrng(19)=' allow bbar annihilaion'
      CTOption(20)=0
      CTOStrng(20)=' dont generate e+e- instead of bbar'
      CTOption(21)=0
      CTOStrng(21)=' use field feynman frgm. function'
      CTOption(22)=1
      CTOStrng(22)=' use lund excitation function'
      CTOption(23)=0
      CTOStrng(23)=' lorentz contraction of projectile & targed'
      CTOption(24)=2      ! 1 is default    2 means fast method
      CTOStrng(24)=' Wood-Saxon initialization'
      CTOption(25)=0
      CTOStrng(25)=' phase space corrections for resonance mass'
      CTOption(26)=0
      CTOStrng(26)=' use z -> 1-z for diquark-pairs'
      CTOption(27)=1
      CTOStrng(27)=' reference frame (1=target, 2=projectile, else=cms)'
      CTOption(28)=0 ! A.F. 2020-04-01 This is different from CORSIKA8
      CTOStrng(28)=' propagate spectators also '
      CTOption(29)=2
      CTOStrng(29)=' no transverse momentum in clustr '
      CTOption(30)=0
      CTOStrng(30)=' frozen fermi motion '
      CTOption(31)=0
      CTOStrng(31)=' reduced mass spectrum in string'
      CTOption(32)=0
      CTOStrng(32)=' masses are distributed acc. to m-dep. widths'
      CTOption(33)=0
      CTOStrng(33)=' use tables & m-dep. for pmean in fprwdt & fwidth'
      CTOption(34)=1
      CTOStrng(34)=' lifetme according to m-dep. width'
      CTOption(35)=1
      CTOStrng(35)=' generate high precision tables'
      CTOption(36)=0
      CTOStrng(36)=' normalize Breit-Wigners with m.dep. widths '
      CTOption(37)=0
      CTOStrng(37)=' heavy quarks form di-quark clusters'
      CTOption(38)=0
      CTOStrng(38)=' scale p-pbar to b-bbar with equal p_lab '
      CTOption(39)=0
      CTOStrng(39)=' dont call pauliblocker'
      CTOption(40)=0
      CTOStrng(40)=' read old fort.14 file '
      CTOption(41)=0
      CTOStrng(41)=' generate extended output for cto40'
      CTOption(42)=0
      CTOStrng(42)=' hadrons now have color fluctuations'
      CTOption(43)=0
      CTOStrng(43)=" don't generate dimuon intead of dielectron output"
      CTOption(44)=1
      CTOStrng(44)=' call PYTHIA for hard scatterings'
      CTOption(45)=0
      CTOStrng(45)=' hydro mode'
      CTOption(46)=0
      CTOStrng(46)=' calculate quark density instead of baryon density'
      CTOption(47)=5
      CTOStrng(47)=' flag for equation of state for hydro'
      CTOption(48)=0
      CTOStrng(48)=' propagate only N timesteps of hydro evolution'
      CTOption(49)=0
      CTOStrng(49)=' propagate also spectators with hydrodynamics'
      CTOption(50)=0
      CTOStrng(50)=' (additional) f14/f19 output after hydro phase'
      CTOption(52)=0
      CTOStrng(52)=' Freezeout procedure changed'
      CTOption(53)=0
      CTOStrng(53)=' efficient momentum generation in Cooperfrye'
      CTOption(54)=0
      CTOStrng(54)=' OSCAR-Output during hydro evolution'
      CTOPtion(55)=0
      CTOStrng(55)=' f19 output adjusted for visualization'
      CTOPtion(56)=0
      CTOStrng(56)=' f15 output has unique particle id'
      CTOPtion(57)=1
      CTOStrng(57)=' legacy event header w/ missing cto and ctp'
      CTOPtion(58)=0
      CTOStrng(58)=' standard event header in collision file (file15)'
      CTOption(59)=1
      CTOStrng(59)=' activate Baryon-Baryon strangeness exchange'

C  INITIALIZE ARRAYS FOR SPECIAL PRO/TAR COMBINATIONS
      do  i = 1, 2
         spityp(i) = 0
         spiso3(i) = 0
      enddo

C  INITIALIZE ARRAYS FOR SPECIAL PARTICLES
      EoS = 0

C  INITIALIZES SOME ARRAYS
      call strini      ! initialize mixing angles for meson-multipletts
      call loginit

      IF ( CTOption(33) .EQ. 0  .OR.  CTOption(9) .EQ. 0 ) THEN
        call loadwtab(io)
      ENDIF

C  CALCULATE NORMALIZATION OF RESONANCES DISTRIBUTION...
      call norm_init

      RETURN
      END
