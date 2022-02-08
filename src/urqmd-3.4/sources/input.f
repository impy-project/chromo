c $Id: input.f,v 1.41 2007/01/30 14:50:25 bleicher Exp $
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine input(io)
c
c     Revision : 1.0
c
c     This subroutine reads the UQMD input file (unit=9)
C
c input : for ({\\tt io=0} input-file will be processed, else default values assumed
c output: information in common-block coms.f
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      implicit none

      include 'coms.f'
      include 'options.f'
      include 'comres.f'
      include 'inputs.f'
      include 'boxinc.f'

      character*3 flag
      character*77 inputstr,file9,fheader,file14,file15,file16,file17
      character*77 file13,file10,file19,file20
      integer line,proflg,tarflg,impflg,beamflg,index,ival,partid
      integer eosflg,i,io
      real*8 rval,caltim,outtim
      logical dtflag,bret

      character CTOStrng(numcto)*60
      character CTPStrg(numctp)*60

c setting of internal parameters values:
c SEH: useless apparently, can be removed in the python
      real*8 valint(1)
      common /values/ valint

      save

      valint(1)=0.d0

c SEH: bret is called after setting the CTOs (return if bret)
      bret=io.ne.0
      goto 108

c SEH: useless too
      entry inpini
c  called by some test programs

c SEH: bret is called after setting the CTOs (return if bret)
      bret=.true.
 108  continue

c initialize counters
c SEH: to be removed, used to read input file
      line=0
c  Line 135688 onwards from corsika
      boxflag=0
      mbflag=0
      edens=0.d0
      para=0
      solid=0
      mbox=0

c the following flags check, wether all necessary input is given
c projectile
c SEH: pro, same principle as beam
      proflg=0
      prspflg=0
c target
c SEH: tar, same principle as beam
      tarflg=0
      trspflg=0
c impact parameter
      impflg=0
c incident beam energy
c SEH beamflg is to check whether the beam energy is defined
c SEH and whether there are multiple definitions
      beamflg=0
      srtflag=0
      firstev=0
c equation of state
c SEH: eos, same principle as beam
      eosflg=0
c excitation function
      nsrt=1
         npb=1
      efuncflag=0
c default number of events
      nevents=1
c default seed for random number generator
      ranseed=0
c default number of timesteps
      nsteps=1000
c use standard time-step
      dtflag=.false.
c skip conditions on unit 14, 15, 16 & 18
      bf13=.false.
      bf14=.false.
      bf15=.false.
      bf16=.false.
      bf18=.false.
      bf19=.false.
      bf20=.false.
      do 111 i=1,numcto
         CTOdc(i)='  '
 111  continue
      do 112 i=1,numctp
         CTPdc(i)='  '
 112  continue
      do 113 i=1,maxstables
         stabvec(i)=0
 113  continue
      nstable = 0

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
c SEH: need to add these parameters wrt to corsika
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
      CTOption(1)=0
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
      CTOption(7)=0
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
      CTOption(24)=1
      CTOStrng(24)=' Wood-Saxon initialization'
      CTOption(25)=0
      CTOStrng(25)=' phase space corrections for resonance mass'
      CTOption(26)=0
      CTOStrng(26)=' use z -> 1-z for diquark-pairs'
      CTOption(27)=1
      CTOStrng(27)=' reference frame (1=target, 2=projectile, else=cms)'
      CTOption(28)=0
      CTOStrng(28)=' propagate spectators also '
      CTOption(29)=2
      CTOStrng(29)=' no transverse momentum in clustr '
      CTOption(30)=1
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
c SEH: parameters updated wrt to corsika
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
cbb Note: If you add more CTOptions, please make sure that all options
c   are included in the standard event header output in output.f.
c   Currently, 60 CTOs are written.
c SEH: remove bret

      if(bret)return


c initialize arrays for special PRO/TAR combinations
      do 10 i=1,2
         spityp(i)=0
         spiso3(i)=0
 10   continue

c SEH: what follows is all input/output stuff that we do not have to
care about. Also, URQINI in corsika stops there.
c header for output files
      fheader=' this is the default uqmd-fileheader'

c open fortran-unit 9 for input
c and units 14, 15 for output
      call getenv('ftn09',file9)
      call getenv('ftn10',file10)
      call getenv('ftn13',file13)
      call getenv('ftn14',file14)
      call getenv('ftn15',file15)
      call getenv('ftn16',file16)
      call getenv('ftn17',file17)
      call getenv('ftn19',file19)
      call getenv('ftn20',file20)
      if (file9(1:4).ne.'    ') then
         OPEN(UNIT=9,FILE=file9,STATUS='OLD',FORM='FORMATTED')
C peak into input file to find out if cto40 is specified.
         do
           read(9,99,end=9) flag,inputstr
           if(flag.eq.'cto') then
             read(inputstr,fmt=*,err=88,end=88) index,ival
             if(index.eq.40) then
               CTOption(40) = ival
               goto 9
             endif
           endif
         enddo
 9       continue
C re-open input file so that everything is fresh.
         rewind(UNIT=9)
      endif
      if (file10(1:4).ne.'    ') then
         OPEN(UNIT=10,FILE=file10,STATUS='OLD',FORM='FORMATTED')
c if CTO(40) is not yet set, but a file10 is given, then do set CTO(40).
         if (CTOption(40).eq.0) then
           CTOption(40)=2
           CTOdc(40)=' *'
C Don't do stuff like that silently.
           write(*,*) 'Old event detected: Assuming CTOption(40) = 2.'
         endif
C We only want as many events as are present in the input file. Yet, the
C number of events must be fixed before we read the input file. So, we
C set it to a ridiculously high number and then see how we can stop the
C evolution gracefully.
         nevents=100000
      endif
      if (file13(1:4).ne.'    ') then
         OPEN(UNIT=13,FILE=file13,STATUS='unknown',FORM='FORMATTED')
      endif
      if (file14(1:4).ne.'    ') then
         OPEN(UNIT=14,FILE=file14,STATUS='unknown',FORM='FORMATTED')
      endif
      if (file15(1:4).ne.'    ') then
         OPEN(UNIT=15,FILE=file15,STATUS='unknown',FORM='FORMATTED')
      endif
      if (file16(1:4).ne.'    ') then
         OPEN(UNIT=16,FILE=file16,STATUS='unknown',FORM='FORMATTED')
      endif
      if (file17(1:4).ne.'    ') then
         OPEN(UNIT=17,FILE=file17,STATUS='unknown',FORM='FORMATTED')
      endif
      if (file19(1:4).ne.'    ') then
         OPEN(UNIT=19,FILE=file19,STATUS='unknown',FORM='FORMATTED')
      endif
      if (file20(1:4).ne.'    ') then
         OPEN(UNIT=20,FILE=file20,STATUS='unknown',FORM='FORMATTED')
      endif
c
 99   format(1A3,1A77)


c stop input if old event is read in
      if(CTOption(40).ne.0) return


c this entry is used to read cto,ctp and tim statements
c in case of old event readin
      entry getparams

      rewind(UNIT=9)

c read input lines
 1    continue
      line=line+1
      read(9,99) flag,inputstr
 3    continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c select action according to flag:
c #  : treat line as a comment
      if(flag(1:1).eq.'#') goto 1
c blanks are comments, too
      if(flag(1:1).eq.' ') goto 1
c xxx: treat line as end of input marker
      if(flag.eq.'xxx'.or.flag.eq.'end') then
        Ap=Ap*int(CTParam(67))
        Zp=Zp*int(CTParam(67))
        At=At*int(CTParam(67))
        Zt=Zt*int(CTParam(67))
         goto 2
cc cal: header for output-files
c      if(flag.eq.'cal') then
c         fheader=inputstr
c pro: define projectile
      elseif(flag.eq.'pro') then
         proflg=proflg+1
         read(inputstr,fmt=*,err=88,end=88) Ap,Zp
         if(proflg.gt.1) then
            write(6,*)'multiple definitions for projectile system:'
            write(6,*)'-> last entry will be used'
         endif
c PRO: define special projectile
      elseif(flag.eq.'PRO') then
         proflg=proflg+1
         prspflg=1
         read(inputstr,fmt=*,err=88,end=88) spityp(1),spiso3(1)
         Ap=1
         if(proflg.gt.1) then
            write(6,*)'multiple definitions for projectile system:'
            write(6,*)'-> last entry will be used'
         endif
c tar: define target
      elseif(flag.eq.'tar') then
         tarflg=tarflg+1
         read(inputstr,fmt=*,err=88,end=88) At,Zt
         if(tarflg.gt.1) then
            write(6,*)'multiple definitions for target system:'
            write(6,*)'-> last entry will be used'
         endif
c TAR: define special target
      elseif(flag.eq.'TAR') then
         tarflg=tarflg+1
         trspflg=1
         read(inputstr,fmt=*,err=88,end=88) spityp(2),spiso3(2)
         At=1
         if(tarflg.gt.1) then
            write(6,*)'multiple definitions for target system:'
            write(6,*)'-> last entry will be used'
         endif
c box: define a box with a length in fm
c       parameters: 2: energie
c                   3: 1 =solid
c                   4: 1 = walls
        elseif(flag.eq.'box') then
          ! don't increase boxflag if we read old data.
          if(CTOption(40).eq.0) then
           boxflag=boxflag+1
          endif
          read(inputstr,fmt=*,err=88,end=88) lbox,edens,solid,para
           if (edens.gt.0.d0) edensflag=1

           if (lbox.le.0) then
              write(6,*) 'Error, length<=0'
              stop
           endif
           lboxhalbe=lbox/2.d0
           lboxd=lbox*2.d0

           if (edens.lt.0.d0) then
              write(6,*) 'Error, a negativ energy '
              stop
           endif

           if(boxflag.gt.1) then
            write(6,*)'multiple boxes are defined'
            stop
        endif
cbb: don't read particle definitions from input file if we have read an
cbb: old event: The particle list in the event header needs to be taken
cbb: from old event instead of the input file, so that the output is
cbb: consistent.
        elseif((flag.eq.'bpe'.or.flag.eq.'bpt').and.CTOption(40).ne.0)
     &   then
         write(6,*) 'ignoring box particle definition in inputfile, '
         write(6,*) 'because we have an old event read in.'
c bpt: define particles in the box
c parameters: ityp, iso3, mpart, pmax
        elseif(flag.eq.'bpt') then
           if (edens.gt.0.d0) then
              write(6,*) 'Error, energy is already defined'
              stop
           endif
           mbox=mbox+1
cbb: how do we know we didn't already add the maximum allowed number of
cbb: bpe lines? Turns out that maximum is not even documented. Bad.
cbb: Shame on you.
           if (mbox.gt.bptmax) then
              write(6,*) 'Error: Only ',bptmax,' particle definitions ',
     &                   'for box allowed. Increase paramater "bptmax"',
     &                   ' in boxinc.f .'
             stop 137
           endif
           read(inputstr,fmt=*,err=88,end=88)
     &     bptityp(mbox),bptiso3(mbox),bptpart(mbox),bptpmax(mbox)
           edensflag=0
           if (bptpart(mbox).le.0) then
              write(6,*) 'Error, a negativ particle number'
              stop
           endif
           if(boxflag.lt.1) then
            write(6,*)'no box is defined'
            stop
        endif
c bpe: define particles in the box with a given energy
c parameters: ityp, iso3, mpart,
        elseif(flag.eq.'bpe') then
           if (edens.le.0) then
              write(6,*) 'Error, no energy is defined'
              stop
           endif
           mbox=mbox+1
cbb: how do we know we didn't already add the maximum allowed number of
cbb: bpe lines? Turns out that maximum is not even documented. Bad.
cbb: Shame on you.
           if (mbox.gt.bptmax) then
              write(6,*) 'Error: Only ',bptmax,' particle definitions ',
     &                   'for box allowed. Increase paramater "bptmax"',
     &                   ' in boxinc.f .'
             stop 137
           endif
           read(inputstr,fmt=*,err=88,end=88)
     &     bptityp(mbox),bptiso3(mbox),bptpart(mbox)
           if(boxflag.lt.1) then
            write(6,*)'no box is defined'
            stop
        endif
c ene: beam energy (lab-system)
      elseif(flag.eq.'ene'.or.flag.eq.'elb') then
         beamflg=beamflg+1
         read(inputstr,fmt=*,err=88,end=88) ebeam
         if(beamflg.gt.1) then
            write(6,*)'multiple definitions for beam-energy:'
            write(6,*)'-> last entry will be used'
         endif
         if (ebeam.le.200) then
            write(6,*)'Calculation at ebeam.le.200 A GeV:'
            write(6,*)'parameter nmax in coms.f may be decreased!'
         endif
c plb: beam momentum (lab-system)
      elseif(flag.eq.'plb') then
         beamflg=beamflg+1
         srtflag=2
         read(inputstr,fmt=*,err=88,end=88) pbeam
         if(beamflg.gt.1) then
            write(6,*)'multiple definitions for beam-energy:'
            write(6,*)'-> last entry will be used'
         endif
       if (pbeam.le.200) then
            write(6,*)'Calculation at pbeam.le.200 A GeV:'
            write(6,*)'parameter nmax in coms.f may be decreased!'
       endif
c PLB: beam momentum ( LAb-system, excitation function possible)
      elseif(flag.eq.'PLB'.or.flag.eq.'PLG') then
         beamflg=beamflg+1
         srtflag=2
         read(inputstr,fmt=*,err=88,end=88) pbmin,pbmax,npb
         pbeam=pbmin
         if(beamflg.gt.1) then
            write(6,*)'multiple definitions for beam-energy:'
            write(6,*)'-> last entry will be used'
         endif
         if(npb.gt.1.and.flag.eq.'PLB') efuncflag=1
         if(npb.gt.1.and.flag.eq.'PLG') efuncflag=2
         if(abs(pbmax-pbmin).le.1.d-6) then
            npb=1
            efuncflag=0
         endif
         if (pbmax.le.200) then
            write(6,*)'Calculations at pbmax.le.200 A GeV:'
            write(6,*)'parameter nmax in coms.f may be decreased!'
         endif
c ecm:  c.m.energy
      elseif(flag.eq.'ecm') then
         beamflg=beamflg+1
         srtflag=1
         read(inputstr,fmt=*,err=88,end=88) ecm
         srtmin=ecm
         srtmax=ecm
         nsrt=1
         efuncflag=0
         if(beamflg.gt.1) then
            write(6,*)'multiple definitions for beam-energy:'
            write(6,*)'-> last entry will be used'
         endif
         if (ecm.le.20) then
            write(6,*)'Calculation at sroot.le.20 A GeV:'
            write(6,*)'parameter nmax in coms.f may be decreased!'
         endif
c ENE: beam energy (sqrt(s): CM-system, excitation function possible)
      elseif(flag.eq.'ENE'.or.flag.eq.'ELG') then
         beamflg=beamflg+1
         srtflag=1
         read(inputstr,fmt=*,err=88,end=88) srtmin,srtmax,nsrt
         ecm=srtmin
c        if(flag.eq.'ELG')ecm=1d1**dlog10(srtmin)
         if(beamflg.gt.1) then
            write(6,*)'multiple definitions for beam-energy:'
            write(6,*)'-> last entry will be used'
         endif
         if(nsrt.gt.1.and.flag.eq.'ENE') efuncflag=1
         if(nsrt.gt.1.and.flag.eq.'ELG') efuncflag=2
         if(abs(srtmax-srtmin).le.1.d-6) then
            nsrt=1
            efuncflag=0
         endif
         if (srtmax.le.20) then
            write(6,*)'Calculations at srootmax.le.20 A GeV:'
            write(6,*)'parameter nmax in coms.f may be decreased!'
         endif
c imp: impact parameter
      elseif(flag.eq.'imp') then
         bmin=0.d0
         impflg=impflg+1
         read(inputstr,fmt=*,err=88,end=88) bdist
         if(bdist.lt.0d0)then
           CTOption(5)=1
           bdist=abs(bdist)
           write(6,*)'randomly choosen impact parameter:',
     ,             ' CTOption(5) is set to 1'
         end if
         if(impflg.gt.1) then
            write(6,*)'multiple definitions for impact parameter:'
            write(6,*)'-> last entry will be used'
         endif
c IMP: impact parameter
      elseif(flag.eq.'IMP') then
         impflg=impflg+1
         read(inputstr,fmt=*,err=88,end=88) bmin,bdist
         CTOption(5)=1
         if(impflg.gt.1) then
            write(6,*)'multiple definitions for impact parameter:'
            write(6,*)'-> last entry will be used'
         endif
c eos: impact parameter
      elseif(flag.eq.'eos') then
         eosflg=eosflg+1
         read(inputstr,fmt=*,err=88,end=88) eos
         if(eosflg.gt.1) then
            write(6,*)'multiple definitions for equation of state:'
            write(6,*)'-> last entry will be used'
         endif
         if (eos.ne.0) then
            CTOption(24)=0
         endif
c nev: number of events
      elseif(flag.eq.'nev') then
         read(inputstr,fmt=*,err=88,end=88) nevents
c rsd:
      elseif(flag.eq.'rsd') then
         read(inputstr,fmt=*,err=88,end=88) ranseed
c cdt: collision time step
      elseif(flag.eq.'cdt') then
         read(inputstr,fmt=*,err=88,end=88) dtimestep
         dtflag=.true.
c tim: time of propatation
      elseif(flag.eq.'tim') then
         read(inputstr,fmt=*,err=88,end=88) caltim, outtim
c stb: keep particle stable
      elseif(flag.eq.'stb') then
         read(inputstr,fmt=*,err=88,end=88) partid
         if (nstable.lt.maxstables) then
            nstable = nstable + 1
            stabvec(nstable) = partid
         else
            write(6,*) 'Warning: too many stable particles defined!'
         endif
c cto: collision term options
      elseif(flag.eq.'cto') then
         read(inputstr,fmt=*,err=88,end=88) index,ival
         write(6,'(a9,i3,a2,x,i3," is changed to ",x,i3," (",a60,")")')
     &        'CTOption(',index,')=',CTOption(index),ival
     &        ,CTOStrng(index)
         CTOption(index)=ival
         CTOdc(index)=' *'
c ctp: collision term parameter
      elseif(flag.eq.'ctp') then
         read(inputstr,fmt=*,err=88,end=88) index,rval
         write(6,'(a8,i3,a2,x,e10.4,a15,x,e10.4," (",a60,")")')
     ,           'CTParam(',index,')=',CTParam(index)," is changed to "
     ,           ,rval,CTPStrg(index)
         CTParam(index)=rval
         CTPdc(index)=' *'
      elseif (flag.eq.'f13') then
         bf13=.true.
         if (info) write(6,*)'(info) no output on unit 13'
      elseif (flag.eq.'f14') then
         bf14=.true.
         if (info) write(6,*)'(info) no output on unit 14'
      elseif (flag.eq.'f15') then
         bf15=.true.
         if (info) write(6,*)'(info) no output on unit 15'
      elseif (flag.eq.'iou') then
         read(inputstr,fmt=*,err=88,end=88) index,ival
         call uounit(index,ival)
         write(6,*)'file',index,'will be written on unit',ival
      elseif (flag.eq.'f16') then
         bf16=.true.
         if (info) write(6,*)'(info) no output on unit 16'
      elseif (flag.eq.'f18') then
          bf18=.true.
          if (info) write(6,*)'(info) no output on unit 18'
      elseif (flag.eq.'f19') then
          bf19=.true.
          if (info) write(6,*)'(info) no output on unit 19'
      elseif (flag.eq.'f20') then
          bf20=.true.
          if (info) write(6,*)'(info) no output on unit 20'
      else
         write(6,*)'undefined opcode in input-file on line',line
         stop
      endif
      goto 1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 2    continue

c SEH: in URQLNK in corsika -- start block
c fast CASCADE mode
      if(.not.dtflag.and.eos.eq.0) dtimestep=outtim
c
      nsteps=int(0.01+caltim/dtimestep)
      outsteps=int(0.01+outtim/dtimestep)


c stop input if old event is read in
      if(CTOption(40).ne.0) return
c SEH: in URQLNK in corsika -- end block


c here some validity checks of the input should be performed
        if (boxflag.eq.1.and.mbox.eq.0) then
            Write(6,*) 'Error: no particles in the box.'
            stop 137
        elseif (boxflag.eq.1.and.mbox.gt.bptmax) then
            Write(6,*) 'Error: too many particles in the box, only ',
     &                 bptmax, ' allowed'
            stop 137
        ElseIf (boxflag.eq.0) then
      if(proflg.eq.0) then
         write(6,*)'Error: no projectile specified in input.'
         stop
      elseif(tarflg.eq.0) then
         write(6,*)'Error: no target specified in input.'
         stop
      elseif((impflg.eq.0)) then
         write(6,*)'Error: no impact parameter in input.'
         stop
      elseif(beamflg.eq.0.and.prspflg.eq.0) then
         write(6,*)'Error: no incident beam energy specified.'
         stop
      endif
c EndIf for the Box
        EndIf
      if (efuncflag.ne.0.and.
     &    mod(nevents,max(nsrt,npb)).ne.0) then
         write(6,*)'INPUT: the number of events divided by the ',
     ,   'number of energies requested is no integer.'
      end if
c SEH: end of the input/output stuff
c
c constraints for skyrme pots:
c SEH: in URQLNK in corsika
      if(eos.ne.0.and.((srtflag.eq.0.and.ebeam.gt.4d0)
     &             .or.(srtflag.eq.1.and.srtmax.gt.3.3d0)
     &             .or.(srtflag.eq.2.and.pbeam.gt.4.9))) then
         write(6,*)'***(W) I switched off the potentials'
         eos=0
      end if
      if(eos.ne.0) then
         CTOption(11)=1
         CTOption(28)=0
         CTOption(30)=0
      endif
c

c SEH: noch mal some input/output stuff
c now print the selected analysis

c...some input combinations should be avoided and/or commented
      if(CTOption(7).ne.0.and.At*Ap.ne.1)then
        write(6,*)'Warning: CTOption(7)=',CTOption(7),
     ,  ' no elastic collisions in NN',
     ,  ' should not be used for serious calculations!'
      end if

      if(CTOption(18).ne.0)then
        write(6,*)'Warning: CTOption(18)=',CToption(18),': ',
     ,  'unstable particles will not decay after propagation.'
      end if


      if(CTOption(31).ne.0)then
        write(6,*)'Warning: CTOption(31)=',CToption(31),': ',
     ,  "Not yet completly implemented. Don't use for serious",
     ,  'calculations (not yet..).'
      end if

      if(CTParam(28).lt.0d0.or.CTParam(28).gt.1d0)then
        write(6,*)'Warning: CTParam(28)=',CTParam(28),': ',
     ,  'should be between 0 and 1. it will be corrected.'
        CTParam(28)=min(1d0,max(0d0,CTParam(28)))
      end if

      return

 88   write(6,*) 'syntax-error in input-file on line ',line,flag
      write(6,*)inputstr
      stop
      end
