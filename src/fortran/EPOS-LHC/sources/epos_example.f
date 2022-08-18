c 15.01.2009 Simplified Main program and random number generator for epos

c-----------------------------------------------------------------------
      program aamain
c-----------------------------------------------------------------------

c HEP common as defined in epos.inc
      parameter (nmxhep=9990)       !max nr of particles in hep ptl list
      double precision phep,vhep
      common/hepevt/nevhep,nhep,isthep(nmxhep),idhep(nmxhep),
     &jmohep(2,nmxhep),jdahep(2,nmxhep),phep(5,nmxhep),vhep(4,nmxhep)
   
      common/cevt/phievt,nevt,bimevt,kolevt,koievt,pmxevt,egyevt,npjevt
     *,ntgevt,npnevt,nppevt,ntnevt,ntpevt,jpnevt,jppevt,jtnevt,jtpevt
     *,xbjevt,qsqevt,nglevt,zppevt,zptevt,minfra,maxfra
     *,typevt
  
c Set parameters to default value
      call aaset(0)

c Update some parameters value to run correctly for LHC
      call IniEpos(nevent,ish)

c The parameters can be changed optionnaly by reading a file (example.param) using the following subroutine
      ! call EposInput(nevent,ish)  !(it can be commented)
c if you put what is in input.optns in example.param, you can even run exactly the same way (coded parameters are overwritten). Don't forget the command :
c "EndEposInput"
c at the end of example.param, otherwise it will not run.
      
c initialization for the given energy
      call ainit

      do  n=1,nevent
c Calculate an inelastic event
        call aepos(-1)
c Fix final particles and some event parameters
        call afinal
c Fill HEP common
        call hepmcstore   !complete filling of hepevt arrays compatible with hepmc library
c        call ustore !user defined in epos-bas-*.f
c Print out (same as the one defined in input.optns)
        if(nhep.gt.nmxhep)then
          print *,'Warning : produced number of particles is too high'
          print *,'          increase nmxhep : ',nhep,' > ',nmxhep
          stop
        endif
        write(*,*)nevhep,nhep,typevt !typevt is type of collision : 1=Non Diffractive, 2=Double Diffractive, 3=Single Diffractive
        do i=1,nhep
          write(*,'(i5,3x,i2,2x,2i5,2x,2i5)')i,isthep(i)
     *               ,jmohep(1,i),jmohep(2,i),jdahep(1,i),jdahep(2,i)
          write(*,'(i10,1x,4(e12.6,1x))')idhep(i),(phep(k,i),k=1,4)
        enddo
      enddo
c optional Statistic information (only with debug level ish=1)
      call astati
      if(ish.ge.1)call bfinal

      end


c-----------------------------------------------------------------------
      subroutine IniEpos(nevto,isho)
c-----------------------------------------------------------------------
c Update some parameters and define path to tab files
c here can be set what particle is stable or not
c transfer number of event and debug level to main program (to avoid
c epos.inc in main)
c-----------------------------------------------------------------------
      include "epos.inc"
      
      seedi=1.d0   !seed for random number generator: at start program
      seedj=1.d0   !seed for random number generator: for first event

c Initialize decay of particles
      nrnody=0       !number of particle types without decay (if 0 (default) : all unstable particles decay (at the end only (anti)nucleons, (anti)electrons and muons)
c Particle code is given as
c     id=+/-ijkl
c
c          mesons--
c          i=0, j<=k, +/- is sign for j
c          id=110 for pi0, id=220 for eta, etc.
c
c          baryons--
c          i<=j<=k in general
c          j<i<k for second state antisymmetric in (i,j), eg. l = 2130
c
c          other--
c          id=1,...,6 for quarks
c          id=9 for gluon
c          id=10 for photon
c          id=11,...,16 for leptons
c          i=17 for deuteron
c          i=18 for triton
c          i=19 for alpha
c          id=20 for ks, id=-20 for kl
c
c          i=21...26 for scalar quarks
c          i=29 for gluino
c
c          i=30 for h-dibaryon
c
c          i=31...36 for scalar leptons
c          i=39 for wino
c          i=40 for zino
c
c          id=80 for w+
c          id=81,...,83 for higgs mesons (h0, H0, A0, H+)
c          id=84,...,87 for excited bosons (Z'0, Z''0, W'+)
c          id=90 for z0
c
c          diquarks--
c          id=+/-ij00, i<j for diquark composed of i,j.
c
c Examples : 2130 = lambda, 1330=xi0, 2330=xi-, 3331=omega
c
c Conversion from epos to  pdg code can be done using
c      id_pdg=idtrafo('nxs','pdg',id_epos)

      nrnody=nrnody+1
      nody(nrnody)=1220    !neutron
      nrnody=nrnody+1
      nody(nrnody)=-1220   !aneutron
      nrnody=nrnody+1
      nody(nrnody)=120     !pi+
      nrnody=nrnody+1
      nody(nrnody)=-120    !pi-
      nrnody=nrnody+1
      nody(nrnody)=130     !K+
      nrnody=nrnody+1
      nody(nrnody)=-130    !K-
      nrnody=nrnody+1
      nody(nrnody)=-20     !Kl
      nrnody=nrnody+1
      nody(nrnody)=-14     !mu+
      nrnody=nrnody+1
      nody(nrnody)=14      !mu-
c      nrnody=nrnody+1
c      nody(nrnody)=idtrafo('pdg','nxs',3122)    !lambda using pdg code

c      ctaumin=1.           !or min decay lenght c*tau can be defined here in cm
                           !(all particles with c*tau>ctaumin are stable)

      call LHCparameters        !LHC tune for EPOS
      isigma=1                  !use analytic cross section for nuclear xs
      ionudi=1

c      isigma=0              !do not print out the cross section on screen
c      ionudi=3              !count diffraction without excitation as elastic

      iframe=11                 !nucleon-nucleon frame (12=target)
      iecho=0                     !"silent" reading mode

      nfnnx=2   
      fnnx="./"                    ! path to main epos subdirectory
      nfnii=10                     ! epos tab file name lenght
      fnii="epos.initl"            ! epos tab file name
      nfnid=10
      fnid="epos.inidi"
      nfnie=10
      fnie="epos.iniev"
      nfnrj=10
      fnrj="epos.inirj"       !'.lhc' is added a the end of the file name in ainit if LHCparameters is called
      nfncs=10
      fncs="epos.inics"       !'.lhc' is added a the end of the file name in ainit if LHCparameters is called

c Debug
      ish=0       !debug level
      ifch=6      !debug output (screen)
c      ifch=31    !debug output (file)
c      fnch="epos.debug"
c      nfnch=index(fnch,' ')-1
c      open(ifch,file=fnch(1:nfnch),status='unknown')


      nevent = 1000  !number of events
      modsho = 1  !printout every modsho events

      ecms=7000  !center of mass energy in GeV/c2
      
      idproj = 1120   !proton
      laproj = 1      !proj Z
      maproj = 1      !proj A
      idtarg = 1120   !proton
      latarg = 1      !targ Z
      matarg = 1      !targ A

      istmax = 0      !only final particles (istmax=1 includes mother particles)

c for main program
      nevto  = nevent
      isho   = ish

      end



c-----------------------------------------------------------------------     
      subroutine EposInput(nevto,isho)
c-----------------------------------------------------------------------     
c Read informations (new options or parameter change) in the file
c "epos.param". The unit "ifop" is used in aread. If not used, it will
c use the default value of all parameters.
c-----------------------------------------------------------------------     
      include "epos.inc"   
      nopen=0
      ifop=35
      open(unit=ifop,file='example.param',status='old')
      call aread
      close(ifop)
c for main program
      nevto  = nevent
      isho   = ish
      end

c-----------------------------------------------------------------------
      subroutine xsection(xsigtot,xsigine,xsigela,xsigdd,xsigsd
     &                          ,xsloela,xsigtotaa,xsigineaa,xsigelaaa)
c-----------------------------------------------------------------------
c     cross section function
c-----------------------------------------------------------------------
      implicit none
      include 'epos.inc'
      double precision xsigtot,xsigine,xsigela,xsigdd,xsigsd
     &                ,xsloela,xsigtotaa,xsigineaa,xsigelaaa

      xsigtot   = dble( sigtot   )
      xsigine   = dble( sigine   )
      xsigela   = dble( sigela   )
      xsigdd    = dble( sigdd    )
      xsigsd    = dble( sigsd    )
      xsloela   = dble( sloela   )
c Nuclear cross section only if needed
      xsigtotaa = 0d0
      xsigineaa = 0d0
      xsigelaaa = 0d0
      if(maproj.gt.1.or.matarg.gt.1)then
        if(model.eq.1)then
          call crseaaEpos(sigtotaa,sigineaa,sigcutaa,sigelaaa)
        else
          call crseaaModel(sigtotaa,sigineaa,sigcutaa,sigelaaa)
        endif
        xsigtotaa = dble( sigtotaa )
        xsigineaa = dble( sigineaa )
        xsigelaaa = dble( sigelaaa )
      endif

      return
      end
