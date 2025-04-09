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
      subroutine InitEpos(emax, datpath, lpath, idbg, iou)
c-----------------------------------------------------------------------
c General initialization of EPOS
c anfe: accepts nuclear PDG id instead of A, Z combos
c-----------------------------------------------------------------------      
      include "epos.inc"
      real emax
      integer lpath, idbg, iou
      character(*) datpath
      
      ! dummy values
      seedi=1   !seed for random number generator: at start program
      seedj=2   !seed for random number generator: for first event

c Initialize decay of particles (all unstable decay)
      nrnody=0

      call LHCparameters        !LHC tune for EPOS
      isigma=2                  !use analytic cross section for nuclear xs
      ionudi=1

c      isigma=0              !do not print out the cross section on screen
c      ionudi=3              !count diffraction without excitation as elastic

      iecho=0                     !"silent" reading mode

      nfnnx=lpath   
      fnnx=datpath                        ! path to main epos subdirectory
      nfnii=lpath + 10                    ! epos tab file name lenght
      fnii=fnnx(1:nfnnx) // "epos.initl"  ! epos tab file name
      nfnid=lpath + 10
      fnid=fnnx(1:nfnnx) // "epos.inidi"
      nfnie=lpath + 10
      fnie=fnnx(1:nfnnx) // "epos.iniev"
      nfnrj=lpath + 10
      fnrj=fnnx(1:nfnnx) // "epos.inirj"       !'.lhc' is added a the end of the file name in ainit if LHCparameters is called
      nfncs=lpath + 10
      fncs=fnnx(1:nfnnx) // "epos.inics"       !'.lhc' is added a the end of the file name in ainit if LHCparameters is called

c Debug
      ish=idbg       !debug level
      ifch=iou      !debug output (screen)
c      ifch=31    !debug output (file)
c      fnch="epos.debug"
c      nfnch=index(fnch,' ')-1
c      open(ifch,file=fnch(1:nfnch),status='unknown')

      nevent = 1  !number of events
      modsho = 1  !printout every modsho events

      ecms=emax  !center of mass energy in GeV/c2
      
      istmax = 1      !only final particles (istmax=1 includes mother particles)

      End

c-----------------------------------------------------------------------
      subroutine InitEposEvt(ecm, ela, ippdg, itpdg, ifram)
c-----------------------------------------------------------------------
c Initialization to be called after changing the energy or beam
c configuration
c define either ecm < 0 and ela > 0 or ecm > 0 and ela < 0
c anfe: accepts nuclear PDG id instead of A, Z combos
c-----------------------------------------------------------------------      
      include "epos.inc"
      integer ifram

      engy = -1.
      ecms = -1.
      elab = -1.
      ekin = -1.
      pnll = -1.

      ecms=ecm  !center of mass energy in GeV/c2
      elab=ela  ! lab energy
      
      if (ippdg.ge.1000000000) then
         izpro = mod(ippdg, 1000) / 10
         iapro = mod(ippdg, 1000000) / 10000
         ippdg = 2212
      else
         izpro = 1
         iapro = 1
      endif
      if (itpdg.ge.1000000000) then
         iztar = mod(itpdg, 1000) / 10
         iatar = mod(itpdg, 1000000) / 10000
         itpdg = 2212
      else
         iztar = 1
         iatar = 1
      endif

      idprojin = idtrafo("pdg","nxs",ippdg)
      if (idpro.ne.2212) izpro = -1
      laproj = izpro      !proj Z
      maproj = iapro      !proj A
      idtargin = idtrafo("pdg","nxs",itpdg)
      latarg = iztar      !targ Z
      matarg = iatar      !targ A

      iframe=10 + ifram           !nucleon-nucleon frame (12=target)

c anfe why here? istmax = 1      !only final particles (istmax=1 includes mother particles)
      call ainit()
      End

c-----------------------------------------------------------------------
      real function GetCharge(idpdg)
c-----------------------------------------------------------------------
c Returns charge for particle with PDG ID
c-----------------------------------------------------------------------
      integer idpdg
C anfe This is a workaround for nuclear fragments in particle history
      if (idpdg.eq.90) then
         GetCharge=0
      else
         call idchrg(idtrafo("pdg","nxs",idpdg),GetCharge)
      endif
      return 

      End

c-----------------------------------------------------------------------
      subroutine xsection(xsigtot,xsigine,xsigela,xsigdd,xsigsd
     &    ,xsloela)
c-----------------------------------------------------------------------
c     cross section function
c-----------------------------------------------------------------------

      implicit none
      include 'epos.inc'
      double precision xsigtot,xsigine,xsigela,xsigdd,xsigsd
     &                ,xsloela

Cf2py intent(out) xsigtot,xsigine,xsigela,xsigdd,xsigsd,xsloela

      xsigtot   = dble( sigtot   )
      xsigine   = dble( sigine   )
      xsigela   = dble( sigela   )
      xsigdd    = dble( sigdd    )
      xsigsd    = dble( sigsd    )
      xsloela   = dble( sloela   )
c Nuclear cross section only if needed
      ! xsigtot = 0d0
      ! xsigine = 0d0
      ! xsigela = 0d0
      if(maproj.gt.1.or.matarg.gt.1)then
        if(model.eq.1)then
          call crseaaEpos(sigtotaa,sigineaa,sigcutaa,sigelaaa)
        else
          call crseaaModel(sigtotaa,sigineaa,sigcutaa,sigelaaa)
        endif
        xsigtot = dble( sigtotaa )
        xsigine = dble( sigineaa )
        xsigela = dble( sigelaaa )
      endif

      return
      end
