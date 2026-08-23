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
      subroutine InitEpos(emax, lrescat, datpath, lpath, idbg, iou,
     &isigma0)   
c-----------------------------------------------------------------------
c General initialization of EPOS
c anfe: accepts nuclear PDG id instead of A, Z combos
c-----------------------------------------------------------------------
         include "epos.inc"
         real emax
c        lrescat dummy parameter to be consistent with Epos-LHC-R       
         logical lrescat
         integer lpath, idbg, iou, isigma0
         character(*) datpath

         ! dummy values
         seedi=1   !seed for random number generator: at start program
         seedj=2   !seed for random number generator: for first event

c Initialize decay of particles (all unstable decay)
         nrnody=0

         call LHCparameters        !LHC tune for EPOS
c isigma=1: cross-section is calculated by a numerical method
c           which is valid only for h-p or h-A (h being pion, kaon or nucleon)
c           but not A-B (nucleus-nucleus)
c           (not good for ionudi=2)
c
c isigma=0: same as isigma=1 but do not print the cross section on screen
c
c isigma=2: all the nuclear cross-sections are calculated by AA pseudo simulations
c           but it takes several minutes to compute

         isigma=isigma0 
         ionudi=1

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
         iwseed=0    !print seed to file
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
**FR fix PDG conversion for nuclei
         if (ippdg.ge.1000000000) then
            iapro = mod(ippdg, 10000) / 10
            izpro = mod(ippdg, 1000000) / 10000            
            iproton = 2212
            idprojin = idtrafo('pdg','nxs',iproton)
         else
            idprojin = idtrafo('pdg','nxs',ippdg)
            izpro = -1
            iapro = 1
         endif

         if (itpdg.ge.1000000000) then
            iatar = mod(itpdg, 10000) / 10
            iztar = mod(itpdg, 1000000) / 10000
            iproton = 2212
            idtargin = idtrafo('pdg','nxs',iproton)
         else
            idtargin = idtrafo('pdg','nxs',itpdg)         
            iztar = -1
            iatar = 1
         endif
         
         laproj = izpro      !proj Z
         maproj = iapro      !proj A         
         latarg = iztar      !targ Z
         matarg = iatar      !targ A

         iframe=10 + ifram           !nucleon-nucleon frame (12=target)

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
     &    ,xsloela,xsigqel)
c-----------------------------------------------------------------------
c     cross section function
c
c     xsigqel = quasi-elastic (non-excited projectile diffraction)
c     part of xsigine; the production cross section is
c     xsigine - xsigqel. 0 for hadron-proton (non-excited projectile
c     diffraction off a proton is elastic scattering); -1 if
c     unavailable (crseaaModel path).
c-----------------------------------------------------------------------

         implicit none
         include 'epos.inc'
         double precision xsigtot,xsigine,xsigela,xsigdd,xsigsd
     &                   ,xsloela,xsigqel
         real sigqelaa

Cf2py intent(out) xsigtot,xsigine,xsigela,xsigdd,xsigsd,xsloela,xsigqel

         xsigtot   = dble( sigtot   )
         xsigine   = dble( sigine   )
         xsigela   = dble( sigela   )
         xsigdd    = dble( sigdd    )
         xsigsd    = dble( sigsd    )
         xsloela   = dble( sloela   )
         xsigqel   = 0d0
c Nuclear cross section only if needed
         ! xsigtot = 0d0
         ! xsigine = 0d0
         ! xsigela = 0d0
         if(maproj.gt.1.or.matarg.gt.1)then
            if(model.eq.1)then
               call crseaaeposqel(sigtotaa,sigineaa,sigcutaa,sigelaaa
     &                           ,sigqelaa)
               xsigqel = dble( sigqelaa )
            else
               call crseaaModel(sigtotaa,sigineaa,sigcutaa,sigelaaa)
               xsigqel = -1d0
            endif
            xsigtot = dble( sigtotaa )
            xsigine = dble( sigineaa )
            xsigela = dble( sigelaaa )
         endif

         return
      end

c-----------------------------------------------------------------------
      subroutine crseaaeposqel(sigt,sigi,sigc,sige,sigql)
c-----------------------------------------------------------------------
c nucleus-nucleus (hadron) cross sections from epocrossc Glauber
c simulations (air-weighted if the target is air), including the
c quasi-elastic component
c  sigt  = sig tot
c  sigi  = sig inelastic (cut + all projectile diffraction)
c  sigc  = sig cut
c  sige  = sig elastic (includes target diffraction)
c  sigql = quasi-elastic part of sigi: projectile emerges non-excited
c          (0 if ionudi.ne.1, where it is counted as elastic)
c The production cross section (projectile destroyed) is sigi - sigql.
c
c epogcr separates gqel from the diffractive part only under ionudi=2,
c so ionudi=2 is used for the duration of the Glauber MC. The
c gqel/gdd split is pure arithmetic after gprod/gabs/gcoh are formed
c and draws no random numbers, so sigt/sigi/sigc/sige do not depend
c on this setting.
c-----------------------------------------------------------------------
      include 'epos.inc'
      niter=20000
      matarg0=matarg
      ionudi0=ionudi
      ionudi=2
      if(idtarg.eq.0)then
        sigt=0.
        sigc=0.
        sigi=0.
        sige=0.
        sigd=0.
        sigql=0.
        do k=1,3
          matarg=int(airanxs(k))
          call epocrossc(niter,xsigt,xsigi,xsigc,xsige,xsigql,xsigd)
          sigt=sigt+airwnxs(k)*xsigt
          sigi=sigi+airwnxs(k)*xsigi
          sigc=sigc+airwnxs(k)*xsigc
          sige=sige+airwnxs(k)*xsige
          sigd=sigd+airwnxs(k)*xsigd
          sigql=sigql+airwnxs(k)*xsigql
        enddo
        matarg=matarg0
      else
        call epocrossc(niter,sigt,sigi,sigc,sige,sigql,sigd)
      endif
      ionudi=ionudi0
      if(ionudi.ne.1)then
        sige=sige+sigql      !add non-excited diffractive projectile to elastic
        sigi=sigi-sigql      !do not count non-excited diffractive projectile in inelastic
        sigql=0.             !no quasi-elastic part left in sigi
        if(maproj+matarg.gt.2)then
          sigc=sigc+sigd*0.95   !for absorbtion cross section remove 5% of the
                                !excited projectile diffractive cross section
                                !which "looks like" non excited (approximation)
        endif
      endif
      end
