c -*- mode: Fortran -*-
c---------------------------------------------------------------------------
c                         dimensions
c---------------------------------------------------------------------------

      integer   mmry,mxptl,nmxhep,myptl,nzeta,nflav,mxstr,mystr,mxtau
     *         ,mxtrig,mxpri,mxbins,matau,mxnucl,mxhisarg,idxD0,idxD1
     *         ,idxD,nclha,nclegy,mamxx,mxjerr,mxvol,mxeps,mxidco,mxcoox
     *         ,mxcooy
      parameter (mmry=1)   !memory saving factor

      parameter (mxptl=300000/mmry) !max nr of particles in epos ptl list
      parameter (nmxhep=9990)    !max nr of particles in hep ptl list
      parameter (myptl=1000)        !max nr of droplets in epos ptl list
      parameter (nzeta=60)          !max nr of zeta bins for droplets
      parameter (nflav=6)           !max nr of flavors
      parameter (mxstr=20000/mmry)  !max nr of strings in epos string list 
      parameter (mystr=20000/mmry)    
      parameter (mxtau=4,mxvol=10,mxeps=16)
      parameter (mxtrig=99,mxidco=99)
      parameter (mxpri=200) 
      parameter (mxbins=10000)
      parameter (matau=10,mxcoox=40,mxcooy=10)
      parameter (mxnucl=20)
      parameter (mxhisarg=100)
      parameter (idxD0=0,idxD1=2,idxD=1,nclha=4,nclegy=100) 
      parameter (mamxx=250)
      parameter (mxjerr=10)
           
c---------------------------------------------------------------------------
c                   epos event common block
c---------------------------------------------------------------------------

      real        phievt,bimevt,pmxevt,egyevt
     *,xbjevt,qsqevt,zppevt,zptevt
      integer     nevt,kolevt,koievt,kohevt,npjevt
     *,ntgevt,npnevt,nppevt,ntnevt,ntpevt,jpnevt,jppevt,jtnevt,jtpevt
     *,nglevt,minfra,maxfra,npglb,ntglb
      common/cevt/phievt,nevt,bimevt,kolevt,koievt,pmxevt,egyevt,npjevt
     *,ntgevt,npnevt,nppevt,ntnevt,ntpevt,jpnevt,jppevt,jtnevt,jtpevt
     *,xbjevt,qsqevt,nglevt,zppevt,zptevt,minfra,maxfra,kohevt
     *,npglb,ntglb
      real         rglevt,sglevt,eglevt,fglevt,typevt,segevt,ctrevt
      integer      ng1evt,ng11evt,ng12evt,ngspecp,ngspecn,ng2evt
     *,ikoevt,ikhevt
      common/c2evt/rglevt,sglevt,eglevt,fglevt,typevt,segevt,ctrevt
     *,ng1evt,ng11evt,ng12evt,ngspecp,ngspecn,ng2evt,ikoevt,ikhevt

c     nevt .......... error code. 1=valid event, 0=invalid event
c     bimevt ........ absolute value of impact parameter
c     phievt ........ angle of impact parameter
c     kolevt ........ number of collisions
c     koievt ........ number of inelastic collisions
c     kohevt ........ number of hard collisions
c     pmxevt ........ reference momentum
c     egyevt ........ pp cm energy (hadron) or string energy (lepton)
c     npjevt ........ number of primary projectile participants
c     ntgevt ........ number of primary target participants
c     npnevt ........ number of primary projectile neutron spectators
c     nppevt ........ number of primary projectile proton spectators
c     ntnevt ........ number of primary target neutron spectators
c     ntpevt ........ number of primary target proton spectators
c     jpnevt ........ number of absolute projectile neutron spectators
c     jppevt ........ number of absolute projectile proton spectators
c     jtnevt ........ number of absolute target neutron spectators
c     jtpevt ........ number of absolute target proton spectators
c     xbjevt ........ bjorken x for dis     
c     qsqevt ........ q**2 for dis  
c     sigtot ........ total cross section
c     nglevt ........ number of collisions acc to  Glauber 
c     zppevt ........ average Z-parton-proj 
c     zptevt ........ average Z-parton-targ
c     ng1evt ........ number of Glauber participants with at least one IAs 
c     ng2evt ........ number of Glauber participants with at least two IAs 
c     ikoevt ........ number of elementary parton-parton scatterings
c     ikoevt ........ number of elementary hard parton-parton scatterings
c     typevt ........ type of event (1=Non Diff, 2=Double Diff, 3=Central Diff, 4=AB->XB, -4=AB->AX)

c---------------------------------------------------------------------------
c                   epos particle list common block
c---------------------------------------------------------------------------
      
      real        pptl,tivptl,xorptl
      integer     nptl,iorptl,idptl,istptl,ifrptl,jorptl,ibptl,ityptl
      real         ekievt
      integer      itsptl
      common/cptl/nptl,pptl(5,mxptl),iorptl(mxptl),idptl(mxptl)
     *,istptl(mxptl),tivptl(2,mxptl),ifrptl(2,mxptl)
     *,jorptl(mxptl)
     *,xorptl(4,mxptl),ibptl(4,mxptl),ityptl(mxptl)
      common/c1ptl/ekievt,itsptl(mxptl)

c     nptl .......... current particle index (=number of ptls stored)
c     idptl(i) ...... particle id
c     pptl(1,i) ..... x-component of particle momentum 
c     pptl(2,i) ..... y-component of particle momentum 
c     pptl(3,i) ..... z-component of particle momentum 
c     pptl(4,i) ..... particle energy 
c     pptl(5,i) ..... particle mass 
c     iorptl(i) ..... particle number of father (if .le. 0 : no father) 
c     jorptl(i) ..... particle number of mother (if .le. 0 : no mother)
c     istptl(i) ..... status: 40 and 41 : Remnant
c                             30 and 31 : Pomeron
c                             20 and 21 : Parton
c                             10 and 11 : Droplet
c                             00 and 01 : Particle
c                            last digit = 0 : last generation
c                            last digit = 1 : not last generation
c     xorptl(1,i) ... x-component of formation point
c     xorptl(2,i) ... y-component of formation point
c     xorptl(3,i) ... z-component of formation point
c     xorptl(4,i) ... formation time
c     tivptl(1,i) ... formation time (always in the pp-cms!)
c     tivptl(2,i) ... destruction time (always in the pp-cms!)
c     ityptl(i)  .... type of particles origin:
c                         10-19: target
c                         20-29: soft Pom
c                         30-39: hard Pom 
c                         40-49: projectile 
c                         50: string, droplet
c     itsptl(i) ..... string type of particles origin (if string)  

      real         radptl
      integer      iaaptl
      real         desptl,dezptl
      common/c2ptl/iaaptl(mxptl),radptl(mxptl)
      common/c3ptl/desptl(mxptl),dezptl(mxptl)
      integer      nptlbd,nptlpt
      common/c4ptl/nptlbd,nptlpt
      real         rinptl,vrad
      integer      inbxxx,icccptl
      real         qsqptl,zpaptl
      common/c6ptl/rinptl(mxptl),vrad,inbxxx
      common/c8ptl/qsqptl(mxptl),zpaptl(2,mxptl)
      common/c9ptl/icccptl


c---------------------------------------------------------------------------
c                   hep standard event commonblock.
c---------------------------------------------------------------------------

      double precision phep,vhep
      integer       nevhep,nhep,isthep,idhep,jmohep,jdahep

      common/hepevt/nevhep,nhep,isthep(nmxhep),idhep(nmxhep),
     &jmohep(2,nmxhep),jdahep(2,nmxhep),phep(5,nmxhep),vhep(4,nmxhep)
     
c---------------------------------------------------------------------------
c
c         nevhep      -   event number
c         nhep        -   number of entries in the event record
c
c         isthep(i)   -   status code
c         idhep(i)    -   particle id (particle data group standard)
c
c         jmohep(1,i) -   position of mother particle in list
c         jmohep(2,i) -   position of second mother particle in list
c         jdahep(1,i) -   position of first daughter in list
c         jdahep(2,i) -   position of first daughter in list
c
c         phep(1,i)   -   p_x momentum in gev/c
c         phep(2,i)   -   p_y momentum in gev/c
c         phep(3,i)   -   p_z momentum in gev/c
c         phep(4,i)   -   energy in gev
c         phep(5,i)   -   mass in gev/c**2
c
c         vhep(1,i)   -   x position of production vertex in mm 
c         vhep(2,i)   -   y position of production vertex in mm 
c         vhep(3,i)   -   z position of production vertex in mm 
c         vhep(4,i)   -   time of production  in mm/c 
c
c          (note:  1 mm = 10^-12 fm = 5.07 10^-12 1/gev)

c------------------------------------------------------------------------
c  Parameters set in sr aaset and variables to communicate between moduls
c------------------------------------------------------------------------

      integer      ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      common/files/ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      integer       ifin
      common/files2/ifin
      character*500  fnch,fnhi,fndt,fnii,fnid,fnie,fnrj,fnmt
     * ,fngrv,fncp,fnnx,fncs,fndr,fnhpf
      common/fname/  fnch, fnhi, fndt, fnii, fnid, fnie, fnrj, fnmt
     * ,fngrv,fncp,fnnx,fncs,fndr,fnhpf
      integer       nfnch,nfnhi,nfndt,nfnii,nfnid,nfnie,nfnrj,nfnmt
     *,nfngrv,nfncp,nfnnx,nfncs,nfndr,nfnhpf
      common/nfname/nfnch,nfnhi,nfndt,nfnii,nfnid,nfnie,nfnrj,nfnmt
     *,nfngrv,nfncp,nfnnx,nfncs,nfndr,nfnhpf
      character*500  fnin
      common/fname2/ fnin
      integer        nfnin
      common/nfname2/nfnin
      real         delvol,deleps,dlzeta,etafac,facnuc,taurea,epscri
     *,taustr
      common/resc2/delvol,deleps,dlzeta,etafac,facnuc,taurea,epscri(3)
     *,taustr
      real         epsdfm,tauhac,radeft,radeft1,radeft2
      integer      iostro
      common/resc4/epsdfm,tauhac,radeft,iostro,radeft1,radeft2
      integer leadcore
      common/resc5/leadcore
      integer nsegmincore,ratiomaxcore
      common/core3/nsegmincore,ratiomaxcore 
      real         pud,pmqu,pmqd,pmqs,pmqc,pmqq
      integer      ndecay,maxres
      common/frag1/ndecay,maxres,pud,pmqu,pmqd,pmqs,pmqc,pmqq
      character*3  hydt
      common/hydr1/hydt
      integer      izmode,maxsurf,iofrout,jzmode
      real         zclass
      common/hydr4/zclass(5,100),izmode,maxsurf,iofrout,jzmode(7)
      real         zlimit
      integer      jtable
      common/hydr5/zlimit(0:100),jtable(7)
      integer      mrclass,ihdim,kclass
      real         rclass
      parameter(mrclass=9)
      common/hydr7/rclass(mrclass,3,100),ihdim(3),kclass(mrclass,3,100)
      integer      nrclass,irclass
      common/hydr8/nrclass(mrclass),irclass(mrclass)
      real         pdiqua,delrex,ptfraqq,ptfra,ptfrasr
      integer      ioptf
      common/frag2/pdiqua,delrex,ptfraqq,ptfra,ptfrasr,ioptf
      real         aouni,pbreak,fkappa,strcut,diqcut,fkappag,pbreakg
     *,zetacut
      integer      itflav
      common/frag3/aouni,pbreak,fkappa,itflav,strcut,diqcut
     *,fkappag,pbreakg,zetacut
      real         difud,difus,difuc,pudd,puds,pudc,difuuu,difuud
     *,difuus,difuuc,difudd,difuds,difudc,difuss,difusc,difucc
      integer      nrflav
      common/frag4/difud,difus,difuc,pudd,puds,pudc,difuuu,difuud
     *,difuus,difuuc,difudd,difuds,difudc,difuss,difusc,difucc,nrflav
      real         qmass
      integer      isospin
      common/frag5/qmass(0:6),isospin(0:6)
      real         pnll,ptq,exmass,cutmss,wproj,wtarg
      common/hadr1/pnll,ptq,exmass,cutmss,wproj,wtarg
      real          rstrau,rstrad,rstras,rstrac,rstrasi
      common/hadr10/rstrau(4),rstrad(4),rstras(4),rstrac(4),rstrasi
      real          wgtqqq,wgtval,wgtsea,wgtdiq
      common/wgtqrk/wgtqqq(4),wgtval,wgtsea,wgtdiq
      double precision timeini,timefin
      common/time1/timeini,timefin
      integer       iotst1,iotst2,iotst3,iotst4
      common/ciotst/iotst1,iotst2,iotst3,iotst4
      integer mxparam
      parameter(mxparam=99)
      real vparam
      common/cparam/vparam(mxparam)
      real         tauzer,tauone,tautwo,deltau,numtau,amsiac,amprif
     *,tauup,tauthree
      common/resc1/tauzer,tauone,tautwo,deltau,numtau,amsiac,amprif
     *,tauup,tauthree
      real         dscale,cepara,delamf,deuamf
      integer      iceopt
      common/resc3/dscale,cepara,iceopt,delamf,deuamf
      integer      ispherio,icotabm,icotabr,icocore
      common/sprio/ispherio,icotabm,icotabr,icocore
      integer       ihacas
      common/chacas/ihacas
      integer      iaddplot,maddplot,mxaddplot1,mxaddplot2
      parameter (mxaddplot1=100,mxaddplot2=10)
      character*20 caddplot
      common/root7/iaddplot,maddplot(mxaddplot1)
     .             ,caddplot(mxaddplot1,1:mxaddplot2)
       real         cutico,dssico
      common/incon/cutico,dssico
      real         gaumx
      integer      istore,istmax,irescl,ntrymx,nclean,iopdg,ioidch
      common/othe1/istore,istmax,gaumx,irescl,ntrymx,nclean,iopdg,ioidch
      integer      ifrade,iframe,idecay,ihdecay,jdecay,iremn
      common/othe2/ifrade,iframe,idecay,ihdecay,jdecay,iremn
      integer      jframe,kframe,iposi,istfor
      common/othe3/jframe,kframe,iposi,istfor
      integer      iselect
      common/othe4/iselect
      real           facposz,facposf,tauzer1,tauzer2
      common/ccorcor/facposz,facposf,tauzer1,tauzer2
      real           fitermet,felamet
      common/cfacmet/fitermet,felamet
      real         taumx,sigj
      integer      jpsi,jpsifi,nsttau,ijphis,ijtauan
      common/jpsif/jpsi,jpsifi,taumx,nsttau,sigj,ijphis,ijtauan
      real         themas
      integer      iopenu
      common/strlt/iopenu,themas
      integer      iappl,model
      common/appli/iappl,model
      integer      nevent,nfull,nfreeze,ninicon
      common/events/nevent,nfull,nfreeze,ninicon
      real         egymin,egymax,elab,ecms,ekin
      common/enrgy/egymin,egymax,elab,ecms,ekin
      integer      iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi
      common/prnt1/iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi
      real         engy,elepti,elepto,angmue
      integer      icinpu
      common/lept1/engy,elepti,elepto,angmue,icinpu
      real         core,fctrmx
      integer       laproj,maproj,latarg,matarg
      common/nucl1/laproj,maproj,latarg,matarg,core,fctrmx
      real         bmaxim,bminim,phimax,phimin
      common/nucl2/bmaxim,bminim,phimax,phimin
      real         ymximi,wtimet,wtimei,wtimea
      integer      imihis,iclhis,iwtime
      common/wana1/ymximi,imihis,iclhis,iwtime,wtimet,wtimei,wtimea
      real         wtmini,wtstep
      integer      isphis,ispall,iwcent,iana,nbdky
      common/wana2/isphis,ispall,wtmini,wtstep,iwcent,iana,nbdky
      real         asuhax,asuhay
      common/drop4/asuhax(7),asuhay(7)
      real         grigam,grirsq,gridel,grislo,gricel,sigppi,sigppd
      common/gribo/grigam,grirsq,gridel,grislo,gricel,sigppi,sigppd
      real         bag4rt,dezzer,amuseg,taurem,yradmx,facts,factb,factq
      common/drop3/bag4rt,dezzer,amuseg,taurem,yradmx,facts,factb,factq
      real         rcoll,ylongmx,facecc,yradpp,yradmi,yrmaxi,fradflii
     &            ,yfgeo,ydslg,ydsrd,dsegce,aminclu,ycori
      integer      nsegsu
      common/drop2/rcoll,ylongmx,facecc,yradpp,yradmi,dsegce,aminclu
     &            ,yrmaxi,fradflii,ycori,yfgeo,ydslg,ydsrd,nsegsu
      real         ptlow,yradpi,yradpx
      integer      ioclude,iocluin,ioquen,kigrid
      common/drop7/ptlow,yradpi,yradpx,ioclude,iocluin,ioquen,kigrid
      real         fzcell,ptupp,qufac,fludiq,quexpo,fxcell
      common/drop8/fzcell,ptupp,qufac,fludiq,quexpo,fxcell
      real         avgdens,avgmass,ratiodens,ratiomass
      common/drop9/avgdens,avgmass,ratiodens,ratiomass
      real         facmicro
      common/dropa2/facmicro 
      integer      iospec,iocova,iopair,iozero,ioflac,iomom
      common/metr1/iospec,iocova,iopair,iozero,ioflac,iomom
      integer      nadd,iograc,iocite,ioceau,iociau
      common/metr2/nadd,iograc,iocite,ioceau,iociau
      integer      iomodl,idproj,idtarg
      real         wexcit
      common/hadr2/iomodl,idproj,idtarg,wexcit
      real          rexdifi,rexndii
      integer       idprojin,idtargin,irdmpr,isoproj,isotarg
      common/hadr25/idprojin,idtargin,rexdifi(4),rexndii(4),irdmpr,
     *              isoproj,isotarg
      integer      iostat,ioinco,ionlat,ioobsv,iosngl,iorejz,iompar
      common/metr3/iostat,ioinco,ionlat,ioobsv,iosngl,iorejz,iompar
      integer      ioinfl,ioinct,iowidn
      real         epsgc
      common/metr4/ioinfl,ioinct,iowidn,epsgc
      real         prob
      integer      nstmax,icbac,icfor
      common/lept2/nstmax,prob(99),icbac(99,2),icfor(99,2)
      integer      iolept,igampr,idisco
      common/lept3/iolept,igampr,idisco
      real        engmin,engmax
      integer     noebin,nrebin,iologe,iologl
      common/ebin/noebin,engmin,engmax,nrebin,iologe,iologl
      real         pi,pii,hquer,prom,piom,ainfin 
      common/cnsta/pi,pii,hquer,prom,piom,ainfin 
      integer      iversn,iverso 
      common/versn/iversn,iverso 
      integer      imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl
      common/accum/imsg,jerr(mxjerr),ntevt,nrevt,naevt,nrstr,nrptl
      integer       nglacc
      common/accum2/nglacc
      integer      nptlu,nrclu
      common/cptlu/nptlu
      common/cnrclu/nrclu
      real         tecm,volu
      common/drop6/tecm,volu
      integer      iterma,iternc,iterpr,iterpl,iozinc,iozevt
      common/metr5/iterma,iternc,iterpr,iterpl,iozinc,iozevt
      real         epsr
      integer      keepr 
      common/metr6/epsr,keepr 
      integer      keu,ked,kes,kec,keb,ket
      common/drop5/keu,ked,kes,kec,keb,ket
      integer      ikolmn,ikolmx,nglmin,nglmax
      real         segmin,segmax
      common/musct/ikolmn,ikolmx,nglmin,nglmax,segmin,segmax
      real          ptrmin
      common/highpt/ptrmin
      integer       itrigg,ioecc
      real          valecc
      common/ctrigg/itrigg,ioecc,valecc
      double precision seedi,seedj,seedj2,seedc
      integer      iseqini,iseqsim
      common/cseed/seedi,seedj,seedj2,seedc,iseqini,iseqsim  
      real          clust
      common/cjintc/clust(mxtau,mxvol,mxeps)
      real          volsum,vo2sum
      integer       nclsum
      common/cjintd/volsum(mxtau),vo2sum(mxtau),nclsum(mxtau)
      integer       iutotc,iutote
      common/ciutot/iutotc,iutote
      integer      nopen,nopenr
      common/copen/nopen,nopenr
      integer      kchopen,khiopen,kdtopen,kcpopen,klgopen,knxopen
      common/kopen/kchopen,khiopen,kdtopen,kcpopen,klgopen,knxopen
      character*6  xvaria,yvaria
      real         xminim,xmaxim,hisfac
      integer      normal,nrbins
      common/vana1/xvaria,yvaria,normal,xminim,xmaxim,nrbins,hisfac
      integer      iologb,iocnxb
      common/vana3/iologb,iocnxb
      integer      mxnody,nrnody,nody
      parameter(mxnody=200)
      common/nodcy/nrnody,nody(mxnody)
      real         ctaumin
      common/ctdcy/ctaumin
      integer ifoele
      common/cifoele/ifoele
      character*20 subpri
      integer      nrpri,ishpri 
      common/prnt2/nrpri,subpri(mxpri),ishpri(mxpri) 
      integer      ishevt,ixtau,iwseed,jwseed,ixgeometry
      common/prnt3/ishevt,ixtau,iwseed,jwseed,ixgeometry
      integer      jprint
      common/prnt4/jprint
      integer      mxcnt,ionoerr
      real         ar,ary,ardy,histowy
      parameter (mxcnt=20)
      common/vana4/ar(mxbins,5),ary(mxbins,mxcnt),ardy(mxbins,mxcnt)
     *,histowy(mxcnt),ionoerr
      real xpar1,xpar2,xpar3,xpar4,xpar5,xpar6,xpar7,xpar8
     *,xpar9,xpar10,xpar11,xpar12,xpar13,xpar14,xpar15,xpar16,xpar17
     *,xpar98,xpar99
      common/xpars/xpar1,xpar2,xpar3,xpar4,xpar5,xpar6,xpar7,xpar8
     *,xpar9,xpar10,xpar11,xpar12,xpar13,xpar14,xpar15,xpar16,xpar17
     *,xpar98,xpar99
      integer      khisto
      common/khist/khisto
      integer      nctcor,ncttim
      common/ctcor/nctcor
      common/ccttim/ncttim
      integer      kdensi
      real         tauv
      common/densi/kdensi(matau,nzeta,mxcoox,mxcooy),tauv(matau)
      integer       iorsce,iorsdf,iorshh,ionudi,kexit
      common/cjinti/iorsce,iorsdf,iorshh,ionudi,kexit
      real         amimfs,amimel
      common/camim/amimfs,amimel
      real          scr,scs,hacore
      common/craddf/scr,scs,hacore
      integer      iokoll
      common/ckoll/iokoll
      integer      ncnt,inicnt,nemsi
      common/cncnt/ncnt
      common/cicnt/inicnt
      common/cnemsi/nemsi
      real        gfactor,gwidth
      integer     iemspl,iemsct
      common/ems1/iemspl,iemsct,gfactor,gwidth
      real           amproj,amtarg,ypjtl,yhaha,pnullx
      common/chadron/amproj,amtarg,ypjtl,yhaha,pnullx
      real         xshift,etacut
      common/vana5/xshift,etacut
      double precision rnucl
      real         bnucl,xbtot
      integer      ixbDens
      common/nucl5/rnucl(mxnucl,2),bnucl(mxnucl,4),xbtot(4),ixbDens
      integer      infragm,ibreit
      common/nucl6/infragm,ibreit
      real         drnucl,rnuclo
      integer      nrnucl
      common/nucl4/nrnucl(2),drnucl(2),rnuclo(mxnucl,2)
      real         rapcms
      common/nucl7/rapcms
      real       xsig,xpom
      common/sig/xsig(7),xpom(7)
      integer      ktnbod
      common/metr7/ktnbod
      integer      iregge,isopom,ishpom,iscreen,nprmax,inueff,irmdrop
      common/hadr3/iregge,isopom,ishpom,iscreen,nprmax,inueff,irmdrop
      real         sigtot,sigcut,sigela,sloela,sigsd,sigine,sigdif
     *,sigineaa,sigtotaa,sigelaaa,sigcutaa,sigdd,sigcd
      common/hadr5/sigtot,sigcut,sigela,sloela,sigsd,sigine,sigdif
     *,sigineaa,sigtotaa,sigelaaa,sigcutaa,sigdd,sigcd
      integer      intpol,isigma,iomega,isetcs
      common/hadr6/intpol,isigma,iomega,isetcs
      real         alppom,slopom,gamhad,r2had,chad,wdiff
     *,gamtil,facdif,facmc,r2hads,gamhads,slopoms
      integer      isplit,iq2sat
      common/hadr4/alppom,slopom,gamhad(4),r2had(4),chad(4),wdiff(4)
     & ,gamtil,facdif,facmc,r2hads(4),gamhads(4),slopoms,isplit,iq2sat
      real          gamhadsi
      common/hadr42/gamhadsi(4)
      real         alpreg,sloreg,gamreg,r2reg,ptdiff,ptsend,xminremn
     &,xmindiff,ptsecu
      common/hadr7/alpreg,sloreg,gamreg,r2reg,ptdiff,ptsend,xminremn
     &,xmindiff,ptsecu
      real         alpqua,alppar,alpsea,alpval,alpdiq,alplea,alpdif
      common/hadr8/alpqua,alppar,alpsea,alpval,alpdiq,alplea(4),alpdif
      real          alpndi,alpdi,ptsendi,zdrinc,zmsinc,ptsems
      integer       irzptn
      common/hadr14/alpndi(2),alpdi(2),ptsendi,zdrinc,zmsinc,ptsems
     &             ,irzptn
      real          zbcut,zopinc,zipinc,zoeinc,xmxrem
      common/hadr15/zbcut,zopinc,zipinc,zoeinc,xmxrem
      real          fkainc,fkamax,zodinc,zbrmax,zdfinc,xzcut,ptvpom
      common/hadr16/fkainc,fkamax,zodinc,zbrmax(4),zdfinc,xzcut,ptvpom
      real          edmaxi,epmaxi
      common/hadr17/edmaxi,epmaxi
      real          ammsqq,ammsqd,ammsdd,cumpom,rexndi,rexdif
     *             ,reminv,rexpdif,rexres,zrminc,rexndf
      common/hadr9/ammsqq,ammsqd,ammsdd,cumpom,rexndi(4),rexdif(4)
     *             ,reminv,rexpdif(4),rexres(4),zrminc,rexndf
      integer      iclpro,icltar,iclegy
      common/had10/iclpro,icltar,iclegy
      integer      iclpro1,iclpro2,icltar1,icltar2,iclegy1,iclegy2
      common/had11/iclpro1,iclpro2,icltar1,icltar2,iclegy1,iclegy2
      real         egylow,egyfac
      common/had12/egylow,egyfac
      real         amdrmax,amdrmin,alpdro
      common/had13/amdrmax,amdrmin,alpdro(3)
      real         alpcoso,alpcose,betcoso,betcose
      common/had14/alpcoso,alpcose,betcoso,betcose
      real         accept,reject
      common/emsx1/accept,reject
      integer      iemspr,iemspm,iemspx,iemsrx,iemspu,iemsi2,iemspbx
      common/xems1/iemspr,iemspm,iemspx,iemsrx,iemspu,iemsi2,iemspbx
      integer      iemsse,iemsi1,iemsb,iemsbg,ioems,iemsdr
      common/xems2/iemsse,iemsi1,iemsb,iemsbg,ioems,iemsdr
      integer        ispacetime
      common/xspatim/ispacetime
      character*500 hisarg
      integer        ihisarg
      common/chisarg/ihisarg,hisarg(2,mxhisarg)
      real           difnuc,radnuc
      common /psar10/difnuc(mamxx),radnuc(mamxx)
      real diproj,ditarg,disize
      common/codipole/diproj(mamxx,2),ditarg(mamxx,2),disize
      integer        mxbarray,nbarray
      real           barray
      parameter (mxbarray=100)
      common/cbarray/barray(mxbarray),nbarray
      real          airznxs,airanxs,airwnxs
     *             ,airavznxs,airavanxs
      common/nxsair/airznxs(3),airanxs(3),airwnxs(3)
     *             ,airavznxs,airavanxs
      real            qgsincs
      common/mod2incs/qgsincs
      real            gheincs
      common/mod3incs/gheincs
      real            pytincs
      common/mod4incs/pytincs
      real            hijincs
      common/mod5incs/hijincs
      real            sibincs
      common/mod6incs/sibincs
      real            qgsIIincs
      common/mod7incs/qgsIIincs
      real            phoincs
      common/mod8incs/phoincs
      real            fluincs
      common/mod9incs/fluincs
      real            urqincs
      common/mod10incs/urqincs
      real            dpmincs
      common/mod12incs/dpmincs
      real            qgsIIIincs
      common/mod13incs/qgsIIIincs
      real           antot,ansh,ansf,pp4max,pp4ini,andropl,anstrg0
     *,anshf,ansff,antotf
     *,anstrg1,anreso0,anreso1,anghadr,antotre
      common/testpom/antot,ansh,ansf,pp4max,pp4ini,andropl,anstrg0
     *,anshf,ansff,antotf
     *,anstrg1,anreso0,anreso1,anghadr,antotre
      real         anintdiff,anintsdif,anintine,anintddif,anintcdif
     *,sigineex,sigdifex,sigsdex,sigddex,sigcdex
      integer      nglndif,nglhard
      common/cdiff/anintdiff,anintsdif,anintine,anintddif,anintcdif
     *,sigineex,sigdifex,sigsdex,sigddex,sigcdex,nglndif,nglhard
      real xnpomcol,xnhardcol,xnsatcol,xfrachard
      common/cnhard/xnpomcol,xnhardcol,xnsatcol,xfrachard
      real           epszero,alpff,betff
      common/cepszer/epszero,alpff(nclha),betff(2)
      real        tgss,wgss
      common/cgss/tgss(7,7),wgss(7,7)
      integer       iopcnt,jcentrality
      common/copcnt/iopcnt,jcentrality
      integer       isyst
      common/cisyst/isyst
      integer            irootcproot,iboein
      common/crootcproot/irootcproot,iboein
      integer      ihepmc,ihepframe
      common/chepmc/ihepmc,ihepframe
      integer          ihyskip
      common /cihyskip/ihyskip

      integer    maxzzsoi
      parameter (maxzzsoi=1)
      real        zzsoival
      character*2 zzsoich
      common/csoi/zzsoival(maxzzsoi),zzsoich(maxzzsoi)
      integer         istatom
      common/cistatom/istatom

      real alpDs,alpDps,alpDpps,betDs,betDps,betDpps,gamDs,delDs
     *    ,alpD,alpDp,alpDpp,betD,betDp,betDpp,gamD,delD
      common/Dparams/alpDs(  idxD0:idxD, nclegy, nclha,nclha),
     *              alpDps( idxD0:idxD, nclegy, nclha,nclha),
     *              alpDpps(idxD0:idxD, nclegy, nclha,nclha),
     *              betDs(  idxD0:idxD, nclegy, nclha,nclha),
     *              betDps( idxD0:idxD, nclegy, nclha,nclha),
     *              betDpps(idxD0:idxD, nclegy, nclha,nclha),
     *              gamDs(  idxD0:idxD, nclegy, nclha,nclha),
     *              delDs(  idxD0:idxD, nclegy, nclha,nclha)
      real          bmxdif,bkmxndif
      integer       idxDmin
      common/Dparam/alpD(  idxD0:idxD1, nclha, nclha),
     *              alpDp( idxD0:idxD1, nclha, nclha),
     *              alpDpp(idxD0:idxD1, nclha, nclha),
     *              betD(  idxD0:idxD1, nclha, nclha),
     *              betDp( idxD0:idxD1, nclha, nclha),
     *              betDpp(idxD0:idxD1, nclha, nclha),
     *              gamD(  idxD0:idxD1, nclha, nclha),
     *              delD(  idxD0:idxD1, nclha, nclha),
     *     idxDmin,bmxdif(nclha, nclha),bkmxndif
      double precision alpUni,betUni,betpUni,fpUni,ftUni,betfom
     *                 ,epspUni,epstUni,zzpUni,zztUni
      real           alpfom,alpfomi,gamfom
      common/DparUni/alpUni(  idxD0:idxD1,2),
     *               betUni(  idxD0:idxD1,2),
     *               betpUni(idxD0:idxD1,2),fpUni,ftUni,
     *               epspUni(idxD0:idxD1),zzpUni,
     *               epstUni(idxD0:idxD1),zztUni,
     *               betfom,alpfom,alpfomi,gamfom
      integer      idlead,ilprtg,mfragmaxp,mfragmaxt
      common/crvar/idlead,ilprtg,mfragmaxp,mfragmaxt
      integer        iLHC,ipytune
      common/LHCtune/iLHC,ipytune
      logical chargex
      integer ichargexin,ichargex,kchargex
      real xschargex(0:3),rexchrg(4)
      double precision pnullx0
      common/chrgexch/pnullx0,xschargex,rexchrg,ichargex,kchargex
     *               ,ichargexin,chargex

      integer nrcode,nidtmax
      parameter (nrcode=5,nidtmax=600)
      integer idtbl,ifl1tbl,mutbl
      integer ifl2tbl,ifl3tbl,jspintbl
      integer nlidtbl,nidtmxx,nidtbl,idegtbl
      real chrgtbl,amtbl,witbl
      character*50 nametbl
      character*1 statbl
      common/cidt/idtbl(nrcode,nidtmax),ifl1tbl(nidtmax),
     *ifl2tbl(nidtmax),ifl3tbl(nidtmax),jspintbl(nidtmax),
     *chrgtbl(nidtmax),amtbl(nidtmax),witbl(nidtmax),
     *nlidtbl(-9990:9900),nametbl(nidtmax),nidtmxx(nrcode),nidtbl,
     *mutbl(nidtmax),idegtbl(nidtmax),statbl(nidtmax)

      real qkonia1,qkonia2,qkonia3,qkonia4
      common/cqkonia/qkonia1,qkonia2,qkonia3,qkonia4

      integer mxdefine,ndefine,l1define,l2define
      character w1define*100,w2define*100
      parameter(mxdefine=500)
      common/cdefine/ndefine,l1define(mxdefine),l2define(mxdefine)
     &               ,w1define(mxdefine),w2define(mxdefine)

      integer k1domax, k2domax
      parameter(k1domax=12, k2domax=12)
      character *160 dooptnsargtext(k1domax,k2domax)
      common/cdooptns/dooptnsargtext

      integer ibeginwrite,jselect,kselect
      common/cbeginwrite/ibeginwrite,jselect,kselect

      integer macronr, macromax, macrolines, macroargmax
      parameter (macromax=2, macrolines=400,macroargmax=20)
      integer macroargs(macromax)
      character *160 macroname(macromax)
      character *160 macrotext(macromax,macrolines)
      integer macroargtextlength
      parameter (macroargtextlength=30)
      character *30 macroargtext(macromax,macroargmax)
      common/cmacros/macronr,macroname,macrotext,macroargs,macroargtext
      integer macrotextlines(macromax)
      common/cmacros2/macrotextlines
      integer macrocurrent,macrocurrentline
      common/cmacros3/macrocurrent,macrocurrentline

      integer mxtxt80,ntxt80
      parameter (mxtxt80=10)
      character*80 txt80(mxtxt80)
      common/ctxt80/ntxt80,txt80 


