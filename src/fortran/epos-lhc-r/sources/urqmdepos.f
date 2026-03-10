C
C  This file (originaly KW/u.f) is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

!--------------------------------------------------------------------------------
! structure of /cptl/ after calling hacas:
!
!     --------------
!      1           \    (A)  ptls before hacas 
!      ...          >           ist = 0 -> 3 when going into UrQMD
!      nptlbur     /            ist = 0 -> 4 when unknown 
!     --------------
!      1           \    (B)  return from hacas,  ist = 0
!      ...          >        copy ist=4 particles -> ist = 0
!      nptl        /        then decayall(1,99),     ist=0,1
!     --------------
!      ....        \    (C) decayall(1,3)
!      ...          >          copies (A),        ist = 8 (last gen)
!      nptl        /           decays,            ist = 6,8
!     --------------
!      ....        \    (D)  utrescxx
!      ...          >      copies (C),  corrects
!      nptl        /              E,p conservation   ist = -2,-1
!     --------------
!
!  nptl is at each stage the maximum particle index
!-------------------------------------------------------------------------------


! 2012 July: new def of ityptl in SR printReaction
!            acticate string suppression via CTOption(12), later comment this
! 2012 Aug:  option to suppress strings via "set iunostring 1"
!            option to skip channels via "hacas chaskip"
!              precise definitions in make22 in "if(iuchaskip.ge.1)"  block
!            print options activated via lpr=.true.
!              use getbranxx instead of getbran to activate prints
! 2012 Sep: unactivate all the above options, unless skip io=10

      subroutine hacas(n,iret)
      include "epos.inc"
      character*500 edir
      integer kcollmax
      common/ckcollmax/kcollmax
      common/cdir/edir
      common/cishuuu/ishuuu
      ishuuu=ish
      call utpri('hacas ',ish,ishini,4)
      edir=fnnx(1:nfnnx)//' '
      kcollmax=0

      call checkengy('before hacas     ')
      call uinitial(np)
      if(np.gt.0)then
        call uepos
        call uexit
        call utrescxx(iret,0)   !before decay because decay conserve energy and not to spoil invariant masses
        if(iret.ne.0)then
        if(ish.ge.1)call utmsg("Pb in utrescxx in uexit, redo event !&")
          return
        endif
      endif
c     Form nuclear fragments
      if(model.eq.1.or.model.eq.12
     *     .and.(maproj.gt.1.or.matarg.gt.1))then
        call emsfrag(3,iret)
        if(iret.ne.0)then
          call utmsg("Pb in emsfrag in uexit, redo event !&")
          return
        endif
      endif
      call decayall(1,99)
      if(ihacas.gt.0.and.np.gt.0)call decayall(1,3)
      call xf15
      if(n.eq.nevent)call uplot
      call checkengy('after hacas     ')
      if(ish.ge.2)then
        write(ifch,'(a,i8,i3)')
     .'end of hacas - nptl icccptl = ',nptl,icccptl
c      call checktime('exit hacas;') 
        write(ifch,'(a,i10)')'maximum collision number',kcollmax
      endif
      call utprix('hacas ',ish,ishini,4)
      end

  !---------------------------------------------------------------------
  ! Resonance special codes, which allow to indentify resonances which 
  ! decay via particular channels (like phis which decay into KK)
  ! in the cascade part and which are not in the usual particle list. 
  ! id and ity is needed. 
  ! ity gives information whether the daughters rescatter or not.
  !    id:  33100 ... phi -> K+ K- (BR 0.5)
  !         33101 ... phi -> K0 K0bar (BR 0.35)
  !         23100 ... K*0 -> pi- K+ (BR 0.667)  (includes K*0bar)
  !         23101 ... K*0 -> pi0 K0 (BR 0.333)  (includes K*0bar)
  !         etc (see printReaction)
  !    ity: 80 ...  no rescattering  (KK detectable)
  !         81 ...  rescattering   (not KK detectable)
  !    istptl = 9 (not really needed, id is unique)
  !   Am ende phi yield und spektra durch 0.85 teilen,
  !   das ist der branching ratio in K+K
  !---------------------------------------------------------------------
  !   phis decaying via KK
  !
  !Technical info:
  !
  !Es werden bei jedem phi zerfall in K+K die phi informationen gespeichert,
  !d.h. 4er imuls, ort und masse:
  !
  !   phip0(phinum),phipx(phinum),phipy(phinum),phipz(phinum)
  !   phir0(phinum),phirx(phinum),phiry(phinum),phirz(phinum)
  !   phimass(phinum)
  !
  !   phiscatt(phinum) :  0 wenn keines der Kaonen rescattert
  !                       1 falls mind. eines rescattert.
  !
  !   phinum : Index phi Liste
  !
  !Wird in addpart.f und delpart.f hoffentlich korrekt weitergegeben.
  !Wenn nicht sollten wir das merken wenn pl?tzlich teilchen ausser Kaonen ein phinum
  !.ne. 0 haben (wird abgefragt).
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  !es werden keine schwachen zerfaelle gemacht. (d.h. Lambdas sind
  !stabil, pionen sind stabil etc., Sigma's sind uebrigens auch stabil,
  !etas ebenfalls).
  !die vollstaendige liste aller erlaubten zerfaelle findest sich im file
  !blockres.f
  !---------------------------------------------------------------------
  !die freeze-out koordinaten haben ein fr davor
  !also z.b. frr0(i), frrx(i), frpx(i),...)
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  ! in u directory remove input.f
  ! in u directory remove tabinit.f
  ! in u directory remove gnuranf.f + replace ranf by ranff
  !---------------------------------------------------------------------
  ! modifs of routines:
  !   upmerge.f:
  !      comment the line "write(ifmt,*) '(Info) pdg2ityp: ..."
  !   output.f:
  !      add in subroutine output(iunit):
  !         integer ninit
  !         common /cninit/ninit
  !      add in entry file15out
  !        ninit=nin
  !   make22.f, SUBROUTINE STREXCT
  !        skip "go to 100" when ntry > 1000
  !        to avoid endless loop
  !   scatter.f, subroutine scatFinal:
  !             if(NBlColl.lt.1000)then                           !***add***
  !             do 306 i=1,nb1+nb2
  !                   ...
  !        306  continue
  !             else                                                !***add***
  !             print*,"Too many attempts => ignore Pauli blocking" !***add***
  !             NBlColl=0                                          !***add***
  !             endif                                               !***add***
  !
  !
  !
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  !  uinitial ... transfer from EPOS
  !  uepos ... main u routine
  !  uexit ... transfer to EPOS
  !  file14outx .. analysis, essentially copied from file14out (output.f)
  !---------------------------------------------------------------------
  ! crash because of too small resonance mass; message:
  !        write(ifmt,*)'anndex(dec): no final state found for:',i1,m1,iz1
  !        write(ifmt,*)'please check minimal masses: m1,m1min,m2min'
  !        write(ifmt,*)'and iso3 of decaying particle'
  !        write(ifmt,*)(prob(j),j=0,maxbr)
  !        stop
  ! avoided by calling "checkdecay", and using epos mass in case of problem
  !---------------------------------------------------------------------
  ! crash because of unknown particle type:
  !      idpdg = 40323 for example
  ! avoided by ignoring (idepos=0)
  !---------------------------------------------------------------------
  ! "odd" times (like 3.96) used instead of 4. Any good reason?????
  !---------------------------------------------------------------------

c#####################################################################
c#####################################################################

                        subroutine uinitial(np)

c#####################################################################
c#####################################################################

      implicit none
c     debug and validity range
      integer nmax, nspl
      real*8 hit_sphere
      parameter (nmax = 100000) ! maximum number of particles
      parameter (nspl = 500)  ! dimension of spline arrays
      parameter (hit_sphere = 8.d0)  ! hard collision cutoff: 251 mbarn
      integer Ap, At, Zp, Zt, npart, nbar, nmes, ctag
      integer nsteps,ranseed,event,eos,dectag,uid_cnt
      integer NHardRes,NSoftRes,NDecRes,NElColl,NBlColl
      real*8  time,  acttime, bdist, ebeam, bimp,bmin,ecm
      logical success
c 7 integer
      common /sys/ npart, nbar, nmes, ctag,nsteps,uid_cnt,
     +             ranseed,event,Ap,At,Zp,Zt,eos,dectag,
     +             NHardRes,NSoftRes,NDecRes,NElColl,NBlColl
     +             ,success
      common /rsys/ time,acttime,bdist,bimp,bmin,ebeam,ecm
      real*8
     +     gw, sgw, delr, fdel, dt,
     +     da, db,
     +     Cb0, Yuk0, Pau0, Sky20, Sky30, gamSky, gamYuk, drPau, dpPau,
     +     dtimestep
c 19 real*8

      real*8 cutmax, cutPau, cutCb, cutYuk, cutSky, cutdww
      common /cuts/ cutmax, cutPau, cutCb, cutYuk, cutSky, cutdww
      real*8 spx(nspl), spPauy(nspl), outPau(nspl),
     +                spCby(nspl),  outCb(nspl),
     +                spYuky(nspl), outYuk(nspl),
     +                spSkyy(nspl), outSky(nspl),
     +                spdwwy(nspl), outdww(nspl)
      common /spdata/ spx, spPauy, outPau, spCby,  outCb,
     +                     spYuky, outYuk, spSkyy, outSky,
     +                     spdwwy, outdww
      real*8
     +     r0(nmax), rx(nmax), ry(nmax), rz(nmax),
     +     p0(nmax), px(nmax), py(nmax), pz(nmax),
     +     airx(nmax), airy(nmax), airz(nmax),
     +     aipx(nmax), aipy(nmax), aipz(nmax),
     +     aorx(nmax,4), aory(nmax,4), aorz(nmax,4),
     +     aopx(nmax,4), aopy(nmax,4), aopz(nmax,4),
     +     fmass(nmax), rww(nmax),
     +     dectime(nmax), tform(nmax), xtotfac(nmax)


      integer spin(nmax),ncoll(nmax),charge(nmax),strid(nmax),
     +        ityp(nmax),lstcoll(nmax),iso3(nmax),origin(nmax),uid(nmax)
      common/isys/spin,ncoll,charge,ityp,lstcoll,iso3,origin,strid,
     +            uid

      common /coor/ r0, rx, ry, rz, p0, px, py, pz, fmass, rww, dectime
      common /frag/ tform, xtotfac
      common /aios/ airx, airy, airz, aipx, aipy, aipz,
     +              aorx, aory, aorz, aopx, aopy, aopz
      common /pots/ Cb0, Yuk0, Pau0, Sky20, Sky30, gamSky,
     +              gamYuk, drPau, dpPau, gw, sgw, delr, fdel,
     +              dt,da, db,dtimestep
c spectator arrays
        integer smax
        parameter(smax=500)  ! maximum number of spectators
        real*8 r0s(smax), rxs(smax), rys(smax), rzs(smax),
     +         p0s(smax), pxs(smax), pys(smax), pzs(smax),
     +         sfmass(smax)

        integer sspin(smax), scharge(smax), sityp(smax), siso3(smax),
     +          suid(smax)
        integer nspec
        common /scoor/ r0s, rxs, rys, rzs, p0s, pxs ,pys, pzs, sfmass
        common /sisys/ sspin, scharge, sityp, siso3, suid
        common /ssys/ nspec
        real*8 p0td(2,nmax),pxtd(2,nmax),pytd(2,nmax),pztd(2,nmax),
     +         fmasstd(2,nmax)
        integer ityptd(2,nmax),iso3td(2,nmax)
        integer itypt(2),uidt(2),origint(2),iso3t(2)
        common /rtdelay/p0td,pxtd,pytd,pztd,fmasstd
        common /itdelay/ityptd,iso3td
        common /svinfo/itypt,uidt,origint,iso3t
        real*8 ffermpx(nmax), ffermpy(nmax), ffermpz(nmax)
        real*8 peq1, peq2
        common /ffermi/ ffermpx, ffermpy, ffermpz
        common /peq/ peq1,peq2
        integer numcto,numctp,maxstables
        parameter(numcto=400) ! maximum number of options
        parameter(numctp=400) ! maximum number of parameters
        parameter(maxstables=20) ! maximum number of stable particles
        integer   CTOption(numcto)
        real*8    CTParam(numctp)
        common /options/CTOption,CTParam
        real*8 frr0(nmax), frrx(nmax), frry(nmax), frrz(nmax),
     +     frp0(nmax), frpx(nmax), frpy(nmax), frpz(nmax)
      common /frcoor/ frr0, frrx, frry, frrz, frp0, frpx, frpy, frpz

      include 'urqmd34/comres.f'

      !----------------------------------------------------------------------------
      !epos common blocks for particle list
      !----------------------------------------------------------------------------
      integer      ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      common/files/ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      integer iappl,model
      common/appli/iappl,model
      integer mmry,mxptl
      parameter (mmry=1)   !memory saving factor
      parameter (mxptl=300000/mmry) !max nr of particles in epos ptl list
      integer iorptl(mxptl),idptl(mxptl),istptl(mxptl),
     *  ifrptl(2,mxptl),jorptl(mxptl),ibptl(4,mxptl)
     * ,ityptl(mxptl)
      real pptl(5,mxptl),tivptl(2,mxptl),xorptl(4,mxptl)
      common/cptl/nptl,pptl,iorptl,idptl
     *,istptl,tivptl,ifrptl,jorptl
     *,xorptl,ibptl,ityptl
      integer nptl
      real         qsqptl,zpaptl
      common/c8ptl/qsqptl(mxptl),zpaptl(2,mxptl)
      real        phievt,bimevt,pmxevt,egyevt
     *,xbjevt,qsqevt,zppevt,zptevt
      integer     nevt,kolevt,koievt,kohevt,npjevt
     *,ntgevt,npnevt,nppevt,ntnevt,ntpevt,jpnevt,jppevt,jtnevt,jtpevt
     *,nglevt,minfra,maxfra,npglb,ntglb
      common/cevt/phievt,nevt,bimevt,kolevt,koievt,pmxevt,egyevt,npjevt
     *,ntgevt,npnevt,nppevt,ntnevt,ntpevt,jpnevt,jppevt,jtnevt,jtpevt
     *,xbjevt,qsqevt,nglevt,zppevt,zptevt,minfra,maxfra,kohevt
     *,npglb,ntglb
      integer laproj,maproj,latarg,matarg,nptlbd,nptlpt
      real core,fctrmx
      common/nucl1/laproj,maproj,latarg,matarg,core,fctrmx
      common/c4ptl/nptlbd,nptlpt
      integer istore,istmax,irescl,ntrymx,nclean,iopdg,ioidch
      real gaumx
      common/othe1/istore,istmax,gaumx,irescl,ntrymx,nclean,iopdg,ioidch
      integer      iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi
      common/prnt1/iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi
      integer ihacas
      common/chacas/ihacas
      integer icinpu
      real engy,elepti,elepto,angmue
      common/lept1/engy,elepti,elepto,angmue,icinpu
      integer mxjerr,imsg,jerr,ntevt,nrevt,naevt,nrstr,nrptl
      parameter (mxjerr=10)
      common/accum/imsg,jerr(mxjerr),ntevt,nrevt,naevt,nrstr,nrptl
      real rminfrg,emaxfrg,facgrey,p3grey
      common /cfragm/rminfrg,emaxfrg,facgrey,p3grey
      !---------------------------------------------------------------
      integer i,nn,idtmp,ityptmp,iso3tmp,itmp,np
      integer idtrafo
      external idtrafo
      integer fchg
      external fchg
      real*8 dectim
      external dectim
      real*8 mintime,eb
      integer j,k,icount,npold
      integer strcount
      common /inewpartxx/ strcount
      real*8 etot
      integer igeteost,igethyt,ireadhyt
      common/hydr2/igeteost,igethyt,ireadhyt
      integer nrap,ntau,nrad,id
      real*8 dcoll(-10:10,-1:40),rcoll(0:20,-1:40),zevents
      real*8 dpart(-10:10,-1:40)
      common /ccoll/dcoll,rcoll,dpart,zevents
      integer nupa(0:10)
      common/cnupa/nupa
      real*8 deluutau
      integer nuutaumx
      common /cf15/deluutau,nuutaumx
      integer nui
      data nui/0/
      save nui
      !integer nptest,ia,jca(3),nflav
      !parameter (nflav=6)
      !integer jc(nflav,2),ic(2)
      integer itypart(nmax), itrace(2,nmax), iorpart(nmax)
      common /city/itypart
      common /ctrace/itrace
      common /corpart/iorpart
      integer kjpsi
      real a,b,c,rangen
      external rangen
      integer idyy,iutest

      real*8 betatx,betaty,betatz,aamt,rrap
      integer phiid(nmax)
      common /phitrack/ phiid

      integer nptlbur
      common/cnptlbur/nptlbur

      integer maxphi
      parameter(maxphi=140000)
      real*8 phip0(maxphi),phipx(maxphi),phipy(maxphi),phipz(maxphi)
      real*8 phir0(maxphi),phirx(maxphi),phiry(maxphi),phirz(maxphi)
      real*8 phimass(maxphi)
      integer phinum,phiscatt(maxphi),phiidres(maxphi)

      common /phidata/ phip0,phipx,phipy,phipz,
     $      phir0,phirx,phiry,phirz,phimass,
     $      phinum,phiscatt,phiidres

      integer ishini,ishuuu
      common/cishuuu/ishuuu

c#####################################################################
c##############################     subroutine uinitial
c#####################################################################

      call utpri('uiniti',ishuuu,ishini,4)

      if(ishuuu.ge.2)write(ifch,'(a)')'******** start uinitial ********'

      np=0
      iutest=1  !0 or 1

      nui=nui+1
      if(nui.eq.1)then
        write(ifmt,'(/a,a)')' (info)',
     .     '   Hadronic rescattering in EPOS done using :'
c        if(ish.ge.1)write(ifch,'(a,50a1)')' (info)',('u',i=1,50)
        do j=-1,40
          do i=-10,10
            dcoll(i,j)=0d0
          enddo
          do i=-10,10
            dpart(i,j)=0d0
          enddo
          do i=0,20
            rcoll(i,j)=0
          enddo
        enddo
        zevents=0
        do i=0,10
        nupa(i)=0
        enddo
      endif
      zevents=zevents+1

      do nn=1,nmax
        itypart(nn)=0
        itrace(1,nn)=0
        itrace(2,nn)=0
        iorpart(nn)=0
      enddo

      call uinit(0)
      call osc_header
      call osc99_header

      npart = 0
      npold = 0
      nbar=0
      nmes=0
      uid_cnt=0
c reset counters
c     all collisions/decays
      ctag  = 0
c     all decays
      dectag = 0
c     number of prod. hard resonances
      NHardRes=0
c     number of prod. soft resonances
      NSoftRes=0
c     number of prod. resonances via decay
      NDecRes=0
c     number of blocked collisions
      NBlColl=0
c     number of elastic collisions
      NElColl=0
c     number of strings
      strcount=1
c
      eb=0D0
c icount is the number of EXTRAordinary pro/tar combinations (i.e. pion ...)
      icount = 0
c reset particle vectors
      do 20 j=1,nmax
        spin(j)  = 0
        ncoll(j) = 0
        lstcoll(j)=0
        r0(j) = 0.0
        rx(j)         = 0.0
        ry(j)         = 0.0
        rz(j)         = 0.0
        p0(j)         = 0.0
        px(j)         = 0.0
        py(j)         = 0.0
        pz(j)         = 0.0
        frr0(j) = 0.0
        frrx(j)    = 0.0
        frry(j)    = 0.0
        frrz(j)    = 0.0
        frp0(j)    = 0.0
        frpx(j)    = 0.0
        frpy(j)    = 0.0
        frpz(j)    = 0.0
        fmass(j) = 0.0
        charge(j)= 0
        iso3(j)  = 0
        ityp(j)  = 0
        dectime(j)= 0.0
        origin(j)=0
        tform(j)=0.0
        xtotfac(j)=1.0
        strid(j)=0
        uid(j)=0
         ffermpx(j) = 0.0
         ffermpy(j) = 0.0
         ffermpz(j) = 0.0

        phiid(j)=0

         do 21 k=1,2
            p0td(k,j)=0.d0
            pxtd(k,j)=0.d0
            pytd(k,j)=0.d0
            pztd(k,j)=0.d0
            fmasstd(k,j)=0.d0
            ityptd(k,j)=0
           iso3td(k,j)=0
 21         continue
 20   continue

       phinum=0


c epos event info to u event info
      bimp = bimevt

c initialise
      npart=0  !number of particles transferred to u
      mintime = 1d2 !the minimum formation time

c energy of spectators
      etot=0
      do nn=1,nptlpt
        if(istptl(nn).eq.0)etot=etot+pptl(4,nn)
      enddo

c ptls from testFile

      !-----------------------------------
      ! 1st run with optns:  set ihacas 2
      ! 2nd run with optns:  set ihacas 3
      !-----------------------------------
      if(ihacas.eq.3)then
        nptlpt=0
        call openOldTestFile
        call readTestFile(nn)
        nptlbd=nn
      endif

c decay unknown particles

c      write(ifch,*)' ~~~~~~~ decay'
      nptlbur=nptlbd
      !nw=nptlbd
      !nwi=nw+1
      do nn=nptlpt+1,nptlbd
        if(istptl(nn).eq.0)then
          call idcorrect(ityptl(nn),idptl(nn),idtmp)
          call pdg2id(ityptmp,iso3tmp,idtmp)
          if(abs(ityptmp).gt.1000)then  !unknown hadrons
            istptl(nn)=4
            !~~~~~~~~~~~~~~~~~~~~~~~~~~
            !nw=nw+1
            !!write(ifch,*)' ~~~~~~~ unknown hadron',idtmp
            !!.      ,ityptmp,iso3tmp,nn,nw
            !istptl(nw)   = istptl(nn)
            !idptl(nw)    = idptl(nn)
            !ityptl(nw)   = ityptl(nn)
            !do i=1,4
            !xorptl(i,nw) = xorptl(i,nn)
            !enddo
            !do i=1,5
            !pptl(i,nw)   = pptl(i,nn)
            !enddo
            !ibptl(1,nw)  = ibptl(1,nn)
            !iorptl(nw)   = nn
            !do i=1,2
            !tivptl(i,nw)  = tivptl(i,nn)
            !enddo
            !zpaptl(1,nw)=zpaptl(1,nn)
            !zpaptl(2,nw)=zpaptl(2,nn)
            !ifrptl(1,nn)=nw
            !ifrptl(2,nn)=nw
            !~~~~~~~~~~~~~~~~~~~~
          elseif(ityptl(nn).ge.70)then  !not to lose particles not counted later
            istptl(nn)=4
c          elseif(ityptl(nn)/10.eq.4)then
c              istptl(nn)=4
          elseif(ityptl(nn).eq.4)then
c           use only part of spectators to limit intra nuclear cascade (otherwise too many grey nucleons)
c            if(rangen().gt.facgrey*7)then   !factor 7 to use the same parameter with and without hadcas
              istptl(nn)=4   !if this is changed, please check utrescxxx !!!
c            endif
          elseif(ityptl(nn).eq.5)then
c           use only part of spectators to limit intra nuclear cascade (otherwise too many grey nucleons)
c            if(rangen().gt.facgrey*7)then
              istptl(nn)=4   !if this is changed, please check utrescxxx !!!
c            endif
          endif
        endif
      enddo
      !~~~~~~~~~~~~~~~~~~~~
      !if(nw.ge.nwi)then
      !  nptl=nw
      !  call decayall(nwi,999)
      !  nptlbd=nptl
      !  !call alist(' ~~~~~~~ list from unknown hadrons&',nwi,nptl)
      !  do nn=nwi,nptl
      !    if(istptl(nn).eq.0)then
      !      call idcorrect(ityptl(nn),idptl(nn),idtmp)
      !    endif
      !    if(istptl(nn).eq.1)istptl(nn)=4
      !  enddo
      !endif
      !~~~~~~~~~~~~~~~~~~~

c fill in the baryons first
      !write(ifmt,'(a/)')' '
      !nptest=0
      !do i=1,3
      !jca(i)=0
      !enddo
      if(ihacas.eq.2)call openNewTestFile
      if(ish.ge.2)call alist('(u) initial baryons&',0,0)
      if(ish.ge.2)write(ifch,*)'check between ',nptlpt+1,' and ',nptlbd
      if(ish.ge.2)write(ifch,*)' '
      nbar = 0
      do nn=nptlpt+1,nptlbd
        if(istptl(nn).eq.0.and.ityptl(nn).lt.70)then
          if(ihacas.eq.2)call writeParticleTestFile(nn)
          if(ityptl(nn).eq.61)then
           idtmp=idptl(nn)
          else
           idtmp=idtrafo('nxs','pdg',idptl(nn))
          endif
          if(abs(idtmp).lt.11.or.abs(idtmp).gt.16)then
            call pdg2id(ityptmp,iso3tmp,idtmp)
            if(abs(ityptmp).ge.1000)then
              jerr(10)=jerr(10)+1
              if(ish.ge.2)then
                write(ifch,*)' ~~~~~~~ ERROR: unknown hadron'
     .          ,ityptl(nn),idtmp,ityptmp,iso3tmp,nn
              endif
            endif
            if(abs(ityptmp).le.maxbar) then
              if(ish.ge.2)call alistc('&',nn,nn)
              !ia=abs(idptl(nn))
              !if(ia.eq.2130.or.ia.eq.1230.or.ia.eq.1130.or.ia.eq.2230
              !.         .or.ia.eq.1131.or.ia.eq.1231.or.ia.eq.2231)
              !.         !nptest=nptest+1
              !call idtr4(idptl(nn),ic)
              !call iddeco(ic,jc)
              !do i=1,3
              !jca(i)=jca(i)+jc(i,1)-jc(i,2)
              !enddo
              !if(jc(3,1)-jc(3,2).ne.0)write(ifmt,'(i6,$)')idptl(nn)
              etot=etot+pptl(4,nn)
              nbar=nbar+1
              !#########################
              if(iutest.eq.1)then
                c=pptl(4,nn)
                if(.not.(c.le.0..or.c.ge.0.))then
                  write(ifmt,*)'p=',(pptl(j,nn),j=1,5)
                  write(ifmt,*)'ERROR NaN catch uinitial (1)'
                  call utstop('ERROR NaN catch uinitial (1)&')
                endif
                if(pptl(4,nn).lt.2.)then
                 a=pptl(4,nn)**2-pptl(3,nn)**2
                 b=pptl(5,nn)**2+pptl(1,nn)**2+pptl(2,nn)**2
                 idyy=idptl(nn)
                 if(abs(a-b).gt.0.1)write(ifch,*)
     .           'mass shell pb  ini :',idyy,a,b,ityptl(nn),nn
                endif
              endif
              !#########################
              call DelayFormation(nn)
              r0(nbar)=xorptl(4,nn)
              rx(nbar)=xorptl(1,nn)
              ry(nbar)=xorptl(2,nn)
              rz(nbar)=xorptl(3,nn)
              p0(nbar)=pptl(4,nn)
              px(nbar)=pptl(1,nn)
              py(nbar)=pptl(2,nn)
              pz(nbar)=pptl(3,nn)
              fmass(nbar)=pptl(5,nn)
              aamt=sqrt(fmass(nbar)**2+px(nbar)**2+py(nbar)**2)
              rrap=sign(1d0,pz(nbar))*log((p0(nbar)+abs(pz(nbar)))/aamt)
              pz(nbar)=aamt*sinh(rrap)
              p0(nbar)=aamt*cosh(rrap)
              betatx=px(nbar)/p0(nbar)
              betaty=py(nbar)/p0(nbar)
              betatz=pz(nbar)/p0(nbar)
              !#########################
              if((betatx**2+betaty**2+betatz**2).gt.1.0d0)then
                write(*,*)'NaN ; beta>1 in uinitial (1)'
                write(*,*)'NaN   ',ityptl(nn),idptl(nn)
                write(*,*)'NaN   ', betatx,betaty,betatz
                write(*,*)'NaN   ',pptl(1,nn),pptl(2,nn),pptl(3,nn)
     .                            ,pptl(4,nn)
                write(*,*)'NaN   ',px(nbar),py(nbar),pz(nbar)
     .                             ,p0(nbar)
              endif
              !#########################
              ityp(nbar)=ityptmp
              iso3(nbar)=iso3tmp
              if(ityptl(nn).eq.61)then
              charge(nbar)=ibptl(1,nn)
              else
              charge(nbar)=fchg(iso3(nbar),ityp(nbar))
              endif
              itypart(nbar)=ityptl(nn)
              itrace(1,nbar)=nbar
              itrace(2,nbar)=nbar
              !?????????????????????????
              !if(px(nbar)**2+py(nbar)**2.gt.20.**2)then
              !write(ifch,'(a,i8,i3,f9.2,f6.2)')
              !.          ' HHHHHHH in  ',nbar,itypart(nbar),r0(nbar)
              !.        ,sqrt(px(nbar)**2+py(nbar)**2)
              !endif
              !????????????????????????????
              iorpart(nbar)=0
              lstcoll(nbar)=0
              ncoll(nbar)=0
              origin(nbar)=0
              tform(nbar)=r0(nbar)
              dectime(nbar)=dectim(nbar,1)+tform(nbar)
              xtotfac(nbar)=0d0
              if(r0(nbar).lt.mintime) mintime = r0(nbar)
              !print*,'uInput:',nbar,nn,idptl(nn),ityptl(nn)
c              write(ifmt,*)'uInput:',nbar,nn,idptl(nn),idtmp
c     *        ,rx(nbar),ry(nbar),rz(nbar),px(nbar),py(nbar),pz(nbar)
            endif
          endif
        endif
      enddo

c then fill in the mesons
      if(ish.ge.2)call alist('(u) initial mesons&',0,0)
      if(ish.ge.2)write(ifch,*)'check between ',nptlpt+1,' and ',nptlbd
      if(ish.ge.2)write(ifch,*)' '
      nmes = 0
      do nn=nptlpt+1,nptlbd
        if(istptl(nn).eq.0.and.ityptl(nn).lt.70)then
          if(ityptl(nn).eq.61)then
            idtmp=idptl(nn)
          else
            idtmp=idtrafo('nxs','pdg',idptl(nn))
          endif
          if(abs(idtmp).lt.11.or.abs(idtmp).gt.16)then
            call pdg2id(ityptmp,iso3tmp,idtmp)
            if(abs(ityptmp).ge.1000)then
              jerr(10)=jerr(10)+1
              if(ishuuu.ge.2)then
               write(ifch,*)' ~~~~~~~ ERROR: unknown hadron'
     .         ,ityptl(nn),idtmp,ityptmp,iso3tmp,nn
              endif
            endif
            if(abs(ityptmp).ge.minmes.and.abs(ityptmp).lt.1000) then
              !if(ishuuu.ge.3)call alistc('&',nn,nn)
              if(ish.ge.2)call alistc('&',nn,nn)
              nmes=nmes+1
              etot=etot+pptl(4,nn)
              itmp=nbar+nmes
              !#########################
              if(iutest.eq.1)then
                c=pptl(4,nn)
                if(.not.(c.le.0..or.c.ge.0.))then
                  write(ifmt,*)'p=',(pptl(j,nn),j=1,5)
                  write(ifmt,*)'ERROR NaN catch uinitial (2)'
                  call utstop('ERROR NaN catch uinitial (2)&')
                endif
                if(pptl(4,nn).lt.2.)then
                  a=pptl(4,nn)**2-pptl(3,nn)**2
                  b=pptl(5,nn)**2+pptl(1,nn)**2+pptl(2,nn)**2
                  idyy=idptl(nn)
                  if(abs(a-b).gt.0.1)write(ifch,*)
     .            'mass shell pb  ini :',idyy,a,b,ityptl(nn),nn
                endif
              endif
              !#########################
              call DelayFormation(nn)
              r0(itmp)=xorptl(4,nn)
              rx(itmp)=xorptl(1,nn)
              ry(itmp)=xorptl(2,nn)
              rz(itmp)=xorptl(3,nn)
              p0(itmp)=pptl(4,nn)
              px(itmp)=pptl(1,nn)
              py(itmp)=pptl(2,nn)
              pz(itmp)=pptl(3,nn)
              fmass(itmp)=pptl(5,nn)
              if(idptl(nn).ne.10)then
                aamt=sqrt(fmass(itmp)**2+px(itmp)**2+py(itmp)**2)
                rrap=sign(1d0,pz(itmp))
     .          *log((p0(itmp)+abs(pz(itmp)))/aamt)
                pz(itmp)=aamt*sinh(rrap)
                p0(itmp)=aamt*cosh(rrap)
                !#########################
                betatx=px(itmp)/p0(itmp)
                betaty=py(itmp)/p0(itmp)
                betatz=pz(itmp)/p0(itmp)
                if((betatx**2+betaty**2+betatz**2).gt.1.0d0)then
                  write(*,*)'NaN ; beta>1 in uinitial (2)'
                  write(*,*)'NaN   ',ityptl(nn),idptl(nn)
                  write(*,*)'NaN   ', betatx,betaty,betatz
                  write(*,*)'NaN   ',pptl(1,nn),pptl(2,nn),pptl(3,nn)
     .                             ,pptl(4,nn)
                  write(*,*)'NaN   ',px(itmp),py(itmp),pz(itmp)
     .                             ,p0(itmp)
                endif
                !#########################
              endif
              ityp(itmp)=ityptmp
              iso3(itmp)=iso3tmp
              if(ityptl(nn).eq.61)then
              charge(itmp)=ibptl(1,nn)
              else
              charge(itmp)=fchg(iso3(itmp),ityp(itmp))
              endif
              if(ishuuu.ge.2)then
                if(abs(idptl(nn))/10.eq.14.or.abs(idptl(nn))/10.eq.24
     .          .or.abs(ityptmp).eq.133.or.abs(ityptmp).eq.134)
     .          write(ifch,*)'CHARM uinitial'
     .          ,idptl(nn),idtmp,ityptmp,iso3tmp
              endif
              itypart(itmp)=ityptl(nn)
              itrace(1,itmp)=itmp
              itrace(2,itmp)=itmp
              iorpart(itmp)=0
              lstcoll(itmp)=0
              ncoll(itmp)=0
              origin(itmp)=0
              tform(itmp)=r0(itmp)
              dectime(itmp)=dectim(itmp,1)+tform(itmp)
              xtotfac(itmp)=0d0
              if(r0(itmp).lt.mintime) mintime = r0(itmp)
              !print*,'uInput:',itmp,nn,idptl(nn)
c              write(ifmt,*)'uInput:',itmp,nn,idptl(nn),idtmp
c     *        ,rx(itmp),ry(itmp),rz(itmp),px(itmp),py(itmp),pz(itmp)
            endif
          endif
        endif
      enddo
      if(kjpsi().eq.1)then
      !jpsi see epos_uj.f
      endif

      if(ihacas.eq.2)call writeEndEventTestFile

      !print*,'+++++ etot/eini=',etot/(197*200.)

      !write(ifmt,'(a,i5,3x,$)')'+++++ Nptest:',nptest
      !write(ifmt,'(a,3i5,3x,$)')'+++++ Flavor:',jca

      npart = nbar + nmes
      if(ish.ge.2)then
        write(ifch,*)' '
        write(ifch,*)'hacas nr of ptls entering: ',npart
      endif

      if(npart.eq.0)then
c        if(maproj.gt.0.and.matarg.gt.0)then
c          write(ifmt,*)'********** uinitial: no particles'
c          write(ifmt,*)'check between ',nptlpt+1,' and ',nptlbd
c          do nn=nptlpt+1,nptlbd
c            if(ityptl(nn).eq.61)then
c              idtmp=idptl(nn)
c            else
c              idtmp=idtrafo('nxs','pdg',idptl(nn))
c            endif
c            call pdg2id(ityptmp,iso3tmp,idtmp)
c            write(ifmt,*)abs(ityptmp).ge.minmes,
c     .           idptl(nn),istptl(nn) ,ityptl(nn)
c          enddo
cc         stop'091130'
c        endif
        return
      else
        np=npart
      endif

c back to the same starting time
      do i = 1, npart
         call getind(p0(i),pz(i),r0(i),rx(i),ry(i),rz(i)
     .    ,nrap,ntau,nrad,ityp(i),iso3(i),id)
          !if(ish.ge.2)write(ifch,'(4x,i6,i7,4x,3i3,3x,2f7.2,3x,2f7.2)')
          !.    i,id,ntau,nrad,nrap,sqrt(r0(i)**2-rz(i)**2)
          !.      ,sqrt(rx(i)**2+ry(i)**2)
         dpart(nrap,ntau)=dpart(nrap,ntau)+1
         dpart(nrap,nuutaumx+1)=dpart(nrap,nuutaumx+1)+1
         !save freeze-out configuration, in case of no further
         !rescatterings
         frr0(i) = r0(i)
         frrx(i) = rx(i)
         frry(i) = ry(i)
         frrz(i) = rz(i)
         frp0(i) = p0(i)
         frpx(i) = px(i)
         frpy(i) = py(i)
         frpz(i) = pz(i)
         rx(i)=rx(i)-px(i)/p0(i)*(r0(i)-mintime)
         ry(i)=ry(i)-py(i)/p0(i)*(r0(i)-mintime)
         rz(i)=rz(i)-pz(i)/p0(i)*(r0(i)-mintime)
         r0(i)=mintime
      enddo

      acttime=mintime

      !call affini
      !call traceity(0)

c      write(*,*)'DEBUG INFO (epos.f): ',mintime,npart,istmax,nbar,nmes

c keep  epos ptl list
      nptl=nptlpt
      if(ish.ge.2)call alist('(u) kept EPOS list&',0,0)
      do nn=nptlpt+1,nptlbd
        if(istptl(nn).eq.0.and.ityptl(nn).lt.70)istptl(nn)=3
        nptl=nptl+1
        if(ish.ge.2)call alist('&',nn,nn)
      enddo
      if(ish.ge.2)call alist('(u) copied unknowns&',0,0)
      do nn=nptlpt+1,nptlbd
        if(istptl(nn).eq.4)then
          nptl=nptl+1
          call utrepl(nptl,nn)
          istptl(nptl)=0
          iorptl(nptl)=nn
          if(ish.ge.2)call alist('&',nptl,nptl)
        endif 
      enddo
      
c      if(nui.eq.1.and.ish.ge.1)
c     &write(ifch,'(a,50a1)')' (info)',('u',i=1,50)

      call utprix('uiniti',ishuuu,ishini,4)

      return
      end

      subroutine DelayFormation(j)
      include "epos.inc"
      if(tauhac.le.0.)return
      if(ityptl(j).ge.60)return
      r=rangen()
      tauran=-tauhac*alog(r)*pptl(4,j)/pptl(5,j)
      do i=1,4
        xorptl(i,j)=xorptl(i,j)
     &   + pptl(i,j) / pptl(4,j) * tauran
      enddo
      end
      subroutine idcorrect(ity,id,idtmp)
      if(ity.eq.61)then
        if(id.eq.310)id=311
        if(id.eq.130)id=-311
        idtmp=id
      else
        if(id.eq.20)id=230
        if(id.eq.-20)id=-230
        idtmp=idtrafo('nxs','pdg',id)
      endif
      end !~~~~~~~~~~~~~~~~


c#####################################################################
c#####################################################################

                subroutine uexit      !transfer -> EPOS

c#####################################################################
c#####################################################################

      implicit none
c     debug and validity range
      integer nmax, nspl
      real*8 hit_sphere
      parameter (nmax = 100000) ! maximum number of particles
      parameter (nspl = 500)  ! dimension of spline arrays
      parameter (hit_sphere = 8.d0)  ! hard collision cutoff: 251 mbarn
      integer Ap, At, Zp, Zt, npart, nbar, nmes, ctag
      integer nsteps,ranseed,event,eos,dectag,uid_cnt
      integer NHardRes,NSoftRes,NDecRes,NElColl,NBlColl
      real*8  time,  acttime, bdist, ebeam, bimp,bmin,ecm
      logical success
c 7 integer
      common /sys/ npart, nbar, nmes, ctag,nsteps,uid_cnt,
     +             ranseed,event,Ap,At,Zp,Zt,eos,dectag,
     +             NHardRes,NSoftRes,NDecRes,NElColl,NBlColl
     +             ,success
      common /rsys/ time,acttime,bdist,bimp,bmin,ebeam,ecm
      real*8
     +     gw, sgw, delr, fdel, dt,
     +     da, db,
     +     Cb0, Yuk0, Pau0, Sky20, Sky30, gamSky, gamYuk, drPau, dpPau,
     +     dtimestep
c 19 real*8
      real*8 cutmax, cutPau, cutCb, cutYuk, cutSky, cutdww
      common /cuts/ cutmax, cutPau, cutCb, cutYuk, cutSky, cutdww
      real*8 spx(nspl), spPauy(nspl), outPau(nspl),
     +                spCby(nspl),  outCb(nspl),
     +                spYuky(nspl), outYuk(nspl),
     +                spSkyy(nspl), outSky(nspl),
     +                spdwwy(nspl), outdww(nspl)
      common /spdata/ spx, spPauy, outPau, spCby,  outCb,
     +                     spYuky, outYuk, spSkyy, outSky,
     +                     spdwwy, outdww
      real*8
     +     r0(nmax), rx(nmax), ry(nmax), rz(nmax),
     +     p0(nmax), px(nmax), py(nmax), pz(nmax),
     +     airx(nmax), airy(nmax), airz(nmax),
     +     aipx(nmax), aipy(nmax), aipz(nmax),
     +     aorx(nmax,4), aory(nmax,4), aorz(nmax,4),
     +     aopx(nmax,4), aopy(nmax,4), aopz(nmax,4),
     +     fmass(nmax), rww(nmax),
     +     dectime(nmax), tform(nmax), xtotfac(nmax)
      integer spin(nmax),ncoll(nmax),charge(nmax),strid(nmax),
     +        ityp(nmax),lstcoll(nmax),iso3(nmax),origin(nmax),uid(nmax)
      common/isys/spin,ncoll,charge,ityp,lstcoll,iso3,origin,strid,
     +            uid
      common /coor/ r0, rx, ry, rz, p0, px, py, pz, fmass, rww, dectime
      common /frag/ tform, xtotfac
      common /aios/ airx, airy, airz, aipx, aipy, aipz,
     +              aorx, aory, aorz, aopx, aopy, aopz
      common /pots/ Cb0, Yuk0, Pau0, Sky20, Sky30, gamSky,
     +              gamYuk, drPau, dpPau, gw, sgw, delr, fdel,
     +              dt,da, db,dtimestep
c spectator arrays
        integer smax
        parameter(smax=500)  ! maximum number of spectators
        real*8 r0s(smax), rxs(smax), rys(smax), rzs(smax),
     +               p0s(smax), pxs(smax), pys(smax), pzs(smax),
     +               sfmass(smax)
       integer sspin(smax), scharge(smax), sityp(smax), siso3(smax),
     +          suid(smax)
        integer nspec
        common /scoor/ r0s, rxs, rys, rzs, p0s, pxs ,pys, pzs, sfmass
        common /sisys/ sspin, scharge, sityp, siso3, suid
        common /ssys/ nspec
        real*8 p0td(2,nmax),pxtd(2,nmax),pytd(2,nmax),pztd(2,nmax),
     +         fmasstd(2,nmax)
        integer ityptd(2,nmax),iso3td(2,nmax)
        integer itypt(2),uidt(2),origint(2),iso3t(2)
        common /rtdelay/p0td,pxtd,pytd,pztd,fmasstd
        common /itdelay/ityptd,iso3td
        common /svinfo/itypt,uidt,origint,iso3t
        real*8 ffermpx(nmax), ffermpy(nmax), ffermpz(nmax)
        real*8 peq1, peq2
        common /ffermi/ ffermpx, ffermpy, ffermpz
        common /peq/ peq1,peq2
        integer numcto,numctp,maxstables
        parameter(numcto=400) ! maximum number of options
        parameter(numctp=400) ! maximum number of parameters
        parameter(maxstables=20) ! maximum number of stable particles
        integer   CTOption(numcto)
        real*8    CTParam(numctp)
        common /options/CTOption,CTParam
        real*8 frr0(nmax), frrx(nmax), frry(nmax), frrz(nmax),
     +     frp0(nmax), frpx(nmax), frpy(nmax), frpz(nmax)
      common /frcoor/ frr0, frrx, frry, frrz, frp0, frpx, frpy, frpz

      include 'urqmd34/comres.f'

      real taugm,rangen
      integer nptl
      real        phievt,bimevt,pmxevt,egyevt
     *,xbjevt,qsqevt,zppevt,zptevt
      integer     nevt,kolevt,koievt,kohevt,npjevt
     *,ntgevt,npnevt,nppevt,ntnevt,ntpevt,jpnevt,jppevt,jtnevt,jtpevt
     *,nglevt,minfra,maxfra,npglb,ntglb
      common/cevt/phievt,nevt,bimevt,kolevt,koievt,pmxevt,egyevt,npjevt
     *,ntgevt,npnevt,nppevt,ntnevt,ntpevt,jpnevt,jppevt,jtnevt,jtpevt
     *,xbjevt,qsqevt,nglevt,zppevt,zptevt,minfra,maxfra,kohevt
     *,npglb,ntglb
      integer istore,istmax,irescl,ntrymx,nclean,iopdg,ioidch
      real gaumx
      common/othe1/istore,istmax,gaumx,irescl,ntrymx,nclean,iopdg,ioidch
      integer nn,idpdgg,idepos  ,nptlep
      integer idtrafo
      external idtrafo
      integer fchg
      external fchg
      integer pdgid
      external pdgid
      real*8 dectim
      external dectim
      integer      ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      common/files/ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      integer      iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi
      common/prnt1/iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi
      !integer nptest, ia
      !integer jc(nflav,2),ic(2),jca(3),i,nflav
      !parameter (nflav=6)
      integer itypart(nmax), itrace(2,nmax), iorpart(nmax)
      common /city/itypart
      common /ctrace/itrace
      common /corpart/iorpart
      real a,b,c
      integer idyy,iutest,j

      integer maxphi
      parameter(maxphi=140000)
      real*8 phip0(maxphi),phipx(maxphi),phipy(maxphi),phipz(maxphi)
      real*8 phir0(maxphi),phirx(maxphi),phiry(maxphi),phirz(maxphi)
      real*8 phimass(maxphi)
      integer phinum,phiscatt(maxphi),phiidres(maxphi)

      common /phidata/ phip0,phipx,phipy,phipz,
     $      phir0,phirx,phiry,phirz,phimass,
     $      phinum,phiscatt,phiidres

      integer ishini,ishuuu
      common/cishuuu/ishuuu
      real p1,p2,p3,p4,p5,x1,x2,x3,x4,t1

c####################################################################################
c######################################   subroutine uexit      !transfer -> EPOS
c####################################################################################

      call utpri('uexit ',ishuuu,ishini,4)

      iutest=  1  !0

      do j=1,npart
         if(frp0(j).ne.frp0(j))then
            write(*,*)'NaN in uexit (1): ',j,frp0(j),p0(j)
            call utstop('NaN in uexit (1)&')
         endif
      enddo

      if(ish.ge.2)call alist('(u) final ptls&',0,0)
      if(ish.ge.2)write(ifch,*)' '

      do j=1,npart
            if(frp0(j).ne.frp0(j))then
               write(*,*)'NaN in uexit (2): ',j,frp0(j),p0(j)
               call utstop('NaN in uexit (2)&')
            endif
      enddo

      !write(ifmt,'(a)')' '
      !nptest=0
      !do i=1,3
      !jca(i)=0
      !enddo
      call getnptl(nptl)
      nptlep=nptl

      if(ish.ge.2)write(ifch,'(a,i8,a,i8,a)')
     . 'hacas: return ptls',nptlep+1,' -',nptlep+npart+phinum,' to epos'

      do nn=1,npart
        idpdgg=pdgid(ityp(nn),iso3(nn))
        idepos=idtrafo('pdg','nxs',idpdgg)
        nptl=nptl+1
        call setnptl(nptl)
        if(idepos.eq.0)then
          write(ifch,*)'pdg -->',idpdgg,ityp(nn),iso3(nn)
          call utstop('uexit: id not recognized&')
        endif
        call setiorptl(nptl,0)
        call setjorptl(nptl,0)
        call setistptl(nptl,0)
        call setifrptl(nptl,0,0)
        call setxorptl(nptl,
     .  sngl(frrx(nn)),sngl(frry(nn)),sngl(frrz(nn)),sngl(frr0(nn)))
        call setpptl(nptl,
     .  sngl(frpx(nn)),sngl(frpy(nn)),sngl(frpz(nn)),sngl(frp0(nn))
     . ,sngl(fmass(nn)))
        call setidptl(nptl,idepos)
        !#########################
        if(iutest.eq.1)then
          c=sngl(frp0(nn)) 
          if(.not.(c.le.0..or.c.ge.0.))then
            write(ifmt,*)'p=',sngl(frpx(nn)),sngl(frpy(nn))
     .      ,sngl(frpz(nn)),sngl(frp0(nn))
     .      ,sngl(fmass(nn))
            write(ifmt,*)'ERROR NaN catch uexit'
            call utstop('ERROR NaN catch uexit&')
          endif
          if(c.lt.2.)then
            a=sngl(frp0(nn))**2-sngl(frpz(nn))**2                    ! p4^2 - p3^2
            b=sngl(fmass(nn))**2+sngl(frpx(nn))**2+sngl(frpy(nn))**2 ! m^2 + p1^2 + p2^2
            if(abs(a-b).gt.0.1)write(ifch,*)
     .      'ERROR mass shell pb  uexit :',idepos,a,b
          endif
        endif
        !#########################
        call setistptl(nptl,0)
        call setityptl(nptl,itypart(nn))
        !?????????????????????????
        !f(frpx(nn)**2+frpy(nn)**2.gt.20.**2)then
        !write(ifch,'(a,i8,i3,f9.2,f6.2)')
        !.          ' HHHHHHH out ',nn,itypart(nn),frr0(nn)
        !.        ,sqrt(frpx(nn)**2+frpy(nn)**2)
        !endif
        !????????????????????????????
        call setiorptl(nptl,iorpart(nn))
        if(itypart(nn).eq.61)call setityptl(nptl,60)
        if(itypart(nn).eq.62)call setityptl(nptl,60)
        call setrinptl(nptl,-9999)
        call setzpaptl(nptl,0.,0.)
        call getpptl(nptl,p1,p2,p3,p4,p5)
        call getxorptl(nptl,x1,x2,x3,x4)
        t1=x4
        call idtau(idepos,p4,p5,taugm)
        call settivptl(nptl,t1,t1+taugm*(-alog(rangen())))
        if(ish.ge.2)call alist('&',nptl,nptl)
        if(ish.ge.2.and.abs(idepos).eq.230)
     .   write(ifch,*)'KKK',nptl,nn
c        if(idepos.eq.441)
c     .    write(ifmt,*)'JPSI out: ',ityp(nn),iso3(nn),frpz(nn)
        if(ishuuu.ge.2)then
          if(abs(idepos)/10.eq.14.or.abs(idepos)/10.eq.24
     .    .or.abs(ityp(nn)).eq.133.or.abs(ityp(nn)).eq.134)
     .    write(ifch,*)'CHARM uexit'
     .    ,idepos,idpdgg,ityp(nn),iso3(nn)
        endif
      enddo

      do nn=1,phinum
        nptl=nptl+1
        call setnptl(nptl)
        call setiorptl(nptl,0)
        call setjorptl(nptl,0)
        call setistptl(nptl,9)
        call setifrptl(nptl,0,0)
        call setxorptl(nptl,
     .  sngl(phirx(nn)),sngl(phiry(nn)),sngl(phirz(nn)),sngl(phir0(nn)))
        call setpptl(nptl,
     .  sngl(phipx(nn)),sngl(phipy(nn)),sngl(phipz(nn)),sngl(phip0(nn))
     . ,sngl(phimass(nn)))
        call setidptl(nptl,phiidres(nn))
        call setityptl(nptl,80+phiscatt(nn)) !rescattered:  80=no  81=yes
        call setrinptl(nptl,-9999)
        call setzpaptl(nptl,0.,0.)
        call getpptl(nptl,p1,p2,p3,p4,p5)
        call getxorptl(nptl,x1,x2,x3,x4)
        t1=x4
        call getidptl(nptl,idyy)
        call idtau(sign( abs(idyy)/100 , idyy ) ,p4,p5,taugm)
        call settivptl(nptl,t1,t1+taugm*(-alog(rangen())))
      enddo
      !????????????????????????????
      !nn=itrace(2,230)
      !write(ifch,'(a,i8,i3,f9.2,f6.2)')
      !. ' HHHHHHH xxx ',nn,itypart(nn),frr0(nn)
      !.  ,sqrt(frpx(nn)**2+frpy(nn)**2)
      !????????????????????????????
      !write(ifmt,'(i5)')nptest
      !write(ifmt,'(3i5)')jca

      call utprix('uexit ',ishuuu,ishini,4)

      end

c####################################################################################
c####################################################################################

                     subroutine  uepos

c####################################################################################
c####################################################################################

cc---------------------------------------------------------------------
c   main modul
cc---------------------------------------------------------------------

      implicit none
      integer idpdgg,idepos,pdgid,idtrafo
      real amepos
      include 'urqmd34/coms.f'
      include 'urqmd34/comres.f'
      include 'urqmd34/options.f'
      include 'urqmd34/colltab.f'
      include 'urqmd34/inputs.f'
      include 'urqmd34/newpart.f'
      include 'urqmd34/boxinc.f'
      integer i,j,k,steps,ii,ocharge,ncharge, mc, mp, noc, it1,it2
      real*8 sqrts,otime,xdummy,st
      logical isstable
      integer stidx,iret
      real*8 Ekinbar, Ekinmes, ESky2, ESky3, EYuk, ECb, EPau
      common /energies/ Ekinbar, Ekinmes, ESky2, ESky3, EYuk, ECb, EPau
      integer cti1sav,cti2sav
      integer      ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      common/files/ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      integer       laproj,maproj,latarg,matarg
      real         core,fctrmx
      common/nucl1/laproj,maproj,latarg,matarg,core,fctrmx
      integer      iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi
      common/prnt1/iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi
      integer nuj
      data nuj/0/
      save nuj
      integer iuskip
      common /ciuskip/iuskip
      double precision iutime
      integer kk
      double precision ddt
      COMMON /NPTLC/ddt(ncollmax,4), iutime(5)
      logical ktime


c####################################################################################
c#######################################    subroutine  uepos
c####################################################################################

      if(iuskip.eq.1)goto 777
      if(iuskip.eq.2)goto 778

      nuj=nuj+1

      mc=0
      mp=0
      noc=0
      kk=0
      ktime=.false.

      time = 0.0  !time is the system time at the BEGINNING of every timestep

      !initialize random number generator
      !call auto-seed generator only for first event and if no seed was fixed
      if(.not.firstseed.and.(.not.fixedseed)) then
         ranseed=-(1*abs(ranseed))
         call sseed(ranseed)
      else
         firstseed=.false.
      endif

      !old time if an old fort.14 is used
      if(CTOption(40).eq.1)time=acttime
      if(CTOption(40).eq.3)time=acttime

      !write headers to file
      call output(13)
      !call output(14)
      call output(15)
      call output(16)
      !if(event.eq.1)call output(17)
      call osc99_event(-1)

      !for CTOption(4)=1 : output of initialization configuration
      if(CTOption(4).eq.1)call file14out(0)
      !participant/spectator model:
      if(CTOption(28).ne.0) call rmspec(0.5d0*bimp,-(0.5d0*bimp))

      otime = outsteps*dtimestep  !compute time of output

      steps = 0  !reset time step counter

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! loop over all timesteps
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      do 20  steps=1,nsteps

         if (eos.ne.0) then
           do j=1,npart
               r0_t(j) = r0(j)
               rx_t(j) = rx(j)
               ry_t(j) = ry(j)
               rz_t(j) = rz(j)
           enddo
         end if


         do j=1,npart
            if(p0(j).ne.p0(j))then
               write(*,*)'NaN in uepos (1): ',j, p0(j),steps,time
                        call utstop('NaN in uepos (1)&')
            endif
         enddo
         !we are at the beginning of the timestep, set current time (acttime)
         acttime = time

         if(CTOption(16).ne.0) goto 103  !option for MD without collision term
c         print *,'urqmd',steps,nsteps,outsteps,dtimestep,nct

         call colload  ! Load collision table with next collisions in current timestep

         ! check for collisions in time-step, nct = # of collisions in table
         if (nct.gt.0) then
 101        continue              !entry-point for collision loop in case
            k = 0                 !      of full colload after every coll
 100        continue              !normal entry-point for collision loop
            kk=kk+1
            !call traceity(1)
            !call traceori(1,230)
            if(ktime)call timer(iutime)
            if(ktime)ddt(kk,1)=iutime(3)+0.001*iutime(4)
            call getnext(k)       !get next collision
            if (k.eq.0) goto 102  !exit collision loop if no collisions are left

            !propagate all particles to next collision time
            !store actual time in acttime, propagation time st=cttime(k)-acttime
            st=cttime(k)-acttime
            if(ktime)call timer(iutime)
            if(ktime)ddt(kk,2)=iutime(3)+0.001*iutime(4)
            call cascstep(acttime,st)
            acttime = cttime(k)   !new actual time (for upcoming collision)

            !perform collision
            !call traceity(2)
            !call traceori(2,230)

            if(cti2(k).gt.0.and.
     .           abs(sqrts(cti1(k),cti2(k))-ctsqrts(k)).gt.1d-3)then
               write(ifch,*)' ***(E) wrong collision update (col) ***'
               write(ifch,*)cti1(k),cti2(k),
     .              ctsqrts(k),sqrts(cti1(k),cti2(k))
            else if(cti2(k).eq.0.and.
     .              abs(fmass(cti1(k))-ctsqrts(k)).gt.1d-3) then
              write(ifch,*)' *** main(W) wrong collision update (decay)'
              write(ifch,*)ctag,cti1(k),ityp(cti1(k)),dectime(cti1(k)),
     .              fmass(cti1(k)),ctsqrts(k)
            endif

            ocharge=charge(cti1(k))
            if(cti2(k).gt.0) ocharge=ocharge+charge(cti2(k))
            !call traceity(3)
            !call traceori(3,230)

            !store quantities in local variables for charge conservation check
            it1= ityp(cti1(k))
            if(cti2(k).gt.0)it2= ityp(cti2(k))
            if(cti2(k).eq.0)then !~~~~~~~to avoid rare crash~~~~~~~~
              call checkdecay(ctsqrts(k),it1,iso3(cti1(k)),iret)
              if(iret.eq.1)then
                idpdgg=pdgid(ityp(cti1(k)),iso3(cti1(k)))
                idepos=idtrafo('pdg','nxs',idpdgg)
                call idmass(idepos,amepos)
                ctsqrts(k)=amepos
                fmass(cti1(k))=amepos
                !print*,'(info) checkdecay: iret=1 --> take epos mass'
                call checkdecay(ctsqrts(k),it1,iso3(cti1(k)),iret)
                if(iret.eq.1)call utstop('checkdecay: still iret=1&')
              endif
            endif !~~~~~~~~~~~~~~~~
            !call traceity(4)
            !call traceori(4,230)
            !increment "dirty" collision counter
            if(cti2(k).gt.0)then !scatter
               mc=mc+1
            endif
            !perform scattering/decay
            cti1sav = cti1(k)
            cti2sav = cti2(k)
            call scatter(cti1(k),cti2(k),ctsigtot(k),ctsqrts(k),
     &                   ctcolfluc(k))

            !update collision table
            !call traceity(5)
            !call traceori(5,230)

            !normal update mode
            if(CTOption(17).eq.0) then
               if(nexit.eq.0) then
                 !new collision partners for pauli-blocked states (nexit=0)
                 if (cti1(k).ne.cti1sav.or.cti2(k).ne.cti2sav) then
                   cti1(k) = cti1sav
                   cti2(k) = cti2sav
                 endif

                 if(ktime)call timer(iutime)
                 if(ktime)ddt(kk,3)=iutime(3)+0.001*iutime(4)
                 call collupd(cti1(k),1)
                 if(cti2(k).gt.0) call collupd(cti2(k),1)
               else
                 ncharge=0
                 !new collision partners for scattered/produced particles (nexit><0)
                 if(ktime)call timer(iutime)
                 if(ktime)ddt(kk,3)=iutime(3)+0.001*iutime(4)
                 do i=1,nexit
                   !ncharge is used for charge conservation check
                   ncharge=ncharge+charge(inew(i))
                   call collupd(inew(i),1)
                 enddo
               endif
               !update collisions for partners of annihilated particles
               do ii=1,nsav
                  call collupd(ctsav(ii),1)
               enddo
               nsav=0
            else ! (CTOption(17).ne.0)
              !full collision load
              call colload
            endif
            if(ktime)call timer(iutime)
            if(ktime)ddt(kk,4)=iutime(3)+0.001*iutime(4)
            if(ktime)print*,'*****',kk
     .       ,nint((ddt(kk,2)-ddt(kk,1))*1000)
     .       ,nint((ddt(kk,3)-ddt(kk,2))*1000)
     .       ,nint((ddt(kk,4)-ddt(kk,3))*1000)
            !call traceity(6)
            !call traceori(6,230)

            if (CTOption(17).eq.0) goto 100
            goto 101

            !this is the point to jump to after all collisions in the timestep
            !have been taken care of
 102        continue

         endif ! (nct.gt.0)

         do j=1,npart
            if(p0(j).ne.p0(j))then
               write(*,*)'NaN in uepos (2): ',j, p0(j),steps,time
                call utstop('NaN in uepos (2)&')
            endif
         enddo

         !After all collisions in the timestep are done, propagate to end of
         !the timestep.

         !point to jump to in case of MD without collision term
 103     continue

         time = time+dtimestep  !increment timestep

         !After all collisions in the timestep are done, propagate to end of
         !the timestep.
         call cascstep(acttime,time-acttime)

         !in case of potential interaction, do MD propagation step
         if (eos.ne.0) then
            ! set initial conditions for MD propagation-step
            do j=1,npart
               r0(j) = r0_t(j)
               rx(j) = rx_t(j)
               ry(j) = ry_t(j)
               rz(j) = rz_t(j)
            enddo
            !now molecular dynamics trajectories
            call proprk(time,dtimestep)
         endif ! (eos.ne.0)

         !perform output if desired
         if(mod(steps,outsteps).eq.0.and.steps.lt.nsteps)then
            if(CTOption(28).eq.2)call spectrans(otime)
            call file14outx(steps/outsteps)
         endif ! output handling

 20   continue ! time step loop

      acttime=time

      if(maproj.gt.0.and.matarg.gt.0)then
        if(npart.eq.0)call utstop(
     .'\n \n   ***** STOP in uepos: no particles (2) *****  \n\n&)')
      else
        return
      endif

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! optional decay of all unstable
  !  particles before final output
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  777 continue
           do j=1,npart
            if(p0(j).ne.p0(j))then
               write(*,*)'NaN in uepos (3): ',j, p0(j),steps,time
                call utstop('NaN in uepos (3)&')
            endif
         enddo
      !DANGER: pauli-blocked decays are not performed !!!
      if(CTOption(18).eq.0.and.CTOption(51).eq.0
     .   .or.iuskip.eq.1 ) then
         !print*,'(info) npart=',npart,' before final decay' !----------
         !no do-loop is used because npart changes in loop-structure
         i=0
         nct=0
         actcol=0
         CTOption(10)=1  !disable Pauli-Blocker for final decays
 40      continue  !decay loop structure starts here
         i=i+1
         if(dectime(i).lt.1.d30) then !if particle unstable
 41         continue
            isstable = .false.
            do stidx=1,nstable
               if (ityp(i).eq.stabvec(stidx)) then
                  !write (6,*) 'no decay of particle ',ityp(i)
                  isstable = .true.
               endif
            enddo
            if (.not.isstable) then
               !~~~~~~~to avoid rare crash~~~~~~~~
               call checkdecay(fmass(i),ityp(i),iso3(i),iret)
               if(iret.eq.1)then
                 idpdgg=pdgid(ityp(i),iso3(i))
                 idepos=idtrafo('pdg','nxs',idpdgg)
                 call idmass(idepos,amepos)
                 fmass(i)=amepos
                 print*,'(info) checkdecay: iret=1 --> '
     .           ,'take epos mass (2)'
                 call checkdecay(fmass(i),ityp(i),iso3(i),iret)
                 if(iret.eq.1)
     .           call utstop('checkdecay: still iret=1 (2)&')
               endif
               !~~~~~~~~~~~~~~~~
               call scatter(i,0,0.d0,fmass(i),xdummy) !perform decay
               !backtracing if decay-product is unstable itself
               if(dectime(i).lt.1.d30) goto 41
            endif
         endif
         !check next particle
         if(i.lt.npart) goto 40
         !print*,'(info) npart=',npart,' after final decay'  !---------------
      endif ! final decay

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !     final output
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         do j=1,npart
            if(p0(j).ne.p0(j))then
               write(*,*)'NaN in uepos (4): ',j, p0(j),steps,time
                call utstop('NaN in uepos (4)&')
            endif
         enddo
      if(CTOption(28).eq.2)call spectrans(otime)

      call file13out(nsteps)
      !call file14out(nsteps)
      call file16out
      call osc_event
      call osc99_event(1)
      call osc99_eoe

      mp=mp+npart
      if(ctag.eq.0)then
         !!!!!write(*,*)'(W) No collision in event ',event
         noc=noc+1
      endif

c      if(nuj.eq.1.and.ish.ge.1)
c     &write(ifch,'(a,50a1)')' (info)',('u',i=1,50)

  778 continue

         do j=1,npart
            if(p0(j).ne.p0(j))then
               write(*,*)'NaN in uepos (5): ',j, p0(j),steps,time
            call utstop('NaN in uepos (5)&')
            endif
         enddo

      end

      subroutine traceori(ii,io)
      implicit none
      include 'urqmd34/coms.f'
      include 'urqmd34/newpart.f'
      include 'urqmd34/freezeout.f'
      integer itrace(2,nmax)  , ii, io, i
      common /ctrace/itrace
      integer      ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      common/files/ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      i=itrace(2,io)
      write(ifch,'(a,i2,a,2i8,f10.2)')' HHHHHHH     traceori',ii,'    '
     .     ,io,i,sqrt(px(i)**2+py(i)**2)
      end

      subroutine traceity(ii)
      implicit none
      include 'urqmd34/coms.f'
      include 'urqmd34/newpart.f'
      include 'urqmd34/freezeout.f'
      integer itrace(2,nmax)  , ii, io, i
      common /ctrace/itrace
      integer      ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      common/files/ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      integer itypart(nmax)
      common /city/itypart
      integer ncnttraceity
      data ncnttraceity /0/
      save ncnttraceity
      if(ncnttraceity.gt.2000)return

      do i = 1, npart
        io=itrace(1,i)
        !if(io.ne.0)then
        if(itypart(i).eq.0.or.tform(i).lt.-0.5)then
          ncnttraceity=ncnttraceity+1
          write(ifch,'(a,i2,a,2i8,2f10.2)')
     .     ' HHHHHHH     traceity',ii,'    '
     .     ,io,i,sqrt(px(i)**2+py(i)**2), tform(i)
        endif
        !endif
      enddo
      end

c      subroutine affini
c      include "epos.inc"
c      do nn=nptlpt+1,nptlbd
c        if(istptl(nn).eq.0)then
c          if(abs(idtmp).lt.11.or.abs(idtmp).gt.16)then
c            if(xorptl(4,nn).lt.-0.5)then
c              write(ifch,'(a,2i8,i4,5f10.2)')
c     .       ' HHHHHHH   affini    ',nn,idptl(nn),ityptl(nn)
c     .       ,sqrt(pptl(1,nn)**2+pptl(2,nn)**2)
c     .       , xorptl(4,nn),xorptl(3,nn),xorptl(1,nn),xorptl(2,nn)
c             endif
c           endif
c        endif
c      enddo
c      end

c###################################################################################
c###################################################################################

                       subroutine input(io)

c###################################################################################
c###################################################################################

c-----------------------------------------------------------------------
c     This subroutine reads the UQMD input file (unit=9)
c
c input : for ({\\tt io=0} input-file will be processed,
c         else default values assumed
c output: information in common-block coms.f
c
c-----------------------------------------------------------------------

      implicit none

      include 'urqmd34/coms.f'
      include 'urqmd34/options.f'
      include 'urqmd34/comres.f'
      include 'urqmd34/inputs.f'
      include 'urqmd34/boxinc.f'

      integer laproj,maproj,latarg,matarg
      real core,fctrmx
      common/nucl1/laproj,maproj,latarg,matarg,core,fctrmx
      integer icinpu
      real engy,elepti,elepto,angmue
      common/lept1/engy,elepti,elepto,angmue,icinpu

      character*3 flag
      character*77 inputstr,file9,fheader,file14,file15,file16,file17
      character*77 file13,file10,file19,file20
      integer line,proflg,tarflg,impflg,beamflg,inx,ival,partid
      integer eosflg,i,io
      real*8 rval,caltim,outtim
      logical dtflag,bret

      character CTOStrng(numcto)*60
      character CTPStrg(numctp)*60

c setting of internal parameters values:
      real*8 valint(1)
      common /values/ valint

      logical infu

      integer      ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      common/files/ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      integer      iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi
      common/prnt1/iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi

      integer ncnt
      data ncnt /0/

      save

c##############################################################################
c######################################  subroutine input(io)
c##############################################################################

      ncnt=ncnt+1
      infu=info
      if(ncnt.gt.1.or.ish.eq.0)infu=.false.


      valint(1)=0.d0

      bret=io.ne.0
      goto 108

      entry inpini
c  called by some test programs

      bret=.true.
 108  continue

c initialize counters
      line=0
      boxflag=0
      mbflag=0
      edens=0.d0
      para=0
      solid=0
      mbox=0

c the following flags check, wether all necessary input is given
c projectile
      proflg=0
      prspflg=0
c target
      tarflg=0
      trspflg=0
c impact parameter
      impflg=0
c incident beam energy
      beamflg=0
      srtflag=0
      firstev=0
c equation of state
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
      CTOption(27)=0
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


      if(bret)return

c initialize arrays for special PRO/TAR combinations
      do 10 i=1,2
         spityp(i)=0
         spiso3(i)=0
 10   continue
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
ckw      if (file9(1:4).ne.'    ') then
ckw         OPEN(UNIT=9,FILE=file9,STATUS='OLD',FORM='FORMATTED')
ckw      endif
      if (file10(1:4).ne.'    ') then
         OPEN(UNIT=10,FILE=file10,STATUS='OLD',FORM='FORMATTED')
         CTOption(40)=1
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
c stop input if old event is read in
      if(CTOption(40).ne.0) return


c this entry is used to read cto,ctp and tim statements
c in case of old event readin
      entry getparams

ckw      close(9)
ckw      OPEN(UNIT=9,FILE=file9,STATUS='OLD',FORM='FORMATTED')

c read input lines
 1    continue
      line=line+1
ckw      read(9,99) flag,inputstr

      inputstr( 1:40)='                                        '
      inputstr(41:77)='                                     '
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !    settings
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(line.eq.1)then
        flag='pro'                         !~~~~projectile (Ap, Zp)
        write(inputstr(1:8),'(2i4)')maproj,laproj
      elseif(line.eq.2)then
        flag='tar'                         !~~~~target (Ap, Zp)
        write(inputstr(1:8),'(2i4)')matarg,latarg
      elseif(line.eq.3)then
        flag='nev'                         !~~~~number of events
        inputstr(1:2)=' 1'
      elseif(line.eq.4)then
        flag='tim'    !~~~~propagation time, output time step (in fm/c)
        inputstr(1:14)=' 150000 150000'
      elseif(line.eq.5)then
        flag='ecm'                         !~~~~cms energy in AGeV
        write(inputstr(1:8),'(i8)')nint(engy)
      elseif(line.eq.6)then
        flag='imp'                         !~~~~impact parameter (in fm)
        inputstr(1:2)=' 0'
      elseif(line.eq.7)then
        flag='   '                         !~~~~random number seed
        !inputstr(1:11)=' 1134570653'
      elseif(line.eq.8)then
        flag='eos'                         !~~~~equation of state
        inputstr(1:2)=' 0'
      elseif(line.eq.9)then
        flag='cto'
        inputstr(1:4)=' 5 1'
      elseif(line.eq.10)then
        flag='cto'
        inputstr(1:4)=' 6 1'
      elseif(line.eq.11)then
        flag='cto'
        inputstr(1:5)=' 40 3'
      elseif(line.eq.12)then
        flag='f13'
        inputstr(1:1)=' '
      elseif(line.eq.13)then
        flag='   '   !'f14'             !~~~~~~~~~suppress output
        inputstr(1:1)=' '
      elseif(line.eq.14)then
        flag='f15'
        inputstr(1:1)=' '
      elseif(line.eq.15)then
        flag='f16'
        inputstr(1:1)=' '
      elseif(line.eq.16)then
        flag='f19'
        inputstr(1:1)=' '
      elseif(line.eq.17)then
        flag='f20'
c      elseif(line.eq.18)then        !~~~~no strings
c        flag='cto'
c        inputstr(1:5)=' 12 1'
c      elseif(line.eq.??)then !~~~~Pauli blocker
c        flag='cto'
c        inputstr(1:5)=' 10 1'
c      elseif(line.eq.??)then  !~~~~~for fast run comment the following flags!!!!
c        flag='cdt'                 !~~~~~~~~~  timestep
c        inputstr(1:1)='1'
c      elseif(line.eq.??)then
c        flag='tim'                 !~~~~~~~~~ caltim, outtim
c        inputstr(1:6)='400 20'
      else
        goto2
      endif
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c 3    continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c select action according to flag:
c #  : treat line as a comment
      if(flag(1:1).eq.'#') goto 1
c blanks are comments, too
      if(flag(1:1).eq.' ') goto 1
c xxx: treat line as end of input marker
      if(flag.eq.'xxx'.or.flag.eq.'end') then
         goto 2
cc cal: header for output-files
c      if(flag.eq.'cal') then
c         fheader=inputstr
c pro: define projectile
      elseif(flag.eq.'pro') then
         proflg=proflg+1
         read(inputstr,fmt=*,err=88,end=88) Ap,Zp
         if(proflg.gt.1) then
            write(ifmt,*)'multiple definitions for projectile system:'
            write(ifmt,*)'-> last entry will be used'
         endif
c PRO: define special projectile
      elseif(flag.eq.'PRO') then
         proflg=proflg+1
         prspflg=1
         read(inputstr,fmt=*,err=88,end=88) spityp(1),spiso3(1)
         Ap=1
         if(proflg.gt.1) then
            write(ifmt,*)'multiple definitions for projectile system:'
            write(ifmt,*)'-> last entry will be used'
         endif
c tar: define target
      elseif(flag.eq.'tar') then
         tarflg=tarflg+1
         read(inputstr,fmt=*,err=88,end=88) At,Zt
         if(tarflg.gt.1) then
            write(ifmt,*)'multiple definitions for target system:'
            write(ifmt,*)'-> last entry will be used'
         endif
c TAR: define special target
      elseif(flag.eq.'TAR') then
         tarflg=tarflg+1
         trspflg=1
         read(inputstr,fmt=*,err=88,end=88) spityp(2),spiso3(2)
         At=1
         if(tarflg.gt.1) then
            write(ifmt,*)'multiple definitions for target system:'
            write(ifmt,*)'-> last entry will be used'
         endif
c box: define a box with a length in fm
c        parameters: 2: energie
c                    3: 1 =solid
c                    4: 1 = walls
        elseif(flag.eq.'box') then
           boxflag=boxflag+1
          read(inputstr,fmt=*,err=88,end=88) lbox,edens,solid,para
           if (edens.gt.0.d0) edensflag=1

           if (lbox.le.0) then
              write(ifmt,*) 'Error, lenght<=0'
              call utstop('Error, lenght<=0&')
           endif
           lboxhalbe=lbox/2.d0
           lboxd=lbox*2.d0

            if (edens.lt.0.d0) then
              write(ifmt,*) 'Error, a negativ energy '
              call utstop('Error, a negativ energy &')
           endif

           if(boxflag.gt.1) then
            write(ifmt,*)'multiple boxes are defined'
            call utstop('multiple boxes are defined &')
        endif
c bpt: define particles in the box
c parameters: ityp, iso3, mpart, pmax
        elseif(flag.eq.'bpt') then
           if (edens.gt.0.d0) then
              write(ifmt,*) 'Error, energie is already defined'
              call utstop('Error, energie is already defined &')
           endif
           mbox=mbox+1
           read(inputstr,fmt=*,err=88,end=88)
     &     bptityp(mbox),bptiso3(mbox),bptpart(mbox),bptpmax(mbox)
           edensflag=0
            if (bptpart(mbox).le.0) then
              write(ifmt,*) 'Error, a negativ particle number'
              call utstop('Error, a negativ particle number&')
           endif
           if(boxflag.lt.1) then
            write(ifmt,*)'no box is defined'
            call utstop('no box is defined&')
        endif
c bpe: define particles in the box with a given energy
c parameters: ityp, iso3, mpart,
        elseif(flag.eq.'bpe') then
           if (edens.le.0) then
              write(ifmt,*) 'Error, no energie is defined'
              call utstop('Error, no energie is defined&')
           endif
           mbox=mbox+1
           read(inputstr,fmt=*,err=88,end=88)
     &     bptityp(mbox),bptiso3(mbox),bptpart(mbox)
           if(boxflag.lt.1) then
            write(ifmt,*)'no box is defined'
            call utstop('no box is defined&')
        endif
c ene: beam energy (lab-system)
      elseif(flag.eq.'ene'.or.flag.eq.'elb') then
         beamflg=beamflg+1
         read(inputstr,fmt=*,err=88,end=88) ebeam
         if(beamflg.gt.1) then
            write(ifmt,*)'multiple definitions for beam-energy:'
            write(ifmt,*)'-> last entry will be used'
         endif
         if (ebeam.le.200) then
           write(ifmt,*)'Calculation at ebeam.le.200 A GeV:'
           write(ifmt,*)'parameter nmax in coms.f may be decreased!'
         endif
c plb: beam momentum (lab-system)
      elseif(flag.eq.'plb') then
         beamflg=beamflg+1
         srtflag=2
         read(inputstr,fmt=*,err=88,end=88) pbeam
         if(beamflg.gt.1) then
            write(ifmt,*)'multiple definitions for beam-energy:'
            write(ifmt,*)'-> last entry will be used'
         endif
       if (pbeam.le.200) then
            write(ifmt,*)'Calculation at pbeam.le.200 A GeV:'
            write(ifmt,*)'parameter nmax in coms.f may be decreased!'
       endif
c PLB: beam momentum ( LAb-system, excitation function possible)
      elseif(flag.eq.'PLB'.or.flag.eq.'PLG') then
         beamflg=beamflg+1
         srtflag=2
         read(inputstr,fmt=*,err=88,end=88) pbmin,pbmax,npb
         pbeam=pbmin
         if(beamflg.gt.1) then
            write(ifmt,*)'multiple definitions for beam-energy:'
            write(ifmt,*)'-> last entry will be used'
         endif
         if(npb.gt.1.and.flag.eq.'PLB') efuncflag=1
         if(npb.gt.1.and.flag.eq.'PLG') efuncflag=2
         if(abs(pbmax-pbmin).le.1.d-6) then
            npb=1
            efuncflag=0
         endif
         if (pbmax.le.200) then
            write(ifmt,*)'Calculations at pbmax.le.200 A GeV:'
            write(ifmt,*)'parameter nmax in coms.f may be decreased!'
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
            write(ifmt,*)'multiple definitions for beam-energy:'
            write(ifmt,*)'-> last entry will be used'
         endif
         if (ecm.le.20) then
          if(infu)then
        write(ifmt,*)'(info) Calculation at sroot.le.20 A GeV:'
        write(ifmt,*)'(info) parameter nmax in coms.f may be decreased!'
          endif
         endif
c ENE: beam energy (sqrt(s): CM-system, excitation function possible)
      elseif(flag.eq.'ENE'.or.flag.eq.'ELG') then
         beamflg=beamflg+1
         srtflag=1
         read(inputstr,fmt=*,err=88,end=88) srtmin,srtmax,nsrt
         ecm=srtmin
c        if(flag.eq.'ELG')ecm=1d1**dlog10(srtmin)
         if(beamflg.gt.1) then
            write(ifmt,*)'multiple definitions for beam-energy:'
            write(ifmt,*)'-> last entry will be used'
         endif
         if(nsrt.gt.1.and.flag.eq.'ENE') efuncflag=1
         if(nsrt.gt.1.and.flag.eq.'ELG') efuncflag=2
         if(abs(srtmax-srtmin).le.1.d-6) then
            nsrt=1
            efuncflag=0
         endif
         if (srtmax.le.20) then
            write(ifmt,*)'Calculations at srootmax.le.20 A GeV:'
            write(ifmt,*)'parameter nmax in coms.f may be decreased!'
         endif
c imp: impact parameter
      elseif(flag.eq.'imp') then
         bmin=0.d0
         impflg=impflg+1
         read(inputstr,fmt=*,err=88,end=88) bdist
         if(bdist.lt.0d0)then
           CTOption(5)=1
           bdist=abs(bdist)
           write(ifmt,*)'randomly choosen impact parameter:',
     ,             ' CTOption(5) is set to 1'
         end if
         if(impflg.gt.1) then
            write(ifmt,*)'multiple definitions for impact parameter:'
            write(ifmt,*)'-> last entry will be used'
         endif
c IMP: impact parameter
      elseif(flag.eq.'IMP') then
         impflg=impflg+1
         read(inputstr,fmt=*,err=88,end=88) bmin,bdist
         CTOption(5)=1
         if(impflg.gt.1) then
            write(ifmt,*)'multiple definitions for impact parameter:'
            write(ifmt,*)'-> last entry will be used'
         endif
c eos: impact parameter
      elseif(flag.eq.'eos') then
         eosflg=eosflg+1
         read(inputstr,fmt=*,err=88,end=88) eos
         if(eosflg.gt.1) then
            write(ifmt,*)'multiple definitions for equation of state:'
            write(ifmt,*)'-> last entry will be used'
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
            write(ifmt,*) 'Warning: too many stable particles defined!'
         endif
c cto: collision term options
      elseif(flag.eq.'cto') then
         read(inputstr,fmt=*,err=88,end=88) inx,ival
         if(ncnt.eq.1.and.ish.ge.1)
     &    write(ifmt,*)'(info) CTOption(',inx,')=',CTOption(inx)
     &        ,CTOStrng(inx)(1:index(CTOStrng(inx),' '))
     &      ,'is changed to',ival
         CTOption(inx)=ival
         CTOdc(inx)=' *'
c ctp: collision term parameter
      elseif(flag.eq.'ctp') then
         read(inputstr,fmt=*,err=88,end=88) inx,rval
         CTParam(inx)=rval
         CTPdc(inx)=' *'
         write(ifmt,*)'CTParam(',inx,'):   ',CTPStrg(inx)
     ,             ,'is changed to',rval
      elseif (flag.eq.'f13') then
         bf13=.true.
         if (infu) write(ifmt,*)'(info) no output on unit 13'
      elseif (flag.eq.'f14') then
         bf14=.true.
         if (infu) write(ifmt,*)'(info) no output on unit 14'
      elseif (flag.eq.'f15') then
         bf15=.true.
         if (infu) write(ifmt,*)'(info) no output on unit 15'
      elseif (flag.eq.'iou') then
         read(inputstr,fmt=*,err=88,end=88) inx,ival
         call uounit(inx,ival)
         write(ifmt,*)'file',inx,'will be written on unit',ival
      elseif (flag.eq.'f16') then
         bf16=.true.
         if (infu) write(ifmt,*)'(info) no output on unit 16'
      elseif (flag.eq.'f18') then
          bf18=.true.
          if (infu) write(ifmt,*)'(info) no output on unit 18'
      elseif (flag.eq.'f19') then
          bf19=.true.
          if (infu) write(ifmt,*)'(info) no output on unit 19'
      elseif (flag.eq.'f20') then
          bf20=.true.
          if (infu) write(ifmt,*)'(info) no output on unit 20'
      else
         write(ifmt,*)'undefined opcode in input-file on line',line
         stop
      endif
      goto 1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 2    continue


c fast CASCADE mode
      if(.not.dtflag.and.eos.eq.0) dtimestep=outtim
c
      nsteps=int(0.01+caltim/dtimestep)
      outsteps=int(0.01+outtim/dtimestep)
      if(infu)write(ifmt,*)'(info) nsteps, outsteps :'
     . ,nsteps, outsteps
      if(infu)write(ifmt,*)'(info) dtimestep :', dtimestep

c stop input if old event is read in
      if(CTOption(40).eq.1) return


c here some validity checks of the input should be performed
        if (boxflag.eq.1.and.mbox.eq.0) then
            write(ifmt,*) 'Error: no particles in the box.'
            stop
        ElseIf (boxflag.eq.0) then
      if(proflg.eq.0) then
         write(ifmt,*)'Error: no projectile specified in input.'
         stop
      elseif(tarflg.eq.0) then
         write(ifmt,*)'Error: no target specified in input.'
         stop
      elseif((impflg.eq.0)) then
         write(ifmt,*)'Error: no impact parameter in input.'
         stop
      elseif(beamflg.eq.0.and.prspflg.eq.0) then
         write(ifmt,*)'Error: no incident beam energy specified.'
         stop
      endif
c EndIf for the Box
        EndIf
      if (efuncflag.ne.0.and.
     &    mod(nevents,max(nsrt,npb)).ne.0) then
         write(ifmt,*)'INPUT: the number of events divided by the ',
     ,   'number of energies requested is no integer.'
      end if
c
c constraints for skyrme pots:
      if(eos.ne.0.and.((srtflag.eq.0.and.ebeam.gt.4d0)
     &             .or.(srtflag.eq.1.and.srtmax.gt.3.3d0)
     &             .or.(srtflag.eq.2.and.pbeam.gt.4.9))) then
         write(ifmt,*)'***(W) I switched off the potentials'
         eos=0
      end if
      if(eos.ne.0) then
         CTOption(11)=1
         CTOption(28)=0
         CTOption(30)=0
      endif
c

c now print the selected analysis

c...some input combinations should be avoided and/or commented
      if(CTOption(7).ne.0.and.At*Ap.ne.1)then
        write(ifmt,*)'Warning: CTOption(7)=',CTOption(7),
     ,  ' no elastic collisions in NN',
     ,  ' should not be used for serious calculations!'
      end if

      if(CTOption(18).ne.0)then
        write(ifmt,*)'Warning: CTOption(18)=',CToption(18),': ',
     ,  'unstable particles will not decay after propagation.'
      end if


      if(CTOption(31).ne.0)then
        write(ifmt,*)'Warning: CTOption(31)=',CToption(31),': ',
     ,  "Not yet completly implemented. Don't use for serious",
     ,  'calculations (not yet..).'
      end if

      if(CTParam(28).lt.0d0.or.CTParam(28).gt.1d0)then
        write(ifmt,*)'Warning: CTParam(28)=',CTParam(28),': ',
     ,  'should be between 0 and 1. it will be corrected.'
        CTParam(28)=min(1d0,max(0d0,CTParam(28)))
      end if

      return

 88   write(ifmt,*) 'syntax-error in input-file on line ',line
     .   ,'   flag ',flag
      write(ifmt,*)inputstr
      stop
      end


c######################################################################
c######################################################################
c##################
c##################        output
c##################
c######################################################################
c######################################################################


c--------------------------------------------------------------------------------
                  subroutine file14outx(itime)
c--------------------------------------------------------------------------------
      implicit none
      include 'urqmd34/comres.f'
      include 'urqmd34/coms.f'
      include 'urqmd34/options.f'
      include 'urqmd34/inputs.f'
      include 'urqmd34/newpart.f'
      include 'urqmd34/freezeout.f'
      include 'urqmd34/boxinc.f'
      integer i,itotcoll,iinelcoll,ii,ix,iy,itime,ncnt
      integer nexpart,idpdgg,idepos,idtrafo,pdgid,neta,ij,nrap,nrapid
      integer zeta(2,-10:10,-10:10,40),zetasum(-10:10,40)
     . ,zrapsum(-10:10,40)
      common /czeta/zeta,zetasum,zrapsum
      real*8 sigmatot,t,tf,z,x,y,p1,p2,p3,p4
      logical go
      common /outco2/sigmatot
      include 'urqmd34/outcom.f'
      data ncnt /0/
      save

      if(bf14)return
      itotcoll=ctag-dectag
      iinelcoll=itotcoll-NBlColl-NElColl
      ! print*,'(file14outx)',ttime,npart
      !@ ,itotcoll,NElColl,iinelcoll,NBlColl,dectag,
      !@     NHardRes,NSoftRes,NDecRes
       !---------------------------------------------------
       !  r0(i), rx(i), ry(i), rz(i)   ................... x4
       !  p0(i),px(i)+ffermpx(i),py(i)+ffermpy(i)
       !                    ,pz(i)+ffermpz(i),fmass(i) ... p5
       !  ityp(i)  ..... particle id
       !  iso3(i) ...... 2 times the isospin of a particle
       !  charge(i) .... charge of the particle
       !  lstcoll(i) ... index of last collision partner
       !  ncoll(i) ..... number of collisions
       !  origin(i) ....
       !  dectime(i) ...
       !  tform(i) ..... formation time
       !  xtotfac(i) ... cross section
       !                 (zero if the particle is not yet formed)
       !  uid(i) ......
       !-----------------------------------------------------
       nexpart=0
       ncnt=ncnt+1
       if(ncnt.eq.1)then
       do neta=-6,6
        zetasum(neta,itime)=0
       enddo
       do nrap=-10,10
        zrapsum(nrap,itime)=0
       enddo
       do ii=1,2
       do ij=-10,10
       do neta=-10,10
        zeta(ii,ij,neta,itime)=0
       enddo
       enddo
       enddo
       endif
       do i=1,npart
         t=r0(i)
         tf=tform(i)
         if(t.gt.tf)then
          nexpart=nexpart+1
          idpdgg=pdgid(ityp(i),iso3(i))
          idepos=idtrafo('pdg','nxs',idpdgg)
          ! some codes like 40323 are not in list -> idepos=0
          z=rz(i)
          x=rx(i)
          y=ry(i)
          ix=-10+(x+10.5)
          iy=-10+(y+10.5)
          p4=p0(i)
          p1=px(i)+ffermpx(i)
          p2=py(i)+ffermpy(i)
          p3=pz(i)+ffermpz(i)
          nrap=nrapid(p3,p4)
          zrapsum(nrap,itime)=zrapsum(nrap,itime)+1
          neta=nrapid(z,t)
          zetasum(neta,itime)=zetasum(neta,itime)+1
          do ii=1,2
           go=.false.
           if(ii.eq.1.and.abs(y).le.1.)go=.true.
           if(ii.eq.2.and.abs(x).le.1.)go=.true.
           if(go)then
            if(ii.eq.1)ij=ix
            if(ii.eq.2)ij=iy
            if(ij.ge.-10.and.ij.le.10)then
             zeta(ii,ij,neta,itime)=zeta(ii,ij,neta,itime)+1
            endif
           endif
          enddo
         endif
       enddo

       !write(ifmt,'(a,2i3,i6)')'(file14outx)',itime,nsteps,nexpart

       !do ii=1,2
       !print*,' '
       !do neta=-6,6
       !write(ifmt,'(3x,13i4)')(zeta(ii,ij,neta,itime),ij=-5,5)
       !enddo
       !enddo

      end


c----------------------------------------------------------------------
            subroutine uplot
c----------------------------------------------------------------------

      include 'urqmd34/coms.f'
      !if(nsteps.gt.1)call uplot3
      end

      !------------------------
      subroutine uplot3
      !------------------------
      include "epos.inc"
      integer zeta(2,-10:10,-10:10,40),zetasum(-10:10,40)
     . ,zrapsum(-10:10,40)
      common /czeta/zeta,zetasum,zrapsum
      character*4 ch
      write(ifhi,'(a)')    '!-------------------'
      write(ifhi,'(a,i3)') '!   uplot2     '
      write(ifhi,'(a)')    '!-------------------'
      print*,'nevent=',nevent
      do itime=1,31,2
       call getchar(itime,ch)
       write(ifhi,'(a)')       '!newpage'
       write(ifhi,'(a)')'openhisto htyp his name u2-'//ch
       write(ifhi,'(a,f4.1)')'xmod lin xrange -5 5'
       write(ifhi,'(a)')    'txt  "xaxis y "'
       write(ifhi,'(a)') 'ymod lin yrange auto auto '
       write(ifhi,'(a,i2,a)')'text 0.6 0.9 "  [t]=',itime,'"'
       write(ifhi,'(a)')'txt "yaxis dn/dy "'
       write(ifhi,'(a)')'array 2'
       do nrap=-5,5
        x=nrap
        y=zrapsum(nrap,itime)
        write(ifhi,'(2e13.5)')x,y/nevent
       enddo
       write(ifhi,'(a)') 'endarray closehisto plot 0'
      enddo
      end

      !------------------------
      subroutine uplot2
      !------------------------
      include "epos.inc"
      integer zeta(2,-10:10,-10:10,40),zetasum(-10:10,40)
     . ,zrapsum(-10:10,40)
      common /czeta/zeta,zetasum,zrapsum
      character*4 ch
      write(ifhi,'(a)')    '!--------------------------'
      write(ifhi,'(a,i3)') '!   uplot2     '
      write(ifhi,'(a)')    '!--------------------------'
      print*,'nevent=',nevent
      do itime=1,31,2
       call getchar(itime,ch)
       write(ifhi,'(a)')       '!newpage'
       write(ifhi,'(a)')'openhisto htyp his name u2-'//ch
       write(ifhi,'(a,f4.1)')'xmod lin xrange -5 5'
       write(ifhi,'(a)')    'txt  "xaxis [c] "'
       write(ifhi,'(a)') 'ymod lin yrange auto auto '
       write(ifhi,'(a,i2,a)')'text 0.6 0.9 "  [t]=',itime,'"'
       write(ifhi,'(a)')'txt "yaxis dn/d[c] "'
       write(ifhi,'(a)')'array 2'
       do neta=-5,5
        x=neta
        y=zetasum(neta,itime)
        write(ifhi,'(2e13.5)')x,y/nevent
       enddo
       write(ifhi,'(a)') 'endarray closehisto plot 0'
      enddo
      end

      !------------------------
      subroutine uplot1
      !------------------------
      include "epos.inc"
      integer zeta(2,-10:10,-10:10,40),zetasum(-10:10,40)
     . ,zrapsum(-10:10,40)
      common /czeta/zeta,zetasum,zrapsum
      character*4 ch
      write(ifhi,'(a)')    '!-----------------------'
      write(ifhi,'(a,i3)') '!   uplot1    '
      write(ifhi,'(a)')    '!-----------------------'
      np=0
      do neta=-4,4,2
      do ii=1,2
      do itime=1,31,2
       np=np+1
       call getchar(np,ch)
       write(ifhi,'(a)')       '!newpage'
       write(ifhi,'(a)')'openhisto htyp his name u2-'//ch
       write(ifhi,'(a,f4.1)')'xmod lin xrange -10 10'
       if(ii.eq.1)write(ifhi,'(a)')    'txt  "xaxis x (fm)"'
       if(ii.eq.2)write(ifhi,'(a)')    'txt  "xaxis y (fm)"'
       write(ifhi,'(a)') 'ymod lin yrange auto auto '
       write(ifhi,'(a,i2,a)')'text 0.1 0.9 "  [c]=',neta,'"'
       write(ifhi,'(a,i2,a)')'text 0.6 0.9 "  [t]=',itime,'"'
       write(ifhi,'(a)')'txt "yaxis ptl density "'
       write(ifhi,'(a)')'array 2'
       do ij=-10,10
        x=ij
        y=zeta(ii,ij,neta,itime)
        write(ifhi,'(2e13.5)')x,y/nevent/2.
       enddo
       write(ifhi,'(a)') 'endarray closehisto plot 0'
      enddo
      enddo
      enddo
      end


c---------------------------------------------------------------------------------
      subroutine checkdecay(mm1,ii1,iiz1,iret)
c---------------------------------------------------------------------------------

cinput mm1   : mass of  particle
cinput ii1   : ID of  particle
cinput iiz1  : $2\cdot I_3$ of particle


      implicit none
      integer i1,iz1,ii1,iiz1,is,iret
      real*8 m1,mm1
      include 'urqmd34/comres.f'
      include 'urqmd34/options.f'
      integer strit

      integer      ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      common/files/ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio

      iret=0
      i1=ii1
      iz1=iiz1
      m1=mm1

      is=iabs(strit(i1))

      if(iabs(i1).ge.minmes)then ! meson dec.
         call anndexx(m1,i1,iz1,
     .        maxbrm ,minmes+1,maxmes,bmtype,branmes,iret)
      else if(is.eq.0)then   ! n*,d,d*
         call anndexx(m1,i1,iz1,
     .        maxbra,minnuc+1,maxdel,brtype,branres,iret)
      else if(is.eq.1)then   !
         call anndexx(m1,i1,iz1,
     .         maxbrs1,minlam+1,maxsig,bs1type,branbs1,iret)
      else if(is.eq.2)then
         call anndexx(m1,i1,iz1,
     .        maxbrs2,mincas+1,maxcas,bs2type,branbs2,iret)
      else
         write(ifmt,*)'make22(anndexx): s=',is,'not included'
         call utstop('make22(anndexx): s not included&')
      end if
      return
      end

      subroutine anndexx(m1,i1,iz1,
     &            maxbr,mini,maxi,btype,branch,iret)

      implicit none

      include 'urqmd34/comres.f'
      include 'urqmd34/comwid.f'
      include 'urqmd34/newpart.f'
      include 'urqmd34/options.f'

      real*8 pi,cc
      parameter(pi=3.1415927,cc=0.38937966)
      integer maxbr,mini,maxi,btype(4,0:maxbr)
      real*8 branch(0:maxbr,mini:maxi)
      integer icnt,iret
      integer i,i1,iz1    !,j
      real*8 m1,prob(0:100),sum
      real*8 mminit,fbrancx
      integer isoit

      do 3 i=0,maxbr
         if(isoit(btype(1,i))+isoit(btype(2,i))+isoit(btype(3,i))+
     &      isoit(btype(4,i)).lt.iabs(iz1).or.
     &        m1.lt.mminit(btype(1,i))+mminit(btype(2,i))
     &             +mminit(btype(3,i))+mminit(btype(4,i)) )then
            prob(i)=0.d0
         else
            prob(i)=fbrancx(i,iabs(i1),iz1,m1,branch(i,iabs(i1)),
     &           btype(1,i),btype(2,i),btype(3,i),btype(4,i))
         endif
 3    continue

      icnt=0
      call getbran(prob,0,100,sum,0,maxbr,i)

      if(i.gt.maxbr)then
         iret=1
         !write(ifmt,*)'anndexx(dec): no final state found for:',i1,m1,iz1
         !write(ifmt,*)'please check minimal masses: m1,m1min,m2min'
         !write(ifmt,*)'and iso3 of decaying particle'
         !write(ifmt,*)(prob(j),j=0,maxbr)
         !stop
      end if


      return
      end

c----------------------------------------------------------------------------------
                   subroutine f15outchxx(colldens)
c----------------------------------------------------------------------------------

      implicit none
      include 'urqmd34/comres.f'
      include 'urqmd34/coms.f'
      include 'urqmd34/options.f'
      include 'urqmd34/inputs.f'
      include 'urqmd34/newpart.f'
      include 'urqmd34/freezeout.f'
      include 'urqmd34/boxinc.f'
      include 'urqmd34/outcom.f'
      integer      ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      common/files/ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      integer ninit,id,nupa(0:10)
      common/cnupa/nupa
      common /cninit/ninit
      integer nrap,ntau,i,ii
      integer nuutaumx,nrad,nin
      real*8 sigmatot,colldens,deluutau
      common /outco2/sigmatot
      real*8 dcoll(-10:10,-1:40),rcoll(0:20,-1:40),zevents
      real*8 dpart(-10:10,-1:40)
      common /ccoll/dcoll,rcoll,dpart,zevents
      common /cf15/deluutau,nuutaumx
      integer      iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi
      common/prnt1/iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi
      real*8 colldensdummy
      colldensdummy=colldens

      nin=ninit
ccc  if (bf15) return
      !if(ish.ge.2)write(ifch,*)'-----------'
      do i=1,nin
      call getind(tp0(i),tpz(i),tr0(i),trx(i),try(i),trz(i)
     . ,nrap,ntau,nrad,tityp(i),tiso3(i),id)
      !if(ish.ge.2)write(ifch,'(4x,i6,i5,3x,3i3)')tind(i),id,ntau,nrad,nrap
      !write(ifch,'(i6,$)')id
       !         istr=strit(tityp(i))
       !         ich = fchg(tiso3(i),tityp(i))
       !
       !         write(15,*) tind(i),tr0(i),trx(i),try(i),trz(i),
       !     @                   tp0(i),tpx(i),tpy(i),tpz(i),tm(i),
       !     @                   tityp(i),tiso3(i),ich,tlcoll(i),
       !     @                   tcoll(i),istr,torigin(i)
      enddo
      !if(ish.ge.2)write(ifch,*)'     --------->'
      !write(ifch,'(a,$)')'     --------->  '
      if(nexit.le.10)nupa(nexit)=nupa(nexit)+1
      do ii=1,nexit
       i=inew(ii)
       call getind(p0(i),pz(i),r0(i),rx(i),ry(i),rz(i)
     .  ,nrap,ntau,nrad,ityp(i),iso3(i),id)
       dcoll(nrap,ntau)=dcoll(nrap,ntau)+1
       rcoll(nrad,ntau)=rcoll(nrad,ntau)+1
       dcoll(nrap,nuutaumx+1)=dcoll(nrap,nuutaumx+1)+1
       rcoll(nrad,nuutaumx+1)=rcoll(nrad,nuutaumx+1)+1
      !if(ish.ge.2)write(ifch,'(4x,i6,i5,3x,3i3)')i,id,ntau,nrad,nrap
      !write(ifch,'(i6,$)')id
      enddo
      !write(ifch,'(a)')' '
      end

      subroutine xf15
      integer nupa(0:10)
      common/cnupa/nupa
      integer      ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      common/files/ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      integer nuutaumx
      real*8 deluutau,avpa,sum
      real*8 dcoll(-10:10,-1:40),rcoll(0:20,-1:40),zevents
      real*8 dpart(-10:10,-1:40)
      common /ccoll/dcoll,rcoll,dpart,zevents
      common /cf15/deluutau,nuutaumx
      integer      iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi
      common/prnt1/iprmpt,ish,ishsub,irandm,irewch,iecho,modsho,idensi
      if(ish.ge.3)write(ifch,*)'Initial particles(y,tau)'
      do j=-1,nuutaumx/2
      if(ish.ge.3)write(ifch,'(5x,11i5)')
     . (nint(dpart(i,j)/zevents),i=-5,5)
      enddo
      if(ish.ge.3)write(ifch,'(5x,11i5)')
     . (nint(dpart(i,nuutaumx+1)/zevents),i=-5,5)
      if(ish.ge.3)write(ifch,*)'ParticlesFromCollisions(y,tau)'
      do j=-1,nuutaumx+1
      if(ish.ge.3)write(ifch,'(5x,11i5)')
     . (nint(dcoll(i,j)/zevents),i=-5,5)
      enddo
      avpa=0
      sum=0
      do i=0,10
        sum=sum+nupa(i)
        avpa=avpa+nupa(i)*i
      enddo
      if(sum.gt.0.)then
        avpa=avpa/sum
      else
        avpa=0
      endif 
      if(ish.ge.3)write(ifch,*)'nexit distribution:'
      if(ish.ge.3)write(ifch,'(5x,7i5,3x,a,f5.2)')
     . (nint(nupa(j)/zevents),j=0,6)
     . ,'mean=',avpa
      if(ish.ge.3)write(ifch,*)
     . 'ParticlesFromCollisions(y,tau)/mean_nexit'
      if(ish.ge.3)write(ifch,'(5x,11i5)')
     . (nint(dcoll(i,nuutaumx+1)/zevents/avpa),i=-5,5)
      if(ish.ge.3)write(ifch,*)
     . 'Collisions/Initial_particle(y,tau)'
      if(ish.ge.3)write(ifch,'(5x,11f5.2)')
     . (dcoll(i,nuutaumx+1)/avpa/max(0.0001,dpart(i,nuutaumx+1)),i=-5,5)
      !if(ish.ge.3)write(ifch,*)'Collisions(rad,tau)'
      !do j=-1,nuutaumx
      !if(ish.ge.3)write(ifch,'(16i4)')(rcoll(i,j),i=1,16)
      !enddo
      end

      subroutine getind(p0,pz,r0,rx,ry,rz,nrap,ntau,nrad,ityp,iso3,id)
      implicit none
      real*8 pz,p0,r0,rx,ry,rz,tau,deluutau,rad
      integer nrap,nrapid,nuutaumx,ntau,nrad,ityp,iso3
     . ,id,idpdgg,pdgid
      common /cf15/deluutau,nuutaumx
      deluutau=1
      nuutaumx=20
      nrap=nrapid(pz,p0)
      ntau=-1
      tau=r0**2-rz**2
      if(tau.gt.0d0)then
       tau=sqrt(tau)
       if(tau+deluutau/2.gt.(nuutaumx+1)*deluutau)then
        ntau=nuutaumx+1
       else
        ntau=(tau+deluutau/2)/deluutau
        ntau=min(ntau,nuutaumx+1)
       endif
      endif
      rad=-1
      rad=rx**2+ry**2
      nrad=1
      if(rad.gt.0d0)then
       rad=sqrt(rad)
       if(rad.gt.20)then
         nrad=20
       else
        nrad=1+rad
        nrad=min(nrad,20)
       endif
      endif
      idpdgg=pdgid(ityp,iso3)
      id=idpdgg
      end

c------------------------------------------------------------------------
      subroutine getResoId(ik1,ik2,ik3,ik4,ka)
c------------------------------------------------------------------------
      !  ik1 ... decaying particle ityp code
      !  ik2 ... decaying particle iso3 code
      !  ik3 ... second daughter ityp code
      !  ik4 ... second daughter iso3 code
      !  ka .... resonance id + decay mode
      !-------------------------------------------------

       jk1=ik1
       jk2=ik2
       jk3=ik3
       jk4=ik4

       if(ik1.lt.0)then
       ! if the mother has negative ityp, convert everything to anti-particles
         jk1=-ik1
         jk2=-ik2
         jk3=-ik3
         jk4=-ik4
       endif

       if(jk3.eq.-100) jk3=100 !photon
       if(jk3.eq.-101) jk3=101 !pion

       !the following numbers are obtained from function pdgid in file UR/ityp2pdg.f
       !relevant codes:   id_UrQMD   id_   id_
       !                 ityp,iso3   PDG   EPOS
       !     photon        100, 0,    22    10
       !     pi-           101,-2,  -211   -120
       !     pi0           101, 0,   111    110
       !     pi+           101, 2,   211    120
       !     eta           102, 0,   221    220
       !     omega         103, 0,   223    221
       !     rho-          104,-2,  -213   -121
       !     rho0          104, 0,   113    111
       !     rho+          104, 2,   213    121
       !     K0            106,-1,   311    230
       !     K0bar        -106, 1,  -311   -230
       !     K+            106, 1,   321    130
       !     K-           -106,-1,  -321   -130
       !     eta'          107, 0,   331    330
       !     K*0           108,-1,   313    231
       !     K*0bar       -108, 1,  -313   -231
       !     K*+           108, 1,   323    131
       !     K*-          -108,-1,  -323   -131
       !     phi           109, 0,   333    331
       !     n               1,-1,  2112   1220
       !     p               1, 1,  2212   1120
       !     Lambda         27, 0,  3122   2130
       !     Sigma-         40,-2,  3112   2230
       !     Sigma0         40, 0,  3212   1230
       !     Sigma+         40, 2,  3222   1130
       !     Xi-            49,-1,  3312   2330
       !     Xi0            49, 1,  3322   1330
       !     Lambda(1520)   29, 0,  3124   1234
       !     Delta(1232)-   17,-3,  1114   2221
       !     Delta(1232)0   17,-1,  2114   1221
       !     Delta(1232)+   17, 1,  2214   1121
       !     Delta(1232)++  17, 3,  2224   1111
       !     Sigma(1385)-   41,-2,  3114   2231
       !     Sigma(1385)0   41, 0,  3214   1231
       !     Sigma(1385)+   41, 2,  3224   1131
       !     Xi(1530)-      50,-1,  3314   2331
       !     Xi(1530)0      50, 1,  3324   1331

       ka=0

       if(jk1.eq.103)then
         if(jk3.eq.101.and.jk4.ne.0)ka=22100 !omega -> pi+ pi- (0.016)
         if(jk3.eq.101.and.jk4.eq.0)ka=22101 !omega -> pi0 pi0 (0.008)
       endif

       if(jk1.eq.104.and.jk2.eq.-2)then
         if(jk3.eq.101.and.jk4.eq.-2)ka=-12100 !rho- -> pi- pi0 (1.0)
         if(jk3.eq.101.and.jk4.eq. 0)ka=-12100 !rho- -> pi0 pi- (1.0)
       endif

       if(jk1.eq.104.and.jk2.eq. 0)then
         if(jk3.eq.101.and.jk4.eq.-2)ka=11100 !rho0 -> pi- pi+ (1.0)
         if(jk3.eq.101.and.jk4.eq. 2)ka=11100 !rho0 -> pi+ pi- (1.0)
       endif

       if(jk1.eq.104.and.jk2.eq. 2)then
         if(jk3.eq.101.and.jk4.eq. 2)ka=12100 !rho+ -> pi+ pi0 (1.0)
         if(jk3.eq.101.and.jk4.eq. 0)ka=12100 !rho+ -> pi0 pi+ (1.0)
       endif

       if(jk1.eq.107.and.jk2.eq.0)then
         if(jk3.eq.104.and.jk4.eq.0)ka=33000 !eta' -> photon rho (0.295)
         if(jk3.eq.103.and.jk4.eq.0)ka=33003 !eta' -> photon omega  (0.025)
         if(jk3.eq.100.and.jk4.eq.0)ka=33004 !eta' -> photon photon (0.025)
       endif

       if(jk1.eq.108.and.jk2.eq.-1)then
         if(jk3.eq.101.and.jk4.eq.-2)ka=23100 !K*0 -> pi- K+ (0.667)
         if(jk3.eq.106.and.jk4.eq. 1)ka=23100 !K*0 -> K+ pi- (0.667)
         if(jk3.eq.101.and.jk4.eq. 0)ka=23101 !K*0 -> pi0 K0 (0.333)
         if(jk3.eq.106.and.jk4.eq.-1)ka=23101 !K*0 -> K0 pi0 (0.333)
       endif

       if(jk1.eq.108.and.jk2.eq. 1)then
         if(jk3.eq.101.and.jk4.eq. 2)ka=13100 !K*+ -> pi+ K0 (0.667)
         if(jk3.eq.106.and.jk4.eq.-1)ka=13100 !K*+ -> K0 pi+ (0.667)
         if(jk3.eq.101.and.jk4.eq. 0)ka=13101 !K*+ -> pi0 K+ (0.333)
         if(jk3.eq.106.and.jk4.eq. 1)ka=13101 !K*+ -> K+ pi0 (0.333)
       endif

       if(jk1.eq.109)then
         if(jk3.eq.106.and.jk4.eq. 1)ka=33100 !phi -> K+ K- (0.5)
         if(jk3.eq.-106.and.jk4.eq.-1)ka=33100 !phi -> K+ K- (0.5)
         if(jk3.eq.106.and.jk4.eq.-1)ka=33101 !phi -> K0 K0bar (0.35)
         if(jk3.eq.-106.and.jk4.eq. 1)ka=33101 !phi -> K0bar K0 (0.35)
         if(jk3.eq.104)ka=33102 !phi -> pi rho (0.13)
         if(jk3.eq.101)ka=33102 !phi -> pi rho (0.13)
       endif

       if(jk1.eq. 29.and.jk2.eq. 0)then
         if(jk3.eq.-106.and.jk4.eq.-1)ka=123400 !Lambda(1520) -> K- proton (0.225)
         if(jk3.eq.  1.and.jk4.eq. 1)ka=123400 !Lambda(1520) -> proton K- (0.225)
         if(jk3.eq.-106.and.jk4.eq. 1)ka=123401 !Lambda(1520) -> anti-K0 n (0.225)
         if(jk3.eq.  1.and.jk4.eq.-1)ka=123401 !Lambda(1520) -> n anti-K0 (0.225)
         if(jk3.eq.101.and.jk4.eq.-2)ka=123402 !Lambda(1520) -> pi- Sigma+ (0.14)
         if(jk3.eq. 40.and.jk4.eq. 2)ka=123402 !Lambda(1520) -> Sigma+  pi- (0.14)
         if(jk3.eq.101.and.jk4.eq. 2)ka=123403 !Lambda(1520) -> pi+ Sigma- (0.14)
         if(jk3.eq. 40.and.jk4.eq.-2)ka=123403 !Lambda(1520) -> Sigma- pi+ (0.14)
       endif

       if(jk1.eq. 17.and.jk2.eq.-3)then
         if(jk3.eq.101.and.jk4.eq.-2)ka=222100 !Delta(1232)- -> pi- neutron (1.0)
         if(jk3.eq.  1.and.jk4.eq.-1)ka=222100 !Delta(1232)- -> neutron pi- (1.0)
       endif

       if(jk1.eq. 17.and.jk2.eq.-1)then
         if(jk3.eq.101.and.jk4.eq.-2)ka=122100 !Delta(1232)0 -> pi- proton (0.333)
         if(jk3.eq.  1.and.jk4.eq. 1)ka=122100 !Delta(1232)0 -> proton pi- (0.333)
         if(jk3.eq.101.and.jk4.eq. 0)ka=122101 !Delta(1232)0 -> pi0 neutron (0.666)
         if(jk3.eq.  1.and.jk4.eq.-1)ka=122101 !Delta(1232)0 -> neutron pi0 (0.666)
       endif

       if(jk1.eq. 17.and.jk2.eq. 1)then
         if(jk3.eq.101.and.jk4.eq. 2)ka=112100 !Delta(1232)+ -> pi+ neutron (0.333)
         if(jk3.eq.  1.and.jk4.eq.-1)ka=112100 !Delta(1232)+ -> neutron pi+ (0.333)
         if(jk3.eq.101.and.jk4.eq. 0)ka=112101 !Delta(1232)+ -> pi0 proton (0.666)
         if(jk3.eq.  1.and.jk4.eq. 1)ka=112101 !Delta(1232)+ -> proton pi0 (0.666)
       endif

       if(jk1.eq. 17.and.jk2.eq. 3)then
         if(jk3.eq.101.and.jk4.eq. 2)ka=111100 !Delta(1232)++ -> pi+ proton (1.0)
         if(jk3.eq.  1.and.jk4.eq. 1)ka=111100 !Delta(1232)++ -> proton pi+ (1.0)
       endif

       if(jk1.eq. 41.and.jk2.eq.-2)then
         if(jk3.eq.101.and.jk4.eq.-2)ka=223100 !Sigma(1385)- -> pi- Lambda (0.87)
         if(jk3.eq. 27.and.jk4.eq. 0)ka=223100 !Sigma(1385)- -> Lambda pi- (0.87)
       endif

       if(jk1.eq. 41.and.jk2.eq. 0)then
         if(jk3.eq.101.and.jk4.eq. 0)ka=123100 !Sigma(1385)0 -> pi0 Lambda (0.87)
         if(jk3.eq. 27.and.jk4.eq. 0)ka=123100 !Sigma(1385)0 -> Lambda pi0 (0.87)
       endif

       if(jk1.eq. 41.and.jk2.eq. 2)then
         if(jk3.eq.101.and.jk4.eq. 2)ka=113100 !Sigma(1385)+ -> pi+ Lambda (0.87)
         if(jk3.eq. 27.and.jk4.eq. 0)ka=113100 !Sigma(1385)+ -> Lambda pi+ (0.87)
       endif

       if(jk1.eq. 50.and.jk2.eq.-1)then
         if(jk3.eq.101.and.jk4.eq.-2)ka=233100 !Xi(1530)- -> pi- Xi0 (0.667)
         if(jk3.eq. 49.and.jk4.eq. 1)ka=233100 !Xi(1530)- -> Xi0 pi- (0.667)
         if(jk3.eq.101.and.jk4.eq. 0)ka=233101 !Xi(1530)- -> pi0 Xi- (0.333)
         if(jk3.eq. 49.and.jk4.eq.-1)ka=233101 !Xi(1530)- -> Xi- pi0 (0.333)
       endif

       if(jk1.eq. 50.and.jk2.eq. 1)then
         if(jk3.eq.101.and.jk4.eq. 2)ka=133100 !Xi(1530)0 -> pi+ Xi- (0.667)
         if(jk3.eq. 49.and.jk4.eq.-1)ka=133100 !Xi(1530)0 -> Xi- pi+ (0.667)
         if(jk3.eq.101.and.jk4.eq. 0)ka=133101 !Xi(1530)0 -> pi0 Xi0 (0.333)
         if(jk3.eq. 49.and.jk4.eq. 1)ka=133101 !Xi(1530)0 -> Xi0 pi0 (0.333)
       endif

       if(ik1.lt.0) ka=-ka

       end

c------------------------------------------------------------------------
      subroutine printReaction
c------------------------------------------------------------------------
      !id of out-particles determined in subroutine anndex
      !branch i from call getbran(prob,0,100,sum,0,maxbr,i)
      !prob contains the branching ratios
      !-------------------------------------------------------
      include 'urqmd34/newpart.f'
      include 'urqmd34/outcom.f'
      include 'urqmd34/coms.f'
      include 'urqmd34/comres.f'
      integer ninit,idpdgg,pdgid,idtrafo,id,ii,kk
      common /cninit/ninit
      !logical lstrange
      logical lstr
      integer      ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      common/files/ifop,ifmt,ifch,ifcx,ifhi,ifdt,ifcp,ifdr,ifio
      integer itypart(nmax),iorpart(nmax)
      common /city/itypart
      common /corpart/iorpart
        real*8 frr0(nmax), frrx(nmax), frry(nmax), frrz(nmax),
     +     frp0(nmax), frpx(nmax), frpy(nmax), frpz(nmax)
      common /frcoor/ frr0, frrx, frry, frrz, frp0, frpx, frpy, frpz
      integer iochan,mm
      common/ciochan/iochan
      real*8  tformii(2),pxii(2),pyii(2)
      integer indii(2),iso3ii(2),itypii(2),ityii(2)
      common /saveii/ tformii,pxii,pyii,ityii,indii,iso3ii,itypii
      real*8  pzii(2)
      real*8  etaii(2)
      common /saveii2/pzii
      integer phiidii(2)                             !### kw ###
      common /saveii3/phiidii                        !### kw ###
      integer itrace(2,nmax)
      integer ityhard, itysoft, ityexot, idxxx(2)
      integer ilinexxx                               !### kw ###
      common /cilinexxx/ilinexxx                     !### kw ###
      common /ctrace/itrace
      common /cishuuu/ishuuu
      real*8 th, ctheta1
      real*8 phi1,phi2
      real*8 rpott(2),rpotx(2),rpoty(2),rpotz(2)
      real*8 pzi1,pzi2,pxi1,pxi2,pyi1,pyi2
      real*8 rstringx(2),rstringy(2),rstringz(2),tstring(2),tformold(2)
      common /scatcomr/rstringx,rstringy,rstringz,tstring,
     &                rpott,rpotx,rpoty,rpotz,
     &                pzi1,pzi2,pxi1,pxi2,pyi1,pyi2,
     &                ctheta1,phi1,th,phi2,tformold
      integer maxphi
      parameter(maxphi=140000)
      real*8 phip0(maxphi),phipx(maxphi),phipy(maxphi),phipz(maxphi)
      real*8 phir0(maxphi),phirx(maxphi),phiry(maxphi),phirz(maxphi)
      real*8 phimass(maxphi)
      real*8 massit
      integer phinum,phiscatt(maxphi),phiidres(maxphi)
      common /phidata/ phip0,phipx,phipy,phipz,
     $      phir0,phirx,phiry,phirz,phimass,
     $      phinum,phiscatt,phiidres

      call utpri('prirea',ishuuu,ishini,4)
      mm=0
      ityhard=0
      itysoft=0   !soft non-fluid
      ityexot=0   !exotic stuff
      itybulk=0   !bulk
      etaii(1)=0
      etaii(2)=0
      idxxx(1)=0
      idxxx(2)=0
      ka=0

      kmin=1
      kmax= 1 !2
      if(ishuuu.ge.2)kmax=2
      ipdg=0   !ipdg may  be set 0 or 1
      lstr=.false.

      if(ninit.eq.1)then
        ik1=itypii(1)           !decaying particle ityp code
        if(nexit.eq.2)then      !consider only two ptl decays
          ik2=iso3ii(1)         !decaying particle iso3 code
          ik3=ityp(inew(2))     !second daughter ityp code   (not first because some BR have the same 1st part and different 2nd)
          ik4=iso3(inew(2))     !second daughter iso3 code
          call getResoId(ik1,ik2,ik3,ik4,ka)
c          if(itypii(1).eq.103)
c     .      print *,ityp(inew(1)),iso3(inew(1)),ityp(inew(2))
c     .             ,iso3(inew(2))
        elseif(nexit.eq.3)then
c          if(itypii(1).eq.107)
c           print *,itypii(1),ityp(inew(1)),iso3(inew(1)),ityp(inew(2))
c     .             ,iso3(inew(2)),ityp(inew(3)),iso3(inew(3))
          if(ik1.eq.107.and.ityp(inew(1)).eq.102)then !count 3 body decay for etap, omega and phi (there is 2 body channels too)
            if(iso3(inew(2)).ne.0)then
              ka=33001          !eta' -> eta-(pi+)-(pi-) (0.435)
            else
              ka=33002          !eta' -> eta-pi0-pi0 (0.215)
            endif
          elseif(ik1.eq.103.and.ityp(inew(1)).eq.101)then
            ka=22102    !omega -> 3 pi  (0.9)
          elseif(ik1.eq.109.and.ityp(inew(1)).eq.101)then
            ka=33103    !phi -> 3 pi  (0.02)
          endif
        endif
        if(ka.ne.0)then
        phinum=phinum+1
        if(phinum.gt.maxphi)call utstop('#####  ERROR 25092014 #####&')
        phiidres(phinum)=ka
        phipx(phinum)=pxii(1)
        phipy(phinum)=pyii(1)
        phipz(phinum)=pzii(1)
        phip0(phinum)=
     .   sqrt(massit(ik1)**2+pxii(1)**2+pyii(1)**2+pzii(1)**2)
        phir0(phinum)=tformii(1)   !1e-30  !unknown
        phirx(phinum)=1e-30
        phiry(phinum)=1e-30
        phirz(phinum)=1e-30
        phimass(phinum)=massit(ik1)
        phiscatt(phinum)=0
        do ii=1,nexit
           i=inew(ii)
           phiid(i)=phinum
        enddo
       endif
      endif

      !check decay products from marked resonances
      if(ninit.eq.2)then
        do i=1,2
          if(phiidii(i).gt.0)then
c            if(phiidres(phiidii(i))/100.eq.330)then
c            print*,'decay from marked:',indii(i)
c     . ,itypii(i),iso3ii(i),phiidii(i),phiidres(phiidii(i))
c            endif
            phiscatt(phiidii(i))=1
            do ii=1,nexit
              ixx=inew(ii)
              phiid(ixx)=0
              !write(*,*)'INFO 04062014 Phi decay prod scatt'
              !.         ,ixx,ityp(ixx),phiid(ixx)
            enddo
          endif
        enddo
      endif

      do kk=kmin,kmax
        if(kk.eq.2.and.lstr)write(ifch,'(a,$)')'    '
        do i=1,ninit
          ii=indii(i)
          io=itrace(1,ii)
          idpdgg=pdgid(itypii(i),iso3ii(i))
          if(ipdg.eq.1)then
          id=idpdgg
          else
          id=idtrafo('pdg','nxs',idpdgg)
          endif
          idxxx(i)=id
          if(ityii(i)/10.eq.3)ityhard=ityii(i)
          if(ityii(i)/10.lt.3.or.ityii(i)/10.eq.4
     .   .or.ityii(i)/10.eq.5)itysoft=ityii(i)
          if(ityii(i)/10.eq.6)itybulk=ityii(i)
          if(ityii(i)/10.eq.7)ityexot=ityii(i)
          etaii(i)=psrap(pxii(i),pyii(i),pzii(i))
          if(kk.eq.1)then
            !if(pxii(i)**2+pyii(i)**2.gt.20.**2)lstr=.true.
            !if(ityii(i).eq.0)lstr=.true.
            !if(tformii(i).lt.-0.5)lstr=.true.
            !if(abs(itypii(i)).eq.133)lstr=.true.
            !if(abs(itypii(i)).eq.134)lstr=.true.
            if(abs(itypii(i)).eq.49)lstr=.true.
          else
            !write(ifch,'(a,i6,i6,i3,f9.2,f6.2,$)')
            !.    ' HHHH',ii,io,ityii(i), tformii(i),sqrt(pxii(i)**2+pyii(i)**2)
            !if(lstr)write(ifch,'(i6,$)')id  !,ii,ityii(i)
            if(lstr)write(ifch,'(i6,i2,$)')itypii(i),iso3ii(i)
          endif
          !~~~~~~~~~~~~~~~~~~~~~~~
        enddo
        if(kk.eq.2.and.lstr)write(ifch,'(a,$)')'  ----->  '
        do ii=1,nexit
          i=inew(ii)
          idpdgg=pdgid(ityp(i),iso3(i))
          if(ipdg.eq.1)then
          id=idpdgg
          else
          id=idtrafo('pdg','nxs',idpdgg)
          endif
          if(ninit.eq.1)then
            itypart(i)=ityii(1)
            iorpart(i)=-999
          else
            iorpart(i)=0
            if(ityexot.gt.0)then
             itypart(i)=ityexot
            elseif(itybulk.gt.0)then
             itypart(i)=62
            elseif(ityhard.gt.0)then
             itypart(i)=ityhard
            elseif(itysoft.gt.0)then
             itypart(i)=itysoft
            else
             call utstop('***** ERROR 01072014 &')
            endif
          endif
          if(kk.eq.1)then
            !if(lstrange(id))lstr=.true.
            !if(px(i)**2+py(i)**2.gt.20.**2)lstr=.true.
            !if(itypart(i).eq.0)lstr=.true.
            !if(ishuuu.ge.2.and.abs(id).eq.230)lstr=.true.
            !if(ninit.ne.1.and.itypart(i).ge.60)lstr=.true.
            !if(abs(ityp(i)).eq.133)lstr=.true.
            !if(abs(ityp(i)).eq.134)lstr=.true.
            if(abs(ityp(i)).eq.49)lstr=.true. !Xi-
          else
            !if(ii.ge.3.and.mod(ii-1,2).eq.0)
            !.    write(ifch,'(/81a,$)')(' ',mm=1,81)
            !write(ifch,'(a,2i6,f9.2,f6.2,$)')' HHHH',i
            !.    ,itypart(i),r0(i),sqrt(px(i)**2+py(i)**2)
            !write(ifch,'(i6,2(e9.3,1x),$)')id,rx(i),ry(i)
            !write(ifch,'(2i6,$)')id,itypart(i)
            !if(lstr)write(ifch,'(i6,$)')id
            if(lstr)write(ifch,'(i6,i2,$)')ityp(i),iso3(i)
          endif
        enddo
        if(kk.eq.2.and.lstr)write(ifch,'(a)')'   '
      enddo
      call utprix('prirea',ishuuu,ishini,4)
      end
      function psrap(p1,p2,p3)
      real*8  pt,p1,p2,p3
      pt=sqrt(p2**2+p1**2)
      psrap=sign(999.d0,p3)
      if(p3.ne.0..and.pt.ne.0.)
     * psrap=sign(1.d0,p3)*
     * log((sqrt(p3**2+pt**2)+abs(p3))/pt)
      end

c----------------------------------------------------------------
      logical function lstrange(id)
c----------------------------------------------------------------
      include 'urqmd34/coms.f'
      integer id,ia
      logical lstr
      lstr=.false.
      ia=abs(id)
      if(ia.eq.2130.or.ia.eq.1230.or.ia.eq.1130.or.ia.eq.2230
     .         .or.ia.eq.1131.or.ia.eq.1231.or.ia.eq.2231)
     . lstr=.true.
      lstrange=lstr
      end

cc------------------------------------------------------------------------
      subroutine getchar(np,ch)
cc------------------------------------------------------------------------
      character*4 ch
      ch='    '
      if(np.le.9)then
       write(ch,'(a,i1)')np
      elseif(np.le.99)then
       write(ch,'(a,i2)')np
      elseif(np.le.999)then
       write(ch,'(a,i3)')np
      else
       ch='????'
      endif
      end

cc---------------------------------------------------------------------
          integer function nrapid(p3,p4)
cc---------------------------------------------------------------------
      real*8 p3,p4
      if(p4-p3.le.0.)then
      nrapid=1000
      elseif(p4+p3.le.0.)then
      nrapid=-1000
      else
      rap=0.5*log((p4+p3)/(p4-p3))
      nrapid=-19+int(rap+19.5)
      endif
      nrapid=max(nrapid,-10)
      nrapid=min(nrapid, 10)
      end

cc---------------------------------------------------------------------
          subroutine writeParticleTestFile(nn)
cc---------------------------------------------------------------------
      implicit none
      include "epos.inc"
      integer nn
      open(unit=89,file='uuu.data',status='old',access='append')
      write(89,*)'ptl'
      write(89,*)
     .   istptl(nn)
     .  ,ityptl(nn)
     .  ,idptl(nn)
     .  ,xorptl(1,nn)
     .  ,xorptl(2,nn)
     .  ,xorptl(3,nn)
     .  ,xorptl(4,nn)
     .  ,pptl(1,nn)
     .  ,pptl(2,nn)
     .  ,pptl(3,nn)
     .  ,pptl(4,nn)
     .  ,pptl(5,nn)
      close(89)
      end

      subroutine readTestFile(nn)
      implicit none
      include "epos.inc"
      character*80 line
      integer nn
      nn=0
   1  continue
      read(89,'(a)',end=2)line
      if(line(2:4).eq.'end')return
      nn=nn+1
      if(line(2:4).ne.'ptl')stop'\n\n STOP in readTestFile \n\n'
      read(89,*)
     .   istptl(nn)
     .  ,ityptl(nn)
     .  ,idptl(nn)
     .  ,xorptl(1,nn)
     .  ,xorptl(2,nn)
     .  ,xorptl(3,nn)
     .  ,xorptl(4,nn)
     .  ,pptl(1,nn)
     .  ,pptl(2,nn)
     .  ,pptl(3,nn)
     .  ,pptl(4,nn)
     .  ,pptl(5,nn)
      goto 1
   2  stop'\n\n STOP in readTestFile: EoF \n\n'
      end

      subroutine openNewTestFile
      data ncnttf/0/
      save ncnttf
      ncnttf=ncnttf+1
      if(ncnttf.eq.1)then
      open(unit=89,file='uuu.data',status='new')
      write(89,*)'test'
      close(89)
      endif
      end

      subroutine openOldTestFile
      character*80 line
      data ncnttf/0/
      save ncnttf
      ncnttf=ncnttf+1
      if(ncnttf.eq.1)then
      open(unit=89,file='uuu.data',status='old')
      else
      rewind(89)
      endif
      read(89,'(a)')line
      if(line(2:5).ne.'test')stop'\n\n STOP in opeOldTestFile\n\n'
c      endif
      end

      subroutine writeEndEventTestFile
      open(unit=89,file='uuu.data',status='old',access='append')
      write(89,*)'end'
      close(89)
      end


c###################################################################################
c###################################################################################
c#######
c#######
c#######       random number stuff
c#######
c#######
c###################################################################################
c###################################################################################


      double precision function ranff(idummy)
      jdummy=idummy
      ranff = rangen()
      return
      end
      subroutine sseed(ranseed)
      integer ranseed
      dummyseed=ranseed
      return
      end


c###################################################################################
c###################################################################################
c#######
c#######
c#######       tables.dat stuff
c#######
c#######
c###################################################################################
c###################################################################################

c----------------------------------------------------------
      subroutine loadwtab (io)
c----------------------------------------------------------
c   output   : information in common-block comwid.f
c
c     load the tabulated branching ratios from disk
c
c----------------------------------------------------------

      implicit none

      include 'urqmd34/comres.f'
      include 'urqmd34/comwid.f'

      integer mxho
      integer       nfnch,nfnhi,nfndt,nfnii,nfnid,nfnie,nfnrj,nfnmt
     *,nfngrv,nfncp,nfnnx,nfncs,nfndr,nfnhpf

      parameter(mxho=10)
      character*500 fnch,fnhi,fndt,fnii,fnid,fnie,fnrj,fnmt
     *,fngrv,fncp,fnnx,fncs,fndr,fnhpf
      common/fname/  fnch, fnhi, fndt, fnii, fnid, fnie, fnrj, fnmt
     *,fngrv,fncp,fnnx,fncs,fndr,fnhpf
      common/nfname/nfnch,nfnhi,nfndt,nfnii,nfnid,nfnie,nfnrj,nfnmt
     *,nfngrv,nfncp,nfnnx,nfncs,nfndr,nfnhpf



c      integer ios, nsp, i, io, ver
      character*40 pwdcmd
      character*50 deftab
      character*8 defexe
      logical b
      integer ios, nsp,i, io, ver
c      character*50 deftab
c      character*8 defexe
c      logical b
      character*500 edir
      common/cdir/edir
      logical fex

      parameter (defexe='uqmd.exe')

      deftab(1:nfnhpf+1)=fnhpf(1:nfnhpf)//' '  !tables.dat including path

      b=io.eq.1

c set the name of the table
      tabname=edir(1:index(edir,' ')-1)
     .     //deftab(1:index(deftab,' '))
      inquire(file=tabname(1:index(tabname,' ')-1),exist=fex)
      if(.not.fex)then
        write (6,*) 'No file ',tabname(1:index(tabname,' ')-1)
        write(6,*) '    ==>   make it ... '
        call mkwtab
        write (6,'(a)') '      +---------------------------+ '
        write (6,'(a)') '      !  Table making finished.   ! '
        write (6,'(a)') '      !  Program will stop!!!     ! '
        write (6,'(a)') '      !  To continue: restart it. ! '
        write (6,'(a)') '      +---------------------------+ '
        stop
      endif

c
      if(b)write (6,*) 'Looking for the tabulated decay width...'
      open (unit=75,iostat=ios,file=tabname(1:index(tabname,' ')-1),
     .     status='old')



c set the defaultname of the file, containing the table
c      parameter (deftab='tables.dat', defexe='uqmd.exe')

c      b=io.eq.1


c get the name of the table from the environment variable
c      call getenv('URQMD_TAB',tabname)
c if it is empty, use the default name
c      if (tabname(1:4).eq.'    ') then
c         tabname=deftab
c      endif

c      if(b)write (6,*) 'Looking for the tabulated decay width...'
c open the table
c      open (unit=75,iostat=ios,file=tabname,form='unformatted',
c     .      status='old')
c if it fails ...
      if (ios.ne.0) then
         write(6,*)"error 1 in reading tables.dat"
         stop
c     close the file handle
c     close (unit=75, status='delete')
         if(b)write (6,*) 'No file:',tabname,'in this directory'
c     get the full path of the executable, ...
         call getenv('_',pwdcmd)
         write (6,*) 'pwd:',pwdcmd
c     extract the path
         i=max(index(pwdcmd,defexe),2)
         write (tabname,*) pwdcmd (1:i-1),deftab
         tabname=tabname(2:)
         if(b)write (6,*) 'Looking for ',tabname,'...'
c     and look for a table in the directory of the executable
         open (unit=75,iostat=ios,file=tabname,
     .        form='unformatted',status='old')
      endif
c     if the last 'open' command succeeds read the file
      if (ios.eq.0) then
         if(b)write (6,*) 'O.K.'
         if(b)write (6,*) 'reading...'
c     read all tables
         read (75,*) ver, nsp, tabx, fbtaby, pbtaby,
     .        fmtaby, pmtaby, bwbarnorm, bwmesnorm,
     .        tabxnd, frrtaby
      !--------------------------------------------------
      !  KW : I changed to formatted !!!!!!!!!!!!!!!
      !--------------------------------------------------
         if(b)write (6,*) 'version=',ver
c     if no errors occur ...
         if (ios.eq.0) then
            if(b)write (6,*) 'O.K.'
            wtabflg=3
c     check, if the version number is correct
            if (ver.eq.tabver) then
               if(b)write (6,*) 'tabver=',ver,'  O.K.'
            else
               write (6,*) 'wrong table!'
               write (6,*) 'tabver should be',tabver,',instead of',ver
               wtabflg=0
            endif
c     check, if the table has the correct 'widnsp'
            if (nsp.eq.widnsp) then
               if(b)write (6,*) 'widnsp=',nsp,'  O.K.'
            else
               write (6,*) 'wrong table!'
               write (6,*) 'widnsp should be',widnsp,', instead of',nsp
               wtabflg=0
            endif
c     if table is O.K. close file
            if (wtabflg.eq.3) then
               close (unit=75, status='keep')
c     otherwise ...
            else
c     delete the present table
               close (unit=75, status='delete')
               write(6,*)"error 2 in reading tables.dat"
               stop
               tabname=deftab
c     and calculate a new one
               call mkwtab
            endif
c     in case of read errors ...
         else
c     delete the present table
            close (unit=75, status='delete')
            write(6,*)"error 3 in reading tables.dat", ios
            stop
            write (6,*) 'Error while reading ',tabname

            tabname=deftab
c     and calculate a new one
            call mkwtab
         endif
c     in any other case ...
      else
         tabname=deftab
c calculate an new table
         call mkwtab
      endif

      return
      end





c      integer ios, nsp, io, ver
c      character*50 deftab
c      character*8 defexe
c      logical b
c      character*500 edir
c      common/cdir/edir

c set the defaultname of the file, containing the table
c      parameter (defexe='uqmd.exe')

c      deftab(1:nfnhpf+1)=fnhpf(1:nfnhpf)//' '  !tables.dat including path

c      b=io.eq.1

c set the name of the table

c      tabname=edir(1:index(edir,' ')-1)
c     .//deftab(1:index(deftab,' '))
c
c      if(b)write (6,*) 'Looking for the tabulated decay width...'
c      open (unit=75,iostat=ios,file=tabname(1:index(tabname,' ')-1),
c     .      status='old')
c      if (ios.eq.0) then
c         if(b)write (6,'(a)') ' (info) reading tables.dat'
c         read (75,*,err=99) ver, nsp, tabx, fbtaby, pbtaby,
c     .   fmtaby, pmtaby, bwbarnorm, bwmesnorm,
c     .         tabxnd, frrtaby
c        if(b)write (6,*) 'O.K.'
c        wtabflg=3
c         if (ver.eq.tabver) then
c            if(b)write (6,*) 'tabver=',ver,'  O.K.'
c         else
c            write (6,*) '(info) wrong tables.dat!'
c            write (6,*) '(info) tabver should be',tabver,'
c     .        ,instead of',ver
c            stop
c         endif
c         if (nsp.eq.widnsp) then
c            if(b)write (6,*) 'widnsp=',nsp,'  O.K.'
c         else
c            write (6,*) '(info) wrong table!'
c            write (6,*) '(info) widnsp should be',widnsp
c     .       ,', instead of',nsp
c               stop
c         endif
c         close (unit=75, status='keep')
c      else
c        write (6,*) 'No file ',tabname(1:index(tabname,' ')-1)
c        close (unit=75, status='delete')
c        call mkwtab
c        stop'\n\n mkwtab finished \n\n'
c      endif
c
c      return
c
c  99  continue
c      write (6,*) '(info) read error, reading from '
c     .,tabname(1:index(tabname,' ')-1)
c      stop
c      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine savewtab
c
c     Revision : 1.0
c
c     save the tabulated branching ratios to disk
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      include 'urqmd34/comres.f'
      include 'urqmd34/comwid.f'

      integer ios

      write (6,*) 'Writing new table...'
c try to generate a new file
      open (unit=75,iostat=ios,file=tabname(1:index(tabname,' ')-1),
     .      status='new')
c if it succedds ...
      if (ios.eq.0) then
c write the tables into the file
         write (75,*) tabver, widnsp, tabx, fbtaby,
     .        pbtaby, fmtaby, pmtaby, bwbarnorm, bwmesnorm,
     .      tabxnd, frrtaby
c otherwise complain
      else
         write (6,*) 'Error: ',tabname(1:index(tabname,' ')-1)
     .    ,'exists!'
      endif
c close the file
      close (unit=75, status='keep')

      return
      end


c      write (6,*) 'Writing new table...'
c try to generate a new file
c      open (unit=75,iostat=ios,file=tabname(1:index(tabname,' ')-1),
c     .      status='new')
c if it succedds ...
c      if (ios.eq.0) then
c write the tables into the file
c         write (75,*) tabver, widnsp, tabx, fbtaby,
c     .        pbtaby, fmtaby, pmtaby, bwbarnorm, bwmesnorm,
c     .            tabxnd, frrtaby
c otherwise complain
c      else
c         write (6,*) 'Error: ',tabname(1:index(tabname,' ')-1)
c     .    ,'exists!'
c      endif
c close the file
c      close (unit=75, status='keep')
c
c      return
c      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mkwtab
c
c     Revision : 1.0
c
coutput   : information in common-block comwid.f
c
c     tabulate the mass dependent branching ratios
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      include 'urqmd34/comres.f'
      include 'urqmd34/comwid.f'


      real*8 fwidth,m,first,last,delta,abl0,abln,mir,mminit,fbrancx
      real*8 massit,bran,smass,bwnorm,fppfit
      integer i,bchan,itp,isoit,cmin,cmax,i1,i2,i3,i4,ii1

      write (6,*) 'Generating table...'
c this indicates, that all tables are still empty
      wtabflg=0

c high precision splines from mintab to maxtab1
c lower precicision between maxtab1 and maxtab2

c now fill the x-values
c start with 'mintab'
      first=mintab
c 66 % of all fixpoints between mintab and maxtab1
c calculate the steps
      delta=(maxtab1-mintab)/((widnsp-1d0)*2d0/3d0)
      if (delta.le.0d0) then
         write(*,*)'(E) Please allow maxtab1>mintab in comwid'
         stop 137
      endif
c store the values into 'tabx'
      do 10 i=1,int(widnsp*2./3.)
         m=first+(i-1)*delta
         tabx(i)=m
 10   continue
c 33 % of all fixpoints with larger delta between maxtab1 and maxtab2
      delta=(maxtab2-maxtab1)/((widnsp-1d0)*1d0/3d0)
      if (delta.le.0d0) then
         write(*,*)'(E) Please allow maxtab2>maxtab1 in comwid'
         stop 137
      endif
c store the values into 'tabx'
        do 11 i=int(widnsp*2./3.)+1,widnsp
         m=maxtab1+(i-1-int(widnsp*2./3.))*delta
         tabx(i)=m
 11   continue

c now fill the y-values of the full branching ratios

c these are the first derivatives at the first an the last point
c of the interpolating function. a value greater than 1E30 signals the
c 'spline' routine to set the boundary condition for a natural spline
c with zero second derivative
      abl0=2D30
      abln=2D30

c loop over all baryons
      do 20 itp=minbar,maxbar
c loop over all x-values
         do 21 i=1,widnsp
c store the values ...
            fbtaby (i,itp,1)=fwidth(itp,isoit(itp),tabx(i))
 21      continue
c calculate the second derivate and store it in 'fbtaby(,,2)'
         call spline (tabx(1),fbtaby(1,itp,1),widnsp,abl0,abln,
     .                fbtaby(1,itp,2))
 20   continue
      write (6,*) '(1/7) ready.'

c loop over all mesons
      do 30 itp=minmes,maxmes
c loop over all x-values
         do 31 i=1,widnsp
c store the values ...
            fmtaby (i,itp,1)=fwidth(itp,isoit(itp),tabx(i))
 31      continue
c calculate the second derivate and store it in 'fmtaby(,,2)'
         call spline (tabx(1),fmtaby(1,itp,1),widnsp,abl0,abln,
     .                fmtaby(1,itp,2))
 30   continue
      write (6,*) '(2/7) ready.'

c the flag indicates, that now all full widths are tabulated
      wtabflg=1

c now fill the y-values of the partial branching ratios

c loop over all baryons
      do 40 itp=minbar,maxbar
c get the mass of this particle
         mir=massit(itp)
c get the range of possible decay channels
         call brange (itp, cmin, cmax)
c check, if there are any decay channels
         if (cmax.gt.0) then
c loop over all decay channels
            do 41 bchan=cmin,cmax
c now get the outgoing particles 'i1' and 'i2' for the channel 'j'
c 'bran' is the mass independent branching ratio (tabulated in blockres)
c 'bflag' indicates, if 'i1', 'i2' or both are broad
               call b3type (itp,bchan,bran,i1,i2,i3,i4)
c check, if decay is allowed

               smass=mminit(i2)
               if(i3.ne.0) smass=smass+mminit(i3)
               if(i4.ne.0) smass=smass+mminit(i4)

               if (bran.gt.1d-9.and.mir.gt.mminit(i1)+smass) then
c loop over all x-values
                  do 42 i=1,widnsp
c store the values
                     pbtaby(i,1,itp,bchan)=
     .                    fbrancx (bchan,itp,isoit(itp),tabx(i),
     .                    bran,i1,i2,i3,i4)
 42               continue
c calculate the second derivate and store it in 'pbtaby(,2,,)'
                  call spline (tabx(1),pbtaby(1,1,itp,bchan),widnsp,
     .                         abl0,abln,pbtaby(1,2,itp,bchan))
               end if
 41         continue
         end if
 40   continue
      write (6,*) '(3/7) ready.'

c loop over all mesons
      do 50 itp=minmes,maxmes
c get the mass of this particle
         mir=massit(itp)
c get the range of possible decay channels
         call brange (itp, cmin, cmax)
c check, if there are any decay channels
         if (cmax.gt.0) then
            do 51 bchan=cmin,cmax
c now get the outgoing particles 'i1' and 'i2' for the channel 'j'
c 'bran' is the mass independent branching ratio (tabulated in blockres)
c 'bflag' indicates, if 'i1', 'i2' or both are broad
               call b3type(itp,bchan,bran,i1,i2,i3,i4)
c!!!
               smass=mminit(i2)
               if(i3.ne.0) smass=smass+mminit(i3)
               if(i4.ne.0) smass=smass+mminit(i4)

               if (bran.gt.1d-9.and.mir.gt.mminit(i1)+smass) then
c loop over all x-values
                  do 52 i=1,widnsp
                     pmtaby(i,1,itp,bchan)=
     .                    fbrancx (bchan,itp,isoit(itp),tabx(i),
     .                    bran,i1,i2,i3,i4)
 52               continue
c calculate the second derivate and store it in 'pmtaby(,2,,)'
                  call spline (tabx(1),pmtaby(1,1,itp,bchan),widnsp,
     .                         abl0,abln,pmtaby(1,2,itp,bchan))
               end if
 51         continue
         end if
 50   continue

      write (6,*) '(4/7) ready.'


c calculate the norm integral of the Breit-Wigner functions
c   with mass dependent widths

c..baryons
        do 60 i=minbar,maxbar
           bwbarnorm(i)=bwnorm(i)
60      continue
      write (6,*) '(5/7) ready.'

c.. mesons
        do 61 i=minmes,maxmes
           bwmesnorm(i)=bwnorm(i)
61      continue
      write (6,*) '(6/7) ready.'

c now all branching ratios and B W - integrals are tabulated
      wtabflg=2

ce tabulate fppfit
c fill the x-values
c range of tabulated cross sections
      first=2d0*massit(nucleon)+massit(pimeson)
        last=maxtab1
c calculate the steps
c the energies are weighted quadratically
      delta=(last-first)/((widnsp-1)*2./3.)**2
c store the values into 'tabx'
c 66 % of all fixpoints between mintab and maxtab1
      do 69 i=1,int(widnsp*2./3.)
         m=first+(i-1)**2*delta
         tabxnd(i)=m
 69   continue
c 33 % of all fixpoints with larger, constant delta between maxtab1 and maxtab2
        delta=(maxtab2-last)/((widnsp-1)*1./3.)
        do 70 i=int(widnsp*2./3.)+1,widnsp
         m=maxtab1+(i-1-int(widnsp*2./3.))*delta
         tabxnd(i)=m
 70   continue


c.. all pp-exit channels
c loop over first out-particle N & D
        do 81 ii1=1,2
          if(ii1.eq.1)i1=minnuc
          if(ii1.eq.2)i1=mindel
c loop over second out-particle N(1440)..maxdel
          do 82 i2=minnuc+1,maxdel
c loop over all x-values
          do 83 i=1,widnsp
c store the values ...
              frrtaby(i,1,ii1,i2)=fppfit(99,tabxnd(i),i1,i2)
83          continue
c calculate the second derivate and store it in 'frrtaby(,,2)'
          call spline (tabxnd(1),frrtaby(1,1,ii1,i2),widnsp,abl0,abln,
     .                frrtaby(1,2,ii1,i2))
82        continue
81      continue


      write (6,*) '(7/7) ready.'

c pp cross sections are now tabulated
        wtabflg=3

c save the table on disk
      call savewtab

      return
      end



c######################################################################
c######################################################################
c#########
c#########                    utilities
c#########
c######################################################################
c######################################################################


ccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function splint (xa,ya,y2a,n,x)
c
c     Unit     : general infrastructure
c     Author   : (C) Copr. 1986-92 Numerical Recipes Software
c     Date     : 03/07/96
c     Revision : 1.1
c
ccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      include 'urqmd34/comres.f'
      include 'urqmd34/comwid.f'

      integer n
      integer k,khi,klo
      real*8 x,y,xa(n),y2a(n),ya(n)
      real*8 a,b,h

      klo=1
      khi=n
1     if (khi-klo.gt.1) then
         k=(khi+klo)/2d0
         if(xa(k).gt.x)then
            khi=k
         else
            klo=k
         endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) call utstop('\n\n STOP in splint: bad xa input\n\n')
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a*a*a-a)*y2a(klo)+
     .            (b*b*b-b)*y2a(khi))*(h*h)/6d0
      splint=y

      return
      end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function splintth (xa,ya,y2a,n,x,th)
c
c     Unit     : general infrastructure
c     Author   : (C) Copr. 1986-92 Numerical Recipes Software
c                modified my H. Weber
c     Date     : 03/07/96
c     Revision : 1.1
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     split routine with nice threshold behaviour for cross sections
c

      implicit none

      include 'urqmd34/comres.f'
      include 'urqmd34/comwid.f'

      integer n
      integer k,khi,klo
      real*8 x,y,xa(n),y2a(n),ya(n)
      real*8 a,b,h,th

      klo=1
      khi=n
1     if (khi-klo.gt.1) then
         k=(khi+klo)/2d0
         if(xa(k).gt.x)then
            khi=k
         else
            klo=k
         endif
         goto 1
      endif
      h=xa(khi)-xa(klo)
      if(h.eq.0.)call utstop('\n\n STOP in splint: bad xa input \n\n')
      if (xa(khi).lt.(th+2*h)) then
c linear approximation close to threshold (within 2h)
         splintth=ya(khi)*(x-th)/(xa(khi)-th)
      else
         a=(xa(khi)-x)/h
         b=(x-xa(klo))/h
         y=a*ya(klo)+b*ya(khi)+((a*a*a-a)*y2a(klo)+
     .        (b*b*b-b)*y2a(khi))*(h*h)/6d0
         splintth=y
      endif

      return
      end

c######################################################################
c######################################################################
c
c    JPSI
c
c######################################################################
c######################################################################

      subroutine getJpsi(engy,px,py,pz,en,am)
      !jpsi see epos_uj.f
      engy=0
      px=0
      py=0
      pz=0
      en=0
      am=0
      end


