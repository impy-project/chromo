c $Id: coms.f,v 1.22 2007/01/30 14:50:24 bleicher Exp $
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c standard common block for uQMD
c
cdes This file contains the standard commom blocks of UrQMD
c
      real*8 emnuc
      parameter (emnuc = 0.938)

      integer nmax, nspl
      real*8 hit_sphere
      parameter (nmax = 40000)
      parameter (nspl = 500)
      parameter (hit_sphere = 8.d0)

c     debug and validity range
      logical check, info, warn
      parameter (check=.true., info=.true., warn=.false.)

      integer Ap, At, Zp, Zt, npart, nbar, nmes, ctag
      integer nsteps,ranseed,event,eos,dectag,uid_cnt
      integer NHardRes,NSoftRes,NDecRes,NElColl,NBlColl
      logical success
      real*8  time,  acttime, bdist, ebeam, bimp,bmin,ecm
c 7 integer

      common /sys/ npart, nbar, nmes, ctag,nsteps,uid_cnt,
     +             ranseed,event,Ap,At,Zp,Zt,eos,dectag,
     +             NHardRes,NSoftRes,NDecRes,NElColl,NBlColl,
     +             success
      common /rsys/ time,acttime,bdist,bimp,bmin,ebeam,ecm

c     Ap       : projectile mass
c     Zp       : projectile charge
c     At       : target mass
c     Zt       : target charge
c     npart    : total number of particles
c     nbar     : number of baryons AND antibaryons
c     nmes     : number of mesons
c     ctag     : counter of All interactions (coll. and dec.)
c     nsteps   : number of timesteps
c     uid_cnt  : counter for assigning unique particle ID-tags
c     ranseed  : random number generator seed of event
c     event    : event counter
c     eos      : flag for the EoS chosen
c     dectag   : counter for decays
c     NHardRes : counter for resonance excitation in BB collisions
c     NSoftRes : counter for resonance excitation in MB collisions
c     NElColl  : counter for elastic collisions
c     NBlColl  : counter for Pauli-blocked collisions
c     time     : system time at beginning of timestep
c     acttime  : current system time
c     bdist    : maximum impact parameter (of event sample)
c     bimp     : actual impact parameter
c     bmin     : minimum impact parameter (of event sample)
c     ebeam    : incident beam energy (lab frame)
c     ecm      : initial projectile-hadron target-hadron c.m. energy 

      logical firstseed
      common /comseed/firstseed

      logical lsct(nmax), 
     +        logSky, logYuk, logCb, logPau
      common /logic/ lsct, logSky, logYuk, logCb, logPau
c 2*nmax*nmax logical

      integer spin(nmax),ncoll(nmax),charge(nmax),strid(nmax),
     +        ityp(nmax),lstcoll(nmax),iso3(nmax),origin(nmax),uid(nmax)
c 6*nmax integer



      real*8 eps, er0, pi, rho0,hqc
      parameter (eps  = 1.0E-12)
      parameter (er0  = 1.128379167)
      parameter (pi   = 3.1415926535)
      parameter (rho0 = 0.16)
      parameter (hqc  = 0.197327)

c IMPORTANT: when you change the version number please change also
c            the versiontag in blockres.f !
      integer version, laires
      character*8 versiontxt

      parameter ( version = 30400)
      parameter ( laires  = 30400)
      parameter ( versiontxt = '3.4' )

c MD temporary arrays
       real*8 r0_t(nmax), rx_t(nmax), ry_t(nmax), rz_t(nmax)

       common/mdprop/ r0_t, rx_t, ry_t, rz_t

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

      common/isys/spin,ncoll,charge,ityp,lstcoll,iso3,origin,
     +            uid
      common /coor/ r0, rx, ry, rz, p0, px, py, pz, fmass, rww, dectime
      common /frag/ tform, xtotfac

c     spin    : particle spin
c     ncoll   : particle number of collisions
c     charge  : particle charge
c     ityp    : particle ID
c     lstcoll : tag of last interaction of particle
c     iso3    : $2 \cdot I_3$ of particle
c     origin  : ID of last interaction of particle
c     uid     : unique particle ID tag
c     r0      : particle time
c     rx      : $x$ coordinate
c     ry      : $y$ coordinate
c     rz      : $z$ coordinate
c     p0      : particle energy
c     px      : $p_x$ momentum
c     py      : $p_y$ momentum
c     pz      : $p_z$ momentum
c     fmass   : mass of particle
c     rww     : ??
c     dectime : decay time of particle
c     tform   : formation time of particle
c     xtotfac : cross section scaling factor during formation time
c
      common /aios/ airx, airy, airz, aipx, aipy, aipz,
     +              aorx, aory, aorz, aopx, aopy, aopz

      common /pots/ Cb0, Yuk0, Pau0, Sky20, Sky30, gamSky, 
     +              gamYuk, drPau, dpPau, gw, sgw, delr, fdel,
     +              dt,da, db,dtimestep


c spectator arrays
        integer smax
        parameter(smax=500)
        real*8 r0s(smax), rxs(smax), rys(smax), rzs(smax),
     +         p0s(smax), pxs(smax), pys(smax), pzs(smax),
     +         sfmass(smax)
        

        integer sspin(smax), scharge(smax), sityp(smax), siso3(smax),
     +          suid(smax)

        integer nspec

        common /scoor/ r0s, rxs, rys, rzs, p0s, pxs ,pys, pzs, sfmass

        common /sisys/ sspin, scharge, sityp, siso3, suid

        common /ssys/ nspec

c     sspin    : spectator particle spin
c     scharge  : spectator particle charge
c     sityp    : spectator particle ID
c     siso3    : $2 \cdot I_3$ of spectator particle
c     r0s      : spectator particle time
c     rxs      : spectator $x$ coordinate
c     rys      : spectator $y$ coordinate
c     rzs      : spectator $z$ coordinate
c     p0s      : spectator particle energy
c     pxs      : spectator $p_x$ momentum
c     pys      : spectator $p_y$ momentum
c     pzs      : spectator $p_z$ momentum
c     sfmass   : mass of spectator particle


c
        real*8 p0td(2,nmax),pxtd(2,nmax),pytd(2,nmax),pztd(2,nmax),
     +         fmasstd(2,nmax)
        integer ityptd(2,nmax),iso3td(2,nmax)
        integer itypt(2),uidt(2),origint(2),iso3t(2)

        common /rtdelay/p0td,pxtd,pytd,pztd,fmasstd
        common /itdelay/ityptd,iso3td
        common /svinfo/itypt,uidt,origint,iso3t

c     p0td    : energy of parent particles of resonace (DP formalism)
c     pxtd    : $p_x$ of parent particles of resonace (DP formalism)
c     pytd    : $p_y$ of parent particles of resonace (DP formalism)
c     pztd    : $p_z$ of parent particles of resonace (DP formalism)
c     fmasstd : mass of parent particles of resonace (DP formalism)
c     ityptd  : ID of parent particles of resonace (DP formalism)
c     iso3td  : $2\cdot I_3$ of parent particles of resonace (DP formalism)

        real*8 ffermpx(nmax), ffermpy(nmax), ffermpz(nmax)
        real*8 peq1, peq2
        common /ffermi/ ffermpx, ffermpy, ffermpz
        common /peq/ peq1,peq2
        
c     ffermpx  : fermi momentum in $x$ direction
c     ffermpy  : fermi momentum in $y$ direction
c     ffermpz  : fermi momentum in $z$ direction

c..type of ccbar
        integer cross(60)
        common cross
