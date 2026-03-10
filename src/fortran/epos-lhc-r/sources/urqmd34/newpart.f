c $Id: newpart.f,v 1.8 2007/01/30 14:50:25 bleicher Exp $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     include-file newpart
c
cdes     this file contains arrays for scattered or new created particles
cdes     and is used to communicate between the different routines which
cdes     are involved in hadling the kinematics
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer mprt,oprt
c maximum number of new particles:
      parameter(mprt=1000) ! maximum number of produced particles
      parameter(oprt=2)   ! maximum number of incoming particles
c pslot : slots of incoming particles
c itypnew: ityps of new particles
c i3new: $2*I_3$ of new particles
c inew: array-indices of new particles (will be assigned in scatter-routines)
c nexit: number of particles in exit-channel
c iline: tag for out-channel process
c nstring1 : number of particles in string 1 
c nstring2 : number of particles in string 2
c pslot : slots of incoming particles
c itot: $2*I_{tot}$ of new particles (will be assigned in scatter-routines)
c pnew(5,mprt) : momenta, energy and mass of produced particles
c     pnew(1,*) = px
c     pnew(2,*) = py
c     pnew(3,*) = pz
c     pnew(4,*) = e
c     pnew(5,*) = mass
c xnew(4,mprt) : locations and time of produced particles
c     xnew(1,*) = rx
c     xnew(2,*) = ry
c     xnew(3,*) = rz
c     xnew(4,*) = r0
c pold(5,oprt) : momenta, defined like pnew
c itypold : incoming itypes
c iso3old : incoming iso3's
c xtotfacold: xtotfacs of incoming particles
c mstring() masses of strings 1 and 2 (or particles 1 and 2)
c leadfac(mprt): (1-leadfac) is the factor, by which the total cross 
c                section of a hadron is multiplied within its formation 
c                time (.ne.1 only for leading hadrons)
      integer itypnew(mprt),i3new(mprt),itot(mprt),inew(mprt),nexit
      integer nstring1, nstring2,iline,itypold(oprt),iso3old(oprt)
      integer pslot(oprt)
      real*8 pnew(5,mprt),xnew(4,mprt),mstring(2),leadfac(mprt)
      real*8 pold(5,oprt),xtotfacold(oprt)


c relative velocity/between comp. frame and two particle rest frame
c is betax, betay, betaz (needed for lotrans)
c momentum vector in two particle restframe is p0nn,pxnn,pynn,pznn

      real*8 betax,betay,betaz,p0nn,pxnn,pynn,pznn,pnn,pnnout
      
      common /inewpart/ itypnew,i3new,itot,inew,nexit,iline,
     &                  pslot,nstring1,nstring2,itypold,iso3old
      common /rnewpart/ pnew,xnew,betax,betay,betaz,pold,
     &                  p0nn,pxnn,pynn,pznn,pnn,mstring,pnnout,
     &                  xtotfacold
      common /fnewpart/ leadfac
