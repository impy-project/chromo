c $Id: colltab.f,v 1.5 1999/11/24 19:47:48 ssoff Exp $
c
cdes  This file contains the  uqmd collision tables
c
      integer ncollmax
      parameter (ncollmax = 20000)
      integer nct,actcol,nsav,apt
      real*8 cttime(0:ncollmax),ctsqrts(ncollmax),ctsigtot(ncollmax)
      real*8 ctcolfluc(ncollmax)
      logical ctvalid(ncollmax)
      real*8 tmin
      integer cti1(ncollmax),cti2(ncollmax),ctsav(ncollmax)
c      integer updi1(ncollmax),updi2(ncollmax)
c
c     cttime  : collision time
c     ctsqrts : $sqrt{s}$ of collision
c     ctsigtot: total cross section in mbarn
c     tmin    : paramteter for {\tt collupd}
c     cti1    : index of particle 1
c     cti2    : index of particle 2
c     nct     : number of collisions in the table
c     actcol  : current collision
c     ctvalid : tag whether collision is {\em true} or {\em false}
c     ctsav   : list of particles which lost their collision partner
c     nsav    : number of entries in {\tt ctsav}
c     apt     : mass of first particle/composite in the part. arrays 

      common /colltab/cttime,ctsqrts,ctsigtot,tmin,cti1,cti2,nct,actcol,
     &     ctvalid,ctsav,nsav,apt,ctcolfluc
