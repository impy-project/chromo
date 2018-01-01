c $Id: delpart.f,v 1.5 2000/01/12 16:02:34 bass Exp $
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine delpart(index)
c
c     Revision : 1.0 
c
cinput index : index of particle to delete
c
c     This subroutine deletes the entry of particle {\tt index} in all 
c     particle arrays.
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      include 'coms.f'
      include 'comres.f'
      include 'newpart.f'
      include 'freezeout.f'
      integer index,i,j,ind

      ind=index
      if(ind.lt.1.or.ind.gt.npart) then
         write(6,*)'***(E) delpart: ind out of bounds ',ind,npart
         stop 137
      endif

c     update nbar and nmes counters
      if(iabs(ityp(ind)).le.maxbar) then
          nbar=nbar-1
      elseif(iabs(ityp(ind)).ge.minmes)then
          nmes=nmes-1
      endif

c     delete slot
      do 10 i=ind+1,npart
         r0(i-1)=r0(i)
         rx(i-1)=rx(i)
         ry(i-1)=ry(i)
         rz(i-1)=rz(i)
         p0(i-1)=p0(i)
         px(i-1)=px(i)
         py(i-1)=py(i)
         pz(i-1)=pz(i)
         fmass(i-1)=fmass(i)
         ityp(i-1)=ityp(i)
         iso3(i-1)=iso3(i)
         ncoll(i-1)=ncoll(i)
         lstcoll(i-1)=lstcoll(i)
         charge(i-1)=charge(i)
         spin(i-1)=spin(i)
         dectime(i-1)=dectime(i)
         tform(i-1)=tform(i)
         xtotfac(i-1)=xtotfac(i)
         origin(i-1)=origin(i)
         uid(i-1)=uid(i)
         frr0(i-1)=frr0(i)
         frrx(i-1)=frrx(i)
         frry(i-1)=frry(i)
         frrz(i-1)=frrz(i)
         frp0(i-1)=frp0(i)
         frpx(i-1)=frpx(i)
         frpy(i-1)=frpy(i)
         frpz(i-1)=frpz(i)
         ffermpx(i-1)=ffermpx(i)
         ffermpy(i-1)=ffermpy(i)
         ffermpz(i-1)=ffermpz(i)
         r0_t(i-1)=r0_t(i)
         rx_t(i-1)=rx_t(i)
         ry_t(i-1)=ry_t(i)
         rz_t(i-1)=rz_t(i)
ctd
         do 11 j=1,2
            p0td(j,i-1)=p0td(j,i)
            pxtd(j,i-1)=pxtd(j,i)
            pytd(j,i-1)=pytd(j,i)
            pztd(j,i-1)=pztd(j,i)
            fmasstd(j,i-1)=fmasstd(j,i)
            ityptd(j,i-1)=ityptd(j,i)
            iso3td(j,i-1)=iso3td(j,i)
 11      continue

c            ...
 10      continue
        npart=npart-1

c     update collision tables
c     and scan for particles which have lost their collision partner
      call scantab(ind,-1)

c     update pointer array for new/scattered particles
      do 20 i=1,nexit
         if(inew(i).gt.ind) then
            inew(i)=inew(i)-1
         elseif(inew(i).eq.ind) then
            inew(i)=0
         endif
 20   continue
      return
      end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine adspec(index)
c
c     Revision : 1.0
c
cinput index : index of particle to delete
c
c     This subroutine deletes the entry of particle {\tt index} in all 
c     particle arrays and writes it to file 14 and 16 by the call 
c     of specout
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      include 'coms.f'
      include 'comres.f'
      include 'newpart.f'
      integer index,i,ind

      ind=index
      if(ind.lt.1.or.ind.gt.npart) then
         write(6,*)'***(E) adspec: ind out of bounds ',ind,npart
         stop 137
      endif

cernst fill spectator-arrays
      nspec=nspec+1
      r0s(nspec)=r0(ind)
      rxs(nspec)=rx(ind)
      rys(nspec)=ry(ind)
      rzs(nspec)=rz(ind)
      p0s(nspec)=p0(ind)
      pxs(nspec)=px(ind)
      pys(nspec)=py(ind)
      pzs(nspec)=pz(ind)
      sfmass(nspec)=fmass(ind)
      sspin(nspec)=spin(ind)
      scharge(nspec)=charge(ind)
      sityp(nspec)=ityp(ind)
      siso3(nspec)=iso3(ind)
      suid(nspec)=uid(ind)

c     update nbar and nmes counters
      if(iabs(ityp(ind)).le.maxbar) then
          nbar=nbar-1
      elseif(iabs(ityp(ind)).ge.minmes)then
          nmes=nmes-1
      endif

      i=ind
      call specout(i,14)
      call specout(i,16)

      do 10 i=ind+1,npart
         r0(i-1)=r0(i)
         rx(i-1)=rx(i)
         ry(i-1)=ry(i)
         rz(i-1)=rz(i)
         p0(i-1)=p0(i)
         px(i-1)=px(i)
         py(i-1)=py(i)
         pz(i-1)=pz(i)
         fmass(i-1)=fmass(i)
         ityp(i-1)=ityp(i)
         iso3(i-1)=iso3(i)
         ncoll(i-1)=ncoll(i)
         lstcoll(i-1)=lstcoll(i)
         charge(i-1)=charge(i)
         spin(i-1)=spin(i)
         dectime(i-1)=dectime(i)
         tform(i-1)=tform(i)
         xtotfac(i-1)=xtotfac(i)
         uid(i-1)=uid(i)
c            ...
 10      continue
        npart=npart-1

c     update collision tables
      call scantab(ind,-1)

c     update pointer array for new/scattered particles
      do 20 i=1,nexit
         if(inew(i).gt.ind) then
            inew(i)=inew(i)-1
         elseif(inew(i).eq.ind) then
            inew(i)=0
         endif
 20   continue
      return
      end


C####C##1#########2#########3#########4#########5#########6#########7##
      subroutine rmspec(bpro,btar)
cccccCcc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      implicit none
      include 'coms.f'
      include 'options.f'
      integer i,n,nspc
      real*8 sr1,sr2,srp,srt,bpro,btar,roff
      real*8 nucrad

      n=npart
      nspc=0
      roff=CTParam(35)
      srp=(roff+nucrad(Ap))**2
      srt=(roff+nucrad(At))**2
      if(n.le.2)return
      i=n
 108  continue

        sr1=ry(i)**2+(rx(i)-bpro)**2
        sr2=ry(i)**2+(rx(i)-btar)**2
        if((sr1.gt.srp.or.sr2.gt.srt).and.ityp(i).eq.1)then
           write(6,*)'rmspec: ',rx(i),ry(i),i
           call adspec(i)
        end if
        i=i-1
        if(i.ge.1)goto 108

      return
      end


