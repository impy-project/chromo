c $Id: saveinfo.f,v 1.7 2002/05/03 00:31:19 weber Exp $
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      subroutine saveinfo(ind,itag)
c
c     Revision : 1.0
c
cinput ind:  particle ID
cinput itag: slot for particle to stored (>0) or extracted (<0) from
c
c     This subroutine stores the information of the {\tt ind} slot in
c     the particle arrays (necessary in case of pauli-blocked collisions)
c     The absolute value of {\tt itag}  indicates the storage slot to be used.
c     For positive values of {\tt itag} the information is stored, for negative
c     values it is restored. Currently two slots are available.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
      implicit none
      include 'coms.f'
c
      integer ind,itag,islot,j
      integer lstcollt(2),ncollt(2)
      integer charget(2),spint(2)
c     origint(2),uidt(2)
      real*8 r0t(2),rxt(2),ryt(2),rzt(2),
     p      r0tt(2),rxtt(2),rytt(2),rztt(2),
     @     p0t(2),pxt(2),pyt(2),pzt(2),fmasst(2),tdectime(2),
     @     xtotfact(2),tformt(2),p0tdt(2,2),pxtdt(2,2),pytdt(2,2),
     @     pztdt(2,2),fmasstdt(2,2)
      integer ityptdt(2,2),iso3tdt(2,2)
      save
c
      if(ind.eq.0) return
c 
      islot=abs(itag)

c     save particle
      if(itag.gt.0) then
         r0t(islot)=r0(ind)
         rxt(islot)=rx(ind)
         ryt(islot)=ry(ind)
         rzt(islot)=rz(ind)
cpot
         r0tt(islot)=r0_t(ind)
         rxtt(islot)=rx_t(ind)
         rytt(islot)=ry_t(ind)
         rztt(islot)=rz_t(ind)

         p0t(islot)=p0(ind)
         pxt(islot)=px(ind)
         pyt(islot)=py(ind)
         pzt(islot)=pz(ind)
         fmasst(islot)=fmass(ind)
         itypt(islot)=ityp(ind)
         iso3t(islot)=iso3(ind)
         ncollt(islot)=ncoll(ind)
         lstcollt(islot)=lstcoll(ind)
         origint(islot)=origin(ind)
         charget(islot)=charge(ind)
         spint(islot)=spin(ind)
         tdectime(islot)=dectime(ind)
         uidt(islot)=uid(ind)
         xtotfact(islot)=xtotfac(ind)
         tformt(islot)=tform(ind)
ctd
         do 11 j=1,2
            p0tdt(j,islot)=p0td(j,ind)
            pxtdt(j,islot)=pxtd(j,ind)
            pytdt(j,islot)=pytd(j,ind)
            pztdt(j,islot)=pztd(j,ind)
            fmasstdt(j,islot)=fmasstd(j,ind)
            ityptdt(j,islot)=ityptd(j,ind)
            iso3tdt(j,islot)=iso3td(j,ind)
 11      continue

c         ...
      elseif(itag.lt.0)then
c     restore particle
         r0(ind)=r0t(islot)
         rx(ind)=rxt(islot)
         ry(ind)=ryt(islot)
         rz(ind)=rzt(islot)
cpot
         r0_t(ind)=r0tt(islot)
         rx_t(ind)=rxtt(islot)
         ry_t(ind)=rytt(islot)
         rz_t(ind)=rztt(islot)

         p0(ind)=p0t(islot)
         px(ind)=pxt(islot)
         py(ind)=pyt(islot)
         pz(ind)=pzt(islot)
         fmass(ind)=fmasst(islot)
         ityp(ind)=itypt(islot)
         iso3(ind)=iso3t(islot)
         ncoll(ind)=ncollt(islot)
         lstcoll(ind)=lstcollt(islot)
         origin(ind)=origint(islot)
         charge(ind)=charget(islot)
         spin(ind)=spint(islot)
         dectime(ind)=tdectime(islot)
         uid(ind)=uidt(islot)
         xtotfac(ind)=xtotfact(islot)
         tform(ind)=tformt(islot)
ctd
         do 12 j=1,2
            p0td(j,ind)=p0tdt(j,islot)
            pxtd(j,ind)=pxtdt(j,islot)
            pytd(j,ind)=pytdt(j,islot)
            pztd(j,ind)=pztdt(j,islot)
            fmasstd(j,ind)=fmasstdt(j,islot)
            ityptd(j,ind)=ityptdt(j,islot)
            iso3td(j,ind)=iso3tdt(j,islot)
 12      continue

      endif
      return
      end
