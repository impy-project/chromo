c $Id: addpart.f,v 1.4 2000/01/12 16:02:32 bass Exp $
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      subroutine addpart(index)
c
c     Revision : 1.0
c
cinput index : index for slot to create in particle arrays
c
c This subroutine creates an entry for a particle with index {\tt index} in all 
c particle arrays. 
c \danger{The {\tt nbar} and {\tt nmes} counters must be updated by the 
c calling routine !}
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      include 'coms.f'
      include 'newpart.f'
      include 'freezeout.f'
      integer ind,i,j,index
      ind=index

c     now shift vectors upwards
      do 10 i=npart,ind,-1
         r0(i+1)=r0(i)
         rx(i+1)=rx(i)
         ry(i+1)=ry(i)
         rz(i+1)=rz(i)

         p0(i+1)=p0(i)
         px(i+1)=px(i)
         py(i+1)=py(i)
         pz(i+1)=pz(i)
         fmass(i+1)=fmass(i)
         ityp(i+1)=ityp(i)
         iso3(i+1)=iso3(i)
         ncoll(i+1)=ncoll(i)
         lstcoll(i+1)=lstcoll(i)
         charge(i+1)=charge(i)
         spin(i+1)=spin(i)
         dectime(i+1)=dectime(i)
         tform(i+1)=tform(i)
         xtotfac(i+1)=xtotfac(i)
         origin(i+1)=origin(i)
         uid(i+1)=uid(i)
         frr0(i+1)=frr0(i)
         frrx(i+1)=frrx(i)
         frry(i+1)=frry(i)
         frrz(i+1)=frrz(i)
         frp0(i+1)=frp0(i)
         frpx(i+1)=frpx(i)
         frpy(i+1)=frpy(i)
         frpz(i+1)=frpz(i)
         ffermpx(i+1)=ffermpx(i)
         ffermpy(i+1)=ffermpy(i)
         ffermpz(i+1)=ffermpz(i)

         r0_t(i+1)=r0_t(i)
         rx_t(i+1)=rx_t(i)
         ry_t(i+1)=ry_t(i)
         rz_t(i+1)=rz_t(i)

         do 11 j=1,2
            p0td(j,i+1)=p0td(j,i)
            pxtd(j,i+1)=pxtd(j,i)
            pytd(j,i+1)=pytd(j,i)
            pztd(j,i+1)=pztd(j,i)
            fmasstd(j,i+1)=fmasstd(j,i)
            ityptd(j,i+1)=ityptd(j,i)
            iso3td(j,i+1)=iso3td(j,i)
 11      continue


 10      continue

c     increment npart
         npart=npart+1
c     
         if(npart.ge.nmax) then  
            write(6,*)'*** error in addpart:too much particles>',nmax
            write(6,*)' -> increase nmax in coms.f '
            call file14out(999)
            stop 137
         endif

c     initialize new slot
         r0(ind)=0.0
         rx(ind)=0.0
         ry(ind)=0.0
         rz(ind)=0.0

         p0(ind)=0.0
         px(ind)=0.0
         py(ind)=0.0
         pz(ind)=0.0
         fmass(ind)=0.0
         ityp(ind)=0
         iso3(ind)=0
         ncoll(ind)=0
         lstcoll(ind)=0
         charge(ind)=0
         spin(ind)=-1
         dectime(ind)=0.d0
         tform(ind)=0.0d0
         xtotfac(ind)=1.0d0
         origin(ind)=0
         uid(ind)=0
         frr0(ind)=0.d0
         frrx(ind)=0.d0
         frry(ind)=0.d0
         frrz(ind)=0.d0
         frp0(ind)=0.d0
         frpx(ind)=0.d0
         frpy(ind)=0.d0
         frpz(ind)=0.d0
         ffermpx(ind)=0.d0
         ffermpy(ind)=0.d0
         ffermpz(ind)=0.d0
cpot
         r0_t(ind)=0.d0
         rx_t(ind)=0.d0
         ry_t(ind)=0.d0
         rz_t(ind)=0.d0
ctd
         do 12 j=1,2
            p0td(j,ind)=0.d0
            pxtd(j,ind)=0.d0
            pytd(j,ind)=0.d0
            pztd(j,ind)=0.d0
            fmasstd(j,ind)=0.d0
            ityptd(j,ind)=0
            iso3td(j,ind)=0
 12      continue

c     rescan collision table - all particle indices which have been
c     shifted upwards must be modified in the collision tables, too. 
      call scantab(ind,1)


      return
      end
