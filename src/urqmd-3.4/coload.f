c $Id: coload.f,v 1.12 2007/01/30 14:50:24 bleicher Exp $
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function sigtot(ind1,ind2,sqrts)
c
c     Revision : 1.0 
C
cinput ind1 : index of particle 1
cinput ind2 : index of particle 2
cinput sqrts: $\sqrt{s}$ of collision between part. 1 and 2
c
c     {\tt sigtot} returns the total cross section (in mbarn) for the collision
c     between the particles with the indices {\tt ind1} and {\tt ind2}.
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       implicit none
       include 'coms.f'
       include 'comres.f'
       include 'newpart.f'
c
       integer isigline
       integer ind1,ind2,collclass
       integer ityp1,ityp2,iso31,iso32
       real*8 sqrts,mminit
c for detailed-balance
       integer nCh,ii
       real*8 e1,e2,sigma


c     reset sigtot, store often needed values in scalars
       sigtot=0.d0 
       ityp1=ityp(ind1)
       ityp2=ityp(ind2)
       iso31=iso3(ind1)
       iso32=iso3(ind2)

c     now get collision class (line-number for sigmaLN array in blockres.f)
c     isigline classifies the collision (pp,pn,Delta N, Meson-Baryon etc)
       isigline=collclass(ityp1,iso31,ityp2,iso32)

       if(isigline.eq.0) return ! zero cross section for collclass=0

c     get pointer for total cross section (column #2 in SigmaLN array)
       iline=SigmaLn(2,1,isigline)

c     if zero, cross section is zero, return
       if(iline.eq.0) return 

       if(iline.gt.0) then
c
c     if gt zero we have a tabulated or parameterized total cross section,
c     all table-lookups and parametrizations are accessed via crossx
c
             call crossx(iline,sqrts,ityp1,iso31,
     &                         max(fmass(ind1),mminit(ityp1)),
     &                         ityp2,iso32,
     &                         max(fmass(ind2),mminit(ityp2)),
     &                         sigtot)

       else  
c
c     total cross section via sum of partial cross sections
c
c     get number of exit-channels:
          if (isigline.gt.maxreac) then 
             write (6,*) '4isigline: ',isigline
          endif
          nCh=SigmaLn(1,1,isigline)
c
c transformation quantities  into NN system for proper kinematics 
c (necessary for detailed balance cross sections)
c     first compute transformation betas
          e1=sqrt(fmass(ind1)**2+px(ind1)**2
     &         +py(ind1)**2+pz(ind1)**2)
          e2=sqrt(fmass(ind2)**2+px(ind2)**2
     &         +py(ind2)**2+pz(ind2)**2)
          betax=(px(ind1)+px(ind2))/(e1+e2)
          betay=(py(ind1)+py(ind2))/(e1+e2)
          betaz=(pz(ind1)+pz(ind2))/(e1+e2)
c     now transform momenta
          pxnn=px(ind1)
          pynn=py(ind1)
          pznn=pz(ind1)
          p0nn=e1
c     call to Lorentz transformation
          call rotbos(0d0,0d0,-betax,-betay,-betaz,
     &         pxnn,pynn,pznn,p0nn)
          pnn=sqrt(pxnn*pxnn+pynn*pynn+pznn*pznn)
c     end of transform part
c
c     loop over exit channels for sum of partial cross sections
c     partial cross sections start in column #3 of SigmaLN in blockres.f
          do 10 ii=3,nCh+2  
c     get pointer for partial cross section
             iline=SigmaLn(ii,1,isigline)
c     normal partial cross sections
             if(iline.gt.0) then
                   call crossx(iline,sqrts,ityp1,iso31,fmass(ind1),
     &                         ityp2,iso32,fmass(ind2),sigma)
             else 
c     detailed balance partial cross section (must be computed now)
c     (crossz delivers inelastic channel for SINGLE resonance)
                call crossz(iline,sqrts,ityp1,iso31,fmass(ind1),
     &                      ityp2,iso32,fmass(ind2),sigma)
             endif
c     perform summation of partial cross sections
             sigtot=sigtot+sigma
c     end of loop for partial cross sections
 10       continue
c     end of total / sum of partial cross sections if
       endif

c     return to caller
       return

       end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getnext(k)
c
c     Revision : 1.0
C
coutput k : index of next collison
c
c {\tt getnext} returns the index of the next collision or decay 
c to be performed.
c If no further collisons occur in the timestep, {\tt getnext} returns
c {\tt k}=0 .
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      implicit none
      integer k
      include 'coms.f'
      include 'options.f'
      include 'colltab.f'

c     shrink collision tables in case 80% load to counter possible overflow
       if(dble(nct)/ncollmax.ge.0.8) call collshrink

c
c     set k to current collision/decay (which has already been performed and
c     and is outdated by the time of this call)
      k = actcol

 1    continue
c     increment counter, next entry in table is now the current one
      k = k+1
c     if the end of the collision table has been reached, return with k=0
      if (k.gt.nct) then
         k = 0
         actcol = k
         return
      endif

c     if the current entry in the collision table is marked "F" - false -
c     (due to previous interaction of one of the collision partners)
c     then find new collision partners for the particle(s) via calls
c     to collupd
      if (.not.ctvalid(k)) then
         call collupd(cti1(k),1)
c     second call only if not a decay entry
         if(cti2(k).gt.0) call collupd(cti2(k),1)
c     if the current collision is now marked "T" - true - return
         if(ctvalid(k)) then
            actcol = k
            return
         endif
c     the current collision is still marked false, go to top of loop
c     (and increment counter)
         goto 1
      endif
c
c     the current entry is marked "T" - true - this is the next
c     collision to be perfomed
      actcol = k
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine collshrink
c
c     Revision : 1.0
c
c     This subroutine deletes all entries in the collision tables between
c     1 and {\tt actcol}-1. It's purpose is to counter a possible overflow
c     of the collision tables.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer i
      include 'colltab.f'

      do 104 i=actcol,nct
         cttime(1+i-actcol) = cttime(i)
         ctsqrts(1+i-actcol) = ctsqrts(i)
         ctsigtot(1+i-actcol) = ctsigtot(i)
         cti1(1+i-actcol) = cti1(i)
         cti2(1+i-actcol) = cti2(i)
         ctvalid(1+i-actcol) = ctvalid(i)
         ctcolfluc(1+i-actcol) = ctcolfluc(i)
 104  continue
c     recalculate number of collisions in tables
      nct=1+nct-actcol
c     reset pointer to current collision
      actcol=1

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine colload
c
c     Revision : 1.0 
c
c
c This routine fills the collision tables with all collisions and decays
c to be performed in the current timestep. Within the timestep, 
c particle propagation is assumed on straight lines. This routine
c actually only performs the outer of the double particle loop and
c calls {\tt collupd} for the inner loop.
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      integer i

      include 'coms.f'
      include 'colltab.f'

c     reset number of collisions in table and current collision pointer
      nct = 0
      actcol = 0

c reset all collision arrays
      do 10 i=1,ncollmax
         cttime(i)=0.d0
         ctsqrts(i)=0.d0
         ctsigtot(i)=0.d0
         cti1(i)=0
         cti2(i)=0
         ctvalid(i)=.false.
         ctcolfluc(i)=1.d0
 10   continue

c     outer loop over all particles
      do 20 i=1,npart
c     call collupd for inner loop, -1: inner loop only from i+1 to npart
         call collupd(i,-1)
 20   continue
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine collupd(i,all)
c
c     Revision : 1.0
C
cinput i   : particle to be checked for collison or decay
cinput all : flag for update mode
c
c {\tt collupd} checks whether particle {\tt i} will collide or decay
c in the time interval between the current time {\tt acttime} 
c and the end of the time step.
c {\tt collupd} uses the variable {\tt tmin} to find the {\bf earliest}
c interaction/decay of particle {\tt i} and store it in the
c collision arrays via a call to {\tt ctupdate}.
c For {\tt all}$>0$ all other particles from 1 to {\tt npart} are checked
c (necessary for update after a collision/decay), whereas for 
c {\tt all}$<0$ only the particles with the indices from {\tt i+1} to
c {\tt npart} are checked (for calls via {\tt colload}).
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none

      include 'coms.f'
      include 'options.f'
      include 'colltab.f'
      include 'boxinc.f'

      real*8 dst2
      integer i,all
      integer j,imin,jmin,j0,A
      integer wall
      integer stidx
      logical isstable
      real*8  tn, wallcoll
cc
      real*8 tcoll,sqrs,sigt,sqrts,sigtot
      real*8 smin,sigmin,sigfac

      real*8 colfaci,colfacj,colfac,colfacmin

c     number of "initial" particles in event
      A = At+Ap
c     initially tmin is set to the time-interval until the end of the timestep
c     tmin is then minimized to the time of the first interaction of 
c     particle i
c
c     acttime: current time
c     time: time at beginning of timestep
c     dtimestep: length of timestep
      tmin=dtimestep-acttime+time
chp... test times
c       write(*,*)'dtimestep,tmin,acttime,time',
c     & dtimestep,tmin,acttime,time
c
c     other information to be stored together with tmin 
      smin = 0 ! sqrt(s) of interaction
      sigmin = 0 ! total cross section of interaction
      imin = 0 ! index of first particle
      jmin = 0 ! index of second particle
c
c  check, if particle i is a resonance, that might decay within the remaining
c  part of the timestep. 
c  If so, than treat decay as collision (with particle 2= 0)
      if (dectime(i)-acttime.lt.tmin) then
         isstable = .false.
         do 132 stidx=1,nstable
            if (ityp(i).eq.stabvec(stidx)) then
c               write (6,*) 'no decay of particle ',ityp(i)
               isstable = .true.
            endif
 132     enddo
         if (.not.isstable) then 
            tmin = dectime(i) - acttime
            smin=fmass(i)
            imin = i
            jmin = 0
         endif
      endif


ccccccccccccccccccccccccccccccccc
c new walls selected if mbflag equals 2
      if ((mbox.gt.0).and.(mbflag.eq.2)) then

c acttime is the ACTUAL time 
c time is the time at the BEGINNING of the timestep
c tmin is being minimized RELATIVE to the beginning of the timestep

c tn is the relativ time to the next wall colision
        tn=wallcoll(i,wall)
        
c comparison wheather a particle decays before a wall colision
            if (tn.lt.tmin) then            
c set the time
                    tmin=tn
c set the particle number
                    imin=i
c set the wall
                    jmin=wall                   
            endif
        endif
cc


c default setting: loop does only go from i+1 to npart
         j0 = i+1
c in case of "update mode" let loop run starting from 1
         if (all.gt.0) j0 = 1
c
c  Now check, which is the earliest collision of particle i
c
         do 101 j=j0,npart

c check for some exclusion cases
            if (
c 1. avoid "self-interaction"
     &           i.ne.j
c 2. particles which have interacted with each other in the past
c    are only allowed to interact with each other if at least one
c    of them has had an interaction in between. For this, check which
c    was their last collision. Initial state particles have a lstcoll 
c    of 0. This matters only if the lstcoll entry of both particles
c    is 0. Thus it is enough to check for one of the entries being 0.
     &          .AND.
     &          ((lstcoll(i).ne.lstcoll(j)).or.(lstcoll(i).eq.0))
c 3. Particles within the projectile or target respectively are per
c    default only allowed to interact with each other in case they
c    have already had an interaction with a particle of the target
c    or projectile respectively. This can be turned off by setting
c    CTOption(6) to a nonzero value.
c    This rule of course must not not apply to produced particles.
c    In the case of a meson nucleus collision, projectile and
c    target may be swapped in the particle vectors (therefore the
c    use of Apt instead of Ap, because Apt has been then swapped
c    accordingly)
     &          .AND.
     &          (CToption(6).ne.0.or.i.gt.A.or.j.gt.A.or.
     &           ncoll(i)+ncoll(j).gt.0.or.
     &          (i.le.Apt.and.j.gt.Apt).or.(j.le.Apt.and.i.gt.Apt))
     &         ) then

c
c  determine time of minimal approach of particles i and j 
c  relative to current time
                  call nxtcoll(i,j,dst2,tcoll)
c
c  are the particles close enough - check cut off 
c (default: 250 mbarn, defined in coms.f)
               if (dst2.lt.hit_sphere) then
c does the collision occur in the current time step?
                  if (tcoll.gt.0.d0.and.tcoll.lt.tmin) then
c     get sqrt(s) and total cross section
                     sqrs = sqrts(i,j)
c
c reduced cross section for leading hadrons of string fragmentation
c within their formation time 
c the scaling factor is sigfac which is determined by the
c individual particles scaling factors xtotfac
                     if(tform(i).le.acttime+tcoll
     &                       .and.tform(j).le.acttime+tcoll) then
                       sigfac=1.d0
                     else if(tform(i).le.acttime+tcoll
     &                       .and.tform(j).gt.acttime+tcoll) then
                       sigfac=xtotfac(j)
                     else if(tform(j).le.acttime+tcoll
     &                       .and.tform(i).gt.acttime+tcoll) then
                       sigfac=xtotfac(i)
                     else 
                       sigfac=xtotfac(i)*xtotfac(j)
                     endif
c     get total cross section via call to sigtot and rescale
                     sigt = sigfac*sigtot(i,j,sqrs)
c     rescale sigtot due to color fluctuations
                     call colorfluc(ityp(i),ityp(j),sqrs,
     &                              colfaci,colfacj)
                     colfac=colfaci*colfacj
                     sigt = sigt*colfac
c     are we within the geometrical cross section
           if ((sigt/CTParam(67)).gt.max(1.d-8,(10.d0*pi*dst2))) then
c     this collision is to beat, now
                        tmin = tcoll 
                        smin = sqrs
                        sigmin = sigt
                        imin = i
                        jmin = j
                        colfacmin=colfac
                     endif
                  endif
               endif
            endif
 101     continue
         if (imin.gt.0) then
c  if we found something, then update table via call to ctupdate
c  (keep in mind: tmin is relative to actual time!)
            call ctupdate(imin,jmin,acttime+tmin,smin,sigmin,colfacmin)
            
            
chp test times
c       write(*,*) 'after colload: acttime, tmin',acttime,tmin            
c     in case of "full load mode" after every collision
c     only the first entry in the collision table is relevant
            if(CTOption(17).ne.0) nct=1
         endif
         
       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine ctupdate(i,j,t,s,sig,cfac)
c
cinput i : index of 1st colliding particle
cinput j : index of 2nd colliding particle (0 for decay)
cinput t : (absolut) time  of collision/decay
cinput s : $sqrt{s}$ (GeV) of collison
cinput sig : total cross section (mbarn) of collision
cinput cfac :  scaling factor for color fluctuation 
c
c This subroutine updates the collision arrays.     
c It determines the (chonologically) correct position for the new
c entry in the collision arrays, creates the respective slot and
c inserts the new entry (via a call to {\tt ctset}. 
c Then the arrays are scanned to tag the chonologically first collision 
c of particle {\tt i} or {\tt j} {\em true} and all subsequent ones
c as {\em false}.
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       implicit none

       integer i,j
       real*8 t,s,sig,cfac
       include 'colltab.f'
       integer k,tfound,ncoll1

       ncoll1 = nct + 1
       tfound = ncoll1
c loop over all collisions in table
       do 101 k=actcol+1,nct
c     get the correct position for the new entry (chronologically sorted)
          if (tfound.eq.ncoll1.and.t.le.cttime(k)) tfound = k
c     the above construct works as follows:
c     as long as t > cttime(k), the first term is true and teh
c     the second is false. At the correct position for the new
c     entry BOTH terms are true, for all later times the first
c     term is false in order not to overwrite the value of tfound
 101   continue

c make sure that collision is not already listed in the table
       if (.not.(t.eq.cttime(tfound).and.(i.eq.cti1(tfound).and.
     &      j.eq.cti2(tfound).or.i.eq.cti2(tfound).and.
     &      j.eq.cti1(tfound)))) then
c then create slot for  new entry 
             do 102 k=nct,tfound,-1
c.. take care of array bounds
                if ((k+1).ge.ncollmax) then
                   write(*,*)'(E) Collision table too small.'
                   write(*,*)'Increase ncollmax in colltab.f'
                   stop 137
                end if
                cttime(k+1) = cttime(k)
                ctsqrts(k+1) = ctsqrts(k)
                ctsigtot(k+1) = ctsigtot(k)
                cti1(k+1) = cti1(k)
                cti2(k+1) = cti2(k)
                ctvalid(k+1) = ctvalid(k)
                ctcolfluc(k+1) = ctcolfluc(k)
 102         continue
c     increment number of collisions/decays in table
          nct = ncoll1
       endif
c insert new entry into the created slot via call to ctset
       call ctset(tfound,i,j,t,s,sig,cfac)
c
c     only the chonologically FIRST collision of particle i is set to true
c
       do 103 k=actcol+1,nct,1
c     the newly found collision must be omitted in the following sequence
          if (k.eq.tfound) goto 103
c  is there already a collision with i or j ?
          if (cti1(k).eq.i.or.cti2(k).eq.i) then
c     if the other collision is at an earlier time (the table is time-ordered,
c     therefore k < tfound corresponds to an earlier time)
c     then set the new one to false or vice versa
             if (k.lt.tfound.and.ctvalid(k)) then
                ctvalid(tfound) = .false.
             else
                ctvalid(k) = .false.
             endif
          endif
c     do likewise for the second particle
          if (j.gt.0.and.(cti1(k).eq.j.or.cti2(k).eq.j)) then
             if (k.lt.tfound.and.ctvalid(k)) then
                ctvalid(tfound) = .false.
             else
                ctvalid(k) = .false.
             endif
          endif
 103   continue

       return
       end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctset(tfound,i,j,t,s,sig,cfac)
c
c     Revision : 1.0
C
cinput tfound : index in coll. array for coll. to be inserted 
cinput i : index of first colliding particle
cinput j : index of second colliing particle (0 for decay, negative: wall)
cinput t : (absolute) time of collision
cinput s : $sqrt{s}$ (GeV) of collison
cinput sig : total cross section (mbarn) of collsion 
cinput cfac :  scaling factor for color fluctuation
c
c {\tt ctset} enters the collision of particles {\tt i} and {\tt j}
c into the collision arrays at index {\tt tfound}.
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      implicit none

      integer tfound,i,j

      real*8 t,s,sig,cfac
      include 'colltab.f'
c      
      cttime(tfound) = t
      ctsqrts(tfound) = s
      ctsigtot(tfound) = sig
      cti1(tfound) = i
      cti2(tfound) = j
      ctvalid(tfound) = .true.
      ctcolfluc(tfound) = cfac
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       real*8 function sqrts(i,j)
c
c     Revision : 1.0
c
c     input: i,j : numbers of colliding particles 
c     output: $\sqrt{s}$ of collision as return value
c
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       implicit none

       include 'coms.f'
       integer i,j
       real*8 p10,p20
       p10 = sqrt((px(i)+ffermpx(i))**2
     +           +(py(i)+ffermpy(i))**2
     +           +(pz(i)+ffermpz(i))**2+fmass(i)**2)
       p20 = sqrt((px(j)+ffermpx(j))**2
     +           +(py(j)+ffermpy(j))**2
     +           +(pz(j)+ffermpz(j))**2+fmass(j)**2)
       sqrts = sqrt((p10+p20)**2
     +               -(px(i)+ffermpx(i)+px(j)+ffermpx(j))**2
     +               -(py(i)+ffermpy(i)+py(j)+ffermpy(j))**2
     +               -(pz(i)+ffermpz(i)+pz(j)+ffermpz(j))**2)
       return
       end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine nxtcoll(j,k,dst,dtauc)
c
c     Revision : 1.0
c
c input:  j,k   : indices of colliding particles 
coutput dst   : impact parameter squared
coutput dtauc : collisiontime in the computational system
c
c     {\tt nxtcoll}  is the heart  of the collision term. It determines
c     the time  in the  computional system, when the  collision between
c     j and k took or will take place ({\tt dtauc}). The squared impact 
c     parameter of the collision is returned in {\tt dst}. {\tt dst} is
c     independent of the computational system  in which the coordinates 
c     of j and k are given.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      implicit none
      integer j,k, i
      include 'coms.f'
      real*8 u1(0:3), u2(0:3), p1(0:3), p2(0:3)
      real*8 dp(0:3), du(0:3), bt_com(3), bt(3)
      real*8 du2, dp2, dudp, bt_du, bt_dp, btdr, bt2
      real*8 dst, gam_com, gam_com2, bt_com2, dtauc

c
c do some assignments first
c
      u1(0) = r0(j)
      u1(1) = rx(j)
      u1(2) = ry(j)
      u1(3) = rz(j)
      u2(0) = r0(k)
      u2(1) = rx(k)
      u2(2) = ry(k)
      u2(3) = rz(k)
      
      p1(0) = sqrt(px(j)**2+py(j)**2+pz(j)**2+fmass(j)**2)
      p1(1) = px(j)
      p1(2) = py(j)
      p1(3) = pz(j)
      p2(0) = sqrt(px(k)**2+py(k)**2+pz(k)**2+fmass(k)**2)
      p2(1) = px(k)
      p2(2) = py(k)
      p2(3) = pz(k)
c   
c -velocity and gamma-factor of the two particle center of momentum 
c frame (com) measured in the computational system
c
      bt_com2 = 0.0d0
      do 1 i=1,3 
        bt_com(i) = -((p1(i)+p2(i))/(p1(0)+p2(0)))
        bt_com2 = bt_com2 + bt_com(i)**2
    1 continue
      gam_com = 1.0d0/sqrt(1.0d0-bt_com2)
      gam_com2 = gam_com**2/(1.0d0+gam_com)
c
c calculate some numbers which are needed for the Lorentz-transformation
c
      bt_du = 0.0d0
      bt_dp = 0.0d0
      do 2 i=1,3 
        bt_du = bt_du + bt_com(i)*(u1(i) - u2(i))
        bt_dp = bt_dp + bt_com(i)*(p1(i) - p2(i))
    2 continue
c
c calculate bt_com square and the dotproduct bt_com*(r(j)-r(k)),
c where the r's are given in the computational frame
c 
c perform Lorentz-transformation of relative distance and relative
c momentum vectors of particles j and k into com-frame 
c
c use the resulting 3-vectors du and dp to obtain du and dp squared
c and the dotproduct du*dp
c
      du2  = 0.0d0
      dp2  = 0.0d0
      dudp = 0.0d0
      bt2  = 0.0d0
      btdr = 0.0d0
      do 3 i=1,3 
        bt(i) = p1(i)/p1(0) - p2(i)/p2(0)
        bt2 = bt2 + bt(i)**2
        btdr = btdr + bt(i)*(u1(i)-u2(i))
        du(i) = u1(i)-u2(i) + bt_com(i)*
     *      (gam_com2*bt_du+gam_com*(u1(0)-u2(0)))
        dp(i) = p1(i)-p2(i) + bt_com(i)*
     *      (gam_com2*bt_dp+gam_com*(p1(0)-p2(0)))
        du2  = du2  + du(i)*du(i)
        dp2  = dp2  + dp(i)*dp(i)
        dudp = dudp + du(i)*dp(i)
    3 continue
c
c  obtain collision time and impact parameter squared
c      
      dtauc = -(btdr/bt2)
      dst   = du2 - dudp*dudp/dp2

      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine scantab(ind,offs)
c
c     Revision : 1.1
c
c     [Changes: Removed wrong counter update for nsav, introduced update of
c      indices from previous call of {\tt scantab}]
c
cinput ind : index of particle
cinput offs: offset
c
c     {\tt scantab} adjusts the collision table to changed particle
c     indices due to calls to {\tt addpart} or {\tt delpart}.
c
c     In case of an annihilation a list is created of those particles
c     which have lost their collision partner and have to be rechecked
c     for possible collisions/decays.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      include 'colltab.f'
      include 'coms.f'
      integer ind,offs,i,k,x,z,control
      logical rescan

c     nsav = counter for the rechecking of particles due to lost collision
c     partners. If nsav is still > 0 from previous call of scantab, the relevant
c     indices are adjusted or deleted.

cccccccccc Bugfix January 2012, Wrong table update ccccccccccc

      if (nsav.gt.0) then
        do 711 z=1,nsav
 713       continue
           if ((ctsav(z).eq.ind).and.(offs.lt.0)) then
              if (z.le.nsav) then
                do 712 x=z+1,nsav
                   ctsav(x-1) = ctsav(x)
 712            continue
                nsav = nsav-1
                goto 713
              endif
           endif
           if ((ctsav(z).gt.ind).or.
     &        ((ctsav(z).eq.ind).and.(offs.gt.0))) then
                 ctsav(z) = ctsav(z) + offs
           endif
 711    continue
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      rescan=.false.

c     call from addpart
      if (offs.gt.0) then
c     shift upwards, if necessary
         do 1 i=1,nct
            if (cti1(i).ge.ind) cti1(i) = cti1(i) + offs
            if (cti2(i).ge.ind) cti2(i) = cti2(i) + offs
 1       continue
c
c     call from delpart
      else

c     omit scan, if last collision in table
         if(actcol.eq.nct) return
c     start loop with next collision in table
         i=actcol+1
 2       continue
         if((cti1(i).eq.ind).or.(cti2(i).eq.ind)) then
c     a dubious entry has been found, now
c     look for particles with 'lost' collision partners
            if(cti1(i).eq.ind.and.cti2(i).gt.0) then
c     save particles with 'lost' partners in the ctsav array
               nsav=nsav+1
               ctsav(nsav)=cti2(i)
               if(ctsav(nsav).gt.ind) ctsav(nsav)=ctsav(nsav) + offs
            elseif(cti2(i).eq.ind.and.cti1(i).gt.0) then
               nsav=nsav+1
               ctsav(nsav)=cti1(i)
               if(ctsav(nsav).gt.ind) ctsav(nsav)=ctsav(nsav) + offs
            endif
c     delete obsolete collision
            do 4 k=i+1,nct
               cttime(k-1) = cttime(k)
               ctsqrts(k-1) = ctsqrts(k)
               ctsigtot(k-1) = ctsigtot(k)
               cti1(k-1) = cti1(k)
               cti2(k-1) = cti2(k)
               ctvalid(k-1) = ctvalid(k)
               ctcolfluc(k-1) = ctcolfluc(k)
 4          continue
c     decrement collision counter
            nct = nct-1
            rescan=.true.
         else
c     else entry is OK
            rescan=.false.
         endif
c     shift rest of table, in case entry has not been deleted
         if(.not.rescan) then
            if (cti1(i).gt.ind) cti1(i) = cti1(i) + offs
            if (cti2(i).gt.ind) cti2(i) = cti2(i) + offs
            i=i+1
         endif
c     condinue procedure until all collisions have been scanned
         if(i.le.nct) goto 2
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine printtab
c
c     Revision: 1.0
c
c     {\tt printtab} prints the contents of the collision arrays on
c     unit 6 and marks the current collision with an *. This subroutine
c     is used for debugging purposes only.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        implicit none

        include 'colltab.f'
        integer i
        character*1 c   
c
        write(6,*) 'colltab:'
        do 101 i=1,nct
           c = ' '
           if (i.eq.actcol) c = '*'
           write(6,'(i4,1x,L1,A1,2(i4,1x),4(f6.3,1x))') 
     &            i,ctvalid(i),c,cti1(i),cti2(i),cttime(i),ctsqrts(i)
 101    continue
        return
        end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine colorfluc(it1,it2,ws,fac1,fac2)
c
c     Revision : 1.0
c
cinput it1 : ityp particle 1
cinput it2 : ityp particle 2
cinput ws  : $\sqrt{s}$ of collision
coutput fac1 : x-section scaling factor of particle 1
coutput fac2 : x-section scaling factor of particle 2
c     
c     Modifies the hadron cross section due to color fluctuations,
c     ref. L. Frankfurt et al.: Ann. Rev. Nucl. Part. Sc. 44 (1994)501
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               
      implicit none
      include 'options.f'

      real*8 ws,fac1,fac2
      real*8 x,y,Pmes,Pbar,ranf
      integer it1,it2


      fac1=1d0
      fac2=1d0

c switched on ?
      if (ctoption(42).eq.0) return


c particle 1 is baryon
      if (abs(it1).lt.100) then
 1     x=ranf(0)*3d0
       y=ranf(0)

c probability disrtibution for x-section factor
       Pbar=(1.19+32.11*x-15.65*x**2-1.24*x**3+0.94*x**4)/18.d0
       Pbar=max(0d0,Pbar)
       if(y.gt.Pbar) goto 1
       fac1=x
c particle 1 is meson
      else
 2     x=ranf(0)*5d0
       y=ranf(0)

c probability disrtibution for x-section factor
       Pmes=(21.76+4.41*x-3-79*x**2+0.40*x**3)/25.d0
       Pmes=max(0d0,Pmes)
       if(y.gt.Pmes) goto 2
       fac1=x
      endif

c particle 2 is baryon
      if (abs(it2).lt.100) then
 3     x=ranf(0)*3d0
       y=ranf(0)

c probability disrtibution for x-section factor
       Pbar=(1.19+32.11*x-15.65*x**2-1.24*x**3+0.94*x**4)/18.d0
       Pbar=max(0d0,Pbar)
       if(y.gt.Pbar) goto 3
       fac2=x
c particle 2 is meson
      else
 4     x=ranf(0)*5d0
       y=ranf(0)

c probability disrtibution for x-section factor
       Pmes=(21.76+4.41*x-3-79*x**2+0.40*x**3)/25.d0
       Pmes=max(0d0,Pmes)
       if(y.gt.Pmes) goto 4
       fac2=x
      endif
       
      return
      end

