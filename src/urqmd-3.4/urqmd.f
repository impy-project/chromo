c $Id: urqmd.f 5115 2016-01-04 19:07:31Z darko $
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
cdh   program UrQMD
      subroutine UrQMD(iflbmax)
c
c      Authors : The UrQMD collaboration 
c                S.A. Bass, M. Belkacem, M. Bleicher, M. Brandstetter,
c                L. Bravina, C. Ernst, L. Gerland, M. Hofmann, 
c                S. Hofmann, J. Konopka, G. Mao, L. Neise, S. Soff,
c                C. Spieles, H. Weber, L.A. Winckelmann, H. Stoecker
c                and W. Greiner
c
c     Revision: 1.2
c
cc    contact address:
cc
cc                     urqmd@th.physik.uni-frankfurt.de
cc
c               
c This is the main module of {\tt urqmd}. It servers as a connection between
c the initialization, the propagation (including the real part of the 
c optical potential) and the collision term (imaginary part of the optical
c potential).
c
c  iflbmax: flag for retrying interaction after non-interaction 
c     0 = do retry until interaction happens
c     1 = do not ( propagate particle, retry then) 
c
C  modifications for use in connection with CORSIKA by D. Heck
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      include 'coms.f'
      include 'comres.f'
      include 'options.f'
      include 'colltab.f'
      include 'inputs.f'
      include 'newpart.f'
      include 'boxinc.f'
c
      integer i,j,k,steps,ii,ocharge,ncharge, mc, mp, noc, it1,it2
      real*8 sqrts,otime,xdummy,st
      logical isstable
      integer stidx,iflbmax
      real*8 Ekinbar, Ekinmes, ESky2, ESky3, EYuk, ECb, EPau
      common /energies/ Ekinbar, Ekinmes, ESky2, ESky3, EYuk, ECb, EPau
      integer cti1sav,cti2sav
      common/cncc/ncc
      integer*8 ncc

c
c     numerical/technical initialisation
c
cdh   call uinit(0)
cdh   call osc_header
cdh   call osc99_header
c
c  Main program
c
      noc=0

 1    mc=0
      mp=0
c
c loop over all events
c
cdh   do 10 event=1,nevents

c     start event here
c    

c     time is the system time at the BEGINNING of every timestep
      time = 0.0

c     initialize random number generator
c     call auto-seed generator only for first event and if no seed was fixed
cdh   if(.not.firstseed.and.(.not.fixedseed)) then
cdh      ranseed=-(1*abs(ranseed))
cdh      call sseed(ranseed)
cdh   else
cdh      firstseed=.false.
cdh   endif
cdh
cdh   write(6,*)'event# ',event,ranseed

c
c     initialisation of physics quantities 
c
      call init

c old time if an old fort.14 is used 
      if(CTOption(40).eq.1)time=acttime

c output preparation

c write headers to file
      call output(13)
      call output(14)
      call output(15)
cdh   call output(16)
cdh   if(event.eq.1)call output(17)
cdh   call osc99_event(-1)


c for CTOption(4)=1 : output of initialization configuration
      if(CTOption(4).eq.1)call file14out(0)

c     participant/spectator model:
cdh   if(CTOption(28).ne.0) call rmspec(0.5d0*bimp,-(0.5d0*bimp))

c     compute time of output
      otime = outsteps*dtimestep

c reset time step counter
      steps = 0

c  loop over all timesteps

      do 20  steps=1,nsteps
c store coordinates in arrays with *_t
c this is needed for MD type propagation
         if (eos.ne.0) then
            do 23 j=1,npart
               r0_t(j) = r0(j)
               rx_t(j) = rx(j)
               ry_t(j) = ry(j)
               rz_t(j) = rz(j)
 23         continue
         end if

c we are at the beginning of the timestep, set current time (acttime) 
         acttime = time

c  option for MD without collision term
         if(CTOption(16).ne.0) goto 103

c  Load collision table with next collisions in current timestep
         call colload
c     check for collisions in time-step, nct = # of collisions in table
         if (nct.gt.0) then

c     entry-point for collision loop in case of full colload after every coll. 
 101        continue
            k = 0
c     normal entry-point for collision loop 
 100        continue
c     get next collision
            call getnext(k)

c     exit collision loop if no collisions are left
            if (k.eq.0) goto 102

c  propagate all particles to next collision time
c  store actual time in acttime, propagation time st=cttime(k)-acttime
	    st=cttime(k)-acttime
            call cascstep(acttime,st)
c  new actual time (for upcoming collision)
            acttime = cttime(k)

c  perform collision 

            if(cti2(k).gt.0.)then
             if(abs(sqrts(cti1(k),cti2(k))-ctsqrts(k)).gt.1d-3)then
               write(6,*)' ***(E) wrong collision update (col) ***'
               write(6,*)cti1(k),cti2(k),
     &              ctsqrts(k),sqrts(cti1(k),cti2(k))
             endif
            else if(cti2(k).eq.0.and.
     &              abs(fmass(cti1(k))-ctsqrts(k)).gt.1d-3) then
               write(6,*)' *** main(W) wrong collision update (decay)'
               write(6,*)ctag,cti1(k),ityp(cti1(k)),dectime(cti1(k)),
     &              fmass(cti1(k)),ctsqrts(k)
            endif

            ocharge=charge(cti1(k))
            if(cti2(k).gt.0) ocharge=ocharge+charge(cti2(k))

c     store quantities in local variables for charge conservation check
            it1= ityp(cti1(k))
            if(cti2(k).gt.0)it2= ityp(cti2(k))

c increment "dirty" collision counter
            if(cti2(k).gt.0)then !scatter      !scatter
               mc=mc+1
            endif
c     perform scattering/decay
            ncc     = ncc+1
            cti1sav = cti1(k)                  ! hjd
            cti2sav = cti2(k)                  ! hjd
            call scatter(cti1(k),cti2(k),ctsigtot(k),ctsqrts(k),
     &                   ctcolfluc(k))

c
c  update collision table 
c
c     normal update mode
            if(CTOption(17).eq.0) then
               if(nexit.eq.0) then
c     new collision partners for pauli-blocked states (nexit=0)
c hjd1  and cdh
                  if (cti1(k).ne.cti1sav.or.cti2(k).ne.cti2sav) then
                    goto 1
c                    cti1(k) = cti1sav !!!  hjd1
c                    cti2(k) = cti2sav !!!  hjd1
                  endif
c hjd1  and cdh
                  call collupd(cti1(k),1)
                  if(cti2(k).gt.0) call collupd(cti2(k),1)
               else
                  ncharge=0
c     new collision partners for scattered/produced particles (nexit><0)
                  do 30 i=1,nexit
c     ncharge is used for charge conservation check
                     ncharge=ncharge+charge(inew(i))
                     call collupd(inew(i),1)
 30               continue

c     charge conservation check
                  if(ocharge.ne.ncharge) then
                     write(6,*)'ch-conservation error coll/dec ',ctag
                     write(6,*)'   it1:',it1,'   it2:',it2
                     write(6,*)'   ch:',ocharge,ncharge
                     write(6,*)'cti1(k),cti2(k),ctsigtot(k),ctsqrts(k)'
                     write(6,*)cti1(k),cti2(k),ctsigtot(k),ctsqrts(k)
                  endif
               endif

c     update collisions for partners of annihilated particles
               do 55 ii=1,nsav
                  call collupd(ctsav(ii),1)
 55            continue
               nsav=0

            else               ! (CTOption(17).ne.0)
c     full collision load
               call colload
            endif

            if (CTOption(17).eq.0) goto 100
            goto 101

c     this is the point to jump to after all collisions in the timestep
c     have been taken care of
 102        continue

         endif                  ! (nct.gt.0)

c  After all collisions in the timestep are done, propagate to end of 
c  the timestep.

c     point to jump to in case of MD without collision term
 103     continue

c     increment timestep
         time = time+dtimestep

c  After all collisions in the timestep are done, propagate to end of 
c  the timestep.
         call cascstep(acttime,time-acttime)

c     in case of potential interaction, do MD propagation step
         if (eos.ne.0) then

c set initial conditions for MD propagation-step
            do 24 j=1,npart
               r0(j) = r0_t(j)
               rx(j) = rx_t(j)
               ry(j) = ry_t(j)
               rz(j) = rz_t(j)
 24         continue

c now molecular dynamics trajectories
            call proprk(time,dtimestep)

         end if                 ! (eos.ne.0)

c     perform output if desired
cdh      if(mod(steps,outsteps).eq.0.and.steps.lt.nsteps)then 
cdh         if(CTOption(28).eq.2)call spectrans(otime)
cdh         call file14out(steps)
cdh      endif                  ! output handling

 20   continue                  ! time step loop

c
	acttime=time
c optional decay of all unstable particles before final output
c DANGER: pauli-blocked decays are not performed !!!
         if(CTOption(18).eq.0) then
c no do-loop is used because npart changes in loop-structure
            i=0
            nct=0
            actcol=0
c disable Pauli-Blocker for final decays
            CTOption(10)=1
c decay loop structure starts here
 40         continue
            i=i+1

c is particle unstable
            if(dectime(i).lt.1.d30) then
 41            continue
               isstable = .false.
               do 44 stidx=1,nstable
                  if (ityp(i).eq.stabvec(stidx)) then
cdh                   write (6,*) 'no decay of particle ',ityp(i)
                     isstable = .true.
                  endif
 44            enddo
               if (.not.isstable) then
c     perform decay
                  call scatter(i,0,0.d0,fmass(i),xdummy)
c     backtracing if decay-product is unstable itself
                  if(dectime(i).lt.1.d30) goto 41
               endif
            endif
c     check next particle
            if(i.lt.npart) goto 40
         endif                       ! final decay
c final output

cdh     if(CTOption(28).eq.2)call spectrans(otime)

         call file13out(nsteps)
         call file14out(nsteps)
cdh      call file16out
cdh      call osc_event
cdh      call osc99_event(1)
cdh      call osc99_eoe
      
cdh           print *,"noc",noc,bdist,ctag,bimp,iflbmax
         mp=mp+npart
         if(iflbmax.eq.0.and.ctag.eq.0)then
cdh        write(*,*)'(W) No collision in event ',event
           noc=noc+1
           if(noc.gt.1000)then
             print *,'no collision problem in UrQMD'
             stop
           endif
           goto 1 
         endif

c     end of event loop
cdh 10   continue

cdh   write(6,*)'no. of collisions = ',mc/dble(nevents), ' (per event)'
cdh   write(6,*)'final particles   = ',mp/dble(nevents), ' (per event)'
cdh   write(6,*)'empty events      : ', noc, ' = ', 
cdh  +      noc*1d2/dble(nevents), '%'
      return
      end
